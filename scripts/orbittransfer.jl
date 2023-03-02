using DrWatson
@quickactivate "SCvxQuadrotor"
using PyPlot

# Here you may include files from the source directory
include(srcdir("includeSource.jl"))

function main()
    # Set number of descrete points in temporal discretization
    N   = 120

    # Define SCvx parameters
    ps  = SCvxParams(N; ϵ = 0.0, ϵr = 1e-4, λ = 1e3)

    # Define problem parameters, constants, and units
    m0      = 22.6              # [kg]
    tMax    = 2.2519e-3         # [N]
    Isp     = 3067.0            # [s]
    tf      = 700.0*86400.0     # [s]
    μ       = 1.3271244e11      # [km^3/s^2]
    g0      = 9.80665e-3        # [km/s^2]
    LU      = 1.495978707e8
    VU      = sqrt(μ / LU)
    TU      = LU / VU
    MU      = m0

    # Scale problem parameters to problem units
    m0s     = m0 / MU
    tMaxs   = tMax * TU^2 / (1000.0 * LU * MU)
    Isps    = Isp / TU
    tfs     = tf / TU 
    μs      = μ * TU^2 / LU^3
    g0s     = g0 * TU^2 / LU

    # Define optimal control problem parameters
    r0      = [-0.70186065, 0.70623244, -3.5115e-5] 
    v0      = [-0.73296949, -0.71590485, 4.40245e-5]
    z0      = log(m0s)
    rf      = [0.41806795, 0.82897114, -0.00143382]
    vf      = [-0.96990332, 0.43630220, -0.00123381]

    # Variable scaling
    Sx      = [3.0, 4.0, 0.004, 3.0, 4.0, 0.004, 0.5]
    Su      = [0.0004, 0.0004, 0.0001, 0.00075]
    Sp      = [13.0]
    cx      = [-1.5, -2.0, -0.002, -1.5, -2.0, -0.002, -0.25]
    cu      = [-0.0002, -0.0002, -5e-5, -0.0001]
    cp      = [0.0]
    scaling = ScaleParams(Sx,Su,Sp,cx,cu,cp)

    # Allocate storage for reference trajectory and update
    xk      = zeros(7,N)
    uk      = zeros(4,N)
    pk      = tfs
    νk      = zeros(7,N - 1)
    Jk      = 0.0
    xkp1    = zeros(7,N)
    ukp1    = zeros(4,N)
    pkp1    = 0.0
    νkp1    = zeros(7,N - 1)

    # Compute initial scaled guess
    ts      = range(start = 0.0, stop = 1.0, length = N)
    xs      = generateGuess(r0, v0, z0, tfs, μs, ts)
    for i in eachindex(ts)
        xk[:,i]     .= scaleState(xs[i], scaling)
    end
    for i in axes(uk, 2)
        uk[:,i]     = scaleControl(uk[:,i], scaling)
    end
    pk      = scaleParam(pk, scaling)

    # Allocate storage for problem matricies
    H0      = diagm(ones(7))
    l0      = vcat(-r0, -v0, -z0)
    Hf      = zeros(6,7)
    for i in 1:6; Hf[i,i] = 1.0; end
    lf      = vcat(-rf, -vf)
    δk      = zeros(7,N - 1)
    δkp1    = zeros(7,N - 1)

    # Compute defects for initial guess
    for i = 1:N - 1
        # Compute times 
        dt      = 1.0 / (N - 1)
        tk      = dt*(i - 1)
        tkp1    = dt*i

        # Compute flow map
        fm       = flowMap(
                    unscaleState(xk[:,i], scaling),
                    unscaleControl(uk[:,i], scaling),
                    unscaleControl(uk[:,i + 1], scaling),
                    unscaleParam(pk, scaling),
                    tk,
                    tkp1,
                    μs,
                    g0s,
                    Isps)

        # Compute defect
        δk[:,i] .= unscaleState(xk[:,i + 1], scaling)  .- fm
    end

    # Initialize convex problem variables
    x   = Variable(7, N)
    u   = Variable(4, N)
    p   = Variable()
    ν   = Variable(7, N - 1)
    νic = Variable(7)
    νfc = Variable(6)

    # Begin trajectory optimization loop
    ηk   = ps.η
    done = false
    iter = 0
    while !done
        # ===== Increment iteration counter
        iter   += 1

        # ===== Begin forming cost expression
        cost    = -x[7,end] + ps.λ * ps.P(0.0,νic) + ps.λ * ps.P(0.0,νfc)

        # ===== Add initial time state constraint
        cons    = Constraint[H0*unscaleState(x[:,1], scaling) + l0 + νic == zeros(7)]

        # ===== Add final time state constraint
        push!(cons, Hf*unscaleState(x[:,end], scaling) + lf + νfc == zeros(6))

        # ===== Convex parameter constraint
        push!(cons, unscaleParam(p, scaling) == tfs)

        # ===== Add constraints at each discrete point
        for i = 1:N
            # Control convex path constraints
            dt = 1.0 / (N - 1)
            t  = dt*(i - 1)
            z0 = log(m0s - tMaxs*unscaleParam(pk, scaling)*t/(g0s*Isps)) # If switching to variable tf, will need to update

            # Constraint on Γ upper bound
            ρ = tMaxs*exp(-z0)
            push!(cons, Su[4]*u[4,i] + cu[4] <= ρ*(1.0 - Sx[7]*x[7,i] + cx[7] + z0))

            # Constraints on z
            push!(cons, Sx[7]*x[7,i] + cx[7] <= log(m0))
            push!(cons, Sx[7]*x[7,i] + cx[7] >= z0)

            # Constraint on thrust magintude
            push!(cons, norm(Su[1:3].*u[1:3,i] + cu[1:3], 2) <= Su[4]*u[4,i] + cu[4])

            # Add trust region constraint
            if ps.q == one
                push!(cons, norm(x[:,i] - xk[:,i],1) + norm(u[:,i] - uk[:,i],1) + norm(p - pk,1) <= ηk)
            elseif ps.q == two
                push!(cons, norm(x[:,i] - xk[:,i],2) + norm(u[:,i] - uk[:,i],2) + norm(p - pk,2) <= ηk)
            elseif ps.q == twop
                push!(cons, sumsquares(x[:,i] - xk[:,i]) + sumsquares(u[:,i] - uk[:,i]) + sumsquares(p - pk) <= ηk)
            else
                push!(cons, maximum(x[:,i] - xk[:,i]) + maximum(u[:,i] - uk[:,i]) + maximum(p - pk) <= ηk)
            end

            # Dynamics constraints only added for i > 1
            if i != 1
                # Compute times 
                tk      = dt*(i - 2)
                tkp1    = dt*(i - 1)

                # Compute linearized continuous time matricies
                Ak, Bkm, Bkp, Fk, rk, Ek = computeDiscreteMatricies(
                                                unscaleState(xk[:,i - 1], scaling),
                                                unscaleControl(uk[:,i - 1], scaling),
                                                unscaleControl(uk[:,i], scaling),
                                                unscaleParam(pk, scaling),
                                                tk,
                                                tkp1,
                                                μs,
                                                g0s,
                                                Isps)

                # Add dynamics constraint
                push!(cons, Ak*unscaleState(x[:,i - 1], scaling) + Bkm*unscaleControl(u[:,i - 1], scaling) + 
                            Bkp*unscaleControl(u[:,i], scaling) + Fk*unscaleParam(p, scaling) + 
                            rk + Ek*ν[:,i - 1] == unscaleState(x[:,i], scaling))

                # Update cost expression
                Γ     = ps.λ * ps.P(Ek*ν[:,i-1], 0.0)
                if i == 2
                    cost += dt*Γ / 2.0
                else
                    cost += dt*Γ
                end
            end 
        end

        # Form convex problem and solve
        problem = minimize(cost, cons)
        Convex.solve!(problem, SCS.Optimizer; silent_solver = true)

        # Grab cost function value
        Lλ    = evaluate(cost)

        # Grab new solution
        xkp1 .= evaluate(x)
        ukp1 .= evaluate(u)
        pkp1  = evaluate(p)
        νkp1 .= evaluate(ν)

        # Compute defects
        for i = 1:N - 1
            # Compute times 
            dt      = 1.0 / (N - 1)
            tk      = dt*(i - 1)
            tkp1    = dt*i

            # Compute flow map
            fm       = flowMap(
                        unscaleState(xkp1[:,i], scaling),
                        unscaleControl(ukp1[:,i], scaling),
                        unscaleControl(ukp1[:,i + 1], scaling),
                        unscaleParam(pkp1, scaling),
                        tk,
                        tkp1,
                        μs,
                        g0s,
                        Isps)

            # Compute defect
            δkp1[:,i] .= unscaleState(xkp1[:,i + 1], scaling) .- fm
        end

        # Compute cost difference
        if Jk == 0.0
            ΔJ = NaN
        else
            ΔJ = Lλ - Jk
        end

        # Compute defect cost 
        Jλk     = defectCost(xk,δk,H0,l0,Hf,lf,scaling,N,ps)
        Jλkp1   = defectCost(xkp1,δkp1,H0,l0,Hf,lf,scaling,N,ps)

        # Compute convergence criterion
        absCriterion = 0.0
        norms        = zeros(N)
        if ps.qh == one
            absCriterion += norm(pkp1 - pk,1) 
            for i in eachindex(norm)
                norms[i] = norm(xkp1[:,i] - xk[:,i],1)
            end
        elseif ps.qh == two
            absCriterion += norm(pkp1 - pk,2)
            for i in eachindex(norms)
                norms[i] = norm(xkp1[:,i] - xk[:,i],2)
            end
        elseif ps.qh == twop 
            absCriterion += sumsquares(pkp1 - pk)
            for i in eachindex(norms)
                norms[i] = sumsquares(xkp1[:,i] - xk[:,i])
            end
        else
            absCriterion += maximum(pkp1 - pk)
            for i in eachindex(norms)
                norms[i] = maximum(xkp1[:,i] - xk[:,i])
            end
        end
        absCriterion += maximum(norms) 
        relCriterion  = (Jλk - Lλ) / abs(Jλk)

        # Print status
        printStatus(iter, maximum(νkp1), max(maximum(evaluate(νic)),maximum(evaluate(νfc))),
            Lλ, abs(ΔJ), maximum(abs.(xkp1 - xk)), maximum(abs.(ukp1 - uk)), maximum(abs.(pkp1 - pk)), 
            maximum(abs.(δkp1)), ηk, absCriterion, relCriterion)

        # Check for convergence
        if absCriterion <= ps.ϵ || relCriterion <= ps.ϵr
            if absCriterion <= ps.ϵ
                println("Absolute criterion satisfied. " * string(absCriterion))
            else
                println("Relative criterion satisfied.")
            end

            # Update solution
            xk .= xkp1
            uk .= ukp1
            pk  = pkp1
            νk .= νkp1

            # Set done flag
            done = true

        else # Did not converge, update trust region
            # Compute convexification acuraccy metric
            ρ       = (Jλk - Jλkp1) / (Jλk - Lλ) 

            # Trust region update rule
            if ρ < ps.ρ0
                ηk = max(ps.η0, ηk / ps.βsh)
            elseif ρ >= ps.ρ0 && ρ < ps.ρ1
                ηk = max(ps.η0, ηk / ps.βsh)
                xk .= xkp1
                uk .= ukp1
                pk  = pkp1
                νk .= νkp1
                δk .= δkp1
                Jk  = Lλ
            elseif ρ >= ps.ρ1 && ρ < ps.ρ1 
                xk .= xkp1
                uk .= ukp1
                pk  = pkp1
                νk .= νkp1
                δk .= δkp1           
                Jk  = Lλ
            else
                ηk = min(ps.η1, ηk * ps.βgr)
                xk .= xkp1
                uk .= ukp1
                pk  = pkp1
                νk .= νkp1
                δk .= δkp1
                Jk  = Lλ
            end
        end
        if iter == 10
            done = true
        end
    end

    # Plot trajectory
    plot3D(Sx[1]*xkp1[1,:] .+ cx[1],Sx[2]*xkp1[2,:] .+ cx[2],Sx[3]*xkp1[3,:] .+ cx[3])
    scatter3D([r0[1],rf[1]],[r0[2],rf[2]],[r0[3],rf[3]])

    display(gcf())

    return xkp1, ukp1, pkp1
end

xkp1, ukp1, pkp1 = main()

println("Max rx: " * string(maximum(xkp1[1,:])))
println("Max ry: " * string(maximum(xkp1[2,:])))
println("Max rz: " * string(maximum(xkp1[3,:])))
println("Max vx: " * string(maximum(xkp1[4,:])))
println("Max vy: " * string(maximum(xkp1[5,:])))
println("Max vz: " * string(maximum(xkp1[6,:])))
println("Max z:  " * string(maximum(xkp1[7,:])))
println("Min rx: " * string(minimum(xkp1[1,:])))
println("Min ry: " * string(minimum(xkp1[2,:])))
println("Min rz: " * string(minimum(xkp1[3,:])))
println("Min vx: " * string(minimum(xkp1[4,:])))
println("Min vy: " * string(minimum(xkp1[5,:])))
println("Min vz: " * string(minimum(xkp1[6,:])))
println("Min z:  " * string(minimum(xkp1[7,:])))
println("Max τx: " * string(maximum(ukp1[1,:])))
println("Max τy: " * string(maximum(ukp1[2,:])))
println("Max τz: " * string(maximum(ukp1[3,:])))
println("Max Γ:  " * string(maximum(ukp1[4,:])))
println("Min τx: " * string(minimum(ukp1[1,:])))
println("Min τy: " * string(minimum(ukp1[2,:])))
println("Min τz: " * string(minimum(ukp1[3,:])))
println("Min Γ:  " * string(minimum(ukp1[4,:])))
println("P: " * string(pkp1))
