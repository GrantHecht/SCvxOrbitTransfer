
# Linearized continuous time equations of motion
# Also compute the state transition matrix and convolution 
# matricies
function linconteoms!(du,u,ps,t)
    # Grab parameters
    uk      = ps[1]
    ukp1    = ps[2]
    tk      = ps[3]
    tkp1    = ps[4]
    p       = ps[5]
    μ       = ps[6]
    g0      = ps[7]
    Isp     = ps[8]

    # Get parts of full state vector and state vector derivative
    x       = view(u,   1:7)
    dx      = view(du,  1:7)
    stm     = reshape(view(u,  8:56),    7, 7)
    dstm    = reshape(view(du, 8:56),    7, 7)
    dBkm    = reshape(view(du, 57:84),   7, 4)
    dBkp    = reshape(view(du, 85:112),  7, 4)
    dFk     = view(du, 113:119)
    drk     = view(du, 120:126)
    dEk     = reshape(view(du, 127:175), 7, 7)

    # Compute control
    σm      = (tkp1 - t) / (tkp1 - tk)
    σp      = (t - tk) / (tkp1 - tk)
    ut      = σm*uk + σp*ukp1

    # Compute continuous time matricies
    A       = zeros(7,7)
    B       = zeros(7,4)

    # Trivial parts of A and B
    for i in 1:3
        A[i, 3 + i] = p
        B[3 + i, i] = p
    end
    B[7, 4] = -p / (g0*Isp)

    # Partial of dv w.r.t. r
    I3      = diagm(ones(3))
    r       = norm(view(x, 1:3))
    invr3   = 1.0 / r^3
    invr5   = 1.0 / r^5
    G       = -μ*invr3*I3 + 3.0*μ*invr5*view(x, 1:3)*transpose(view(x, 1:3))
    for i in 1:3
        for j in 1:3
            A[3 + j,i] = p*G[j,i]
        end
    end

    F       = [x[4], x[5], x[6], 
              -μ*x[1]*invr3 + ut[1],
              -μ*x[2]*invr3 + ut[2],
              -μ*x[3]*invr3 + ut[3],
              -ut[4] / (g0*Isp)]
    r       = -A*x - B*ut
    E       = diagm(ones(7))

    # State dynamics
    dx     .= p*F

    # State transition matrix dynamics
    mul!(dstm, A, stm)

    # Convolution matrix dynamics
    stmInv  = inv(stm)
    mul!(dFk,  stmInv, F)
    mul!(drk,  stmInv, r)
    mul!(dEk,  stmInv, E)
    mul!(dBkm, stmInv, B, σm, 0.0)
    mul!(dBkp, stmInv, B, σp, 0.0)
    return nothing
end

function computeDiscreteMatricies(xk,uk,ukp1,p,tk,tkp1,μ,g0,Isp)
    # Form initial state vector
    u0          = zeros(175)
    u0[1:7]    .= xk
    for i in 1:7
        u0[7 + i + 7*(i - 1)] = 1.0
    end

    # Construct parameter tuple
    ps  = (uk,ukp1,tk,tkp1,p,μ,g0,Isp)

    # Construct time span
    tspan = (tk, tkp1)

    # Perform numerical integration
    prob  = ODEProblem(linconteoms!, u0, tspan, ps)
    sol   = solve(prob, Vern9(); reltol = 1e-12, abstol = 1e-12, save_everystep = false, save_start = false)

    # Grab solutions
    @views Bkmt = reshape(sol[end][57:84],  7,4)
    @views Bkpt = reshape(sol[end][85:112], 7,4)
    @views Fkt  = sol[end][113:119]
    @views rkt  = sol[end][120:126]
    @views Ekt  = reshape(sol[end][127:175], 7,7)
    Ak          = reshape(sol[end][8:56], 7,7)

    # Compute convolution matricies
    Bkm         = Ak*Bkmt
    Bkp         = Ak*Bkpt
    Fk          = Ak*Fkt
    rk          = Ak*rkt
    Ek          = Ak*Ekt
    return (Ak, Bkm, Bkp, Fk, rk, Ek)
end