
function flowMapEOMs(u,ps,t)
     # Grab parameters
    uk      = ps[1]
    ukp1    = ps[2]
    tk      = ps[3]
    tkp1    = ps[4]
    p       = ps[5]
    μ       = ps[6]   
    g0      = ps[7]
    Isp     = ps[8]

    # Compute control
    σm      = (tkp1 - t) / (tkp1 - tk)
    σp      = (t - tk) / (tkp1 - tk)
    ut      = σm*uk + σp*ukp1

    # Compute dynamics
    r       = norm(view(u,1:3))
    invr3   = 1.0 / r^2
    du      = SVector(p*u[4], p*u[5], p*u[6],
                     -μ*u[1]*invr3 + ut[1],
                     -μ*u[2]*invr3 + ut[2],
                     -μ*u[3]*invr3 + ut[3],
                     -ut[4] / (g0*Isp))
    return du
end

function flowMap(xk, uk, ukp1, p, tk, tkp1, μ, g0, Isp)
    # Setup parameters and time span
    ps  = (uk,ukp1,tk,tkp1,p,μ,g0,Isp)
    ts  = (tk, tkp1)

    # Perform integration
    prob = ODEProblem(flowMapEOMs, SVector(xk...), ts, ps)
    sol  = solve(prob, Vern9(); reltol = 1e-12, abstol = 1e-12, 
            save_everystep = false, save_start = false)

    return sol[end]
end