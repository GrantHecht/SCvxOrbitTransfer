
function guessGenEOMs(u,ps,t)
    # Grab parameters
    μ = ps[1]
    p = ps[2]

    # Orbital dynamics
    r       = norm(view(u, 1:3))
    invr3   = 1.0 / r^3
    du      = SVector(p*u[4], p*u[5], p*u[6],
                     -p*μ*u[1]*invr3,
                     -p*μ*u[2]*invr3,
                     -p*μ*u[3]*invr3,
                      p*0.0)
    return du
end

function generateGuess(r0,v0,z0,p,μ,ts)
    prob = ODEProblem(guessGenEOMs,SVector(r0[1],r0[2],r0[3],v0[1],v0[2],v0[3],z0),(ts[1], ts[end]),(μ,p))
    sol  = solve(prob, Vern9(); abstol = 1e-12, reltol = 1e-12, saveat = ts)
    return sol.u
end