# Struct of SCvx Parameters
mutable struct SCvxParams{PenaltyFunction <: Function}
    # Density of temporal discretization
    N::Int

    # Trust region norm
    q::NormType

    # Stopping criterion norm
    qh::NormType

    # Positive definite penalty function
    P::PenaltyFunction

    # Penalty weight
    λ::Float64

    # Initial trust region radius
    η::Float64

    # Minimum and maximum trust region radius
    η0::Float64
    η1::Float64

    # Trust region update parameters
    ρ0::Float64
    ρ1::Float64
    ρ2::Float64

    # Trust region shrinking and growth factors
    βsh::Float64
    βgr::Float64

    # Absolute trajectory change tolerance
    ϵ::Float64

    # Relative cost improvement tolerance
    ϵr::Float64
end

# Constructor function
function SCvxParams(N::Int; 
        q       = two,
        qh      = two,
        P       = (a,b) -> norm(a,1) + norm(b,1),
        λ       = 30.0,
        η       = 1.0,
        η0      = 1e-3,
        η1      = 10.0,
        ρ0      = 0.0,
        ρ1      = 0.1,
        ρ2      = 0.7,
        βsh     = 2.0,
        βgr     = 2.0,
        ϵ       = 0.0,
        ϵr      = 0.0)
    SCvxParams(N,q,qh,P,λ,η,η0,η1,ρ0,ρ1,ρ2,βsh,βgr,ϵ,ϵr)
end