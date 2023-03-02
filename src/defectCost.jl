# Functions for grabing only positive part of array or number
positivePart(x::AbstractArray) = positivePart.(x)
positivePart(x::Number) = x >= 0.0 ? x : 0.0

function defectCost(x,δ,H0,l0,Hf,lf,scaling,N,ps)
    # Compute time info
    dt      = 1.0 / (N - 1)

    # Compute terminal state componants
    gic     = H0*unscaleState(x[:,1], scaling) + l0
    gtc     = Hf*unscaleState(x[:,end], scaling) + lf
    cost    = -x[7,end] + ps.λ * ps.P(0.0, gic) + ps.λ * ps.P(0.0, gtc)

    # Compute trapezoidal rule for running cost
    for i = 2:N
        Γ   = ps.λ * ps.P(δ[:,i - 1], 0.0)
        if i == 2
            cost += 0.5*dt*Γ
        else
            cost += dt*Γ
        end
    end

    # Return cost
    return cost
end