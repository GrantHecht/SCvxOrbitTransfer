
struct ScaleParams 
    Sx::Vector{Float64}
    Su::Vector{Float64}
    Sp::Vector{Float64}
    cx::Vector{Float64}
    cu::Vector{Float64}
    cp::Vector{Float64}
end

scaleState(ux, p::ScaleParams)      = (ux - p.cx) ./ p.Sx
scaleControl(uc, p::ScaleParams)    = (uc - p.cu) ./ p.Su
scaleParam(up, p::ScaleParams)      = (up - p.cp[1]) / p.Sp[1]
scaleParam(up::AbstractArray, p::ScaleParams) = (up - p.cp) ./ p.Sp

unscaleState(sx, p::ScaleParams)    = p.Sx.*sx + p.cx
unscaleControl(su, p::ScaleParams)  = p.Su.*su + p.cu
unscaleParam(sp, p::ScaleParams)    = p.Sp[1]*sp + p.cp[1]
unscaleParam(sp::AbstractArray, p::ScaleParams) = p.Sp.*sp + p.cp