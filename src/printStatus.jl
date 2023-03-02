
function printStatus(iter,ν,νbc,J,ΔJ,Δx,Δu,Δp,δ,η,abs,rel)
    if iter == 1
    println("k  | ν        | νbc      | J        | ΔJ       | Δx       | Δu       | Δp       | δ        | η        | abs.     | rel.    ")
    println("---+----------+----------+----------+----------+----------+----------+----------+----------+----------+----------+---------")
    end
    if isnan(ΔJ)
        @printf("%02i | %.2e | %.2e | %.2e |    -     | %.2e | %.2e | %.2e | %.2e | %.2e | %.2e | %.2e\n",
            iter, ν, νbc, J, Δx, Δu, Δp, δ, η, abs, rel)
    else
        @printf("%02i | %.2e | %.2e | %.2e | %.2e | %.2e | %.2e | %.2e | %.2e | %.2e | %.2e | %.2e\n",
            iter, ν, νbc, J, ΔJ, Δx, Δu, Δp, δ, η, abs, rel)
    end
end