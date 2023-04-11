using LinearAlgebra
using ITensors
using QuadGK
using Dates
using MAT

const βc = 0.5 * log(√2 + 1)

@time let
    include("miscellaneous.jl")
    include("mainTRG.jl")

    println("-----------------------------------------")
    println(Dates.now())
    timeStamp = Dates.format(now(), "DyyyymmddTHHMMSS")

    J = 1.0
    h = 0.0

    nβ = 2
    lsβ = range(.9*βc, βc, length=nβ)
    lsfe = zeros(nβ)

    n_notconv = 0
    βls_notconv = []

    for (idx, β) in enumerate(lsβ)
        feTRG, convergenceQ = mainTRG(β, J, h; nstep=42, maxdim=64, cutoff=1e-12, stop_ratio=1e-12)
        lsfe[idx] = feTRG
        feExact = ising_free_energy(β, J)
        @show abs((feTRG-feExact)/feExact)

        if !convergenceQ
            n_notconv += 1
            append!(βls_notconv)
        end
    end

    file = matopen(timeStamp*".mat", "w")
    write(file, "lsfe", lsfe)
    write(file, "lsb", collect(lsβ))
    close(file)

    println("------------------------------")
    println("# of cases without convergence = $n_notconv")
    println("------------------------------")
    @pt βls_notconv
end

return nothing
