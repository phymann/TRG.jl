"""
    initialize T tensor
"""
function setupT(β, J, h, lattType::String="original")
    if lattType == "original"
        d = 2
        i = Index(d, "up, site")
        j = Index(d, "left, site")
        k = Index(d, "down, site")
        l = Index(d, "right, site")

        Qmat = [
            exp(-β * (-J-2h)) exp(-β * J)
            exp(-β * J) exp(-β * (-J + 2h))
        ]

        # QW = WD
        data = eigen(Qmat)
        sqrtDmat = Diagonal(sqrt.(data.values))
        Wmat = data.vectors
        sqrtQmat = Wmat * sqrtDmat
        #=
                    i'
                    |
                    [X]
                    [X†]
                    |
                    i
                    |
       j'-[X X†]-j-[T]-l-[X X†]-l'
                    |
                    k
                    |
                    [X]
                    [X†]
                    |
                    k'
        =#
        Qi = ITensor(sqrtQmat', i', i)
        Qj = ITensor(sqrtQmat', j', j)
        Qk = ITensor(sqrtQmat, k, k')
        Ql = ITensor(sqrtQmat, l, l')
        delta4 = delta(i, j, k, l)
        T = Qi * Qj * Qk * Ql * delta4
        # replaceinds!(T, [i', j', k', l'], [i, j, k, l])
        return T
    elseif lattType == "dual"
        error("to be implemented")
    else
        error("lattType should be \"original\" or \"dual\"")
    end
end

"""
    update T tensor
"""
function updateT(T; kwargs...)
    # maxdim = get(kwargs, :maxdim, 16)
    # cutoff = get(kwargs, :cutoff, 1e-8)
    # @show maxdim

    l = commonind(T, T, tags="right, site")
    j = commonind(T, T, tags="left, site")
    i = commonind(T, T, tags="up, site")
    k = commonind(T, T, tags="down, site")

    #=
            i
            |
            |
       j----T----l
            |
            |
            k

            i
            |
            ld-l
           /
       j-ru
          |
          k

          i
          |
       j-rd
           \
            lu-l
            |
            k

          i          l
           \   u    /
            lu - ru
           l |   |  r
            ld - rd
           /   d   \
          j         k
    =#

    # ru, ld
    U, S, V = svd(T, j, k; kwargs...)
    # up index = j
    # right index = k
    lnew = commonind(U, S)
    # left index = i
    # down index = l
    jnew = commonind(V, S)

    # remenber to absorb √S into U and V!!!
    sqrtS = ITensor(sqrt.(array(S, lnew, jnew)), lnew, jnew)
    ru = U * sqrtS
    replaceind!(ru, jnew, lnew)
    ld = V * sqrtS
    replaceind!(ld, lnew, jnew)

    # rd, lu
    U, S, V = svd(T, i, j; kwargs...)
    knew = commonind(U, S)
    inew = commonind(V, S)

    sqrtS = ITensor(sqrt.(array(S, knew, inew)), knew, inew)
    rd = U * sqrtS
    replaceind!(rd, inew, knew)
    replaceinds!(rd, [i, j], [k, l])
    lu = V * sqrtS
    replaceind!(lu, knew, inew)
    replaceinds!(lu, [l, k], [j, i])

    Tnew = lu * ru * ld * rd
    nrm = abs(scalar(Tnew * δ(inew, knew) * δ(jnew, lnew)))

    settags!(Tnew, "right, site", lnew)
    settags!(Tnew, "left, site", jnew)
    settags!(Tnew, "down, site", knew)
    settags!(Tnew, "up, site", inew)

    Tnew /= nrm
    return Tnew, nrm
end

function mainTRG(β::Real, J::Real, h::Real; kwargs...)
    convergenceQ = true
    nstep = get(kwargs, :nstep, 32)
    stop_ratio = get(kwargs, :stop_ratio, 1e-18)

    T = setupT(β, J, h)

    fe = 0
    for stepi = 1:nstep
        T, nrm = updateT(T; kwargs...)
        @show size(T)
        #=
            say, after step 1,
                number of local tensor becomes N/2
                each of N/2 tensors are divided by a factor nrm
                i.e., the full partital function is divided by nrm^(N/2)
        =#
        if 2.0^(-stepi) * log(nrm)/fe < stop_ratio
            println("for β = $β, converged at $stepi/$nstep")
            break
        end
        fe -= 2.0^(-stepi) * log(nrm)
        if stepi == nstep
            convergenceQ = false
        end
    end

    fe /= β

    return fe, convergenceQ
end
