##############################################################################
##
## Solve fokker-planck equation
##
##############################################################################
function solve_fp!(fd::FDExp, sol::SolutionExp; maxit::Int64 = 1500, tolFK::Float64 = 1e-7)

    #== sizes ==#
    b, a, z = fd.b, fd.a, fd.z
    bn  = length(b)
    an  = length(a)
    zn  = length(fd.z)
    abn = length(b)*length(a)

    ΔTab = fd.ΔTab
    ΔTagrid = fd.ΔTagrid
    ΔTbgrid = fd.ΔTbgrid

    #== Put density in 3dim matrix ==#
    # gvec_old = vcat(sol.g...)
    # gvec_out = copy(gvec_old)
    #
    # ΔTab = zeros(gvec_old)
    # bvec = zeros(gvec_old)
    # avec = zeros(gvec_old)
    # ijk = zero(Int)
    #
    # for zi = 1:zn, abi = 1:abn
    #     ijk += 1
    #     ΔTab[ijk] = ΔTagrid[afromab(abi)] * ΔTbgrid[bfromab(abi)]
    #     bvec[ijk] = b[bfromab(abi)]
    #     avec[ijk] = a[afromab(abi)]
    # end

    #== build LHS ==#
    fillB!(fd, sol)
    B = fd.B

    #== build RHS ==#
    λoff  = fd.λoff
    invΔᴷ = fd.invΔᴷ
    C = invΔᴷ*eye(zn) + λoff

    distance = 1.0
    gout = hcat(sol.g...)
    gold = copy(gout)

    it = 1
    while distance>tolFK && it<=maxit
        g̃ = gold * C
        updateg!(gout, g̃, B)

        ## Aggregating household decisions ##
        # copy!(gvec_old, gold)
        # copy!(gvec_out, gout)

        # Eaold = sum(avec .* gvec_out .* ΔTab)
        # Eaout = sum(avec .* gvec_out .* ΔTab)
        # Eb = sum(bvec .* gvec .* ΔTab)

        distance = chebyshev(vec(gout), vec(gold))
        # distance = abs(Eaold - Eaout)

        if distance < tolFK
            println("kfe solved : $(it) iterations")
            # println(out_text,"kfe solved : $(it) iterations")
            break
        else
            (it % 50 == 0) && @printf("  density iteration %d, distance %.4f \n", it, distance)
            (it % 50 == 0) && @printf("   sum %.4f --> %.4f\n", sum(gout'*ΔTab), sum(gold'*ΔTab) )
            # (it % 50 == 0) && @printf(out_text,"  density iteration %d, distance %.4f \n", it, distance)
            # (it % 50 == 0) && @printf(out_text,"   sum %.4f --> %.4f\n", sum(gout'*ΔTab), sum(gold'*ΔTab) )
        end

        copy!(gold, gout)
        it += 1
    end

    #== Update the converged density on solution ==#
    for zi=1:zn
        sol.g[zi] = gout[:,zi]
    end

    return Void
end

function updateg!(gout, g̃, B)

    zn = size(gout,2)
    for zi=1:zn
        #== Solve the system B(k) gⁿ⁺¹ (k) = C[] ==#
        gout[:,zi] = B[zi] \ g̃[:,zi]
    end
    return Void
end

function fillB!(fd::FDExp, sol::SolutionExp)

    abn = length(fd.a)*length(fd.b)
    zn  = length(fd.z)
    #= Construct a vector nb*na with the terms Δb̃[bi] * Δã[ai] =#
    ΔTab = fd.ΔTab

    #= D matrix =#
    D = ΔTab
    Dinv = 1./ΔTab
    # =========================================================================

    ### Allocate LHS/RHS matrices ###
    ## ........................................................................

    #== Back out information from fd ==#
    λdiag = fd.λdiag
    invΔᴷ = fd.invΔᴷ

    #== matrices to be filled ==#
    Au = sol.Au          # intensity matrix Au
    A  = sol.A           # use A from HJB as storage
    B  = fd.B            # storage matrix

    for zi = 1:zn
        Avals = nonzeros(A[zi])
        fill!(Avals, zero(Float64))

        Auvals = nonzeros(Au[zi])
        Arows = rowvals(Au[zi])
        for abi in 1:abn
            for k in nzrange(Au[zi], abi)
                row = Arows[k]
                current   = Dinv[row]*Auvals[k]*D[abi]
                Avals[k] = current
            end
        end
    end

    @inbounds for zi in 1:zn

        Avals = nonzeros(A[zi])
        Bvals = nonzeros(B[zi])

        #== erase Bvals ==#
        fill!(Bvals, zero(Float64))

        Brows = rowvals(B[zi])
        for abi in 1:abn
            for k in nzrange(B[zi], abi)
                # loop over elements in the column ij
                row      = Brows[k]
                Bvals[k] = -Avals[k]
                row == abi && (Bvals[k] += invΔᴷ - λdiag[zi])
                #NOTE:110 `λdiag` is entering only in this part of the code
            end
        end
    end
    # =========================================================================

    return Void
end




function solve_fp!(fd::FDImp; maxit::Int64 = 400, tolFK::Float64 = 1e-5)

    #== Back out information ==#
    invΔᴷ = fd.invΔᴷ
    Au = fd.Au          # use Au
    B = fd.B            # storage matrix

    ΔTbgrid = fd.ΔTbgrid
    ΔTagrid = fd.ΔTagrid

    abn = length(ΔTbgrid) * length(ΔTagrid)
    zn  = length(fd.z)

    #== Construct a vector nb*na with the terms Δb̃[bi] * Δã[ai] ==#
    ΔTab = zeros(fd.g⁰)
    ijk = zero(Int)
    for zi = 1:zn

        for abi = 1:abn
            ijk += 1
            ΔTab[ijk] = ΔTagrid[afromab(abi)] * ΔTbgrid[bfromab(abi)]
        end
    end
    D = ΔTab
    Dinv = 1./ΔTab

    Avals = nonzeros(Au)
    Arows = rowvals(Au)
    ijk = zero(Int)
    for zi in 1:zn, abi in 1:abn
        ijk += 1

        for k in nzrange(Au, ijk)
            row = Arows[k]
            current  = Dinv[row]*Avals[k]*D[ijk]
            Avals[k] = current
        end
    end

    Bvals = nonzeros(B); fill!(Bvals, zero(Float64))
    Brows = rowvals(B)
    ijk = zero(Int)
    @inbounds for zi in 1:zn, abi in 1:abn
        ijk += 1
        for k in nzrange(B, ijk)
            # loop over elements in the column ij
            row = Brows[k]
            Bvals[k] = - Avals[k]
            row == ijk && (Bvals[k] += invΔᴷ)
        end
    end

    gold = copy(fd.g⁰)
    gout = zeros(gold)
    it = 1
    distance = 1
    while distance>tolFK && it<=maxit
        g̃ = gold * invΔᴷ
        gout = B \ g̃

        distance = chebyshev(vec(gout), vec(gold))

        if it % 5 == 0
            @printf("Density iteration %d, distance %.4f \n", it, distance)

            #== Check if the summation is held constant across iterations ==#
            @printf("sum gold %.4f \n", sum(gold'*ΔTab))
            @printf("sum gout %.4f \n", sum(gout'*ΔTab))
        end
        copy!(gold, gout)
        it += 1
    end

    fd.g = gout

    return Void
end

# function solve_fp!(fd::TwoAssetsFD)
#
#     g⁰ = fd.g⁰
#     Δa = fd.Δa; Δb = fd.Δb
#
#     # == A' in the notes == #
#     A = fd.A
#
#     i_fix = 1
#     Aurows = rowvals(A)
#     Auvals = nonzeros(A)
#     @inbounds for ij in 1:size(A, 2)
#         for k in nzrange(A, ij)
#             Arows[k] == i_fix && ( Avals[k] = zero(Float64) )
#         end
#     end
#     A[i_fix,i_fix] = one(Float64)
#
#     # == solve system Ag = b == #
#     gg = A \ g⁰
#
#     total = 1/( Δb*Δa*sum(gg) )
#     fd.gg = total * gg
#
#     return Void
# end
