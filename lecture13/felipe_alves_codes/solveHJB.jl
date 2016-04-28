##############################################################################
##
## Solve hamilton-jacobi-bellman
##
##############################################################################


##############################
##
## FDSpec
##
##############################
"""
Use finite difference method to solve for HJB equation.

### INPUTS
- `fd`    : Finite Difference scheme
- `sol`   : holds solution information
- `twoap` : household problem description
"""
function solve_hjb!(twoap::TwoAssetsProb, fd::FDSpec, sol::SolutionExp;
                    maxit::Int = 400, crit::Float64 = 1e-8, monotonicity = true, verbose::Bool = true)

    bn, an, zn = length(fd.b), length(fd.a), length(fd.z)

    #== Initialize Vold ==#
    Vold = zeros(bn,an,zn)

    #== Construct functions ==#
    compute_pol!(solut, V_old) = comp_pol!(twoap, fd, solut, V_old)
    update_V!(solut, V_old)  = updateV!(twoap, fd, solut, V_old, test_mon = monotonicity)

    for it in 1:maxit
        copy!(Vold, vcat(sol.V...))
        # vcat is important in the FDExp case where fd.V is vector of vectors

        #== Compute the Optimal Policy ==#
        compute_pol!(sol, Vold)

        #== Update vⁿ ==#
        update_V!(sol, Vold)

        # check convergence
        distance = chebyshev(vcat(sol.V...), vec(Vold))
        # distance = norm(vec(fd.V)-vec(Vold), Inf)
        if distance < crit
            println("hjb solved : $(it) iterations")
            # println(out_text,"hjb solved : $(it) iterations")
            break
        else
            # update V using newV
            verbose && (it % 10 == 0) && @printf("  value function iteration %d, distance %.4f \n", it,distance)
            # verbose && (it % 10 == 0) && @printf(out_text,"  value function iteration %d, distance %.4f \n", it,distance)
        end
    end
end

function optconsump(twoap::TwoAssetsProb, ∂Vb::Float64, liqinc::Float64, ℓutil::Float64)

    invγ   = 1.0/twoap.γ

    # ∂Vb>=0 ? c = ∂Vb^(-invγ) : error("negative value not allowed")
    # REVIEW:110 trick to not get stuck in case of non-monotonicity
    ∂Vb>=0 ? c = ∂Vb^(-invγ) + ℓutil : c = liqinc
    sc = liqinc - c
    dV = utilfn(twoap, c, ℓutil) + sc*∂Vb

    return c, sc, dV
end

function comp_pol!(twoap::TwoAssetsProb, fd::FDSpec, sol::SolutionExp, V::Array{Float64,3})

    #== Construct fns ==#
    dχinv1(∂ratio, a) = dχinv(twoap, ∂ratio, a)
    χ1(d,a) = χ(twoap, d, a)
    utilfn1(c, ℓutil) = utilfn(twoap, c, ℓutil)

    #== Check if function is changing ==#
    # println(dχinv1(1.1, 1.0))
    # println(χ1(1e-5, 100.0))

    #== Unpack from FD ==#
    z = fd.z; a = fd.a ; b = fd.b;
    Δbgrid = fd.Δbgrid
    Δagrid = fd.Δagrid
    netbinc = fd.netbinc

    ℓutilgrid = fd.ℓutilgrid

    # precompute the lengths
    bn = length(b); an = length(a); zn = length(z)

    invdb = 1.0./Δbgrid
    invda = 1.0./Δagrid

    for zi in 1:zn, ai in 1:an, bi in 1:bn

        # == Derivative w.r.t B == #
        bi<bn && ( ∂Vbᶠ = ( V[bi+1, ai, zi] - V[bi, ai, zi] ) * invdb[bi] )
        bi>1  && ( ∂Vbᴮ = ( V[bi, ai, zi] - V[bi-1, ai, zi] ) * invdb[bi-1] )

        # == Derivative w.r.t A == #
        ai<an && ( ∂Vaᶠ = ( V[bi, ai+1, zi] - V[bi, ai, zi] ) * invda[ai] )
        ai>1  && ( ∂Vaᴮ = ( V[bi, ai, zi] - V[bi, ai-1, zi] ) * invda[ai-1] )

        ### CONSUMPTION decision ###
        liqinc = netbinc[bi,zi]
        ℓutil  = ℓutilgrid[zi]

        # == Foward == #
        if bi<bn
            cᶠ, scᶠ, Hcᶠ  = optconsump(twoap, ∂Vbᶠ, liqinc, ℓutil)
        else
            scᶠ, dVᶠ = 0.0, -1.0e12
        end
        scᶠ>0.0 ? (validᶠ = 1 ) : (validᶠ = 0 )

        # == Backward == #
        if bi>1
            cᴮ, scᴮ, Hcᴮ  = optconsump(twoap, ∂Vbᴮ, liqinc, ℓutil)
        else
            scᴮ, Hcᴮ = 0.0, -1.0e12
        end
        scᴮ<0.0 ? (validᴮ = 1) : (validᴮ = 0)

        #== Catch Other cases ==#
        c⁰  = liqinc
        sc⁰ = 0.0
        Hc⁰ = utilfn1(c⁰, ℓutil)

        if validᶠ==1   && ( validᴮ==0 || Hcᶠ>=Hcᴮ ) && ( Hcᶠ>=Hc⁰ )
            sol.c[bi, ai, zi]  = cᶠ
            sol.sc[bi, ai, zi] = scᶠ
        elseif validᴮ==1 && ( validᶠ==0 || Hcᴮ>=Hcᶠ ) && ( Hcᴮ>=Hc⁰ )
            sol.c[bi, ai, zi]  = cᴮ
            sol.sc[bi, ai, zi] = scᴮ
        else
            sol.c[bi, ai, zi]  = c⁰
            sol.sc[bi, ai, zi] = sc⁰
        end

        #== DEPOSIT decision ==#
        #== b backward, a forward ==#
        if bi>1 && ai<an
            dᴮᶠ = dχinv1(∂Vaᶠ/∂Vbᴮ, a[ai])
            # change in utility
            Hdᴮᶠ = ∂Vaᶠ * dᴮᶠ - ∂Vbᴮ*( dᴮᶠ + χ1(dᴮᶠ,a[ai]) )

            (dᴮᶠ>0.0 && Hdᴮᶠ>=0.0) ? ( validᴮᶠ=1 ) : ( validᴮᶠ=0 )
        else
            validᴮᶠ = 0; Hdᴮᶠ = -1.0e-12
        end

        #== b forward, a backward ==#
        if bi < bn && ai > 1
            dᶠᴮ = dχinv1(∂Vaᴮ/∂Vbᶠ, a[ai])
            # change in utility
            Hdᶠᴮ = ∂Vaᴮ * dᶠᴮ - ∂Vbᶠ*( dᶠᴮ + χ1(dᶠᴮ,a[ai]) )

            (dᶠᴮ<-χ1(dᶠᴮ,a[ai]) && Hdᶠᴮ>=0.0) ? ( validᶠᴮ=1 ) : ( validᶠᴮ=0 )
        else
            validᶠᴮ = 0; Hdᶠᴮ = -1.0e-12
        end

        #== b backward, a backward ==#
        if bi>1 && ai>1
            dᴮᴮ = dχinv1(∂Vaᴮ/∂Vbᴮ, a[ai])
            # change in utility
            Hdᴮᴮ = ∂Vaᴮ * dᴮᴮ - ∂Vbᴮ*( dᴮᴮ + χ1(dᴮᴮ,a[ai]) )

            (dᴮᴮ>-χ1(dᴮᴮ, a[ai]) && dᴮᴮ<=0 && Hdᴮᴮ>=0.0) ? ( validᴮᴮ=1 ) : ( validᴮᴮ=0 )

            # REVIEW:80 BOUNDARY Adjustment 01
            #        force the use of ∂Vbᴮ in case bi==bn
            # bi == bn && (dᴮᴮ<=0 && Hdᴮᴮ>=0.0) && (validᴮᴮ = 1)
        else
            validᴮᴮ = 0; Hdᴮᴮ = -1.0e-12
        end

        #== Check and assign deposit ==#
        if    validᴮᶠ==1                 && (validᶠᴮ==0 || Hdᴮᶠ>=Hdᶠᴮ)   && (validᴮᴮ==0 || Hdᴮᶠ>=Hdᴮᴮ)

            sol.d[bi, ai, zi] = dᴮᶠ
        # ------------------------------------------------------------------------------------------ #
        elseif (validᴮᶠ==0 || Hdᶠᴮ>=Hdᴮᶠ) &&  validᶠᴮ==1                  && (validᴮᴮ==0 || Hdᶠᴮ>=Hdᴮᴮ)

            sol.d[bi, ai, zi] = dᶠᴮ
        # ------------------------------------------------------------------------------------------ #
        elseif (validᴮᶠ==0 || Hdᴮᴮ>=Hdᴮᶠ) && (validᶠᴮ==0 || Hdᴮᴮ>=Hdᶠᴮ)   &&  validᴮᴮ==1

            sol.d[bi, ai, zi] = dᴮᴮ
        # ------------------------------------------------------------------------------------------ #
        elseif validᴮᶠ==0 && validᶠᴮ==0 && validᴮᴮ==0

            sol.d[bi, ai, zi] = 0.0
        end

        # u[bi, ai, zi] = utilfn1( c[bi, ai, zi] )
        # bdot[bi, ai, zi] = sc[bi, ai, zi] - d[bi, ai, zi] - χ1(d, a[ai])
    end

    Void
end


function updateV!(twoap::TwoAssetsProb, fd::FDExp, sol::SolutionExp, V::Array{Float64,3}; test_mon = true)

    #== Construct function ==#
    χ1(d,a) = χ(twoap, d, a)
    utilfn1(c, ℓutil) = utilfn(twoap, c, ℓutil)
    #== Idiosyncratic component ==#
    λdiag = fd.λdiag
    λtoff = fd.λoff'

    #== HOUSEHOLD Parameters ==#
    γ, ρ, _, σ, ψ = _unpackparams(twoap)

    #== Optimal policies ==#
    c = sol.c; sc = sol.sc; d = sol.d

    #== Extract GRIDS ==#
    b = fd.b; a = fd.a; z = fd.z
    Δbgrid, Δagrid = fd.Δbgrid, fd.Δagrid
    netainc = fd.netainc;
    ℓutilgrid = fd.ℓutilgrid

    invΔ = fd.invΔ

    #== SOLUTION Matrices to fill ==#
    Vnew = sol.V
    A    = sol.A
    Au   = sol.Au

    #== Storage matrices ==#
    B  = fd.B
    b̃  = fd.b̃

    # precompute the lengths
    bn = length(b); an = length(a); zn = length(z); abn = bn*an

    invdb = 1.0./Δbgrid
    invda = 1.0./Δagrid
    inv1γ = 1.0/(1.0-γ)

    for zi in 1:zn
        #== set values to zero ==#
        fill!(nonzeros(A[zi]), zero(Float64))
        fill!(nonzeros(Au[zi]), zero(Float64))

        for abi in 1:abn
        # NOTE:90 Positions of A that I need to fill up
            # A[zi][abi-bn, abi] - v(i  ,j-1,k) --> - invda*[ d⁻ ]                  (SHOULD BE POSITIVE)
            # A[zi][abi-1 , abi] - v(i+1,j  ,k) --> - invdb*[ (scᴮ)⁻ + (sdᴮ)⁻ ]     (SHOULD BE POSITIVE)
            # A[zi][abi   , abi] - v(i  ,j  ,k) --> + invdb*[ (scᴮ)⁻ + (sdᴮ)⁻ - ( (scᶠ)⁺ - (sdᶠ)⁺ ) ]
            #                                       + invda*[ d⁻ - (d⁺ + ξwz + rᵃa) ]
            # A[zi][abi+1 , abi] - v(i-1,j  ,k) --> + invdb*[ (scᶠ)⁺ + (sdᶠ)⁺ ]     (SHOULD BE POSITIVE)
            # A[zi][abi+bn, abi] - v(i  .j+1,k) --> + invda*[ d⁺ + ξwx + rᵃa ]      (SHOULD BE POSITIVE)

            ai = afromab(abi)
            bi = bfromab(abi)
            ℓutil  = ℓutilgrid[zi]

            #== RHS of the Bellman equation==#
            b̃[zi][abi] = utilfn1(c[bi, ai, zi], ℓutil) + V[bi, ai, zi] * invΔ + dot(λtoff[:,zi] ,V[bi,ai,:][:] )

            ## COMPUTE THE ENTRIES ##
            # ==================================================================
            #== bdrift ==#
            bdriftᶠ = max( sc[bi, ai, zi], 0.0 ) + max( -d[bi, ai, zi] -χ1(d[bi, ai, zi], a[ai]) , 0.0 )
            bdriftᴮ = min( sc[bi, ai, zi], 0.0 ) + min( -d[bi, ai, zi] -χ1(d[bi, ai, zi], a[ai]) , 0.0 )

            bi <  bn && ( budriftᶠ = max( sc[bi, ai, zi] - d[bi, ai, zi] - χ1(d[bi, ai, zi], a[ai]) , 0.0 ) )
            bi <  bn && ( budriftᴮ = min( sc[bi, ai, zi] - d[bi, ai, zi] - χ1(d[bi, ai, zi], a[ai]) , 0.0 ) )
            bi == bn && ( budriftᶠ = max( sc[bi, ai, zi] - d[bi-1, ai, zi] - χ1(d[bi-1, ai, zi], a[ai]) , 0.0 ) )
            bi == bn && ( budriftᴮ = min( sc[bi, ai, zi] - d[bi-1, ai, zi] - χ1(d[bi-1, ai, zi], a[ai]) , 0.0 ) )

            #== adrift ==#
            adriftᶠ = max( d[bi, ai, zi], 0.0 ) + netainc[ai,zi]
            adriftᴮ = min( d[bi, ai, zi], 0.0 )

            bi <  bn && ( audriftᶠ = max( d[bi, ai, zi] + netainc[ai,zi], 0.0 ) )
            bi <  bn && ( audriftᴮ = min( d[bi, ai, zi] + netainc[ai,zi], 0.0 ) )
            bi == bn && ( audriftᶠ = max( d[bi-1, ai, zi] + netainc[ai,zi], 0.0 ) )
            bi == bn && ( audriftᴮ = min( d[bi-1, ai, zi] + netainc[ai,zi], 0.0 ) )

            #== BOUNDARY Adjustment 02 ==#
            # bi <  bn && ( bdriftᴮ = min( sc[bi, ai, zi], 0.0 ) + min( -d[bi, ai, zi] -χ1(d[bi, ai, zi], a[ai]) , 0.0 ) )
            # bi == bn && ( bdriftᴮ = min( sc[bi, ai, zi], 0.0 ) -d[bi, ai, zi] -χ1(d[bi, ai, zi], a[ai]) )
            #== BOUNDARY Adjustment 03 ==#
            # ai <  an && ( adriftᴮ = min( d[bi, ai, zi], 0.0 ) )
            # ai == an && ( adriftᴮ = d[bi, ai, zi] + netainc[ai,zi] )

            # v[i, j- 1, k]
            if adriftᴮ != 0.0 && ai>1
                row = abfromab(ai-1,bi)
                A[zi][row, abi] = -adriftᴮ * invda[ai-1]
            end

            # v[i, j- 1, k] KF
            if audriftᴮ != 0.0 && ai>1
                row = abfromab(ai-1,bi)
                Au[zi][row, abi] = -audriftᴮ * invda[ai-1]
            end

            # v[i-1, j, k]
            if bdriftᴮ != 0.0 && bi>1
                row = abfromab(ai,bi-1)
                A[zi][row, abi] = -bdriftᴮ * invdb[bi-1]
            end

            # v[i-1, j, k] KF
            if budriftᴮ != 0.0 && bi>1
                row = abfromab(ai,bi-1)
                Au[zi][row, abi] = -budriftᴮ * invdb[bi-1]
            end

            # v[i, j, k]
            val = 0.0
            ( bdriftᴮ != 0.0 && bi > 1 ) && ( val += bdriftᴮ * invdb[bi-1] )
            ( bdriftᶠ != 0.0 && bi < bn) && ( val -= bdriftᶠ * invdb[bi]   )
            ( adriftᴮ != 0.0 && ai > 1)  && ( val += adriftᴮ * invda[ai-1] )
            ( adriftᶠ != 0.0 && ai < an) && ( val -= adriftᶠ * invda[ai]   )
            val != 0.0 && ( A[zi][abi, abi] = val )

            # v[i, j, k] KF
            val = 0.0
            ( budriftᴮ != 0.0 && bi > 1 ) && ( val += budriftᴮ * invdb[bi-1] )
            ( budriftᶠ != 0.0 && bi < bn) && ( val -= budriftᶠ * invdb[bi]   )
            ( audriftᴮ != 0.0 && ai > 1)  && ( val += audriftᴮ * invda[ai-1] )
            ( audriftᶠ != 0.0 && ai < an) && ( val -= audriftᶠ * invda[ai]   )
            val != 0.0 && ( Au[zi][abi, abi] = val )

            # v[i-1, j, k]
            if bdriftᶠ != 0.0 && bi<bn
                row = abfromab(ai,bi+1)
                A[zi][row, abi] = bdriftᶠ * invdb[bi]
            end

            # v[i-1, j, k] KF
            if budriftᶠ != 0.0 && bi<bn
                row = abfromab(ai,bi+1)
                Au[zi][row, abi] = budriftᶠ * invdb[bi]
            end

            # v[i, j+1, k]
            if adriftᶠ != 0.0 && ai<an
                row = abfromab(ai+1,bi)
                A[zi][row, abi] = adriftᶠ * invda[ai]
            end

            # v[i, j+1, k]
            if audriftᶠ != 0.0 && ai<an
                row = abfromab(ai+1,bi)
                Au[zi][row, abi] = audriftᶠ * invda[ai]
            end

        end
    end

    ## Set B = diag(invΔ + ρ) - A
    # TODO: think of a better way to do this part
    @inbounds for zi in 1:zn
        Avals = nonzeros(A[zi])
        Bvals = nonzeros(B[zi])
        #== clean the sotorage matrix before assigning values ==#
        fill!(Bvals, zero(Float64))
        Brows = rowvals(B[zi])
        for abi in 1:abn
            for k in nzrange(B[zi], abi)
                # loop over elements in the column ij
                row = Brows[k]
                Bvals[k] = - Avals[k]
                row == abi && (Bvals[k] += invΔ + ρ - λdiag[zi])
                #NOTE:120 `λdiag` is entering only in this ste
            end
        end
        Vnew[zi][:] = B[zi]' \ b̃[zi]
    end

    #== Check for monotonicity ==#
    if test_mon
        for zi in 1:zn, abi in 1:abn
            ai = afromab(abi)
            bi = bfromab(abi)
            #= Check mon in the b dimension =#
            if bi>1
                Vnew[zi][abi]<Vnew[zi][abfromab(ai,bi-1)] &&
                @printf("   non-mon on b dimension at = (%d,%d)  \n", bi,ai)
            end

            #= Check mon in the a dimension =#
            if ai>1
                Vnew[zi][abi]<Vnew[zi][abfromab(ai-1,bi)] &&
                @printf("   non-mon on a dimension at = (%d,%d)  \n", bi,ai)
            end
        end
    end
    Void
end

"""
Update value function based on Fortran code
"""
function updateV!(fd::FDImp, sol::SolutionImp, twoap::TwoAssetsProb, V::Array{Float64,3}; test_mon = true)

    #== HOUSEHOLD Parameters ==#
    γ, ρ, ξ, _, _ = _unpackparams(twoap)

    #== Optimal policies ==#
    c = sol.c; sc = sol.sc; d = sol.d

    #== Extract grids and matrices ==#
    b = fd.b; a = fd.a; z = fd.z
    Δbgrid, Δagrid = fd.Δbgrid, fd.Δagrid
    netainc = fd.netainc;
    ℓutilgrid = fd.ℓutilgrid
    invΔ = fd.invΔ

    #== Solution matrices to fill ==#
    V  = sol.V
    A  = sol.A
    Au = sol.Au

    #== Storage matrices to fill ==#
    b̃ = fd.b̃
    Λ = fd.Λ
    B = fd.B

    # precompute the lengths
    bn = length(b); an = length(a); zn = length(z); abn = bn*an

    invdb = 1.0./Δbgrid
    invda = 1.0./Δagrid
    inv1γ = 1.0/(1.0-γ)

    ijk = zero(Int)

    #== set past values to zero ==#
    fill!(nonzeros(Au), zero(Float64))
    fill!(nonzeros(A), zero(Float64))

    for zi in 1:zn, abi in 1:abn
    # NOTE:190 Positions of A that I need to fill
        # A[ijk-bn, ijk] - v(i  ,j-1,k) --> - invda*[ d⁻ ]
        # A[ijk- 1, ijk] - v(i+1,j  ,k) --> - invdb*[ (scᴮ)⁻ + (sdᴮ)⁻ ]
        # A[ijk   , ijk] - v(i  ,j  ,k) -->   invdb*[ (scᴮ)⁻ + (sdᴮ)⁻ -(scᶠ)⁺ -(sdᶠ)⁺] + invda*[ d⁻ -(d⁺ + ξwz + rᵃa) ]
        # A[ijk+ 1, ijk] - v(i-1,j  ,k) -->   invdb*[ (scᶠ)⁺ + (sdᶠ)⁺ ]
        # A[ijk+bn, ijk] - v(i  .j+1,k) -->   invda*[ d⁺ + ξwx + rᵃa ]

        ijk +=1

        ai = afromab(abi)
        bi = bfromab(abi)
        ℓutil  = ℓutilgrid[zi]

        #== RHS of the Bellman equation==#
        b̃[ijk] =  utilfn1(c[bi, ai, zi], ℓutil)  + V[bi, ai, zi] * invΔ
        # ==================================================================

        ### COMPUTE THE ENTRIES ###

        #== bdrift ==#
        bdriftᶠ = max( sc[bi, ai, zi], 0.0 ) + max( -d[bi, ai, zi] -χ1(d[bi, ai, zi], a[ai]) , 0.0 )
        bdriftᴮ = min( sc[bi, ai, zi], 0.0 ) + min( -d[bi, ai, zi] -χ1(d[bi, ai, zi], a[ai]) , 0.0 )
        #== bdrift for KF ==#
        bi <  bn && ( budriftᶠ = max( sc[bi, ai, zi] - d[bi, ai, zi] - χ1(d[bi, ai, zi], a[ai]) , 0.0 ) )
        bi <  bn && ( budriftᴮ = min( sc[bi, ai, zi] - d[bi, ai, zi] - χ1(d[bi, ai, zi], a[ai]) , 0.0 ) )
        bi == bn && ( budriftᶠ = max( sc[bi, ai, zi] - d[bi-1, ai, zi] - χ1(d[bi-1, ai, zi], a[ai]) , 0.0 ) )
        bi == bn && ( budriftᴮ = min( sc[bi, ai, zi] - d[bi-1, ai, zi] - χ1(d[bi-1, ai, zi], a[ai]) , 0.0 ) )

        #== adrift ==#
        adriftᶠ = max( d[bi, ai, zi], 0.0 ) + netainc[ai,zi]
        adriftᴮ = min( d[bi, ai, zi], 0.0 )
        # for KF
        bi <  bn && ( audriftᶠ = max( d[bi, ai, zi] + netainc[ai,zi], 0.0 ) )
        bi <  bn && ( audriftᴮ = min( d[bi, ai, zi] + netainc[ai,zi], 0.0 ) )
        bi == bn && ( audriftᶠ = max( d[bi-1, ai, zi] + netainc[ai,zi], 0.0 ) )
        bi == bn && ( audriftᴮ = min( d[bi-1, ai, zi] + netainc[ai,zi], 0.0 ) )

        #== REVIEW:100 BOUNDARY ADJUSTMENT 02 ==#
        # bi <  bn && ( bdriftᴮ = min( sc[bi, ai, zi], 0.0 ) + min( -d[bi, ai, zi] -χ1(d[bi, ai, zi], a[ai]) , 0.0 ) )
        # bi == bn && ( bdriftᴮ = min( sc[bi, ai, zi], 0.0 ) -d[bi, ai, zi] -χ1(d[bi, ai, zi], a[ai]) )
        #== REVIEW:100 BOUNDARY ADJUSTMENT 03 ==#
        # ai <  an && ( adriftᴮ = min( d[bi, ai, zi], 0.0 ) )
        # ai == an && ( adriftᴮ = d[bi, ai, zi] + netainc[ai,zi] )

        # v[i, j- 1, k]
        if adriftᴮ != 0.0 && ai>1
            row = abzfromabz(ai-1,bi,zi)
            A[row, ijk] = -adriftᴮ * invda[ai-1]
        end

        # v[i, j- 1, k] KF
        if audriftᴮ != 0.0 && ai>1
            row = abzfromabz(ai-1,bi,zi)
            Au[row, ijk] = -audriftᴮ * invda[ai-1]
        end

        # v[i-1, j, k]
        if bdriftᴮ != 0.0 && bi>1
            row = abzfromabz(ai,bi-1,zi)
            A[row,ijk] = -bdriftᴮ*invdb[bi-1]
        end

        # v[i-1, j, k] KF
        if budriftᴮ != 0.0 && bi>1
            row = abzfromabz(ai,bi-1,zi)
            Au[row,ijk] = -budriftᴮ*invdb[bi-1]
        end

        # v[i, j, k]
        val = 0.0
        ( bdriftᴮ != 0.0 && bi > 1 ) && ( val += bdriftᴮ * invdb[bi-1] )
        ( bdriftᶠ != 0.0 && bi < bn) && ( val -= bdriftᶠ * invdb[bi]   )
        ( adriftᴮ != 0.0 && ai > 1)  && ( val += adriftᴮ * invda[ai-1] )
        ( adriftᶠ != 0.0 && ai < an) && ( val -= adriftᶠ * invda[ai]   )
        val != 0.0 && ( A[ijk, ijk] = val )

        # v[i, j, k] KF
        val = 0.0
        ( budriftᴮ != 0.0 && bi > 1 ) && ( val += budriftᴮ * invdb[bi-1] )
        ( budriftᶠ != 0.0 && bi < bn) && ( val -= budriftᶠ * invdb[bi]   )
        ( audriftᴮ != 0.0 && ai > 1)  && ( val += audriftᴮ * invda[ai-1] )
        ( audriftᶠ != 0.0 && ai < an) && ( val -= audriftᶠ * invda[ai]   )
        val != 0.0 && ( Au[ijk, ijk] = val )

        # v[i-1, j, k]
        if bdriftᶠ != 0.0 && bi<bn
            row = abzfromabz(ai,bi+1,zi)
            A[row, ijk] = bdriftᶠ*invdb[bi]
        end

        # v[i-1, j, k] KF
        if budriftᶠ != 0.0 && bi<bn
            row = abzfromabz(ai,bi+1,zi)
            Au[row, ijk] = budriftᶠ*invdb[bi]
        end

        # v[i, j+1, k]
        if adriftᶠ != 0.0 && ai<an
            row = abzfromabz(ai+1,bi,zi)
            A[row, ijk] = adriftᶠ * invda[ai]
        end

        # v[i, j+1, k] KF
        if audriftᶠ != 0.0 && ai<an
            row = abzfromabz(ai+1,bi,zi)
            Au[row, ijk] = audriftᶠ * invda[ai]
        end

    end
    #== set A = A + Λ ==#
    # NOTE:70 Notice that all matrices are initialized transposed
    #        they are going to be transposed at the END
    Λvals = nonzeros(Λ)
    Avals = nonzeros(A)
    copy!(Avals, Avals + Λvals)

    #== Set B = diag(invΔ + ρ) - A ==#
    # TODO: think of a better way to do this part
    Bvals = nonzeros(B)
    fill!(Bvals, zero(Float64))
    Brows = rowvals(B)
    ijk = zero(Int)
    @inbounds for zi in 1:zn, ai in 1:an, bi in 1:bn
        ijk += 1
        for k in nzrange(B, ijk)
            # loop over elements in the column ijk
            row = Brows[k]
            Bvals[k] = - Avals[k]
            if row == ijk
                Bvals[k] += invΔ + ρ
            end
        end
    end

    # solve for new V (NOTE:280 the B transposed)
    V[:] = B' \ b̃

    return Void
end



##############################
##
## TwoAssetsFD
##
# soldated
##############################
# """
# Use finite difference method to solve for HJB equation.
# """
# function solve_hjb!(fd::TwoAssetsFD;
#                     maxit::Int = 100,
#                     crit::Float64 = 1e-5,
#                     verbose::Bool = true)
#
#     for iter in 1:maxit
#
#         ## UPDATE newV
#         updateV!(fd)
#
#         # check convergence
#         # distance = chebyshev(vec(fd.newV), vec(fd.V))
#         distance = norm(vec(fd.newV)-vec(fd.V),Inf)
#         if distance < crit
#             if verbose
#                 println("hjb solved : $(iter) iterations")
#             end
#             break
#         else
#             # update V using newV
#             fd.V = fd.newV
#             @printf("Value function iteration %d, distance %.4f \n", iter,distance)
#         end
#     end
# end
#
# """
#
# Update value function based on MATLAB code
# """
# function updateV!(fd::TwoAssetsFD)
#
#     #== Parameters ==#
#     γ = fd.twoap.γ; ρ = fd.twoap.ρ; ξ = fd.twoap.ξ; rᴬ = fd.twoap.rᴬ;
#     rᴮ⁺ = fd.twoap.rᴮ; w = fd.twoap.w; χ₀ = fd.twoap.χ₀; χ₁ = fd.twoap.χ₁
#     rᴮ⁻ = fd.twoap.rᴮ + fd.twoap.wedge
#
#     #== FD method ==#
#     invΔ = fd.invΔ;
#     z = fd.z; a = fd.a ; b = fd.b
#     Δb = fd.Δb; Δa = fd.Δa
#     u = fd.u;
#     B = fd.B; A = fd.A; Λ = fd.Λ; V = fd.V
#
#     # precompute the lengths
#     bn = length(b)
#     an = length(a)
#     zn = length(z)
#     N = an*bn*zn
#
#     invdb = 1.0/(Δb)
#     invda = 1.0/(Δa)
#     invγ  = 1/γ
#     inv1γ = 1.0/(1.0-γ)
#
#     V = reshape(V, bn, an, zn)
#
#     # set A = Λ
#     # NOTE:20 Notice that all matrices are initialized transposed
#     #        they are going to be transposed at the END
#     Λvals = nonzeros(Λ)
#     Avals = nonzeros(A)
#     copy!(Avals, Λvals)
#
#     # update A and u
#     ij = zero(Int)
#     @inbounds for zi in 1:zn, ai in 1:an, bi in 1:bn
#         ij += 1
#
#         rᴮ = (b[bi]>0)*rᴮ⁺ + ( 1-(b[bi]>0) )*rᴮ⁻ # set the interest rate
#
#         ## == Derivative w.r.t B == ##
#         bi<bn ? ( ∂Vbᶠ = ( V[bi+1, ai, zi] - V[bi, ai, zi] ) * invdb ) : ( ∂Vbᶠ = ( ( 1-ξ ) * w * z[zi] + rᴮ⁺ * b[end] )^(-γ) )
#         bi>1  ? ( ∂Vbᴮ = ( V[bi, ai, zi] - V[bi-1, ai, zi] ) * invdb ) : ( ∂Vbᴮ = ( ( 1-ξ ) * w * z[zi] + rᴮ⁻ * b[1]   )^(-γ) )
#
#         ## == Derivative w.r.t A == ##
#         ai<an ? ( ∂Vaᶠ = ( V[bi, ai+1, zi] - V[bi, ai, zi] ) * invda ) : ∂Vaᶠ = zero(Float64) # make sure below that don't use ∂Vaᶠ
#         ai>1  ? ( ∂Vaᴮ = ( V[bi, ai, zi] - V[bi, ai-1, zi] ) * invda ) : ∂Vaᴮ = zero(Float64)
#
#         ## CONSUMPTION decision ##
#         cᶠ = ∂Vbᶠ^(-invγ); cᴮ = ∂Vbᴮ^(-invγ)
#         scᴮ = (1-ξ) * z[zi] * w + b[bi] * rᴮ - cᴮ
#         scᶠ = (1-ξ) * z[zi] * w + b[bi] * rᴮ - cᶠ
#
#         ## DEPOSIT decision ##
#         dᴮᴮ = foc_dep( ∂Vbᴮ, ∂Vaᴮ, a[ai] )
#         dᴮᶠ = foc_dep( ∂Vbᴮ, ∂Vaᶠ, a[ai] )
#         dᶠᴮ = foc_dep( ∂Vbᶠ, ∂Vaᴮ, a[ai] )
#         dᶠᶠ = foc_dep( ∂Vbᶠ, ∂Vaᶠ, a[ai] )
#
#     # DOUBT: having a hard time to understand these constraints...
#         #    Although the deposits are computed for bi=1,bn, later on
#         #    these are set to zero when doing upwind for b
#
#         #== eq (13),(14) + set the boundaries in a ==#
#         dᴮ = min( dᴮᴮ, 0.0 ) + max(dᴮᶠ, 0.0)     # eq (13)/ note dᴮᴮ > dᴮᶠ
#         ai==1  && ( dᴮ = (dᴮᶠ>0)*dᴮᶠ )              # Moll: make sure d>=0 at amin, don't use VaB(:,1,:)
#         ai==an && ( dᴮ = (dᴮᴮ<0)*dᴮᴮ )              # Moll: make sure d<=0 at amax, don't use VaF(:,J,:)
#
#         dᶠ = min( dᶠᴮ,0.0 ) + max(dᶠᶠ, 0.0)         # eq (14)
#         ai==1  && ( dᶠ = (dᶠᶠ>0)*dᶠᶠ )
#         ai==an && ( dᶠ = (dᶠᴮ<0)*dᶠᴮ )
#
#         sdᴮ = - dᴮ - χ(dᴮ,a[ai])
#         sdᶠ = - dᶠ - χ(dᶠ,a[ai])
#         bi==bn && ( sdᶠ = min(sdᶠ,0.0) )
#
#         # == consumption UPWIND SCHEME == #
#         c₀ = (1-ξ) * z[zi] * w + b[bi] * rᴮ
#         Icᴮ = (scᴮ<0) ; Icᶠ = (scᶠ>0 && scᴮ>=0); Ic⁰ = 1 - Icᴮ - Icᶠ
#
#         c  = Icᴮ*cᴮ + Icᶠ*cᶠ + Ic⁰*c₀
#         u[ij] = inv1γ * c^(1-γ) + invΔ * V[ij]  # notice that V[ij] still valid
#
#         # == dep UPWIND SCHEME for Vb == #
#         Idᶠ = (sdᶠ>0); Idᴮ = (sdᴮ<0 && sdᶠ<=0)
#
#         # NOTE:30 Boundary constraints of b
#         bi==1  && ( Idᴮ =0  )               # make sure don't use ∂Vbᴮ  if ib=1  for deposit decision
#         bi==bn && ((Idᴮ,Idᶠ)=(1,0))         # make sure don't use ∂Vbᶠ  if ib=bn for deposit decision
#                                             # and FORCE ∂Vbᴮ is used at bn for deposit decision
#         Id⁰ = 1 - Idᶠ - Idᴮ                 # INACTION Region
#
#         #== dep UPWINDE SCHEMEq for Va- eq (17) ==#
#         d⁻ = min(Idᴮ * dᴮᴮ + Idᶠ * dᶠᴮ,0.0)       # use backward for ∂Va
#         d⁺ = max(Idᴮ * dᴮᶠ + Idᶠ * dᶠᶠ,0.0)       # use forward for ∂Va
#         mf = d⁺ + ξ*w*z[zi] + rᴬ*a[ai]
#
#         # NOTE:40 Boundary constraints of a
#         #   Aggregate all terms in ∂Vaᴮ at an
#         ai==an && (mf = 0; d⁻ = ξ*w*z[zi] + rᴬ*a[an] + d⁻)
#
#         # NOTE:50 Positions of A that I need to fill
#             # A[ij-bn, ij] - v(i  ,j-1,k) --> - invda*[ d⁻ ]
#             # A[ij- 1, ij] - v(i+1,j  ,k) --> - invdb*[(scᴮ)⁻ + (sdᴮ)⁻ ]
#             # A[ij   , ij] - v(i  ,j  ,k) --> invdb*[ (scᴮ)⁻ - (scᶠ)⁺ + (sdᴮ)⁻ - (sdᶠ)⁺ ] + invda*[ d⁻ - (d⁺ + ξwz + rᵃa) ]
#             # A[ij+ 1, ij] - v(i-1,j  ,k) --> invdb*[ (scᶠ)⁺ + (sdᶠ)⁺ ]
#             # A[ij+bn, ij] - v(i  .j+1,k) --> invda*[ d⁺ + ξwx + rᵃa ]
#             #
#             # the matrix a is changing its entrances... Note that code from aigary only changed when
#             # it was nonzero...
#
#         current = Array(Float64,0)
#         zi == 2   && push!(current, 0.0)
#
#         ij>bn     && ( ai>1  ? push!(current, - invda * d⁻ ) : push!(current, 0.0) )
#         ij>1      && ( bi>1  ? push!(current, - invdb *( Icᴮ*scᴮ + Idᴮ*sdᴮ ) ) : push!(current, 0.0) )
#         push!(current, invdb *( Icᴮ*scᴮ + Idᴮ*sdᴮ - Icᶠ*scᶠ - Idᶠ*sdᶠ ) + invda *( d⁻ - mf ) )
#         ij<N      && ( bi<bn ? push!(current, invdb *( Icᶠ*scᶠ + Idᶠ*sdᶠ ) ) : push!(current,0.0) )
#         ij<=N-bn  && ( ai<an ? push!(current, invda * mf ) : push!(current,0.0) )
#
#         zi == 1   && push!(current, 0.0)
#
#         Avals[nzrange(A, ij)] += current
#     end
#
#     Bvals = nonzeros(B)
#     Brows = rowvals(B)
#     ij = zero(Int)
#     @inbounds for zi in 1:zn, ai in 1:an, bi in 1:bn
#         ij += 1
#         for k in nzrange(B, ij)
#             # loop over elements in the column ij
#             row = Brows[k]
#             Bvals[k] = - Avals[k]
#             if row == ij
#                 Bvals[k] += invΔ + ρ
#             end
#         end
#     end
#
#     #== Solve for new V ==#
#     #NOTE:60 the B transposed
#     fd.newV = B' \ u
#
#     return Void
# end
