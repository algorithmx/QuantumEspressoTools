
function floodfill(
    BAND::Vector{Float64}, 
    OCC::Vector{Bool}, 
    ENREF::Float64
    )
    ord = sortperm(BAND)
    (v0, p0) = findmin(abs.(BAND.-ENREF))
    if !OCC[p0]
        @error "floodfill() has got wrong energy $ENREF"
        return ENREF-v0-0.1, ENREF+v0+0.1, 0
    else
        fmax = ENREF
        fmin = ENREF
        i0 = findfirst(x->x==p0, ord)
        pmax = findfirst(i->(!OCC[ord[i0+i]]), 1:length(ord)-i0)
        if (pmax===nothing) 
            pmax = (length(ord)-i0)
        end
        fmax = BAND[ord[i0+pmax-1]]
        pmin = findfirst(i->(!OCC[ord[i0-i]]), 1:i0-1)
        if (pmin===nothing) 
            pmin = i0-1
        end
        fmin = BAND[ord[i0-pmin+1]]
        return fmin, fmax, pmin+pmax-1
    end
end

count_wann(OCC::Array{Bool,2}) = minimum([sum(OCC[ikpoint,:]) for ikpoint=1:size(OCC,1)])

count_wann(OCC::BitArray{2}) = count_wann(Array{Bool,2}(OCC))


#! not good
function find_window0__version_I(
    projwfcx_output_lines::Vector{T0},
    nbnd::Int,
    nwann::Int,
    projectors::Vector{PROJWFC_PROJ_TYPE},
    band_lines::Vector{T1}, 
    proj_lines::Vector{T2};
    thr = 0.5,
    E_ref = 0.0
    ) where {T0 <: AbstractString, T1 <: AbstractString, T2 <: AbstractString}
    bands =  hcat(last.(values.(pw_bands(band_lines))) ...)'
    OCC = Array{Bool,2}(
        calculate_occupation(
            projwfcx_output_lines,
            projectors,
            proj_lines
        ) .> thr
    )
    @assert size(bands) == size(OCC)
    outer_min = minimum(bands[OCC])
    outer_max = maximum(bands[OCC])

    #* ----------------------------------------
    #* sweep method to find best fmin, fmax of the frozen window
    fmin_fmax_list0 =  [ floodfill(bands[:], OCC[:], En) 
                        for En ∈ outer_min.+((outer_max-outer_min)/100).*collect(3:97) ]
    fmin_fmax_list  = [(a,b,c) for (a,b,c) in fmin_fmax_list0 if c > nwann+1]
    fmin_fmax_list1 = sort(fmin_fmax_list, by=x->(x[2]-x[1]) * max(1,x[3]-10))
    @info "XXXXXXXXXXXXXXXXXXXXXXXXXXX"
    @show fmin_fmax_list0
    @show fmin_fmax_list
    @show fmin_fmax_list1
    fmin, fmax = fmin_fmax_list1[end]
    #* ----------------------------------------

    nwann1 = count_wann(OCC)
    if nwann != nwann1
        @warn "find_window0() found, ***at thr=$(thr)***, count_wann(OCC)=$(nwann1) but got nwann=$(nwann) ."
    end
    ret =  Dict( "num_bands"   => nbnd,
                 "num_wann"    => nwann,
                 "dis_win_min" => outer_min-1e-3,
                 "dis_win_max" => outer_max+1e-3,
                 "dis_froz_min"=> fmin+1e-3,
                 "dis_froz_max"=> fmax-1e-3 
    )
    return ret
end


#* OKay for Mo2N
function find_window0__version_II(
    projwfcx_output_lines::Vector{T0},
    nbnd::Int,
    nwann::Int,
    projectors::Vector{PROJWFC_PROJ_TYPE},
    band_lines::Vector{T1}, 
    proj_lines::Vector{T2};
    thr = (0.9, 0.05),
    E_ref = 0.0
    ) where {T0 <: AbstractString, T1 <: AbstractString, T2 <: AbstractString}
    bands =  hcat(last.(values.(pw_bands(band_lines))) ...)'
    OCC_ABS, OCC_REL = calculate_occupation(
            projwfcx_output_lines,
            projectors,
            proj_lines
    )
    @assert size(bands) == size(OCC_ABS)

    #* ----------------------------------------
    #* outer window use OCC_ABS
    #* select the bands with "significant" weights : abs > thr[2]
    #* thr[2] should be low such that the outer window contains 
    #* enough (not all!) bands to disentangle
    #TODO how to increase `thr[2]` to tighten the outer window?
    OCC01 = Array{Bool,2}(OCC_ABS .> thr[2])
    outer_min = minimum(bands[OCC01])
    outer_max = maximum(bands[OCC01])

    #* ----------------------------------------
    #* inner window uses OCC_REL and OCC_ABS 
    #* (select ONLY WITHIN the bands with significant weight)
    #* thr[1] should be high ~0.5 such that the inner frozen window is tight
    #* and the bands inside the window contains bands with large weight percentage
    OCC = Array{Bool,2}((OCC_REL.>thr[1]) .* (OCC_ABS .> thr[2]))
    #TODO automatically relax thr[1] to increase `nwann1`, for example bisection
    nwann1 = count_wann(OCC)
    if nwann != nwann1
        @warn "find_window0() found, ***at thr=$(thr)***,\n count_wann(OCC)=$(nwann1) but got nwann=$(nwann) \nin the frozen window."
    end
    fmin, fmax, num_en = floodfill(bands[:], OCC[:], E_ref)
    if num_en < nwann
        @error "find_window0() : too few states in the frozen window."
    end
    #* ----------------------------------------
    
    ret =  Dict( "num_bands"   => nbnd,
                 "num_wann"    => nwann,
                 "dis_win_min" => outer_min-1e-3,
                 "dis_win_max" => outer_max+1e-3,
                 "dis_froz_min"=> fmin+1e-3,
                 "dis_froz_max"=> fmax-1e-3 
    )
    return ret
end
