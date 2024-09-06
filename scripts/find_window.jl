function dict__projwfc_result(result_lines)
    @inline to_atom_symbol(x) = replace(join(split(x,"(",keepempty=false)[[2,1]].|>strip, ""), 
                                        r"\s*\)\s*"=>"")
    @inline ι(x) = parse(Int,x)
    @inline organize(x) = (ι(x[1])=>(atom=x[2],l=ι(x[4]),m=ι(x[5]),wfc=ι(x[3])))
    @inline process_state_line(l) = [
        replace(l[1], "state #"=>""),
        replace(l[2], "atom"=>"") |> to_atom_symbol,
        replace(l[3], "wfc"=>""),
        replace(l[4], r"l\s*="=>""),
        replace(l[5], r"m\s*="=>""),
    ] .|> strip |> organize 
    @inline sortpair(pl) = sort(pl, by=first)
    @inline cont_list(lst) = [((@assert i==p); v) for (i,(p,v)) in enumerate(lst)]
    state_line_reg = r"state\s*#\s*\d+:\s*atom\s*\d+\s*\(\w+\s*\),\s*wfc\s*\d+\s*\(l=\d+\s*m=\s*\d+\)"
    state_reg = r"state\s*#\s*\d+"
    atom_reg  = r"atom\s*\d+\s*\(\w+\s*\)"
    wfc_reg   = r"wfc\s*\d+"
    l_reg     = r"l\s*=\s*\d+"
    m_reg     = r"m\s*=\s*\d+"
    #TODO meaning of "wfc 2" ?
    info = QuantumEspressoTools.extract_all(
        result_lines, 
        state_line_reg, 
        [state_reg, atom_reg, wfc_reg, l_reg, m_reg]
    ) .|> process_state_line |> sortpair |> cont_list
    return Dict(
        "states" => info,
        "Lowdin Charges" => []
    )
end



function floodfill(
    BAND::Vector{Float64}, 
    OCC::Vector{Bool}, 
    ENREF0::Float64,
    dEmax = 0.1
    )
    # correct E_ref
    (_, k0) = findmin(abs.(BAND[OCC].-ENREF0))
    ENREF = BAND[OCC][k0]
    if abs(ENREF-ENREF0)>dEmax
        @warn "floodfill() has got wrong reference energy $ENREF0 and uses $ENREF."
    end
    # flood fill
    fmax = ENREF
    fmin = ENREF
    ord = sortperm(BAND)
    (v0, p0) = findmin(abs.(BAND.-ENREF))
    i0 = findfirst(x->x==p0, ord)
    pmax = findfirst(i->(!OCC[ord[i0+i]]), 1:length(ord)-i0)
    (pmax===nothing && (pmax = (length(ord)-i0) ;))
    fmax = BAND[ord[i0+pmax-1]]
    pmin = findfirst(i->(!OCC[ord[i0-i]]), 1:i0-1)
    (pmin===nothing && (pmin = i0-1 ;))
    fmin = BAND[ord[i0-pmin]]
    return fmin, fmax, pmin+pmax-1
end



#TODO 2D findgap for energy-weight plot
#TODO see test_proj.jl

# find energy interval in BAND = BAND0[OCC]
# which contains 
function findgap_core(
    BAND0::Vector{Float64}, 
    OCC::Vector{Bool}, 
    ENREF::Float64,
    δE::Float64
    )
    BAND = BAND0[OCC]
    ord = sortperm(BAND)
    (v0, p0) = findmin(abs.(BAND.-ENREF))
    fmax = ENREF
    fmin = ENREF
    i0 = findfirst(x->x==p0, ord)
    pmax = findfirst(i->abs(BAND[ord[i0+i]]-BAND[ord[i0+i+1]])>δE, 1:length(ord)-i0-1)
    (pmax===nothing && (pmax = (length(ord)-i0) ;))
    fmax = BAND[ord[i0+pmax]]
    pmin = findfirst(i->abs(BAND[ord[i0-i]]-BAND[ord[i0-i-1]])>δE, 1:i0-2)
    (pmin===nothing && (pmin = i0 ;))
    fmin = BAND[ord[i0-pmin+1]]
    return fmin, fmax, pmin+pmax-1
end


#TODO validate
function findgap_upscaling_dE(
    BAND0::Vector{Float64}, 
    OCC::Vector{Bool}, 
    ENREF::Float64,
    δE0::Float64
    )
    @inline isclose(t1,t2) = abs(t1[1]-t2[1])<1e-3 && abs(t1[2]-t2[2])<1e-3 && abs(t1[3]-t2[3])<3
    t0 = findgap_core(BAND0, OCC, ENREF, δE0)
    close_count = 0
    for λ ∈ 1.5:0.5:12.0
        #* based on the observation that 
        #* upscaling dE leads to stable output of findgap_core()
        t1 = findgap_core(BAND0, OCC, ENREF, λ*δE0) 
        if isclose(t0,t1)
            close_count += 1
            if close_count>2
                break
            end
        else
            t0 = t1
            close_count = 0
        end
    end
    return t0
end


findgap(
    BAND0::Vector{Float64}, 
    OCC::Vector{Bool}, 
    ENREF::Float64,
    δE0::Float64
    ) = findgap_upscaling_dE(BAND0, OCC, ENREF, δE0)


const PROJWFC_PROJ_TYPE = NamedTuple{(:atom, :l, :m, :orbit), Tuple{String, Int64, Int64, String}}


function calculate_occupation(
    projwfcx_output_lines::Vector{T0}, # 
    projectors::Vector{PROJWFC_PROJ_TYPE},
    proj_lines::Vector{T2}
    ) where {T0 <: AbstractString, T2 <: AbstractString}
    d = dict__projwfc_result(projwfcx_output_lines)
    P, all_states = projwfcx_output_projwfc_up(proj_lines) #TODO slooooooooooow
    @inline match_proj(p) = findall(x->(x.l==p.l&&x.m==p.m&&x.atom==p.atom&&x.orbit==p.orbit), all_states)
    states_id_to_proj0 = match_proj.(projectors)
    states_id_to_proj  = vcat([x for x in states_id_to_proj0 if x!==nothing]...) |> sort
    # @inline match_proj1(p) = findall(x->(x.l==p.l&&x.m==p.m&&x.atom==p.atom), d["states"])
    # states_id_to_proj1 = match_proj1.(projectors)
    # states_id_to_proja = vcat([x for x in states_id_to_proj1 if x!==nothing]...) |> sort
    #@assert states_id_to_proja==states_id_to_proj
    @inline occ(i,j) = sum(P[i,j,states_id_to_proj])
    @inline tot(i,j) = sum(P[i,j,:])
    OCC_ABS = [occ(i,j) for i=1:size(P,1), j=1:size(P,2)]
    OCC_REL = OCC_ABS ./ [tot(i,j) for i=1:size(P,1), j=1:size(P,2)]
    return OCC_ABS, OCC_REL
end

##

count_wann(OCC::Array{Bool,2}) = minimum([sum(OCC[ikpoint,:]) for ikpoint=1:size(OCC,1)])

count_wann(OCC::BitArray{2}) = count_wann(Array{Bool,2}(OCC))

#* OKay for Mo2N
function find_window0(
    projwfc_output_lines::Vector{T0},
    nbnd::Int,
    nwann::Int,
    projectors::Vector{PROJWFC_PROJ_TYPE},
    band_lines::Vector{T1}, 
    projwfc_up_lines::Vector{T2};
    thr = (0.99, 0.2),
    E_ref = 0.0,
    dE = 0.3
    ) where {T0 <: AbstractString, T1 <: AbstractString, T2 <: AbstractString}

    @assert length(projectors)==nwann
    BANDS =  hcat(last.(values.(pw_bands(band_lines))) ...)'

    OCC_ABS, OCC_REL = calculate_occupation(
            projwfc_output_lines,
            projectors,
            projwfc_up_lines
    )
    @assert size(BANDS) == size(OCC_ABS) == size(OCC_REL)

    #* ----------------------------------------
    #* outer window use OCC_ABS
    #* select the bands with "significant" weights : abs > thr[2]
    #* thr[2] should be low such that the outer window contains 
    #* enough (not all!) bands to disentangle
    #TODO how to INCREASE `thr[2]` to tighten the outer window?
    thr_2 = thr[2]
    OCCUPATION_POSSIBLE = Array{Bool,2}(OCC_ABS.>thr_2)
    outer_min, outer_max, outer_num_en = findgap(BANDS[:], OCCUPATION_POSSIBLE[:], E_ref, dE)
    num_try = 3
    for itry = 1:num_try
        if count_wann(OCCUPATION_POSSIBLE) < nwann
            @error "find_window0() : too few states of significant weight. \nRetry with reduced threshold from $(thr_2) to $(0.5*thr_2)."
            thr_2 = 0.5 * thr_2
        else
            break
        end
        OCCUPATION_POSSIBLE = Array{Bool,2}(OCC_ABS.>thr_2)
        outer_min, outer_max, outer_num_en = findgap(BANDS[:], OCCUPATION_POSSIBLE[:], E_ref, dE)
    end

    #* ----------------------------------------
    #* inner window uses OCC_REL and OCC_ABS 
    #* (select ONLY WITHIN the bands with significant weight)
    #* thr[1] should be high ~0.5 such that the inner frozen window is tight
    #* and the bands inside the window contains bands with large weight percentage
    #TODO automatically relax thr[1] to increase `nwann1`, for example bisection
    thr_1 = thr[1]
    OCCUPATION_HIGH = Array{Bool,2}((OCC_ABS.>thr_2) .* (OCC_REL.>thr_1))
    inner_min, inner_max, inner_num_en = floodfill(BANDS[:], OCCUPATION_HIGH[:], E_ref)
    @info "find_window0() : found $inner_num_en states in the frozen window. "

    #* ----------------------------------------
    ret =  Dict( "num_bands"   => nbnd,
                 "num_wann"    => nwann,
                 "dis_win_min" => outer_min-1e-3,
                 "dis_win_max" => outer_max+1e-3,
                 "dis_froz_min"=> inner_min+1e-3,
                 "dis_froz_max"=> inner_max-1e-3 
    )
    @info "find_window0() result : $(showdict(ret))"
    return ret
end

##

function set_window0(
    projwfcx_output_lines::Vector{T0}, # 
    nbnd::Int,
    nwann::Int,
    inner_window::Tuple{Float64,Float64},
    outer_window::Tuple{Float64,Float64}
    ) where {T0 <: AbstractString}
    ret =  Dict( "num_bands"   => nbnd,
                 "num_wann"    => nwann,
                 "dis_win_min" => outer_window[1]-1e-3,
                 "dis_win_max" => outer_window[2]+1e-3,
                 "dis_froz_min"=> inner_window[1]+1e-3,
                 "dis_froz_max"=> inner_window[2]-1e-3
    )
    @info "set_window0() result : $(showdict(ret))"
    return ret
end

##
# test
# 
# using Pkg
# Pkg.activate("/home/dabajabaza/jianguoyun/Nutstore/QuantumEspressoTools")
# using QuantumEspressoTools
# 
# ##
# ##
# 
# global const WKSPC = "/tmp/band_0.9_0.9/ps___PSL_USPP_PBESOL_SR___C-n-1.0.0___kp_20,16,1,0,0,0_cut_60.0,480.0"
# global const bandfn = "$WKSPC/strain_rx_0.9_ry_0.9_nscf/strain_rx_0.9_ry_0.9.pw.x.out"
# global const projfn = "$WKSPC/strain_rx_0.9_ry_0.9_projwfc/strain_rx_0.9_ry_0.9.proj.projwfc_up"
# global const projretfn = "$WKSPC/strain_rx_0.9_ry_0.9_projwfc/strain_rx_0.9_ry_0.9.projwfc.x.out"
# 
# #dd = dict__projwfc_result(readlines(projretfn))
# #bb = pw_bands(readlines(bandfn))
# 
# ##
# 
# wd = find_window0(
#     readlines(projretfn), 
#     32,
#     12,
#     [(atom="C$k",l=1,m=1) for k=1:6],
#     readlines(bandfn), 
#     readlines(projfn);
#     E_ref = pw_fermi_energy_eV(readlines(bandfn))
# )
# 
##
