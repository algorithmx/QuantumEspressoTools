#TODO () => []

function make_PSINFO(ps)::PSEUDO
    @inline _tk_(x) = strip(x, ['"',']','['])
    @inline trim_key_dic(x) = Dict{Tuple{String,String},Tuple{String,Float64,Float64}}(
                                    (_tk_.(SPLTA(k))...,)=>(v...,) for (k,v) ∈ x
                              )
    ps_DATA = Dict(  (startswith(k,"PSL_") 
                            ? k=>Dict(elem=>trim_key_dic(d) for (elem,d) ∈ v) 
                            : k=>v)
                        for (k,v) ∈ ps["DATA"]  )
    PSEUDO(ps["ROOT"], ps["PATHS"], ps_DATA)
end

function make_PSINFO_change_tuple_to_vec(ps::PSEUDO)::PSEUDO
    @inline tuple_to_vec(x) = Dict{Vector{String},Vector{Any}}([k...,]=>[v...,] for (k,v) ∈ x)
    ps_DATA = Dict(  (startswith(k,"PSL_") 
                            ? k=>Dict(elem=>tuple_to_vec(d) for (elem,d) ∈ v) 
                            : k=>v)
                        for (k,v) ∈ ps.DATA  )
    PSEUDO(ps.ROOT, ps.PATHS, ps_DATA)
end


#function build_load_pseudo(pseudo_path, pseudo_fn="pseudopotential_info.jld2")::PSEUDO
function build_load_pseudo(pseudo_path, pseudo_fn="pseudopotential_info.json")::PSEUDO
    #* banner
    show_banner()
    _ROOT_   = ENV["PSEUDO_ROOT"]
    _PSINFO_ = PSEUDO( _ROOT_, 
                        pseudo_path,
                        Dict(k=>Dict() for k in keys(pseudo_path)) )

    #
    @inline extract_upf(psm,f) = (last(SPLTX(f,"/")), Find_suggested_cutoff(psm,f)...)
    @inline interprete_sssp(v) = (v["filename"], v["cutoff"], v["cutoff"]*v["dual"])

    if isfile("$(_ROOT_)/$pseudo_fn")
        #@load ("$(_ROOT_)/$pseudo_fn") _PSINFO_
        _PSINFO_ = JSON.parse(read("$(_ROOT_)/$pseudo_fn", String)) |> make_PSINFO
        if __VERIFY__PS__ 
            if verify_ps(_PSINFO_, PSEUDO_PATHS)
                @info "Found PSEUDO INFO FILE at $(_ROOT_)/$pseudo_fn ."
            else
                @info "Found PSEUDO INFO FILE at $(_ROOT_)/$pseudo_fn but it is inconsistent."
            end
        else
            @info "Found PSEUDO INFO FILE at $(_ROOT_)/$pseudo_fn without verification."
        end
        return _PSINFO_
    else
        #*  SSSP 
        # https://www.materialscloud.org/discover/sssp/table/precision
        # https://www.materialscloud.org/discover/sssp/table/efficiency
        j_sssp_prec = JSON.parsefile("$(_ROOT_)/sssp_precision.json")
        j_sssp_eff  = JSON.parsefile("$(_ROOT_)/sssp_efficiency.json")
        _PSINFO_.DATA["SSSP_Precision"]  = Dict( k=>interprete_sssp(v) for (k,v) ∈ j_sssp_prec )
        _PSINFO_.DATA["SSSP_Efficiency"] = Dict( k=>interprete_sssp(v) for (k,v) ∈ j_sssp_eff  )

        j_sssp_prec1 = JSON.parsefile("$(_ROOT_)/sssp_precision_pbesol.json")
        j_sssp_eff1  = JSON.parsefile("$(_ROOT_)/sssp_efficiency_pbesol.json")
        _PSINFO_.DATA["SSSP_Precision_PBEsol"]  = Dict( k=>interprete_sssp(v) for (k,v) ∈ j_sssp_prec1 )
        _PSINFO_.DATA["SSSP_Efficiency_PBEsol"] = Dict( k=>interprete_sssp(v) for (k,v) ∈ j_sssp_eff1  )

        #*  SG15
        # http://www.quantum-simulation.org/potentials/sg15_oncv/
        _PSINFO_.DATA["SG15"] = Dict(   k => ("$(k)_ONCV_PBE-1.2.upf",1.2*v[2],1.2*v[3])
                                        for (k,v) ∈ _PSINFO_.DATA["SSSP_Precision"]   )

        #*  OTHER 
        #*  Single file in the folder
        for _p_ in [    "SG15_FR", "MT_PBE", "MT_PW_LDA",
                        "GHH_PBE","GHH_PBE_SP","GHH_PZ","GHH_PZ_SP","GHH_BLYP","GHH_BLYP_SP",
                        "GBRV_LDA","GBRV_PBE","GBRV_PBESOL"   ]
            _PSINFO_.DATA[_p_] = Dict(  Find_element_name("$(PSEUDO_PATHS[_p_])/$f") => extract_upf(_p_,"$(PSEUDO_PATHS[_p_])/$f") 
                                        for f ∈ readdir(PSEUDO_PATHS[_p_])  )
        end
        #*  Multiple files in the folder
        for _p_ in [    "PSL_USPP_LDA_SR", "PSL_USPP_LDA_FR", "PSL_USPP_PBE_SR", "PSL_USPP_PBE_FR", "PSL_USPP_PBESOL_SR", "PSL_USPP_PBESOL_FR",
                        "PSL_PAW_PZ_SR"  , "PSL_PAW_PZ_FR"  , "PSL_PAW_PBE_SR" , "PSL_PAW_PBE_FR" , "PSL_PAW_PBESOL_SR" , "PSL_PAW_PBESOL_FR"   ]
            all_fn       = readdir(PSEUDO_PATHS[_p_])
            # map : file name => element name
            file2elem    = Dict(f=>Find_element_name("$(PSEUDO_PATHS[_p_])/$f") for f ∈ all_fn)
            # pair (f,x) = (file name, settings)
            all_settings = [(f,x) for f ∈ all_fn for x ∈ ALL_PSL_SETTINGS if (x[1] ⊂ f) && (x[2] ⊂ f)]
            # organize fn by elements
            elms         = Dict(el=>[f for f in all_fn if file2elem[f]==el] for el ∈ unique(values(file2elem)))
            # structure: PSL_XXX  >>>  element  >>>  setting  >>>  extract_upf
            _PSINFO_.DATA[_p_] = Dict(  el => Dict( x => extract_upf(_p_,"$(PSEUDO_PATHS[_p_])/$f") 
                                                    for (f,x) ∈ all_settings  if f ∈ fns )
                                        for (el,fns) in elms  )
            # deal with H
            if "H" ∈ keys(elms)
                _PSINFO_.DATA[_p_]["H"]  = Dict( ("",v)=>extract_upf(_p_,"$(PSEUDO_PATHS[_p_])/$f") 
                                                    for f ∈ elms["H"]  
                                                        for v ∈ __ALL_VERS if (v ⊂ f) )
            end
            # deal with He
            if "He" ∈ keys(elms)
                _PSINFO_.DATA[_p_]["He"] = Dict( ("",v)=>extract_upf(_p_,"$(PSEUDO_PATHS[_p_])/$f") 
                                                    for f ∈ elms["He"]
                                                        for v ∈ __ALL_VERS if (v ⊂ f) )
            end
            # deal with Li
            if "Li" ∈ keys(elms)
                nothing
                #TODO
            end
        end
        if __MODIFY__PS__
            #@save ("$(_ROOT_)/$pseudo_fn")  _PSINFO_
            #@info "Saving PSEUDO INFO FILE at $(_ROOT_)/pseudopotential_info.jld2."
            JSON.json(make_PSINFO_change_tuple_to_vec(_PSINFO_)) >> ("$(_ROOT_)/$pseudo_fn")
            @info "Saving PSEUDO INFO FILE at \"$(_ROOT_)/pseudopotential_info.json\"."
        end
        return _PSINFO_
    end
end


function verify_ps(PS::PSEUDO, pseudo_path)::Bool
    consis = true

    @info "pseudopotential.jl verifies paths ..."
    if length(keys(PS.PATHS)) != length(keys(pseudo_path))
        AB = setdiff(keys(PS.PATHS),    keys(pseudo_path))
        BA = setdiff(keys(pseudo_path), keys(PS.PATHS))
        println("Inconsistency : paths in the jld2 file not found in pseudo root :" *join(AB, "\n"))
        println("Inconsistency : paths in the pseudo root not found in jld2 file :" *join(BA, "\n"))
        consis = false
    end

    @info "pseudopotential.jl verifies PS Library ..."
    for __p__ in [  "PSL_USPP_LDA_SR", "PSL_USPP_LDA_FR", "PSL_USPP_PBE_SR", "PSL_USPP_PBE_FR", "PSL_USPP_PBESOL_SR", "PSL_USPP_PBESOL_FR",
                    "PSL_PAW_PZ_SR"  , "PSL_PAW_PZ_FR"  , "PSL_PAW_PBE_SR" , "PSL_PAW_PBE_FR" , "PSL_PAW_PBESOL_SR" , "PSL_PAW_PBESOL_FR"   ]
        for e in keys(PS.DATA[__p__])
            A = length([k for k in readdir(PS.PATHS[__p__]) if startswith(k,"$(e).")])
            B = length(PS.DATA[__p__][e])
            if A!=B
                print("Inconsistency : ")
                println((__p__,e,A,B))
                consis = false
            end
        end
    end

    return consis
end

