#* -------------------------------------------------------------------

@inline j_s(x) = join(x,"  ")
@inline p_p(x,t=0) = repeat(" ",t)*"$x"
@inline j_n(x,t=0) = join(pp.(x,t),"\n")


function adaptg(x)
    if x isa AbstractString
        return "$(x)"
    elseif x isa Integer
        return string(x)
    elseif x isa Float64
        return strip((@sprintf "%20.15f" x))
    else
        return string(x)
    end
end

#* -------------------------------------------------------------------

function gulp_title(config_gulp)::Vector{String}
    title = config_gulp[:title]
    return ["title", ((title isa Vector) ? j_n(title) : title), "end"]
end


function is_valid(config_gulp, mode)
    #TODO
    if mode==:fit
        return occursin("fit",config_gulp[:control_line])
    elseif mode==:phonon
        return occursin("phon",config_gulp[:control_line])
    else
        return true
    end
end


function gulp_first_line(config_gulp, mode)::Vector{String}
    @assert is_valid(config_gulp, mode)
    return [config_gulp[:control_line],]
end


function gulp_structure_info(config_gulp)::Vector{String}

    # helper functions
    @inline spln(s) = split(s,"\n",keepempty=false)
    @inline make_frac_line(f) = (   (f[1] isa String) 
                                    ? f[1]*"  "*j_s(f[2])   # "Sc"=>[0,0,0]
                                    : ((f[1][2] isa String) 
                                        ? "$(f[1][1])  $(f[1][2])  "*j_s(f[2])  # ("Sc","core")=>[0,0,0] or ("Sc","shel")=>[0,0,0]
                                        : j_n(["$(f[1][1])  $t  "*j_s(f[2]) for t ∈ f[1][2]]))  ) # ("Sc",["core","shel"])=>[0,0,0]
    #
    @inline make_spec_line(f) = (   (f[1] isa String) 
                                    ? f[1]*"  "*j_s(f[2])   # "Sc"=>3.0
                                    : ((f[1][2] isa String) 
                                        ? "$(f[1][1])  $(f[1][2])  "*j_s(f[2])  # ("Sc","core")=>3.0 or ("Sc","shel")=>3.0
                                        : j_n(["$(f[1][1])  $t  "*j_s(f[2]) for t ∈ f[1][2]]))  ) # ("Sc",["core","shel"])=>3.0
    # cell
    cell       = "cell\n    "* j_s(config_gulp[:cell])
    space      = "space\n    $(get(config_gulp,:IT,1))"
    
    # fractional
    frac_lines = spln(j_n(make_frac_line.(config_gulp[:fractional])))
    nfrac      = length(frac_lines)
    fractional = "fractional   $nfrac\n" * j_n(frac_lines, 4)

    # fractispeciesonal
    spec_lines = spln(j_n(make_spec_line.(config_gulp[:species])))
    nspec      = length(spec_lines)
    species    = "species  $nspec\n" * j_n(spec_lines, 4)

    # library
    lib        = ((:library ∈ keys(config_gulp)) ? "library\n$(config_gulp[:library])" : "")

    return [cell, space, fractional, species, lib]

end


VVF64 = Vector{Vector{Float64}}
VTupIFI = Vector{Tuple{Int,Float64,Int}}
function select_fitting_frequencies(
    all_modes::Vector{Pair{Vector{Float64},Vector{TupF64VC64}}}, 
    selection::Dict{Vector{Float64},Vector{Int}};
    )::Tuple{VVF64, VTupIFI}

    @inline searchQ(x) = findfirst(y->norm(x-y[1])<1e-5, all_modes)

    SELECTED_QP_FREQ = []
    cnt_q = 0
    for (q,mode_i) ∈ selection
        iq = searchQ(q)
        if iq !== nothing
            cnt_q += 1
            tmp = Tuple{Int64,Float64,Int64}[]
            for ii ∈ mode_i
                (f,m) = last(all_modes[iq])[ii]
                push!(tmp, (ii, f, cnt_q))
            end
            push!(SELECTED_QP_FREQ, q=>tmp)
        end
    end
    return first.(SELECTED_QP_FREQ), vcat(last.(SELECTED_QP_FREQ)...)
end


TupF64VC64  = Tuple{Float64,Vector{ComplexF64}}


function select_fitting_modes(
    natoms::Int,
    all_modes::Vector{Pair{Vector{Float64},Vector{TupF64VC64}}}, 
    selection::Dict{Vector{Float64},Vector{Int}};
    atom_order = [],
    IMCHAR = "+ i"
    )::Vector{Pair{String,Vector{String}}}

    # re-ordering of the loaded mode to be aligned to the GULP convention
    AORD = length(unique(atom_order))==natoms ? atom_order : collect(1:natoms)

    @inline searchQ(x) = findfirst(y->norm(x-y[1])<1e-5, all_modes)
    #@inline c2s(x) = (abs(imag(x))<1e-5 ? "$(real(x))" : "$(real(x)) $(IMCHAR) $(imag(x))")
    @inline c2s(x) = " $(real(x)) $(imag(x)) "
    SELECTED_MODES = []
    for (q,mode_i) ∈ selection
        iq = searchQ(q)
        if iq !== nothing
            for (f,m) ∈ last(all_modes[iq])[mode_i]
                #@assert sum(abs.(imag.(m))) < 1e-5  "q = $q ; f = $f ; m = $m"
                if sum(abs.(imag.(m))) < 1e-5
                    md = String["$f", ]
                    for i ∈ AORD
                        push!(md,"$(j_s(real.(m[3(i-1)+1:3i])))")
                        #push!(md,"$(j_s(c2s.(m[3(i-1)+1:3i])))")
                    end
                    push!(SELECTED_MODES, "mode"=>md)
                end
            end
        end
    end
    return SELECTED_MODES

end


function gulp_observables(config_gulp)::Vector{String}
    lines = ["observables",]
    kw = ["mode", "energy", "stress", "frequency"]
    for i  = 1:length(config_gulp[:observables])
        key_obs  = first(config_gulp[:observables][i])
        val_obs0 = last(config_gulp[:observables][i])
        val_obs  = String[]
        if eltype(val_obs0) <: String
            val_obs = val_obs0
        elseif eltype(val_obs0) <: Tuple || eltype(val_obs0) <: Vector
            @inline _ltostr_(x) = j_s(adaptg.(x))
            val_obs = _ltostr_.(val_obs0)
        end
        if any([occursin(k,key_obs) for k in kw])
            lines = [lines; [p_p(key_obs,4),]; pp.(val_obs,8)]
        end
    end
    return [lines; "end"]
end


function gulp_specific(config_gulp)::Vector{String}
    #> stepmax opt 1.0
    #> shrink 8 8 8
    @inline _ltostr_(x) = ((x isa String) ? x : j_s(adaptg.(x)))
    return _ltostr_.(config_gulp[:specific])
end


function gulp_potential(
    pot_name::String,
    elem_tuple,
    specs...
    )::Vector{String}
    @inline trzo(x) = (occursin(r"\d+\.[0-9]+0+$",x) ? replace(x,r"0+$"=>"")*"0" : x)
    nelem = length(elem_tuple)
    @inline make_elem_core_shell_word(elem) = ((elem isa String) 
                                                ? (@sprintf "%2s %5s"  elem  "core")
                                                : ( (elem[2] isa String) 
                                                    ? (@sprintf "%2s %5s"  elem[1]  elem[2])
                                                    : j_s([(@sprintf "%2s %5s"  elem[1]  t) for t ∈ elem[2]]) ))
    elem_word = j_s(make_elem_core_shell_word.(elem_tuple))
    return [
        pot_name,
        repeat(" ",4) * (elem_word * "  " * j_s(trzo.(adaptg.(specs))))
    ]
end


function gulp_potential_pairs(config_gulp)::Vector{String}
    pdic = config_gulp[:potentials]
    return vcat([gulp_potential(pdic[a1_a2][1], a1_a2, pdic[a1_a2][2:end]...) for a1_a2 ∈ keys(pdic)]...)
end


function gulp_qpoints(config_gulp)::Vector{String}
    @assert (config_gulp[:qpoints] isa Vector) && (eltype(config_gulp[:qpoints])<:Pair)
    n = sum(last.(config_gulp[:qpoints]))
    return [
        "dispersion 1 $n", 
        repeat(" ",4) * join(js.(first.(config_gulp[:qpoints]))," to ")
    ]
end


##* ==============================================


function gulp_fit_phonon_input(config)::Vector{String}

    lines = [   gulp_first_line(config, :fit); 
                gulp_title(config); 
                gulp_structure_info(config);
                gulp_specific(config);
                gulp_observables(config);
                gulp_potential_pairs(config);
                gulp_qpoints(config);
            ]

end


function gulp_phonon_input(config)::Vector{String}

    lines = [   gulp_first_line(config, :phonon); 
                gulp_title(config); 
                gulp_structure_info(config);
                gulp_specific(config);
                # gulp_observables(config); #! no fitting
                gulp_potential_pairs(config);
                gulp_qpoints(config);
            ]

end

##* ==============================================
