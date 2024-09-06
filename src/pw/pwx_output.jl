
#+ =======================

pw_convergence_achieved(d::Dict) = parse(Bool,d["output"]["convergence_info"]["scf_conv"]["convergence_achieved"])

pw_info(d::Dict) = d["general_info"]

#+ =======================

function scf_sections(result_lines::Vector{String})
    p = find_all_lines(result_lines, r"^\s+Self\-consistent\s+Calculation\s*")
    pb = [1, p...]
    pe = [(p.-1)..., length(result_lines)]
    return Vector{String}[result_lines[a:b] for (a,b) in zip(pb,pe)]
end

function final_scf_section(lines0::Vector{String}; quiet=true)  # ... if any
    p = find_line(lines0, r"final\s+scf\s+calculation", quiet=quiet)
    if p===nothing
        if !quiet  @warn "final_scf_section():\n Final scf section not found." end
        return lines0
    else
        return lines0[p:end]
    end
end

#+ =======================

pw_alat_bohr(d::Dict) = d["output"]["atomic_structure"][:alat]
function pw_alat_bohr(result_lines::Vector{String})
    a = extract_last(   result_lines,
                        r"lattice\s+parameter\s+\(alat\)\s+\=\s+\d+(\.\d*)?\s+a\.u\.",
                        r"\d+(\.\d*)?"  )
    return (a=="" ? -1000.0 : strtonum(a))
end

pw_alat_Ang(d::Dict) = pw_alat_bohr(d)*_BOHR_RADIUS_
pw_alat_Ang(result_lines::Vector{String}) = pw_alat_bohr(result_lines)*_BOHR_RADIUS_

#+ =======================

function pw_ibrav(d::Dict)
    atm_strct = d["output"]["atomic_structure"]
    return (:bravais_index ∈ keys(atm_strct)) ? atm_strct[:bravais_index] : 0
end 
pw_ibrav(result_lines::Vector{String}) = strtonum(extract_last(result_lines, r"bravais\-lattice\s+index\s*\=\s*\d+", r"\d+"))

#+ =======================
#+ unit in Ry

pw_charge_cutoff_Ry(result_lines::Vector{String}) = parsef(extract_last(result_lines, r"charge\s+density\s+cutoff\s*\=\s*\d+\.\d+\s+Ry", r"\d+\.\d+"))
pw_charge_cutoff_Ry(d::Dict) = 2*d["input"]["basis"]["ecutrho"]

#+ =======================
#+ unit in Ry

pw_wave_function_cutoff_Ry(result_lines::Vector{String}) = parsef(extract_last(result_lines, r"kinetic\-energy\s+cutoff\s*\=\s*\d+\.\d+\s+Ry", r"\d+\.\d+"))
pw_wave_function_cutoff_Ry(d::Dict) = 2*d["input"]["basis"]["ecutwfc"]

#+ =======================

#TODO spin?
pw_num_bands(result_lines::Vector{String}) = strtonum(extract_last(result_lines, r"number\s*of\s*Kohn\-Sham\s*states\s*\=\s*\d+",r"\d+"))
pw_num_bands(d::Dict) = d["output"]["band_structure"]["nbnd"]

#+ =======================

function pw_bands(result_lines::Vector{String})
    @inline all_numbers(x) = !any(isnan.(parsenf(x)))
    #!
    ##@inline extrakt(x) = strtonum_fort.(SPLTS((SPLTX(SPLTEQ(x)[end],"(")[1])))
    @inline extrakt(x) = strtonum_fort.(getindex.(string(rstrip(split(x,r"[=\(]",keepempty=false)[2])), [1:7,8:14,15:21]))
    @inline skip_empty(j) = (result_lines[j]!="" || j==length(result_lines)) ? j : skip_empty(j+1)
    band_pos = find_all_lines(result_lines, "bands (ev):") |> sort
    bands = []
    for p ∈ band_pos
        kpoint = extrakt(result_lines[p])
        i = p+1
        i = skip_empty(i)
        bnd = []
        while all_numbers(result_lines[i])
            bnd = [bnd; parsenf(result_lines[i])]
            i += 1
            i = skip_empty(i)
        end
        push!(bands, kpoint=>bnd)
    end
    return bands
end

#+ =======================

pw_total_energy(d::Dict) = 2.0*d["output"]["total_energy"]["etot"]   #! unit Ry !!!!!!
pw_total_energy(result_lines::Vector{String}) = parsef(extract(result_lines, r"\!\s*total\s+energy\s*\=\s*" * num_f_rstr * r"\s+Ry", num_f_rstr))

#+ =======================

#! reciprocal axes: ( cart. coord. in units 2 pi / Bohr )
function pw_reciprocal_axes(d::Dict)
    alat = pw_alat_bohr(d)
    B    = d["output"]["basis_set"]["reciprocal_lattice"]
    return (1.0/alat) .* parse3f.([B["b1"], B["b2"], B["b3"]])
end
function pw_reciprocal_axes(result_lines::Vector{String})
    lines_f = final_scf_section(result_lines; quiet=true)
    alat    = pw_alat_bohr(lines_f)
    p       = find_last_line(lines_f, r"reciprocal\s+axes\s*\:")
    if p!==nothing
        # a(1) = (   0.442196  -0.255302   0.859816 )
        b1_rstr = r"b\(\d\)\s*\=\s*\(\s*"*num_f_rstr3*r"\s*\)"
        b = extract_all(lines_f[p+1:p+3], b1_rstr, num_f_rstr3)
        # reciprocal axes: (cart. coord. in units 2 pi/alat)
        return (1.0/alat) .* parse3f.(b)
    else
        @error "pw_reciprocal_axes():\ncrystal axes not found."
        return [zeros(Float64,3),zeros(Float64,3),zeros(Float64,3)]
    end
end

##+ =======================


function pw_crystal_axes_Bohr(d::Dict) 
    A = d["output"]["atomic_structure"]["cell"]
    return parse3f.([A["a1"], A["a2"], A["a3"]])
end


function pw_crystal_axes_Bohr(result_lines::Vector{String}; quiet=false)
    lines_f   = final_scf_section(result_lines; quiet=true)
    alat_Bohr = pw_alat_bohr(lines_f)  #!
    p         = find_last_line(lines_f, r"crystal\s+axes\s*\:"; quiet=quiet)
    if p!==nothing
        # a(1) = (   0.442196  -0.255302   0.859816 )
        a1_rstr = r"a\(\d\)\s*\=\s*\(\s*"*num_f_rstr3*r"\s*\)"
        a = extract_all(lines_f[p+1:p+3],  a1_rstr,  num_f_rstr3)
        return alat_Bohr .* parse3f.(a)
    else
        @error "pw_crystal_axes_Angstrom():\ncrystal axes not found."
        return [zeros(Float64,3),zeros(Float64,3),zeros(Float64,3)]
    end
end


pw_crystal_axes_Angstrom(x; quiet=false) = (_BOHR_RADIUS_ .* pw_crystal_axes_Bohr(x; quiet=quiet))


pw_lattice_parameters(x) = basis_to_lattice_parameters(pw_crystal_axes_Angstrom(x)).p

#+ =======================


pw_num_kpoints(result_lines::Vector{String}) = extract_last(result_lines, r"\s*number\s+of\s+k\s+points\s*\=\s*\d+", r"\d+")
pw_num_kpoints(d::Dict) = d["output"]["band_structure"]["nks"]

#+ =======================
#IMPORTANT  prefer Cartesian coordinates, because we have to 
#IMPORTANT  rescale crystal axes in the relax results 
#IMPORTANT  before return fractional coordinates.

function pw_atom_list_XXX(d::Dict, XXX="input")::Vector{Tuple{String,Float64,Float64,Float64}}
    Ainv = inv(hcat(pw_crystal_axes_Angstrom(d)...))
    atom_list_sorted_by_index = sort(d[XXX]["atomic_structure"]["atomic_positions"]["atom"], by=x->x[:index])
    return [(x[:name], (Ainv*parse3f(x[""]))...) for x in atom_list_sorted_by_index]
end
pw_atom_list_intput(d::Dict) = pw_atom_list_XXX(d, "input")
pw_atom_list_output(d::Dict) = pw_atom_list_XXX(d, "output")
pw_atom_list(d::Dict) = pw_atom_list_output(d::Dict)


#* extract tau() directly 
function pw_extract_tau(result_lines0::Vector{String}, cart_or_frac)
    #NOTE  simply extract numbers inside the brackets
    # 1           In  tau(   1) = (   0.0000000   0.0000000  -0.0000000  )
    result_lines = final_scf_section(result_lines0, quiet=true)
    p1 = find_line(result_lines, cart_or_frac, quiet=true)
    if p1!==nothing
        pstart = p1     + delta_find_until_not(result_lines[p1:end],     x->!occursin(r"\s*\d+\s+\w+\s+tau.?",x))
        pend   = pstart + delta_find_until_not(result_lines[pstart:end], x-> occursin(r"\s*\d+\s+\w+\s+tau.?",x))
        pos_rstr = r"\w+\s+tau\(\s*\d+\s*\)\s*\=\s*\(\s*" * num_f_rstr3 * r"\s*\)"
        parse_pair(a,b) = (replace(a,r"\s*tau\s*\(\s*"=>"(") => parse3f(replace(b,r"[\(\)]"=>"")))
        line2pair(l) = parse_pair(strip.(split(l,"=",keepempty=false))...)
        atm_lines_cart = extract_all(result_lines[pstart:pend], pos_rstr, pos_rstr)
        return Dict(line2pair.(atm_lines_cart))
    else
        return Dict()
    end
end


function pw_atom_list(
    result_lines0::Vector{String}
    )::Vector{Tuple{String,Float64,Float64,Float64}}

    result_lines = final_scf_section(result_lines0, quiet=true)
    alat_Bohr    = pw_alat_bohr(result_lines)
    cryst_axes_a0= (1.0/alat_Bohr) .* pw_crystal_axes_Bohr(result_lines)

    @inline rm_num(x)    = strip(replace(x,r"\(\s*\d+\s*\)"=>""))
    @inline rescale_frac(frac) = (1.0/(alat_Bohr*norm(cryst_axes_a0[1]))) .* frac  #TODO

    pos_cart     = pw_extract_tau(result_lines, "Cartesian axes")  # positions (alat units)
    pos_frac     = pw_extract_tau(result_lines, "Crystallographic axes")

    atom_list_frac = []
    if length(pos_cart)>0
        ainv = inv(hcat(cryst_axes_a0...))
        atom_list_frac = [(rm_num(k), (ainv*v)...) for (k,v) in pos_cart] |> sort
    elseif length(pos_frac) > 0
        atom_list_frac = [(rm_num(k), rescale_frac(v)...) for (k,v) in pos_frac] |> sort
    else
        @error "pw_atom_list():\nNo fract coordinates found."
    end

    return atom_list_frac
end

#+ =======================

function pw_mixing_beta(result_lines::Vector{String})
    beta = extract_last(result_lines, r"beta\s*\=\s*"*num_f_rstr, num_f_rstr)
    return parsef(beta)
end

function pw_diag_style(result_lines::Vector{String})
    david = extract_last(result_lines, r"Davidson\s+diagonalization\s+with\s+overlap", r"Davidson")
    CG = extract_last(result_lines, r"CG\s+style\s+diagonalization", r"CG")
    @assert david != "" || CG != ""
    return david=="" ? "cg" : "david"
end

#+ =======================

global const Ry_in_eV = 13.605662285137
pw_fermi_energy_eV(d::Dict) = 2Ry_in_eV*d["output"]["band_structure"]["fermi_energy"]
pw_fermi_energy_eV(result_lines::Vector{String}) = parsef(extract_last(result_lines, r"Fermi\s+energy\s+is\s+\-?\d+.\d+\s+ev", r"\-?\d+.\d+"))

pw_degauss(result_lines::Vector{String}) = parsef(extract_last(result_lines, r"smearing,\s*width\s*\(Ry\)", r"0.\d+"))

#+ =======================

pw_pseudo_folder(d::Dict) = d["input"]["control_variables"]["pseudo_dir"]
function pw_pseudo_folder(result_lines::Vector{String})
    pos = find_all_lines(result_lines, r"PseudoPot\.\s+\#\s+\d+\s+for\s+[A-Z][a-z]?\s+read\s+from\s+file")
    fds = [strip(dirname(f)) for f in result_lines[pos.+1]]
    @assert length(unique(fds))==1
    return first(fds)
end

#+ =======================

pw_pseudo_files(d::Dict) = Dict(p[:name] => basename(p["pseudo_file"]) for p in d["input"]["atomic_species"]["species"])
function pw_pseudo_files(result_lines::Vector{String})
    pos = find_all_lines(result_lines, r"PseudoPot\.\s+\#\s+\d+\s+for\s+[A-Z][a-z]?\s+read\s+from\s+file")
    _get_elm_(l) = first(SPLTS(last(SPLTX(l,"for"))))
    ps_dic = Dict(_get_elm_(result_lines[p]) => basename(result_lines[p+1]) for p in pos )
    return ps_dic
end

#+ =======================

function findsym_pw_result(
    SG_convention_used_to_get_result,
    result_lines::Vector{String};
    findsym_settings = (latticeTolerance = 0.00001, atomicPositionTolerance = 0.00001, occupationTolerance = 0.0001)
    )
    IT_num   = 1
    latt_p   = pw_lattice_parameters(result_lines)
    atm_list = pw_atom_list(result_lines)
    try
        IT_num = findsym_input( 
                    "Called by findsym_pw_result()",
                    latt_p,
                    [al[1:4] for al ∈ atm_list];
                    SG_setting = convert_setting_to_findsym(SG_convention_used_to_get_result),
                    latticeTolerance=findsym_settings[:latticeTolerance],
                    atomicPositionTolerance=findsym_settings[:atomicPositionTolerance],
                    occupationTolerance=findsym_settings[:occupationTolerance], 
        ) |> SPLTN |> findsym |> extract_space_group
    catch __e__
        @error "findsym_pw_result() : error $(__e__). Use IT_num=1."
    end

    return IT_num
end


#+ ======= ENERGY ========

function pw_energy(result_lines00::Vector{String})
    total_energy_line_rstr = r"\!\s*total\s+energy\s*\=\s*" * num_f_rstr * r"\s+Ry"
    total_energy_line1_rstr = r"\s*total\s+energy\s*\=\s*" * num_f_rstr * r"\s+Ry"
    Harris_Foulkes_estimate_line_rstr = r"\s*Harris-Foulkes\s+estimate\s*\=\s*" * num_f_rstr * r"\s+Ry"
    estimated_scf_accuracy_line_rstr = r"\s*estimated\s+scf\s+accuracy\s*\<\s*" * num_f_rstr * r"\s+Ry"
    one_electron_contribution_line_rstr = r"\s*one-electron\s+contribution\s*\=\s*" * num_f_rstr * r"\s+Ry"
    hartree_contribution_line_rstr = r"\s*hartree\s+contribution\s*\=\s*" * num_f_rstr * r"\s+Ry"
    xc_contribution_line_rstr = r"\s*xc\s+contribution\s*\=\s*" * num_f_rstr * r"\s+Ry"
    ewald_contribution_line_rstr = r"\s*ewald\s+contribution\s*\=\s*" * num_f_rstr * r"\s+Ry"

    result_lines0 = final_scf_section(result_lines00; quiet=true)
    endp = find_last_line(result_lines0, "End of self-consistent calculation")
    result_lines  = result_lines0[endp+1:end]
    en = Dict(
        "total_energy"                  => parsef(extract(result_lines, total_energy_line_rstr,              num_f_rstr)),
        "Harris_Foulkes_estimate"       => parsef(extract(result_lines, Harris_Foulkes_estimate_line_rstr,   num_f_rstr)),
        "estimated_scf_accuracy"        => parsef(extract(result_lines, estimated_scf_accuracy_line_rstr,    num_f_rstr)),
        "one_electron_contribution"     => parsef(extract(result_lines, one_electron_contribution_line_rstr, num_f_rstr)),
        "hartree_contribution"          => parsef(extract(result_lines, hartree_contribution_line_rstr,      num_f_rstr)),
        "xc_contribution"               => parsef(extract(result_lines, xc_contribution_line_rstr,           num_f_rstr)),
        "ewald_contribution"            => parsef(extract(result_lines, ewald_contribution_line_rstr,        num_f_rstr)),
    )
    return en
end

pw_total_energy_history(result_lines::Vector{String}) = parsef.(extract_all(result_lines, total_energy_line1_rstr, num_f_rstr))

##* ================================ TEST ================================

