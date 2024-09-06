# pwx_output_relax.jl

##+ =======================

#* search for "Begin final coordinates"
#* after "End of BFGS Geometry Optimization"
#* do nothing if anything wrong

function pw_relax_final_section(result_lines::Vector{String}; quiet=false)
    
    t = find_line(result_lines, "End of BFGS Geometry Optimization", quiet=quiet)

    if t === nothing
        if !quiet  @error  "pw_relax_final_section():\nLine \"End of BFGS Geometry Optimization\" not found."  end
        return result_lines  # do nothing
    else
        result_lines1 = result_lines[(t+1):end]
        p = find_line(result_lines1, "Begin final coordinates")
        if p === nothing
            if !quiet  @error  "pw_relax_final_section():\nThe output file does not contain a line \"Begin final coordinates\"."  end
            return result_lines  # do nothing
        end

        result_lines2 = result_lines1[(p+1):end]
        q = find_line(result_lines2, "End final coordinates")
        if q === nothing
            if !quiet  @error  "pw_relax_final_section():\nThe output file does not contain a line \"End final coordinates\"."  end
            return result_lines  # do nothing
        end

        result_lines3 = result_lines2[1:q-1]
        if length(result_lines3)<=2
            if !quiet  @error   "pw_relax_final_section():\nThe output file contains nothing in the final result section."  end
            return result_lines  # do nothing
        end
    
        return result_lines3

    end
end

##+ =======================

function pw_relax_convergence_steps(res::Vector{String})
    n =  extract_last(  res, 
                        r"bfgs\s+converged\s+in\s+\d+\s+scf\s+cycles\s+and\s+\d+\s+bfgs\s+steps", 
                        r"\d+"  )
    return (n=="" ? -1 : strtonum(n))
end

##+ =======================

function pw_relax_final_enthalpy(result_lines::Vector{String}; quiet=true)
    en = extract(   pw_relax_final_section(result_lines; quiet=quiet), 
                    r"Final\s+enthalpy\s*\=\s*\-?\d+\.\d+\s*Ry", r"-?\d+\.\d+"  )
    return length(en)>0 ? parsef(en) : -10000.0
end

##+ =======================


function pw_relax_cell_parameters(result_lines::Vector{String}; quiet=false)
    basis = pw_crystal_axes_Angstrom(result_lines; quiet=quiet)
    lines_f = pw_relax_final_section(result_lines,quiet=true)
    t = find_last_line(lines_f, "CELL_PARAMETERS"; quiet=true)
    if t===nothing  
        # "relax" final result does not contain "CELL_PARAMETERS"
        # use input
        return basis
    else  # "vc-relax" 
        f3_rstr = r"\s*-?\d+\.\d+\s+-?\d+\.\d+\s+-?\d+\.\d+"
        lines   = delete_empty_lines(lines_f[t:end])
        if !(occursin("CELL_PARAMETERS", lines[1])
            && (f3_rstr ⊂ lines[2])
            && (f3_rstr ⊂ lines[3])
            && (f3_rstr ⊂ lines[4]))
            @error "pw_relax_cell_parameters():\nInvalid block of CELL_PARAMETERS.\n"*join(lines,"\n")
            return basis
        end
        # find alat
        alat    = 1.0
        if "bohr" ⊂ lines[1]
            alat = _BOHR_RADIUS_
        elseif "angstrom" ⊂ lines[1]
            alat = 1.0
        elseif "alat" ⊂ lines[1]
            alat = parsef(extract(lines, r"\(\s*alat\s*\=\s*\d+\.\d+\s*\)", num_f_rstr)) * _BOHR_RADIUS_  #???
        else
            @error "pw_relax_cell_parameters():\nInvalid block of CELL_PARAMETERS.\n"*join(lines,"\n")
            return basis
        end
        # return
        return alat .* parse3f.(lines[2:4])
    end
end


##+ ========== ATOMIC POSITIONS ===========

function pw_relax_atomic_positions(result_lines::Vector{String}; quiet=true)
    final_section         = pw_relax_final_section(result_lines, quiet=quiet)
    one_str_three_numbers = r"\s*\w+\s+-?\d+\.\d+\s+-?\d+\.\d+\s+-?\d+\.\d+"
    header_r              = r"\s*ATOMIC\_POSITIONS"

    p = find_last_line(final_section, header_r, quiet=quiet)
    @assert p!==nothing  "pw_relax_atomic_positions(): ATOMIC_POSITIONS section not found."
    unit = extract(final_section, header_r * r"\s+\(?\w+\)?$", r"\(?\w+\)?$")

    p += 1
    pos = []
    while p <= length(final_section)
        L =  strip(final_section[p])
        if L==""
            p += 1
        elseif one_str_three_numbers ⊂ L
            push!(pos, parse1s3f(L))
            p += 1
        else
            break
        end
    end
    return unit, pos
end

##+ =======================

function pw_relax_atom_list(BFGS_result_lines::Vector{String})
    (unit,atm)      = pw_relax_atomic_positions(BFGS_result_lines)
    @assert "crystal" ⊂ unit
    return [(a,p...) for (a,p) ∈ atm]
end

##+ =======================

pw_relax_to_lattice_parameters(BFGS_result_lines::Vector{String}) =  basis_to_lattice_parameters(pw_relax_cell_parameters(BFGS_result_lines)).p

##+ =======================

function findsym_pw_relax_result(
    SG_convention_used_to_get_result,
    BFGS_result_lines::Vector{String};
    findsym_settings = (latticeTolerance = 0.00001, atomicPositionTolerance = 0.00001, occupationTolerance = 0.0001),
    )
    IT_num   = 1
    latt_p   = pw_relax_to_lattice_parameters(BFGS_result_lines)
    atm_list = pw_relax_atom_list(BFGS_result_lines)
    try
        IT_num = findsym_input( 
                    "Called by findsym_pw_relax_result()",
                    latt_p,
                    [al[1:4] for al ∈ atm_list];
                    SG_setting = convert_setting_to_findsym(SG_convention_used_to_get_result),
                    latticeTolerance=findsym_settings[:latticeTolerance],
                    atomicPositionTolerance=findsym_settings[:atomicPositionTolerance],
                    occupationTolerance=findsym_settings[:occupationTolerance], 
        ) |> SPLTN |> findsym |> extract_space_group
    catch __e__
        @error "findsym_pw_relax_result() : error $(__e__). Use IT_num=1."
    end

    return IT_num
end

##+ =======================

function pw_relax_result_to_cif_lines(
    SG_convention_used_to_get_result,
    BFGS_result_lines::Vector{String};
    findsym_settings = (latticeTolerance = 0.00001, atomicPositionTolerance = 0.00001, occupationTolerance = 0.0001),
    )::Vector{String}
    cif_lines = []
    latt_p   = pw_relax_to_lattice_parameters(BFGS_result_lines)
    atm_list = pw_relax_atom_list(BFGS_result_lines)
    try
        cif_lines = findsym_input( 
                    "Called by findsym_pw_relax_result()",
                    latt_p,
                    [al[1:4] for al ∈ atm_list];
                    SG_setting = convert_setting_to_findsym(SG_convention_used_to_get_result),
                    latticeTolerance=findsym_settings[:latticeTolerance],
                    atomicPositionTolerance=findsym_settings[:atomicPositionTolerance],
                    occupationTolerance=findsym_settings[:occupationTolerance], 
        ) |> SPLTN |> findsym |> extract_cif |> STRPRM
    catch __e__
        @error "pw_relax_result_to_cif_file() : error $(__e__)."
    end
    return cif_lines
end


pw_relax_result_to_cif_lines(BFGS_result_lines::Vector{String}) = pw_relax_result_to_cif_lines(QE_default_symmetry_group_convention, BFGS_result_lines)

##+ =======================

function pw_relax_result_to_cell(
    SG_convention,
    BFGS_result_lines::Vector{String};
    )::UNITCELL

end


function pw_prepare_next_from_relax_result(
    SG_convention,
    BFGS_result_lines::Vector{String};
    findsym_settings = (latticeTolerance = 0.00001, atomicPositionTolerance = 0.00001, occupationTolerance = 0.0001),
    )

    latt_p   = pw_relax_to_lattice_parameters(BFGS_result_lines)
    atm_list = pw_relax_atom_list(BFGS_result_lines)
    IT_num   = findsym_pw_relax_result(SG_convention, BFGS_result_lines, findsym_settings=findsym_settings)
    return Dict( "positions"          => atm_list,
                 :cif                 => (latt_p..., 1), 
                 :do_not_use_symmetry => true  )
    @inline tovec(x) = [x...,]
    @inline rd(x) = round(x,digits=ntol)
    atm = tovec.(atom_list_ext)
    return Dict("positions"          => atm,
                :cif                 => (a, b, c, α, β, γ, 1),
                :cell_parameters     => [rd.(basis[1]),rd.(basis[2]),rd.(basis[3])],
                :do_not_use_symmetry => true)
end
