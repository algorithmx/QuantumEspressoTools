function vasp_energy(result_lines00)
    en_begin = find_last_line(result_lines00, "Free energy of the ion-electron system (eV)")
    result_lines0 = result_lines00[en_begin:end]
    en_end   = find_line(result_lines0, "energy without entropy")
    result_lines  = result_lines0[1:en_end]
    en = Dict(
        "total_energy"                  => parsef(extract(result_lines, r"TOTEN\s*\=\s*" * num_f_rstr * r"\s*eV",           num_f_rstr)),
        "energy_without_entropy"        => parsef(extract(result_lines, r"energy\s+without\s+entropy\s*\=\s*" * num_f_rstr, num_f_rstr)),
        "alpha_Z"                       => parsef(extract(result_lines, r"alpha\s+Z\s*PSCENC\s*\=\s*" * num_f_rstr,         num_f_rstr)),
        "ewald_contribution"            => parsef(extract(result_lines, r"Ewald\s+ener?g?y?\s*TEWEN\s*=\s*"*num_f_rstr,     num_f_rstr)),
        "hartree_contribution"          => parsef(extract(result_lines, r"Hartree\s+ener?g?y?\s*DENC\s*=\s*"*num_f_rstr,    num_f_rstr)),
        "xc_contribution"               => parsef(extract(result_lines, r"V\(xc\)\+E\(xc\)\s*XCENC\s*=\s*"*num_f_rstr,      num_f_rstr)),
        "exchange"                      => parsef(extract(result_lines, r"exchange\s*EXHF\s*=\s*"*num_f_rstr,               num_f_rstr)),
        "PAW_double_counting"           => parse2f(extract(result_lines,r"PAW\s+double\s+counting\s*=\s*"*num_f_rstr*r"\s+"*num_f_rstr, num_f_rstr*r"\s+"*num_f_rstr)),
        "entropy"                       => parsef(extract(result_lines, r"entropy\s+T\*S\s*EENTRO\s*=\s*"*num_f_rstr,       num_f_rstr)),
        "bands"                         => parsef(extract(result_lines, r"eigenvalues\s*EBANDS\s*=\s*"*num_f_rstr,          num_f_rstr)),
        "atomic_energy"                 => parsef(extract(result_lines, r"atomic\s*ener?g?y?\s*EATOM\s*=\s*"*num_f_rstr,    num_f_rstr)),
        "Solvation_Ediel_sol"           => parsef(extract(result_lines, r"Solvation\s*Ediel\_sol\s*=\s*"*num_f_rstr,        num_f_rstr)),
    )
    return en
end


function vasp_pressure_volume(result_lines)
    return Dict(
        "pressure"    => parsef(extract_last(result_lines, r"external\s+pressure\s*\=\s*" * num_f_rstr * r"\s*kB",   num_f_rstr)),
        "cell_volume" => parsef(extract_last(result_lines, r"volume\s*of\s*cell\s*\:\s*"  * num_f_rstr,   num_f_rstr)),
    )
end


function vasp_contcar_to_minimal_cif(
    contcar_lines::Vector{S};
    axis=[1,2,3],
    autoswap=false,
    origin=[0,0,0],
    ) where {S<:AbstractString}
    title   = replace(contcar_lines[1]," "=>"")
    alat    = parsef(contcar_lines[2])
    @assert alat > 0
    (a0,b0,c0) = parse3f.((contcar_lines[3], contcar_lines[4], contcar_lines[5]))
    @assert norm(a0)>0 && norm(b0)>0 && norm(c0)>0
    lattp = basis_to_lattice_parameters((alat.*a0,alat.*b0,alat.*c0))
    (a,b,c,alpha,beta,gamma) = lattp.p
    @assert axis==[1,2,3] || (alpha≈90 && beta≈90 && gamma≈90)
    # autoswap
    axis1 = copy(axis)
    if autoswap && (alpha≈90 && beta≈90 && gamma≈90)
        axis1 = sortperm([a,b,c])
        @info "vasp_contcar_to_minimal_cif():\nUses autoswap for crystal with \n(a,b,c)=($a,$b,$c)"
    end
    v3 = [a0,b0,c0]
    # 乾坤大挪移 axis1
    lattp1 = basis_to_lattice_parameters((alat.*v3[axis1[1]][axis1],alat.*v3[axis1[2]][axis1],alat.*v3[axis1[3]][axis1]))
    latt_params = lattp1.p
    atoms   = SPLTS(contcar_lines[6])
    N_atoms = strtonum.(SPLTS(contcar_lines[7]))
    @assert all(N_atoms .> 0)
    @assert length(N_atoms)==length(atoms)
    @assert strip(contcar_lines[8])=="Direct"
    Tot_atom = sum(N_atoms)
    atom_pos = parse3f.(contcar_lines[8+1:8+Tot_atom])
    shift_frac(f3,orig) = [mod(f3[1]+orig[1]+2,1), mod(f3[2]+orig[2]+2,1), mod(f3[3]+orig[3]+2,1)]
    atom_frac_positions = vcat([    [ (atoms[ni], shift_frac(atom_pos[j+sum(N_atoms[1:ni-1])][axis1],origin)...) 
                                        for j=1:N_atoms[ni] ]
                                for ni=1:length(N_atoms) ]... )
    return minimal_cif( title, latt_params, atom_frac_positions )
end

