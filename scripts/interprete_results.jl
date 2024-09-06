global const __PROGRAM_DO_NOTHING__ = [
    "pp.x", 
    "bands.x", 
    "plotband.x", 
    "q2r.x", 
    "matdyn.x", 
    "ld1.x", 
    "wannier90.x", 
    "projwfc.x", 
    "pw2wannier90.x"
]


function interprete_results(p, output_lines)
    mpi, prog, mode, np, maby = analyze_prog(p)
    ret = Dict()
    if prog=="pwdry.x"
        try 
            ret = dict__pw_dryrun_result(output_lines)
        catch _e_
            @error _e_
        end
    elseif prog=="phdry.x"
        try 
            ret = dict__ph_dryrun_result(output_lines)
        catch _e_
            @error _e_
        end
    elseif prog=="pw.x"
        try 
            ret = dict__pw_result(output_lines)
        catch _e_
            @error _e_
        end
    elseif prog=="ph.x"
        try 
            ret = dict__ph_result(output_lines)
        catch _e_
            @error _e_
        end
    elseif prog ∈ __PROGRAM_DO_NOTHING__
        nothing
    else
        throw(error("interprete_results(): Unknown program $prog ."))
    end

    return ret
end


function dict__ph_dryrun_result(output_lines)
    return Dict()
end


function dict__pw_dryrun_result(output_lines)
    return Dict()
end


function dict__ph_result(ph_out)
    return Dict("frequency"=>ph_frequency(ph_out))
end


@inline correct_cell_parameters(basis_Ang_final) = 
    QuantumEspressoTools.lattice_parameters_to_basis(
        QuantumEspressoTools.basis_to_lattice_parameters(basis_Ang_final)
    )

@inline calculate_reciprocal_basis(avec_Ang) = (2π/norm(avec_Ang[1])).*inv(hcat(avec_Ang...))'

function dict__pw_relax_result(
    pw_out;
    findsym_settings = (latticeTolerance = 0.00001, 
                        atomicPositionTolerance = 0.00001, 
                        occupationTolerance = 0.0001),
    SG_convention    = QE_default_symmetry_group_convention,
    )
    # for completed task with a good output file
    config_seed  = Dict()
    try 
        n_bands      = pw_num_bands(pw_out)  #TODO spin?
        ecut         = pw_charge_cutoff_Ry(pw_out)
        wcut         = pw_wave_function_cutoff_Ry(pw_out)
        ps_files     = pw_pseudo_files(pw_out)
        ps_folder    = pw_pseudo_folder(pw_out)
        ibrav        = pw_ibrav(pw_out)

        #: -------------------------------------------
        IT           = findsym_pw_relax_result(SG_convention, pw_out, findsym_settings=findsym_settings)
        latt_params  = (pw_relax_to_lattice_parameters(pw_out)..., IT)
        atom_pos     = pw_relax_atom_list(pw_out)

        #: -------------------------------------------
        # CHANGE CELL_PARAMETERS -> calculate reciprocal basis
        avec_Ang_0   = pw_relax_cell_parameters(pw_out)
        avec_Ang     = correct_cell_parameters(avec_Ang_0)
        bvec_2pi_a0  = calculate_reciprocal_basis(avec_Ang)

        config_seed  = Dict( :cif               => latt_params, 
                            :cell_parameters    => avec_Ang,
                            :reciprocal_basis   => bvec_2pi_a0,
                            "positions"         => atom_pos,
                            "nbnd"              => n_bands,
                            "ecutwfc"           => wcut, 
                            "ecutrho"           => ecut,
                            "pseudo_files"      => ps_files,
                            "pseudo_dir"        => ps_folder,
                            "ibrav"             => ibrav,
        )
    catch _e_
        @error "dict__pw_relax_result() error : $(_e_)"
    end
    return config_seed
end


function dict__pw_result(
    pw_out;
    findsym_settings = (latticeTolerance = 0.00001, atomicPositionTolerance = 0.00001, occupationTolerance = 0.0001),
    SG_convention    = QE_default_symmetry_group_convention,
    )

    IT           = findsym_pw_result( SG_convention, pw_out, findsym_settings=findsym_settings )
    latt_params  = (pw_lattice_parameters(pw_out)..., IT)
    atom_pos     = pw_atom_list(pw_out)

    ibrav        = pw_ibrav(pw_out)
    alat_Bohr    = pw_alat_bohr(pw_out)
    avec_Bohr    = pw_crystal_axes_Bohr(pw_out)
    bvec_2pi_bohr= pw_reciprocal_axes(pw_out)
    bvec_2pi_a0  = [alat_Bohr.*bvec_2pi_bohr[1], alat_Bohr.*bvec_2pi_bohr[2], alat_Bohr.*bvec_2pi_bohr[3]]
    diff_iden    = sum(abs.([dot(avec_Bohr[i], bvec_2pi_bohr[j]) for i=1:3, j=1:3] .- diagm(0=>[1,1,1])))
    @assert  diff_iden < 2e-5   "diff_iden = $(diff_iden)" #! increase the fucking precision in the fucking output file !!!!!!!!

    n_bands      = pw_num_bands(pw_out)  #TODO spin?
    ecut         = pw_charge_cutoff_Ry(pw_out)
    e_fermi      = pw_fermi_energy_eV(pw_out)
    wcut         = pw_wave_function_cutoff_Ry(pw_out)
    ps_files     = pw_pseudo_files(pw_out)
    ps_folder    = pw_pseudo_folder(pw_out)
    bands        = pw_bands(pw_out)
    degauss      = pw_degauss(pw_out)
    config_seed  = Dict( :cif                => latt_params, 
                         :cell_parameters    =>  _BOHR_RADIUS_ .* avec_Bohr,
                         :reciprocal_basis   => bvec_2pi_a0,
                         "positions"         => atom_pos,
                         "nbnd"              => n_bands,
                         "ecutwfc"           => wcut, 
                         "ecutrho"           => ecut,
                         :E_fermi            => e_fermi,
                         "pseudo_files"      => ps_files,
                         "pseudo_dir"        => ps_folder,
                         "ibrav"             => ibrav,
                         "bands"             => bands,
                         "degauss"           => degauss, 
    )
    return config_seed
end

