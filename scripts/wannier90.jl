global const wannier90_woof = (
    woof_per_x_min=1/3, 
    max_watch=10000, 
    max_silence_woofs=1200, 
    tail_length=40, 
    quiet=false
)

global const w90scf_woof = (
    woof_per_x_min=1/10, 
    max_watch=80000, 
    max_silence_woofs=800, 
    tail_length=20, 
    quiet=false
)


##* ==========
##*  updaters
##* ==========

wannier90_upd_nothing(x) = Dict()

wannier90_upd_scf_nscf(x) = Dict()

function wannier90_upd_nscf_projwfc(x)
    d = dict__pw_result(x)
    bands = hcat(last.(values.(d["bands"])) ...)'
    bmin, bmax = (minimum(bands), maximum(bands))
    ret =  Dict(
            "Emax"         => bmax+0.1,  
            "Emin"         => bmin-0.01, 
            "DeltaE"       => 0.005,
            #"ngauss"       => 1,   #TODO => d["smearing"]
            #"degauss"      => 0.02, #d["degauss"]
        )
    @info "\t\t\twannier90_upd_nscf_projwfc() updates, ret = \n\t\t\t$(ret)"
    return ret
end 


##* ======================
##*  exteranl controllers
##* ======================

wannier90_common(title) = Dict(
    "prefix"           => title, 
    "num_iter"         => 32000,
    "dis_num_iter"     => 32000,
    "num_guide_cycles" => 1,
    :updater           => wannier90_upd_nothing,
    :watchdog_setting  => (
        woof_per_x_min=1/10, 
        max_watch=80000, 
        max_silence_woofs=800, 
        tail_length=20, 
        quiet=false
    ),
)


pw2w90_common(title) = Dict(
    "prefix"           => title, 
    "seedname"         => title, 
    :updater           => wannier90_upd_nothing,
    :watchdog_setting  => (
        woof_per_x_min=1/10, 
        max_watch=80000, 
        max_silence_woofs=800, 
        tail_length=20, 
        quiet=false
    )
)


wannier90_scf_common(title, psmode) = Dict(
    "title"           => title, 
    "prefix"          => title, 
    "conv_thr"        => 1e-12,
    "etot_conv_thr"   => 1e-7,
    "forc_conv_thr"   => 1e-6,
    "disk_io"         => "medium",
    :pseudo_mode      => psmode,
    :watchdog_setting => (
        woof_per_x_min=1/2, 
        max_watch=10000, 
        max_silence_woofs=100, 
        tail_length=40, 
        quiet=false
    )
)


wannier90_scf_external()          = Dict()
wannier90_nscf_external()         = Dict()
wannier90_wannier90_external()    = Dict()
wannier90_wannier90pp_external()  = Dict()
wannier90_pw2wannier90_external() = Dict()
wannier90_projwfc_external()  = Dict()


global const wannier90_additional_settings_default = Dict(
        "diag" => ["david",],
        "beta" => [0.8, 0.4],
        "cutoff_scales" => [1.0,],
        "atoms" => [],
        "wann_occ_thr" => (0.9,0.01),
        "downscale" => 20.0,
)


#! ========
#!   MAIN
#! ========


##* ==========================================================


function wannier90_kgrid(
    workspace::String,
    common_title::String, 
    cif0_fn::String,  # management
    ps_mode, 
    kgrid_scf_nscf::NTuple{6,Int64},
    cutoff::Tuple{Float64,Float64},
    nbnd::Int,
    nwann::Int,
    projections_for_wannier90_prog::Vector{String},
    projections_for_projwfc::Vector{PROJWFC_PROJ_TYPE};
    PROG_PWX    = `mpiexec -np 64 --map-by=slot pw.x -npool 8`,
    PROG_PW2W90 = `pw2wannier90.x`,
    PROG_W90    = `wannier90.x`,
    PROG_W90PP  = `wannier90.x -pp`,
    PROG_PROJ   = `mpiexec -np 64 --map-by=slot projwfc.x`,
    additional_settings = wannier90_additional_settings_default
    )
    wannier90_kgrid(
        workspace,
        common_title, 
        cif0_fn,
        ps_mode, 
        kgrid_scf_nscf,
        kgrid_scf_nscf, #! same setting
        cutoff,
        nbnd,
        nwann,
        projections_for_wannier90_prog,
        projections_for_projwfc;
        PROG_PWX    = PROG_PWX,
        PROG_PW2W90 = PROG_PW2W90,
        PROG_W90    = PROG_W90,
        PROG_W90PP  = PROG_W90PP,
        PROG_PROJ   = PROG_PROJ,
        additional_settings = additional_settings
    )
end


##* ==========================================================


function wannier90_kgrid(
    workspace0::String,
    common_title::String, 
    cif0_fn::String,  # management
    ps_mode, 
    kgrid_scf::NTuple{6,Int64},
    kgrid_nscf::NTuple{6,Int64},
    cutoff::Tuple{Float64,Float64},
    nbnd::Int,
    nwann::Int,
    projections_for_wannier90_prog::Vector{String},
    projections_for_projwfc::Vector{PROJWFC_PROJ_TYPE};
    cleanup     = false,
    PROG_PWX    = `mpiexec -np 64 --map-by=slot pw.x -npool 8`,
    PROG_PW2W90 = `pw2wannier90.x`,
    PROG_W90    = `wannier90.x`,
    PROG_W90PP  = `wannier90.x -pp`,
    PROG_PROJ   = `mpiexec -np 64 --map-by=slot projwfc.x`,
    additional_settings = wannier90_additional_settings_default
    )

    # management
    workspace = try_mkdir(workspace0, ".", "wannier90_kgrid()")
    cd(workspace)
    WORKSPACE = pwd()
    ENV["LOGGER_FILE_FULL_PATH"] = "$WORKSPACE/wannier90_kgrid.log"

    # customized instruction
    diag = additional_settings["diag"]
    beta = additional_settings["beta"]
    cutoff_scales = additional_settings["cutoff_scales"]
    atoms = additional_settings["atoms"]
    #* nscf_conv_thr_downscale can be crucial to the stability of calculations
    #* sometimes the nscf calculation requires looser conv_thr than scf
    nscf_conv_thr_downscale = get(additional_settings, "downscale", 20.0)
    findwin_thr = get(additional_settings,"wann_occ_thr", (0.9,0.01))

    function wannier90_instructor(struct0, FD, psmode, kp, kp_nscf, ct, atms=atoms)
        # specific parts of the configuration
        cutoffs  = Dict("ecutwfc" => ct[1], "ecutrho" => ct[2])
        spin_orb = Dict("lspinorb" => occursin("_FR",pseudo_mode_name(psmode)))
        scf_settings = Dict(
            "nbnd"         => nbnd,
            :beta          => beta,
            :diag_choices  => diag,
            :cutoff_scales => cutoff_scales,
        )
        proj_w90 = Dict("projections" =>  projections_for_wannier90_prog)
        ## struct0 is the result of get_structure_from_cif()
        cryst_struct_w09 = Dict(
            "unit_cell_cart" => struct0[:cell_parameters],
            "atoms_frac" => struct0["positions"],
        )
        kp_w90 = Dict(
            "mp_grid" => (kp_nscf[1:3]),
            "kpoints" => kmesh_pl(kp_nscf[1:3]...,1)
        )
        # to update the "seeds" (default settings)
        pw_comm     = w90scf_common(common_title, psmode) ∪ struct0 ∪ cutoffs ∪ spin_orb ∪ scf_settings
        wg_comm     = wannier90_common(common_title)      ∪ proj_w90 ∪ cryst_struct_w09 ∪ kp_w90
        proj_comm   = Dict()
        pw2w90_comm = pw2w90_common(common_title)
        # INSTRUCTION
        projwfc_updater_local(x) = (
            ("inner_window" ∈ keys(additional_settings) && 
             "outer_window" ∈ keys(additional_settings)) 
            ? set_window0(
                x,
                nbnd,
                nwann,
                additional_settings["inner_window"],
                additional_settings["outer_window"]  )
            : find_window0(
                x, 
                nbnd,
                nwann,
                projections_for_projwfc,
                readlines("$WORKSPACE/$FD/$(common_title)_nscf/$(common_title).pw.x.out"),
                readlines("$WORKSPACE/$FD/$(common_title)_projwfc/$(common_title).proj.projwfc_up");
                E_ref = pw_fermi_energy_eV(readlines("$WORKSPACE/$FD/$(common_title)_scf/$(common_title).pw.x.out")),
                thr = findwin_thr )
        )
        INSTR    = INSTRUCTION(
            # DEFAULT PART
            Dict(
                PROG_PWX    => (_electron_scf_seed_conf_     ⬱  pw_comm     ),
                PROG_PW2W90 => (_electron_pw2w90_seed_conf_  ⬱  pw2w90_comm ),
                PROG_W90PP  => (_wannier90_seed_conf_        ⬱  wg_comm     ),
                PROG_W90    => (_wannier90_seed_conf_        ⬱  wg_comm     ),
                PROG_PROJ   => (_projwfc_seed_conf_          ⬱  proj_comm   )
            ),
            # STEPS
            [   (   "$(common_title)_scf",
                    PROG_PWX, 
                    Dict(   "outdir"       => "$WORKSPACE/$FD/$(common_title)_scf",
                            :calc          => "scf",
                            :pw_mode       => "energy",
                            :kpoint_mode   => "automatic", 
                            :kpoints       => kp,
                            :updater       => wannier90_upd_scf_nscf,
                    ) ∪ wannier90_scf_external()
                ),
                (   "$(common_title)_nscf",
                    PROG_PWX, 
                    Dict(   "outdir"       => "$WORKSPACE/$FD/$(common_title)_scf",
                            :calc          => "nscf",
                            "conv_thr"     => nscf_conv_thr_downscale*get(pw_comm, "conv_thr", 1e-11),
                            #! must set to "file"
                            "startingpot"  => "file", 
                            #! if kpoints settings are different in scf and nscf, cannot use "file"
                            #! kpoint settings are determined by "occupation" in scf but free to choose in nscf
                            #! they have to be consistent in order to use startingwfc = 'file'
                            #"startingwfc"  => "file", 
                            :kpoint_mode   => "crystal",
                            :kpoints       => kmesh_pl(kp_nscf[1:3]...,0),
                            :updater       => wannier90_upd_nscf_projwfc,
                    ) ∪ wannier90_nscf_external()
                ),
                (   "$(common_title)_projwfc", 
                    PROG_PROJ,
                    Dict(   "prefix"       => common_title,
                            "outdir"       => "$WORKSPACE/$FD/$(common_title)_scf",
                            "filproj"      => common_title * ".proj",
                            "Emax"         => 10.0,    #* to be updated 
                            "Emin"         => -10.0,   #* to be updated 
                            "DeltaE"       => 0.005,   #* to be updated 
                            "ngauss"       => 1,       #* to be updated 
                            "degauss"      => 0.02,    #* to be updated 
                            :updater       => projwfc_updater_local,
                    )  ∪ wannier90_projwfc_external()
                ),
                (   "$(common_title)_wannier90",
                    PROG_W90PP, 
                    Dict(
                    ) ∪ wannier90_wannier90pp_external()
                ),
                (   "$(common_title)_wannier90",
                    PROG_PW2W90, 
                    Dict(   "outdir"       => "$WORKSPACE/$FD/$(common_title)_scf",
                    ) ∪ wannier90_pw2wannier90_external()
                ),
                (   "$(common_title)_wannier90",
                    PROG_W90, 
                    Dict(
                    ) ∪ wannier90_wannier90_external()
                ),
            ]
        )
        return INSTR
    end

    res = grid_test(
        WORKSPACE, 
        common_title,
        cif0_fn, 
        wannier90_instructor,
        _FDNM__ps_kp_kp_ct_, # settings_to_folder_name
        #% the grid (below)
        [[ps_mode,], [kgrid_scf,], [kgrid_nscf,], [cutoff,]],
        cleanup = cleanup
    )

    return res

end


##* ==========================================================


function wannier90_kgrid__REPORT(
    workspace::String,
    common_title::String, 
    cif0_fn::String,  # management
    ps_mode, 
    kgrid_scf_nscf::NTuple{6,Int64},
    cutoff::Tuple{Float64,Float64},
    nbnd::Int,
    nwann::Int,
    projections_for_wannier90_prog::Vector{String},
    projections_for_projwfc::Vector{PROJWFC_PROJ_TYPE};
    PROG_PWX    = `mpiexec -np 64 --map-by=slot pw.x -npool 8`,
    PROG_PW2W90 = `pw2wannier90.x`,
    PROG_W90    = `wannier90.x`,
    PROG_W90PP  = `wannier90.x -pp`,
    PROG_PROJ   = `mpiexec -np 64 --map-by=slot projwfc.x`,
    additional_settings = wannier90_additional_settings_default
    )
    wannier90_kgrid__REPORT(
        workspace,
        common_title, 
        cif0_fn,
        ps_mode, 
        kgrid_scf_nscf,
        kgrid_scf_nscf,
        cutoff,
        nbnd,
        nwann,
        projections_for_wannier90_prog,
        projections_for_projwfc;
        PROG_PWX    = PROG_PWX,
        PROG_PW2W90 = PROG_PW2W90,
        PROG_W90    = PROG_W90,
        PROG_W90PP  = PROG_W90PP,
        PROG_PROJ   = PROG_PROJ,
        additional_settings = additional_settings
    )
end

##* ==========================================================


function wannier90_kgrid__REPORT(
    workspace0::String,
    common_title::String, 
    cif0_fn::String,  # management
    ps_mode, 
    kgrid_scf::NTuple{6,Int64},
    kgrid_nscf::NTuple{6,Int64},
    cutoff::Tuple{Float64,Float64},
    nbnd::Int,
    nwann::Int,
    projections_for_wannier90_prog::Vector{String},
    projections_for_projwfc::Vector{PROJWFC_PROJ_TYPE};
    PROG_PWX    = `mpiexec -np 64 --map-by=slot pw.x -npool 8`,
    PROG_PW2W90 = `pw2wannier90.x`,
    PROG_W90    = `wannier90.x`,
    PROG_W90PP  = `wannier90.x -pp`,
    PROG_PROJ   = `mpiexec -np 64 --map-by=slot projwfc.x`,
    additional_settings = wannier90_additional_settings_default
    )

    # management
    workspace = try_mkdir(workspace0, ".", "wannier90_kgrid()")
    cd(workspace)
    WORKSPACE = pwd()
    ENV["LOGGER_FILE_FULL_PATH"] = "$WORKSPACE/wannier90_kgrid.log"

    # customized instruction
    diag = additional_settings["diag"]
    beta = additional_settings["beta"]
    cutoff_scales = additional_settings["cutoff_scales"]
    atoms = additional_settings["atoms"]
    findwin_thr = get(additional_settings,"wann_occ_thr",0.4)

    function wannier90_instructor(struct0, FD, psmode, kp, kp_nscf, ct, atms=atoms)
        # specific parts of the configuration
        cutoffs  = Dict("ecutwfc" => ct[1], "ecutrho" => ct[2])
        spin_orb = Dict("lspinorb" => occursin("_FR",pseudo_mode_name(psmode)))
        scf_settings = Dict(
            "nbnd"         => nbnd,
            :beta          => beta,
            :diag_choices  => diag,
            :cutoff_scales => cutoff_scales,
        )
        proj_w90 = Dict("projections" =>  projections_for_wannier90_prog)
        ## struct0 is the result of get_structure_from_cif()
        cryst_struct_w09 = Dict(
            "unit_cell_cart" => struct0[:cell_parameters],
            "atoms_frac" => struct0["positions"],
        )
        kp_w90 = Dict(
            "mp_grid" => (kp_nscf[1:3]),
            "kpoints" => kmesh_pl(kp_nscf[1:3]...,1)
        )
        # to update the "seeds" (default settings)
        pw_comm     = w90scf_common(common_title, psmode) ∪ struct0 ∪ cutoffs ∪ spin_orb ∪ scf_settings
        wg_comm     = wannier90_common(common_title)      ∪ proj_w90 ∪ cryst_struct_w09 ∪ kp_w90
        proj_comm   = Dict()
        pw2w90_comm = pw2w90_common(common_title)
        # INSTRUCTION
        projwfc_updater_local(x) = (
            ("inner_window" ∈ keys(additional_settings) && 
             "outer_window" ∈ keys(additional_settings)) 
            ? set_window0(
                x,
                nbnd,
                nwann,
                additional_settings["inner_window"],
                additional_settings["outer_window"]  )
            : find_window0(
                x, 
                nbnd,
                nwann,
                projections_for_projwfc,
                readlines("$WORKSPACE/$FD/$(common_title)_nscf/$(common_title).pw.x.out"),
                readlines("$WORKSPACE/$FD/$(common_title)_projwfc/$(common_title).proj.projwfc_up");
                E_ref = pw_fermi_energy_eV(readlines("$WORKSPACE/$FD/$(common_title)_nscf/$(common_title).pw.x.out")),
                thr = findwin_thr )
        )
        INSTR    = INSTRUCTION(
            # DEFAULT PART
            Dict(
                PROG_PWX    => (_electron_scf_seed_conf_     ⬱  pw_comm     ),
                PROG_PW2W90 => (_electron_pw2w90_seed_conf_  ⬱  pw2w90_comm ),
                PROG_W90PP  => (_wannier90_seed_conf_        ⬱  wg_comm     ),
                PROG_W90    => (_wannier90_seed_conf_        ⬱  wg_comm     ),
                PROG_PROJ   => (_projwfc_seed_conf_          ⬱  proj_comm   )
            ),
            # STEPS
            [   (   "$(common_title)_scf",
                    PROG_PWX, 
                    Dict(   "outdir"       => "$WORKSPACE/$FD/$(common_title)_scf",
                            :calc          => "scf",
                            :pw_mode       => "energy",
                            :kpoint_mode   => "automatic", 
                            :kpoints       => kp,
                            :updater       => wannier90_upd_scf_nscf,
                    ) ∪ wannier90_scf_external()
                ),
                (   "$(common_title)_nscf",
                    PROG_PWX, 
                    Dict(   "outdir"       => "$WORKSPACE/$FD/$(common_title)_scf",
                            :calc          => "nscf",
                            "conv_thr"     => 10.0*get(pw_comm, "conv_thr", 1e-11),
                            #! must set to "file"
                            "startingpot"  => "file", 
                            #! if kpoints settings are different in scf and nscf, cannot use "file"
                            #! kpoint settings are determined by "occupation" in scf but free to choose in nscf
                            #! they have to be consistent in order to use startingwfc = 'file'
                            #"startingwfc"  => "file", 
                            :kpoint_mode   => "crystal",
                            :kpoints       => kmesh_pl(kp_nscf[1:3]...,0),
                            :updater       => wannier90_upd_nscf_projwfc,
                    ) ∪ wannier90_nscf_external()
                ),
                (   "$(common_title)_projwfc", 
                    PROG_PROJ,
                    Dict(   "prefix"       => common_title,
                            "outdir"       => "$WORKSPACE/$FD/$(common_title)_scf",
                            "filproj"      => common_title * ".proj",
                            "Emax"         => 10.0,    #* to be updated 
                            "Emin"         => -10.0,   #* to be updated 
                            "DeltaE"       => 0.005,   #* to be updated 
                            "ngauss"       => 1,       #* to be updated 
                            "degauss"      => 0.02,    #* to be updated 
                            :updater       => projwfc_updater_local,
                    )  ∪ wannier90_project_wfc_external()
                ),
                (   "$(common_title)_wannier90",
                    PROG_W90PP, 
                    Dict(
                    ) ∪ wannier90_wannier90pp_external()
                ),
                (   "$(common_title)_wannier90",
                    PROG_PW2W90, 
                    Dict(   "outdir"       => "$WORKSPACE/$FD/$(common_title)_scf",
                    ) ∪ wannier90_pw2wannier90_external()
                ),
                (   "$(common_title)_wannier90",
                    PROG_W90, 
                    Dict(
                    ) ∪ wannier90_wannier90_external()
                ),
            ]
        )
        return INSTR
    end

    res = grid_test__REPORT(
        WORKSPACE, 
        common_title,
        cif0_fn, 
        wannier90_instructor,
        _FDNM__ps_kp_kp_ct_, # settings_to_folder_name
        #% the grid (below)
        [[ps_mode,], [kgrid_scf,], [kgrid_nscf,], [cutoff,]]
    )

    return res

end
