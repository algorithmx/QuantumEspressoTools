##* ==========
##*  updaters
##* ==========

band_wannier90_upd_nothing(x) = Dict()


function upd_copy_XXX(x, FD, title, XXX)
    isfile("$FD/$(title).$(XXX).xml")  && (return Dict())  #! DO NOT OVERWRITE
    # copy
    try
        cp("$FD/$(title).xml", "$FD/$(title).$(XXX).xml")
    catch _e_
        @error "upd_copy_$(XXX)() : \n$(_e_)\nFILE $FD/$(title).xml NOT COPIED."
    end
    return Dict()
end

upd_copy_scf( x, FD, title) = upd_copy_XXX(x, FD, title, "scf")

upd_copy_band(x, FD, title) = upd_copy_XXX(x, FD, title, "band")

upd_copy_nscf(x, FD, title) = upd_copy_XXX(x, FD, title, "nscf")



function band_wannier90_upd_nscf_projwfc(x)
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
    @info "wannier90_upd_nscf_projwfc() updates  ret = $(showdict(ret))"
    return ret
end 


##* ======================
##*  exteranl controllers
##* ======================

band_wannier90_common(title) = Dict(
    "prefix"           => title, 
    "num_iter"         => 32000,
    "dis_num_iter"     => 32000,
    "num_guide_cycles" => 1,
    :updater           => band_wannier90_upd_nothing,
    :watchdog_setting  => (
        woof_per_x_min=1/10, 
        max_watch=80000, 
        max_silence_woofs=800, 
        tail_length=20, 
        quiet=false
    ),
)


band_wannier90_pw2w90_common(title) = Dict(
    "prefix"           => title, 
    "seedname"         => title, 
    :updater           => band_wannier90_upd_nothing,
    :watchdog_setting  => (
        woof_per_x_min=1/10, 
        max_watch=80000, 
        max_silence_woofs=800, 
        tail_length=20, 
        quiet=false
    )
)


band_wannier90_scf_common(title, psmode) = Dict(
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


band_wannier90_scf_external() = Dict()
band_wannier90_nscf_external() = Dict()
band_wannier90_wannier90_external() = Dict()
band_wannier90_wannier90pp_external() = Dict()
band_wannier90_pw2wannier90_external() = Dict()
band_wannier90_projwfc_external() = Dict()

function band_wannier90_band_external0(
    kpoint_mode, 
    kpoint_setting
    )
    if kpoint_mode=="crystal_b" # kpoint_setting = high_symm_line
        # kLINE(rx, ry) = [(0.0, 0.0, 0.0) => 100,
        #                  (0.5, 0.0, 0.0) => 100,
        #                  (0.5, 0.5, 0.0) => 100,
        #                  (0.0, 0.0, 0.0) => 100,
        #                  (0.0, 0.5, 0.0) => 1,]
        return Dict(:kpoint_mode=>"crystal_b", :kpoints=>kpoint_setting)
    elseif kpoint_mode=="tpiba_b" # kpoint_setting = kpoint_line
        # kLINE(rx, ry) = [(0.0, 0.0, 0.0) => 40, 
        #                  (0.0, 0.5*bfrac(rx,ry), 0.0) => 1,]
        return Dict(:kpoint_mode=>"tpiba_b", :kpoints=>kpoint_setting)
    elseif kpoint_mode=="tpiba_c" # kpoint_setting = kpoint_rect
        #kRECT(rx, ry) = [(0.0, 0.0, 0.0) => 1, 
        #                 (1.0, 0.0, 0.0) => 40,
        #                 (0.0, bfrac(rx,ry), 0.0)=>40,]
        return Dict(:kpoint_mode=>"tpiba_c", :kpoints=>kpoint_setting)
    else
        @error "band_wannier90_band_external($(kpoint_mode), $(kpoint_setting)) : unknown kpoint_mode $(kpoint_mode)."
        return Dict()
    end
end


band_wannier90_band_external(kpoint_mode, kpoint_setting) = band_wannier90_band_external0(kpoint_mode, kpoint_setting)


#! ========
#!   MAIN
#! ========


##* ==========================================================


global const band_wannier90_additional_settings_default = Dict(
        "diag" => ["david",],
        "beta" => [0.8, 0.4],
        "cutoff_scales" => [1.0,],
        "atoms" => [],
        "wann_occ_thr" => (0.95,0.1),
        "downscale" => 20.0,
)


function band_wannier90(
    workspace0::String,
    common_title::String, 
    cif0_fn::String,
    ps_mode, 
    kgrid_scf::NTuple{6,Int64},
    (kpoint_mode_band, kpoint_band),
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
    additional_settings = band_wannier90_additional_settings_default
    )

    # management
    workspace = try_mkdir(workspace0, ".", "band_wannier90()")
    cd(workspace)
    WORKSPACE = pwd()
    ENV["LOGGER_FILE_FULL_PATH"] = "$WORKSPACE/band_wannier90.log"

    # customized instruction
    #* nscf_conv_thr_downscale can be crucial to the stability of calculations
    #* sometimes the nscf calculation requires looser conv_thr than scf
    @inline getf(X) = get(additional_settings, X, band_wannier90_additional_settings_default[X])
    nscf_conv_thr_downscale = getf("downscale")
    findwin_thr = getf("wann_occ_thr")
    beta = getf("beta")
    diag = getf("diag")
    cts  = getf("cutoff_scales")
    atms = getf("atoms")

    function band_wannier90_instructor(struct0, FD, psmode, kp, kp_nscf, ct, atms=atms)
        #% find window
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
        # to update the "seeds" (default settings)
        pw_comm     = band_wannier90_scf_common(common_title, psmode) ∪ struct0 ∪ Dict(
            "nbnd"            => nbnd,
            "ecutwfc"         => ct[1], 
            "ecutrho"         => ct[2],
            :beta             => beta,
            :diag_choices     => diag,
            :cutoff_scales    => cts,
        )
        w90_comm    = band_wannier90_common(common_title) ∪ Dict(
            "projections"     => projections_for_wannier90_prog,
            "unit_cell_cart"  => struct0[:cell_parameters],
            "atoms_frac"      => struct0["positions"],
            "mp_grid"         => (kp_nscf[1:3]),
            "kpoints"         => kmesh_pl(kp_nscf[1:3]...,1)
        )
        proj_comm   = Dict()
        pw2w90_comm = band_wannier90_pw2w90_common(common_title)

        # INSTRUCTION
        SCF_OUTDIR  = "$WORKSPACE/$FD/$(common_title)_scf"
        INSTR    = INSTRUCTION(
            # DEFAULT PART
            Dict(
                PROG_PWX    => (_electron_scf_seed_conf_     ⬱  pw_comm     ),
                PROG_PW2W90 => (_electron_pw2w90_seed_conf_  ⬱  pw2w90_comm ),
                PROG_W90PP  => (_wannier90_seed_conf_        ⬱  w90_comm    ),
                PROG_W90    => (_wannier90_seed_conf_        ⬱  w90_comm    ),
                PROG_PROJ   => (_projwfc_seed_conf_          ⬱  proj_comm   )
            ),
            # STEPS
            [   (   "$(common_title)_scf",
                    PROG_PWX, 
                    Dict(   "outdir"       => SCF_OUTDIR,
                            :calc          => "scf",
                            :pw_mode       => "energy",
                            :kpoint_mode   => "automatic", 
                            :kpoints       => kp,
                            :updater       => x->upd_copy_scf(x,SCF_OUTDIR,common_title),
                    ) ∪ band_wannier90_scf_external()
                ),
                (   "$(common_title)_band",
                    PROG_PWX, 
                    Dict(   "outdir"       => SCF_OUTDIR,
                            :calc          => "bands",
                            :pw_mode       => "bands",
                            :beta          => 0.8.*beta,
                            :diag_choices  => diag,
                            "conv_thr"     => nscf_conv_thr_downscale*get(pw_comm, "conv_thr", 1e-11),
                            "startingpot"  => "file", 
                            :updater       => x->upd_copy_band(x,SCF_OUTDIR,common_title) ∪ band_wannier90_upd_nscf_projwfc(x),
                    ) ∪ band_wannier90_band_external(
                            kpoint_mode_band, 
                            kpoint_band
                    )
                ),
                (   "$(common_title)_band_projwfc", 
                    PROG_PROJ,
                    Dict(   "prefix"       => common_title,
                            "outdir"       => SCF_OUTDIR,
                            "filproj"      => common_title * ".proj",
                            "Emax"         => 10.0,    #* to be updated 
                            "Emin"         => -10.0,   #* to be updated 
                            "DeltaE"       => 0.005,   #* to be updated 
                            #"ngauss"       => 1,       #* to be updated 
                            #"degauss"      => 0.02,    #* to be updated 
                            :updater       => band_wannier90_upd_nothing,
                    )  ∪ band_wannier90_projwfc_external()
                ),
                (   "$(common_title)_nscf",
                    PROG_PWX, 
                    Dict(   "outdir"       => SCF_OUTDIR,
                            :calc          => "nscf",
                            ##:pw_mode       => "bands",
                            :beta          => 0.8.*beta,
                            "conv_thr"     => nscf_conv_thr_downscale*get(pw_comm, "conv_thr", 1e-11),
                            "startingpot"  => "file",  # must set to "file"
                            #* if kpoints settings are different in scf and nscf, cannot use "file"
                            #* kpoint settings are determined by "occupation" in scf but free to choose in nscf
                            #* they have to be consistent in order to use startingwfc = 'file'
                            #"startingwfc"  => "file", 
                            :kpoint_mode   => "crystal",
                            :kpoints       => kmesh_pl(kp_nscf[1:3]...,0),
                            :updater       => x->upd_copy_nscf(x,SCF_OUTDIR,common_title) ∪ band_wannier90_upd_nscf_projwfc(x),
                    ) ∪ band_wannier90_nscf_external()
                ),
                (   "$(common_title)_projwfc", 
                    PROG_PROJ,
                    Dict(   "prefix"       => common_title,
                            "outdir"       => SCF_OUTDIR,
                            "filproj"      => common_title * ".proj",
                            "Emax"         => 10.0,    #* to be updated 
                            "Emin"         => -10.0,   #* to be updated 
                            "DeltaE"       => 0.005,   #* to be updated 
                            #"ngauss"       => 1,       #* to be updated 
                            #"degauss"      => 0.02,    #* to be updated 
                            :updater       => projwfc_updater_local,
                    )  ∪ band_wannier90_projwfc_external()
                ),
                (   "$(common_title)_wannier90",
                    PROG_W90PP, 
                    Dict(
                    ) ∪ (band_wannier90_wannier90pp_external() ↓ [
                            "kpoint_path", "bands_plot", "bands_plot_mode", "bands_num_points", "bands_plot_project"
                        ])
                ),
                (   "$(common_title)_wannier90",
                    PROG_PW2W90, 
                    Dict( "outdir" => SCF_OUTDIR, 
                    ) ∪ band_wannier90_pw2wannier90_external()
                ),
                (   "$(common_title)_wannier90",
                    PROG_W90, 
                    Dict(
                    ) ∪ band_wannier90_wannier90_external()
                ),
            ]
        )
        return INSTR
    end

    res = grid_test(
        WORKSPACE, 
        common_title,
        cif0_fn, 
        band_wannier90_instructor,
        _FDNM__ps_kp_kp_ct_, # settings_to_folder_name
        #% the grid (below)
        [[ps_mode,], [kgrid_scf,], [kgrid_nscf,], [cutoff,]],
        cleanup = cleanup
    )

    return res

end


##* ==========================================================


global const kgrid_wannier90_additional_settings_default = Dict(
        "diag" => ["david",],
        "beta" => [0.8, 0.4],
        "cutoff_scales" => [1.0,],
        "atoms" => [],
        "wann_occ_thr" => (0.95,0.1),
        "downscale" => 20.0,
)


function kgrid_wannier90(
    workspace0::String,
    common_title::String, 
    cif0_fn::String,
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
    additional_settings = kgrid_wannier90_additional_settings_default
    )

    # management
    workspace = try_mkdir(workspace0, ".", "kgrid_wannier90()")
    cd(workspace)
    WORKSPACE = pwd()
    ENV["LOGGER_FILE_FULL_PATH"] = "$WORKSPACE/kgrid_wannier90.log"

    # customized instruction
    #* nscf_conv_thr_downscale can be crucial to the stability of calculations
    #* sometimes the nscf calculation requires looser conv_thr than scf
    @inline getf(X) = get(additional_settings, X, kgrid_wannier90_additional_settings_default[X])
    nscf_conv_thr_downscale = getf("downscale")
    findwin_thr = getf("wann_occ_thr")
    beta = getf("beta")
    diag = getf("diag")
    cts  = getf("cutoff_scales")
    atms = getf("atoms")

    function kgrid_wannier90_instructor(struct0, FD, psmode, kp, kp_nscf, ct, atms=atms)
        #% find window
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
        # to update the "seeds" (default settings)
        pw_comm     = band_wannier90_scf_common(common_title, psmode) ∪ struct0 ∪ Dict(
            "nbnd"            => nbnd,
            "ecutwfc"         => ct[1], 
            "ecutrho"         => ct[2],
            :beta             => beta,
            :diag_choices     => diag,
            :cutoff_scales    => cts,
        )
        w90_comm    = band_wannier90_common(common_title) ∪ Dict(
            "projections"     => projections_for_wannier90_prog,
            "unit_cell_cart"  => struct0[:cell_parameters],
            "atoms_frac"      => struct0["positions"],
            "mp_grid"         => (kp_nscf[1:3]),
            "kpoints"         => kmesh_pl(kp_nscf[1:3]...,1)
        )
        proj_comm   = Dict()
        pw2w90_comm = band_wannier90_pw2w90_common(common_title)

        # INSTRUCTION
        SCF_OUTDIR  = "$WORKSPACE/$FD/$(common_title)_scf"
        INSTR    = INSTRUCTION(
            # DEFAULT PART
            Dict(
                PROG_PWX    => (_electron_scf_seed_conf_     ⬱  pw_comm     ),
                PROG_PW2W90 => (_electron_pw2w90_seed_conf_  ⬱  pw2w90_comm ),
                PROG_W90PP  => (_wannier90_seed_conf_        ⬱  w90_comm    ),
                PROG_W90    => (_wannier90_seed_conf_        ⬱  w90_comm    ),
                PROG_PROJ   => (_projwfc_seed_conf_          ⬱  proj_comm   )
            ),
            # STEPS
            [   (   "$(common_title)_scf",
                    PROG_PWX, 
                    Dict(   "outdir"       => SCF_OUTDIR,
                            :calc          => "scf",
                            :pw_mode       => "energy",
                            :kpoint_mode   => "automatic", 
                            :kpoints       => kp,
                            :updater       => x->upd_copy_scf(x,SCF_OUTDIR,common_title),
                    ) ∪ band_wannier90_scf_external()
                ),
                #// (   "$(common_title)_band",
                #//     PROG_PWX, 
                #//     Dict(   "outdir"       => SCF_OUTDIR,
                #//             :calc          => "bands",
                #//             :pw_mode       => "bands",
                #//             :beta          => 0.8.*beta,
                #//             :diag_choices  => diag,
                #//             "conv_thr"     => nscf_conv_thr_downscale*get(pw_comm, "conv_thr", 1e-11),
                #//             "startingpot"  => "file", 
                #//             :updater       => x->upd_copy_band(x,SCF_OUTDIR,common_title) ∪ band_wannier90_upd_nscf_projwfc(x),
                #//     ) ∪ band_wannier90_band_external(
                #//             kpoint_mode_band, 
                #//             kpoint_band
                #//     )
                #// ),
                #// (   "$(common_title)_band_projwfc", 
                #//     PROG_PROJ,
                #//     Dict(   "prefix"       => common_title,
                #//             "outdir"       => SCF_OUTDIR,
                #//             "filproj"      => common_title * ".proj",
                #//             "Emax"         => 10.0,    #* to be updated 
                #//             "Emin"         => -10.0,   #* to be updated 
                #//             "DeltaE"       => 0.005,   #* to be updated 
                #//             #"ngauss"       => 1,       #* to be updated 
                #//             #"degauss"      => 0.02,    #* to be updated 
                #//             :updater       => band_wannier90_upd_nothing,
                #//     )  ∪ band_wannier90_projwfc_external()
                #// ),
                (   "$(common_title)_nscf",
                    PROG_PWX, 
                    Dict(   "outdir"       => SCF_OUTDIR,
                            :calc          => "nscf",
                            ##:pw_mode       => "bands",
                            :beta          => 0.8.*beta,
                            "conv_thr"     => nscf_conv_thr_downscale*get(pw_comm, "conv_thr", 1e-11),
                            "startingpot"  => "file",  # must set to "file"
                            #* if kpoints settings are different in scf and nscf, cannot use "file"
                            #* kpoint settings are determined by "occupation" in scf but free to choose in nscf
                            #* they have to be consistent in order to use startingwfc = 'file'
                            #"startingwfc"  => "file", 
                            :kpoint_mode   => "crystal",
                            :kpoints       => kmesh_pl(kp_nscf[1:3]...,0),
                            :updater       => x->upd_copy_nscf(x,SCF_OUTDIR,common_title) ∪ band_wannier90_upd_nscf_projwfc(x),
                    ) ∪ band_wannier90_nscf_external()
                ),
                (   "$(common_title)_projwfc", 
                    PROG_PROJ,
                    Dict(   "prefix"       => common_title,
                            "outdir"       => SCF_OUTDIR,
                            "filproj"      => common_title * ".proj",
                            "Emax"         => 10.0,    #* to be updated 
                            "Emin"         => -10.0,   #* to be updated 
                            "DeltaE"       => 0.005,   #* to be updated 
                            #"ngauss"       => 1,       #* to be updated 
                            #"degauss"      => 0.02,    #* to be updated 
                            :updater       => projwfc_updater_local,
                    )  ∪ band_wannier90_projwfc_external()
                ),
                (   "$(common_title)_wannier90",
                    PROG_W90PP, 
                    Dict(
                    ) ∪ band_wannier90_wannier90pp_external()
                ),
                (   "$(common_title)_wannier90",
                    PROG_PW2W90, 
                    Dict( "outdir" => SCF_OUTDIR, 
                    ) ∪ band_wannier90_pw2wannier90_external()
                ),
                (   "$(common_title)_wannier90",
                    PROG_W90, 
                    Dict(
                    ) ∪ band_wannier90_wannier90_external()
                ),
            ]
        )
        return INSTR
    end

    res = grid_test(
        WORKSPACE, 
        common_title,
        cif0_fn, 
        kgrid_wannier90_instructor,
        _FDNM__ps_kp_kp_ct_, # settings_to_folder_name
        #% the grid (below)
        [[ps_mode,], [kgrid_scf,], [kgrid_nscf,], [cutoff,]],
        cleanup = cleanup
    )

    return res

end


##* ==========================================================


global const band_NO_wannier90_additional_settings_default = Dict(
        "diag" => ["david",],
        "beta" => [0.8, 0.4],
        "cutoff_scales" => [1.0,],
        "atoms" => [],
        #//"wann_occ_thr" => (0.95,0.1),
        "downscale" => 20.0,
)


function band_NO_wannier90(
    workspace0::String,
    common_title::String, 
    cif0_fn::String,
    ps_mode, 
    kgrid_scf::NTuple{6,Int64},
    (kpoint_mode_band, kpoint_band),
    #//kgrid_nscf::NTuple{6,Int64},
    cutoff::Tuple{Float64,Float64},
    nbnd::Int;
    #//nwann::Int,
    #//projections_for_wannier90_prog::Vector{String},
    #//projections_for_projwfc::Vector{PROJWFC_PROJ_TYPE};
    cleanup     = false,
    PROG_PWX    = `mpiexec -np 64 --map-by=slot pw.x -npool 8`,
    #//PROG_PW2W90 = `pw2wannier90.x`,
    #//PROG_W90    = `wannier90.x`,
    #//PROG_W90PP  = `wannier90.x -pp`,
    PROG_PROJ   = `mpiexec -np 64 --map-by=slot projwfc.x`,
    additional_settings = band_NO_wannier90_additional_settings_default
    )

    # management
    workspace = try_mkdir(workspace0, ".", "band_NO_wannier90()")
    cd(workspace)
    WORKSPACE = pwd()
    ENV["LOGGER_FILE_FULL_PATH"] = "$WORKSPACE/band_NO_wannier90.log"

    # customized instruction
    #* nscf_conv_thr_downscale can be crucial to the stability of calculations
    #* sometimes the nscf calculation requires looser conv_thr than scf
    @inline getf(X) = get(additional_settings, X, band_NO_wannier90_additional_settings_default[X])
    nscf_conv_thr_downscale = getf("downscale")
    beta = getf("beta")
    diag = getf("diag")
    cts  = getf("cutoff_scales")
    atms = getf("atoms")

    function band_NO_wannier90_instructor(struct0, FD, psmode, kp, kp_nscf, ct, atms=atms)
        #% find window
        #//projwfc_updater_local(x) = (
        #//    ("inner_window" ∈ keys(additional_settings) && 
        #//     "outer_window" ∈ keys(additional_settings)) 
        #//    ? set_window0(
        #//        x,
        #//        nbnd,
        #//        nwann,
        #//        additional_settings["inner_window"],
        #//        additional_settings["outer_window"]  )
        #//    : find_window0(
        #//        x, 
        #//        nbnd,
        #//        nwann,
        #//        projections_for_projwfc,
        #//        readlines("$WORKSPACE/$FD/$(common_title)_nscf/$(common_title).pw.x.out"),
        #//        readlines("$WORKSPACE/$FD/$(common_title)_projwfc/$(common_title).proj.projwfc_up");
        #//        E_ref = pw_fermi_energy_eV(readlines("$WORKSPACE/$FD/$(common_title)_scf/$(common_title).pw.x.out")),
        #//        thr = findwin_thr )
        #//)
        # to update the "seeds" (default settings)
        pw_comm     = band_wannier90_scf_common(common_title, psmode) ∪ struct0 ∪ Dict(
            "nbnd"            => nbnd,
            "ecutwfc"         => ct[1], 
            "ecutrho"         => ct[2],
            :beta             => beta,
            :diag_choices     => diag,
            :cutoff_scales    => cts,
        )
        #//w90_comm    = band_wannier90_common(common_title) ∪ Dict(
        #//    "projections"     => projections_for_wannier90_prog,
        #//    "unit_cell_cart"  => struct0[:cell_parameters],
        #//    "atoms_frac"      => struct0["positions"],
        #//    "mp_grid"         => (kp_nscf[1:3]),
        #//    "kpoints"         => kmesh_pl(kp_nscf[1:3]...,1)
        #//)
        #//pw2w90_comm = band_wannier90_pw2w90_common(common_title)
        proj_comm   = Dict()

        # INSTRUCTION
        SCF_OUTDIR  = "$WORKSPACE/$FD/$(common_title)_scf"
        INSTR    = INSTRUCTION(
            # DEFAULT PART
            Dict(
                PROG_PWX    => (_electron_scf_seed_conf_     ⬱  pw_comm     ),
                #//PROG_PW2W90 => (_electron_pw2w90_seed_conf_  ⬱  pw2w90_comm ),
                #//PROG_W90PP  => (_wannier90_seed_conf_        ⬱  w90_comm    ),
                #//PROG_W90    => (_wannier90_seed_conf_        ⬱  w90_comm    ),
                PROG_PROJ   => (_projwfc_seed_conf_          ⬱  proj_comm   )
            ),
            # STEPS
            [   (   "$(common_title)_scf",
                    PROG_PWX, 
                    Dict(   "outdir"       => SCF_OUTDIR,
                            :calc          => "scf",
                            :pw_mode       => "energy",
                            :kpoint_mode   => "automatic", 
                            :kpoints       => kp,
                            :updater       => x->upd_copy_scf(x,SCF_OUTDIR,common_title),
                    ) ∪ band_wannier90_scf_external()
                ),
                (   "$(common_title)_band",
                    PROG_PWX, 
                    Dict(   "outdir"       => SCF_OUTDIR,
                            :calc          => "bands",
                            :pw_mode       => "bands",
                            :beta          => 0.8.*beta,
                            :diag_choices  => diag,
                            "conv_thr"     => nscf_conv_thr_downscale*get(pw_comm, "conv_thr", 1e-11),
                            "startingpot"  => "file", 
                            :updater       => x->upd_copy_band(x,SCF_OUTDIR,common_title) ∪ band_wannier90_upd_nscf_projwfc(x),
                    ) ∪ band_wannier90_band_external(
                            kpoint_mode_band, 
                            kpoint_band
                    )
                ),
                (   "$(common_title)_band_projwfc", 
                    PROG_PROJ,
                    Dict(   "prefix"       => common_title,
                            "outdir"       => SCF_OUTDIR,
                            "filproj"      => common_title * ".proj",
                            "Emax"         => 10.0,    #* to be updated 
                            "Emin"         => -10.0,   #* to be updated 
                            "DeltaE"       => 0.005,   #* to be updated 
                            #"ngauss"       => 1,       #* to be updated 
                            #"degauss"      => 0.02,    #* to be updated 
                            :updater       => band_wannier90_upd_nothing,
                    )  ∪ band_wannier90_projwfc_external()
                ),
                #// (   "$(common_title)_nscf",
                #//     PROG_PWX, 
                #//     Dict(   "outdir"       => SCF_OUTDIR,
                #//             :calc          => "nscf",
                #//             ##:pw_mode       => "bands",
                #//             :beta          => 0.8.*beta,
                #//             "conv_thr"     => nscf_conv_thr_downscale*get(pw_comm, "conv_thr", 1e-11),
                #//             "startingpot"  => "file",  # must set to "file"
                #//             #* if kpoints settings are different in scf and nscf, cannot use "file"
                #//             #* kpoint settings are determined by "occupation" in scf but free to choose in nscf
                #//             #* they have to be consistent in order to use startingwfc = 'file'
                #//             #"startingwfc"  => "file", 
                #//             :kpoint_mode   => "crystal",
                #//             :kpoints       => kmesh_pl(kp_nscf[1:3]...,0),
                #//             :updater       => x->upd_copy_nscf(x,SCF_OUTDIR,common_title) ∪ band_wannier90_upd_nscf_projwfc(x),
                #//     ) ∪ band_wannier90_nscf_external()
                #// ),
                #// (   "$(common_title)_projwfc", 
                #//     PROG_PROJ,
                #//     Dict(   "prefix"       => common_title,
                #//             "outdir"       => SCF_OUTDIR,
                #//             "filproj"      => common_title * ".proj",
                #//             "Emax"         => 10.0,    #* to be updated 
                #//             "Emin"         => -10.0,   #* to be updated 
                #//             "DeltaE"       => 0.005,   #* to be updated 
                #//             #"ngauss"       => 1,       #* to be updated 
                #//             #"degauss"      => 0.02,    #* to be updated 
                #//             :updater       => projwfc_updater_local,
                #//     )  ∪ band_wannier90_projwfc_external()
                #// ),
                #// (   "$(common_title)_wannier90",
                #//     PROG_W90PP, 
                #//     Dict(
                #//     ) ∪ (band_wannier90_wannier90pp_external() ↓ [
                #//             "kpoint_path", "bands_plot", "bands_plot_mode", "bands_num_points", "bands_plot_project"
                #//         ])
                #// ),
                #// (   "$(common_title)_wannier90",
                #//     PROG_PW2W90, 
                #//     Dict( "outdir" => SCF_OUTDIR, 
                #//     ) ∪ band_wannier90_pw2wannier90_external()
                #// ),
                #// (   "$(common_title)_wannier90",
                #//     PROG_W90, 
                #//     Dict(
                #//     ) ∪ band_wannier90_wannier90_external()
                #// ),
            ]
        )
        return INSTR
    end

    res = grid_test(
        WORKSPACE, 
        common_title,
        cif0_fn, 
        band_NO_wannier90_instructor,
        _FDNM__ps_kp_kp_ct_, # settings_to_folder_name
        #% the grid (below)
        [[ps_mode,], [kgrid_scf,], [kgrid_nscf,], [cutoff,]],
        cleanup = cleanup
    )

    return res

end


##* ==========================================================
