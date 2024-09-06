global const phonon_freq_woof = (
    woof_per_x_min=1/3, 
    max_watch=60000, 
    max_silence_woofs=600, 
    tail_length=40, 
    quiet=false
)

global const phonon_freq_scf_woof = (
    woof_per_x_min=1/3, 
    max_watch=10000, 
    max_silence_woofs=300, 
    tail_length=40, 
    quiet=false
)

##* ==========
##*  updaters
##* ==========

#TODO
function phonon_freq_updater_after_prerelax(x)
    d = dict__pw_relax_result(x)  
    return (d ↓ ["ibrav", :cif, :reciprocal_basis]) ⬱ Dict("ibrav"=>0,:do_not_use_symmetry=>true)
end


function phonon_freq_updater_after_relax(x)
    d = dict__pw_relax_result(x)
    return (d ↓ ["ibrav", :cif, :reciprocal_basis]) ⬱ Dict("ibrav"=>0,:do_not_use_symmetry=>true)
end


# for scf calculation, lock the mixing beta and diag method
function phonon_freq_updater_scf_ph(x)
    d = dict__pw_result(x)
    beta = pw_mixing_beta(x)
    diag = pw_diag_style(x)
    return Dict(:reciprocal_basis=>d[:reciprocal_basis], :beta=>[beta,], :diag_choices=>[diag,])
end


##* ========
##*  common
##* ========

phonon_freq_common(title, psmode, kp) = Dict(
        "title"           => title, 
        "prefix"          => title, 
        :pseudo_mode      => psmode,
        :pw_mode          => "energy",
        :kpoint_mode      => "automatic", 
        :kpoints          => kp,
        :watchdog_setting => phonon_freq_scf_woof
)


phonon_freq_relax(kBar, dg) = Dict(
        "occupations"     => "smearing",
        "smearing"        => "mp",
        "degauss"         => dg,
        "cell_dofree"     => "all",
        "press"           => Float64(kBar),
        "press_conv_thr"  => 0.0001, 
        "cell_dynamics"   => "bfgs", 
        "bfgs_ndim"       => 3,
        "ion_dynamics"    => "bfgs", 
        "upscale"         => 100.0,
)

phonon_freq_phx(title::String, kp, qp) = Dict(
        :ph_mode          => :multiple,
        "prefix"          => title, 
        "title"           => title,
        :kpoints          => kp,
        :qpoints          => qp,
        :watchdog_setting => phonon_freq_woof
)

##* ======================
##*  exteranl controllers
##* ======================

phonon_freq_prerelax_external() = Dict()
phonon_freq_relax_external() = Dict()
phonon_freq_scf_external() = Dict()
phonon_freq_ph_external() = Dict()

#! ========
#!   MAIN
#! ========

global const phonon_freq_additional_settings = Dict(
        "diag" => ["david",],
        "beta" => [0.8, 0.5, 0.3],
        "cutoff_scales" => [1.0,],
        "atoms" => [],
)

function phonon_freq(
    workspace0::String, 
    common_title::String, 
    cif0_fn::String, 
    ps_modes::Vector,
    kpoint_configs,
    cutoff_list,
    qpoints, 
    kBar_list::Vector{Float64},
    degauss_list::Vector{Float64};
    cleanup            = false,
    tr2_ph             = 1e-14,
    conv_thresholds    = Dict("scf"=>1e-11),
    scf_cutoff_upscale = 1.1,
    PROG_PWX           = `mpiexec -np 48 --map-by=slot pw.x -npool 8`,
    PROG_PHX           = `mpiexec -np 48 --map-by=slot ph.x -npool 8`,
    additional_settings= phonon_freq_additional_settings
    )

    workspace = try_mkdir(workspace0, ".", "phonon_freq()")
    cd(workspace)
    WORKSPACE = pwd()
    ENV["LOGGER_FILE_FULL_PATH"] = "$WORKSPACE/phonon_freq.log"

    #> ---------------------------------
    # customized instruction
    diag = additional_settings["diag"]
    beta = additional_settings["beta"]
    atoms = additional_settings["atoms"]
    cutoff_scales = additional_settings["cutoff_scales"]
    conv_thr_scf = get(conv_thresholds, "scf",     1e-11)
    conv_thr_pre = get(conv_thresholds, "prerelax", 1e-7)
    conv_thr_rlx = get(conv_thresholds, "relax",    1e-9)

    function ph_instructor(i0, FD, psmode, kp, ct, kBar, dg, atms=atoms)
        scf_settings = Dict(
            :beta          => beta,
            :diag_choices  => diag,
            :cutoff_scales => cutoff_scales,
            "lspinorb"     => occursin("_FR",pseudo_mode_name(psmode))
        )
        ph_comm  = phonon_freq_common(common_title, psmode, kp) ∪ phonon_freq_relax(kBar, dg) ∪ scf_settings
        ph_phx   = phonon_freq_phx(common_title, kp, qpoints)
        INSTR    =  INSTRUCTION(
            # DEFAULT PART
            Dict(   PROG_PWX => ( _electron_scf_seed_conf_  ⬱ ph_comm),
                    PROG_PHX => ((_phonon_seed_conf_        ⬱ ph_comm) ⬱ ph_phx),
            ),
            [   ("$(common_title)_prerelax", PROG_PWX,
                                    Dict(   :calc             => "vc-relax", 
                                            "conv_thr"        => conv_thr_pre,
                                            "etot_conv_thr"   => 1e-5,
                                            "forc_conv_thr"   => 1e-4,
                                            "ecutwfc"         => 0.8*ct[1], 
                                            "ecutrho"         => 0.8*ct[2],
                                            "disk_io"         => "none",
                                            "lspinorb"        => false,
                                            "noncolin"        => false,
                                            :updater          => phonon_freq_updater_after_prerelax,
                                    )  ∪ i0 ∪ phonon_freq_prerelax_external()),
                ("$(common_title)_relax", PROG_PWX,
                                    Dict(   :calc             => "vc-relax", 
                                            "conv_thr"        => conv_thr_rlx,
                                            "etot_conv_thr"   => 1e-7,
                                            "forc_conv_thr"   => 1e-6,
                                            "ecutwfc"         => ct[1], 
                                            "ecutrho"         => ct[2],
                                            "disk_io"         => "none",
                                            "lspinorb"        => false,
                                            "noncolin"        => false,
                                            :updater          => phonon_freq_updater_after_relax,
                                    )  ∪ phonon_freq_relax_external()),
                ("$(common_title)_scf",   PROG_PWX, 
                                    Dict(   :calc             => "scf",
                                            "conv_thr"        => conv_thr_scf,
                                            "ecutwfc"         => scf_cutoff_upscale*ct[1], 
                                            "ecutrho"         => scf_cutoff_upscale*ct[2],
                                            "disk_io"         => "medium",
                                            :updater          => phonon_freq_updater_scf_ph,
                                    )  ∪ phonon_freq_scf_external()),
                ("$(common_title)_ph",    PROG_PHX, 
                                    Dict(   "outdir"          => "$WORKSPACE/$FD/$(common_title)_scf",
                                            "tr2_ph"          => tr2_ph,
                                    )  ∪ phonon_freq_ph_external()),
            ]
        )
        return INSTR
    end

    #> ---------------------------------
    #> serial excution
    res = grid_test(    WORKSPACE, 
                        common_title,
                        cif0_fn, 
                        ph_instructor,
                        _FDNM__ps_kp_ct_kBar_dg_,  # settings_to_folder_name,
                        [ps_modes, kpoint_configs, cutoff_list, kBar_list, degauss_list];
                        cleanup = cleanup
    )

    return res
end


function phonon_freq__REPORT(
    workspace0::String, 
    common_title::String, 
    cif0_fn::String, 
    ps_modes::Vector,
    kpoint_configs,
    cutoff_list,
    qpoints, 
    kBar_list::Vector{Float64},
    degauss_list::Vector{Float64};
    tr2_ph             = 1e-14,
    conv_thresholds    = Dict("scf"=>1e-11),
    scf_cutoff_upscale = 1.1,
    PROG_PWX           = `mpiexec -np 48 --map-by=slot pw.x -npool 8`,
    PROG_PHX           = `mpiexec -np 48 --map-by=slot ph.x -npool 8`,
    additional_settings= phonon_freq_additional_settings
    )

    workspace = try_mkdir(workspace0, ".", "phonon_freq()")
    cd(workspace)
    WORKSPACE = pwd()
    ENV["LOGGER_FILE_FULL_PATH"] = "$WORKSPACE/phonon_freq.log"

    #> ---------------------------------
    # customized instruction
    diag = additional_settings["diag"]
    beta = additional_settings["beta"]
    atoms = additional_settings["atoms"]
    cutoff_scales = additional_settings["cutoff_scales"]
    conv_thr_scf = get(conv_thresholds, "scf",     1e-11)
    conv_thr_pre = get(conv_thresholds, "prerelax", 1e-7)
    conv_thr_rlx = get(conv_thresholds, "relax",    1e-9)

    function ph_instructor(i0, FD, psmode, kp, ct, kBar, dg, atms=atoms)
        scf_settings = Dict(
            :beta          => beta,
            :diag_choices  => diag,
            :cutoff_scales => cutoff_scales,
            "lspinorb"     => occursin("_FR",pseudo_mode_name(psmode))
        )
        ph_comm  = phonon_freq_common(common_title, psmode, kp) ∪ phonon_freq_relax(kBar, dg) ∪ scf_settings
        ph_phx   = phonon_freq_phx(common_title, kp, qpoints)
        INSTR    =  INSTRUCTION(
            # DEFAULT PART
            Dict(   PROG_PWX => ( _electron_scf_seed_conf_  ⬱ ph_comm),
                    PROG_PHX => ((_phonon_seed_conf_        ⬱ ph_comm) ⬱ ph_phx),
            ),
            [   ("$(common_title)_prerelax", PROG_PWX,
                                    Dict(   :calc             => "vc-relax", 
                                            "conv_thr"        => conv_thr_pre,
                                            "etot_conv_thr"   => 1e-5,
                                            "forc_conv_thr"   => 1e-4,
                                            "ecutwfc"         => 0.8*ct[1], 
                                            "ecutrho"         => 0.8*ct[2],
                                            "disk_io"         => "none",
                                            "lspinorb"        => false,
                                            "noncolin"        => false,
                                            :updater          => phonon_freq_updater_after_prerelax,
                                    )  ∪ i0 ∪ phonon_freq_prerelax_external()),
                ("$(common_title)_relax", PROG_PWX,
                                    Dict(   :calc             => "vc-relax", 
                                            "conv_thr"        => conv_thr_rlx,
                                            "etot_conv_thr"   => 1e-7,
                                            "forc_conv_thr"   => 1e-6,
                                            "ecutwfc"         => ct[1], 
                                            "ecutrho"         => ct[2],
                                            "disk_io"         => "none",
                                            "lspinorb"        => false,
                                            "noncolin"        => false,
                                            :updater          => phonon_freq_updater_after_relax,
                                    )  ∪ phonon_freq_relax_external()),
                ("$(common_title)_scf",   PROG_PWX, 
                                    Dict(   :calc             => "scf",
                                            "conv_thr"        => conv_thr_scf,
                                            "ecutwfc"         => scf_cutoff_upscale*ct[1], 
                                            "ecutrho"         => scf_cutoff_upscale*ct[2],
                                            "disk_io"         => "medium",
                                            :updater          => phonon_freq_updater_scf_ph,
                                    )  ∪ phonon_freq_scf_external()),
                ("$(common_title)_ph",    PROG_PHX, 
                                    Dict(   "outdir"          => "$WORKSPACE/$FD/$(common_title)_scf",
                                            "tr2_ph"          => tr2_ph,
                                    )  ∪ phonon_freq_ph_external()),
            ]
        )
        return INSTR
    end

    #> ---------------------------------
    #> serial excution
    res = grid_test__REPORT(
                WORKSPACE, 
                common_title,
                cif0_fn, 
                ph_instructor,
                _FDNM__ps_kp_ct_kBar_dg_, # settings_to_folder_name,
                [ps_modes, kpoint_configs, cutoff_list, kBar_list, degauss_list];
                status="completed",
                additional_info = x -> check_ph_negative_freq(last(x))
    )

    return res
end

