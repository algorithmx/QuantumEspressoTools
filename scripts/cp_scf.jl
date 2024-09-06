global const cp_scf_woof = (woof_per_x_min=1/6, max_watch=20000, max_silence_woofs=10, tail_length=20, quiet=false)

cp_scf_common(title, psmode, dt, EDYN, ECONF, nr, kBar, dg) = Dict(
    "title"           => title, 
    "prefix"           => title, 
    :pseudo_mode      => psmode, 

    :calc             => "scf", 
    "verbosity"       => "low", 
    "restart_mode"    => "from_scratch",
    "disk_io"         => "default",

    "nbnd"            => :default,
    "conv_thr"        => 1e-6,
    "etot_conv_thr"   => 1e-6,
    "ekin_conv_thr"   => 1e-5,
    "forc_conv_thr"   => 1e-4,

    "dt"              => dt,
    "degauss"         => dg,
    "press"           => Float64(kBar),

    #! electron
    "electron_temperature" => "not_controlled",
    "electron_dynamics"    => EDYN,
    "ekincw"               => ECONF[1], #! electron
    "emass"                => ECONF[2], #! electron
    "emass_cutoff"         => ECONF[3], #! electron

    #
    "nr1"                  => nr[1], 
    "nr2"                  => nr[2], 
    "nr3"                  => nr[3],
    #
    "nr1s"                 => nr[4], 
    "nr2s"                 => nr[5], 
    "nr3s"                 => nr[6],
    #
    "nr1b"                 => nr[7], 
    "nr2b"                 => nr[8], 
    "nr3b"                 => nr[9],

    :watchdog_setting      => cp_scf_woof
)

#: ==================================================================


function cp_scf(
    workspace0::String,
    common_title::String,
    cif0_fn::String,
    ps_modes::Vector,
    dt_list::Vector,
    elec_dynamics_list::Vector,
    elec_conf_list::Vector,
    nrXb_list::Vector,
    cutoff_list::Vector{Tuple{Float64,Float64}},
    kBar_list::Vector{Float64},
    degauss_list::Vector{Float64};
    # 
    PROG_CPX = `mpiexec -np 64 --map-by=slot cp.x -npool 8`,
    atoms = [],
    # :pilot_rules, example AUTOPILOT = [(31, "dt", 5.0), (91, "iprint", 100), ...]
    AUTOPILOT = [],
    NSTEP     = 2500
    )

    workspace = try_mkdir(workspace0, ".", "cp_scf()")
    cd(workspace)
    WORKSPACE = pwd()
    ENV["LOGGER_FILE_FULL_PATH"] = "$WORKSPACE/cp_scf.log"

    #> ---------------------------------
    #> customized instruction
    function instructor(init_struct, FD, psmode, dt, EDYN, ECONF, nrXb, cutoffs, kBar, dg, atms=atoms)
        wcut, ecut = cutoffs
        cp_scf_comm  = cp_scf_common(common_title, psmode, dt, EDYN, ECONF, nrXb, kBar, dg)
        INSTR = INSTRUCTION(
            Dict(PROG_CPX => (_cp_scf_seed_conf_  ⬱ cp_scf_comm)),
            [
                (   "$(common_title)_cp_scf", 
                    PROG_CPX,
                    Dict(   :calc             => "scf", 
                            "nstep"           => NSTEP,
                            :cutoff_scales    => [1.0, ], 
                            "ecutwfc"         => wcut,
                            "ecutrho"         => ecut,
                            :updater          => x->Dict(),
                            :pilot_rules      => AUTOPILOT
                    ) ∪ init_struct ) 
            ]
        )
        return INSTR
    end

    #> ---------------------------------
    #> folders and settings
    cp_nt(ps,dt,EDYN,ECONF,nrXb,cc,kBar,dg) = join([
        "ps___$(pseudo_mode_name(ps))___dt_$(dt)", 
        "DYN_$(EDYN)",
        "ECONF_$(ECONF[1]),$(ECONF[2]),$(ECONF[3])",
        "nrx_$(nrXb[1]),$(nrXb[2]),$(nrXb[3])", 
        "nrxs_$(nrXb[4]),$(nrXb[5]),$(nrXb[6])",
        "nrxb_$(nrXb[7]),$(nrXb[8]),$(nrXb[9])",
        "cutoff_$(cc[1]),$(cc[2])_kBar_$(kBar)_dg_$(dg)"
    ], "_")

    #> ---------------------------------
    #> serial excution
    res = grid_test(    WORKSPACE, 
                        common_title,
                        cif0_fn, 
                        instructor,
                        cp_nt,
                        [ps_modes, dt_list, elec_dynamics_list, elec_conf_list, nrXb_list, cutoff_list, kBar_list, degauss_list]
    )
    return res
end
