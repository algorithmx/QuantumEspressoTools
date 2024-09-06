test_dry_run_common(title, psmode, kp) = Dict(
        "title"           => title, 
        "prefix"           => title, 
        :pseudo_mode      => psmode,
        "mixing_beta"     => 0.8,
        "electron_maxstep"=> 500,
        :pw_mode          => "energy",
        :kpoint_mode      => "automatic", 
        :kpoints          => kp,
        :watchdog_setting => default_intercepter_setting
)


function test_dry_run(
    workspace0::String, 
    common_title::String, 
    cif0_fn::String, 
    ps_modes::Vector,
    kpoint_configs;
    PROG_PWX = `mpiexec -np 64 --map-by=slot pwdry.x -npool 8`,
    atoms = []
    )

    workspace = try_mkdir(workspace0, ".", "grid_test()")
    cd(workspace)
    WORKSPACE = pwd()

    #> ---------------------------------
    #> customized instruction
    function instructor(i0, FD, psmode, kp, atms=atoms)
        wcut = maximum([100.0, 
                        1.2*make_pseudo_wcut(psmode,atms), 
                        0.3*make_pseudo_ecut(psmode,atms)])
        INSTR =  INSTRUCTION(
            Dict(
                PROG_PWX => _electron_scf_seed_conf_  ⬱ test_dry_run_common(common_title, psmode, kp),
            ),
            [  (   "$(common_title)_pwdry",
                    PROG_PWX,
                    Dict(   :calc       => "scf", 
                            "nstep"     => 0,
                            "nbnd"      => :default,
                            "electron_maxstep" => 1,
                            "conv_thr"  => 1e-8,
                            "ecutwfc"   => wcut,
                            "ecutrho"   => 4wcut,
                            "verbosity" => "high",
                            "disk_io"   => "none",
                    ) ∪ i0 ),
            ]
        )
        return INSTR
    end

    #> ---------------------------------
    #> folders and settings
    nt(ps,kp) = "ps_$(pseudo_mode_name(ps))_kp_$(kp[1]),$(kp[2]),$(kp[3]),$(kp[4]),$(kp[5]),$(kp[6])"

    #> ---------------------------------
    #> serial excution
    res = grid_test(  WORKSPACE, 
                        common_title,
                        cif0_fn, 
                        instructor,
                        nt,
                        [ps_modes, kpoint_configs]
    )

end


##* ===========================================================


