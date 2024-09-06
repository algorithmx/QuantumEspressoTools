global const woof = (woof_per_x_min=1/2, max_watch=5000, max_silence_woofs=12, tail_length=30, quiet=false)


scf_energy_common(title, psmode, kp, nbnd, dg) = Dict(
    "title"           => title, 
    "prefix"           => title, 

    :pseudo_mode      => psmode,
    :kpoint_mode      => "automatic",
    :kpoints          => kp,
    "nbnd"            => nbnd,     

    "degauss"         => dg,

    :watchdog_setting => woof,
)


#TODO abstract `instructor`
function scf_energy(
    workspace0::String, 
    common_title::String, 
    cif0_fn::String, 
    ps_modes::Vector,
    cutoffs::Vector{Tuple{Float64,Float64}},
    kpoint_configs::Vector,
    nband_list::Vector,
    degauss_list::Vector{Float64};
    beta = [0.8, 0.5, 0.3],
    cutoff_scales = [1.0, 1.25, 1.5],
    PROG_PWX = `mpiexec -np 64 --map-by=slot pw.x -npool 8`,
    atoms = []
    )

    workspace = try_mkdir(workspace0, ".", "grid_test()")
    cd(workspace)
    WORKSPACE = pwd()
    ENV["LOGGER_FILE_FULL_PATH"] = "$WORKSPACE/scf_energy.log"

    #> ---------------------------------
    #> customized instruction
    #> ins = instructor(init_struct, fd, cc...)
    function instructor(structure0, FD, psmode, cuts, kp, nbnd, dg, atms=atoms)
        wcut,ecut = cuts
        scf_comm  = scf_energy_common(common_title, psmode, kp, nbnd, dg)
        spin_orb  = Dict("lspinorb" => occursin("_FR",pseudo_mode_name(psmode)))
        INSTR     = INSTRUCTION(
            Dict(   PROG_PWX => _electron_scf_seed_conf_  ⬱ scf_comm ),
            [   (   common_title,   PROG_PWX, 
                    Dict(   :calc          => "scf",
                            "outdir"       => "$WORKSPACE/$FD/$(common_title)",
                            :beta          => beta,
                            :cutoff_scales => cutoff_scales,
                            "ecutwfc"      => wcut,
                            "ecutrho"      => max(ecut,4wcut),
                    ) ∪ structure0 ∪ spin_orb )
            ]
        )
        return INSTR
    end

    #> ---------------------------------
    #> folders and settings
    nt(ps,ct,kp,nbd,dg) = "ps___$(pseudo_mode_name(ps))___cut_$(ct[1])_$(ct[2])_kp_$(kp[1]),$(kp[2]),$(kp[3]),$(kp[4]),$(kp[5]),$(kp[6])_nbd_$(nbd)_dg_$(dg)"

    #> ---------------------------------
    res = grid_test(
            WORKSPACE, 
            common_title,
            cif0_fn, 
            instructor,
            nt,
            [ps_modes, cutoffs, kpoint_configs, nband_list, degauss_list]
    )

    return res

end

##* ===========================================================


