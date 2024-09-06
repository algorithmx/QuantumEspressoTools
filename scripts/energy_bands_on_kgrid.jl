global const kgrid_woof = (woof_per_x_min=1/10, max_watch=80000, max_silence_woofs=800, tail_length=20, quiet=false)


en_kgrid_common(title, psmode) = Dict(
        "title"           => title, 
        "prefix"          => title, 
        "conv_thr"        => 1e-11,
        "etot_conv_thr"   => 1e-7,
        "forc_conv_thr"   => 1e-6,
        "disk_io"         => "default",
        :pseudo_mode      => psmode,
        :watchdog_setting => kgrid_woof
)


en_kgrid_scf() = Dict()


en_kgrid_band(k_rect) = Dict(:kpoints => k_rect)


function en_kgrid_upd1(x)
    d = dict__pw_relax_result(x)
    return (d ↓ ["ibrav", :cif, :reciprocal_basis]) ⬱ Dict("ibrav"=>0,:do_not_use_symmetry=>true)
end


# for scf calculation, lock the mixing beta and diag method
function en_kgrid_upd2(x)
    d = dict__pw_result(x)
    #beta = pw_mixing_beta(x)
    #diag = pw_diag_style(x)
    ret = Dict(
            :reciprocal_basis => d[:reciprocal_basis], 
            #:beta             => 0.8.*[beta,], 
            #:diag_choices     => [diag,],
            #"ecutwfc"         => d["ecutwfc"], 
            #"ecutrho"         => d["ecutrho"],
            ## "positions" => d["positions"],
    )
    return ret
end


en_kgrid_upd_nothing(x) = Dict()


global const en_kgrid_additional_settings_default = Dict(
        "diag" => ["david",],
        "beta" => [0.8, 0.4],
        "cutoff_scales" => [1.0,],
        "atoms" => [],
)


function en_kgrid_kline(
    kpoint_mode_tpiba__,  # "tpiba_c" or "tpiba_b", OR "crystal_b"
    workspace0::String, 
    common_title::String, 
    cif0_fn::String,  # management
    kpoint_rect_line, #! important
    ps_modes, 
    kpoint_configs,
    cutoff_list;  # grid
    PROG_PWX = `mpiexec -np 64 --map-by=slot pw.x -npool 8`,
    additional_settings = en_kgrid_additional_settings_default
    )

    # management
    workspace = try_mkdir(workspace0, ".", "en_kgrid()")
    cd(workspace)
    WORKSPACE = pwd()
    ENV["LOGGER_FILE_FULL_PATH"] = "$WORKSPACE/en_kgrid.log"

    # customized instruction
    diag = additional_settings["diag"]
    beta = additional_settings["beta"]
    cutoff_scales = additional_settings["cutoff_scales"]
    atoms = additional_settings["atoms"]
    function en_kgrid_instructor(i0, FD, psmode, kp, ct, atms=atoms)
        cutoffs = Dict("ecutwfc" => ct[1], "ecutrho" => ct[2])
        spin_orb = Dict("lspinorb" => occursin("_FR",pseudo_mode_name(psmode)))
        kg_comm  = en_kgrid_common(common_title, psmode) ∪ i0 ∪ cutoffs ∪ spin_orb
        INSTR    = INSTRUCTION(
            Dict(PROG_PWX => (_electron_scf_seed_conf_  ⬱  kg_comm)),
            [   (   "$(common_title)_scf",   
                    PROG_PWX, 
                    Dict(   :calc          => "scf",
                            :beta          => beta,
                            :diag_choices  => diag,
                            :cutoff_scales => cutoff_scales,
                            :pw_mode       => "energy",
                            :kpoint_mode   => "automatic", 
                            :kpoints       => kp,
                            :updater       => en_kgrid_upd_nothing,
                    ) ∪ en_kgrid_scf()
                ),
                (   "$(common_title)_band",
                    PROG_PWX, 
                    Dict(   "outdir"     => "$WORKSPACE/$FD/$(common_title)_scf",
                            :calc        => "bands",
                            :pw_mode     => "bands",
                            :kpoint_mode => kpoint_mode_tpiba__, 
                            :beta          => 0.8.*beta,
                            :diag_choices  => diag,
                    ) ∪ en_kgrid_band(kpoint_rect_line)
                ),
            ]
        )
        return INSTR
    end

    settings_to_folder_name(ps,kp,ct) = 
        "ps___$(pseudo_mode_name(ps))___kp_$(kp[1]),$(kp[2]),$(kp[3]),$(kp[4]),$(kp[5]),$(kp[6])_cut_$(ct[1]),$(ct[2])"

    res = grid_test(
        WORKSPACE, 
        common_title,
        cif0_fn, 
        en_kgrid_instructor,
        settings_to_folder_name,
        #% the grid (below)
        [ps_modes, kpoint_configs, cutoff_list]
    )

    return res
end


en_kgrid(
    workspace0::String, 
    common_title::String, 
    cif0_fn::String,  # management
    kpoint_rect, #! important
    ps_modes, 
    kpoint_configs,
    cutoff_list;  # grid
    PROG_PWX = `mpiexec -np 64 --map-by=slot pw.x -npool 8`,
    additional_settings = en_kgrid_additional_settings_default
    ) =  en_kgrid_kline(
    "tpiba_c",
    workspace0, 
    common_title, 
    cif0_fn,
    kpoint_rect,
    ps_modes, 
    kpoint_configs,
    cutoff_list;  # grid
    PROG_PWX = PROG_PWX,
    additional_settings = additional_settings
)


en_kline(
    workspace0::String, 
    common_title::String, 
    cif0_fn::String,  # management
    kpoint_line, #! important
    ps_modes, 
    kpoint_configs,
    cutoff_list;  # grid
    PROG_PWX = `mpiexec -np 64 --map-by=slot pw.x -npool 8`,
    additional_settings = en_kgrid_additional_settings_default
    ) =  en_kgrid_kline(
    "tpiba_b",
    workspace0, 
    common_title, 
    cif0_fn,
    kpoint_line,
    ps_modes, 
    kpoint_configs,
    cutoff_list;  # grid
    PROG_PWX = PROG_PWX,
    additional_settings = additional_settings
)


en_band(
    workspace0::String, 
    common_title::String, 
    cif0_fn::String,  # management
    #! --------------
    #! important
    high_symm_line,
    #! --------------
    ps_modes, 
    kpoint_configs,
    cutoff_list;  # grid
    PROG_PWX = `mpiexec -np 64 --map-by=slot pw.x -npool 8`,
    additional_settings = en_kgrid_additional_settings_default
    ) =  en_kgrid_kline(
    "crystal_b",
    workspace0, 
    common_title, 
    cif0_fn,
    #! --------------
    #! important
    high_symm_line,   
    #! --------------
    ps_modes, 
    kpoint_configs,
    cutoff_list;  # grid
    PROG_PWX = PROG_PWX,
    additional_settings = additional_settings
)


## ================================================================


function en_kgrid_line__STATUS(
    kpoint_mode_tpiba__,
    workspace0::String, 
    common_title::String, 
    cif0_fn::String,  # management
    kpoint_rect, #! important
    ps_modes, 
    kpoint_configs;  # grid
    PROG_PWX = `mpiexec -np 64 --map-by=slot pw.x -npool 8`,
    additional_settings = en_kgrid_additional_settings_default
    )

    cd(workspace0)
    WORKSPACE = pwd()

    #> -------- IDENTCAL PART  ---------
    diag = additional_settings["diag"]
    beta = additional_settings["beta"]
    cutoff_scales = additional_settings["cutoff_scales"]
    atoms = additional_settings["atoms"]
    function en_kgrid_instructor(i0, FD, psmode, kp, ct, atms=atoms)
        cutoffs = Dict("ecutwfc" => ct[1], "ecutrho" => ct[2])
        spin_orb = Dict("lspinorb" => occursin("_FR",pseudo_mode_name(psmode)))
        kg_comm  = en_kgrid_common(common_title, psmode) ∪ i0 ∪ cutoffs ∪ spin_orb
        INSTR    = INSTRUCTION(
            Dict(PROG_PWX => (_electron_scf_seed_conf_  ⬱  kg_comm)),
            [   (   "$(common_title)_scf",   
                    PROG_PWX, 
                    Dict(   :calc          => "scf",
                            :beta          => beta,
                            :diag_choices  => diag,
                            :cutoff_scales => cutoff_scales,
                            :pw_mode       => "energy",
                            :kpoint_mode   => "automatic", 
                            :kpoints       => kp,
                            :updater       => en_kgrid_upd2,
                    ) ∪ en_kgrid_scf()
                ),
                (   "$(common_title)_band",
                    PROG_PWX, 
                    Dict(   "outdir"     => "$WORKSPACE/$FD/$(common_title)_scf",
                            :calc        => "bands",
                            :pw_mode     => "bands",
                            :kpoint_mode => kpoint_mode_tpiba__, #"tpiba_c", 
                    ) ∪ en_kgrid_band(kpoint_rect)
                ),
            ]
        )
        return INSTR
    end

    settings_to_folder_name(ps,kp,ct) = 
        "ps___$(pseudo_mode_name(ps))___kp_$(kp[1]),$(kp[2]),$(kp[3]),$(kp[4]),$(kp[5]),$(kp[6])_cut_$(ct[1]),$(ct[2])"

    #> -------- DIFFERNT PART  ---------
    res = grid_test__STATUS(  
        WORKSPACE, 
        common_title,
        cif0_fn, 
        en_kgrid_instructor,
        settings_to_folder_name,
        #% the grid (below)
        [ps_modes, kpoint_configs, cutoff_list]
    )

    return res
end


function en_kgrid_line__STATUS_print(status0)
    @inline spltend(x) = split(x,"_",keepempty=false)[end]
    status = sort(status0, by=first)
    for (fd, content) ∈ status
        (config_i, status_list) = content
        L2 = []
        for (A,B) in status_list
            (plan_name,program) = A
            (status,FD) = B
            push!(L2, spltend(plan_name)=>status)    
        end
        println(join(string.(last.(L2)), " ") * "    " * fd)
    end
    return
end
