global const relax_woof = (
                woof_per_x_min=1/3, 
                max_watch=20000, 
                max_silence_woofs=200, 
                tail_length=20, 
                quiet=false
              )


relax_common(title, psmode, kp, kBar, dg) = Dict(
        # common
        "title"           => title, 
        "prefix"           => title, 
        "disk_io"         => "medium",
        :pseudo_mode      => psmode,

        # k-points
        :kpoint_mode      => "automatic", 
        :kpoints          => kp,

        "conv_thr"        => 1e-9,
        "etot_conv_thr"   => 1e-7,
        "forc_conv_thr"   => 1e-6,

        # spin-orbit
        "lspinorb"        => occursin("_FR",pseudo_mode_name(psmode)),

        # relax
        "smearing"        => "mp",
        "degauss"         => dg,
        "press"           => Float64(kBar),

        # watch dog
        :watchdog_setting => (
            woof_per_x_min=1/3, 
            max_watch=20000, 
            max_silence_woofs=200, 
            tail_length=20, 
            quiet=false
        )
)


function pw_to_cif_old(result_lines, fn)
    try
        (result_lines |> pw_relax_result_to_cif_lines) ⇶ fn
    catch _e_
        @error "pw_to_cif() failed : error = $(_e_)"
    end
    return Dict()
end


function pw_to_cif(result_lines, fn)
    try
        latt_p   = QuantumEspressoTools.pw_relax_to_lattice_parameters(result_lines)
        atm_list = QuantumEspressoTools.pw_relax_atom_list(result_lines)
        cif = QuantumEspressoTools.minimal_cif("generated_by_pw_to_cif", latt_p, atm_list)
        QuantumEspressoTools.write_to_file(cif,fn)
    catch _e_
        @error "pw_to_cif() failed : error = $(_e_)"
    end
    return Dict()
end


function diff_cif(cif1,cif2)
    lat_diff  = [ maximum(abs.([p1...,].-[p2...,]))
                  for (p1,p2) ∈ zip(get_cell_params(cif1),get_cell_params(cif2)) ]
    pos_diff  = [ maximum(abs.([p1[2]-p2[2], p1[3]-p2[3], p1[4]-p2[4]]))
                  for (p1,p2) ∈ zip(get_atom_frac_pos(cif1),get_atom_frac_pos(cif2)) ]
    return maximum([lat_diff..., pos_diff...])
end


#: ==================================================================

global const vc_relax_additional_settings_default = Dict(
    "diag" => ["david",],  #! no cg choice 
    "beta" => [0.8, 0.4],
    "cutoff_scales" => [1.0,],
    "atoms" => [],
    "assume_isolated" => "none",
    "folder_number_shift" => 0,
    "stop_crit" => 1e-6,
)


function vc_relax(
    cell_dofree__conv_thr__kick_amp__list::Vector{Tuple{String,Float64,Float64}},  #! difference
    workspace0::String, 
    common_title::String, 
    cif0_fn::String, 
    ps_modes::Vector,
    kpoint_configs::Vector{NTuple{6,Int64}},
    cutoff_list::Vector{Tuple{Float64,Float64}},
    kBar_list::Vector{Float64},
    degauss_list::Vector{Float64};
    cleanup = true,
    kicker = relax_update_random_kick,  #! difference
    PROG_PWX = `mpiexec -np 48 pw.x -npool 6`,
    additional_settings = Dict()
    )

    workspace = try_mkdir(workspace0, ".", "vc_relax()")
    cd(workspace)
    WORKSPACE = pwd()
    ENV["LOGGER_FILE_FULL_PATH"] = "$WORKSPACE/vc_relax.log"  #TODO  thread unsafe

    #> ---------------------------------
    #> customized instruction
    @inline gtdf(X) = get(additional_settings, X, vc_relax_additional_settings_default[X])
    diag  = gtdf("diag")
    beta  = gtdf("beta")
    ctf   = gtdf("cutoff_scales")
    atoms = gtdf("atoms")
    iso  = gtdf("assume_isolated")
    itern_SHIFT = gtdf("folder_number_shift")
    stop_crit_thr = gtdf("stop_crit")

    function vc_relax_instructor(init_struct, FD, psmode, kp, cutoffs, kBar, dg, atms_to_kick=atoms)
        wcut, ecut = cutoffs
        relax_comm  = (relax_common(common_title, psmode, kp, kBar, dg) ↓ ["cell_dofree",])  #! difference
        @inline mkfn(ii) = "$WORKSPACE/$FD/$(common_title)_vc_relax_iter_$(ii).cif"
        function updater_local(x,ii,km,el)
            pw_to_cif(x,mkfn(ii))
            @info "updater_local() : cif  = \n$(join(readlines(mkfn(ii)),"\n"))"
            exit_label = false
            if ii >= 2
                dff = diff_cif(mkfn(ii),mkfn(ii-1))
                if dff < stop_crit_thr
                    @info "\n---------------------------\nvc_relax_instructor() : updater_local() : diff_cif = $(dff), structure is stable.\n---------------------------\n"
                    exit_label = true
                end
            end
            ret = kicker(x,km,el) ∪ (exit_label ? Dict(:exit=>true) : Dict())
            @info "updater_local() : \n$(showdict(ret))"
            return ret
        end
        INSTR = INSTRUCTION(
            # common settings
            Dict(PROG_PWX => (_electron_relax_seed_conf_  ⬱ relax_comm)),
            # steps
            [   (   "$(common_title)_vc_relax_iter_$(itern+itern_SHIFT)", 
                    PROG_PWX,
                    (Dict(  :calc             => "vc-relax", 
                            :beta             => beta,
                            :diag_choices     => diag,  
                            :cutoff_scales    => ctf,
                            "conv_thr"        => 0.01*thr,  #! difference
                            "etot_conv_thr"   => thr,  #! difference
                            "assume_isolated" => iso,
                            "ecutwfc"         => wcut,
                            "ecutrho"         => ecut,
                            :updater          => x->updater_local(x,itern,kickamp,atms_to_kick), #! difference
                    ) ∪ ((itern==1) ? (init_struct) : Dict()) ) ↑ ("cell_dofree"=>dofree) ) 
                for (itern, (dofree,thr,kickamp)) ∈ enumerate(cell_dofree__conv_thr__kick_amp__list) #! difference
            ]
        )
        return INSTR
    end

    #> ---------------------------------
    #> serial excution
    res = grid_test(    WORKSPACE, 
                        common_title,
                        cif0_fn, 
                        vc_relax_instructor,
                        _FDNM__ps_kp_ct_kBar_dg_,
                        [ps_modes, kpoint_configs, cutoff_list, kBar_list, degauss_list];
                        cleanup = cleanup
    )
    return res
end


function vc_relax(
    cell_dofree_list::Vector{String},
    workspace0::String, 
    common_title::String, 
    cif0_fn::String, 
    ps_modes::Vector,
    kpoint_configs::Vector{NTuple{6,Int64}},
    cutoff_list::Vector{Tuple{Float64,Float64}},
    kBar_list::Vector{Float64},
    degauss_list::Vector{Float64};
    cleanup = true,
    relax_iter_updater = x->relax_update_random_kick(x,0.002,[]),
    PROG_PWX = `mpiexec -np 64 --map-by=slot pw.x -npool 8`,
    additional_settings = Dict()
    )

    workspace = try_mkdir(workspace0, ".", "vc_relax()")
    cd(workspace)
    WORKSPACE = pwd()
    ENV["LOGGER_FILE_FULL_PATH"] = "$WORKSPACE/vc_relax.log"

    #> ---------------------------------
    #> customized instruction
    @inline gtdf(X) = get(additional_settings, X, vc_relax_additional_settings_default[X])
    diag  = gtdf("diag")
    beta  = gtdf("beta")
    ctf   = gtdf("cutoff_scales")
    atoms = gtdf("atoms")
    iso  = gtdf("assume_isolated")
    itern_SHIFT = gtdf("folder_number_shift")

    function vc_relax_instructor(init_struct, FD, psmode, kp, cutoffs, kBar, dg, atms=atoms)
        wcut, ecut = cutoffs
        relax_comm  = relax_common(common_title, psmode, kp, kBar, dg)
        @inline mkfn(ii) = "$WORKSPACE/$FD/$(common_title)_vc_relax_iter_$(ii).cif"
        INSTR = INSTRUCTION(
            # common settings
            Dict(PROG_PWX => (_electron_relax_seed_conf_  ⬱ (relax_comm ↓ "cell_dofree"))),
            # steps
            [   (   "$(common_title)_vc_relax_iter_$(itern+itern_SHIFT)", 
                    PROG_PWX,
                    (Dict(  :calc             => "vc-relax", 
                            :beta             => beta,
                            :diag_choices     => diag,  
                            :cutoff_scales    => ctf,
                            "assume_isolated" => iso,
                            "ecutwfc"         => wcut,
                            "ecutrho"         => ecut,
                            :updater          => x->(pw_to_cif(x,mkfn(ii)) ∪ relax_iter_updater(x)),
                    ) ∪ ((itern==1) ? (init_struct) : Dict()) ) ↑ ("cell_dofree"=>dofree) ) 
                for (itern, dofree) ∈ enumerate(cell_dofree_list)
            ]
        )
        return INSTR
    end

    #> ---------------------------------
    #> serial excution
    res = grid_test(    WORKSPACE, 
                        common_title,
                        cif0_fn, 
                        vc_relax_instructor,
                        _FDNM__ps_kp_ct_kBar_dg_,
                        [ps_modes, kpoint_configs, cutoff_list, kBar_list, degauss_list];
                        cleanup = cleanup
    )
    return res
end


function vc_relax(
        dofree_N::Tuple{String,Int},
        workspace0::String, 
        common_title::String, 
        cif0_fn::String, 
        ps_modes::Vector,
        kpoint_configs::Vector{NTuple{6,Int64}},
        cutoff_list::Vector{Tuple{Float64,Float64}},
        kBar_list::Vector{Float64},
        degauss_list::Vector{Float64};
        cleanup = true,
        relax_iter_updater = x->relax_update_random_kick(x,0.002,[]),
        PROG_PWX = `mpiexec -np 64 --map-by=slot pw.x -npool 8`,
        additional_settings = Dict()
    )

    vc_relax(
        String[dofree_N[1] for i=1:dofree_N[2]],
        workspace0, 
        common_title, 
        cif0_fn, 
        ps_modes,
        kpoint_configs,
        cutoff_list,
        kBar_list,
        degauss_list;
        cleanup = cleanup,
        relax_iter_updater = relax_iter_updater,
        PROG_PWX = PROG_PWX,
        additional_settings = additional_settings

    )
end


function vc_relax(
        Nkick::Int,
        workspace0::String, 
        common_title::String, 
        cif0_fn::String, 
        ps_modes::Vector,
        kpoint_configs::Vector{NTuple{6,Int64}},
        cutoff_list::Vector{Tuple{Float64,Float64}},
        kBar_list::Vector{Float64},
        degauss_list::Vector{Float64};
        cleanup = true,
        PROG_PWX = `mpiexec -np 64 --map-by=slot pw.x -npool 8`,
        additional_settings = Dict()
    )

    vc_relax(
        ("all", Nkick),
        workspace0, 
        common_title, 
        cif0_fn, 
        ps_modes,
        kpoint_configs,
        cutoff_list,
        kBar_list,
        degauss_list;
        cleanup = cleanup,
        relax_iter_updater = update_nothing,
        PROG_PWX = PROG_PWX,
        additional_settings = additional_settings
    )
end


function vc_relax(
        workspace0::String, 
        common_title::String, 
        cif0_fn::String, 
        ps_modes::Vector,
        kpoint_configs::Vector{NTuple{6,Int64}},
        cutoff_list::Vector{Tuple{Float64,Float64}},
        kBar_list::Vector{Float64},
        degauss_list::Vector{Float64};
        cleanup = true,
        PROG_PWX = `mpiexec -np 64 --map-by=slot pw.x -npool 8`,
        additional_settings = Dict()
    )

    vc_relax(
        1,
        workspace0, 
        common_title, 
        cif0_fn, 
        ps_modes,
        kpoint_configs,
        cutoff_list,
        kBar_list,
        degauss_list;
        cleanup = cleanup,
        relax_iter_updater = update_nothing,
        PROG_PWX = PROG_PWX,
        additional_settings = additional_settings
    )
end


function vc_relax_2Dxy(
        Nkick::Int,
        workspace0::String, 
        common_title::String, 
        cif0_fn::String, 
        ps_modes::Vector,
        kpoint_configs::Vector{NTuple{6,Int64}},
        cutoff_list::Vector{Tuple{Float64,Float64}},
        kBar_list::Vector{Float64},
        degauss_list::Vector{Float64};
        cleanup = true,
        relax_iter_updater = x->relax_update_random_kick(x,0.002,[]),
        PROG_PWX = `mpiexec -np 64 --map-by=slot pw.x -npool 8`,
        additional_settings = Dict()
    )

    vc_relax(
        ("2Dxy",Nkick), #! not fc_relax(), have to specify cell_dofree
        workspace0, 
        common_title, 
        cif0_fn, 
        ps_modes,
        kpoint_configs,
        cutoff_list,
        kBar_list,
        degauss_list;
        cleanup = cleanup,
        relax_iter_updater = relax_iter_updater,
        PROG_PWX = PROG_PWX,
        additional_settings = (additional_settings ⬱ Dict("assume_isolated"=>"2D"))
    )
end


function vc_relax_2Dxy(
        workspace0::String, 
        common_title::String, 
        cif0_fn::String, 
        ps_modes::Vector,
        kpoint_configs::Vector{NTuple{6,Int64}},
        cutoff_list::Vector{Tuple{Float64,Float64}},
        kBar_list::Vector{Float64},
        degauss_list::Vector{Float64};
        cleanup = true,
        relax_iter_updater = x->relax_update_random_kick(x,0.002,[]),
        PROG_PWX = `mpiexec -np 64 --map-by=slot pw.x -npool 8`,
        additional_settings = Dict()
    )

    vc_relax_2Dxy(
        1,
        workspace0, 
        common_title, 
        cif0_fn, 
        ps_modes,
        kpoint_configs,
        cutoff_list,
        kBar_list,
        degauss_list;
        cleanup = cleanup,
        relax_iter_updater = relax_iter_updater,
        PROG_PWX = PROG_PWX,
        additional_settings = additional_settings

    )
end

#: ==================================================================

global const fc_relax_additional_settings_default = Dict(
    "diag" => ["david",],  #! no cg choice 
    "beta" => [0.8, 0.4],
    "cutoff_scales" => [1.0,],
    "atoms" => [],
    "assume_isolated" => "none",
    "folder_number_shift" => 0,
    "stop_crit" => 1e-6,
)


function fc_relax(
    Nkick::Int,
    workspace0::String, 
    common_title::String, 
    cif0_fn::String, 
    ps_modes::Vector,
    kpoint_configs::Vector{NTuple{6,Int64}},
    cutoff_list::Vector{Tuple{Float64,Float64}},
    kBar_list::Vector{Float64},
    degauss_list::Vector{Float64};
    cleanup = true,
    relax_iter_updater = x->relax_update_random_kick(x,0.002,[]),
    PROG_PWX = `mpiexec -np 64 --map-by=slot pw.x -npool 8`,
    additional_settings = Dict()
    )

    workspace = try_mkdir(workspace0, ".", "fc_relax()")
    cd(workspace)
    WORKSPACE = pwd()
    ENV["LOGGER_FILE_FULL_PATH"] = "$WORKSPACE/fc_relax.log"

    #> ---------------------------------
    #> customized instruction
    @inline gtdf(X) = get(additional_settings, X, fc_relax_additional_settings_default[X])
    diag  = gtdf("diag")
    beta  = gtdf("beta")
    ctf   = gtdf("cutoff_scales")
    atoms = gtdf("atoms")
    iso   = gtdf("assume_isolated")
    itern_SHIFT = gtdf("folder_number_shift")

    function fc_relax_instructor(init_struct, FD, psmode, kp, cutoffs, kBar, dg, atms=atoms)
        wcut, ecut = cutoffs
        relax_comm = relax_common(common_title, psmode, kp, kBar, dg)
        @inline mkfn(ii) = "$WORKSPACE/$FD/$(common_title)_fc_relax_iter_$(ii).cif"
        INSTR  =  INSTRUCTION(
            Dict(PROG_PWX => (_electron_relax_seed_conf_  ⬱ (relax_comm ↓ "cell_dofree"))),
            [   (   "$(common_title)_fc_relax_iter_$(itern+itern_SHIFT)", 
                    PROG_PWX,
                   (Dict(   :calc             => "relax", 
                            :beta             => beta,
                            :diag_choices     => diag,  
                            :cutoff_scales    => ctf,
                            "assume_isolated" => iso,
                            "ecutwfc"   => wcut,
                            "ecutrho"   => ecut,
                            :updater    => x -> ( pw_to_cif(x,mkfn(ii)) ∪ relax_iter_updater(x) ),
                    ) ∪ ((itern==1) ? (init_struct) : Dict())) ) 
                for itern = 1:Nkick
            ]
        )
        return INSTR
    end

    #> ---------------------------------
    #> serial excution
    res = grid_test(  
            WORKSPACE, 
            common_title,
            cif0_fn, 
            fc_relax_instructor,
            _FDNM__ps_kp_ct_kBar_dg_,
            [ps_modes, kpoint_configs, cutoff_list, kBar_list, degauss_list];
            cleanup = cleanup
    )
    return res
end


function fc_relax(
    etot_conv_thr__kick_amp__list::Vector{Tuple{Float64,Float64}},  #! difference
    workspace0::String, 
    common_title::String, 
    cif0_fn::String, 
    ps_modes::Vector,
    kpoint_configs::Vector{NTuple{6,Int64}},
    cutoff_list::Vector{Tuple{Float64,Float64}},
    kBar_list::Vector{Float64},
    degauss_list::Vector{Float64};
    cleanup = true,
    kicker = relax_update_random_kick,  #! difference
    PROG_PWX = `mpiexec -np 48 pw.x -npool 6`,
    additional_settings = Dict()
    )

    workspace = try_mkdir(workspace0, ".", "fc_relax()")
    cd(workspace)
    WORKSPACE = pwd()
    ENV["LOGGER_FILE_FULL_PATH"] = "$WORKSPACE/fc_relax.log"

    #> ---------------------------------
    #> customized instruction
    @inline gtdf(X) = get(additional_settings, X, fc_relax_additional_settings_default[X])
    diag  = gtdf("diag")
    beta  = gtdf("beta")
    ctf   = gtdf("cutoff_scales")
    atoms = gtdf("atoms")
    iso   = gtdf("assume_isolated")
    itern_SHIFT = gtdf("folder_number_shift")
    stop_crit_thr = gtdf("stop_crit")

    function fc_relax_instructor(init_struct, FD, psmode, kp, cutoffs, kBar, dg, atms_to_kick=atoms)
        wcut, ecut = cutoffs
        relax_comm  = (relax_common(common_title, psmode, kp, kBar, dg) ↓ ["cell_dofree"])  #! difference
        @inline mkfn(ii) = "$WORKSPACE/$FD/$(common_title)_fc_relax_iter_$(ii).cif"
        function updater_local(x,ii,km,el)
            pw_to_cif(x,mkfn(ii))
            @info "updater_local() : cif  = \n$(join(readlines(mkfn(ii)),"\n"))"
            exit_label = false
            if ii >= 2
                dff = diff_cif(mkfn(ii),mkfn(ii-1))
                if dff < stop_crit_thr
                    @info "\n---------------------------\nfc_relax_instructor() : updater_local() : diff_cif = $(dff), structure is stable.\n---------------------------\n"
                    exit_label = true
                end
            end
            ret = kicker(x,km,el) ∪ (exit_label ? Dict(:exit=>true) : Dict())
            @info "updater_local() : \n$(showdict(ret))"
            return ret
        end
        INSTR  =  INSTRUCTION(
            Dict(PROG_PWX => (_electron_relax_seed_conf_  ⬱ relax_comm)),
            [   (   "$(common_title)_fc_relax_iter_$(itern+itern_SHIFT)", 
                    PROG_PWX,
                   (Dict(   :calc             => "relax", 
                            :beta             => beta,
                            :diag_choices     => diag,  
                            :cutoff_scales    => ctf,
                            "conv_thr"        => 0.01*thr,  #! difference
                            "etot_conv_thr"   => thr,  #! difference
                            "assume_isolated" => iso,
                            "ecutwfc"         => wcut,
                            "ecutrho"         => ecut,
                            :updater          => x->updater_local(x,itern,kickamp,atms_to_kick), #! difference
                    ) ∪ ((itern==1) ? (init_struct) : Dict())) ) 
                for (itern, (thr,kickamp)) ∈ enumerate(etot_conv_thr__kick_amp__list) #! difference
            ]
        )
        return INSTR
    end

    #> ---------------------------------
    #> serial excution
    res = grid_test(  
            WORKSPACE, 
            common_title,
            cif0_fn, 
            fc_relax_instructor,
            _FDNM__ps_kp_ct_kBar_dg_,
            [ps_modes, kpoint_configs, cutoff_list, kBar_list, degauss_list];
            cleanup = cleanup
    )
    return res
end


function fc_relax(
        workspace0::String, 
        common_title::String, 
        cif0_fn::String, 
        ps_modes::Vector,
        kpoint_configs::Vector{NTuple{6,Int64}},
        cutoff_list::Vector{Tuple{Float64,Float64}},
        kBar_list::Vector{Float64},
        degauss_list::Vector{Float64};
        cleanup = true,
        relax_iter_updater = x->relax_update_random_kick(x,0.002,[]),
        PROG_PWX = `mpiexec -np 64 --map-by=slot pw.x -npool 8`,
        additional_settings = Dict()
    )

    fc_relax(
        1,
        workspace0, 
        common_title, 
        cif0_fn, 
        ps_modes,
        kpoint_configs,
        cutoff_list,
        kBar_list,
        degauss_list;
        cleanup = cleanup,
        relax_iter_updater = relax_iter_updater,
        PROG_PWX = PROG_PWX,
        additional_settings = additional_settings
    )
end
