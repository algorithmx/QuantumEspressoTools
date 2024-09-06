global const projwfc_woof = (
    woof_per_x_min=1/10, 
    max_watch=9600, 
    max_silence_woofs=360, 
    tail_length=40, 
    quiet=false
)

##* ========
##*  common
##* ========

project_wfc_common(title, psmode, kp) = Dict(
        "title"           => title, 
        "prefix"          => title, 
        :watchdog_setting => projwfc_woof
)


#TODO
function pick_weight(pw_out)
    ret = Dict()
    try 
        proj, header = projwfcx_output_projwfc_up(pw_out)
        ret = Dict("weights"=>proj, "orbits"=>header)
    catch _
        nothing
    end
    return ret
end



##* ======================
##*  exteranl controllers
##* ======================

project_wfc_external() = Dict()

#! ========
#!   MAIN
#! ========

global const project_wfc_additional_settings = Dict(
        "ngauss"    => 1,
        "Emax"      => 10.0,
        "Emin"      => -10.0,
        "DeltaE"    => 0.005,
        "degauss"   => 0.02
)


function project_wavefunctions(
    workspace0::String, 
    common_title::String, 
    cif0_fn::String, 
    ps_modes::Vector,
    kpoint_configs,
    cutoff_list,
    kBar_list::Vector{Float64},
    degauss_list::Vector{Float64};
    PROG_PROJ          = `mpiexec -np 48 --map-by=slot projwfc.x -npool 8`,
    additional_settings= project_wfc_additional_settings
    )

    workspace = try_mkdir(workspace0, ".", "project_wavefunctions()")
    cd(workspace)
    WORKSPACE = pwd()
    ENV["LOGGER_FILE_FULL_PATH"] = "$WORKSPACE/project_wavefunctions.log"

    #> ---------------------------------
    # customized instruction
    DeltaE =  get(additional_settings, "DeltaE",  0.005)
    emax =  get(additional_settings, "Emax",  10.0)
    emin =  get(additional_settings, "Emin", -10.0)
    try 
        bands = pw_bands(readlines("$WORKSPACE/$FD/$(common_title)_scf/$(common_title).pw.x.out"))
        allb  = vcat(last.(bands)...)
        emax, emin = (maximum(allb)+0.1, minimum(allb)-0.1)
    catch _
        nothing
    end
    ng = get(additional_settings, "ngauss", 1)

    function proj_instructor(i0, FD, psmode, kp, ct, kBar, dg, atms=atoms)
        proj_comm   = project_wfc_common(common_title, psmode, kp)
        INSTR    =  INSTRUCTION(
            Dict(PROG_PROJ => proj_comm
            ),
            [("$(common_title)_projwfc", PROG_PROJ,
                    Dict(   "prefix"    => common_title,
                            "filproj"   => common_title * ".proj",
                            "outdir"    => "$WORKSPACE/$FD/$(common_title)_scf",
                            "Emax"      => emax,
                            "Emin"      => emin,
                            "DeltaE"    => DeltaE,
                            "ngauss"    => ng,
                            "degauss"   => dg,
                            :updater    => x->Dict(),
                    )  ∪ project_wfc_external()),
            ]
        )
        return INSTR
    end

    #> ---------------------------------
    #> serial excution
    res = grid_test(    WORKSPACE, 
                        common_title,
                        cif0_fn, 
                        proj_instructor,
                        _FDNM__ps_kp_ct_kBar_dg_,  #> folders and settings (same as most cases)
                        [ps_modes, kpoint_configs, cutoff_list, kBar_list, degauss_list];
                        cleanup = false
    )
    
    return res
end


function project_wavefunctions_W90(
    workspace0::String, 
    common_title::String, 
    cif0_fn::String, 
    ps_mode, 
    kgrid_scf_nscf::NTuple{6,Int64},
    cutoff::Tuple{Float64,Float64};
    PROG_PROJ          = `mpiexec -np 48 --map-by=slot projwfc.x -npool 8`,
    additional_settings= project_wfc_additional_settings
    )

    workspace = try_mkdir(workspace0, ".", "project_wavefunctions()")
    cd(workspace)
    WORKSPACE = pwd()
    ENV["LOGGER_FILE_FULL_PATH"] = "$WORKSPACE/project_wavefunctions.log"

    #> ---------------------------------
    # customized instruction
    DeltaE =  get(additional_settings, "DeltaE",  0.005)
    emax =  get(additional_settings, "Emax",  10.0)
    emin =  get(additional_settings, "Emin", -10.0)
    try 
        bands = pw_bands(readlines("$WORKSPACE/$FD/$(common_title)_scf/$(common_title).pw.x.out"))
        allb  = vcat(last.(bands)...)
        emax, emin = (maximum(allb)+0.1, minimum(allb)-0.1)
    catch _
        nothing
    end
    ng = get(additional_settings, "ngauss", 1)
    dg = get(additional_settings, "degauss", 0.02)
    function proj_instructor(i0, FD, psmode, kp, ct, atms=[])
        proj_comm   = project_wfc_common(common_title, psmode, kp)
        INSTR    =  INSTRUCTION(
            Dict(PROG_PROJ => proj_comm
            ),
            [("$(common_title)_projwfc", PROG_PROJ,
                    Dict(   "prefix"    => common_title,
                            "filproj"   => common_title * ".proj",
                            "outdir"    => "$WORKSPACE/$FD/$(common_title)_scf",
                            "Emax"      => emax,
                            "Emin"      => emin,
                            "DeltaE"    => DeltaE,
                            "ngauss"    => ng,
                            "degauss"   => dg,
                            :updater    => x->Dict(),
                    )  ∪ project_wfc_external()),
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
        proj_instructor,
        _FDNM__ps_kp_ct_, # settings_to_folder_name
        #% the grid (below)
        [[ps_mode,], [kgrid_scf_nscf,], [cutoff,]],
        cleanup = false
    )

    return res
end
