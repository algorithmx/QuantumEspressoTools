
function wannier90_run_watchdog(
    b,
    output_file_of_b::String,
    program::Cmd,
    settings
    )
    return watchdog_serial(
            b, 
            output_file_of_b,
            r->verify_result(r,program;final=false);
            woof_per_x_min     = settings[:woof_per_x_min],
            max_watch          = settings[:max_watch],
            max_silence_woofs  = settings[:max_silence_woofs],
            tail_length        = settings[:tail_length],
            quiet              = settings[:quiet],
            grace_time_sec     = get(settings,:grace_time_sec,300) 
        )
end


function wannier90_run_report(prog, terminated, status)
    @info "The program Wannier90 has terminated $(terminated) with status $(status)."
end

#* -------------------------------------------------------------------

function wannier90_run(
    program::Cmd,
    scripts::AbstractString, 
    workspace::String;
    IGNORE_SCRIPTS = false,
    f_root  = "tmp",
    watchdog_setting = _default_watchdog_setting
    )::Tuple{Vector{String},Bool,Symbol}

    FD = rstrip(workspace,'/')*"/"
    cd(FD)

    input_file_loc  = FD*f_root*".win"
    output_file_loc = FD*f_root*".program.out"
    wout_file_loc   = FD*f_root*".wout"
    if !IGNORE_SCRIPTS
        scripts >> input_file_loc
    end

    terminated  = true
    status      = :NotStarted
    success     = false
    program_run = nothing
    results     = String[]

    try
        program_run = run(pipeline(`$program  $f_root`, output_file_loc), wait=false)
        (terminated, status) = wannier90_run_watchdog(program_run, wout_file_loc, program, watchdog_setting)
    catch _e_
        @error "wannier90_run() : Program $(join(program.exec,' ')) lauching error : $(_e_)"
    end

    if !terminated && program_run!==nothing
        kill(program_run)
    end

    wannier90_run_report(program, terminated, status)

    success = terminated && (status ∈ [ :JobDoneGoodFile, :JobKilledAfterGoodFile, 
                                        :JobEndedGoodFile, :JobSilentTerminatedGoodFile ])
    results = isfile(wout_file_loc) ? readlines(wout_file_loc) : String[]

    return results, success, status

end


@inline function wannier90_run(
    program::Cmd, 
    scripts::Vector{S},
    workspace; 
    IGNORE_SCRIPTS = false,
    f_root = "tmp", 
    watchdog_setting = _default_watchdog_setting
    ) where {S<:AbstractString}

    wannier90_run(   
        program,
        join(scripts,"\n"), 
        workspace;
        f_root = f_root,
        IGNORE_SCRIPTS = IGNORE_SCRIPTS,
        watchdog_setting = watchdog_setting    
    )

end


#* -------------------------------------------------------------------


function wannier90(
    scripts;
    workspace="./",
    f_root="wannier90.tmp",
    watchdog_setting = _default_watchdog_setting
    )::Tuple{Vector{String},Bool,Symbol}
    return wannier90_run(  `wannier90.x`, scripts, workspace; 
                            IGNORE_SCRIPTS = true,
                            f_root=f_root, watchdog_setting=watchdog_setting  )
end


function wannier90pp(
    scripts;
    workspace="./",
    f_root="wannier90pp.tmp",
    watchdog_setting = _default_watchdog_setting
    )::Tuple{Vector{String},Bool,Symbol}
    return wannier90_run(  `wannier90.x -pp`, scripts, workspace, 
                            IGNORE_SCRIPTS = false,
                            f_root=f_root, watchdog_setting=watchdog_setting  )
end
