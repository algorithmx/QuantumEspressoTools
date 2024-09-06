function GULP_run_watchdog(
    b,
    output_file_of_b::String,
    program::Cmd,
    settings
    )

    watchdog_serial(
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


function GULP_run_report(prog, terminated, status)
    @info "The program gulp has terminated $(terminated) with status $(status)."
end

#* -------------------------------------------------------------------


function GULP_run(
    program::Cmd,
    scripts::AbstractString, 
    workspace::String;
    mode  = :pipeline, 
    f_in  = "tmp.in",
    f_out = "tmp.out",
    watchdog_setting = _default_watchdog_setting
    )::Tuple{Vector{String},Bool,Symbol}

    FD = rstrip(workspace,'/')*"/"
    cd(FD)

    input_file_loc = FD*f_in
    output_file_loc = FD*f_out
    scripts >> input_file_loc

    terminated  = true
    status      = :NotStarted
    success     = false
    program_run = nothing
    results     = String[]

    if mode == :pipeline
        try
            #> run a pipeline :
            #> cat tmp.in | pw.x | tmp.out
            program_run = run(pipeline(`cat $input_file_loc`, program, output_file_loc), wait=false)
            (terminated, status) = GULP_run_watchdog(program_run, output_file_loc, program, watchdog_setting)
        catch _e_
            @error "GULP_run() : Program $(join(program.exec,' ')) lauching error : $(_e_)"
        end
    elseif mode == :file
        try
            program_run = run(pipeline(`$program  $input_file_loc`, output_file_loc), wait=false)
            (terminated, status) = QE_run_watchdog(program_run, output_file_loc, program, watchdog_setting)
        catch _e_
            @error "GULP_run() : Program $(join(program.exec,' ')) lauching error : $(_e_)"
        end
    else
        @error "GULP_run() : Unknown mode $mode. Program $program not run. Exiting ..."
        return results, success, status
    end

    if !terminated && program_run!==nothing
        kill(program_run)
    end

    GULP_run_report(program, terminated, status)

    success = terminated && (status ∈ [ :JobDoneGoodFile, :JobKilledAfterGoodFile, 
                                        :JobEndedGoodFile, :JobSilentTerminatedGoodFile ])
    results = isfile(output_file_loc) ? readlines(output_file_loc) : String[]

    return results, success, status

end


@inline function GULP_run(
    program::Cmd, 
    scripts::Vector{S},
    workspace; 
    mode=:pipeline,
    f_in = "tmp.in", f_out = "tmp.out",
    watchdog_setting = _default_watchdog_setting
    ) where {S<:AbstractString}

    GULP_run(   program,
                join(scripts,"\n"), 
                workspace;
                mode=mode,
                f_in=f_in, f_out=f_out,
                watchdog_setting=watchdog_setting   )

end


#* -------------------------------------------------------------------

function gulp(
    np::Int,
    map_by::String,
    scripts::Vector{S};
    workspace="./",
    fin="gulp.tmp.in", fout="gulp.tmp.out",
    watchdog_setting = _default_watchdog_setting
    )::Tuple{Vector{String},Bool,Symbol} where {S<:AbstractString}
    if np==1
        return GULP_run(   `gulp`, scripts, workspace, 
                            f_in=fin, f_out=fout, watchdog_setting=watchdog_setting )
    else
        opt = map_by=="" ? [] : [map_by]
        if np > 1
            return GULP_run(   `mpiexec  -np  $np  $opt  gulp`, scripts, workspace,
                                f_in=fin, f_out=fout, watchdog_setting=watchdog_setting )
        else # np == -1, oversubscribe
            return GULP_run(   `mpiexec  --oversubscribe  $opt  gulp`, scripts, workspace,
                                f_in=fin, f_out=fout, watchdog_setting=watchdog_setting )
        end
    end
end
