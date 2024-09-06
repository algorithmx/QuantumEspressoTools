#TODO monitor the SIGABRT
#+ need to launch `pw.x` asynchronously

#TODO "ghost run" by output equivalent bash script >> consider logger

global const _default_watchdog_setting = (
    woof_per_x_min=1, 
    max_watch=60*24*3,     # 3 days
    max_silence_woofs=60*2,  # 2 hour
    tail_length=40, 
    quiet=false
)

global const _default_intercepter_setting = (
    woof_per_x_min=1/120,
    max_watch=120*4, 
    max_silence_woofs=120*3, 
    tail_length=300, 
    quiet=false
)

global const __QE_PROGS__ = [
    "pw.x", 
    "cp.x", 
    "ph.x", 
    "bands.x", 
    "q2r.x", 
    "matdyn.x", 
    "pw2wannier90.x",
    "projwfc.x"
]


#* -------------------------------------------------------------------
#* -------------------------------------------------------------------


function QE_run_watchdog(
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


function QE_dry_run_watchdog(
    b,
    output_file_of_b::String,
    stop_crit::Union{String,Regex},
    settings
    )

    watchdog_intercepter(
        b, 
        output_file_of_b,
        stop_crit;
        woof_per_x_min     = settings[:woof_per_x_min],
        max_watch          = settings[:max_watch],
        max_silence_woofs  = settings[:max_silence_woofs],
        tail_length        = settings[:tail_length],
        quiet              = settings[:quiet],
    )

end


function QE_run_report(program::Cmd, terminated, status)
    @info "The $(join(program.exec, " ")) has terminated $(terminated) with status $(status)."
end


#* -------------------------------------------------------------------
#* -------------------------------------------------------------------


function QE_dry_run(
    program::Cmd,
    scripts::AbstractString, 
    workspace::String,
    stop_crit::Union{String,Regex};
    mode  = :pipeline, 
    f_in  = "tmp.in",
    f_out = "tmp.out",
    watchdog_setting = _default_intercepter_setting
    )::Tuple{Vector{String},Bool,Symbol}

    FD = rstrip(workspace,'/')*"/"
    cd(FD)
    TMP_DIR0 = get(ENV,"ESPRESSO_TMPDIR","./")
    ENV["ESPRESSO_TMPDIR"] = FD

    input_file_loc = FD*f_in
    output_file_loc = FD*f_out
    scripts >> input_file_loc

    terminated  = true
    status      = :NotStarted
    success     = false
    program_run = nothing
    results     = String[]

    if mode == :pipeline
        #> run a pipeline :
        try
            #> cat tmp.in | pw.x | tmp.out
            program_run = run(pipeline(`cat $input_file_loc`, program, output_file_loc), wait=false)
            (terminated, status) = QE_dry_run_watchdog(program_run, output_file_loc, stop_crit, watchdog_setting)
        catch _e_
            @error "QE_dry_run() : Program $(join(program.exec,' ')) lauching error : $(_e_)"
        end
    elseif mode == :file
        try
            program_run = run(pipeline(`$program  $input_file_loc`, output_file_loc), wait=false)
            (terminated, status) = QE_dry_run_watchdog(program_run, output_file_loc, stop_crit, watchdog_setting)
        catch _e_
            @error "QE_dry_run() : Program $(join(program.exec,' ')) lauching error : $(_e_)"
        end
    else
        @error "QE_run() : Unknown mode $mode. Program $program not run. Exiting ..."
        ENV["ESPRESSO_TMPDIR"] = TMP_DIR0
        return results, success, status
    end

    if !terminated && program_run!==nothing
        kill(program_run)
    end

    QE_run_report(program, terminated, status)

    success = terminated && (status ∈ [ :JobDoneCritSat,  :JobKilledAfterCritSat, 
                                        :JobEndedCritSat, :JobSilentTerminatedCritSat ])
    results = isfile(output_file_loc) ? readlines(output_file_loc) : String[]

    ENV["ESPRESSO_TMPDIR"] = TMP_DIR0
    return results, success, status

end



@inline function QE_dry_run(
    program::Cmd, 
    scripts::Vector{S},
    workspace::String,
    stop_crit::Union{String,Regex};
    mode=:pipeline,
    f_in = "tmp.in", f_out = "tmp.out",
    watchdog_setting = _default_intercepter_setting
    ) where {S<:AbstractString}

    QE_dry_run( program,
            join(scripts,"\n"), 
            workspace,
            stop_crit;
            mode=mode,
            f_in=f_in, f_out=f_out,
            watchdog_setting=watchdog_setting )

end


#* -------------------------------------------------------------------
#* -------------------------------------------------------------------


function QE_run(
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
    TMP_DIR0 = get(ENV,"ESPRESSO_TMPDIR","./")
    ENV["ESPRESSO_TMPDIR"] = FD

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
            (terminated, status) = QE_run_watchdog(program_run, output_file_loc, program, watchdog_setting)
        catch _e_
            @error "QE_run() : Program $(join(program.exec,' ')) lauching error : $(_e_)"
        end
    elseif mode == :file
        try
            program_run = run(pipeline(`$program  $input_file_loc`, output_file_loc), wait=false)
            (terminated, status) = QE_run_watchdog(program_run, output_file_loc, program, watchdog_setting)
        catch _e_
            @error "QE_run() : Program $(join(program.exec,' ')) lauching error : $(_e_)"
        end
    else
        @error "QE_run() : Unknown mode $mode. Program $program not run. Exiting ..."
        ENV["ESPRESSO_TMPDIR"] = TMP_DIR0
        return results, success, status
    end

    if !terminated && program_run!==nothing
        kill(program_run)
    end

    QE_run_report(program, terminated, status)

    success = terminated && (status ∈ [ :JobDoneGoodFile, :JobKilledAfterGoodFile, 
                                        :JobEndedGoodFile, :JobSilentTerminatedGoodFile ])
    results = isfile(output_file_loc) ? readlines(output_file_loc) : String[]

    ENV["ESPRESSO_TMPDIR"] = TMP_DIR0
    return results, success, status

end


@inline function QE_run(
    program::Cmd, 
    scripts::Vector{S},
    workspace; 
    mode=:pipeline,
    f_in = "tmp.in", 
    f_out = "tmp.out",
    watchdog_setting = _default_watchdog_setting
    ) where {S<:AbstractString}

    QE_run( program,
            join(scripts,"\n"), 
            workspace;
            mode=mode,
            f_in=f_in, f_out=f_out,
            watchdog_setting=watchdog_setting )

end

#* -------------------------------------------------------------------
#* -------------------------------------------------------------------


# pw.x
function ld1(
    mpi::String,
    np::Int,
    map_by::String,
    npool::Int,
    scripts::Vector{S};
    workspace="./",
    fin="ld1.x.tmp.in", 
    fout="ld1.x.tmp.out",
    watchdog_setting = _default_watchdog_setting
    )::Tuple{Vector{String},Bool,Symbol} where {S<:AbstractString}
    if np==1 && mpi==""
        return QE_run( `ld1.x`, scripts, workspace, 
                        f_in=fin, f_out=fout, watchdog_setting=watchdog_setting )
    elseif mpi=="srun"
        return QE_run( `srun  ld1.x  -npool  $(npool)`, scripts, workspace,
                        f_in=fin, f_out=fout, watchdog_setting=watchdog_setting )
    else
        opt = map_by=="" ? [] : [map_by]
        if np > 1
            return QE_run( `mpiexec  -np  $np  $opt  ld1.x  -npool  $(npool)`, scripts, workspace,
                            f_in=fin, f_out=fout, watchdog_setting=watchdog_setting )
        else # np == -1, oversubscribe
            return QE_run( `mpiexec  --oversubscribe  $opt  ld1.x  -npool  $(npool)`, scripts, workspace,
                            f_in=fin, f_out=fout, watchdog_setting=watchdog_setting )
        end
    end
end


# pw.x dry run
function pw_dry(
    mpi::String,
    np::Int,
    map_by::String,
    npool::Int,
    scripts::Vector{S},
    stop_crit::Union{String,Regex};
    workspace="./",
    fin="pw.x.tmp.in", fout="pw.x.tmp.out",
    watchdog_setting = _default_intercepter_setting
    )::Tuple{Vector{String},Bool,Symbol} where {S<:AbstractString}
    if np==1 && mpi==""
        return QE_dry_run( `pw.x`, scripts, workspace, stop_crit,
                            f_in=fin, f_out=fout, watchdog_setting=watchdog_setting )
    
    elseif mpi=="srun"
        return QE_dry_run( `srun  pw.x  -npool  $(npool)`, scripts, workspace, stop_crit,
                            f_in=fin, f_out=fout, watchdog_setting=watchdog_setting )
    else
        opt = map_by=="" ? [] : [map_by]
        if np > 1
            return QE_dry_run( `mpiexec  -np  $np  $opt  pw.x  -npool  $(npool)`, scripts, workspace, stop_crit,
                                f_in=fin, f_out=fout, watchdog_setting=watchdog_setting )
        else # np == -1, oversubscribe
            return QE_dry_run( `mpiexec  --oversubscribe  $opt  pw.x  -npool  $(npool)`, scripts, workspace, stop_crit,
                                f_in=fin, f_out=fout, watchdog_setting=watchdog_setting )
        end
    end
end


# pw.x
function pw(
    mpi::String,
    np::Int,
    map_by::String,
    npool::Int,
    scripts::Vector{S};
    workspace="./",
    fin="pw.x.tmp.in", fout="pw.x.tmp.out",
    watchdog_setting = _default_watchdog_setting
    )::Tuple{Vector{String},Bool,Symbol} where {S<:AbstractString}
    if np==1 && mpi==""
        return QE_run( `pw.x`, scripts, workspace, 
                        f_in=fin, f_out=fout, watchdog_setting=watchdog_setting )
    elseif mpi=="srun"
        return QE_run( `srun  pw.x  -npool  $(npool)`, scripts, workspace,
                        f_in=fin, f_out=fout, watchdog_setting=watchdog_setting )
    else
        opt = map_by=="" ? [] : [map_by]
        if np > 1
            return QE_run( `mpiexec  -np  $np  $opt  pw.x  -npool  $(npool)`, scripts, workspace,
                            f_in=fin, f_out=fout, watchdog_setting=watchdog_setting )
        else # np == -1, oversubscribe
            return QE_run( `mpiexec  --oversubscribe  $opt  pw.x  -npool  $(npool)`, scripts, workspace,
                            f_in=fin, f_out=fout, watchdog_setting=watchdog_setting )
        end
    end
end


#* -------------------------------------------------------------------
#* -------------------------------------------------------------------


# cp.x dry run
function cp_dry(
    mpi::String,
    np::Int,
    map_by::String,
    npool::Int,
    scripts::Vector{S},
    stop_crit::Union{String,Regex};
    workspace="./",
    fin="cp.x.tmp.in", 
    fout="cp.x.tmp.out",
    watchdog_setting = _default_intercepter_setting
    )::Tuple{Vector{String},Bool,Symbol} where {S<:AbstractString}
    if np==1 && mpi==""
        return QE_dry_run( `cp.x`, scripts, workspace, stop_crit,
                            f_in=fin, f_out=fout, watchdog_setting=watchdog_setting )
    
    elseif mpi=="srun"
        return QE_dry_run( `srun  cp.x  -npool  $(npool)`, scripts, workspace, stop_crit,
                            f_in=fin, f_out=fout, watchdog_setting=watchdog_setting )
    else
        opt = map_by=="" ? [] : [map_by]
        if np > 1
            return QE_dry_run( `mpiexec  -np  $np  $opt  cp.x  -npool  $(npool)`, scripts, workspace, stop_crit,
                                f_in=fin, f_out=fout, watchdog_setting=watchdog_setting )
        else # np == -1, oversubscribe
            return QE_dry_run( `mpiexec  --oversubscribe  $opt  cp.x  -npool  $(npool)`, scripts, workspace, stop_crit,
                                f_in=fin, f_out=fout, watchdog_setting=watchdog_setting )
        end
    end
end


# cp.x
function cpx(
    mpi::String,
    np::Int,
    map_by::String,
    npool::Int,
    scripts::Vector{S};
    workspace="./",
    fin="cp.x.tmp.in", 
    fout="cp.x.tmp.out",
    watchdog_setting = _default_watchdog_setting
    )::Tuple{Vector{String},Bool,Symbol} where {S<:AbstractString}
    if np==1 && mpi==""
        return QE_run( `cp.x`, scripts, workspace, 
                        f_in=fin, f_out=fout, watchdog_setting=watchdog_setting )
    elseif mpi=="srun"
        return QE_run( `srun  cp.x  -npool  $(npool)`, scripts, workspace,
                        f_in=fin, f_out=fout, watchdog_setting=watchdog_setting )
    else
        opt = map_by=="" ? [] : [map_by]
        if np > 1
            return QE_run( `mpiexec  -np  $np  $opt  cp.x  -npool  $(npool)`, scripts, workspace,
                            f_in=fin, f_out=fout, watchdog_setting=watchdog_setting )
        else # np == -1, oversubscribe
            return QE_run( `mpiexec  --oversubscribe  $opt  cp.x  -npool  $(npool)`, scripts, workspace,
                            f_in=fin, f_out=fout, watchdog_setting=watchdog_setting )
        end
    end
end


#* -------------------------------------------------------------------
#* -------------------------------------------------------------------


# pp.x
function pp(
    scripts;
    workspace="./",
    fin="pp.x.tmp.in",
    fout="pp.x.tmp.out",
    watchdog_setting = _default_watchdog_setting
    )::Tuple{Vector{String},Bool,Symbol}
    return QE_run( `pp.x`, scripts, workspace, 
                    f_in=fin, f_out=fout, watchdog_setting=watchdog_setting )
end


# projwfc.x
function projwfc(
    mpi::String,
    np::Int,
    map_by::String,
    npool::Int,
    scripts::Vector{S};
    workspace = "./",
    fin       = "projwfc.x.tmp.in",
    fout      = "projwfc.x.tmp.out",
    watchdog_setting = _default_watchdog_setting
    )::Tuple{Vector{String},Bool,Symbol} where {S<:AbstractString}
    if np==1 && mpi==""
        return QE_run( `projwfc.x`, scripts, workspace, 
                        f_in=fin, f_out=fout, watchdog_setting=watchdog_setting )
    elseif mpi=="srun"
        return QE_run( `srun  projwfc.x  -npool  $(npool)`, scripts, workspace,
                        f_in=fin, f_out=fout, watchdog_setting=watchdog_setting )
    else
        opt = map_by=="" ? [] : [map_by]
        if np > 1
            return QE_run( `mpiexec  -np  $np  $opt  projwfc.x  -npool  $(npool)`, scripts, workspace,
                            f_in=fin, f_out=fout, watchdog_setting=watchdog_setting )
        else # np == -1, oversubscribe
            return QE_run( `mpiexec  --oversubscribe  $opt  projwfc.x  -npool  $(npool)`, scripts, workspace,
                            f_in=fin, f_out=fout, watchdog_setting=watchdog_setting )
        end
    end
end


#* -------------------------------------------------------------------
#* -------------------------------------------------------------------


# ph.x dry run
function ph_dry(
    mpi::String,
    np::Int,
    map_by::String,
    npool::Int,
    scripts::Vector{S},
    stop_crit::Union{String,Regex};
    workspace="./",
    fin="ph.x.tmp.in", fout="ph.x.tmp.out",
    watchdog_setting = _default_intercepter_setting
    )::Tuple{Vector{String},Bool,Symbol} where {S<:AbstractString}
    if np==1 && mpi==""
        return QE_dry_run( `ph.x`, scripts, workspace, stop_crit,
                            f_in=fin, f_out=fout, watchdog_setting=watchdog_setting )
    elseif mpi=="srun"
        return QE_dry_run( `srun  ph.x  -npool  $(npool)`, scripts, workspace, stop_crit,
                            f_in=fin, f_out=fout, watchdog_setting=watchdog_setting )
    else
        opt = map_by=="" ? [] : [map_by]
        if np > 1
            return QE_dry_run( `mpiexec  -np  $np  $opt  ph.x  -npool  $(npool)`, scripts, workspace, stop_crit,
                                f_in=fin, f_out=fout, watchdog_setting=watchdog_setting )
        else # np == -1, oversubscribe
            return QE_dry_run( `mpiexec  --oversubscribe  $opt  ph.x  -npool  $(npool)`, scripts, workspace, stop_crit,
                                f_in=fin, f_out=fout, watchdog_setting=watchdog_setting )
        end
    end
end


# ph.x
function ph(
    mpi::String,
    np::Int,
    map_by::String,
    npool::Int,
    scripts::Vector{S};
    workspace="./",
    fin="ph.x.tmp.in",
    fout="ph.x.tmp.out",
    watchdog_setting = _default_watchdog_setting
    )::Tuple{Vector{String},Bool,Symbol} where {S<:AbstractString}
    if np==1 && mpi==""
        return QE_run( `ph.x`, scripts, workspace,
                        f_in=fin, f_out=fout, watchdog_setting=watchdog_setting )
    elseif mpi=="srun"
        return QE_run( `srun  ph.x  -npool  $(npool)`, scripts, workspace, 
                        f_in=fin, f_out=fout, watchdog_setting=watchdog_setting )
    else
        opt = map_by=="" ? [] : [map_by]
        if np > 1
            return QE_run( `mpiexec  -np  $np  $opt  ph.x  -npool  $(npool)`, scripts, workspace, 
                            f_in=fin, f_out=fout, watchdog_setting=watchdog_setting )
        else # np == -1, oversubscribe
            return QE_run( `mpiexec  --oversubscribe  $opt  ph.x  -npool  $(npool)`, scripts, workspace, 
                            f_in=fin, f_out=fout, watchdog_setting=watchdog_setting )
        end
    end
end

#* -------------------------------------------------------------------
#* -------------------------------------------------------------------

# bands.x
function QE_bands(
    mpi::String,
    np::Int,
    map_by::String,
    scripts;
    workspace="./",
    fin="bands.x.tmp.in",
    fout="bands.x.tmp.out",
    watchdog_setting = _default_watchdog_setting
    )::Tuple{Vector{String},Bool,Symbol}
    if np==1 && mpi==""
        return QE_run(  `bands.x`, scripts, workspace, 
                        f_in=fin, f_out=fout, watchdog_setting=watchdog_setting )
    elseif mpi=="srun"
        return QE_run(  `srun  bands.x`, scripts, workspace, 
                                    f_in=fin, f_out=fout, watchdog_setting=watchdog_setting )
    else
        opt = map_by=="" ? [] : [map_by]
        if np > 1
            return QE_run(  `mpiexec  -np  $np  $opt  bands.x`, scripts, workspace, 
                            f_in=fin, f_out=fout, watchdog_setting=watchdog_setting )
        else # np == -1, oversubscribe
            return QE_run(  `mpiexec  --oversubscribe  $opt  bands.x`, scripts, workspace, 
                            f_in=fin, f_out=fout, watchdog_setting=watchdog_setting )
        end
    end
end


# bands.x
function QE_plotbands(
    scripts;
    workspace="./",
    fin="plotband.x.tmp.in",
    fout="plotband.x.tmp.out",
    watchdog_setting = _default_watchdog_setting
    )::Tuple{Vector{String},Bool,Symbol}
    return QE_run(  `plotband.x`, scripts, workspace, 
                                  f_in=fin, f_out=fout, watchdog_setting=watchdog_setting )
end

#* -------------------------------------------------------------------
#* -------------------------------------------------------------------

# q2r.x
function q2r(
    scripts;
    workspace="./",
    fin="q2r.x.tmp.in",
    fout="q2r.x.tmp.out",
    watchdog_setting = _default_watchdog_setting
    )::Tuple{Vector{String},Bool,Symbol}
    return QE_run(  `q2r.x`, scripts, workspace, 
                             f_in=fin, f_out=fout, watchdog_setting=watchdog_setting )
end


# matdyn.x
#: both matdyn.x and dynmat.x exist !
    #>> matdyn for real-space IFC produced by q2r
    #>> dynmat for ph.x output  
function matdyn(
    scripts;
    workspace=".",
    fin="matdyn.x.tmp.in",
    fout="matdyn.x.tmp.out",
    watchdog_setting = _default_watchdog_setting
    )::Tuple{Vector{String},Bool,Symbol}
    
    R = QE_run( `matdyn.x`, scripts, workspace, 
                            f_in=fin, f_out=fout, 
                            watchdog_setting=watchdog_setting )
    
    freq_fn = first([basename(last(SPLTEQ(l))) for l in scripts if occursin("flfrq",l)])
    
    add_white_space(freq_fn, workspace)
    
    return R
end

#* -------------------------------------------------------------------
#* -------------------------------------------------------------------

# pw2wannier90.x
function pw2wannier90(
    mpi::String,
    np::Int,
    map_by::String,
    scripts::Vector{S};
    workspace="./",
    fin="pw2wannier90.x.tmp.in",
    fout="pw2wannier90.x.tmp.out",
    watchdog_setting = _default_watchdog_setting
    )::Tuple{Vector{String},Bool,Symbol} where {S<:AbstractString}
    if np==1 && mpi==""
        return QE_run( `pw2wannier90.x`, scripts, workspace, 
                        f_in=fin, f_out=fout, watchdog_setting=watchdog_setting )
    elseif mpi=="srun"
        return QE_run( `srun  pw2wannier90.x`, scripts, workspace,
                        f_in=fin, f_out=fout, watchdog_setting=watchdog_setting )
    else
        opt = map_by=="" ? [] : [map_by]
        if np > 1
            return QE_run( `mpiexec  -np  $np  $opt  pw2wannier90.x`, scripts, workspace,
                            f_in=fin, f_out=fout, watchdog_setting=watchdog_setting )
        else # np == -1, oversubscribe
            return QE_run( `mpiexec  --oversubscribe  $opt  pw2wannier90.x`, scripts, workspace,
                            f_in=fin, f_out=fout, watchdog_setting=watchdog_setting )
        end
    end
end
