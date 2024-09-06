
global const default_watchdog_setting = (woof_per_x_min=1/10, max_watch=14400, max_silence_woofs=50, tail_length=30, quiet=false)
global const default_intercepter_setting = (woof_per_x_min=1/1000, max_watch=3000, max_silence_woofs=2000, tail_length=300, quiet=false)


function lauch_program(
    p::Cmd,
    i::Vector{S};
    workspace="",
    fin="",
    fout="",
    watchdog_setting = default_watchdog_setting,
    from_scratch=true
    )::Tuple{Vector{String},Bool,Symbol} where {S<:AbstractString}

    @inline input_filename_root(fn) = (endswith(fn,"wannier90.x.in") 
                                        ? string(replace(basename(fn),".wannier90.x.in"=>"")) 
                                        : fn)
    #input_filename_root("tmp..in")

    # 1. enter workspace
    @info log_info("Now in lauch_program() with workspace=$workspace", level=1)
    cd(workspace)

    # 2.1. analyze program parameters
    mpi, prog, npool, mode, np, maby = analyze_prog(p)

    # 2.2. prepare files 
    fin1  = fin==""  ? "$(prog).tmp.in"  : fin
    fout1 = fout=="" ? "$(prog).tmp.out" : fout

    # 3. check output file, if exist, return without calculation
    if isfile("$workspace/$fout1") && !from_scratch
        lines = readlines("$workspace/$fout1")
        if verify_result(lines, p)
            @info log_info("lauch_program() found good file $workspace/$fout1.", level=2)
            return lines, true, :FoundGoodFile
        end
    end

    # 4. launch the program indicated by `prog`
    __results__ = nothing
    if prog=="pwdry.x"
        @info log_info("lauch_program() starting dry run of pw($(np),$(maby)...,workspace=$workspace)", level=2)
        __results__ =  pw_dry(
                    mpi, np, maby, npool, i, r"Dense\s+grid";
                    workspace=workspace, fin=fin1, fout=fout1,
                    watchdog_setting = watchdog_setting )
    elseif prog=="pw.x"
        @info log_info("lauch_program() starting pw($(np),$(maby)...,workspace=$workspace)", level=2)
        __results__ =  pw(
                    mpi, np, maby, npool, i;
                    workspace=workspace, fin=fin1, fout=fout1,
                    watchdog_setting = watchdog_setting )
    elseif prog=="cpdry.x"
        @info log_info("lauch_program() starting dry run of cpx($(np),$(maby)...,workspace=$workspace)", level=2)
        __results__ =  cp_dry(
                    mpi, np, maby, npool, i, r"Dense\s+grid";  #TODO when to terminate ?
                    workspace=workspace, fin=fin1, fout=fout1,
                    watchdog_setting = watchdog_setting )
    elseif prog=="cp.x"
        @info log_info("lauch_program() starting cpx($(np),$(maby)...,workspace=$workspace)", level=2)
        __results__ =  cpx(
                    mpi, np, maby, npool, i;
                    workspace=workspace, fin=fin1, fout=fout1,
                    watchdog_setting = watchdog_setting )
    elseif prog=="phdry.x"
        @info log_info("lauch_program() starting dry run of ph($(np),$(maby)...,workspace=$workspace)", level=2)
        __results__ =  ph_dry(
                    mpi, np, maby, npool, i, r"Dense\s+grid";
                    workspace=workspace, fin=fin1, fout=fout1,
                    watchdog_setting = watchdog_setting )
    elseif prog=="ph.x"
        @info log_info("lauch_program() starting ph($(np),$(maby)...,workspace=$workspace)", level=2)
        __results__ =  ph(  
                    mpi, np, maby, npool, i;
                    workspace=workspace, fin=fin1, fout=fout1,
                    watchdog_setting = watchdog_setting )
    elseif prog=="pp.x"
        @info log_info("lauch_program() starting pp(...,workspace=$workspace)", level=2)
        __results__ =  pp(  i;
                    workspace=workspace, fin=fin1, fout=fout1,
                    watchdog_setting = watchdog_setting )
    elseif prog=="bands.x"
        @info log_info("lauch_program() starting QE_bands($(np),$(maby)...,workspace=$workspace)", level=2)
        __results__ =  QE_bands(   
                    mpi, np, maby, i;
                    workspace=workspace, fin=fin1, fout=fout1,
                    watchdog_setting = watchdog_setting )
    elseif prog=="plotband.x"
        @info log_info("lauch_program() starting QE_plotbands(...,workspace=$workspace)", level=2)
        __results__ = QE_plotbands(
                    i;
                    workspace=workspace, fin=fin1, fout=fout1,
                    watchdog_setting = watchdog_setting )
    elseif prog=="q2r.x"
        @info log_info("lauch_program() starting q2r(...,workspace=$workspace)", level=2)
        __results__ =  q2r( i;
                    workspace=workspace, fin=fin1, fout=fout1,
                    watchdog_setting = watchdog_setting )
    elseif prog=="matdyn.x"
        @info log_info("lauch_program() starting matdyn(...,workspace=$workspace)", level=2)
        __results__ =  matdyn( i;
                    workspace=workspace, fin=fin1, fout=fout1,
                    watchdog_setting = watchdog_setting )
    elseif prog=="ld1.x"
        @info log_info("lauch_program() starting ld1($(np),$(maby)...,workspace=$workspace)", level=2)
        __results__ =  ld1(
                    mpi, np, maby, npool, i;
                    workspace=workspace, fin=fin1, fout=fout1,
                    watchdog_setting = watchdog_setting )
    elseif prog=="pw2wannier90.x"
        @info log_info("lauch_program() starting pw2wannier90($(np),$(maby)...,workspace=$workspace)", level=2)
        __results__ =  pw2wannier90(
                    mpi, np, maby, i;
                    workspace=workspace, fin=fin1, fout=fout1,
                    watchdog_setting = watchdog_setting )
    elseif prog=="wannier90.x"
        if "-pp" ∈ p.exec
            @info log_info("lauch_program() starting wannier90pp($(np),$(maby)...,workspace=$workspace)", level=2)
            __results__ =  wannier90pp(
                        i;
                        workspace=workspace, 
                        f_root=input_filename_root(fin),
                        watchdog_setting = watchdog_setting )
        else
            @info log_info("lauch_program() starting wannier90($(np),$(maby)...,workspace=$workspace)", level=2)
            __results__ =  wannier90(
                        i;
                        workspace=workspace, 
                        f_root=input_filename_root(fin),
                        watchdog_setting = watchdog_setting )
        end
    elseif prog=="projwfc.x"
        @info log_info("lauch_program() starting projwfc($(np),$(maby)...,workspace=$workspace)", level=2)
        __results__ =  projwfc(
                    mpi, np, maby, npool, i;
                    workspace=workspace, fin=fin1, fout=fout1,
                    watchdog_setting = watchdog_setting )
    else
        @error log_fucked_up("lauch_program(): Unknown program $prog .")
        throw(error("lauch_program(): Unknown program $prog ."))
    end

    # 5. return
    @info log_info("Now exiting lauch_program() with workspace=$workspace", level=1)
    return __results__

end



function lauch_program(R::QEResult) 

    output_lines, success, status = lauch_program(
        R.PROG,
        R.INP, 
        workspace = R.WORKSPACE, 
        watchdog_setting = R.WD, 
        fin  = R.INP_F,
        fout = R.OUTP_F,
        from_scratch = R.SCRATCH
    )

    ##********************* patch *********************
    cnt = 0
    while verify_MPI_ABORT(output_lines) && cnt < 10
        @info log_info("lauch_program(): MPI_ABORT triggered, relaunch, i = $cnt", level=1)
        #% relaunch
        cnt += 1
        sleep(10rand())
        output_lines, success, status = lauch_program(
            R.PROG,
            R.INP, 
            workspace = R.WORKSPACE, 
            watchdog_setting = R.WD, 
            fin  = R.INP_F,
            fout = R.OUTP_F,
            from_scratch = R.SCRATCH
        )
    end
    if verify_MPI_ABORT(output_lines) && cnt >= 10
        @error log_error("lauch_program(): MPI_ABORT triggered, relaunched $cnt times, all failed.")
    end
    ##*************************************************

    return (output_lines, success, status)

end
