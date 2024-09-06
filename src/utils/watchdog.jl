__terminated__(x::Task) = istaskdone(x)

__terminated__(x::Union{Base.Process,Base.ProcessChain}) = process_exited(x)

__abnormal__(f) = false

function ∞_kill(prc::Union{Base.Process,Base.ProcessChain})
    @warn "∞ kill !!!"
    cnt = 0
    while !__terminated__(prc) && cnt<10000
        #@warn "$cnt kill !!!"
        try
            kill(prc)
            kill(prc)
            kill(prc)
            kill(prc)
            kill(prc)
        catch
            nothing
        end
        if __terminated__(prc)
                return
        else
            cnt += 1
            continue
        end
    end
    #TODO deal with catastrophic failure !!
    return
end


function ∞_kill(tsk::Task)
    @warn "∞ kill !!!"
    ex = InterruptException()
    cnt = 0
    while !__terminated__(tsk) && cnt<10000
        #@warn "$cnt kill !!!"
        try
            Base.throwto(tsk, ex)
        catch  _e_
            if (_e_ isa InterruptException) && __terminated__(tsk)
                return
            else
                cnt += 1
                continue
            end
        end
    end
    #TODO deal with catastrophic failure !!
    return
end


@inline tail(file, k=3) = (try lines=readlines(file); lines[max(1,length(lines)-k+1):end] catch; repeat(["",],k) end)


function watchdog_serial(
    b::Union{Task, Base.Process, Base.ProcessChain}, 
    output_file_of_b::String,
    is_job_done::Function;

    #! Hi task b ! I woof every 1 min.
    woof_per_x_min = 1,

    #! Hi task b ! I woof 7200 times.
    #! If you keep working while I'm watching you, 
    #! I'll leave you working FOREVER !!!
    max_watch = 7200,

    #! Hi task b ! 
    #! If I woof 5 times and you don't update the output_file_of_b, 
    #! I'LL KILL YOU !!!
    #! If you survived, you're good to work FOREVER !!!
    max_silence_woofs = 20,

    #! Hi task b ! I check your tail at a length of 50
    tail_length = 50,

    #! I'll wooooooooooooooooooof
    quiet = false,

    #! I allow you to terminate gracefully after 
    #! the output file has been judged "GoodFile"
    grace_time_sec = 600,

    woof_per = 5

    )::Tuple{Bool,Symbol}

    @inline all_same(x::Vector) = length(unique(x))==1

    @inline same_tail(x::Vector, nt::Int) = (length(x)>nt && all_same(x[end-(nt-1):end]))

    if !quiet
        @info join(["", " -------- WATCHDOG V0.2 --------",
                    "",
                    " Author : github.com/algorithmx  ",
                    "",
                    "                     |\\",
                    "            \\`-. _.._| \\",
                    "             |_,'  __`. \\",
                    "             (Θ\\ _/Θ| _  |",
                    "  WOOOF!    ,'      __ \\ |",
                    "          ,'     __/||\\  |",
                    "         (XX   ,/||||||/ |",
                    "            `-'_----    /",
                    "               /`-._.-'/",
                    "               `-.__.-'",
                    "-------------------------------"], "\n")
    end

    #IMPORTANT wait 1min for the file to be created
    cnt = 0
    while !isfile(output_file_of_b) && cnt<60
        sleep(1)  #> wait 20 seconds for the program to create a file
        cnt += 1
    end
    if cnt>=60
        throw(error("File $output_file_of_b not created since the program has started for 20 seconds."))
    end

    @assert isfile(output_file_of_b)
    hist = []
    k = 0
    while k < max_watch

        #> see what's in $file
        t = tail(output_file_of_b, tail_length)
        push!(hist, t)

        #! almost sure that if we have good file, watchdog() will exit here
        if is_job_done(t)
            __w__ = 0
            while (__w__< grace_time_sec) && (! __terminated__(b))
                sleep(1)
                __w__ += 1
            end
            if ! __terminated__(b)
                #: sorry b ....
                ∞_kill(b)
                return (__terminated__(b), :JobKilledAfterGoodFile)
            else
                return (__terminated__(b), :JobDoneGoodFile)
            end
        elseif __terminated__(b)
            if is_job_done(t)
                #! marginally possible because the last time for is_job_done(t)==false is just above
                return (true, :JobEndedGoodFile)
            else
                return (true, :JobEndedBadFile)
            end
        elseif same_tail(hist, max_silence_woofs)
            @info "Trying to terminate the silent process ... "
            ∞_kill(b)
            if __terminated__(b)  @info "The silent process has been terminated."  end

            if is_job_done(t)  
                #! not possible because when you fuck up a job, its output file is never complete
                return (true, :JobSilentTerminatedGoodFile)
            else
                return (true, :JobSilentTerminatedBadFile)
            end
        elseif __abnormal__(output_file_of_b)
            @info "Trying to terminate the abnormal process ... "
            ∞_kill(b)
            if __terminated__(b)  @info "The abnormal process has been terminated."  end
            return (true, :JobAbnormalTerminatedBadFile)
        end

        #% Wxxxooooooooooooof !
        sleep(woof_per_x_min * 60)
        if (!quiet && k%woof_per==0)
            num_woofs = k ÷ woof_per + 1
            num_o = num_woofs % 10 
            num_x = num_woofs ÷ 10 
            @info "W$(repeat('x',num_x))$(repeat('o',num_o))f !" 
        end
        k += 1

    end

    #! max_watch tooo short
    if !__terminated__(b)  ∞_kill(b)  end
    return (__terminated__(b), :RunTimeTooLong)

end


global const _behaviour_tranlator_ = Dict(
    :JobDoneGoodFile             => :JobDoneCritSat,
    :JobKilledAfterGoodFile      => :JobKilledAfterCritSat,
    :JobEndedGoodFile            => :JobEndedCritSat,
    :JobEndedBadFile             => :JobEndedCritNotSat,
    :JobSilentTerminatedGoodFile => :JobSilentTerminatedCritSat,
    :JobSilentTerminatedBadFile  => :JobSilentTerminatedCritNotSat,
    :RunTimeTooLong              => :RunTimeTooLong,
    :JobAbnormalTerminatedBadFile=> :JobAbnormal
)


function watchdog_intercepter(
    b::Union{Task, Base.Process, Base.ProcessChain}, 
    output_file_of_b::String,
    stop_crit::Union{String,Regex};
    woof_per_x_min = 1/120,
    max_watch = 240,
    max_silence_woofs = 120,
    tail_length = 300,
    quiet = true,
    woof_per = 20,
    )::Tuple{Bool,Symbol}

    intercept_result = watchdog_serial( b, output_file_of_b, x->(stop_crit ⊂ x);
                                        woof_per_x_min = woof_per_x_min,
                                        max_watch = max_watch,
                                        max_silence_woofs = max_silence_woofs,
                                        tail_length = tail_length,
                                        quiet = quiet,
                                        grace_time_sec = -1,
                                        woof_per = woof_per  )

    return (intercept_result[1], _behaviour_tranlator_[intercept_result[2]])

end
