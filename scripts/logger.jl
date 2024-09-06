# logger.jl
using Dates


function timestamp()
    return Dates.format(now(), "Y u d HH:MM:SS")
end


function logger(file_full_path::String, title::String, msg::String; level=0)::String

    t = timestamp()
    prep = join((level==0 ? [""] : (level==1 ? ["|---",] : ["|---", ["----" for i=2:level]...])), "")
    to_write = "[$t] <$title> $msg\n"

    log_io = undef
    try 
        log_io = open(file_full_path, "a")
    catch _
        log_io = open(file_full_path, "w")
    end
    for i=1:10
        try 
            write(log_io, prep * to_write)
            sleep(0.1rand())
            flush(log_io)
            break
        catch _
            sleep(0.1rand())
            continue
        end
    end
    for i=1:10
        try 
            close(log_io)
            break
        catch _
            sleep(0.1rand())
            continue
        end
    end

    return msg
end


function logger(title::String, msg::String; level=0)

    if "LOGGER_FILE_FULL_PATH" ∉ keys(ENV)
        @warn "logger.jl : PLEASE SPECIFY THE ENV VAR \"LOGGER_FILE_FULL_PATH\" !!!!!!"
    end

    file_full_path = get(ENV, "LOGGER_FILE_FULL_PATH", "log.txt")

    logger(file_full_path, title, msg; level=level)

end


log_info(msg::String; level=0) = logger("info", msg; level=level)

log_warn(msg::String; level=0) = logger("warn", msg; level=level)

log_error(msg::String; level=0) = logger("error", msg; level=level)

log_fucked_up(msg::String) = logger("FUUUUCKED UUUUP !!!!", msg; level=0)


