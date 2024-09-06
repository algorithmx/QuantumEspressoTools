include("/home/dabajabaza/jianguoyun/Workspace/QuantumEspressoTools/src/watchdog.jl")

job_done(lines) = (@show lines ; findfirst(x->x=="21",lines)!==nothing)

job_done2(lines) = findfirst(x->x=="200",lines)!==nothing

function xxx(file, N=101, s=0.9)
    for i=1:N
        open(file, "a") do io
            write(io, "$i\n")
        end
        sleep(s)
    end
    return
end

function reset(f)
    open(f, "w") do io
        write(io, "")
    end
    return
end

##* ====================================

FN = "/home/dabajabaza/input_tmp.in"
reset(FN)

##* ====================================

for mv in [5,4,3,2]
    a1() = xxx(FN, 10, 1);
    w1(bb) = watchdog_serial(   bb, FN, job_done, 
                                max_watch=100, woof_per_x_min=1/120, max_silence_woofs=mv, tail_length=3  )
    #b1 = @async a1(); x1 = w1(b1); @show x1; reset(FN)
end

##* ====================================

for mv in [5,4,3,2]
    a2() = xxx(FN, 30, 0.5);
    w2(bb) = watchdog_serial(   bb, FN, job_done, 
                                max_watch=100, woof_per_x_min=1/240, max_silence_woofs=mv, tail_length=3  )
    #b2 = @async a2(); x2 = w2(b2); @show x2; reset(FN)
end

##* ====================================

for mw in [5,4,3,2]
    a3() = xxx(FN, 20, 0.1);
    w3(bb) = watchdog_serial( bb, FN, job_done, 
                            max_watch=100, woof_per_x_min=1/1200, max_silence_woofs=mw, tail_length=3 )
    #b3 = @async a3(); x3 = w3(b3); @show x3; reset(FN)
end

##* ====================================

#* sampling rate too low, missed the "21 line"
for tl in [5,4,3,2]
    a4() = xxx(FN, 40, 0.2);
    w4(bb) = watchdog_serial( bb, FN, job_done, 
                              max_watch=100, woof_per_x_min=1/60, max_silence_woofs=5, tail_length=tl )
    b4 = @async a4(); x4 = w4(b4); @show x4; reset(FN)
end

##* ====================================

