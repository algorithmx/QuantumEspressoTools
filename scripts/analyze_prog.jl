function analyze_prog(p::Cmd)::Tuple{String, String, Int, Symbol, Int, String}
    exec = p.exec
    i_npool  = findfirst(x->x=="-npool",exec)
    if (exec[1]=="mpiexec" || exec[1]=="mpirun")
        i_map_by  = findfirst(x->occursin("--map-by",x),exec)
        i_oversub = findfirst(x->occursin("--oversubscribe",x),exec)
        i_np      = findfirst(x->x=="-np",exec)
        if i_np!==nothing
            @inline maxn(a,b) = (a===nothing ? b+2 : (a>b ? a+1 : b+2))
            return (    exec[1],
                        exec[maxn(i_map_by,i_np)],  # prog
                        (i_npool===nothing ? 1 : parse(Int,exec[i_npool+1])),
                        (("-i" ∈ exec) ? :file : :pipeline), 
                        parse(Int,exec[i_np+1]), 
                        (i_map_by===nothing ? "" : exec[i_map_by])   )
        else ## --oversubscribe
            @assert i_oversub!==nothing
            @inline maxn2(a,b) = (a===nothing ? b+1 : (a>b ? a+1 : b+1))
            return (    exec[1],        #* mpi prog
                        exec[maxn2(i_map_by,i_oversub)],   #* prog
                        (i_npool===nothing ? 1 : parse(Int,exec[i_npool+1])),  #* npool
                        (("-i" ∈ exec) ? :file : :pipeline),  #* -i tag
                        -1,    #* np = -1 means oversubscribe !
                        (i_map_by===nothing ? "" : exec[i_map_by])   )
        end
    elseif exec[1]=="srun"
        return (exec[1], exec[2], (i_npool===nothing ? 1 : parse(Int,exec[i_npool+1])), (("-i" ∈ exec) ? :file : :pipeline), 1, "")
    else
        return ("", exec[1], (i_npool===nothing ? 1 : parse(Int,exec[i_npool+1])), (("-i" ∈ exec) ? :file : :pipeline), 1, "")
    end
end


##

#=

PROG_PWX    = `mpiexec -np 64 --map-by=slot pw.x -npool 8`
PROG_PWX    = `mpirun --map-by=slot --oversubscribe  pw.x -npool 8`
PROG_PW2W90 = `mpiexec -np 64 --map-by=slot pw2wannier90.x`
PROG_W90    = `mpiexec -np 64 --map-by=slot wannier90.x`
PROG_W90PP  = `mpiexec -np 64 --map-by=slot wannier90.x -pp`

analyze_prog(PROG_PWX)
analyze_prog(PROG_PW2W90)
analyze_prog(PROG_W90)
analyze_prog(PROG_W90PP)

=#