function verify_GULP_result(
    results::Vector{S}, 
    program::Cmd; 
    final=false
    ) where {S<:AbstractString}
    return true # ("JOB DONE" ⊂ results[max(1,length(results)-6):end])
end


function verify_QE_result(
    results::Vector{S}, 
    program::Cmd; 
    final=false
    ) where {S<:AbstractString}
    if any([(p ⊂ program.exec) for p ∈ __QE_PROGS__])
        status = ("JOB DONE" ⊂ results[max(1,length(results)-6):end])
        if !status && final
            @error join(["QE ERROR:",results[max(1,length(results)-6):end]...], "\n")
        end
        return status
    else
        return true
    end
end

function verify_wannier90_result(
    results::Vector{S}, 
    prog::Cmd; 
    final=false
    ) where {S<:AbstractString}
    status = ((("All done"     ⊂ results[max(1,length(results)-6):end]))
          || (("-pp"∈prog.exec)&&("nnkp written" ⊂ results[max(1,length(results)-2):end])))
    if !status && final
        @error join(["WANNIER90 ERROR:",results[max(1,length(results)-6):end]...], "\n")
    end
    return status
end

function verify_result(
    results::Vector{S}, 
    prog::Cmd; 
    final=false
    ) where {S<:AbstractString}
    pp = prog.exec
    if "gulp" ∈ pp
        return verify_GULP_result(results, prog; final=final)
    elseif "wannier90.x" ∈ pp
        return verify_wannier90_result(results, prog; final=final)
    else
        return verify_QE_result(results, prog; final=final)
    end
end


verify_MPI_ABORT(results::Vector{S}) where {S<:AbstractString} = 
        ("MPI_ABORT" ⊂ results[end-min(length(results)-1,10):end])

export verify_result, verify_MPI_ABORT