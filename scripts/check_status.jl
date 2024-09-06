#TODO

"""
    check_status(workspace0)

return "DONE" / "WORKING" / "FAILED" for calculation in `workspace0`

"""
function check_status(
    S::Union{INSTRUCTION,Nothing},
    WORKSPACE_FULL
    )
    if S === nothing
        @error "check_status() : \nS === nothing in WORKSPACE \n$WORKSPACE_FULL"
        return []
    end
    P  = expand(S)
    RESULTS   = []
    for (plan_name, program, config) ∈ P
        mpi, prog, npool, mode, np, maby = analyze_prog(program)
        FD  = "$WORKSPACE_FULL/$plan_name"
        status = "?"
        if "wannier90.x" ∉ program.exec
            INP_F  = "$(config["prefix"]).$(prog).in"
            OUTP_F = "$(config["prefix"]).$(prog).out"
            if isfile("$FD/$INP_F") && isfile("$FD/$OUTP_F")
                if verify_result(readlines("$FD/$OUTP_F"), program; final=true)
                    status = "O"
                else
                    status = "X"
                end
            else
                status = "_"
            end
            #*[crossref] = grid_test__STATUS_print
            push!( RESULTS, 
                   (plan_name,program) => (status, FD, "$FD/$INP_F", "$FD/$OUTP_F"))
        else  # "wannier90.x" is S.P.E.C.I.A.L.
            INP_F  = "$(config["prefix"]).win"
            OUTP_F = "$(config["prefix"]).wout"
            if isfile("$FD/$INP_F") && isfile("$FD/$OUTP_F")
                if verify_result(readlines("$FD/$OUTP_F"), program; final=true)
                    status = "O"
                else
                    status = "X"
                end
            else
                status = "_"
            end
            #*[crossref] = grid_test__STATUS_print
            push!( RESULTS, 
                   (plan_name,program) => (status, FD, "$FD/$INP_F", "$FD/$OUTP_F"))
        end
    end
    return RESULTS
end


function check_status_all_successful(
    S::Union{INSTRUCTION,Nothing},
    WORKSPACE
    )
    STT = check_status(S, WORKSPACE)
    return length(STT)>0 && all(first.(last.(STT)) .== "O")
end