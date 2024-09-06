#TODO  allow for configurations on file types to be deleted

global const __DELETE_CMDS__ = [ 
    `find . -name '*.wfc*' -delete`,
    `find . -name '*.dat' -delete`,
    `find . -name '*.dwf*' -delete`,
    `find . -name '*.mixd*' -delete`,
    `find . -name '*.mix*' -delete`,
    `find . -name '*.bar*' -delete`,
    `find . -name '*.recover*' -delete`,
]


function grid_test__CLEANUP(
    S::INSTRUCTION;
    WORKSPACE="."
    )

    # 1. enter workspace
    if !isdir(WORKSPACE)
        @error "grid_test_cleanup() : \nWORKSPACE not found : \n$WORKSPACE"
        return 
    end
    cd(WORKSPACE)
    WORKSPACE_FULL = pwd()

    # 2. clean up each plan
    P  = expand(S)
    for (plan_name, program, config) ∈ P
        mpi, prog, npool, mode, np, maby = analyze_prog(program)
        FD  = "$WORKSPACE_FULL/$plan_name"
        if isdir(FD)
            fd_now = pwd()
            cd(FD)
            for c in __DELETE_CMDS__
                try
                    run(c)
                catch 
                    nothing
                end
            end
            cd(fd_now)
        end
    end

    # 3.
    return
end

