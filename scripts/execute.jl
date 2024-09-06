update_nothing = x->Dict()

global const prog_with_no_beta = [
    "wannier90.x", 
    "wannier90.x -pp", 
    "pw2wannier90.x", 
    "bands.x", 
    "plotbands.x",
    "projwfc.x",
    "cp.x",
    "cpdry.x", 
    "ph.x", 
    "phdry.x"
]


function execute_serial(
    S::INSTRUCTION;
    WORKSPACE=".",
    from_scratch=false
    )

    # 0. enter workspace
    WORKSPACE1 = try_mkdir(WORKSPACE, ".", "execute_serial()")
    cd(WORKSPACE1)
    WORKSPACE_FULL = pwd()

    # 1. set logger file path
    LOGGER_FILE_FULL_PATH_PREV = ENV["LOGGER_FILE_FULL_PATH"]
    ENV["LOGGER_FILE_FULL_PATH"] = "$WORKSPACE_FULL/execute_serial.log"

    @info log_info(">>> execute_serial() WORKSPACE=$WORKSPACE_FULL", level=0)

    # 2. compute on each plan
    P  = expand(S)
    RESULTS   = QEResult[]
    REGISTER     = Dict()  #TODO extension to "register"
    for (plan_name, program, config) ∈ P
        mpi, prog, npool, mode, np, maby = analyze_prog(program)

        #! the inner launch grid for each calculation : multiple beta, diag methods and cutoff
        GRID_OF_SETTINGS = []
        cutoff_scales    = get(config, :cutoff_scales, [1.0,])
        cutoffs          = collect(zip(
                                cutoff_scales.*get(config,"ecutwfc",100.0), 
                                cutoff_scales.*get(config,"ecutrho",400.0)
                           ))
        if prog ∈ prog_with_no_beta
            if prog ∈ ["ph.x", "phdry.x"] 
                #* deal with ph.x, which does not require "mixing_beta"
                #* but requires "diagonalization"
                decreasing_beta = [0.0,]
                diag_choices    = get(config, :diag_choices, [get(config,"diagonalization","david"),] )
                GRID_OF_SETTINGS = [Iterators.product(decreasing_beta, diag_choices, cutoffs) ...]
            else
                #* deal with programs that do not require "mixing_beta" nor "diagonalization"
                decreasing_beta = [0.0,]
                diag_choices    = ["dummy",]            
                GRID_OF_SETTINGS = [Iterators.product(decreasing_beta, diag_choices, cutoffs) ...]
            end
        else    #* requires "mixing_beta" ANDS "diagonalization"
            decreasing_beta = get(config, :beta,         [get(config,"mixing_beta",0.7),])
            diag_choices    = get(config, :diag_choices, [get(config,"diagonalization","david"),] )
            GRID_OF_SETTINGS = [Iterators.product(decreasing_beta, diag_choices, cutoffs) ...]
        end

        #! other settings
        FD   = try_mkdir( "$WORKSPACE_FULL/$plan_name", WORKSPACE_FULL, "execute_serial()" )
        #* !!!!!!!!!!!!!!!!! UPDATE !!!!!!!!!!!!!!!!!!!!
        conf = (Dict("outdir"=>FD) ⬱ config) ⬱ REGISTER
        #TODO how does "outdir" propagate ?
        #* !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

        @info log_info("execute_serial() \n\t\t\t\tWORKSPACE=$WORKSPACE_FULL\n\t\t\t\tplan_name=$plan_name\n\t\t\t\tprogram=$program\n\t\t\t\toutdir=$(conf["outdir"])\n\t\t\t\tdecreasing_beta=$(decreasing_beta)\n\t\t\t\tdiag_choices=$(diag_choices)", level=0)

        #! initialize recorder
        R  = QEResult( 
                plan_name, # PLAN::String
                program,   # PROG::Cmd
                conf, # CONF::Dict
                FD, # WORKSPACE::String
                String[], # INP::Vector{String}
                "$(config["prefix"]).$(prog).in", #! INP_F::String
                true,  # SUCC::Bool
                :NotStarted, # STAT::Symbol
                String[],  # OUTP::Vector{String}
                "$(config["prefix"]).$(prog).out", #! OUTP_F::String
                Dict(),  # RES::Dict
                get(config, :updater, update_nothing),  # UPDATER::Function
                get(config, :watchdog_setting, default_watchdog_setting), # WD::NamedTuple
                get(config, :from_scratch, from_scratch)  # SCRATCH::Bool  #TODO simplify control
        )
        #: the updater must be very robust, 
        #: because there's no internal guarantees for it to run correctly
        #: it all depends on the user !!!

        #! ------------ main loop ------------
        output_lines = String[]
        success = false
        status  = :NotStarted   
        for (β, diag_method, (wcut,ecut))  ∈ GRID_OF_SETTINGS
            @info  log_info("execute_serial(): beta=$β , diag_method=$diag_method", level=0)
            new_conf = conf
            #: ---------------
            try
                if prog ∉ prog_with_no_beta
                    new_conf = conf ⬱ Dict("mixing_beta"=>β, "diagonalization"=>diag_method, "ecutwfc"=>wcut, "ecutrho"=>ecut)
                end
                R.INP = generate_input_script(program, new_conf)
                if prog ∉ prog_with_no_beta
                    @info log_info("execute_serial(): finished generate_input_script(): mixing_beta=$(new_conf["mixing_beta"]) , diagonalization=$(new_conf["diagonalization"]) ", level=0)
                else
                    @info log_info("execute_serial(): finished generate_input_script() .", level=0)
                end
            catch _e_
                @error log_error("generate_input_script() gives error $(_e_)")
                rethrow(_e_)
                continue
            end
            #: ---------------
            try
                output_lines, success, status = lauch_program(R)
            catch _e_
                if prog ∉ prog_with_no_beta
                    @error log_error("lauch_program() \n\t\t\tbeta=$(new_conf["mixing_beta"])\n\t\t\tdiagonalization=$(new_conf["diagonalization"])\n\t\t\tgives error $(_e_)")
                else
                    @error log_error("lauch_program() \n\t\t\tgives error $(_e_)")
                end
                continue
            end
            @info log_info("execute_serial() finished from lauch_program() status=$success", level=0)
            if success  break  end
        end
        #! ------------ end of main loop ------------

        R.OUTP      = output_lines
        R.SUCC      = success
        R.STAT      = status
        if success
            R.RES   = interprete_results(program, output_lines)
        end
        push!(RESULTS, R)

        if !(success)
            @error log_fucked_up("Quiting execute_serial() because lauch_program() failed with all combinations in \n\t\t\t\tdecreasing_beta=$(decreasing_beta)\n\t\t\t\tdiag_choices=$(diag_choices)\n\t\t\t\tcutoffs=$(cutoffs)")
            break
        end

        #* !!!!!!!!!!! update !!!!!!!!!!!!
        REGISTER = R.UPDATER(output_lines)
        #* !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

        #* !!!!!!!!!!!! break loop !!!!!!!!!!!!!!!!
        if (:exit ∈ keys(REGISTER)) && REGISTER[:exit]
            @warn "execute_serial() : exit in R.UPDATER() triggered. \nExiting loop with REGISTER =  \n$(showdict(REGISTER))"
            break
        end
        #* !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    end

    # 3. restore and return
    ENV["LOGGER_FILE_FULL_PATH"] = LOGGER_FILE_FULL_PATH_PREV 
    @info log_info("<<< execute_serial() WORKSPACE=$WORKSPACE_FULL", level=0)
    return RESULTS
end

