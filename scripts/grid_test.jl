function grid_test(
    WORKSPACE::String, 
    common_title::String, 
    init_cif_fn::String, 
    instructor::Function,
    folder_name_translator::Function,
    arg_lists::Vector;
    cleanup=true
    )::Dict

    # 0. set log file path
    LOGGER_FILE_FULL_PATH_PREV = ENV["LOGGER_FILE_FULL_PATH"]
    ENV["LOGGER_FILE_FULL_PATH"] = "$WORKSPACE/grid_test.log"

    # 1. entering folder
    cd(WORKSPACE)

    # 2. load initial structure 
    init_struct = get_structure_from_cif(init_cif_fn)
    @info log_info("grid_test() : init_struct\n\t\tinit_struct[:cif] = $(init_struct[:cif])\n\t\tinit_struct[\"positions\"] = $(init_struct["positions"])", level=0)

    # 3. folders and settings
    folders_and_settings = Dict( folder_name_translator(v...)=>(v...,) 
                                 for v ∈ Iterators.product(arg_lists...) )
    Nsettings = length(folders_and_settings)
    dN = max(Nsettings ÷ 3, 2)

    # 4. serial excution
    result_v = []
    _cnt_ = 0
    for (fd, config_i) ∈ folders_and_settings
        @info log_info("grid_test() : \n\t\tfd = $fd\n\t\tconfig_i = $config_i", level=0)
        #% report 
        _cnt_ += 1
        if (_cnt_ % dN) == 0
            grid_test__REPORT(
                WORKSPACE, 
                common_title, 
                init_cif_fn, 
                instructor, 
                folder_name_translator, 
                arg_lists
            )
        end
        #% =================
        #%     main part    
        #% =================
        ins = nothing
        instr_generated = false
        try 
            ins = instructor(init_struct, fd, config_i...)
            instr_generated = true
        catch _e_
            @error log_error("grid_test() : instructor(...) error : \n\t\tfd=$fd \n\t\terror $(_e_)\n\t\tR=execute_serial(...) skipped.", level=1)
            rethrow(_e_)
        end
        # run 
        R = QEResult[]
        if instr_generated
            #TODO add control of `from_scratch`
            #TODO currently is just a safe choice
            R = execute_serial(ins; WORKSPACE="$WORKSPACE/$fd", from_scratch=false)
        end
        # results
        push!(result_v, fd=>(config_i,R))
        # clean up, only if instructed and all calculations are successful
        if cleanup && check_status_all_successful(ins, "$WORKSPACE/$fd")
            grid_test__CLEANUP(ins, WORKSPACE="$WORKSPACE/$fd")
        end
    end
    RESULTS = Dict(result_v)

    # 5. report and save
    #% report 
    grid_test__REPORT(
        WORKSPACE, 
        common_title, 
        init_cif_fn, 
        instructor, 
        folder_name_translator, 
        arg_lists
    )
    #@save "$WORKSPACE/grid_test_$(common_title).all.jld2" RESULTS

    # 6. return
    ENV["LOGGER_FILE_FULL_PATH"] = LOGGER_FILE_FULL_PATH_PREV
    return RESULTS
end


"""
    function grid_test__STATUS(
        WORKSPACE::String, 
        common_title::String,   #! not used
        init_cif_fn::String, 
        instructor::Function,
        folder_name_translator::Function,
        arg_lists::Vector
        )


###NOTE

This function checks each folder in the `WORKSPACE`  
according to all combinations of `arg_lists`, whether 
the instruction `ins` has been completed. 

The function `check_status(ins, "WORKSPACE/fd")` does the 
main job. 
"""
function grid_test__STATUS(
    WORKSPACE::String, 
    common_title::String,   #! not used
    init_cif_fn::String, 
    instructor::Function,
    folder_name_translator::Function,
    arg_lists::Vector
    )::Vector
    cd(WORKSPACE)
    init_struct = get_structure_from_cif(init_cif_fn)

    #> -------- IDENTCAL PART  ---------
    #> folders and settings
    folders_and_settings = 
        Dict( folder_name_translator(v...)=>(v...,) 
              for v ∈ Iterators.product(arg_lists...) )

    #> ---------------------------------
    #> serial excution
    RESULTS = []
    for (fd, config_i) ∈ folders_and_settings
        ins = nothing
        try
            ins = instructor(init_struct, fd, config_i...)
        catch _e_
            nothing
        end
        status_list = check_status(ins, "$WORKSPACE/$fd")
        push!(RESULTS, fd=>(config_i, status_list))
    end

    #> ---------------------------------
    #> save
    return RESULTS
end

"""
    grid_test__STATUS_print(status0)


###NOTE 

This function prints the results `status0` of `grid_test__STATUS()`.

"""
function grid_test__STATUS_print(
    status0;
    additional_info = x->[]
    )
    @inline spltend(x) = split(x,"_",keepempty=false)[end]
    status = sort(status0, by=first)  # sort by key `fd`
    for (fd, content) ∈ status
        (config_i, status_list) = content
        _L_ = []
        _O_ = []
        for (A,B) in status_list
            (plan_name,program) = A
            #*[crossref] = grid_test__STATUS_print
            (status,FD,INP_F,OUTP_F) = B
            N = spltend(plan_name)
            push!(_L_, N => status)
            push!(_O_, N => OUTP_F)
        end
        # print
        println(join(string.(last.(_L_)), " ") * "    " * fd)
        println(join(additional_info(last.(_O_)), "\n"))
    end
    return
end


#* ======================================================

# example of additional_info
# additional_info = x -> check_ph_negative_freq(last(x))
function check_ph_negative_freq(ph_OUTP_F)
    neg_freq = []
    if !isfile(ph_OUTP_F)
        return []
    end
    try 
        q_f      = ph_frequency_THz(readlines(ph_OUTP_F))
        neg_freq = [q=>f for (q,f) ∈ q_f if minimum(f)<0]
    catch _e_
        return ["check_ph_negative_freq() error : \n $(_e_)",]
    end
    return ( length(neg_freq)==0 
             ? []
             : [ "Negative frequencies : ",
                 ["q = $q \nf = $nf" for (q,nf) ∈ neg_freq]... ] )
end

#* ======================================================


function grid_test__REPORT(
    WORKSPACE::String, 
    common_title::String,
    init_cif_fn::String, 
    instructor::Function,
    folder_name_translator::Function,
    arg_lists::Vector;
    status="completed",
    additional_info = x->[]
    )
    println("\n---------------------------------")
    println("Status Report ($status)")
    println("WORKSPACE : $WORKSPACE")
    println("common_title : $common_title")
    println("---------------------------------\n")
    grid_test__STATUS_print(
        grid_test__STATUS(
            WORKSPACE, 
            common_title,
            init_cif_fn,
            instructor,
            folder_name_translator,
            arg_lists
        );
        additional_info = additional_info
    )
    println("---------------------------------\n")
    return
end