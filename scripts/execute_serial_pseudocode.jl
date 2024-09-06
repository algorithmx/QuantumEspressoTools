
function execute_serial(
    S::INSTRUCTION;
    WORKSPACE=".",
    from_scratch=false
    )
    RESULTS   = []
    REGISTER  = Dict()
    for (title, program, config) ∈ expand(S)
        folder = create_folder_in_workspace(title)
        confoguration = (Dict("outdir"=>folder) << config) << REGISTER
        updater = confoguration[:updater]
        input_scripts = generate_input_script(program, new_conf)
        output_lines, success, status = 
            lauch_program(program, input_scripts, from_scratch)
        results  = interprete_results(program, output_lines)
        push!(RESULTS, 
            (title, program, config, 
             input_scripts, output_lines, 
             results))
        if !(success)
            break
        end
        REGISTER = updater(output_lines)
        if (:exit ∈ keys(REGISTER)) && REGISTER[:exit]
            break
        end
    end
    return RESULTS
end

