# ld1x_input.jl

function INPUTLD1(config)
    NL = Namelist("INPUT", Dict(), Dict(), Dict())

    NL.req = Dict(
        "iswitch"  => config["iswitch"],
        "title"    => config["title"],
        "prefix"   => config["prefix"],
        "zed"      => config["zed"],
        "rel"      => config["rel"],
        "lsd"      => config["lsd"],
        "config"   => config["config"],
        "dft"      => config["dft"],
    )
    NL.default = Dict(
        "verbosity"=> "high",
        "max_out_wfc" => 10,
        "noscf"    => false,
        "tr2"      => 1e-14,
        #* Radial grid parameters:
        "xmin"     => ((config["iswitch"]>1 || config["rel"]==0) ? -7.0 : -8.0),
        "dx"       => ((config["iswitch"]>1) ? 0.0125 : 0.008),
        "rmax"     => 100.0,
        #* Parameters for logarithmic derivatives:
        "nld"      => nothing,
        "rlderiv"  => nothing,
        "eminld"   => nothing,
        "emaxld"   => nothing,
        "deld"     => nothing,
        "rpwe"     => nothing,
    )
    NL.opt = NL.default ← config
    return NL
end


function AllElectron(config)
    CD = Card{Union{String,Float64,Int}}("", "", Vector{Union{String,Int,Float64}}[])

    rel = config["rel"]
    econf = config[:electronic_configuration]
    if rel<2
        #:  nl(1) n(1) l(1) oc(1) isw(1) 
        typ1(x) = collect(Tuple{String,Int,Int,Float64,Int}((x...,))) |> Vector{Union{String,Int,Float64}}
        CD.contents = Vector{Union{String,Int,Float64}}[[length(econf)], typ1.(econf)...]
    elseif rel==2
        #:  nl(1) n(1) l(1) oc(1) jj(1) 
        typ2(x) = collect(Tuple{String,Int,Int,Float64,Float64}((x...,))) |> Vector{Union{String,Int,Float64}}
        CD.contents = Vector{Union{String,Int,Float64}}[[length(econf)], typ2.(econf)...]
    end
    return CD
end



function INPUTP(config)
    NL = Namelist("INPUTP", Dict(), Dict(), Dict())

    if "PP generation" ∈ config[:ld1_mode]
        NL.req = Dict(
            "pseudotype"     => config["pseudotype"],
            "tm"             => true,
            "nlcc"           => config["nlcc"],
            "file_pseudopw"  => config["file_pseudopw"],  #: example 'Sc.pbesol-spn-rrkjus_psl.1.0.0.UPF',
        )
        NL.default = Dict(
            "author"         => "ld1.x",
            "lloc"           => -1,
            "lpaw"           => false,
            "file_recon"     => "",
            "which_augfun"   => "PSQ",
            "rmatch_augfun_nc" => false,
            "rmatch_augfun"  => 0.5,
            "rcore"          => nothing,
            "rcloc"          => nothing,
            "use_xsd"        => false,
            "new_core_ps"    => true,
            "lsave_wfc"      => true,
            "file_chi"       => nothing,   
            "file_beta"      => nothing,   
            "file_qvan"      => nothing,   
            "file_screen"    => nothing,   
            "file_core"      => nothing,   
            "file_wfcaegen"  => nothing,   
            "file_wfcncgen"  => nothing,   
            "file_wfcusgen"  => nothing,   
        )
        NL.opt = NL.default ← config
    end
    return NL
end


function PseudoPotentialGeneration(config)
    CD = Card{Union{String,Float64,Int}}("", "", Vector{Union{String,Int,Float64}}[])

    rel = config["rel"]
    if "PP generation" ∈ config[:ld1_mode]
        wfc_spec = config[:wavefunctions_spec]
        if rel==0 || rel==2
            #: nls(1) nns(1) lls(1) ocs(1) ener(1) rcut(1) rcutus(1) jjs(1) 
            typ3(x) = collect(Tuple{String,Int,Int,Float64,Float64,Float64,Float64,Float64}((x[1:8]...,))) |> Vector{Union{String,Int,Float64}}
            CD.contents = Vector{Union{String,Int,Float64}}[[length(wfc_spec)], typ3.(wfc_spec)...]
        else
            #: nls(1) nns(1) lls(1) ocs(1) ener(1) rcut(1) rcutus(1) 
            typ4(x) = collect(Tuple{String,Int,Int,Float64,Float64,Float64,Float64}((x[1:7]...,))) |> Vector{Union{String,Int,Float64}}
            CD.contents = Vector{Union{String,Int,Float64}}[[length(wfc_spec)], typ4.(wfc_spec)...]
        end
    end
    return CD
end


function check_conf_ld1x(config)
    @assert "iswitch" ∈ keys(config)
    @assert "rel"     ∈ keys(config)

    if config["iswitch"]==4
        @assert "nconf" ∈ keys(config)
        @assert config["nconf"]==2
    elseif config["iswitch"]==2
        @assert "nconf" ∈ keys(config)
    elseif config["iswitch"]==3
        @assert "lsd" ∉ keys(config)
    end

    if config["rel"]==2
        @assert "lsd" ∉ keys(config)
    end

    if "test" ∈ config[:ld1_mode]
        if "configts(1)" ∈ keys(config)
            @assert :test_spec ∉ keys(config)
        end
        if :test_spec ∈ keys(config)
            @assert "configts(1)" ∉ keys(config)
        end
    end

    return true
end


function TEST(config)
    NL = Namelist("TEST", Dict(), Dict(), Dict())

    if "test" ∈ config[:ld1_mode]
        check_conf_ld1x(config)
        configts = Dict("configts($i)"=>config["configts($i)"] for i=1:config["nconf"] if "configts($i)" ∈ keys(config))
        lsdts = Dict("lsdts($i)"=>config["lsdts($i)"] for i=1:config["nconf"] if "lsdts($i)" ∈ keys(config))
        NL.req = Dict(
            "nconf"       => config["nconf"],
            "file_pseudo" => config["file_pseudo"]
        ) ∪ configts
        NL.default = Dict(
            "frozen_core" => false,
            "rcutv"       => -1.0,
            "rm"          => 30.0,
            "ecutmin"     => 0.0,
            "ecutmax"     => 0.0,
            "decut"       => 5.0,
        )
        NL.opt = NL.default ← config
    end

    return NL
end


function PseudoPotentialTest(config)
    CD = Card{Union{String,Int,Float64}}("", "", Vector{Union{String,Int,Float64}}[])

    if "test" ∈ config[:ld1_mode]
        rel = config["rel"]
        lsd = config["lsd"]
        test_spec = config[:test_spec]
        if lsd==1
            #: elts(1) nnts(1) llts(1) octs(1) enerts(1) rcutts(1) rcutusts(1) iswts(1) 
            typ5(x) = collect(Tuple{String,Int,Int,Float64,Float64,Float64,Float64,Int}((x...,))) |> Vector{Union{String,Int,Float64}}
            CD.contents = Vector{Union{String,Int,Float64}}[[length(test_spec)], typ5.(test_spec)...]
        elseif rel==2
            #: elts(1) nnts(1) llts(1) octs(1) enerts(1) rcutts(1) rcutusts(1) jjts(1) 
            typ6(x) = collect(Tuple{String,Int,Int,Float64,Float64,Float64,Float64,Float64}((x...,))) |> Vector{Union{String,Int,Float64}}
            CD.contents = Vector{Union{String,Int,Float64}}[[length(test_spec)], typ6.(test_spec)...]
        else
            #: elts(1) nnts(1) llts(1) octs(1) enerts(1) rcutts(1) rcutusts(1) 
            typ7(x) = collect(Tuple{String,Int,Int,Float64,Float64,Float64,Float64}((x...,))) |> Vector{Union{String,Int,Float64}}
            CD.contents = Vector{Union{String,Int,Float64}}[[length(test_spec)], typ7.(test_spec)...]
        end
    end
    return CD
end


##* ===============================================================


function electron_configuration_report(config::Dict)
    rep       = []
    push!(rep, "----------- LD1.X CONFIGURATION -----------")

    econf     = config[:electronic_configuration]
    wfc_spec  = ("PP generation" ∈ config[:ld1_mode]) ? config[:wavefunctions_spec] : []
    test_spec = ("test" ∈ config[:ld1_mode]) ? config[:test_spec] : []
    rel       = config["rel"]

    #: ----------------------------------
    push!(rep, "All-electron calculation configuration : ")
    if rel<2
        #:  nl(1) n(1) l(1) oc(1) isw(1) 
        parse_join1(line,i) = replace(join(map(x->x[1]*strip(x[2]),zip(
                                        ["nl($i)=","n($i)=","l($i)=","oc($i)=","isw($i)="],
                                        parse_elem.(line))),
                                        "  "), r"(\s\s\s\s\s|0000000000)"=>"")
        for (i,c) ∈ enumerate(econf)
            push!(rep, parse_join1(c,i))
        end
        if config["config"]!=""
            push!(rep, config["config"])
        end
    elseif rel==2
        #:  nl(1) n(1) l(1) oc(1) jj(1) 
        parse_join2(line,i) = replace(join(map(x->x[1]*strip(x[2]),zip(
                                        ["nl($i)=","n($i)=","l($i)=","oc($i)=","jj($i)="],
                                        parse_elem.(line))),
                                        "  "), r"(\s\s\s\s\s|0000000000)"=>"")
        for (i,c) ∈ enumerate(econf)
            push!(rep, parse_join2(c,i))
        end
        if config["config"]!=""
            push!(rep, strip(config["config"]))
        end
    end

    #: ----------------------------------
    if "PP generation" ∈ config[:ld1_mode]
        push!(rep, "PP generation configuration : ")
        if rel==0 || rel==2
            #: nls(1) nns(1) lls(1) ocs(1) ener(1) rcut(1) rcutus(1) jjs(1) 
            parse_join3(line,i) = replace(join(map(x->x[1]*strip(x[2]),zip(
                                            ["nls($i)=","nns($i)=","lls($i)=","ocs($i)=","ener($i)=","rcut($i)=","rcutus($i)=","jjs($i)="],
                                            parse_elem.(line[1:8]))),
                                            "  "), r"(\s\s\s\s\s|0000000000)"=>"")
            for (i,c) ∈ enumerate(wfc_spec)
                push!(rep, parse_join3(c,i))
            end
        else
            #: nls(1) nns(1) lls(1) ocs(1) ener(1) rcut(1) rcutus(1) 
            parse_join4(line,i) = replace(join(map(x->x[1]*strip(x[2]),zip(
                                            ["nls($i)=","nns($i)=","lls($i)=","ocs($i)=","ener($i)=","rcut($i)=","rcutus($i)="],
                                            parse_elem.(line[1:7]))),
                                            "  "), r"(\s\s\s\s\s|0000000000)"=>"")
            for (i,c) ∈ enumerate(wfc_spec)
                push!(rep, parse_join4(c,i))
            end
        end
    end

    #: ----------------------------------
    if "test" ∈ config[:ld1_mode]
        lsd = config["lsd"]
        push!(rep, "Test configuration : ")
        if lsd==1
            #: elts(1) nnts(1) llts(1) octs(1) enerts(1) rcutts(1) rcutusts(1) iswts(1) 
            parse_join5(line,i) = replace(join(map(x->x[1]*strip(x[2]),zip(
                                            ["elts($i)=","nnts($i)=","llts($i)=","octs($i)=","enerts($i)=","rcutts($i)=","rcutusts($i)=","iswts($i)="],
                                            parse_elem.(line))),
                                            "  "), r"(\s\s\s\s\s|0000000000)"=>"")
            for (i,c) ∈ enumerate(test_spec)
                push!(rep, parse_join5(c,i))
            end
        elseif rel==2
            #: elts(1) nnts(1) llts(1) octs(1) enerts(1) rcutts(1) rcutusts(1) jjts(1) 
            parse_join6(line,i) = replace(join(map(x->x[1]*strip(x[2]),zip(
                                            ["elts($i)=","nnts($i)=","llts($i)=","octs($i)=","enerts($i)=","rcutts($i)=","rcutusts($i)=","jjts($i)="],
                                            parse_elem.(line))),
                                            "  "), r"(\s\s\s\s\s|0000000000)"=>"")
            for (i,c) ∈ enumerate(test_spec)
                push!(rep, parse_join6(c,i))
            end
        else
            #: elts(1) nnts(1) llts(1) octs(1) enerts(1) rcutts(1) rcutusts(1) 
            parse_join7(line,i) = replace(join(map(x->x[1]*strip(x[2]),zip(
                                            ["elts($i)=","nnts($i)=","llts($i)=","octs($i)=","enerts($i)=","rcutts($i)=","rcutusts($i)=",],
                                            parse_elem.(line))),
                                            "  "), r"(\s\s\s\s\s|0000000000)"=>"")
            for (i,c) ∈ enumerate(test_spec)
                push!(rep, parse_join7(c,i))
            end
        end
    end

    #: ----------------------------------
    push!(rep, "---------------------------------------------------")
    return join(rep, "\n")
end


function ld1x_inupt(config::Dict)
    @info electron_configuration_report(config)
    return INPUT_LD1X(
                config["title"], 
                map( x->x(config), 
                     [INPUTLD1, AllElectron, 
                      INPUTP, PseudoPotentialGeneration, 
                      TEST, PseudoPotentialTest]
                )...) |> build_input_file
end
