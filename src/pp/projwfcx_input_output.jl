function PROJWFC(config::Dict)
    NL = Namelist("PROJWFC",Dict(),Dict(),Dict())
    NL.req  = Dict(
            "outdir"    => config["outdir"],
            "prefix"    => config["prefix"],
            "filproj"   => config["filproj"],
    )
    NL.default = Dict(
            "ngauss"    => 0,
            "DeltaE"    => 0.005,
            "degauss"   => 0.0,
            "Emin"      => -10.0,
            "Emax"      =>  10.0,
    )
    NL.opt = NL.default ← config
    return NL
end


function projwfcx_input(config::Dict)
    return INPUT_PROJWFC("", PROJWFC(config)) |> build_input_file
end


#% parse projwfc_up file (no soc)
function projwfcx_output_projwfc_up(result_lines)
    three_int_regex = r"^\s*\d{1,3}\s+\d{1,3}\s+\d{1,3}\s*$"
    f3p = find_last_line(result_lines, three_int_regex; quiet=true)
    @assert f3p !== nothing
    (norb, nkp, nbnd) = strtonum.(SPLTS(result_lines[f3p]))
    #@inline parse_block_header(l) = SPLTS(l)[[3,2,4]] => strtonum.(SPLTS(l)[[5,6,7]])
    @inline ι(x) = parse(Int, x)
    @inline organize(x) = ι(x[3]) => (
        atom="$(x[1][1])$(x[2][1])", 
        l=(x[2][3]), 
        m=(x[2][4]), 
        wfc=(x[2][2]), 
        orbit=(x[1][2]), 
    )
    @inline parse_block_header(l) = (
        SPLTS(l)[[3,4]], 
        strtonum.(SPLTS(l)[[2,5,6,7]]),
        SPLTS(l)[1]
    ) |> organize

    block_first_line_regex = r"^\s*\d+\s+\d+\s+\w+\s+[1-6][spdfSPDF]\s+\d+\s+\d+\s+\d+\s*$"
    pos = find_all_lines(result_lines, block_first_line_regex)
    @assert length(pos)==norb
    proj = zeros(Float64, nkp, nbnd, norb)
    header = []
    iik = 0
    iib = 0
    weight = 0.0
    #@time begin
        #TODO parser is really slow here but also file is big 
        #TODO (~30k lines, 10s processing time)
        for (iorb,p) ∈ enumerate(pos)
            @assert strtonum(SPLTS(result_lines[p])[1]) == iorb
            push!(header, parse_block_header(result_lines[p]))
            k = 1
            for ik=1:nkp
                for ib=1:nbnd
                    (iik,iib,weight) = strtonum.(SPLTS(result_lines[p+k]))
                    @assert iik==ik && iib==ib
                    proj[ik,ib,iorb] = weight
                    k += 1
                end
            end
            @assert k==nkp*nbnd+1
        end
    #end
    @inline cont_list(lst) = [((@assert i==p); v) for (i,(p,v)) in enumerate(lst)]
    states = header |> cont_list
    return proj, states
end
