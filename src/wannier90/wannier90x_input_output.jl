# wannier90_input.jl

function kmesh_pl(n1,n2,n3,n4=1)
    omit_weight = (n4 > 0)
    weight      = 1.0/prod([n1,n2,n3])
    return [ omit_weight ? Float64[i//n1, j//n2, k//n3] : Float64[i//n1, j//n2, k//n3, weight]
                for i=0:n1-1 for j=0:n2-1 for k=0:n3-1 ]
end


function adapt_w90(x)
    if x isa Bool
        return (x ? ".true." : ".false.")
    elseif x isa AbstractString
        return x
    elseif x isa Integer
        return string(x)
    elseif x isa Float64
        return strip((@sprintf "%20.15f" x))
    elseif (x isa Tuple)
        return join(adapt_w90.(x), "  ")
    elseif (x isa Vector)
        return join(adapt_w90.(x), "  ")
    else
        return string(x)
    end
end


wannier90x_kpath_from_QE(kpath_QE) = [
    "A$i  $(kpath_QE[i][1][1])  $(kpath_QE[i][1][2])  $(kpath_QE[i][1][3])  A$(i+1)  $(kpath_QE[i+1][1][1])  $(kpath_QE[i+1][1][2])  $(kpath_QE[i+1][1][3])"
    for i=1:length(kpath_QE)-1
]


function wannier90x_block2lines(block_name, block_data)
    lines = []
    push!(lines, "begin $block_name")
    for data ∈ block_data
        push!(lines, adapt_w90(data))
    end
    push!(lines, "end $block_name")
    return lines
end


const wannier90_config_default = Dict(
    "num_iter"        => 200,
    "dis_num_iter"    => 16000,
    "dis_mix_ratio"   => 0.9,
    "wannier_plot"    => false,
    "bands_plot"      => false,
    "postproc_setup"  => false,
    "write_bvec"      => true,
    "guiding_centres" => false,
    "write_hr"        => true,
    "spinors"         => false,
    "postproc_setup"  => false,
)


function wannier90x_inupt(config::Dict)
    config_lines = []
    for (k,v) ∈ config
        if k isa AbstractString
            if v isa Vector
                config_lines = [ config_lines;
                                 wannier90x_block2lines(k,v) ]
            else
                push!(config_lines, "$k = $(adapt_w90(v))")
            end
        end
    end
    push!(config_lines, "\n")
    return config_lines
end

global const W90_spred_i_reg   = r"WF\s+centre\s+and\s+spread\s+\d+\s*\(\s*-?\d+\.\d+,\s*-?\d+\.\d+,\s*-?\d+\.\d+\s*\)\s*\d+\.\d+"

global const W90_spred_sum_reg = r"Sum\s*of\s*centres\s*and\s*spreads\s*\(\s*-?\d+\.\d+,\s*-?\d+\.\d+,\s*-?\d+\.\d+\s*\)\s*\d+\.\d+"

function wannier90_get_spread(
    W90_lines::Vector{S}
    ) where {S<:AbstractString}
    @inline ι(x) = parse(Int,x)
    @inline φ(x) = parse(Float64,x)
    @inline _get_id(x) = ι(split(strip(x)," ",keepempty=false)[2])
    @inline _get_pos(x) = φ.(strip.((split(replace(x,r"[\(\)]"=>""),",",keepempty=false))))
    @inline _get_spread(x) = φ(strip((split(x," ",keepempty=false)[2])))
    @inline extract_spread(x) = (_get_id(x[1]))=>(pos=_get_pos(x[2]),spread=_get_spread(x[3]))
    pos_Final_State = find_line(W90_lines, "Final State"; quiet=true)
    if pos_Final_State !== nothing
        return QuantumEspressoTools.extract_all(
                W90_lines[pos_Final_State:end], 
                W90_spred_i_reg, 
                [r"spread\s+\d+\s*\(", 
                r"\(\s*-?\d+\.\d+,\s*-?\d+\.\d+,\s*-?\d+\.\d+\s*\)",
                r"\)\s*\d+\.\d+"]) .|> extract_spread
    else
        @error "wannier90_get_spread() did not find \"Final State\" section."
        return []
    end
    # Wspred_sum = QuantumEspressoTools.extract_all(
    #     W90_lines[pos_Final_State:end], 
    #     W90_spred_sum_reg, 
    #     [r"\)\s*\d+\.\d+",])
end

global const W90_outer_reg = r"Outer:\s*-?\d+.\d+\s*to\s*-?\d+.\d+\s*\(eV\)"
global const W90_inner_reg = r"Inner:\s*-?\d+.\d+\s*to\s*-?\d+.\d+\s*\(eV\)"

function wannier90_get_window(
    W90_lines::Vector{S}
    ) where {S<:AbstractString}
    @inline ι(x) = parse(Int,x)
    @inline φ(x) = parse(Float64,x)
    @inline extract_window(x) = φ.(x[1][1])
    pos_Energy_Windows = find_line(W90_lines, r"Energy\s+Windows"; quiet=true)
    if pos_Energy_Windows !== nothing
        try
            outer = QuantumEspressoTools.extract_all2(
                    W90_lines[pos_Energy_Windows:pos_Energy_Windows+4], 
                    W90_outer_reg, [r"-?\d+.\d+",]) |> extract_window
            inner = QuantumEspressoTools.extract_all2(
                    W90_lines[pos_Energy_Windows:pos_Energy_Windows+4], 
                    W90_inner_reg, [r"-?\d+.\d+",]) |> extract_window
            return Dict("Inner"=>inner, "Outer"=>outer)
        catch _e_
            @error "wannier90_get_window() find \"Energy  Windows\" error : \n$(_e_)\n\n"
            return Dict("Inner"=>[], "Outer"=>[])
        end
    else
        @error "wannier90_get_window() did not find \"Energy  Windows\" section."
        return Dict("Inner"=>[], "Outer"=>[])
    end
end
