cd("/data/phonopy_database/")

using JLD2

@load "mat.jld2"

mp_list = map(x->"mp-$(x)-20180417", last.(mat))

##

f(y) = ("mp-$(replace(y[1],"Materials id "=>""))-20180417" => (replace(y[2],r"\s+"=>""), y[3]))

mp_names = Dict( f(strip.(split(k,"/",keepempty=false))) for k in first.(mat) )

##

@load "all_freqs.jld2"

all_freqs


##? --------------------------------------------------------------

function first_gap_in_freq(freqs::Vector{Float64},df=0.01)
    if length(freqs)==0
        return -1.0
    else
        f = freqs[freqs.>df]
        return length(f)>0 ? minimum(f) : -1.0
    end
end


function imag_f_in_freq(freqs::Vector{Float64},df=0.01)
    if length(freqs)==0
        return false
    else
        return minimum(freqs) < -df
    end
end


##? --------------------------------------------------------------

freq_table = Dict(k=>(mp_names[k][1],mp_names[k][2],imag_f_in_freq(v),first_gap_in_freq(v)) for (k,v) in all_freqs)


for (k,v) ∈ freq_table
    
end

##? --------------------------------------------------------------

using Distributed

#addprocs(100)

##? --------------------------------------------------------------

@everywhere using YAML


##? --------------------------------------------------------------

@everywhere function get_lzma(mp, workspace="/data/phonopy_database", phonopy_root="/home/dabajabaza/anaconda3/bin")
    if !isfile("$workspace/$(mp).tar.lzma")
        download("http://phonondb.mtl.kyoto-u.ac.jp/_downloads/$(mp).tar.lzma", "$workspace/$(mp).tar.lzma")
    end
    fd = "$workspace/$(mp)"
    if !isdir(fd)
        try 
            run(`tar --lzma -xvf $workspace/$(mp).tar.lzma --directory $workspace/`)
        catch _e_
            @error "Error in run tar :"
            @error _e_
            return Dict()
        end
    end
    cd(fd)
    fn = readdir(fd)
    if "disp.yaml" ∈ fn
        if "BORN" ∈ fn
            try 
                run(`$phonopy_root/phonopy -c disp.yaml --include-all --nac --qpoints 0 0 0`)
            catch _e_
                @error "Error in run phonopy :"
                @error _e_
                return Dict()
            end
        else
            try 
                run(`$phonopy_root/phonopy -c disp.yaml --include-all --qpoints 0 0 0`)
            catch _e_
                @error "Error in run phonopy :"
                @error _e_
                return Dict()
            end
        end
    else
        @error "disp.yaml not found in $fd ."
        return Dict()
    end

    qpoints = YAML.load_file("$fd/qpoints.yaml"; dicttype=Dict{String,Any})
    
    return qpoints

end


@everywhere function qpoint_parser(qp::Dict)::Vector{Float64}
    if "phonon" ∈ keys(qp)
        f0 = qp["phonon"]
        if length(f0)>1 && "band"∈keys(f0[1])
            f = f0[1]["band"]
            if (f isa Vector) && length(f)>0 && "frequency"∈keys(f[1])
                return map(x->x["frequency"], f)
            end
        end
    end
    return Float64[]
end


@everywhere get_gap(mp) = first_gap_in_freq(qpoint_parser(get_lzma(mp)))

##

results = pmap(get_gap, mp_list)

##







#=

## parse ph20180417_index.html

using Gumbo

page_source = join(readlines("/home/dabajabaza/Downloads/ph20180417_index.html"),"\n") ;

function parse_ph20180417_index(s)
    html = parsehtml(s)
    r = html.root
    useful_lines = r.children[2].children[2].children[1].children
    materials = []
    for l in useful_lines
        if length(l.children[2].children) >= 5
            try
                href = l.children[2].children[2].children[8].attributes["href"]
                txt = l.children[2].children[3].text
                if !occursin("index.html",href)
                    push!(materials, (txt, href, replace(basename(href),r"(\.html|mp\-)"=>"")))
                end
            catch 
                nothing
            end
        end
    end
    return materials
end

=#
