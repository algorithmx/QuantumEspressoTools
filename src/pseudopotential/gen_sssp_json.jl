using JSON

FD = "/home/dabajabaza/abinitio/pseudopotentials"

##

sssp_efficiency = JSON.parsefile("$FD/sssp_efficiency.json") ;

sssp_efficiency_pbesol_files = readdir("$FD/SSSP_efficiency_PBEsol")

findel(el) = findfirst(x->startswith(lowercase(x),Regex("$(lowercase(el))[_.]")),sssp_efficiency_pbesol_files)

function convertv(el,dic)
    ret = []
    for (k,v) ∈ dic
        if k=="filename"
            p = findel(el)
            if p !== nothing
                push!(ret, k=>sssp_efficiency_pbesol_files[p])
            else
                return Dict()
            end
        else
            push!(ret, k=>v)
        end
    end
    return Dict(ret)
end

sssp_efficiency_pbesol = Dict(el=>convertv(el,v) for (el,v) ∈ sssp_efficiency if length(convertv(el,v))>0)

##

open("$FD/sssp_efficiency_pbesol.json","w") do file
    js = JSON.json(sssp_efficiency_pbesol)
    write(file, js)
end


##

##

sssp_precision = JSON.parsefile("$FD/sssp_precision.json") ;

sssp_precision_pbesol_files = readdir("$FD/SSSP_precision_PBEsol")

findel(el) = findfirst(x->startswith(lowercase(x),Regex("$(lowercase(el))[_.]")),sssp_precision_pbesol_files)

function convertv(el,dic)
    ret = []
    for (k,v) ∈ dic
        if k=="filename"
            p = findel(el)
            if p !== nothing
                push!(ret, k=>sssp_precision_pbesol_files[p])
            else
                return Dict()
            end
        else
            push!(ret, k=>v)
        end
    end
    return Dict(ret)
end

sssp_precision_pbesol = Dict(el=>convertv(el,v) for (el,v) ∈ sssp_precision if length(convertv(el,v))>0)

##

open("$FD/sssp_precision_pbesol.json","w") do file
    js = JSON.json(sssp_precision_pbesol)
    write(file, js)
end
