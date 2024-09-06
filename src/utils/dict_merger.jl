#> difference between ∪, ⬱ and ← : 
#> ∪ is simple naive merge
#> ← pick from dic1 ∩ dic2, prefers values in dic2 and defaults on dic1
#> ⬱ take all from dic1


@inline ∪(dic1::Dict, dic2::Dict) = merge(Dict{Any,Any}(dic1),Dict{Any,Any}(dic2))

@inline ∪(dic1::Dict, p::Pair) = merge(Dict{Any,Any}(dic1),Dict{Any,Any}(p))


function ⬱(dic1::Dict, dic2::Dict)
    #: firstly, ⬱ is a merge
    dicm = copy(dic1 ∪ dic2)  #! FIXME
    #: then, contrtol the behavior from dic2
    for (k2,v2) in dic2
        if k2 ∈ keys(dic1)
            if v2 === nothing
                delete!(dicm, k2)
            else # v2 !== nothing
                if dic1[k2]!==nothing  &&  typeof(v2)!=typeof(dic1[k2])
                    @warn "The operator ⬱ found incompatible value types for k2=$(k2) : \ndic1[k2]=$(dic1[k2]) of type $(typeof(dic1[k2]))\nv2=$(v2) of type $(typeof(v2))."
                    dicm[k2] = v2  #! to be safe ...
                else
                    dicm[k2] = v2  #* update
                end
            end
        else 
            #: for k2 ∉ keys(dic1), ⬱ is just a merge
            nothing
        end
    end
    return dicm
end


@inline ⬱(dic1::Dict, p::Pair) = ⬱(dic1, Dict{Any,Any}(p))

⬱(dic1::Dict, P::Tuple) = length(P)==1 ? (dic1 ⬱ P[1]) : ((dic1 ⬱ P[1]) ⬱ P[2:end])

⬱(dic1::Dict, P::Vector) = length(P)==1 ? (dic1 ⬱ P[1]) : ((dic1 ⬱ P[1]) ⬱ P[2:end])


@inline use_default(v) = (v===nothing || v==:default || v==:use_default)

@inline ←(dic1::Dict, dic2::Dict) = Dict(   k=>(use_default(dic2[k]) ? dic1[k] : dic2[k]) 
                                            for k ∈ keys(dic1) if (k ∈ keys(dic2))  )

@inline ←(dic1::Dict, p::Pair) = ←(dic1, Dict{Any,Any}(p)) 


function purge(dic::Dict, k::Union{String,Symbol})

    if length(dic)==0  return copy(dic)  end
    if k===nothing return copy(dic)  end

    dic1 = copy(dic)
    if k ∉ keys(dic)
        return dic1
    else
        return delete!(dic1, k)
    end

end

function purge(dic::Dict, kk::Vector)
    d0 = copy(dic)
    for k ∈ kk
        d0 = purge(d0,k)
    end
    return d0
end

↓(dic::Dict, k) = purge(dic,k)


#: force update
↑(dic::Dict, k_v::Pair) = (dic ↓ k_v[1]) ∪ Dict(k_v)

