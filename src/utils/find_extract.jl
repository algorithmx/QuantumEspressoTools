function find_line(res, cont_str; quiet=false)
    p = findfirst(x->occursin(cont_str,x), res)
    if p===nothing && !quiet
        @error "find_line():\nString ($cont_str) not found."
    end
    return p
end


function find_last_line(res, cont_str; quiet=false)
    p = findlast(x->occursin(cont_str,x), res)
    if p===nothing
        if !quiet  @error "find_last_line():\nString ($cont_str) not found." end
    end
    return p
end


function find_all_lines(res, cont_str)
    p = findall(x->occursin(cont_str,x), res)
    if length(p)==0
        @error "find_all_lines():\nString ($cont_str) not found."
    end
    return p
end


function find_empty_line(res; quiet=true)
    find_line(res, r"^\s*$"; quiet=quiet)
end

# return num of steps to go for cond(line) to be false
function delta_find_until_not(lines, cond)
    delta = 0
    for l in lines
        if cond(l)
            delta += 1
        else
            break
        end
    end
    return delta
end


function extract(res, reg1, reg2)
    pos  = findfirst(x->occursin(reg1,x), res)
    m1   = (pos !== nothing ? match(reg1,res[pos]) : nothing)
    m2   = (m1  !== nothing ? m1.match : "")
    m3   = match(reg2,m2)
    m4   = (m3  !== nothing ? m3.match : "")
    m4
end


function extract_all(
    res, 
    reg1::Union{Regex,String}, 
    reg2::Union{Regex,String}
    )
    pos  = findall(x->occursin(reg1,x), res)
    m1   = [match(reg1,res[p]) for p in pos]
    m2   = String[m.match for m in m1 if m!==nothing]
    m3   = [match(reg2,m) for m in m2]
    m4   = String[m.match for m in m3 if m!==nothing]
    m4
end


function extract_all(res, reg1::Union{Regex,String}, reg2L::Vector{Regex})
    pos  = findall(x->occursin(reg1,x), res)
    m1   = [match(reg1,res[p]) for p in pos]
    m2   = String[m.match for m in m1 if m!==nothing]
    m3   = [[match(reg2,m) for reg2 in reg2L] for m in m2]
    Vector{String}[String[x.match for x in m] for m in m3]
end


function extract_all2(res, reg1::Union{Regex,String}, reg2L::Vector{Regex})
    pos  = findall(x->occursin(reg1,x), res)
    m1   = [match(reg1,res[p]) for p in pos]
    m2   = String[m.match for m in m1 if m!==nothing]
    m3   = [[String[m[i] for i in findall(reg2,m)] for reg2 in reg2L] for m in m2]
    Vector{Vector{String}}[Vector{String}[m...] for m in m3]
end


function extract_last(res, reg1, reg2)
    pos  = findlast(x->occursin(reg1,x), res)
    m1   = (pos !== nothing ? match(reg1,res[pos]) : nothing)
    m2   = (m1  !== nothing ? m1.match : "")
    m3   = match(reg2,m2)
    m4   = (m3  !== nothing ? m3.match : "")
    m4
end


⊂(subs::AbstractString, s::AbstractString) = occursin(subs, s)

⊂(subs::Regex, s::AbstractString) = occursin(subs, s)

⊂(s::AbstractString, L::Vector{S}) where {S<:AbstractString} = (findfirst(x->(s ⊂ x),L)!==nothing)

⊂(s::Regex, L::Vector{S}) where {S<:AbstractString} = (findfirst(x->(s ⊂ x),L)!==nothing)
