
function Find_XXX(XXX, pseudo_fn)
    lines = readlines(pseudo_fn)
    p  = findfirst(x->occursin(XXX * r"\s*\:",x), lines)
    if p!==nothing
        return last(SPLTC(lines[p]))
    end
    p2  = findfirst(x->occursin(r"\s+"*XXX,x), lines)
    if p2!==nothing
        return first(SPLTS(lines[p2]))
    end
    p1  = findfirst(x->occursin(lowercase(XXX) * r"\s*\=\s*",x), lines)
    if p1!==nothing
        return strip(last(SPLTEQ(lines[p1])),'"')
    else
        throw(error("Find_$(XXX)(): Pseudo file $pseudo_fn has unsupported format."))
    end
end


Find_element_name(pseudo_fn) = strip(Find_XXX("Element", pseudo_fn))


Find_functional(pseudo_fn) = Find_XXX("Functional", pseudo_fn)

