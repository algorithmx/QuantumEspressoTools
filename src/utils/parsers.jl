
global const num_f_rstr = r"-?\d+(\.\d*)?((e|E)-?\d+)?"

global const num_f_rstr3 = num_f_rstr*r"\s+"*num_f_rstr*r"\s+"*num_f_rstr

global const num_f_rstr6 = num_f_rstr3 * r"\s+" * num_f_rstr3

##* ================================================================

@inline check_nothing(a...) = all(a.!==nothing)


delete_empty_lines(lines::Vector{String}) = [s for s in lines if length(strip(s))>0]




@inline strtonum(s) = (s isa Number) ? s : (try parse(Int, s) catch; parse(Float64, s) end)


@inline function strtonum_fort(s)
    r = 0
    if (s isa Number)
        return s
    elseif (s isa Vector)
        return strtonum_fort.(s)
    else
        try
            r = parse(Int, s)
        catch
            try r = parse(Float64, s) catch; return s end
        end
    end
    return r
end


parsef(x) = (try parse(Float64,strip(x)) catch; NaN end)

parse2f(x) = (try parse.(Float64, SPLTS(x)) catch; [NaN,NaN] end)

parse3f(x) = (try parse.(Float64, SPLTS(x)) catch; [NaN,NaN,NaN] end)

parsenf(x) = (try parse.(Float64, SPLTS(x)) catch; repeat([NaN,], length(SPLTS(x))) end)


function parse_1_string_3_float64(x)
    x1 = split(x, " ", keepempty=false, limit=2)
    return (string(x1[1]) => parse3f(x1[2]))
end

parse1s3f(x) = parse_1_string_3_float64(x)

evalt(t) = eval(Meta.parse(t))


##* ================================================================


function parse_elem(a)
    if a isa Integer
        return (@sprintf "%8d" a)
    elseif a isa Float64
        return (@sprintf "%20.15f" a)
    elseif a isa String
        return (@sprintf "%-8s" a)
    else
        return string(a)
    end
end


parse_join(line) = join(parse_elem.(line), "  ")


showdict(dic::Dict) = (
    "\n" 
    * join([("$k\t = $((v isa Dict) 
                      ? showdict(v) 
                      : ((v isa Vector) 
                            ? "\n"*join(string.(v),"\n")*"\n" 
                            : v))") 
             for (k,v) ∈ dic], "\n") 
    * "\n"
)
