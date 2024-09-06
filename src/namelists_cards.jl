using Printf

#* ============================ Namelist  ==============================

mutable struct Namelist

    #* title of the namelist, for example, "SYSTEM", "CONTROL", ...
    title::String

    #* MUST BE SPECIFIED IN EACH CALCULATION
    req::Dict{Union{String,Symbol},Any}

    #* these "opt" settings distinguishes the types of calculations
    opt::Dict{Union{String,Symbol},Any}

    #* for default settings we prefer to leave it unspecified 
    #> however, some settings have big impact on the quality and efficiency 
    #> of the calculation
    #> the following "default" field is a summary of them
    #+ it depends on the experience on QE and may be expanded later
    default::Dict{Union{String,Symbol},Any}
end


#! required is required, optional is optional !!!

is_valid(nl::Namelist) = (length(intersect(keys(nl.req),keys(nl.opt)))==0)

empty_Namelist(nl::Namelist) = (length(nl.req)==0 && length(nl.opt)==0)


#* ============================ Card{T}  ==============================

mutable struct Card{T}

    #* title of the card, for example, "ATOM_POSITIONS", ...
    title::String

    #* options in the curly bracket
    option::String

    #* data 
    contents::Vector{Vector{T}}

end

empty_Card(c::Card) = length(c.contents)==0


##* ============================ BUILD  ==============================


function build_namelist(nl::Namelist)

    function adapt(x)
        if x isa Bool
            return (x ? ".true." : ".false.")
        elseif x isa AbstractString
            return "'$(x)'"
        elseif x isa Integer
            return string(x)
        elseif x isa Float64
            return strip((@sprintf "%20.15f" x))
        else
            return string(x)
        end
    end

    if empty_Namelist(nl) return String[] end
    @assert is_valid(nl)
    lines = ["&" * strip(nl.title),]
    dict = nl.req  ⬱  nl.opt    #+ FIXME
    for r ∈ sort(collect(keys(dict))) 
        #* ignore symbol keys 
        #* symbol keys are for control purposes
        if r isa String  
            v = dict[r]
            if v!==nothing
                push!(lines, (@sprintf "    %18s = %s" string(r) adapt(v)))
            end
        end
    end
    push!(lines, "/")
    return lines
end



function build_card(c::Card{T}) where T
    if empty_Card(c) return String[] end
    first_line = strip(c.title * "  " * c.option)
    return [ length(first_line)==0 ? [] : [first_line,]; 
             parse_join.(c.contents) ]
end

#* ======================  FOR DIFFERENT PROGRAM ========================


mutable struct INPUT_PWX

    TITLE::String
    
    CONTROL::Namelist
    SYSTEM::Namelist
    ELECTRONS::Namelist
    IONS::Namelist
    CELL::Namelist

    ATOMIC_SPECIES::Card{Union{String,Float64,Int}}
    ATOMIC_POSITIONS::Card{Union{String,Float64,Int}}
    K_POINTS::Card{Union{Float64,Int}}
    CELL_PARAMETERS::Card{Float64}
    CONSTRAINTS::Card{Union{String,Float64,Int}}
    OCCUPATIONS::Card{Union{Float64,Int}}
    ATOMIC_VELOCITIES::Card{Union{Float64,Int}}
    ATOMIC_FORCES::Card{Union{Float64,Int}}

end


mutable struct INPUT_PHX

    TITLE::String

    INPUTPH::Namelist

    Q::Card{Union{Float64,Int}}

end


mutable struct INPUT_PPX

    TITLE::String

    INPUTPP::Namelist

    PLOT::Namelist

end


mutable struct INPUT_LD1X

    TITLE::String

    INPUT::Namelist
    AllElectron::Card{Union{String,Float64,Int}}

    INPUTP::Namelist
    PseudoPotentialGeneration::Card{Union{String,Float64,Int}}

    TEST::Namelist
    PseudoPotentialTest::Card{Union{String,Float64,Int}}

end


mutable struct INPUT_CPX

    TITLE::String

    CONTROL::Namelist
    SYSTEM::Namelist
    ELECTRONS::Namelist
    IONS::Namelist
    CELL::Namelist
    WANNIER::Namelist

    ATOMIC_SPECIES::Card{Union{String,Float64,Int}}
    ATOMIC_POSITIONS::Card{Union{String,Float64,Int}}
    CELL_PARAMETERS::Card{Float64}
    OCCUPATIONS::Card{Union{Float64,Int}}
    CONSTRAINTS::Card{Union{String,Float64,Int}}
    ATOMIC_VELOCITIES::Card{Union{Float64,Int}}
    ATOMIC_FORCES::Card{Union{Float64,Int}}
    AUTOPILOT::Card{String}

end


mutable struct INPUT_MATDYN

    TITLE::String

    INPUT::Namelist

    Q::Card{Union{Float64,Int}}

end

mutable struct INPUT_PW2WANNIER90

    TITLE::String

    INPUT::Namelist

end


mutable struct INPUT_Q2R

    TITLE::String

    INPUT::Namelist

end

mutable struct INPUT_PROJWFC

    TITLE::String

    PROJWFC::Namelist

end


mutable struct INPUT_BANDS

    TITLE::String

    BANDS::Namelist

end


function build_input_file(IN)::Vector{String}
    lines = []
    for f ∈ fieldnames(typeof(IN))
        x = getfield(IN,f)
        if x isa AbstractString
            lines = [lines; [x,]]
        elseif x isa Namelist
            lines = [lines; build_namelist(x)]
        elseif x isa Card
            lines = [lines; build_card(x)]
        else
            throw(error("unknown type of field $f."))
        end
    end
    return [lines; ["\n"]]
end
