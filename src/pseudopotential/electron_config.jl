#NOTE
#> fractional occupations : example GBRV/all_pbesol_UPF_v1.5/re_pbesol_v1.2.uspp.F.UPF
#> mixed pseudo files : example SSSP_Efficiency uses SG15 and others


mutable struct ELECTRON_CONFIG
    s1::Int
    p1::Int

    s2::Int
    p2::Int

    s3::Int
    p3::Int
    d3::Int

    s4::Int
    p4::Int
    d4::Int
    f4::Int

    s5::Int
    p5::Int
    d5::Int
    f5::Int

    s6::Int
    p6::Int
    d6::Int
    f6::Int

    s7::Int
    p7::Int
    d7::Int
    f7::Int
end


ELECTRON_CONFIG() = ELECTRON_CONFIG(0,0, 0,0, 0,0,0, 0,0,0,0, 0,0,0,0, 0,0,0,0, 0,0,0,0)


function ELECTRON_CONFIG(orbits::NTuple)
    elc = ELECTRON_CONFIG()
    for orb in orbits
        for nl in orb
            o = Symbol(lowercase(nl[2:-1:1]))
            i = parse(Int,nl[3:end])
            setfield!(elc, o, i)
        end
    end
    return elc
end


global const shell_max = ELECTRON_CONFIG(2,6, 2,6, 2,6,10, 2,6,10,14, 2,6,10,14, 2,6,10,14, 2,6,10,14)


function ELECTRON_CONFIG(a, pseudo_lines::Vector{String}, psmode::String)
    elc = conventional_valence_electron_config[Atoms[a]]
    if any([startswith(psmode, x) for x in ["MT_","GHH_","GBRV_","PSL_","SSSP_", "SG15"]])
        p   = find_line(pseudo_lines, r"nl\s+pn\s+l\s+occ\s+Rcut\s+Rcut\s+US\s+E\s+pseu",  quiet=true)
        if p===nothing
            @warn "ELECTRON_CONFIG() found no valence configuration for element $a with $psmode ."
        else
            p = p+1
            L = strip(pseudo_lines[p])
            while (L!="") && !occursin(r"[Gg]eneration\s+configuration\s*\:",L) && !occursin("</PP_INFO>",L)
                (nl,pn,l,occ,Rcut,RcutUS,Epseu) = split(L," ",keepempty=false)
                orbit = Symbol(lowercase(nl[2:-1:1]))
                occ0  = getfield(elc, orbit)
                #> fractional occupations : example GBRV/all_pbesol_UPF_v1.5/re_pbesol_v1.2.uspp.F.UPF
                occf = parse(Float64,occ)
                setfield!(elc, orbit, occ0+Int(ceil(occf)))
                p = p+1
                L = pseudo_lines[p]
            end
        end
    end
    return elc
end


ELECTRON_CONFIG(a, pseudo_fn::String, psmode::String) = ELECTRON_CONFIG(a, readlines(pseudo_fn), psmode)


##
global const conventional_valence_electron_config =  ELECTRON_CONFIG.(parse_electron_config.(last.(ELEC_CONF)))


function num_electrons(c::ELECTRON_CONFIG)
    return sum(Int[getfield(c,f) for f in fieldnames(ELECTRON_CONFIG)])
end


function max_electrons(c::ELECTRON_CONFIG, spin=0)
    #TODO SPIN
    return sum(Int[getfield(shell_max,f) for f in fieldnames(ELECTRON_CONFIG) if getfield(c,f)>0])
end
