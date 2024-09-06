#PSLibrary.jl

#> don't try to remove "-" or "psl.". They help to identify the file names.
global const __ALL_SPECS = ["-n-", "-nl-", "-dn-", "-dnl-", "-s-", "-sl-", "-spn-", "-spdn-", "-spfn-", "-spnl-"]

global const __ALL_VERS  = ["psl.0.1.UPF",
                            "psl.0.2.UPF",   "psl.0.2.1.UPF", "psl.0.2.2.UPF", "psl.0.2.3.UPF", 
                            "psl.0.3.0.UPF", "psl.0.3.1.UPF",
                            "psl.1.0.0.UPF", "psl.1.0.1.UPF",]
global const ALL_PSL_SETTINGS = reshape(Iterators.product(__ALL_SPECS,__ALL_VERS)|>collect,:)


"""
    psl_all_possible_configs(atom_list, additional_filter)

returns a list of tuples
`(k ∈ keys(PSEUDO_DATA), Dict(elem=>(t1(c[1]),t2(c[2])) for (el,c) in zip(elemlist,conf)))`
for the possibel configurations of a PSL pseudo spec `k`. The dictionary has key => val paris of the following format:
`"F" => ("dnl", "1.0.0")`


The function `additional_filter(PSL_pseudo, spec)` selects pseudo settings of `PSL_pseudo`
with specification `spec`. For example: \n
`F_filter(PSL,spec) = (PSL=="PSL_PAW_PBESOL_SR" && spec[1][2]=="psl.1.0.0.UPF")`

"""
function psl_all_possible_configs(atom_list, additional_filter)
    elemlist = unique(atom_list)
    psl_configs = []
    # trim off the "-" or "psl." in version and spec notions
    @inline t1(x) = replace(x,"-"=>"")
    @inline t2(x) = replace(replace(x,".UPF"=>""),"psl."=>"")
    @inline all_in(k) = all([(el∈keys(PSEUDO_DATA[k])) for el in elemlist])
    for k in keys(PSEUDO_DATA)
        if startswith(k,"PSL") && all_in(k)
            all_conf = reshape(Iterators.product([keys(PSEUDO_DATA[k][el]) for el in elemlist]...)|>collect,:)
            for conf in all_conf
                if additional_filter(k, conf)
                    push!(  psl_configs, 
                            (k, Dict(el=>(t1(c[1]),t2(c[2])) for (el,c) in zip(elemlist,conf))) )
                end
            end
        end
    end
    return psl_configs
end


##


function pseudo_mode_name(ps)
    if ps isa String
        return ps
    else
        return "$(ps[1])___"*join(["$(a)-$(c[1])-$(c[2])" for (a,c) ∈ ps[2]],"_")
    end
end


take_part_from_str(str, delim, parts) = join(split(str,delim,keepempty=false)[parts],delim)


function PSL_name_pseudo_mode(name; delim1="___", delim2="_", delim3="-")
    @inline splt(s,d) = split(s,d,keepempty=false)
    @assert startswith(name,"PSL")
    parts = splt(name,delim1)
    part2 = splt(parts[2],delim2)
    ps1   = string(parts[1])
    ps2   = []
    for x in part2
        (a,c1,c2) = splt(x,delim3)
        push!(ps2, string(a)=>(string(c1),string(c2)))
    end
    return (ps1,Dict(ps2))
end



##

#nt(ps,ct,kp,nbd,dg) = "ps___$(pseudo_mode_name(ps))___cut_$(ct[1])_$(ct[2])_kp_$(kp[1]),$(kp[2]),$(kp[3]),$(kp[4]),$(kp[5]),$(kp[6])_nbd_$(nbd)_dg_$(dg)"

#S  = psl_all_possible_configs(["Sc","F"], x->true)

#ps = S[1]

#n = nt(ps, (60.0,480.0), (12,12,12,1,1,1), 44, 0.02)

#n1 = take_part_from_str(n,"___",[2,3])

#p = PSL_name_pseudo_mode(n1)

#make_pseudofile_dict(p,["Sc","F"])
