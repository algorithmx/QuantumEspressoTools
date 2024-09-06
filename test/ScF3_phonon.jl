#include("/public3/home/sc55341/softwares/QuantumEspressoTools/scripts/scripts.jl")
include("/home/dabajabaza/jianguoyun/Workspace/QuantumEspressoTools/scripts/scripts.jl")
include("/home/dabajabaza/jianguoyun/Workspace/QuantumEspressoTools/scripts/test_dry_run.jl")

minimal_cif("ScF3", (4.06, 4.06, 4.06, 90, 90, 90, 221), 
            [("Sc", 0.0, 0.0, 0.0), ("F", 1/2, 0, 0), ("F", 0, 1/2, 0), ("F", 0, 0, 1/2)]
            )  ⇶ "/data/ScF3.0.cif"

##

psl_psmodes = []

@inline t1(x) = replace(x,"-"=>"")
@inline t2(x) = replace(replace(x,".UPF"=>""),"psl."=>"")


for k in keys(PSEUDO_DATA)
    if startswith(k,"PSL") && endswith(k,"SR")
        for xf in keys(PSEUDO_DATA[k]["F"])
            if "Sc" ∈ keys(PSEUDO_DATA[k])
                for xsc in keys(PSEUDO_DATA[k]["Sc"])
                    push!(psl_psmodes, (k,Dict("F"=>(t1(xf[1]),t2(xf[2])),"Sc"=>(t1(xsc[1]),t2(xsc[2])))))
                end
            end
        end
    end
end

##

test_dry_run(
    "/data", 
    "ScF3", 
    "/data/ScF3.0.cif", 
    psl_psmodes,
    [(8,8,8,1,1,1), (12,12,12,1,1,1)];
    PROG_PWX=`pw.x`,
    atoms = ["Sc", "F"]
)
