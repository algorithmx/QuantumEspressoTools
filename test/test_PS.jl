include("/home/dabajabaza/jianguoyun/Workspace/QuantumEspressoTools/scripts/scripts.jl")

##

num_KS_states(("PSL_PAW_PZ_SR",Dict("O"=>("n","0.1"),"Re"=>("spn","1.0.0"))), ["O","O","O","Re"])

num_KS_states(("GHH_PZ"), ["F","F","F","Sc"])


num_KS_states(("SSSP_Precision"), ["O","O","O","Re"])

num_KS_states(("SSSP_Precision"), ["F","F","F","Sc"])

##

Find_suggested_cutoff(PSEUDO_PATHS["SG15"]*"/F_ONCV_PBE-1.2.upf")
Find_suggested_cutoff(PSEUDO_PATHS["SG15"]*"/Sc_ONCV_PBE-1.2.upf")
Find_suggested_cutoff(PSEUDO_PATHS["SG15"]*"/O_ONCV_PBE-1.2.upf")


Find_suggested_cutoff("PSL_PAW_PZ_SR", PSEUDO_PATHS["PSL_PAW_PZ_SR"]*"/F.pz-n-kjpaw_psl.0.1.UPF")
Find_suggested_cutoff("PSL_PAW_PZ_SR", PSEUDO_PATHS["PSL_PAW_PZ_SR"]*"/Sc.pz-spn-kjpaw_psl.1.0.0.UPF")


PSEUDO_DATA["SSSP_Precision"]["F"]
PSEUDO_DATA["SSSP_Precision"]["Mn"]

make_pseudo_ecut(("PSL_PAW_PBE_SR",Dict("Re"=>("spn","1.0.0"),"O"=>("n","1.0.0"))), ["O","Re"])
