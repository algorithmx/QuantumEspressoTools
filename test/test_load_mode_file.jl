include("/home/dabajabaza/jianguoyun/Workspace/QuantumEspressoTools/scripts/scripts.jl")


#

fn_freq = "/home/dabajabaza/Downloads/ScF3_PSL_PAW_PBE_SR.frq"
ff = parse_energy_file(readlines(fn_freq))

##

fn_eig  = "/home/dabajabaza/Downloads/ScF3_PSL_PAW_PBE_SR.eig"
qq = parse_mode_file(readlines(fn_eig))
