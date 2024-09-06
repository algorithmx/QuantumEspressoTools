global const WKSPC = "/mnt/dg_hpc/MXene/Mo2N/ps___SG15___kp_16,16,1,0,0,0_kp_16,16,1,0,0,0_cut_100.0,400.0"
ENV["PSEUDO_ROOT"] = "/home/dabajabaza/abinitio/pseudopotentials"
global const SFT = "/home/dabajabaza/jianguoyun/Nutstore"
include("$SFT/QuantumEspressoTools/scripts/scripts.jl")

##

pwfn = "$WKSPC/Mo2N_scf/Mo2N.pw.x.out"

##

QuantumEspressoTools.pw_fermi_energy_eV(
    readlines(pwfn)
)
