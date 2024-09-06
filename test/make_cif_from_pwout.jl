
ENV["JULIA_PKG_SERVER"] = "https://mirrors.tuna.tsinghua.edu.cn/julia"
ENV["PSEUDO_ROOT"] = "/home/dabajabaza/abinitio/pseudopotentials"

global const SFT = "/home/dabajabaza/jianguoyun/Nutstore"
include("$SFT/QuantumEspressoTools/scripts/scripts.jl")

##

function make_cif(title, lines)
    positions = QuantumEspressoTools.pw_atom_list(lines)
    cellparams = QuantumEspressoTools.pw_crystal_axes_Angstrom(lines,quiet=true)
    lattp = QuantumEspressoTools.basis_to_lattice_parameters(cellparams)
    pp = (lattp.p[1:6]...,)
    return QuantumEspressoTools.minimal_cif(title, pp, positions)
end

##

WS2_PSL_PAW_PBE_SR_pwout = readlines("/mnt/dg_hpc/WS2/ps___PSL_PAW_PBE_SR___S-n-1.0.0_W-spn-1.0.0___kp_6,6,8,0,0,0_cut_80.0,640.0_kBar_0.0_dg_0.02/WS2_scf/WS2.pw.x.out")

WS2_PSL_PAW_PBE_SR_cif = make_cif("WS2_PSL_PAW_PBE_SR", WS2_PSL_PAW_PBE_SR_pwout) ;

WS2_PSL_PAW_PBE_SR_cif ⇶ "WS2.PSL_PAW_PBE_SR.cif"

##

WS2_SSSP_Precision_pwout = readlines("/mnt/dg_hpc/WS2/ps___SSSP_Precision___kp_8,8,10,0,0,0_cut_80.0,640.0_kBar_0.0_dg_0.02/WS2_scf/WS2.pw.x.out")

WS2_SSSP_Precision_cif = make_cif("WS2_SSSP_Precision", WS2_SSSP_Precision_pwout) ;

WS2_SSSP_Precision_cif ⇶ "WS2.SSSP_Precision.cif"

##

WS2_SSSP_Precision_pwout = readlines("/mnt/dg_hpc/WS2/ps___SSSP_Precision_PBEsol___kp_8,8,10,0,0,0_cut_80.0,640.0_kBar_0.0_dg_0.02/WS2_scf/WS2.pw.x.out")

WS2_SSSP_Precision_cif = make_cif("SSSP_Precision_PBEsol", WS2_SSSP_Precision_pwout) ;

WS2_SSSP_Precision_cif ⇶ "WS2.SSSP_Precision_PBEsol.cif"
