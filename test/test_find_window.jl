using Plots
gr()

ENV["PSEUDO_ROOT"] = "/home/dabajabaza/abinitio/pseudopotentials"
global const SFT = "/home/dabajabaza/jianguoyun/Nutstore"
include("$SFT/QuantumEspressoTools/scripts/scripts.jl")

##* ====================================

global const WKSPC = "/home/dabajabaza/jianguoyun/Nutstore/Biphenylene"
include("$WKSPC/include/Biphenylene.jl")

global const FD_STRAIN = "/mnt/dg_hpc/STRAIN2"
global const ps1 = "ps___PSL_USPP_PBESOL_SR___C-n-1.0.0___kp_20,16,1,0,0,0_cutoff_60.0,480.0_kBar_0.0_dg_0.02"
global const ps2 = "ps___PSL_USPP_PBESOL_SR___C-n-1.0.0___kp_20,16,1,0,0,0_kp_20,16,1,0,0,0_cutoff_60.0,480.0"
hr_file(rx,ry)   = "$FD_STRAIN/relaxed_$(rx)_$(ry)/$ps2/strain_rx_$(rx)_ry_$(ry)_wannier90/strain_rx_$(rx)_ry_$(ry)_hr.dat"
cif_file(rx,ry,i=3)  = "$FD_STRAIN/relaxed_$(rx)_$(ry)/$ps1/strain_rx_$(rx)_ry_$(ry)_fc_relax_iter_$(i).cif"
xml_file(rx,ry)  = "$FD_STRAIN/relaxed_$(rx)_$(ry)/$ps2/strain_rx_$(rx)_ry_$(ry)_scf/strain_rx_$(rx)_ry_$(ry).band.xml"
scf_file(rx,ry)  = "$FD_STRAIN/relaxed_$(rx)_$(ry)/$ps2/strain_rx_$(rx)_ry_$(ry)_band/strain_rx_$(rx)_ry_$(ry).pw.x.out"
scf0_file(rx,ry)  = "$FD_STRAIN/relaxed_$(rx)_$(ry)/$ps2/strain_rx_$(rx)_ry_$(ry)_scf/strain_rx_$(rx)_ry_$(ry).pw.x.out"
nscf_file(rx,ry)  = "$FD_STRAIN/relaxed_$(rx)_$(ry)/$ps2/strain_rx_$(rx)_ry_$(ry)_nscf/strain_rx_$(rx)_ry_$(ry).pw.x.out"
wout_file(rx,ry) = "$FD_STRAIN/relaxed_$(rx)_$(ry)/$ps2/strain_rx_$(rx)_ry_$(ry)_wannier90/strain_rx_$(rx)_ry_$(ry).wout"
proj_file(rx,ry) = "$FD_STRAIN/relaxed_$(rx)_$(ry)/$ps2/strain_rx_$(rx)_ry_$(ry)_band_projwfc/strain_rx_$(rx)_ry_$(ry).projwfc.x.out"
up_file(rx,ry)   = "$FD_STRAIN/relaxed_$(rx)_$(ry)/$ps2/strain_rx_$(rx)_ry_$(ry)_band_projwfc/strain_rx_$(rx)_ry_$(ry).proj.projwfc_up"
proj_file1(rx,ry) = "$FD_STRAIN/relaxed_$(rx)_$(ry)/$ps2/strain_rx_$(rx)_ry_$(ry)_projwfc/strain_rx_$(rx)_ry_$(ry).projwfc.x.out"
up_file1(rx,ry)   = "$FD_STRAIN/relaxed_$(rx)_$(ry)/$ps2/strain_rx_$(rx)_ry_$(ry)_projwfc/strain_rx_$(rx)_ry_$(ry).proj.projwfc_up"

##* ====================================

E_fermi = pw_fermi_energy_eV(readlines(scf0_file(0.9,1.02)))

OCC_ABS, OCC_REL = calculate_occupation(
    readlines(proj_file1(0.9,1.02)),
    [(atom="C$(k)",l=1,m=1,orbit="2P") for k=1:6],
    readlines(up_file1(0.9,1.02))
) ;

BANDS =  hcat(last.(values.(pw_bands(readlines(nscf_file(0.9,1.02))))) ...)' ;

OCCUPATION_POSSIBLE = Array{Bool,2}(OCC_ABS.>0.1) ;

for ee = 0.05 : 0.1 : 5.0
    @info findgap(BANDS[:], OCCUPATION_POSSIBLE[:], E_fermi, ee)
end


outer_min, outer_max, outer_num_en = findgap(BANDS[:], OCCUPATION_POSSIBLE[:], E_fermi, 0.3)