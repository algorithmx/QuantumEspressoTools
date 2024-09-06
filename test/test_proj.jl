using Plots
gr()

#global const WKSPC = "/tmp/band_0.9_0.9/ps___PSL_USPP_PBESOL_SR___C-n-1.0.0___kp_20,16,1,0,0,0_cut_60.0,480.0"
global const WKSPC = "/mnt/dg_hpc/MXene/Mo2N/ps___SG15___kp_16,16,1,0,0,0_kp_16,16,1,0,0,0_cut_100.0,400.0"
#global const WKSPC = "/mnt/dg_hpc/BANDWU/band_1.0_1.0/ps___PSL_USPP_PBESOL_SR___C-n-1.0.0___kp_20,16,1,0,0,0_kp_20,16,1,0,0,0_cut_60.0,480.0"
ENV["PSEUDO_ROOT"] = "/home/dabajabaza/abinitio/pseudopotentials"
global const SFT = "/home/dabajabaza/jianguoyun/Nutstore"
include("$SFT/QuantumEspressoTools/scripts/scripts.jl")

##

# projfn = "$WKSPC/strain_rx_0.9_ry_0.9_projwfc/strain_rx_0.9_ry_0.9.proj.projwfc_up"
# bandfn = "$WKSPC/strain_rx_0.9_ry_0.9_nscf/strain_rx_0.9_ry_0.9.pw.x.out" 

projfn = "$WKSPC/Mo2N_projwfc/Mo2N.proj.projwfc_up"
proj_out_fn = "$WKSPC/Mo2N_projwfc/Mo2N.projwfc.x.out"
bandfn = "$WKSPC/Mo2N_nscf/Mo2N.pw.x.out"
band_highsymm_fn = "$WKSPC/Mo2N_band/Mo2N.pw.x.out"
scffn = "$WKSPC/Mo2N_scf/Mo2N.pw.x.out"
band_w90_fn = "$WKSPC/Mo2N_wannier90/Mo2N_band.dat"

# projfn = "$WKSPC/strain_rx_1.0_ry_1.0_projwfc/strain_rx_1.0_ry_1.0.proj.projwfc_up"
# proj_out_fn = "$WKSPC/strain_rx_1.0_ry_1.0_projwfc/strain_rx_1.0_ry_1.0.projwfc.x.out"
# bandfn = "$WKSPC/strain_rx_1.0_ry_1.0_nscf/strain_rx_1.0_ry_1.0.pw.x.out"


##% ===========================================================================

#% check weights at fermi energy

##% ===========================================================================

E_fermi = QuantumEspressoTools.pw_fermi_energy_eV(readlines(scffn))

P, all_states = projwfcx_output_projwfc_up(readlines(projfn)) ;

bands = hcat(last.(pw_bands(readlines(bandfn)))...)' ;

close_to_Fermi = findall((bands .> E_fermi-0.1) .* (bands .< E_fermi+0.1))

weight_close_to_Fermi = hcat([P[:,:,i][close_to_Fermi] for i=1:length(all_states)]...) ;

@inline snd(x) = x[2]

states_close_to_Fermi = findall(weight_close_to_Fermi.>0.02) .|> snd |> unique

##% ===========================================================================

all_states[states_close_to_Fermi]

##% ===========================================================================


bands1 = hcat(last.(pw_bands(readlines(band_highsymm_fn)))...)' ;


##

l2f(x) = parse.(Float64,split(x," ",keepempty=false)[1:2])

L = readlines(band_w90_fn) ;

SECTIONS = []
B = Float64[]
for l in L   
    if strip(l)==""
        push!(SECTIONS, B)
        global B = Float64[]
    else
        push!(B, l2f(l)[2])
    end
end

##

plt = plot()

for i=6:30
    #plt = plot!(bands[:,i],color="black",legends=nothing)
end

for s in SECTIONS
    plt = plot!(s,color="red",legends=nothing)
end

ylims!(-22.0,8.0)

plt

##

projector_proj = [
    [(atom="Mo1",l=2,m=i) for i=1:5]...,
    [(atom="Mo2",l=2,m=i) for i=1:5]...,
    [(atom="N1",l=1,m=i)  for i=1:3]...   ### pz only
]

OCC_ABS, OCC_REL = calculate_occupation(
    readlines(proj_out_fn), # 
    projector_proj,
    readlines(projfn)
) ;

##

close_to_Fermi = findall((bands .> E_fermi-0.1) .* (bands .< E_fermi+0.1))

##



##

find_window0(
    readlines(proj_out_fn),
    60,
    13,
    projector_proj,
    readlines(bandfn), 
    readlines(projfn); 
    thr = (0.6, 0.2),
    E_ref = -3.8617,
    dE = 0.5
) ;

##

#! good idea ~
MASK = OCC_ABS.>0.05
scatter(OCC_ABS[MASK], bands[MASK],c="black")
scatter!(OCC_REL[MASK], bands[MASK],c="black")

##

#! good idea ~
MASK1 = (OCC_REL.>0.99) .& MASK
scatter!(OCC_ABS[MASK1], bands[MASK1],c="red")
scatter!(OCC_REL[MASK1], bands[MASK1],c="red", legends=nothing)

##

sum(MASK1)

##

#! good idea ~
MASK1 = (OCC_REL.>0.5) .& MASK
scatter!(OCC_ABS[MASK1], bands[MASK1],c="blue")
scatter!(OCC_REL[MASK1], bands[MASK1],c="blue")

##

plot!(0.01.*collect(0:100), 0.0.*collect(0:100) .- 8.4)

##

histogram(OCC_REL[OCC_ABS.<0.2])

##

occ_pz(i,j) = sum(proj[i,j,2:4:24])
occ_pzi(i,j,z) = sum(proj[i,j,[z...,]])

PZ_ord = [sortperm([-occ_pzi(i,j,2:4:24) for j=1:32]) for i=1:320] ;
PZ_occ = [[occ_pzi(i,j,2:4:24) for j=1:32] for i=1:320] ;


PZ_occ[22][PZ_ord[22][1:16]]

PZ_ord[9][1:12] |> sort

bands[PZ_ord[22][1:16], 22]
bands[PZ_ord[23][1:16], 23] |> sort
