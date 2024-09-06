using Pkg
Pkg.activate("/csns_workspace/CSG/lianyl/softwares/QuantumEspressoTools")
using QuantumEspressoTools
QT = QuantumEspressoTools

vf(rx,ry) = (
try 
    QT.verify_result(
"./relaxed_$(rx)_$(ry)/ps___PSL_USPP_PBESOL_SR___C-n-1.0.0___kp_30,24,1,0,0,0_cutoff_60.0,480.0_kBar_0.0_dg_0.02/strain_rx_$(rx)_ry_$(ry)_fc_relax_iter_1/strain_rx_$(rx)_ry_$(ry).pw.x.out" |> readlines, 
`pw.x`)
catch _; 
    false 
end )

cd("/dg_hpc/CSG/lianyl/STRAIN2")
P = [(i,j) for i=0.9:0.01:1.1 for j=0.9:0.01:1.1] ;
STR2_DONE = Dict((rx,ry) => vf(rx,ry) for (rx,ry) in P) ;

undone = length(P) - sum(values(STR2_DONE))

@info "undone = $undone"
for k in keys(STR2_DONE)
    if !STR2_DONE[k]
        @info k
    end
end

using BSON
try rm("STR2_DONE.bson") catch _; nothing end
BSON.bson("STR2_DONE.bson",STR2_DONE=STR2_DONE)
exit()
