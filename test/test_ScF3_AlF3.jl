#include("/public3/home/sc55341/softwares/QuantumEspressoTools/scripts/scripts.jl")

include("/home/dabajabaza/jianguoyun/Workspace/QuantumEspressoTools/scripts/scripts.jl")
include("/home/dabajabaza/jianguoyun/Workspace/QuantumEspressoTools/scripts/test_phonon_freq.jl")

#minimal_cif("ScF3", (3.8,3.8,3.8, 90,90,90, 221), 
#            [("Sc", 0.0, 0.0, 0.0), ("F", 1/2, 0, 0), ("F", 0, 1/2, 0), ("F", 0, 0, 1/2)]) ⇶ "/data/ScF3/ScF3.0.cif"


##

test_phonon_freq(
    "/data/ScF3", 
    "ScF3", 
    "/data/ScF3/ScF3.0.cif", 
    ["GBRV_LDA"],
    [(4,4,4,1,1,1)],
    [[0,0,0],[1/2,0,0],], 
    [0.0],
    [0.02];
    PROG_PWX=`mpiexec -np 4 pw.x -npool 2`,
    PROG_PHX=`mpiexec -np 4 ph.x -npool 2`
)

##


