## header

global const WKSPC = "/dg_hpc/CSG/lianyl/MoO3_relax"
cd(WKSPC)

global const iX = parse(Int,ARGS[1])



PS_LIST = [
    "SSSP_Precision",
    "SSSP_Efficiency",
    "GBRV_PBE",
    "GBRV_PBESOL",
    "SG15",
    ("PSL_USPP_PBE_SR",    Dict("Mo"=>("spn","1.0.0"), "O"=>("n","1.0.1"))),
    ("PSL_USPP_PBESOL_SR", Dict("Mo"=>("spn","1.0.0"), "O"=>("n","1.0.0"))),
    ("PSL_PAW_PBE_SR",     Dict("Mo"=>("spn","1.0.0"), "O"=>("n","1.0.0"))),
    ("PSL_PAW_PBESOL_SR",  Dict("Mo"=>("spn","1.0.0"), "O"=>("n","1.0.0"))),
]

CT_LIST = [
    (50.0,400.0), 
    (60.0,480.0), 
    (80.0,640.0)
]

KP_LIST = [
    (4,12,4,0,0,0),
    (6,18,6,0,0,0),
]


(i, j, k) = (0, 0, 0)

if iX > length(PS_LIST) * length(CT_LIST) * length(KP_LIST)
    exit()
else
    (i, j, k) = ( ((iX-1) % length(KP_LIST)) + 1,
                 (((iX-1)÷length(KP_LIST)) % length(CT_LIST)) + 1,
                 (((iX-1)÷(length(KP_LIST)*length(CT_LIST))) % length(PS_LIST)) + 1, )
    @info (KP_LIST[i], CT_LIST[j], PS_LIST[k])
end

ENV["JULIA_PKG_SERVER"] = "https://mirrors.tuna.tsinghua.edu.cn/julia"
ENV["PSEUDO_ROOT"] = "/dg_hpc/CSG/lianyl/pseudopotentials"
global const CSNS = "/csns_workspace/CSG/lianyl/softwares"

sleep(30*rand()+(iX%60))
include("$CSNS/QuantumEspressoTools/scripts/scripts.jl")
include("$CSNS/QuantumEspressoTools/scripts/phonon_freq.jl")

## configs 

MoO3_settings = copy(phonon_freq_additional_settings)
MoO3_settings["atoms"] = ["Mo","O"]

## structure

cif_fn = "$WKSPC/Mo4O12.VASP.cif"

if !isfile(cif_fn)
    minimal_cif(
        "Mo4O12", 
        (3.9624, 13.86, 3.6971, 90, 90, 90, 0), 
        [("Mo", 0.08502, 0.10132, 0.25000),
        ("Mo", 0.91496, 0.89867, 0.75000),
        ("Mo", 0.58503, 0.39867, 0.75000),
        ("Mo", 0.41497, 0.60132, 0.25000),
        ("O",  0.03480, 0.22120, 0.25000),
        ("O",  0.96520, 0.77880, 0.75000),
        ("O",  0.53479, 0.27880, 0.75000),
        ("O",  0.46520, 0.72119, 0.25000),
        ("O",  0.52109, 0.08806, 0.25000),
        ("O",  0.47890, 0.91193, 0.75000),
        ("O",  0.02109, 0.41192, 0.75000),
        ("O",  0.97890, 0.58806, 0.25000),
        ("O",  0.50190, 0.43610, 0.25000),
        ("O",  0.49809, 0.56389, 0.75000),
        ("O",  0.00190, 0.06389, 0.75000),
        ("O",  0.99809, 0.93610, 0.25000),]
    )  ⇶ "$WKSPC/Mo4O12.VASP.cif"
end

## common settings + DFT-D3

phonon_freq_common(title, psmode, kp) = Dict(
        "title"           => title, 
        "prefix"          => title, 
        "vdw_corr"        => "DFT-D3", 
        :pseudo_mode      => psmode,
        :pw_mode          => "energy",
        :kpoint_mode      => "automatic", 
        :kpoints          => kp,
        :watchdog_setting => phonon_freq_scf_woof
)

## MAIN


phonon_freq(
    WKSPC,               #*  workspace0::String, 
    "Mo4O12",            #*  common_title::String, 
    cif_fn,              #*  cif0_fn::String, 
    [PS_LIST[k],],       #*  ps_modes::Vector,
    [KP_LIST[i],],       #*  kpoint_configs,
    [CT_LIST[j],],       #*  cutoff_list,,
    [[0,0,0],
        [1/2,0,0],[0,1/2,0],[0,0,1/2],
        [1/2,1/2,0],[1/2,0,1/2],[0,1/2,1/2],
        [1/2,1/2,1/2],],    #*  qpoints
    [0.0, 2.0],          #*  kBar_list::Vector{Float64},
    [0.02];              #*  degauss_list::Vector{Float64};
    tr2_ph              = 1e-14,
    conv_thr_scf        = 1e-11,
    scf_cutoff_upscale  = 1.1,
    PROG_PWX            = `mpiexec -np 48  pw.x -npool 6`,
    PROG_PHX            = `mpiexec -np 48  ph.x -npool 6`,
    additional_settings = MoO3_settings
)

