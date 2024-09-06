ENV["JULIA_PKG_SERVER"] = "https://mirrors.tuna.tsinghua.edu.cn/julia"
ENV["PSEUDO_ROOT"] = "/home/dabajabaza/abinitio/pseudopotentials"

using BSON
global const WKSPC = "/tmp"
global const iX = 1
global const STRAIN_LIST = [Float64.((1+(dx//100), 1+(dy//100))) for dx=-10:10 for dy=-10:10]
global const C_frac_positions = BSON.load("$WKSPC/C_frac_positions2.bson")[:pos]

global const SFT = "/home/dabajabaza/jianguoyun/Nutstore"
include("$SFT/QuantumEspressoTools/scripts/scripts.jl")
include("$SFT/QuantumEspressoTools/scripts/wannier90.jl")

## * ===============================================================

strained_CELL_PARAMETERS(rx, ry) = [rx*3.745747057, ry*4.506578341, 20.0]

make_strained_relaxed_cif(title,cell,pos) =
    QuantumEspressoTools.minimal_cif(title, (cell...,90.0,90.0,90.0), [("C",k...) for k in pos])

en_kgrid_common(title, psmode) = Dict(
        "assume_isolated" => "2D",
        "title"           => title, 
        "prefix"          => title, 
        "disk_io"         => "medium",
        "conv_thr"        => 1e-11,
        "etot_conv_thr"   => 1e-7,
        "forc_conv_thr"   => 1e-6,
        :pseudo_mode      => psmode,
        :watchdog_setting => kgrid_woof
)

wannier90_wannier90pp_external() = Dict(
    "guiding_centres" => true,
    "wannier_plot" => true,
    #"wannier_plot" => true,
    "dis_num_iter" => 16000
)
wannier90_wannier90_external() = Dict(
    "guiding_centres" => true,
    "wannier_plot" => true,
    #"wannier_plot" => true,
    "dis_num_iter" => 16000
)

## * ===============================================================

const PSM = ("PSL_USPP_PBESOL_SR", Dict("C"=>("n", "1.0.0")))

if iX > length(STRAIN_LIST)
    exit()
else
    (rx,ry) = STRAIN_LIST[iX]
    cell    = strained_CELL_PARAMETERS(rx,ry)
    pos     = C_frac_positions[[rx,ry]]
    ttl     = "strain_rx_$(rx)_ry_$(ry)"
    fn_rlxd = "$(WKSPC)/relaxed_structures/$(ttl).cif"
    make_strained_relaxed_cif(ttl,cell,pos)   ⇶   fn_rlxd
    if !isfile(fn_rlxd)  exit()  end

    @info "\n\n\n---------------------------------\n\nComputation for strain ( $rx, $ry )\n\n---------------------------------\n\n"
    wkspc0 = try_mkdir("$(WKSPC)/band_$(rx)_$(ry)/")

    wannier90_kgrid(
        wkspc0,                          #*  workspace0::String, 
        ttl,                             #*  common_title::String, 
        fn_rlxd,                         #*  cif0_fn::String,
        PSM,                             #*  ps_mode
        (20,16,1,0,0,0),                 #*  kpoint_config::NTuple{6,Int64},
        (60.0,480.0),                    #*  cutoff::Tuple{Float64,Float64},
        32,                              #*  nbnd
        6,                               #*  nwann 
        ["f=$(pos[k][1]),$(pos[k][2]),$(pos[k][3]):pz" for k=1:6],
        [(atom="C$k",l=1,m=1) for k=1:6];
        cleanup     = false,
        PROG_PWX    = `mpiexec -np 8 pw.x -npool 2`,
        PROG_PW2W90 = `mpiexec -np 8 pw2wannier90.x`,
        PROG_W90    = `wannier90.x`,
        PROG_W90PP  = `wannier90.x -pp`,
        PROG_PROJ   = `mpiexec -np 8 projwfc.x`,
    )

end


##


##

