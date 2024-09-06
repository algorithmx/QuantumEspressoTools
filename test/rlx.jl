global const WKSPC = "/dg_hpc/CSG/lianyl/STRAIN2"

global const iX = parse(Int,ARGS[1])

ENV["JULIA_PKG_SERVER"] = "https://mirrors.tuna.tsinghua.edu.cn/julia"

ENV["PSEUDO_ROOT"] = "/dg_hpc/CSG/lianyl/pseudopotentials"

global const CSNS = "/csns_workspace/CSG/lianyl/softwares"

sleep(30*rand()+(iX%60))
include("$CSNS/QuantumEspressoTools/scripts/scripts.jl")
include("$CSNS/QuantumEspressoTools/scripts/relax.jl")

# * ===============================================================

global const a0 = 3.745747057
global const b0 = 4.506578341
global const c0 = 20.0

# * ===============================================================

global const STRAIN_LIST = [Float64.((1+(dx//800), 1+(dy//800))) 
                              for dx=-8:8 for dy=-8:8]

# * ===============================================================

relax_common(title, psmode, kp, kBar, dg) = Dict(
        # common
        "assume_isolated" => "2D",
        "title"           => title, 
        "prefix"          => title, 
        :pseudo_mode      => psmode,
        # k-points
        :kpoint_mode      => "automatic", 
        :kpoints          => kp,
        "conv_thr"        => 1e-9,
        "etot_conv_thr"   => 1e-7,
        "forc_conv_thr"   => 1e-6,
        "lspinorb"        => false,
        # relax
        "smearing"        => "mp",
        "degauss"         => dg,
        "press"           => Float64(kBar),
        # watch dog
        :watchdog_setting => relax_woof
)

# * ===============================================================

pos_1_x_k(k) = [
  0.1936036265283593,
  0.010408359838014497,
  0.02748106408330539,
  0.8841096947236762,
 -0.18968454050429218,
 -0.11137417276031691,
]' * [1.0,k[1],k[2],k[1]^2,k[2]^2,k[1]*k[2]]
pos_1_x(rx,ry) = round(pos_1_x_k([rx-1,ry-1]), digits=5)

pos_1_y_k(k) = [
  0.3387588461972236,
  0.015223611804431339,
  0.07088088201377858,
  0.6752142122265306,
 -0.1753490323154585,
  0.11712643815402868,
]' * [1.0,k[1],k[2],k[1]^2,k[2]^2,k[1]*k[2]]
pos_1_y(rx,ry) = round(pos_1_y_k([rx-1,ry-1]), digits=5)

pos_3_y_k(k) = [
  0.16001417344658708,
  0.08132767536064613,
 -0.018263070233655823,
 -0.029510020757950395,
  0.34118606861087974,
  0.1264237827571706,
]' * [1.0,k[1],k[2],k[1]^2,k[2]^2,k[1]*k[2]]
pos_3_y(rx,ry) = round(pos_3_y_k([rx-1,ry-1]), digits=5)

C_frac_pos(rx,ry) = [
    ["C",pos_1_x(rx,ry),   pos_1_y(rx,ry), 0.5],
    ["C",pos_1_x(rx,ry), 1-pos_1_y(rx,ry), 0.5],
    ["C",0.5,   pos_3_y(rx,ry), 0.5],
    ["C",0.5, 1-pos_3_y(rx,ry), 0.5],
    ["C",1-pos_1_x(rx,ry),   pos_1_y(rx,ry), 0.5],
    ["C",1-pos_1_x(rx,ry), 1-pos_1_y(rx,ry), 0.5],
]

make_strained_cif(title, rx, ry) = 
    QuantumEspressoTools.minimal_cif(
        title, 
        (a0*rx, b0*ry, c0, 90.0, 90.0, 90.0), 
        C_frac_pos(rx,ry)
    )

## * ===============================================================


if iX > length(STRAIN_LIST)
    exit()
else
    (rx,ry) = STRAIN_LIST[iX]
    ttl = "strain_rx_$(rx)_ry_$(ry)"
    fn_unrelaxed = "$(WKSPC)/unrelaxed_structures/$(ttl).cif"
    make_strained_cif(ttl, rx, ry)   ⇶   fn_unrelaxed
    if !isfile(fn_unrelaxed)  exit()  end
    @info "\n\n\n---------------------------------\n\nComputation for strain ( $rx, $ry )\n\n---------------------------------\n\n"
    wkspc0 = "$(WKSPC)/relaxed_$(rx)_$(ry)/"
    wkspc0 = try_mkdir(wkspc0)
    sleep(30*rand())
    fc_relax(
        wkspc0,                  #*  workspace0::String, 
        ttl,                     #*  common_title::String, 
        fn_unrelaxed,            #*  cif0_fn::String, 
        [("PSL_USPP_PBESOL_SR", Dict("C" => ("n", "1.0.0"))),],
                                 #*  ps_modes::Vector,
        [(20,16,1,0,0,0)],       #*  kpoint_configs::Vector{NTuple{6,Int64}},
        [(60.0,480.0),],         #*  cutoff_list::Vector{Tuple{Float64,Float64}},
        [0.0,],                  #*  kBar_list::Vector{Float64},
        [0.02,];                 #*  degauss_list::Vector{Float64};
        atoms = ["C",],    
        beta = [0.8, 0.5],
        PROG_PWX = `mpirun -np 48 pw.x -npool 8`,
        relax_iter_updater = x->Dict(), 
    )
end
