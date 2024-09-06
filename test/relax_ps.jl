global const WKSPC = "/dg_hpc/CSG/lianyl/SHEAR1"
ENV["PSEUDO_ROOT"] = "/dg_hpc/CSG/lianyl/pseudopotentials"
global const SFT = "/csns_workspace/CSG/lianyl/softwares"


global const PS_LIST = [
    ("PSL_USPP_PBESOL_SR", Dict("C" => ("n", "1.0.0"))),
    "SSSP_Efficiency",
    "SSSP_Precision",
    "SG15",
    "GBRV_PBESOL"
]
global const CTFF_LIST = [
    (60.0,480.0),
    (45.0,360.0),
    (45.0,360.0),
    (80.0,640.0),
    (45.0,360.0),
]

global const iX = parse(Int,ARGS[1])


global const a0 = 3.745747057
global const b0 = 4.506578341
global const c0 = 20.0


sleep(30*rand()+(iX%57))
ENV["JULIA_PKG_SERVER"] = "https://mirrors.tuna.tsinghua.edu.cn/julia"
include("$SFT/QuantumEspressoTools/scripts/scripts.jl")
include("$SFT/QuantumEspressoTools/scripts/relax.jl")

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

        # spin-orbit
        #"lspinorb"        => occursin("_FR",pseudo_mode_name(psmode)),
        "lspinorb"        => false,

        # relax
        "smearing"        => "mp",
        "degauss"         => dg,
        "press"           => Float64(kBar),

        # watch dog
        :watchdog_setting => relax_woof
)


# * ===============================================================

global const C_frac_pos0 = [
    ["C", 0.8065498337, 0.3385975908, 0.5000000000],
    ["C", 0.1934501663, 0.6614023092, 0.5000000000],
    ["C", 0.5000000000, 0.8401238552, 0.5000000000],
    ["C", 0.5000000000, 0.1598761448, 0.5000000000],
    ["C", 0.8065498337, 0.6614023092, 0.5000000000],
    ["C", 0.1934501663, 0.3385975908, 0.5000000000],
]

make_sheared_cif(title) = 
    QuantumEspressoTools.minimal_cif(title, (a0, b0, 20.0, 90.0, 90.0, 90.0), C_frac_pos0)


## * ===============================================================

if iX > length(DA_LIST)
    exit()
else
    ps  = PS_LIST[iX]
    ct  = CTFF_LIST[iX]
    ttl = "shear_rx_1.0_ry_1.0"
    fn_unrelaxed = "$(WKSPC)/unrelaxed_structures/$(ttl).cif"
    make_sheared_cif(ttl)  ⇶  fn_unrelaxed
    if !isfile(fn_unrelaxed)  exit()  end

    @info "\n\n\n---------------------------------\n\nComputation for pseudopotential and cutoff ( $ps, $ct )\n\n---------------------------------\n\n"
    wkspc0 = "$(WKSPC)/relaxed_$(rx)_$(ry)/"
    wkspc0 = try_mkdir(wkspc0)
    vc_relax(
        "",
        wkspc0,                  #*  workspace0::String, 
        ttl,                     #*  common_title::String, 
        fn_unrelaxed,            #*  cif0_fn::String, 
        [ps,],                   #*  ps_modes::Vector,
        [(10, 8,1,0,0,0),
         (20,16,1,0,0,0)],       #*  kpoint_configs::Vector{NTuple{6,Int64}},
        [ct,],                   #*  cutoff_list::Vector{Tuple{Float64,Float64}},
        [0.0,],                  #*  kBar_list::Vector{Float64},
        [0.02,];                 #*  degauss_list::Vector{Float64};
        atoms = ["C",],    
        beta = [0.8, 0.5],
        PROG_PWX = `mpirun -np 48 pw.x -npool 8`,
        relax_iter_updater = x->Dict(),
    )
end
