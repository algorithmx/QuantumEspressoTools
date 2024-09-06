global const WKSPC = "/dg_hpc/CSG/lianyl/STRAINE"
global const iX = parse(Int,ARGS[1])

ENV["JULIA_PKG_SERVER"] = "https://mirrors.tuna.tsinghua.edu.cn/julia"
ENV["PSEUDO_ROOT"] = "/dg_hpc/CSG/lianyl/pseudopotentials"
global const CSNS = "/csns_workspace/CSG/lianyl/softwares"

sleep(40*rand()+(iX%67))
include("$CSNS/QuantumEspressoTools/scripts/scripts.jl")
include("$CSNS/QuantumEspressoTools/scripts/relax.jl")

# * ===============================================================

global const A0_LIST = collect(3.68:0.002:3.8)
global const ryx = 3.0/(1+sqrt(3))
global const C_frac_pos0 = [
    ["C", 0.5/(1+sqrt(3)), 1/3, 0.5],
    ["C", 0.5/(1+sqrt(3)), 2/3, 0.5],
    ["C", 0.5, 1/6, 0.5],
    ["C", 0.5, 5/6, 0.5],
    ["C", (0.5+sqrt(3))/(1+sqrt(3)), 1/3, 0.5],
    ["C", (0.5+sqrt(3))/(1+sqrt(3)), 2/3, 0.5],
]

make_strained_cif(title, a0) = 
    QuantumEspressoTools.minimal_cif(
        title, 
        (a0, a0*ryx, 20.0, 90.0, 90.0, 90.0), 
        C_frac_pos0
    )

## * ===============================================================

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

## * ===============================================================


if iX > length(A0_LIST)
    exit()
else
    a0 = A0_LIST[iX]
    ttl = "strain_a0_$(a0)"
    fn_unrelaxed = "$(WKSPC)/unrelaxed_structures/$(ttl).cif"
    make_strained_cif(ttl, rx, ry)   ⇶   fn_unrelaxed
    if !isfile(fn_unrelaxed)  exit()  end
    @info "\n\n\n---------------------------------\n\nComputation for a0 ( $a0 )\n\n---------------------------------\n\n"
    wkspc0 = "$(WKSPC)/relaxed_$(a0)/"
    wkspc0 = try_mkdir(wkspc0)
    fc_relax(
        wkspc0,                  #*  workspace0::String, 
        ttl,                     #*  common_title::String, 
        fn_unrelaxed,            #*  cif0_fn::String, 
        [("PSL_USPP_PBESOL_SR", Dict("C" => ("n", "1.0.0"))),
          "SSSP_Efficiency"],
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
