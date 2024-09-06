global const WKSPC = "/dg_hpc/CSG/lianyl/SHEAR1"
ENV["PSEUDO_ROOT"] = "/dg_hpc/CSG/lianyl/pseudopotentials"
global const SFT = "/csns_workspace/CSG/lianyl/softwares"

global const psmode = ("PSL_USPP_PBESOL_SR", Dict("C" => ("n", "1.0.0")))
global const DA_LIST = 0.005collect(-10:10)
global const iX = parse(Int,ARGS[1])

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


function make_sheared_cif(title, lines, δa)
    pospairs = QuantumEspressoTools.pw_relax_atomic_positions(lines,quiet=true)
    positions = [(a[1],a[2]...) for a in pospairs[2]]
    cellparams = QuantumEspressoTools.pw_relax_cell_parameters(lines,quiet=true)
    lattp = QuantumEspressoTools.basis_to_lattice_parameters(cellparams)
    pp = lattp.p
    @assert pp[4] ≈ pp[5] ≈ pp[6] ≈ 90.0
    pp1 = (pp[1], sqrt(pp[2]^2+(δa*pp[1])^2), pp[3], 90.0, 90.0, 90-(180atan(δa*pp[1],pp[2])/π))
    return QuantumEspressoTools.minimal_cif(title, pp1, positions)
end


# * ===============================================================

global const STRUCT0 = """
CELL_PARAMETERS (angstrom)
   3.745747057   0.000000000   0.000000000
   0.000000000   4.506578341   0.000000000
   0.000000000   0.000000000  20.000000000

ATOMIC_POSITIONS (crystal)
C             0.8065498337        0.3385975908        0.5000000000
C             0.1934501663        0.6614023092        0.5000000000
C             0.5000000000        0.8401238552        0.5000000000
C             0.5000000000        0.1598761448        0.5000000000
C             0.8065498337        0.6614023092        0.5000000000
C             0.1934501663        0.3385975908        0.5000000000
""" |> SPLTN

## * ===============================================================

if iX > length(DA_LIST)
    exit()
else
    δa  = DA_LIST[iX]
    ttl = "shear_rx_$(rx)_ry_$(ry)"
    fn_unrelaxed = "$(WKSPC)/unrelaxed_structures/$(ttl).cif"
    make_sheared_cif(ttl, STRUCT0, δa)  ⇶  fn_unrelaxed
    if !isfile(fn_unrelaxed)  exit()  end

    @info "\n\n\n---------------------------------\n\nComputation for shear ( $alpha )\n\n---------------------------------\n\n"
    wkspc0 = "$(WKSPC)/relaxed_$(rx)_$(ry)/"
    wkspc0 = try_mkdir(wkspc0)
    fc_relax(
        wkspc0,                  #*  workspace0::String, 
        ttl,                     #*  common_title::String, 
        fn_unrelaxed,            #*  cif0_fn::String, 
        [psmode,],               #*  ps_modes::Vector,
        [(30,24,1,0,0,0)],       #*  kpoint_configs::Vector{NTuple{6,Int64}},
        [(60.0,480.0),],         #*  cutoff_list::Vector{Tuple{Float64,Float64}},
        [0.0,],                  #*  kBar_list::Vector{Float64},
        [0.02,];                 #*  degauss_list::Vector{Float64};
        atoms = ["C",],    
        beta = [0.8, 0.5],
        PROG_PWX = `mpirun -np 48 pw.x -npool 8`,
        relax_iter_updater = x->Dict(),
    )
end
