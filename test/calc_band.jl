global const WKSPC = "/dg_hpc/CSG/lianyl/SHEAR2"
ENV["PSEUDO_ROOT"] = "/dg_hpc/CSG/lianyl/pseudopotentials"
global const SFT = "/csns_workspace/CSG/lianyl/softwares"
global const DA_LIST = 0.005collect(-10:10)
global const iX = parse(Int,ARGS[1])

sleep(30*rand()+(iX%57))
include("$SFT/QuantumEspressoTools/scripts/scripts.jl")
include("$SFT/QuantumEspressoTools/scripts/band_wannier90.jl")

using BSON
global const C_frac_positions = BSON.load("$WKSPC/C_frac_positions_shear.bson")[:pos]

## * ===============================================================

function make_sheared_cif(title, δa, positions)
    pp = [3.745747057, 4.506578341, 20.000000000]
    pp1 = (pp[1], sqrt(pp[2]^2+(δa*pp[1])^2), pp[3], 90.0, 90.0, 90-(180atan(δa*pp[1],pp[2])/π))
    @inline addC(x) = ["C",x...]
    return QuantumEspressoTools.minimal_cif(title, pp1, addC.(positions))
end

## * ===============================================================

println("iX = $iX")

if iX > length(DA_LIST)
    exit()
else
    δa  = DA_LIST[iX]
    pos = C_frac_positions[δa]
    ttl = "shear_da_$(δa)"
    fn_rlxd = "$(WKSPC)/relaxed_structures/$(ttl).cif"
    make_sheared_cif(ttl, δa, pos)  ⇶  fn_rlxd
    if !isfile(fn_rlxd)  exit()  end

    @info "\n\n\n---------------------------------\n\nComputation for shear ( $(δa) )\n\n---------------------------------\n\n"
    wkspc0 = "$(WKSPC)/band_$(δa)/"
    wkspc0 = try_mkdir(wkspc0)

    band_w90_kgrid(
        wkspc0,                          #*  workspace0::String, 
        ttl,                             #*  common_title::String, 
        fn_rlxd,                         #*  cif0_fn::String,
        ("PSL_USPP_PBESOL_SR", 
         Dict("C"=>("n", "1.0.0"))),     #*  ps_mode
        (20,16,1,0,0,0),                 #*  kpoint_config::NTuple{6,Int64},
        ("tpiba_c", 
         [(-0.2, 0.0, 0.0) => 1, 
          ( 0.2, 0.0, 0.0) => 40,
          (-0.2, 0.4, 0.0) => 40,]
        ),
        (60.0,480.0),                    #*  cutoff::Tuple{Float64,Float64},
        32,                              #*  nbnd
        (10,15),
        ["C:pz",];
        cleanup     = false,
        PROG_PWX    = `mpiexec -np 48 pw.x -npool 6`,
        PROG_PW2W90 = `mpiexec -np 48 pw2wannier90.x`,
        PROG_W90    = `wannier90.x`,
        PROG_W90PP  = `wannier90.x -pp`
    )

end


