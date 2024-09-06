global const WKSPC = "/dg_hpc/CSG/lianyl/BANDR1"
global const iX = parse(Int,ARGS[1])
global const NX_LIST = collect(3:20)

@info Base.read(`hostname`,String)

sleep(30*rand())
using BSON

ENV["JULIA_PKG_SERVER"] = "https://mirrors.tuna.tsinghua.edu.cn/julia"
global const CSNS = "/csns_workspace/CSG/lianyl/softwares"
ENV["PSEUDO_ROOT"] = "$CSNS/pseudopotentials"

sleep(20*rand()+2*(iX%47))
include("$CSNS/QuantumEspressoTools/scripts/scripts.jl")
include("$CSNS/QuantumEspressoTools/scripts/energy_bands_on_kgrid.jl")
include("$CSNS/QuantumEspressoTools/scripts/band_wannier90.jl")

sleep(20*rand()+2*(iX%17))
Pkg.activate("$CSNS/cifio")
using cifio

sleep(20*rand()+2*(iX%17))
Pkg.activate("$CSNS/LatticeLab")
using LatticeLab

sleep(3*rand())
global const C_frac_positions = BSON.load("$WKSPC/C_frac_positions2.bson")[:pos]

## * ===============================================================

function cellp2a(cellp;digits=12)
    @inline rd(x) = round(x,digits=digits)
    (a_Ang, b_Ang, c_Ang, alpha, beta, gamma) = cellp[1:6]
    alphar, betar, gammar = map(x->Float64((π/big(180))*x), [alpha, beta, gamma])
    cell_px = a_Ang.*[1.0,         0.0,                                              0.0 ]
    cell_py = b_Ang.*[cos(gammar), sin(gammar),                                      0.0 ]
    cell_pz = c_Ang.*[cos(betar),  (cos(alphar)-cos(betar)*cos(gammar))/sin(gammar), 
                        sqrt(1.0 - cos(alphar)^2 - cos(betar)^2 - cos(gammar)^2 + 2*cos(alphar)*cos(betar)*cos(gammar))/sin(gammar)]
    return hcat(cell_px, cell_py, cell_pz) .|> rd
end

function cif2UC(cif::CIF, orbs)
    natom = length(cif.frac)
    mm  = [Symbol("$c$i") for (i,c) in enumerate(LatticeLab.Masses(Symbol.(first.(cif.frac))))]
    aa  = LatticeLab.Coordinates(cellp2a(cif.cellp))
    ff  = hcat([f[2:4] for f ∈ cif.frac]...)
    dd  = aa * ff |> LatticeLab.Coordinates  # Cartesian
    ξξ  = LatticeLab.Orbits.(orbs)
    UC = LatticeLab.UnitCell(3, natom, aa, dd, mm, ξξ)
    return UC
end

function UC2CIF(uc::LatticeLab.UnitCell, tit="")
    IT = 0
    cellp = cifio.basis_to_lattice_parameters(uc.a, IT)
    frac0  = inv(uc.a) * uc.δ
    frac = [[replace(string(uc.m[i]),r"\d+"=>""), frac0[:,i]...,] for i=1:size(frac0,2)]
    misc = nothing
    symm_ops = ["x,y,z"]
    return CIF(tit, cellp, frac, IT, symm_ops, misc)
end

function make_ring_along_x(cif::CIF, Nuc_x::Int; padding_x=2.0, padding_z=2.0)
    ξξ      = [[:pz],[:pz],[:pz],[:pz],[:pz],[:pz]] .|> LatticeLab.Orbits
    UC      = cif2UC(cif, ξξ)
    a1      = UC.a[1,1]
    R       = Nuc_x*a1/(2π)
    @inline x2x(x) = R*sin(x/R)
    @inline x2z(x) = R*cos(x/R)
    δδ      = hcat([UC.δ .+ k.*UC.a[:,1] for k=0:Nuc_x-1]...)
    δδ1      = copy(δδ[1,:])
    δδ[1,:] .= x2x.(δδ1)
    δδ[3,:] .= x2z.(δδ1)
    shift_δ  = [R+padding_x, 0.0, R+padding_z]
    natom    = UC.nsubl * Nuc_x
    mm       = vcat([UC.m for i=1:Nuc_x]...)
    aa       = hcat([2(R+padding_x),0,0], UC.a[:,2], [0,0,2(R+padding_z)])
    dd       = (δδ.+shift_δ) |> LatticeLab.Coordinates  # Cartesian
    ξξ       = vcat([UC.ξ for i=1:Nuc_x]...) .|> LatticeLab.Orbits
    return LatticeLab.UnitCell(3, natom, aa, dd, mm, ξξ) |> UC2CIF |> cifio.minimal_cif |> cifio.SPLTN
end

## * ===============================================================

strained_CELL_PARAMETERS(rx, ry) = [rx*3.745747057, ry*4.506578341, 20.0]

make_strained_relaxed_cif(title,cell,pos) =
    QuantumEspressoTools.minimal_cif(title, (cell...,90.0,90.0,90.0), [("C",k...) for k in pos])

## * ===============================================================


function find_window0(
    x,
    nbnd,
    large_window_band_id; 
    margin=5.0, 
    dE = 2.0
    )
    d = dict__pw_result(x)
    bnd = hcat(last.(values.(d["bands"])) ...)'  #
    E_fermi = d[:E_fermi]
    @assert size(bnd,2)==nbnd
    band_id_start    = first(large_window_band_id)
    band_id_end      = last( large_window_band_id)
    large_window_min = minimum(bnd[:,band_id_start])
    large_window_max = maximum(bnd[:,band_id_end  ])
    small_window_min = E_fermi - dE
    small_window_max = E_fermi + dE
    ret =  Dict("num_bands"   =>  nbnd,
                "num_wann"    =>  band_id_end - band_id_start + 1,
                "dis_win_min" =>  large_window_min-margin,
                "dis_win_max" =>  large_window_max+margin,
                "dis_froz_min"=>  small_window_min,
                "dis_froz_max"=>  small_window_max )
    return ret
end


find_window(x,nbnd,large_window_band_id) = 
    find_window0( x, nbnd, large_window_band_id; 
                  margin=5.0, dE = 2.0 )


## * ===============================================================

println("iX = $iX")

if iX > length(NX_LIST)
    exit()
else
    Nux     = NX_LIST[iX]
    w90scf_common(title, psmode) = Dict(
        "assume_isolated" => "2D",
        "title"           => title, 
        "prefix"          => title, 
        "disk_io"         => "medium",
        "conv_thr"        => Nux*1e-10,
        "etot_conv_thr"   => 1e-6,
        "forc_conv_thr"   => 1e-5,
        :pseudo_mode      => psmode,
        :watchdog_setting => w90scf_woof
    )
    cell    = strained_CELL_PARAMETERS(1.0,1.0)
    pos     = C_frac_positions[[1.0,1.0]]
    ttl1    = "strain00_$(Nux)"
    fn_rlxd = "$(WKSPC)/relaxed_structures/$(ttl1).cif"
    cifio.write_to_file(make_strained_relaxed_cif(ttl1,cell,pos), fn_rlxd)
    if !isfile(fn_rlxd)  exit()  end
    ttl     = "ring_N_$(Nux)"
    fn_ring = "$(WKSPC)/ring_structures/$(ttl).cif"
    cifio.write_to_file(make_ring_along_x(CIF(fn_rlxd; is_sort=false), Nux), fn_ring)
    if !isfile(fn_ring)  exit()  end

    @info "\n\n\n---------------------------------\n\nComputation for ring size ( $Nux )\n\n---------------------------------\n\n"

    wkspc0 = try_mkdir("$(WKSPC)/ring_$(Nux)/")

    band_w90_kgrid(
        wkspc0,                          #*  workspace0::String, 
        ttl,                             #*  common_title::String, 
        fn_ring,                         #*  cif0_fn::String,
        ("PSL_USPP_PBESOL_SR", 
         Dict("C"=>("n", "1.0.0"))),     #*  ps_mode
        (1,16,1,0,0,0),                   #*  kpoint_config::NTuple{6,Int64},
        ("crystal_b", 
         [(0.0, 0.0, 0.0) => 100, 
          (0.0, 1.0, 0.0) => 1,]),
        (60.0,480.0),                    #*  cutoff::Tuple{Float64,Float64},
        32*Nux,                          #*  nbnd
        (10*Nux,15*Nux),
        ["C:pz",];
        cleanup     = false,
        PROG_PWX    = `mpiexec -np 48 pw.x -npool 4`,
        PROG_PW2W90 = `mpiexec -np 48 pw2wannier90.x`,
        PROG_W90    = `wannier90.x`,
        PROG_W90PP  = `wannier90.x -pp`
    )
end

exit()
