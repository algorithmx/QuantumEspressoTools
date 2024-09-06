include("/public3/home/sc55341/softwares/QuantumEspressoTools/scripts/scripts.jl")

global const TEST_FOLDER = "/public3/home/sc55341/ScF3"

global const FD1 = "/public3/home/sc55341/ScF3/old"

##* ===================================================

function grab_alat(PS)
    ff = "$FD1/ScF3_$(PS)/pw.x.tmp.$(PS).k8x8.out"
    if isfile (ff)
        return pw_alat_bohr(readlines(ff)) * _BOHR_RADIUS_
    else
        return 7.1 * _BOHR_RADIUS_
    end
end


d2(PS) = Dict(
        "title"           => "ScF3_$PS", 
        "prefix"          => "ScF3_$PS", 
        :pseudo_mode      => PS,
        "mixing_beta"     => 0.8,
        "electron_maxstep"=> 500,
        "nbnd"            => 42,
        :pw_mode          => "energy",
        :kpoint_mode      => "automatic", 
        :kpoints          => (8,8,8,1,1,1),
        :watchdog_setting => (woof_per_x_min=1, max_watch=1440, max_silence_woofs=10, tail_length=20, quiet=false),
)

init_struct(PS,a0) = Dict(
        "positions"       => [("Sc", 0.0, 0.0, 0.0), ("F", 1/2, 0, 0), ("F", 0, 1/2, 0), ("F", 0, 0, 1/2)],
        :cif              => (a0,a0,a0, 90,90,90, 221),
        :do_not_use_symmetry => true,
)

dict_relax = Dict(
        "occupations"     => "smearing",
        "smearing"        => "mp",
        "degauss"         => 0.1,

        "cell_dofree"     => "all", 
        "press_conv_thr"  => 0.01, 
        "cell_dynamics"   => "bfgs", 
        "bfgs_ndim"       => 3,

        "ion_dynamics"    => "bfgs", 
        "upscale"         => 100.0,
)

dict_phx(PS) = Dict(
        :ph_mode    => :grid,
        "prefix"    => "ScF3_$(PS)", 
        "title"     => "ScF3_$(PS)",
        :kpoints=>(8,8,8,1,1,1),
        :qpoints    => (4,4,4),
)


function upd1(x)
    d = dict__pw_relax_result(x)
    return (d ↓ ["ibrav", :cif, :reciprocal_basis]) ⬱ Dict("ibrav"=>0,:do_not_use_symmetry=>true)
end


function upd2(x)
    d = dict__pw_result(x)
    return Dict(:reciprocal_basis => d[:reciprocal_basis])
end


function upd3(x,t="")
    freq_min = ph_min_frequency_cm_1(readlines(t))
    freq_max = ph_max_frequency_cm_1(readlines(t))
    return Dict(:freq_min_max_cm_1 => (min(0,1.05*freq_min), 1.05*freq_max))
end


dict_matdyn(PS) = Dict( "prefix" => "ScF3_$(PS)",
                        "title"  => "ScF3_$(PS)",
                        "flfrc"  => "ScF3_$(PS).fc",    #* IFC_filename
                        "flfrq"  => "ScF3_$(PS).frq",   #* frequency file
                        "flvec"  => "ScF3_$(PS).vec",   #* eigen vector 
                        "fleig"  => "ScF3_$(PS).eig",   #* eigen value
)

dict_q2r(PS) = Dict(    "prefix" => "ScF3_$(PS)",
                        "title"  => "ScF3_$(PS)",
                        "acoustic_sum_rule" => "simple",
                        "fildyn" => "ScF3_$(PS).dynmat.",  #* 
                        "flfrc"  => "ScF3_$(PS).fc",       #* IFC_filename
)


dict_plotband(PS) = Dict("prefix" => "ScF3_$(PS)",
                        :freq     => "ScF3_$(PS).frq",
                        :plot     => "ScF3_$(PS).phonon.band.plot",
                        :ps       => "ScF3_$(PS).phonon.band.ps"
)


instr(PS) = INSTRUCTION(
    Dict(
        `mpiexec -np 64 --map-by=slot pw.x`       => _electron_scf_seed_conf_  ⬱ d2(PS),
        `mpiexec -np 64 --map-by=slot ph.x`       => _phonon_seed_conf_  ⬱ d2(PS),
        `q2r.x`      => _phonon_q2r_seed_conf_,
        `matdyn.x`   => _phonon_matdyn_seed_conf_,
        `plotband.x` => _phonon_plotband_seed_conf_,
    ),
    [
        ("ScF3_relax",   `mpiexec -np 64 --map-by=slot pw.x`,  Dict(   
                                                :calc       => "vc-relax", 
                                                "conv_thr"  => 1e-8,
                                                "ecutwfc"   => 100.0,
                                                "ecutrho"   => 400.0,
                                                :updater    => upd1,
                                        ) ∪ init_struct(PS,grab_alat(PS)) ∪ dict_relax),
        ("ScF3_scf",     `mpiexec -np 64 --map-by=slot pw.x`, Dict(
                                                :calc       => "scf",
                                                "conv_thr"  => 1e-11,
                                                "ecutwfc"   => 150.0, 
                                                "ecutrho"   => 600.0,
                                                :updater    => upd2,
                                        )),
        ("ScF3_ph_grid", `mpiexec -np 64 --map-by=slot ph.x`, Dict(   
                                                "outdir"=>"$TEST_FOLDER/$PS/ScF3_scf",
                                        ) ∪ dict_phx(PS)),
        ("ScF3_ph_grid", `q2r.x`,       dict_q2r(PS) ),
        ("ScF3_ph_grid", `matdyn.x`,    Dict(   :qpath_rel=>[(0,0,0)=>30,
                                                             (0.5,0,0)=>30,
                                                             (0.5,0.5,0)=>30,
                                                             (0.5,0.5,0.5)=>45,
                                                             (0,0,0)=>1],
                                                :updater  => x->upd3(x,"ScF3_$(PS).ph.x.out")
                                        ) ∪ dict_matdyn(PS) ),
        ("ScF3_ph_grid", `plotband.x`,  dict_plotband(PS) ),
    ]
)

#* ===================================================

const PS_TO_TEST = ["GHH_PBE", "GHH_PZ", "MT_PBE", "MT_PW_LDA"]
const PS_TO_TEST = ["GBRV_LDA", "GBRV_PBE", "GBRV_PBESOL", "GHH_PBE", "GHH_PZ", "MT_PBE", "MT_PW_LDA", "SG15", "SSSP_Efficiency", "SSSP_Precision"]

##* ===================================================

i = parse(Int,ENV["SLURM_ARRAY_TASK_ID"])

if i > length(PS_TO_TEST)
    @error "i = $i > length(PS_TO_TEST) = $(length(PS_TO_TEST)). Exciting ..."
    exit()
end

PSi = PS_TO_TEST[i]

R = execute_serial(instr(PSi), WORKSPACE="$TEST_FOLDER/$PSi", from_scratch=false)

@save "ScF3_ALL_RESULTS_$(PSi).jld2" R

##* ===================================================
