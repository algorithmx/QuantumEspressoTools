include("/public3/home/sc55341/softwares/QuantumEspressoTools/scripts/scripts.jl")

global const TEST_FOLDER = "/public3/home/sc55341/ScF3_pressure"

global const PS_TO_TEST = [ 
"PSL_USPP_LDA_SR", "PSL_USPP_PBE_SR", "PSL_USPP_PBESOL_SR", 
"PSL_PAW_PZ_SR"  , "PSL_PAW_PBE_SR" , "PSL_PAW_PBESOL_SR" , 
]

i = parse(Int,ENV["SLURM_ARRAY_TASK_ID"])

if i > length(PS_TO_TEST)
    @error "i = $i > length(PS_TO_TEST) = $(length(PS_TO_TEST)). Exciting ..."
    exit()
end

PSi = PS_TO_TEST[i]

##* ===================================================

function grab_alat(PS)
    ff = "$TEST_FOLDER/$(PS)_0bar/ScF3_scf/ScF3_$(PS).pw.x.out"
    if isfile(ff) && verify_QE_result(readlines(ff), `pw.x`)
        return pw_alat_bohr(readlines(ff)) * _BOHR_RADIUS_
    else
        return 3.8
    end
end

d2(PS) = Dict(
        "title"           => "ScF3_$PS", 
        "prefix"          => "ScF3_$PS", 
        "verbosity"       => "high",
        "tstress"         => true,
        "tprnfor"         => true,
        :pseudo_mode      => PS,
        "mixing_beta"     => 0.8,
        "electron_maxstep"=> 500,
        "nbnd"            => 42,
        :pw_mode          => "energy",
        :kpoint_mode      => "automatic", 
        :kpoints          => (8,8,8,1,1,1),
        :watchdog_setting => (woof_per_x_min=1, max_watch=1440, max_silence_woofs=10, tail_length=20, quiet=false),
)

init_struct(a0) = Dict(
        "positions"       => [("Sc", 0.0, 0.0, 0.0), ("F", 1/2, 0, 0), ("F", 0, 1/2, 0), ("F", 0, 0, 1/2)],
        :cif              => (a0,a0,a0, 90,90,90, 221),
        :do_not_use_symmetry => true,
)

dict_relax(kBar) = Dict(
        "etot_conv_thr"   => 1e-7,
        "forc_conv_thr"   => 1e-6,
        "occupations"     => "smearing",
        "smearing"        => "mp",
        "degauss"         => 0.02,

        "cell_dofree"     => "all", 
        "press"           => Float64(kBar),
        "press_conv_thr"  => 0.0001, 
        "cell_dynamics"   => "bfgs", 
        "bfgs_ndim"       => 3,

        "ion_dynamics"    => "bfgs", 
        "upscale"         => 100.0,
)

dict_phx(PS) = Dict(
        :ph_mode    => :grid,
        "prefix"    => "ScF3_$(PS)", 
        "title"     => "ScF3_$(PS)",
        :kpoints    => (8,8,8,1,1,1),
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


instr(PS,kBar) = INSTRUCTION(
    Dict(
        `mpiexec -np 64 --map-by=slot pwdry.x -npool 8`       => _electron_scf_seed_conf_  ⬱ d2(PS),
        `mpiexec -np 64 --map-by=slot pw.x -npool 8`       => _electron_scf_seed_conf_  ⬱ d2(PS),
        `mpiexec -np 64 --map-by=slot ph.x -npool 8`       => _phonon_seed_conf_  ⬱ d2(PS),
        `q2r.x`      => _phonon_q2r_seed_conf_,
        `matdyn.x`   => _phonon_matdyn_seed_conf_,
        `plotband.x` => _phonon_plotband_seed_conf_,
    ),
    [
        ("ScF3_scf_dryrun",   `mpiexec -np 64 --map-by=slot pwdry.x -npool 8`,  
                                Dict(   :calc       => "scf", 
                                        "nstep"     => 0,
                                        "electron_maxstep" => 0,
                                        "disk_io"   => "none",
                                        "conv_thr"  => 1e-11,
                                        "ecutwfc"   => 150.0, 
                                        "ecutrho"   => 600.0,
                                        :from_scratch => true,
                                ) ∪ init_struct(grab_alat(PS)) ),
        ("ScF3_relax",   `mpiexec -np 64 --map-by=slot pw.x -npool 8`, 
                                Dict(   :calc       => "vc-relax", 
                                        "disk_io"   => "none",
                                        "conv_thr"  => 1e-8,
                                        "ecutwfc"   => 150.0, 
                                        "ecutrho"   => 600.0,
                                        :updater    => upd1,
                                        :from_scratch => true,
                                ) ∪ init_struct(grab_alat(PS)) ∪ dict_relax(kBar) ),
        ("ScF3_scf",     `mpiexec -np 64 --map-by=slot pw.x -npool 8`, Dict(
                                                :calc       => "scf",
                                                "disk_io"   => "high",
                                                "conv_thr"  => 1e-11,
                                                "ecutwfc"   => 150.0, 
                                                "ecutrho"   => 600.0,
                                                "wfcdir"    => "$TEST_FOLDER/$(PS)_$(Int(1000kBar))bar/ScF3_scf/wfc",
                                                :updater    => upd2,
                                        )),
        ("ScF3_ph_grid", `mpiexec -np 64 --map-by=slot ph.x -npool 8`, Dict( 
                                                "outdir"=>"$TEST_FOLDER/$(PS)_$(Int(1000kBar))bar/ScF3_scf",
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


##* ===================================================

R = [   execute_serial(instr(PSi,(1//20)*ikBar), WORKSPACE="$TEST_FOLDER/$(PSi)_$(Int((1000//20)*ikBar))bar", from_scratch=false) 
        for ikBar ∈ [0, 1, 3, 5, 9, 21]   ]

@save "ScF3_ALL_RESULTS_$(PSi)_pressure.jld2" R

##* ===================================================
