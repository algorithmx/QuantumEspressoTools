include("/home/dabajabaza/jianguoyun/Workspace/QuantumEspressoTools/scripts/scripts.jl")

global const TEST_FOLDER = "/data/test"

global const PSi = "GBRV_LDA"

ff = "/home/dabajabaza/Downloads/ScF3_PSL_PAW_PBE_SR.pw.x.out"
pw_alat_bohr(readlines(ff))

##* ===================================================

function grab_alat(PS)
    return 3.796
end

d2(PS) = Dict(
        "title"           => "ReO3_$PS", 
        "prefix"          => "ReO3_$PS", 
        :pseudo_mode      => PS,
        "mixing_beta"     => 0.8,
        "verbosity"       => "high",
        :calc             => "scf",
        "conv_thr"        => 1e-11,
        "ecutwfc"         => 150.0, 
        "ecutrho"         => 600.0,
        "electron_maxstep"=> 500,
        "nbnd"            => 42,
        :pw_mode          => "energy",
        :kpoint_mode      => "automatic", 
        :kpoints          => (4,4,4,1,1,1),
)

init_struct(a0) = Dict(
        "positions"       => [("Re", 0.0, 0.0, 0.0), ("O", 1/2, 0, 0), ("O", 0, 1/2, 0), ("O", 0, 0, 1/2)],
        :cif              => (a0,a0,a0, 90,90,90, 221),
        :do_not_use_symmetry => true,
)

instr(PS) = INSTRUCTION(
    Dict(
        `mpiexec -np 4 pw.x -npool 2`       => _electron_scf_seed_conf_  ⬱ d2(PS),
        `mpiexec -np 4 pwdry.x -npool 2`    => _electron_scf_seed_conf_  ⬱ d2(PS),
    ),
    [
        ("ReO3_scf",  `mpiexec -np 4 pw.x -npool 2`, Dict( 
                                        :watchdog_setting  => default_watchdog_setting,
                              ) ∪ init_struct(grab_alat(PS))),
        ("ReO3_scf",  `mpiexec -np 4 pw.x -npool 2`, Dict( 
                                        :watchdog_setting  => default_watchdog_setting,
                              ) ∪ init_struct(grab_alat(PS))),
        ("ReO3_scf",  `mpiexec -np 4 pwdry.x -npool 2`, Dict( 
                                        "nstep"           => 0,
                                        "electron_maxstep"=> 1,
                                        :watchdog_setting  => default_intercepter_setting,
                              ) ∪ init_struct(grab_alat(PS))),
    ]
)

##* ===================================================

R = execute_serial(instr(PSi), WORKSPACE="$TEST_FOLDER/$PSi")

##* ===================================================
