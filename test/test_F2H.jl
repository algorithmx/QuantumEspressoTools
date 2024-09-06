include("/home/dabajabaza/jianguoyun/Workspace/QuantumEspressoTools/scripts/scripts.jl")

global const TEST_FOLDER = "/data/test/F2H"

global const PSi = "SG15_FR"

##* ===================================================

d2(PS) = Dict(
        "title"           => "F2H_$PS", 
        "prefix"          => "F2H_$PS", 
        :pseudo_mode      => PS,
        "mixing_beta"     => 0.8,
        "verbosity"       => "high",
        :calc             => "scf",
        "conv_thr"        => 1e-12,
        "ecutwfc"         => 200.0, 
        "ecutrho"         => 800.0,
        "electron_maxstep"=> 1500,
        "nbnd"            => 18,
        :pw_mode          => "energy",
        :kpoint_mode      => "automatic", 
        :kpoints          => (1,1,1,0,0,0),
)

init_struct(d,a0=20) = Dict(
        "positions"       => [("F", 1/2, 1/2, 1/2-d/a0), ("F", 1/2, 1/2, 1/2+d/a0), ("H", 1/2, 1/2, 1/2)],
        :cif              => (a0, a0, a0,  90, 90, 90, 221),
        :do_not_use_symmetry => true,
)

instr(PS,d,a0=20) = INSTRUCTION(
    Dict(
        `mpiexec -np 4 pwdry.x -npool 2`    => _electron_scf_seed_conf_  ⬱ d2(PS),
        `mpiexec -np 4 pw.x -npool 2`       => _electron_scf_seed_conf_  ⬱ d2(PS),
    ),
    [
        ("F2H_scf",  `mpiexec -np 4 pwdry.x -npool 2`, 
                        Dict( :watchdog_setting  => default_intercepter_setting,
                        ) ∪ init_struct(d,a0) ),
        ("F2H_scf",  `mpiexec -np 4 pw.x -npool 2`,
                        Dict( :watchdog_setting  => default_watchdog_setting,
                        ) ∪ init_struct(d,a0) ),
    ]
)

##* ===================================================

R = execute_serial(instr(PSi, 2.1), WORKSPACE="$TEST_FOLDER/$PSi")

##* ===================================================
