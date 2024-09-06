include("/home/dabajabaza/jianguoyun/Workspace/QuantumEspressoTools/scripts/scripts.jl")

global const TEST_FOLDER = "/data/ld1"

##* ===================================================

ld1_config = Dict(
        "title"     => "Sc",
        "prefix"    => "Sc",
        "zed"       => 21.0,
        #"config"    => "[Ar] 4s2 4p0 3d1",
        "config"    => "[Ne] 3s2 3p6 4s2 4p0 3d1",
        "iswitch"   => 3,
        "nlcc"      => true,
        "rel"       => 1,
        "lsd"       => 1,
        "dft"       => "PBE",
        "pseudotype"=> 3,
        "tm"        => true,
        "use_xsd"   => false,
        "file_pseudopw" => "Sc.pbe-spn-kjpaw_psl.1.0.0.Ne.UPF",
        "lloc"      => -1,
        "rcloc"     => 1.2,
        "which_augfun" => "PSQ",
        "rmatch_augfun_nc" => true,
        "beta"      => 0.2,
        "new_core_ps"  => true,
        "rcore"     => 0.8,
        :ld1_mode   => ["All electron", "PP generation"],
        :watchdog_setting => (woof_per_x_min=1/10, max_watch=1440, max_silence_woofs=10, tail_length=20, quiet=false),
)

c = INSTRUCTION(
    Dict(
        `mpiexec -np 4 ld1.x -npool 2`  => _electron_scf_seed_conf_  ⬱ ld1_config,
    ),
    [   ("ld1",   `mpiexec -np 4 ld1.x -npool 2`,  
                    Dict(  :electronic_configuration => [
                            ],
                           :wavefunctions_spec => [ 
                                #: nls(1) nns(1) lls(1) ocs(1) ener(1) rcut(1) rcutus(1) 
                                #: {String,Int,Int,Float64,Float64,Float64,Float64}
                                ["3S",  1,  0,  2.00,  0.00,  0.90,  1.30,   0.0],
                                ["4S",  2,  0,  2.00,  0.00,  0.90,  1.30,  0.0],
                                ["3P",  2,  1,  6.00,  0.00,  0.90,  1.45,  0.0],
                                ["4P",  3,  1,  0.00,  7.00,  0.90,  1.45,  0.0],
                                ["3D",  3,  2,  1.00,  0.00,  1.00,  1.50,  0.0],
                                ["3D",  3,  2,  0.00,  6.00,  1.00,  1.50,  0.0],
                            ]
                    )
        ),
    ]
) ;

##* ===================================================

R = execute_serial(c, WORKSPACE="$TEST_FOLDER")

