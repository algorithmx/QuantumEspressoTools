global const _electron_scf_seed_conf_ = Dict(
        :calc             => "scf", 
        "verbosity"       => "high",
        "disk_io"         => "medium",

        "nbnd"            => :default,
        "conv_thr"        => 1.0e-12,
        "mixing_beta"     => 0.8,
        "electron_maxstep"=> 500,
        "ecutwfc"         => 150.0, 
        "ecutrho"         => 600.0,
        :diag_choices     => ["david", "cg"],
        
        "occupations"     => "smearing",
        "smearing"        => "mp",
        "degauss"         => nothing,

        :pw_mode          => "energy",
        :cif              => nothing,
        :do_not_use_symmetry=> nothing,
)


global const _electron_pw2w90_seed_conf_ = Dict(
        "write_mmn"       => true,
        "write_amn"       => true,
        "write_unk"       => false,
        "write_uHu"       => false,
        "spin_component"  => "none",
)


global const _wannier90_seed_conf_ = Dict(
        "num_iter"        => 0,
        "dis_num_iter"    => 8000,
        "dis_mix_ratio"   => 0.9,

        "wannier_plot"    => false,
        "bands_plot"      => false,

        "postproc_setup"  => false,

        "write_bvec"      => true,
        "write_hr"        => true,
        "guiding_centres" => true,

        "spinors"         => false,
        "postproc_setup"  => false,
)


global const _projwfc_seed_conf_ = Dict(
        #"ngauss"          => 1,
        #"degauss"         => 0.02,
        "Emin"            => -10.0,
        "Emax"            =>  10.0,
        "DeltaE"          => 0.01,
) 


global const _cp_scf_seed_conf_ = Dict(
        "verbosity"       => "low",
        "restart_mode"    => "from_scratch",
        "nbnd"            => :default,
        "conv_thr"        => 1e-9,
        "etot_conv_thr"   => 1e-9,
        "ekin_conv_thr"   => 1e-6,
        "forc_conv_thr"   => 1e-4,
        "electron_maxstep"=> 2500,
        "ecutwfc"         => 60.0, 
        "ecutrho"         => 480.0,
        "tstress"         => true,
        "tprnfor"         => true,
        "degauss"         => 0.02,
        "press"           => 0.0,
        #! electron
        "electron_temperature" => "not_controlled",
        "electron_dynamics" => "sd",
        "emass"           => 400.0,
        "ekincw"          => 0.001,  
        "emass_cutoff"    => 2.5,
        "orthogonalization"    => "ortho",
        "ortho_eps"            => 1e-8,
        "ortho_max"            => 40,
        #! ion
        "ion_dynamics"    => "none",
        #! cell
        "cell_dynamics"   => "none",
        #
        "nr1"              => nothing, 
        "nr2"              => nothing, 
        "nr3"              => nothing,
        #
        "nr1s"              => nothing, 
        "nr2s"              => nothing, 
        "nr3s"              => nothing,
        #
        "nr1b"              => nothing, 
        "nr2b"              => nothing, 
        "nr3b"              => nothing,
)


global const _electron_relax_seed_conf_ = Dict(
        :calc             => "vc-relax", 
        "verbosity"       => "high",
        "restart_mode"    => "from_scratch",

        "nbnd"            => :default,
        "conv_thr"        => 1.0e-8,
        "mixing_beta"     => 0.8,
        "electron_maxstep"=> 500,
        "ecutwfc"         => 150.0, 
        "ecutrho"         => 600.0,

        "cell_dofree"     => "all",
        "press_conv_thr"  => 0.0001, 
        "cell_dynamics"   => "bfgs", 
        "bfgs_ndim"       => 3,
        "ion_dynamics"    => "bfgs", 
        "upscale"         => 100.0,
        "occupations"     => "smearing",

        :pw_mode          => "energy",
)


global const _cp_relax_seed_conf_ = Dict(
        "verbosity"       => "high",
        "restart_mode"    => "from_scratch",

        "nbnd"            => :default,
        "conv_thr"        => 1.0e-8,
        "etot_conv_thr"   => 1e-6,
        "forc_conv_thr"   => 1e-5,
        "electron_maxstep"=> 500,
        "ecutwfc"         => 60.0, 
        "ecutrho"         => 240.0,


        "degauss"         => 0.02,
        "press"           => 0.0,

        #
        "emass"           => nothing, #! electron
        "emass_cutoff"    => nothing, #! electron
        "wmass"           => nothing, #! cell
 
        #
        "ekincw"          => nothing,  #! electron
        "tempw"           => nothing,  #! ion
        "temph"           => nothing,  #! cell

        #
        "nr1"              => nothing, 
        "nr2"              => nothing, 
        "nr3"              => nothing,
        #
        "nr1s"              => nothing, 
        "nr2s"              => nothing, 
        "nr3s"              => nothing,
        #
        "nr1b"              => nothing, 
        "nr2b"              => nothing, 
        "nr3b"              => nothing,

        "cell_temperature"     => "not_controlled",
        "ion_temperature"      => "nose",
        "electron_temperature" => "nose",

        "cell_dynamics"        => nothing,
        "ion_dynamics"         => "verlet",
        "electron_dynamics"    => "verlet",

        "occupations"     => nothing,
        "smearing"        => nothing,

        "cell_dofree"     => "all",
        "press_conv_thr"  => 0.001, 
)


global const _md_seed_conf_ = Dict(
)


global const _phonon_seed_conf_ = Dict(
        "verbosity"     => "high",
        "mixing_beta"     => 0.8,
        "electron_maxstep"=> 500,

        "alpha_mix"     =>  0.5,
        "conv_thr_ph"   => 1e-12,
        "search_sym"    =>  false,
        "tr2_ph"        =>  1e-14,
)


global const _phonon_q2r_seed_conf_ = Dict(
        "acoustic_sum_rule"   => "simple"
)


global const _phonon_matdyn_seed_conf_ = Dict(
        "acoustic_sum_rule"    => "simple"
)


global const _phonon_plotband_seed_conf_ = Dict(
        :E_Fermi => 0.0,
        :dE      => 50.0,
)
