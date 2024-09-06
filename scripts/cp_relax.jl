global const cp_relax_woof = (woof_per_x_min=1/6, max_watch=20000, max_silence_woofs=10, tail_length=20, quiet=false)

cp_relax_common(title, psmode, dt, TL, ML, nr, kBar, dg) = Dict(
        :calc             => "vc-relax", 
        "title"           => title, 
        "prefix"           => title, 
        :pseudo_mode      => psmode, 
        "disk_io"         => "default", 
        "verbosity"       => "low", 
        "conv_thr"        => 1e-9,

        "dt"              => dt,

        "degauss"         => dg,
        "press"           => Float64(kBar),

        #
        "electron_temperature" => "nose",
        "electron_dynamics"    => "verlet",
        "ekincw"               => TL[1],  #! electron
        "emass"                => ML[1], #! electron
        "emass_cutoff"         => ML[2], #! electron

        "ion_dynamics"         => "verlet",
        "ion_temperature"      => "nose",
        "tempw"                => TL[2],  #! ion
        
        "cell_dynamics"        => "pr",
        "cell_temperature"     => "nose",
        "temph"                => TL[3],  #! cell
        "wmass"                => ML[3], #! cell

        #
        "nr1"                  => nr[1], 
        "nr2"                  => nr[2], 
        "nr3"                  => nr[3],
        #
        "nr1s"                 => nr[4], 
        "nr2s"                 => nr[5], 
        "nr3s"                 => nr[6],
        #
        "nr1b"                 => nr[7], 
        "nr2b"                 => nr[8], 
        "nr3b"                 => nr[9],

        :watchdog_setting      => cp_relax_woof
)


#TODO test
function cp_relax_update_random_kick(x,delta=0.01,elems=[])
    d = dict__cp_relax_result(x)
    d = random_kick_positions(d, delta, elems)
    d = random_kick_lattp(d, delta)
    d1 = (d ↓ ["ibrav", :cif, :reciprocal_basis]) ⬱ Dict("ibrav"=>0,:do_not_use_symmetry=>true)
    return d1
end


#: ==================================================================


function cp_vc_relax(
    cell_dofree_list::Vector{String},
    workspace0::String,
    common_title::String,
    cif0_fn::String,
    ps_modes::Vector,
    dt_list::Vector,
    temp_list::Vector,
    m_list::Vector,
    nrXb_list::Vector,
    cutoff_list::Vector{Tuple{Float64,Float64}},
    kBar_list::Vector{Float64},
    degauss_list::Vector{Float64};
    # 
    relax_iter_updater = x->cp_relax_update_random_kick(x,0.02,[]),
    PROG_CPX = `mpiexec -np 64 --map-by=slot cp.x -npool 8`,
    atoms = []
    )

    workspace = try_mkdir(workspace0, ".", "cp_vc_relax()")
    cd(workspace)
    WORKSPACE = pwd()
    ENV["LOGGER_FILE_FULL_PATH"] = "$WORKSPACE/cp_vc_relax.log"

    #> ---------------------------------
    #> customized instruction
    function instructor(init_struct, FD, psmode, dt, TL, ML, nrXb, cutoffs, kBar, dg, atms=atoms)
        wcut, ecut = cutoffs
        cp_relax_comm  = cp_relax_common(common_title, psmode, dt, TL, ML, nrXb, kBar, dg)
        INSTR = INSTRUCTION(
            Dict(PROG_CPX => (_cp_relax_seed_conf_  ⬱ cp_relax_comm)),
            [   (   "$(common_title)_cp_vc_relax_iter_$(itern)", 
                    PROG_CPX,
                    (Dict(  :calc             => "vc-relax", 
                            :cutoff_scales    => [1.0, ],     #! fixed
                            "ecutwfc"   => wcut,
                            "ecutrho"   => ecut,
                            :updater    => relax_iter_updater,
                    ) ∪ ((itern==1) ? init_struct : Dict()) ) ↑ ("cell_dofree"=>dofree) ) 
                for (itern, dofree) ∈ enumerate(cell_dofree_list)
            ]
        )
        return INSTR
    end

    #> ---------------------------------
    #> folders and settings
    cp_nt(ps,dt,T,M,nrXb,cc,kBar,dg) = join([
        "ps___$(pseudo_mode_name(ps))___dt_$(dt)", 
        "T_$(T[1]),$(T[2]),$(T[3])",
        "M_$(M[1]),$(M[2]),$(M[3])",
        "nrx_$(nrXb[1]),$(nrXb[2]),$(nrXb[3])", 
        "nrxs_$(nrXb[4]),$(nrXb[5]),$(nrXb[6])",
        "nrxb_$(nrXb[7]),$(nrXb[8]),$(nrXb[9])",
        "cutoff_$(cc[1]),$(cc[2])_kBar_$(kBar)_dg_$(dg)"
    ], "_")

    #> ---------------------------------
    #> serial excution
    res = grid_test(    WORKSPACE, 
                        common_title,
                        cif0_fn, 
                        instructor,
                        cp_nt,
                        [ps_modes, dt_list, temp_list, m_list, nrXb_list, cutoff_list, kBar_list, degauss_list]
    )
    return res
end


function cp_vc_relax(
    workspace0::String,
    common_title::String,
    cif0_fn::String,
    ps_modes::Vector,
    dt_list::Vector,
    temp_list::Vector,
    m_list::Vector,
    nrXb_list::Vector,
    cutoff_list::Vector{Tuple{Float64,Float64}},
    kBar_list::Vector{Float64},
    degauss_list::Vector{Float64};
    # 
    relax_iter_updater = x->cp_relax_update_random_kick(x,0.02,[]),
    PROG_CPX = `mpiexec -np 64 --map-by=slot cp.x -npool 8`,
    atoms = []
    )
    cp_vc_relax(
        String["all",],
        workspace0,
        common_title,
        cif0_fn,
        ps_modes,
        dt_list,
        temp_list,
        m_list,
        nrXb_list,
        cutoff_list,
        kBar_list,
        degauss_list;
        # 
        relax_iter_updater = relax_iter_updater,
        PROG_CPX = PROG_CPX,
        atoms = atoms
    )
end


function cp_vc_relax(
    dofree_N::Tuple{String,Int},
    workspace0::String,
    common_title::String,
    cif0_fn::String,
    ps_modes::Vector,
    dt_list::Vector,
    temp_list::Vector,
    m_list::Vector,
    nrXb_list::Vector,
    cutoff_list::Vector{Tuple{Float64,Float64}},
    kBar_list::Vector{Float64},
    degauss_list::Vector{Float64};
    # 
    relax_iter_updater = x->cp_relax_update_random_kick(x,0.02,[]),
    PROG_CPX = `mpiexec -np 64 --map-by=slot cp.x -npool 8`,
    atoms = []
    )
    cp_vc_relax(
        String[dofree_N[1] for i=1:dofree_N[2]],
        workspace0,
        common_title,
        cif0_fn,
        ps_modes,
        dt_list,
        temp_list,
        m_list,
        nrXb_list,
        cutoff_list,
        kBar_list,
        degauss_list;
        # 
        relax_iter_updater = relax_iter_updater,
        PROG_CPX = PROG_CPX,
        atoms = atoms
    )
end
