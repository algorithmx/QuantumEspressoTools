# pwx_input.jl
# ref : https://www.quantum-espresso.org/Doc/INPUT_PW.html

PW_Float_Int = Union{Float64, Int64}

#* ============================================================


function CONTROL(config::Dict)
    #> two types of pseudo-mode
    @inline psn(ps) = (ps isa String) ? ps : ps[1]
    NL = Namelist("CONTROL",Dict(),Dict(),Dict())
    NL.req = Dict(
        "calculation" => config[:calc],
        "title"       => config["title"] * "_" * config[:calc],
        "prefix"      => config["prefix"],
        "pseudo_dir"  => PSEUDO_PATHS[psn(config[:pseudo_mode])],
    )
    NL.default = Dict(
        "outdir"        => "./",
        "restart_mode"  => "from_scratch",
        "verbosity"     => "default",
        "disk_io"       => "default",
        "wf_collect"    => false,
        "nstep"         => (config[:calc]∈["scf","nscf","bands"] ? 1 : 250),
        "max_seconds"   => 1e6,
        #"wfcdir"        => rstrip(config["outdir"],'/') * "/wfc",
        #> I had experimented the option 'wfcdir'. 
        #> specifying this folder with disk_io=high
        #> just got a copy of .wfc1 files inside that folder. 
        #> The other copy is already saved in the workspace (outdir).

        "tstress"       => true,
        "tprnfor"       => true,
        "etot_conv_thr" => 1e-5,
        "forc_conv_thr" => 1e-4,

        #> for cp.x only
        # time step for molecular dynamics, in Hartree atomic units
        # (1 a.u.=2.4189 * 10^-17 s) : beware, PW code use
        # Rydberg atomic units, twice that much!!!)
        "dt"            => 1.0,
        "isave"         => 5,
        "iprint"        => 5,
        "saverho"       => false,
        "ekin_conv_thr" => 1e-6,
        "ndw"           => 50,
        "ndr"           => 50,
        "tabps"         => false,
    )
    NL.opt = NL.default ← config
    return NL
end


function ELECTRONS(config::Dict)
    NL = Namelist("ELECTRONS",Dict(),Dict(),Dict())
    NL.req = Dict(
    )
    NL.default = Dict(
        "electron_maxstep"  => 500,
        "conv_thr"          => 1e-10,
        "startingwfc"       => "atomic+random",     #* set to 'file' when MC 
        "ampre"             => 0.0,
        "startingpot"       => "atomic",            #* set to 'file' when MC 

        # pw.x
        "diagonalization"   => "david",
        "mixing_mode"       => "plain",
        "diago_thr_init"    => 1e-4,
        "diago_david_ndim"  => 2,
        "diago_cg_maxiter"  => 200,
        "mixing_beta"       => 0.7,
        "mixing_ndim"       => 8,

        #> for cp.x only
        "electron_dynamics"    => "none",
        "emass"                => 400.0,
        "emass_cutoff"         => 2.5,
        "orthogonalization"    => "ortho",
        "ortho_eps"            => 1e-8,
        "ortho_max"            => 30,
        "electron_damping"     => 0.1,
        "electron_velocities"  => "default",
        "electron_temperature" => "not_controlled",
        "ekincw"               => 1e-3, 
        "fnosee"               => 1.0,
        # cg
        "niter_cg_restart"     => 20,
        "tcg"                  => false,
        "maxiter"              => 100,
        "passop"               => 0.3,
        "n_inner"              => 2,
        "ninter_cold_restart"  => 1,
        "lambda_cold"          => 0.03,
        "grease"               => 1.0,
        "ampre"                => 0.0,
    )
    NL.opt = NL.default ← config
    return NL
end


function natom(positions)
    @inline wyck_in_pos(l) = typeof(l[2])<:AbstractString
    sg  = all(wyck_in_pos.(positions))
    sg1 = any(wyck_in_pos.(positions))
    @assert (sg && sg1) || (!sg && !sg1)
    if sg
        return sum(strtonum(replace(l[2],r"[a-z]"=>"")) for l in positions)
    else
        return length(positions)
    end
end


function nband(positions)
    @inline wyck_in_pos(l) = typeof(l[2])<:AbstractString
    sg  = all(wyck_in_pos.(positions))
    sg1 = any(wyck_in_pos.(positions))
    @assert (sg && sg1) || (!sg && !sg1)
    if sg
        _nbnds = sum(N_SHELL_ELECTRONS[l[1]]*strtonum(replace(l[2],r"[a-z]"=>""))  for l in positions)
        return _nbnds
    else
        _nbnds = sum(N_SHELL_ELECTRONS[a]  for a in first.(positions))
        return _nbnds
    end
end


#IMPORTANT  priority of the keys in config : celldm(1) > :cell_parameters > :cif
function SYSTEM(config::Dict)

    NL = Namelist("SYSTEM",Dict(),Dict(),Dict())

    #* ---------------- PREPARE -----------------
    uniqueb = get(config,"uniqueb",false)  #! must be specified here because we use it later

    atoms_al= first.(config["positions"])
    atoms   = unique(atoms_al)

    #: -------------------------
    #: ad-hoc method to update a few required keys 
    ntyp    = length(atoms)
    nat     = natom(config["positions"])
    nbnd    = num_KS_states( config[:pseudo_mode], atoms_al) #TODO wyckoff positions ?
    wcut0   = 0.0 # max(1.2*make_pseudo_wcut(config[:pseudo_mode], atoms), 50.0) # "ecutwfc"
    ecut0   = 0.0 # max(1.2*make_pseudo_ecut(config[:pseudo_mode], atoms), 200.0) # "ecutrho"
    ecut    = 0.0 # max(ecut0, 4wcut0)
    dic1    = Dict( "ecutwfc"        => 0.25*ecut,
                    "ecutrho"        => ecut,
                    "nbnd"           => nbnd, 
              ) ← config   #! ad-hoc method to update a few required keys  
    #: -------------------------

    #: -------------------------
    # filter the space-group related settings via struct  `UNITCELL` 
    cell_inp= config_to_cell_to_QE_pwx_input(config)
    k2del   = (get(config,:do_not_use_symmetry,false) ? "space_group" : nothing)
    cell_cnf= ((cell_inp ↓ k2del) ↓ :cell_parameters)
    #: -------------------------

    NL.req = Dict(
        "nat"   => nat,
        "ntyp"  => ntyp,
    ) ∪ dic1 ∪ cell_cnf
    NL.default = Dict(
        "tot_charge"        => 0.0,
        "tot_magnetization" => -1,
        "assume_isolated"   => "none",
        "occupations"       => "smearing",
        "smearing"          => "mp",
        "degauss"           => 0.01,

        "lspinorb"          => false,
        "noncolin"          => false,
        "noinv"             => false,
        "nosym"             => false,
        "origin_choice"     => 1,
        "rhombohedral"      => true,
        
        "nr1"               => 60, 
        "nr2"               => 60, 
        "nr3"               => 60,
        "nr1s"              => 60, 
        "nr2s"              => 60, 
        "nr3s"              => 60,

        "lda_plus_u"        => false,
       ["Hubbard_U($i)"     => 0.0    for i=1:ntyp]...,

        "vdw_corr"          => "none",
        "xdm_a1"            => 1.2153,
        "xdm_a2"            => 2.3704,
        
        #> cp.x only
        "nr1b"              => 40, 
        "nr2b"              => 40, 
        "nr3b"              => 40,
        "nspin"             => 1,
        "london_s6"         => 0.75,
        "london_rcut"       => 200.0,
        "ts_vdw_econv_thr"  => 1e-6,  #! too safe ?
        "ts_vdw_isolated"   => false,
    )
    NL.opt = (NL.default ← config)
    return NL
end


function IONS(config::Dict)
    NL = Namelist("IONS",Dict(),Dict(),Dict())
    if "ion_dynamics" ∉ keys(config)
        return NL
    else
        ion_dynamics = lowercase(config["ion_dynamics"])
        if  ion_dynamics=="bfgs"
            #* BFGS specific
            NL.req     = Dict(
                "ion_dynamics"          => ion_dynamics,
            )
            NL.default = Dict(
                "bfgs_ndim"             => 3, # Number of old forces and displacements vectors used in the PULAY mixing
                "upscale"               => 100.0,  # Max reduction factor for conv_thr during structural optimization
                "trust_radius_max"      => 0.8,
                "trust_radius_min"      => 0.001,
                "trust_radius_ini"      => 0.5,
                "pot_extrapolation"     => "none",
                "wfc_extrapolation"     => "none",
            )
        elseif config[:calc] ∈ ["md","vc-md"]
            #* variables used for molecular dynamics
            NL.req     = Dict(
                "ion_dynamics"          => ion_dynamics,
            )
            NL.default = Dict(
                "ion_temperature"       => "not_controlled",
                "nraise"                => 1,
            )
        else
            NL.req = Dict("ion_dynamics" => ion_dynamics)
        end
        NL.opt = NL.default ← config
        return NL
    end
end


function CELL(config::Dict)

    NL = Namelist("CELL",Dict(),Dict(),Dict())

    if (config[:calc]!="vc-md" && config[:calc]!="vc-relax") return NL end
    cell_dynamics = (("cell_dynamics" ∈ keys(config)) ? lowercase(config["cell_dynamics"]) : (config[:calc]=="vc-relax" ? "bfgs" : "none"))

    if config[:calc]=="vc-relax"
        @assert cell_dynamics ∈ ["none", "sd", "damp-pr", "damp-w", "bfgs"]  "$cell_dynamics not allowed for vc-relax."
    elseif config[:calc]=="vc-md"
        @assert cell_dynamics ∈ ["none", "pr", "w"]  "$cell_dynamics not allowed for vc-md."
    else
        #@info "CELL card not necessary for calculation $(config[:calc])."
        return NL
    end

    NL.default = Dict(
        "cell_dofree" => "all",
        "press" => 0.0, # [KBar]
        "press_conv_thr" => 0.1, 
    ) ∪ Dict("cell_dynamics"=>cell_dynamics)

    NL.opt = NL.default ← config
    return NL

end



function K_POINTS(config::Dict)

    CD = Card{PW_Float_Int}("K_POINTS", config[:kpoint_mode], Vector{PW_Float_Int}[])

    if config[:kpoint_mode] == "gamma"
        nothing
    elseif  config[:kpoint_mode] == "automatic"
        @assert (config[:kpoints] isa NTuple{6,Int64})
        CD.contents = [Vector{PW_Float_Int}([config[:kpoints]...])]
    elseif config[:kpoint_mode] ∈ ["crystal",]
        #* example of "crystal"
        # [[0.0,0.5,0.0,0.125], [0.0,0.5,0.5,0.125], ...]
        @assert (config[:kpoints] isa Vector) && eltype(config[:kpoints][1])<:Real
        CD.contents = [ Vector{PW_Float_Int}([length(config[:kpoints]),]),
                        [Vector{PW_Float_Int}(p) for p in config[:kpoints]]... ]
    elseif config[:kpoint_mode] ∈ ["tpiba_b","crystal_b"]
        #* example of "crystal_b"
        # [(0.0,0.0,0.0)=>5,(0.5,0,0)=>5,(0.5,0.5,0)=>5,(0,0,0)=>1]
        #* example of "tpiba_b"
        # IN UNITS OF 2π/alat
        # [(0.0,0.0,0.0)=>5,(0.5,0,0)=>5,(0.5,0.23,0)=>5,(0,0,0)=>1]
        @assert (config[:kpoints] isa Vector) && eltype(config[:kpoints])<:Pair
        CD.contents = [ Vector{PW_Float_Int}([length(config[:kpoints]),]),
                        [Vector{PW_Float_Int}([p...,n]) 
                         for (p,n) in config[:kpoints]]... ]
    elseif config[:kpoint_mode] == "tpiba_c"
        @assert length(config[:kpoints])==3
        # in unit of 2 pi/a.
        # k_0_x, k_0_y, k_0_z, n_0  ## origin
        # k_1_x, k_1_y, k_1_z, n_1  ## base vec x
        # k_2_x, k_2_y, k_2_z, n_2  ## base vec y
        # range : k0 + α (k1-k0) + β (k2-k0)
        CD.contents = [ Vector{PW_Float_Int}([3,]),
                        [Vector{PW_Float_Int}([p...,n]) 
                         for (p,n) in config[:kpoints]]... ]
    end

    return CD

end


#* =========================================================

#IMPORTANT  priority of the keys in config : celldm(1) > :cell_parameters > :cif

#> pw.x , cp.x
function CELL_PARAMETERS(config::Dict)
    #: cross-ref : dict__pw_result(pw_out_fn), qrel_to_qabs()
    CD = Card{Float64}("CELL_PARAMETERS", "bohr", Vector{Float64}[]) #TODO bohr -> angstrom
    if :cell_parameters in keys(config)
        CD.contents = config[:cell_parameters]./_BOHR_RADIUS_  # config[:cell_parameters] given in Angstrom
    else
        # filter the space-group related settings via struct  `UNITCELL` 
        cell_cnf = config_to_cell_to_QE_pwx_input(config)
        if :cell_parameters in keys(cell_cnf)
            CD.contents = cell_cnf[:cell_parameters]./_BOHR_RADIUS_  # config[:cell_parameters] given in Angstrom
        end
    end
    return CD
end


#> pw.x , cp.x
function ATOMIC_SPECIES(config::Dict)
    CD           = Card{Union{Float64, Int64, String}}("ATOMIC_SPECIES", "", Vector{Union{Float64, Int64, String}}[])
    cnt          = Vector{Union{Float64, Int64, String}}[]
    atoms        = unique(first.(config["positions"]))
    pseudo_files = get(config, "pseudo_files", make_pseudofile_dict(config[:pseudo_mode], atoms))
    for (atm, pseudo)  in pseudo_files
        @assert (atm in AtomSymb)
        s = Union{Float64, Int64, String}[atm, AtomMasses[atm], pseudo]
        push!(cnt, s)
    end
    CD.contents = copy(cnt)
    return CD
end


#> pw.x , cp.x
function ATOMIC_POSITIONS(config::Dict)
    CD = Card{Union{Float64, Int64, String}}(
            "ATOMIC_POSITIONS", 
            "crystal", 
            Vector{Union{Float64, Int64, String}}[]
         )
    #if config[:calc]=="bands"   return CD   end
    
    #TODO
    positions = config["positions"]
    @inline wyck_in_pos(l) = typeof(l[2])<:AbstractString
    sg = all(wyck_in_pos.(positions))
    CD.option = sg ? "crystal_sg" : "crystal"

    if sg  #TODO
        #> When atomic positions are of type crystal_sg coordinates can be given
        #>      in the following four forms (Wyckoff positions):
        #>         C  1a
        #>         C  8g   x
        #>         C  24m  x y
        #>         C  48n  x y z
        #>      The first form must be used when the Wyckoff letter determines uniquely
        #>      all three coordinates, forms 2,3,4 when the Wyckoff letter and 1,2,3
        #>      coordinates respectively are needed.
        CD.contents = Vector{Union{Float64, Int64, String}}.(collect.(positions))  # (atm,x,y,z)
    else
        CD.contents = Vector{Union{Float64, Int64, String}}.(collect.(positions))  # (atm,x,y,z)
    end

    return CD

end


#> pw.x , cp.x
function CONSTRAINTS(config::Dict)

    CD = Card{Union{Float64, Int64, String}}("CONSTRAINTS", "", Vector{Union{Float64, Int64, String}}[])

    return CD

end


#> pw.x , cp.x
function OCCUPATIONS(config::Dict)

    CD = Card{PW_Float_Int}("OCCUPATIONS", "", Vector{PW_Float_Int}[])

    return CD

end


#> pw.x , cp.x
function ATOMIC_VELOCITIES(config::Dict)

    CD = Card{PW_Float_Int}("ATOMIC_VELOCITIES", "", Vector{PW_Float_Int}[])

    return CD

end


#> pw.x , cp.x
function ATOMIC_FORCES(config::Dict)

    CD = Card{PW_Float_Int}("ATOMIC_FORCES", "", Vector{PW_Float_Int}[])

    return CD

end


#* =========================================================


# dry run for debug
function pwx_inupt_dry_run(config0::Dict)
    config = copy(config0)
    config[:calc] = "scf"
    return  INPUT_PWX(  config["title"], 
                        map( x->x(config ⬱ ("nstep"=>0)),   #!!!!!!
                              [ CONTROL, SYSTEM, ELECTRONS, IONS, CELL,
                                ATOMIC_SPECIES, ATOMIC_POSITIONS, K_POINTS, CELL_PARAMETERS,
                                CONSTRAINTS, OCCUPATIONS, ATOMIC_VELOCITIES, ATOMIC_FORCES, ] ) ...) |> build_input_file
end


# this is a condensation of experience
function pwx_inupt_compute_energy_or_bands(
    mode_energy_or_bands,
    config::Dict
    )

    CALC_TYPE = config[:calc]
    @assert CALC_TYPE ∈ ["scf", "bands"]

    _KPMODE_  = config[:kpoint_mode]
    _PWMODE_  = config[:pw_mode]
    @assert _PWMODE_ == mode_energy_or_bands
    @assert ( (mode_energy_or_bands=="energy" && (_KPMODE_=="gamma" || _KPMODE_=="automatic")) 
            ||(mode_energy_or_bands=="bands" && (_KPMODE_=="crystal_b" 
                                              || _KPMODE_=="tpiba_b" 
                                              || _KPMODE_=="tpiba_c")) 
            )  "config[\"kpoint_mode\"] = $(config[:kpoint_mode]) is not recognized !!!"

    return  INPUT_PWX(  config["title"],
                         map( x->x(config),
                            [   CONTROL, SYSTEM, ELECTRONS, IONS, CELL,
                                ATOMIC_SPECIES, ATOMIC_POSITIONS, K_POINTS, CELL_PARAMETERS,
                                CONSTRAINTS, OCCUPATIONS, ATOMIC_VELOCITIES, ATOMIC_FORCES,   ] ) ...) |> build_input_file

end


pwx_inupt_compute_energy(config::Dict) = pwx_inupt_compute_energy_or_bands( "energy", config )

pwx_inupt_compute_bands(config::Dict) =  pwx_inupt_compute_energy_or_bands( "bands", config )


function pwx_inupt_nscf(
    config::Dict
    )

    @assert config[:calc] == "nscf"
    @assert config["startingpot"] == "file"

    return  INPUT_PWX(  config["title"],
                         map( x->x(config),
                            [   CONTROL, SYSTEM, ELECTRONS, IONS, CELL,
                                ATOMIC_SPECIES, ATOMIC_POSITIONS, K_POINTS, CELL_PARAMETERS,
                                CONSTRAINTS, OCCUPATIONS, ATOMIC_VELOCITIES, ATOMIC_FORCES,   ] ) ...) |> build_input_file

end


# this is a condensation of experience
function pwx_inupt_relax(config::Dict)
    CALC_TYPE = config[:calc]
    @assert (CALC_TYPE == "relax") || (CALC_TYPE == "vc-relax")

    _KPMODE_  = config[:kpoint_mode]
    @assert (_KPMODE_=="gamma" || _KPMODE_=="automatic")

    return  INPUT_PWX(  config["title"],
                        map( x->x(config),
                            [   CONTROL, SYSTEM, ELECTRONS, IONS, CELL,
                                ATOMIC_SPECIES, ATOMIC_POSITIONS, K_POINTS, CELL_PARAMETERS,
                                CONSTRAINTS, OCCUPATIONS, ATOMIC_VELOCITIES, ATOMIC_FORCES, ] ) ...) |> build_input_file

end

