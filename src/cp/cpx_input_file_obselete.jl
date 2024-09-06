# cpx_input_file.jl
# TODO finish the cards and the vc-md input generator
# ref : https://www.quantum-espresso.org/Doc/INPUT_CP.html

## ----------------------------------------------------------------
## SECTIONS
## ----------------------------------------------------------------


function CELL_PARAMETERS(config)
    cell_px = [0.0,0.0,0.0]
    cell_py = [0.0,0.0,0.0]
    cell_pz = [0.0,0.0,0.0]
    if :cell_parameters in keys(config)
        cell_px = config[:cell_parameters][1]  # in Angstrom
        cell_py = config[:cell_parameters][2]  # in Angstrom
        cell_pz = config[:cell_parameters][3]  # in Angstrom
    else
        CIF                   = (config[:cif] isa Dict) ? config[:cif] : to_cif(config[:cif]...)
        a_Ang, b_Ang, c_Ang   = strtonum.([CIF["_cell_length_a"], CIF["_cell_length_b"], CIF["_cell_length_c"]])
        alpha, beta, gamma    = strtonum.([CIF["_cell_angle_alpha"], CIF["_cell_angle_beta"], CIF["_cell_angle_gamma"]])
        alphar, betar, gammar = map(x->Float64((π/big(180))*x), [alpha, beta, gamma])
        IT_num                = (config[:cif] isa Dict) ? get(CIF,"_symmetry_Int_Tables_number",0) : (length(config[:cif])==7 ? config[:cif][7] : 0)
        cell_px = [a_Ang,   b_Ang*cos(gammar),  c_Ang*cos(betar)]
        cell_py = [0.0,     b_Ang*sin(gammar),  c_Ang*(cos(alphar)-cos(betar)*cos(gammar))/sin(gammar)]
        cell_pz = [0.0,     0.0,                c_Ang*sqrt(1.0 - cos(alphar)^2 - cos(betar)^2 - cos(gammar)^2 + 2*cos(alphar)*cos(betar)*cos(gammar))/sin(gammar)]
    end

    lines = ["\nCELL_PARAMETERS {bohr}", ]
    for i=1:3
        s = (@sprintf "  %19.14f   %19.14f   %19.14f"  cell_px[i]/_BOHR_RADIUS_  cell_py[i]/_BOHR_RADIUS_  cell_pz[i]/_BOHR_RADIUS_)
        push!(lines, s)
    end
    push!(lines, "\n")
    return lines
end


function CONTROL(config)

    CALC_TYPE       = config["calc"]
    TITLE           = config["title"] * "_" * CALC_TYPE
    FILE_PREFIX     = config["prefix"]

    TSTRESS         = get(config,"tstress",true)
    TPRNFOR         = get(config,"tprnfor",true)
    PSEUDO_FD       = get(config,"pseudo_folder",   "./"        )
    OUT_FD          = get(config,"outdir",   "./"        )
    DISK_IO         = get(config,"disk_io",         "default"   )
    VERBOSITY       = get(config,"verbosity",       "high"      )
    EN_THR          = get(config,"etot_conv_thr",   0.0001      )
    F_THR           = get(config,"forc_conv_thr",   0.001       )

    lines = [   "&CONTROL" ,
                "                       title = '" * TITLE * "'" ,
                "                 calculation = '" * CALC_TYPE * "'" ,
                "                restart_mode = 'from_scratch'" ,  #! in rare cases use 'restart'
                "                      outdir = '$OUT_FD'" ,
                "                  pseudo_dir = '$PSEUDO_FD'" ,
                "                      prefix = '$FILE_PREFIX'" ,
                "                     disk_io = '$DISK_IO'" ,
                "                   verbosity = '$VERBOSITY'" ,
                "                     tstress = ." * (TSTRESS ? "true" : "false") * "." ,
                "                     tprnfor = ." * (TPRNFOR ? "true" : "false") * "." ,
                "               etot_conv_thr = $EN_THR" ,     # Convergence threshold on total energy (a.u) for ionic minimization
                "               forc_conv_thr = $F_THR" ,      # Convergence threshold on forces (a.u) for ionic minimization
                "/", ]
    
    return lines
end


function ELECTRONS(config)

    ESTEP           = get(config, "electron_maxstep",  500    )
    CONVERGENCE     = get(config, "conv_thr", get(config, "convergence_threshold", 1e-8))
    DIAG_THR        = get(config, "diago_thr_init",    1e-4   )
    MIX_BETA        = get(config, "mixing_beta",       1e-8   )
    MIX_NDIM        = get(config, "mixing_ndim",       8      )
    DIAG_METHOD     = get(config, "diagonalization",   "david")

    lines = [   "&ELECTRONS",
                "            electron_maxstep = $ESTEP",
                "                    conv_thr = $CONVERGENCE",
                "              diago_thr_init = $DIAG_THR",
                "                 startingpot = 'atomic'",
                "                 startingwfc = 'atomic'",
                "                 mixing_mode = 'plain'",
                "                 mixing_beta = $MIX_BETA",
                "                 mixing_ndim = $MIX_NDIM",
                "             diagonalization = '$DIAG_METHOD'",
                "/",   ]

    return lines
end


function SYSTEM(config; compute_ibrav=true, use_IT=false, find_sugg=false)
    #! if compute_ibrav=false, ibrav = 0 => must present celldm or CELL_PARAMETERS
    #! IT_num and ibrav are independent!

    PSEUDO_FILES    = config["pseudo_files"]
    PSEUDO_FD       = get(config, "pseudo_folder", "./")
    FUNCTIONAL      = get(config, "input_dft", nothing )
    sugg_ecuts      = find_sugg ? [Find_suggested_cutoff(PSEUDO_FD*pf) for pf in values(PSEUDO_FILES)] : [(-1.0,-1.0) for pf in values(PSEUDO_FILES)]
    
    ecutwfc         = max(config["ecut_wavefunction"], maximum(first.(sugg_ecuts)))
    ecutrho         = max(config["ecut_charge"], maximum(last.(sugg_ecuts)))
    
    #* unit cell conventions 
    uniqueb         = get(config,   "uniqueb",        false )
    origin_choice   = get(config,   "origin_choice",  1     )
    rhombohedral    = get(config,   "rhombohedral",   true  )

    DIM             = get(config,   "dimension",      3     )
    NATOM           = length(config["positions"])
    NTYPE           = length(config["pseudo_files"])

    OCC             = config["occupation"]
    SMEARING        = get(config, "smearing",   "gaussian"  )
    DEGAUSS         = get(config, "degauss",    0.0         )

    #TODO more complete settings for van der Waals correction
    VdW_CORR_MODE   = get(config, "vdw_corr",    "none"     )
    xdm_a1          = get(config, "xdm_a1",      1.2153     )
    xdm_a2          = get(config, "xdm_a2",      2.3704     )

    lines0 = [      "                     uniqueb = .$uniqueb.",
                    "               origin_choice = $origin_choice",
                    "                rhombohedral = .$rhombohedral.",
                    DIM==2 ? "             assume_isolated = '2D'" : "",
                    (FUNCTIONAL!==nothing && FUNCTIONAL!=="") ? "                   input_dft = '$FUNCTIONAL'" : "",
                    "                 occupations = '$OCC'",
                    "                    smearing = '$SMEARING'",
                    "                     degauss = $DEGAUSS",  ## 1Ry = 13.6056923(12)eV
                    "                         nat = $NATOM",    ## number of atoms in the unit cell (ALL atoms, except if space_group is set, in which case, INEQUIVALENT atoms),
                    "                        ntyp = $NTYPE",    ## number of types of atoms in the unit cell,
                    "                     ecutwfc = $ecutwfc",  # kinetic energy cutoff (Ry) for wavefunctions,
                    "                     ecutrho = $ecutrho",  # Kinetic energy cutoff (Ry) for charge density and potential,
                    "                    vdw_corr = '$VdW_CORR_MODE'", #
                    "                      xdm_a1 = $xdm_a1",
                    "                      xdm_a2 = $xdm_a2",
                    "/",   
    ]

    if (:cif ∈ keys(config)) && (:cell_parameters ∉ keys(config))
        CIF             = (config[:cif] isa Dict) ? config[:cif] : to_cif(config[:cif]...)
        a_Ang, b_Ang, c_Ang  = strtonum.([CIF["_cell_length_a"], CIF["_cell_length_b"], CIF["_cell_length_c"]])
        alpha, beta, gamma   = strtonum.([CIF["_cell_angle_alpha"], CIF["_cell_angle_beta"], CIF["_cell_angle_gamma"]])
        alphar, betar, gammar = map(x->Float64((π/big(180))*x), [alpha, beta, gamma])
        IT_num          = (config[:cif] isa Dict) ? get(CIF,"_symmetry_Int_Tables_number",0) : (length(config[:cif])==7 ? config[:cif][7] : 0)

        #* depends on unit cell conventions
        latt            = Find_Lattice(a_Ang,b_Ang,c_Ang, alpha,beta,gamma, uniqueb)
        ibrav           = ((compute_ibrav && IT_num>0) ? Find_ibrav(Int_Tables[IT_num], latt, uniqueb) : 0)  #! if compute_ibrav=false, ibrav = 0
        @info "lattice type : $latt  ibrav = $ibrav"
        cd              = celldm(ibrav, (a_Ang, b_Ang, c_Ang, alphar, betar, gammar))
        lines = [  ["&SYSTEM",
                    "                       ibrav = $ibrav",
                    cd,
                    (use_IT && IT_num>1) ? "                 space_group = $IT_num" : ""];
                    lines0;
        ]
        return lines
    elseif (:cell_parameters ∈ keys(config))  #! supercedes :cif
        lines = [  ["&SYSTEM",
                    "                       ibrav = 0" ];
                    lines0;
        ]
        return lines
    end
end


function IONS(config)
    ION_DYN = get(config, "ion_dynamics", "")
    #* BFGS specific
    # Max reduction factor for conv_thr during structural optimization
    MAX_REDUC = get(config, "upscale", 100.0)
    # Number of old forces and displacements vectors used in the PULAY mixing
    NDIM = get(config,"bfgs_ndim",1)
    # displacement max, min and init
    R_max = get(config,"trust_radius_max",0.8)
    R_min = get(config,"trust_radius_min",0.001)
    R_0   = get(config,"trust_radius_ini",0.5)
    #* extrapolations
    pot_xtrpl = get(config,"pot_extrapolation","none")
    wfc_xtrpl = get(config,"wfc_extrapolation","none")
    #* variables used for molecular dynamics
    ion_T = get(config,"ion_temperature","not_controlled")
    N_raise = get(config,"nraise",1)

    lines = [   "&IONS",
                (length(ION_DYN)>0) ? "                ion_dynamics = '$ION_DYN'" : "",
                (lowercase(ION_DYN)=="bfgs") ? "                       upscale = $MAX_REDUC" : "",
                (lowercase(ION_DYN)=="bfgs") ? "                     bfgs_ndim = $NDIM" : "",
                (lowercase(ION_DYN)=="bfgs") ? "              trust_radius_max = $R_max" : "",
                (lowercase(ION_DYN)=="bfgs") ? "              trust_radius_min = $R_min" : "",
                (lowercase(ION_DYN)=="bfgs") ? "              trust_radius_ini = $R_0" : "",
                "             pot_extrapolation = '$pot_xtrpl'",
                "             wfc_extrapolation = '$wfc_xtrpl'",
                "               ion_temperature = '$ion_T'",
                "                        nraise = $N_raise",
                "/",   ]

    return lines
end 


function CELL(config)
    cell_dynamics           = get(config, "cell_dynamics",  "none"  )
    if config["calc"]=="vc-relax"
        @assert lowercase(cell_dynamics) in ["none", "sd", "damp-pr", "damp-w", "bfgs"]  "$(cell_dynamics) not allowed for vc-relax."
    elseif config["calc"]=="vc-md"
        @assert lowercase(cell_dynamics) in ["none", "pr", "w"]  "$(cell_dynamics) not allowed for vc-md."
    else
        @info "CELL card not necessary for calculation $(config["calc"])."
        return String[]
    end

    cell_dofree             = get(config, "cell_dofree",    "all" )
    target_pressure_KBar    = get(config, "press",          0.0   ) # [KBar]
    press_conv_thr          = get(config, "press_conv_thr", 0.1   )

    lines = [   "&CELL",
                "               cell_dynamics = '$cell_dynamics'",
                "                 cell_dofree = '$cell_dofree'",
                (@sprintf "                       press = %12.8f" target_pressure_KBar),
                (@sprintf "              press_conv_thr = %12.8f" press_conv_thr),
                "/",   ]

    return lines
end


function ATOMIC_SPECIES(config)
    lines = ["\nATOMIC_SPECIES",]
    for (atm, pseudo)  in config["pseudo_files"]
        @assert (atm in AtomSymb)
        s = (@sprintf "  %3s  %14.10f  %s"  atm  Atoms[atm]  pseudo)
        push!(lines, s)
    end
    return lines
end


function ATOMIC_POSITIONS(config)
    lines = ["\nATOMIC_POSITIONS crystal", ]
    for (atm,x,y,z) in config["positions"]
        s = ( @sprintf "%3s   %19.14f   %19.14f   %19.14f"  atm  x  y  z )
        push!(lines, s)
    end
    return lines
end

function ATOMIC_SPECIES(config)

end

# When atomic positions are of type crystal_sg coordinates can be given
#      in the following four forms (Wyckoff positions):
#         C  1a
#         C  8g   x
#         C  24m  x y
#         C  48n  x y z
#      The first form must be used when the Wyckoff letter determines uniquely
#      all three coordinates, forms 2,3,4 when the Wyckoff letter and 1,2,3
#      coordinates respectively are needed.


function ATOMIC_POSITIONS_SG(config)
    lines = ["\nATOMIC_POSITIONS crystal_sg"]
    # @info "ATOMIC_POSITIONS_SG():\n" * (join(string.(config["positions"]),"\n"))
    for tup in config["positions"]
        #! MODIFIED !!!
        (atm, nwyck) = tup[1:2]
        p1 = length(tup)>2 ? (@sprintf "%19.14f   " tup[3]) : ("")
        p2 = length(tup)>3 ? (@sprintf "%19.14f   " tup[4]) : ("")
        p3 = length(tup)>4 ? (@sprintf "%19.14f   " tup[5]) : ("")
        s = (@sprintf "%3s   %4s   %s"  atm  nwyck  p1*p2*p3)
        push!(lines, s)
    end
    return lines
end

## =======================================================================


# this is a condensation of experience
function cpx_inupt_CP(
    config::Dict;
    compute_ibrav = true, # very helpful for the calculation 
    use_IT = false, # easy to fail due to different conventions
    find_sugg = false
    )
    @assert config["calc"] == "cp"
    IT_num  = (config[:cif] isa Dict) ? get(CIF,"_symmetry_Int_Tables_number",0) : (length(config[:cif])==7 ? config[:cif][7] : 0)
    return  join(
                    vcat(

                        CONTROL(config),

                        SYSTEM(config, compute_ibrav=compute_ibrav, use_IT=use_IT, find_sugg=find_sugg),

                        ELECTRONS(config),

                        IONS(config), # necessary only for 'relax' or 'vc-relax'

                        CELL(config),

                        ATOMIC_SPECIES(config),

                        #! IT = 1 supercedes use_IT
                        (use_IT && IT_num>1) ? ATOMIC_POSITIONS_SG(config) : ATOMIC_POSITIONS(config),

                        compute_ibrav ? String[] : CELL_PARAMETERS(config)
                    ),
            "\n")

end


# this is a condensation of experience
function cpx_inupt_relax(
    config::Dict;
    compute_ibrav = true,       # very helpful for the calculation 
    use_IT = false,             # easy to fail due to different conventions
    find_sugg = false
    )
    CALC_TYPE = config["calc"]
    @assert (CALC_TYPE == "relax") || (CALC_TYPE == "vc-relax")
    IT_num  = (config[:cif] isa Dict) ? get(CIF,"_symmetry_Int_Tables_number",0) : (length(config[:cif])==7 ? config[:cif][7] : 0)
    return  join(
                    delete_empty_lines(vcat(

                        CONTROL(config),

                        SYSTEM(config, compute_ibrav=compute_ibrav, use_IT=use_IT, find_sugg=find_sugg),

                        ELECTRONS(config),

                        IONS(config), 

                        (CALC_TYPE=="vc-relax") ? CELL(config) : "",

                        ATOMIC_SPECIES(config),

                        #! IT = 1 supercedes use_IT
                        (use_IT && IT_num>1) ? ATOMIC_POSITIONS_SG(config) : ATOMIC_POSITIONS(config),

                        compute_ibrav ? String[] : CELL_PARAMETERS(config)
                    )),  
            "\n")

end

