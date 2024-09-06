#pseudopotential.jl


export ELECTRON_CONFIG
export num_electrons, max_electrons
export conventional_valence_electron_config
export Find_suggested_cutoff
export PSEUDO_PATHS, PSEUDO_DATA
export make_pseudo_wcut, make_pseudo_ecut, make_pseudofile_dict
export num_valence_electrons, num_KS_states
export psl_all_possible_configs
export pseudo_mode_name, take_part_from_str, PSL_name_pseudo_mode


#: --------------------------------------------------------------



if "PSEUDO_ROOT" ∉ keys(ENV)
    @warn "Please specify ENV[\"PSEUDO_ROOT\"] before using QuantumEspressoTools."
    ENV["PSEUDO_ROOT"] = "/home/dabajabaza/abinitio/pseudopotentials"
end



mutable struct PSEUDO
    ROOT::String
    
    PATHS::Dict{String,String}

    # general format of pseudopotential dict
    # (  "file name",  wave function cutoff,  charge cutoff  )
    DATA::Dict{String,Any}
end


include("pseudopotential/electron_config.data.jl")

include("pseudopotential/electron_config.jl")

include("pseudopotential/paths.jl")

include("pseudopotential/build_load_verify.jl")

include("pseudopotential/PSLibrary.jl")

include("pseudopotential/find_in_ps_file.jl")

include("pseudopotential/Find_suggested_cutoff.jl")



global const PSINFO = build_load_pseudo(PSEUDO_PATHS)

global const PSEUDO_DATA = copy(PSINFO.DATA)


#* USER SPECIFIED PSEUDOPOTENTIALS

include("pseudopotential/user_defined.jl")

global const User_defined_pseudo = copy(user_pseudo(ENV))


include("pseudopotential/make_pseudo.jl")



#* num_valence_electrons

include("pseudopotential/nbnd.jl")
