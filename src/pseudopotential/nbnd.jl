
#TODO wyckoff positions ?

function num____electrons(
    psmode::String, 
    atom_list_all::Vector{String}, 
    fff
    )::Int
    @assert all([a ∈ keys(PSEUDO_DATA[psmode]) for a in unique(atom_list_all)])    "num____electrons() : $psmode does not contain pseudopotential files for all atoms in $(unique(atom_list_all)) ."
    s1 = sum([fff(ELECTRON_CONFIG(a, PSEUDO_PATHS[psmode]*"/"*PSEUDO_DATA[psmode][a][1],psmode)) for a in atom_list_all])
    if s1 != 0
        return s1
    else
        s2 = sum([fff(conventional_valence_electron_config[Atoms[a]]) for a in atom_list_all])
        return s2
    end
end


function num____electrons(
    psmode_t::PSL_PS_MODE{T}, 
    atom_list_all::Vector{String}, 
    fff
    )::Int where {T<:Union{String,Symbol}}
    (psmode, spec0) = psmode_t
    spec = Dict(a=>("-$(s[1])-","psl.$(s[2]).UPF") for (a,s) in spec0)
    @assert all([a ∈ keys(PSEUDO_DATA[psmode]) for a in unique(atom_list_all)])  "num____electrons() : $psmode does not contain pseudopotential files for all atoms in $(unique(atom_list_all)) ."
    @assert all([a ∈ keys(spec)                for a in unique(atom_list_all)])  "num____electrons() : $spec0 does not fit pseudopotential files for all atoms in $(unique(atom_list_all)) ."
    @assert all([f ∈ ALL_PSL_SETTINGS for f in values(spec)])                    "num____electrons() : $spec0 contains unsupported specifications."
    s1 = sum([fff(ELECTRON_CONFIG(a, PSEUDO_PATHS[psmode]*"/"*PSEUDO_DATA[psmode][a][spec[a]][1],psmode)) for a in atom_list_all])
    if s1 != 0
        return s1
    else
        s2 = sum([fff(conventional_valence_electron_config[Atoms[a]]) for a in atom_list_all])
        return s2
    end
end


# fractional occupatiosn rounded up !
num_valence_electrons(psmode, atom_list_all::Vector{String})::Int = num____electrons(psmode, atom_list_all, num_electrons)


num_KS_states(psmode, atom_list_all::Vector{String})::Int = num____electrons(psmode, atom_list_all, max_electrons)
