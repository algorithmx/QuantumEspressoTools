

function make_pseudo_(psmode::String, atom_list::Vector{String}, i::Int)
    if psmode=="USER_DEFINED"
        if length(User_defined_pseudo)==0
            throw(error("ENV[\"USER_PSEUDO\"] not valid or not specified !!!"))
        else
            @assert all([a ∈ keys(User_defined_pseudo) for a in atom_list])  "make_pseudo_() : User_defined_pseudo does not contain pseudopotential files for all atoms in $atom_list ."
            return Dict(a=>User_defined_pseudo[a][i] for a in atom_list)
        end
    else
        @assert all([a ∈ keys(PSEUDO_DATA[psmode]) for a in atom_list])  "make_pseudo_() : $psmode does not contain pseudopotential files for all atoms in $atom_list ."
        return Dict(a=>PSEUDO_DATA[psmode][a][i] for a in atom_list)
    end
end


PSL_PS_MODE{T} = Tuple{String, Dict{String,Tuple{T,String}}}


function make_pseudo_(psmode_t::PSL_PS_MODE{T}, atom_list::Vector{String}, i::Int) where {T<:Union{String,Symbol}}
    (psmode, spec0) = psmode_t
    spec = Dict(a=>("-$(s[1])-","psl.$(s[2]).UPF") for (a,s) in spec0)
    @assert all([a ∈ keys(PSEUDO_DATA[psmode]) for a in atom_list])  "make_pseudo_() : $psmode does not contain pseudopotential files for all atoms in $atom_list ."
    @assert all([a ∈ keys(spec)                for a in atom_list])  "make_pseudo_() : $spec0 does not fit pseudopotential files for all atoms in $atom_list ."
    @assert all([f ∈ ALL_PSL_SETTINGS for f in values(spec)])        "make_pseudo_() : $spec0 contains unsupported specifications."
    return Dict(a=>PSEUDO_DATA[psmode][a][spec[a]][i] for a in atom_list)
end


make_pseudofile_dict(psmode::Union{String,PSL_PS_MODE}, atom_list::Vector{String}) = make_pseudo_(psmode, atom_list, 1)

make_pseudofile_dict(psmode::PSL_PS_MODE) = make_pseudofile_dict(psmode, unique(keys(psmode)))

make_pseudo_wcut(psmode::Union{String,PSL_PS_MODE}, atom_list::Vector{String})     = maximum(values(make_pseudo_(psmode,atom_list,2)))

make_pseudo_ecut(psmode::Union{String,PSL_PS_MODE}, atom_list::Vector{String})     = maximum(values(make_pseudo_(psmode,atom_list,3)))
