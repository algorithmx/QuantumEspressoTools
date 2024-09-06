export UNITCELL, CELLDM, LATTPARAM
export basis_to_lattice_parameters, lattice_parameters_to_basis, latt_param_h2r
export check_valid_celldm, celldm_to_basis_Bohr, celldm_to_basis_Angstrom
export celldm_vec, celldm_dic
export Find_ibrav_for_IT, consistent
export convert_unit_to_Ang, convert_unit_to_Bohr
export cell_to_basis_Angstrom, cell_to_basis_Bohr
export check_compat_QE_ibrav_convention, check_valid_cell
export basis_to_cell, cell_to_lattparam, lattparam_to_cell
export QE_celldm_SG_to_cell, cell_to_QE_celldm_ibrav, QE_celldm_ibrav_to_cell, config_to_cell
export cell_to_QE_pwx_input


##--------------------------------------------------------------------------------

#> NOTE: upto rotation,  

#> A1,A2,A3 ---- basis_to_lattice_parameters ---> (a,b,c,alpha,beta,gamma) 
#> A1,A2,A3 <--- lattice_parameters_to_basis ---- (a,b,c,alpha,beta,gamma) 
#> SG>2 + conventions ---> ibrav

#> SG  <--->  (A1,A2,A3) are not 1:1 correspondance, but
#> SG>2 + conventions ---> ibrav ----QE----> (A1,A2,A3)
#> reverse is not true

#: lattp ----> celldm  and  ibrav -----> basis follows DIFFERENT cell conventions
#: lattp <---- basis  <---- celldm  and  ibrav <--\-- basis 
#: cannot get celldm from lattp directly, ibrav is necessary
#! a in lattice parameter is NOT alat !!!!!!! (consider ibrav=2)

#>NOTE
#> It is usually convenient to start with LATTPARAM since
#> a cif file usually starts with (a,b,c,α,β,γ)
#> the danger is, cif file may follow an UNKNOWN convention
#> or, even the convention is specified, it is INCOMPATIBLE 
#> with the QE convention. Moreover, the QE convention itself is wierd.
#> Any of the above factors may FUCK UP the calculation. 
#: The program is NOT RESPONSIBLE for conventions used in the cif file.

#--------------------------------------------------------------------------------


mutable struct UNITCELL
    UNIT::Symbol                #> :Angstrom or :Bohr
    SG::Int                     #> 1 --230
    ALAT::Float64               #> alat scale parameter, should be combined with UNIT
    A1_ALAT::Vector{Float64}    #> three basis in unit of alat
    A2_ALAT::Vector{Float64}
    A3_ALAT::Vector{Float64}
end


mutable struct CELLDM
    ibrav::Int
    v::Vector{Float64}
end


mutable struct LATTPARAM
    SG::Int
    p::NTuple{6,Real}
end

import Base.isapprox
import Base.copy
import Base.isequal
import Base.==

isapprox(a::NTuple{6,Float64}, b::NTuple{6,Float64}) = all((x≈y) for (x,y) ∈ zip(a,b))

isapprox(a::NTuple{6,Float64}, b::Tuple) = (length(a)==length(b) && all((x≈y) for (x,y) ∈ zip(a,b)))

isapprox(b::Tuple, a::NTuple{6,Float64}) = isapprox(b,a)

include("crystal_structure/basis_operations.jl")

include("crystal_structure/ibrav.jl")

include("crystal_structure/consistent.jl")

include("crystal_structure/UNITCELL.jl")

include("crystal_structure/LATTPARAM.jl")

include("crystal_structure/CELLDM.jl")

include("crystal_structure/check_adapt_lattice_params_to_QE.jl")

include("crystal_structure/QE_pwx_input.jl")

#include("crystal_structure/test.jl")

##

# Find_ibrav_for_IT.(1:230, false) |> print