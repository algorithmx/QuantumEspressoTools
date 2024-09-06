__precompile__()

module QuantumEspressoTools

__VERIFY__PS__ = false 
__MODIFY__PS__ = false


# using JLD2  #TODO (to test) replace everything by JSON !!!
using LinearAlgebra
using Printf
using JSON
using Isosuite


#* DATA
export _BOHR_RADIUS_, AtomMasses, Atoms
export IT_NUMBER_CRYSTAL_CLASS, CRYSTAL_CLASS, Int_Tables, Int_Tables_No
export KVECS_BY_CRYSTAL_CLASS_QE_DEFAULT
include("data.jl")

export convert_setting_to_findsym
export QE_default_equivalent_settings_findsym
export QE_default_symmetry_group_convention
include("default_settings.jl")

#*  Utility tools 
include("utils.jl")


#* CRYSTAL STRUCTURE
include("crystal_structure.jl")


#* pseudopotentials
include("pseudopotential.jl")
#TODO (to test) replace everything by JSON !!!


#* RUN 
include("verify_result.jl")
include("run.jl")


#* Data Structures 
include("namelists_cards.jl")


#* pw.x
include("pw.jl")


#* cp.x
include("cp.jl")


#* pp.x
include("pp.jl")


#* ph.x
include("ph.jl")


#* ld1.x
include("ld1.jl")


#* VASP
include("vasp.jl")


#* GULP
include("GULP.jl")


#* MaterialsProject 
include("MaterialsProject.jl")


include("wannier90.jl")




end # module
