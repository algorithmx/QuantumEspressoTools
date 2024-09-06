"""
    get_structure_from_cif(cif0_fn::String)

IMPORTANT ASSUMPTION ON CIF FILE CONVENTION
TODO

input  : cif file name
output : initial configuration, contains key
         "posiitons", :cif (with SG = 1), :cell_parameters and 
         :do_not_use_symmetry => true
note   : uses 
         QuantumEspressoTools.lattice_parameters_to_basis
         Isosuite.symmetry_operators
         Isosuite.get_atom_frac_pos
         Isosuite.extend_positions

"""
function get_structure_from_cif(cif0_fn::String; ntol=9)
    (a, b, c, α, β, γ) = get_cell_params(cif0_fn)
    basis = QuantumEspressoTools.lattice_parameters_to_basis(QuantumEspressoTools.LATTPARAM((a, b, c, α, β, γ)))
    symm_ops = Isosuite.symmetry_operators(cif0_fn)
    atom_list = Isosuite.get_atom_frac_pos(cif0_fn)
    atom_list_ext = (length(symm_ops)>0) ? Isosuite.extend_positions(atom_list, symm_ops) : atom_list
    @inline tovec(x) = [x...,]
    @inline rd(x) = round(x,digits=ntol)
    atm = tovec.(atom_list_ext)
    return Dict("positions"          => atm,
                :cif                 => (a, b, c, α, β, γ, 1),
                :cell_parameters     => [rd.(basis[1]),rd.(basis[2]),rd.(basis[3])],
                :do_not_use_symmetry => true)
end


#TODO test it on ALL materials in the Materials Project library

##

#=
using PyPlot
using BSON
using LsqFit
using LightXML
using XMLDict
using LinearAlgebra
using SparseArrays
using PyCall
optimize = pyimport("scipy.optimize")

SFT = "/home/dabajabaza/jianguoyun/Nutstore"
ENV["JULIA_PKG_SERVER"] = "https://mirrors.tuna.tsinghua.edu.cn/julia"

using Pkg

Pkg.activate("$SFT/LatticeLab")
using LatticeLab

Pkg.activate("$SFT/BandStructures")
using BandStructures

Pkg.activate("$SFT/QuantumEspressoTools")
using QuantumEspressoTools
include("$SFT/QuantumEspressoTools/scripts/scripts.jl")
include("$SFT/QuantumEspressoTools/scripts/energy_bands_on_kgrid.jl")

Pkg.activate("$SFT/cifio")
using cifio


##

get_structure_from_cif("/home/dabajabaza/jianguoyun/Workspace/ScF3/materials/ABO3/mp-22540.cif")

=#