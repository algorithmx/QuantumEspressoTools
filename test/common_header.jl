ENV["ISOSUITE_FOLDER"] = "/home/dabajabaza/abinitio/iso/"

ENV["QE_DEV"] = "/home/dabajabaza/jianguoyun/Workspace/QuantumEspressoTools"

using LinearAlgebra

using JLD2

if @isdefined QuantumEspressoTools
    nothing
else
    using Pkg
    Pkg.activate(ENV["QE_DEV"])
    using QuantumEspressoTools
end

import QuantumEspressoTools: ∪

if @isdefined Isosuite
    nothing
else
    if (@isdefined use_dev) && use_dev
        using Pkg
        Pkg.activate("/home/dabajabaza/jianguoyun/Workspace/Isosuite")
        using Isosuite
    else
        using Isosuite
    end
end

I_am_here = Base.@__DIR__
