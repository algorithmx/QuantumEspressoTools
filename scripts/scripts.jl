ENV["ISOSUITE_FOLDER"] = "/home/dabajabaza/abinitio/iso/"

ENV["QE_DEV"] = "/home/dabajabaza/jianguoyun/Nutstore/QuantumEspressoTools"

using Distributed

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
import QuantumEspressoTools.↓
import QuantumEspressoTools.↑

if @isdefined Isosuite
    nothing
else
    if (@isdefined use_dev) && use_dev
        using Pkg
        Pkg.activate("/home/dabajabaza/jianguoyun/Nutstore/Isosuite")
        using Isosuite
    else
        using Isosuite
    end
end

I_am_here = Base.@__DIR__
include("$I_am_here/logger.jl")
include("$I_am_here/instruction.jl")

include("$I_am_here/QEResult.jl")

include("$I_am_here/analyze_prog.jl")
include("$I_am_here/check_status.jl")
include("$I_am_here/get_structure_from_cif.jl")
include("$I_am_here/random_kick.jl")
include("$I_am_here/generate_input_script.jl")
include("$I_am_here/launch_program.jl")
include("$I_am_here/execute.jl")
include("$I_am_here/seeds.jl")
include("$I_am_here/interprete_results.jl")
include("$I_am_here/grid_test.jl")
include("$I_am_here/common_folder_names.jl")
include("$I_am_here/find_window.jl")
include("$I_am_here/cleanup.jl")