# cpx_input_file.jl
# ref : https://www.quantum-espresso.org/Doc/INPUT_CP.html


CP_Float_Int = Union{Float64, Int64}

#* ============================================================


function IONS_CP(config::Dict)
    NL = Namelist("IONS",Dict(),Dict(),Dict())
    if "ion_dynamics" ∉ keys(config)
        NL.req = Dict("ion_dynamics"=> "none")
        return NL
    else
        atoms_al= first.(config["positions"])
        atoms   = unique(atoms_al)
        ntyp    = length(atoms)
        NL.req = Dict(
        )
        NL.default = Dict(
            "ion_temperature"   => "not_controlled",
            "ion_dynamics"      => "none",
            "tempw"             => 300.0,
            "fnosep"            => 1.0,
            "ion_nstepe"        => 20,
            "ion_damping"       => 0.1,
            "ion_velocities"    => "default", # ["change_step", "random", "zero", "from_input"]
           ["ion_radius($i)"    => 0.5 for i=1:ntyp]...,
        )
        NL.opt = NL.default ← config
        return NL
    end
end


function CELL_CP(config::Dict)
    NL = Namelist("CELL",Dict(),Dict(),Dict())
    # input this namelist only if calculation = 'vc-relax', 'vc-cp', 'vc-cp-wf'
    if !(config[:calc] ∈ ["vc-relax", "vc-cp", "vc-cp-wf"])
        NL.req = Dict("cell_dynamics"=> "none")
        return NL
    end
    NL.req = Dict(
        "cell_temperature" => config["cell_temperature"],
    )
    NL.default = Dict(
        "cell_dynamics"    => "none",
        "temph"             => 0.0,
        "fnoseh"            => 1.0,

        "wmass"            => nothing,  # unable to calculate that value ...
        "cell_parameters"  => "default",
        "cell_dofree"      => "all",
        "cell_factor"      => 1.2,
        "press"            => 0.0, # [KBar]
    )
    NL.opt = NL.default ← config
    return NL
end


function AUTOPILOT_CP(config::Dict)
    function adapt(x)
        if x isa Bool
            return (x ? ".true." : ".false.")
        elseif x isa AbstractString
            return "'$(x)'"
        elseif x isa Integer
            return string(x)
        elseif x isa Float64
            return strip((@sprintf "%20.15f" x))
        else
            return string(x)
        end
    end

    CD = Card{String}("AUTOPILOT", "", Vector{String}[])
    # example of rules : [(31, "dt", 5.0), (91, "iprint", 100), ...]
    if :pilot_rules ∈ keys(config)
        CD.contents = Vector{String}[ 
                            [   String["on_step = ", string(rule[1]), " : ", string(rule[2]), " = ", adapt(rule[3])] 
                                for rule ∈ config[:pilot_rules]    ]...,
                            ["ENDRULES"]
        ]
    end

    return CD

end


function WANNIER_CP(config::Dict)
    NL = Namelist("WANNIER",Dict(),Dict(),Dict())
    return NL
end


#* =========================================================
#> same as pw.x


SYSTEM_CP(config::Dict) = SYSTEM(config)


ELECTRONS_CP(config::Dict) = ELECTRONS(config)


CONTROL_CP(config::Dict) = CONTROL(config) 


CELL_PARAMETERS_CP(config::Dict) = CELL_PARAMETERS(config)


ATOMIC_SPECIES_CP(config::Dict) = ATOMIC_SPECIES(config)


ATOMIC_POSITIONS_CP(config::Dict) = ATOMIC_POSITIONS(config)


CONSTRAINTS_CP(config::Dict) = CONSTRAINTS(config)


OCCUPATIONS_CP(config::Dict) = OCCUPATIONS(config)


ATOMIC_VELOCITIES_CP(config::Dict) = ATOMIC_VELOCITIES(config)


ATOMIC_FORCES_CP(config::Dict) = ATOMIC_FORCES(config)


#* =========================================================


#TODO more careful verification of input configs
function cpx_inupt_dry_run(config0::Dict)
    config = copy(config0)
    config[:calc] = "cp"
    return  INPUT_CPX(  config["title"], 
                        map( x->x(config ⬱ ("nstep"=>0)),   #!!!!!!
                             [ CONTROL_CP, SYSTEM_CP, ELECTRONS_CP, IONS_CP, CELL_CP, WANNIER_CP,
                               ATOMIC_SPECIES_CP, ATOMIC_POSITIONS_CP, CELL_PARAMETERS_CP,
                               OCCUPATIONS_CP, CONSTRAINTS_CP, ATOMIC_VELOCITIES_CP, ATOMIC_FORCES_CP,
                               AUTOPILOT_CP, ] ) ...) |> build_input_file
end


#TODO more careful verification of input configs
function cpx_inupt_MD(
    config::Dict
    )
    CALC_TYPE = config[:calc]
    @assert CALC_TYPE ∈ ["cp", "vc-cp"]

    return  INPUT_CPX(  config["title"],
                        map( x->x(config),
                            [  CONTROL_CP, SYSTEM_CP, ELECTRONS_CP, IONS_CP, CELL_CP, WANNIER_CP,
                               ATOMIC_SPECIES_CP, ATOMIC_POSITIONS_CP, CELL_PARAMETERS_CP,
                               OCCUPATIONS_CP, CONSTRAINTS_CP, ATOMIC_VELOCITIES_CP, ATOMIC_FORCES_CP,
                               AUTOPILOT_CP, ] ) ...) |> build_input_file

end


#TODO more careful verification of input configs
function cpx_inupt_energy(
    config::Dict
    )
    CALC_TYPE = config[:calc]
    @assert CALC_TYPE ∈ ["scf", "nscf"]

    return  INPUT_CPX(  config["title"],
                        map( x->x(config),
                            [  CONTROL_CP, SYSTEM_CP, ELECTRONS_CP, IONS_CP, CELL_CP, WANNIER_CP,
                               ATOMIC_SPECIES_CP, ATOMIC_POSITIONS_CP, CELL_PARAMETERS_CP,
                               OCCUPATIONS_CP, CONSTRAINTS_CP, ATOMIC_VELOCITIES_CP, ATOMIC_FORCES_CP,
                               AUTOPILOT_CP, ] ) ...) |> build_input_file

end


#TODO more careful verification of input configs
function cpx_inupt_relax(config::Dict)
    CALC_TYPE = config[:calc]
    @assert CALC_TYPE ∈ ["relax", "vc-relax"]

    return  INPUT_CPX(  config["title"],
                        map( x->x(config),
                            [  CONTROL_CP, SYSTEM_CP, ELECTRONS_CP, IONS_CP, CELL_CP, WANNIER_CP,
                               ATOMIC_SPECIES_CP, ATOMIC_POSITIONS_CP, CELL_PARAMETERS_CP,
                               OCCUPATIONS_CP, CONSTRAINTS_CP, ATOMIC_VELOCITIES_CP, ATOMIC_FORCES_CP,
                               AUTOPILOT_CP, ] ) ...) |> build_input_file

end

