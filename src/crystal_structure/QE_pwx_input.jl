"""
    config_to_cell(config::Dict)::UNITCELL

input  : `config` is the configuration, may contain keys "celldm(1)", :cell_parameters or :cif
output : internal structure `UNITCELL`
note   : priority of the keys in config : celldm(1) > :cell_parameters > :cif

"""
function config_to_cell(config::Dict)::UNITCELL
    @assert  check_valid_config(config)
    SG     = get(config,"space_group",1)
    if "celldm(1)" ∈ keys(config)
        return UNITCELL((config, config["ibrav"]), SG)
    elseif :cell_parameters ∈ keys(config)
        basis3 = [  config[:cell_parameters][1],
                    config[:cell_parameters][2],
                    config[:cell_parameters][3]  ]
        return UNITCELL(:Angstrom, SG, norm_basis(basis3)...)
    elseif :cif ∈ keys(config)
        return UNITCELL(config_to_lattp(config), get_sign_of_ibrav(config))
    else
        throw(error("config_to_cell(): no celldm(1) no :cif no :cell_parameters. What do you want to do ???"))
    end
end



"""
    cell_to_QE_pwx_input(uc::UNITCELL, ibrv_ngtv_sgn::Bool)

input  :  `uc`, `ibrv_ngtv_sgn`
output :  configuration for the function `pw_input()`
note   : space group = 1 is the special case 

"""
function cell_to_QE_pwx_input(uc::UNITCELL, ibrv_ngtv_sgn::Bool)::Dict
    if uc.SG<2
        return Dict(:cell_parameters=>cell_to_basis_Angstrom(uc), "ibrav"=>0)
    else
        cd = CELLDM(uc, ibrv_ngtv_sgn)
        uniqueb_dic = ((abs(cd.ibrav)==12 || abs(cd.ibrav)==13) ? Dict("uniqueb"=>ibrv_ngtv_sgn) : Dict())
        return celldm_dic(cd)  ∪  uniqueb_dic  ∪  Dict("space_group"=>uc.SG, "ibrav"=>cd.ibrav)
    end
end



"""
    config_to_cell_to_QE_pwx_input(config::Dict)

= cell_to_QE_pwx_input(config_to_cell(config), get_sign_of_ibrav(config))

input  : `config` is the configuration, may contain keys "celldm(1)", :cell_parameters or :cif
output :  configuration for the function `pw_input()`

"""
config_to_cell_to_QE_pwx_input(config::Dict) = cell_to_QE_pwx_input(config_to_cell(config), get_sign_of_ibrav(config))
