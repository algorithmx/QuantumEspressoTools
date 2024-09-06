function check_valid_config(config::Dict)
    @inline no_celldm_in_config(c) = all(["celldm($i)" ∉ keys(c) for i=1:6])
    # ibrav v.s :cell_parameters
    if "ibrav" ∈ keys(config)
        if :cell_parameters ∈ keys(config)
            @assert  no_celldm_in_config(config)  "check_valid_config(): [1.1] both keys :cell_parameters and celldm(1) appear in config !!!!!"
            @assert  config["ibrav"]==0  "check_valid_config(): [1.2] config has :cell_parameters but ibrav!=0  !!!!!"
        elseif  "celldm(1)" ∈ keys(config)
            @assert  config["ibrav"]> 0  "check_valid_config(): [1.3] config doesn't have :cell_parameters but ibrav=0  !!!!!"
            @assert  is_valid(CELLDM(config,config["ibrav"]))  "check_valid_config(): [1.4] incompatible ibrav anc celldm !!!!!"
            if "space_group" ∈ keys(config)
                @assert consistent(config["ibrav"], config["space_group"])  "check_valid_config(): [1.5] incompatible config[\"ibrav\"] and config[\"space_group\"] !!!!!"
            end
        elseif          :cif ∈ keys(config)
            LP = LATTPARAM(config[:cif])
            if "space_group" ∈ keys(config)
                @assert  config["space_group"]==LP.SG  "check_valid_config(): [1.6] incompatible config[\"space_group\"] and LP.SG !!!!!"
            end
            @assert  consistent(config["ibrav"],LP.SG)  "check_valid_config(): [1.7] incompatible LP and config[\"ibrav\"] !!!!!"
            @assert  is_valid(CELLDM(LP,config["ibrav"]))  "check_valid_config(): [1.8] config[:cif] gives invalid celldm !!!!!"
        else
            throw(error("check_valid_config(): [1.9] no :cell_parameters no celldm(1) !!!!!"))
        end

        if "uniqueb" ∈ keys(config)
            @assert  (config["ibrav"]<0) == config["uniqueb"]   "check_valid_config(): [1.10] sgn(ibrav) and uniqueb inconsistent !!!!!"
        elseif :ibrv_ngtv_sgn ∈ keys(config)
            @assert  (config["ibrav"]<0) == config[:ibrv_ngtv_sgn]  "check_valid_config(): [1.11] sgn(ibrav) and :ibrv_ngtv_sgn inconsistent !!!!!"
        end
    else  # "ibrav" ∉ keys(config)
        if  :cell_parameters ∈ keys(config)
            @assert no_celldm_in_config(config)  "check_valid_config(): [2.1] both keys :cell_parameters and celldm(1) appear in config !!!!!"
        elseif         :cif  ∈ keys(config)
            LP = LATTPARAM(config[:cif])
            #: no ibrav, no :cell_parameters, no IT_num, only (a,b,c, alpha,beta,gamma)
            @assert  (length(config[:cif])==7 && config[:cif][7]>1) "check_valid_config(): [2.2] key :cell_parameters, ibrav are absent in config, and\nconfig[:cif] does not tell SG.\nCan't determine anything !!!!!"
            @assert  is_valid(CELLDM(LP, get(config,"uniqueb",false)))  "check_valid_config(): [2.3] config[:cif] gives invalid celldm !!!!!"
        else
            throw(error("check_valid_config(): [2.4] no :cell_parameters no :cif no celldm(1) !!!!!"))
        end

        @assert :ibrv_ngtv_sgn ∉ keys(config)  "check_valid_config(): [2.5] meaningless :ibrv_ngtv_sgn in keys of config !!!!!"
        @assert "uniqueb"      ∉ keys(config)  "check_valid_config(): [2.6] meaningless uniqueb in keys of config !!!!!"
    end

    #* ALLES GUT !!!
    return true
end

