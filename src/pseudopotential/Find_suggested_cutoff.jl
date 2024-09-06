function Find_suggested_cutoff_KW(pseudo_lines, wcut_kw, ecut_kw)
    @inline extract_number(res) = match(r"\d+\.(\d+([EeDd]\+\d+)?)?",res).match

    ecut_wavefunctions  = 100.0
    p_wavefunctions  = findfirst(x->occursin(wcut_kw,x), pseudo_lines)
    if p_wavefunctions!==nothing
        a  = extract_number(pseudo_lines[p_wavefunctions])  |> strtonum_fort
        if a > 1.0  ecut_wavefunctions = a end
    end
    
    ecut_charge_density = 400.0
    p_charge_density = findfirst(x->occursin(ecut_kw,x), pseudo_lines)
    if p_charge_density!==nothing
        a = extract_number(pseudo_lines[p_charge_density]) |> strtonum
        if a > 1.0  ecut_charge_density = a end
    end

    wfc_cut = max(ecut_wavefunctions, 0.25*ecut_charge_density)
    return wfc_cut, 4wfc_cut
end


Find_suggested_cutoff_PSL(pseudo_lines) = Find_suggested_cutoff_KW(
                                            pseudo_lines, 
                                            "minimum cutoff for wavefunctions",
                                            "minimum cutoff for charge density")


Find_suggested_cutoff_XXX(pseudo_lines) = Find_suggested_cutoff_KW(
                                            pseudo_lines, 
                                            "wfc_cutoff",
                                            "rho_cutoff")


function Find_suggested_cutoff_SG15(pseudo_lines)
    return 100.0, 400.0
end


function Find_suggested_cutoff_GBRV(pseudo_lines)
    # http://www.physics.rutgers.edu/gbrv/
    # GBRV v1.5 + Elemental Testing
    # http://www.physics.rutgers.edu/gbrv/delta_update_v15.pdf
    # GBRV phonon update and JTH/PSlibrary testing
    # http://www.physics.rutgers.edu/gbrv/gbrv_phonon_update2.pdf
    return 40.0, 200.0
end


function Find_suggested_cutoff_GHH(pseudo_lines)
    return 100.0, 400.0
end



function Find_suggested_cutoff_MT(pseudo_lines)
    return 100.0, 400.0
end


function Find_suggested_cutoff_other(pseudo_lines)
    return 100.0, 400.0
end


function Find_suggested_cutoff(psmode,pseudo_fn)
    pseudo_lines = readlines(pseudo_fn)
    if startswith(psmode, "PSL")
        return Find_suggested_cutoff_PSL(pseudo_lines)
    elseif startswith(psmode, "GBRV")
        return Find_suggested_cutoff_GBRV(pseudo_lines)
    elseif startswith(psmode, "GHH")
        return Find_suggested_cutoff_GHH(pseudo_lines)
    elseif startswith(psmode, "MT")
        return Find_suggested_cutoff_MT(pseudo_lines)
    elseif startswith(psmode, "SG15")
        return Find_suggested_cutoff_SG15(pseudo_lines)
    else
        return Find_suggested_cutoff_other(pseudo_lines)
    end
end
