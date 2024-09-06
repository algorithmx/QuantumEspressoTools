function LATTPARAM(cif::Tuple)
    #@info "[L.1] LATTPARAM(cif::Tuple=$(cif))"
    LATTPARAM((length(cif)==7 ? cif[7] : 1), cif[1:6])
end


function LATTPARAM(c::UNITCELL)
    #@info "[L.2] LATTPARAM(c::UNITCELL=$(c))"
    basis_to_lattice_parameters(cell_to_basis_Angstrom(c), c.SG)
end


function copy(lp::LATTPARAM)
    return LATTPARAM(lp.SG, copy(lp.p))
end


isequal(lp1::LATTPARAM, lp2::LATTPARAM) = (lp1.SG==lp2.SG && lp1.p==lp2.p)

isapprox(lp1::LATTPARAM, lp2::LATTPARAM) = (lp1.SG==lp2.SG  &&  lp1.p≈lp2.p)

==(lp1::LATTPARAM, lp2::LATTPARAM) = isequal(lp1, lp2)


function basis_to_lattice_parameters(basis_Ang_final, sg=1)::LATTPARAM
    basis = basis_Ang_final
    (a,b,c) = (norm(basis[1]), norm(basis[2]), norm(basis[3]))
    d(x,y) = 180acos(dot(normalize(basis[x]),normalize(basis[y])))/π
    (ang_bc, ang_ac, ang_ab) = (d(2,3), d(1,3), d(1,2))
    return LATTPARAM(sg, (a, b, c, ang_bc, ang_ac, ang_ab))
end


function lattice_parameters_to_basis(lp::LATTPARAM, ibrv_ngtv_sgn::Bool=false)
    if lp.SG==1
        (a_Ang, b_Ang, c_Ang, alpha, beta, gamma) = lp.p
        alphar, betar, gammar = map(x->Float64((π/big(180))*x), [alpha, beta, gamma])
        cell_px = a_Ang.*[1.0,         0.0,                                              0.0 ]
        cell_py = b_Ang.*[cos(gammar), sin(gammar),                                      0.0 ]
        cell_pz = c_Ang.*[cos(betar),  (cos(alphar)-cos(betar)*cos(gammar))/sin(gammar), 
                            sqrt(1.0 - cos(alphar)^2 - cos(betar)^2 - cos(gammar)^2 + 2*cos(alphar)*cos(betar)*cos(gammar))/sin(gammar)]
        return [cell_px, cell_py, cell_pz]
    else
        return celldm_to_basis_Angstrom(CELLDM(lp, ibrv_ngtv_sgn))
    end    
end


lattice_parameters_to_basis((a,b,c,al,be,ga)) = lattice_parameters_to_basis(LATTPARAM(1,(a,b,c,al,be,ga)))


function config_to_lattp(config::Dict)::LATTPARAM
    @assert  check_valid_config(config)
    @assert :cif ∈ keys(config)
    LP0 = ((config[:cif] isa LATTPARAM) ? config[:cif] : LATTPARAM(config[:cif]))
    neg_sgn = get_sign_of_ibrav(config)
    return check_adapt_lattice_params_to_QE(LP0, neg_sgn)
end


function get_sign_of_ibrav(config::Dict)
    if "ibrav" ∈ keys(config)
        if "uniqueb" ∈ keys(config)
            @assert  (config["ibrav"]<0) == config["uniqueb"]
        elseif :ibrv_ngtv_sgn ∈ keys(config)
            @assert  (config["ibrav"]<0) == config[:ibrv_ngtv_sgn]
        else
            nothing
        end
        return (config["ibrav"]<0)
    else
        @assert :ibrv_ngtv_sgn ∉ keys(config)
        @assert "uniqueb"      ∉ keys(config)
        return false
    end
end


