#NOTE 1: upto rotation,  
#> A1,A2,A3 ---- basis_to_lattice_parameters ---> (a,b,c,alpha,beta,gamma) 
#> A1,A2,A3 <--- lattice_parameters_to_basis ---- (a,b,c,alpha,beta,gamma) 
#> SG>2 and conventions ---> ibrav
#> SG  <--->  (A1,A2,A3) are not 1:1 correspondance, but
#: SG>2 and conventions ---> ibrav ----QE----> (A1,A2,A3)
#> reverse is not true
#: lattp ----> celldm  and  ibrav -----> basis follows DIFFERENT cell conventions
#: lattp <---- basis  <---- celldm  and  ibrav <--\-- basis 
#: cannot get celldm from lattp directly, ibrav is necessary
#! a in lattice parameter is NOT alat !!!!!!! (consider ibrav=2)




function copy(uc::UNITCELL)
    return UNITCELL(uc.UNIT, uc.SG, uc.ALAT, copy(uc.A1_ALAT), copy(uc.A2_ALAT), copy(uc.A3_ALAT))
end



function isapprox(uc1::UNITCELL, uc2::UNITCELL)
    if uc1.UNIT==uc2.UNIT
        return (uc1.SG==uc2.SG && uc1.ALAT≈uc2.ALAT 
             && uc1.A1_ALAT≈uc2.A1_ALAT && uc1.A2_ALAT≈uc2.A2_ALAT && uc1.A3_ALAT≈uc2.A3_ALAT)
    else
        return isapprox(UNITCELL(uc1,uc2.UNIT), uc2)
    end
end



function isequal(uc1::UNITCELL, uc2::UNITCELL)
    if uc1.UNIT==uc2.UNIT
        return (uc1.SG==uc2.SG && uc1.ALAT==uc2.ALAT 
            && uc1.A1_ALAT==uc2.A1_ALAT && uc1.A2_ALAT==uc2.A2_ALAT && uc1.A3_ALAT==uc2.A3_ALAT)
    else
        return isequal(UNITCELL(uc1,uc2.UNIT), uc2)
    end
end


==(uc1::UNITCELL, uc2::UNITCELL) = isequal(uc1, uc2)


function check_valid_cell(c::UNITCELL, ibrv_ngtv_sgn::Bool)
    # try to find another basis set 
    # from SG and ibrv_ngtv_sgn, according to 
    #: SG>2 and conventions ---> ibrav ----QE----> (A1,A2,A3)
    # and compare the three axes
    if c.SG<2
        return true
    else
        A123_Ang   = cell_to_basis_Angstrom(c)
        lattp1     = basis_to_lattice_parameters(A123_Ang, s.SG)
        ibrav2     = Find_ibrav_for_IT(c.SG, ibrv_ngtv_sgn)
        A123_Ang_2 = celldm_to_basis_Angstrom(CELLDM(ibrav2,lattp1)) |> regulate_a1
        A123_Ang_1 = A123_Ang |> regulate_a1
    end
    return all([A123_Ang_1[i] ≈ A123_Ang_2[i] for i=1:3])
end


# convert units
function UNITCELL(uc::UNITCELL, unit::Symbol)
    #@info "[U.1] UNITCELL(uc::UNITCELL=$(uc), unit::Symbol=$(unit))"
    if unit==uc.UNIT
        return copy(uc)
    elseif uc.UNIT==:Angstrom && unit==:Bohr
        return UNITCELL(:Bohr,     uc.SG, uc.ALAT/_BOHR_RADIUS_, uc.A1_ALAT, uc.A2_ALAT, uc.A3_ALAT)
    elseif uc.UNIT==:Bohr && unit==:Angstrom
        return UNITCELL(:Angstrom, uc.SG, uc.ALAT*_BOHR_RADIUS_, uc.A1_ALAT, uc.A2_ALAT, uc.A3_ALAT)
    else
        throw(error("Unknown unit $(uc.UNIT) and $(unit) !!!"))
    end
end



# lattice parameter to cell, via check_adapt_lattice_params_to_QE
# for both conventional and primitive settings
function UNITCELL(lp::LATTPARAM, ibrv_ngtv_sgn::Bool)
    #@info "[U.2] UNITCELL(lp::LATTPARAM=$(lp), ibrv_ngtv_sgn::Bool=$(ibrv_ngtv_sgn))"
    return UNITCELL(CELLDM(lp,ibrv_ngtv_sgn), lp.SG)
end



function UNITCELL(basis3::Vector{Vector{T}}, SG::Int, ibrv_ngtv_sgn::Bool) where {T<:Real}
    #@info "[U.3] UNITCELL(basis3=$(basis3), SG::Int=$(SG), ibrv_ngtv_sgn::Bool=$(ibrv_ngtv_sgn))"
    LP0 = basis_to_lattice_parameters(basis3,SG)
    return UNITCELL(LP0, ibrv_ngtv_sgn)

end



# celldm + SG to cell
# automatically consistent because CELLDM is only available via QE
function UNITCELL(cd::CELLDM, SG::Int)
    #@info "[U.4] UNITCELL(cd::CELLDM=$(cd), SG::Int=$(SG))"
    @assert SG>1         "QE_celldm_ibrav_to_cell() get sg=1, celldm is meaningless in this case !!!"
    @assert cd.ibrav!=0  "QE_celldm_ibrav_to_cell() get ibrav=0, celldm is meaningless in this case !!!"
    @assert consistent(cd, SG)
    (v1,v2,v3) = celldm_to_basis_Bohr(cd)
    return UNITCELL(:Bohr, SG, norm_basis([v1,v2,v3])...)
end


# celldm_dic + ibrav + SG to cell
function UNITCELL(cd_dic_ibrav::Tuple, SG::Int)
    #@info "[U.5] UNITCELL(cd_dic_ibrav::Tuple=$(cd_dic_ibrav), SG::Int=$(SG))"
    (cd_dic::Dict, ibrav::Int) = cd_dic_ibrav
    @assert consistent(ibrav, SG)
    return UNITCELL(CELLDM(cd_dic,ibrav), SG)
end

##* ============================================

function cell_to_basis_Angstrom(c::UNITCELL)
    c_Ang = UNITCELL(c, :Angstrom)
    return  [ c_Ang.ALAT .* c_Ang.A1_ALAT,  
              c_Ang.ALAT .* c_Ang.A2_ALAT,  
              c_Ang.ALAT .* c_Ang.A3_ALAT ]
end


function cell_to_basis_Bohr(c::UNITCELL)
    c_Bohr = UNITCELL(c, :Bohr)
    return  [ c_Bohr.ALAT .* c_Bohr.A1_ALAT,  
              c_Bohr.ALAT .* c_Bohr.A2_ALAT,  
              c_Bohr.ALAT .* c_Bohr.A3_ALAT ]
end
