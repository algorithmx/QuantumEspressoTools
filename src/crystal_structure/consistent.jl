#> the case for ibrav=0 or space_group=1 is always valid !
consistent(ib::Int,SG::Int,u::Bool) = (ib==0 || SG==1 || (ib!=0 && SG>1 && Find_ibrav_for_IT(SG,u)==ib))

consistent(ib::Int,SG::Int) = consistent(ib,SG,true) || consistent(ib,SG,false)

consistent(cd::CELLDM, SG::Int) = consistent(cd.ibrav, SG)

consistent(cd::CELLDM, lp::LATTPARAM, ibrv_ngtv_sgn::Bool) = (consistent(cd.ibrav, lp.SG, ibrv_ngtv_sgn) && UNITCELL(cd,lp.SG) ≈ UNITCELL(lp,ibrv_ngtv_sgn))

consistent(lp::LATTPARAM, cd::CELLDM, ibrv_ngtv_sgn::Bool) = consistent(cd, lp, ibrv_ngtv_sgn)

