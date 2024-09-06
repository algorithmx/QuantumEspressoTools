#TODO inverse of the function check_adapt_lattice_params_to_QE_tup6
#+ restore the conventional setting from lattice parameters following QE convention 
#+ ibrav = 3

"""
    check_adapt_lattice_params_to_QE_tup6(ibrav::Int, latt_params_conv_or_prim::Tuple)

input  : `ibrav`; `latt_params_conv_or_prim` lattice parameters, conventional or primitive 
return : lattice parameters following the QE convention 
note   : 1. the unit of `latt_params_conv_or_prim` is not converted
         2. the criteria for "conventional" or "primitive" in the if-else branch may be 
         INCOMPATIBLE with the input. USE WITH CAUTION !!! 

"""
function check_adapt_lattice_params_to_QE_tup6(ibrav::Int, latt_params_conv_or_prim::Tuple)::Tuple
    #TODO what if I use the output as latt_params_conv_or_prim and run the function again?
    @inline r2d(x) = 180x/pi
    @inline d2r(x) = pi*x/180
    @inline got_inconsistent(lp,ibrav) = "check_adapt_lattice_params_to_QE_tup6() : inconsistent \nlatt_params=$lp and ibrav=$ibrav ."
    (a,b,c,alpha0,beta0,gamma0) = latt_params_conv_or_prim[1:6]
    (alpha, beta, gamma) = d2r.((alpha0,beta0,gamma0))

    if ibrav in [0,14]
        #: the case ibrav=14 is awkward : 
        #: cos(v1,v2) should be cos(ab) but in fact is celldm6[4]=cos(bc)
        return latt_params_conv_or_prim[1:6]
    elseif ibrav==1 # cubic P
        if alpha ≈ beta ≈ gamma ≈ pi/2 && a ≈ b ≈ c
            #: conventional
            return (a, b, c, 90, 90, 90)
        else
            @error  got_inconsistent(latt_params_conv_or_prim,ibrav)
            return (a, b, c, r2d(alpha), r2d(beta), r2d(gamma))
        end
    elseif ibrav==2 # cubic F (fcc)
        if alpha ≈ beta ≈ gamma ≈ pi/2 && a ≈ b ≈ c
            #: conventional
            return (a, b, c, 90, 90, 90)
        elseif alpha ≈ beta ≈ gamma ≈ pi/3 && a ≈ b ≈ c
            #: primitive
            return (sqrt(2)*a, sqrt(2)*a, sqrt(2)*a, 90.0, 90.0, 90.0)
        else
            @error  got_inconsistent(latt_params_conv_or_prim,ibrav)
            return (a, b, c, r2d(alpha), r2d(beta), r2d(gamma))
        end
    elseif ibrav==3 # cubic I (bcc)
        if alpha ≈ beta ≈ gamma ≈ pi/2 && a ≈ b ≈ c
            #: conventional
            return (a, b, c, 90, 90, 90)
        elseif cos(alpha) ≈ -cos(beta) ≈ cos(gamma) ≈ 1/3 && a ≈ b ≈ c
            #: primitive
            COEFF_3 = Float64(2/sqrt(big(3)))
            return (COEFF_3*a, COEFF_3*a, COEFF_3*a, 90, 90, 90)
        else
            @error  got_inconsistent(latt_params_conv_or_prim,ibrav)
            return (a, b, c, r2d(alpha), r2d(beta), r2d(gamma))
        end
    elseif ibrav==-3 # cubic I (bcc), more symmetric axis:
        if alpha ≈ beta ≈ gamma ≈ pi/2 && a ≈ b ≈ c
            #: conventional
            return (a, b, c, 90, 90, 90)
        elseif cos(alpha) ≈ cos(beta) ≈ cos(gamma) ≈ -1/3 && a ≈ b ≈ c
            #: primitive
            COEFF_3 = Float64(2/sqrt(big(3)))
            return (COEFF_3*a, COEFF_3*a, COEFF_3*a, 90, 90, 90)
        else
            @error  got_inconsistent(latt_params_conv_or_prim,ibrav)
            return (a, b, c, r2d(alpha), r2d(beta), r2d(gamma))
        end
    elseif ibrav==4  # Hexagonal and Trigonal P 
        if alpha ≈ beta ≈ pi/2 && gamma ≈ 2pi/3 && a ≈ b
            #: conventional
            return (a, b, c, 90, 90, 120)
        else
            @error  got_inconsistent(latt_params_conv_or_prim,ibrav)
            return (a, b, c, r2d(alpha), r2d(beta), r2d(gamma))
        end
    elseif ibrav in [5,-5]   # Trigonal R, 3fold axis c # Trigonal R, 3fold axis <111> 
        if alpha ≈ beta ≈ gamma && a ≈ b ≈ c
            #: rhombohedral
            return (a, b, c, r2d(alpha), r2d(beta), r2d(gamma))
        elseif alpha ≈ beta ≈ pi/2 && gamma ≈ 2pi/3 && a ≈ b
            #: hexagonal
            cd4 = 1 - 3/(2*(1+(c/a)^2))
            g = acos(cd4)
            ar  = sqrt(a^2+c^2)/3
            return (ar, ar, ar, r2d(g), r2d(g), r2d(g))
        else
            @error  got_inconsistent(latt_params_conv_or_prim,ibrav)
            return (a, b, c, r2d(alpha), r2d(beta), r2d(gamma))
        end
    elseif ibrav==6
        if alpha ≈ beta ≈ gamma ≈ pi/2 && a ≈ b
            #: conventional
            return (a, b, c, 90, 90, 90)
        else
            @error  got_inconsistent(latt_params_conv_or_prim,ibrav)
            return (a, b, c, r2d(alpha), r2d(beta), r2d(gamma))
        end
    elseif ibrav==7 # Tetragonal I (bct)              celldm(3)=c/a
        if alpha ≈ beta ≈ gamma ≈ pi/2 && a ≈ b
            #: conventional
            return (a, b, c, 90, 90, 90)
        elseif beta ≈ gamma && a ≈ b ≈ c && -cos(alpha)+cos(beta)+cos(gamma) ≈ 1.0
            #: QE primitive
            cd3 = sqrt(2cos(beta)/(1-cos(beta)))
            l   = a * sqrt(2-2cos(beta))
            return (l, l, cd3*l, 90, 90, 90)
        else
            @error  got_inconsistent(latt_params_conv_or_prim,ibrav)
            return (a, b, c, r2d(alpha), r2d(beta), r2d(gamma))
        end
    elseif ibrav==8
        if alpha ≈ beta ≈ gamma ≈ pi/2
            #: conventional
            return (a, b, c, 90, 90, 90)
        else
            @error  got_inconsistent(latt_params_conv_or_prim,ibrav)
            return (a, b, c, r2d(alpha), r2d(beta), r2d(gamma))
        end
    elseif ibrav==9 # Orthorhombic base-centered(bco) celldm(2)=b/a
        if alpha ≈ beta ≈ gamma ≈ pi/2
            #: conventional
            return (a, b, c, 90, 90, 90)
        elseif alpha ≈ beta ≈ pi/2 && a ≈ b
            #: QE primitive
            return (a*sqrt(2-2cos(gamma)), a*sqrt(2+2cos(gamma)), c, 90, 90, 90)
        else
            @error  got_inconsistent(latt_params_conv_or_prim,ibrav)
            return (a, b, c, r2d(alpha), r2d(beta), r2d(gamma))
        end
    elseif ibrav==-9 #  as 9, alternate description
        if alpha ≈ beta ≈ gamma ≈ pi/2
            #: conventional
            return (a, b, c, 90, 90, 90)
        elseif alpha ≈ beta ≈ pi/2 && a ≈ b
            #: primitive
            return (a*sqrt(2+2cos(gamma)), a*sqrt(2-2cos(gamma)), c, 90, 90, 90)
        else
            @error  got_inconsistent(latt_params_conv_or_prim,ibrav)
            return (a, b, c, r2d(alpha), r2d(beta), r2d(gamma))
        end
    elseif ibrav==91 # Orthorhombic one-face base-centered A-type
        if alpha ≈ beta ≈ gamma ≈ pi/2
            #: conventional
            return (a, b, c, 90, 90, 90)
        elseif gamma ≈ beta ≈ pi/2 && c ≈ b
            #: primitive
            return (a, c*sqrt(2+2cos(alpha)), c*sqrt(2-2cos(alpha)), 90, 90, 90)
        else
            @error  got_inconsistent(latt_params_conv_or_prim,ibrav)
            return (a, b, c, r2d(alpha), r2d(beta), r2d(gamma))
        end
    elseif ibrav==10 # Orthorhombic face-centered      celldm(2)=b/a
        if alpha ≈ beta ≈ gamma ≈ pi/2
            #: conventional
            return (a, b, c, 90, 90, 90)
        elseif cos(alpha)≈(-a^2+b^2+c^2)/(2b*c) && cos(beta)≈(a^2-b^2+c^2)/(2a*c) && cos(gamma)≈(a^2+b^2-c^2)/(2a*b)
            #: QE primitive
            return (sqrt(2(a^2+b^2-c^2)), sqrt(2(-a^2+b^2+c^2)), sqrt(2(a^2-b^2+c^2)), 90, 90, 90)
        else
            @error  got_inconsistent(latt_params_conv_or_prim,ibrav)
            return (a, b, c, r2d(alpha), r2d(beta), r2d(gamma))
        end
    elseif ibrav==11 # Orthorhombic body-centered      celldm(2)=b/a
        if alpha ≈ beta ≈ gamma ≈ pi/2
            #: conventional
            return (a, b, c, 90, 90, 90)
        elseif a ≈ b ≈ c && cos(alpha)-cos(beta)+cos(gamma) ≈ 1.0
            #: primitive, transformation rules identical to ibrav=7
            l = a
            return (sqrt(2(1 - cos(gamma)))*l,
                    sqrt(2(1 - cos(alpha)))*l,
                    sqrt(2(1 + cos(beta )))*l,  #!!!!!!!!! 
                    90, 90, 90)
        else
            @error  got_inconsistent(latt_params_conv_or_prim,ibrav)
            return (a, b, c, r2d(alpha), r2d(beta), r2d(gamma))
        end
    elseif ibrav==12 # Monoclinic P, unique axis c     celldm(2)=b/a
        if alpha ≈ beta ≈ pi/2
            #: conventional
            return (a, b, c, 90, 90, r2d(gamma))
        else
            @error  got_inconsistent(latt_params_conv_or_prim,ibrav)
            return (a, b, c, r2d(alpha), r2d(beta), r2d(gamma))
        end
    elseif ibrav==-12 # Monoclinic P, unique axis b     celldm(2)=b/a
        if alpha ≈ gamma ≈ pi/2
            #: conventional
            return (a, b, c, 90, r2d(beta), 90)
        else
            @error  got_inconsistent(latt_params_conv_or_prim,ibrav)
            return (a, b, c, r2d(alpha), r2d(beta), r2d(gamma))
        end
    elseif ibrav==13 #  13          Monoclinic base-centered        celldm(2)=b/a
        #                           (unique axis c)                 celldm(3)=c/a,
        #                                                           celldm(4)=cos(gamma)
        if alpha ≈ beta ≈ pi/2
            #: conventional
            return (a, b, c, 90, 90, r2d(gamma))
        elseif  alpha ≈ gamma  &&  a ≈ c
            #: primitive
            return ( a*sqrt(2(1+cos(beta))),  b,  a*sqrt(2(1-cos(beta))), 90, 90, r2d(acos(sqrt(2)*cos(gamma)/sqrt(1+cos(beta)))) )
        else
            @error  got_inconsistent(latt_params_conv_or_prim,ibrav)
            return (a, b, c, r2d(alpha), r2d(beta), r2d(gamma))
        end
    elseif ibrav==-13 # Monoclinic base-centered        celldm(2)=b/a
        #                           (unique axis b)     celldm(3)=c/a,
        #                                               celldm(5)=cos(beta)
        if gamma ≈ alpha ≈ pi/2
            #: conventional
            return (a, b, c, 90, r2d(beta), 90)
        elseif  alpha + beta ≈ pi  &&  a ≈ b
            #: primitive
            return ( a*sqrt(2(1-cos(gamma))), a*sqrt(2(1+cos(gamma))), c, 90, r2d(acos(sqrt(2)*cos(beta)/sqrt(1-cos(gamma)))), 90 )
        else
            @error  got_inconsistent(latt_params_conv_or_prim,ibrav)
            return (a, b, c, r2d(alpha), r2d(beta), r2d(gamma))
        end
    else
        @error  "check_adapt_lattice_params_to_QE_tup6() got unsupported ibrav=$ibrav ."
        return (a, b, c, r2d(alpha), r2d(beta), r2d(gamma))
    end
end


function check_adapt_lattice_params_to_QE(lp::LATTPARAM, ibrv_ngtv_sgn::Bool)::LATTPARAM
    SG      = lp.SG
    ibrav   = Find_ibrav_for_IT(SG, ibrv_ngtv_sgn)
    lattp6  = check_adapt_lattice_params_to_QE_tup6(ibrav, lp.p)
    return LATTPARAM(SG, lattp6)
end