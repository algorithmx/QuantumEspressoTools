function copy(cd::CELLDM)
    return CELLDM(cd.ibrav, copy(cd.v))
end


function isapprox(cd1::CELLDM, cd2::CELLDM)
    return cd1.ibrav==cd2.ibrav  &&  cd1.v ≈ cd2.v
end


function isequal(cd1::CELLDM, cd2::CELLDM)
    return cd1.ibrav==cd2.ibrav  &&  cd1.v==cd2.v
end


==(cd1::CELLDM, cd2::CELLDM) = isequal(cd1, cd2)

function CELLDM(c::UNITCELL, ibrv_ngtv_sgn::Bool)
    @info "[C.1] CELLDM(c::UNITCELL=$(c), ibrv_ngtv_sgn::Bool=$(ibrv_ngtv_sgn))"
    CELLDM(LATTPARAM(c), Find_ibrav_for_IT(c.SG, ibrv_ngtv_sgn))
end


# ibrav + celldm_dic
function CELLDM(celldm_dic::Dict, ibrav::Int)
    @info "[C.2] CELLDM(celldm_dic::Dict=$(celldm_dic), ibrav::Int=$(ibrav))"
    if ibrav==0
        return CELLDM(0, zeros(Float64,6))
    else
        return CELLDM(ibrav, [get(celldm_dic,"celldm($k)",0.0) for k=1:6])
    end
end


# ibrav + LATTPARAM
function CELLDM(lp_conv::LATTPARAM, ibrav::Int)
    @info "[C.3] CELLDM(lp_conv::LATTPARAM=$(lp_conv), ibrav::Int=$(ibrav))"
    @assert consistent(ibrav, lp_conv.SG)
    lattp_adatped = check_adapt_lattice_params_to_QE_tup6(ibrav, lp_conv.p)
    return CELLDM(ibrav, celldm_vec(ibrav, (lattp_adatped..., lp_conv.SG)))
end


# ibrav + LATTPARAM
function CELLDM(lp::LATTPARAM, ibrv_ngtv_sgn::Bool)
    @info "[C.4] CELLDM(lp::LATTPARAM=$(lp), ibrv_ngtv_sgn::Bool=$(ibrv_ngtv_sgn))"
    return CELLDM(lp, Find_ibrav_for_IT(lp.SG, ibrv_ngtv_sgn))
end


#* ====================================================


function is_valid(cd::CELLDM)
    ibrav, celldm6 = (cd.ibrav, cd.v)
    @inline in_dic(iii::Vector{Int}) = all([ ( ((k in iii) && abs(celldm6[k])>1e-8) 
                                            || (!(k in iii) && !(abs(celldm6[k])>1e-8)) )
                                            for k=1:6 ])
    # https://www.quantum-espresso.org/Doc/INPUT_PW.html#idm218
    if ibrav in [4,6,7]
        return in_dic([1,3])
        #return (@sprintf  "                   celldm(1) = %19.14f, celldm(3) = %19.14f"  (a_Ang/_BOHR_RADIUS_)  (c_Ang/a_Ang))
    elseif ibrav in [5,-5]
        return in_dic([1,4])
        #return (@sprintf  "                   celldm(1) = %19.14f, celldm(4) = %19.14f"  (a_Ang/_BOHR_RADIUS_)  cos(alphar))
    elseif ibrav in [14,]
        return in_dic([1,2,3,4,5,6])
        #return (@sprintf  "                   celldm(1) = %19.14f, celldm(2) = %19.14f, celldm(3) = %19.14f, celldm(4) = %19.14f, celldm(5) = %19.14f, celldm(6) = %19.14f" (a_Ang/_BOHR_RADIUS_) (b_Ang/a_Ang) (c_Ang/a_Ang) cos(alphar) cos(betar) cos(gammar))
    elseif ibrav in [-12,-13]
        return in_dic([1,2,3,5])
        #return (@sprintf  "                   celldm(1) = %19.14f, celldm(2) = %19.14f, celldm(3) = %19.14f, celldm(5) = %19.14f" (a_Ang/_BOHR_RADIUS_) (b_Ang/a_Ang) (c_Ang/a_Ang) cos(betar))
    elseif ibrav in [12,13]
        return in_dic([1,2,3,4])
        #return (@sprintf  "                   celldm(1) = %19.14f, celldm(2) = %19.14f, celldm(3) = %19.14f, celldm(4) = %19.14f" (a_Ang/_BOHR_RADIUS_) (b_Ang/a_Ang) (c_Ang/a_Ang) cos(gammar))
    elseif ibrav in [8,9,-9,91,10,11]
        return in_dic([1,2,3])
        #return (@sprintf  "                   celldm(1) = %19.14f, celldm(2) = %19.14f, celldm(3) = %19.14f" (a_Ang/_BOHR_RADIUS_) (b_Ang/a_Ang) (c_Ang/a_Ang))
    elseif ibrav in [1,2,3,-3]
        return in_dic([1])
        #return (@sprintf  "                   celldm(1) = %19.14f" (a_Ang/_BOHR_RADIUS_))
    elseif ibrav in [0,]  #! use cell_parameter card
        @warn  "is_valid(cd::CELLDM) cd.ibrav=0, return true !!!!!"
        return true
    end
end


"""
    celldm_vec(ibrav::Int, latt_params_conv_Ang::Tuple)

input  : ibrav parameter (with sign) and 
         lattice parameter (conventional, unit Angstrom) 
return : a vector [celldm(1), ..., celldm(6)]
note   : alat (Bohr) = celldm(1) is the lattice parameter "a (Ang) / a0"

Implements https://www.quantum-espresso.org/Doc/INPUT_PW.html#idm218

"""
function celldm_vec(ibrav::Int, latt_params_conv_Ang::Tuple)::Vector{Float64}
    @inline d2r(x) = pi*x/180
    (a_Ang,b_Ang,c_Ang, alpha0,beta0,gamma0) = latt_params_conv_Ang[1:6]
    (alphar,betar,gammar) = d2r.((alpha0,beta0,gamma0))  #!!!!!!!!!
    c1 = a_Ang/_BOHR_RADIUS_
    # https://www.quantum-espresso.org/Doc/INPUT_PW.html#idm218
    if ibrav in [4,6,7]
        c3 = c_Ang/a_Ang
        return [c1, 0, c3, 0, 0, 0]
    elseif ibrav in [5,-5]
        c4 = cos(gammar)
        return [c1, 0, 0, c4, 0, 0]
    elseif ibrav in [14,]
        c2 = (b_Ang/a_Ang)
        c3 = (c_Ang/a_Ang)
        c4 = cos(alphar)
        c5 = cos(betar)
        c6 = cos(gammar)
        return [c1, c2, c3, c4, c5, c6]
    elseif ibrav in [-12,-13]
        c2 = (b_Ang/a_Ang)
        c3 = (c_Ang/a_Ang)
        c5 = cos(betar)
        return [c1,c2,c3,0,c5,0]
    elseif ibrav in [12,13]
        c2 = b_Ang/a_Ang
        c3 = c_Ang/a_Ang 
        c4 = cos(gammar)
        return [c1,c2,c3,c4,0,0]
    elseif ibrav in [8,9,-9,91,10,11]
        c2 = b_Ang/a_Ang
        c3 = c_Ang/a_Ang
        return [c1,c2,c3,0,0,0]
    elseif ibrav in [1,2,3,-3]
        return [c1,0,0,0,0,0]
    elseif ibrav in [0,]
        @warn  "celldm_vec() ibrav=0, return [0,0,0,0,0,0] !!!!!"
        return [0,0,0,0,0,0]
    end
end

celldm_vec(cd::CELLDM) = (cd.ibrav>0 ? cd.v : zeros(Float64,6))

celldm_dic(cd::CELLDM) = Dict{String,Float64}("celldm($i)"=>f for (i,f) in enumerate(cd.v) if (abs(f)>1e-8 && cd.ibrav>0))

celldm_dic(ibrav::Int, latt_params::Tuple) = Dict{String,Float64}("celldm($i)"=>f for (i,f) in enumerate(celldm_vec(ibrav,latt_params)))


function celldm_to_basis_Bohr(cd::CELLDM)::Vector{Vector{Float64}}
    celldm6, ibrav = (cd.v, cd.ibrav)
    v1 = zeros(3)
    v2 = zeros(3)
    v3 = zeros(3)
    if ibrav==0 # free
        return [v1,v2,v3]
    elseif ibrav==1 # cubic P (sc)
        v1 = [1,0,0]
        v2 = [0,1,0]
        v3 = [0,0,1]
    elseif ibrav==2 # cubic F (fcc)
        v1 = (1/2).*[-1,0,1]
        v2 = (1/2).*[0,1,1]
        v3 = (1/2).*[-1,1,0]
    elseif ibrav==3 # cubic I (bcc)
        v1 = (1/2).*[1,1,1]
        v2 = (1/2).*[-1,1,1]
        v3 = (1/2).*[-1,-1,1]
    elseif ibrav==-3 # cubic I (bcc), more symmetric axis:
        v1 = (1/2).*[-1,1,1]
        v2 = (1/2).*[1,-1,1]
        v3 = (1/2).*[1,1,-1]
    elseif ibrav==4 # Hexagonal and Trigonal P        celldm(3)=c/a
        v1 =  [1,0,0]
        v2 =  [-1/2,sqrt(3)/2,0]
        v3 =  [0,0,celldm6[3]]
    elseif ibrav==5 # Trigonal R, 3fold axis c        celldm(4)=cos(gamma)
        #The crystallographic vectors form a three-fold star around
        #the z-axis, the primitive cell is a simple rhombohedron:
        tx = sqrt((1-celldm6[4])/2)
        ty = sqrt((1-celldm6[4])/6)
        tz = sqrt((1+2celldm6[4])/3)
        v1 = [tx,-ty,tz]
        v2 = [0,2ty,tz]
        v3 = [-tx,-ty,tz]
        #where c=cos(gamma) is the cosine of the angle gamma between
        #any pair of crystallographic vectors
    elseif ibrav == -5 #  Trigonal R, 3fold axis <111>    celldm(4)=cos(gamma)
        # tx, ty, tz as for case ibrav=5
        tx = sqrt((1-celldm6[4])/2)
        ty = sqrt((1-celldm6[4])/6)
        tz = sqrt((1+2celldm6[4])/3)
        #The crystallographic vectors form a three-fold star around
        #<111>. Defining a' = a/sqrt(3) :
        u  = tz - 2*sqrt(2)*ty
        v  = tz + sqrt(2)*ty
        ap = 1/sqrt(3)
        v1 = ap .* [u,v,v]
        v2 = ap .* [v,u,v]
        v3 = ap .* [v,v,u]
        # Note: if you prefer x,y,z as axis in the cubic limit,
        #  set  u = tz + 2*sqrt(2)*ty,  v = tz - sqrt(2)*ty
        #  See also the note in Modules/latgen.f90
    elseif ibrav==6 # Tetragonal P (st)               celldm(3)=c/a
        v1 = [1,0,0]
        v2 = [0,1,0]
        v3 = [0,0,celldm6[3]]
    elseif ibrav==7 # Tetragonal I (bct)              celldm(3)=c/a
        v1=(1/2).*[ 1,-1,celldm6[3]]
        v2=(1/2).*[ 1, 1,celldm6[3]]
        v3=(1/2).*[-1,-1,celldm6[3]]
    elseif ibrav==8 # Orthorhombic P                  celldm(2)=b/a
        #                                            celldm(3)=c/a
        v1 = [1,0,0]
        v2 = [0,celldm6[2],0]
        v3 = [0,0,celldm6[3]]
    elseif ibrav==9 # Orthorhombic base-centered(bco) celldm(2)=b/a
        #                                             celldm(3)=c/a
        v1 = [1/2, celldm6[2]/2,0]
        v2 = [-1/2,celldm6[2]/2,0]
        v3 = [0,0,celldm6[3]]
    elseif ibrav==-9 #  as 9, alternate description
        v1 = [1/2, -celldm6[2]/2,0]
        v2 = [1/2,  celldm6[2]/2,0]
        v3 = [0,0,celldm6[3]]
    elseif ibrav==91 # Orthorhombic one-face base-centered A-type
        #                                           celldm(2)=b/a
        #                                           celldm(3)=c/a
        v1 = [1, 0, 0]
        v2 = [0,celldm6[2]/2,-celldm6[3]/2]
        v3 = [0,celldm6[2]/2, celldm6[3]/2]
    elseif ibrav==10 # Orthorhombic face-centered      celldm(2)=b/a
        #                                              celldm(3)=c/a
        v1 = [1/2, 0, celldm6[3]/2]
        v2 = [1/2, celldm6[2]/2, 0]
        v3 = [0, celldm6[2]/2, celldm6[3]/2]
    elseif ibrav==11 # Orthorhombic body-centered      celldm(2)=b/a
        #                                           celldm(3)=c/a
        v1=[ 1/2, celldm6[2]/2, celldm6[3]/2]
        v2=[-1/2, celldm6[2]/2, celldm6[3]/2]
        v3=[-1/2,-celldm6[2]/2, celldm6[3]/2]
    elseif ibrav==12 # Monoclinic P, unique axis c     celldm(2)=b/a
        #                                           celldm(3)=c/a,
        #                                           celldm(4)=cos(ab)
        v1 = [1,0,0]
        v2 = [celldm6[2]*celldm6[4], celldm6[2]*sqrt(1-celldm6[4]^2), 0]
        v3 = [0,0,celldm6[3]]
        # where gamma is the angle between axis a and b.
    elseif ibrav==-12 # Monoclinic P, unique axis b     celldm(2)=b/a
        #                                               celldm(3)=c/a,
        #                                               celldm(5)=cos(ac)
        v1 = [1,0,0]
        v2 = [0,celldm6[2],0]
        v3 = celldm6[3].*[celldm6[5], 0, sqrt(1-celldm6[5]^2)]
        # where beta is the angle between axis a and c
    elseif ibrav==13 # Monoclinic base-centered        celldm(2)=b/a
        # (unique axis c)                 celldm(3)=c/a,
        #                                celldm(4)=cos(gamma)
        v1 = [  1/2,         0,          -celldm6[3]/2]
        v2 = [celldm6[2]*celldm6[4],celldm6[2]*sqrt(1-celldm6[4]^2),0]
        v3 = [  1/2,         0,           celldm6[3]/2]
        # where gamma=angle between axis a and b projected on xy plane
    elseif ibrav==-13 # Monoclinic base-centered        celldm(2)=b/a
        # (unique axis b)                 celldm(3)=c/a,
        #                                 celldm(5)=cos(beta)
        v1 = [  1/2, celldm6[2]/2, 0]
        v2 = [ -1/2, celldm6[2]/2, 0]
        v3 = [celldm6[3]*celldm6[5], 0, celldm6[3]*sqrt(1-celldm6[5]^2)]
        #where beta=angle between axis a and c projected on xz plane
        #! IMPORTANT NOTICE: until QE v.6.4.1, axis for ibrav=-13 had a
        #! different definition: v1(old) =-v2(now), v2(old) = v1(now)
    elseif ibrav==14 # Triclinic                    celldm(2)= b/a,
        #                                           celldm(3)= c/a,
        #                                           celldm(4)= cos(bc),
        #                                           celldm(5)= cos(ac),
        #                                           celldm(6)= cos(ab)
        #: this case is awkward : 
        #: cos(v1,v2) should be cos(ab) but in fact is celldm6[4]=cos(bc)
        v1 = [1, 0, 0]
        v2 = [celldm6[2]*celldm6[4],celldm6[2]*sqrt(1-celldm6[4]^2),0]
        v3 = [celldm6[3]*celldm6[5],  celldm6[3]*(celldm6[4]-celldm6[5]*celldm6[6])/sqrt(1-celldm6[6]^2), 
                        ( celldm6[3]*sqrt( 1 + 2*celldm6[4]*celldm6[5]*celldm6[6]
                        - celldm6[4]^2 - celldm6[5]^2 - celldm6[6]^2 ) / sqrt(1-celldm6[6]^2)) ]
        #where alpha is the angle between axis b and c
        #       beta is the angle between axis a and c
        #    gamma is the angle between axis a and b
    end
    return [celldm6[1].*v1,  celldm6[1].*v2,  celldm6[1].*v3]
end

celldm_to_basis_Angstrom(cd::CELLDM) = _BOHR_RADIUS_ .* celldm_to_basis_Bohr(cd)



function config_to_celldm(config::Dict)::CELLDM
    @assert check_valid_config(config)
    @assert "ibrav"     ∈ keys(config)
    @assert "celldm(1)" ∈ keys(config)
    SG = get(config,"space_group",1)
    return UNITCELL((config, config["ibrav"]), SG)
end