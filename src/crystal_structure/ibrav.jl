#> for space group 1 P1, 
#> we never specify space_group in the input of pw.x
#> and never use crystal_sg
#> and never use celldm

function Find_ibrav_for_IT(
    IT_num::Int,
    ibrv_ngtv_sgn::Bool
    )

    sgn = (ibrv_ngtv_sgn ? -1 : 1)

    @assert IT_num>=1 && IT_num<=230

    if IT_num==1  return 0  end

    spacegroup   = Int_Tables[IT_num]
    primitive    = startswith(spacegroup,"P")
    facecentered = startswith(spacegroup,"F")
    bodycentered = startswith(spacegroup,"I")
    basecentered = startswith(spacegroup,"C") ||  startswith(spacegroup,"A")
    rhombo       = startswith(spacegroup,"R")

    ibrav = 0
    class = CRYSTAL_CLASS(IT_num)
    if class=="cubic"
        if (primitive)    ibrav=1  end
        if (facecentered) ibrav=2  end
        if (bodycentered) ibrav=3sgn  end  #TODO -3 ?
    elseif class=="hexagonal"
        if (primitive)
            ibrav=4
        else
            @warn "Find_ibrav_for_IT($spacegroup,$ibrv_ngtv_sgn):\nNo P no R."
            ibrav=0
        end
    elseif class=="trigonal"
        if (rhombo) 
            ibrav=5sgn
        elseif (primitive)
            ibrav=4
        else
            @warn "Find_ibrav_for_IT($spacegroup,$ibrv_ngtv_sgn):\nNo P no R."
            ibrav=0
        end
    elseif class=="tetragonal"
        if (primitive)    ibrav=6  end
        if (bodycentered) ibrav=7  end
    elseif class=="orthorhombic"
        if (primitive)    ibrav=8  end
        if (basecentered) ibrav=9sgn  end #TODO 91
        if (facecentered) ibrav=10 end
        if (bodycentered) ibrav=11 end
    elseif class=="monoclinic"
        #! unique axis b
        if (primitive)
            ibrav = 12sgn
        elseif (basecentered)
            ibrav = 13sgn
        end
    elseif class=="triclinic"
        ibrav=14
    else
        ibrav=0
    end
    return ibrav
end

## Find_ibrav_for_IT.(1:230, false) |> print

global const SG2ibrav = [
    0, 14, 12, 12, 13, 12, 12, 13, 13, 12, 12, 13, 12, 12, 13, 
    8, 8, 8, 8, 9, 9, 
    10, 11, 11, 
    8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 9, 9, 9, 9, 9, 9, 9, 10, 10, 
    11, 11, 11, 
    8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 9, 9, 9, 9, 9, 9, 10, 10, 
    11, 11, 11, 11, 
    6, 6, 6, 6, 7, 7, 6, 7, 6, 6, 6, 6, 7, 7, 6, 6, 6, 6, 
    6, 6, 6, 6, 7, 7, 6, 6, 6, 6, 6, 6, 6, 6, 7, 
    7, 7, 7, 6, 6, 6, 6, 6, 6, 6, 6, 7, 7, 7, 
    7, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 
    6, 6, 6, 7, 7, 7, 7, 
    4, 4, 4, 5, 4, 5, 4, 4, 4, 4, 4, 4, 5, 4, 
    4, 4, 4, 5, 5, 4, 4, 4, 4, 5, 5, 4, 4, 4, 
    4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 
    4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 1, 2, 3, 1, 
    3, 1, 1, 2, 2, 3, 1, 3, 1, 1, 2, 2, 3, 1, 
    1, 3, 1, 2, 3, 1, 2, 3, 1, 1, 1, 1, 2, 2, 
    2, 2, 3, 3
]

