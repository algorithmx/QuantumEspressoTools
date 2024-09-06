using Test

using LinearAlgebra

include("/home/dabajabaza/jianguoyun/Workspace/QuantumEspressoTools/src/data.jl")

##

function test_compatibility_for_ibrav(ib, lattp)
    for SG=2:230
        if consistent(ib,SG)
            @info "\n$SG\n"
            LP = LATTPARAM(SG, lattp)
            @info "--------------------"
            CD = CELLDM(LP,ib)
            @info "--------------------"
            UC = UNITCELL(CD,SG)
            @info "--------------------"
            t1 = (CD ≈ CELLDM(UC,false))
            @info "--------------------"
            LP1 = LATTPARAM(UC)
            t2 = (LP ≈ LATTPARAM(SG, check_adapt_lattice_params_to_QE_tup6(ib,LP1.p)))
            @info "--------------------"
            if t1 && t2
                nothing
            else
                return false
            end
        end
    end
    return true
end


@testset "crystal structure interpreters" begin
    ## 1. cubic
    lattp = (1.0,1.0,1.0,90.0,90.0,90.0)
    @test test_compatibility_for_ibrav(1, lattp)
    @test test_compatibility_for_ibrav(2, lattp)
    @test test_compatibility_for_ibrav(3, lattp)
    ## 2. Hexagonal and Trigonal P 
    lattp = (1.0,1.0,1.3,90.0,90.0,120.0)
    @test test_compatibility_for_ibrav(4, lattp)
    ## 3. 
    lattp = (1.0,1.0,1.0,88.0,88.0,88.0)
    @test test_compatibility_for_ibrav(5, lattp)
    ## 4. 
    lattp = (1.0,1.0,1.8,90.0,90.0,90.0)
    @test test_compatibility_for_ibrav(6, lattp)
    @test test_compatibility_for_ibrav(7, lattp)
    ## 5. 
    lattp = (1.0,1.2,1.8,90.0,90.0,90.0)
    @test test_compatibility_for_ibrav(8, lattp)
    ## 6. 
    lattp = (1.0,1.2,1.8,90.0,90.0,90.0)
    @test test_compatibility_for_ibrav(8, lattp)
    @test test_compatibility_for_ibrav(9, lattp)
    # @test test_compatibility_for_ibrav(-9, lattp)
    @test test_compatibility_for_ibrav(10, lattp)
    @test test_compatibility_for_ibrav(11, lattp)
    ## 7. 
    lattp = (1.0,1.2,1.8,90.0,90.0,88.0)
    @test test_compatibility_for_ibrav(12, lattp)
    @test test_compatibility_for_ibrav(13, lattp)
    ## 8. 
    lattp = (1.0,1.2,1.8,79.0,57.0,79.0)
    @test test_compatibility_for_ibrav(14, lattp)
end;

##
