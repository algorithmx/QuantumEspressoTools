using LinearAlgebra
using Isosuite
include("../src/crystallography.jl")
include("../src/check_adapt_lattice_params_to_QE.jl")
include("../src/lattice_parameters.jl")
include("../src/celldm.jl")
include("../src/data.jl")
include("../src/ibrav.jl")
include("../src/crystal.jl")

##

function test_crystal(a1, a2, a3, SG::Int, uniqueb::Bool)
    UC0 = UNITCELL(:Angstrom, SG, norm_basis([a1,a2,a3])...)
    IB  = Find_ibrav_for_IT(SG,uniqueb)
    @show SG,uniqueb,IB
    LP0 = check_adapt_lattice_params_to_QE(basis_to_lattice_parameters([a1,a2,a3],SG), uniqueb)
    CD0 = CELLDM(IB, celldm_vec(IB, LP0.p[1:6]))

    UC1 = UNITCELL([a1, a2, a3], SG, uniqueb)
    UC2 = UNITCELL(LP0, uniqueb)
    UC3 = UNITCELL(CD0, SG)

    t1  = UC0 ≈ UC1
    t2  = UC0 ≈ UC2
    t3  = UC0 ≈ UC3
    t12  = UC1 ≈ UC2
    t23  = UC2 ≈ UC3
    t31  = UC3 ≈ UC1

    if ! all((t1,t2,t3,t12,t23,t31))
        @warn  "component"
        @show (t1,t2,t3,t12,t23,t31)
        @show (a1,a2,a3)
        @show UC0
        @show UC1
        @show UC2
        @show UC3
    end
    return all((t1,t2,t3,t12,t23,t31))
end

##

for SG = 2:230
    sleep(0.2)
    println()
    println("------------------")
    println("SG=$SG")
    @info "SG=$SG"
    sleep(0.2)

    ibrav0 = Find_ibrav_for_IT(SG,false)
    r6 = rand(6)
    if ibrav0 ∈ [12,13]
        r6[5]=0
        r6[6]=0
    end
    CD_rand = CELLDM(ibrav0, r6)
    (a1,a2,a3) = celldm_to_basis_Angstrom(CD_rand)
    test_crystal(a1, a2, a3, SG, false)

    ibrav1 = Find_ibrav_for_IT(SG,true)
    if ibrav1 == -ibrav0
        if ibrav1 ∈ [-12,-13]
            r6[5]=0
            r6[6]=0
        end
        CD_rand = CELLDM(ibrav1, rand(6))
        (a1,a2,a3) = celldm_to_basis_Angstrom(CD_rand)
        test_crystal(a1, a2, a3, SG, true)
    end
end


##

(a1, a2, a3) = ([0.2469166181700511, 0.0, 0.0], [0.0077459509458863655, 0.060434548288843246, 0.0], [0.00012448256765456707, 2.1053132333292574e-6, 0.0005533350335722911])
UC0 = UNITCELL(:Angstrom, 2, 0.2469166181700511, [1.0, 0.0, 0.0], [0.03137071535846057, 0.2447569091814714, 0.0], [0.0005041481961689437, 8.526413689496312e-6, 0.0022409793138800003])
UC1 = UNITCELL(:Bohr, 2, 0.46660478396019167, [1.0, 0.0, 0.0], [0.007793800052846888, 0.24663601329673016, 0.0], [0.0005041481961689439, 8.526413689496284e-6, 0.0022409793138800003])

##

cell_to_basis_Angstrom(UC1)