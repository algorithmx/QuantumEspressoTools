
## -------------------------------------------------
## DATA
## -------------------------------------------------

## source ???
# global const _BOHR_RADIUS_ = 0.52917720859 
## source Mathematica : Quantity[1, "BohrRadius"]/Quantity[1, "Angstroms"] // UnitConvert
#global const _BOHR_RADIUS_ = 0.52917721067
#!! cat constants.f90 | grep "BOHR_"
global const _BOHR_RADIUS_ = 0.529177210903  


global const AtomSymb = split("H He Li Be B C N O F Ne Na Mg Al Si P S Cl Ar K Ca Sc Ti V Cr Mn Fe Co Ni Cu Zn Ga Ge As Se Br Kr Rb Sr Y Zr Nb Mo Tc Ru Rh Pd Ag Cd In Sn Sb Te I Xe Cs Ba La Ce Pr Nd Pm Sm Eu Gd Tb Dy Ho Er Tm Yb Lu Hf Ta W Re Os Ir Pt Au Hg Tl Pb Bi Po At Rn Fr Ra Ac Th Pa U Np Pu Am Cm Bk Cf Es Fm Md No Lr Rf Db Sg Bh Hs Mt", " ", keepempty=false)

function nshell(i)
    if i<=2
        return 2
    elseif i<=10
        return 8
    elseif i<=18
        return 8
    elseif i<=20
        return 2+6
    elseif i<=30
        return 2+10
    elseif i<=36
        return 2+6
    elseif i<=38
        return 2+6
    elseif i<=48
        return 2+10
    elseif i<=54
        return 2+6+10
    elseif i<=56
        return 2+6+10
    elseif i<=71
        return 2+10+14
    elseif i<=80
        return 2+6+10
    elseif i<=86
        return 2+6+10+14
    elseif i<=88
        return 2+6+10+14
    elseif i<=103
        return 2+6+10+14
    elseif i<=112
        return 2+6+10+14
    else
        return 0
    end
end

global const N_SHELL_ELECTRONS = Dict(a=>nshell(i) for (i,a) in enumerate(AtomSymb))


#TODO use NIST data
global const AtomMass = Float64[
1.0079,4.0026,6.941,9.0122,10.811,12.0107,14.0067,15.9994,18.9984,20.1797,22.9897,24.305,26.9815,28.0855,30.9738,32.065,35.453,
39.948,39.0983,40.078,44.9559,47.867,50.9415,51.9961,54.938,55.845,58.9332,58.6934,63.546,65.39,69.723,72.64,74.9216,78.96,79.904,
83.8,85.4678,87.62,88.9059,91.224,92.9064,95.94,98,101.07,102.906,106.42,107.868,112.411,114.818,118.71,121.76,127.6,126.904,131.293,
132.905,137.327,138.905,140.116,140.908,144.24,145,150.36,151.964,157.25,158.925,162.5,164.93,167.259,168.934,173.04,174.967,178.49,
180.948,183.84,186.207,190.23,192.217,195.078,196.966,200.59,204.383,207.2,208.98,209,210,222,223,226,227,232.038,231.036,238.029,
237,244,243,247,247,251,252,257,258,259,262,261,262,266,264,277,268
]


global const Atoms = Dict(AtomSymb[i]=>i for i=1:length(AtomSymb))

global const AtomMasses = Dict(AtomSymb[i]=>AtomMass[i] for i=1:length(AtomMass))


#> International Tables
# 1-2 Triclinic, 3-15 Monoclinic, 16-74 Orthorhombic, 75-142 Tetragonal, 143-167 Trigonal, 168-194 Hexagonal, 195-230 Cubic
global const IT_NUMBER_CRYSTAL_CLASS = [1:2, 3:15, 16:74, 75:142, 143:167, 168:194, 195:230]

function CRYSTAL_CLASS(SG)
    if     1  <=SG && SG<=2
        return "triclinic" 
    elseif 3  <=SG && SG<=15 
        return "monoclinic" 
    elseif 16 <=SG && SG<=74 
        return "orthorhombic" 
    elseif 75 <=SG && SG<=142 
        return "tetragonal" 
    elseif 143<=SG && SG<=167 
        return "trigonal" 
    elseif 168<=SG && SG<=194 
        return "hexagonal" 
    elseif 195<=SG && SG<=230
        return "cubic" 
    else
        return "wtf"
    end
end

global const Int_Tables = split("P1 P-1 P2 P2(1) C2 Pm Pc Cm Cc P2/m P2(1)/m C2/m P2/c P2(1)/c C2/c P222 P222(1) P2(1)2(1)2 P2(1)2(1)2(1) C222(1) C222 F222 I222 I2(1)2(1)2(1) Pmm2 Pmc2(1) Pcc2 Pma2 Pca2(1) Pnc2 Pmn2(1) Pba2 Pna2(1) Pnn2 Cmm2 Cmc2(1) Ccc2 Amm2 Abm2 Ama2 Aba2 Fmm2 Fdd2 Imm2 Iba2 Ima2 Pmmm Pnnn Pccm Pban Pmma Pnna Pmna Pcca Pbam Pccn Pbcm Pnnm Pmmn Pbcn Pbca Pnma Cmcm Cmca Cmmm Cccm Cmma Ccca Fmmm Fddd Immm Ibam Ibca Imma P4 P4(1) P4(2) P4(3) I4 I4(1) P-4 I-4 P4/m P4(2)/m P4/n P4(2)/n I4/m I4(1)/a P422 P42(1)2 P4(1)22 P4(1)2(1)2 P4(2)22 P4(2)2(1)2 P4(3)22 P4(3)2(1)2 I422 I4(1)22 P4mm P4bm P4(2)cm P4(2)nm P4cc P4nc P4(2)mc P4(2)bc I4mm I4cm I4(1)md I4(1)cd P-42m P-42c P-42(1)m P-42(1)c P-4m2 P-4c2 P-4b2 P-4n2 I-4m2 I-4c2 I-42m I-42d P4/mmm P4/mcc P4/nbm P4/nnc P4/mbm P4/mnc P4/nmm P4/ncc P4(2)/mmc P4(2)/mcm P4(2)/nbc P4(2)/nnm P4(2)/mbc P4(2)/mnm P4(2)/nmc P4(2)/ncm I4/mmm I4/mcm I4(1)/amd I4(1)/acd P3 P3(1) P3(2) R3 P-3 R-3 P312 P321 P3(1)12 P3(1)21 P3(2)12 P3(2)21 R32 P3m1 P31m P3c1 P31c R3m R3c P-31m P-31c P-3m1 P-3c1 R-3m R-3c P6 P6(1) P6(5) P6(2) P6(4) P6(3) P-6 P6/m P6(3)/m P622 P6(1)22 P6(5)22 P6(2)22 P6(4)22 P6(3)22 P6mm P6cc P6(3)cm P6(3)mc P-6m2 P-6c2 P-62m P-62c P6/mmm P6/mcc P6(3)/mcm P6(3)/mmc P23 F23 I23 P2(1)3 I2(1)3 Pm-3 Pn-3 Fm-3 Fd-3 Im-3 Pa-3 Ia-3 P432 P4(2)32 F432 F4(1)32 I432 P4(3)32 P4(1)32 I4(1)32 P-43m F4-3m I-43m P-43n F-43c I-43d Pm-3m Pn-3n Pm-3n Pn-3m Fm-3m Fm-3c Fd-3m Fd-3c Im-3m Ia-3d", " ", keepempty=false) .|> string

global const Int_Tables_No = Dict(Int_Tables[i]=>i for i=1:length(Int_Tables))

# ---------------------------------------------------------

#@inline take0(dic) = [v[1] for (k,v) in dic if v[2]==0] |> sort
#KVECS_BY_CRYSTAL_CLASS = [  union(unique([take0(kvecs[i]) for i in IT_NUMBER_CRYSTAL_CLASS[k]])...)
#                            for k=1:length(IT_NUMBER_CRYSTAL_CLASS)  ]

global const KVECS_BY_CRYSTAL_CLASS_QE_DEFAULT = [
    ["(0,0,0)", "(0,0,1/2)", "(0,1/2,0)", "(0,1/2,1/2)", "(1/2,0,0)", "(1/2,0,1/2)", "(1/2,1/2,0)", "(1/2,1/2,1/2)"],
    ["(0,0,0)", "(0,0,1/2)", "(0,1/2,0)", "(0,1/2,1/2)", "(1/2,0,0)", "(1/2,0,1/2)", "(1/2,1/2,0)", "(1/2,1/2,1/2)", "(0,0,1)", "(1/2,0,1)"],
    ["(0,0,0)", "(0,0,1/2)", "(0,1/2,0)", "(0,1/2,1/2)", "(1/2,0,0)", "(1/2,0,1/2)", "(1/2,1/2,0)", "(1/2,1/2,1/2)", "(1,0,0)", "(1,0,1/2)", "(0,0,1)", "(0,1,0)", "(1,1,1)", "(0,0,-1)", "(0,1/2,-1/2)", "(1/2,0,-1)", "(1/2,1/2,-1/2)"],
    ["(0,0,0)", "(0,0,1/2)", "(0,1/2,0)", "(0,1/2,1/2)", "(1/2,1/2,0)", "(1/2,1/2,1/2)", "(1,1,1)", "(1/2,0,1/2)"],
    ["(0,0,0)", "(0,0,1/2)", "(1/2,0,0)", "(1/2,0,1/2)", "(1/3,1/3,0)", "(1/3,1/3,1/2)", "(0,1/2,0)", "(1/2,1/2,0)", "(1/2,1/2,1/2)"],
    ["(0,0,0)", "(0,0,1/2)", "(1/2,0,0)", "(1/2,0,1/2)", "(1/3,1/3,0)", "(1/3,1/3,1/2)"],
    ["(0,0,0)", "(0,1/2,0)", "(1/2,1/2,0)", "(1/2,1/2,1/2)", "(0,1,0)", "(1/2,1,0)", "(1,1,1)"],
]


