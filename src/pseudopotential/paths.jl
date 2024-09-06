global const PSEUDO_PATHS = Dict(
    "SG15"                => "$(ENV["PSEUDO_ROOT"])/SG15_SR/V1_2",
    "SG15_FR"             => "$(ENV["PSEUDO_ROOT"])/SG15_FR",

    "SSSP_Precision"         => "$(ENV["PSEUDO_ROOT"])/SSSP_precision_pseudos",
    "SSSP_Efficiency"        => "$(ENV["PSEUDO_ROOT"])/SSSP_efficiency_pseudos",
    "SSSP_Precision_PBEsol"  => "$(ENV["PSEUDO_ROOT"])/SSSP_precision_PBEsol",
    "SSSP_Efficiency_PBEsol" => "$(ENV["PSEUDO_ROOT"])/SSSP_efficiency_PBEsol",

    "PSL_USPP_LDA_SR"     => "$(ENV["PSEUDO_ROOT"])/PSLibrary/USPP_LDA_SR",
    "PSL_USPP_LDA_FR"     => "$(ENV["PSEUDO_ROOT"])/PSLibrary/USPP_LDA_FR",
    "PSL_USPP_PBE_SR"     => "$(ENV["PSEUDO_ROOT"])/PSLibrary/USPP_PBE_SR",
    "PSL_USPP_PBE_FR"     => "$(ENV["PSEUDO_ROOT"])/PSLibrary/USPP_PBE_FR",
    "PSL_USPP_PBESOL_SR"  => "$(ENV["PSEUDO_ROOT"])/PSLibrary/USPP_PBESOL_SR",
    "PSL_USPP_PBESOL_FR"  => "$(ENV["PSEUDO_ROOT"])/PSLibrary/USPP_PBESOL_FR",

    "PSL_PAW_PZ_SR"       => "$(ENV["PSEUDO_ROOT"])/PSLibrary/PAW_PZ_SR",
    "PSL_PAW_PZ_FR"       => "$(ENV["PSEUDO_ROOT"])/PSLibrary/PAW_PZ_FR",
    "PSL_PAW_PBE_SR"      => "$(ENV["PSEUDO_ROOT"])/PSLibrary/PAW_PBE_SR",
    "PSL_PAW_PBE_FR"      => "$(ENV["PSEUDO_ROOT"])/PSLibrary/PAW_PBE_FR",
    "PSL_PAW_PBESOL_SR"   => "$(ENV["PSEUDO_ROOT"])/PSLibrary/PAW_PBESOL_SR",
    "PSL_PAW_PBESOL_FR"   => "$(ENV["PSEUDO_ROOT"])/PSLibrary/PAW_PBESOL_FR",

    "MT_PBE"              => "$(ENV["PSEUDO_ROOT"])/Martins-Troullier/PBE",
    "MT_PW_LDA"           => "$(ENV["PSEUDO_ROOT"])/Martins-Troullier/PW_LDA",

    "GHH_PBE"             => "$(ENV["PSEUDO_ROOT"])/Goedecker-Hartwigsen-Hutter-Teter/PBE/pbe",
    "GHH_PBE_SP"          => "$(ENV["PSEUDO_ROOT"])/Goedecker-Hartwigsen-Hutter-Teter/PBE/sp",
    "GHH_PZ"              => "$(ENV["PSEUDO_ROOT"])/Goedecker-Hartwigsen-Hutter-Teter/PZ/pz",
    "GHH_PZ_SP"           => "$(ENV["PSEUDO_ROOT"])/Goedecker-Hartwigsen-Hutter-Teter/PZ/sp",
    "GHH_BLYP"            => "$(ENV["PSEUDO_ROOT"])/Goedecker-Hartwigsen-Hutter-Teter/BLYP/blyp",
    "GHH_BLYP_SP"         => "$(ENV["PSEUDO_ROOT"])/Goedecker-Hartwigsen-Hutter-Teter/BLYP/sp",

    "GBRV_LDA"            => "$(ENV["PSEUDO_ROOT"])/GBRV/all_lda_UPF_v1.5",
    "GBRV_PBE"            => "$(ENV["PSEUDO_ROOT"])/GBRV/all_pbe_UPF_v1.5",
    "GBRV_PBESOL"         => "$(ENV["PSEUDO_ROOT"])/GBRV/all_pbesol_UPF_v1.5",
)

