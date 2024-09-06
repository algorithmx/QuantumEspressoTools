# https://www.quantum-espresso.org/pseudopotentials/naming-convention

##* ==============================================================

n(t,f,r) = r"Pseudopotential\s+type\:\s+$(t)\s*\nFunctional\s+type\:\s+$(f)\s*\n(Non\s+Linear\s+Core\s+Correction\s*\n)?$(r)\s+relativistic"

PS = Dict(
    "USPP_LDA_SR"    => n("USPP", "LDA", "Scalar"), 
    "USPP_LDA_FR"    => n("USPP", "LDA", "Full"), 
    "USPP_PBE_SR"    => n("USPP", "PBE", "Scalar"),
    "USPP_PBE_FR"    => n("USPP", "PBE", "Full"),
    "USPP_PBESOL_SR" => n("USPP", "PBESOL", "Scalar"),
    "USPP_PBESOL_FR" => n("USPP", "PBESOL", "Full"),
    "PAW_PZ_SR"      => n("PAW", r"PERDEW\-ZUNGER\s+\(LDA\)\s+exch\-corr", "Scalar"), 
    "PAW_PZ_FR"      => n("PAW", r"PERDEW\-ZUNGER\s+\(LDA\)\s+exch\-corr", "Full"),
    "PAW_PBE_SR"     => n("PAW",  "PBE", "Scalar"),
    "PAW_PBE_FR"     => n("PAW",  "PBE", "Full"),
    "PAW_PBESOL_SR"  => n("PAW",  "PBESOL", "Scalar"),
    "PAW_PBESOL_FR"  => n("PAW",  "PBESOL", "Full"),
)

PS_ROOT = "/data"

PS_FILES = Dict(ps=>readdir("$PS_ROOT/$ps/") for ps in keys(PS))

##* ==============================================================

get_el(fn) = first(split(fn,".",keepempty=false))

PS_FILES_REFINED = Dict((el,ps)=>[f for f in PS_FILES[ps] if startswith(f,el)] 
                        for ps in keys(PS_FILES) 
                            for el in unique(get_el.(PS_FILES[ps])))

##* ==============================================================

Sc = Dict((k,x)=>readlines("/data/$k/$x")[5:40] for k in keys(PS) for x in PS_FILES_REFINED[("Sc",k)])

