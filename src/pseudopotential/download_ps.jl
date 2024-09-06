using Distributed 

addprocs(40)

@everywhere using HTTP
@everywhere using Gumbo

@everywhere const AtomSymb = split("H He Li Be B C N O F Ne Na Mg Al Si P S Cl Ar K Ca Sc Ti V Cr Mn Fe Co Ni Cu Zn Ga Ge As Se Br Kr Rb Sr Y Zr Nb Mo Tc Ru Rh Pd Ag Cd In Sn Sb Te I Xe Cs Ba La Ce Pr Nd Pm Sm Eu Gd Tb Dy Ho Er Tm Yb Lu Hf Ta W Re Os Ir Pt Au Hg Tl Pb Bi Po At Rn Fr Ra Ac Th Pa U Np Pu Am Cm Bk Cf Es Fm Md No Lr Rf Db Sg Bh Hs Mt", " ", keepempty=false)

@everywhere function search_in_page(h::HTMLElement, level::Int, kw::Union{String,Regex})
    if level > 9
        return [n for n ∈ h.children if length(search_in_page(n,level+1,kw))>0]
    else
        vcat([search_in_page(n, level+1, kw) for n ∈ h.children]...)
    end
end

@everywhere function search_in_page(t::HTMLText, level::Int, kw::Union{String,Regex})
    return occursin(kw, t.text) ? [t] : []
end

@everywhere function download_ps(el, spec, FD)
    @info "inside download_ps($el, spec, $FD)"
    url(el) = "https://www.quantum-espresso.org/pseudopotentials/ps-library/$el"
    r = undef
    html = undef
    try
        r  = HTTP.get(url(el))
        html = parsehtml(String(r.body))
    catch _
        @warn "error happend to $el , $FD"
        return
    end
    
    for k ∈ search_in_page(html.root, 0, spec)
        href = k.children[1].attributes["href"]
        if !isfile("$FD/$(basename(href))")
            try
                @show "Downloading $href"
                download("https://www.quantum-espresso.org$href", "$FD/$(basename(href))")
            catch _
                @warn "error happend to $el , $FD"
                return
            end
        end
    end
    return
end

n(t,f,r) = r"Pseudopotential\s+type\:\s+$(t)\s*\nFunctional\s+type\:\s+$(f)\s*\n(Non\s+Linear\s+Core\s+Correction\s*\n)?$(r)\s+relativistic"

@everywhere PS = Dict(
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

try mkdir("/data/PS") catch ; nothing end

for ps in keys(PS)
    try mkdir("/data/PS/$ps") catch ; nothing end
    try mkdir("/data/$ps") catch ; nothing end
end

@everywhere g(x,el) = download_ps(el, PS[x], "/data/$x")

##

@everywhere f(x) = g(x[1],lowercase(x[2]))

@everywhere XXX = [(ps,el) for el in AtomSymb for ps in collect(keys(PS))]

pmap(f, XXX)

##

#check_tail() = (cd("/home/dabajabaza/abinitio/pseudopotentials/PSLibrary"); run("for X in \$(ls); do diff <(ls \$X | wc | awk '{print $1}') <(tail \$X/*.UPF | grep \"</UPF>\" | wc | awk '{print \$1}'); done"))
