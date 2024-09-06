function show_banner()
    banner = [
        "",
        "****************** WELCOME ******************",
        "* ",
        "*   QuantumEspressoTools.jl ",
        "*   Developers : Yunlong Lian",
        "*   Version    : 0.2",
        "* ",
        "*   Hints ",
        "*   (1) Pseudopotentials: ",
        "*             "*join(collect(keys(PSEUDO_PATHS)),"\n*             "),
        "*   (2) ",
        "*********************************************",
    ]
    print(join(banner,"\n"))
    println()
end

##