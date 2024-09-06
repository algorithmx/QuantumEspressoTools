function BANDS(config::Dict)
    NL = Namelist("BANDS",Dict(),Dict(),Dict())
    NL.req  = Dict(
            "prefix"        => config["prefix"],
            "outdir"        => config["outdir"],
            "filband"       => "$(prefix).bands"
    )
    NL.default = Dict(
            "lsym"          => true,
            "no_overlap"    => true,
    )
    NL.opt = NL.default ← config
    return NL
end


# this is a condensation of experience
function bandsx_input(config::Dict)
    return  INPUT_BANDS( config["title"], BANDS(config) ) |> build_input_file
end

