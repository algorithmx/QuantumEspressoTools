# pw2wannier90x_input.jl
# ref : ???

#* ============================================================

function INPUTPPW90(config::Dict)
    NL = Namelist("INPUTPP",Dict(),Dict(),Dict())
    prefix   = config["prefix"]
    seedname = config["seedname"]
    NL.req   = Dict(
            "prefix"        => prefix,
            "seedname"      => seedname,
            "write_mmn"     => true,
            "write_amn"     => true,
            )
    NL.default = Dict(
            "outdir"        => "./",
            "write_unk"     => false,
            "write_uHu"     => false,
            "spin_component"=> "none",
    )
    NL.opt = NL.default ← config
    return NL
end


#* =========================================================

function pw2wannier90x_inupt(config::Dict)
    return  INPUT_PW2WANNIER90("", INPUTPPW90(config)) |> build_input_file
end

