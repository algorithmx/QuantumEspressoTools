function cleanup_config_for_wannier90(c)
    c ↓ ["title", "prefix", "outdir"]
end


function generate_input_script(p::Cmd, c::Dict)
    @info log_info("Now entering generate_input_script() for p=$p", level=1)
    mpi, prog, npool, mode, np, maby = analyze_prog(p) #TODO npool
    i = String[]
    if prog=="pwdry.x"
        c[:calc]    = "scf"
        c[:pw_mode] ="energy"
        i = pwx_inupt_compute_energy(c)
    elseif prog=="pw.x"
        if c[:calc] ∈ ["scf","bands"]
            i = (c[:pw_mode]=="energy" ? pwx_inupt_compute_energy(c) : pwx_inupt_compute_bands(c))
        elseif c[:calc] == "nscf"
            i = pwx_inupt_nscf(c)
        elseif c[:calc] ∈ ["relax", "vc-relax"]
            i = pwx_inupt_relax(c)
        else
            @error  log_error("generate_input_script() : calculation mode $(c[:calc]) and kpoint mode $(c[:kpoint_mode]) is not supported.", level=2)
        end
    elseif prog=="cpdry.x"
        c[:calc] = "cp"
        i = cpx_inupt_dry_run(c ↓ ["mixing_beta", "diagonalization"])
    elseif prog=="cp.x"
        if c[:calc] ∈ ["scf", "nscf"]
            i = cpx_inupt_energy(c ↓ ["mixing_beta", "diagonalization"])
        elseif c[:calc] ∈ ["cp", "vc-cp"]
            i = cpx_inupt_MD(c ↓ ["mixing_beta", "diagonalization"])
        elseif c[:calc] ∈ ["relax", "vc-relax"]
            i = cpx_inupt_relax(c ↓ ["mixing_beta", "diagonalization"])
        else
            @error  log_error("generate_input_script() : cpx calculation mode $(c[:calc]) unsupported.", level=2)
        end
    elseif prog=="phdry.x"
        @assert c[:ph_mode] == :single
        c["verbosity"] = "high" 
        c["only_init"] = true 
        c["reduce_io"] = true 
        c["lqdir"]     = true 
        i = phx_inupt_single(c ↓ "diagonalization")
    elseif prog=="ph.x"
        if c[:ph_mode] == :single
            i = phx_inupt_single(c ↓ "diagonalization")
        elseif c[:ph_mode] == :multiple
            i = phx_inupt_multiple(c ↓ "diagonalization")
        elseif c[:ph_mode] == :grid
            i = phx_inupt_grid(c ↓ "diagonalization")
        elseif c[:ph_mode] == :bands
            i = phx_inupt_bands(c ↓ "diagonalization")
        else
            @error  log_error("generate_input_script() : ph_mode mode $(c[:ph_mode]) is not supported.", level=2)
        end
    elseif prog=="pp.x"
        if c[:pp_mode] == :total_charge
            i = pp_input_plot_total_charge(c)
        else
            @error  log_error("generate_input_script() : pp_mode mode $(c[:pp_mode]) is not supported.", level=2)
        end
    elseif prog=="bands.x"
        #TODO
        i = String[]
    elseif prog=="plotband.x"
        i = plotband_input(c)
    elseif prog=="q2r.x"
        i = q2r_input(c)
    elseif prog=="matdyn.x"
        i = matdyn_band_input(c)
    elseif prog=="ld1.x"
        i = ld1x_inupt(c)
    elseif prog=="wannier90.x"
        i = wannier90x_inupt(cleanup_config_for_wannier90(c))
    elseif prog=="pw2wannier90.x"
        i = pw2wannier90x_inupt(c)
    elseif prog=="projwfc.x"
        i = projwfcx_input(c)
    else
        @error log_fucked_up("generate_input_script(): Unknown program $prog .")
        throw(error("generate_input_script(): Unknown program $prog ."))
    end

    return i
end
