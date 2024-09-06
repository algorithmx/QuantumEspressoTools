# pwx_output_file.jl

## -------------------------------------------------------------
## process pw.x results
## -------------------------------------------------------------



function ph_qpoints(res)
    _f_rstr = r"-?\d+\.\d*((e|E)-?\d+)?"
    _3f_rstr = _f_rstr*r"\s+"*_f_rstr*r"\s+"*_f_rstr
    _1d3f_rstr = r"\d+\s+" * _3f_rstr

    p = find_line(res, r"Dynamical\s+matrices\s+for\s*")
    if p===nothing
        @error "ph_qpoints():\nNo q-points found."
        return Vector{Float64}[]
    else
        nq = extract(res[p:end], r"\(\s*\d+\s*q\-points\s*\)\s*\:", r"\d+") |> strtonum
        qpoints = extract_all(res[p:p+nq+4], _1d3f_rstr, _3f_rstr) .|> parse3f
        @assert length(qpoints)==nq
        return qpoints
    end
end


function ph_qpoint_section(res)
    qpos = find_all_lines(res, r"Calculation\s+of\s+q\s+=\s*")
    if qpos===nothing
        @error "ph_qpoint_section():\nNo q-point sections found."
        return Vector{String}[]
    else
        qpos2 = [qpos[2:end].-1; length(res)]
        return Vector{String}[res[a:b] for (a,b) ∈ zip(qpos,qpos2)]
    end
end


function ph_qpoint_coords(res_section)
    _f_rstr = r"-?\d+\.\d*"
    _3f_rstr = _f_rstr * r"\s+" * _f_rstr * r"\s+" * _f_rstr
    extract(res_section, r"Calculation\s+of\s+q\s*\=\s*"*_3f_rstr, _3f_rstr) |> parse3f
end


function ph_qpoint_freqs(res_section)
    # freq (    1) =       2.536750 [THz] =      84.616876 [cm-1]
    _f_rstr = r"-?\d+\.\d*((e|E)-?\d+)?"
    freq_rstr = r"freq\s*\(\s*\d+\s*\)\s*\=\s*" * _f_rstr
    fpos = find_all_lines(res_section, freq_rstr)
    if fpos===nothing || length(fpos)==0
        @error "ph_qpoint_freqs():\nNo frequencies found."
        return Float64[]
    else
        fvals_THz  = extract_all(res_section[fpos], _f_rstr*r"\s*\[THz\]",   _f_rstr) .|> strtonum
        fvals_cm_1 = extract_all(res_section[fpos], _f_rstr*r"\s*\[cm\-1\]", _f_rstr) .|> strtonum
        return fvals_THz, fvals_cm_1
    end
end

ph_frequency(result_lines) = Dict(  ph_qpoint_coords(sect)=>ph_qpoint_freqs(sect)
                                    for sect in ph_qpoint_section(result_lines)  )