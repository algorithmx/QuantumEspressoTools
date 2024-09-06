function pw_decreasing_beta(
    np::Int,
    config0::Dict,
    pw_script_generator::Function,
    beta_diag_pairs=[];
    workspace="./",
    fin="pw.x.tmp.in",
    fout="pw.x.tmp.out",
    watchdog_setting = (woof_per_x_min=6, max_watch=720, max_silence_woofs=5, tail_length=10, quiet=false)
    )::Tuple{Vector{String},Bool,Float64,String}
    
    if length(beta_diag_pairs)==0
        (r,s) = pw( np, pw_script_generator(config0); 
                    workspace=workspace, fin=fin, fout=fout,
                    watchdog_setting=watchdog_setting )
        return (r, s, 0.0, "")
    else
        r = String[]
        s = false
        beta = 0.0
        diag = ""
        for (beta,diag) ∈ beta_diag_pairs
            config = config0 ⬱ Dict("diagonalization"=>diag, "mixing_beta"=>beta)  #* haha
            (r,s) = pw( np, pw_script_generator(config); 
                        workspace=workspace, fin=fin, fout=fout,
                        watchdog_setting=watchdog_setting )
            if s return (r, s, beta, diag) end
        end

        @warn "pw_decreasing_beta():\nAfter trying out \n$(beta_diag_pairs)\nAll calculations failed."
        return (r, s, beta, diag)
    end
end
