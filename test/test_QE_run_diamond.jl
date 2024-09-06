include("test/common_header.jl")

##* ===================================================

using JLD2

function test(PS)
    dict_relax = Dict("title"  => "graphene", 
                    "prefix" => "graphene", 
                    "calc"   => "vc-relax", 
                    "assume_isolated"    => "2D",
                    :do_not_use_symmetry => true,

                    "outdir"=>"/data/test",
                    "pseudo_folder"=>PSEUDO_PATHS[PS], 
                    "pseudo_files" => make_pseudofile_dict(PS,["C"]),
                    "positions" => [("C", 1/3, 2/3, 1/2), ("C", 2/3, 1/3, 1/2)],
                    :cif=>(2.46,2.46,2.46*12, 90,90,120, 1),

                    "uniqueb" => true, "origin_choice"=>1, "rhombohedral"=>false,

                    "conv_thr"        => 1.0e-8,
                    "nbnd"            => 8,
                    "mixing_beta"     => 0.8,
                    "ecutwfc"         => 100.0, 
                    "ecutrho"         => 400.0,
                    "electron_maxstep"=> 200,
                    "mixing_mode"     => "plain",
                    "diagonalization" => "david",

                    "cell_dofree"     => "2Dxy", 
                    "press_conv_thr"  => 0.01, 
                    "cell_dynamics"   => "bfgs", 
                    "bfgs_ndim"       => 3,

                    "ion_dynamics"    => "bfgs", 
                    "upscale"         => 100.0,

                    :kpoint_mode=>"automatic", 
                    :kpoints=>(8,8,1,0,0,0), 
    )

    c = pwx_inupt_relax(dict_relax) ;

    R = pw( 4, c;
        workspace="/data/test",
        fin="pw.x.tmp.$(PS).k8x8.in",
        fout="pw.x.tmp.$(PS).k8x8.out",
        watchdog_setting = (woof_per_x_min=1/20, max_watch=7200, max_silence_woofs=5, tail_length=10, quiet=false)
    )

    return R
end

##* ===================================================

results_k8x8 = []

for PS ∈ ["GHH_PBE", "SSSP_Efficiency", "MT_PW_LDA", "GHH_PZ", "SG15", "MT_PBE", "SSSP_Precision"]
    @info "Now testing $PS ..."
    R = test(PS)
    push!(results_k8x8, R)
end

@save "/data/test/k8x8.jld2" results_k8x8

##* ===================================================

for R ∈ results_k8x8
    @show pw_relax_cell_parameters(R[1])
end
