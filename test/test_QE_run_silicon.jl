include("test/common_header.jl")

##* ===================================================

using JLD2

alat_au         = 10.20
a0              = alat_au * _BOHR_RADIUS_
basis_Ang_final = a0 .* celldm_to_basis([alat_au,0,0,0,0,0], 2)
latt_p          = basis_to_lattice_parameters(basis_Ang_final)

function test(PS)
    dict_relax = Dict("title"  => "silicon", 
                    "prefix" => "silicon", 
                    "calc"   => "vc-relax", 

                    "outdir"=>"/data/test",
                    "pseudo_folder"=>PSEUDO_PATHS[PS], 
                    "pseudo_files" => make_pseudofile_dict(PS,["Si"]),
                    "positions" => [("Si", 0.0, 0.0, 0.0), ("Si", 1/4, 1/4, 1/4)],
                    :cif=>(latt_p..., 227),
                    :do_not_use_symmetry => true,

                    "uniqueb" => true, "origin_choice"=>1, "rhombohedral"=>false,

                    "conv_thr"        => 1.0e-8,
                    "nbnd"            => 8,
                    "mixing_beta"     => 0.8,
                    "ecutwfc"         => 100.0, 
                    "ecutrho"         => 400.0,
                    "electron_maxstep"=> 2000,
                    "mixing_mode"     => "plain",
                    "diagonalization" => "david",

                    "cell_dofree"     => "all", 
                    "press_conv_thr"  => 0.01, 
                    "cell_dynamics"   => "bfgs", 
                    "bfgs_ndim"       => 3,

                    "ion_dynamics"    => "bfgs", 
                    "upscale"         => 100.0,

                    :kpoint_mode=>"automatic", 
                    :kpoints=>(8,8,8,0,0,0), 
    )

    c = pwx_inupt_relax(dict_relax) ;

    if !(isfile("/data/test/pw.x.tmp.$(PS).k8x8.out") && verify_QE_result(readlines("/data/test/pw.x.tmp.$(PS).k8x8.out"), `pw.x`))
        R = pw( 4, c;
                workspace="/data/test",
                fin="pw.x.tmp.$(PS).k8x8.in",
                fout="pw.x.tmp.$(PS).k8x8.out",
                watchdog_setting = (woof_per_x_min=1/10, max_watch=7200, max_silence_woofs=5, tail_length=10, quiet=false)
        )
        return R
    else
        return (readlines("/data/test/pw.x.tmp.$(PS).k8x8.out"), true)
    end

end

##* ===================================================

results_k8x8 = []

for PS ∈ ["SSSP_Efficiency", "MT_PW_LDA", "GHH_PZ", "SG15", "MT_PBE", "SSSP_Precision", "GHH_PBE", "GBRV_LDA", "GBRV_PBE", "GBRV_PBESOL"]
    @info "Now testing $PS ..."
    R = test(PS)
    push!(results_k8x8, R)
end

@save "/data/test/silicon_k8x8.jld2" results_k8x8

##* ===================================================

for R ∈ results_k8x8
    @show pw_relax_cell_parameters(R[1])
end
