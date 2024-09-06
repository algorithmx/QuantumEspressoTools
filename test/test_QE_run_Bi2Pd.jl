include("test/common_header.jl")

##* ===================================================

using JLD2

alat_au         = 6.3341
a0              = alat_au * _BOHR_RADIUS_
cd3             = 4.4510
basis_Ang_final = a0 .* celldm_to_basis([alat_au, 0, cd3, 0,0,0], 6)
latt_p          = basis_to_lattice_parameters(basis_Ang_final)

function test(PS)
    dict_relax = Dict("title"  => "Bi2Pd", 
                    "prefix" => "Bi2Pd", 
                    "calc"   => "vc-relax", 

                    "outdir"=>"/data/test",
                    "positions" => [("Bi", 1/2, 1/2, 0.3768), ("Bi", 1/2, 1/2, 1-0.3768), ("Pd", 0.0, 0.0, 1/2)],
                    :cif=>(latt_p..., 1),
                    :do_not_use_symmetry => true,
                    :pseudo_mode => PS,

                    "uniqueb" => true, "origin_choice"=>1, "rhombohedral"=>false,

                    "conv_thr"        => 1.0e-8,
                    "mixing_beta"     => 0.8,
                    "ecutwfc"         => 100.0, 
                    "ecutrho"         => 400.0,
                    "electron_maxstep"=> 2000,
                    "mixing_mode"     => "plain",
                    "diagonalization" => "david",

                    "occupations"       => "smearing",
                    "smearing"          => "mp",
                    "degauss"           => 0.02,

                    "cell_dofree"     => "all", 
                    "press_conv_thr"  => 0.01, 
                    "cell_dynamics"   => "bfgs", 
                    "bfgs_ndim"       => 3,

                    "ion_dynamics"    => "bfgs", 
                    "upscale"         => 100.0,

                    :kpoint_mode=>"automatic", 
                    :kpoints=>(16,16,1,0,0,0), 
    )

    c = pwx_inupt_relax(dict_relax) ;

    if !(isfile("/data/test/pw.x.tmp.$(PS).k16x16.out") && verify_QE_result(readlines("/data/test/pw.x.tmp.$(PS).k16x16.out"), `pw.x`))
        R = pw( 4, c;
                workspace="/data/test",
                fin="pw.x.tmp.$(PS).k16x16.in",
                fout="pw.x.tmp.$(PS).k16x16.out",
                watchdog_setting = (woof_per_x_min=1/5, max_watch=7200, max_silence_woofs=20, tail_length=10, quiet=false)
        )
        return R
    else
        return (readlines("/data/test/pw.x.tmp.$(PS).k16x16.out"), true)
    end

end

##* ===================================================

results_k16x16 = []

for PS ∈ ["SSSP_Efficiency", "MT_PW_LDA", "GHH_PZ", "SG15", "MT_PBE", "SSSP_Precision", "GHH_PBE", "GBRV_LDA", "GBRV_PBE", "GBRV_PBESOL"]
    @info "Now testing $PS ..."
    R = test(PS)
    push!(results_k16x16, R)
end

@save "/data/test/Bi2Pd_k16x16.jld2" results_k16x16

##* ===================================================

for R ∈ results_k16x16
    @show pw_relax_cell_parameters(R[1])
end
