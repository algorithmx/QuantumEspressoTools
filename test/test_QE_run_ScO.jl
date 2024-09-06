include("test/common_header.jl")

global const TEST_FOLDER = ""

##* ===================================================

using JLD2

function test(PS)
    dict_relax = Dict(
        :calc             => "vc-relax", 
        "title"           => "ScO", 
        "prefix"          => "ScO", 
        "outdir"   => TEST_FOLDER,

        :pseudo_mode      => PS,

        "positions"       => [("Sc", 0.0, 0.0, 0.0), ("O", 1/2, 1/2, 1/2)],
        :cif             => (4.47983381807588,4.47983381807588,4.47983381807588, 90,90,90, 225),
        :do_not_use_symmetry => true,

        "nbnd"            => 20,
        "conv_thr"        => 1.0e-8,
        "mixing_beta"     => 0.8,
        "ecutwfc"         => 100.0, 
        "ecutrho"         => 400.0,
        "electron_maxstep"=> 2000,

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
        :kpoints=>(8,8,8,0,0,0), 
    )

    c = pwx_inupt_relax(dict_relax) ;

    if !(isfile("$TEST_FOLDER/pw.x.tmp.$(PS).k8x8.out") && verify_QE_result(readlines("$TEST_FOLDER/pw.x.tmp.$(PS).k8x8.out"), `pw.x`))
        R = pw( 4, c;
                workspace="$TEST_FOLDER",
                fin="pw.x.tmp.$(PS).k8x8.in",
                fout="pw.x.tmp.$(PS).k8x8.out",
                watchdog_setting = (woof_per_x_min=1/3, max_watch=7200, max_silence_woofs=20, tail_length=10, quiet=false)
        )
        return R
    else
        return (readlines("$TEST_FOLDER/pw.x.tmp.$(PS).k8x8.out"), true)
    end
end


#* ===================================================

function test1(PS::String, cell::Dict)
    cd("$TEST_FOLDER")
    try mkdir("ScO_$(PS)") catch; nothing end
    dict_scf = Dict(
        :calc             => "scf", 
        "title"           => "ScO_$(PS)", 
        "prefix"          => "ScO_$(PS)", 
        "outdir"   => "$TEST_FOLDER/ScO_$(PS)",

        :pseudo_mode      => PS,

        "nbnd"            => 20,
        "conv_thr"        => 1.0e-11,
        "mixing_beta"     => 0.8,
        "ecutwfc"         => 100.0, 
        "ecutrho"         => 400.0,
        "electron_maxstep"=> 2000,

        "occupations"       => "smearing",
        "smearing"          => "mp",
        "degauss"           => 0.02,

        :kpoint_mode=>"automatic", 
        :kpoints=>(8,8,8,0,0,0), 
    ) ∪ cell

    c = pwx_inupt_compute_energy(dict_scf) ;

    if !(isfile("$TEST_FOLDER/ScO_$(PS)/pw.x.tmp.$(PS).k8x8.out") && verify_QE_result(readlines("$TEST_FOLDER/ScO_$(PS)/pw.x.tmp.$(PS).k8x8.out"), `pw.x`))
        R = pw( 4, c;
                workspace="$TEST_FOLDER/ScO_$(PS)",
                fin="pw.x.tmp.$(PS).k8x8.in",
                fout="pw.x.tmp.$(PS).k8x8.out",
                watchdog_setting = (woof_per_x_min=1/3, max_watch=7200, max_silence_woofs=20, tail_length=10, quiet=false)
        )
        return R
    else
        return (readlines("$TEST_FOLDER/ScO_$(PS)/pw.x.tmp.$(PS).k8x8.out"), true)
    end

end

#* ===================================================

function test2(PS::String)
    if !isdir("$TEST_FOLDER/ScO_$(PS)") return end
    cd("$TEST_FOLDER/ScO_$(PS)")

    scf_result = readlines("pw.x.tmp.$(PS).k8x8.out")

    dict_phx = Dict(
        :ph_mode             => :single,
        "prefix"          => "ScO_$(PS)", 
        "title"           => "ScO_$(PS)",
        "outdir"   => "$TEST_FOLDER/ScO_$(PS)",
        :reciprocal_basis=>pw_reciprocal_axes(scf_result),
        :qpoints         =>[[0.0,0.0,0.0],],
    )

    c = phx_inupt_single(dict_phx) ;
    @info @show c

    if  !( isfile("$$TEST_FOLDER/ScO_$(PS)/ph.x.tmp.$(PS).k8x8.out") 
        && verify_QE_result(readlines("$$TEST_FOLDER/ScO_$(PS)/ph.x.tmp.$(PS).k8x8.out"), `ph.x`) )
        R = ph( 4, c;
                workspace="$$TEST_FOLDER/ScO_$(PS)",
                fin="ph.x.tmp.$(PS).k8x8.in",
                fout="ph.x.tmp.$(PS).k8x8.out",
                watchdog_setting = (woof_per_x_min=1/5, max_watch=7200, max_silence_woofs=20, tail_length=10, quiet=false)
        )
        return R
    else
        return (readlines("$TEST_FOLDER/ScO_$(PS)/ph.x.tmp.$(PS).k8x8.out"), true)
    end

end

#* ===================================================

const PS_TO_TEST = ["SSSP_Efficiency", "MT_PW_LDA", "GHH_PZ", "SG15", "MT_PBE", "SSSP_Precision", "GBRV_LDA", "GBRV_PBE", "GBRV_PBESOL"]

##* ===================================================

results_k8x8 = []

for PS ∈ PS_TO_TEST
    @info "Now computing $PS ..."
    R0 = test(PS)
    R1 = test1(PS, pw_prepare_next_from_BFGS_result(QE_default_symmetry_group_convention, R0[1]))
    push!(results_k8x8, R1)
end

@save "$TEST_FOLDER/ScO_k8x8.energy.jld2" results_k8x8

##* ===================================================

results_k8x8 = []

for PS ∈ ["SSSP_Efficiency", "MT_PW_LDA", "GHH_PZ", "SG15", "MT_PBE", "SSSP_Precision", "GBRV_LDA", "GBRV_PBE", "GBRV_PBESOL"]
    @info "Now running ph.x on $PS ..."
    R = test2(PS)
    push!(results_k8x8, R)
end

@save "$TEST_FOLDER/ScO_k8x8.phonon.jld2" results_k8x8

##* ===================================================
