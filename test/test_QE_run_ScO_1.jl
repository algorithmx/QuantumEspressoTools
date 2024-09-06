include("/home/dabajabaza/jianguoyun/Workspace/QuantumEspressoTools/scripts/scripts.jl")

global const TEST_FOLDER = "/data/test"

#* ===================================================

d1 = Dict(
        "mixing_beta"     => 0.8,
        "ecutwfc"         => 150.0, 
        "ecutrho"         => 600.0,
        "electron_maxstep"=> 500,
        "nbnd"            => 26,
        :kpoint_mode      => "automatic", 
        :kpoints          => (4,4,4,1,1,1),

        :watchdog_setting => (woof_per_x_min=1/20, max_watch=1440, max_silence_woofs=80, tail_length=20, quiet=false),
)

d2(PS) = Dict(
        "title"           => "ScO_$PS", 
        "prefix"          => "ScO_$PS", 
        :pseudo_mode      => PS,
)

init_struct = Dict(
        "positions"       => [("Sc", 0.0, 0.0, 0.0), ("O", 1/2, 1/2, 1/2)],
        :cif              => (4.47983381807588,4.47983381807588,4.47983381807588, 90,90,90, 225),
        :do_not_use_symmetry => true,
)

dict_relax = Dict(
        "conv_thr"=>1e-8, 
        "ecutwfc"         => 100.0,
        "ecutrho"         => 400.0,

        "occupations"     => "smearing",
        "smearing"        => "mp",
        "degauss"         => 0.1,

        "cell_dofree"     => "all", 
        "press_conv_thr"  => 0.01, 
        "cell_dynamics"   => "bfgs", 
        "bfgs_ndim"       => 3,

        "ion_dynamics"    => "bfgs", 
        "upscale"         => 100.0,
)

dict_phx(PS) = Dict(
        :ph_mode    => :single,
        "prefix"    => "ScO_$(PS)", 
        "title"     => "ScO_$(PS)",
        :qpoints    => [[0.0,0.0,0.0],],
)


function upd1(x)
    d = dict__pw_relax_result(x)
    return (d ↓ ["ibrav", :cif, :reciprocal_basis]) ⬱ Dict("ibrav"=>0,:do_not_use_symmetry=>true)
end

function upd2(x)
    d = dict__pw_result(x)
    return Dict(:reciprocal_basis => d[:reciprocal_basis])
end

instr(PS) = INSTRUCTION(
    Dict(
        `pw.x`=> _electron_scf_seed_conf_ ⬱ (d1 ∪ d2(PS)),
        `ph.x`=> _phonon_seed_conf_ ⬱ (d1 ∪ d2(PS)),
    ),
    [
        ("ScO_relax", `pw.x`, Dict(   :calc=>"vc-relax", 
                                        :updater => upd1,
                                        :pw_mode => "energy",
                                ) ∪ init_struct ∪ dict_relax),
        ("ScO_scf",   `pw.x`, Dict(   :calc=>"scf",
                                        "conv_thr"=>1e-11, 
                                        :updater => upd2,
                                        :pw_mode => "energy",
                                )),
        ("ScO_ph",    `ph.x`, Dict("outdir"=>"$TEST_FOLDER/ScO_scf") ∪ dict_phx(PS)),
    ]
)

##

instr1(PS) = INSTRUCTION(
    Dict(
        `pw.x`=> _electron_scf_seed_conf_ ⬱ (d1 ∪ d2(PS)),
        `ph.x`=> _phonon_seed_conf_ ⬱ (d1 ∪ d2(PS)),
    ),
    [
        ("ScO_relax", `pw.x`, Dict(   :calc=>"vc-relax", 
                                        :updater => upd1,
                                        :pw_mode => "energy",
                                ) ∪ init_struct ∪ dict_relax),
        ("ScO_scf1",   `pw.x`, Dict(   :calc=>"scf",
                                        "conv_thr"=>1e-11, 
                                        :updater => upd2,
                                        :pw_mode => "energy",
                                        :kpoints => (4,4,4,0,0,0),
                                )),
        ("ScO_ph1",    `ph.x`, Dict(    "outdir"=>"$TEST_FOLDER/ScO_scf1",
                                        :kpoints=>(4,4,4,0,0,0)
                                ) ∪ dict_phx(PS)),
    ]
)


dict_phx2(PS) = Dict(
        :ph_mode    => :grid,
        "prefix"    => "ScO_$(PS)", 
        "title"     => "ScO_$(PS)",
        :qpoints    => (2,2,2),
)

dict_matdyn(PS) = Dict( "prefix" => "ScO_$(PS)",
                        "title"  => "ScO_$(PS)",
                        "flfrc"  => "ScO_$(PS).fc",   #* IFC_filename
                        "flfrq"  => "ScO_$(PS).frq",   #* frequency file
                        "flvec"  => "ScO_$(PS).vec",   #* eigen vector 
                        "fleig"  => "ScO_$(PS).eig",   #* eigen value
)

dict_q2r(PS) = Dict(    "prefix" => "ScO_$(PS)",
                        "title"  => "ScO_$(PS)",
                        "acoustic_sum_rule"   => "simple",
                        "fildyn" => "ScO_$(PS).dynmat.",  #* 
                        "flfrc"  => "ScO_$(PS).fc",    #* IFC_filename
)

instr2(PS) = INSTRUCTION(
    Dict(
        `pw.x`     => _electron_scf_seed_conf_ ⬱ (d1 ∪ d2(PS)),
        `ph.x`     => _phonon_seed_conf_ ⬱ (d1 ∪ d2(PS)),
        `q2r.x`    => _phonon_q2r_seed_conf_,
        `matdyn.x` => _phonon_matdyn_seed_conf_,
    ),
    [
        ("ScO_relax",   `pw.x`,         Dict(   :calc=>"vc-relax", 
                                                :updater => upd1,
                                                :pw_mode => "energy",
                                        ) ∪ init_struct ∪ dict_relax),
        ("ScO_scf1",    `pw.x`,         Dict(   :calc=>"scf",
                                                "conv_thr"=>1e-11, 
                                                :updater => upd2,
                                                :pw_mode => "energy",
                                                :kpoints => (4,4,4,0,0,0),
                                        )),
        ("ScO_ph1_grid",`ph.x`,         Dict(   "outdir"=>"$TEST_FOLDER/ScO_scf1",
                                                :kpoints=>(4,4,4,0,0,0)
                                        ) ∪ dict_phx2(PS)),
        ("ScO_ph1_grid",`q2r.x`,        dict_q2r(PS) ),
        ("ScO_ph1_grid",`matdyn.x`,     Dict(   :qpath_rel=>[(0,0,0)=>10,(0.5,0,0)=>10,(0.5,0.5,0)=>10,(0.5,0.5,0.5)=>10,(0,0,0)=>1],
                                        ) ∪ dict_matdyn(PS) ),
    ]
)


const PS_TO_TEST = ["SSSP_Efficiency", "MT_PW_LDA", "GHH_PZ", "SG15", "MT_PBE", "SSSP_Precision", "GBRV_LDA", "GBRV_PBE", "GBRV_PBESOL"]

##* ===================================================

#execute_serial(instr("SSSP_Efficiency"), WORKSPACE=TEST_FOLDER, from_scratch=false)

#execute_serial(instr1("SSSP_Efficiency"), WORKSPACE=TEST_FOLDER, from_scratch=false)

execute_serial(instr2("SSSP_Efficiency"), WORKSPACE=TEST_FOLDER, from_scratch=false)

##* ===================================================

f1 = "/data/test/ScO_relax/ScO_MT_PW_LDA.pw.x.out"

dict__pw_relax_result(readlines(f1))

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
