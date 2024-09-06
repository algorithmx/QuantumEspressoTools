include("/home/dabajabaza/jianguoyun/Workspace/QuantumEspressoTools/scripts/scripts.jl")

global const TEST_FOLDER = "/data/ReO3"

global const PS_TO_TEST = [ "GBRV_LDA", "GBRV_PBE", "GBRV_PBESOL", 
                            "GHH_PBE", "GHH_PZ", 
                            "MT_PBE", "MT_PW_LDA", 
                            "SG15", "SSSP_Efficiency", "SSSP_Precision" ]

i = 3  #// parse(Int,ENV["SLURM_ARRAY_TASK_ID"])

if i > length(PS_TO_TEST)
    @error "i = $i > length(PS_TO_TEST) = $(length(PS_TO_TEST)). Exciting ..."
    exit()
end

PSi = PS_TO_TEST[i]

##* ===================================================

function grab_alat(PS)
    ff = "$TEST_FOLDER/ReO3_$(PS).pw.x.out"
    a0 = -1.0
    if isfile(ff)
        a0 = pw_alat_bohr(readlines(ff)) * _BOHR_RADIUS_
        if a0 > 100.0 || a0 < 0
            return 7.090709 * _BOHR_RADIUS_
        else
            @info "ici"
            return a0
        end
    else
        return 7.090709 * _BOHR_RADIUS_
    end
end

d2(PS,a0) = Dict(
        "title"           => "ReO3_$PS", 
        "prefix"          => "ReO3_$PS", 
        :calc             => "scf",
        :pseudo_mode      => PS,
        "mixing_beta"     => 0.8,
        "electron_maxstep"=> 500,
        "conv_thr"        => 1e-11,
        "ecutwfc"         => 100.0, 
        "ecutrho"         => 400.0,
        "nbnd"            => 42,
        "disk_io"         => "high",
        
        :cif              => (a0,a0,a0, 90,90,90, 221),
        :do_not_use_symmetry => true,
        :pw_mode          => "energy",
        :kpoint_mode      => "automatic", 
        :kpoints          => (2,2,2,1,1,1),
        :updater          => upd1,
        :watchdog_setting => (woof_per_x_min=1/20, max_watch=24400, max_silence_woofs=30, tail_length=20, quiet=false),
)


function upd1(x)
    d = dict__pw_result(x)
    return (d ↓ ["ibrav", :cif, :reciprocal_basis]) ⬱ Dict("ibrav"=>0,:do_not_use_symmetry=>true)
end


instr(PS,dx,dy,dz) = INSTRUCTION(
    Dict( `pw.x` => _electron_scf_seed_conf_  ⬱ d2(PS,grab_alat(PS)) ),
    [
        (   ".", 
            `pw.x`, 
            Dict( "positions"  => [("Re", 0.0, 0.0, 0.0), ("O", 1/2+dx, dy, dz), ("O", 0, 1/2, 0), ("O", 0, 0, 1/2)], )
        ),
    ]
)

##


instr1(PS,dx,dy,dz) = INSTRUCTION(
    Dict( `pw.x` => _electron_scf_seed_conf_  ⬱ d2(PS,grab_alat(PS)) ),
    [
        (   ".", 
            `pw.x`, 
            Dict("startingpot" => "file",     #* set to 'file' when MC 
                 #"startingwfc" => "file",     #* set to 'file' when MC 
                 "positions"   => [("Re", 0.0, 0.0, 0.0), ("O", 1/2+dx, dy, dz), ("O", 0, 1/2, 0), ("O", 0, 0, 1/2)], )
        ),
    ]
)


##* ===================================================

R = execute_serial(instr(PSi,0.0,0.0,0.0), WORKSPACE=TEST_FOLDER)
R = execute_serial(instr1(PSi,0.01,0.005,0.002), WORKSPACE=TEST_FOLDER)

##* ===================================================
