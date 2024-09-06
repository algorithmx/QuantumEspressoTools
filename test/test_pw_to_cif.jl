using Plots
gr()

#global const WKSPC = "/mnt/dg_hpc/STRAIN2/relaxed_1.0_1.0/ps___PSL_USPP_PBESOL_SR___C-n-1.0.0___kp_30,24,1,0,0,0_cutoff_60.0,480.0_kBar_0.0_dg_0.02"
global const WKSPC = "/mnt/dg_hpc/STRAIN2/relaxed_1.05_1.0/ps___PSL_USPP_PBESOL_SR___C-n-1.0.0___kp_20,16,1,0,0,0_cutoff_60.0,480.0_kBar_0.0_dg_0.02"
ENV["PSEUDO_ROOT"] = "/home/dabajabaza/abinitio/pseudopotentials"
global const SFT = "/home/dabajabaza/jianguoyun/Nutstore"
include("$SFT/QuantumEspressoTools/scripts/scripts.jl")

##

@inline check_fd(fd) = (length([f for f in readdir(fd) if occursin("fc_relax_iter",f)])>=2
                    &&  length([f for f in readdir(fd) if endswith(f,".cif")])>=2)

@inline last_cif(fd) = last(sort([f for f in readdir(fd) if endswith(f,".cif")]))

##

ciflast = last_cif(WKSPC)

STRCT  = get_structure_from_cif("$WKSPC/$ciflast")

POS0   = sort(STRCT["positions"])

regularize_structure(POS0)

##

function regularize_structure_old(pos)
    POS0   = sort(pos)
    x1  = 0.5*(POS0[1][2]+POS0[2][2])
    y1  = 0.5*(POS0[1][3]+POS0[5][3])
    y3  = 0.5*(POS0[3][3]+1-POS0[4][3])
    POS = [
        ("C", x1,     y1, 0.5),
        ("C", x1,   1-y1, 0.5),
        ("C", 0.5,    y3, 0.5),
        ("C", 0.5,  1-y3, 0.5),
        ("C", 1-x1,   y1, 0.5),
        ("C", 1-x1, 1-y1, 0.5),
    ]
    return POS
end

function regularize_structure(pos)
    POS0 = sort(pos)
    @assert POS0[3][2]≈0.5
    @assert POS0[4][2]≈0.5
    @inline sm(x) = (x<0.5) ? x : 1-x
    @inline rd(x) = round(x,digits=7)
    x1  = 0.25*sum(sm.([POS0[1][2],POS0[2][2],POS0[5][2],POS0[6][2]])) |> rd
    y1  = 0.25*sum(sm.([POS0[1][3],POS0[2][3],POS0[5][3],POS0[6][3]])) |> rd
    y3  = 0.5*sum(sm.([POS0[3][3],POS0[4][3]])) |> rd
    POS = [
        ("C", x1,     y1, 0.5),
        ("C", x1,   1-y1, 0.5),
        ("C", 0.5,    y3, 0.5),
        ("C", 0.5,  1-y3, 0.5),
        ("C", 1-x1,   y1, 0.5),
        ("C", 1-x1, 1-y1, 0.5),
    ]
    return POS
end

##


##


##

function pw_to_cif(result_lines, fn)
    try
        latt_p   = QuantumEspressoTools.pw_relax_to_lattice_parameters(result_lines)
        atm_list = QuantumEspressoTools.pw_relax_atom_list(result_lines)
        cif = QuantumEspressoTools.minimal_cif("generated_by_pw_to_cif", latt_p, atm_list)
        QuantumEspressoTools.write_to_file(cif,fn)
    catch _e_
        @error "pw_to_cif() failed : error = $(_e_)"
    end
    return Dict()
end

##

pw_out = readlines("$(WKSPC)/strain_rx_1.0_ry_1.0_fc_relax_iter_1/strain_rx_1.0_ry_1.0.pw.x.out")

pw_to_cif(pw_out,"test.cif")


pw_out1 = readlines("$(WKSPC)/strain_00_fc_relax_iter_1/strain_00.pw.x.out")
pw_out2 = readlines("$(WKSPC)/strain_00_fc_relax_iter_2/strain_00.pw.x.out")
pw_out3 = readlines("$(WKSPC)/strain_00_fc_relax_iter_3/strain_00.pw.x.out")
pw_out4 = readlines("$(WKSPC)/strain_00_fc_relax_iter_4/strain_00.pw.x.out")


##

pw_relax_result_to_cif_lines(pw_out3)

##



##


global const WKSPC = "/mnt/dg_hpc/WS2"
ENV["PSEUDO_ROOT"] = "/home/dabajabaza/abinitio/pseudopotentials"
global const SFT = "/home/dabajabaza/jianguoyun/Nutstore"
include("$SFT/QuantumEspressoTools/scripts/scripts.jl")

##

get_structure_from_cif("$WKSPC/WS2.0.cif")


