function random_kick_positions(dic0::Dict, delta::Float64, elems=[])
    L_is_0 = (length(elems)==0)
    @inline rd3(x) = (delta*Int(L_is_0 || x)).*(2rand(3).-1)
    dic1 = copy(dic0)
    dic1["positions"] = [(atm, ([x,y,z].+rd3(atm∈elems))...) 
                         for (atm,x,y,z) ∈ dic0["positions"]]
    return dic1
end


function random_kick_positions_masked(
    dic0::Dict, 
    delta::Float64, 
    mask::Vector{Vector{Int}}, 
    elems=[]
    )
    @assert length(dic0["positions"]) == length(mask)
    @assert all(length.(mask).==3)
    # elems=[] ==> L_is_0, which means update all elements
    L_is_0 = (length(elems)==0) 
    @inline rd3(x) = (delta*Int(L_is_0 || x)).*(2rand(3).-1)
    dic1 = copy(dic0)
    dic1["positions"] = [(atm, ([x,y,z] .+ mv .* rd3(atm∈elems))..., mv...) 
                         for ((atm,x,y,z),mv) ∈ zip(dic0["positions"],mask)]
    return dic1
end


function random_kick_positions_masked(
    dic0::Dict, 
    delta::Float64, 
    mask_all::Vector{Int}, 
    elems=[]
    )
    random_kick_positions_masked(
        dic0, 
        delta, 
        [mask_all for i=1:length(dic0["positions"])], 
        elems
    )
end


random_kick_positions_2D(dic0::Dict, delta::Float64, elems=[]) = 
    random_kick_positions_masked(dic0, delta, Int[1,1,0], elems)
    

# enforce [1,1,0]
random_kick_positions_2D_masked(dic0::Dict, delta::Float64, mask::Vector{Vector{Int}}, elems=[]) = 
    random_kick_positions_masked(
        dic0, 
        delta, 
        Vector{Int}[Int[1,1,0].*mv for mv ∈ mask], 
        elems
    )


# enforce [1,1,0]
random_kick_positions_2D_masked(dic0::Dict, delta::Float64, mask_all::Vector{Int}, elems=[]) = 
    random_kick_positions_masked(dic0, delta, Int[1,1,0] .* mask_all, elems)
    


#: ==================================================================


function random_kick_lattp(dic0::Dict, delta::Float64, w6=[1,1,1,0,0,0])
    (a, b, c, α, β, γ) = dic0[:cif][1:6]
    @inline rd(x) = round(x,digits=9)
    @inline RRR() = delta*(2rand().-1)
    dic1 = copy(dic0)
    dic1[:cif] = (a+w6[1]*RRR(), b+w6[2]*RRR(), c+w6[3]*RRR(), 
                  α+w6[4]*RRR(), β+w6[5]*RRR(), γ+w6[6]*RRR(), 1)
    basis1 = lattice_parameters_to_basis((dic1[:cif][1:6]...,))
    dic1[:cell_parameters] = [rd.(basis1[1]),rd.(basis1[2]),rd.(basis1[3])]
    return dic1
end


function relax_update_random_kick(x,delta=0.01,elems=[])
    d = dict__pw_relax_result(x)
    d = random_kick_positions(d, delta, elems)
    d = random_kick_lattp(d, delta)
    d1 = (d ↓ ["ibrav", :cif, :reciprocal_basis]) ⬱ Dict("ibrav"=>0,:do_not_use_symmetry=>true)
    return d1
end


function relax_update_random_kick_2D(x,delta=0.01,elems=[])
    d = dict__pw_relax_result(x)
    d = random_kick_positions_2D(d, delta, elems)
    d = random_kick_lattp(d, delta, [1,1,0,0,0,0])
    d1 = (d ↓ ["ibrav", :cif, :reciprocal_basis]) ⬱ Dict("ibrav"=>0,:do_not_use_symmetry=>true)
    return d1
end


function relax_update_random_kick_pos(x,delta=0.01,elems=[])
    d = dict__pw_relax_result(x)
    d = random_kick_positions(d, delta, elems)
    d1 = (d ↓ ["ibrav", :cif, :reciprocal_basis]) ⬱ Dict("ibrav"=>0,:do_not_use_symmetry=>true)
    return d1
end


function relax_update_random_kick_pos_2D(x,delta=0.01,elems=[])
    d = dict__pw_relax_result(x)
    d = random_kick_positions_2D(d, delta, elems)
    d1 = (d ↓ ["ibrav", :cif, :reciprocal_basis]) ⬱ Dict("ibrav"=>0,:do_not_use_symmetry=>true)
    return d1
end


function relax_update_random_kick_pos_2D_masked(
    x::Vector{String},
    delta::Float64,
    mask::Vector{Vector{Int}},
    elems::Vector{String} = String[]
    )
    d = dict__pw_relax_result(x)
    d = random_kick_positions_2D_masked(d, delta, mask, elems)
    d1 = (d ↓ ["ibrav", :cif, :reciprocal_basis]) ⬱ Dict("ibrav"=>0,:do_not_use_symmetry=>true)
    return d1
end

