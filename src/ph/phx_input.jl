# phx_input_file.jl

#IMPORTANT overall convention for dictionary keys : Symbol for internal control, String to pass over

PH_Float_Int = Union{Float64, Int64}

function k_mode(mode)
    # TODO mode == "kpath"
    #+ kpath driver: (1) Python seekkpath ; (2) ISOTROPY

    ldisp=true
    qplot=true
    q_in_band_form=false
    q2d=false

    if mode == :multiple
        ldisp=true
        qplot=true
        q_in_band_form=false
        q2d=false
    elseif mode == :single
        ldisp=false
        qplot=false
        q_in_band_form=false
        q2d=false
    elseif mode == :grid
        ldisp=true
        qplot=false
        q_in_band_form=false
        q2d=false
    elseif mode == :bands
        ldisp=true
        qplot=true
        q_in_band_form=true
        q2d=false
    else # default to  mode == "multiple"
        ldisp=true
        qplot=true
        q_in_band_form=false
        q2d=false
    end

    return Dict(
        "ldisp"              => ldisp,
        "qplot"              => qplot,
        "q_in_band_form"     => q_in_band_form,
        "q2d"                => q2d
    )

end


function INPUTPH(config)
    NL = Namelist("INPUTPH", Dict(), Dict(), Dict())
    NL.req = Dict(
        "prefix"            => config["prefix"],
        "outdir"            => config["outdir"],
        #//"tmp_dir"           => config["outdir"],
        "fildyn"            => "$(config["prefix"]).dynmat.",  #* the dot makes life easier
    ) ⬱ k_mode(config[:ph_mode]) #* (ldisp, qplot, q_in_band_form, q2d)
    NL.default = Dict(
        "search_sym"        => false,
        "verbosity"         => "high",
        "tr2_ph"            => 1e-14,
        "alpha_mix"         => 0.7,
        "nmix_ph"           => 4,
        "diagonalization"   => "david",
        "nq1"               => 1,
        "nq2"               => 1,
        "nq3"               => 1,
        "mass"              => nothing,
    )
    # config[k] can be "nothing" or ":default"
    NL.opt = NL.default ← config
    return NL
end


function QSECTION(config)

    CD = Card{PH_Float_Int}("", "", Vector{PH_Float_Int}[])

    if config[:ph_mode] == :single
        CD.contents = Vector{PH_Float_Int}[PH_Float_Int[config[:qvecs][1]...],]
    elseif config[:ph_mode] == :multiple
        qv = [PH_Float_Int[q..., 1] for q in config[:qvecs]]
        CD.contents = Vector{PH_Float_Int}[[PH_Float_Int[length(config[:qvecs]),],]; qv]
    elseif config[:ph_mode] == :bands
        #* difference compared to "multiple" is that q is a pair: (qx,qy,qz) => nq
        qv = [PH_Float_Int[q[1]..., q[2]] for q in config[:qvecs]]
        CD.contents = Vector{PH_Float_Int}[[PH_Float_Int[length(config[:qvecs]),],]; qv]
    elseif config[:ph_mode] == :grid
        nothing
    else
        throw(error("QSECTION():\nmode=$(config[:ph_mode]) unsupported.")) 
    end
    return CD
end

_check_input_X(X::Symbol, config::Dict) = ( (X==:single   && true) ||
                                            (X==:multiple && true) ||
                                            (X==:grid     && true) ||
                                            (X==:bands    && true)    )

function phx_inupt_X(X::Symbol, config::Dict)
    @assert  X == config[:ph_mode]
    #TODO check input
    @assert  _check_input_X(X, config)
    return  INPUT_PHX(config["title"], map(x->x(config), [INPUTPH, QSECTION])...) |> build_input_file
end

##* ================  qrel_to_qabs  ================

#: ph.x accept qabs in units of 2pi/a0 (a0 = lattice parameter in a.u.)
#: xq1, xq2, xq3 are the q-point coordinates; used only with ldisp=.true. and qplot=.true.
#: it is easier to work with qrel
#: cross-ref : dict__pw_result(pw_out_fn), CELL_PARAMETERS(config::Dict)

@inline qrel_to_qabs(qrel, b1_2pi_over_a0::Vector, b2_2pi_over_a0::Vector, b3_2pi_over_a0::Vector) = (
    qrel[1].*b1_2pi_over_a0 .+ qrel[2].*b2_2pi_over_a0 .+ qrel[3].*b3_2pi_over_a0
)

@inline qrel_to_qabs(qrel, b::Vector{Vector{Float64}}) = qrel_to_qabs(qrel,b[1],b[2],b[3])

@inline qrel_to_qabs(qrel, b::Array{Float64,2}) = qrel_to_qabs(qrel,b[:,1],b[:,2],b[:,3])

##* ================================================


@inline phx_inupt_multiple_single_check(c::Dict) = (
    c[:qpoints] isa Vector) && (eltype(c[:qpoints])<:Vector) && all([eltype(v)<:Real for v ∈ c[:qpoints]]
)

function phx_inupt_multiple(config::Dict)
    @assert config[:ph_mode] == :multiple
    @assert phx_inupt_multiple_single_check(config)
    config1 = config ⬱ ( :qvecs => [qrel_to_qabs(qrel,config[:reciprocal_basis]) for qrel ∈ config[:qpoints]] )
    phx_inupt_X(:multiple, config1)
end


function phx_inupt_single(config::Dict)
    @assert config[:ph_mode] == :single
    @assert phx_inupt_multiple_single_check(config)
    config1 = config ⬱ ( :qvecs => [qrel_to_qabs(qrel,config[:reciprocal_basis]) for qrel ∈ config[:qpoints]] )
    phx_inupt_X(:single, config1)
end


function phx_inupt_grid(config::Dict)
    @assert config[:ph_mode] == :grid
    @assert (config[:qpoints] isa Tuple) && (eltype(config[:qpoints])<:Integer)
    config1 = config ⬱ Dict( "nq1" => config[:qpoints][1], 
                              "nq2" => config[:qpoints][2], 
                              "nq3" => config[:qpoints][3],
                              :qvecs => []  )
    phx_inupt_X(:grid, config1)
end


@inline phx_inupt_bands_check(c::Dict) = (
    c[:qpoints] isa Vector) && (eltype(c[:qpoints])<:Pair) && all([eltype(v[1])<:Real && (v[2] isa Integer) for v ∈ c[:qpoints]]
)


function phx_inupt_bands(config::Dict)
    @assert config[:ph_mode] == :bands
    @assert phx_inupt_bands_check(config)
    config1 = config  ⬱ ( :qvecs => [qrel_to_qabs(v,config[:reciprocal_basis])=>n for (v,n) in config[:qpoints]] )
    phx_inupt_X(:bands, config1)
end


##* ======================== TEST ========================

#=

include("src/namelists_cards.jl")
include("src/parsers.jl")
include("src/dict_merger.jl")
include("src/data.jl")
include("src/crystallography.jl")

##

dict_test = Dict( "title"             => "test", 
                  "prefix"            => "test", 
                  :ph_mode               => :bands, 
                  "outdir"     => "./", 
                  :reciprocal_basis  => [[1.0,0.0,0.0], [0.0,1.0,0.0], [0.0,0.0,1.0]], 
                  :qpoints           => [(0.0,0.0,0.0)=>10], 
)

c = phx_inupt_bands(dict_test) ;
println.(c) ;

##

dict_test = Dict( "title"             => "test", 
                  "prefix"            => "test", 
                  :ph_mode               => :grid, 
                  "outdir"     => "./", 
                  :reciprocal_basis  => [[1.0,0.0,0.0], [0.0,1.0,0.0], [0.0,0.0,1.0]], 
                  :qpoints           => (2,2,2), 
)

c = phx_inupt_grid(dict_test) ;
println.(c) ;

##

dict_test = Dict( "title"             => "test", 
                  "prefix"            => "test", 
                  :ph_mode               => :multiple, 
                  "outdir"     => "./", 
                  :reciprocal_basis  => [[1.0,0.0,0.0], [0.0,1.0,0.0], [0.0,0.0,1.0]], 
                  :qpoints           => [[0.0,0.0,0.0], [0.5,0.5,0.5]], 
)

c = phx_inupt_multiple(dict_test) ;
println.(c) ;

##

dict_test = Dict( "title"             => "test", 
                  "prefix"            => "test", 
                  :ph_mode               => :single, 
                  "outdir"     => "./", 
                  :reciprocal_basis  => [[1.0,0.0,0.0], [0.0,1.0,0.0], [0.0,0.0,1.0]], 
                  :qpoints           => [[0.0,0.0,0.0]], 
)

c = phx_inupt_single(dict_test) ;
println.(c) ;

=#