#>  this program calculates the phonon frequencies for a list of generic
#>  q vectors starting from the interatomic force constants generated
#>  from the dynamical matrices as written by DFPT phonon code through
#>  the companion program q2r
#>
#>  matdyn can generate a supercell of the original cell for mass
#>  approximation calculation. If supercell data are not specified
#>  in input, the unit cell, lattice vectors, atom types and positions
#>  are read from the force constant file
#>
#>  Input cards: namelist &input
#>     flfrc     file produced by q2r containing force constants (needed)
#>               It is the same as in the input of q2r.x (+ the .xml extension
#>               if the dynamical matrices produced by ph.x were in xml
#>               format). No default value: must be specified.
#>      asr      (character) indicates the type of Acoustic Sum Rule imposed
#>               - 'no': no Acoustic Sum Rules imposed (default)
#>               - 'simple':  previous implementation of the asr used
#>                  (3 translational asr imposed by correction of
#>                  the diagonal elements of the force constants matrix)
#>               - 'crystal': 3 translational asr imposed by optimized
#>                  correction of the force constants (projection).
#>               - 'one-dim': 3 translational asr + 1 rotational asr
#>                  imposed by optimized correction of the force constants
#>                  (the rotation axis is the direction of periodicity;
#>                   it will work only if this axis considered is one of
#>                   the cartesian axis).
#>               - 'zero-dim': 3 translational asr + 3 rotational asr
#>                  imposed by optimized correction of the force constants
#>               Note that in certain cases, not all the rotational asr
#>               can be applied (e.g. if there are only 2 atoms in a
#>               molecule or if all the atoms are aligned, etc.).
#>               In these cases the supplementary asr are cancelled
#>               during the orthonormalization procedure (see below).
#>     dos       if .true. calculate phonon Density of States (DOS)
#>               using tetrahedra and a uniform q-point grid (see below)
#>               NB: may not work properly in noncubic materials
#>               if .false. calculate phonon bands from the list of q-points
#>               supplied in input (default)
#>     nk1,nk2,nk3  uniform q-point grid for DOS calculation (includes q=0)
#>                  (must be specified if dos=.true., ignored otherwise)
#>     deltaE    energy step, in cm^(-1), for DOS calculation: from min
#>               to max phonon energy (default: 1 cm^(-1) if ndos, see
#>               below, is not specified)
#>     ndos      number of energy steps for DOS calculations
#>               (default: calculated from deltaE if not specified)
#>     fldos     output file for dos (default: 'matdyn.dos')
#>               the dos is in states/cm(-1) plotted vs omega in cm(-1)
#>               and is normalised to 3*nat, i.e. the number of phonons
#>     flfrq     output file for frequencies (default: 'matdyn.freq')
#>     flvec     output file for normalized phonon displacements 
#>               (default: 'matdyn.modes'). The normalized phonon displacements
#>               are the eigenvectors divided by the square root of the mass,
#>               then normalized. As such they are not orthogonal.
#>     fleig     output file for phonon eigenvectors (default: 'matdyn.eig')
#>               The phonon eigenvectors are the eigenvectors of the dynamical
#>               matrix. They are orthogonal.
#>     fldyn     output file for dynamical matrix (default: ' ' i.e. not written)
#>     at        supercell lattice vectors - must form a superlattice of the
#>               original lattice (default: use original cell)
#>     l1,l2,l3  supercell lattice vectors are original cell vectors times
#>               l1, l2, l3 respectively (default: 1, ignored if at specified)
#>     ntyp      number of atom types in the supercell (default: ntyp of the
#>               original cell)
#>     amass     masses of atoms in the supercell (a.m.u.), one per atom type
#>               (default: use masses read from file flfrc)
#>     readtau   read  atomic positions of the supercell from input
#>               (used to specify different masses) (default: .false.)
#>     fltau     write atomic positions of the supercell to file "fltau"
#>               (default: fltau=' ', do not write)
#>     la2F      if .true. interpolates also the el-ph coefficients.
#>     q_in_band_form if .true. the q points are given in band form:
#>               Only the first and last point of one or more lines 
#>               are given. See below. (default: .false.).
#>     q_in_cryst_coord if .true. input q points are in crystalline 
#>              coordinates (default: .false.)
#>     eigen_similarity: use similarity of the displacements to order 
#>                       frequencies  (default: .false.)
#>                NB: You cannot use this option with the symmetry
#>                analysis of the modes.
#>     fd         (logical) if .t. the ifc come from the finite displacement calculation
#>     na_ifc     (logical) add non analitic contributions to the interatomic force 
#>                constants if finite displacement method is used (as in Wang et al.
#>                Phys. Rev. B 85, 224303 (2012)) [to be used in conjunction with fd.x]
#>     nosym      if .true., no symmetry and no time reversal are imposed
#>     loto_2d    set to .true. to activate two-dimensional treatment of LO-TO splitting.
#>     loto_disable (logical) if .true. do not apply LO-TO splitting for q=0
#>                  (default: .false.)
#>
#>  if (readtau) atom types and positions in the supercell follow:
#>     (tau(i,na),i=1,3), ityp(na)
#>  IF (q_in_band_form.and..not.dos) THEN
#>     nq   #> number of q points
#>     (q(i,n),i=1,3), nptq   nptq is the number of points between this point
#>                            and the next. These points are automatically
#>                            generated. the q points are given in Cartesian
#>                            coordinates, 2pi/a units (a=lattice parameters)
#>  ELSE, if (.not.dos) :
#>     nq         number of q-points
#>     (q(i,n), i=1,3)    nq q-points in cartesian coordinates, 2pi/a units
#>  If q = 0, the direction qhat (q=>0) for the non-analytic part
#>  is extracted from the sequence of q-points as follows:
#>     qhat = q(n) - q(n-1)  or   qhat = q(n) - q(n+1)
#>  depending on which one is available and nonzero.
#>  For low-symmetry crystals, specify twice q = 0 in the list
#>  if you want to have q = 0 results for two different directions
#>
#> ---------------------------------------------------------------------
#> example
#>  &input
#>    asr='simple'  
#>    amass(1)=26.98
#>    flfrc='Al444.fc'
#>    flfrq='al.freq'
#>    q_in_band_form=.true.
#> /
#> 6
#>  0.000 0.0 0.0   40
#>  1.0   0.0 0.0   20
#>  1.0   0.5 0.0   20
#>  1.0   1.0 0.0   40
#>  0.00  0.0 0.0   40
#>  0.5   0.5 0.5    1
#>

function INPUT_M(config::Dict)
    prefix = config["prefix"]
    NL = Namelist("INPUT",Dict(),Dict(),Dict())
    NL.req = Dict(
        "asr"              => config["acoustic_sum_rule"],          #* sum_rule
        "q_in_cryst_coord" => true,
        "q_in_band_form"   => true,
    )
    NL.default = Dict(
        "flfrc"       => "$prefix.flfrc",   #* IFC_filename
        "flfrq"       => "$prefix.flfrq",   #* frequency file
        "flvec"       => "$prefix.flvec",   #* eigen vector 
        "fleig"       => "$prefix.fleig",   #* eigen value
        "fldyn"       => "$prefix.fldyn",   #* dynamical matrix
    )
    # config[k] can be "nothing"
    NL.opt = NL.default ← config
    return NL
end


function check_matdyn_config(config)
    return ((config[:qpath_rel] isa Vector) 
            && (eltype(config[:qpath_rel])<:Pair) 
            && all([eltype(v[1])<:Real && (v[2] isa Integer) for v ∈ config[:qpath_rel]]))
end


function matdyn_band_input(config::Dict)
    @assert check_matdyn_config(config)
    config1 = config ⬱ Dict(:qvecs=>config[:qpath_rel], :ph_mode=>:bands)
    return  INPUT_MATDYN(config["title"], map(x->x(config1), [INPUT_M, QSECTION])...) |> build_input_file
end


function add_white_space(fn, wkspc)
    if isfile("$wkspc/$fn")
        lines = readlines("$wkspc/$fn")
        ret = []
        nbnd = strtonum( extract([lines[1]], r"nbnd\s*\=\s*\d+,", r"\d+") )
        nks  = strtonum( extract([lines[1]], r"nks\s*\=\s*\d+", r"\d+") )
        if mod(length(lines)-1,nks)!=0
            @warn "add_white_space():\nFile $fn in $(wkspc) has $(length(lines)) lines and $(nks) qpoints. Incompatible!"
            return
        end
        nfrq = (length(lines)-1)÷nks -1
        push!(ret, lines[1])
        for i=1:nks
            push!(ret, lines[1+(i-1)*(nfrq+1)+1])
            for j=1:nfrq
                L = rstrip(lines[1+(i-1)*(nfrq+1)+1+j])
                N = length(L)
                if mod(N,10)!=0
                    @warn "add_white_space():\nFile $fn in $(wkspc) line $(1+(i-1)*(nfrq+1)+1+j) format incompatible!"
                    return
                else
                    L1 = [string(L[10a+1:10a+10]) for a=0:(N÷10-1)]
                    push!(ret, join(L1," "))
                end
            end
        end
        open("$wkspc/$fn", "w") do fout
            write(fout, join(ret,"\n"))
        end
    else
        @warn "add_white_space():\nFile $fn not in $(wkspc)."
        return
    end
end


##* ============================================

TupVF64VF64 = Tuple{Vector{Float64},Vector{Float64}}

TupF64VC64  = Tuple{Float64,Vector{ComplexF64}}

PairVF64VF64 = Pair{Vector{Float64},Vector{Float64}}

function parse_energy_file(
    lines::Vector{String}; 
    SEG=10
    )::Vector{PairVF64VF64}

    #nbnd  = strtonum( extract([lines[1]], r"nbnd\s*\=\s*\d+,", r"\d+") )
    nks   = strtonum( extract([lines[1]], r"nks\s*\=\s*\d+",   r"\d+") )
    if mod(length(lines)-1,nks)!=0
        @warn "parse_energy_file():\nFile $en_file has $(length(lines)) lines and $(nks) kpoints. Incompatible!"
        return []
    end
    nfrq = (length(lines)-1) ÷ nks - 1
    RESULT = PairVF64VF64[]
    for i=1:nks
        q = strtonum.(split(lines[1+(i-1)*(nfrq+1)+1],keepempty=false))
        freqs = Float64[]
        for j=1:nfrq
            L = lines[1+(i-1)*(nfrq+1)+1+j]
            N = length(L)
            #: a signed float in flfrq has SEG=10 char width (including whitespace and sign)
            if mod(N,SEG)!=0
                @warn "parse_energy_file():\nFile $en_file line $(1+(i-1)*(nfrq+1)+1+j) format incompatible!"
                return []
            else
                L1    = [strip(string(L[SEG*a+1:SEG*a+SEG])) for a=0:(N÷SEG-1)]
                freqs = Float64[freqs; Float64[strtonum(l) for l ∈ L1 if l!=""]]
            end
        end
        push!(RESULT, (q=>sort(freqs)))
    end
    return RESULT
end


function parse_mode_file(
    MODELINES::Vector{String}
    )::Vector{Pair{Vector{Float64},Vector{TupF64VC64}}}

    #: step 0. helper functions
    @inline _spln_(s)     = String[strip(x) for x ∈ split(s,"\n",keepempty=false) if length(strip(x))>0 && !startswith(strip(x),"*****")]
    @inline _sections_(s) = _spln_.( split(s,r"\s*\*\s*\n", keepempty=false) )
    @inline extract_q(qblock)  = parsenf(extract(qblock, r"q\s*\=\s*"*num_f_rstr3, num_f_rstr3))
    @inline extract_freq(line) = parsenf(extract([line], r"\s*\[THz\]\s*\=\s*" * num_f_rstr * r"\s*\[cm-1\]", num_f_rstr))
    @inline extract_xyz(line)  = parsenf(extract([line], r"\(\s*" * num_f_rstr6 * r"\s*\)", num_f_rstr6))

    #: step 1. load and section    
    sections = _sections_(join(MODELINES,"\n"))
    @assert length(unique(length.(sections)))==2 "unique(length.(sections)) = $(unique(length.(sections)))"

    #: step 2. parse blocks
    @inline searchQ(x,lst) = findfirst(y->norm(x-y[1])<1e-3, lst)
    MODES = []
    for ib = 1:length(sections)÷2
        qblck = sections[2ib-1]
        fblck = sections[2ib]
        q = extract_q(qblck)
        f_curr = 0.0
        m_curr = Vector{ComplexF64}[]
        freq_mode_tuples = TupF64VC64[]
        for l ∈ fblck
            if startswith(l,"freq")
                if length(m_curr)>0
                    # (f1, [a1 a2 a3 ...])
                    # (f2, [b1 b2 b3 ...])
                    # (f3, [c1 c2 c3 ...])
                    # ...
                    push!( freq_mode_tuples, (f_curr, vcat(m_curr...)) ) #! no transpose!
                end
                f_curr = extract_freq(l) |> first
                m_curr = Vector{ComplexF64}[]
            elseif startswith(l,"(")
                (xr,xi,yr,yi,zr,zi) = extract_xyz(l)
                push!(m_curr, ComplexF64[xr+im*xi,yr+im*yi,zr+im*zi])
            else
                throw(error("Corrupted file $flvec_or_fleig !!!!!"))
            end
        end
        if length(m_curr)>0
            # (f1, [a1 a2 a3 ...])
            # (f2, [b1 b2 b3 ...])
            # (f3, [c1 c2 c3 ...])
            # ... 
            push!( freq_mode_tuples, (f_curr, vcat(m_curr...)) ) #! no transpose!
        end

        ASS = sortperm(first.(freq_mode_tuples))
        push!(MODES , q=>freq_mode_tuples[ASS])
    end

    #: step 3. return
    return MODES
end
