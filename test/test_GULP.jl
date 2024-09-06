include("/home/dabajabaza/jianguoyun/Workspace/QuantumEspressoTools/scripts/scripts.jl")

#

fn_eig  = "/home/dabajabaza/Downloads/ScF3_GBRV_LDA.eig"
#fn_eig  = "/home/dabajabaza/Downloads/ScF3_PSL_PAW_PBE_SR.eig"

avail_modes = parse_mode_file(readlines(fn_eig));

avail_q = unique(first.(avail_modes)) ;

fit_q = Vector{Float64}[
    [0,0,0], [0.1,0,0], [0.2,0,0], [0.3,0,0], [0.4,0,0],
    [1/2,0,0], [1/2,0.1,0], [1/2,0.2,0], [1/2,0.3,0], [1/2,0.4,0],
    [1/2,1/2,0], [1/2,1/2,0.1], [1/2,1/2,0.2], [1/2,1/2,0.3], [1/2,1/2,0.4],
    [1/2,1/2,1/2], [0.4,0.4,0.4], [0.3,0.3,0.3], [0.2,0.2,0.2], [0.1,0.1,0.1] 
]

mode_sel = Dict(q=>[1,2,] for q ∈ fit_q)

qq, ff = select_fitting_frequencies(avail_modes, mode_sel) ;

#fit_q1 = Vector{Float64}[
#    [0.2,0,0],
#    [1/2,0.3,0], 
#    [1/2,1/2,0.2],]

#mode_sel1 = Dict(q=>[1,2,] for q ∈ fit_q1)

#m = select_fitting_modes( 4, parse_mode_file(readlines(fn_eig)), mode_sel; atom_order = [1,2,3,4])

##

conf(D0, a0, r0, K90, K180) = Dict(
    :title       => "ScF3 fit phonon frequencies",
    :control_line=> "fit conv phon include_imaginary",
    :cell        => (3.9647544, 3.9647544, 3.9647544, 90.00, 90.00, 90.00, 0, 0, 0),
    :IT          => 221,
    :fractional  => ["Sc"=>[0,0,0], "F"=>[1/2,0,0],],  # only symmetry inequivalent points
    :species     => ["Sc"=> 3.0,    "F"=> -1.0,    ],
    :observables => ["frequency  $(length(ff))"=>ff],
    :specific    => [("shrink", 10, 10, 10) , "kpoints", qq...],
    :potentials  => Dict(("Sc","F")    => ["morse",  D0,  a0, r0,  2.1,  1, 1, 1],
                         ("Sc","F","F")  => ["three",  K90,   90.0,  2.1,  2.1,  2.1*sqrt(2), 1, 0],
                         ("F","Sc","Sc") => ["lin3",   K180,   1,  1,   2.1,  2.1,  4.2,  1] 
                        ),
    :qpoints     => [[0,0,0]=>20, [1/2,0,0]=>20, [1/2,1/2,0]=>20, [1/2,1/2,1/2]=>30,[0,0,0]=>1]
)

##

params = [  conf(D0, a0, r0, K90, K180) 
            for D0 ∈ [0.1,10.0,100.0]
            for a0 ∈ [0.1,2.0,]
            for r0 ∈ [1.5,2.5]
            for K90 ∈ [0.1,10.0]
            for K180 ∈ [0.1,10.0]
] ;

##

for (i,c) ∈ enumerate(params)
    gulp_fit_phonon_input(c) ⇶ "/data/ScF3_gulp/c/ScF3.fitmode.$(i).gin" ;
end
N
##
#=
--------------------------------------------------------------------------------
 Parameter No.    Parameter        Parameter      Parameter Type  Species
                  Original         Final                                 
--------------------------------------------------------------------------------
         1              0.001000         0.001000 Three-body cnst     1
         2              0.015000         0.901569 Three-body cnst     2
         3              3.000000         2.615992 Morse  De      
         4              0.900000         0.861038 Morse  a0      
         5              1.000000         2.623883 Morse  r0      
--------------------------------------------------------------------------------

--------------------------------------------------------------------------------
 Parameter No.    Parameter        Parameter      Parameter Type  Species
                  Original         Final                                 
--------------------------------------------------------------------------------
         1              1.354801         1.381638 Three-body cnst     1
         2              4.169320         4.151799 Morse  De      
         3              1.155725         1.153397 Morse  a0      
--------------------------------------------------------------------------------


--------------------------------------------------------------------------------
 Parameter No.    Parameter        Parameter      Parameter Type  Species
                  Original         Final                                 
--------------------------------------------------------------------------------
         1              1.354801         1.586537 Three-body cnst     1
         2              4.169320         3.862526 Morse  De      
--------------------------------------------------------------------------------


--------------------------------------------------------------------------------
 Parameter No.    Parameter        Parameter      Parameter Type  Species
                  Original         Final                                 
--------------------------------------------------------------------------------
         1              1.354801         1.777698 Three-body cnst     1
         2              4.169320         3.611470 Morse  De      
--------------------------------------------------------------------------------

--------------------------------------------------------------------------------
 Parameter No.    Parameter        Parameter      Parameter Type  Species
                  Original         Final                                 
--------------------------------------------------------------------------------
         1              1.354801         1.162070 Three-body cnst     1
         2              4.169320         4.431529 Morse  De      
--------------------------------------------------------------------------------

--------------------------------------------------------------------------------
 Parameter No.    Parameter        Parameter      Parameter Type  Species
                  Original         Final                                 
--------------------------------------------------------------------------------
         1              1.354801         1.147452 Three-body cnst     1
         2              4.169320         4.428382 Morse  De      
--------------------------------------------------------------------------------



=#