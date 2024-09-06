# pp_input_file.jl

#: -----------------------------------------------------------------------------
    #:Variable:       output_format
    #Type:           INTEGER
    #Description:    (ignored on 1D plot)
    #                
    #                0  = format suitable for gnuplot   (1D)
    #                1  = obsolete format no longer supported
    #                2  = format suitable for plotrho   (2D)
    #                3  = format suitable for XCRYSDEN  (2D or user-supplied 3D region)
    #                4  = obsolete format no longer supported
    #                5  = format suitable for XCRYSDEN  (3D, using entire FFT grid)
    #                6  = format as gaussian cube file  (3D)
    #                     (can# be read by many programs)
    #                7  = format suitable for gnuplot   (2D) x, y, f(x,y)
#: -----------------------------------------------------------------------------



# -----------------------------------------------------------------------------
#:Variable:   plot_num    INTEGER
    # Selects what to save in filplot:

    #>   0  = electron (pseudo-)charge density

    #   1  = total potential V_bare + V_H + V_xc

    #   2  = local ionic potential V_bare

    #   3  = local density of states at specific energy or grid of energies
    #        (number of states per volume, in bohr^3, per energy unit, in Ry)

    #   4  = local density of electronic entropy

    #   5  = STM images
    #        Tersoff and Hamann, PRB 31, 805 (1985)

    #   6  = spin polarization (rho(up)-rho(down))

    #   7  = contribution of selected wavefunction(s) to the
    #        (pseudo-)charge density. For norm-conserving PPs,
    #        |psi|^2 (psi=selected wavefunction). Noncollinear case:
    #        contribution of the given state to the charge or
    #        to the magnetization along the direction indicated
    #        by spin_component (0 = charge, 1 = x, 2 = y, 3 = z )

    #   8  = electron localization function (ELF)

    #   9  = charge density minus superposition of atomic densities

    #   10 = integrated local density of states (ILDOS)
    #        from emin to emax (emin, emax in eV)
    #        if emax is not specified, emax=E_fermi

    #   11 = the V_bare + V_H potential

    #   12 = the sawtooth electric field potential (if present)

    #   13 = the noncollinear magnetization.

    #   17 = all-electron valence charge density
    #        can be performed for PAW calculations only
    #        requires a very dense real-space grid!

    #   18 = The exchange and correlation magnetic field in the noncollinear case

    #   19 = Reduced density gradient
    #        ( J. Chem. Theory Comput. 7, 625 (2011), doi:10.1021/ct100641a )
    #        Set the isosurface between 0.3 and 0.6 to plot the
    #        non-covalent interactions (see also plot_num = 20)

    #   20 = Product of the electron density (charge) and the second
    #        eigenvalue of the electron-density Hessian matrix;
    #        used to colorize the RDG plot (plot_num = 19)

    #   21 = all-electron charge density (valence+core).
    #        For PAW calculations only; requires a very dense real-space grid.

    #   22 = kinetic energy density (for meta-GGA and XDM only)
# -----------------------------------------------------------------------------

global const _pp_plot_num_ = Dict(
    :total_charge => 0,
)

global const _pp_output_format_ = Dict(
    :total_charge => 7,
)

function INPUTPP(config::Dict)
    NL = Namelist("INPUTPP",Dict(),Dict(),Dict())
    if config[:pp_mode]==:total_charge
        prefix = config["prefix"]
        NL.req  = Dict(
                "prefix"        => prefix,
                "plot_num"      => _pp_plot_num_[:total_charge],
                "filplot"       => "$(prefix).charge"
        )
        NL.default = Dict(
                "outdir"   => "./",
        )
        NL.opt = NL.default ← config
        return NL
    else
        return NL
    end
end


function PLOT(config::Dict)
    #TODO
    #cmds2 = [
    #    join( [ "e1(1) = $(basis[1][1]), e1(2) = $(basis[1][2]), e1(3) = $(basis[1][3])"*(dimension==1 ? "" : ","),
    #            (dimension>1) ? "e2(1) = $(basis[2][1]), e2(2) = $(basis[2][2]), e2(3) = $(basis[2][3])"*(dimension==2 ? "" : ",") : "",
    #            (dimension>2) ? "e3(1) = $(basis[3][1]), e3(2) = $(basis[3][2]), e3(3) = $(basis[3][3])"*(dimension==3 ? "" : ",") : "" 
    #          ], " "),
    #    "x0(1) = $(origin[1]), x0(2) = $(origin[2]), x0(3) = $(origin[3])",
    #    join( [ "nx = $(Ngrid[1])"*(dimension==1 ? "" : ","),
    #            (dimension>1) ? "ny = $(Ngrid[2])"*(dimension==2 ? "" : ",") : "",
    #            (dimension>2) ? "nz = $(Ngrid[3])"*(dimension==3 ? "" : ",") : "" 
    #          ], " "),
    #    "/",
    #]
    NL = Namelist("PLOT",Dict(),Dict(),Dict())
    if config[:pp_mode]==:total_charge
        NL.req  = Dict(
                "iflag"    => config["dimension"],
                "fileout"  => config["fileout"],
                "output_format" => _pp_output_format_[:total_charge],
                )
        NL.default = Dict(
        )
        NL.opt = NL.default ← config
        return NL
    else
        return NL
    end

end


# this is a condensation of experience
function pp_input_plot_total_charge(config::Dict)
    @assert config[:pp_mode]==:total_charge
    return  INPUT_PPX(  config["title"],
                        map(x->x(config), [INPUTPP, PLOT]) ... ) |> build_input_file

end
