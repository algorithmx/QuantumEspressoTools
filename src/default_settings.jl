
#TODO finish all defaults

global const findsym_equivalent_settings_QE_default = []

convert_findsym_setting_to_QE(findsym_setting::Vector) = Dict()

convert_setting_to_findsym(QE_setting::Dict) = [
    "monoclinicAxes" => (QE_setting["uniqueb"] ? "a(b)c" : "ab(c)"),
    #* uniqueb default is FALSE, the two fold axis or the mirror normal is parallel to the c axis.
    #Axes of a monoclinic space group.  Enter a(b)c, c(-b)a, ab(c), ba(-c), (a)bc, 
    #or (-a)cb.  The axis in parentheses is the unique axis. 
    #The default value is a(b)c.
    "monoclinicCell" => "1",
    # Cell choice of a monoclinic space group.  Enter 1, 2, or 3.  The default value is 1.
    "orthorhombicAxes" => "abc",
    #Axes of an orthorhombic space group.  Enter abc, ba-c, cab, -cba, bca, or a-cb.
    #The default value is abc.
    "originChoice" => string(QE_setting["origin_choice"]),
    #* origin_choice default is 1
    #Origin choice for a space group.  Enter 1 or 2.  The default value is 2.
    #"hexagonalAxes" => "",
    #Use hexagonal axes for R-centered hexagonal space groups.  This is the default value.
    (QE_setting["rhombohedral"] ? ("rhombohedralAxes" => "") : ("hexagonalAxes" => "")),
    #* rhombohedral default is TRUE
]

convert_setting_to_findsym(QE_setting::NT) where {NT<:NamedTuple} = convert_setting_to_findsym(
    Dict( "uniqueb"=>QE_setting[:uniqueb],
          "rhombohedral"=>QE_setting[:rhombohedral], 
          "origin_choice"=>QE_setting[:origin_choice]  )
)

convert_setting_to_findsym(QE_setting::Vector{P}) where {P<:Pair} = convert_setting_to_findsym(Dict(QE_setting))

global const QE_default_symmetry_group_convention = (uniqueb=false,rhombohedral=true,origin_choice=1)

global const QE_default_equivalent_settings_findsym = convert_setting_to_findsym(QE_default_symmetry_group_convention)


#* ============================================================
