#export conv_xml_dic, conv_xml
#include("utils/xml.jl")


export strtonum, strtonum_fort, parsef, parse2f, parse3f, parse1s3f, evalt
export parse_join, parse_elem
export showdict
include("utils/parsers.jl")


export ↓  ,  ↑  ,  purge  ,  ⬱  ,  ←  ,  ∪
include("utils/dict_merger.jl")


export find_line, find_all_lines, find_last_line, extract, extract_last, ⊂
include("utils/find_extract.jl")


export plane_basis_from_normal_vector, householder_make_v2_perp_v1
export extend_positions
include("utils/geomtry_tools.jl")


export try_mkdir
include("utils/file_io.jl")


include("utils/cif_op.jl")


export watchdog_serial, ∞_kill
include("utils/watchdog.jl")


include("utils/banner.jl")