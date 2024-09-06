householder_make_v2_perp_v1(v1,v2) = (v2 .- dot(v2,normalize(v1)).*normalize(v1))

function plane_basis_from_normal_vector(n::Vector, xref::Vector)
    z = normalize(n)
    x = normalize(householder_make_v2_perp_v1(z,xref))
    y = normalize(cross(z,x))
    return x, y
end

