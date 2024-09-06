function norm_basis(basis3::Vector)
    #@assert eltype(basis3[1]) <: Real && eltype(basis3[2]) <: Real && eltype(basis3[3]) <: Real
    alat = norm(basis3[1])
    return alat, (1/alat).*basis3[1], (1/alat).*basis3[2], (1/alat).*basis3[3]
end


function regulate_a1_M(a1::Vector{Float64})
    (a,b,c) = BigFloat.(a1)
    if abs(b^2+c^2) < 1e-10  return Float64[1 0 0; 0 1 0; 0 0 1]  end
    n = norm([a,b,c])
    #=
    # Mathematica
    M = FullSimplify[RotationMatrix[-ArcCos[Dot[{1,0,0},{a,b,c}]/Sqrt[Dot[{a,b,c},{a,b,c}]]],Normalize[Cross[{1,0,0},{a,b,c}]]],Assumptions->a>0&&b\[Element]Reals&&c\[Element]Reals];
    M/.{Sqrt[a^2+b^2+c^2]->n,1/Sqrt[a^2+b^2+c^2]->1/n}
    FullSimplify[M.{a, b, c}, Assumptions -> a > 0 && b \[Element] Reals && c \[Element] Reals]
    =#
    M = [ a/n  b/n  c/n;  
          (-b/n)  ((c^2+(a * b^2)/n)/(b^2+c^2))  ((b * c * (-1 + a/n))/(b^2 + c^2));  
          (-c/n)  ((b * c * (-1 + a/n))/(b^2 + c^2))  ((b^2 + (a * c^2)/n)/(b^2 + c^2)) ]
    return M
end


function regulate_a1(basis::Vector{Vector{Float64}})
    M = regulate_a1_M(basis[1])
    return Vector{Float64}.((M*basis[1], M*basis[2], M*basis[3]))
end
# test
# basis1=regulate_a1((basis = [rand(3), rand(3), rand(3)]; basis)); det(hcat(basis...))-det(hcat(basis1...))
