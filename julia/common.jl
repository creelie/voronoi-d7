# Shared helpers: exact rational polytopes through cddlib (GMP mode).
using Polyhedra, CDDLib
const Q = Rational{BigInt}
const LIB = CDDLib.Library(:exact)

"The 240 roots of E8 in the standard coordinates."
function e8roots()
    R = Vector{Vector{Q}}()
    for i in 1:8, j in i+1:8, a in (1, -1), b in (1, -1)
        v = zeros(Q, 8); v[i] = a; v[j] = b; push!(R, v)
    end
    for bits in 0:255
        s = [((bits >> k) & 1) == 1 ? -1//2 : 1//2 for k in 0:7]
        iseven(count(<(0), s)) && push!(R, Q.(s))
    end
    R
end

"""
A lattice in a rational subspace U of R^N given by a basis matrix B (rows span U).
A point w = B' p of U has chart coordinates p; <w, x> = p ⋅ (B x).
"""
struct Chart
    B::Matrix{Q}
end
grad(C::Chart, x) = C.B * x                       # p ↦ <B'p, x>
lift(C::Chart, p) = transpose(C.B) * p
nrm(C::Chart, p) = (w = lift(C, p); sum(w .^ 2))

"Polytope {p : <B'p, r> ≤ 1 for r in R} ∩ extra half-spaces (x, b): <B'p, x> ≤ b."
function cell(C::Chart, R, extra = Tuple{Vector{Q},Q}[])
    hs = HalfSpace{Q,Vector{Q}}[]
    for r in R
        push!(hs, HalfSpace(grad(C, r), one(Q)))
    end
    for (x, b) in extra
        push!(hs, HalfSpace(grad(C, x), Q(b)))
    end
    polyhedron(hrep(hs), LIB)
end
vertices(P) = collect(points(vrep(P)))

"""
Volume of the cap {z ∈ P : <z, x> ≥ c}, in chart coordinates: the convex hull of the
vertices of P beyond the hyperplane and of the vertices of P ∩ {<z,x> ≤ c} on it.
Returns (cap volume, vertices of P ∩ {<z,x> ≤ c}, section vertices, vertices beyond).
"""
function cap(C::Chart, R, V, x, c)
    Pc = cell(C, R, [(x, Q(c))])
    W = vertices(Pc)
    g = grad(C, x)
    sec = [w for w in W if g ⋅ w == c]
    top = [p for p in V if g ⋅ p > c]
    K = polyhedron(vrep(vcat(top, sec)), LIB)
    volume(K), W, sec, top
end
