# Exact checks for the Voronoi cell P of sqrt(2) E7 and the 56 deep-hole caps.
# L = {x ∈ E8 : x1 + x2 = 0} is E7 with roots of squared length 2 (the L-picture);
# Lambda = sqrt(2) L.  7-volumes in the Lambda-picture are 2^(7/2) * sqrt(2) = 16 times
# the volumes in the chart (x1, x3, ..., x8), whose Gram determinant is 2.
include("common.jl")
using LinearAlgebra

const C7 = Chart(Q[1 -1 0 0 0 0 0 0;
                   0 0 1 0 0 0 0 0;
                   0 0 0 1 0 0 0 0;
                   0 0 0 0 1 0 0 0;
                   0 0 0 0 0 1 0 0;
                   0 0 0 0 0 0 1 0;
                   0 0 0 0 0 0 0 1] .// 1)
# with this basis, w = B'p gives w = (p1, -p1, p2..p7) and the chart volume is the
# coordinate volume in (x1, x3..x8) (the basis vector e1 - e2 is not normalised).
const RL = filter(r -> r[1] + r[2] == 0, e8roots())
@assert length(RL) == 126

P = cell(C7, RL)
V = vertices(P)
norms = [nrm(C7, p) for p in V]
println("P: ", length(V), " vertices; squared norms (L-picture) ",
        [(q, count(==(q), norms)) for q in sort(unique(norms))])
volP = volume(P)
@assert volP^2 * det(C7.B * transpose(C7.B)) == 2      # covolume of E7 is sqrt 2
println("vol P = 16 * ", volP, " = ", 16volP, "   (Lambda-picture)")
@assert 16volP == 16

deep = [p for (p, q) in zip(V, norms) if q == 3//2]
shallow = [p for (p, q) in zip(V, norms) if q == 7//8]
@assert length(deep) == 56 && length(shallow) == 576
for v in deep
    @assert count(r -> grad(C7, r) ⋅ v == 1, RL) == 27
end
println("each of the 56 deep holes lies on exactly 27 facets")

# One cut <w, v_L> <= 1, i.e. the bisector of the extra centre (4/3) v_L.
vL = lift(C7, deep[1])
k1, W1, sec1, top1 = cap(C7, RL, V, vL, 1)
println("one cut: ", length(W1), " vertices, cross-section ", length(sec1),
        " vertices, vertices beyond the cut ", length(top1))
println("  cap volume (hull of apex and cross-section) = 16 * ", k1, " = ", 16k1)
@assert 16k1 == 1//105
v1 = volume(cell(C7, RL, [(vL, Q(1))]))
println("  volume of P ∩ {<w,v_L> <= 1} = 16 * ", v1, " = ", 16v1)
@assert 16v1 == 16 - 1//105

# The seventh-power law at other thresholds c (cap = (1/105)(3 - 2c)^7).
for c in (11//10, 5//4, 7//5)
    kc, _, _, _ = cap(C7, RL, V, vL, c)
    @assert 16kc == (1//105) * (3 - 2c)^7
    println("  threshold ", c, ": cap ", 16kc, " = (1/105)(3-2c)^7")
end

# The 56 caps have disjoint interiors: every vertex of cap j satisfies <w, v_i> <= 1, i ≠ j.
for (j, vj) in enumerate(deep)
    kj, _, secj, topj = cap(C7, RL, V, lift(C7, vj), 1)
    @assert 16kj == 1//105 && topj == [vj]
    for (i, vi) in enumerate(deep)
        i == j && continue
        gi = grad(C7, lift(C7, vi))
        @assert all(gi ⋅ w <= 1 for w in vcat(topj, secj))
    end
end
println("the 56 caps are cones of volume 1/105 with disjoint interiors; vol P56 = ",
        16 - 56//105)

# P56 = P ∩ conv(M).  conv(M) = {x : <x, v> <= 1 for all vertices v of P} (L-picture).
P56 = cell(C7, RL, [(lift(C7, v), Q(1)) for v in deep])
V56 = vertices(P56)
G = [grad(C7, lift(C7, v)) for v in V]
println("P56: ", length(V56), " vertices")
@assert length(V56) == 5364
@assert all(all(g ⋅ x <= 1 for g in G) for x in V56)
println("every vertex of P56 satisfies all 632 facet inequalities of conv(M): P56 ⊆ conv(M)")
println("vertices of P outside conv(M): ", count(x -> any(g ⋅ x > 1 for g in G), V),
        " (the deep holes)")
shallow = [p for (p, q) in zip(V, norms) if q == 7//8]
Gs = [grad(C7, lift(C7, s)) for s in shallow]
mx = maximum(g ⋅ x for g in Gs for x in V56)
@assert mx == 7//8
@assert all(count(x -> g ⋅ x == mx, V56) == 1 for g in Gs)
println("max of <z, s> over P56 and the 576 shallow holes s: ", mx,
        " = |s|^2 < 1, attained only at z = s")
