# The one-point criterion for the irreducible root lattices (root normalisation).
# For each lattice: the cell cut out by the roots, its covolume check, the covering radius,
# and the fraction of the cell removed by the half-space <w, v> <= 1 at a deep hole v.
include("common.jl")
using LinearAlgebra

function Aroots(n)            # roots e_i - e_j in R^(n+1)
    R = Vector{Vector{Q}}()
    for i in 1:n+1, j in 1:n+1
        i == j && continue
        v = zeros(Q, n + 1); v[i] = 1; v[j] = -1; push!(R, v)
    end
    R
end
function Droots(n)
    R = Vector{Vector{Q}}()
    for i in 1:n, j in i+1:n, a in (1, -1), b in (1, -1)
        v = zeros(Q, n); v[i] = a; v[j] = b; push!(R, v)
    end
    R
end
Achart(n) = Chart(hcat(Matrix{Q}(I, n, n), -ones(Q, n)))      # rows e_i - e_{n+1}
Dchart(n) = Chart(Matrix{Q}(I, n, n))
E6roots() = filter(r -> r[1] == r[2] == r[3], e8roots())     # ⟂ e1-e2, e2-e3
E6chart() = Chart(vcat(Q[1 1 1 0 0 0 0 0], hcat(zeros(Q, 5, 3), Matrix{Q}(I, 5, 5))))
E7roots() = filter(r -> r[1] + r[2] == 0, e8roots())
E7chart() = Chart(vcat(Q[1 -1 0 0 0 0 0 0], hcat(zeros(Q, 6, 2), Matrix{Q}(I, 6, 6))))
E8chart() = Chart(Matrix{Q}(I, 8, 8))

lattices = [
    ("A3", Aroots(3), Achart(3), 4), ("D4", Droots(4), Dchart(4), 4),
    ("A4", Aroots(4), Achart(4), 5), ("A5", Aroots(5), Achart(5), 6),
    ("D5", Droots(5), Dchart(5), 4), ("D6", Droots(6), Dchart(6), 4),
    ("D7", Droots(7), Dchart(7), 4), ("E6", E6roots(), E6chart(), 3),
    ("E7", E7roots(), E7chart(), 2),
]
for (name, R, C, d) in lattices
    P = cell(C, R)
    V = vertices(P)
    vol = volume(P)
    @assert vol^2 * det(C.B * transpose(C.B)) == d          # covolume = sqrt(det)
    norms = [nrm(C, p) for p in V]
    rho2 = maximum(norms)
    v = V[argmax(norms)]
    if rho2 <= 1
        println(name, ": ", length(V), " vertices, rho0^2 = ", rho2, " <= 1: P ⊆ conv(M)")
        continue
    end
    k, _, _, _ = cap(C, R, V, lift(C, v), 1)
    println(name, ": ", length(V), " vertices, rho0^2 = ", rho2,
            ", fraction removed by the deep-hole centre = ", k / vol)
end
# E8: vertex enumeration only (19440 vertices); covolume 1 is classical.
let C = E8chart(), V = vertices(cell(C, e8roots()))
    norms = [nrm(C, p) for p in V]
    println("E8: ", length(V), " vertices, squared norms ",
            [(q, count(==(q), norms)) for q in sort(unique(norms))],
            "; rho0^2 = ", maximum(norms), " <= 1: P ⊆ conv(M)")
end
