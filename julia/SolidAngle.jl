# The solid-angle fraction Omega_n at a vertex of the regular n-simplex, as the orthant
# probability of an equicorrelated Gaussian vector with correlation -1/n:
#   Omega_n = 2 ∫_0^∞ φ(t) Re[(1/2 + (i/2) erfi(a t/√2))^n] dt,  a = (n+1)^(-1/2),
# and the Rogers–Böröczky cell bound vol(T_n)/((n+1) Omega_n) for unit balls.
using QuadGK, Printf
setprecision(BigFloat, 256)

function erfi_big(x::BigFloat)            # 2/√π Σ x^(2k+1) / (k! (2k+1))
    s = zero(BigFloat); term = x; k = 0
    while true
        add = term / (2k + 1)
        s += add
        abs(add) < eps(BigFloat) * abs(s) && k > 5 && break
        k += 1; term *= x^2 / k
    end
    2s / sqrt(BigFloat(pi))
end
function omega(n)
    a = 1 / sqrt(BigFloat(n + 1))
    f(t) = exp(-t^2 / 2) / sqrt(2 * BigFloat(pi)) *
           real((BigFloat(1) / 2 + im * erfi_big(a * t / sqrt(BigFloat(2))) / 2)^n)
    I, err = quadgk(f, BigFloat(0), BigFloat(4), BigFloat(12), BigFloat(30), BigFloat(80);
                    rtol = BigFloat(10)^-45)
    2I
end
volT(n) = BigFloat(2)^n / factorial(big(n)) * sqrt(BigFloat(n + 1) / BigFloat(2)^n)
for n in 2:8
    w = omega(n)
    @printf("%d  Omega = %.22e   cell bound = %.15f\n", n, w, volT(n) / ((n + 1) * w))
end
w7 = omega(7)
@printf("1/(1260 Omega_7) = %.18f\n", 1 / (1260w7))
@assert 1 / (1260w7) > BigFloat("14.3234")
