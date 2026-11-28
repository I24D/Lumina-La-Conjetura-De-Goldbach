#!/usr/bin/env python3
"""Doc 94d: perfil de la diagonalizacion positiva gcd para Q_s.
Procedencia: LUMINA, GPT-6.2; 2026-11-27.
Diagnostico numerico: compara V+ truncado con el divisor completo.
No demuestra cotas uniformes ni Goldbach. Solo Python standard library.
"""
import math
from collections import defaultdict

def primes_upto(n):
    a = bytearray(b"\x01") * (n + 1)
    a[:2] = b"\x00\x00"
    for p in range(2, math.isqrt(n) + 1):
        if a[p]:
            j = p * p
            a[j:n + 1:p] = b"\x00" * ((n - j) // p + 1)
    return [p for p in range(2, n + 1) if a[p]]

def vplus(fs_desc, level):
    product = 1
    for j, p in enumerate(fs_desc, 1):
        product *= p
        if j % 2 == 1 and p * p * product >= level:
            return False
    return True

def phi_sf(fs):
    return math.prod(p - 1 for p in fs)

def divisors_with_factors(fs):
    out = [()]
    for p in fs:
        out += [t + (p,) for t in out]
    return out

def run_one(N, A, s, primes):
    X = N // A
    D = math.isqrt(X)
    z = int(X ** 0.1)
    L = D // s
    small = [p for p in primes if p < z and math.gcd(p, A*N) == 1]
    T = defaultdict(float)
    diag_de = 0.0
    accepted = 0

    for mask in range(1 << len(small)):
        fs = [small[i] for i in range(len(small)) if (mask >> i) & 1]
        d = math.prod(fs)
        if d >= L or not vplus(sorted(fs, reverse=True), L):
            continue
        accepted += 1
        phid = phi_sf(fs)
        mu = -1.0 if len(fs) % 2 else 1.0
        diag_de += 1.0 / phid
        for rs in divisors_with_factors(fs):
            r = math.prod(rs)
            T[r] += mu / phid

    # Exact positive gcd diagonalization: Q = sum_r phi(r) C_s(r)^2.
    q_terms = {}
    for r, value in T.items():
        # T only contains squarefree r made from small primes.
        factors = [p for p in small if r % p == 0]
        q_terms[r] = phi_sf(factors) * value * value
    Q = sum(q_terms.values())

    # The unrestricted squarefree divisor sum has a closed Euler product.
    P = math.prod(1 - 2/(p-1) + p/(p-1)**2 for p in small)
    all_divisors = divisors_with_factors(small)
    Qfull_sum = 0.0
    ratios = []
    for rs in all_divisors:
        r = math.prod(rs)
        phir = phi_sf(rs)
        mu_r = -1.0 if len(rs) % 2 else 1.0
        local = math.prod(1 - 1/(p-1) for p in small if p not in rs)
        cfull = mu_r * local / phir
        Qfull_sum += phir * cfull * cfull
        if r in T and cfull != 0:
            ratios.append(abs(T[r] / cfull))

    cut = L / z
    q_low = sum(v for r, v in q_terms.items() if r <= cut)
    q_high = sum(v for r, v in q_terms.items() if r > cut)
    q_high_sqrt = sum(v for r, v in q_terms.items() if r > math.sqrt(L))
    max_r = max(T, default=1)
    return {
        "z": z, "L": L, "u": math.log(L)/math.log(z),
        "small": len(small), "accepted": accepted, "Q": Q, "P": P,
        "Qlog": Q*math.log(z), "QoverP": Q/P,
        "Qfull_sum": Qfull_sum, "full_err": abs(Qfull_sum-P),
        "diag_de": diag_de, "diag_over_Q": diag_de/Q if Q else float("inf"),
        "q_low": q_low, "q_high": q_high, "q_high_share": q_high/Q if Q else 0,
        "q_high_sqrt_share": q_high_sqrt/Q if Q else 0,
        "max_r": max_r, "max_logr_logz": math.log(max_r)/math.log(z),
        "ratio_C_median": sorted(ratios)[len(ratios)//2] if ratios else float("nan"),
        "ratio_C_max": max(ratios, default=float("nan")),
        "ratio_C_min": min(ratios, default=float("nan")),
    }

def main():
    results = []
    for e in (40, 45, 50, 55, 60):
        N = 2**e
        A = 3
        X = N // A
        D = math.isqrt(X)
        z = int(X ** 0.1)
        y = int(X ** (1/3))
        primes = primes_upto(y)
        candidates = [p for p in primes if z <= p < y and math.gcd(p, A*N) == 1]
        for u0 in (5/3, 2, 3, 4):
            target = D / (z ** u0)
            s = min(candidates, key=lambda p: abs(p-target))
            r = run_one(N, A, s, primes)
            r.update({"e": e, "s": s, "u_target": u0})
            results.append(r)
            print(
                f"e={e} u0={u0:.3f} u={r['u']:.3f} z={r['z']} s={s} "
                f"Qlog={r['Qlog']:.5f} Q/P={r['QoverP']:.5f} "
                f"d=e/Q={r['diag_over_Q']:.3e} "
                f"r>L/z share={r['q_high_share']:.4f} "
                f"r>sqrt(L) share={r['q_high_sqrt_share']:.4f} "
                f"max|C/Cfull|={r['ratio_C_max']:.4f}"
            )
    print("FULL_PRODUCT_IDENTITY_MAX_ERROR",
          max(r["full_err"] for r in results))
    print("QLOG_RANGE",
          min(r["Qlog"] for r in results), max(r["Qlog"] for r in results))
    print("HIGH_R_SHARE_RANGE",
          min(r["q_high_share"] for r in results),
          max(r["q_high_share"] for r in results))
    print("DIAGONAL_d_equals_e_OVER_Q_RANGE",
          min(r["diag_over_Q"] for r in results),
          max(r["diag_over_Q"] for r in results))

if __name__ == "__main__":
    main()