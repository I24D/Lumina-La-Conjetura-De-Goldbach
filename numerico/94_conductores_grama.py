#!/usr/bin/env python3
"""Doc 94: identidad exacta de Gram/gcd para los coeficientes por conductor.
Procedencia: LUMINA, GPT-6.2; 2026-11-27. Método: enumeración exacta del soporte y
comparación numérica de la identidad Gram/gcd. Estado [NUM], no prueba asintótica.
Fuentes: docs 86-94; Nathanson, Additive Number Theory, caps. 9-10.
Control FINITO de identidades; no prueba una cota asintotica.
Parametros proxy: N=10^12, A=3, X=N/A, z=floor(X^0.1),
y=floor(X^(1/3)), D=floor(sqrt(X)) (sin perdida logaritmica).
"""
import math
from collections import defaultdict

N = 10**12
A = 3
X = N // A
D = math.isqrt(X)
z = int(X ** 0.1)
y = int(X ** (1.0 / 3.0))


def primes_upto(n):
    is_p = bytearray(b"\x01") * (n + 1)
    is_p[:2] = b"\x00\x00"
    for p in range(2, math.isqrt(n) + 1):
        if is_p[p]:
            start = p * p
            is_p[start:n + 1:p] = b"\x00" * ((n - start) // p + 1)
    return [p for p in range(2, n + 1) if is_p[p]]


def phi_sieve(n):
    ph = list(range(n + 1))
    for p in range(2, n + 1):
        if ph[p] == p:
            for j in range(p, n + 1, p):
                ph[j] -= ph[j] // p
    return ph


def vplus(factors_desc, level):
    product = 1
    for j, p in enumerate(factors_desc, start=1):
        product *= p
        if j % 2 == 1 and p * p * product >= level:
            return False
    return True


def squarefree_divisors(factors):
    out = [1]
    for p in factors:
        out += [d * p for d in out]
    return out


primes = primes_upto(y - 1)
small = [p for p in primes if p < z and math.gcd(p, A * N) == 1]
phi = phi_sieve(D)
rows = []
for s in primes:
    if s < z or s >= y or math.gcd(s, A * N) != 1:
        continue
    level = D // s
    for mask in range(1 << len(small)):
        fs = [small[i] for i in range(len(small)) if (mask >> i) & 1]
        d = math.prod(fs)
        if d < level and vplus(sorted(fs, reverse=True), level):
            mu = -1 if len(fs) % 2 else 1
            rows.append((s, d, mu, s * d, fs))

# K_A(h) = 1/phi(A) * sum_{m:h|m} b(m)/phi(m).
K = defaultdict(float)
K_abs = defaultdict(float)
h_factors = {}
for s, d, b, m, fs in rows:
    factors = [s] + fs
    denom = (A - 1) * phi[m]  # phi(A)=2 for A=3
    coeff = b / denom
    for h in squarefree_divisors(factors):
        K[h] += coeff
        K_abs[h] += abs(coeff)
        h_factors[h] = [p for p in factors if h % p == 0]

def K_of_conductor(r):
    a = math.gcd(r, A)
    h = r // a
    # Admissible only if a|A, h squarefree, (h,A*N)=1.
    if A % a != 0 or math.gcd(h, A * N) != 1:
        return 0.0
    x = h
    p = 2
    while p * p <= x:
        if x % p == 0:
            x //= p
            if x % p == 0:
                return 0.0
            while x % p == 0:
                x //= p
        p += 1 if p == 2 else 2
    return K.get(h, 0.0)

norm = sum(phi[h] * v * v for h, v in K.items())
norm_abs = sum(phi[h] * v * v for h, v in K_abs.items())

# Exact Gram identity:
# sum_h phi(h) K(h)^2 = 1/phi(A)^2 *
# sum_{m,n} b(m)b(n) gcd(m,n)/(phi(m)phi(n)).
gram = 0.0
for s, d, b, m, fs in rows:
    for s2, d2, b2, m2, fs2 in rows:
        gram += (b * b2 * math.gcd(m, m2) /
                 ((A - 1) ** 2 * phi[m] * phi[m2]))

# Exact support restriction: h has no prime >= y and at most one prime >= z.
bad_support = []
for h, value in K.items():
    if abs(value) < 1e-15:
        continue
    fs = h_factors[h]
    if any(p >= y for p in fs) or sum(p >= z for p in fs) > 1:
        bad_support.append((h, fs, value))

examples = {r: K_of_conductor(r) for r in (3, 7, 11, 13, 21, 17)}

# Cauchy exponent using #M <= D ~ N^((1-c)/2):
# sqrt(D*N)/[N^(1-c)] = N^((-1+3c)/4).
thresholds = []
for c in (0.039884, 2.0 / 17.0, 0.30, 1.0 / 3.0):
    thresholds.append((c, (-1.0 + 3.0 * c) / 4.0))

print(f"PARAMS N={N} A={A} X={X} z={z} y={y} D={D} small_primes={small}")
print(f"SUPPORT rows={len(rows)} distinct_h={len(K)} bad_support={len(bad_support)}")
print("K_A(r) examples: " + " ".join(f"K({r})={v:.12f}" for r, v in examples.items()))
print(f"GRAM norm={norm:.12f} direct_gram={gram:.12f} abs_majorant={norm_abs:.12f} signed/abs={norm/norm_abs:.6f}")
print("CAUCHY exponents: " + "; ".join(f"c={c:.6f}: exponent={e:.6f}" for c, e in thresholds))
ok = (len(rows) == 1444 and not bad_support and abs(norm - gram) < 1e-10 and norm >= 0
      and abs(K_of_conductor(21) - K.get(7, 0.0)) < 1e-12)
print("PASS" if ok else "FAIL")
