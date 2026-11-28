#!/usr/bin/env python3
"""Doc 94c: auditoria finita de la particion exacta de la norma Gram.
Procedencia: LUMINA, GPT-6.2; 2026-11-27.
Metodo: reconstruir K_A(h), T_s(h0), G_smooth y G_switch por rutas separadas.
Estado [NUM]; no prueba cotas asintoticas. Python standard library.
Parametros proxy N=10^12, A=3; se omite el factor log^{-C0} en D.
"""
import math
from collections import defaultdict

N = 10**12
A = 3
PHIA = A - 1
X = N // A
D = math.isqrt(X)
Z = int(X ** 0.1)
Y = int(X ** (1.0 / 3.0))


def primes_upto(n):
    sieve = bytearray(b"\x01") * (n + 1)
    sieve[:2] = b"\x00\x00"
    for p in range(2, math.isqrt(n) + 1):
        if sieve[p]:
            start = p * p
            sieve[start:n + 1:p] = b"\x00" * ((n - start) // p + 1)
    return [p for p in range(2, n + 1) if sieve[p]]


def phi_sieve(n):
    phi = list(range(n + 1))
    for p in range(2, n + 1):
        if phi[p] == p:
            for m in range(p, n + 1, p):
                phi[m] -= phi[m] // p
    return phi


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


primes = primes_upto(Y - 1)
small = [p for p in primes if p < Z and math.gcd(p, A * N) == 1]
phi = phi_sieve(D)

# One row for each admissible pair (s,d), with m=s*d and b=mu(d).
by_s = defaultdict(list)
rows = []
for s in primes:
    if not (Z <= s < Y) or math.gcd(s, A * N) != 1:
        continue
    level = D // s
    for mask in range(1 << len(small)):
        factors = [small[i] for i in range(len(small)) if (mask >> i) & 1]
        d = math.prod(factors)
        if d < level and vplus(sorted(factors, reverse=True), level):
            b = -1 if len(factors) % 2 else 1
            by_s[s].append((d, b, factors))
            rows.append((s, d, b, s * d))

# Build K_A(h) from the original m-sum, and T_s(h0) independently.
K = defaultdict(float)
K_abs = defaultdict(float)
h_factors = {}
T_by_s = {}
Q_by_s = {}
for s, entries in by_s.items():
    T = defaultdict(float)
    for d, b, factors in entries:
        m = s * d
        coeff = b / (PHIA * phi[m])
        all_factors = [s] + factors
        for h in squarefree_divisors(all_factors):
            K[h] += coeff
            K_abs[h] += abs(coeff)
            h_factors[h] = [p for p in all_factors if h % p == 0]
        for h0 in squarefree_divisors(factors):
            T[h0] += b / phi[d]
    T_by_s[s] = T
    Q_by_s[s] = sum(phi[h0] * value * value for h0, value in T.items())

# Route 1: direct Gram norm from K_A(h).
G_direct = sum(phi[h] * value * value for h, value in K.items())
G_abs = sum(phi[h] * value * value for h, value in K_abs.items())

# Route 2a: h with only small prime factors, reconstructed from sum_s T_s/phi(s).
all_h0 = set()
for T in T_by_s.values():
    all_h0.update(T)
G0_formula = 0.0
for h0 in all_h0:
    coeff = sum(T_by_s[s].get(h0, 0.0) / phi[s] for s in T_by_s) / PHIA
    G0_formula += phi[h0] * coeff * coeff

# Route 2b: h=s*h0, with its unique switching prime s.
G1_formula = sum(Q_by_s[s] / phi[s] for s in Q_by_s) / (PHIA * PHIA)

# Split Route 1 directly by the prime support of h.
G0_direct = 0.0
G1_direct = 0.0
bad_support = []
for h, value in K.items():
    factors = h_factors[h]
    switching = [p for p in factors if p >= Z]
    if any(p >= Y for p in factors) or len(switching) > 1:
        bad_support.append((h, factors))
    if switching:
        G1_direct += phi[h] * value * value
    else:
        G0_direct += phi[h] * value * value

sum_inv_phi_s = sum(1.0 / phi[s] for s in by_s)
ratio = G_direct / G_abs if G_abs else float("nan")
errors = (
    abs(G_direct - G0_formula - G1_formula),
    abs(G0_direct - G0_formula),
    abs(G1_direct - G1_formula),
)

print(f"PARAMS N={N} A={A} X={X} z={Z} y={Y} D={D} small={small}")
print(f"SUPPORT pairs={len(rows)} distinct_h={len(K)} bad_support={len(bad_support)}")
print(f"GRAM direct={G_direct:.12f} abs_majorant={G_abs:.12f} signed/abs={ratio:.6f}")
print(f"SMOOTH direct={G0_direct:.12f} formula={G0_formula:.12f}")
print(f"SWITCH direct={G1_direct:.12f} formula={G1_formula:.12f}")
print(f"PARTITION errors total/smooth/switch={errors[0]:.3e}/{errors[1]:.3e}/{errors[2]:.3e}")
print(f"SWITCH blocks={len(Q_by_s)} sum_1/phi(s)={sum_inv_phi_s:.9f} "
      f"Qs_min={min(Q_by_s.values()):.9f} Qs_max={max(Q_by_s.values()):.9f}")

ok = (
    len(rows) == 1444
    and len(K) == 1451
    and not bad_support
    and max(errors) < 1e-10
)
print("PASS" if ok else "FAIL")