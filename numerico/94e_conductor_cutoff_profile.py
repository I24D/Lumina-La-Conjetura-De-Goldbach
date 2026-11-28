#!/usr/bin/env python3
"""Doc 94e: perfil de masa Gram por corte de conductor completo ah <= R0.
Procedencia: LUMINA, GPT-6.2; 2026-11-27.
Control finito del proxy N=10^12, A=3; no prueba uniformidad asintotica.
"""
import math
import runpy
from pathlib import Path

DOC94 = str(Path(__file__).resolve().with_name("94_conductores_grama.py"))
d = runpy.run_path(DOC94)
K, phi, A, N = d["K"], d["phi"], d["A"], d["N"]

def divisors_with_phi(n):
    out = []
    for a in range(1, n + 1):
        if n % a == 0:
            # A=3 in this proxy, but compute phi(a) directly.
            x, ph = a, a
            p = 2
            while p * p <= x:
                if x % p == 0:
                    ph -= ph // p
                    while x % p == 0:
                        x //= p
                p += 1
            if x > 1:
                ph -= ph // x
            out.append((a, ph))
    return out

divsA = divisors_with_phi(A)
G = sum(phi[h] * v * v for h, v in K.items())
total = A * G
print(f"PARAMS N={N} A={A} G={G:.12f} total_conductor_mass={total:.12f}")
for kappa in (1.0, 2.0, 2.5, 2.7, 2.8, 2.9, 3.0):
    R0 = math.log(N) ** kappa
    low = 0.0
    low_count = 0
    high_count = 0
    for h, value in K.items():
        for a, phia in divsA:
            term = phia * phi[h] * value * value
            if a * h <= R0:
                low += term
                low_count += 1
            else:
                high_count += 1
    high = total - low
    print(
        f"kappa={kappa:.1f} R0={R0:.3f} "
        f"low_mass={low:.12f} low_share={low/total:.6%} "
        f"high_mass={high:.12f} high_share={high/total:.6%} "
        f"low_pairs={low_count} high_pairs={high_count}"
    )