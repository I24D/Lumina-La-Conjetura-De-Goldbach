"""Verificacion exacta de (35.1)--(35.4); no es evidencia de Goldbach."""

from collections import defaultdict
from fractions import Fraction
from math import isqrt


def prime_factors_squarefree(n: int):
    factors = []
    p = 2
    while p * p <= n:
        if n % p == 0:
            n //= p
            if n % p == 0:
                return None
            factors.append(p)
        p += 1 if p == 2 else 2
    if n > 1:
        factors.append(n)
    return factors


def dyadic_cell(n: int) -> int:
    return n.bit_length() - 1


checked = 0
energy_by_cell = defaultdict(Fraction)
energy_formula = defaultdict(Fraction)

for d in range(1, 20001):
    factors = prime_factors_squarefree(d)
    if factors is None:
        continue
    omega = len(factors)
    mu = -1 if omega % 2 else 1
    coeff = Fraction(mu, 2**omega)
    counts = defaultdict(int)

    for mask in range(1 << omega):
        r = 1
        for j, p in enumerate(factors):
            if mask & (1 << j):
                r *= p
        e = d // r
        counts[(dyadic_cell(r), dyadic_cell(e))] += 1

    c = {cell: coeff * multiplicity for cell, multiplicity in counts.items()}
    assert sum(c.values(), Fraction()) == mu
    assert sum((abs(v) for v in c.values()), Fraction()) == 1

    # Sustituto determinista no nulo de log(d) Delta_d(F).
    a = Fraction((d * d + 3 * d + 7) % 101 + 1, (d % 17) + 1)
    for cell, value in c.items():
        energy_by_cell[cell] += value * value * a * a
        m = counts[cell]
        energy_formula[cell] += Fraction(m * m, 4**omega) * a * a
    checked += 1

assert energy_by_cell == energy_formula
assert checked > 10000
print(f"PASS: {checked} squarefree d; identidades (35.1)-(35.4) exactas")
