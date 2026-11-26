#!/usr/bin/env python3
"""Verificacion racional de las identidades layer-cake del doc 60."""
from fractions import Fraction as F
import random


def pos(x):
    return x if x > 0 else F(0)


def layer_integral(A, u):
    # longitud exacta de {t in [0,A] : u<t}
    return pos(A - u)


def discrepancy(values, residue, d, c):
    """L_d(f)=sum_{n=residue(d)}f(n)-c*sum_n f(n)."""
    return sum(v for n, v in values if n % d == residue % d) - c * sum(
        v for _, v in values
    )


def integrated_threshold_discrepancy(items, A, residue, d, c):
    """Integra L_d(1_{u<t}) por intervalos de quiebre, exactamente."""
    cuts = sorted({F(0), A, *(u for _, u in items if 0 < u < A)})
    total = F(0)
    for left, right in zip(cuts, cuts[1:]):
        mid = (left + right) / 2
        vals = [(n, F(1) if u < mid else F(0)) for n, u in items]
        total += (right - left) * discrepancy(vals, residue, d, c)
    return total


def main():
    rng = random.Random(600044)
    for _ in range(5000):
        A = F(rng.randrange(0, 41), rng.randrange(1, 17))
        u = sum((F(rng.randrange(0, 31), rng.randrange(1, 19))
                 for _ in range(rng.randrange(0, 8))), F(0))
        assert pos(A - u) == layer_integral(A, u)

    for _ in range(500):
        A = F(rng.randrange(1, 31), rng.randrange(1, 13))
        items = []
        for n in range(1, rng.randrange(4, 25)):
            u = sum((F(rng.randrange(0, 15), rng.randrange(1, 13))
                     for _ in range(rng.randrange(1, 6))), F(0))
            items.append((n, u))
        d = rng.randrange(1, 9)
        residue = rng.randrange(0, d)
        c = F(rng.randrange(0, 7), rng.randrange(1, 9))
        weighted = [(n, pos(A - u)) for n, u in items]
        lhs = discrepancy(weighted, residue, d, c)
        rhs = integrated_threshold_discrepancy(items, A, residue, d, c)
        assert lhs == rhs, (A, d, residue, c, lhs, rhs)

    print("PASS: (60.2) 5000 tuplas; (60.3) 500 funcionales exactos; Fraction")


if __name__ == "__main__":
    main()
