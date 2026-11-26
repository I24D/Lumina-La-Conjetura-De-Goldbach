#!/usr/bin/env python3
"""Verifica la atomizacion finita del peso Kuhn del doc 62."""
from fractions import Fraction as F
import random

LO, HI = F(1, 20), F(1, 6)

def inside(x):
    return LO <= x < HI

def direct(alpha, factors):
    j = int(inside(alpha))
    k = sum(inside(x) for x in factors)
    return max(F(0), F(1) - F(j + k, 2))

def atom(alpha, factors):
    j = int(inside(alpha))
    k = sum(inside(x) for x in factors)
    if j:
        return F(1, 2) if k == 0 else F(0)
    if k == 0:
        return F(1)
    if k == 1:
        return F(1, 2)
    return F(0)

def main():
    for j in (0, 1):
        for k in range(6):
            lhs = max(F(0), F(1) - F(j + k, 2))
            rhs = (F(1, 2) if k == 0 else F(0)) if j else (
                F(1) if k == 0 else F(1, 2) if k == 1 else F(0)
            )
            assert lhs == rhs, (j, k, lhs, rhs)

    rng = random.Random(620267)
    checked = 0
    while checked < 20000:
        alpha = F(rng.randrange(6, 31), 120)
        total_factors = rng.choice((4, 5))
        fs = sorted(F(rng.randrange(max(7, int(alpha * 120) + 1), 61), 120)
                    for _ in range(total_factors - 1))
        if alpha + sum(fs) > 1:
            continue
        assert all(x > alpha for x in fs)
        assert direct(alpha, fs) == atom(alpha, fs)
        checked += 1
    print('PASS: 12 estados exhaustivos; 20000 configuraciones racionales de 4/5 factores')

if __name__ == '__main__':
    main()
