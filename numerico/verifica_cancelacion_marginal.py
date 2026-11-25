#!/usr/bin/env python3
"""Verificacion racional de las identidades del documento 46."""

from fractions import Fraction
from random import Random


def zero_margin_matrix(rng: Random, m: int, n: int):
    c = [[0 for _ in range(n)] for _ in range(m)]
    for i in range(m - 1):
        for j in range(n - 1):
            c[i][j] = rng.randint(-20, 20)
    for i in range(m - 1):
        c[i][n - 1] = -sum(c[i][j] for j in range(n - 1))
    for j in range(n - 1):
        c[m - 1][j] = -sum(c[i][j] for i in range(m - 1))
    c[m - 1][n - 1] = -sum(c[m - 1][j] for j in range(n - 1))
    assert all(sum(row) == 0 for row in c)
    assert all(sum(c[i][j] for i in range(m)) == 0 for j in range(n))
    return c


def check_matrix(c, qs, rs):
    a = [Fraction(1, q - 2) for q in qs]
    b = [Fraction(1, r - 2) for r in rs]
    u = [1 + x for x in a]
    v = [1 + x for x in b]

    weighted = sum(Fraction(c[i][j]) * u[i] * v[j]
                   for i in range(len(qs)) for j in range(len(rs)))
    mixed = sum(Fraction(c[i][j]) * a[i] * b[j]
                for i in range(len(qs)) for j in range(len(rs)))
    assert weighted == mixed

    l1 = sum(abs(x) for row in c for x in row)
    bound = Fraction(l1, (min(qs) - 2) * (min(rs) - 2))
    assert abs(weighted) <= bound


def check_rectangle(q1, q2, r1, r2):
    aq1, aq2 = Fraction(1, q1 - 2), Fraction(1, q2 - 2)
    br1, br2 = Fraction(1, r1 - 2), Fraction(1, r2 - 2)
    lhs = (1 + aq1) * (1 + br1) - (1 + aq1) * (1 + br2)
    lhs -= (1 + aq2) * (1 + br1)
    lhs += (1 + aq2) * (1 + br2)
    rhs = Fraction((q2 - q1) * (r2 - r1),
                   (q1 - 2) * (q2 - 2) * (r1 - 2) * (r2 - 2))
    assert lhs == (aq1 - aq2) * (br1 - br2) == rhs


def main():
    rng = Random(20261125)
    primes = [5, 7, 11, 13, 17, 19, 23, 29, 31, 37, 41, 43, 47, 53]
    for _ in range(1000):
        m = rng.randint(2, 6)
        n = rng.randint(2, 6)
        qs = sorted(rng.sample(primes, m))
        rs = sorted(rng.sample(primes, n))
        c = zero_margin_matrix(rng, m, n)
        check_matrix(c, qs, rs)
        q1, q2 = sorted(rng.sample(primes, 2))
        r1, r2 = sorted(rng.sample(primes, 2))
        check_rectangle(q1, q2, r1, r2)
    print("PASS: 1000 matrices de margenes nulos; identidad, cota de segundo orden y atomos rectangulares")


if __name__ == "__main__":
    main()
