#!/usr/bin/env python3
"""Chequeos racionales para el documento 45; no es evidencia de Goldbach."""

from fractions import Fraction
import random


PRIMES = [5, 7, 11, 13, 17, 19, 23, 29, 31, 37, 41, 43, 47]


def u(p: int) -> Fraction:
    return Fraction(p - 1, p - 2)


def check_once(rng: random.Random) -> None:
    qs = sorted(rng.sample(PRIMES, rng.randint(2, 6)))
    rs = sorted(rng.sample(PRIMES, rng.randint(2, 6)))
    coeff = {(q, r): rng.randint(-9, 9) for q in qs for r in rs}

    weighted = sum((Fraction(c) * u(q) * u(r) for (q, r), c in coeff.items()), Fraction())
    base = sum((Fraction(c) for c in coeff.values()), Fraction())
    correction = sum(
        (
            Fraction(c)
            * (
                Fraction(1, q - 2)
                + Fraction(1, r - 2)
                + Fraction(1, (q - 2) * (r - 2))
            )
            for (q, r), c in coeff.items()
        ),
        Fraction(),
    )
    assert weighted == base + correction

    eps = (
        Fraction(1, qs[0] - 2)
        + Fraction(1, rs[0] - 2)
        + Fraction(1, (qs[0] - 2) * (rs[0] - 2))
    )
    total_variation = sum((abs(Fraction(c)) for c in coeff.values()), Fraction())
    assert abs(weighted - base) <= eps * total_variation


def main() -> None:
    rng = random.Random(451144)
    for _ in range(1000):
        check_once(rng)

    q1, q2, r = 11, 17, 19
    base = Fraction(1) - Fraction(1)
    weighted = u(r) * (u(q1) - u(q2))
    expected = u(r) * (Fraction(1, q1 - 2) - Fraction(1, q2 - 2))
    assert base == 0
    assert weighted == expected and weighted > 0
    print("PASS: 1000 matrices racionales; identidad, cota L1 y contraejemplo firmados")


if __name__ == "__main__":
    main()
