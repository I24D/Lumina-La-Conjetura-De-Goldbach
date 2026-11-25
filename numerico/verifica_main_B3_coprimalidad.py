"""Verifica la identidad finita (28.3) con pesos racionales arbitrarios."""

from fractions import Fraction
from math import gcd
import random


def mobius(n: int) -> int:
    if n == 1:
        return 1
    out = 1
    p = 2
    while p * p <= n:
        if n % p == 0:
            n //= p
            out = -out
            if n % p == 0:
                return 0
            while n % p == 0:
                n //= p
        p += 1
    if n > 1:
        out = -out
    return out


def phi(n: int) -> int:
    return sum(1 for a in range(1, n + 1) if gcd(a, n) == 1)


def divisors(n: int):
    return [r for r in range(1, n + 1) if n % r == 0]


def in_interval(n: int, N: int, C: int, d: int) -> bool:
    return max(1, N - 2 * C * d) <= n < N - C * d


def verify_case(N: int, C: int, K: int, values: dict[int, int]) -> None:
    # L(d) is deliberately non-logarithmic: the algebra holds for any weight.
    L = lambda d: d * d + 1
    lhs = Fraction(0)
    rhs = Fraction(0)

    for d in range(K + 1, 2 * K + 1):
        if gcd(d, N) != 1:
            continue
        wd = sum(
            values[n]
            for n in range(1, N + 1)
            if in_interval(n, N, C, d) and gcd(n, N) == gcd(n, d) == 1
        )
        lhs += Fraction(mobius(d) * L(d) * wd, phi(d))

        for r in divisors(d):
            e = d // r
            if gcd(r, e) != 1 or mobius(r) == 0:
                continue
            wr = sum(
                values[n]
                for n in range(1, N + 1)
                if in_interval(n, N, C, d) and gcd(n, N) == 1 and n % r == 0
            )
            rhs += Fraction(mobius(e) * L(d) * wr, phi(r) * phi(e))

    assert lhs == rhs, (N, C, K, lhs, rhs)


def main() -> None:
    rng = random.Random(20261125)
    for _ in range(500):
        N = 2 * rng.randint(20, 80)
        C = rng.randint(1, 6)
        K = rng.randint(2, 20)
        values = {n: rng.randint(-7, 7) for n in range(1, N + 1)}
        verify_case(N, C, K, values)
    print("PASS: 500 exact checks of (28.3)")


if __name__ == "__main__":
    main()
