"""Verifica las identidades finitas del doc 34 con aritmetica racional."""

from fractions import Fraction
from math import gcd
from random import Random


def mobius(n: int) -> int:
    ans = 1
    p = 2
    while p * p <= n:
        if n % p == 0:
            n //= p
            ans = -ans
            if n % p == 0:
                return 0
            while n % p == 0:
                n //= p
        p += 1
    return -ans if n > 1 else ans


def omega(n: int) -> int:
    out = 0
    p = 2
    while p * p <= n:
        if n % p == 0:
            out += 1
            while n % p == 0:
                n //= p
        p += 1
    return out + (n > 1)


def phi(n: int) -> int:
    out = n
    p = 2
    while p * p <= n:
        if n % p == 0:
            out -= out // p
            while n % p == 0:
                n //= p
        p += 1
    return out - out // n if n > 1 else out


def rho(n: int) -> Fraction:
    mu = mobius(n)
    return Fraction(mu, 1 << omega(n)) if mu else Fraction(0)


def divisors(n: int):
    return [d for d in range(1, n + 1) if n % d == 0]


def one_case(N: int, C: int, K: int, F: list[int], L) -> None:
    direct = Fraction(0)
    main = Fraction(0)
    error = Fraction(0)
    split = Fraction(0)
    deltas = {}

    for d in range(K + 1, 2 * K + 1):
        mu = mobius(d)
        if not mu or gcd(d, N) != 1:
            continue
        lo = max(1, N - 2 * C * d)
        hi = N - C * d
        A = sum(F[n] for n in range(lo, hi) if n % d == N % d and gcd(n, N) == 1)
        W = sum(F[n] for n in range(lo, hi) if gcd(n, N) == 1 and gcd(n, d) == 1)
        delta = Fraction(A) - Fraction(W, phi(d))
        deltas[d] = delta
        direct += mu * L(d) * A
        main += Fraction(mu * L(d) * W, phi(d))
        error += mu * L(d) * delta

    for d, delta in deltas.items():
        for r in divisors(d):
            e = d // r
            if gcd(r, e) == 1:
                split += rho(r) * rho(e) * L(d) * delta

    assert direct == main + error
    assert error == split


def main() -> None:
    rng = Random(20261125)
    for _ in range(500):
        N = 2 * rng.randint(12, 70)
        C = rng.randint(1, 7)
        K = rng.randint(2, 14)
        F = [rng.randint(-9, 11) for _ in range(N + 1)]
        salt = rng.randint(1, 9)
        one_case(N, C, K, F, lambda d, s=salt: d * d + 3 * d + s)
    print("PASS: 500 configuraciones; B=M+E y split d=re exactos")


if __name__ == "__main__":
    main()
