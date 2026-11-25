"""Verifica las densidades locales y la formula singular truncada del doc 44."""

from fractions import Fraction
from math import gcd
import random


def primes_upto(n: int) -> list[int]:
    out: list[int] = []
    for x in range(2, n + 1):
        if all(x % p for p in out if p * p <= x):
            out.append(x)
    return out


PRIMES = primes_upto(97)
ODD_PRIMES = [p for p in PRIMES if p > 2]


def local_roots(N: int, q: int, r: int, p: int) -> int:
    return sum(1 for x in range(p) if (x * (N - q * r * x)) % p == 0)


def local_factor(nu: int, p: int) -> Fraction:
    return Fraction(p * (p - nu), (p - 1) ** 2)


def check_case(N: int, q: int, r: int) -> int:
    assert N % 2 == 0 and q > 2 and r > 2 and gcd(N, q * r) == 1
    actual = Fraction(1)
    generic = Fraction(1)
    correction = Fraction(1)
    checks = 0
    for p in PRIMES:
        nu = local_roots(N, q, r, p)
        expected = 1 if (N * q * r) % p == 0 else 2
        assert nu == expected, (N, q, r, p, nu, expected)
        actual *= local_factor(nu, p)
        if p == 2:
            generic *= Fraction(2)
        else:
            generic *= local_factor(2, p)
            if (N * q * r) % p == 0:
                correction *= Fraction(p - 1, p - 2)
        checks += 1
    assert actual == generic * correction, (N, q, r)
    return checks


def main() -> None:
    rng = random.Random(20261125)
    checks = 0
    for _ in range(400):
        q = rng.choice(ODD_PRIMES)
        r = rng.choice([p for p in ODD_PRIMES if p != q])
        while True:
            N = 2 * rng.randint(2, 5000)
            if gcd(N, q * r) == 1:
                break
        checks += check_case(N, q, r)
    print(
        f"PASS: 400 ternas; {checks} factores locales; "
        "formula singular truncada exacta"
    )


if __name__ == "__main__":
    main()

