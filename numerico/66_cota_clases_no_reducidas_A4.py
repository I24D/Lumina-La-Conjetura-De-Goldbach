#!/usr/bin/env python3
"""Chequeo exacto de (66.1)-(66.2); no prueba una estimacion asintotica."""

from fractions import Fraction
from math import gcd


def prime_factors(n: int) -> list[int]:
    out = []
    p = 2
    while p * p <= n:
        if n % p == 0:
            out.append(p)
            while n % p == 0:
                n //= p
        p += 1
    if n > 1:
        out.append(n)
    return out


def squarefree(n: int) -> bool:
    p = 2
    while p * p <= n:
        if n % (p * p) == 0:
            return False
        p += 1
    return True


def tau(n: int) -> int:
    ans = 1
    p = 2
    while p * p <= n:
        e = 0
        while n % p == 0:
            n //= p
            e += 1
        if e:
            ans *= e + 1
        p += 1
    return ans * (2 if n > 1 else 1)


def rough(t: int, z: int) -> bool:
    return t >= 2 and min(prime_factors(t)) >= z


def check(N: int, Q: int, z: int) -> None:
    E = {
        t: Fraction((7 * t * t + 3 * N + z) % 17, 16)
        for t in range(2, N)
        if rough(t, z)
    }
    # Suma original por modulos.
    lhs = sum(
        E.get(t, Fraction(0))
        for d in range(1, Q + 1)
        if squarefree(d) and gcd(d, N) > 1
        for t in range(2, N)
        if (N - t) % d == 0
    )
    # Mismo doble conteo, ahora por t.
    swapped = sum(
        w
        * sum(
            1
            for d in range(1, Q + 1)
            if squarefree(d) and gcd(d, N) > 1 and (N - t) % d == 0
        )
        for t, w in E.items()
    )
    assert lhs == swapped

    delta = max(tau(m) for m in range(1, N))
    large_common_primes = [p for p in prime_factors(N) if p >= z]
    rhs = Fraction(delta) * sum((N - 1) // p for p in large_common_primes)
    assert lhs <= rhs


def main() -> None:
    cases = 0
    for N in range(40, 200, 2):
        for z in (2, 3, 5, 7):
            Q = 1 + ((11 * N + 7 * z) % (N - 1))
            check(N, Q, z)
            cases += 1
    assert cases == 320
    print(f"PASS: {cases} configuraciones; (66.1) exacta y (66.2) verificada")


if __name__ == "__main__":
    main()
