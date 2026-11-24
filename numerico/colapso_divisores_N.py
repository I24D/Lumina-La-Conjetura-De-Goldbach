#!/usr/bin/env python3
"""Verifica el colapso exacto de (12.8) a divisores a=1 o a primo."""

from __future__ import annotations

import argparse
import json
import math
import random
from pathlib import Path


def sieve_primes(n: int) -> list[int]:
    mark = bytearray(b"\x01") * (n + 1)
    if n >= 0:
        mark[0] = 0
    if n >= 1:
        mark[1] = 0
    for p in range(2, math.isqrt(n) + 1):
        if mark[p]:
            mark[p * p : n + 1 : p] = b"\x00" * (((n - p * p) // p) + 1)
    return [i for i in range(2, n + 1) if mark[i]]


def von_mangoldt_table(n: int, primes: list[int]) -> list[float]:
    lam = [0.0] * (n + 1)
    for p in primes:
        x = p
        lp = math.log(p)
        while x <= n:
            lam[x] = lp
            if x > n // p:
                break
            x *= p
    return lam


def distinct_prime_factors(n: int, primes: list[int]) -> list[int]:
    out = []
    for p in primes:
        if p * p > n:
            break
        if n % p == 0:
            out.append(p)
            while n % p == 0:
                n //= p
        if n == 1:
            break
    if n > 1:
        out.append(n)
    return out


def squarefree_divisors(factors: list[int]):
    vals = [(1, 1)]
    for p in factors:
        vals += [(a * p, -mu) for a, mu in vals]
    return vals


def g_a(N: int, a: int, primes: list[int], lam: list[float]) -> float:
    total = 0.0
    for p in primes:
        if p >= N or N % p == 0:
            continue
        x = p
        lp = math.log(p)
        while a * x < N:
            total += lp * lam[N - a * x]
            if x > (N - 1) // (a * p):
                break
            x *= p
    return total


def check_one(N: int, primes: list[int], lam: list[float], tol: float = 1e-9):
    fac = distinct_prime_factors(N, primes)
    full_first = 0.0
    full_second = 0.0
    high_support = 0.0
    for a, mu in squarefree_divisors(fac):
        if a < N:
            full_first += mu * math.log(a) * lam[N - a]
        ga = g_a(N, a, primes, lam)
        full_second += mu * ga
        if len(distinct_prime_factors(a, primes)) >= 2:
            high_support = max(high_support, abs(ga), abs(lam[N - a]) if a < N else 0.0)
    lhs = full_first - full_second
    rhs = -g_a(N, 1, primes, lam)
    for ell in fac:
        rhs += g_a(N, ell, primes, lam) - math.log(ell) * lam[N - ell]
    return abs(lhs - rhs) <= tol, abs(lhs - rhs), high_support


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--limit", type=int, default=5000)
    ap.add_argument("--random", type=int, default=300)
    ap.add_argument("--seed", type=int, default=20261124)
    ap.add_argument("--output", type=Path)
    args = ap.parse_args()
    primes = sieve_primes(args.limit)
    lam = von_mangoldt_table(args.limit, primes)
    rng = random.Random(args.seed)
    cases = list(range(4, min(args.limit, 500) + 1, 2))
    cases += [2 * rng.randint(2, args.limit // 2) for _ in range(args.random)]
    failures = []
    max_error = 0.0
    max_high_support = 0.0
    for N in cases:
        ok, err, high = check_one(N, primes, lam)
        max_error = max(max_error, err)
        max_high_support = max(max_high_support, high)
        if not ok:
            failures.append({"N": N, "error": err})
    result = {
        "claim": "Eq. (14.1): only a=1 or a prime divisor of N survives in (12.8)",
        "limit": args.limit,
        "random_cases": args.random,
        "total_cases": len(cases),
        "seed": args.seed,
        "max_abs_error": max_error,
        "max_forbidden_support": max_high_support,
        "failures": failures,
        "status": "PASS" if not failures and max_high_support == 0.0 else "FAIL",
    }
    print(json.dumps(result, indent=2))
    if args.output:
        args.output.parent.mkdir(parents=True, exist_ok=True)
        args.output.write_text(json.dumps(result, indent=2) + "\n", encoding="utf-8")
    raise SystemExit(0 if result["status"] == "PASS" else 1)


if __name__ == "__main__":
    main()
