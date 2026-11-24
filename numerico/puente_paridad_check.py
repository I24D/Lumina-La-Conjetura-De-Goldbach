#!/usr/bin/env python3
"""Checks for the Goldbach parity-bridge notes; finite evidence, not a proof.

Agent: Codex. Model: exact model identifier not exposed by this interface.
Date: 2026-11-23. Status: [P] implementation; outputs [N].
Method: exact integer/Fraction identities, independent direct counting, and a
small floating-point DFT check. Sources: mathematical notes 03 and 04;
Sawin-Shusterman, arXiv:1808.04001v2, sections 2-3 and Corollary 6.1.
No third-party dependencies, model calls, or shared-service changes.
"""

from __future__ import annotations

import argparse
import cmath
from collections import Counter
from datetime import datetime, timezone
from fractions import Fraction
import hashlib
import itertools
import json
import math
from pathlib import Path
import platform
import time


def factors(n: int) -> dict[int, int]:
    out: dict[int, int] = {}
    p = 2
    while p * p <= n:
        while n % p == 0:
            out[p] = out.get(p, 0) + 1
            n //= p
        p += 1
    if n > 1:
        out[n] = out.get(n, 0) + 1
    return out


def mobius(f: dict[int, int]) -> int:
    return 0 if any(e > 1 for e in f.values()) else (-1) ** len(f)


def integer_checks(limit: int) -> dict:
    fs = [{}] + [factors(n) for n in range(1, limit + 1)]
    mu = [0] + [mobius(fs[n]) for n in range(1, limit + 1)]
    # A nonzero Mangoldt value is represented by its prime label log(p).
    mangoldt = [None] + [next(iter(fs[n])) if len(fs[n]) == 1 else None
                            for n in range(1, limit + 1)]
    prime = [False] + [len(fs[n]) == 1 and fs[n].get(n) == 1
                       for n in range(1, limit + 1)]

    def add_product(poly: Counter, d: int, x: int, sign: int) -> None:
        q = mangoldt[x]
        if sign and q is not None:
            for p, exponent in fs[d].items():
                poly[tuple(sorted((p, q)))] += sign * exponent

    identities = 0
    prime_power_checks = 0
    for n in range(4, limit + 1, 2):
        direct: Counter = Counter()
        prime_only: Counter = Counter()
        for a in range(1, n):
            if mangoldt[a] is not None and mangoldt[n - a] is not None:
                key = tuple(sorted((mangoldt[a], mangoldt[n - a])))
                direct[key] += 1
                if prime[a] and prime[n - a]:
                    prime_only[key] += 1
        assert all(direct[k] - prime_only[k] >= 0 for k in direct)
        prime_power_checks += 1
        cube_root = 1
        while (cube_root + 1) ** 3 <= n:
            cube_root += 1
        for cutoff in sorted({1, cube_root, math.isqrt(n), n - 1}):
            small: Counter = Counter()
            tail: Counter = Counter()
            for d in range(1, cutoff + 1):
                if mu[d]:
                    for c in range(1, (n - 1) // d + 1):
                        add_product(small, d, n - c * d, -mu[d])
            for c in range(1, (n - 1) // (cutoff + 1) + 1):
                for d in range(cutoff + 1, (n - 1) // c + 1):
                    add_product(tail, d, n - c * d, mu[d])
            residual = Counter(direct)
            residual.subtract(small)
            residual.update(tail)  # R = small - tail.
            assert all(v == 0 for v in residual.values()), (n, cutoff, residual)
            identities += 1

    convolution_checks = 0
    for parameter in (2, 4, 6, 30, 210):
        for n in range(1, limit + 1):
            phi = n
            for p in fs[n]:
                phi = phi // p * (p - 1)

            def h(b: int) -> Fraction:
                value = Fraction(1)
                for p, k in fs[b].items():
                    value *= (Fraction(1, p**k) if parameter % p == 0
                              else Fraction(-1, (p - 1) * p**k))
                return value

            got = sum((Fraction(mu[n // b], n // b) * h(b)
                       for b in range(1, n + 1) if n % b == 0), Fraction())
            expected = Fraction(mu[n], phi) if math.gcd(n, parameter) == 1 else Fraction()
            assert got == expected, (parameter, n, got, expected)
            convolution_checks += 1
    budget_checks = 0
    for n in range(4, limit + 1, 2):
        cutoff = 1
        while (cutoff + 1) ** 5 <= n ** 3:
            cutoff += 1
        for c in range(1, (n - 1) // (cutoff + 1) + 1):
            for d in range(cutoff + 1, (n - 1) // c + 1):
                assert c ** 3 < d ** 2  # sigma=1/10, alpha=2/3; exact integers.
                budget_checks += 1
    size_obstruction = None
    if limit >= 600:
        assert math.gcd(31, 600) == 1 and factors(197) == {197: 1}
        assert 600 - 31 * 13 == 197 and 13 > 8 and 31 > 13
        size_obstruction = {"N": 600, "D": 8, "c": 31, "d": 13,
                            "N_minus_cd_prime": 197,
                            "meaning": "admissible primitive tail term with coefficient c larger than variable d"}
    diagnostics = []
    noncoprime_terms_checked = 0
    for n in sorted({v for v in (30, 64, 100, 210, 256, 512, limit) if v <= limit}):
        cutoff = 1
        while (cutoff + 1) ** 3 <= n:
            cutoff += 1
        primitive_u = 0.0
        signed_tail = 0.0
        bad_absolute = 0.0
        for c in range(1, (n - 1) // (cutoff + 1) + 1):
            running = 0.0
            supremum = 0.0
            for d in range(cutoff + 1, (n - 1) // c + 1):
                q = mangoldt[n - c * d]
                term = mu[d] * math.log(q) if q is not None else 0.0
                running += term
                supremum = max(supremum, abs(running))
                signed_tail += term * math.log(d)
                if math.gcd(c, n) > 1 and term:
                    assert n % q == 0 and c % q == 0
                    noncoprime_terms_checked += 1
                    bad_absolute += abs(term * math.log(d))
            if math.gcd(c, n) == 1:
                primitive_u += supremum
        elementary_bad_bound = math.log(n) * sum(
            math.log(p) * sum(math.prod(e + 1 for e in fs[n - p**k].values())
                              for k in range(1, n.bit_length()) if p**k < n)
            for p in fs[n])
        assert bad_absolute <= elementary_bad_bound + 1e-8
        diagnostics.append({"N": n, "D": cutoff,
                            "primitive_U": primitive_u,
                            "U_times_logN_over_N": primitive_u * math.log(n) / n,
                            "signed_tail": signed_tail,
                            "noncoprime_absolute_tail": bad_absolute,
                            "elementary_noncomp_bound": elementary_bad_bound})
    return {"even_N_through": limit,
            "hyperbola_symbolic_log_polynomial_identities": identities,
            "prime_power_nonnegative_remainder_checks": prime_power_checks,
            "rational_euler_coefficient_checks": convolution_checks,
            "beyond_half_level_integer_budget_checks": budget_checks,
            "coefficient_size_obstruction_example": size_obstruction,
            "arithmetic": "identities: exact integer/Fraction; tail diagnostics: floating logs",
            "noncoprime_tail_terms_checked": noncoprime_terms_checked,
            "tail_diagnostics": diagnostics}


def paley_checks(limit: int = 127) -> dict:
    cases = []
    max_error = 0.0
    for p in range(3, limit + 1, 4):
        if factors(p) != {p: 1}:
            continue
        support = [x for x in range(1, p) if pow(x, (p - 1) // 2, p) == p - 1]
        support_set = set(support)
        count = sum((-x) % p in support_set for x in support)
        assert len(support) == (p - 1) // 2 and count == 0
        spectrum = [sum((cmath.exp(2j * math.pi * k * x / p)
                         for x in support), 0j) for k in range(p)]
        expected_norm2 = (p + 1) / 4
        error = max(abs(abs(z) ** 2 - expected_norm2) for z in spectrum[1:])
        square_sum_error = abs(sum((z * z for z in spectrum), 0j))
        assert error < 1e-7 and square_sum_error < 1e-7
        max_error = max(max_error, error, square_sum_error)
        odd_support = [b for b in range(1, 2 * p, 2) if b % p in support_set]
        odd_set = set(odd_support)
        assert len(odd_support) == len(support)
        assert not any((-b) % (2 * p) in odd_set for b in odd_support)
        odd_spectrum = [sum((cmath.exp(2j * math.pi * k * b / (2 * p))
                             for b in odd_support), 0j) for k in range(2 * p)]
        parity_error = max(abs(odd_spectrum[k + p] + odd_spectrum[k]) for k in range(p))
        odd_norm_error = max(abs(abs(odd_spectrum[k]) ** 2 - expected_norm2)
                             for k in range(2 * p) if k not in (0, p))
        assert parity_error < 1e-7 and odd_norm_error < 1e-7
        max_error = max(max_error, parity_error, odd_norm_error)
        cases.append({"prime": p, "support_size": len(support), "self_sum_at_zero": count,
                      "odd_lift_modulus": 2 * p, "odd_lift_sum_at_even_target": 0})
    return {"cases": cases, "max_floating_error": max_error,
            "scope": "dense finite-field surrogate, not a counterexample about primes in Z"}


def trim(a: list[int], p: int) -> list[int]:
    a = [x % p for x in a]
    while len(a) > 1 and a[-1] == 0:
        a.pop()
    return a


def divmod_poly(a: list[int], b: list[int], p: int) -> tuple[list[int], list[int]]:
    a, b = trim(a, p), trim(b, p)
    assert b != [0]
    quotient = [0] * max(1, len(a) - len(b) + 1)
    inverse = pow(b[-1], -1, p)
    while a != [0] and len(a) >= len(b):
        shift = len(a) - len(b)
        coefficient = a[-1] * inverse % p
        quotient[shift] = coefficient
        for i, value in enumerate(b):
            a[i + shift] = (a[i + shift] - coefficient * value) % p
        a = trim(a, p)
    return trim(quotient, p), a


def derivative(a: list[int], p: int) -> list[int]:
    return trim([i * a[i] for i in range(1, len(a))] or [0], p)


def monic_polynomials(degree: int, p: int):
    for coefficients in itertools.product(range(p), repeat=degree):
        yield list(coefficients) + [1]


def polynomial_mu(a: list[int], p: int) -> int:
    remaining = a[:]
    count = 0
    for degree in range(1, (len(a) - 1) // 2 + 1):
        for divisor in monic_polynomials(degree, p):
            exponent = 0
            while len(remaining) >= len(divisor):
                q, r = divmod_poly(remaining, divisor, p)
                if r != [0]:
                    break
                exponent += 1
                remaining = q
            if exponent >= 2:
                return 0
            count += exponent
    if len(remaining) > 1:
        count += 1
    return (-1) ** count


def determinant_mod(matrix: list[list[int]], p: int) -> int:
    a = [row[:] for row in matrix]
    value = 1
    for col in range(len(a)):
        pivot = next((r for r in range(col, len(a)) if a[r][col] % p), None)
        if pivot is None:
            return 0
        if pivot != col:
            a[pivot], a[col] = a[col], a[pivot]
            value = -value
        entry = a[col][col] % p
        value = value * entry % p
        inverse = pow(entry, -1, p)
        for r in range(col + 1, len(a)):
            multiple = a[r][col] * inverse % p
            for j in range(col, len(a)):
                a[r][j] = (a[r][j] - multiple * a[col][j]) % p
    return value % p


def resultant(a: list[int], b: list[int], p: int) -> int:
    b = trim(b, p)
    if b == [0]:
        return 0
    m, n = len(a) - 1, len(b) - 1
    rows = []
    for i in range(n):
        rows.append([0] * i + a[::-1] + [0] * (n - i - 1))
    for i in range(m):
        rows.append([0] * i + b[::-1] + [0] * (m - i - 1))
    return determinant_mod(rows, p)


def polynomial_checks() -> dict:
    p = 3
    checked = 0
    for degree in range(1, 7):
        for f in monic_polynomials(degree, p):
            disc = (-1) ** (degree * (degree - 1) // 2) * resultant(f, derivative(f, p), p) % p
            chi = 0 if disc == 0 else (1 if pow(disc, (p - 1) // 2, p) == 1 else -1)
            assert polynomial_mu(f, p) == (-1) ** degree * chi, (degree, f)
            checked += 1
    fiber_sums = []
    fiber_cases = 0
    for a1, a2 in itertools.product(range(p), repeat=2):
        r = [0, a1, a2, 0, 1]
        total = 0
        for b0, b1 in itertools.product(range(p), repeat=2):
            f = [b0, a1, a2, b1, 1]  # r + (b0+b1*T)^3 over F_3.
            assert derivative(f, p) == derivative(r, p)
            total += polynomial_mu(f, p)
            fiber_cases += 1
        fiber_sums.append(total)
    return {"field": "F_3", "max_degree": 6,
            "pellet_via_factorization_vs_sylvester_determinant": checked,
            "fixed_derivative_fiber_cases": fiber_cases,
            "fiber_mobius_sums": fiber_sums,
            "scope": "algebraic identities only; F_3 does not meet the large-q twin-prime theorem hypothesis"}


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--hasta", type=int, default=600)
    args = parser.parse_args()
    if not 4 <= args.hasta <= 2000:
        parser.error("Keep this lightweight: 4 <= --hasta <= 2000")
    started = time.perf_counter()
    result = {
        "provenance": {"agent": "Codex", "model": "exact identifier not exposed",
                       "project_date": "2026-11-23", "run_utc": datetime.now(timezone.utc).isoformat(),
                       "status": "[N] finite checks, not Goldbach proof or Lean verification",
                       "script_sha256": hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
                       "python": platform.python_version(),
                       "sources": ["../matematicas/03_sawin_shusterman_alpoge_puente.md",
                                   "../matematicas/04_reduccion_mobius_goldbach.md"]},
        "integer_checks": integer_checks(args.hasta),
        "dense_paley_counterexample": paley_checks(),
        "polynomial_checks": polynomial_checks(),
    }
    result["elapsed_seconds"] = round(time.perf_counter() - started, 4)
    result["pass"] = True
    target = Path(__file__).parent / "resultados" / (
        datetime.now(timezone.utc).strftime("%Y%m%d_%H%M%S_%f") + "_puente_paridad.json")
    with target.open("x", encoding="utf-8") as output:
        output.write(json.dumps(result, indent=2, ensure_ascii=False) + "\n")
    print(json.dumps({"pass": True, "result_file": str(target),
                      "elapsed_seconds": result["elapsed_seconds"],
                      "integer_checks": result["integer_checks"],
                      "paley_cases": len(result["dense_paley_counterexample"]["cases"]),
                      "polynomial_checks": result["polynomial_checks"]}, ensure_ascii=True))


if __name__ == "__main__":
    main()
