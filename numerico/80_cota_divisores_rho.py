#!/usr/bin/env python3
"""Checks for doc 80. Exact finite sanity checks; not a proof of Goldbach."""

from math import comb, isqrt


def tau_k_prime_power(a: int, k: int) -> int:
    return comb(a + k - 1, k - 1)


def divisors(n: int):
    out = []
    for d in range(1, isqrt(n) + 1):
        if n % d == 0:
            out.append(d)
            if d * d != n:
                out.append(n // d)
    return out


def tau2(n: int) -> int:
    return len(divisors(n))


def tau4_from_factorization(n: int) -> int:
    value = 1
    p = 2
    m = n
    while p * p <= m:
        if m % p == 0:
            a = 0
            while m % p == 0:
                m //= p
                a += 1
            value *= tau_k_prime_power(a, 4)
        p += 1
    if m > 1:
        value *= 4
    return value


checks = 0
for a in range(501):
    d4 = tau_k_prime_power(a, 4)
    d6 = tau_k_prime_power(a, 6)
    assert d4 <= d6
    checks += 1

for n in range(1, 5001):
    convolution = sum(tau2(d) * tau2(n // d) for d in divisors(n))
    assert convolution == tau4_from_factorization(n)
    checks += 1

print(f"PASS {checks}/{checks}")
print("local exponent for sum tau_6(n)^2: 6^2-1 = 35")
