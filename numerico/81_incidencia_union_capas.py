"""Checks exact incidence multiplicities behind C_k+C_{k+1}.

Pre-registered claim: choosing an unordered (k-1)-subset of the distinct
prime divisors of m and requiring the remaining radical to have 1 or 2
prime factors counts layer k exactly k times and layer k+1 exactly
binom(k+1, 2) times, and no other layer.  Refutation criterion: one mismatch
for any 2 <= m <= 20000 and 1 <= k <= 8.
"""

from itertools import combinations
from math import comb


def prime_factors(n: int) -> list[int]:
    out = []
    p = 2
    while p * p <= n:
        if n % p == 0:
            out.append(p)
            while n % p == 0:
                n //= p
        p += 1 if p == 2 else 2
    if n > 1:
        out.append(n)
    return out


def incidence_count(m: int, k: int) -> int:
    ps = prime_factors(m)
    if k < 1 or k - 1 > len(ps):
        return 0
    total = 0
    for chosen in combinations(ps, k - 1):
        chosen_set = set(chosen)
        remaining = sum(p not in chosen_set for p in ps)
        total += remaining in (1, 2)
    return total


def predicted_count(m: int, k: int) -> int:
    s = len(prime_factors(m))
    if s == k:
        return k
    if s == k + 1:
        return comb(k + 1, 2)
    return 0


def ordered_count(m: int, k: int) -> int:
    """The same selection with ordered distinct tuples."""
    return incidence_count(m, k) * factorial(k - 1)


def factorial(n: int) -> int:
    out = 1
    for j in range(2, n + 1):
        out *= j
    return out


def main() -> None:
    checks = 0
    for m in range(2, 20001):
        for k in range(1, 9):
            got = incidence_count(m, k)
            want = predicted_count(m, k)
            assert got == want, (m, k, got, want)
            assert ordered_count(m, k) == factorial(k - 1) * want
            checks += 2

    # Weighted identity: arbitrary integer weights show that no analytic
    # property of Lambda or log(rad m) is being smuggled into the count.
    for k in range(1, 8):
        lhs = rhs = 0
        for m in range(2, 5001):
            weight = (37 * m + 11) % 101
            lhs += weight * incidence_count(m, k)
            s = len(prime_factors(m))
            rhs += weight * (
                (k if s == k else 0)
                + (comb(k + 1, 2) if s == k + 1 else 0)
            )
        assert lhs == rhs, (k, lhs, rhs)
        checks += 1

    print(f"PASS {checks}/{checks}: exact incidence identity")


if __name__ == "__main__":
    main()
