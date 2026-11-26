"""Chequeos finitos del doc 73; no constituyen una prueba de Goldbach."""

from fractions import Fraction
from math import isqrt, log


def primes_upto(n: int) -> list[int]:
    sieve = bytearray(b"\x01") * (n + 1)
    sieve[:2] = b"\x00\x00"
    for p in range(2, isqrt(n) + 1):
        if sieve[p]:
            sieve[p * p : n + 1 : p] = b"\x00" * (((n - p * p) // p) + 1)
    return [i for i in range(2, n + 1) if sieve[i]]


def check_a2_polytope() -> None:
    # Si siete exponentes son >=1/20 y a lo sumo uno es <1/6,
    # su suma minima es 1/20+6/6>1.
    assert Fraction(1, 20) + 6 * Fraction(1, 6) > 1
    for band_count in (0, 1):
        lower = band_count * Fraction(1, 20) + (7 - band_count) * Fraction(1, 6)
        assert lower > 1


def check_a3_scale() -> None:
    ratios = []
    for k in (20, 40, 80, 160, 320, 400, 500):
        # N=exp(k): N^(19/20)/(N/log^3 N)=k^3 exp(-k/20).
        ratios.append(k**3 * __import__("math").exp(-k / 20))
    # La sucesion no es monotona al principio, pero ya decrece en esta cola.
    assert all(a > b for a, b in zip(ratios[2:], ratios[3:]))
    assert ratios[-1] < 1


def check_nonreduced_classes() -> None:
    for N in (1000, 2310, 10000, 19998):
        ps = primes_upto(N - 1)
        prime_set = set(ps)
        D = max(2, int(N**0.267))
        for d in range(2, D + 1):
            if __import__("math").gcd(d, N) == 1:
                continue
            represented = [q for q in ps if (N - q) % d == 0]
            assert all(N % q == 0 for q in represented)
            assert len(represented) <= len([q for q in prime_set if N % q == 0])


if __name__ == "__main__":
    check_a2_polytope()
    check_a3_scale()
    check_nonreduced_classes()
    print("PASS doc 73: A2 exacta, escala A3 y clases no reducidas")
