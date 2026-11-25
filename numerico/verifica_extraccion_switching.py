#!/usr/bin/env python3
from fractions import Fraction


def primes_upto(n):
    sieve = bytearray(b"\x01") * (n + 1)
    sieve[:2] = b"\x00\x00"
    for p in range(2, int(n**0.5) + 1):
        if sieve[p]:
            sieve[p*p:n+1:p] = b"\x00" * (((n-p*p)//p)+1)
    return [i for i in range(2, n + 1) if sieve[i]]


def pminus(n):
    if n == 1:
        return 10**18
    p = 2
    while p*p <= n:
        if n % p == 0:
            return p
        p += 1
    return n

pr = primes_upto(5000)
small_r = [p for p in pr if p <= 97]
checked = 0
for r in small_r:
    for ell in range(2, 5001):
        lhs = pminus(ell) > r
        s = pminus(ell)
        t = ell // s
        rhs = (s in pr) and (s > r) and (pminus(t) >= s)
        assert lhs == rhs, (r, ell, s, t)
        checked += 1

D = 240
boxes = 0
for di in range(1, D//2):
    delta = Fraction(di, D)
    lower = Fraction(1, 2) - delta
    upper = (1 - delta) / 2
    composite_upper = (1 - delta) / 3
    for ri in range(1, D):
        rho = Fraction(ri, D)
        if not (lower < rho <= upper):
            continue
        lam = 1 - delta - rho
        can_be_composite = rho <= lam/2
        assert can_be_composite == (rho <= composite_upper)
        if delta <= Fraction(1, 4):
            assert not can_be_composite
        boxes += 1
    if delta > Fraction(1, 4):
        assert composite_upper - lower == (4*delta - 1)/6
    else:
        assert composite_upper <= lower

print(f"PASS: identidad Buchstab en {checked} pares; geometria en {boxes} cajas racionales")
