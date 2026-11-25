"""Verifica identidades finitas (53.5)-(53.6); no prueba distribucion."""

from math import isqrt


def primes_upto(n):
    sieve = bytearray(b"\x01") * (n + 1)
    sieve[:2] = b"\x00\x00"
    for p in range(2, isqrt(n) + 1):
        if sieve[p]:
            sieve[p * p : n + 1 : p] = b"\x00" * (((n - p * p) // p) + 1)
    return [i for i in range(2, n + 1) if sieve[i]]


def factors_with_multiplicity(n, primes):
    out = []
    for p in primes:
        if p * p > n:
            break
        while n % p == 0:
            out.append(p)
            n //= p
    if n > 1:
        out.append(n)
    return out


def eligible_weight(t, lo, hi, S, R, primes):
    fs = factors_with_multiplicity(t, primes)
    if not fs or not (lo < fs[0] <= hi):
        return 0
    mfs = fs[1:]
    if not (S <= len(mfs) <= R - 1):
        return 0
    if any(p <= fs[0] for p in mfs) or len(set(mfs)) != len(mfs):
        return 0
    # A positive synthetic exact weight; identities do not use its formula.
    return 1 + len(mfs) + sum(mfs) % 7


def check_case(N, lo, hi, S, R):
    primes = primes_upto(N)
    E = [eligible_weight(t, lo, hi, S, R, primes) for t in range(N + 1)]
    b = [0] * (N + 1)
    for ell in range(1, N):
        b[ell] = E[N - ell]
        assert b[ell] == E[N - ell]
    for d in range(1, min(N, 31)):
        lhs = sum(b[ell] for ell in range(d, N, d))
        rhs = sum(E[t] for t in range(1, N) if (N - t) % d == 0)
        assert lhs == rhs, (N, lo, hi, S, R, d, lhs, rhs)


def main():
    cases = 0
    for N in range(120, 361, 12):
        for lo, hi in ((2, 7), (3, 11), (5, 17)):
            for S, R in ((1, 4), (2, 5), (3, 6), (1, 6)):
                check_case(N, lo, hi, S, R)
                cases += 1
    assert cases == 252
    print(f"PASS: {cases} configuraciones; identidades (53.5)-(53.6) exactas")


if __name__ == "__main__":
    main()
