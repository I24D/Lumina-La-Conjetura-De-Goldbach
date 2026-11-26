#!/usr/bin/env python3
"""Control finito del selector polilogaritmico del doc 83.

PRE-REGISTRO (escrito antes de ejecutar):
H1: para cada representacion seleccionada, m=N-p=Qv tiene exactamente
    k o k+1 primos distintos y Q aparece con exponente uno.
H2: para cada N par en [1000, 50000] y k=2,3,4 existe al menos una
    representacion de la subincidencia u=1.
Refutacion: una sola violacion de H1 o un par (N,k) sin representacion.
La comprobacion de H2 es evidencia numerica, no una demostracion.
"""

from math import gcd, log

LIMIT = 50_000


def sieve(n: int):
    prime = bytearray(b"\x01") * (n + 1)
    prime[0:2] = b"\x00\x00"
    for p in range(2, int(n**0.5) + 1):
        if prime[p]:
            prime[p * p : n + 1 : p] = b"\x00" * (((n - p * p) // p) + 1)
    primes = [p for p in range(2, n + 1) if prime[p]]
    omega = bytearray(n + 1)
    for p in primes:
        for m in range(p, n + 1, p):
            omega[m] += 1
    return prime, primes, omega


def selector(N: int, r: int, odd_primes):
    qs = []
    for q in odd_primes:
        if N % q:
            qs.append(q)
            if len(qs) == r:
                break
    Q = 1
    for q in qs:
        Q *= q
    return qs, Q


def main():
    is_prime, primes, omega = sieve(LIMIT)
    odd_primes = primes[1:]
    h1_failures = []
    h2_count = {2: 0, 3: 0, 4: 0}
    h2_first = {2: None, 3: None, 4: None}
    h2_last = {2: None, 3: None, 4: None}
    checked = 0
    max_ratio = (0.0, None)

    for N in range(1000, LIMIT + 1, 2):
        for k in (2, 3, 4):
            qs, Q = selector(N, k - 1, odd_primes)
            ratio = qs[-1] / log(N)
            if ratio > max_ratio[0]:
                max_ratio = (ratio, (N, k, tuple(qs), Q))
            found = False
            for p in primes:
                if p >= N:
                    break
                if gcd(p, N) != 1:
                    continue
                m = N - p
                if m % Q:
                    continue
                v = m // Q
                if gcd(v, Q) != 1 or omega[v] not in (1, 2):
                    continue
                # H1: Q es squarefree, disjunto de v; por tanto omega(m)=k o k+1.
                checked += 1
                if omega[m] != (k - 1) + omega[v]:
                    h1_failures.append((N, k, p, Q, v, omega[m], omega[v]))
                    break
                found = True
                break
            if not found:
                h2_count[k] += 1
                item = (N, tuple(qs), Q)
                h2_first[k] = h2_first[k] or item
                h2_last[k] = item

    if h1_failures:
        print("FAIL H1", h1_failures[:20])
        raise SystemExit(1)
    print(
        "PASS H1; H2 REFUTED" if sum(h2_count.values()) else "PASS H1+H2",
        f"pares={(LIMIT - 1000)//2 + 1}",
        f"casos={3*((LIMIT - 1000)//2 + 1)}",
        f"representaciones_controladas={checked}",
        f"max_qr/logN={max_ratio[0]:.6f}",
        f"argmax={max_ratio[1]}",
        f"h2_count={h2_count}",
        f"h2_first={h2_first}",
        f"h2_last={h2_last}",
    )


if __name__ == "__main__":
    main()
