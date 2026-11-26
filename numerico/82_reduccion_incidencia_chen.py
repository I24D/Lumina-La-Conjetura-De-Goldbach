"""Pre-registro (antes de ejecutar), bucle Goldbach 2026-11-26 02:07 ET.

Hipotesis exactas que se intentan refutar:
H1. Para m y Q cuadrado-libre con Q|rad(m), si r=m/Q y u es la parte
    Q-lisa maxima de r, entonces v=r/u es coprimo con Q y
    omega(rad(m)/Q)=omega(v).
H2. Para N par y (Q,N)=1, las dos formas r y N-Qr excluyen modulo ell
    una clase si ell|NQ, y dos clases si ell no divide NQ.

Criterio de refutacion fijado antes de correr: cualquier discrepancia en
2<=m<=100000, 1<=k<=7 para H1, o en N par <=500, Q cuadrado-libre <=300
coprimo con N y primos ell<=97 para H2. La prueba es algebraica; el calculo
solo controla errores de formulacion/implementacion.
"""

from math import gcd


def distinct_primes(n: int) -> list[int]:
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
    return all(n % (p * p) for p in distinct_primes(n))


def primes_upto(n: int) -> list[int]:
    return [x for x in range(2, n + 1) if len(distinct_primes(x)) == 1 and distinct_primes(x)[0] == x]


def subsets_of_size(xs: list[int], size: int):
    def rec(i: int, chosen: list[int]):
        if len(chosen) == size:
            q = 1
            for x in chosen:
                q *= x
            yield q
            return
        if len(chosen) + len(xs) - i < size:
            return
        for j in range(i, len(xs)):
            yield from rec(j + 1, chosen + [xs[j]])

    yield from rec(0, [])


h1 = 0
for m in range(2, 100001):
    ps = distinct_primes(m)
    for k in range(1, 8):
        if k - 1 > len(ps):
            continue
        direct = 0
        decomposed = 0
        for q in subsets_of_size(ps, k - 1):
            remainder_distinct = len(ps) - (k - 1)
            direct += remainder_distinct in (1, 2)

            r = m // q
            u = 1
            for p in distinct_primes(q):
                while r % p == 0:
                    u *= p
                    r //= p
            v = r
            assert gcd(v, q) == 1
            decomposed += len(distinct_primes(v)) in (1, 2)
        assert direct == decomposed, (m, k, direct, decomposed)
        h1 += 1

h2 = 0
small_primes = primes_upto(97)
for n in range(4, 501, 2):
    for q in range(1, 301):
        if not squarefree(q) or gcd(n, q) != 1:
            continue
        for ell in small_primes:
            roots = sum(1 for r in range(ell) if (r * (n - q * r)) % ell == 0)
            expected = 1 if (n * q) % ell == 0 else 2
            assert roots == expected, (n, q, ell, roots, expected)
            h2 += 1

print(f"PASS H1={h1} H2={h2}")
