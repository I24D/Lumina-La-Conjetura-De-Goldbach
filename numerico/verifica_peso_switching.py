#!/usr/bin/env python3
"""Verifica la obstruccion de rango del peso K(r,l)=1_{P^-(l)>r}."""

from fractions import Fraction


def primes(count: int):
    out = []
    n = 3
    while len(out) < count:
        if all(n % p for p in range(2, int(n**0.5) + 1)):
            out.append(n)
        n += 2
    return out


def rank_q(matrix):
    a = [[Fraction(x) for x in row] for row in matrix]
    rows = len(a)
    cols = len(a[0]) if rows else 0
    rank = 0
    for col in range(cols):
        pivot = next((i for i in range(rank, rows) if a[i][col]), None)
        if pivot is None:
            continue
        a[rank], a[pivot] = a[pivot], a[rank]
        scale = a[rank][col]
        a[rank] = [x / scale for x in a[rank]]
        for i in range(rows):
            if i != rank and a[i][col]:
                factor = a[i][col]
                a[i] = [x - factor * y for x, y in zip(a[i], a[rank])]
        rank += 1
    return rank


def main():
    for h in range(1, 26):
        ps = primes(2 * h)
        r = ps[0::2]
        s = ps[1::2]
        k = [[int(r[i] < s[j]) for j in range(h)] for i in range(h)]
        expected = [[int(i <= j) for j in range(h)] for i in range(h)]
        assert k == expected
        assert rank_q(k) == h
    print("PASS: rangos completos para h=1..25; nucleo no separable")


if __name__ == "__main__":
    main()
