#!/usr/bin/env python3
"""Verifica las identidades racionales del documento 47."""

from fractions import Fraction
from random import Random


def center(c):
    m, n = len(c), len(c[0])
    rows = [sum(row, Fraction()) for row in c]
    cols = [sum(c[i][j] for i in range(m)) for j in range(n)]
    total = sum(rows, Fraction())
    z = [
        [c[i][j] - rows[i] / n - cols[j] / m + total / (m * n)
         for j in range(n)]
        for i in range(m)
    ]
    return z, rows, cols, total


def weighted(c, a, b):
    return sum(
        (c[i][j] * (1 + a[i]) * (1 + b[j])
         for i in range(len(c)) for j in range(len(c[0]))),
        Fraction(),
    )


def main():
    rng = Random(471125)
    for _ in range(1000):
        m, n = rng.randint(2, 9), rng.randint(2, 9)
        c = [[Fraction(rng.randint(-9, 9), rng.randint(1, 9))
              for _ in range(n)] for _ in range(m)]
        z, rows, cols, total = center(c)
        assert all(sum(row, Fraction()) == 0 for row in z)
        assert all(sum(z[i][j] for i in range(m)) == 0 for j in range(n))
        for i in range(m):
            for j in range(n):
                additive = rows[i] / n + cols[j] / m - total / (m * n)
                assert z[i][j] + additive == c[i][j]

        a = [Fraction(rng.randint(-3, 3), rng.randint(2, 11)) for _ in range(m)]
        b = [Fraction(rng.randint(-3, 3), rng.randint(2, 11)) for _ in range(n)]
        expanded = (total
                    + sum(a[i] * rows[i] for i in range(m))
                    + sum(b[j] * cols[j] for j in range(n))
                    + sum(c[i][j] * a[i] * b[j]
                          for i in range(m) for j in range(n)))
        assert weighted(c, a, b) == expanded
        assert weighted(z, a, b) == sum(
            z[i][j] * a[i] * b[j]
            for i in range(m) for j in range(n)
        )

    for h in range(2, 61):
        k = [[Fraction(int(i < j)) for j in range(h)] for i in range(h)]
        _, rows, cols, total = center(k)
        assert rows == [Fraction(h - 1 - i) for i in range(h)]
        assert cols == [Fraction(j) for j in range(h)]
        assert total == Fraction(h * (h - 1), 2)
        assert total == sum(abs(x) for row in k for x in row)

    print("PASS: 1000 matrices racionales; nucleo triangular h=2..60")


if __name__ == "__main__":
    main()
