#!/usr/bin/env python3
"""Verifica las identidades exactas del documento 48."""

from fractions import Fraction


def formulas(h):
    c = [[Fraction(int(i < j) * (2 + i - j)) for j in range(h)]
         for i in range(h)]
    rows = [sum(row, Fraction()) for row in c]
    cols = [sum(c[i][j] for i in range(h)) for j in range(h)]
    total = sum(rows, Fraction())
    l1 = sum((abs(x) for row in c for x in row), Fraction())
    return c, rows, cols, total, l1


def main():
    for h in range(2, 301):
        c, rows, cols, total, l1 = formulas(h)
        for i in range(h):
            for j in range(h):
                k2 = -(j - i - 1) if i < j else 0
                assert c[i][j] == int(i < j) + k2

                # Enumeracion literal independiente en tamanos pequenos.
                if h <= 40:
                    literal = -sum(1 for k in range(h) if i < k < j)
                    assert k2 == literal

        expected_rows = [
            Fraction((h - 1 - i) * (3 - (h - 1 - i)), 2)
            for i in range(h)
        ]
        expected_cols = [Fraction(j * (3 - j), 2) for j in range(h)]
        assert rows == expected_rows
        assert cols == expected_cols
        assert total == Fraction(h * (h - 1) * (5 - h), 6)
        assert l1 == Fraction(h - 1) * (
            1 + Fraction((h - 3) * (h - 2), 6)
        )
        if h >= 5:
            row_l1 = sum(abs(x) for x in rows)
            expected = Fraction(2) + sum(
                (Fraction(m * (m - 3), 2) for m in range(4, h)),
                Fraction(),
            )
            assert row_l1 == expected

    _, rows, _, _, l1 = formulas(300)
    ratio = float(sum(abs(x) for x in rows) / l1)
    assert 0.97 < ratio < 1.03
    print(f"PASS: identidades h=2..300; cociente marginal/L1 en h=300 = {ratio:.8f}")


if __name__ == "__main__":
    main()
