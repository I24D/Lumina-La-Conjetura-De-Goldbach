#!/usr/bin/env python3
"""Verifica la geometria racional del doc 40; no prueba Goldbach."""

from fractions import Fraction


def in_critical_wedge(delta: Fraction, rho: Fraction) -> bool:
    lam = 1 - delta - rho
    return Fraction(1, 2) - delta < rho <= (1 - delta) / 2 and rho <= lam


def main() -> None:
    den = 240
    checked = 0
    rescued_by_swap = 0

    for di in range(1, den // 2):
        delta = Fraction(di, den)
        lower = Fraction(1, 2) - delta
        upper = (1 - delta) / 2
        assert upper - lower == delta / 2

        for ri in range(1, den):
            rho = Fraction(ri, den)
            lam = 1 - delta - rho
            critical = in_critical_wedge(delta, rho)
            formula = lower < rho <= upper
            assert critical == formula
            if not critical:
                continue

            checked += 1
            assert rho <= lam
            assert delta + rho > Fraction(1, 2)
            assert delta + lam > Fraction(1, 2)
            assert rho < Fraction(1, 2)
            assert lam < Fraction(1, 2)
            assert rho + lam == 1 - delta > Fraction(1, 2)
            assert (1 - delta) / 2 <= lam < Fraction(1, 2)

            if delta + min(rho, lam) <= Fraction(1, 2):
                rescued_by_swap += 1

    examples = {
        Fraction(1, 12): (Fraction(5, 12), Fraction(11, 24), Fraction(1, 24)),
        Fraction(1, 8): (Fraction(3, 8), Fraction(7, 16), Fraction(1, 16)),
        Fraction(1, 6): (Fraction(1, 3), Fraction(5, 12), Fraction(1, 12)),
        Fraction(1, 4): (Fraction(1, 4), Fraction(3, 8), Fraction(1, 8)),
    }
    for delta, (lower, upper, width) in examples.items():
        assert Fraction(1, 2) - delta == lower
        assert (1 - delta) / 2 == upper
        assert upper - lower == width == delta / 2

    assert checked == 3540
    assert rescued_by_swap == 0
    print(
        f"PASS: {checked} cajas racionales; {len(examples)} ejemplos; "
        f"{rescued_by_swap} rescates por intercambio"
    )


if __name__ == "__main__":
    main()
