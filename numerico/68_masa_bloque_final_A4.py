"""Chequeos algebraicos del doc 68; no verifica el PNT ni Goldbach."""
from fractions import Fraction


def windows(Y: int, eta: Fraction):
    P = Fraction(Y, 1) / (1 + eta)
    D = Fraction(Y, 1) - P
    J = [(P + Fraction(2*i-1, 10)*D, P + Fraction(2*i, 10)*D)
         for i in range(1, 5)]
    return P, D, J


def main():
    checks = 0
    for Y in range(100, 1100):
        eta = Fraction((37*Y) % 997 + 1, 1000)
        if eta >= 1:
            continue
        P, D, J = windows(Y, eta)
        assert D == Fraction(Y)*eta/(1+eta)
        assert D >= Fraction(Y)*eta/2
        assert P < J[0][0]
        assert J[-1][1] < Y
        for i, (a, b) in enumerate(J):
            assert b-a == D/10
            if i:
                assert J[i-1][1] < a
        # Todo punto escogido en ventanas sucesivas queda ordenado.
        mids = [(a+b)/2 for a,b in J]
        assert mids == sorted(mids) and len(set(mids)) == 4
        product = Fraction(1)
        for x in mids:
            product *= x
        assert product < Y**4
        checks += 1
    assert checks == 1000
    print(f"PASS: {checks} configuraciones racionales; geometria (68.2)-(68.3).")


if __name__ == "__main__":
    main()
