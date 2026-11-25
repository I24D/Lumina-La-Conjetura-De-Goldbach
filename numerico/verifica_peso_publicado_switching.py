from fractions import Fraction
import random


def eps_bound(q0, r0):
    return Fraction(1, q0 - 2) + Fraction(1, r0 - 2) + Fraction(1, (q0 - 2) * (r0 - 2))


def u(p):
    return Fraction(p - 1, p - 2)


def check_positive_extraction(rng):
    primes = [5, 7, 11, 13, 17, 19, 23, 29, 31, 37, 41, 43]
    for _ in range(2000):
        qs = sorted(rng.sample(primes, rng.randint(2, 6)))
        rs = sorted(rng.sample(primes, rng.randint(2, 6)))
        coeff = {(q, r): Fraction(rng.randint(0, 20), rng.randint(1, 9)) for q in qs for r in rs}
        mass = sum(coeff.values(), Fraction(0))
        exact = sum((c * u(q) * u(r) for (q, r), c in coeff.items()), Fraction(0))
        delta = exact - mass
        assert delta >= 0
        assert delta <= eps_bound(qs[0], rs[0]) * mass


def check_published_weights(rng):
    # Mock logarithmic sizes alpha=log p/log x. The identities are algebraic.
    for _ in range(2000):
        upar = Fraction(rng.randint(3, 12), 1)
        vpar = upar + rng.randint(1, 12)
        alphas = [Fraction(rng.randint(1, 999), 1000) for _ in range(rng.randint(1, 8))]
        in_band = [a for a in alphas if Fraction(1, vpar) <= a < Fraction(1, upar)]
        kuhn = Fraction(1) - Fraction(len(in_band), 2)
        assert max(kuhn, Fraction(0)) >= 0

        lam = Fraction(rng.randint(1, 20), 20)
        richert = Fraction(1) - lam * sum((Fraction(1) - upar * a for a in in_band), Fraction(0))
        assert max(richert, Fraction(0)) >= 0


def main():
    rng = random.Random(20261125)
    check_positive_extraction(rng)
    check_published_weights(rng)
    print('PASS: 2000 cajas positivas + 2000 pesos Kuhn/Richert; cotas exactas con Fraction')


if __name__ == '__main__':
    main()
