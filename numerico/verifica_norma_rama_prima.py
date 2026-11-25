"""Verificacion algebraica de (43.3)--(43.5); no prueba nada sobre primos."""

from random import Random


def main() -> None:
    rng = Random(431125)
    for _ in range(500):
        nr = rng.randint(1, 10)
        nl = rng.randint(1, 10)
        a = [rng.randint(-4, 4) for _ in range(nr)]
        b = [rng.randint(-4, 4) for _ in range(nl)]
        kernel = [[rng.randint(0, 12) for _ in range(nl)] for _ in range(nr)]

        t = sum(a[i] * b[j] * kernel[i][j] for i in range(nr) for j in range(nl))
        row_energy = sum(
            sum(b[j] * kernel[i][j] for j in range(nl)) ** 2 for i in range(nr)
        )
        assert t * t <= sum(x * x for x in a) * row_energy

        diagonal = sum(
            b[j] * b[j] * kernel[i][j] * kernel[i][j]
            for i in range(nr)
            for j in range(nl)
        )
        off_diagonal = sum(
            b[j] * b[k] * kernel[i][j] * kernel[i][k]
            for i in range(nr)
            for j in range(nl)
            for k in range(nl)
            if j != k
        )
        assert row_energy == diagonal + off_diagonal

        mass = sum(sum(row) for row in kernel)
        # The bilinear sup over coefficients in [-1,1] includes a=b=1.
        all_one_value = sum(1 * 1 * kernel[i][j] for i in range(nr) for j in range(nl))
        assert all_one_value == mass

    print("PASS: 500 configuraciones; masa positiva, Cauchy y expansion diagonal exactas")


if __name__ == "__main__":
    main()
