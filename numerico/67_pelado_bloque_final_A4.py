"""Verifica la geometria exacta del ultimo bloque diadico de doc 67."""
from fractions import Fraction
from math import floor


def last_block(N: int):
    # Evitamos floats: P_j=2^j*N^(1/20), Y=N^(1/4) se compara
    # elevando a la potencia 20: (2^j)^20*N <= N^5.
    j = 0
    while (2 ** (j + 1)) ** 20 * N <= N**5:
        j += 1
    # Para la verificacion finita usamos los extremos reales solo al enumerar.
    P = (2**j) * (N ** (1 / 20))
    Y = N ** (1 / 4)
    return j, P, Y


def main():
    checks = 0
    for k in range(100, 500):
        N = 10_000 + k * k * k
        j, P, Y = last_block(N)
        assert P <= Y + 1e-10
        assert Y < 2 * P + 1e-10
        if j:
            prev = P / 2
            assert 2 * prev <= Y + 1e-10
        # Residuo entero, con todos los enteros en lugar de solo primos.
        lo, hi = floor(P) + 1, floor(Y)
        exact = sum(N // n for n in range(lo, hi + 1))
        count = max(0, hi - lo + 1)
        coarse = (N / P) * count if count else 0.0
        assert exact <= coarse + count + 1e-7
        # Version enteramente discreta usada como respaldo de redondeo.
        assert exact <= count * (N // lo) if count else exact == 0
        checks += 1
    print(f"PASS: {checks} valores; particion y cota del residuo verificadas")


if __name__ == "__main__":
    main()
