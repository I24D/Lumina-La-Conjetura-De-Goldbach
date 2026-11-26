#!/usr/bin/env python3
"""balance_peso_nuevo_doc54.py -- Doc 55, Etapa 2.
Balance del peso nuevo del doc 54: B = sum c(q,r) a_q b_r (46.5),
la suma que multiplica a la serie singular S_N (y al termino
S(N)*N/log^2 N). Marco: docs 46 (46.1)-(46.6), 52 (constante de
Chen 0.6706 = 2.6408 - 3.9404/2, [L]), 54 (peso por atomos).
Modelo: N = 10^12 (doc 52, modelo A); primos representativos
q_i = menor primo >= N^{alpha_i}; a_q = 1/(q-2) (46.1).
Control de escala: N = 10^7 (modelo B del doc 52) para el
decaimiento de B con N. v = 10 y control v = 8.
La diagonal q=r entra una vez (pares ordenados, doc 54 sec. 1):
se reporta su parte D y se valida con la identidad de atomos
(Prop 46.2), que es exacta y la incluye automaticamente.
"""
from fractions import Fraction
import bisect
import math

CHEN = Fraction(26408, 10000) - Fraction(39404, 20000)  # 0.6706


def primes_upto(n):
    bs = bytearray(b"\x01") * (n + 1)
    bs[0:2] = b"\x00\x00"
    for i in range(2, int(n ** 0.5) + 1):
        if bs[i]:
            bs[i * i:n + 1:i] = b"\x00" * ((n - i * i) // i + 1)
    return [i for i in range(n + 1) if bs[i]]


PRIMES = primes_upto(10 ** 6 + 5000)


def rep_prime(alpha, E):
    """Menor primo >= (10**E)**alpha (float solo para elegirlo)."""
    tgt = 10 ** (E * float(alpha))
    return PRIMES[bisect.bisect_left(PRIMES, math.ceil(tgt))]


def mesh(v, h):
    lo = Fraction(1, v)
    hi = Fraction(1, 2)
    step = (hi - lo) / (h + 1)
    return [lo + (t + 1) * step for t in range(h)], step


def psi(a, b):
    return (a * b) ** 2 + (a * b) ** 3


def build_c(v, h):
    alpha, s = mesh(v, h)
    c = [[Fraction(0)] * h for _ in range(h)]
    for i in range(h - 1):
        for j in range(h - 1):
            w = psi(alpha[i], alpha[j])
            c[i][j] += w
            c[i][j + 1] -= w
            c[i + 1][j] -= w
            c[i + 1][j + 1] += w
    return alpha, c


def balance(v, h, E=12):
    alpha, c = build_c(v, h)
    q = [rep_prime(a, E) for a in alpha]
    # La malla fina puede mapear dos alpha seguidos al mismo primo
    # (huecos entre primos pequenos): esos pasos aportan 0 al balance,
    # no lo cambian de signo. Se cuentan y se reportan.
    assert all(q[i] <= q[i + 1] for i in range(h - 1)), "primos no monotonos"
    n_degen = sum(1 for i in range(h - 1) if q[i] == q[i + 1])
    assert any(q[i] < q[i + 1] for i in range(h - 1)), "sin pasos estrictos"
    a = [Fraction(1, qq - 2) for qq in q]
    # (i) margenes cero (46.4) y suma total cero
    R = [sum(c[i][j] for j in range(h)) for i in range(h)]
    Cc = [sum(c[i][j] for i in range(h)) for j in range(h)]
    assert all(x == 0 for x in R) and all(x == 0 for x in Cc)
    assert sum(c[i][j] for i in range(h) for j in range(h)) == 0
    # (ii) B directo: pares ordenados, diagonal (i,i) una vez
    B = sum(c[i][j] * a[i] * a[j] for i in range(h) for j in range(h))
    # (iii) B por atomos (Prop 46.2): cada atomo aporta > 0
    Ba = sum(psi(alpha[i], alpha[j]) * (a[i] - a[i + 1]) * (a[j] - a[j + 1])
             for i in range(h - 1) for j in range(h - 1))
    assert B == Ba, "directo != atomos"
    assert B > 0, "B no positivo"
    # (iv) diagonal y normas
    D = sum(c[i][i] * a[i] * a[i] for i in range(h))
    M1 = sum(abs(c[i][j]) for i in range(h) for j in range(h))
    Q0 = q[0]
    bound = M1 / (Fraction(Q0 - 2) * Fraction(Q0 - 2))  # cota (46.6)
    return {"B": B, "D": D, "M1": M1, "bound": bound, "Q0": Q0,
            "n_degen": n_degen,
            "r_M1": B / M1, "r_bound": B / bound, "r_diag": D / B,
            "r_chen": B / CHEN}


def main():
    print("B = sum c(q,r) a_q b_r (46.5); Chen ref = 0.6706")
    for E in (12, 7):
        print(f"-- N = 10^{E} --")
        for v in (10, 8):
            for h in (12, 24, 48, 96):
                r = balance(v, h, E)
                print(f"  v={v} h={h}: B={float(r['B']):.6e} "
                      f"(exacto {r['B'].numerator}/{r['B'].denominator}) "
                      f"B/||c||_1={float(r['r_M1']):.4e} "
                      f"B/cota46.6={float(r['r_bound']):.4e} "
                      f"diag/B={float(r['r_diag']):.4f} "
                      f"B/0.67={float(r['r_chen']):.4e} "
                      f"Q0={r['Q0']} degen={r['n_degen']}")
    print("-- control de convergencia h=192, N = 10^12 --")
    for v in (10, 8):
        r = balance(v, 192, 12)
        print(f"  v={v} h=192: B={float(r['B']):.6e} "
              f"B/||c||_1={float(r['r_M1']):.4e} "
              f"B/cota46.6={float(r['r_bound']):.4e} "
              f"diag/B={float(r['r_diag']):.4f} "
              f"degen={r['n_degen']} Q0={r['Q0']}")
    print("PASS: B > 0 exacto en los 16 casos + 2 de convergencia; "
          "B directo == B por atomos; margenes cero; suma total cero; "
          "diagonal una vez.")


if __name__ == "__main__":
    main()
