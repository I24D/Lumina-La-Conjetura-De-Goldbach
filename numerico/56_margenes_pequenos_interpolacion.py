#!/usr/bin/env python3
"""56_margenes_pequenos_interpolacion.py -- Doc 56.
Relajacion de (46.4): margenes PEQUENOS pero no cero.
Familia: c_t = (1-t)*c0/||c0||_1 + t*c1/||c1||_1, con
c0 = peso de atomos del doc 54 (psi0 = (ab)^2+(ab)^3, o
psi1 = (ab)+(ab)^2, variante del doc 55) y c1 = peso de
primer orden (constante 1, o rama favorable de Chen c1 del
doc 52). Marco exacto (Fraction): malla del doc 52,
N = 10^12, a_q = 1/(q-2), primos representativos segun el
doc 55 (sin forzado estricto; degenerados contados).
Medidas por (46.3): D1, D2, sum c, P1 = sum a_q R_q +
sum b_r C_r, B = sum c a_q b_r; cocientes |P1|/sum c,
B/sum c, B/D, sum c/D con D = max(D1,D2).
Lema 56.1 [T]: para todo peso con sum c > 0,
sum c <= D1*||c||_1 y sum c <= D2*||c||_1 (desigualdad
triangular: sum|R_i| >= |sum R_i| = |sum c|). Normalizado
a ||c||_1 = 1, sum c <= min(D1,D2): el trade-off lineal
es OPTIMO y la interpolacion lo alcanza (sum c = D1 = D2
= t, salvo la correccion por ||c_t||_1 < 1).
"""
from fractions import Fraction
import bisect
import math

CHEN = Fraction(6706, 10000)  # 0.6706, constante de Chen (doc 52)


def primes_upto(n):
    bs = bytearray(b"\x01") * (n + 1)
    bs[0:2] = b"\x00\x00"
    for i in range(2, int(n ** 0.5) + 1):
        if bs[i]:
            bs[i * i:n + 1:i] = b"\x00" * ((n - i * i) // i + 1)
    return [i for i in range(n + 1) if bs[i]]


PRIMES = primes_upto(10 ** 6 + 5000)


def rep_prime(alpha, E=12):
    """Menor primo >= (10**E)**alpha (float solo para elegirlo)."""
    tgt = 10 ** (E * float(alpha))
    return PRIMES[bisect.bisect_left(PRIMES, math.ceil(tgt))]


def mesh(v, h):
    lo = Fraction(1, v)
    hi = Fraction(1, 2)
    step = (hi - lo) / (h + 1)
    return [lo + (t + 1) * step for t in range(h)], step


def psi0(a, b):
    return (a * b) ** 2 + (a * b) ** 3


def psi1(a, b):
    return (a * b) + (a * b) ** 2


def build_atoms(alpha, psi):
    """Peso del doc 54: suma de atomos rectangulares A[i][j] con
    coeficiente w[i][j] = psi(alpha_i, alpha_j)."""
    h = len(alpha)
    c = [[Fraction(0)] * h for _ in range(h)]
    for i in range(h - 1):
        for j in range(h - 1):
            w = psi(alpha[i], alpha[j])
            c[i][j] += w
            c[i][j + 1] -= w
            c[i + 1][j] -= w
            c[i + 1][j + 1] += w
    return c


def n3p_bisect(alpha, i, j):
    """#{k>j : alpha_i+alpha_j+alpha_k <= 1} (bisect; alpha creciente)."""
    lim = Fraction(1) - alpha[i] - alpha[j]
    K = bisect.bisect_right(alpha, lim) - 1
    return max(0, K - j)


def n3p_directo(alpha, i, j):
    """Misma cuenta por triple bucle (solo autocheck h=12)."""
    return sum(1 for k in range(j + 1, len(alpha))
               if alpha[i] + alpha[j] + alpha[k] <= 1)


def build_chen1(alpha, b):
    """Rama favorable de Chen (doc 52): c1[i][j] = (1+N3p) si
    i<j, j>=b; 0 en el resto. Entradas >= 0."""
    h = len(alpha)
    c = [[Fraction(0)] * h for _ in range(h)]
    for i in range(h):
        for j in range(i + 1, h):
            if j >= b:
                n = n3p_bisect(alpha, i, j)
                if i >= b:
                    assert n == 0, (i, j, n)  # ec. (3) de [L]
                c[i][j] = Fraction(1 + n)
    return c


def build_const(h):
    return [[Fraction(1)] * h for _ in range(h)]


def l1(M):
    return sum(abs(M[i][j]) for i in range(len(M)) for j in range(len(M)))


def normalize(M):
    m = l1(M)
    assert m > 0
    return [[M[i][j] / m for j in range(len(M))] for i in range(len(M))], m


def measure(c, a):
    """Medidas exactas por (46.3)."""
    h = len(c)
    R = [sum(c[i][j] for j in range(h)) for i in range(h)]
    C = [sum(c[i][j] for i in range(h)) for j in range(h)]
    T = sum(R)
    M1 = sum(abs(c[i][j]) for i in range(h) for j in range(h))
    D1 = Fraction(sum(abs(x) for x in R), M1)
    D2 = Fraction(sum(abs(x) for x in C), M1)
    P1 = sum(a[i] * R[i] for i in range(h)) \
        + sum(a[j] * C[j] for j in range(h))
    B = sum(c[i][j] * a[i] * a[j] for i in range(h) for j in range(h))
    # Lema 56.1: sum c <= D*||c||_1 (desigualdad triangular)
    assert T <= D1 * M1 and T <= D2 * M1, "Lema 56.1 falla"
    return dict(T=T, M1=M1, D1=D1, D2=D2, P1=P1, B=B, R=R, C=C)


def rep_primes(alpha, E=12):
    q = [rep_prime(a, E) for a in alpha]
    assert all(q[i] <= q[i + 1] for i in range(len(q) - 1))
    n_degen = sum(1 for i in range(len(q) - 1) if q[i] == q[i + 1])
    return q, [Fraction(1, qq - 2) for qq in q], n_degen


def fmt(x):
    return f"{float(x):.4e}"


def main():
    # --- autocheck: N3p por bisect == triple bucle (h=12) ---
    al12, _ = mesh(10, 12)
    for i in range(12):
        for j in range(12):
            assert n3p_bisect(al12, i, j) == n3p_directo(al12, i, j)
    print("autocheck N3p: bisect == triple bucle (h=12) OK")
    # --- referencias cruzadas docs 54/55 (B y ||c||_1, v=10,8 h=96,12) ---
    for (v, h, Bref, Mref) in [(10, 96, 1.404e-06, 0.2905),
                               (8, 96, 6.474e-07, 0.2918),
                               (10, 12, None, 0.1763),
                               (8, 12, None, 0.1830)]:
        al, _ = mesh(v, h)
        c0 = build_atoms(al, psi0)
        m0 = measure(c0, [Fraction(0)] * h)
        assert m0["D1"] == 0 and m0["D2"] == 0 and m0["T"] == 0  # (46.4)
        assert abs(float(m0["M1"]) - Mref) / Mref < 1e-3, (v, h, m0["M1"])
        if Bref is not None:
            q, a, _ = rep_primes(al)
            mb = measure(c0, a)
            assert abs(float(mb["B"]) - Bref) / Bref < 1e-3, (v, h, mb["B"])
    print("cross-check docs 54/55: ||c0||_1 y B(psi0) reproducidos OK")
    # --- psi1: B(v=10,h=96) ~ 6.570e-05 (doc 55 sec. 5a) ---
    al, _ = mesh(10, 96)
    q, a, _ = rep_primes(al)
    assert abs(float(measure(build_atoms(al, psi1), a)["B"])
               - 6.570e-05) / 6.570e-05 < 1e-3
    print("cross-check doc 55 sec. 5a: B(psi1, v=10, h=96) OK")
    # --- endpoints de primer orden: D1 = D2 = 1 exactos ---
    for (v, h) in [(10, 96), (8, 48)]:
        al, _ = mesh(v, h)
        b = sum(1 for x in al if x < Fraction(1, 3))
        q, a, _ = rep_primes(al)
        for tag, M in [("const", build_const(h)), ("chen1", build_chen1(al, b))]:
            m = measure(M, a)
            assert m["D1"] == 1 and m["D2"] == 1, (tag, v, h)
            assert m["T"] == m["M1"], (tag, v, h)  # entradas >= 0
    print("cross-check doc 52: D1 = D2 = 1 en endpoints const y chen1 OK")
    print()
    print("t | sum c | D1 | D2 | |P1| | |P1|/sum c | B | B/sum c | sum c/D")
    print("Chen ref = 0.6706; D = max(D1,D2); v=10 salvo control")
    best = None  # (T, etiqueta) con D<=1/10 y |P1|/T<=1/20, max T
    for (v, hs) in [(10, [12, 24, 48, 96]), (8, [48])]:
        for h in hs:
            al, _ = mesh(v, h)
            b = sum(1 for x in al if x < Fraction(1, 3))
            q, a, nd = rep_primes(al)
            c0a, M0a = normalize(build_atoms(al, psi0))
            c0b, M0b = normalize(build_atoms(al, psi1))
            c1c, M1c = normalize(build_const(h))
            c1k, M1k = normalize(build_chen1(al, b))
            fams = [("A psi0+const", c0a, c1c, measure(c0a, a),
                     measure(c1c, a)),
                    ("B psi1+const", c0b, c1c, measure(c0b, a),
                     measure(c1c, a)),
                    ("C psi0+chen1", c0a, c1k, measure(c0a, a),
                     measure(c1k, a))]
            for t in [Fraction(1, 100), Fraction(1, 20), Fraction(1, 10),
                      Fraction(1, 5), Fraction(1, 2)]:
                for (tag, c0n, c1n, m0, m1) in fams:
                    ct = [[(1 - t) * c0n[i][j] + t * c1n[i][j]
                           for j in range(h)] for i in range(h)]
                    m = measure(ct, a)
                    # linealidad exacta: componentes se suman
                    assert m["T"] == (1 - t) * m0["T"] + t * m1["T"]
                    assert m["P1"] == (1 - t) * m0["P1"] + t * m1["P1"]
                    assert m["B"] == (1 - t) * m0["B"] + t * m1["B"]
                    D = max(m["D1"], m["D2"])
                    rP = abs(m["P1"]) / m["T"] if m["T"] != 0 else Fraction(0)
                    rB = m["B"] / m["T"] if m["T"] != 0 else Fraction(0)
                    rD = m["T"] / D if D != 0 else Fraction(0)
                    print(f"v={v} h={h} {tag} t={float(t):.2f}: "
                          f"sumc={fmt(m['T'])} D1={fmt(m['D1'])} "
                          f"D2={fmt(m['D2'])} |P1|={fmt(abs(m['P1']))} "
                          f"|P1|/sumc={fmt(rP)} B={fmt(m['B'])} "
                          f"B/sumc={fmt(rB)} sumc/D={fmt(rD)} "
                          f"B/D={fmt(m['B'] / D if D != 0 else 0)}")
                    if (m["T"] > 0 and m["D1"] <= Fraction(1, 10)
                            and m["D2"] <= Fraction(1, 10)
                            and rP <= Fraction(1, 20)
                            and (best is None or m["T"] > best[0])):
                        best = (m["T"], f"v={v} h={h} {tag} "
                                f"t={float(t):.2f} D1={fmt(m['D1'])} "
                                f"|P1|/sumc={fmt(rP)}")
    print()
    if best is not None:
        print(f"BEST (max sum c con D1,D2<=1/10 y |P1|/sum c<=1/20): "
              f"sum c={fmt(best[0])} -- {best[1]}")
    else:
        print("BEST: ningun candidato cumple D<=1/10 y |P1|/sum c<=1/20")
    print("PASS: interpolaciones medidas exactamente (Fraction); "
          "linealidad de T, P1, B verificada; Lema 56.1 verificado; "
          "cross-checks docs 52/54/55 OK.")


if __name__ == "__main__":
    main()
