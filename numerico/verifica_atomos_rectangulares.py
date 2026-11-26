#!/usr/bin/env python3
"""verifica_atomos_rectangulares.py -- Doc 54, Etapa 1.
Atom os rectangulares de un peso NUEVO con margenes cero por construccion.
Marco de medicion exacto de los docs 46-48 y 52 (Fraction, malla alpha).
El peso es c = suma_{i,j} w[i][j] * A[i][j], donde A[i][j] es el atomo
rectangular (+1,-1,-1,+1) en las esquinas (i,j),(i,j+1),(i+1,j),(i+1,j+1)
y w[i][j] = psi(alpha_i, alpha_j) con psi(a,b) = (a*b)**2 + (a*b)**3.
Cada atomo tiene margenes cero; luego c tambien (criterio 46.4).
psi es genuinamente bidimensional: su parcial mixta es ab(4+9ab), ni
aditiva ni separable (un solo termino separable daria c de rango 1).
"""
from fractions import Fraction


def mesh(v, h):
    lo = Fraction(1, v)
    hi = Fraction(1, 2)
    step = (hi - lo) / (h + 1)
    return [lo + (t + 1) * step for t in range(h)], step


def psi(a, b):
    return (a * b) ** 2 + (a * b) ** 3


def build(v, h):
    alpha, s = mesh(v, h)
    c = [[Fraction(0)] * h for _ in range(h)]
    for i in range(h - 1):
        for j in range(h - 1):
            w = psi(alpha[i], alpha[j])
            c[i][j] += w
            c[i][j + 1] -= w
            c[i + 1][j] -= w
            c[i + 1][j + 1] += w
    return alpha, s, c


def margins(c, h):
    R = [sum(c[i][j] for j in range(h)) for i in range(h)]
    C = [sum(c[i][j] for i in range(h)) for j in range(h)]
    return R, C


def check(v, h):
    alpha, s, c = build(v, h)
    R, C = margins(c, h)
    M1 = sum(abs(c[i][j]) for i in range(h) for j in range(h))
    # (a) criterio 46.4: todos los margenes cero (exacto)
    assert all(x == 0 for x in R) and all(x == 0 for x in C)
    # (b) no trivial: masa L1 positiva
    assert M1 > 0
    # (c) formula interior exacta de c(i,j)
    for i in range(1, h - 1):
        for j in range(1, h - 1):
            f = (2 * alpha[i] - s) * (2 * alpha[j] - s)
            g = (3 * alpha[i] ** 2 - 3 * alpha[i] * s + s ** 2) * \
                (3 * alpha[j] ** 2 - 3 * alpha[j] * s + s ** 2)
            assert c[i][j] == s * s * (f + g), (i, j)
    # (d) no aditiva: diferencias mixtas interiores adyacentes no nulas
    md = md_tot = 0
    for i in range(1, h - 2):
        for j in range(1, h - 2):
            md_tot += 1
            if c[i][j] - c[i + 1][j] - c[i][j + 1] + c[i + 1][j + 1] != 0:
                md += 1
    assert md == md_tot, (md, md_tot)
    # (e) no rango 1: menores 2x2 interiores adyacentes no nulos
    mn = mn_tot = 0
    for i in range(1, h - 2):
        for j in range(1, h - 2):
            mn_tot += 1
            if c[i][j] * c[i + 1][j + 1] - c[i][j + 1] * c[i + 1][j] != 0:
                mn += 1
    assert mn == mn_tot, (mn, mn_tot)
    npos = sum(1 for i in range(h) for j in range(h) if c[i][j] > 0)
    nneg = sum(1 for i in range(h) for j in range(h) if c[i][j] < 0)
    print(f"  v={v} h={h}: margenes cero (46.4 SATISFECHO), ||c||_1={M1}, "
          f"entradas +:{npos} -:{nneg}, dif.mixtas {md}/{md_tot}, "
          f"menores 2x2 {mn}/{mn_tot}")
    return True


def main():
    for v in (10, 8):
        for h in (12, 24, 48, 96):
            check(v, h)
    print("PASS: atomos rectangulares, margenes cero exactos (46.4), "
          "peso genuinamente bidimensional (no aditivo, no rango 1), "
          "v=10,8; h=12,24,48,96.")


if __name__ == "__main__":
    main()
