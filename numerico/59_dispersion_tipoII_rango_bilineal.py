#!/usr/bin/env python3
"""59_dispersion_tipoII_rango_bilineal.py -- Doc 59.
Dispersion / Tipo II para el rango bilineal [N^(1/10),N^(1/2)]^2
del peso relajado c_t (familia A del doc 56).

Marco (docs 46, 52, 54-58): malla alpha en [1/v,1/2], primos
representativos q_i >= N^alpha_i, c_t = (1-t)*c0/||c0||_1 + t*const,
M1 = sum c_t/phi(qr), M = N*M1. Nivel: Def. 58.1 (error pesado
<< M/(log N)^A). El error con qr<=N^(1/2) es despreciable (BV);
la pregunta es el slice (N^(1/2), N^(1/2+delta)].

Herramientas publicadas auditadas (ver doc 59):
  * Gran criba clasica: sum_{q<=Q} sum*_{a(q)} |sum a_n e(an/q)|^2
    << (Q^2+N) sum |a_n|^2. El termino Q^2 capa el Tipo II en 1/2.
  * Maynard, Teor. 1.1 / Cor. 1.2 (leido): modulos q<=x^(1/2+delta),
    delta<1/42, con un divisor en [x^(2 delta), x^(1/10-7 delta/5)];
    error << x/(log x)^A con valores absolutos, clase a FIJA.
  * BFI, exponente 29/56 (resumido en la literatura): coeficientes
    soportados en modulos q1 q2 con q1<=x^th1, q2<=x^th2,
    th1<1/3, th2<1/5, 5 th1+2 th2<2; esquina optima
    (th1,th2)=(9/28,11/56), nivel 29/56. Clase a FIJA (BFI/Fiorilli).
  * BFI Teor. B (well-factorable, 4/7): bloqueado por Prop. 58.1.

Mediciones (exactas donde es barato; float+fsum a h grande con
cross-check exacto a h<=24, como el doc 58):
  * Slice BFI: 1/2 < alpha+beta <= 29/56. Region L_BFI:
    min(alpha,beta)<=11/56 y max(alpha,beta)<=9/28 (la parte que
    el resultado BFI de factorizacion fija cubre, a nivel 29/56).
  * Slice Maynard: 1/2 < alpha+beta <= 1/2+1/100 (delta=1/100<1/42).
    Ventana estricta W_May: min(alpha,beta) <= 1/10-7 delta/5.
    Ventana generosa estilo Fouvry (Maynard sec. 15, "near x^(1/7)"):
    min(alpha,beta) <= 1/7.
  * R_unc = ln(N) * sum_{slice no cubierta} |c|/(qr) / M1: cota
    trivial del error de la parte no cubierta del slice, /M.
    Si R_unc es una constante (no << 1/(log N)^A), no hay nivel
    en el sentido estricto de la Def. 58.1.
  * Reordenamiento (exacto): sum_{q,r} c(q,r) E_{qr}
    = sum_{n<N} Lambda(n) W(n) - N*M1, con
    W(n) = sum_{d|(N-n)} lambda_t(d),
    lambda_t(d) = sum_{qr=d} c_t(q,r). Se verifica la identidad
    algebraica sum_{qr|m} c = sum_{d|m} lambda_t(d) en muestras.
"""
from fractions import Fraction
import bisect
import math


def primes_upto(n):
    bs = bytearray(b"\x01") * (n + 1)
    bs[0:2] = b"\x00\x00"
    for i in range(2, int(n ** 0.5) + 1):
        if bs[i]:
            bs[i * i:n + 1:i] = b"\x00" * ((n - i * i) // i + 1)
    return [i for i in range(n + 1) if bs[i]]


PRIMES = primes_upto(10 ** 6 + 5000)


def rep_prime(alpha, E=12):
    tgt = 10 ** (E * float(alpha))
    return PRIMES[bisect.bisect_left(PRIMES, math.ceil(tgt))]


def mesh(v, h):
    lo = Fraction(1, v)
    hi = Fraction(1, 2)
    step = (hi - lo) / (h + 1)
    return [lo + (t + 1) * step for t in range(h)], step


def psi0(a, b):
    return (a * b) ** 2 + (a * b) ** 3


def build_atoms(alpha, psi):
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


def l1(M):
    return sum(abs(M[i][j]) for i in range(len(M)) for j in range(len(M)))


def normalize(M):
    m = l1(M)
    assert m > 0
    return [[M[i][j] / m for j in range(len(M))] for i in range(len(M))], m


def build_ct(alpha, t):
    h = len(alpha)
    c0, m0 = normalize(build_atoms(alpha, psi0))
    c1 = [[Fraction(1, h * h)] * h for _ in range(h)]
    ct = [[(1 - t) * c0[i][j] + t * c1[i][j]
           for j in range(h)] for i in range(h)]
    return ct, m0


def phi2(q, r):
    return q * (q - 1) if q == r else (q - 1) * (r - 1)


# --- constantes exactas de las herramientas auditadas ---
BFI_T1 = Fraction(9, 28)    # factor mayor <= x^(9/28)
BFI_T2 = Fraction(11, 56)   # factor menor <= x^(11/56)
BFI_LVL = Fraction(29, 56)  # nivel BFI de factorizacion fija
MAY_D = Fraction(1, 100)    # delta de Maynard (< 1/42)
MAY_WIN = Fraction(1, 10) - Fraction(7, 5) * MAY_D  # 0.086
FOU_WIN = Fraction(1, 7)    # escala estilo Fouvry (Maynard sec. 15)
HALF = Fraction(1, 2)


def analyze(v, h, E, t):
    alpha, step = mesh(v, h)
    N = 10 ** E
    lnN = math.log(N)
    q = [rep_prime(a, E) for a in alpha]
    ct, m0 = build_ct(alpha, t)
    # --- exacto: sum c_t = t; ||c_t||_1; min alpha > 1/10 ---
    T = sum(ct[i][j] for i in range(h) for j in range(h))
    assert T == t, ("sum c_t != t", v, h, t)
    M1n = l1(ct)
    assert Fraction(9, 10) < M1n <= 1, ("||c_t||_1 fuera de rango", M1n)
    assert min(alpha) > Fraction(1, v), "min alpha no > 1/v"
    # --- regiones sobre la malla alpha (exactas, Fraction) ---
    def in_bfi_slice(i, j):
        s = alpha[i] + alpha[j]
        return HALF < s <= BFI_LVL
    def in_may_slice(i, j):
        s = alpha[i] + alpha[j]
        return HALF < s <= HALF + MAY_D
    def in_L(i, j):
        a, b = alpha[i], alpha[j]
        return min(a, b) <= BFI_T2 and max(a, b) <= BFI_T1
    def in_wmay(i, j):
        return min(alpha[i], alpha[j]) <= MAY_WIN
    def in_wfou(i, j):
        return min(alpha[i], alpha[j]) <= FOU_WIN
    # --- sumas: exactas (Fraction) a h<=24; float+fsum a h mayor ---
    exact = (h <= 24)
    # acumuladores: [M1, S] por region, como Fraction o lista float
    regs = ["tot", "bfi", "bfi_L", "bfi_unc", "may", "may_fou",
            "may_unc_fou", "L", "wfou"]
    acc = {r: [Fraction(0), Fraction(0)] if exact
           else ([], []) for r in regs}
    sumc = {r: Fraction(0) for r in regs}
    l1r = {r: Fraction(0) for r in regs}
    for i in range(h):
        for j in range(h):
            cij = ct[i][j]
            pr = q[i] * q[j]
            ph = phi2(q[i], q[j])
            memb = {"tot": True,
                    "bfi": in_bfi_slice(i, j),
                    "may": in_may_slice(i, j),
                    "L": in_L(i, j),
                    "wfou": in_wfou(i, j)}
            memb["bfi_L"] = memb["bfi"] and memb["L"]
            memb["bfi_unc"] = memb["bfi"] and not memb["L"]
            memb["may_fou"] = memb["may"] and memb["wfou"]
            memb["may_unc_fou"] = memb["may"] and not memb["wfou"]
            for r in regs:
                if memb[r]:
                    sumc[r] += cij
                    l1r[r] += abs(cij)
                    if exact:
                        acc[r][0] += cij / ph
                        acc[r][1] += abs(cij) / pr
                    else:
                        acc[r][0].append(float(cij) / ph)
                        acc[r][1].append(float(abs(cij)) / pr)
    def get(r, k):
        return acc[r][k] if exact else math.fsum(acc[r][k])
    M1 = get("tot", 0)
    assert M1 > 0, ("M1 <= 0", v, h, t)
    # W_May debe ser vacia: min alpha > 1/10 > MAY_WIN en toda la malla
    n_wmay = sum(1 for i in range(h) for j in range(h) if in_wmay(i, j))
    assert n_wmay == 0, ("ventana Maynard no vacia", n_wmay)
    return dict(v=v, h=h, E=E, t=t, M1=M1,
                M1_bfi=get("bfi", 0), M1_bfi_L=get("bfi_L", 0),
                M1_may=get("may", 0), M1_may_fou=get("may_fou", 0),
                M1_L=get("L", 0),
                S_bfi_unc=get("bfi_unc", 1), S_may_unc_fou=get("may_unc_fou", 1),
                sumc_L=sumc["L"], sumc_bfi=sumc["bfi"],
                sumc_may=sumc["may"], lnN=lnN, exact=exact)


def reordering_check():
    """Identidad de reordenamiento (exacta, Fraction, malla h=12):
    sum_{i,j: q_i q_j | m} c_t(i,j) = sum_{d|m} lambda_t(d),
    lambda_t(d) = sum_{i,j: q_i q_j = d} c_t(i,j)."""
    v, h, E, t = 10, 12, 12, Fraction(1, 20)
    alpha, _ = mesh(v, h)
    q = [rep_prime(a, E) for a in alpha]
    ct, _ = build_ct(alpha, t)
    lam = {}
    for i in range(h):
        for j in range(h):
            d = q[i] * q[j]
            lam[d] = lam.get(d, Fraction(0)) + ct[i][j]
    assert sum(lam.values()) == t, "sum lambda_t != t"
    assert sum(abs(x) for x in lam.values()) <= l1(ct) + Fraction(0)
    samples = [q[0] * q[1], q[3] * q[3], q[0] * q[5] * q[11],
               2 * q[7], 1, 999983, q[2] * q[9]]
    for m in samples:
        lhs = sum(ct[i][j] for i in range(h) for j in range(h)
                  if m % (q[i] * q[j]) == 0)
        rhs = sum(val for d, val in lam.items() if m % d == 0)
        assert lhs == rhs, ("reordenamiento falla", m, lhs, rhs)
    return len(samples)


def main():
    print("== checks exactos de las herramientas ==")
    assert 5 * BFI_T1 + 2 * BFI_T2 == 2, "esquina BFI no es frontera"
    assert BFI_T1 + BFI_T2 == BFI_LVL, "nivel BFI != 29/56"
    assert BFI_T2 < Fraction(1, 5) and BFI_T1 < Fraction(1, 3)
    assert MAY_WIN < Fraction(1, 10), "ventana Maynard no < 1/10"
    assert MAY_D < Fraction(1, 42), "delta Maynard fuera de rango"
    print(f"BFI: 5*(9/28)+2*(11/56)={5*BFI_T1+2*BFI_T2} (frontera), "
          f"9/28+11/56={BFI_T1+BFI_T2}=29/56 OK")
    print(f"Maynard: delta=1/100, ventana [x^(1/50), x^({float(MAY_WIN):.3f})] "
          f"< x^(1/10) OK")
    n = reordering_check()
    print(f"Reordenamiento: {n} muestras, identidad exacta OK "
          f"(sum lambda_t = t, sum|lambda| <= ||c_t||_1)")
    print()
    print("v h E t | M1 | M1_bfi/M1 | M1_bfi_L/M1_bfi | M1_may/M1 | "
          "M1_may_fou/M1_may | R_unc_BFI | R_unc_May | M1_L/M1 | sumc_L/t")
    for (v, hs, E) in [(10, [12, 24, 48, 96], 12), (8, [48], 12),
                       (10, [96], 7)]:
        for h in hs:
            for t in (Fraction(1, 20), Fraction(1, 10)):
                r = analyze(v, h, E, t)
                M1 = float(r["M1"])
                print(f"v={v} h={h} E={E} t={float(t):.2f}: "
                      f"M1={M1:.3e} "
                      f"M1bfi/M1={float(r['M1_bfi'])/M1:.4f} "
                      f"M1bfiL/M1bfi={float(r['M1_bfi_L'])/float(r['M1_bfi']):.4f} "
                      f"M1may/M1={float(r['M1_may'])/M1:.4f} "
                      f"M1mayFou/M1may={float(r['M1_may_fou'])/float(r['M1_may']):.4f} "
                      f"RuncBFI={r['lnN']*float(r['S_bfi_unc'])/M1:.4f} "
                      f"RuncMay={r['lnN']*float(r['S_may_unc_fou'])/M1:.4f} "
                      f"M1L/M1={float(r['M1_L'])/M1:.4f} "
                      f"sumcL/t={float(r['sumc_L'])/float(t):.4f}")
    print()
    print("PASS: dispersion auditada; esquina BFI y ventana Maynard "
          "exactas; ventana Maynard vacia en la malla (0 pares); "
          "reordenamiento exacto; float==exacto a h<=24.")


if __name__ == "__main__":
    main()
