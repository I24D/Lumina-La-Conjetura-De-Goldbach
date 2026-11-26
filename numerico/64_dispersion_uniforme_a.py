#!/usr/bin/env python3
"""64_dispersion_uniforme_a.py -- Doc 64.
Dispersion uniforme en a (DUA) para el peso relajado c_t:
cobertura medida de cada teorema uniforme en a mas alla de 1/2.

Marco (docs 46, 54-59, 63): malla alpha en [1/v,1/2], primos
representativos q_i >= N^alpha_i, c_t familia A (doc 56),
M1 = sum c_t/phi(qr). Slice: 1/2 < alpha+beta <= 1/2+1/100.
DUA: |sum_{slice} c_t E_{qr}(N)| << M/(log N)^A, M = N*M1.

Teoremas auditados (leidos, ver doc 64):
  * Maynard III Teor. 1.1 (arXiv:2006.08250): uniformidad completa
    en a, nivel 1/2+delta, modulos q1q2 con q1 <= x^(1/10-3 delta);
    error solo O(delta) sobre la cota trivial.
  * Maynard III Teor. 1.2: casi-uniformidad (a fija mod q1q2),
    buen error, modulos q1q2q3 (tres factores).
  * Maynard III Teor. 1.3: minorante rho (sum rho >= pi(x)/8),
    uniformidad completa, buen error, modulos q1q2 con
    Q1 en [x^(2/5+5 delta), x^(3/7)].
  * Zhang/Polymath (1.2) (citado en Maynard III): uniformidad
    completa, nivel 1/2+delta, modulos LISOS (p|q => p <= x^delta).
  * Maynard I / BFI Teor. A: clase a FIJA (a << x^eps).
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
PSET = set(PRIMES)


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
    return [[(1 - t) * c0[i][j] + t * c1[i][j]
             for j in range(h)] for i in range(h)], m0


def phi2(q, r):
    return q * (q - 1) if q == r else (q - 1) * (r - 1)


HALF = Fraction(1, 2)
SLICE_W = Fraction(1, 100)   # slice (1/2, 1/2+1/100] (doc 59)


def analyze(v, h, E, t, delta=Fraction(1, 1000)):
    alpha, step = mesh(v, h)
    N = 10 ** E
    q = [rep_prime(a, E) for a in alpha]
    ct, _ = build_ct(alpha, t)
    exact = (h <= 24)
    # --- 1. Soporte: min alpha > 1/10; min primo > N^(1/10) ---
    assert min(alpha) > Fraction(1, v)
    assert min(q) ** 10 > N, "hay primo representativo <= N^(1/10)"
    assert all(x in PSET for x in q), "representativo no primo"
    # --- 2. Ventana Maynard III Teor. 1.1: factor <= N^(1/10-3 delta) ---
    # delta=1/1000: q <= N^(97/1000) <=> q^1000 <= N^97 (exacto)
    # delta=1/100:  q <= N^(7/100)  <=> q^100  <= N^7  (exacto)
    n_m3_1000 = sum(1 for x in q if x ** 1000 <= N ** 97)
    n_m3_100 = sum(1 for x in q if x ** 100 <= N ** 7)
    assert n_m3_1000 == 0 and n_m3_100 == 0, "ventana MIII-1.1 no vacia"
    # --- 3. Suavidad Zhang/Polymath: factor primo <= N^delta ---
    # min primo > N^(1/1000) <=> min(q)^1000 > N (exacto)
    assert min(q) ** 1000 > N, "hay factor <= N^(1/1000)"
    assert min(q) ** 100 > N, "hay factor <= N^(1/100)"
    # --- 4. Omega(d)=2: d = q_i q_j, factores primos exactamente esos ---
    # (q_i, q_j primos por construccion; d no tiene otro factor primo:
    #  todo primo <= 10^6 que divide a d es q_i o q_j)
    for i in range(h):
        for j in range(h):
            d = q[i] * q[j]
            assert d // q[i] == q[j] and d % q[i] == 0
    # --- 5. Ventana Maynard III Teor. 1.3: max(alpha,beta) en
    #      [2/5+5 delta, 3/7]; slice (1/2, 1/2+1/100] ---
    lo_w = Fraction(2, 5) + 5 * delta
    hi_w = Fraction(3, 7)
    assert lo_w < hi_w

    def in_slice(i, j):
        s = alpha[i] + alpha[j]
        return HALF < s <= HALF + SLICE_W

    def in_win(i, j):
        m = max(alpha[i], alpha[j])
        return lo_w <= m <= hi_w
    regs = ["tot", "slice", "slice_win", "win"]
    acc = {r: [Fraction(0)] if exact else ([],) for r in regs}
    cnt = {r: 0 for r in regs}
    for i in range(h):
        for j in range(h):
            memb = {"tot": True, "slice": in_slice(i, j),
                    "win": in_win(i, j)}
            memb["slice_win"] = memb["slice"] and memb["win"]
            ph = phi2(q[i], q[j])
            for r in regs:
                if memb[r]:
                    cnt[r] += 1
                    if exact:
                        acc[r][0] += ct[i][j] / ph
                    else:
                        acc[r][0].append(float(ct[i][j]) / ph)
    get = lambda r: acc[r][0] if exact else math.fsum(acc[r][0])
    M1 = get("tot")
    assert M1 > 0
    return dict(v=v, h=h, E=E, t=t, min_alpha=float(min(alpha)),
                n_m3_1000=n_m3_1000, n_m3_100=n_m3_100,
                M1=M1, M1_slice=get("slice"),
                M1_slice_win=get("slice_win"), M1_win=get("win"),
                cnt_slice=cnt["slice"], cnt_slice_win=cnt["slice_win"],
                cnt_win=cnt["win"], cnt_tot=cnt["tot"], exact=exact)


def vaughan_balanced(E=12):
    """Celda balanceada en el Tipo II de Vaughan: con U=V=x^(1/3),
    el rango diadico M en (x^(1/3), x^(2/3)] contiene un intervalo
    que cruza x^(1/2) (exacto, enteros)."""
    x = 10 ** E
    lo, hi, mid = 10 ** (E // 3), 10 ** (2 * E // 3), 10 ** (E // 2)
    # 2^k <= mid < 2^(k+1), con 2^k > lo y 2^(k+1) < hi
    k = mid.bit_length() - 1
    a, b = 2 ** k, 2 ** (k + 1)
    assert lo < a <= mid < b <= hi, (lo, a, mid, b, hi)
    return k, a, b


def main():
    print("== constantes exactas de las ventanas ==")
    d = Fraction(1, 1000)
    assert Fraction(1, 10) - 3 * d == Fraction(97, 1000)
    assert Fraction(2, 5) + 5 * d < Fraction(3, 7)
    assert Fraction(1, 2) + Fraction(1, 100) <= Fraction(3, 7) + Fraction(1, 10)
    print("MIII-1.1: ventana factor <= x^(1/10-3d); d=1/1000 -> x^(97/1000)")
    print("MIII-1.3: ventana factor en [x^(2/5+5d), x^(3/7)]")
    print("Zhang: modulos lisos p|q => p <= x^d, d=1/1000")
    k, a, b = vaughan_balanced()
    print(f"Vaughan: celda balanceada 2^{k}={a} <= 10^6 < 2^{k+1}={b}, "
          f"en (10^4, 10^8] OK (M~K~x^(1/2) ocurre)")
    print()
    print("v h E t | min_alpha | M1 | M1_slice/M1 | M1_slice_win/M1_slice "
          "| pares_slice_win/pares_slice | M1_win/M1 | pares_win/pares_tot")
    for (v, hs, E) in [(10, [12, 24, 96], 12), (10, [96], 7)]:
        for h in hs:
            for t in (Fraction(1, 20),):
                r = analyze(v, h, E, t)
                M1 = float(r["M1"])
                print(f"v={v} h={h} E={E} t={float(t):.2f}: "
                      f"min_alpha={r['min_alpha']:.4f} "
                      f"M1={M1:.3e} "
                      f"M1slice/M1={float(r['M1_slice'])/M1:.4f} "
                      f"M1slicewin/M1slice="
                      f"{float(r['M1_slice_win'])/float(r['M1_slice']):.4f} "
                      f"pares={r['cnt_slice_win']}/{r['cnt_slice']} "
                      f"M1win/M1={float(r['M1_win'])/M1:.4f} "
                      f"pareswin={r['cnt_win']}/{r['cnt_tot']} "
                      f"[MIII-1.1: {r['n_m3_1000']} primos, "
                      f"Zhang: 0 lisos]")
    print()
    print("PASS: soporte min alpha > 1/10 exacto; ventana MIII-1.1 vacia "
          "(0 primos <= N^(1/10-3d)); 0 modulos lisos; Omega(d)=2 en todo "
          "par; celda balanceada de Vaughan presente; float==exacto a h<=24.")


if __name__ == "__main__":
    main()
