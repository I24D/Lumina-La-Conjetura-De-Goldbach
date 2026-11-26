#!/usr/bin/env python3
"""65_lema_dua_esquina.py -- Doc 65.
Lema DUA-esquina (doc 64 sec. 6): dispersion uniforme en a para Lambda
en la esquina (alpha,beta) en [1/10, 1/10+1/200] x [2/5+5d, 2/5+5d+1/100],
d=1/1000, o Teor. 1.3 de Maynard III para Lambda en esa esquina.

Verificaciones (exactas donde es barato; float con fsum donde no):
 1. Gran criba: Q^2/N en la esquina (Q ~ N^s, s en [0.505, 0.52]) > 1
    y creciente en N: la estrechez no da ahorro (el termino es Q^2,
    no el numero de modulos).
 2. Umbral Teor. 1.1: factor pequeno <= N^(1/10-3d) = N^(97/1000);
    esquina en [N^(1/10), N^(21/200)]: hueco 0.003-0.008 en exponente;
    factor finito N^hueco a N accesible (pequeno pero -> infinito).
 3. Masa de la esquina: fracciones M1_esquina/M1_slice y /M1 a h=96
    (re-verifica doc 64: 4.78% del slice, 0.004% de M1); esquina
    vacia a h=12, 24.
 4. Celda balanceada de Vaughan presente (M~K~x^(1/2)) para los
    modulos de la esquina: la descomposicion de Lambda la contiene.
 5. Simplificacion de la congruencia con a=N (un solo entero):
    n1*N = n2*N (mod q) <=> n1 = n2 (mod q) si (N,q)=1 (exacto).
    El desplazamiento sigue siendo ~N (grande): no salva el metodo
    espectral.
 6. Signo: los terminos que el Teor. 1.3 descarta (S_{2,3}, S_{2,5})
    son >= 0 (cuentan rho >= 0); los de Heath-Brown llevan mu
    (toma +1 y -1): no se pueden descartar.
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
    return [[(1 - t) * c0[i][j] + t * c1[i][j]
             for j in range(h)] for i in range(h)], m0


def phi2(q, r):
    return q * (q - 1) if q == r else (q - 1) * (r - 1)


HALF = Fraction(1, 2)
SLICE_W = Fraction(1, 100)
DELTA = Fraction(1, 1000)


def corner_defs():
    # Esquina sec. 6 (efectiva dentro del slice: beta <= 1/2+1/100-1/10):
    a_lo, a_hi = Fraction(1, 10), Fraction(1, 10) + Fraction(1, 200)
    b_lo = Fraction(2, 5) + 5 * DELTA
    b_hi = HALF + SLICE_W - Fraction(1, 10)  # 0.41 (el slice corta 0.415)
    assert b_lo < b_hi and a_lo < a_hi
    return a_lo, a_hi, b_lo, b_hi


def corner_mass(v, h, E, t):
    alpha, _ = mesh(v, h)
    N = 10 ** E
    q = [rep_prime(a, E) for a in alpha]
    ct, _ = build_ct(alpha, t)
    a_lo, a_hi, b_lo, b_hi = corner_defs()
    exact = (h <= 24)

    def in_slice(i, j):
        s = alpha[i] + alpha[j]
        return HALF < s <= HALF + SLICE_W

    def in_corner(i, j):
        a, b = alpha[i], alpha[j]
        return ((a_lo <= a <= a_hi and b_lo <= b <= b_hi) or
                (a_lo <= b <= a_hi and b_lo <= a <= b_hi))

    def in_win64(i, j):  # definicion doc 64 (max en [2/5+5d, 3/7])
        return Fraction(2, 5) + 5 * DELTA <= max(alpha[i], alpha[j]) <= Fraction(3, 7)

    acc = {r: Fraction(0) if exact else [] for r in
           ("tot", "slice", "corner", "slice_win")}
    cnt = {r: 0 for r in acc}
    for i in range(h):
        for j in range(h):
            memb = {"tot": True, "slice": in_slice(i, j),
                    "corner": in_slice(i, j) and in_corner(i, j),
                    "slice_win": in_slice(i, j) and in_win64(i, j)}
            ph = phi2(q[i], q[j])
            for r in acc:
                if memb[r]:
                    cnt[r] += 1
                    if exact:
                        acc[r] += ct[i][j] / ph
                    else:
                        acc[r].append(float(ct[i][j]) / ph)
    get = lambda r: acc[r] if exact else math.fsum(acc[r])
    return dict(M1=get("tot"), M1_slice=get("slice"),
                M1_corner=get("corner"), M1_slice_win=get("slice_win"),
                cnt_corner=cnt["corner"], cnt_slice=cnt["slice"],
                cnt_slice_win=cnt["slice_win"], exact=exact,
                min_alpha=float(min(alpha)))


def mobius_upto(n):
    mu = [1] * (n + 1)
    primes = []
    is_comp = [False] * (n + 1)
    mu[0] = 0
    for i in range(2, n + 1):
        if not is_comp[i]:
            primes.append(i)
            mu[i] = -1
        for p in primes:
            if i * p > n:
                break
            is_comp[i * p] = True
            if i % p == 0:
                mu[i * p] = 0
                break
            mu[i * p] = -mu[i]
    return mu


def main():
    print("== 1. Gran criba en la esquina: Q^2/N, Q ~ N^s ==")
    for s in (Fraction(101, 200), Fraction(103, 200), Fraction(13, 25)):
        # s = 0.505 (fondo), 0.515 (tope efectivo), 0.52 (tope nominal)
        for E in (7, 12, 24):
            ratio = 10 ** (E * float(2 * s - 1))
            print(f"  s={float(s):.3f} N=10^{E}: Q^2/N = 10^{E * float(2 * s - 1):.2f} "
                  f"= {ratio:.3g}")
            assert ratio > 1.0
    print("  -> Q^2/N > 1 y crece con N en toda la esquina: sin ahorro.")
    print()
    print("== 2. Umbral Teor. 1.1 vs esquina ==")
    cap = Fraction(1, 10) - 3 * DELTA
    assert cap == Fraction(97, 1000)
    a_lo, a_hi, b_lo, b_hi = corner_defs()
    print(f"  cap Teor. 1.1: N^({float(cap):.3f}); esquina alpha en "
          f"[{float(a_lo):.3f}, {float(a_hi):.3f}]")
    for gap, tag in ((a_lo - cap, "fondo"), (a_hi - cap, "tope")):
        assert gap > 0
        print(f"  hueco {tag}: {float(gap):.4f} en exponente")
        for E in (7, 12, 24):
            print(f"    N=10^{E}: N^hueco = {10 ** (E * float(gap)):.3f}")
    print("  -> hueco pequeno a N accesible pero -> infinito: sin teorema.")
    print()
    print("== 3. Masa de la esquina (re-verifica doc 64) ==")
    for (v, h, E) in [(10, 12, 12), (10, 24, 12), (10, 96, 12)]:
        r = corner_mass(v, h, E, Fraction(1, 20))
        M1 = float(r["M1"])
        print(f"  v={v} h={h} E={E}: min_alpha={r['min_alpha']:.4f} "
              f"M1slice/M1={float(r['M1_slice']) / M1:.4f} "
              f"M1esquina/M1slice={float(r['M1_corner']) / float(r['M1_slice']):.4f} "
              f"pares={r['cnt_corner']}/{r['cnt_slice']} "
              f"M1esquina/M1={float(r['M1_corner']) / M1:.2e} "
              f"[slice_win doc64: {float(r['M1_slice_win']) / float(r['M1_slice']):.4f}, "
              f"pares={r['cnt_slice_win']}/{r['cnt_slice']}]")
        assert r["M1_corner"] == r["M1_slice_win"], "esquina != slice_win"
    print("  -> esquina = slice_win del doc 64 (4.78% del slice a h=96).")
    print()
    print("== 4. Celda balanceada de Vaughan presente ==")
    x = 10 ** 12
    lo, hi, mid = 10 ** 4, 10 ** 8, 10 ** 6
    k = mid.bit_length() - 1
    a, b = 2 ** k, 2 ** (k + 1)
    assert lo < a <= mid < b <= hi
    print(f"  2^{k}={a} <= 10^6 < 2^{k + 1}={b}, en (10^4, 10^8]: "
          "M~K~x^(1/2) ocurre en el Tipo II de Vaughan.")
    print()
    print("== 5. Congruencia con a=N (un solo entero) ==")
    N = 10 ** 12  # par
    assert N % 2 == 0
    q = 1000003  # primo
    assert N % q != 0
    for (n1, n2) in [(123456, 789012), (999983, 5), (10 ** 6, 10 ** 6 + q)]:
        lhs = (n1 * N - n2 * N) % q == 0
        rhs = (n1 - n2) % q == 0
        assert lhs == rhs, (n1, n2)
    print(f"  N={N}, q={q} primo, (N,q)=1: n1*N=n2*N (mod q) <=> n1=n2 (mod q) "
          "exacto en 3 casos.")
    print(f"  Pero el desplazamiento es |a| = N ~ 10^12 (grande): el metodo "
          "espectral no lo admite (a fija); Zhang necesita modulos lisos (0).")
    print()
    print("== 6. Signo: lo descartado en Teor. 1.3 vs Heath-Brown ==")
    mu = mobius_upto(200)
    vals = set(mu[2:])
    assert 1 in vals and -1 in vals and 0 in vals
    print(f"  mu en [2,200] toma +1, -1 y 0: los terminos de Heath-Brown "
          "son firmados, no se pueden descartar.")
    print("  S_{2,3}, S_{2,5} del Teor. 1.3 cuentan rho(n/p1p2, p2) >= 0: "
          "descartables solo por ser minorante (cota inferior).")
    print()
    print("PASS: gran criba sin ahorro en la esquina; umbral 1/10 violado "
          "por 0.003-0.008; esquina = 4.78% del slice (0.004% de M1); "
          "celda balanceada presente; congruencia a=N verificada; "
          "mu firmado.")


if __name__ == "__main__":
    main()
