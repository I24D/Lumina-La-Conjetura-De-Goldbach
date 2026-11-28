#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Verificación de los puntos [P-rev] del doc 97 §8 y del doc 100 §7.

Procedencia: Claude Code (Anthropic, modelo claude-opus-5-5), 2026-11-28.
Método: reimplementación PROPIA a partir de las definiciones escritas en los docs
94, 96, 98 y 99. No importa ni ejecuta código de Lumina Chat ni de Muse. Los valores
publicados por ellos aparecen solo como "afirmación a contrastar".
Aritmética exacta (fractions.Fraction) en las partes A y B; coma flotante (math.fsum)
en las partes C y D, que usan psi reales.
Estado: PROPUESTA (pendiente de contraste por Muse). Control de reproducibilidad [M]
y comprobaciones de identidades finitas.
Fuentes: docs 94, 96, 98, 99 (OneDrive) y 95, 97, 100 (repositorio público, clon
anónimo del 2026-11-28, último commit 3a82328).
No es revisión independiente en el sentido del contrato (README, regla 4): este
agente comparte usuario de Windows con las corridas de Lumina en esta PC.
Solo biblioteca estándar. Uso:  python -B verificacion_97_100_claude_code.py
Nada de esto demuestra Goldbach, que sigue abierta.
"""
import cmath
import math
import sys
import time
from fractions import Fraction

T0 = time.time()
CHECKS = []


def check(name, ok, detail=""):
    CHECKS.append((name, bool(ok)))
    print(f"[{'PASS' if ok else 'FAIL'}] {name}" + (f" | {detail}" if detail else ""))


# ----------------------------------------------------------------------------
# Aritmética básica
# ----------------------------------------------------------------------------
def primes_upto(n):
    if n < 2:
        return []
    s = bytearray([1]) * (n + 1)
    s[0] = s[1] = 0
    for p in range(2, math.isqrt(n) + 1):
        if s[p]:
            s[p * p::p] = bytearray(len(range(p * p, n + 1, p)))
    return [i for i in range(n + 1) if s[i]]


def factor(n):
    f = {}
    p = 2
    while p * p <= n:
        while n % p == 0:
            f[p] = f.get(p, 0) + 1
            n //= p
        p += 1
    if n > 1:
        f[n] = f.get(n, 0) + 1
    return f


def phi(n):
    out = n
    for p in factor(n):
        out = out // p * (p - 1)
    return out


def mobius(n):
    f = factor(n)
    if any(a > 1 for a in f.values()):
        return 0
    return -1 if len(f) % 2 else 1


def divisors(n):
    ds = [1]
    for p, a in factor(n).items():
        ds = [d * p ** k for d in ds for k in range(a + 1)]
    return sorted(ds)


def iroot(x, k):
    """Mayor entero z con z**k <= x."""
    z = int(round(x ** (1.0 / k)))
    while z ** k > x:
        z -= 1
    while (z + 1) ** k <= x:
        z += 1
    return z


def is_prime(n):
    return n >= 2 and factor(n) == {n: 1}


def vplus_set(small, L):
    """Conjunto V^+(L) del doc 94 §0 (beta=2), con la restricción d < L:
    d libre de cuadrados con primos en `small`; si d = p1...pj, p1 > ... > pj,
    entonces p_l^2 * (p1...p_l) < L para todo l impar.
    Enumeración recursiva por primos decrecientes con poda (algoritmo distinto
    de la enumeración por máscaras de bits de los scripts 94/94d/99).
    Devuelve lista de (d, tupla_de_primos)."""
    ps = sorted(small, reverse=True)
    out = []

    def rec(start, prod, ell, chosen):
        out.append((prod, tuple(chosen)))
        for i in range(start, len(ps)):
            p = ps[i]
            nprod = prod * p
            if nprod >= L:
                continue
            l = ell + 1
            if l % 2 == 1 and p * p * nprod >= L:
                continue
            chosen.append(p)
            rec(i + 1, nprod, l, chosen)
            chosen.pop()

    if L > 1:
        rec(0, 1, 0, [])
    return out


def subsets(fs):
    out = [()]
    for p in fs:
        out += [t + (p,) for t in out]
    return out


# ----------------------------------------------------------------------------
# PARTE A — 94d: diagonalización de Q_s, rango Q_s log z, Q_s/P, cortes en r
# ----------------------------------------------------------------------------
# Afirmación a contrastar: salida guardada de Lumina (20261127_94d_diagonal_profile.txt)
LUMINA_94D = {
    (40, 0): (2.63906, 1.70054, 0.0000, 0.0000), (40, 1): (2.14423, 1.38169, 0.0000, 0.0000),
    (40, 2): (1.60213, 1.03237, 0.0000, 0.1396), (40, 3): (1.55603, 1.00266, 0.0000, 0.0165),
    (45, 0): (2.43403, 1.54896, 0.0000, 0.0000), (45, 1): (2.09597, 1.33383, 0.0596, 0.0596),
    (45, 2): (1.72964, 1.10070, 0.0000, 0.1401), (45, 3): (1.57595, 1.00290, 0.0011, 0.0291),
    (50, 0): (2.70742, 1.61921, 0.0000, 0.0000), (50, 1): (2.33139, 1.39432, 0.0596, 0.0596),
    (50, 2): (1.77395, 1.06094, 0.0036, 0.1120), (50, 3): (1.67517, 1.00186, 0.0005, 0.0243),
    (55, 0): (2.58093, 1.53355, 0.0596, 0.0596), (55, 1): (2.37318, 1.41010, 0.0648, 0.0648),
    (55, 2): (1.81175, 1.07651, 0.0136, 0.1468), (55, 3): (1.69068, 1.00457, 0.0010, 0.0308),
    (60, 0): (2.82873, 1.67584, 0.0596, 0.0596), (60, 1): (2.45448, 1.45412, 0.0984, 0.0984),
    (60, 2): (1.82346, 1.08028, 0.0267, 0.1815), (60, 3): (1.69998, 1.00713, 0.0021, 0.0374),
}


def part_A():
    print("\n=== PARTE A — doc 96 §3 / script 94d: Q_s exacto por la diagonal y por Gram ===")
    U0 = (5 / 3, 2.0, 3.0, 4.0)
    res = []
    for e in (40, 45, 50, 55, 60):
        N = 2 ** e
        A = 3
        X = N // A
        D = math.isqrt(X)
        z = iroot(X, 10)
        y = iroot(X, 3)
        same_float = (z == int(X ** 0.1)) and (y == int(X ** (1 / 3)))
        pr = primes_upto(y)
        cand = [p for p in pr if z <= p < y and math.gcd(p, A * N) == 1]
        small = [p for p in pr if p < z and math.gcd(p, A * N) == 1]
        excluded = sorted(set(factor(A * N)))
        for iu, u0 in enumerate(U0):
            target = D / (z ** u0)
            s = min(cand, key=lambda p: abs(p - target))
            L = D // s
            V = vplus_set(small, L)
            # C_s(r) exacto: suma de mu(d)/phi(d) sobre d en V+ con r | d
            C = {}
            phi_of = {}
            for d, fs in V:
                w = Fraction(-1 if len(fs) % 2 else 1, math.prod(p - 1 for p in fs))
                for rs in subsets(fs):
                    r = math.prod(rs)
                    C[r] = C.get(r, Fraction(0)) + w
                    phi_of[r] = math.prod(p - 1 for p in rs)
            q = {r: phi_of[r] * c * c for r, c in C.items()}
            Q = sum(q.values(), Fraction(0))
            P = Fraction(1)
            for p in small:
                P *= Fraction(p * p - 3 * p + 3, (p - 1) ** 2)
            # Ruta 2 (Gram, coma flotante): sum_{d,e} mu mu gcd(d,e)/(phi phi)
            ds = [(d, (-1 if len(fs) % 2 else 1) / math.prod(p - 1 for p in fs)) for d, fs in V]
            Qg = math.fsum(wd * we * math.gcd(d, e2) for d, wd in ds for e2, we in ds)
            share_Lz = float(sum((v for r, v in q.items() if r * z > L), Fraction(0)) / Q)
            share_sqrtL = float(sum((v for r, v in q.items() if r * r > L), Fraction(0)) / Q)
            share_L = float(sum((v for r, v in q.items() if r > L), Fraction(0)) / Q)
            ddiag = math.fsum(1.0 / math.prod(p - 1 for p in fs) for _, fs in V)
            lz = math.log(z)
            row = dict(e=e, iu=iu, u0=u0, z=z, s=s, L=L, nV=len(V), Qlog=float(Q) * lz,
                       QP=float(Q / P), Plog=float(P) * lz, share_Lz=share_Lz,
                       share_sqrtL=share_sqrtL, share_L=share_L, gram_err=abs(Qg - float(Q)),
                       ddiag_Q=ddiag / float(Q), excluded=excluded, same_float=same_float,
                       u=math.log(L) / lz)
            res.append(row)
            lum = LUMINA_94D[(e, iu)]
            print(f"e={e} u0={u0:.3f} u={row['u']:.3f} z={z} s={s} L={L} |V+|={len(V)} "
                  f"Qlogz={row['Qlog']:.5f} (Lumina {lum[0]:.5f}) Q/P={row['QP']:.5f} (Lumina {lum[1]:.5f}) "
                  f"P*logz={row['Plog']:.5f} r>L/z={share_Lz:.4f} r>sqrt(L)={share_sqrtL:.4f} "
                  f"r>L={share_L:.4f} |Gram-diag|={row['gram_err']:.1e} excl={excluded}")
    # Comprobaciones contra la afirmación de Lumina
    md = max(max(abs(r['Qlog'] - LUMINA_94D[(r['e'], r['iu'])][0]),
                 abs(r['QP'] - LUMINA_94D[(r['e'], r['iu'])][1])) for r in res)
    check("A1 94d: Qlogz y Q/P de los 20 escenarios coinciden con la salida de Lumina", md < 6e-6,
          f"max desviación {md:.2e} (redondeo a 5 decimales)")
    ms = max(max(abs(r['share_Lz'] - LUMINA_94D[(r['e'], r['iu'])][2]),
                 abs(r['share_sqrtL'] - LUMINA_94D[(r['e'], r['iu'])][3])) for r in res)
    check("A2 94d: porciones r>L_s/z y r>sqrt(L_s) coinciden escenario a escenario", ms < 6e-5,
          f"max desviación {ms:.2e}")
    ge = max(r['gram_err'] for r in res)
    check("A3 Identidad diagonal Q_s = sum phi(r)C_s(r)^2 = forma de Gram (doble suma)", ge < 1e-12,
          f"max |diag exacta - Gram float| = {ge:.1e}")
    check("A4 z, y por raíz entera coinciden con int(X**0.1), int(X**(1/3))", all(r['same_float'] for r in res))
    qlog = [r['Qlog'] for r in res]
    qp = [(r['QP'], r['e'], r['u0']) for r in res]
    mn = min(qp)
    mx = max(qp)
    print(f"RANGO Qlogz = [{min(qlog):.5f}, {max(qlog):.5f}]")
    print(f"RANGO Q/P = [{mn[0]:.5f} (e={mn[1]}, u0={mn[2]:.3f}), {mx[0]:.5f} (e={mx[1]}, u0={mx[2]:.3f})]")
    check("A5 Mínimo real de Q/P en 94d es 1.00186 (e=50,u0=4), NO 1.00266 como dice doc 96 §3",
          abs(mn[0] - 1.00186) < 6e-6 and mn[1] == 50)
    pl = [r['Plog'] for r in res]
    print(f"RANGO P*logz (parte del producto local completo) = [{min(pl):.5f}, {max(pl):.5f}]")
    excl = {tuple(r['excluded']) for r in res}
    check("A6 Factor local fijo en 94d: los primos excluidos son {2,3} en los 20 escenarios", excl == {(2, 3)},
          f"conjuntos excluidos: {excl}")
    # Por u: el rango de Qlogz viene de Q/P (truncación V+), no del factor local
    for iu, u0 in enumerate(U0):
        sub = [r for r in res if r['iu'] == iu]
        print(f"  u0={u0:.3f}: Q/P en [{min(r['QP'] for r in sub):.4f}, {max(r['QP'] for r in sub):.4f}], "
              f"Qlogz en [{min(r['Qlog'] for r in sub):.4f}, {max(r['Qlog'] for r in sub):.4f}]")
    top = max(res, key=lambda r: r['Qlog'])
    check("A7 El extremo superior 2.82873 aparece en u0=5/3 (máxima truncación V+), con N=2^60",
          abs(top['u0'] - 5 / 3) < 1e-9 and top['e'] == 60, f"Q/P={top['QP']:.5f}, P*logz={top['Plog']:.5f}")
    check("A8 La porción r>L_s es exactamente 0 en los 20 escenarios (r|d<L_s obliga r<L_s)",
          all(r['share_L'] == 0.0 for r in res))
    mx_Lz = max(res, key=lambda r: r['share_Lz'])
    mx_sq = max(res, key=lambda r: r['share_sqrtL'])
    print(f"MAX r>L/z = {mx_Lz['share_Lz']:.4%} (e={mx_Lz['e']}, u0={mx_Lz['u0']:.3f}); "
          f"MAX r>sqrt(L) = {mx_sq['share_sqrtL']:.4%} (e={mx_sq['e']}, u0={mx_sq['u0']:.3f})")
    dd = [r['ddiag_Q'] for r in res]
    print(f"RANGO D_diag/Q_s = [{min(dd):.3f}, {max(dd):.3f}]")


# ----------------------------------------------------------------------------
# PARTE B — 94e: regla exacta de atribución de masa por conductor completo ah
# ----------------------------------------------------------------------------
LUMINA_94E = {1.0: 49.169430, 2.0: 73.278523, 2.5: 87.245823, 2.7: 92.878861,
              2.8: 95.249831, 2.9: 97.320588, 3.0: 99.111943}


def part_B():
    print("\n=== PARTE B — doc 96 §4 / doc 98 §5 / script 94e: masa Gram por corte ah <= R0 ===")
    N = 10 ** 12
    A = 3
    X = N // A
    D = math.isqrt(X)
    z = iroot(X, 10)
    y = iroot(X, 3)
    pr = primes_upto(y)
    small = [p for p in pr if p < z and math.gcd(p, A * N) == 1]
    sw = [p for p in pr if z <= p < y and math.gcd(p, A * N) == 1]
    phiA = phi(A)
    K = {}
    phi_h = {}
    nrows = 0
    for s in sw:
        L = D // s
        for d, fs in vplus_set(small, L):
            nrows += 1
            fac = (s,) + fs
            b = -1 if len(fs) % 2 else 1
            w = Fraction(b, phiA * (s - 1) * math.prod(p - 1 for p in fs))
            for hs in subsets(fac):
                h = math.prod(hs)
                K[h] = K.get(h, Fraction(0)) + w
                phi_h[h] = math.prod(p - 1 for p in hs)
    G = sum((phi_h[h] * v * v for h, v in K.items()), Fraction(0))
    total = A * G
    print(f"PARAMS N=10^12 A=3 X={X} D={D} z={z} y={y} small={small} filas(s,d)={nrows} |H|={len(K)}")
    print(f"G_A = {float(G):.12f}  A*G_A = {float(total):.12f}")
    check("B1 Soporte y G_A reproducen doc 94 §5 (1444 filas, 1451 h, G_A=0.450566952208)",
          nrows == 1444 and len(K) == 1451 and abs(float(G) - 0.450566952208) < 1e-12)
    adivs = divisors(A)
    worst = 0.0
    worst_mix = Fraction(0)
    for kappa in (1.0, 2.0, 2.5, 2.7, 2.8, 2.9, 3.0):
        R0 = math.log(N) ** kappa
        exact = sum((phi(a) * phi_h[h] * v * v for h, v in K.items() for a in adivs
                     if math.gcd(a, h) == 1 and a * h <= R0), Fraction(0))
        hcut = sum((A * phi_h[h] * v * v for h, v in K.items() if h <= R0), Fraction(0))
        ahcut = sum((A * phi_h[h] * v * v for h, v in K.items() if A * h <= R0), Fraction(0))
        mix = (phi(1) * hcut + phi(A) * ahcut) / A
        worst_mix = max(worst_mix, abs(mix - exact))
        sh = float(exact / total)
        worst = max(worst, abs(100 * sh - LUMINA_94E[kappa]))
        hi_exact = 1 - sh
        hi_h = 1 - float(hcut / total)
        hi_ah = 1 - float(ahcut / total)
        amp = lambda f: (1 / math.sqrt(f)) if f > 0 else float('inf')
        print(f"kappa={kappa:.1f} R0={R0:.1f}  regla exacta 94e: bajo={100*sh:.6f}% alto={100*hi_exact:.4f}% "
              f"(energía x{1/hi_exact:.2f}, amplitud x{amp(hi_exact):.2f}) | lectura 'h<=R0': alto={100*hi_h:.4f}% "
              f"(amplitud x{amp(hi_h):.2f}) | lectura 'A*h<=R0': alto={100*hi_ah:.4f}% (amplitud x{amp(hi_ah):.2f})")
    check("B2 Tabla kappa del 94e reproducida con la regla m(a,h)=phi(a)phi(h)K(h)^2, corte ah<=R0", worst < 6e-7,
          f"max desviación {worst:.1e} puntos porcentuales")
    check("B3 La regla exacta es la mezcla (phi(1)*[corte h] + phi(3)*[corte 3h])/3, identidad exacta",
          worst_mix == 0)
    # K(21) = K(7): el coeficiente depende solo de la parte nueva h
    # (en 94 se evalúa K en conductores: K(ah) := K(h))
    check("B4 K_A(h) coincide con los ejemplos del doc 94 (K(7), K(11), K(13), K(17))",
          abs(float(K[7]) + 0.067848416650) < 1e-12 and abs(float(K[11]) + 0.028330369695) < 1e-12
          and abs(float(K[13]) + 0.020169817129) < 1e-12 and abs(float(K[17]) - 0.021527777778) < 1e-12)


# ----------------------------------------------------------------------------
# Psi real y caracteres (partes C y D)
# ----------------------------------------------------------------------------
def lambda_list(Nmax):
    out = []
    for p in primes_upto(Nmax):
        lp = math.log(p)
        x = p
        while x <= Nmax:
            out.append((x, lp))
            x *= p
    return out


class PsiData:
    def __init__(self, lam, N):
        self.N = N
        self.lam = [(n, w) for n, w in lam if n <= N]
        self.psi = math.fsum(w for _, w in self.lam)
        self._cnt = {}

    def counts(self, r):
        """Suma de Lambda(n), n<=N, por clase n mod r."""
        if r not in self._cnt:
            c = [0.0] * r
            for n, w in self.lam:
                c[n % r] += w
            self._cnt[r] = c
        return self._cnt[r]


def S_charfree(P, r):
    """Ruta 1, SIN caracteres. Para (n,r)=1 y (N,r)=1:
    sum*_{chi mod r} conj(chi(N)) chi(n) = sum_{d|r} mu(r/d) phi(d) [n = N mod d]
    (inversión de Möbius sobre conductores). Como chi primitivo mod r se anula
    exactamente en (n,r)>1, S_r = sum_{(n,r)=1} Lambda(n) * esa suma."""
    N = P.N
    cnt = P.counts(r)
    ds = [(d, mobius(r // d) * phi(d)) for d in divisors(r) if mobius(r // d) != 0]
    tot = []
    for a in range(r):
        if cnt[a] == 0.0 or math.gcd(a, r) != 1:
            continue
        g = sum(c for d, c in ds if (a - N) % d == 0)
        if g:
            tot.append(cnt[a] * g)
    return math.fsum(tot)


def M_mobius(P, r):
    """M_r = sum_{d|r} mu(r/d) phi(d) E(N;d,N), psi(N;d,N) sobre TODOS los n=N mod d."""
    N = P.N
    cnt = P.counts(r)
    terms = []
    for d in divisors(r):
        m = mobius(r // d)
        if m == 0:
            continue
        psi_d = math.fsum(cnt[a] for a in range(r) if (a - N) % d == 0)
        terms.append(m * phi(d) * (psi_d - P.psi / phi(d)))
    return math.fsum(terms)


def J_of(p, N):
    j = 0
    x = p
    while x <= N:
        j += 1
        x *= p
    return j


def delta_charfree(N, r):
    """Delta_r por la forma reducida (98.3) del doc 98, evaluada sin caracteres:
    1_{r primo} Theta_r(N) + sum_{p||r, r/p>1} sum_{j: p^j<=N} log p * sum*_{chi mod r/p} conj(chi(N)) chi(p^j),
    con sum*_{chi mod e} conj(chi(N)) chi(m) = sum_{d|e} mu(e/d) phi(d)[m = N mod d], válida porque (p,e)=1."""
    f = factor(r)
    out = []
    if len(f) == 1 and list(f.values())[0] == 1:
        out.append(J_of(r, N) * math.log(r))
    for p, a in f.items():
        e = r // p
        if a != 1 or e == 1:
            continue
        ds = [(d, mobius(e // d) * phi(d)) for d in divisors(e) if mobius(e // d) != 0]
        lp = math.log(p)
        x = p
        while x <= N:
            g = sum(c for d, c in ds if (x - N) % d == 0)
            out.append(lp * g)
            x *= p
    return math.fsum(out)


def primitive_root_pp(p, a):
    m = p ** a
    ph = m // p * (p - 1)
    qs = list(factor(ph))
    for g in range(2, m):
        if math.gcd(g, m) == 1 and all(pow(g, ph // q, m) != 1 for q in qs):
            return g
    raise RuntimeError


class CharGroupOdd:
    """Grupo de caracteres de (Z/rZ)*, r impar, construcción propia por CRT.
    La PRIMITIVIDAD se decide con una prueba directa e independiente de etiquetas:
    chi es primitivo si y solo si, para cada primo p|r, chi no es trivial sobre
    {n : n = 1 mod r/p, (n,r)=1}. El test es exacto (aritmética de Fraction)."""

    def __init__(self, r):
        assert r % 2 == 1 and r > 1
        self.r = r
        self.comps = []
        for p, a in sorted(factor(r).items()):
            m = p ** a
            ph = m // p * (p - 1)
            g = primitive_root_pp(p, a)
            ind = {}
            x = 1
            for t in range(ph):
                ind[x] = t
                x = x * g % m
            self.comps.append((m, ph, ind))
        self.units = [n for n in range(r) if math.gcd(n, r) == 1]
        self.idx = {n: tuple(c[2][n % c[0]] for c in self.comps) for n in self.units}
        self.phs = [c[1] for c in self.comps]
        labels = [()]
        for ph in self.phs:
            labels = [l + (k,) for l in labels for k in range(ph)]
        kers = []
        for p in factor(r):
            q = r // p
            kers.append([n for n in self.units if (n - 1) % q == 0])
        self.primitive = []
        for k in labels:
            if all(any(self._frac(k, n) != 0 for n in ker) for ker in kers):
                self.primitive.append(k)

    def _frac(self, k, n):
        t = self.idx[n]
        return sum((Fraction(ki * ti, ph) for ki, ti, ph in zip(k, t, self.phs)), Fraction(0)) % 1

    def value(self, k, n):
        n %= self.r
        if math.gcd(n, self.r) != 1:
            return 0j
        t = self.idx[n]
        ang = 2 * math.pi * math.fsum(ki * ti / ph for ki, ti, ph in zip(k, t, self.phs))
        return cmath.exp(1j * ang)


def nprim_formula(r):
    """Número de caracteres primitivos mod r (fórmula multiplicativa conocida)."""
    out = 1
    for p, a in factor(r).items():
        if a == 1:
            out *= p - 2
        else:
            out *= p ** a * (1 - Fraction(1, p)) ** 2
    return int(out)


def S_chars(P, G):
    """Ruta 2: suma directa sobre caracteres primitivos de mi construcción."""
    cnt = P.counts(G.r)
    tot = 0j
    for k in G.primitive:
        psi_chi = sum(cnt[a] * G.value(k, a) for a in G.units)
        tot += G.value(k, P.N).conjugate() * psi_chi
    return tot


def delta_chars(N, r, groups):
    """Delta_r por (98.3) con caracteres explícitos de mi construcción (r impar)."""
    f = factor(r)
    tot = 0j
    if len(f) == 1 and list(f.values())[0] == 1:
        tot += J_of(r, N) * math.log(r)
    for p, a in f.items():
        e = r // p
        if a != 1 or e == 1:
            continue
        G = groups.setdefault(e, CharGroupOdd(e))
        lp = math.log(p)
        J = J_of(p, N)
        for k in G.primitive:
            cp = G.value(k, p)
            pp = sum(cp ** j for j in range(1, J + 1)) * lp
            tot += G.value(k, N).conjugate() * pp
    return tot


# ----------------------------------------------------------------------------
# PARTE C — Delta_r (doc 98 §2-3) con psi reales, N = 999983
# ----------------------------------------------------------------------------
LUMINA16 = (3, 5, 7, 11, 15, 21, 35, 55, 105, 385, 9, 25, 45, 75, 225, 315)
DOC100_16 = (3, 5, 7, 9, 11, 13, 15, 21, 25, 27, 35, 45, 49, 63, 75, 99)
DOC97_10 = (3, 5, 6, 7, 9, 15, 21, 25, 27, 35)
# Afirmación a contrastar: salida guardada del script 97 de Lumina (|Delta|/|S| en %)
LUMINA97_RATIO = {3: 1.843815, 5: 1.229959, 7: 2.588591, 11: 0.472407, 15: 0.0, 21: 0.884261,
                  35: 0.351828, 55: 0.531024, 105: 0.722364, 385: 1.192253, 9: 0.0, 25: 0.0,
                  45: 0.050155, 75: 0.104442, 225: 0.0, 315: 0.632908}
# Afirmación a contrastar: doc 97 §3.2 (S_r, RHS sin Delta)
DOC97_TABLE = {3: (-715.0038, -728.1872), 5: (1046.8240, 1033.9485), 6: (0.0, -0.6931),
               7: (526.2080, 512.5867), 9: (674.8847, 674.8847), 15: (1563.0200, 1563.0200),
               21: (1540.4249, 1554.0463), 25: (357.5449, 357.5449), 27: (791.9846, 791.9846),
               35: (2382.8836, 2374.4999)}
# Afirmación a contrastar: doc 100 §2.2, columna "espurio"
DOC100_SPUR = {9: -13.1833, 25: -12.8755, 45: 13.1833, 49: -13.6214, 63: 13.1833, 75: 12.8755,
               99: 13.1833, 121: -11.9895}


def part_C(lam):
    print("\n=== PARTE C — doc 98 §2-3 / script 97: Delta_r con psi reales, N=999983 ===")
    N = 999983
    P = PsiData(lam, N)
    print(f"psi(N) = {P.psi:.4f}  (doc 97: 999586.5975)")
    check("C0 psi(999983) coincide con doc 97", abs(P.psi - 999586.5975) < 6e-5)
    groups = {}
    allr = sorted(set(LUMINA16) | set(DOC100_16) | set(DOC97_10))
    rows = {}
    for r in allr:
        S1 = S_charfree(P, r)
        M = M_mobius(P, r)
        Dcf = delta_charfree(N, r)
        row = dict(S1=S1, M=M, Dcf=Dcf, Dpure=S1 - M)
        if r % 2 == 1:
            G = groups.setdefault(r, CharGroupOdd(r))
            S2 = S_chars(P, G)
            Dch = delta_chars(N, r, groups)
            row.update(S2=S2, Dch=Dch, nprim=len(G.primitive), nprim_f=nprim_formula(r))
        rows[r] = row
    print(" r  #chi*  S_r(ruta1 sin caract.)  |S1-S2(caract.)|  M_r(Möbius sin Delta)  Delta(98.3)  |S-(M+Delta)|  |Delta|/|S|")
    for r in allr:
        x = rows[r]
        s2e = abs(x['S1'] - x['S2']) if 'S2' in x else float('nan')
        npr = x.get('nprim', 0)
        ratio = abs(x['Dcf']) / abs(x['S1']) if abs(x['S1']) > 1e-6 else float('nan')
        print(f"{r:4d} {npr:5d} {x['S1']:20.6f} {s2e:16.2e} {x['M']:20.6f} {x['Dcf']:12.6f} "
              f"{abs(x['S1'] - (x['M'] + x['Dcf'])):14.2e} {ratio:10.6%}")
    odd = [r for r in allr if r % 2 == 1]
    e12 = max(abs(rows[r]['S1'] - rows[r]['S2']) / max(1, abs(rows[r]['S1'])) for r in odd)
    check("C1 S_r: ruta sin caracteres = ruta con caracteres primitivos propios", e12 < 1e-10,
          f"max error relativo {e12:.1e} en {len(odd)} conductores impares")
    check("C2 Número de primitivos (test directo) = fórmula multiplicativa",
          all(rows[r]['nprim'] == rows[r]['nprim_f'] for r in odd))
    im = max(abs(rows[r]['S2'].imag) for r in odd)
    check("C3 S_r es real (parte imaginaria de la ruta por caracteres)", im < 1e-7, f"max |Im| {im:.1e}")
    eid = max(abs(rows[r]['S1'] - (rows[r]['M'] + rows[r]['Dcf'])) / max(1, abs(rows[r]['S1'])) for r in allr)
    check("C4 Identidad corregida S_r = M_r + Delta_r (98.3), con S_r por ruta sin caracteres", eid < 1e-10,
          f"max error relativo {eid:.1e} en {len(allr)} conductores")
    edd = max(abs(rows[r]['Dch'] - rows[r]['Dcf']) for r in odd)
    check("C5 Delta_r (98.3): evaluación con caracteres = evaluación sin caracteres", edd < 1e-8,
          f"max |dif| {edd:.1e}")
    # Doc 98 §3 / script 97 de Lumina: los 16 conductores y sus |Delta|/|S|
    dev = max(abs(100 * abs(rows[r]['Dcf']) / abs(rows[r]['S1']) - LUMINA97_RATIO[r]) for r in LUMINA16)
    check("C6 |Delta|/|S| de los 16 conductores del script 97 de Lumina reproducidos", dev < 2e-6,
          f"max desviación {dev:.1e} puntos porcentuales")
    for name, st in (("Lumina (doc 98)", LUMINA16), ("doc 100", DOC100_16)):
        best = max(st, key=lambda r: abs(rows[r]['Dcf']) / abs(rows[r]['S1']))
        print(f"  conjunto {name}: max |Delta|/|S| = {abs(rows[best]['Dcf'])/abs(rows[best]['S1']):.6%} en r={best}")
    common = sorted(set(LUMINA16) & set(DOC100_16))
    print(f"  conductores comunes a ambos conjuntos ({len(common)}): {common}")
    print(f"  solo Lumina: {sorted(set(LUMINA16) - set(DOC100_16))}; solo doc 100: {sorted(set(DOC100_16) - set(LUMINA16))}")
    # Doc 97 §3.2: tabla de 10 conductores
    d97 = max(max(abs(rows[r]['S1'] - DOC97_TABLE[r][0]), abs(rows[r]['M'] - DOC97_TABLE[r][1])) for r in DOC97_10)
    check("C7 Tabla del doc 97 §3.2 (S_r y RHS sin Delta, 10 conductores, incluye r=6) reproducida",
          d97 < 6e-5, f"max desviación {d97:.1e}")
    check("C8 r=6: no hay primitivos mod 6, S_6=0 y Delta_6 = log 2",
          abs(rows[6]['S1']) < 1e-8 and abs(rows[6]['Dcf'] - math.log(2)) < 1e-9)
    # p || r: en la forma con caracteres, los términos con p^2 | r se anulan solos (chi*(p)=0 si p | r/p).
    # El "espurio" del doc 100 §2.2 equivale a sustituir chi*(p)^j por 1:
    sp_ok = True
    for r, claim in DOC100_SPUR.items():
        tot = 0.0
        for p, a in factor(r).items():
            if a >= 2:
                e = r // p
                s_star = sum(mobius(e // d) * phi(d) for d in divisors(e) if (N - 1) % d == 0)  # sum* conj(chi(N))
                tot += math.log(p) * J_of(p, N) * s_star
        sp_ok &= abs(tot - claim) < 6e-5
        print(f"  r={r}: 'espurio' doc 100 = log p * floor(log_p N) * sum*_(chi mod r/p) conj chi(N) = {tot:+.4f} (doc 100: {claim:+.4f}); "
              f"término real con chi*(p)=0: 0")
    check("C9 El 'espurio' del doc 100 se reproduce SOLO si se da a chi*(p) el valor 1 con p | r/p", sp_ok)
    # Contraste: la forma ortogonal del propio doc 100 (P_e = sum_{d|e} mu(e/d) A_d,
    # A_d = phi(d) sum_{j: p^j = N mod d} log p) aplicada SIN la condición p||r.
    differ = []
    for r, claim in DOC100_SPUR.items():
        tot = 0.0
        for p, a in factor(r).items():
            if a >= 2:
                e = r // p
                x = p
                while x <= N:
                    tot += math.log(p) * sum(mobius(e // d) * phi(d) for d in divisors(e) if (x - N) % d == 0)
                    x *= p
        if abs(tot - claim) > 6e-5:
            differ.append(r)
        print(f"  r={r}: forma ortogonal del doc 100 sin p||r = {tot:+.4f} (columna 'espurio' del doc 100: {claim:+.4f})")
    check("C10 La forma ortogonal del doc 100 sin p||r NO reproduce su columna 'espurio' en todos los r "
          "(fuera de su hipótesis (p,r/p)=1 da valores sin significado)", len(differ) > 0, f"difiere en r={differ}")


# ----------------------------------------------------------------------------
# PARTE D — segundo factor (doc 99): soporte, B_M, B_cross, B_Delta
# ----------------------------------------------------------------------------
DOC99_DISK = {1.0: (1.981992172e7, 1.513555593e4, 2.784367244e2, 7.777020e-4, 29),
              1.5: (1.073276864e7, -1.393005151e4, 9.808629168e1, -1.288760e-3, 16),
              2.0: (1.295385198e5, -3.566649987e1, 2.455059728e-3, -2.753161e-4, 1)}
CHAT_DOC100 = (3.9222e7, 5.4023e3, 466.7, 1.4966e-4, 69)   # transcripción del chat en doc 100 §1/§3


def doc99_model(lam, N, A=3):
    X = N // A
    D = math.isqrt(X)
    z = iroot(X, 10)
    y = iroot(X, 3)
    pr = primes_upto(y)
    small = [p for p in pr if p < z and math.gcd(p, A * N) == 1]
    sw = [p for p in pr if z <= p < y and math.gcd(p, A * N) == 1]
    H = set()
    nrows = 0
    for s in sw:
        for d, fs in vplus_set(small, D // s):
            nrows += 1
            for hs in subsets((s,) + fs):
                H.add(math.prod(hs))
    cond = sorted({a * h for a in divisors(A) for h in H if math.gcd(a, h) == 1 and a * h > 1})
    P = PsiData(lam, N)
    data = {}
    for r in cond:
        S = S_charfree(P, r)
        M = M_mobius(P, r)
        data[r] = (S, M, S - M, delta_charfree(N, r))
    return dict(X=X, D=D, z=z, y=y, small=small, nrows=nrows, H=len(H), cond=cond, data=data)


def B_parts(model, N, kappa):
    R0 = math.log(N) ** kappa
    hi = [r for r in model['cond'] if r > R0]
    BM = math.fsum(model['data'][r][1] ** 2 / phi(r) for r in hi)
    BC = math.fsum(2 * model['data'][r][1] * model['data'][r][2] / phi(r) for r in hi)
    BD = math.fsum(model['data'][r][2] ** 2 / phi(r) for r in hi)
    Bdir = math.fsum(model['data'][r][0] ** 2 / phi(r) for r in hi)
    return len(hi), BM, BC, BD, Bdir


def part_D(lam):
    print("\n=== PARTE D — doc 99 / script 99: segundo factor corregido ===")
    out = {}
    for N in (1000000, 999983):
        m = doc99_model(lam, N)
        out[N] = m
        print(f"N={N}: X={m['X']} D={m['D']} z={m['z']} y={m['y']} primos pequeños={m['small']} "
              f"filas(s,d)={m['nrows']} |H|={m['H']} conductores={len(m['cond'])}")
        eid = max(abs(S - (M + Dcf)) / max(1, abs(S)) for S, M, Dp, Dcf in m['data'].values())
        print(f"  identidad S=M+Delta(98.3) en el soporte: max error relativo {eid:.1e}")
        for kappa in (1.0, 1.5, 2.0):
            nh, BM, BC, BD, Bdir = B_parts(m, N, kappa)
            print(f"  kappa={kappa}: altos={nh} B_M={BM:.9e} B_cross={BC:.9e} B_Delta={BD:.9e} "
                  f"cambio rel={(BC + BD) / BM:.6e} |B_corr-B_directo|={abs(BM + BC + BD - Bdir):.1e}")
    m6 = out[1000000]
    dev = 0.0
    for kappa, (bm, bc, bd, rel, nh) in DOC99_DISK.items():
        nh2, BM, BC, BD, _ = B_parts(m6, 1000000, kappa)
        dev = max(dev, abs(BM / bm - 1), abs(BC / bc - 1), abs(BD / bd - 1))
        if nh2 != nh:
            dev = float('inf')
    check("D1 Tabla del doc 99 (N=10^6, 33 conductores) reproducida por ruta sin caracteres",
          len(m6['cond']) == 33 and dev < 1e-6, f"max desviación relativa {dev:.1e}")
    m9 = out[999983]
    nh, BM, BC, BD, _ = B_parts(m9, 999983, 1.0)
    ok = (len(m9['cond']) == 69 and abs(BM / CHAT_DOC100[0] - 1) < 2e-5 and abs(BC / CHAT_DOC100[1] - 1) < 2e-5
          and abs((BC + BD) / BM / CHAT_DOC100[3] - 1) < 5e-5)
    check("D2 Las cifras del chat que audita el doc 100 (69 conductores, B_M=3.9222e7, B_cross=5.4023e3, "
          "cambio 1.4966e-4) salen del MISMO modelo con N=999983", ok,
          f"B_M={BM:.6e} B_cross={BC:.6e} B_Delta={BD:.4f} cambio={(BC+BD)/BM:.6e} altos={nh}")
    check("D3 B_Delta real con N=999983 es 467.66, no 466.7 (doc 100 lo transcribe 466,7)",
          abs(BD - 467.664) < 0.01 and abs(BD - 466.7) > 0.5, f"B_cross/B_Delta = {BC/BD:.3f}")


def main():
    try:
        sys.stdout.reconfigure(encoding="utf-8")
    except (AttributeError, ValueError):
        pass
    print(f"Python {sys.version.split()[0]}")
    part_A()
    part_B()
    lam = lambda_list(1000000)
    part_C(lam)
    part_D(lam)
    npass = sum(ok for _, ok in CHECKS)
    print(f"\nRESUMEN: {npass}/{len(CHECKS)} PASS, {len(CHECKS) - npass} FAIL. Tiempo {time.time() - T0:.1f} s")
    return 0 if npass == len(CHECKS) else 1


if __name__ == "__main__":
    sys.exit(main())
