#!/usr/bin/env python3
# 101_experimento_e100_1_conteo_a100_1.py  (Doc 101)
# Experimento E100.1 + conteo A100.1 propuestos en el doc 100.
# Solo stdlib.
#
# E100.1: crecimiento de W_max(N) y S_BV(N) con N par creciente.
#   - Parte A (soporte fijo R2 = soporte del doc 94): S_BV(N), B_M(N),
#     W_max(N) medidos con psi REAL (criba segmentada, SegPsi) para
#     N par en tres familias: 2 x primo, 2^e, primoriales.
#   - Parte B (soporte escalado): familia de tendencia de los docs
#     94/95/100: D = z^5, y = z^(10/3), A = 3, con z = N^(1/10)
#     (marco SW-Q, docs 87/88: D = N^(1/2), y = N^(1/3)).
#     W_max y |H| EXACTOS por la formula de conteo de A100.1.
# A100.1: conteo semi-analitico del soporte:
#   |H|  = sum_s |V^+(D/s)| + |V^+(D/s_min)|,
#   W_1  = sum_{h en V^+(D/s_min)} h/phi(h)
#          + sum_s (s/(s-1)) * sum_{h en V^+(D/s)} h/phi(h),
#   y W_max = W_1 (suma sobre subconjunto, terminos positivos).
#   La formula usa: (i) V^+ cerrado por divisores [T, prueba en el doc
#   101 sec. 2.1 + verificacion numerica J1]; (ii) monotonía de V^+ en
#   el nivel; (iii) la descomposicion de H por el unico primo >= z.
#
# Chequeos:
#   J0  Validacion: Model97 reproduce R2 (|H|=880); psi(999983) del
#       doc 97; SegPsi vs PsiEngine en E_d (N=999982); bv_decomposition
#       reproduce W_max=1011.28 y S_BV=4.3593e9 del doc 100 (N=999983).
#   J1  V^+ cerrado por divisores y monotono en el nivel (R1/R2/R3 y
#       tablas escaladas z=13,17).
#   J2  Formula de conteo vs Model97 (R1/R2/R3): |H| exacto y W_1 =
#       W_max de la descomposicion BV (N=999983, filtro vacuo).
#   J3  Familia escalada: formula vs enumeracion explicita (z=13,17,19,23).
#   J4  Regresion de alpha: ln W_1 y ln |H| contra ln z (z=11..29);
#       alpha = pendiente/10 (N = z^10). Intervalos honestos (IC 95%,
#       leave-one-out, pendientes locales). Envolvente Nsf sin V^+.
#   J5  E100.1 parte A: tabla de mediciones (8 valores de N par),
#       regresiones de S_BV y B_M, delta(N) = S_BV/(#divs * N ln N),
#       comparacion con N^(3/2) ln N, colapso primorial de W_max.
#   J6  W_max = W_1 (argmax en d=1) en cada N de la parte A.
# Uso: python3 101_experimento_e100_1_conteo_a100_1.py
import math
import sys
import time
from array import array


# ------------------------------------------------------------------ utiles
# (copiados de los scripts de los docs 97/100, misma convencion)
def primes_upto(n):
    if n < 2:
        return []
    bs = bytearray(b"\x01") * (n + 1)
    bs[0:2] = b"\x00\x00"
    for i in range(2, int(n ** 0.5) + 1):
        if bs[i]:
            bs[i * i:n + 1:i] = b"\x00" * ((n - i * i) // i + 1)
    return [i for i in range(2, n + 1) if bs[i]]


def mu_spf(n):
    mu = [0] * (n + 1)
    spf = [0] * (n + 1)
    primes = []
    mu[1] = 1
    for i in range(2, n + 1):
        if spf[i] == 0:
            spf[i] = i
            primes.append(i)
            mu[i] = -1
        for p in primes:
            if i * p > n:
                break
            spf[i * p] = p
            if i % p == 0:
                mu[i * p] = 0
                break
            mu[i * p] = -mu[i]
    spf[1] = 1
    return mu, spf, primes


_phi_cache = {}


def phi_of(n, primes):
    if n in _phi_cache:
        return _phi_cache[n]
    r = n
    x = n
    for p in primes:
        if p * p > x:
            break
        if x % p == 0:
            r -= r // p
            while x % p == 0:
                x //= p
    if x > 1:
        r -= r // x
    _phi_cache[n] = r
    return r


def pf_desc(d, spf):
    out = []
    while d > 1:
        p = spf[d]
        out.append(p)
        while d % p == 0:
            d //= p
    return sorted(out, reverse=True)


def pf_trial(d, primes):
    """Factores primos distintos de d, decrecientes, por division trial."""
    out = []
    x = d
    for p in primes:
        if p * p > x:
            break
        if x % p == 0:
            out.append(p)
            while x % p == 0:
                x //= p
    if x > 1:
        out.append(x)
    return sorted(out, reverse=True)


def in_Vplus(d, D, pdesc):
    if d == 1:
        return True
    if not pdesc:
        return False
    if pdesc[0] >= D:
        return False
    prod = 1
    for i, p in enumerate(pdesc):
        m = i + 1
        prod *= p
        if m % 2 == 1:
            if p * p * prod >= D:
                return False
    return True


def divisors_from_pf(pdesc):
    divs = [1]
    for p in pdesc:
        divs += [d * p for d in list(divs)]
    return divs


def divisors_of(n, spf):
    if n == 1:
        return [1]
    out = [1]
    x = n
    facs = []
    while x > 1:
        p = spf[x]
        e = 0
        while x % p == 0:
            x //= p
            e += 1
        facs.append((p, e))
    for (p, e) in facs:
        cur = list(out)
        pw = 1
        for _ in range(e):
            pw *= p
            out += [d * pw for d in cur]
    return sorted(out)


SMALL_PRIMES = primes_upto(10000)

# ------------------------------------------------------------------ modelo
# (Model97 copiado del script del doc 97/100 sin cambios)


class Model97:
    def __init__(self, A, z, y, D, mode="fixed", excl=(), full=True,
                 s_list=None):
        t0 = time.time()
        self.A, self.z, self.y, self.D, self.mode = A, z, y, D, mode
        self.excl = set(excl)
        self.full = full
        self.phiA = phi_of(A, SMALL_PRIMES)
        allpr = primes_upto(max(y, 100))
        self.sprimes = [p for p in allpr
                        if p >= z and p < y and math.gcd(p, A) == 1
                        and p not in self.excl]
        if s_list is not None:
            self.sprimes = [s for s in self.sprimes if s in set(s_list)]
        self.small_primes = SMALL_PRIMES
        maxd = D if full else D // z
        mu, spf, _ = mu_spf(maxd)
        self.mu, self.spf = mu, spf
        self.Vs = {}
        for s in self.sprimes:
            Ds = D // s
            Plim = z if mode == "fixed" else s
            Vs = []
            if mode == "fixed":
                allow = [p for p in allpr
                         if 2 < p < z and math.gcd(p, A) == 1
                         and p not in self.excl]
                divs = [1]
                for p in allow:
                    divs += [d * p for d in list(divs)]
                for d in divs:
                    if d >= Ds:
                        continue
                    if d == 1 or in_Vplus(d, Ds, pf_desc(d, spf)):
                        Vs.append(d)
            else:
                for d in range(1, Ds):
                    if mu[d] == 0 or d % 2 == 0 or math.gcd(d, A) != 1:
                        continue
                    if any(d % p == 0 for p in self.excl):
                        continue
                    pf = pf_desc(d, spf)
                    if pf and pf[0] >= Plim:
                        continue
                    if in_Vplus(d, Ds, pf):
                        Vs.append(d)
            self.Vs[s] = sorted(Vs)
        self.build_t = time.time() - t0
        if not full:
            return
        self.b = {}
        self.pairs_by_m = {}
        for s in self.sprimes:
            for d in self.Vs[s]:
                m = s * d
                self.b[m] = self.b.get(m, 0) + mu[d]
                self.pairs_by_m.setdefault(m, []).append((s, d))
        self.M = sorted(self.b)
        H = {1}
        for m in self.M:
            pf = pf_desc(m, spf)
            for h in divisors_from_pf(pf):
                H.add(h)
        self.H = sorted(H)
        self.K = {}
        for h in self.H:
            ssum = 0.0
            for m in self.M:
                if m % h == 0:
                    ssum += self.b[m] / phi_of(m, self.small_primes)
            self.K[h] = ssum / self.phiA
        self.Hsmooth = []
        for h in self.H:
            pf = pf_desc(h, spf) if h > 1 else []
            if all(p < z for p in pf):
                self.Hsmooth.append(h)
        self.T = {}
        self.Q = {}
        for s in self.sprimes:
            for h0 in self.Hsmooth:
                t = 0.0
                for d in self.Vs[s]:
                    if d % h0 == 0:
                        t += mu[d] / phi_of(d, self.small_primes)
                self.T[(s, h0)] = t
            q = 0.0
            for h0 in self.Hsmooth:
                q += phi_of(h0, self.small_primes) * self.T[(s, h0)] ** 2
            self.Q[s] = q
        self.build_t = time.time() - t0

    def G_A(self, Hset=None):
        Hset = self.H if Hset is None else Hset
        return sum(phi_of(h, self.small_primes) * self.K[h] ** 2
                   for h in Hset)


# ------------------------------------------------------------------ psi
def lambda_sieve(N):
    lam = [0.0] * (N + 1)
    is_p = [True] * (N + 1)
    is_p[0] = is_p[1] = False
    for i in range(2, int(N ** 0.5) + 1):
        if is_p[i]:
            for j in range(i * i, N + 1, i):
                is_p[j] = False
    for i in range(2, N + 1):
        if is_p[i]:
            lam[i] = math.log(i)
            pw = i * i
            while pw <= N:
                lam[pw] = math.log(i)
                pw *= i
    return lam


class PsiEngine:
    """psi(N;d,a) real con cache (doc 100, sin cambios)."""

    def __init__(self, N):
        self.N = N
        self.lam = lambda_sieve(N)
        self.psiN = sum(self.lam)
        self.cache = {}

    def psi_mod(self, d, a):
        a %= d
        key = (d, a)
        if key in self.cache:
            return self.cache[key]
        lam = self.lam
        N = self.N
        s = 0.0
        n = a
        while n < 2:
            n += d
        while n <= N:
            s += lam[n]
            n += d
        self.cache[key] = s
        return s

    def E(self, d):
        if d == 1:
            return 0.0
        return self.psi_mod(d, self.N) - self.psiN / phi_of(d, SMALL_PRIMES)


def prime_powers_j2(N):
    """Lista de (q, log p) con q = p^j <= N, j >= 2."""
    out = []
    for p in primes_upto(int(N ** 0.5) + 1):
        q = p * p
        while q <= N:
            out.append((q, math.log(p)))
            if q > N // p:
                break
            q *= p
    return out


class SegPsi:
    """psi(N; d, N mod d) para muchos modulos, por criba segmentada.

    Equivalente a PsiEngine pero sin materializar Lambda sobre [0,N]:
    procesa bloques, acumula sum_{n<=N, n=N (mod d)} log p sobre primos
    (via rebanadas de array) y anade las potencias de primos j>=2 aparte.
    Devuelve .E (dict d -> E(N;d,N)) y .psiN."""

    def __init__(self, N, moduli, block=2_000_000):
        self.N = N
        mods = sorted(set(m for m in moduli if m > 1))
        self.mods = mods
        acc = {d: 0.0 for d in mods}
        psiN = 0.0
        sprimes = primes_upto(int(N ** 0.5) + 1)
        lo = 1
        while lo <= N:
            hi = min(N + 1, lo + block)
            size = hi - lo
            is_p = bytearray(b"\x01") * size
            if lo == 1:
                is_p[0] = 0
            for p in sprimes:
                if p * p >= hi:
                    break
                start = max(p * p, ((lo + p - 1) // p) * p)
                if start < hi:
                    is_p[start - lo:hi - lo:p] = \
                        b"\x00" * ((hi - 1 - start) // p + 1)
            lam = array('d', bytes(8 * size))
            idx = is_p.find(1)
            while idx != -1:
                lam[idx] = math.log(lo + idx)
                idx = is_p.find(1, idx + 1)
            psiN += sum(lam)
            for d in mods:
                a = N % d
                i0 = (a - lo) % d
                if i0 < size:
                    acc[d] += sum(lam[i0::d])
            lo = hi
        pp2 = prime_powers_j2(N)
        psiN += sum(lp for (_, lp) in pp2)
        for (q, lp) in pp2:
            qm = None
            for d in mods:
                if q % d == N % d:
                    acc[d] += lp
        self.psiN = psiN
        self.E = {1: 0.0}
        for d in mods:
            self.E[d] = acc[d] - psiN / phi_of(d, SMALL_PRIMES)


# ------------------------------------------------------- segundo factor / BV
def M_r_E(r, E, mu, spf):
    """M_r con una tabla E (dict): sum_{d|r} mu(r/d) phi(d) E_d."""
    tot = 0.0
    for d in divisors_of(r, spf):
        tot += mu[r // d] * phi_of(d, SMALL_PRIMES) * E.get(d, 0.0)
    return tot


def bv_from_E(mod, E, N, mu, spf):
    """Version de bv_decomposition (doc 100) con tabla E dada:
    divs = {d : d|r, r en H, (r,N)=1}; W_d = sum_{r ni d} r/phi(r);
    S_BV = sum_d phi(d) |E_d|^2. Devuelve dict con W_max, argmax,
    S_BV, #divs, |H_N| y la tabla W."""
    divs = set()
    HN = []
    for r in mod.H:
        if math.gcd(r, N) != 1:
            continue
        HN.append(r)
        for d in divisors_of(r, spf):
            divs.add(d)
    divs = sorted(divs)
    S_BV = 0.0
    for d in divs:
        if d > 1:
            e = E.get(d, 0.0)
            S_BV += phi_of(d, SMALL_PRIMES) * e * e
    W = {}
    for d in divs:
        w = 0.0
        for r in HN:
            if r % d == 0:
                w += r / phi_of(r, SMALL_PRIMES)
        W[d] = w
    W_max = max(W.values()) if W else 0.0
    argmax = max(W, key=W.get) if W else None
    return {"W_max": W_max, "argmax": argmax, "S_BV": S_BV,
            "ndivs": len(divs), "nHN": len(HN), "divs": divs, "W": W}


def second_factor_E(mod, E, N, kappa, mu, spf):
    """B_M sobre r en H, r > (ln N)^kappa, (r,N)=1 (def. doc 100)."""
    R0 = math.log(N) ** kappa
    BM = 0.0
    n = 0
    for r in mod.H:
        if r <= R0 or math.gcd(r, N) != 1:
            continue
        Mr = M_r_E(r, E, mu, spf)
        BM += Mr * Mr / phi_of(r, SMALL_PRIMES)
        n += 1
    return BM, n


def bv_decomposition(mod, eng, N, mu, spf):
    """bv_decomposition del doc 100 (con PsiEngine), para J0/J2."""
    divs = set()
    for r in mod.H:
        if math.gcd(r, N) != 1:
            continue
        for d in divisors_of(r, spf):
            divs.add(d)
    divs = sorted(divs)
    E2 = {}
    for d in divs:
        E2[d] = 0.0 if d == 1 else eng.E(d) ** 2
    S_BV = sum(phi_of(d, SMALL_PRIMES) * E2[d] for d in divs if d > 1)
    W = {}
    for d in divs:
        w = 0.0
        for r in mod.H:
            if r % d == 0 and math.gcd(r, N) == 1:
                w += r / phi_of(r, SMALL_PRIMES)
        W[d] = w
    W_max = max(W.values()) if W else 0.0
    return W_max, S_BV, len(divs)


# ------------------------------------------------------------------ conteo
def subset_table(z, A):
    """Tabla de subconjuntos: para cada d libre de cuadrados con primos
    en (2, z) y (d, A) = 1: (d, pf_desc, d/phi(d)). Se calcula UNA vez
    por z; los niveles L solo filtran."""
    allow = [p for p in primes_upto(z + 1)
             if 2 < p < z and math.gcd(p, A) == 1]
    divs = [1]
    for p in allow:
        divs += [d * p for d in list(divs)]
    tab = []
    for d in divs:
        pf = pf_trial(d, SMALL_PRIMES) if d > 1 else []
        ratio = 1.0
        for p in pf:
            ratio *= p / (p - 1.0)
        tab.append((d, pf, ratio))
    return tab


def vplus_from_table(tab, L):
    """( |V^+(L)| , sum_{h en V^+(L)} h/phi(h) , #{d: d<L} sin V^+ )
    filtrando la tabla de subconjuntos con d < L y la condicion V^+."""
    c = 0
    g = 0.0
    nsf = 0
    for (d, pf, ratio) in tab:
        if d >= L:
            continue
        nsf += 1
        if d == 1 or in_Vplus(d, L, pf):
            c += 1
            g += ratio
    return c, g, nsf


def sprimes_of(A, z, y):
    return [p for p in primes_upto(y + 1)
            if p >= z and p < y and math.gcd(p, A) == 1]


def formula_counts(A, z, y, D, tab):
    """A100.1: |H| y W_1 por la formula de descomposicion.
    Devuelve (|H|, W_1, #sprimes, |H|_env, W_env) donde _env es la
    envolvente sin condiciones V^+ (solo d < L)."""
    sprimes = sprimes_of(A, z, y)
    smin = sprimes[0]
    tot_H = 0
    W1 = 0.0
    tot_env = 0
    cache = {}
    for s in sprimes:
        L = D // s
        if L not in cache:
            cache[L] = vplus_from_table(tab, L)
        c, g, nsf = cache[L]
        tot_H += c
        W1 += (s / (s - 1.0)) * g
        tot_env += nsf
    c0, g0, nsf0 = cache[D // smin]
    tot_H += c0
    W1 += g0
    tot_env += nsf0
    return tot_H, W1, len(sprimes), tot_env


def explicit_H(A, z, y, D, tab):
    """Enumeracion explicita del soporte: pares (s,d) -> M -> H.
    Devuelve (|H|, W_1 = sum_{r en H} r/phi(r)). Independiente de la
    formula: construye H como conjunto de divisores de los m = s*d."""
    sprimes = sprimes_of(A, z, y)
    H = {}
    for s in sprimes:
        L = D // s
        for (d, pf, ratio) in tab:
            if d >= L:
                continue
            if not (d == 1 or in_Vplus(d, L, pf)):
                continue
            subs = [(1, 1.0)]
            for p in [s] + pf:
                subs += [(v * p, q * p / (p - 1.0))
                         for (v, q) in list(subs)]
            for (v, q) in subs:
                if v not in H:
                    H[v] = q
    return len(H), sum(H.values())


# ------------------------------------------------------------------ regresion
T_CRIT = {1: 12.706, 2: 4.303, 3: 3.182, 4: 2.776, 5: 2.571,
          6: 2.447, 7: 2.365, 8: 2.306, 9: 2.262, 10: 2.228}


def ols(xs, ys):
    """Minimos cuadrados y = a + b x. Devuelve (b, a, se_b, ci95_b,
    loo_min, loo_max) con loo = rango de pendientes leave-one-out."""
    n = len(xs)
    mx = sum(xs) / n
    my = sum(ys) / n
    sxx = sum((x - mx) ** 2 for x in xs)
    sxy = sum((x - mx) * (y - my) for (x, y) in zip(xs, ys))
    b = sxy / sxx
    a = my - b * mx
    resid = [y - (a + b * x) for (x, y) in zip(xs, ys)]
    if n > 2:
        s2 = sum(r * r for r in resid) / (n - 2)
        se = math.sqrt(s2 / sxx)
        ci = T_CRIT.get(n - 2, 2.0) * se
    else:
        se = float("nan")
        ci = float("nan")
    loo = []
    for i in range(n):
        xs2 = xs[:i] + xs[i + 1:]
        ys2 = ys[:i] + ys[i + 1:]
        m2x = sum(xs2) / (n - 1)
        m2y = sum(ys2) / (n - 1)
        s2xx = sum((x - m2x) ** 2 for x in xs2)
        s2xy = sum((x - m2x) * (y - m2y) for (x, y) in zip(xs2, ys2))
        loo.append(s2xy / s2xx)
    return b, a, se, ci, min(loo), max(loo)


MU, SPF, _ = mu_spf(600000)


# ------------------------------------------------------------------ chequeos
def check_J0(m2):
    print("J0 validacion de la maquinaria:")
    ok = True
    if len(m2.H) != 880:
        ok = False
    print("     Model97 R2: |H|=%d (doc 97/100: 880) %s"
          % (len(m2.H), "OK" if len(m2.H) == 880 else "FAIL"))
    eng = PsiEngine(999983)
    rel = abs(eng.psiN - 999586.5975) / 999586.5975
    print("     PsiEngine psi(999983)=%.4f (doc 97: 999586.5975, "
          "rel=%.2e) %s" % (eng.psiN, rel, "OK" if rel < 1e-9 else "FAIL"))
    ok &= rel < 1e-9
    W_max, S_BV, nd = bv_decomposition(m2, eng, 999983, MU, SPF)
    print("     bv_decomposition R2 (N=999983): W_max=%.4f (doc 100: "
          "1011.28), S_BV=%.6e (doc 100: 4.3593e9), #divs=%d"
          % (W_max, S_BV, nd))
    okw = abs(W_max - 1011.28) < 0.01
    oks = abs(S_BV - 4.3593e9) / 4.3593e9 < 1e-3
    print("     W_max %s | S_BV %s" % ("OK" if okw else "FAIL",
                                       "OK" if oks else "FAIL"))
    ok &= okw and oks
    # SegPsi vs PsiEngine en E_d, N par 999982
    eng2 = PsiEngine(999982)
    divs = set()
    for r in m2.H:
        if math.gcd(r, 999982) == 1:
            for d in divisors_of(r, SPF):
                divs.add(d)
    seg = SegPsi(999982, sorted(divs))
    maxrel = 0.0
    worst = None
    for d in divs:
        if d == 1:
            continue
        e1 = eng2.E(d)
        e2 = seg.E[d]
        rr = abs(e1 - e2) / max(1.0, abs(e1))
        if rr > maxrel:
            maxrel, worst = rr, d
    print("     SegPsi vs PsiEngine (N=999982, %d modulos): max rel "
          "diff en E_d = %.2e (peor d=%s) %s"
          % (len(divs), maxrel, worst,
             "OK" if maxrel < 1e-9 else "FAIL"))
    ok &= maxrel < 1e-9
    print("     J0 %s" % ("PASS" if ok else "FAIL"))
    return ok


def check_J1(models, tabs_scaled):
    print("J1 V^+ cerrado por divisores y monotono en el nivel:")
    ok = True
    for (mod, tag) in models:
        smax = max(mod.Vs, key=lambda s: mod.D // s)
        base = set(mod.Vs[smax])
        for s, vs in mod.Vs.items():
            S = set(vs)
            for d in vs:
                for h in divisors_of(d, SPF):
                    if h not in S:
                        ok = False
                        print("     [%s] NO cerrado: s=%d d=%d, divisor "
                              "%d fuera" % (tag, s, d, h))
            if not S <= base and s != smax:
                # monotonia: V^+(D/s) ⊆ V^+(D/s_min) (nivel mayor)
                if not S.issubset(set(mod.Vs[mod.sprimes[0]])):
                    ok = False
                    print("     [%s] NO monotono: s=%d" % (tag, s))
        print("     [%s] cerrado+monotono en %d niveles: OK"
              % (tag, len(mod.Vs)))
    for (z, tab) in tabs_scaled:
        D = z ** 5
        sprimes = sprimes_of(3, z, int(round(z ** (10.0 / 3.0))))
        ref = None
        for s in sprimes:
            L = D // s
            S = set()
            for (d, pf, ratio) in tab:
                if d < L and (d == 1 or in_Vplus(d, L, pf)):
                    S.add(d)
            for d in list(S):
                for h in divisors_from_pf(pf_trial(d, SMALL_PRIMES)):
                    if h not in S:
                        ok = False
                        print("     [escalado z=%d] NO cerrado: L=%d "
                              "d=%d, divisor %d fuera" % (z, L, d, h))
            if ref is None:
                ref = S
            elif not S <= ref:
                ok = False
                print("     [escalado z=%d] NO monotono en s=%d" % (z, s))
        print("     [escalado z=%d] cerrado+monotono en %d niveles: OK"
              % (z, len(sprimes)))
    print("     J1 %s" % ("PASS" if ok else "FAIL"))
    return ok


def check_J2(models):
    print("J2 formula de conteo A100.1 vs enumeracion (Model97 + BV):")
    eng = PsiEngine(999983)
    ok = True
    for (mod, tag, A, z, y, D) in models:
        tab = subset_table(z, A)
        fH, fW, ns, env = formula_counts(A, z, y, D, tab)
        W_max, S_BV, nd = bv_decomposition(mod, eng, 999983, MU, SPF)
        match_H = (fH == len(mod.H))
        relW = abs(fW - W_max) / W_max
        ok &= match_H and relW < 1e-9
        print("     [%s] |H| formula=%d vs Model97=%d %s | W_1 formula="
              "%.6f vs W_max(BV)=%.6f (rel=%.2e) %s"
              % (tag, fH, len(mod.H), "OK" if match_H else "FAIL",
                 fW, W_max, relW, "OK" if relW < 1e-9 else "FAIL"))
    print("     J2 %s" % ("PASS" if ok else "FAIL"))
    return ok


def check_J3_J4():
    print("J3/J4 familia escalada D=z^5, y=z^(10/3), A=3 (z = N^(1/10)):")
    zs = [11, 13, 17, 19, 23, 29]
    rows = []
    ok = True
    for z in zs:
        D = z ** 5
        y = int(round(z ** (10.0 / 3.0)))
        tab = subset_table(z, 3)
        fH, fW, ns, env = formula_counts(3, z, y, D, tab)
        extra = ""
        if z in (13, 17, 19, 23):
            eH, eW = explicit_H(3, z, y, D, tab)
            good = (eH == fH and abs(eW - fW) / fW < 1e-9)
            ok &= good
            extra = " | explicito |H|=%d W_1=%.4f %s" % (
                eH, eW, "OK" if good else "DISCREPA")
        rows.append((z, D, y, ns, fH, fW, env))
        print("     z=%2d D=%9d y=%6d #s=%5d |H|=%7d W_1=%12.2f "
              "env_sinV+=%8d%s" % (z, D, y, ns, fH, fW, env, extra))
    xs = [math.log(z) for (z, *_rest) in rows]
    bW, aW, seW, ciW, looW0, looW1 = ols(xs, [math.log(r[5]) for r in rows])
    bH, aH, seH, ciH, looH0, looH1 = ols(xs, [math.log(r[4]) for r in rows])
    bE, aE, seE, ciE, looE0, looE1 = ols(xs, [math.log(r[6]) for r in rows])
    print("     Regresion ln W_1 ~ ln z: pendiente=%.4f +- %.4f (SE), "
          "IC95%% +-%.4f, LOO [%.4f, %.4f]" % (bW, seW, ciW, looW0, looW1))
    print("       => alpha_W = pendiente/10 = %.4f, IC95%% [%.4f, %.4f], "
          "LOO [%.4f, %.4f]"
          % (bW / 10, (bW - ciW) / 10, (bW + ciW) / 10,
             looW0 / 10, looW1 / 10))
    print("     Regresion ln |H| ~ ln z: pendiente=%.4f +- %.4f (SE), "
          "IC95%% +-%.4f, LOO [%.4f, %.4f]"
          % (bH, seH, ciH, looH0, looH1))
    print("       => alpha_H = %.4f, IC95%% [%.4f, %.4f], LOO [%.4f, %.4f]"
          % (bH / 10, (bH - ciH) / 10, (bH + ciH) / 10,
             looH0 / 10, looH1 / 10))
    print("     Envolvente sin V^+ (solo d<L): pendiente=%.4f => "
          "alpha_env=%.4f (la envolvente NO decide: ver doc sec. 2.4)"
          % (bE, bE / 10))
    loc = []
    for i in range(len(rows) - 1):
        loc.append((math.log(rows[i + 1][5]) - math.log(rows[i][5]))
                   / (math.log(rows[i + 1][0]) - math.log(rows[i][0])) / 10)
    print("     alpha locales (W_1, entre z consecutivos): %s"
          % ", ".join("%.3f" % v for v in loc))
    print("     J3 %s | J4 (reporte) hecho" % ("PASS" if ok else "FAIL"))
    return ok, {"alpha_W": bW / 10, "ci_W": ciW / 10,
                "loo_W": (looW0 / 10, looW1 / 10),
                "alpha_H": bH / 10, "ci_H": ciH / 10,
                "alpha_env": bE / 10, "rows": rows}


def check_J5_J6(m2):
    print("J5/J6 E100.1 parte A: soporte fijo R2 (doc 94), N par "
          "creciente, psi segmentado:")
    Ns = [(510510, "primorial 2*3*5*7*11*13*17"),
          (999982, "2*79*6329 (79=sprime)"),
          (1048576, "2^20"),
          (9699690, "primorial ...*19"),
          (16777216, "2^24"),
          (67108864, "2^26"),
          (77597520, "2^4 x primorial 19"),
          (99999998, "2*7*23*310559")]
    rows = []
    ok6 = True
    base_W = None
    for (N, fam) in Ns:
        t0 = time.time()
        divs = set()
        for r in m2.H:
            if math.gcd(r, N) == 1:
                for d in divisors_of(r, SPF):
                    divs.add(d)
        seg = SegPsi(N, sorted(divs))
        bv = bv_from_E(m2, seg.E, N, MU, SPF)
        BM, nBM = second_factor_E(m2, seg.E, N, 1.0, MU, SPF)
        delta = bv["S_BV"] / (bv["ndivs"] * N * math.log(N))
        pred32 = N ** 1.5 * math.log(N)
        rows.append((N, fam, bv, BM, delta, seg.psiN))
        if bv["argmax"] != 1:
            ok6 = False
        if N == 999982:
            base_W = bv["W_max"]
        print("     N=%9d (%s): |H_N|=%d #divs=%d W_max=%.2f (argmax "
              "d=%s) S_BV=%.6e B_M=%.6e delta=%.4f S_BV/(N^1.5 lnN)="
              "%.3e psi/N=%.4f (%.0fs)"
              % (N, fam, bv["nHN"], bv["ndivs"], bv["W_max"],
                 bv["argmax"], bv["S_BV"], BM, delta,
                 bv["S_BV"] / pred32, seg.psiN / N, time.time() - t0))
    xs = [math.log(r[0]) for r in rows]
    bS, aS, seS, ciS, looS0, looS1 = ols(
        xs, [math.log(r[2]["S_BV"]) for r in rows])
    bB, aB, seB, ciB, looB0, looB1 = ols(
        xs, [math.log(r[3]) for r in rows])
    print("     Regresion ln S_BV ~ ln N (soporte fijo): beta=%.4f "
          "+- %.4f (SE), IC95%% +-%.4f, LOO [%.4f, %.4f]"
          % (bS, seS, ciS, looS0, looS1))
    print("     Regresion ln B_M ~ ln N (soporte fijo): beta=%.4f "
          "+- %.4f (SE), IC95%% +-%.4f, LOO [%.4f, %.4f]"
          % (bB, seB, ciB, looB0, looB1))
    # subfamilias puras (por conjunto de N, sin ambiguedad de etiquetas)
    for keyset, tag in [({1048576, 16777216, 67108864}, "2^e puras"),
                        ({510510, 9699690, 77597520}, "primoriales")]:
        sub = [r for r in rows if r[0] in keyset]
        if len(sub) >= 3:
            bb, aa, see, cii, l0, l1 = ols(
                [math.log(r[0]) for r in sub],
                [math.log(r[2]["S_BV"]) for r in sub])
            print("     Subfamilia %s (n=%d): beta_S_BV=%.4f, LOO "
                  "[%.4f, %.4f] (IC amplio, n pequeno)"
                  % (tag, len(sub), bb, l0, l1))
    wprim = [r[2]["W_max"] for r in rows if "primorial" in r[1]]
    print("     Colapso primorial de W_max: base (N=999982) "
          "=%.2f; valores primoriales=%s"
          % (base_W, ", ".join("%.2f" % w for w in wprim)))
    print("     J5 (reporte) hecho | J6 W_max=W_1 (argmax d=1 en todo N)"
          " %s" % ("PASS" if ok6 else "FAIL"))
    return ok6, {"beta_S": bS, "ci_S": ciS, "beta_B": bB, "rows": rows}


# ------------------------------------------------------------------ main
def main():
    allok = True
    t0 = time.time()
    print("== DOC 101: experimento E100.1 + conteo A100.1 ==")
    print("Solo stdlib. MEDICION, no prueba. Goldbach sigue abierta.")
    print()
    m1 = Model97(3, 11, 60, 20000, "fixed")
    m2 = Model97(3, 11, 3000, 200000, "fixed")
    m3 = Model97(3, 3, 349, 200000, "fixed")
    print("Modelos: R1 |H|=%d | R2 |H|=%d (doc 94) | R3 |H|=%d"
          % (len(m1.H), len(m2.H), len(m3.H)))
    print()
    allok &= check_J0(m2)
    print()
    tabs = [(13, subset_table(13, 3)), (17, subset_table(17, 3))]
    allok &= check_J1([(m1, "R1"), (m2, "R2"), (m3, "R3")], tabs)
    print()
    allok &= check_J2([(m1, "R1", 3, 11, 60, 20000),
                       (m2, "R2", 3, 11, 3000, 200000),
                       (m3, "R3", 3, 3, 349, 200000)])
    print()
    ok34, res34 = check_J3_J4()
    allok &= ok34
    print()
    ok56, res56 = check_J5_J6(m2)
    allok &= ok56
    print()
    print("RESUMEN: alpha_W=%.4f (IC95%% +-%.4f) | alpha_H=%.4f | "
          "beta_S_BV(soporte fijo)=%.4f (IC95%% +-%.4f)"
          % (res34["alpha_W"], res34["ci_W"], res34["alpha_H"],
             res56["beta_S"], res56["ci_S"]))
    print("PASS" if allok else "FAIL (ver lineas marcadas)")
    print("Tiempo total: %.1f s" % (time.time() - t0))
    return 0 if allok else 1


if __name__ == "__main__":
    sys.exit(main())
