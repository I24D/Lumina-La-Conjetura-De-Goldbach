#!/usr/bin/env python3
# 102_experimentos_e101_1_e101_2.py  (Doc 102)
# Experimentos E101.1 y E101.2 propuestos en el doc 101.
# Solo stdlib.
#
# E101.1 (holgura de Cauchy): a N fijo, sobre soportes escalados,
#   se mide la holgura de la separacion de Cauchy
#       slack = W_max * S_BV / B_M
#   y su descomposicion en la cadena (cada eslabon >= 1)
#       B_M <= BV1* <= BV1_full <= W_max * S_BV,
#   donde BV1* = sum_{r en el conjunto de B_M} (r/phi(r)) *
#                sum_{d|r} phi(d) |E_d|^2   (Cauchy por conductor),
#         BV1_full = sum_d phi(d) |E_d|^2 W_d  (mismos terminos,
#                todos los r de H_N, sin corte kappa).
#   slack = S_Cauchy * S_cut * S_max con
#       S_Cauchy = BV1*/B_M, S_cut = BV1_full/BV1*,
#       S_max = W_max*S_BV/BV1_full.
#   Descomposicion primo/compuesto: B_M = B_M^P + B_M^C (r primo /
#   compuesto), BV1* analogo, Diag_r = (1/phi(r)) sum_{d|r} phi(d)^2
#   |E_d|^2 (sin terminos cruzados; Diag^P = B_M^P exacto [T]).
#   Familias de N: 999982 (generico par), 1048576 = 2^20 (soporte
#   intacto), 720720 = 2^4*3^2*5*7*11*13 (muchos factores pequenos:
#   control de colapso primorial, B_M^C = 0 y B_M = S_BV [T]).
#   Soportes: R1, R3, R2 (docs 94/100) y escalados z=11,13,17
#   (familia D=z^5, y=z^(10/3) del doc 101).
#
# E101.2 (alpha en z mayor): formula de conteo A100.1 en z=31 y z=37;
#   alpha_W y alpha_H con IC95% y pendientes locales, comparando el
#   ajuste z=11..29 (doc 101) con z=11..37 y con la mitad superior.
#
# Chequeos:
#   K0  Validacion: Model97 R2 |H|=880; build_H == Model97 (R2 y
#       escalado z=11); SegPsi reproduce el doc 101 en R2/N=999982
#       (W_max=1006.34, S_BV=4.139e9, B_M=3.777e9); SegPsi vs
#       PsiEngine en una muestra de E_d; formula de conteo
#       reproduce la tabla escalada del doc 101 (z=11..29).
#   K1  Identidades: Diag^P = B_M^P; en N=720720: B_M^C = 0 y
#       B_M = S_BV (colapso primorial) en todos los soportes.
#   K2  Cadena B_M <= BV1* <= BV1_full <= W_max*S_BV en todo (N, sop).
#   K3  argmax W_d = 1 (W_max = W_1 sobre H_N) en todo (N, soporte).
#   K4  E101.2: alpha_W/alpha_H en z=11..29 reproducen el doc 101
#       (|dal| < 0.002); variante ordenada de vplus == formula del
#       doc 101 en z=29 (|H| exacto, W_1 rel < 1e-9).
#   K5  E101.1: tablas completas y regresiones theta (reporte).
# Uso: python3 102_experimentos_e101_1_e101_2.py
import math
import sys
import time
from array import array


# ------------------------------------------------------------------ utiles
# (copiados de los scripts de los docs 97/100/101, misma convencion)
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
# (Model97 copiado del script del doc 97/100/101 sin cambios)


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


def build_H(A, z, y, D, mu, spf):
    """Version ligera de Model97 (modo fixed): devuelve (H, sprimes,
    Vs) replicando exactamente la construccion de H de Model97
    (pares (s,d) -> m -> cierre por divisores), sin K/T/Q."""
    allpr = primes_upto(max(y, 100))
    sprimes = [p for p in allpr
               if p >= z and p < y and math.gcd(p, A) == 1]
    allow = [p for p in allpr if 2 < p < z and math.gcd(p, A) == 1]
    divs = [1]
    for p in allow:
        divs += [d * p for d in list(divs)]
    Vs = {}
    Mset = set()
    for s in sprimes:
        Ds = D // s
        Vs_s = []
        for d in divs:
            if d >= Ds:
                continue
            if d == 1 or in_Vplus(d, Ds, pf_desc(d, spf)):
                Vs_s.append(d)
                Mset.add(s * d)
        Vs[s] = sorted(Vs_s)
    H = {1}
    for m in Mset:
        for h in divisors_from_pf(pf_desc(m, spf)):
            H.add(h)
    return sorted(H), sprimes, Vs


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
    """psi(N;d,a) real con cache (docs 97/100, sin cambios)."""

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
    """psi(N; d, N mod d) para muchos modulos, criba segmentada
    (doc 101, sin cambios). Devuelve .E (dict) y .psiN."""

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
            for d in mods:
                if q % d == N % d:
                    acc[d] += lp
        self.psiN = psiN
        self.E = {1: 0.0}
        for d in mods:
            self.E[d] = acc[d] - psiN / phi_of(d, SMALL_PRIMES)


# ------------------------------------------------------- segundo factor / BV
def M_r_E(r, E, mu, spf):
    tot = 0.0
    for d in divisors_of(r, spf):
        tot += mu[r // d] * phi_of(d, SMALL_PRIMES) * E.get(d, 0.0)
    return tot


# ------------------------------------------------------------------ conteo
def subset_table(z, A):
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
    """A100.1 (doc 101): (|H|, W_1, #sprimes, |H|_env)."""
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


def vplus_sorted(tab_sorted, L):
    """Como vplus_from_table pero con la tabla ordenada por d y
    corte temprano. Misma semantica; se valida en K4."""
    c = 0
    g = 0.0
    nsf = 0
    for (d, pf, ratio) in tab_sorted:
        if d >= L:
            break
        nsf += 1
        if d == 1 or in_Vplus(d, L, pf):
            c += 1
            g += ratio
    return c, g, nsf


def formula_counts_fast(A, z, y, D, tab_sorted):
    sprimes = sprimes_of(A, z, y)
    smin = sprimes[0]
    tot_H = 0
    W1 = 0.0
    cache = {}
    for s in sprimes:
        L = D // s
        if L not in cache:
            cache[L] = vplus_sorted(tab_sorted, L)
        c, g, nsf = cache[L]
        tot_H += c
        W1 += (s / (s - 1.0)) * g
    c0, g0, nsf0 = cache[D // smin]
    tot_H += c0
    W1 += g0
    return tot_H, W1, len(sprimes)


# ------------------------------------------------------------------ regresion
T_CRIT = {1: 12.706, 2: 4.303, 3: 3.182, 4: 2.776, 5: 2.571,
          6: 2.447, 7: 2.365, 8: 2.306, 9: 2.262, 10: 2.228}


def ols(xs, ys):
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


MU, SPF, _ = mu_spf(1_420_000)


def is_prime_spf(n):
    return n > 1 and SPF[n] == n


# ------------------------------------------------------------------ medida
def measure(H, N, kappa=1.0):
    """Mide la cadena BV completa sobre el soporte H a N dado.
    Devuelve un dict con todas las descomposiciones de E101.1."""
    R0 = math.log(N) ** kappa
    HN = [r for r in H if math.gcd(r, N) == 1]
    divs_set = set()
    for r in HN:
        for d in divisors_of(r, SPF):
            divs_set.add(d)
    divs = sorted(divs_set)
    seg = SegPsi(N, divs)
    E = seg.E
    phi = {v: phi_of(v, SMALL_PRIMES) for v in set(divs) | set(HN)}
    # W_d por acumulacion sobre divisores de cada r
    W = {d: 0.0 for d in divs}
    for r in HN:
        rr = r / phi[r]
        for d in divisors_of(r, SPF):
            W[d] += rr
    W_max = max(W.values())
    argmax = max(W, key=W.get)
    S_BV = 0.0
    S_BV_P = 0.0
    S_BV_C = 0.0
    BV1_full = 0.0
    for d in divs:
        if d == 1:
            continue
        t = phi[d] * E[d] * E[d]
        S_BV += t
        if is_prime_spf(d):
            S_BV_P += t
        else:
            S_BV_C += t
        BV1_full += t * W[d]
    B_M = B_M_P = B_M_C = 0.0
    BV1s = BV1s_P = BV1s_C = 0.0
    Diag_P = Diag_C = 0.0
    nBM = nBM_P = nBM_C = 0
    for r in HN:
        if r <= R0:
            continue
        Mr = M_r_E(r, E, MU, SPF)
        term = Mr * Mr / phi[r]
        s1 = 0.0   # sum_{d|r} phi(d) |E_d|^2
        s2 = 0.0   # sum_{d|r} phi(d)^2 |E_d|^2
        for d in divisors_of(r, SPF):
            if d == 1:
                continue
            e2 = E[d] * E[d]
            s1 += phi[d] * e2
            s2 += phi[d] * phi[d] * e2
        bound = (r / phi[r]) * s1
        diag = s2 / phi[r]
        if is_prime_spf(r):
            B_M_P += term
            BV1s_P += bound
            Diag_P += diag
            nBM_P += 1
        else:
            B_M_C += term
            BV1s_C += bound
            Diag_C += diag
            nBM_C += 1
        B_M += term
        BV1s += bound
        nBM += 1
    slack = W_max * S_BV / B_M if B_M > 0 else float("inf")
    return {
        "N": N, "nH": len(H), "nHN": len(HN), "ndivs": len(divs),
        "W_max": W_max, "argmax": argmax, "S_BV": S_BV,
        "S_BV_P": S_BV_P, "S_BV_C": S_BV_C, "BV1_full": BV1_full,
        "B_M": B_M, "B_M_P": B_M_P, "B_M_C": B_M_C,
        "BV1s": BV1s, "BV1s_P": BV1s_P, "BV1s_C": BV1s_C,
        "Diag_P": Diag_P, "Diag_C": Diag_C,
        "nBM": nBM, "nBM_P": nBM_P, "nBM_C": nBM_C,
        "slack": slack,
        "S_Cauchy": BV1s / B_M if B_M > 0 else float("inf"),
        "S_cut": BV1_full / BV1s if BV1s > 0 else float("inf"),
        "S_max": (W_max * S_BV / BV1_full) if BV1_full > 0
                 else float("inf"),
        "psiN": seg.psiN,
    }


def fmt_row(tag, m):
    exc = m["BV1s"] - m["B_M"]
    excC = m["BV1s_C"] - m["B_M_C"]
    fracC = excC / exc if exc > 0 else float("nan")
    return ("     [%s] N=%d |H_N|=%d #divs=%d W_max=%.2f S_BV=%.6e "
            "B_M=%.6e slack=%.1f | S_Cauchy=%.2f S_cut=%.3f S_max=%.1f "
            "| B_M^C/B_M=%.3f BV1exc_C=%.3f Diag_C/B_M^C=%s"
            % (tag, m["N"], m["nHN"], m["ndivs"], m["W_max"],
               m["S_BV"], m["B_M"], m["slack"], m["S_Cauchy"],
               m["S_cut"], m["S_max"],
               m["B_M_C"] / m["B_M"] if m["B_M"] > 0 else float("nan"),
               fracC,
               "%.3f" % (m["Diag_C"] / m["B_M_C"])
               if m["B_M_C"] > 0 else "n/a (B_M^C=0)"))


# ------------------------------------------------------------------ chequeos
def check_K0(m2, H_s11_model, H_s11_light):
    print("K0 validacion de la maquinaria:")
    ok = True
    c1 = len(m2.H) == 880
    ok &= c1
    print("     Model97 R2: |H|=%d (doc 94/101: 880) %s"
          % (len(m2.H), "OK" if c1 else "FAIL"))
    c2 = set(H_s11_light) == set(H_s11_model)
    ok &= c2
    print("     build_H escalado z=11 == Model97: %s (|H|=%d)"
          % ("OK" if c2 else "FAIL", len(H_s11_light)))
    # SegPsi reproduce el doc 101 en R2, N=999982
    m = measure(m2.H, 999982)
    eW = abs(m["W_max"] - 1006.34) / 1006.34
    eS = abs(m["S_BV"] - 4.139e9) / 4.139e9
    eB = abs(m["B_M"] - 3.777e9) / 3.777e9
    c3 = eW < 1e-3 and eS < 1e-3 and eB < 1e-3
    ok &= c3
    print("     SegPsi R2/N=999982: W_max=%.2f (doc 101: 1006.34), "
          "S_BV=%.6e (4.139e9), B_M=%.6e (3.777e9) %s"
          % (m["W_max"], m["S_BV"], m["B_M"], "OK" if c3 else "FAIL"))
    # SegPsi vs PsiEngine en una muestra de modulos
    eng = PsiEngine(999982)
    divs = set()
    for r in m2.H:
        if math.gcd(r, 999982) == 1:
            for d in divisors_of(r, SPF):
                divs.add(d)
    seg = SegPsi(999982, sorted(divs))
    sample = sorted(divs)[::29]
    maxrel = 0.0
    for d in sample:
        if d == 1:
            continue
        e1 = eng.E(d)
        e2 = seg.E[d]
        maxrel = max(maxrel, abs(e1 - e2) / max(1.0, abs(e1)))
    c4 = maxrel < 1e-9
    ok &= c4
    print("     SegPsi vs PsiEngine (muestra de %d modulos): max rel "
          "diff = %.2e %s" % (len(sample), maxrel,
                              "OK" if c4 else "FAIL"))
    # formula de conteo reproduce la tabla escalada del doc 101
    doc101 = {11: (802, 912.79), 13: (1628, 1898.02),
              17: (5033, 6015.88), 19: (8245, 9943.43),
              23: (17999, 21933.30), 29: (47178, 58126.02)}
    c5 = True
    for z, (eH, eW) in doc101.items():
        D = z ** 5
        y = int(round(z ** (10.0 / 3.0)))
        tab = subset_table(z, 3)
        fH, fW, ns, env = formula_counts(3, z, y, D, tab)
        good = (fH == eH and abs(fW - eW) < 0.02)
        c5 &= good
        if not good:
            print("     [conteo z=%d] |H|=%d (doc 101: %d), W_1=%.2f "
                  "(doc 101: %.2f) DISCREPA" % (z, fH, eH, fW, eW))
    ok &= c5
    print("     Formula de conteo vs tabla del doc 101 (z=11..29): %s"
          % ("OK" if c5 else "FAIL"))
    print("     K0 %s" % ("PASS" if ok else "FAIL"))
    return ok, m


def check_K1_K2_K3(results):
    print("K1/K2/K3 identidades, cadena y argmax:")
    ok = True
    for (tag, m) in results:
        # K1a: Diag^P = B_M^P
        if m["B_M_P"] > 0:
            rel = abs(m["Diag_P"] - m["B_M_P"]) / m["B_M_P"]
            if rel > 1e-9:
                ok = False
                print("     [%s] Diag^P != B_M^P (rel=%.2e)"
                      % (tag, rel))
        # K2: cadena
        tol = 1e-9
        if not (m["B_M"] * (1 - tol) <= m["BV1s"] + 1e-6
                and m["BV1s"] * (1 - tol) <= m["BV1_full"] + 1e-6
                and m["BV1_full"] * (1 - tol)
                <= m["W_max"] * m["S_BV"] + 1e-6):
            ok = False
            print("     [%s] CADENA ROTA: B_M=%.6e BV1*=%.6e "
                  "BV1_full=%.6e W_max*S_BV=%.6e"
                  % (tag, m["B_M"], m["BV1s"], m["BV1_full"],
                     m["W_max"] * m["S_BV"]))
        # K3: argmax
        if m["argmax"] != 1:
            ok = False
            print("     [%s] argmax W_d = %s != 1" % (tag, m["argmax"]))
        # K1b: colapso primorial en N=720720
        if m["N"] == 720720:
            relB = abs(m["B_M"] - m["S_BV"]) / m["S_BV"]
            if m["B_M_C"] != 0.0 or relB > 1e-9:
                ok = False
                print("     [%s] COLAPSO FALLO: B_M^C=%.6e, "
                      "|B_M-S_BV|/S_BV=%.2e"
                      % (tag, m["B_M_C"], relB))
    print("     Diag^P=B_M^P en todo punto con B_M^P>0; cadena "
          "B_M<=BV1*<=BV1_full<=W_max*S_BV en todo punto; argmax=1 "
          "en todo punto; colapso primorial (B_M^C=0, B_M=S_BV) en "
          "N=720720: %s" % ("PASS" if ok else "FAIL (ver arriba)"))
    return ok


def theta_fit(points, tag):
    """points: lista de (|H_N|, slack, S_Cauchy, S_max). Ajusta
    ln slack ~ ln |H_N| y devuelve el reporte."""
    xs = [math.log(p[0]) for p in points]
    out = {}
    for idx, name in [(1, "slack"), (2, "S_Cauchy"), (3, "S_max")]:
        ys = [math.log(p[idx]) for p in points]
        b, a, se, ci, l0, l1 = ols(xs, ys)
        out[name] = (b, se, ci, l0, l1)
        print("     [%s] ln %s ~ ln |H_N|: theta=%.4f +- %.4f (SE), "
              "IC95%% +-%.4f, LOO [%.4f, %.4f]"
              % (tag, name, b, se, ci, l0, l1))
    return out


def e101_2():
    print("E101.2 alpha en z mayor (formula de conteo A100.1):")
    zs_old = [11, 13, 17, 19, 23, 29]
    zs_new = [31, 37]
    rows = {}
    ok = True
    for z in zs_old:
        D = z ** 5
        y = int(round(z ** (10.0 / 3.0)))
        tab = subset_table(z, 3)
        fH, fW, ns, env = formula_counts(3, z, y, D, tab)
        rows[z] = (fH, fW, ns)
        print("     z=%2d D=%10d y=%7d #s=%6d |H|=%8d W_1=%14.2f"
              % (z, D, y, ns, fH, fW))
    # validacion de la variante ordenada en z=29 (K4, segunda parte)
    z = 29
    D = z ** 5
    y = int(round(z ** (10.0 / 3.0)))
    tab = subset_table(z, 3)
    tab_sorted = sorted(tab, key=lambda t: t[0])
    fH2, fW2, ns2 = formula_counts_fast(3, z, y, D, tab_sorted)
    goodv = (fH2 == rows[z][0]
             and abs(fW2 - rows[z][1]) / rows[z][1] < 1e-9)
    ok &= goodv
    print("     Variante ordenada en z=29: |H|=%d, W_1=%.4f %s"
          % (fH2, fW2, "OK" if goodv else "DISCREPA"))
    for z in zs_new:
        t0 = time.time()
        D = z ** 5
        y = int(round(z ** (10.0 / 3.0)))
        tab = subset_table(z, 3)
        tab_sorted = sorted(tab, key=lambda t: t[0])
        fH, fW, ns = formula_counts_fast(3, z, y, D, tab_sorted)
        rows[z] = (fH, fW, ns)
        print("     z=%2d D=%10d y=%7d #s=%6d |H|=%8d W_1=%14.2f "
              "(%.0fs)" % (z, D, y, ns, fH, fW, time.time() - t0))
    zs_all = zs_old + zs_new

    def fit(zs, col):
        xs = [math.log(z) for z in zs]
        ys = [math.log(rows[z][col]) for z in zs]
        return ols(xs, ys)

    res = {}
    for col, name in [(1, "W"), (0, "H")]:
        for tagset, zs in [("z=11..29 (doc 101)", zs_old),
                           ("z=11..37 (extendida)", zs_all),
                           ("z=19..37 (mitad superior)",
                            [19, 23, 29, 31, 37])]:
            b, a, se, ci, l0, l1 = fit(zs, col)
            res[(name, tagset)] = (b / 10, ci / 10, l0 / 10, l1 / 10)
            print("     alpha_%s [%s]: %.4f, IC95%% [%.4f, %.4f], "
                  "LOO [%.4f, %.4f]"
                  % (name, tagset, b / 10, (b - ci) / 10,
                     (b + ci) / 10, l0 / 10, l1 / 10))
    locW = []
    locH = []
    for i in range(len(zs_all) - 1):
        z0, z1 = zs_all[i], zs_all[i + 1]
        dl = math.log(z1) - math.log(z0)
        locW.append((z0, z1, (math.log(rows[z1][1])
                              - math.log(rows[z0][1])) / dl / 10))
        locH.append((z0, z1, (math.log(rows[z1][0])
                              - math.log(rows[z0][0])) / dl / 10))
    print("     alpha locales W_1: %s"
          % ", ".join("%d->%d: %.3f" % t for t in locW))
    print("     alpha locales |H|: %s"
          % ", ".join("%d->%d: %.3f" % t for t in locH))
    # K4 primera parte: reproduccion del doc 101 en z=11..29
    aW = res[("W", "z=11..29 (doc 101)")][0]
    aH = res[("H", "z=11..29 (doc 101)")][0]
    good4 = abs(aW - 0.4289) < 0.002 and abs(aH - 0.4208) < 0.002
    ok &= good4
    print("     Reproduccion doc 101: alpha_W=%.4f (doc 101: 0.4289), "
          "alpha_H=%.4f (doc 101: 0.4208) %s"
          % (aW, aH, "OK" if good4 else "DISCREPA"))
    print("     E101.2 %s" % ("PASS" if ok else "FAIL"))
    return ok, rows, res, locW, locH


# ------------------------------------------------------------------ main
def main():
    allok = True
    t00 = time.time()
    print("== DOC 102: experimentos E101.1 + E101.2 ==")
    print("Solo stdlib. MEDICION, no prueba. Goldbach sigue abierta.")
    print()
    m1 = Model97(3, 11, 60, 20000, "fixed")
    m2 = Model97(3, 11, 3000, 200000, "fixed")
    m3 = Model97(3, 3, 349, 200000, "fixed")
    print("Modelos base: R1 |H|=%d | R2 |H|=%d | R3 |H|=%d"
          % (len(m1.H), len(m2.H), len(m3.H)))
    # soportes escalados: Model97 solo para z=11 (validacion de build_H)
    D11 = 11 ** 5
    y11 = int(round(11 ** (10.0 / 3.0)))
    m_s11 = Model97(3, 11, y11, D11, "fixed")
    supports = {}
    supports["R1"] = m1.H
    supports["R3"] = m3.H
    supports["R2"] = m2.H
    for z in (11, 13, 17):
        D = z ** 5
        y = int(round(z ** (10.0 / 3.0)))
        Hz, _, _ = build_H(3, z, y, D, MU, SPF)
        supports["S%d" % z] = Hz
        print("Soporte escalado z=%d: |H|=%d (doc 101: %s)"
              % (z, len(Hz),
                 {11: 802, 13: 1628, 17: 5033}[z]))
    print()
    ok0, _ = check_K0(m2, m_s11.H, supports["S11"])
    allok &= ok0
    print()
    print("E101.1 holgura de Cauchy por (familia de N, soporte):")
    plan = [(999982, ["R1", "R3", "S11", "R2", "S13", "S17"],
             "generico"),
            (1048576, ["R2", "S11", "S13", "S17"], "2^e"),
            (720720, ["R2", "S11", "S13", "S17"], "colapso")]
    results = []
    fam_points = {}
    for (N, tags, fam) in plan:
        pts = []
        for tag in tags:
            t0 = time.time()
            m = measure(supports[tag], N)
            results.append(("%s/%s" % (fam, tag), m))
            pts.append((m["nHN"], m["slack"], m["S_Cauchy"],
                        m["S_max"]))
            print(fmt_row("%s/%s" % (fam, tag), m)
                  + " (%.0fs)" % (time.time() - t0))
        fam_points[fam] = pts
        print()
    allok &= check_K1_K2_K3(results)
    print()
    print("Crecimiento de la holgura con |H_N| (regresion log-log):")
    theta = {}
    for fam in ("generico", "2^e"):
        theta[fam] = theta_fit(fam_points[fam], fam)
    print("     [colapso N=720720] slack frente a W_max (B_M=S_BV "
          "implica slack = W_max):")
    for (tag, m) in results:
        if m["N"] == 720720:
            print("       [%s] slack=%.2f W_max=%.2f cociente=%.6f"
                  % (tag, m["slack"], m["W_max"],
                     m["slack"] / m["W_max"]))
    print()
    okE2, rowsE2, resE2, locW, locH = e101_2()
    allok &= okE2
    print()
    print("RESUMEN E101.1: theta_slack(generico)=%.4f, "
          "theta_slack(2^e)=%.4f"
          % (theta["generico"]["slack"][0], theta["2^e"]["slack"][0]))
    print("RESUMEN E101.2: alpha_W(11..37)=%.4f IC95%% [%.4f, %.4f] | "
          "alpha_W(11..29)=%.4f"
          % (resE2[("W", "z=11..37 (extendida)")][0],
             resE2[("W", "z=11..37 (extendida)")][0]
             - resE2[("W", "z=11..37 (extendida)")][1],
             resE2[("W", "z=11..37 (extendida)")][0]
             + resE2[("W", "z=11..37 (extendida)")][1],
             resE2[("W", "z=11..29 (doc 101)")][0]))
    print("PASS" if allok else "FAIL (ver lineas marcadas)")
    print("Tiempo total: %.1f s" % (time.time() - t00))
    return 0 if allok else 1


if __name__ == "__main__":
    sys.exit(main())
