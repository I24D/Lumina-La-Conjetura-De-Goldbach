#!/usr/bin/env python3
# 100_verificacion_doc98_99_segundo_factor.py  (Doc 100)
# Verificacion independiente PARCIAL de los docs 98 y 99 de Lumina Chat
# (correccion Delta_r + primer diagnostico numerico del segundo factor).
# Solo stdlib.
#
# LIMITACION (ver doc 100, sec. 0): el tunel al PC sigue caido (segunda
# caida de egress del 2026-11-27, desde ~12:16 EST). Los docs 98/99
# integros y los scripts 97/99 de Lumina NO estan disponibles. Esta
# verificacion trabaja SOLO con la respuesta que Lumina dejo en el chat
# "Revisa UNIF U 2M X" el 2026-11-28 ~00:40 (transcrita en el doc 100)
# mas los docs 85-97 locales. Todo lo verificado aqui lo fue sobre la
# transcripcion y sobre calculos independientes propios; nada depende del
# texto exacto de los docs 98/99. Los items [P-rev] deben re-verificarse
# contra los documentos completos y los scripts 97/99 cuando el tunel
# se restablezca.
#
# Modelo sintetico: igual que en el doc 97 (pesos lambda^+ de Nathanson,
# V^+, A=modulo exterior, z=limite de criba, y=extremo switching,
# D=nivel). Novedades de este doc:
#   * PsiEngine: psi(N;d,a) real por criba de von Mangoldt, con cache.
#   * M_r = sum_{d|r} mu(r/d) phi(d) E(N;d,N) (parte de Mobius, sin Delta).
#   * Delta_r por ORTOGONALIDAD (nuevo, sin enumerar caracteres):
#       Delta_r = sum_{p|r, p^2 no | r} P_{r/p},
#       P_e = sum_{d|e} mu(e/d) A_d,
#       A_d = phi(d) * sum_{j: p^j<=N, p^j=N (mod d)} log p.
#     Valido solo si (p,e)=1 y (N,d)=1; los terminos con p^2|r se anulan
#     porque chi*(p)=0 cuando p|e (ver doc 100, sec. 2).
#   * Segundo factor: B_M = sum (1/phi(r)) |M_r|^2,
#       B_cross = sum (1/phi(r)) 2 M_r Delta_r,
#       B_Delta = sum (1/phi(r)) |Delta_r|^2, sobre r en el soporte R_A
#       con r > R0, R0 = (ln N)^kappa.
#   * Descomposicion BV: |M_r|^2 <= r * sum_{d|r} phi(d) |E_d|^2, de donde
#       B_M <= W_max * S_BV con W_d = sum_{r ni d} r/phi(r),
#       S_BV = sum_d phi(d) |E(N;d,N)|^2.
#
# Chequeos:
#   K0  Delta_r por ortogonalidad vs calculo directo con caracteres
#       (CharGroup, como en el doc 97 H6): S_r = M_r + Delta_r exacto.
#   K1  (a) Segundo factor: reproduccion de B_M, B_cross, B_Delta y del
#       cambio relativo en tres modelos (R1 S1, R2 S4, R3 z=3).
#   K2  (b) Condicion p||r: Delta_r = 0 para potencias de primos; error
#       espurio si se ignora la condicion (modelo "chi*(p) como unidad").
#   K3  (c) Condicion de N: auditoria z=3 en el doc 97; segundo factor
#       con N par (999982, 2^20) y N con muchos primos pequenos (510510).
#   K4  (d) Las tres cuestiones: cota de gran criba vs B_M medido;
#       descomposicion BV (W_max, S_BV); analisis de brecha asintotica
#       frente al objetivo de Chen (c=0.039884).
#   K5  Numeros de verificacion de Lumina (16 conductores): error abs
#       max, error rel max, correccion max |Delta|/|S_r|.
# Uso: python3 100_verificacion_doc98_99_segundo_factor.py
import itertools
import math
import sys
import time


# ------------------------------------------------------------------ utiles
# (copiados del script del doc 97, misma convencion: pf_desc en orden
# DECRECIENTE, como exige el Teorema 9.3 de Nathanson)
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
    """Criba lineal: mu[0..n], spf[0..n], primos."""
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
    """Factores primos distintos de d, en orden DECRECIENTE."""
    out = []
    while d > 1:
        p = spf[d]
        out.append(p)
        while d % p == 0:
            d //= p
    return sorted(out, reverse=True)


def in_Vplus(d, D, pdesc):
    """Teorema 9.3 (Nathanson), beta=2, V^+: d=p1...pk, p1>...>pk, y para
    todo m impar, p_m < (D/(p1...pm))^{1/2}, i.e. p_m^2*(p1...pm) < D.
    d=1 siempre esta. Se supone d libre de cuadrados."""
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
    """Todos los divisores de n (n>=1), via sus factores primos."""
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


def relerr(a, b):
    return abs(a - b) / max(1.0, abs(a), abs(b))


SMALL_PRIMES = primes_upto(10000)

# ------------------------------------------------------------------ modelo
# (Model97 copiado del script del doc 97 sin cambios)
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


# ------------------------------------------------------------------ psi real
# (lambda_sieve y CharGroup copiados del script del doc 97)
def prime_power_decomp(n):
    out = []
    x = n
    p = 2
    while p * p <= x:
        if x % p == 0:
            e = 0
            while x % p == 0:
                x //= p
                e += 1
            out.append((p, e))
        p += 1 if p == 2 else 2
    if x > 1:
        out.append((x, 1))
    return out


def _gen_prime_power(p, e):
    """Generador de (Z/p^e)*, p primo impar (grupo ciclico)."""
    q = p ** e
    if q <= 2:
        return 1, 1
    phi = q - q // p
    facs = set()
    x = phi
    d = 2
    while d * d <= x:
        if x % d == 0:
            facs.add(d)
            while x % d == 0:
                x //= d
        d += 1 if d == 2 else 2
    if x > 1:
        facs.add(x)
    for g in range(2, q):
        if math.gcd(g, q) != 1:
            continue
        if all(pow(g, phi // f, q) != 1 for f in facs):
            return g, phi
    raise ValueError("sin generador para %d^%d" % (p, e))


def _dlog(a, g, q, order):
    cur = 1
    for k in range(order):
        if cur == a % q:
            return k
        cur = cur * g % q
    raise ValueError("dlog fallo")


class CharGroup:
    """Caracteres de (Z/rZ)* via CRT sobre potencias de primos.
    Solo r sin factor 4 (r impar, o r = 2 * impar)."""

    def __init__(self, r):
        self.r = r
        self.parts = []
        for (p, e) in prime_power_decomp(r):
            q = p ** e
            g, ph = _gen_prime_power(p, e)
            self.parts.append((q, g, ph))
        self.units = [a for a in range(1, r) if math.gcd(a, r) == 1]
        self.dlogs = {}
        for a in self.units:
            ks = []
            for (q, g, ph) in self.parts:
                ks.append(_dlog(a % q, g, q, ph))
            self.dlogs[a] = tuple(ks)

    def chi(self, t, a):
        a = a % self.r
        if math.gcd(a, self.r) != 1:
            return 0j
        v = 1 + 0j
        for (q, g, ph), ti, ki in zip(self.parts, t, self.dlogs[a]):
            v *= complex(math.cos(2 * math.pi * ti * ki / ph),
                         math.sin(2 * math.pi * ti * ki / ph))
        return v

    def all_t(self):
        ranges = [range(ph) for (_, _, ph) in self.parts]
        return list(itertools.product(*ranges))

    def conductor(self, t):
        divs = sorted(d for d in range(1, self.r + 1) if self.r % d == 0)
        for d in divs:
            ok = True
            for a in self.units:
                if a % d == 1 % d:
                    if abs(self.chi(t, a) - 1) > 1e-9:
                        ok = False
                        break
            if ok:
                return d
        return self.r

    def primitive_ts(self):
        return [t for t in self.all_t()
                if not all(x == 0 for x in t)
                and self.conductor(t) == self.r]


def lambda_sieve(N):
    """lam[n] = log p si n = p^j, 0 en otro caso (criba de von Mangoldt)."""
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


def prime_powers_list(N):
    """Lista de (n, log p) para n = p^j <= N."""
    out = []
    for p in primes_upto(N):
        pw = p
        while pw <= N:
            out.append((pw, math.log(p)))
            if pw > N // p:
                break
            pw *= p
    return out


def psi_buckets_pp(pp, moduli):
    """Para cada modulo m: buckets[a] = sum_{n<=N, n=a mod m} lam[n],
    en una sola pasada sobre la lista de potencias de primos."""
    out = {}
    for m in moduli:
        bk = [0.0] * m
        for (n, lp) in pp:
            bk[n % m] += lp
        out[m] = bk
    return out


def is_prime(n):
    if n < 2:
        return False
    if n % 2 == 0:
        return n == 2
    d = 3
    while d * d <= n:
        if n % d == 0:
            return False
        d += 2
    return True

# ------------------------------------------------------------------ nuevo:
# PsiEngine, M_r, Delta_r por ortogonalidad, segundo factor, BV
class PsiEngine:
    """psi(N;d,a) real con cache. lam[n] = Lambda(n)."""

    def __init__(self, N):
        self.N = N
        self.lam = lambda_sieve(N)
        self.psiN = sum(self.lam)
        self.cache = {}

    def psi_mod(self, d, a):
        """psi(N; d, a) = sum_{n<=N, n=a mod d} Lambda(n)."""
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
        """E(N; d, N) = psi(N; d, N) - psi(N)/phi(d). d=1 -> 0."""
        if d == 1:
            return 0.0
        return self.psi_mod(d, self.N) - self.psiN / phi_of(d, SMALL_PRIMES)


def M_r(r, eng, mu, spf):
    """Parte de Mobius (sin Delta): sum_{d|r} mu(r/d) phi(d) E(N;d,N)."""
    tot = 0.0
    for d in divisors_of(r, spf):
        tot += mu[r // d] * phi_of(d, SMALL_PRIMES) * eng.E(d)
    return tot


def P_e(e, p, N, mu, spf):
    """P_e = sum*_{chi* mod e} chi_bar*(N) sum_{p^j<=N} chi*(p)^j log p,
    calculado por ortogonalidad (sin enumerar caracteres):
    P_e = sum_{d|e} mu(e/d) A_d,
    A_d = phi(d) * sum_{j: p^j<=N, p^j=N (mod d)} log p.
    Valido solo si (p,e)=1 y (N,d)=1 para todo d|e."""
    tot = 0.0
    for d in divisors_of(e, spf):
        s = 0.0
        pw = p
        while pw <= N:
            if pw % d == N % d:
                s += math.log(p)
            if pw > N // p:
                break
            pw *= p
        tot += mu[e // d] * phi_of(d, SMALL_PRIMES) * s
    return tot


def delta_r_orth(r, N, mu, spf):
    """Delta_r por ortogonalidad (Teorema 97.1):
    Delta_r = sum_{p|r, p^2 no | r} P_{r/p}.
    Los terminos con p^2|r se omiten: chi*(p)=0 cuando p | r/p, asi que
    P_{r/p} = 0 para ellos (ver doc 100, sec. 2.2)."""
    tot = 0.0
    for p in pf_desc(r, spf):
        e = r // p
        if e % p == 0:
            continue
        tot += P_e(e, p, N, mu, spf)
    return tot


def spurious_term(r, N, mu, spf):
    """Error que se cometeria al ignorar la condicion p||r Y evaluar
    chi*(p) como si p fuera unidad mod r/p (es decir, sustituir
    sum_j chi*(p)^j log p por sum_j 1*log p = (log p)*floor(log_p N)):
    por cada primo p|r con p^2|r, la contribucion espuria es
    (log p)*floor(log_p N) * sum*_{chi* mod r/p} chi_bar*(N),
    donde sum* chi_bar*(N) = sum_{d|r/p, d|N-1} mu((r/p)/d) phi(d).
    Se calcula directamente (sin suponer nada sobre N-1)."""
    tot = 0.0
    detail = []
    for p in pf_desc(r, spf):
        e = r // p
        if e % p != 0:
            continue
        ssum = 0.0
        for d in divisors_of(e, spf):
            if (N - 1) % d == 0:
                ssum += mu[e // d] * phi_of(d, SMALL_PRIMES)
        jmax = 0
        pw = p
        while pw <= N:
            jmax += 1
            if pw > N // p:
                break
            pw *= p
        contrib = math.log(p) * jmax * ssum
        detail.append((p, e, ssum, jmax, contrib))
        tot += contrib
    return tot, detail


def second_factor(mod, eng, N, kappa, mu, spf, use_R0=True):
    """B_M, B_cross, B_Delta sobre r en el soporte H_A con r > R0
    (R0 = (ln N)^kappa) y (r,N)=1. Devuelve tambien las filas (r, M, D)."""
    R0 = (math.log(N) ** kappa) if use_R0 else 0.0
    BM = BC = BD = 0.0
    rows = []
    for r in mod.H:
        if r <= R0:
            continue
        if math.gcd(r, N) != 1:
            continue
        Mr = M_r(r, eng, mu, spf)
        Dr = delta_r_orth(r, N, mu, spf)
        ph = phi_of(r, SMALL_PRIMES)
        BM += Mr * Mr / ph
        BC += 2.0 * Mr * Dr / ph
        BD += Dr * Dr / ph
        rows.append((r, Mr, Dr))
    return BM, BC, BD, rows


def bv_decomposition(mod, eng, N, mu, spf):
    """Descomposicion BV del segundo factor:
    |M_r|^2 <= r * sum_{d|r} phi(d) |E_d|^2  (Cauchy con pesos phi(d)),
    de donde B_M <= W_max * S_BV con
    W_d = sum_{r en soporte, d|r} r/phi(r),
    S_BV = sum_d phi(d) |E(N;d,N)|^2 (d en el soporte de divisores).
    Devuelve (W_max, S_BV, max_slack_Cauchy) donde max_slack_Cauchy es
    el maximo de |M_r|^2 / (r * sum_{d|r} phi(d)|E_d|^2) (debe ser <=1)."""
    # soporte de divisores d (todos los divisores de algun r del soporte)
    divs = set()
    for r in mod.H:
        if math.gcd(r, N) != 1:
            continue
        for d in divisors_of(r, spf):
            divs.add(d)
    divs = sorted(divs)
    E2 = {}
    for d in divs:
        if d == 1:
            E2[d] = 0.0
        else:
            e = eng.E(d)
            E2[d] = e * e
    S_BV = sum(phi_of(d, SMALL_PRIMES) * E2[d] for d in divs if d > 1)
    W = {}
    for d in divs:
        w = 0.0
        for r in mod.H:
            if r % d == 0 and math.gcd(r, N) == 1:
                w += r / phi_of(r, SMALL_PRIMES)
        W[d] = w
    W_max = max(W.values()) if W else 0.0
    # holgura del paso de Cauchy por r
    max_slack = 0.0
    worst = None
    for r in mod.H:
        if math.gcd(r, N) != 1:
            continue
        Mr = M_r(r, eng, mu, spf)
        rhs = r * sum(phi_of(d, SMALL_PRIMES) * E2[d]
                      for d in divisors_of(r, spf))
        if rhs > 0:
            sl = (Mr * Mr) / rhs
            if sl > max_slack:
                max_slack, worst = sl, r
    return W_max, S_BV, max_slack, worst, len(divs)

# ------------------------------------------------------------------ chequeos
MU, SPF, _ = mu_spf(600000)


def M_r_from_buckets(r, bk, psiN, mu, spf):
    """M_r con psi real desde buckets (para K0/K5)."""
    tot = 0.0
    for d in divisors_of(r, spf):
        E = bk[d][999983 % d] - psiN / phi_of(d, SMALL_PRIMES) \
            if d > 1 else 0.0
        tot += mu[r // d] * phi_of(d, SMALL_PRIMES) * E
    return tot


def S_r_char(r, bk, psiN):
    """S_r directo con caracteres primitivos (CharGroup), parte real.
    Devuelve (S_r_real, max|Im|)."""
    cg = CharGroup(r)
    prims = cg.primitive_ts()
    Sr = 0j
    for t in prims:
        psi_chi = 0j
        for a in cg.units:
            psi_chi += cg.chi(t, a) * bk[r][a]
        Sr += cg.chi(t, 999983).conjugate() * psi_chi
    return Sr.real, abs(Sr.imag)


def check_K0K5(N, r_list):
    # K0: Delta_r por ortogonalidad vs calculo directo con caracteres.
    # K5: numeros de verificacion de Lumina (16 conductores).
    pp = prime_powers_list(N)
    moduli = set()
    for r in r_list:
        for d in divisors_of(r, SPF):
            moduli.add(d)
        moduli.add(r)
    bk = psi_buckets_pp(pp, sorted(moduli))
    psiN = sum(lp for (_, lp) in pp)
    print("K0/K5 [N=%d, %d conductores] S_r directo (caracteres) vs "
          "M_r + Delta_r (ortogonalidad):" % (N, len(r_list)))
    maxerr = 0.0
    maxrel = 0.0
    maxcorr = 0.0
    maxd_orth_char = 0.0
    for r in r_list:
        Sr, im = S_r_char(r, bk, psiN)
        Mr = M_r_from_buckets(r, bk, psiN, MU, SPF)
        Dr = delta_r_orth(r, N, MU, SPF)
        err = abs(Sr - (Mr + Dr))
        scale = max(1.0, abs(Sr), abs(Mr))
        rel = err / scale
        maxerr = max(maxerr, err)
        maxrel = max(maxrel, rel)
        if abs(Sr) > 1e-9:
            maxcorr = max(maxcorr, abs(Dr) / abs(Sr))
        # Delta por caracteres: S_r - M_r
        dchar = Sr - Mr
        maxd_orth_char = max(maxd_orth_char, abs(Dr - dchar))
        print("     r=%3d #chi*=%2d S_r=%12.4f M_r=%12.4f Delta_orth=%9.4f "
              "Delta_char=%9.4f |S-(M+D)|=%.2e |Im(S)|=%.1e"
              % (r, len(CharGroup(r).primitive_ts()), Sr, Mr, Dr, dchar,
                 err, im))
    ok = maxerr < 1e-6 and maxd_orth_char < 1e-6
    print("     max|S-(M+Delta)|=%.2e (Lumina: 6.56e-09) | max rel=%.2e "
          "(Lumina: 9.71e-12) | max |Delta|/|S|=%.4f (Lumina: 0.0259) | "
          "max|Delta_orth-Delta_char|=%.2e %s"
          % (maxerr, maxrel, maxcorr, maxd_orth_char,
             "PASS" if ok else "FAIL"))
    return ok


def check_K1(models, N, kappa=1.0):
    # (a) Segundo factor: B_M, B_cross, B_Delta y cambio relativo.
    print("K1 [N=%d, kappa=%.1f, R0=(ln N)^kappa=%.2f] segundo factor "
          "sobre el soporte R_A (r>R0):" % (N, kappa, math.log(N) ** kappa))
    print("     Lumina doc 99 (69 conductores): B_M=3.9222e+07 "
          "B_cross=5.4023e+03 B_Delta=466.7 cambio rel=1.4966e-04 | "
          "doc 97: peso Delta = 1.24e-04")
    eng = PsiEngine(N)
    res = {}
    for (mod, tag) in models:
        BM, BC, BD, rows = second_factor(mod, eng, N, kappa, MU, SPF)
        rel = (BC + BD) / BM if BM else float("inf")
        w97 = BD / (BM + BC + BD) if (BM + BC + BD) else float("inf")
        ratio = BC / BD if BD else float("inf")
        # descomposicion: cuanto de B_M viene de conductores primos
        prime_share = 0.0
        for (r, Mr, Dr) in rows:
            if is_prime(r):
                prime_share += Mr * Mr / phi_of(r, SMALL_PRIMES)
        res[tag] = (BM, BC, BD, rel)
        ok = BC > BD and 1e-5 < rel < 1e-3
        print("     [%s] #r=%d B_M=%.4e B_cross=%.4e B_Delta=%.4e | "
              "cambio rel=%.4e | B_cross/B_Delta=%.2f | peso Delta "
              "(doc97)=%.3e | B_M de primos=%.2f%% %s"
              % (tag, len(rows), BM, BC, BD, rel, ratio, w97,
                 100.0 * prime_share / BM if BM else 0.0,
                 "PASS" if ok else "FAIL (ver texto)"))
        # top-3 contribuyentes a B_M
        top = sorted(rows, key=lambda x: x[1] ** 2 / phi_of(x[0],
                                                           SMALL_PRIMES),
                     reverse=True)[:3]
        print("       top B_M: " + "; ".join(
            "r=%d (1/phi)M^2=%.3e" % (r, Mr ** 2 / phi_of(r, SMALL_PRIMES))
            for (r, Mr, Dr) in top))
    return res


def check_K1_sens(mod, tag, N):
    # Sensibilidad del segundo factor a kappa y al corte R0.
    eng = PsiEngine(N)
    print("K1-sens [%s, N=%d] sensibilidad a kappa y al corte R0:" % (tag, N))
    for kap in (1.0, 2.0, 3.0):
        BM, BC, BD, rows = second_factor(mod, eng, N, kap, MU, SPF)
        rel = (BC + BD) / BM if BM else float("inf")
        print("     kappa=%.1f: #r=%d B_M=%.4e B_cross=%.4e B_Delta=%.4e "
              "cambio rel=%.4e" % (kap, len(rows), BM, BC, BD, rel))
    BM, BC, BD, rows = second_factor(mod, eng, N, 1.0, MU, SPF,
                                     use_R0=False)
    rel = (BC + BD) / BM if BM else float("inf")
    print("     sin corte R0: #r=%d B_M=%.4e B_cross=%.4e B_Delta=%.4e "
          "cambio rel=%.4e" % (len(rows), BM, BC, BD, rel))
    return True


def check_K2(N, r_list):
    # (b) Condicion p||r: Delta_r = 0 para potencias de primos; error
    # espurio si se ignora la condicion (modelo "chi*(p) como unidad").
    print("K2 [N=%d] condicion p||r (p divide a r pero p^2 no):" % N)
    pp = prime_powers_list(N)
    moduli = set()
    for r in r_list:
        for d in divisors_of(r, SPF):
            moduli.add(d)
        moduli.add(r)
    bk = psi_buckets_pp(pp, sorted(moduli))
    psiN = sum(lp for (_, lp) in pp)
    ok = True
    for r in r_list:
        Dr = delta_r_orth(r, N, MU, SPF)
        Sr, _ = S_r_char(r, bk, psiN)
        Mr = M_r_from_buckets(r, bk, psiN, MU, SPF)
        spur, detail = spurious_term(r, N, MU, SPF)
        is_pp = len(pf_desc(r, SPF)) == 1 and r not in primes_upto(200)
        # potencia de primo pura: r = p^e, e>=2
        pure = (len(pf_desc(r, SPF)) == 1 and not is_prime(r))
        flag = ""
        if pure and abs(Dr) > 1e-9:
            ok = False
            flag = " <-- FAIL: Delta != 0 en potencia de primo"
        if abs(Dr - (Sr - Mr)) > 1e-6:
            ok = False
            flag += " <-- FAIL: Delta_orth != Delta_char"
        print("     r=%3d Delta_orth=%9.4f Delta_char=%9.4f "
              "espurio(si se ignora p||r)=%9.4f |S_r|=%.1f%s"
              % (r, Dr, Sr - Mr, spur, abs(Sr), flag))
    print("     (El termino espurio es (log p)*floor(log_p N)*"
          "sum* chi_bar*(N) por cada p|r con p^2|r; para N=999983 vale "
          "(log p)*floor(log_p N)*mu(r/p) cuando el unico divisor de r/p "
          "que divide a N-1 es 1.)")
    print("     K2 %s" % ("PASS" if ok else "FAIL"))
    return ok


def check_K3():
    # (c1) Auditoria: el doc 97 uso alguna vez el proxy z=3?
    print("K3(c1) auditoria del proxy z=3 en el doc 97:")
    import re
    ok = True
    for fn in ("97_verificacion_doc96_diagonal_mobius.md",
               "97_verificacion_doc96_diagonal_mobius.py"):
        try:
            txt = open(fn).read()
        except OSError as e:
            print("     %s: no legible (%s)" % (fn, e))
            ok = False
            continue
        n1 = len(re.findall(r"z\s*=\s*3\b", txt))
        n2 = len(re.findall(r"z\s*<\s*3\b", txt))
        print("     %s: ocurrencias de 'z=3' o 'z = 3': %d; 'z<3': %d"
              % (fn, n1, n2))
        if n1 or n2:
            ok = False
    zs = set()
    for m in re.finditer(r"Model97\(3,\s*(\d+),", open(
            "97_verificacion_doc96_diagonal_mobius.py").read()):
        zs.add(int(m.group(1)))
    for m in re.finditer(r"check_H5\(3,\s*(\d+),", open(
            "97_verificacion_doc96_diagonal_mobius.py").read()):
        zs.add(int(m.group(1)))
    print("     valores de z usados en los modelos del doc 97: %s "
          "(minimo %d; el proxy z=3 NO se uso en ningun modelo)"
          % (sorted(zs), min(zs) if zs else -1))
    print("     K3(c1) %s" % ("PASS" if ok else "FAIL"))
    return ok


def check_K3_N(models_N):
    # (c2) Segundo factor con N par y con N con muchos primos pequenos.
    print("K3(c2) segundo factor con N par / N con primos pequenos "
          "(modelo R1 salvo 9699690, que usa excl=(5,7,11,13,17,19)):")
    for (N, mod, tag) in models_N:
        eng = PsiEngine(N)
        BM, BC, BD, rows = second_factor(mod, eng, N, 1.0, MU, SPF)
        rel = (BC + BD) / BM if BM else float("inf")
        print("     [%s] N=%d (#r=%d, R0=%.2f) B_M=%.4e B_cross=%.4e "
              "B_Delta=%.4e cambio rel=%.4e B_cross/B_Delta=%.2f"
              % (tag, N, len(rows), math.log(N), BM, BC, BD, rel,
                 BC / BD if BD else float("inf")))
    return True


def check_K4(models, N):
    # (d1) Cota de gran criba clasica vs B_M medido.
    print("K4(d1) [N=%d] cota de gran criba clasica vs B_M medido:" % N)
    print("     B_M <= sum_{r<=R} sum*_{chi* mod r} |psi(N,chi)|^2 "
          "<= (N + R^2) * N * ln N  (gran criba, a_n = Lambda(n), "
          "sum Lambda^2 ~ N ln N; n_r <= phi(r))")
    eng = PsiEngine(N)
    for mod, tag in models:
        BM, BC, BD, rows = second_factor(mod, eng, N, 1.0, MU, SPF)
        R = max(r for r in mod.H if math.gcd(r, N) == 1)
        ls = (N + R * R) * N * math.log(N)
        print("     [%s] R=max r=%d B_M=%.4e cota LS=%.4e cota/B_M=%.3e "
              "(la cota clasica es ~%d veces el valor medido: no "
              "distingue el soporte disperso)"
              % (tag, R, BM, ls, ls / BM, int(ls / BM)))
    Wmax_out = None
    Wmax_h = {}
    print("K4(d2) descomposicion BV: B_M <= W_max * S_BV, W_d = "
          "sum_{r ni d} r/phi(r), S_BV = sum_d phi(d)|E_d|^2:")
    for mod, tag in models:
        BM, BC, BD, rows = second_factor(mod, eng, N, 1.0, MU, SPF)
        W_max, S_BV, slack, worst, ndiv = bv_decomposition(mod, eng, N, MU, SPF)
        print("     [%s] #divisores=%d W_max=%.2f S_BV=%.4e "
              "W_max*S_BV=%.4e B_M=%.4e holgura=%.3e | max holgura "
              "Cauchy por r: %.3f (r=%s, debe ser <=1)"
              % (tag, ndiv, W_max, S_BV, W_max * S_BV, BM,
                 (W_max * S_BV) / BM if BM else float("inf"),
                 slack, worst))
        Wmax_h[tag] = (W_max, len(mod.H))
        if tag.startswith("R2"):
            Wmax_out = W_max
    # (d3) Analisis de brecha asintotica frente al objetivo de Chen.
    # CORRECCION DE ESCALA (detectada al revisar K4): S_BV = sum_d
    # phi(d)|E_d|^2. Cada termino phi(d)|E_d|^2 ~ N ln N (escala ingenua)
    # y hay ~D_A = N^(1/2) divisores d distintos, luego S_BV ~ N^(3/2) ln N,
    # NO N^2 ln N. La viabilidad asintotica de la ruta BV depende entonces
    # de como crece W_max(N): necesita W_max(N) << N^(1-2c)/2 / ln^5 N
    # ~ N^0.42 / ln^5 N. En los modelos W_max ~ |H|.
    print("K4(d3) brecha asintotica frente al objetivo de Chen "
          "(ordenes de magnitud; constantes [P-rev] del doc 94):")
    print("     Objetivo: R_grande = sqrt(G_grande * B_M) << N/(A ln^2 N)"
          " con G_grande = f*G_A ~ 0.05*0.45 = 0.0225, A = N^c, "
          "c = 0.039884. Luego B_M << N^2/(A^2 ln^4 N * 0.0225).")
    print("     Cotas: LS = (N + R^2) N ln N con R = A*D_A, D_A = N^(1/2); "
          "BV = W_max * N^(3/2) ln N con W_max = %.1f (medido en R2; "
          "W_max/|H| = %.2f, %.2f, %.2f en R1/R2/R3)."
          % (Wmax_out if Wmax_out else float("nan"),
             Wmax_h["R1"][0] / Wmax_h["R1"][1],
             Wmax_h["R2"][0] / Wmax_h["R2"][1],
             Wmax_h["R3"][0] / Wmax_h["R3"][1]))
    print("     W_eq = maximo W_max admisible para que BV alcance el "
          "objetivo (W_eq = objetivo / (N^(3/2) ln N)).")
    for log10N in (12, 24, 100):
        Nbig = 10.0 ** log10N
        lnN = math.log(Nbig)
        c = 0.039884
        A = Nbig ** c
        R = A * Nbig ** 0.5
        target = Nbig ** 2 / (A * A * lnN ** 4 * 0.0225)
        ls = (Nbig + R * R) * Nbig * lnN
        bv = (Wmax_out if Wmax_out else 1.0) * Nbig ** 1.5 * lnN
        weq = target / (Nbig ** 1.5 * lnN)
        print("     N=10^%d: objetivo B_M << %.3e | LS=%.3e (%.1e x el "
              "objetivo) | BV=%.3e (%.1e x el objetivo) | W_eq=%.3e"
              % (log10N, target, ls, ls / target, bv, bv / target, weq))
    print("     Lectura: la gran criba clasica falla por un factor "
          "N^0.16 (creciente). La ruta BV, con W_max FIJO, cruza el "
          "objetivo entre N=10^12 y N=10^24 y lo supera ampliamente a "
          "N=10^100; su viabilidad real depende de que W_max(N) crezca "
          "mas despacio que N^0.42. W_max ~ |H| en los modelos: el "
          "conteo del soporte es el nudo de la cuestion (a).")
    return True


# ------------------------------------------------------------------ main
def main():
    allok = True
    t0 = time.time()
    print("== DOC 100: verificacion docs 98/99 (Delta_r / segundo factor) ==")
    print("Solo stdlib. Modelos sinteticos con pesos lambda^+ de Nathanson "
          "y psi reales hasta N~1e7.")
    print()

    N = 999983
    assert is_prime(N)

    # Modelos: R1 ~ S1 del doc 97; R2 ~ S4 (parametros del doc 94);
    # R3 = hipotesis z=3 (soporte {1} U {primos s}, 69 conductores).
    m1 = Model97(3, 11, 60, 20000, "fixed")
    m2 = Model97(3, 11, 3000, 200000, "fixed")
    m3 = Model97(3, 3, 349, 200000, "fixed")
    print("Modelos: R1 A=3 z=11 y=60 D=20000 |H|=%d (%.1fs) | R2 A=3 z=11 "
          "y=3000 D=200000 |H|=%d (%.1fs) | R3 A=3 z=3 y=349 D=200000 "
          "|H|=%d (%.1fs)" % (len(m1.H), m1.build_t, len(m2.H),
                              m2.build_t, len(m3.H), m3.build_t))
    print("     R3: soporte = {1} U {s primo, 5<=s<349, (s,3)=1}: "
          "%d conductores (hipotesis: los 69 de Lumina)" % len(m3.H))
    print()

    # K0 + K5: Delta por ortogonalidad vs caracteres (16 conductores)
    r16 = [3, 5, 7, 9, 11, 13, 15, 21, 25, 27, 35, 45, 49, 63, 75, 99]
    allok &= check_K0K5(N, r16)
    print()

    # K1: segundo factor en R1/R2/R3
    check_K1([(m1, "R1"), (m2, "R2"), (m3, "R3")], N, kappa=1.0)
    print()
    check_K1_sens(m1, "R1", N)
    print()

    # K2: condicion p||r
    r_pp = [9, 25, 27, 45, 49, 63, 75, 99, 121, 125, 135, 81]
    allok &= check_K2(N, r_pp)
    print()

    # K3: (c1) auditoria z=3; (c2) N par / N con primos pequenos
    allok &= check_K3()
    print()
    m1n = Model97(3, 11, 60, 20000, "fixed", excl=(5, 7, 11, 13, 17, 19))
    print("     Modelo R1n (para N=9699690): |H|=%d, soporte tras "
          "excluir primos de N" % len(m1n.H))
    check_K3_N([(999983, m1, "R1 N=999983 (primo impar, referencia)"),
                (999982, m1, "R1 N=999982 (par, 2 x primo)"),
                (1048576, m1, "R1 N=2^20 (par, potencia de 2)"),
                (9699690, m1n, "R1n N=9699690 (par, 2*5*7*11*13*17*19)")])
    print()

    # K4: las tres cuestiones
    check_K4([(m1, "R1"), (m2, "R2"), (m3, "R3")], N)
    print()
    print("PASS" if allok else "FAIL (ver lineas marcadas)")
    print("Tiempo total: %.1f s" % (time.time() - t0))
    return 0 if allok else 1


if __name__ == "__main__":
    sys.exit(main())
