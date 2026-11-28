#!/usr/bin/env python3
# 97_verificacion_doc96_diagonal_mobius.py  (Doc 97)
# Verificacion independiente PARCIAL del doc 96 de Lumina Chat (ruta
# diagonal + identidad de Mobius + experimentos 94d/94e). Solo stdlib.
#
# LIMITACION (ver doc 97, sec. 0): el doc 96 integro y los scripts 94d/94e
# viven en OneDrive en el PC y NO estan disponibles (tunel caido, segunda
# caida de egress del 2026-11-27, desde ~12:16 EST). Esta verificacion
# trabaja SOLO con la respuesta que Lumina dejo en el chat "Revisa UNIF U
# 2M X" el 2026-11-27 (transcrita en el doc 97) mas los docs 85-95
# locales. Todo lo verificado aqui lo fue sobre la transcripcion y sobre
# calculos independientes propios; nada depende del texto exacto del
# doc 96. Los items [P-rev] deben re-verificarse contra el doc 96
# completo y los scripts 94d/94e cuando el tunel se restablezca.
#
# Modelo sintetico. Parametros (A, z, y, D, modo, excl):
#   A    = modulo exterior (impar, (A,N)=1; juega el papel de Qu).
#   z    = limite inferior de criba / extremo inferior del rango switching.
#   y    = extremo superior del rango switching.
#   D    = nivel (D_A); proxy D = X_A^{1/2} sin log^C, como en docs 85-88.
#   modo = 'fixed'    : d | P(z)  (limite FIJO z, version certificada doc 87)
#          'variable' : d | P(s)  (limite VARIABLE s, version doc 85)
#   excl = primos extra excluidos del soporte (ademas de 2 y de los
#          primos que dividen a A): modelan los primos que dividen a N.
# Para cada s primo de switching en [z,y), (s,A)=1, s no en excl:
#   V_s = {d : d < D/s, mu(d)!=0, (d,A)=1, d impar, primos(d) < Plim(s),
#          d en V^+(D/s), primos(d) no en excl} con Plim(s)=z (fixed) o s
#          (variable).
#   b(m) = sum_{(s,d): sd=m} mu(d)   (peso de Nathanson).
#   T_s(h0) = sum_{d en V_s, h0|d} mu(d)/phi(d).
#   Q_s = sum_{h0} phi(h0) |T_s(h0)|^2  (h0 suaves, primos < z).
#   C_s(r) = sum_{d en V_s, r|d} mu(d)/phi(d)   (nueva: doc 96).
#   Q_s^diag = sum_r phi(r) C_s(r)^2   (representacion diagonal).
#   G_A = sum_{h en H_A} phi(h) |K_A(h)|^2.
#
# Chequeos:
#   H0  hipotesis estructurales y conteos.
#   H1  identidad diagonal: Q_s directo vs sum_r phi(r) C_s(r)^2 (exacta);
#       y los terminos d=e solos NO bastan (diag_de / Q_s ~ 2-4).
#   H2  producto completo: (a) identidad exacta en el caso sin criba
#       (Q_s sin restriccion V^+ == Q_completo(z)); (b) cociente
#       Q_s / Q_completo con la restriccion V^+ (medicion, orden 1);
#       (c) Q_completo(z) ~ cte / ln z (orden de Mertens).
#   H3  perfil kappa (94e): fraccion de G_A con h <= R0, R0=(log N)^kappa.
#   H4  mejora del primer factor de Cauchy al restringir a r>R0 (punto 4).
#   H5  precaucion de N (punto 5): variante E0 (N=2^e) vs E1 (N con muchos
#       primos pequenos); cociente Q_s(E1)/Q_s(E0) vs factor producto.
#   H6  identidad de Mobius con psi reales: S_r vs sum mu(r/d)phi(d)E,
#       correccion imprimitiva Delta_r (formula derivada en el doc).
#   H7  perfil de masa en r de la representacion diagonal (cola en r).
# Uso: python3 97_verificacion_doc96_diagonal_mobius.py
import itertools
import math
import sys
import time


# ------------------------------------------------------------------ utiles
# (copiados del script del doc 95, con la misma convencion: pf_desc en
# orden DECRECIENTE, como exige el Teorema 9.3 de Nathanson)
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


def relerr(a, b):
    return abs(a - b) / max(1.0, abs(a), abs(b))


# ------------------------------------------------------------------ modelo
class Model97:
    def __init__(self, A, z, y, D, mode="fixed", excl=(), full=True,
                 s_list=None):
        t0 = time.time()
        self.A, self.z, self.y, self.D, self.mode = A, z, y, D, mode
        self.excl = set(excl)
        self.full = full
        self.phiA = phi_of(A, primes_upto(10000))
        allpr = primes_upto(max(y, 100))
        self.sprimes = [p for p in allpr
                        if p >= z and p < y and math.gcd(p, A) == 1
                        and p not in self.excl]
        if s_list is not None:
            self.sprimes = [s for s in self.sprimes if s in set(s_list)]
        self.small_primes = primes_upto(10000)
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
        # --- pares (s,d), b(m) ---
        self.b = {}
        self.pairs_by_m = {}
        for s in self.sprimes:
            for d in self.Vs[s]:
                m = s * d
                self.b[m] = self.b.get(m, 0) + mu[d]
                self.pairs_by_m.setdefault(m, []).append((s, d))
        self.M = sorted(self.b)
        # --- H_A: divisores libres de cuadrados de algun m ---
        H = {1}
        for m in self.M:
            pf = pf_desc(m, spf)
            for h in divisors_from_pf(pf):
                H.add(h)
        self.H = sorted(H)
        # --- K_A(h) por definicion directa ---
        self.K = {}
        for h in self.H:
            ssum = 0.0
            for m in self.M:
                if m % h == 0:
                    ssum += self.b[m] / phi_of(m, self.small_primes)
            self.K[h] = ssum / self.phiA
        # --- h0 suaves (primos < z) ---
        self.Hsmooth = []
        for h in self.H:
            pf = pf_desc(h, spf) if h > 1 else []
            if all(p < z for p in pf):
                self.Hsmooth.append(h)
        # --- T_s(h0) y Q_s (suave) ---
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

    # --- cantidades derivadas ---
    def G_A(self, Hset=None):
        Hset = self.H if Hset is None else Hset
        return sum(phi_of(h, self.small_primes) * self.K[h] ** 2
                   for h in Hset)

    def T_of(self, s, h0):
        return sum(self.mu[d] / phi_of(d, self.small_primes)
                   for d in self.Vs[s] if d % h0 == 0)

    def C_of(self, s, r):
        """C_s(r) = sum_{d en V_s, r|d} mu(d)/phi(d)  (doc 96)."""
        return sum(self.mu[d] / phi_of(d, self.small_primes)
                   for d in self.Vs[s] if d % r == 0)

    def R_of(self, s):
        """Todos los divisores de algun d en V_s (el rango de r)."""
        R = {1}
        for d in self.Vs[s]:
            for h in divisors_from_pf(pf_desc(d, self.spf)):
                R.add(h)
        return sorted(R)

    def Q_diag(self, s):
        """sum_r phi(r) C_s(r)^2 sobre todos los divisores de V_s."""
        tot = 0.0
        for r in self.R_of(s):
            c = self.C_of(s, r)
            tot += phi_of(r, self.small_primes) * c * c
        return tot

    def Q_of_h0set(self, s, h0set):
        return sum(phi_of(h0, self.small_primes) * self.T_of(s, h0) ** 2
                   for h0 in h0set)

    def diag_de(self, s):
        """Terminos d=e de la doble suma original: sum_d d/phi(d)^2."""
        return sum(d / phi_of(d, self.small_primes) ** 2
                   for d in self.Vs[s])

    def smooth_part(self, h):
        """Descomposicion de h en (s, h0): s = unico primo en [z,y) que
        divide a h (o None), h0 = resto. Devuelve (s, h0) o (None, h)."""
        pf = pf_desc(h, self.spf) if h > 1 else []
        sw = [p for p in pf if self.z <= p < self.y]
        if not sw:
            return (None, h)
        if len(sw) == 1:
            s = sw[0]
            return (s, h // s)
        return ("multi", h)


def Q_completo(z, A, excl=()):
    """Producto completo: prod_{p<z, p no | AN} (1 - 1/(p-1) + 1/(p-1)^2),
    con N modelado por excl (primos que dividen a N) y p=2 siempre
    excluido (N par). Igual a prod (1 - 2/(p-1) + p/(p-1)^2) (doc 96)."""
    tot = 1.0
    for p in primes_upto(z - 1):
        if p == 2 or A % p == 0 or p in excl:
            continue
        tot *= (1.0 - 1.0 / (p - 1) + 1.0 / (p - 1) ** 2)
    return tot


def Q_unrestricted(A, z, excl=()):
    """Q_s SIN restriccion V^+ ni cota D/s: d recorre TODOS los libres de
    cuadrados con primos < z, (d,AN)=1. Debe coincidir exactamente con
    Q_completo(z) (identidad del producto completo, doc 97 sec. 2.2)."""
    allow = [p for p in primes_upto(z - 1)
             if p > 2 and A % p != 0 and p not in excl]
    divs = [1]
    for p in allow:
        divs += [d * p for d in list(divs)]
    mu, spf, _ = mu_spf(max(divs) if divs else 1)
    small = primes_upto(10000)
    tot = 0.0
    for r in divs:
        c = 0.0
        for d in divs:
            if d % r == 0:
                c += mu[d] / phi_of(d, small)
        tot += phi_of(r, small) * c * c
    return tot


# ------------------------------------------------------------------ H6:
# caracteres primitivos y psi reales (identidad de Mobius)
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
    """Generador de (Z/p^e)*, p primo impar (grupo ciclico).
    Para (p,e)=(2,1) devuelve (1,1) (grupo trivial)."""
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
        self.parts = []  # (q, g, phi_q) por potencia de primo
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
        """t: tupla de exponentes (uno por parte); a: residuo mod r."""
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
        """Minimo d|r tal que chi(a)=1 para todo a unidad con a=1 mod d."""
        divs = sorted(d for d in range(1, self.r + 1) if self.r % d == 0)
        for d in divs:
            ok = True
            for a in self.units:
                if a % d == 1 % d:  # a == 1 (mod d); para d=1 es a%1==0
                    if abs(self.chi(t, a) - 1) > 1e-9:
                        ok = False
                        break
            if ok:
                return d
        return self.r

    def primitive_ts(self):
        # Caracteres primitivos NO principales (conductor == r).
        # El caracter principal (t todo cero) tiene conductor 1: para
        # r=1 no hay ningun caracter primitivo no principal.
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


def psi_buckets(lam, N, moduli):
    """Para cada modulo m: buckets[a] = sum_{n<=N, n=a mod m} lam[n]."""
    out = {}
    for m in moduli:
        bk = [0.0] * m
        for n in range(2, N + 1):
            if lam[n]:
                bk[n % m] += lam[n]
        out[m] = bk
    return out


# ------------------------------------------------------------------ H0
def check_H0(mod, tag):
    ok = True
    msgs = []
    for m in mod.M:
        if math.gcd(m, mod.A) != 1 or m % 2 == 0 or mod.mu[m] == 0:
            ok = False
            msgs.append("m=%d viola hipotesis" % m)
            break
    maxmult = max(len(v) for v in mod.pairs_by_m.values())
    npairs = sum(len(v) for v in mod.pairs_by_m.values())
    if mod.mode == "fixed" and maxmult != 1:
        ok = False
        msgs.append("multiplicidad fixed=%d" % maxmult)
    print("H0 [%s] |M|_pares=%d #m_distintos=%d max_mult=%d |H|=%d "
          "build=%.1fs %s" % (tag, npairs, len(mod.M), maxmult, len(mod.H),
                              mod.build_t, "PASS" if ok else "FAIL " +
                              ";".join(msgs)))
    return ok


# ------------------------------------------------------------------ H1
def check_H1(mod, tag):
    # Q_s directo (suave) vs representacion diagonal sum_r phi(r)C_s(r)^2.
    # En modo fixed ambas coinciden; en modo variable la diagonal necesita
    # el conjunto COMPLETO de h0 (todos los divisores de V_s).
    maxe = 0.0
    worst = None
    maxe_full = 0.0
    ratios = []
    for s in mod.sprimes:
        qd = mod.Q_diag(s)
        q = mod.Q[s]
        e = abs(qd - q)
        if e > maxe:
            maxe, worst = e, s
        qfull = mod.Q_of_h0set(s, mod.R_of(s))
        ef = abs(qd - qfull)
        if ef > maxe_full:
            maxe_full = ef
        dd = mod.diag_de(s)
        ratios.append(dd / q if q else float("inf"))
    ok = maxe < 1e-9 if mod.mode == "fixed" else maxe_full < 1e-9
    print("H1 [%s] diagonal: max|Q_diag - Q_s(suave)|=%.2e (s=%s) | "
          "max|Q_diag - Q_full|=%.2e | diag(d=e)/Q_s en [%.2f, %.2f] %s"
          % (tag, maxe, worst, maxe_full, min(ratios), max(ratios),
             "PASS" if ok else "FAIL"))
    if mod.mode == "variable":
        print("     (variable: Q_s^suave != diagonal; Q_full si coincide "
              "-- la diagonal necesita h0 con primos en [z,s))")
    return ok


# ------------------------------------------------------------------ H2
def check_H2(A, z):
    # (a) identidad exacta sin criba: Q sin restriccion == Q_completo(z)
    qu = Q_unrestricted(A, z)
    qc = Q_completo(z, A)
    e = relerr(qu, qc)
    ok_a = e < 1e-12
    print("H2(a) [A=%d z=%d] Q_sin_criba=%.12f Q_completo=%.12f "
          "relerr=%.2e %s" % (A, z, qu, qc, e,
                              "PASS" if ok_a else "FAIL"))
    return ok_a


def check_H2b(mod, tag):
    # (b) con la restriccion V^+: cociente Q_s / Q_completo por s.
    qc = Q_completo(mod.z, mod.A, mod.excl)
    ratios = []
    for s in mod.sprimes:
        ratios.append(mod.Q[s] / qc if qc else float("inf"))
    print("H2(b) [%s] Q_completo=%.6f | Q_s/Q_completo en [%.3f, %.3f] "
          "(media %.3f) -- orden 1, V^+ conserva la cota como medicion"
          % (tag, qc, min(ratios), max(ratios),
             sum(ratios) / len(ratios)))
    return True


def check_H2c(A):
    # (c) Q_completo(z) ~ cte / ln z (orden de Mertens).
    print("H2(c) [A=%d] Q_completo(z)*ln z:" % A)
    for z in (11, 17, 29, 100, 1000, 100000):
        qc = Q_completo(z, A)
        print("     z=%6d Q_completo=%.6f Q_c*ln z=%.4f" % (z, qc,
                                                           qc * math.log(z)))
    # limite teorico: e^-gamma * prod_{p|AN}(1-1/p)^-1 * C2, con
    # C2 = prod_{p no |AN} (1-1/(p-1)+1/(p-1)^2)/(1-1/p) (convergente)
    import math as _m
    euler_gamma = 0.5772156649015329
    c2 = 1.0
    for p in primes_upto(100000):
        if p == 2 or A % p == 0:
            continue
        c2 *= ((1.0 - 1.0 / (p - 1) + 1.0 / (p - 1) ** 2)
               / (1.0 - 1.0 / p))
    loc = 1.0
    for p in primes_upto(100):
        if p == 2 or A % p == 0:
            loc *= 1.0 / (1.0 - 1.0 / p)
    print("     limite teorico e^-gamma * prod_{p|AN}(1-1/p)^-1 * C2 "
          "= %.4f (C2=%.4f, factor local=%.4f; convergencia lenta, "
          "error O(1/ln z))" % (_m.exp(-euler_gamma) * loc * c2, c2, loc))
    return True


# ------------------------------------------------------------------ H3
def check_H3(mod, tag, lnN_list):
    # Perfil kappa (94e): fraccion de G_A con h <= R0, R0=(ln N)^kappa.
    # Lectura primaria: corte en h (parte nueva). Sensibilidad: corte en
    # A*h (lectura "ah<=R0" con a=A, el mayor a posible).
    gtot = mod.G_A()
    print("H3 [%s] perfil kappa, G_A=%.6f (fraccion de masa con h<=R0):"
          % (tag, gtot))
    for lnN, lab in lnN_list:
        row = []
        rowA = []
        for kap in (1.0, 2.0, 2.5, 2.7, 2.8, 3.0):
            R0 = lnN ** kap
            f = mod.G_A(Hset=[h for h in mod.H if h <= R0]) / gtot
            fA = mod.G_A(Hset=[h for h in mod.H
                               if h * mod.A <= R0]) / gtot
            row.append("%.3f" % f)
            rowA.append("%.3f" % fA)
        print("     %s (lnN=%.2f): kappa 1.0/2.0/2.5/2.7/2.8/3.0 -> %s"
              % (lab, lnN, row))
        print("       sensibilidad corte A*h<=R0:                 -> %s"
              % rowA)
    print("     Lumina 94e (N=1e12, A=3): 0.492/0.733/0.873/0.929/0.953/"
          "0.991")
    return True


# ------------------------------------------------------------------ H4
def check_H4(mod, tag, lnN):
    # Mejora del primer factor de Cauchy al restringir a r>R0.
    # f = fraccion de energia en el bloque grande (h>R0, kappa=2.7).
    gtot = mod.G_A()
    R0 = lnN ** 2.7
    f_small = mod.G_A(Hset=[h for h in mod.H if h <= R0]) / gtot
    f = 1.0 - f_small
    A, z = mod.A, mod.z
    N = A * (z ** 10)
    dA = sum(1 for a in range(1, A + 1) if A % a == 0)
    loss = dA * (A ** 2) * (math.log(N) ** 5) / math.log(z)
    if f <= 0:
        print("H4 [%s] kappa=2.7 R0=%.0f: bloque grande VACIO en este "
              "modelo (R0 > max h=%.0f); mejora no definida aqui"
              % (tag, R0, max(mod.H)))
        return True
    imp_sqrt = 1.0 / math.sqrt(f)
    imp_lin = 1.0 / f
    print("H4 [%s] kappa=2.7 R0=%.0f: masa bloque grande f=%.4f | mejora "
          "primer factor Cauchy: 1/sqrt(f)=%.2fx, 1/f=%.1fx | perdida "
          "ingenua doc95 %.2e -> tras mejora %.2e (sigue enorme)"
          % (tag, R0, f, imp_sqrt, imp_lin, loss, loss / imp_sqrt))
    return True


# ------------------------------------------------------------------ H5
def check_H5(A, z, y, D):
    # Precaucion de N: variante E0 (N=2^e: excluye solo 2 y primos|A) vs
    # E1 (N con muchos primos pequenos: excluye ademas 5,7,11,13).
    # Prediccion: Q_s(E1)/Q_s(E0) ~= prod_{p en {5,7,11,13}, p<z, p no |A}
    #   (1 - 1/(p-1) + 1/(p-1)^2)^{-1}.
    # Muestreo de 10 primos s log-uniformes en [z,y).
    allpr = primes_upto(y)
    spr = [p for p in allpr if p >= z and math.gcd(p, A) == 1]
    ns = len(spr)
    idx = sorted({round(i * (ns - 1) / 9.0) for i in range(10)})
    s_list = [spr[i] for i in idx]
    m0 = Model97(A, z, y, D, "fixed", excl=(), full=False, s_list=s_list)
    m1 = Model97(A, z, y, D, "fixed", excl=(5, 7, 11, 13), full=False,
                 s_list=s_list)
    pred = 1.0
    for p in (5, 7, 11, 13):
        if p < z and A % p != 0:
            pred *= 1.0 / (1.0 - 1.0 / (p - 1) + 1.0 / (p - 1) ** 2)
    # factor local tipo serie singular: prod_{p|N} p/(p-1)
    sing = 1.0
    for p in (5, 7, 11, 13):
        if p < z and A % p != 0:
            sing *= p / (p - 1.0)
    print("H5 [A=%d z=%d y=%d D=%d] N-factor: E0 (N=2^e) vs E1 (N con "
          "primos 5,7,11,13)" % (A, z, y, D))
    print("     prediccion producto Q(E1)/Q(E0) = %.4f | factor local "
          "tipo serie singular = %.4f" % (pred, sing))
    ratios = []
    for s in s_list:
        q0 = sum(phi_of(h0, m0.small_primes) * m0.T_of(s, h0) ** 2
                 for h0 in m0.R_of(s))
        q1 = sum(phi_of(h0, m1.small_primes) * m1.T_of(s, h0) ** 2
                 for h0 in m1.R_of(s))
        ratios.append(q1 / q0 if q0 else float("inf"))
        print("     s=%4d Q_s(E0)=%.6f Q_s(E1)=%.6f cociente=%.4f "
              "Q_s*ln z: E0=%.3f E1=%.3f"
              % (s, q0, q1, q1 / q0 if q0 else float("inf"),
                 q0 * math.log(z), q1 * math.log(z)))
    print("     cociente E1/E0 en [%.3f, %.3f], media %.3f (prediccion "
          "%.3f)" % (min(ratios), max(ratios),
                     sum(ratios) / len(ratios), pred))
    return True


# ------------------------------------------------------------------ H6
def check_H6(N, r_list):
    # Identidad de Mobius con psi reales: S_r vs sum mu(r/d)phi(d)E,
    # correccion imprimitiva Delta_r (formula derivada en el doc).
    lam = lambda_sieve(N)
    psiN = sum(lam)
    moduli = set()
    for r in r_list:
        for d in range(1, r + 1):
            if r % d == 0:
                moduli.add(d)
    bk = psi_buckets(lam, N, sorted(moduli))
    print("H6 [N=%d primo, psi(N)=%.4f] identidad de Mobius S_r vs "
          "sum mu(r/d)phi(d)E(N;d,N):" % (N, psiN))
    num = 0.0
    den = 0.0
    for r in r_list:
        cg = CharGroup(r)
        prims = cg.primitive_ts()
        # S_r = sum*_{chi* mod r} chi_bar*(N) psi(N,chi*)
        Sr = 0j
        for t in prims:
            psi_chi = 0j
            for a in cg.units:
                psi_chi += cg.chi(t, a) * bk[r][a]
            Sr += (cg.chi(t, N).conjugate()) * psi_chi
        # RHS = sum_{d|r} mu(r/d) phi(d) E(N;d,N)
        mu, _, _ = mu_spf(r)
        rhs = 0.0
        for d in range(1, r + 1):
            if r % d != 0:
                continue
            E = bk[d][N % d] - psiN / phi_of(d, primes_upto(100))
            rhs += mu[r // d] * phi_of(d, primes_upto(100)) * E
        # Delta_r predicho = sum_{p|r, p no | r/p} sum*_{chi* mod r/p}
        #   chi_bar*(N) * sum_{p^j<=N} chi*(p^j) log p
        delta = 0j
        for p in primes_upto(r):
            if r % p != 0:
                continue
            rp = r // p
            if rp % p == 0:
                continue  # p | r/p: el termino se anula
            if rp == 1:
                # Caracter primitivo mod 1 = caracter principal mod 1:
                # chi*(p^j)=1, chi_bar*(N)=1. Contribucion exacta
                # sum_{p^j<=N} log p = (log p)*floor(log_p N).
                j = 1
                pw = p
                while pw <= N:
                    delta += math.log(p)
                    j += 1
                    if pw > N // p:
                        break
                    pw *= p
                continue
            cg2 = CharGroup(rp)
            for t in cg2.primitive_ts():
                chip = cg2.chi(t, p)
                # suma directa sobre j: sum_{p^j<=N} chi*(p)^j log p
                inner = 0j
                j = 1
                pw = p
                while pw <= N:
                    inner += (chip ** j) * math.log(p)
                    j += 1
                    if pw > N // p:
                        break
                    pw *= p
                delta += cg2.chi(t, N).conjugate() * inner
        err = abs(Sr - (rhs + delta))
        scale = max(1.0, abs(Sr), abs(rhs))
        print("     r=%2d #chi*=%d S_r=%12.4f RHS=%12.4f Delta=%9.4f "
              "|S_r-(RHS+Delta)|=%.2e (rel %.1e) |Delta|/|S_r|=%s"
              % (r, len(prims), Sr.real, rhs, delta.real, err,
                 err / scale,
                 ("%.3f" % (abs(delta) / abs(Sr)) if abs(Sr) > 1e-9
                  else "S_r=0")))
        if abs(Sr) > 1e-9:
            num += (abs(delta) ** 2) / phi_of(r, primes_upto(100))
            den += (abs(Sr) ** 2) / phi_of(r, primes_upto(100))
    # cancelacion exacta del termino principal: sum_{d|r} mu(r/d) = 0
    for r in (6, 15):
        mu, _, _ = mu_spf(r)
        s = sum(mu[r // d] for d in range(1, r + 1) if r % d == 0)
        print("     cancelacion termino principal r=%d: sum mu(r/d)=%d "
              "(exacto, [T])" % (r, s))
    print("     peso relativo de la correccion en el segundo factor: "
          "sum (1/phi)|Delta|^2 / sum (1/phi)|S_r|^2 = %.3e"
          % (num / den if den else float("inf")))
    return True


# ------------------------------------------------------------------ H7
def check_H7(mod, tag):
    # Perfil de masa en r de la representacion diagonal: fraccion de Q_s
    # con r > t. (La notacion L_s del 94d queda [P-rev]: vease doc 97.)
    print("H7 [%s] masa de Q_s (forma diagonal) con r > t:" % tag)
    for s in mod.sprimes[:4]:
        R = mod.R_of(s)
        tot = mod.Q_diag(s)
        row = []
        for t in (mod.z, mod.z ** 2, mod.D // (s * mod.z), mod.D // s):
            m = sum(phi_of(r, mod.small_primes) * mod.C_of(s, r) ** 2
                    for r in R if r > t)
            row.append("r>%d: %.3f" % (t, m / tot if tot else 0.0))
        print("     s=%d Q_s=%.6f | %s" % (s, tot, " | ".join(row)))
    return True


# ------------------------------------------------------------------ main
def main():
    allok = True
    print("== DOC 97: verificacion doc 96 (diagonal / Mobius) ==")
    print("Solo stdlib. Modelos sinteticos con pesos lambda^+ de Nathanson "
          "y psi reales hasta N~1e6.")
    print()

    # S1 fixed: pequeno, todo por fuerza bruta
    m1 = Model97(3, 11, 60, 20000, "fixed")
    tag1 = "S1-fixed A=3 z=11 y=60 D=20000"
    allok &= check_H0(m1, tag1)
    allok &= check_H1(m1, tag1)
    allok &= check_H2b(m1, tag1)
    check_H7(m1, tag1)
    print()

    # S2 fixed: mediano
    m2 = Model97(3, 17, 120, 200000, "fixed")
    tag2 = "S2-fixed A=3 z=17 y=120 D=200000"
    allok &= check_H0(m2, tag2)
    allok &= check_H1(m2, tag2)
    allok &= check_H2b(m2, tag2)
    check_H7(m2, tag2)
    print()

    # S1 variable: la diagonal necesita el conjunto COMPLETO de h0
    m1v = Model97(3, 11, 60, 20000, "variable")
    tag1v = "S1-variable A=3 z=11 y=60 D=20000"
    allok &= check_H1(m1v, tag1v)
    print()

    # S4 fixed: parametros aproximados a los de Lumina (doc 94/95)
    m4 = Model97(3, 11, 3000, 200000, "fixed")
    tag4 = "S4-fixed A=3 z=11 y=3000 D=200000"
    allok &= check_H0(m4, tag4)
    allok &= check_H1(m4, tag4)
    allok &= check_H2b(m4, tag4)
    print()

    # H2(a): identidad exacta del producto completo (sin criba)
    for (A, z) in ((3, 11), (3, 17), (15, 19)):
        allok &= check_H2(A, z)
    # H2(c): orden de Mertens del producto completo
    check_H2c(3)
    print()

    # H3: perfil kappa (94e) en S1/S2/S4, con N modelo y proxy 1e12
    lnN_1e12 = math.log(10 ** 12)
    for (mod, tag) in ((m1, tag1), (m2, tag2), (m4, tag4)):
        lnN_mod = math.log(mod.A * (mod.z ** 10))
        check_H3(mod, tag, [(lnN_mod, "N modelo A*z^10"),
                            (lnN_1e12, "N proxy 1e12")])
    print()

    # H4: mejora del primer factor de Cauchy (kappa=2.7, proxy 1e12)
    for (mod, tag) in ((m1, tag1), (m2, tag2), (m4, tag4)):
        check_H4(mod, tag, lnN_1e12)
    print()

    # H5: precaucion de N (E0 vs E1), z=29
    check_H5(3, 29, 200, 29 ** 5)
    print()

    # H6: identidad de Mobius con psi reales (N primo)
    N = 999983
    assert is_prime(N)
    check_H6(N, [3, 5, 6, 7, 9, 15, 21, 25, 27, 35])
    print()
    print("PASS" if allok else "FAIL")
    return 0 if allok else 1


if __name__ == "__main__":
    sys.exit(main())
