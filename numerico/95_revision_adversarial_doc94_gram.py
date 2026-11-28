#!/usr/bin/env python3
# 95_revision_adversarial_doc94_gram.py  (Doc 95)
# Revision adversarial INDEPENDIENTE de la derivacion del doc 94 de Lumina
# Chat (norma de Gram, UNIF-U-COND). Solo stdlib.
#
# LIMITACION (ver doc 95, sec. 0): el doc 94 completo vive en OneDrive en
# el PC y NO esta disponible (tunel caido). Esta revision trabaja SOLO con
# la derivacion que Lumina compartio en el chat "Revisa UNIF U 2M X" el
# 2026-11-27 (transcrita en el doc 95) mas los docs 85-93 locales.
# Nada de lo aqui verificado depende del texto exacto del doc 94: cada
# identidad se comprueba de forma autonoma en modelos sinteticos con los
# pesos lambda^+ de Nathanson (doc 86, Teorema 9.3, beta=2).
#
# Modelo sintetico. Parametros (A, z, y, D, modo):
#   A  = modulo exterior (impar, (A,N)=1; juega el papel de Qu).
#   z  = limite inferior de criba / extremo inferior del rango switching.
#   y  = extremo superior del rango switching.
#   D  = nivel (D_A); proxy D = X_A^{1/2} sin log^C, como en docs 85-88.
#   modo = 'fixed'    : d | P(z)  (limite FIJO z, version certificada doc 87)
#          'variable' : d | P(s)  (limite VARIABLE s, version doc 85)
# Primos de switching s en [z,y), (s,A)=1. Para cada s:
#   V_s = {d : d < D/s, mu(d)!=0, (d,A)=1, d impar, primos(d) < Plim(s),
#          d en V^+(D/s)} con Plim(s)=z (fixed) o s (variable).
#   b(m) = sum_{(s,d): sd=m} mu(d)   (peso de Nathanson; multiplicidad 1
#          en fixed, <=5 en variable).
#   K_A(h) = (1/phi(A)) sum_{m: h|m} b(m)/phi(m).
#   T_s(h0) = sum_{d en V_s, h0|d} mu(d)/phi(d),  h0 suave (primos < z).
#   Q_s = sum_{h0} phi(h0) |T_s(h0)|^2.
#   G_A = sum_{h en H_A} phi(h) |K_A(h)|^2.
#
# Chequeos:
#   H0  hipotesis estructurales y conteos.
#   H1  identidad de Gram por dos calculos independientes.
#   H2  segunda reduccion en version de limite FIJO (debe valer exacta).
#   H3  segunda reduccion en version de limite VARIABLE (contraejemplo).
#   H4  Q_s: forma de Gram, diagonal vs total (cancelacion), Q_s*log z.
#   H5  objetivo sup Q_s y perdida de Cauchy en G_smooth.
#   H6  inclusion de h=1 en H_A.
#   H7  obstruccion dual: la cota ingenua de G_A no implica R_A (analitico).
#   H8  las tres correcciones de Lumina (d1 conteo O(D_A); d2 ah>R0;
#       d3 implicacion GRH en un solo sentido).
#   H9  norma signada vs mayorante absoluto (numeros de Lumina) y conteos.
# Uso: python3 95_revision_adversarial_doc94_gram.py
import math
import sys
import time


# ------------------------------------------------------------------ utiles
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
    """Factores primos distintos de d, en orden DECRECIENTE (p1 > ... > pk,
    como exige el Teorema 9.3: la condicion de cadena se lee sobre los
    primos en orden decreciente)."""
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


# ------------------------------------------------------------------ modelo
class Model:
    def __init__(self, A, z, y, D, mode, full=True):
        t0 = time.time()
        self.A, self.z, self.y, self.D, self.mode = A, z, y, D, mode
        self.full = full
        self.phiA = phi_of(A, primes_upto(10000))
        allpr = primes_upto(max(y, 100))
        self.sprimes = [p for p in allpr
                        if p >= z and p < y and math.gcd(p, A) == 1]
        # primos pequenos para phi (trial division)
        self.small_primes = primes_upto(10000)
        # --- V_s ---
        maxd = D if full else D // z
        mu, spf, _ = mu_spf(maxd)
        self.mu, self.spf = mu, spf
        self.Vs = {}
        for s in self.sprimes:
            Ds = D // s
            Plim = z if mode == "fixed" else s
            Vs = []
            if mode == "fixed":
                # enumeracion por subconjuntos de primos < z, (p,A)=1, p impar
                allow = [p for p in allpr if 2 < p < z and math.gcd(p, A) == 1]
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
        # --- T_s(h0) y Q_s ---
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
        return sum(phi_of(h, self.small_primes) * self.K[h] ** 2 for h in Hset)

    def G_A_bruteforce(self):
        """(1/phi(A)^2) sum_{m,n} b(m)b(n)(m,n)/(phi(m)phi(n)), doble suma
        directa con mcd: NO usa la agrupacion por h."""
        tot = 0.0
        ph = {m: phi_of(m, self.small_primes) for m in self.M}
        for m in self.M:
            bm = self.b[m] / ph[m]
            for n in self.M:
                tot += bm * (self.b[n] / ph[n]) * math.gcd(m, n)
        return tot / (self.phiA ** 2)

    def G_A_abs(self):
        """Mayorante absoluto: b(m) reemplazado por |b(m)| por pares."""
        babs = {}
        for m, lst in self.pairs_by_m.items():
            babs[m] = sum(abs(self.mu[d]) for (_s, d) in lst)
        Kabs = {}
        for h in self.H:
            ssum = 0.0
            for m in self.M:
                if m % h == 0:
                    ssum += babs[m] / phi_of(m, self.small_primes)
            Kabs[h] = ssum / self.phiA
        return sum(phi_of(h, self.small_primes) * Kabs[h] ** 2 for h in self.H)

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


def relerr(a, b):
    return abs(a - b) / max(1.0, abs(a), abs(b))


# ------------------------------------------------------------------ H0
def check_H0(mod, tag):
    ok = True
    msgs = []
    # (A,m)=1, m impar, m libre de cuadrados
    for m in mod.M:
        if math.gcd(m, mod.A) != 1 or m % 2 == 0 or mod.mu[m] == 0:
            ok = False
            msgs.append("m=%d viola hipotesis" % m)
            break
    maxmult = max(len(v) for v in mod.pairs_by_m.values())
    npairs = sum(len(v) for v in mod.pairs_by_m.values())
    if mod.mode == "fixed" and maxmult != 1:
        ok = False
        msgs.append("multiplicidad fixed=%d, se esperaba 1" % maxmult)
    if mod.mode == "variable" and maxmult > 5:
        ok = False
        msgs.append("multiplicidad variable=%d > 5" % maxmult)
    # h: mu^2=1, (h,A)=1
    for h in mod.H:
        if mod.mu[h] == 0 or math.gcd(h, mod.A) != 1:
            ok = False
            msgs.append("h=%d viola hipotesis" % h)
            break
    # particion en fixed
    nmulti = 0
    for h in mod.H:
        s, _h0 = mod.smooth_part(h)
        if s == "multi":
            nmulti += 1
    if mod.mode == "fixed" and nmulti != 0:
        ok = False
        msgs.append("fixed: %d h con >=2 primos switching" % nmulti)
    print("H0 [%s] |M|_pares=%d #m_distintos=%d max_mult=%d |H|=%d "
          "h_multi=%d build=%.1fs %s" % (
              tag, npairs, len(mod.M), maxmult, len(mod.H), nmulti,
              mod.build_t, "PASS" if ok else "FAIL " + ";".join(msgs)))
    return ok


# ------------------------------------------------------------------ H1
def check_H1(mod, tag):
    g1 = mod.G_A()
    g2 = mod.G_A_bruteforce()
    re = relerr(g1, g2)
    # lema sum_{h|d} phi(h) = d en muestra de pares
    lem_ok = True
    for i in range(0, len(mod.M), max(1, len(mod.M) // 200)):
        for j in range(0, len(mod.M), max(1, len(mod.M) // 200)):
            d = math.gcd(mod.M[i], mod.M[j])
            s = 0
            for h in mod.H:
                if d % h == 0:
                    s += phi_of(h, mod.small_primes)
            if s != d:
                lem_ok = False
    ok = re < 1e-9 and lem_ok
    print("H1 [%s] G_A(K)=%.12f G_A(bruto)=%.12f relerr=%.2e lema_div=%s %s"
          % (tag, g1, g2, re, lem_ok, "PASS" if ok else "FAIL"))
    return ok, g1, g2


# ------------------------------------------------------------------ H2
def check_H2(mod, tag):
    # K_A(h0) = (1/phiA) sum_s T_s(h0)/phi(s) para h0 suave
    maxe1 = 0.0
    worst1 = None
    for h0 in mod.Hsmooth:
        lhs = mod.K[h0]
        rhs = sum(mod.T[(s, h0)] / phi_of(s, mod.small_primes)
                  for s in mod.sprimes) / mod.phiA
        e = abs(lhs - rhs)
        if e > maxe1:
            maxe1, worst1 = e, h0
    # K_A(s*h0) = T_s(h0)/(phiA*phi(s))
    maxe2 = 0.0
    worst2 = None
    checked2 = 0
    for h0 in mod.Hsmooth:
        for s in mod.sprimes:
            h = s * h0
            if h not in mod.K:
                continue
            checked2 += 1
            lhs = mod.K[h]
            rhs = mod.T[(s, h0)] / (mod.phiA * phi_of(s, mod.small_primes))
            e = abs(lhs - rhs)
            if e > maxe2:
                maxe2, worst2 = e, (s, h0)
    # particion G = G_smooth + G_switch
    gs = sum(phi_of(h, mod.small_primes) * mod.K[h] ** 2
             for h in mod.H if mod.smooth_part(h)[0] is None)
    gw = sum(phi_of(h, mod.small_primes) * mod.K[h] ** 2
             for h in mod.H if mod.smooth_part(h)[0] not in (None, "multi"))
    gtot = mod.G_A()
    epart = abs(gs + gw - gtot)
    ok = maxe1 < 1e-9 and maxe2 < 1e-9 and epart < 1e-9
    print("H2 [%s] K(h0) maxerr=%.2e (h0=%s) | K(s*h0) maxerr=%.2e (%s, "
          "n=%d) | G_smooth=%.10f G_switch=%.10f G_tot=%.10f epart=%.2e %s"
          % (tag, maxe1, worst1, maxe2, worst2, checked2, gs, gw, gtot,
             epart, "PASS" if ok else "FAIL"))
    return ok


# ------------------------------------------------------------------ H3
def check_H3(mod, tag):
    # En modo variable: K_A(s*h0) debe DIFERIR de T_s(h0)/(phiA*phi(s)).
    # Descomposicion exacta con los terminos s'>s:
    #   K_A(s*h0)*phiA = T_s(h0)/phi(s) + sum_{s'>s} U_{s'}(s*h0)/phi(s')
    # con U_{s'}(t) = sum_{d en V_{s'}, t|d} mu(d)/phi(d).
    def U(sp, t):
        tot = 0.0
        for d in mod.Vs[sp]:
            if d % t == 0:
                tot += mod.mu[d] / phi_of(d, mod.small_primes)
        return tot
    maxdisc = 0.0
    worst = None
    maxdecomp = 0.0
    nmulti = 0
    gmulti = 0.0
    for h in mod.H:
        s, h0 = mod.smooth_part(h)
        if s == "multi":
            nmulti += 1
            gmulti += phi_of(h, mod.small_primes) * mod.K[h] ** 2
    for h0 in mod.Hsmooth:
        for s in mod.sprimes:
            h = s * h0
            if h not in mod.K:
                continue
            claimed = mod.T[(s, h0)] / (mod.phiA * phi_of(s, mod.small_primes))
            true = mod.K[h]
            disc = abs(true - claimed)
            if disc > maxdisc:
                maxdisc, worst = disc, (s, h0, true, claimed)
            # descomposicion con terminos s'>s
            miss = 0.0
            for sp in mod.sprimes:
                if sp > s:
                    miss += U(sp, h) / phi_of(sp, mod.small_primes)
            lhs = true * mod.phiA
            rhs = mod.T[(s, h0)] / phi_of(s, mod.small_primes) + miss
            maxdecomp = max(maxdecomp, abs(lhs - rhs))
    gtot = mod.G_A()
    frac = gmulti / gtot if gtot else 0.0
    ok = maxdisc > 1e-6 and maxdecomp < 1e-9 and nmulti > 0
    print("H3 [%s] CONTRAEJEMPLO variable: max|K(s*h0)-T_s/(phiA*phi(s))|"
          "=%.6f en %s | descomposicion con s'>s maxerr=%.2e | "
          "h_multi=%d G_multi/G_A=%.4f %s"
          % (tag, maxdisc, worst, maxdecomp, nmulti, frac,
             "PASS(contraejemplo confirmado)" if ok else "FAIL"))
    if worst:
        s, h0, true, claimed = worst
        print("     detalle: s=%d h0=%d h=%d K_A(h)=%.10f formula=%.10f "
              "terminos s'>s presentes" % (s, h0, s * h0, true, claimed))
    return ok


# ------------------------------------------------------------------ H4
def check_H4(mod, tag):
    # Q_s directo vs forma de Gram sum_{d,e} mu(d)mu(e)(d,e)/(phi(d)phi(e)).
    # Dos versiones de Q_s: 'smooth' (h0 con primos < z, la de la
    # derivacion transcrita) y 'full' (h0 sobre TODOS los divisores de
    # algun d en V_s). En modo fixed coinciden; en modo variable difieren:
    # la expansion de Gram necesita todos los divisores de (d,e).
    maxe_s = 0.0
    maxe_f = 0.0
    worst_s = None
    worst_f = None
    diag_ratio = []
    qlogz = []
    for s in mod.sprimes:
        Vs = mod.Vs[s]
        tot = 0.0
        for d in Vs:
            md = mod.mu[d] / phi_of(d, mod.small_primes)
            for e in Vs:
                tot += md * (mod.mu[e] / phi_of(e, mod.small_primes)) * math.gcd(d, e)
        # Q_s full: h0 sobre todos los divisores de d en Vs
        Hfull = {1}
        for d in Vs:
            for h in divisors_from_pf(pf_desc(d, mod.spf)):
                Hfull.add(h)
        qfull = 0.0
        for h0 in sorted(Hfull):
            t = 0.0
            for d in Vs:
                if d % h0 == 0:
                    t += mod.mu[d] / phi_of(d, mod.small_primes)
            qfull += phi_of(h0, mod.small_primes) * t * t
        es = abs(tot - mod.Q[s])
        ef = abs(tot - qfull)
        if es > maxe_s:
            maxe_s, worst_s = es, s
        if ef > maxe_f:
            maxe_f, worst_f = ef, s
        diag = sum(d / phi_of(d, mod.small_primes) ** 2 for d in Vs)
        diag_ratio.append((s, diag, mod.Q[s], diag / mod.Q[s] if mod.Q[s] else float("inf")))
        qlogz.append((s, mod.Q[s] * math.log(mod.z)))
    qmax = max(qlogz, key=lambda t: t[1])
    qmin = min(qlogz, key=lambda t: t[1])
    ok = maxe_f < 1e-9 and (mod.mode == "fixed" or True)
    ok_smooth = maxe_s < 1e-9
    print("H4 [%s] Q_s vs Gram: smooth maxerr=%.2e (s=%s) full maxerr=%.2e "
          "(s=%s) | max_s Q_s*ln z=%.4f (s=%d) min=%.4f (s=%d) %s"
          % (tag, maxe_s, worst_s, maxe_f, worst_f, qmax[1], qmax[0],
             qmin[1], qmin[0],
             "PASS" if ok else "FAIL"))
    if mod.mode == "variable" and not ok_smooth:
        print("     (variable: Q_s^smooth != forma de Gram -- la expansion "
              "necesita h0 con primos en [z,s); en fixed coinciden)")
    for (s, dg, q, r) in diag_ratio[:6]:
        print("     s=%d diagonal=%.4f Q_s=%.6f diag/Q_s=%.1f" % (s, dg, q, r))
    return ok, qlogz


# ------------------------------------------------------------------ H5
def check_H5(mod, tag):
    S = sum(1.0 / phi_of(s, mod.small_primes) for s in mod.sprimes)
    gswitch = sum(mod.Q[s] / phi_of(s, mod.small_primes)
                  for s in mod.sprimes) / (mod.phiA ** 2)
    gsmooth = sum(phi_of(h, mod.small_primes) * mod.K[h] ** 2
                  for h in mod.H if mod.smooth_part(h)[0] is None)
    cauchy = S * gswitch
    loss = cauchy / gsmooth if gsmooth > 0 else float("inf")
    cmax = max(mod.Q[s] * math.log(mod.z) for s in mod.sprimes)
    bound = (cmax / (mod.phiA ** 2 * math.log(mod.z))) * S * (1.0 + S)
    gtot = mod.G_A()
    # La cota implicada supone la particion completa (fixed). En variable
    # hay ademas G_multi, asi que ahi solo se comprueba la desigualdad de
    # Cauchy de la parte suave.
    if mod.mode == "fixed":
        ok = gtot <= bound * (1 + 1e-9) and gsmooth <= cauchy * (1 + 1e-9)
    else:
        gmulti = sum(phi_of(h, mod.small_primes) * mod.K[h] ** 2
                     for h in mod.H if mod.smooth_part(h)[0] == "multi")
        ok = gsmooth <= cauchy * (1 + 1e-9)
        print("     (variable: G_multi=%.6f explica que G_A=%.6f supere "
              "la cota %.6f)" % (gmulti, gtot, bound))
    print("H5 [%s] S=sum 1/phi(s)=%.4f G_switch=%.6f G_smooth=%.6f "
          "Cauchy(S*G_switch)=%.6f perdida=%.3f | C=max Q_s*ln z=%.4f "
          "cota_implicada=%.6f G_A=%.6f ratio=%.3f %s"
          % (tag, S, gswitch, gsmooth, cauchy, loss, cmax, bound, gtot,
             bound / gtot if gtot else float("inf"),
             "PASS" if ok else "FAIL"))
    return ok


# ------------------------------------------------------------------ H6
def check_H6(mod, tag):
    g_with = mod.G_A()
    g_without = mod.G_A(Hset=[h for h in mod.H if h != 1])
    k1sq = mod.K[1] ** 2
    e1 = abs((g_with - g_without) - k1sq)
    g2 = mod.G_A_bruteforce()
    e2 = relerr(g_with, g2)
    e3 = relerr(g_without, g2)
    ok = e1 < 1e-12 and e2 < 1e-9 and e3 > 1e-6
    print("H6 [%s] |K_A(1)|^2=%.10f (%.2f%% de G_A) G_con1-G_sin1 "
          "dif=%.2e | bruto vs con1 relerr=%.2e vs sin1 relerr=%.2e %s"
          % (tag, k1sq, 100 * k1sq / g_with if g_with else 0, e1, e2, e3,
             "PASS" if ok else "FAIL"))
    return ok


# ------------------------------------------------------------------ H7
def check_H7(mod, tag):
    # Obstruccion dual (analitico, no medicion finita). Con G_A << 1/(phi(A)^2
    # log z) y la cota ingenua del dual D <= d(A)*A*(N+(A*D_A)^2)*N*log N
    # (criba grande), el cociente (cota de |R_A|)^2/(objetivo)^2 es:
    #   d(A) * A^2 * (log N)^5 / log z   (>> 1).
    # Aqui N = A * z^10 (escala del modelo, X_A = z^10).
    A = mod.A
    z = mod.z
    N = A * (z ** 10)
    dA = sum(1 for a in range(1, A + 1) if A % a == 0)
    ratio = dA * (A ** 2) * (math.log(N) ** 5) / math.log(z)
    print("H7 [%s] dual ingenuo: (cota|R_A|)^2/(objetivo)^2 ~= "
          "d(A)*A^2*(ln N)^5/ln z = %.3e (N=A*z^10=%.3e) -- la cota "
          "ingenua pierde por ese factor; el paso G_A -> R_A NO esta en "
          "la derivacion transcrita" % (tag, ratio, float(N)))
    return True


# ------------------------------------------------------------------ H8
def check_H8(mod, tag):
    # (d1) conteo de modulos: pares vs distintos
    npairs = sum(len(v) for v in mod.pairs_by_m.values())
    ndist = len(mod.M)
    ok1 = npairs <= 5 * mod.D and ndist <= mod.D
    print("H8(d1) [%s] pares(s,d)=%d distintos=%d D_A=%d: pares<=5*D_A %s, "
          "distintos<=D_A %s -> el conteo es O(D_A), no necesariamente "
          "<= D_A" % (tag, npairs, ndist, mod.D, ok1, ok1))
    # (d2) fraccion de Q_s proveniente de h0 <= R0
    for R0 in (10, 30, 100):
        fr = []
        for s in mod.sprimes[:6]:
            tot = mod.Q[s]
            small = sum(phi_of(h0, mod.small_primes) * mod.T[(s, h0)] ** 2
                        for h0 in mod.Hsmooth if h0 <= R0)
            fr.append(small / tot if tot else 0.0)
        print("H8(d2) [%s] frac Q_s de h0<=%d (primeros s): %s -- la "
              "estimacion de bloque debe restringirse a ah>R0; el rango "
              "ah<=R0 necesita tratamiento aparte [P-rev en doc 94]"
              % (tag, R0, ["%.3f" % f for f in fr]))
    # (d3) implicacion GRH en un solo sentido: se verifica el mecanismo
    # (ortogonalidad) como en doc 93 H7, r=3, N=1e5.
    N = 100000
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
            p = i
            pw = i * i
            while pw <= N:
                lam[pw] = math.log(i)
                pw *= i
    psiN = sum(lam)
    def psi_mod(m, a):
        s = 0.0
        n = a % m
        if n == 0:
            n = m
        while n <= N:
            s += lam[n]
            n += m
        return s
    # chi mod 3 no principal: chi(1)=1, chi(2)=-1
    chi = {1: 1.0, 2: -1.0}
    direct = sum(lam[n] * (chi[n % 3] if n % 3 in chi else 0.0)
                 for n in range(2, N + 1))
    viaE = 0.0
    for a in (1, 2):
        E = psi_mod(3, a) - psiN / 2.0
        viaE += chi[a] * E
    ok3 = abs(direct - viaE) < 1e-6 * max(1.0, abs(direct))
    print("H8(d3) [%s] ortogonalidad r=3: psi(N,chi) directo=%.4f "
          "via E=%.4f dif=%.2e %s -- UNIF-U-2M-R => GRH [T, doc 93 Prop. "
          "93.6]; GRH =/=> UNIF-U-2M-R NO probado [N]; 'estrictamente mas "
          "fuerte que GRH' queda DEGRADADO a implicacion en un solo sentido"
          % (tag, direct, viaE, abs(direct - viaE),
             "PASS" if ok3 else "FAIL"))
    return ok1 and ok3


# ------------------------------------------------------------------ H9
def check_H9(mod, tag, lumina_counts=None, lumina_quot=None):
    g = mod.G_A()
    gabs = mod.G_A_abs()
    q = g / gabs if gabs else float("inf")
    qs = math.sqrt(q)
    print("H9 [%s] |M|=%d |H|=%d G_A=%.10f G_A^abs=%.10f cociente "
          "G/G^abs=%.6f sqrt=%.6f" % (tag, len(mod.M), len(mod.H), g, gabs,
                                      q, qs))
    if lumina_counts:
        print("     Lumina (doc 94): 1444 pesos, 1451 partes nuevas, "
              "norma signada 0.4505669522, mayorante 1.1723577317, "
              "cociente 0.384325; aqui |M|=%d |H|=%d"
              % (len(mod.M), len(mod.H)))
    return True


# ------------------------------------------------------------------ tendencia Q_s
def trend_Qs(A, zs):
    print("TENDENCIA Q_s (modo fixed, D=z^5, y=z^(10/3), A=%d):" % A)
    for z in zs:
        D = z ** 5
        y = int(z ** (10.0 / 3.0)) + 1
        mod = Model(A, z, y, D, "fixed", full=False)
        # muestreo de 10 primos s log-uniformes en [z,y)
        ns = len(mod.sprimes)
        idx = sorted({round(i * (ns - 1) / 9.0) for i in range(10)})
        vals = []
        for i in idx:
            s = mod.sprimes[i]
            # T_s y Q_s para este s
            Vs = mod.Vs[s]
            # h0: divisores suaves de d en Vs
            H0 = {1}
            for d in Vs:
                for h in divisors_from_pf(pf_desc(d, mod.spf)):
                    H0.add(h)
            q = 0.0
            diag = 0.0
            for h0 in sorted(H0):
                t = 0.0
                for d in Vs:
                    if d % h0 == 0:
                        t += mod.mu[d] / phi_of(d, mod.small_primes)
                q += phi_of(h0, mod.small_primes) * t * t
            for d in Vs:
                diag += d / phi_of(d, mod.small_primes) ** 2
            vals.append((s, q, q * math.log(z), diag))
        qmax = max(v[2] for v in vals)
        qmin = min(v[2] for v in vals)
        dg = max(v[3] for v in vals)
        print("  z=%d D=%d y=%d #s=%d | Q_s*ln z en [%0.4f, %0.4f] "
              "max=%.4f diag_max=%.2f (diag/Q_s~%.0f)"
              % (z, D, y, ns, qmin, qmax, qmax, dg,
                 dg / max(v[1] for v in vals)))
    return True


# ------------------------------------------------------------------ main
def main():
    allok = True
    print("== DOC 95: revision adversarial doc 94 (Gram / UNIF-U-COND) ==")
    print("Solo stdlib. Modelos sinteticos con pesos lambda^+ de Nathanson.")
    print()

    # S1: pequeno, ambos modos, todo por fuerza bruta
    for mode in ("fixed", "variable"):
        m1 = Model(3, 11, 60, 20000, mode)
        tag = "S1-%s A=3 z=11 y=60 D=20000" % mode
        allok &= check_H0(m1, tag)
        ok, _, _ = check_H1(m1, tag)
        allok &= ok
        if mode == "fixed":
            allok &= check_H2(m1, tag)
        else:
            allok &= check_H3(m1, tag)
        ok, _ = check_H4(m1, tag)
        allok &= ok
        allok &= check_H5(m1, tag)
        allok &= check_H6(m1, tag)
        check_H7(m1, tag)
        allok &= check_H8(m1, tag)
        check_H9(m1, tag)
        print()

    # S2: mediano, ambos modos
    for mode in ("fixed", "variable"):
        m2 = Model(3, 17, 120, 200000, mode)
        tag = "S2-%s A=3 z=17 y=120 D=200000" % mode
        allok &= check_H0(m2, tag)
        ok, _, _ = check_H1(m2, tag)
        allok &= ok
        if mode == "fixed":
            allok &= check_H2(m2, tag)
        else:
            allok &= check_H3(m2, tag)
        ok, _ = check_H4(m2, tag)
        allok &= ok
        allok &= check_H5(m2, tag)
        allok &= check_H6(m2, tag)
        check_H7(m2, tag)
        check_H9(m2, tag)
        print()

    # S3: A=15, fixed
    m3 = Model(15, 19, 150, 500000, "fixed")
    tag = "S3-fixed A=15 z=19 y=150 D=500000"
    allok &= check_H0(m3, tag)
    ok, _, _ = check_H1(m3, tag)
    allok &= ok
    allok &= check_H2(m3, tag)
    ok, _ = check_H4(m3, tag)
    allok &= ok
    allok &= check_H5(m3, tag)
    allok &= check_H6(m3, tag)
    check_H7(m3, tag)
    allok &= check_H8(m3, tag)
    check_H9(m3, tag)
    print()

    # S4: hipotesis de parametros de Lumina (1444 pesos / 1451 partes)
    m4 = Model(3, 11, 3000, 200000, "fixed")
    tag = "S4-fixed A=3 z=11 y=3000 D=200000"
    allok &= check_H0(m4, tag)
    ok, _, _ = check_H1(m4, tag)
    allok &= ok
    allok &= check_H2(m4, tag)
    ok, _ = check_H4(m4, tag)
    allok &= ok
    allok &= check_H5(m4, tag)
    allok &= check_H6(m4, tag)
    check_H7(m4, tag)
    check_H9(m4, tag, lumina_counts=True)
    print()

    # Tendencia de Q_s con la escala (D=z^5, y=z^(10/3))
    trend_Qs(3, [11, 13, 17, 19, 23, 29])
    print()
    print("PASS" if allok else "FAIL")
    return 0 if allok else 1


if __name__ == "__main__":
    sys.exit(main())
