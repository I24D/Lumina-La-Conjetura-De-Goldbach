#!/usr/bin/env python3
# 127_reparaciones.py  (Doc 127)
# Reparaciones y verificaciones de la errata del doc 125 tras la
# auditoria adversarial del doc 126 (LUMINA Chat):
#   R1 dominio de la formula PW: (3/2, 9/4), no (3/2, 3).
#   R2 constante de Mertens M': correccion +2/3 (5/6, no 3/2).
#   R3 signo de eps_Pi: suma filtrada = completa - eps_Pi (125a).
#   R4 chequeo K3 del 125b reparado (sqfree_cs) y re-ejecucion
#      real de la enumeracion del estrato 5 (barrido + K4).
#   R5 contraejemplo exacto de la auditora (z=61, L=3926, r=11)
#      verificado con las funciones de la serie (docs 120a/125b).
# Solo stdlib.
import importlib.util
import math
import os
import time
from fractions import Fraction

BASE = os.path.dirname(os.path.abspath(__file__))
os.chdir(BASE)


def load(name, fname):
    spec = importlib.util.spec_from_file_location(name,
                                                  os.path.join(BASE, fname))
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod


def primes_upto(n):
    s = [True] * (n + 1)
    s[0] = s[1] = False
    for i in range(2, int(n ** 0.5) + 1):
        if s[i]:
            for j in range(i * i, n + 1, i):
                s[j] = False
    return [i for i, v in enumerate(s) if v]


def is_prime(n):
    if n < 2:
        return False
    i = 2
    while i * i <= n:
        if n % i == 0:
            return False
        i += 1 if i == 2 else 2
    return True


# ------------------------------------------------------------------ R1
def check_R1():
    print("R1 dominio de la formula PW (correccion 1 del doc 126):")
    ok = True

    def F(u):
        return math.log(3.0 / u) * math.log((u / 3.0) / (u - 1.5))

    f2 = F(2.0)
    print("  F(2) = log(3/2)*log(4/3) = %.10f (doc 126: 0.1166450426)"
          % f2)
    ok &= abs(f2 - 0.1166450426) < 1e-9
    lnz = math.log(10.0 ** 6)
    for u in (1.6, 2.0, 2.2, 2.25, 2.5, 2.9):
        lnL = u * lnz
        nonempty = (lnL - 1.5 * lnz) < (lnL / 3.0)
        lg = math.log((u / 3.0) / (u - 1.5))
        print("  u=%.2f: intervalo (L z^-3/2, L^1/3] no vacio: %s |"
              " log((u/3)/(u-3/2)) = %+.6f" % (u, nonempty, lg))
        if u < 2.25:
            ok &= nonempty and lg > 0
        elif u == 2.25:
            ok &= (not nonempty) and abs(lg) < 1e-12
        else:
            ok &= (not nonempty) and lg < 0
    print("  => dominio correcto (3/2, 9/4); u=2 dentro: %s"
          % ("PASS" if ok else "FAIL"))
    return ok


# ------------------------------------------------------------------ R2
def check_R2():
    print("R2 constante de Mertens M' (correccion 2 del doc 126):")
    ok = True
    sub = Fraction(1, 2) + Fraction(1, 3)
    dif = Fraction(3, 2) - sub
    print("  quitar p=2,3 de sum 1/p resta 1/2+1/3 = %s (no 3/2);"
          " 3/2 - 5/6 = %s" % (sub, dif))
    ok &= sub == Fraction(5, 6) and dif == Fraction(2, 3)
    x = 100000
    pr = primes_upto(x)
    s1 = sum(1.0 / (p - 1) for p in pr if p >= 5)
    s2 = sum(1.0 / p for p in pr if p >= 5)
    s3 = sum(1.0 / (p * (p - 1)) for p in pr if p >= 5)
    sall = sum(1.0 / p for p in pr)
    print("  identidad finita x=1e5: |S_{1/(p-1)} - (S_{1/p} +"
          " S_{1/(p(p-1))})| = %.2e | sum_{p>=5} 1/p = (sum_{p<x}"
          " 1/p) - 5/6, residuo %.2e"
          % (abs(s1 - (s2 + s3)), abs(s2 - (sall - 5.0 / 6.0))))
    ok &= abs(s1 - (s2 + s3)) < 1e-9
    ok &= abs(s2 - (sall - 5.0 / 6.0)) < 1e-12
    X = 10 ** 7
    prX = primes_upto(X)
    serie = sum(1.0 / (p * (p - 1)) for p in prX if p >= 5)
    B1 = 0.26149721284764278
    Mcorr = B1 - 5.0 / 6.0 + serie
    Mimpr = B1 - 1.5 + serie
    SX = sum(1.0 / (p - 1) for p in prX if p >= 5)
    Ecorr = SX - math.log(math.log(X)) - Mcorr
    Eimpr = SX - math.log(math.log(X)) - Mimpr
    print("  sum_{p>=5} 1/(p(p-1)) = %.8f (a 1e7; cola < 1e-7)"
          % serie)
    print("  M'_corr = %.6f | M'_impresa = %.6f | diferencia"
          " = %.6f (2/3 = %.6f)" % (Mcorr, Mimpr, Mcorr - Mimpr, 2.0 / 3.0))
    print("  E_corr(1e7) = %+.4f (debe -> 0) | E_impresa(1e7)"
          " = %+.4f (debe -> 2/3)" % (Ecorr, Eimpr))
    ok &= abs((Mcorr - Mimpr) - 2.0 / 3.0) < 1e-6
    ok &= abs(Ecorr) < 0.15 and abs(Eimpr - 2.0 / 3.0) < 0.15
    print("  => %s" % ("PASS" if ok else "FAIL"))
    return ok


# ------------------------------------------------------------------ R3
def check_R3():
    print("R3 signo de eps_Pi (correccion 3 del doc 126; esta en"
          " el 125a sec. 3, no en el 125b):")
    ok = True
    z, L = 61, 61 ** 2
    lo = L * z ** (-1.5)
    hi = L ** (1.0 / 3.0)
    pr = primes_upto(int(hi) + 1)
    interval = [p for p in pr if lo < p <= hi]
    full = sum(1.0 / (p - 1) for p in interval)
    excl = [p for p in interval if p % 3 == 2]
    eps = sum(1.0 / (p - 1) for p in excl)
    filt = full - eps
    print("  intervalo real (L z^-3/2, L^1/3] con z=61, L=61^2:"
          " primos %s | masa completa %.6f" % (interval, full))
    print("  filtro artificial (p %% 3 == 2): excluidos %s,"
          " eps = %.6f > 0" % (excl, eps))
    print("  filtrada = %.6f | completa - eps = %.6f (exacta) |"
          " completa + eps = %.6f (la impresa: falsa)"
          % (filt, full - eps, full + eps))
    ok &= eps > 0 and abs(filt - (full - eps)) < 1e-12 \
        and abs(filt - (full + eps)) > 1e-6
    print("  => %s" % ("PASS" if ok else "FAIL"))
    return ok


# ------------------------------------------------------------------ R4
def check_R4(d125b):
    print("R4 K3 del 125b reparado y enumeracion del estrato 5"
          " re-ejecutada (correccion 4 del doc 126):")
    ok = True
    pr = d125b.primes_upto(100000)
    b1, b2, b3 = 125 * 7 / 121, 125 * 11 / 169, 625 / 49

    def sqfree_orig(B):
        out = []

        def rec(start, cur):
            if cur >= 5:
                out.append(cur)
            for p in pr:
                if p < start or cur * p >= B:
                    break
                if p >= 5:
                    rec(p + 1, cur * p)
        rec(5, 1)
        return out

    def sqfree_fixed(B):
        out = []

        def rec(start, cur):
            if cur >= 5:
                out.append(cur)
            for p in pr:
                if p < start:
                    continue
                if cur * p >= B:
                    break
                rec(p + 1, cur * p)
        rec(5, 1)
        return out

    def brute(B):
        res = []
        for c in range(5, B):
            m, p, fac = c, 2, []
            while p * p <= m:
                if m % p == 0:
                    fac.append(p)
                    m //= p
                    if m % p == 0:
                        fac = None
                        break
                else:
                    p += 1
            if fac is None:
                continue
            if m > 1:
                fac.append(m)
            if fac and min(fac) >= 5:
                res.append(c)
        return res

    orig = sqfree_orig(125)
    fixed = sqfree_fixed(125)
    bf = brute(125)
    print("  sqfree_cs ORIGINAL del 125b: %s (vacia: el PASS de"
          " K3 era vacuo)" % orig)
    print("  sqfree_cs reparado: %d candidatos | fuerza bruta"
          " independiente: %d | coinciden: %s"
          % (len(fixed), len(bf), sorted(fixed) == bf))
    ok &= orig == [] and len(fixed) == 37 and sorted(fixed) == bf
    alive = []
    for c in fixed:
        if c in (1, 5, 7):
            continue
        fac = [p for p in pr if p <= c and c % p == 0]
        if 5 in fac:
            dead = c > b3
        else:
            dead = c > b1 and c > b2
        if not dead:
            alive.append(c)
    print("  filtro de ramas (POS) sobre los 37: supervivientes"
          " fuera de {1,5,7}: %s" % alive)
    ok &= not alive
    nxt = {pr[i]: pr[i + 1] for i in range(len(pr) - 1)}
    sup_bad = [p for p in pr if p >= 11 and p in nxt
               and 125.0 * p / nxt[p] ** 2 > b2 + 1e-12]
    eq = [p for p in pr if p >= 11 and p in nxt
          and abs(125.0 * p / nxt[p] ** 2 - b2) < 1e-12]
    mono_bad = [p for p in pr if 13 <= p <= 100000
                and 11 * p * p - 125 * p + 44 <= 0]
    print("  aritmetica K3 (no afectada por el bug): 875/121"
          " = %.4f < 7.24 | 1375/169 = %.4f < 8.14 | sup no"
          " excede: %s | igualdad solo p=11: %s | monotonia: %s"
          % (b1, b2, not sup_bad, eq == [11], not mono_bad))
    ok &= b1 < 7.24 and b2 < 8.14 and b3 < 12.76 \
        and not sup_bad and eq == [11] and not mono_bad
    print("  re-ejecucion del K4 original del 125b (enumeracion"
          " extendida, no usa sqfree_cs):")
    ok &= d125b.check_K4()
    print("  => %s" % ("PASS" if ok else "FAIL"))
    return ok


# ------------------------------------------------------------------ R5
def check_R5(d125b):
    print("R5 contraejemplo exacto del doc 126 sec. 2.2,"
          " verificado con las funciones de la serie:")
    ok = True
    z, L, r = 61, 3926, 11
    D, s = z ** 5, 215077
    print("  z=61: D = 61^5 = %d | s = 215077 primo: %s |"
          " floor(D/s) = %d" % (D, is_prime(s), D // s))
    ok &= D // s == L and is_prime(s)
    ur = math.log(L / r) / math.log(z)
    print("  r=11: 11^3 = %d < %d | u_r = %.5f en [1, 3/2) |"
          " L/r = %.3f <= 1331" % (r ** 3, L, ur, L / r))
    ok &= r ** 3 < L and 1.0 <= ur < 1.5 and L / r <= 1331
    lv = d125b.Level110(z, 3, (), L)
    found = None
    for tup in d125b.band_instances(lv, L, math.log(z)):
        if tup[0] == r:
            found = tup
    print("  instancia de banda (maquinaria 120a): presente %s"
          " | u_r = %.5f | Pi(11) = %.6f < 1/2 | nhat = %s > 1"
          % (found is not None, found[3], found[4], found[6]))
    ok &= found is not None and abs(found[3] - ur) < 1e-9 \
        and abs(found[4] - 151963425 / 385875968) < 1e-9 \
        and found[4] < 0.5 and found[6] > 1
    print("  testigo no degenerado: 11*7^3 = %d < %d"
          % (11 * 7 ** 3, L))
    ok &= 11 * 7 ** 3 < L
    ld = d125b.lados(lv, r)
    print("  lados(11): clase de d'=35 -> %s" % ld[35][0])
    ok &= ld[35][0] == "Rem"
    f = d125b.merged_failure([11], [5, 7], L)
    print("  fallo fusionado de 385: (m, j, p_m) = %s |"
          " 5^2*385 = %d >= %d (Nathanson falla en m=3)"
          % (f, 25 * 385, L))
    ok &= f == (3, 2, 5) and 25 * 385 >= L
    n, er, elems = d125b.rem_info(lv, r)
    print("  rem_info(11): |Rem| = %d, elementos = %s,"
          " E_Rem(11) = %.10f" % (n, [e[0] for e in elems], er))
    print("  phi(35) = 24 -> -1/phi(35) = %.10f | E_Rem(11)"
          " <= -1/24 < 0: %s" % (-1.0 / 24.0, er <= -1.0 / 24.0 + 1e-12))
    ok &= n == 1 and elems[0][0] == 35 \
        and abs(er - (-1.0 / 24.0)) < 1e-12 and er < 0
    print("  => anulacion uniforme E_Rem = 0 refutada con"
          " instancia exacta [T-instancia]: %s"
          % ("PASS" if ok else "FAIL"))
    return ok


def main():
    t0 = time.time()
    d125b = load("d125b", "125b_erem_pos.py")
    ok = True
    for chk in (check_R1, check_R2, check_R3):
        try:
            ok &= chk()
        except Exception as e:
            ok = False
            print("  EXCEPCION en %s: %r" % (chk.__name__, e))
    for chk in (check_R4, check_R5):
        try:
            ok &= chk(d125b)
        except Exception as e:
            ok = False
            print("  EXCEPCION en %s: %r" % (chk.__name__, e))
    print("Tiempo total: %.1f s" % (time.time() - t0))
    print("PASS" if ok else "FAIL")
    return 0 if ok else 1


if __name__ == "__main__":
    raise SystemExit(main())
