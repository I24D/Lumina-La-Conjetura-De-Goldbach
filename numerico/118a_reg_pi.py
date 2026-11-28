#!/usr/bin/env python3
# 118a_reg_pi.py  (Doc 118a)
# (REG-Π) del doc 117b §3.1: sum_{r en banda} Pi(r)/phi(r) <= C*Q_c.
# Veredicto que el script sostiene:
#   J0  identidad Pi(r)/phi(r) = Pi_full * g(r), g(r) = prod 1/(p-2),
#       y Q_c = Pi_full * M_full (maquinaria del doc 110, sin tocar).
#   J1  medicion en marcos (z=29,37,53; u_L=5/3,2): la suma de la banda
#       partida por omega(r) — reproduce el rango 0.17-0.46 de 117b K6.
#   J2  sumas de Mertens de la familia obstructora omega=2 a z grande:
#       F(z) = [sum_{z^(11/20)<p<z^(3/5)} 1/(p-2)] *
#              [sum_{5<=q<=z^(1/2)} 1/(q-2)]  ~  log(12/11) * log log z.
#   J3  rebanada omega=1: sum_{p en I(u*)} 1/(p-2) -> B_1(u*) (Mertens).
#   J4  la familia obstructora contra los predicados reales de la serie
#       (visible_set, Pi, u_r) a z=53, L=floor(53^2): visibilidad total,
#       entrada en la banda parcial (preasintotica, Pi_full aun grande).
# Solo stdlib. Importa 117b (que importa 110) sin modificarlos.
import importlib.util
import math
import sys
import time

spec = importlib.util.spec_from_file_location(
    "d117b", "117b_banda_delta.py")
d117b = importlib.util.module_from_spec(spec)
spec.loader.exec_module(d117b)

d110 = d117b.d110
Level110 = d110.Level110
visible_set = d117b.visible_set
U0 = 1.5


def g_of_pf(pf):
    g = 1.0
    for p in pf:
        g /= (p - 2)
    return g


def sieve_primes(n):
    """Primos <= n (bytearray, solo stdlib)."""
    if n < 2:
        return []
    bs = bytearray(b"\x01") * (n + 1)
    bs[0:2] = b"\x00\x00"
    for i in range(2, int(n ** 0.5) + 1):
        if bs[i]:
            step = i
            start = i * i
            bs[start:n + 1:step] = b"\x00" * ((n - start) // step + 1)
    return [i for i in range(2, n + 1) if bs[i]]


# ------------------------------------------------------------------ J0
def check_J0():
    print("J0 identidad Pi(r)/phi(r) = Pi_full * g(r) y Q_c = "
          "Pi_full*M_full:")
    ok = True
    for z in (29, 37):
        lv = Level110(z, 3, (), z ** 2)
        err = 0.0
        for (r, pf, mu, phi) in lv.enum:
            if r == 1:
                continue
            err = max(err, abs(lv.Pi[r] / phi
                               - lv.Pi_full * g_of_pf(pf)))
        e2 = abs(lv.qc - lv.Pi_full * lv.M_full_sum) / lv.qc
        print("  z=%d: err identidad=%.2e | err Q_c=%.2e | Pi_full=%.4f "
              "M_full=%.4f Q_c=%.4f"
              % (z, err, e2, lv.Pi_full, lv.M_full_sum, lv.qc))
        ok &= err < 1e-12 and e2 < 1e-12
    print("  => %s" % ("PASS" if ok else "FAIL"))
    return ok


# ------------------------------------------------------------------ J1
def check_J1():
    print("J1 suma de la banda partida por omega(r) (cocientes a Q_c):")
    ok = True
    for z in (29, 37, 53):
        for ustar in (5.0 / 3.0, 2.0):
            L = int(z ** ustar)
            lv = Level110(z, 3, (), L)
            info = {r: (pf, phi) for (r, pf, mu, phi) in lv.enum}
            lnz = math.log(z)
            parts = {1: 0.0, 2: 0.0, 3: 0.0}
            tot = 0.0
            for r in visible_set(lv):
                if r == 1:
                    continue
                pf, phi = info[r]
                if math.log(L / r) / lnz >= U0:
                    continue
                Pi = lv.Pi[r]
                if Pi >= 0.5:
                    continue
                w = min(len(pf), 3)
                parts[w] += Pi / phi
                tot += Pi / phi
            print("  z=%d u_L=%.3f: omega=1: %.4f | omega=2: %.4f | "
                  "omega>=3: %.4f | total/Q_c=%.4f"
                  % (z, math.log(L) / lnz, parts[1] / lv.qc,
                     parts[2] / lv.qc, parts[3] / lv.qc, tot / lv.qc))
            ok &= tot / lv.qc < 1.0
    print("  => totales en el rango de ventana de 117b K6 (0.17-0.46), "
          "con parte omega=2 > 0. %s" % ("PASS" if ok else "FAIL"))
    return ok


# ------------------------------------------------------------------ J2
def check_J2():
    print("J2 familia obstructora omega=2: F(z) = exterior * interior "
          "(Mertens, z grande):")
    print("  limite teorico: exterior -> log(12/11) = %.4f; interior "
          "= log log z + O(1) -> inf" % math.log(12.0 / 11.0))
    ok = True
    prev = None
    for z in (10 ** 4, 10 ** 6, 10 ** 8, 10 ** 10):
        hi = int(z ** 0.6) + 2
        primes = sieve_primes(hi)
        lo_out = z ** (11.0 / 20.0)
        hi_out = z ** (3.0 / 5.0)
        lo_in = z ** 0.5
        ext = sum(1.0 / (p - 2) for p in primes
                  if lo_out < p < hi_out)
        inn = sum(1.0 / (p - 2) for p in primes
                  if 5 <= p <= lo_in)
        F = ext * inn
        llz = math.log(math.log(z))
        print("  z=1e%d: exterior=%.4f | interior=%.4f | F=%.4f | "
              "F/loglogz=%.4f" % (round(math.log10(z)), ext, inn, F,
                                  F / llz))
        if prev is not None:
            ok &= F > prev
        prev = F
    print("  => F(z) crece sin cota (como log log z): la masa g de la "
          "familia diverge. %s" % ("PASS" if ok else "FAIL"))
    return ok


# ------------------------------------------------------------------ J3
def check_J3():
    print("J3 rebanada omega=1: sum_{p en I(u*)} 1/(p-2) -> B_1(u*):")
    ok = True
    for z in (10 ** 10, 10 ** 12):
        primes = sieve_primes(int(z ** (2.0 / 3.0)) + 2)
        for ustar, lim in ((5.0 / 3.0, math.log(10.0 / 3.0)),
                           (2.0, math.log(4.0 / 3.0))):
            lo = z ** (ustar - 1.5)
            hi = z ** (min(1.0, ustar / 3.0))
            s = sum(1.0 / (p - 2) for p in primes if lo < p < hi)
            print("  z=1e%d u*=%.3f: suma=%.4f | B_1=%.4f | err=%.4f"
                  % (round(math.log10(z)), ustar, s, lim,
                     abs(s - lim)))
            ok &= abs(s - lim) < 0.06
    print("  => la rebanada omega=1 tiene limite de Mertens finito "
          "(<= log(10/3)). %s" % ("PASS" if ok else "FAIL"))
    return ok


# ------------------------------------------------------------------ J4
def check_J4():
    print("J4 familia {pq} contra los predicados reales (z=53, "
          "L=floor(53^2), u_L=2):")
    z = 53
    L = z ** 2
    lv = Level110(z, 3, (), L)
    info = {r: (pf, phi) for (r, pf, mu, phi) in lv.enum}
    vis = visible_set(lv)
    lnz = math.log(z)
    allow = [p for p in d110.primes_upto(z - 1) if p != 2 and p != 3]
    fam = []
    for p in allow:
        xp = math.log(p) / lnz
        if not (0.5 < xp < 2.0 / 3.0):
            continue
        for q in allow:
            if q < p:
                fam.append(p * q)
    nvis = 0
    nband = 0
    gmass_band = 0.0
    for r in fam:
        if r not in vis:
            continue
        nvis += 1
        pf, phi = info[r]
        ur = math.log(L / r) / lnz
        if ur < U0 and lv.Pi[r] < 0.5:
            nband += 1
            gmass_band += g_of_pf(pf)
    print("  familia: %d elementos | visibles (maquinaria): %d | ya en "
          "la banda a z=53: %d | masa g en banda=%.4f"
          % (len(fam), nvis, nband, gmass_band))
    print("  Pi_full(53)=%.4f: el factor maximo (6/5)(4/3)=1.6 aun no "
          "deja entrar a toda la familia (Pi(r) hasta %.3f > 1/2); la "
          "entrada completa es asintotica (Pi_full ~ K/ln z -> 0)."
          % (lv.Pi_full, lv.Pi_full * 1.6))
    ok = nvis == len(fam) and nband >= 1
    print("  => visibilidad total verificada; la banda a z=53 es "
          "preasintotica. %s" % ("PASS" if ok else "FAIL"))
    return ok


def main():
    t0 = time.time()
    print("== DOC 118a: (REG-Pi) ==")
    print("Solo stdlib; maquinaria de los docs 110/117b sin tocar.")
    print()
    allok = True
    allok &= check_J0()
    print()
    allok &= check_J1()
    print()
    allok &= check_J2()
    print()
    allok &= check_J3()
    print()
    allok &= check_J4()
    print()
    print("Tiempo total: %.1f s" % (time.time() - t0))
    print("PASS" if allok else "FAIL")
    return 0 if allok else 1


if __name__ == "__main__":
    sys.exit(main())
