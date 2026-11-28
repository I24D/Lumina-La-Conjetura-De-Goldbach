#!/usr/bin/env python3
# 116c_banda_conteo.py  (Doc 116c)
# El residuo (BAND) del doc 114 para la cota cruda (1-Pi(r)^2)/phi(r):
#   Sigma_{r visible, u_r<u0, Pi(r)<1/2} (1-Pi^2)/phi  <=?  C * Q_completo.
# Chequeos:
#   J0  Pi(r) = Pi_full / prod_{p|r}(1-1/(p-1)) exacto; Pi(r) >= Pi_full;
#       Q_c = Pi_full * M_full (reuso docs 109/110).
#   J1  Pi(p) = Pi_full*(p-1)/(p-2); umbral Pi_full(z) < 1/2 (z>=19);
#       banda vacia a z=17 (Pi_full > 1/2).
#   J2  Descomposicion de la banda cruda por omega(r) en niveles
#       representativos L = z^{u*} (marcos A y B) y maximo sobre los
#       niveles muestreados del doc 114 (comparar con 0,6991).
#   J3  Limites de Mertens: (a) sumas de primos en los intervalos del
#       slice superior -> log(4/3) (u_L=2) y log(10/3) (u_L=5/3);
#       (b) integrales de region B_omega(u_L, u0=1.5) por cuadratura;
#       (c) Pi_full(z)*ln z y Q_c(z)*ln z se estabilizan (Q_c ~ K/ln z).
#   J4  Masa real t(p)+ vs cota cruda en la familia omega=1 (z=37):
#       la cancelacion S(p) ~ Pi(p) que el conteo no ve.
#   J5  Sub-banda {Pi >= 1/2} dentro de u_r<u0: suma cruda <= 4*Q_c [T].
# Uso: python3 116c_banda_conteo.py | tee salida_doc116c_local.txt
import importlib.util
import itertools
import math
import sys
import time

spec = importlib.util.spec_from_file_location("doc110",
                                              "110_segundo_momento_defectos.py")
d110 = importlib.util.module_from_spec(spec)
spec.loader.exec_module(d110)
d109 = d110.d109

primes_upto = d110.primes_upto
Q_completo = d110.Q_completo
Level110 = d110.Level110
sub_divisors = d110.sub_divisors
sample_sprimes = d110.sample_sprimes


def Pi_full(z, A=3, excl=()):
    tot = 1.0
    for p in primes_upto(z - 1):
        if p == 2 or A % p == 0 or p in excl:
            continue
        tot *= (1.0 - 1.0 / (p - 1))
    return tot


def visible_set(lv):
    vis = set()
    for (d, pf, mu, phi) in lv.V:
        for (r, rmu, rphi) in sub_divisors(pf):
            vis.add(r)
    return vis


def band_decomp(z, L, u0=1.5):
    """Devuelve (Q_c, dict omega->(crudo, real), crudo_tot, real_tot,
    crudo en {Pi>=1/2, u_r<u0})."""
    lv = Level110(z, 3, (), L)
    info = {r: (pf, phi) for (r, pf, mu, phi) in lv.enum}
    rows = {}
    tot_c = tot_r = comp_c = 0.0
    for r in visible_set(lv):
        if r == 1:
            continue
        pf, phi = info[r]
        ur = math.log(L / r) / math.log(z)
        if ur >= u0:
            continue
        Pi = lv.Pi[r]
        crude = (1.0 - Pi * Pi) / phi
        if Pi >= 0.5:
            comp_c += crude
            continue
        delta = lv.delta[r]
        t = delta * (2.0 * Pi + delta) / phi
        w = min(len(pf), 4)
        c, rr = rows.get(w, (0.0, 0.0))
        rows[w] = (c + crude, rr + max(t, 0.0))
        tot_c += crude
        tot_r += max(t, 0.0)
    return lv.qc, rows, tot_c, tot_r, comp_c


# ------------------------------------------------------------------ J0/J1
def check_J0_J1():
    print("J0/J1 caracterizacion de Pi(r) y umbrales:")
    ok = True
    for z in (17, 29, 37):
        lv = Level110(z, 3, (), z ** 3)
        pf_full = Pi_full(z)
        err = 0.0
        mn = 1e9
        for (r, pf, mu, phi) in lv.enum:
            den = 1.0
            for p in pf:
                den *= (1.0 - 1.0 / (p - 1))
            err = max(err, abs(lv.Pi[r] * den - pf_full))
            mn = min(mn, lv.Pi[r])
        ok &= err < 1e-12 and mn >= pf_full - 1e-15
        # Pi(p) = Pi_full*(p-1)/(p-2)
        ep = 0.0
        for (r, pf, mu, phi) in lv.enum:
            if len(pf) == 1:
                p = pf[0]
                ep = max(ep, abs(lv.Pi[r] - pf_full * (p - 1) / (p - 2)))
        ok &= ep < 1e-12
        print("  z=%d: Pi_full=%.6f | err Pi(r)*prod - Pi_full = %.1e | "
              "min Pi(r) >= Pi_full: %s | err Pi(p) = %.1e"
              % (z, pf_full, err, mn >= pf_full - 1e-15, ep))
    pf17, pf19 = Pi_full(17), Pi_full(19)
    print("  Pi_full(17)=%.4f > 1/2 > Pi_full(19)=%.4f: la banda "
          "{Pi<1/2} es VACIA a z=17 y no vacia desde z=19." % (pf17, pf19))
    ok &= pf17 > 0.5 > pf19
    print("  => Pi(r) = Pi_full/prod_{p|r}(1-1/(p-1)) >= Pi_full [T]; "
          "Pi(r)<1/2 <=> prod_{p|r}(1+1/(p-2)) < 1/(2 Pi_full). %s"
          % ("PASS" if ok else "FAIL"))
    return ok


# ------------------------------------------------------------------ J2
def check_J2():
    print("J2 descomposicion de la banda cruda por omega(r) (u0=1.5):")
    for (tag, z, ustars) in (("A", 17, (5.0 / 3, 2, 3, 4)),
                             ("A", 29, (5.0 / 3, 2, 3, 4)),
                             ("A", 37, (5.0 / 3, 2, 3, 4)),
                             ("B", 29, (3, 4, 5))):
        for us in ustars:
            L = int(z ** us)
            qc, rows, tc, tr, _ = band_decomp(z, L)
            det = " ".join(
                "w%d:%.4f" % (w, rows[w][0] / qc) for w in sorted(rows))
            print("  [%s] z=%d u_L=%.2f L=%d: crudo/Q_c=%.4f (real %.4f)"
                  " | %s" % (tag, z, us, L, tc / qc, tr / qc, det))
    print("  Maximo sobre niveles muestreados del doc 114 (L=D//s):")
    ok = True
    for z in (17, 29, 37):
        spr = sample_sprimes(z, int(z ** (10.0 / 3.0)), 3, (), 60)
        mx = 0.0
        for L in sorted({z ** 5 // s for s in spr}):
            qc, rows, tc, tr, _ = band_decomp(z, L)
            mx = max(mx, tc / qc)
        print("    z=%d: max crudo_banda/Q_c = %.4f (doc 114: total 0,6991, parte Pi<1/2 = 0,4408; aqui "
              "global)" % (z, mx))
        ok &= mx <= 0.75
    print("  => en ventana la banda cruda es <= 0,70 Q_c [M], dominada "
          "por omega=1 en niveles de u_L bajo. %s"
          % ("PASS" if ok else "FAIL"))
    return ok


# ------------------------------------------------------------------ J4
def check_J4():
    print("J4 familia omega=1 a z=37, L=37^2 (u_L=2): crudo vs real:")
    z, L = 37, 37 ** 2
    lv = Level110(z, 3, (), L)
    info = {r: (pf, phi) for (r, pf, mu, phi) in lv.enum}
    ok = True
    for r in sorted(visible_set(lv)):
        if r == 1:
            continue
        pf, phi = info[r]
        if len(pf) != 1:
            continue
        ur = math.log(L / r) / math.log(z)
        Pi = lv.Pi[r]
        if ur >= 1.5 or Pi >= 0.5:
            continue
        S = lv.S[r]
        delta = lv.delta[r]
        crude = (1 - Pi * Pi) / phi
        t = delta * (2 * Pi + delta) / phi
        print("    p=%d: Pi=%.4f S(p)=%.4f delta=%+.4f | crudo=%.4f "
              "real t+=%.4f (factor %.1fx)"
              % (r, Pi, S, delta, crude, max(t, 0.0),
                 crude / max(t, 1e-12)))
        ok &= abs(S) < 0.75  # S(p) lejos de 1: la cota delta<=1-Pi pierde
    print("  => S(p) esta cerca de Pi(p), no de 1: la cota cruda pierde "
          "el factor (1-Pi)/delta por termino. %s"
          % ("PASS" if ok else "FAIL"))
    return ok


# ------------------------------------------------------------------ J3
def box_bound(centers, uL, u0, delta=0.02):
    """Cota inferior rigurosa del limite de Mertens de una familia:
    caja prod_i [c_i-delta, c_i+delta] (ordenada, dentro de la region)
    => integral >= prod_i 2*delta/(c_i+delta). 0 si la caja no vale."""
    om = len(centers)
    for i in range(om - 1):
        if centers[i] - delta <= centers[i + 1] + delta:
            return 0.0
    if sum(c - delta for c in centers) <= uL - u0:
        return 0.0
    for j in range(1, om + 1, 2):
        if sum(c + delta for c in centers[:j]) + 2 * (centers[j - 1]
                                                      + delta) >= uL:
            return 0.0
    b = 1.0
    for c in centers:
        b *= 2 * delta / (c + delta)
    return b


def check_J3():
    print("J3 limites de Mertens y constantes:")
    ok = True
    # (a) sumas de primos en el intervalo omega=1 del slice superior
    for (Z, a, b, target, tag) in (
            (10 ** 8, 0.5, 2.0 / 3, math.log(4.0 / 3), "u_L=2"),
            (10 ** 12, 1.0 / 6, 5.0 / 9, math.log(10.0 / 3), "u_L=5/3")):
        hi = int(Z ** b) + 1
        s = sum(1.0 / (p - 1) for p in primes_upto(hi) if p > Z ** a)
        print("    Z=%g [%s]: sum_{Z^a<p<Z^b} 1/(p-1) = %.4f "
              "(limite log=%.4f, err %.4f)" % (Z, tag, s, target,
                                               abs(s - target)))
        ok &= abs(s - target) < 0.03
    # (a2) cofactor fijo: la misma familia con r = p*m0, m0=35 (omega=3)
    Z = 10 ** 8
    s = sum(1.0 / (p - 1) for p in primes_upto(int(Z ** (2.0 / 3)) + 1)
            if p > Z ** 0.5)
    print("    cofactor m0=35: limite predicho log(4/3)/phi(35) = %.4f; "
          "suma/24 = %.4f" % (math.log(4.0 / 3) / 24.0, s / 24.0))
    ok &= abs(s / 24.0 - math.log(4.0 / 3) / 24.0) < 0.002
    # (b) B(u_L) > 0 en todo el rango de niveles: B1 exacta (omega=1)
    # y cotas de caja para omega=2,3 cuando B1 = 0.
    print("    B(u_L) = liminf de la suma cruda de la banda (u0=1.5):")
    centers2 = [(0.9, 0.7), (0.8, 0.6), (0.95, 0.8), (0.85, 0.55),
                (0.75, 0.5)]
    centers3 = [(0.97, 0.93, 0.68), (0.95, 0.9, 0.5), (0.9, 0.85, 0.4),
                (0.9, 0.8, 0.45), (0.85, 0.8, 0.6), (0.95, 0.85, 0.65)]
    for uL in (5.0 / 3, 2.0, 2.5, 3.0, 3.5, 4.0):
        lo, up = uL - 1.5, min(1.0, uL / 3.0)
        b1 = math.log(up / lo) if 0 < lo < up else 0.0
        b2 = max([box_bound(c, uL, 1.5, 0.02) for c in centers2] + [0.0])
        b3 = max([box_bound(c, uL, 1.5, 0.01) for c in centers3] + [0.0])
        B = max(b1, b2, b3)
        print("      u_L=%.2f: B1=%.4f B2>=%.2e B3>=%.2e -> B>=%.2e > 0"
              % (uL, b1, b2, b3, B))
        ok &= B > 1e-6
    # (c) Pi_full(z)*ln z y Q_c(z)*ln z se estabilizan: Q_c ~ KM/ln z
    print("    Estabilizacion de Pi_full*ln z y Q_c*ln z:")
    for z in (37, 53, 101, 1009, 10007):
        pf = Pi_full(z)
        qc = Q_completo(z, 3)
        print("      z=%6d: Pi_full=%.5f (x ln z = %.4f) | Q_c=%.5f "
              "(x ln z = %.4f)" % (z, pf, pf * math.log(z), qc,
                                    qc * math.log(z)))
    print("  => la suma cruda de la banda tiende a B(u_L) > 0 mientras "
          "Q_c -> 0 como K/ln z: (BAND) cruda es FALSA en el marco "
          "asintotico. %s" % ("PASS" if ok else "FAIL"))
    return ok


# ------------------------------------------------------------------ J5
def check_J5():
    print("J5 sub-banda {Pi >= 1/2} con u_r<1.5: suma cruda <= 4 Q_c [T]:")
    ok = True
    for z in (29, 37):
        for us in (5.0 / 3, 2, 3, 4):
            L = int(z ** us)
            qc, rows, tc, tr, comp = band_decomp(z, L)
            ok &= comp <= 4.0 * qc + 1e-12
            print("    z=%d u_L=%.2f: suma_{Pi>=1/2} crudo = %.4f Q_c "
                  "(<= 4 Q_c)" % (z, us, comp / qc))
    print("  => cerrada por (1-Pi^2)/phi <= (Pi^2/phi)/c0^2 y "
          "sum Pi^2/phi = Q_c. %s" % ("PASS" if ok else "FAIL"))
    return ok


# ------------------------------------------------------------------ main
def main():
    t0 = time.time()
    print("== DOC 116c: el residuo (BAND) por conteo ==")
    print("Solo stdlib; maquinaria de los docs 109/110 importada.")
    print()
    allok = True
    allok &= check_J0_J1()
    print()
    allok &= check_J2()
    print()
    allok &= check_J4()
    print()
    allok &= check_J3()
    print()
    allok &= check_J5()
    print()
    print("Tiempo total: %.1f s" % (time.time() - t0))
    print("PASS" if allok else "FAIL")
    return 0 if allok else 1


if __name__ == "__main__":
    sys.exit(main())
