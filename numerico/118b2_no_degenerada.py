#!/usr/bin/env python3
# 118b2_no_degenerada.py  (Doc 118b2)
# Parte no degenerada de (BAND-LOW) y reconciliacion con la refutacion
# de (REG-Pi) del doc 118a.
#   J0: aritmetica de Mertens de la reconciliacion (familia de 118a,
#       hermanas omega=3 y escalera, en la capa baja).
#   J1: medicion en marcos A/B (maquinaria del doc 117b reproducida):
#       sumas (Shat-Pi)+/phi y Pi/phi sobre la banda no degenerada,
#       partidas por u_r < 1 y 1 <= u_r < 3/2.
#   J2: la familia de 118a en ventana (niveles representativos).
# Solo stdlib. Las funciones in_Vhat/Shat_of/visible_set reproducen
# verbatim las del doc 117b (que a su vez reproduce in_Vhat del 113).
import importlib.util
import math
import sys
import time

spec = importlib.util.spec_from_file_location(
    "d110", "110_segundo_momento_defectos.py")
d110 = importlib.util.module_from_spec(spec)
spec.loader.exec_module(d110)

Level110 = d110.Level110
sub_divisors = d110.sub_divisors
sample_sprimes = d110.sample_sprimes
primes_upto = d110.primes_upto


def in_Vhat(pf, L, r):
    """Doc 113, verbatim: pertenencia a V^+(L/r; P\\r) en enteros."""
    if not pf:
        return True
    if pf[0] * r >= L:
        return False
    prod = 1
    for i, p in enumerate(pf):
        m = i + 1
        prod *= p
        if m % 2 == 1 and r * p * p * prod >= L:
            return False
    return True


def visible_set(lv):
    vis = set()
    for (d, pf, mu, phi) in lv.V:
        for (r, rmu, rphi) in sub_divisors(pf):
            vis.add(r)
    return vis


def Shat_of(lv, r):
    """Shat(r) = sum_{d' en Vhat^(r)} mu(d')/phi(d'); tambien |Vhat|."""
    L = lv.L
    tot = 0.0
    cnt = 0
    for (dp, dpf, dmu, dphi) in lv.enum:
        if math.gcd(dp, r) != 1:
            continue
        if (dp * r < L) and in_Vhat(dpf, L, r):
            tot += dmu / dphi
            cnt += 1
    return tot, cnt


FRAMES = [("A", 17, 17 ** 5, int(17 ** (10.0 / 3.0))),
          ("A", 29, 29 ** 5, int(29 ** (10.0 / 3.0))),
          ("A", 37, 37 ** 5, int(37 ** (10.0 / 3.0))),
          ("B", 29, 29 ** 6, 29 ** 3)]
U0 = 1.5  # corte de la banda (doc 117b)


def adm_primes(Z):
    """Primos de P (impares, != 3) hasta Z."""
    return [p for p in primes_upto(Z) if p != 2 and p != 3]


# ------------------------------------------------------------------ J0
def check_J0():
    print("J0 reconciliacion [T/M]: donde vive la masa de 118a y las "
          "familias hermanas de la capa baja:")
    ok = True
    # (a) suma de Mertens de la ventana x_p de la familia F
    tgt = math.log(12.0 / 11.0)
    for Z in (10 ** 6, 10 ** 8):
        lz = math.log(Z)
        s = sum(1.0 / (p - 2) for p in adm_primes(Z)
                if 11.0 / 20.0 < math.log(p) / lz < 3.0 / 5.0)
        print("  (a) z=10^%d: sum_{x_p en (11/20,3/5)} 1/(p-2) = %.4f "
              "(log(12/11) = %.4f)" % (round(math.log10(Z)), s, tgt))
        ok &= abs(s - tgt) < 0.02
    # (b) fraccion de la masa interna de F con u_r >= 1 (u_L = 2):
    #     x_q <= 1 - x_p; se evalua en el primo mayor de la ventana.
    for Z in (10 ** 4, 10 ** 6, 10 ** 8):
        lz = math.log(Z)
        win = [p for p in adm_primes(Z)
               if 11.0 / 20.0 < math.log(p) / lz < 3.0 / 5.0]
        p = win[-1]
        xp = math.log(p) / lz
        tot = sum(1.0 / (q - 2) for q in adm_primes(p - 1) if q >= 5)
        hi = sum(1.0 / (q - 2) for q in adm_primes(p - 1)
                 if q >= 5 and math.log(q) / lz <= 1.0 - xp)
        print("  (b) z=10^%d: p=%d x_p=%.3f | fraccion de la masa de "
              "cofactores con u_r>=1: %.3f (tiende a 1)"
              % (round(math.log10(Z)), p, xp, hi / tot))
        ok &= hi / tot > 0.5
    # (c) hermana omega=3 en la capa baja (u_L=2, u_0*=1): subcaja
    #     x_1 en (0.56,0.65), x_2 en (0.50,0.55): x_2 < x_1,
    #     x_1+x_2 > 1.06 > u_L-u_0*, x_1 < 2/3: masa > 0.
    Z = 10 ** 6
    lz = math.log(Z)
    box = [(0.56, 0.65), (0.50, 0.55)]
    masses = []
    for (a, b) in box:
        masses.append(sum(1.0 / (p - 2) for p in adm_primes(Z)
                          if a < math.log(p) / lz < b))
    print("  (c) subcaja omega=3 (u_L=2, u_0*=1): masas de Mertens "
          "por coordenada = %.3f, %.3f (>0: region de medida "
          "positiva; x_2<x_1, x_1+x_2>1.06, x_1<2/3)"
          % (masses[0], masses[1]))
    ok &= all(m > 0.05 for m in masses)
    # (d) escalera en u_L=2 para u_0*=1/2: pares (2/3,2/3),(2/9,2/9)
    #     (omega=5 con el primo pequeno libre): suma central 16/9,
    #     u_r central 2/9 < 1/2; condiciones de Nathanson en el centro
    #     con margen: se usan los valores x*(1-1e-3).
    xs = [2.0 / 3.0, 2.0 / 3.0, 2.0 / 9.0, 2.0 / 9.0]
    f = 1.0 - 1e-3
    xs = [x * f for x in xs]
    nath = True
    for m in (1, 3):
        lhs = sum(xs[:m - 1]) + 3.0 * xs[m - 1]
        nath &= lhs < 2.0
    ssum = sum(xs)
    print("  (d) escalera omega=5 (u_L=2, u_0*=1/2): suma de las "
          "coordenadas fijas = %.4f > u_L-u_0* = 1.5 | u_r = %.4f "
          "< 0.5 | Nathanson m=1,3 estrictas: %s"
          % (ssum, 2.0 - ssum, nath))
    ok &= nath and ssum > 1.5 and (2.0 - ssum) < 0.5
    print("  => %s" % ("PASS" if ok else "FAIL"))
    return ok


# ------------------------------------------------------------------ J1
def frame_levels(z, A, excl, D, y, nsamp=25):
    spr = sample_sprimes(z, y, A, excl, nsamp)
    levels = sorted({D // s for s in spr})
    return [Level110(z, A, excl, L) for L in levels]


def check_J1():
    print("J1 medicion en marcos A/B: banda no degenerada partida en "
          "u_r<1 (capa baja de referencia) y 1<=u_r<3/2:")
    ok = True
    for (tag, z, D, y) in FRAMES:
        levels = frame_levels(z, 3, (), D, y)
        lnz = math.log(z)
        # maximos sobre niveles de las sumas / Q_c, por celda
        best = {"lo_nd": 0.0, "hi_nd": 0.0, "lo_pi": 0.0, "hi_pi": 0.0,
                "deg": 0.0}
        cnt = {"lo_nd": 0, "hi_nd": 0, "deg": 0}
        relmax_lo = 0.0
        char_ok = True
        for lv in levels:
            info = {r: (pf, phi) for (r, pf, mu, phi) in lv.enum}
            sums = {"lo_nd": 0.0, "hi_nd": 0.0, "lo_pi": 0.0,
                    "hi_pi": 0.0, "deg": 0.0}
            for r in visible_set(lv):
                if r == 1:
                    continue
                pf, phi = info[r]
                ur = math.log(lv.L / r) / lnz
                if ur >= U0:
                    continue
                Pi = lv.Pi[r]
                if Pi >= 0.5:
                    continue
                sh, nhat = Shat_of(lv, r)
                q0 = next(p for p in lv.allow if p not in pf)
                char_ok &= (nhat == 1) == (q0 ** 3 * r >= lv.L)
                dif = sh - Pi
                if nhat == 1:
                    sums["deg"] += (1.0 - Pi) / phi
                    cnt["deg"] += 1
                elif ur < 1.0:
                    cnt["lo_nd"] += 1
                    if dif > 0:
                        sums["lo_nd"] += dif / phi
                        relmax_lo = max(relmax_lo, dif / Pi)
                    sums["lo_pi"] += Pi / phi
                else:
                    cnt["hi_nd"] += 1
                    if dif > 0:
                        sums["hi_nd"] += dif / phi
                    sums["hi_pi"] += Pi / phi
            for k in best:
                best[k] = max(best[k], sums[k] / lv.qc)
        print("  [%s] z=%d (%d niveles): max sum/Q_c | no-deg u_r<1: "
              "(Shat-Pi)+ %.4f, Pi %.4f | no-deg 1<=u_r<1.5: (Shat-Pi)+ "
              "%.4f, Pi %.4f | deg: %.4f | instancias lo/hi/deg = "
              "%d/%d/%d | max (Shat-Pi)+/Pi en capa baja no-deg = %.3f "
              "| carac. deg %s"
              % (tag, z, len(levels), best["lo_nd"], best["lo_pi"],
                 best["hi_nd"], best["hi_pi"], best["deg"],
                 cnt["lo_nd"], cnt["hi_nd"], cnt["deg"], relmax_lo,
                 "OK" if char_ok else "FALLA"))
        ok &= char_ok
    print("  => %s" % ("PASS" if ok else "FAIL"))
    return ok


# ------------------------------------------------------------------ J2
def check_J2():
    print("J2 la familia de 118a en ventana (omega=2, x_p en "
          "(11/20, 3/5)):")
    ok = True
    for z in (29, 37, 53):
        allow = [p for p in primes_upto(z - 1) if p != 2 and p != 3]
        lnz = math.log(z)
        for uu in (5.0 / 3.0, 2.0):
            L = int(z ** uu)
            lv = Level110(z, 3, (), L)
            info = {r: (pf, phi) for (r, pf, mu, phi) in lv.enum}
            n = 0
            urs = []
            s_ex = s_pi = 0.0
            for r in visible_set(lv):
                if r == 1:
                    continue
                pf, phi = info[r]
                if len(pf) != 2:
                    continue
                xp = math.log(pf[0]) / lnz
                if not (11.0 / 20.0 < xp < 3.0 / 5.0):
                    continue
                ur = math.log(L / r) / lnz
                if ur >= U0:
                    continue
                Pi = lv.Pi[r]
                if Pi >= 0.5:
                    continue
                sh, nhat = Shat_of(lv, r)
                n += 1
                urs.append(ur)
                if sh - Pi > 0:
                    s_ex += (sh - Pi) / phi
                s_pi += Pi / phi
            rng = ("u_r en [%.3f, %.3f]" % (min(urs), max(urs))
                   if urs else "sin instancias")
            print("  z=%d u_L=%.2f L=%d: %d instancias de la familia "
                  "en banda | %s | sum (Shat-Pi)+/phi/Q_c = %.4f | "
                  "sum Pi/phi/Q_c = %.4f"
                  % (z, uu, L, n, rng, s_ex / lv.qc, s_pi / lv.qc))
    print("  => PASS (reporte)")
    return ok


# ----------------------------------------------------------------- J2b
def check_J2b():
    print("J2b ventana ampliada de 118a (omega=2, x_p en (1/2, 2/3), "
          "como su J4):")
    for z in (53,):
        lnz = math.log(z)
        for uu in (2.0,):
            L = int(z ** uu)
            lv = Level110(z, 3, (), L)
            info = {r: (pf, phi) for (r, pf, mu, phi) in lv.enum}
            fam = []
            for r in visible_set(lv):
                if r == 1:
                    continue
                pf, phi = info[r]
                if len(pf) != 2:
                    continue
                xp = math.log(pf[0]) / lnz
                if 0.5 < xp < 2.0 / 3.0:
                    fam.append(r)
            inband = []
            for r in fam:
                pf, phi = info[r]
                ur = math.log(L / r) / lnz
                Pi = lv.Pi[r]
                if ur < U0 and Pi < 0.5:
                    sh, nhat = Shat_of(lv, r)
                    inband.append((r, ur, nhat))
            print("  z=%d u_L=%.2f L=%d: %d elementos en la ventana "
                  "ampliada, %d en banda: %s"
                  % (z, uu, L, len(fam), len(inband),
                     "; ".join("r=%d u_r=%.3f %s" % (r, ur,
                                 "degenerado" if nh == 1
                                 else "no-deg")
                               for (r, ur, nh) in inband)))
    print("  => PASS (reporte)")
    return True


# ----------------------------------------------------------------- main
def main():
    t0 = time.time()
    print("== DOC 118b2: parte no degenerada de (BAND-LOW) ==")
    print()
    allok = True
    allok &= check_J0()
    print()
    allok &= check_J1()
    print()
    allok &= check_J2()
    print()
    allok &= check_J2b()
    print()
    print("Tiempo total: %.1f s" % (time.time() - t0))
    print("PASS" if allok else "FAIL")
    return 0 if allok else 1


if __name__ == "__main__":
    sys.exit(main())
