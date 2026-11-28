#!/usr/bin/env python3
# 118b1_degenerada.py  (Doc 118b1)
# Subcapa degenerada de (BAND-LOW): r visible en banda con Vhat^(r)={1}
# sii q0^3*r >= L, q0 = min(P \ primos(r)) [T, doc 117b sec. 3.2].
# Chequeos:
#   J0  caracterizacion: degenerado sii Shat cuenta 1 sii q0^3 r >= L,
#       sobre las instancias de banda de los marcos A/B.
#   J1  Teorema A: ningun degenerado visible con omega(r) <= 2;
#       p1(r)^3 < L para todo r visible.
#   J2  suma degenerada por estrato q0 contra Q_c (max por marco).
#   J3  constantes de Euler: E(z), kappa0(z), cola q0 >= q0-daga.
#   J4  obstruccion en ventana: la misma suma SIN visibilidad (q0=5).
# Solo stdlib; importa la maquinaria de los docs 110/117b sin tocarla.
import importlib.util
import math
import sys
import time

spec = importlib.util.spec_from_file_location(
    "d117b", "117b_banda_delta.py")
d = importlib.util.module_from_spec(spec)
spec.loader.exec_module(d)

Level110 = d.Level110
sample_sprimes = d.sample_sprimes
visible_set = d.visible_set
Shat_of = d.Shat_of
sub_divisors = d.sub_divisors

FRAMES = [("A", 17, 17 ** 5, int(17 ** (10.0 / 3.0))),
          ("A", 29, 29 ** 5, int(29 ** (10.0 / 3.0))),
          ("A", 37, 37 ** 5, int(37 ** (10.0 / 3.0))),
          ("B", 29, 29 ** 6, 29 ** 3)]
U0 = 1.5


def frame_levels(z, D, y):
    spr = sample_sprimes(z, y, 3, ())
    return [Level110(z, 3, (), L) for L in sorted({D // s for s in spr})]


def q0_of(lv, pf):
    s = set(pf)
    for p in lv.allow:
        if p not in s:
            return p
    return None


def band_degen(lv, z):
    """(lista de (r, pf, phi, q0) degenerados en banda, n_banda)."""
    info = {r: (pf, phi) for (r, pf, mu, phi) in lv.enum}
    vis = visible_set(lv)
    out = []
    nb = 0
    for r in vis:
        pf, phi = info[r]
        if lv.Pi[r] >= 0.5:
            continue
        if math.log(lv.L / r) / math.log(z) >= U0:
            continue
        nb += 1
        q0 = q0_of(lv, pf)
        if q0 is None or q0 ** 3 * r >= lv.L:
            out.append((r, pf, phi, q0))
    return out, nb


def check_J0_J1_J2():
    print("J0/J1/J2 caracterizacion, Teorema A (condicional) y sumas:")
    ok = True
    for (tag, z, D, y) in FRAMES:
        levels = frame_levels(z, D, y)
        mx_tot = 0.0
        mx_byq0 = {}
        nband = 0
        ndeg = 0
        bad_char = 0
        bad_sharp = 0
        bad_p1 = 0
        hist = {}
        for lv in levels:
            info = {r: (pf, phi) for (r, pf, mu, phi) in lv.enum}
            for r in visible_set(lv):
                pf, _ = info[r]
                if pf and pf[0] ** 3 >= lv.L:
                    bad_p1 += 1
            deg, nb = band_degen(lv, z)
            nband += nb
            ndeg += len(deg)
            tot = 0.0
            per = {}
            for (r, pf, phi, q0) in deg:
                # J0: cruce con Shat (conteo de Vhat) en la banda
                _, cnt = Shat_of(lv, r)
                if cnt != 1:
                    bad_char += 1
                hist[len(pf)] = hist.get(len(pf), 0) + 1
                # J1: omega<=2 degenerado visible => L < q0^9 (agudeza)
                if len(pf) <= 2 and not (q0 and lv.L < q0 ** 9):
                    bad_sharp += 1
                tot += 1.0 / phi
                per[q0] = per.get(q0, 0.0) + 1.0 / phi
            mx_tot = max(mx_tot, tot / lv.qc if lv.qc else 0.0)
            for q0, v in per.items():
                mx_byq0[q0] = max(mx_byq0.get(q0, 0.0), v / lv.qc)
        ok &= bad_char == 0 and bad_sharp == 0 and bad_p1 == 0
        estr = " ".join("q0=%s:%.3f" % (q, v)
                        for q, v in sorted(mx_byq0.items(),
                                           key=lambda kv: str(kv[0])))
        print("  [%s] z=%d: banda=%d deg=%d hist-omega=%s | "
              "max sum_deg/Q_c=%.4f | estratos: %s | viol char=%d "
              "agudeza=%d p1=%d"
              % (tag, z, nband, ndeg, hist, mx_tot, estr,
                 bad_char, bad_sharp, bad_p1))
    print("  => %s" % ("PASS" if ok else "FAIL"))
    return ok


def check_J3():
    print("J3 constantes de Euler y cola q0 >= q0-daga:")
    ok = True
    for z in (17, 29, 37, 53):
        allow = [p for p in range(5, z + 1)
                 if all(p % q for q in range(2, int(p ** 0.5) + 1))]
        E = 1.0
        kap = 1.0
        for p in allow:
            E *= p / (p - 1)
            kap *= (1.0 - 1.0 / (p - 1) ** 2)
        # P(q0) y q0-daga para este z
        P = 1
        dag = None
        for q in allow:
            if P >= math.log(z) ** 3:
                dag = q
                break
            P *= q
        tail = 2.0 * E / (math.log(z) ** 3)
        print("  z=%d: E(z)=%.4f E/log z=%.4f kappa0=%.4f "
              "q0-daga=%s cola<=%.5f" % (z, E, E / math.log(z), kap,
                                         dag, tail))
        ok &= E / math.log(z) <= 0.66 and kap <= 0.90
    print("  => %s" % ("PASS" if ok else "FAIL"))
    return ok


def check_J4():
    print("J4 obstruccion de ruta: total sin visibilidad y cota Euler "
          "del estrato q0=5:")
    ok = True
    for (tag, z, D, y) in FRAMES:
        if z == 17:
            continue
        levels = frame_levels(z, D, y)
        qc = levels[0].qc
        mx_all = 0.0
        for lv in levels:
            tot = 0.0
            for (r, pf, mu, phi) in lv.enum:
                if lv.Pi[r] >= 0.5:
                    continue
                if math.log(lv.L / r) / math.log(z) >= U0:
                    continue
                q0 = q0_of(lv, pf)
                if q0 is None or q0 ** 3 * r < lv.L:
                    continue
                tot += 1.0 / phi
            mx_all = max(mx_all, tot / qc)
        E = 1.0
        for p in levels[0].allow:
            E *= p / (p - 1)
        route5 = E * 0.8 / qc
        print("  [%s] z=%d: max sum_ALLdeg/Q_c=%.4f | cota Euler q0=5 "
              "= %.2f Q_c (la ruta pierde (log z)^2)"
              % (tag, z, mx_all, route5))
        ok &= route5 > 1.0
    print("  => %s" % ("PASS" if ok else "FAIL"))
    return ok


def main():
    t0 = time.time()
    print("== DOC 118b1: subcapa degenerada de (BAND-LOW) ==")
    ok = True
    ok &= check_J0_J1_J2()
    print()
    ok &= check_J3()
    print()
    ok &= check_J4()
    print()
    print("Tiempo total: %.1f s" % (time.time() - t0))
    print("PASS" if ok else "FAIL")
    return 0 if ok else 1


if __name__ == "__main__":
    sys.exit(main())
