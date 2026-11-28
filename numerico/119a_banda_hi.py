#!/usr/bin/env python3
# 119a_banda_hi.py  (Doc 119a)
# (BAND-MOD-HI): la capa alta {u_r >= u0*} de (BAND-MOD) atacada DIRECTO
# sobre Shat, sin el mayorante Pi/phi (muerto en el doc 118).
# Mide, en marcos A/B y estres z=53:
#   sum_{r en banda} (Shat-Pi)^+/phi(r) / Q_c, por rebanadas de u_r y
#   por omega(r); la familia F del doc 118a (r=pq, x_p en (0.55,0.6))
#   con su peso REAL (Shat-Pi)^+ contra su masa Pi; el perfil de la
#   razon (Shat-Pi)/Pi por rebanada de u_r; y los momentos L2
#   sum (Shat-Pi)^2/phi y sum (Shat-Pi)^2/(Pi*phi) contra Q_c.
# Banda: r visible, r>1, u_r < 1.5, Pi(r) < 1/2 (docs 114/117b).
# Solo stdlib; importa la maquinaria del doc 117b sin tocarla.
import importlib.util
import math
import sys
import time

spec = importlib.util.spec_from_file_location(
    "d117b", "117b_banda_delta.py")
d117 = importlib.util.module_from_spec(spec)
spec.loader.exec_module(d117)

Level110 = d117.Level110
Shat_of = d117.Shat_of
That_of = d117.That_of
visible_set = d117.visible_set
sample_sprimes = d117.sample_sprimes

SLICES = [(0.0, 0.5), (0.5, 1.0), (1.0, 1.25), (1.25, 1.4),
          (1.4, 1.45), (1.45, 1.5)]


def slice_idx(u):
    for i, (a, b) in enumerate(SLICES):
        if a <= u < b:
            return i
    return None


def pick_levels(z, D, y, cap, near_u=None):
    spr = sample_sprimes(z, y, 3, (), 60)
    levels = sorted({D // s for s in spr})
    if len(levels) <= cap:
        return levels
    if near_u is not None:
        lnz = math.log(z)
        levels.sort(key=lambda L: abs(math.log(L) / lnz - near_u))
        chosen = levels[:cap]
        return sorted(chosen)
    step = len(levels) / float(cap)
    return sorted({levels[int(i * step)] for i in range(cap)})


# ------------------------------------------------------------------ K0
def check_K0():
    print("K0 identidad Shat + That = Pi (maquinaria 117b):")
    lv = Level110(29, 3, (), 29 ** 5 // 37)
    vis = [r for r in visible_set(lv) if r > 1][:14]
    err = 0.0
    for r in vis:
        sh, _ = Shat_of(lv, r)
        th = That_of(lv, r)
        err = max(err, abs(sh + th - lv.Pi[r]))
    print("  14 r visibles, z=29: max |Shat+That-Pi| = %.2e %s"
          % (err, "PASS" if err < 1e-9 else "FAIL"))
    return err < 1e-9


# ------------------------------------------------------------- barrido
def sweep():
    print("K1 suma real de la banda por rebanada de u_r (max sobre "
          "niveles, / Q_c):")
    frames = [("A", 17, 17 ** 5, int(17 ** (10.0 / 3.0)), 60, None),
              ("A", 29, 29 ** 5, int(29 ** (10.0 / 3.0)), 60, None),
              ("A", 37, 37 ** 5, int(37 ** (10.0 / 3.0)), 60, None),
              ("B", 29, 29 ** 6, 29 ** 3, 60, None),
              ("A-estres", 53, 53 ** 5, int(53 ** (10.0 / 3.0)), 4, 2.0)]
    allok = True
    for (tag, z, D, y, cap, near_u) in frames:
        levels = pick_levels(z, D, y, cap, near_u)
        lnz = math.log(z)
        mx_slice = [0.0] * len(SLICES)
        mx_om = {1: 0.0, 2: 0.0, 3: 0.0, 4: 0.0}
        mx_tot = mx_hi100 = mx_hi125 = mx_hi140 = 0.0
        mx_maj = mx_l2 = mx_l2p = 0.0
        mx_F = mx_Fmaj = 0.0
        ratio_prof = {}
        for L in levels:
            lv = Level110(z, 3, (), L)
            qc = lv.qc
            info = {r: pf for (r, pf, mu, phi) in lv.enum}
            sl = [0.0] * len(SLICES)
            om = {1: 0.0, 2: 0.0, 3: 0.0, 4: 0.0}
            tot = maj = l2 = l2p = 0.0
            fsum = fmaj = fnum = 0.0
            prof = {}
            for r in visible_set(lv):
                if r == 1:
                    continue
                ur = math.log(L / r) / lnz
                if ur >= 1.5:
                    continue
                Pi = lv.Pi[r]
                if Pi >= 0.5:
                    continue
                pf = info[r]
                phi = 1.0
                for p in pf:
                    phi *= (p - 1)
                sh, _ = Shat_of(lv, r)
                d = sh - Pi
                dp = d if d > 0.0 else 0.0
                w = dp / phi
                tot += w
                maj += Pi / phi
                l2 += d * d / phi
                l2p += d * d / (Pi * phi)
                si = slice_idx(ur)
                if si is not None:
                    sl[si] += w
                om[min(len(pf), 4)] += w
                key = si
                a = prof.setdefault(key, [0.0, 0.0, 0.0])
                a[0] += dp / phi
                a[1] += Pi / phi
                if Pi > 0:
                    a[2] = max(a[2], d / Pi)
                if len(pf) == 2:
                    pmax = pf[0]
                    xp = math.log(pmax) / lnz
                    if 0.55 < xp < 0.60 and pf[1] >= 5:
                        fsum += w
                        fmaj += Pi / phi
                        fnum += 1
            for i in range(len(SLICES)):
                mx_slice[i] = max(mx_slice[i], sl[i] / qc)
                a = prof.get(i)
                if a and a[1] > 0:
                    cur = ratio_prof.get(i, (0.0, 0.0))
                    ratio_prof[i] = (max(cur[0], a[0] / a[1]),
                                     max(cur[1], a[2]))
            for k in om:
                mx_om[k] = max(mx_om[k], om[k] / qc)
            mx_tot = max(mx_tot, tot / qc)
            mx_hi100 = max(mx_hi100,
                           sum(sl[i] for i in range(2, 6)) / qc)
            mx_hi125 = max(mx_hi125,
                           sum(sl[i] for i in range(3, 6)) / qc)
            mx_hi140 = max(mx_hi140,
                           sum(sl[i] for i in range(4, 6)) / qc)
            mx_maj = max(mx_maj, maj / qc)
            mx_l2 = max(mx_l2, l2 / qc)
            mx_l2p = max(mx_l2p, l2p / qc)
            mx_F = max(mx_F, fsum / qc)
            mx_Fmaj = max(mx_Fmaj, fmaj / qc)
        print("  [%s] z=%d, %d niveles: banda total=%.4f | capa "
              "u_r>=1.00: %.4f | >=1.25: %.4f | >=1.40: %.4f"
              % (tag, z, len(levels), mx_tot, mx_hi100, mx_hi125,
                 mx_hi140))
        print("      rebanadas u_r: " + " ".join(
            "[%.2f,%.2f)=%.4f" % (a, b, mx_slice[i])
            for i, (a, b) in enumerate(SLICES)))
        print("      omega: 1:%.4f 2:%.4f 3:%.4f >=4:%.4f | mayorante "
              "Pi/phi=%.4f | L2: sum d^2/phi=%.4f, sum d^2/(Pi phi)="
              "%.4f" % (mx_om[1], mx_om[2], mx_om[3], mx_om[4],
                        mx_maj, mx_l2, mx_l2p))
        print("      familia F (pq, x_p en (0.55,0.6)): real=%.4f | "
              "su masa Pi=%.4f | cociente=%.3f"
              % (mx_F, mx_Fmaj, (mx_F / mx_Fmaj) if mx_Fmaj > 0
                 else float("nan")))
        for i in sorted(ratio_prof):
            a, b = SLICES[i]
            wm, mx = ratio_prof[i]
            print("      perfil (Shat-Pi)/Pi en [%.2f,%.2f): media "
                  "ponderada=%.4f max=%.4f" % (a, b, wm, mx))
        allok &= mx_tot < 50.0
    print("  => tablas completas y finitas. %s"
          % ("PASS" if allok else "FAIL"))
    return allok


def main():
    t0 = time.time()
    print("== DOC 119a: (BAND-MOD-HI) directo sobre Shat ==")
    print("Solo stdlib. Banda: u_r<1.5, Pi<1/2, r visible>1.")
    print()
    ok = True
    ok &= check_K0()
    print()
    ok &= sweep()
    print()
    print("Tiempo total: %.1f s" % (time.time() - t0))
    print("PASS" if ok else "FAIL")
    return 0 if ok else 1


if __name__ == "__main__":
    sys.exit(main())
