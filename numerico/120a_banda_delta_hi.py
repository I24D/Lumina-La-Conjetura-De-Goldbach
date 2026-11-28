#!/usr/bin/env python3
# 120a_banda_delta_hi.py  (Doc 120a)
# (BAND-delta-HI): la capa alta {u_r >= u0*} de la banda, directa sobre
# delta(r) = (Shat(r)-Pi(r)) + E(r). Pregunta decisiva: E(r) <= 0
# cancela a (Shat-Pi)^+ en la capa alta, r a r o solo en suma?
#   K0  identidades: delta = (Shat-Pi) + E; E = sum_Add a - sum_Rem a
#       (lados del doc 114); E <= 0 en la capa (DESP-SIGN en ventana).
#   K1  medicion decisiva en marcos A/B y estres z=53: por nivel,
#       sumas de capa alta (u_r>=1) y sub-rebanada (u_r>=1.4) de
#       delta^+/phi, (Shat-Pi)^+/phi, |E|/phi contra Q_c; separando
#       degenerados (Vhat={1}, E=0 exacto) de no degenerados;
#       cociente |E|/(Shat-Pi)^+ por instancia; familia F (r=pq,
#       x_p en (11/20,3/5)) en banda y fuera de banda.
#   K2  estructura: E_Add vs E_Rem sobre la capa alta no degenerada
#       (quien cancela), y exceso no cubierto ((Shat-Pi)^+ - |E|)^+.
# Solo stdlib; importa la maquinaria de los docs 110/117b sin tocarla;
# lados() reproduce la clasificacion del doc 114 con in_Vhat del 117b.
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
visible_set = d117.visible_set
sample_sprimes = d117.sample_sprimes
in_Vhat = d117.in_Vhat

U0 = 1.5          # banda: u_r < 3/2 (docs 114/117b)
U0STAR = 1.0      # capa alta medida: u_r >= 1 (perfil del doc 119a)
U14 = 1.4         # sub-rebanada donde se concentra la masa de F


def lados(lv, r):
    """Doc 114: clasifica cada d' coprimo con r en Add/Rem/Both/Neither.
    Add = V^(r) \\ Vhat^(r); Rem = Vhat^(r) \\ V^(r)."""
    L = lv.L
    Vset = lv.Vset
    out = {}
    for (dp, dpf, dmu, dphi) in lv.enum:
        if math.gcd(dp, r) != 1:
            continue
        inhat = (dp * r < L) and in_Vhat(dpf, L, r)
        inv = (r * dp) in Vset
        out[dp] = (("Add" if inv and not inhat else
                    "Rem" if inhat and not inv else
                    "Both" if inv and inhat else "Neither"),
                   dpf, dmu, dphi)
    return out


def phi_of(pf):
    f = 1
    for p in pf:
        f *= (p - 1)
    return f


def frame_levels(z, D, y, cap=None, near_u=None):
    spr = sample_sprimes(z, y, 3, (), 60)
    levels = sorted({D // s for s in spr})
    if cap is not None and len(levels) > cap:
        lnz = math.log(z)
        levels.sort(key=lambda L: abs(math.log(L) / lnz - near_u))
        levels = sorted(levels[:cap])
    return levels


FRAMES = [("A", 17, 17 ** 5, int(17 ** (10.0 / 3.0)), None, None),
          ("A", 29, 29 ** 5, int(29 ** (10.0 / 3.0)), None, None),
          ("A", 37, 37 ** 5, int(37 ** (10.0 / 3.0)), None, None),
          ("B", 29, 29 ** 6, 29 ** 3, None, None),
          ("A-estres", 53, 53 ** 5, int(53 ** (10.0 / 3.0)), 4, 2.0)]


def band_instances(lv, L, lnz):
    """Genera (r, pf, phi, ur, Pi, Shat, nhat, dif, E, delta) de banda."""
    info = {r: pf for (r, pf, mu, phi) in lv.enum}
    for r in visible_set(lv):
        if r == 1:
            continue
        pf = info[r]
        ur = math.log(L / r) / lnz
        if ur >= U0:
            continue
        Pi = lv.Pi[r]
        if Pi >= 0.5:
            continue
        sh, nhat = Shat_of(lv, r)
        dif = sh - Pi
        E = lv.S[r] - sh
        yield (r, pf, phi_of(pf), ur, Pi, sh, nhat, dif, E, lv.delta[r])


# ------------------------------------------------------------------ K0
def check_K0():
    print("K0 identidades en la capa alta (z=29, dos niveles):")
    ok = True
    lv = None
    for L in (29 ** 5 // 37, 29 ** 5 // 53):
        lv = Level110(29, 3, (), L)
        lnz = math.log(29)
        e_id = e_flow = 0.0
        npos = 0
        for (r, pf, phi, ur, Pi, sh, nhat, dif, E, dl) in \
                band_instances(lv, L, lnz):
            if ur < U0STAR:
                continue
            e_id = max(e_id, abs(dl - (dif + E)))
            ld = lados(lv, r)
            Ef = sum(v[2] / v[3] for v in ld.values() if v[0] == "Add") \
                - sum(v[2] / v[3] for v in ld.values() if v[0] == "Rem")
            e_flow = max(e_flow, abs(E - Ef))
            if E > 1e-12:
                npos += 1
        print("  L=%d: |delta-((Shat-Pi)+E)| max = %.2e | "
              "|E - (sum_Add a - sum_Rem a)| max = %.2e | E>0: %d"
              % (L, e_id, e_flow, npos))
        ok &= e_id < 1e-12 and e_flow < 1e-12 and npos == 0
    print("  => %s" % ("PASS" if ok else "FAIL"))
    return ok


# ------------------------------------------------------------------ K1
def sweep_K1():
    print("K1 medicion decisiva de la capa alta (max sobre niveles, "
          "/ Q_c):")
    ok = True
    ratios = []        # |E|/(Shat-Pi)^+ en capa alta no degenerada
    ratios_all = []    # idem incluyendo degenerados (E=0)
    fail_mass = 0.0    # masa delta^+ de instancias con delta>0 (capa)
    tot_dpos = 0.0
    F_rows = []
    for (tag, z, D, y, cap, near_u) in FRAMES:
        levels = frame_levels(z, D, y, cap, near_u)
        lnz = math.log(z)
        mx = dict(hi_d=0.0, hi_m=0.0, hi_e=0.0, hi_d_nd=0.0,
                  hi14_d=0.0, hi14_m=0.0, n_hi=0, n_hi_nd=0, n_deg=0)
        for L in levels:
            lv = Level110(z, 3, (), L)
            qc = lv.qc
            s = dict(hi_d=0.0, hi_m=0.0, hi_e=0.0, hi_d_nd=0.0,
                     hi14_d=0.0, hi14_m=0.0)
            for (r, pf, phi, ur, Pi, sh, nhat, dif, E, dl) in \
                    band_instances(lv, L, lnz):
                # familia F: pq, primo mayor con x_p en (11/20, 3/5)
                if len(pf) == 2:
                    pmax = max(pf)
                    xp = math.log(pmax) / lnz
                    if 11.0 / 20.0 < xp < 3.0 / 5.0 and min(pf) >= 5:
                        F_rows.append((z, L, r, ur, Pi, dif, E, dl))
                if ur < U0STAR:
                    continue
                deg = (nhat == 1)
                dp = dl / phi if dl > 0 else 0.0
                mp = dif / phi if dif > 0 else 0.0
                s["hi_d"] += dp
                s["hi_m"] += mp
                s["hi_e"] += abs(E) / phi
                if deg:
                    mx["n_deg"] += 1
                else:
                    s["hi_d_nd"] += dp
                    mx["n_hi_nd"] += 1
                    if dif > 0:
                        ratios.append(abs(E) / dif)
                mx["n_hi"] += 1
                if dif > 0:
                    ratios_all.append(abs(E) / dif)
                if dl > 0:
                    fail_mass += dp
                    tot_dpos += dp
                if ur >= U14:
                    s["hi14_d"] += dp
                    s["hi14_m"] += mp
            for k in ("hi_d", "hi_m", "hi_e", "hi_d_nd",
                      "hi14_d", "hi14_m"):
                mx[k] = max(mx[k], s[k] / qc)
        print("  [%s] z=%d (%d niveles): capa u_r>=1: delta^+ = %.4f | "
              "(Shat-Pi)^+ = %.4f | |E| = %.4f | delta^+ no-deg = %.4f || "
              "u_r>=1.4: delta^+ = %.4f | (Shat-Pi)^+ = %.4f || "
              "instancias capa: %d (no-deg %d, deg %d)"
              % (tag, z, len(levels), mx["hi_d"], mx["hi_m"],
                 mx["hi_e"], mx["hi_d_nd"], mx["hi14_d"],
                 mx["hi14_m"], mx["n_hi"], mx["n_hi_nd"], mx["n_deg"]))
    if ratios:
        rs = sorted(ratios)
        n = len(rs)
        print("  cociente |E|/(Shat-Pi)^+ (capa alta, no degenerado, "
              "dif>0; n=%d): min = %.4f | p25 = %.4f | mediana = %.4f "
              "| p75 = %.4f | max = %.4f"
              % (n, rs[0], rs[n // 4], rs[n // 2], rs[(3 * n) // 4],
                 rs[-1]))
        full = sum(1 for x in rs if x >= 1.0 - 1e-12)
        print("  cancelacion puntual completa (|E| >= (Shat-Pi)^+): "
              "%d/%d = %.1f%% de las instancias"
              % (full, n, 100.0 * full / n))
    print("  familia F en banda con peso real (z, L, r, u_r, Pi, "
          "Shat-Pi, E, delta):")
    if F_rows:
        for row in F_rows[:20]:
            print("    z=%d L=%d r=%d u_r=%.3f Pi=%.4f dif=%+.4f "
                  "E=%+.4f delta=%+.4f" % row)
    else:
        print("    (ningun miembro de F en banda en ventana)")
    print("  => %s" % ("PASS" if ok else "FAIL"))
    return ok


# ------------------------------------------------------------------ K2
def sweep_K2(cap_inst=300):
    print("K2 estructura: E_Add vs E_Rem en la capa alta no degenerada:")
    ok = True
    done = 0
    err_id = 0.0
    n_add_pos = n_rem_pos = 0
    add_share = rem_share = 0.0
    cover_add = cover_rem = cover_both = 0
    n_dif = 0
    r_add = []
    r_rem = []
    frames = [("A", 29, 29 ** 5, int(29 ** (10.0 / 3.0)), None, None),
              ("A-estres", 53, 53 ** 5, int(53 ** (10.0 / 3.0)), 4, 2.0),
              ("A", 37, 37 ** 5, int(37 ** (10.0 / 3.0)), None, None)]
    for (tag, z, D, y, cap, near_u) in frames:
        if done >= cap_inst:
            break
        levels = frame_levels(z, D, y, cap, near_u)
        lnz = math.log(z)
        for L in levels:
            if done >= cap_inst:
                break
            lv = Level110(z, 3, (), L)
            for (r, pf, phi, ur, Pi, sh, nhat, dif, E, dl) in \
                    band_instances(lv, L, lnz):
                if ur < U0STAR or nhat == 1:
                    continue
                ld = lados(lv, r)
                EA = sum(v[2] / v[3] for v in ld.values()
                         if v[0] == "Add")
                ER = sum(v[2] / v[3] for v in ld.values()
                         if v[0] == "Rem")
                # ojo: E = sum_Add a - sum_Rem a (doc 114)
                Ef = EA - ER
                err_id = max(err_id, abs(E - Ef))
                if EA > 1e-12:
                    n_add_pos += 1
                if ER < -1e-12:
                    n_rem_pos += 1
                add_share += abs(EA) / phi
                rem_share += abs(ER) / phi
                done += 1
                if dif > 0:
                    n_dif += 1
                    r_add.append(abs(EA) / dif)
                    r_rem.append(abs(ER) / dif)
                    if abs(EA) >= dif - 1e-15:
                        cover_add += 1
                    if abs(ER) >= dif - 1e-15:
                        cover_rem += 1
                    if abs(Ef) >= dif - 1e-15:
                        cover_both += 1
    print("  instancias descompuestas: %d | |E - (E_Add - E_Rem)| max "
          "= %.2e" % (done, err_id))
    print("  sum_Add a > 0 en %d instancias; sum_Rem a < 0 en %d  "
          "[doc 114: E_Add = sum_Add a <= 0 y E_Rem = -sum_Rem a <= 0 "
          "por separado en ventana]" % (n_add_pos, n_rem_pos))
    tot = add_share + rem_share
    print("  masa |.|/phi en la muestra: |E_Add| = %.4f (%.1f%%) | "
          "|E_Rem| = %.4f (%.1f%%)"
          % (add_share, 100.0 * add_share / max(tot, 1e-300),
             rem_share, 100.0 * rem_share / max(tot, 1e-300)))
    if r_add:
        ra = sorted(r_add)
        rr = sorted(r_rem)
        n = len(ra)
        print("  cocientes medianos contra (Shat-Pi)^+: |E_Add| = %.4f "
              "| |E_Rem| = %.4f (n=%d, dif>0)" % (ra[n // 2], rr[n // 2],
                                                  n))
        print("  cobertura puntual (dif>0): solo Add %d/%d | solo Rem "
              "%d/%d | E total %d/%d"
              % (cover_add, n_dif, cover_rem, n_dif, cover_both, n_dif))
    ok &= err_id < 1e-12
    print("  => %s" % ("PASS" if ok else "FAIL"))
    return ok


# ------------------------------------------------------------------ K3
def sweep_K3():
    print("K3 detalle del cociente: E exactamente 0 vs E != 0 en la "
          "capa alta no degenerada (dif>0):")
    n_tot = n_zero = 0
    ratios_nz = []
    conc = {}
    for (tag, z, D, y, cap, near_u) in FRAMES:
        levels = frame_levels(z, D, y, cap, near_u)
        lnz = math.log(z)
        best = (0.0, 0.0)
        for L in levels:
            lv = Level110(z, 3, (), L)
            s_d = s_d14 = 0.0
            for (r, pf, phi, ur, Pi, sh, nhat, dif, E, dl) in \
                    band_instances(lv, L, lnz):
                if ur < U0STAR:
                    continue
                if dl > 0:
                    s_d += dl / phi
                    if ur >= U14:
                        s_d14 += dl / phi
                if nhat == 1 or dif <= 0:
                    continue
                n_tot += 1
                if E == 0.0:
                    n_zero += 1
                else:
                    ratios_nz.append(abs(E) / dif)
            if s_d > best[0]:
                best = (s_d, s_d14)
        conc[(tag, z)] = best
    print("  instancias: %d | E == 0 exacto: %d (%.1f%%) | E != 0: %d"
          % (n_tot, n_zero, 100.0 * n_zero / max(n_tot, 1),
             n_tot - n_zero))
    if ratios_nz:
        rn = sorted(ratios_nz)
        n = len(rn)
        print("  cociente |E|/(Shat-Pi)^+ solo con E != 0 (n=%d): "
              "min = %.4f | mediana = %.4f | max = %.4f"
              % (n, rn[0], rn[n // 2], rn[-1]))
    for k in sorted(conc):
        s_d, s_d14 = conc[k]
        print("  [%s] z=%d: en el nivel de max delta^+ de capa, "
              "fraccion en u_r>=1.4: %.3f"
              % (k[0], k[1], s_d14 / s_d if s_d > 0 else float("nan")))
    print("  => PASS")
    return True


def main():
    t0 = time.time()
    print("== DOC 120a: (BAND-delta-HI), cancelacion de E en la capa "
          "alta ==")
    ok = True
    ok &= check_K0()
    print()
    ok &= sweep_K1()
    print()
    ok &= sweep_K2()
    print()
    ok &= sweep_K3()
    print()
    print("Tiempo total: %.1f s" % (time.time() - t0))
    print("PASS" if ok else "FAIL")
    return 0 if ok else 1


if __name__ == "__main__":
    sys.exit(main())
