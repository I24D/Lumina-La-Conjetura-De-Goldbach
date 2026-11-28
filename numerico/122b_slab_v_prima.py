#!/usr/bin/env python3
# 122b_slab_v_prima.py  (Doc 122b)
# (SLAB-V') del doc 119b sec. 4: P(q0) = sum_c U_c^deg(q0) <=
# C(q0) Q_c para q0 en {5,7,11}. Medicion extendida de la familia
# completa de sumas U_c^deg en marcos A (z=17,29,37,53) y B (z=29):
#   K0  identidad phi(c)*phi(r) = phi(d) en todo par; cobertura
#       N(r) >= 1 en todo r degenerado visible de estrato finito.
#   K1  doble contabilidad: P por cofactores == P por r con
#       multiplicidad N(r); verdadera <= P.
#   K2  tablas por (marco, q0, c): max_nivel U_c^deg/Q_c, P/Q_c,
#       verdadera/Q_c, concentracion (parte de c=1 y top-3),
#       nivel del maximo de P, tendencia en z y en L.
#   K3  Teorema A: si L >= q0^9, ningun degenerado visible del
#       estrato q0 con omega <= 2; conteo de omega<=2 fuera de
#       hipotesis (artefactos de ventana).
# Solo stdlib; importa la maquinaria de los docs 110/117b/118b1.
import importlib.util
import math
import sys
import time

spec = importlib.util.spec_from_file_location(
    "b1", "118b1_degenerada.py")
b1 = importlib.util.module_from_spec(spec)
spec.loader.exec_module(b1)

Level110 = b1.Level110
frame_levels = b1.frame_levels
q0_of = b1.q0_of
band_degen = b1.band_degen
sub_divisors = b1.sub_divisors

FRAMES = list(b1.FRAMES) + [("A", 53, 53 ** 5, int(53 ** (10.0 / 3.0)))]
STRATA = (5, 7, 11)


def level_data(lv, z):
    """Por nivel: U^deg por (q0, c), P, verdadera, N(r), omega min,
    rebanada sum_{d en V+, d >= L/q0^3} 1/phi(d), violaciones."""
    info = {r: (pf, phi) for (r, pf, mu, phi) in lv.enum}
    L, qc = lv.L, lv.qc
    deg, _ = band_degen(lv, z)
    true_sum = {q0: 0.0 for q0 in STRATA}
    deg_set = {}
    omega = {}
    for (r, pf, phi, q0) in deg:
        deg_set[r] = q0
        omega[r] = len(pf)
        if q0 in STRATA:
            true_sum[q0] += 1.0 / phi
    Udeg = {q0: {} for q0 in STRATA}
    mult = {}
    bad_phi = 0
    for (d, dpf, dmu, dphi) in lv.V:
        subs = sub_divisors(dpf)
        subpf = {}
        for (c, cmu, cphi) in subs:
            cc, pp = c, []
            for p in dpf:
                if cc % p == 0:
                    pp.append(p)
                    cc //= p
            subpf[c] = tuple(pp)
        for (c, cmu, cphi) in subs:
            r = d // c
            cpf = set(subpf[c])
            rpf = tuple(p for p in dpf if p not in cpf)
            rphi = info[r][1]
            if cphi * rphi != dphi:
                bad_phi += 1
            if r in deg_set:
                q0r = q0_of(lv, rpf)
                if deg_set[r] == q0r and q0r in STRATA:
                    mult[r] = mult.get(r, 0) + 1
                    Udeg[q0r][c] = Udeg[q0r].get(c, 0.0) + 1.0 / rphi
    P_by_c = {q0: sum(Udeg[q0].values()) for q0 in STRATA}
    P_by_r = {q0: 0.0 for q0 in STRATA}
    for r, n in mult.items():
        P_by_r[deg_set[r]] += n / info[r][1]
    uncovered = sum(1 for r, q0 in deg_set.items()
                    if q0 in STRATA and mult.get(r, 0) == 0)
    slice_sum = {}
    for q0 in STRATA:
        thr = L / (q0 ** 3)
        slice_sum[q0] = sum(1.0 / dphi for (d, dpf, dmu, dphi)
                            in lv.V if d >= thr)
    min_om = {}
    n_om_le2 = {}
    for q0 in STRATA:
        oms = [omega[r] for r, qq in deg_set.items() if qq == q0]
        min_om[q0] = min(oms) if oms else None
        n_om_le2[q0] = sum(1 for o in oms if o <= 2)
    return {"L": L, "qc": qc, "true": true_sum, "Udeg": Udeg,
            "P_by_c": P_by_c, "P_by_r": P_by_r, "mult": mult,
            "bad_phi": bad_phi, "uncovered": uncovered,
            "slice": slice_sum, "min_om": min_om,
            "n_om_le2": n_om_le2,
            "ndeg": {q0: sum(1 for qq in deg_set.values() if qq == q0)
                     for q0 in STRATA},
            "deg_set": deg_set}


def run():
    ok = True
    for (tag, z, D, y) in FRAMES:
        levels = frame_levels(z, D, y)
        agg = {q0: {"true": 0.0, "P": 0.0, "L_at_P": None,
                    "Uc": {}, "Nmax": 0, "slice": 0.0}
               for q0 in STRATA}
        bad_phi = uncovered = 0
        swap_err = 0.0
        thmA_viol = 0
        om_le2_in = om_le2_out = 0
        rows = []
        for lv in levels:
            dat = level_data(lv, z)
            qc = dat["qc"]
            bad_phi += dat["bad_phi"]
            uncovered += dat["uncovered"]
            for q0 in STRATA:
                swap_err = max(swap_err, abs(dat["P_by_c"][q0]
                                             - dat["P_by_r"][q0]))
                if dat["true"][q0] > dat["P_by_c"][q0] + 1e-12:
                    ok = False
                p = dat["P_by_c"][q0] / qc
                if p > agg[q0]["P"]:
                    agg[q0]["P"] = p
                    agg[q0]["L_at_P"] = dat["L"]
                agg[q0]["true"] = max(agg[q0]["true"],
                                      dat["true"][q0] / qc)
                agg[q0]["slice"] = max(agg[q0]["slice"],
                                       dat["slice"][q0] / qc)
                for c, v in dat["Udeg"][q0].items():
                    if v / qc > agg[q0]["Uc"].get(c, 0.0):
                        agg[q0]["Uc"][c] = v / qc
                for r, n in dat["mult"].items():
                    if dat["deg_set"].get(r) == q0:
                        agg[q0]["Nmax"] = max(agg[q0]["Nmax"], n)
                if dat["L"] >= q0 ** 9:
                    if dat["n_om_le2"][q0] > 0:
                        thmA_viol += dat["n_om_le2"][q0]
                    om_le2_in += 0
                else:
                    om_le2_out += dat["n_om_le2"][q0]
            rows.append((dat["L"], dat))
        ok &= bad_phi == 0 and uncovered == 0 and swap_err < 1e-9
        ok &= thmA_viol == 0
        print("[%s] z=%d (%d niveles): viol phi=%d sin-cobertura=%d "
              "err doble-contab=%.1e | TeorA viol=%d, omega<=2 "
              "fuera de hipotesis=%d"
              % (tag, z, len(levels), bad_phi, uncovered, swap_err,
                 thmA_viol, om_le2_out))
        for q0 in STRATA:
            a = agg[q0]
            ucs = sorted(a["Uc"].items(), key=lambda kv: -kv[1])
            tot_uc = sum(a["Uc"].values())
            top3 = sum(v for _, v in ucs[:3])
            print("  q0=%2d: verdadera=%.4f P=%.4f (en L=%d) "
                  "Nmax=%d rebanada=%.4f Q_c | #c=%d c=1: %.4f "
                  "top3/maxU=%.2f"
                  % (q0, a["true"], a["P"], a["L_at_P"] or -1,
                     a["Nmax"], a["slice"], len(ucs),
                     a["Uc"].get(1, 0.0),
                     top3 / tot_uc if tot_uc > 0 else 0.0))
            if ucs:
                print("        U_c^deg/Q_c (max/nivel): %s"
                      % " ".join("c=%d:%.4f" % (c, v)
                                 for c, v in ucs[:8]))
        # tendencia en L dentro del marco (estrato 5)
        trend = [(dat_L, d) for (dat_L, d) in rows]
        if trend:
            lo = min(trend, key=lambda t: t[0])[1]
            hi = max(trend, key=lambda t: t[0])[1]
            print("  tendencia L: P(5)/Q_c en Lmin=%d: %.4f | en "
                  "Lmax=%d: %.4f"
                  % (lo["L"], lo["P_by_c"][5] / lo["qc"],
                     hi["L"], hi["P_by_c"][5] / hi["qc"]))
    print("  => %s" % ("PASS" if ok else "FAIL"))
    return ok


def main():
    t0 = time.time()
    print("== DOC 122b: (SLAB-V') familia U_c^deg, z=17..53 ==")
    ok = run()
    print()
    print("Tiempo total: %.1f s" % (time.time() - t0))
    print("PASS" if ok else "FAIL")
    return 0 if ok else 1


if __name__ == "__main__":
    sys.exit(main())
