#!/usr/bin/env python3
# 119b_slab_v.py  (Doc 119b)
# (SLAB-V) del doc 118b1 sec. 5: la suma degenerada visible por
# estrato q0 se reindexa por pares (d, c), d en V^+(L), c = d/r <
# q0^3, con 1/phi(r) = phi(c)/phi(d) (libre de cuadrados).
# Chequeos:
#   K0  identidad phi(c)*phi(r) = phi(d) en todo par; cofactor de
#       todo testigo de un r degenerado es < q0^3 [T].
#   K1  cobertura y cadena: suma verdadera <= suma de pares P(q0)
#       (con multiplicidad N(r) >= 1) <= W(q0); constantes contra Q_c.
#   K2  intercambio: W(q0) = sum_c U_c exacto; lista de cofactores
#       que ocurren, mayores contribuyentes, prediccion del modelo
#       Pi_full * A(c) * ln(q0^3).
#   K3  estratos: q0 en {5,7,11} por separado; estratos >= 13 y la
#       parte c=1 de W.
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
visible_set = b1.visible_set
sub_divisors = b1.sub_divisors
FRAMES = b1.FRAMES
U0 = b1.U0

STRATA = (5, 7, 11)
PQ0 = {5: 1, 7: 5, 11: 35}   # P(q0) = prod_{p en P, p < q0} p


def divisors_with_phi(pf):
    """(c, pf_c, phi_c) para cada c | d, d el de pf (libre de cuad.)."""
    return [(c, None, cphi) for (c, cmu, cphi) in sub_divisors(pf)]


def level_data(lv, z):
    """Por nivel: sumas verdadera / de pares / W por estrato,
    multiplicidades, U_c, y violaciones de identidad."""
    info = {r: (pf, phi) for (r, pf, mu, phi) in lv.enum}
    vis = visible_set(lv)
    L, qc = lv.L, lv.qc
    logz = math.log(z)
    # degenerados visibles en banda, por estrato (como doc 118b1)
    deg, _ = band_degen(lv, z)
    true_sum = {}
    deg_set = {}
    for (r, pf, phi, q0) in deg:
        key = q0 if q0 in STRATA else ("otro" if q0 else None)
        true_sum[key] = true_sum.get(key, 0.0) + 1.0 / phi
        deg_set[r] = q0
    # pares (d, c): d en V, c | d
    pair_sum = {}          # por estrato de r = d/c
    pair_c = {}            # pair_c[q0][c]: suma de pares por cofactor
    mult = {}              # N(r) para r degenerado visible en banda
    U = {q0: {} for q0 in STRATA}   # U[(q0)][c]
    W = {q0: 0.0 for q0 in STRATA}
    W_c1 = {q0: 0.0 for q0 in STRATA}
    bad_phi = 0
    bad_cof = 0
    for (d, dpf, dmu, dphi) in lv.V:
        subs = sub_divisors(dpf)
        subphi = {c: cphi for (c, cmu, cphi) in subs}
        subpf = {}
        for (c, cmu, cphi) in subs:
            # primos de c: los de dpf cuyo producto parcial... se
            # reconstruyen por division sucesiva (libre de cuadrados)
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
            # contribucion a W(q0) si el estrato de r es q0 y d en
            # rebanada; el peso mayorante solo pide (c, P(q0)) = 1
            q0r = q0_of(lv, rpf)
            if r in deg_set and deg_set[r] == q0r:
                if c >= (q0r ** 3 if q0r else c + 1):
                    bad_cof += 1
                mult[r] = mult.get(r, 0) + 1
                if q0r in STRATA:
                    pair_sum[q0r] = pair_sum.get(q0r, 0.0) + 1.0 / rphi
                    pc = pair_c.setdefault(q0r, {})
                    pc[c] = pc.get(c, 0.0) + 1.0 / rphi
            if q0r in STRATA:
                B = q0r ** 3
                if d >= L / B and c < B and math.gcd(c, PQ0[q0r]) == 1:
                    w = cphi / dphi
                    W[q0r] += w
                    U[q0r][c] = U[q0r].get(c, 0.0) + 1.0 / rphi
                    if c == 1:
                        W_c1[q0r] += w
    return {
        "true": true_sum, "pair": pair_sum, "pair_c": pair_c,
        "W": W, "W_c1": W_c1,
        "U": U, "mult": mult, "deg_set": deg_set,
        "bad_phi": bad_phi, "bad_cof": bad_cof, "qc": qc,
    }


def A_of_c(lv, c):
    """Constante de Mertens del modelo para el cofactor c."""
    a = 1.0
    for p in lv.allow:
        if c % p == 0:
            a *= (1.0 - 1.0 / p)
        else:
            a *= (1.0 + 1.0 / (p * (p - 1)))
    return a


def run():
    print("K0/K1/K2/K3 reindexacion, cobertura, intercambio y "
          "constantes por estrato:")
    ok = True
    for (tag, z, D, y) in FRAMES:
        levels = frame_levels(z, D, y)
        agg = {q0: {"true": 0.0, "pair": 0.0, "W": 0.0, "W_c1": 0.0}
               for q0 in STRATA}
        otro = 0.0
        bad_phi = bad_cof = 0
        n_uncovered = 0
        mult_max = 0
        mult_vals = []
        swap_err = 0.0
        U_tot = {q0: {} for q0 in STRATA}
        U_pred = {q0: {} for q0 in STRATA}
        PC_tot = {q0: {} for q0 in STRATA}
        for lv in levels:
            dat = level_data(lv, z)
            qc = dat["qc"]
            bad_phi += dat["bad_phi"]
            bad_cof += dat["bad_cof"]
            for r, q0 in dat["deg_set"].items():
                n = dat["mult"].get(r, 0)
                if q0 in STRATA:
                    if n == 0:
                        n_uncovered += 1
                    mult_max = max(mult_max, n)
                    mult_vals.append(n)
            for q0 in STRATA:
                agg[q0]["true"] = max(
                    agg[q0]["true"],
                    dat["true"].get(q0, 0.0) / qc)
                agg[q0]["pair"] = max(
                    agg[q0]["pair"],
                    dat["pair"].get(q0, 0.0) / qc)
                agg[q0]["W"] = max(agg[q0]["W"], dat["W"][q0] / qc)
                agg[q0]["W_c1"] = max(
                    agg[q0]["W_c1"], dat["W_c1"][q0] / qc)
                sU = sum(dat["U"][q0].values())
                swap_err = max(swap_err, abs(sU - dat["W"][q0]))
                for c, v in dat["U"][q0].items():
                    if v / qc > U_tot[q0].get(c, (0.0,))[0]:
                        U_tot[q0][c] = (v / qc,)
                        U_pred[q0][c] = (
                            lv.Pi_full * A_of_c(lv, c)
                            * math.log(q0 ** 3) / qc,)
                for c, v in dat["pair_c"].get(q0, {}).items():
                    if v / qc > PC_tot[q0].get(c, 0.0):
                        PC_tot[q0][c] = v / qc
            otro = max(otro, sum(
                v for k, v in dat["true"].items()
                if k not in STRATA) / qc)
        ok &= bad_phi == 0 and bad_cof == 0 and n_uncovered == 0
        ok &= swap_err < 1e-9
        mean_mult = (sum(mult_vals) / len(mult_vals)
                     if mult_vals else 0.0)
        print("  [%s] z=%d: viol phi=%d cof=%d sin-cobertura=%d | "
              "N(r): max=%d media=%.2f | err intercambio=%.1e"
              % (tag, z, bad_phi, bad_cof, n_uncovered, mult_max,
                 mean_mult, swap_err))
        for q0 in STRATA:
            a = agg[q0]
            tops = sorted(U_tot[q0].items(),
                          key=lambda kv: -kv[1][0])[:4]
            tops_s = " ".join(
                "c=%d:%.3f(pred %.3f)" % (c, v[0], U_pred[q0][c][0])
                for c, v in tops)
            print("    q0=%2d: verdadera=%.4f pares=%.4f W=%.4f "
                  "(c=1: %.4f) Q_c | #c=%d | top U_c/Q_c: %s"
                  % (q0, a["true"], a["pair"], a["W"], a["W_c1"],
                     len(U_tot[q0]), tops_s))
            if PC_tot[q0]:
                tops_pc = sorted(PC_tot[q0].items(),
                                 key=lambda kv: -kv[1])[:5]
                print("           pares por cofactor (max/nivel): %s"
                      % " ".join("c=%d:%.4f" % (c, v)
                                 for c, v in tops_pc))
        print("    estratos >=13/None: verdadera max = %.4f Q_c"
              % otro)
    print("  => %s" % ("PASS" if ok else "FAIL"))
    return ok


def main():
    t0 = time.time()
    print("== DOC 119b: (SLAB-V) por reindexacion divisor-cofactor ==")
    ok = run()
    print()
    print("Tiempo total: %.1f s" % (time.time() - t0))
    print("PASS" if ok else "FAIL")
    return 0 if ok else 1


if __name__ == "__main__":
    sys.exit(main())
