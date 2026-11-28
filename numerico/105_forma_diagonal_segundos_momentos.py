#!/usr/bin/env python3
# 105_forma_diagonal_segundos_momentos.py  (Doc 105)
# Redireccion del doc 104: B_M en su forma diagonal con SOLO
# segundos momentos. Solo stdlib.
#
# Forma diagonal [T]: con F_d = phi(d) E_d (d >= 2; F_1 = 0),
#   B_M^(M) = sum_{d,e en D} F_d F_e W_{de},
#   W_{de} = sum_{r en Hcal, d|r, e|r} mu(r/d) mu(r/e) / phi(r),
#   Hcal = {r en H_N : r > R0}, R0 = (ln N)^kappa, kappa = 1,
#   D = {d >= 2 : d | r para algun r en Hcal}.
# Factorizacion [T]: todo r en H es libre de cuadrados; para
# L = mcm(d, e), r/d = (r/L)(L/d) con factores coprimos, luego
# mu(r/d) mu(r/e) = mu(L/d) mu(L/e) (mu(r/L)^2 = 1) y
#   W_{de} = mu(L/d) mu(L/e) * U_L,
#   U_L := sum_{r en Hcal, L|r} 1/phi(r) >= 0,
# en particular |W_{de}| = U_{mcm(d,e)} (0 si ningun r lo cubre).
#
# Candidata Bhat5 (mayorante [T], solo segundos momentos):
#   Omega_d^X := sum_{e en D} |W^X_{de}|
#              = sum_{r en Hcal_X, d|r} (tau(r) - 1) / phi(r),
#   (la igualdad es conteo: sum_e U_{mcm(d,e)} suma sobre r el
#    numero de e >= 2 con e | r, que es tau(r) - 1)
#   Dhat^X := sum_d phi(d)^2 E_d^2 Omega_d^X   (X = P, C),
#   Bhat5^X := ( sqrt(Dhat^X) + sqrt(B_Delta^X) )^2,
#   B_Delta^X := sum_{r en Hcal_X} Delta_r^2 / phi(r),
#   Bhat5 := Bhat5^P + Bhat5^C >= B_M^(S).
# Validez [T]: (i) por pieza, AM-GM sobre pares (d, e) con la
# simetria |W_{de}| = |W_{ed}| da B_M^{X,(M)} <= Dhat^X;
# (ii) Minkowski en l^2 con pesos 1/phi(r) da
# B_M^{X,(S)} <= (sqrt(B_M^{X,(M)}) + sqrt(B_Delta^X))^2.
# W^P es diagonal en D (r primo: unicos divisores >= 2: el mismo
# r), luego Dhat^P = Diag^P = B_M^{P,(M)} exactamente [T].
#
# Criterio previo (doc 104): (i) expresabilidad solo en segundos
# momentos (verificado por construccion: Dhat^X usa E_d^2 y pesos
# del soporte; B_Delta^X es el segundo momento de la correccion
# explicita Delta_r, docs 97-100); (ii) C1+C2 en los seis puntos,
# aplicado en agregado y por piezas (ver doc 105, seccion 2, para
# el argumento de cual lectura de C2 es la honesta).
#
# Chequeos:
#   K0  Maquinaria: Model97 R2 |H|=880; build_H z=11 == Model97.
#   K1  B_M^(S) reproduce los docs 103/104 en los seis puntos
#       (rel < 2e-3) y Diag^C/B_M^C reproduce 1.488-1.674
#       (rel < 3e-2).
#   K2  Forma diagonal en R1 y R3 (directo): la forma cuadratica
#       sum_{d,e} F_d F_e W_{de} = B_M^(M); |W_{de}| =
#       U_{mcm(d,e)} con signo mu(L/d)mu(L/e); sumas por fila
#       sum_e |W_{de}| = Omega_d (formula tau); Diag por la suma
#       en d = Diag por la suma en r.
#   K3  Mayorante valido: Dhat^X >= B_M^{X,(M)} y Bhat5^X >=
#       B_M^{X,(S)} en todo punto; Bhat5 >= B_M^(S).
#   K4  Pieza prima exacta: Dhat^P = B_M^{P,(M)} (rel < 1e-9)
#       en todo punto.
#   K5  Cociente agregado = promedio ponderado de los cocientes
#       por pieza (pesos B_M^{X,(S)}/B_M^(S)); en particular
#       cociente agregado <= max(cocientes por pieza).
#   K6  C1+C2 para Bhat5: agregado y por piezas (reporte).
# Uso: python3 105_forma_diagonal_segundos_momentos.py
import importlib.util
import math
import sys
import time

spec = importlib.util.spec_from_file_location(
    "s104", "104_experimento_e103_1_bhat3.py")
s104 = importlib.util.module_from_spec(spec)
spec.loader.exec_module(s104)

Model97 = s104.Model97
build_H = s104.build_H
SegPsi = s104.SegPsi
M_r_E = s104.M_r_E
ols = s104.ols
divisors_of = s104.divisors_of
phi_of = s104.phi_of
is_prime_spf = s104.is_prime_spf
delta_r_orth = s104.delta_r_orth
MU, SPF, SMALL_PRIMES = s104.MU, s104.SPF, s104.SMALL_PRIMES


def measure105(H, N, kappa=1.0):
    R0 = math.log(N) ** kappa
    HN = [r for r in H if math.gcd(r, N) == 1]
    divs_set = set()
    for r in HN:
        for d in divisors_of(r, SPF):
            divs_set.add(d)
    divs = sorted(divs_set)
    seg = SegPsi(N, divs)
    E = seg.E
    phi = {v: phi_of(v, SMALL_PRIMES) for v in set(divs) | set(HN)}
    Hcal = [r for r in HN if r > R0]
    # pesos Omega^X y U^X por conteo (formula tau)
    Om = {"P": {}, "C": {}}
    Uu = {"P": {}, "C": {}}
    for r in Hcal:
        X = "P" if is_prime_spf(r) else "C"
        tau = len(divisors_of(r, SPF))
        w_om = (tau - 1) / phi[r]
        w_u = 1.0 / phi[r]
        for d in divisors_of(r, SPF):
            if d < 2:
                continue
            Om[X][d] = Om[X].get(d, 0.0) + w_om
            Uu[X][d] = Uu[X].get(d, 0.0) + w_u
    # momentos por pieza
    BMm = {"P": 0.0, "C": 0.0}
    BMs = {"P": 0.0, "C": 0.0}
    BD = {"P": 0.0, "C": 0.0}
    Bcross = {"P": 0.0, "C": 0.0}
    for r in Hcal:
        X = "P" if is_prime_spf(r) else "C"
        Mr = M_r_E(r, E, MU, SPF)
        Dr = delta_r_orth(r, N, MU, SPF)
        ph = phi[r]
        BMm[X] += Mr * Mr / ph
        BMs[X] += (Mr + Dr) * (Mr + Dr) / ph
        BD[X] += Dr * Dr / ph
        Bcross[X] += 2.0 * Mr * Dr / ph
    out = {"N": N, "nHN": len(HN), "nHcal": len(Hcal)}
    for X in ("P", "C"):
        Dhat = 0.0
        Diag = 0.0
        for d, om in Om[X].items():
            f2 = (phi[d] * E.get(d, 0.0)) ** 2
            Dhat += f2 * om
            Diag += f2 * Uu[X][d]
        out["Dhat_" + X] = Dhat
        out["Diag_" + X] = Diag
        out["BMm_" + X] = BMm[X]
        out["BMs_" + X] = BMs[X]
        out["BD_" + X] = BD[X]
        out["Bcross_" + X] = Bcross[X]
        out["Bhat_" + X] = (math.sqrt(Dhat) + math.sqrt(BD[X])) ** 2
    out["BMm"] = BMm["P"] + BMm["C"]
    out["BMs"] = BMs["P"] + BMs["C"]
    out["Bhat"] = out["Bhat_P"] + out["Bhat_C"]
    out["Dhat"] = out["Dhat_P"] + out["Dhat_C"]
    out["Diag"] = out["Diag_P"] + out["Diag_C"]
    return out


def w_direct_check(H, N, kappa=1.0):
    """K2: forma cuadratica directa en soportes pequenos."""
    R0 = math.log(N) ** kappa
    HN = [r for r in H if math.gcd(r, N) == 1]
    divs_set = set()
    for r in HN:
        for d in divisors_of(r, SPF):
            divs_set.add(d)
    seg = SegPsi(N, sorted(divs_set))
    E = seg.E
    phi = {v: phi_of(v, SMALL_PRIMES)
           for v in set(divs_set) | set(HN)}
    Hcal = [r for r in HN if r > R0]
    Dset = sorted({d for r in Hcal for d in divisors_of(r, SPF)
                   if d >= 2})
    idx = {d: i for i, d in enumerate(Dset)}
    n = len(Dset)
    W = [[0.0] * n for _ in range(n)]
    U = {}
    for r in Hcal:
        for L in divisors_of(r, SPF):
            U[L] = U.get(L, 0.0) + 1.0 / phi[r]
        dv = [d for d in divisors_of(r, SPF) if d >= 2]
        for d in dv:
            for e in dv:
                W[idx[d]][idx[e]] += \
                    MU[r // d] * MU[r // e] / phi[r]
    F = {d: phi[d] * E.get(d, 0.0) for d in Dset}
    quad = 0.0
    for d in Dset:
        for e in Dset:
            quad += F[d] * F[e] * W[idx[d]][idx[e]]
    BMm = 0.0
    Diag_r = 0.0
    for r in Hcal:
        Mr = M_r_E(r, E, MU, SPF)
        BMm += Mr * Mr / phi[r]
        s2 = 0.0
        for d in divisors_of(r, SPF):
            if d >= 2:
                s2 += F[d] * F[d]
        Diag_r += s2 / phi[r]
    ok_quad = abs(quad - BMm) / BMm < 1e-9
    ok_fact = True
    worst = 0.0
    for d in Dset:
        for e in Dset:
            L = d * e // math.gcd(d, e)
            u = U.get(L, 0.0)
            pred = MU[L // d] * MU[L // e] * u
            err = abs(W[idx[d]][idx[e]] - pred)
            scale = max(1.0, abs(pred))
            worst = max(worst, err / scale)
            if err / scale > 1e-9:
                ok_fact = False
    # sumas por fila contra la formula tau (Omega total)
    Om = {}
    for r in Hcal:
        tau = len(divisors_of(r, SPF))
        for d in divisors_of(r, SPF):
            if d >= 2:
                Om[d] = Om.get(d, 0.0) + (tau - 1) / phi[r]
    ok_row = True
    for d in Dset:
        rs = sum(abs(W[idx[d]][idx[e]]) for e in Dset)
        if abs(rs - Om.get(d, 0.0)) / max(1.0, Om.get(d, 0.0)) \
                > 1e-9:
            ok_row = False
    # Diag por suma en d (pesos U totales) contra Diag por suma en r
    Utot = {}
    for r in Hcal:
        for d in divisors_of(r, SPF):
            if d >= 2:
                Utot[d] = Utot.get(d, 0.0) + 1.0 / phi[r]
    Diag_d = sum(F[d] * F[d] * Utot[d] for d in Dset)
    ok_diag = abs(Diag_d - Diag_r) / Diag_r < 1e-9
    return {"ok_quad": ok_quad, "ok_fact": ok_fact,
            "ok_row": ok_row, "ok_diag": ok_diag,
            "worst_fact": worst, "nD": n, "quad": quad, "BMm": BMm}


def fmt_pt(tag, m):
    return ("     [%s] N=%d |H_N|=%d |Hcal|=%d\n"
            "        B_M^(S)=%.6e (P %.6e / C %.6e)  "
            "B_M^(M)=%.6e\n"
            "        Bhat5/B_M^(S)=%.3f  (P x%.3f  C x%.3f)\n"
            "        Dhat/B_M^(M): tot x%.3f  P x%.4f  C x%.3f;"
            "  Diag/B_M^(M): P x%.4f  C x%.3f\n"
            "        precio tau (Dhat/Diag): P x%.4f  C x%.3f;"
            "  share_C(B_M^(S))=%.3f  B_Delta/B_M^(S)=%.2e"
            % (tag, m["N"], m["nHN"], m["nHcal"], m["BMs"],
               m["BMs_P"], m["BMs_C"], m["BMm"],
               m["Bhat"] / m["BMs"], m["Bhat_P"] / m["BMs_P"],
               m["Bhat_C"] / m["BMs_C"] if m["BMs_C"] > 0
               else float("nan"),
               m["Dhat"] / m["BMm"], m["Dhat_P"] / m["BMm_P"],
               m["Dhat_C"] / m["BMm_C"] if m["BMm_C"] > 0
               else float("nan"),
               m["Diag_P"] / m["BMm_P"],
               m["Diag_C"] / m["BMm_C"] if m["BMm_C"] > 0
               else float("nan"),
               m["Dhat_P"] / m["Diag_P"],
               m["Dhat_C"] / m["Diag_C"] if m["Diag_C"] > 0
               else float("nan"),
               m["BMs_C"] / m["BMs"],
               (m["BD_P"] + m["BD_C"]) / m["BMs"]))


def main():
    allok = True
    t00 = time.time()
    print("== DOC 105: forma diagonal con solo segundos momentos ==")
    print("Solo stdlib. MEDICION, no prueba. Goldbach sigue "
          "abierta.")
    print("Candidata Bhat5: filas de |W| (AM-GM) por pieza + "
          "Minkowski para Delta_r.")
    print()
    m1 = Model97(3, 11, 60, 20000, "fixed")
    m2 = Model97(3, 11, 3000, 200000, "fixed")
    m3 = Model97(3, 3, 349, 200000, "fixed")
    supports = {"R1": (m1.H, 11), "R2": (m2.H, 11), "R3": (m3.H, 3)}
    for z in (11, 13, 17):
        D = z ** 5
        y = int(round(z ** (10.0 / 3.0)))
        Hz, _, _ = build_H(3, z, y, D, MU, SPF)
        supports["S%d" % z] = (Hz, z)
    print("Soportes: R1=%d R3=%d R2=%d S11=%d S13=%d S17=%d"
          % tuple(len(supports[t][0])
                  for t in ("R1", "R3", "R2", "S11", "S13", "S17")))
    D11 = 11 ** 5
    y11 = int(round(11 ** (10.0 / 3.0)))
    m_s11 = Model97(3, 11, y11, D11, "fixed")
    k0 = len(m2.H) == 880 and set(supports["S11"][0]) == set(m_s11.H)
    allok &= k0
    print("K0 maquinaria (R2 |H|=880; build_H z=11 == Model97): %s"
          % ("PASS" if k0 else "FAIL"))
    print()
    print("Calibracion (N=999982): R1 y R3:")
    allpts = {}
    for tag in ("R1", "R3"):
        H, z = supports[tag]
        t0 = time.time()
        m = measure105(H, 999982)
        allpts["gen/" + tag] = m
        print(fmt_pt("gen/" + tag, m) + " (%.0fs)" % (time.time() - t0))
    print()
    print("Seis puntos del criterio C1+C2 (doc 102):")
    plan = [(999982, ["S11", "S13", "S17"], "gen"),
            (1048576, ["S11", "S13", "S17"], "2^e")]
    six = []
    for (N, tags, fam) in plan:
        for tag in tags:
            H, z = supports[tag]
            t0 = time.time()
            m = measure105(H, N)
            six.append(("%s/%s" % (fam, tag), m))
            allpts["%s/%s" % (fam, tag)] = m
            print(fmt_pt("%s/%s" % (fam, tag), m)
                  + " (%.0fs)" % (time.time() - t0))
    print()
    # K1: reproduccion de los docs 103/104
    doc103_BM = {"gen/S11": 3.3658e9, "gen/S13": 7.6433e9,
                 "gen/S17": 2.7875e10, "2^e/S11": 3.7105e9,
                 "2^e/S13": 8.3332e9, "2^e/S17": 2.9761e10}
    doc103_DiagC = {"gen/S11": 1.495, "gen/S13": 1.488,
                    "gen/S17": 1.533, "2^e/S11": 1.674,
                    "2^e/S13": 1.616, "2^e/S17": 1.576}
    k1 = True
    for tag, ref in doc103_BM.items():
        m = allpts[tag]
        if abs(m["BMs"] - ref) / ref > 2e-3:
            k1 = False
            print("     K1 DISCREPA B_M^(S) [%s]: %.6e vs %.6e"
                  % (tag, m["BMs"], ref))
    for tag, ref in doc103_DiagC.items():
        m = allpts[tag]
        got = m["Diag_C"] / m["BMs_C"]
        if abs(got - ref) / ref > 3e-2:
            k1 = False
            print("     K1 DISCREPA Diag^C/B_M^C(S) [%s]: %.3f vs "
                  "%.3f" % (tag, got, ref))
    allok &= k1
    print("K1 reproduccion docs 103/104 (B_M^(S), Diag^C/B_M^C): %s"
          % ("PASS" if k1 else "FAIL (ver arriba)"))
    # K2: forma diagonal directa en R1 y R3
    k2 = True
    for tag in ("R1", "R3"):
        H, z = supports[tag]
        w = w_direct_check(H, 999982)
        okk = w["ok_quad"] and w["ok_fact"] and w["ok_row"] \
            and w["ok_diag"]
        k2 &= okk
        print("K2 [%s] forma cuadratica=%s factorizacion |W|=%s "
              "(peor rel %.2e) filas tau=%s Diag d=r %s (|D|=%d): %s"
              % (tag, "OK" if w["ok_quad"] else "FALLO",
                 "OK" if w["ok_fact"] else "FALLO", w["worst_fact"],
                 "OK" if w["ok_row"] else "FALLO",
                 "OK" if w["ok_diag"] else "FALLO", w["nD"],
                 "PASS" if okk else "FAIL"))
    allok &= k2
    # K3: mayorante valido en todo punto
    k3 = all(m["Dhat_P"] >= m["BMm_P"] * (1 - 1e-12)
             and m["Dhat_C"] >= m["BMm_C"] * (1 - 1e-12)
             and m["Bhat_P"] >= m["BMs_P"] * (1 - 1e-12)
             and m["Bhat_C"] >= m["BMs_C"] * (1 - 1e-12)
             and m["Bhat"] >= m["BMs"] * (1 - 1e-12)
             for m in allpts.values())
    allok &= k3
    print("K3 mayorante: Dhat^X >= B_M^{X,(M)}, Bhat5^X >= "
          "B_M^{X,(S)}, Bhat5 >= B_M^(S) en todo punto: %s"
          % ("PASS" if k3 else "FAIL"))
    # K4: pieza prima exacta (forma M)
    k4 = all(abs(m["Dhat_P"] - m["BMm_P"]) / m["BMm_P"] < 1e-9
             for m in allpts.values())
    allok &= k4
    print("K4 pieza prima exacta: Dhat^P = B_M^{P,(M)} "
          "(rel < 1e-9) en todo punto: %s"
          % ("PASS" if k4 else "FAIL"))
    # K5: cociente agregado = promedio ponderado de piezas
    k5 = True
    for (tag, m) in six:
        q = m["Bhat"] / m["BMs"]
        qw = (m["Bhat_P"] + m["Bhat_C"]) / m["BMs"]
        qp = m["Bhat_P"] / m["BMs_P"]
        qc = m["Bhat_C"] / m["BMs_C"]
        if abs(q - qw) > 1e-12 or q > max(qp, qc) * (1 + 1e-12):
            k5 = False
            print("     K5 FALLO [%s]" % tag)
    allok &= k5
    print("K5 cociente agregado = promedio ponderado de piezas "
          "y <= max(piezas) en los seis puntos: %s"
          % ("PASS" if k5 else "FAIL"))
    print()
    # Pendientes: agregado (S), piezas (M y S), Diag de referencia,
    # contrafactual de composicion fija
    xs = [math.log(m["nHN"]) for (_, m) in six]
    def slope(ys):
        return ols(xs, ys)
    y_agg = [math.log(m["Bhat"] / m["BMs"]) for (_, m) in six]
    y_PM = [math.log(m["Dhat_P"] / m["BMm_P"]) for (_, m) in six]
    y_CM = [math.log(m["Dhat_C"] / m["BMm_C"]) for (_, m) in six]
    y_PS = [math.log(m["Bhat_P"] / m["BMs_P"]) for (_, m) in six]
    y_CS = [math.log(m["Bhat_C"] / m["BMs_C"]) for (_, m) in six]
    y_DC = [math.log(m["Diag_C"] / m["BMm_C"]) for (_, m) in six]
    b_agg = slope(y_agg)
    b_PM = slope(y_PM)
    b_CM = slope(y_CM)
    b_PS = slope(y_PS)
    b_CS = slope(y_CS)
    b_DC = slope(y_DC)
    # contrafactual: composicion (shares de B_M^(S)) congelada en
    # la media de los seis puntos
    wC = sum(m["BMs_C"] / m["BMs"] for (_, m) in six) / 6.0
    y_cf = []
    for (_, m) in six:
        qp = m["Bhat_P"] / m["BMs_P"]
        qc = m["Bhat_C"] / m["BMs_C"]
        y_cf.append(math.log((1 - wC) * qp + wC * qc))
    b_cf = slope(y_cf)
    def pr(name, t):
        print("     %s pendiente=%.4f +- %.4f (SE), IC95%% "
              "+-%.4f, LOO [%.4f, %.4f]"
              % (name, t[0], t[2], t[3], t[4], t[5]))
    print("Pendientes log-log contra |H_N| (seis puntos):")
    pr("agregado Bhat5/B_M^(S)      ", b_agg)
    pr("pieza P, forma M (Dhat)     ", b_PM)
    pr("pieza C, forma M (Dhat)     ", b_CM)
    pr("pieza P, forma S (Bhat5^P)  ", b_PS)
    pr("pieza C, forma S (Bhat5^C)  ", b_CS)
    pr("Diag^C/B_M^{C,(M)} (ref.)   ", b_DC)
    pr("contrafactual compos. fija  ", b_cf)
    print("     (share_C medio congelado en el contrafactual: "
          "%.3f)" % wC)
    print()
    c1 = all(m["Bhat"] / m["BMs"] <= 10 for (_, m) in six)
    c2_agg = b_agg[0] <= 0.1
    c2_pie = b_PM[0] <= 0.1 and b_CM[0] <= 0.1 \
        and b_PS[0] <= 0.1 and b_CS[0] <= 0.1
    print("Criterio C1+C2 para Bhat5 en los seis puntos:")
    print("     C1 (Bhat5/B_M^(S) <= 10 en los 6): %s "
          "(min=%.3f, max=%.3f)"
          % ("PASS" if c1 else "FAIL",
             min(m["Bhat"] / m["BMs"] for (_, m) in six),
             max(m["Bhat"] / m["BMs"] for (_, m) in six)))
    print("     C2 agregado (pendiente <= 0.1): %.4f -> %s"
          % (b_agg[0], "PASS" if c2_agg else "FAIL"))
    print("     C2 por piezas (P y C, formas M y S, <= 0.1): "
          "P_M=%.4f P_S=%.4f C_M=%.4f C_S=%.4f -> %s"
          % (b_PM[0], b_PS[0], b_CM[0], b_CS[0],
             "PASS" if c2_pie else "FAIL"))
    print()
    print("RESUMEN: Bhat5 C1=%s C2_agregado=%s C2_piezas=%s; "
          "precios: agregado [%.3f, %.3f], C_M [%.3f, %.3f]"
          % ("PASS" if c1 else "FAIL",
             "PASS" if c2_agg else "FAIL",
             "PASS" if c2_pie else "FAIL",
             min(m["Bhat"] / m["BMs"] for (_, m) in six),
             max(m["Bhat"] / m["BMs"] for (_, m) in six),
             min(m["Dhat_C"] / m["BMm_C"] for (_, m) in six),
             max(m["Dhat_C"] / m["BMm_C"] for (_, m) in six)))
    print("PASS" if allok else "FAIL (ver lineas marcadas)")
    print("Tiempo total: %.1f s" % (time.time() - t00))
    return 0 if allok else 1


if __name__ == "__main__":
    sys.exit(main())
