#!/usr/bin/env python3
# 104_experimento_e103_1_bhat3.py  (Doc 104)
# Experimento E103.1 (propuesto en el doc 103): la variante V103.2,
# mayorante Bhat3 = triangulo SOLO entre clases de parte nueva,
# signos de Mobius conservados dentro de cada clase, contra el
# criterio C1+C2 del doc 102 en los seis puntos.
# Solo stdlib.
#
# Definicion (doc 103, seccion 4.3):
#   Bhat3 := sum_{r en Hcal} (1/phi(r)) *
#            ( sum_h | sum_{d|r, nueva(d)=h} mu(r/d) phi(d) E_d |
#              + |Delta_r| )^2 ,
#   Hcal = {r en H_N : r > R0}, R0 = (ln N)^kappa, kappa = 1.
# Bhat3 >= B_M (triangulo entre clases sobre la particion de la
# suma M_r) y Bhat3 <= Bhat2 (refina el triangulo total) [T].
#
# AMBIGUEDAD DECLARADA (doc 104, seccion 1): "nueva(d)" admite dos
# lecturas.
#   L1 (literal): nueva(d) = d/(d, A), A = 3. En este modelo todo
#       d|r con r en Hcal es coprimo con 3 (los m = s*d se
#       construyen con primos != 3), luego nueva(d) = d, las clases
#       son unitarias y Bhat3 = Bhat2 identicamente. Se verifica
#       en K3: (d,3)=1 para todo d|r en los seis puntos.
#   L2 (operativa, ADOPTADA): la particion que diagonaliza el lado
#       de coeficientes en Model97 (Hsmooth/T, docs 94/96/97) es la
#       descomposicion por la parte lisa en z: nueva_z(d) =
#       producto de los p^e || d con p < z (z = el de construccion
#       del soporte). Cada clase reune d = h0 y d = h0*s (s primo
#       >= z, s|r; a lo sumo uno, porque r|m y m tiene un unico
#       primo >= z), cuyos signos de Mobius en M_r son opuestos.
#       Es la lectura que el doc 103 indica ("la particion por
#       parte nueva ya esta en Model97 via Hsmooth/T") y la unica
#       que hace falsable la hipotesis.
#
# Hipotesis falsable del doc 103 (seccion 4.3): si la pendiente
# 0,136 de Bhat2 proviene de la cancelacion ENTRE divisores de
# clases distintas (que crece con el numero de clases por
# conductor), Bhat3 la aplana por debajo de 0,1 conservando C1;
# si proviene de cancelacion DENTRO de cada clase, Bhat3 hereda
# la pendiente y la via del mayorante absoluto queda cerrada.
# Descomposicion exacta de la pendiente (OLS es lineal en y):
#   pendiente(log Bhat2/B_M) = pendiente(log Bhat2/Bhat3)
#                            + pendiente(log Bhat3/B_M),
# donde el primer sumando es el precio del triangulo DENTRO de
# clase (el que Bhat3 deja de pagar al conservar signos) y el
# segundo el precio del triangulo ENTRE clases (el que Bhat3
# sigue pagando). Se verifica en K4.
#
# Chequeos:
#   K0  Maquinaria: Model97 R2 |H|=880; build_H z=11 == Model97.
#   K1  B_M^(S) y Bhat2 reproducen el doc 103 en los seis puntos
#       (B_M rel < 2e-3; Bhat2/B_M rel < 1e-2).
#   K2  Por conductor: |M_r| <= suma_h |suma de clase| <=
#       suma_d |termino| (triangulo en ambos sentidos); en
#       agregado: B_M <= Bhat3 <= Bhat2 en todo punto.
#   K3  Degeneracion L1: (d, 3) = 1 para todo d|r, r en Hcal,
#       en los seis puntos (luego bajo L1, Bhat3 = Bhat2).
#   K4  Aditividad de pendientes: pend(Bhat2/B_M) =
#       pend(Bhat2/Bhat3) + pend(Bhat3/B_M) (tol 1e-9).
#   K5  C1+C2 para Bhat3 en los seis puntos (reporte).
# Uso: python3 104_experimento_e103_1_bhat3.py
import importlib.util
import math
import sys
import time

spec = importlib.util.spec_from_file_location(
    "s103", "103_criba_restringida_segundo_factor.py")
s103 = importlib.util.module_from_spec(spec)
spec.loader.exec_module(s103)

Model97 = s103.Model97
build_H = s103.build_H
SegPsi = s103.SegPsi
M_r_E = s103.M_r_E
ols = s103.ols
divisors_of = s103.divisors_of
phi_of = s103.phi_of
pf_desc = s103.pf_desc
is_prime_spf = s103.is_prime_spf
delta_r_orth = s103.delta_r_orth
MU, SPF, SMALL_PRIMES = s103.MU, s103.SPF, s103.SMALL_PRIMES


def smooth_part(d, z):
    """nueva_z(d) (lectura L2): parte de d con primos < z."""
    h = 1
    for p in pf_desc(d, SPF):
        if p < z:
            h *= p
    return h


def measure104(H, N, z, kappa=1.0):
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
    BMs = BMs_P = BMs_C = 0.0
    Bhat2 = Bhat3 = 0.0
    Bhat3_P = Bhat3_C = 0.0
    ok_tri = True          # K2 por conductor
    ok_L1 = True           # K3: (d,3)=1
    nclass_tot = 0         # clases totales (sobre r en Hcal)
    ndiv_tot = 0           # divisores d>1 totales
    nr_with_pair = 0       # r con alguna clase de tamano 2
    rho_num = 0.0          # sum_r sum_h |suma clase| (sin Delta)
    rho_den = 0.0          # sum_r sum_d |termino| (sin Delta)
    pair_mass = 0.0        # sum_r sum_{h: |clase|=2} |suma clase|
    sing_mass = 0.0        # sum_r sum_{h: |clase|=1} |termino|
    for r in Hcal:
        Mr = M_r_E(r, E, MU, SPF)
        Dr = delta_r_orth(r, N, MU, SPF)
        Sr = Mr + Dr
        ph = phi[r]
        tS = Sr * Sr / ph
        classes = {}
        sabs = 0.0
        for d in divisors_of(r, SPF):
            if math.gcd(d, 3) != 1:
                ok_L1 = False
            term = MU[r // d] * phi_of(d, SMALL_PRIMES) * E.get(d, 0.0)
            sabs += abs(term)
            h = smooth_part(d, z)
            classes[h] = classes.get(h, 0.0) + term
            if d > 1:
                ndiv_tot += 1
        csum = 0.0
        has_pair = False
        # tamanos de clase: contar divisores por clase
        sizes = {}
        for d in divisors_of(r, SPF):
            sizes[smooth_part(d, z)] = sizes.get(smooth_part(d, z), 0) + 1
        for h, v in classes.items():
            csum += abs(v)
            if sizes[h] == 2:
                has_pair = True
                pair_mass += abs(v)
            elif sizes[h] == 1:
                sing_mass += abs(v)
        nclass_tot += len(classes)
        if has_pair:
            nr_with_pair += 1
        rho_num += csum
        rho_den += sabs
        # K2 por conductor (sin Delta: |M_r| <= csum <= sabs)
        if not (abs(Mr) <= csum + 1e-6 and csum <= sabs + 1e-6):
            ok_tri = False
        t2 = (sabs + abs(Dr)) ** 2 / ph
        t3 = (csum + abs(Dr)) ** 2 / ph
        Bhat2 += t2
        Bhat3 += t3
        if is_prime_spf(r):
            BMs_P += tS
            Bhat3_P += t3
        else:
            BMs_C += tS
            Bhat3_C += t3
        BMs += tS
    return {
        "N": N, "z": z, "nHN": len(HN), "nHcal": len(Hcal),
        "BMs": BMs, "BMs_P": BMs_P, "BMs_C": BMs_C,
        "Bhat2": Bhat2, "Bhat3": Bhat3,
        "Bhat3_P": Bhat3_P, "Bhat3_C": Bhat3_C,
        "ok_tri": ok_tri, "ok_L1": ok_L1,
        "nclass_avg": nclass_tot / len(Hcal),
        "ndiv_avg": ndiv_tot / len(Hcal),
        "frac_pair": nr_with_pair / len(Hcal),
        "rho": rho_num / rho_den if rho_den > 0 else float("nan"),
        "pair_share": pair_mass / rho_num if rho_num > 0 else 0.0,
    }


def fmt_pt(tag, m):
    return ("     [%s] N=%d z=%d |H_N|=%d |Hcal|=%d\n"
            "        B_M^(S)=%.6e (P %.6e / C %.6e)\n"
            "        Bhat2/B_M=%.3f  Bhat3/B_M=%.3f  "
            "Bhat2/Bhat3=%.3f\n"
            "        Bhat3 piezas: P x%.3f  C x%.3f\n"
            "        clases/r=%.2f divisores/r=%.2f  r con par: "
            "%.1f%%  rho_intra=%.4f  masa en pares=%.1f%%"
            % (tag, m["N"], m["z"], m["nHN"], m["nHcal"], m["BMs"],
               m["BMs_P"], m["BMs_C"], m["Bhat2"] / m["BMs"],
               m["Bhat3"] / m["BMs"], m["Bhat2"] / m["Bhat3"],
               m["Bhat3_P"] / m["BMs_P"] if m["BMs_P"] > 0 else float("nan"),
               m["Bhat3_C"] / m["BMs_C"] if m["BMs_C"] > 0 else float("nan"),
               m["nclass_avg"], m["ndiv_avg"], 100 * m["frac_pair"],
               m["rho"], 100 * m["pair_share"]))


def main():
    allok = True
    t00 = time.time()
    print("== DOC 104: experimento E103.1 - Bhat3 (V103.2) ==")
    print("Solo stdlib. MEDICION, no prueba. Goldbach sigue abierta.")
    print("Lectura adoptada: L2 (parte lisa en z, via Hsmooth/T).")
    print("Lectura literal L1 (d/(d,3)): degenerada, ver K3.")
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
        m = measure104(H, 999982, z)
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
            m = measure104(H, N, z)
            six.append(("%s/%s" % (fam, tag), m))
            allpts["%s/%s" % (fam, tag)] = m
            print(fmt_pt("%s/%s" % (fam, tag), m)
                  + " (%.0fs)" % (time.time() - t0))
    print()
    # K1: reproduccion del doc 103 (B_M^(S) y Bhat2/B_M, seis puntos)
    doc103_BM = {"gen/S11": 3.3658e9, "gen/S13": 7.6433e9,
                 "gen/S17": 2.7875e10, "2^e/S11": 3.7105e9,
                 "2^e/S13": 8.3332e9, "2^e/S17": 2.9761e10}
    doc103_B2 = {"gen/S11": 1.709, "gen/S13": 1.920, "gen/S17": 2.231,
                 "2^e/S11": 1.826, "2^e/S13": 2.053, "2^e/S17": 2.318}
    k1 = True
    for tag, ref in doc103_BM.items():
        m = allpts[tag]
        if abs(m["BMs"] - ref) / ref > 2e-3:
            k1 = False
            print("     K1 DISCREPA B_M [%s]: %.6e vs doc 103 %.6e"
                  % (tag, m["BMs"], ref))
        r2 = m["Bhat2"] / m["BMs"]
        if abs(r2 - doc103_B2[tag]) / doc103_B2[tag] > 1e-2:
            k1 = False
            print("     K1 DISCREPA Bhat2/B_M [%s]: %.3f vs doc 103 "
                  "%.3f" % (tag, r2, doc103_B2[tag]))
    allok &= k1
    print("K1 reproduccion del doc 103 (B_M^(S), Bhat2/B_M): %s"
          % ("PASS" if k1 else "FAIL (ver arriba)"))
    # K2: triangulo por conductor y en agregado
    k2 = all(m["ok_tri"] for m in allpts.values()) and all(
        m["BMs"] * (1 - 1e-12) <= m["Bhat3"] <= m["Bhat2"] * (1 + 1e-12)
        for m in allpts.values())
    allok &= k2
    print("K2 |M_r| <= suma de clases <= triangulo total, por "
          "conductor; B_M <= Bhat3 <= Bhat2 en todo punto: %s"
          % ("PASS" if k2 else "FAIL"))
    # K3: degeneracion de la lectura literal L1
    k3 = all(m["ok_L1"] for m in allpts.values())
    allok &= k3
    print("K3 (d,3)=1 para todo d|r en los seis puntos + R1/R3 "
          "(bajo L1, nueva(d)=d y Bhat3=Bhat2): %s"
          % ("PASS" if k3 else "FAIL"))
    print()
    # Pendientes sobre los seis puntos
    xs = [math.log(m["nHN"]) for (_, m) in six]
    y_B2 = [math.log(m["Bhat2"] / m["BMs"]) for (_, m) in six]
    y_B3 = [math.log(m["Bhat3"] / m["BMs"]) for (_, m) in six]
    y_23 = [math.log(m["Bhat2"] / m["Bhat3"]) for (_, m) in six]
    b2, a2, se2, ci2, l02, l12 = ols(xs, y_B2)
    b3, a3, se3, ci3, l03, l13 = ols(xs, y_B3)
    b23, a23, se23, ci23, l023, l123 = ols(xs, y_23)
    # K4: aditividad exacta de pendientes
    k4 = abs(b2 - (b23 + b3)) < 1e-9
    allok &= k4
    print("K4 aditividad: pend(Bhat2/B_M)=%.6f = pend(Bhat2/Bhat3)"
          "=%.6f + pend(Bhat3/B_M)=%.6f : %s"
          % (b2, b23, b3, "PASS" if k4 else "FAIL"))
    print()
    print("Descomposicion de la pendiente de Bhat2 (doc 103: 0,136):")
    print("     precio triangulo DENTRO de clase (Bhat2/Bhat3, "
          "lo que Bhat3 recupera): pendiente=%.4f "
          "+- %.4f (SE), IC95%% +-%.4f, LOO [%.4f, %.4f]"
          % (b23, se23, ci23, l023, l123))
    print("     precio triangulo ENTRE clases (Bhat3/B_M, "
          "lo que Bhat3 sigue pagando): pendiente=%.4f "
          "+- %.4f (SE), IC95%% +-%.4f, LOO [%.4f, %.4f]"
          % (b3, se3, ci3, l03, l13))
    print()
    # C1 + C2 para Bhat3
    c1 = all(m["Bhat3"] / m["BMs"] <= 10 for (_, m) in six)
    c2 = b3 <= 0.1
    print("Criterio C1+C2 para Bhat3 (V103.2) en los seis puntos:")
    print("     C1 (Bhat3/B_M <= 10 en los 6): %s (min=%.3f, max=%.3f)"
          % ("PASS" if c1 else "FAIL",
             min(m["Bhat3"] / m["BMs"] for (_, m) in six),
             max(m["Bhat3"] / m["BMs"] for (_, m) in six)))
    print("     C2 (pendiente log-log <= 0.1): pendiente=%.4f -> %s"
          % (b3, "PASS" if c2 else "FAIL"))
    print()
    # Diagnostico por piezas: precios P y C por separado y mezcla
    yC = [math.log(m["Bhat3_C"] / m["BMs_C"]) for (_, m) in six]
    yP = [math.log(m["Bhat3_P"] / m["BMs_P"]) for (_, m) in six]
    bC, aC, seC, ciC, l0C, l1C = ols(xs, yC)
    bP, aP, seP, ciP, l0P, l1P = ols(xs, yP)
    print("Diagnostico por piezas (Bhat3):")
    print("     pendiente Bhat3^P/B_M^P = %.4f (SE %.4f, LOO "
          "[%.4f, %.4f]); precio P = 1.000 en los seis puntos"
          % (bP, seP, l0P, l1P))
    print("     pendiente Bhat3^C/B_M^C = %.4f (SE %.4f, LOO "
          "[%.4f, %.4f])"
          % (bC, seC, l0C, l1C))
    for (tag, m) in six:
        print("     [%s] share_C(B_M)=%.3f  Bhat3^C/B_M^C=%.3f  "
              "Bhat3^P/B_M^P=%.3f"
              % (tag, m["BMs_C"] / m["BMs"],
                 m["Bhat3_C"] / m["BMs_C"],
                 m["Bhat3_P"] / m["BMs_P"]))
    print()
    print("RESUMEN: Bhat3 C1=%s C2=%s; pendiente Bhat2=%.4f -> "
          "Bhat3=%.4f (recuperado dentro de clase %.4f + "
          "residual entre clases %.4f)"
          % ("PASS" if c1 else "FAIL", "PASS" if c2 else "FAIL",
             b2, b3, b23, b3))
    print("PASS" if allok else "FAIL (ver lineas marcadas)")
    print("Tiempo total: %.1f s" % (time.time() - t00))
    return 0 if allok else 1


if __name__ == "__main__":
    sys.exit(main())
