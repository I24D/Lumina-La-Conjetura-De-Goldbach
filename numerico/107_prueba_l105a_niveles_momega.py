#!/usr/bin/env python3
# 107_prueba_l105a_niveles_momega.py  (Doc 107)
# Tarea 1: verificaciones numericas de la prueba escrita de
# L105(a)  (B_Delta <= C_Delta * T(N)):
#   - cota puntual (106.6): |Delta_r| <= phi(r) omega(r) ln N
#     para todo r en Hcal, via delta_r_orth (docs 97/100);
#   - B_Delta <= ln^2 N * M_w,  M_w = sum_{r en Hcal} phi(r)
#     omega(r)^2, en los 12 puntos de los docs 105/106;
#   - M_w por descomposicion en ramas estilo A100.1 (doc 101)
#     contra enumeracion directa;
#   - cadena del marco [T]: B_Delta <= N (ln N)^4/(4 (ln 2)^2),
#     cociente contra T(N) y N_0 explicito para C_Delta dado.
# Tarea 2: descomposicion por niveles de M_Omega:
#   masa de peso sum_d Omega_d^C phi(d)^2 por nivel
#   (tipo liso/conmutado x k = omega(d)), cotas clasicas por
#   nivel (gran criba via Parseval; BDH pagando el peso phi(d))
#   contra T(N_f) en el marco N_f = z^10, y calibracion con psi
#   de la masa real de M_Omega por nivel (Dhat^C por niveles).
# Solo stdlib. Importa sin cambios la maquinaria de los docs
# 102-106. MEDICION de apoyo a una prueba [T]; no es la prueba.
# Uso: python3 107_prueba_l105a_niveles_momega.py
import importlib.util
import math
import sys
import time

spec = importlib.util.spec_from_file_location(
    "s106", "106_precio_tau_conteo_omega.py")
s106 = importlib.util.module_from_spec(spec)
spec.loader.exec_module(s106)
s105 = s106.s105
s102 = s105.s104.s103.s102

Model97 = s106.Model97
build_H = s106.build_H
SegPsi = s106.SegPsi
M_r_E = s106.M_r_E
ols = s106.ols
measure105 = s106.measure105
Counting = s106.Counting
build_H106 = s106.build_H106
T_of_N = s106.T_of_N
tau_phi = s106.tau_phi
delta_r_orth = s105.delta_r_orth
is_prime_spf = s105.is_prime_spf
divisors_of = s105.divisors_of
phi_of = s105.phi_of
MU, SPF, SMALL_PRIMES = s105.MU, s105.SPF, s105.SMALL_PRIMES
pf_trial = s102.pf_trial
divisors_from_pf = s102.divisors_from_pf

C_CHEN = 0.039884
G_GRANDE = 0.0225
LN2 = math.log(2.0)


def level_key(d, z):
    """Nivel de d en D: (tipo, k), tipo 'L' liso (todo primo < z)
    o 'S' conmutado (exactamente un primo >= z, el de conmutacion
    s0); k = omega(d) = |supp(d)| (tau(d) = 2^k)."""
    pf = pf_trial(d, SMALL_PRIMES)
    big = [p for p in pf if p >= z]
    assert len(big) <= 1, "d con dos primos >= z: %d" % d
    return ("S" if len(big) == 1 else "L", len(pf))


def measure107(H, N, z, kappa=1.0):
    """Como measure105 pero con: Delta_r por conductor (ratios de
    la cota puntual), momentos de conteo M_w, M_phi, y Dhat^C
    descompuesto por niveles."""
    R0 = math.log(N) ** kappa
    HN = [r for r in H if math.gcd(r, N) == 1]
    divs_set = set()
    for r in HN:
        for d in divisors_of(r, SPF):
            divs_set.add(d)
    divs = sorted(divs_set)
    seg = SegPsi(N, divs)
    E = seg.E
    phi = {v: phi_of(v, SMALL_PRIMES)
           for v in set(divs) | set(HN)}
    Hcal = [r for r in HN if r > R0]
    lnN = math.log(N)
    # Omega_d^C por enumeracion (divisores de r compuestos)
    OmC = {}
    for r in Hcal:
        if is_prime_spf(r):
            continue
        tau = len(divisors_of(r, SPF))
        w = (tau - 1) / phi[r]
        for d in divisors_of(r, SPF):
            if d >= 2:
                OmC[d] = OmC.get(d, 0.0) + w
    # Delta_r, momentos de conteo, B_Delta
    BD = 0.0
    Mw = 0.0
    Mphi = 0.0
    Mw1 = 0.0
    ratio1 = 0.0
    ratio2 = 0.0
    for r in Hcal:
        Dr = delta_r_orth(r, N, MU, SPF)
        BD += Dr * Dr / phi[r]
        om = len(pf_trial(r, SMALL_PRIMES))
        Mw += phi[r] * om * om
        Mw1 += phi[r] * om
        Mphi += phi[r]
        if om >= 1:
            ratio1 = max(ratio1, abs(Dr) / (phi[r] * om * lnN))
            ratio2 = max(ratio2, abs(Dr) / (om * lnN))
    # Dhat^C por niveles (masa real de M_Omega con psi)
    lev = {}
    DhatC = 0.0
    for d, om in OmC.items():
        f2 = (phi[d] * E.get(d, 0.0)) ** 2
        key = level_key(d, z)
        lev[key] = lev.get(key, 0.0) + f2 * om
        DhatC += f2 * om
    return {"N": N, "nHN": len(HN), "nHcal": len(Hcal),
            "BD": BD, "Mw": Mw, "Mw1": Mw1, "Mphi": Mphi,
            "ratio1": ratio1, "ratio2": ratio2,
            "crude": lnN * lnN * Mw, "lev": lev, "DhatC": DhatC,
            "OmC": OmC, "Hcal": Hcal, "phi": phi}


def branch_Mw(sprimes, Vs, N, kappa=1.0):
    """M_w = sum_{r en Hcal} phi(r) omega(r)^2 por la
    descomposicion en ramas del doc 101 (A100.1):
    H = V0 (rama lisa, V+(D/s_min)) disjunto con
    union_s s*V+(D/s); en la rama s: phi(s h) = (s-1) phi(h),
    omega(s h) = omega(h) + 1 (s no divide a h: los primos de h
    son < z <= s). Filtros de Hcal aplicados por elemento."""
    R0 = math.log(N) ** kappa
    tot = 0.0
    smin = sprimes[0]
    for h in Vs[smin]:
        if h < 2 or h <= R0 or math.gcd(h, N) != 1:
            continue
        om = len(pf_trial(h, SMALL_PRIMES))
        tot += tau_phi(h)[1] * om * om
    for s in sprimes:
        if N % s == 0:
            continue
        for h in Vs[s]:
            r = s * h
            if r <= R0 or math.gcd(r, N) != 1:
                continue
            if h == 1:
                ph, omh = 1, 0
            else:
                ph = tau_phi(h)[1]
                omh = len(pf_trial(h, SMALL_PRIMES))
            tot += (s - 1) * ph * (omh + 1) ** 2
    return tot


def chain_bound(N):
    """Cota [T] del marco: B_Delta <= N (ln N)^4/(4 (ln 2)^2)."""
    return N * math.log(N) ** 4 / (4.0 * LN2 * LN2)


def ratio_chain(N):
    """chain_bound(N)/T(N) = (0.0225/(4 (ln2)^2)) (ln N)^8
    N^{-(1-2c)}, decreciente para ln N > 8/0.920232 (N > 5970)."""
    return (G_GRANDE / (4.0 * LN2 * LN2)) * math.log(N) ** 8 \
        * N ** (-(1.0 - 2.0 * C_CHEN))


def find_N0(target):
    lo, hi = 1.0e4, 1.0e20
    for _ in range(200):
        mid = math.sqrt(lo * hi)
        if ratio_chain(mid) > target:
            lo = mid
        else:
            hi = mid
    return hi


def main():
    allok = True
    t00 = time.time()
    print("== DOC 107: prueba de L105(a) + niveles de M_Omega ==")
    print("Solo stdlib. Las verificaciones son MEDICION de apoyo;")
    print("la prueba de L105(a) es el texto del doc, cadena [T].")
    print("Goldbach sigue abierta.")
    print()
    # ---------------- soportes ----------------
    m1 = Model97(3, 11, 60, 20000, "fixed")
    m2 = Model97(3, 11, 3000, 200000, "fixed")
    m3 = Model97(3, 3, 349, 200000, "fixed")
    supports = {"R1": (m1.H, 11, m1.sprimes, m1.Vs),
                "R3": (m3.H, 3, m3.sprimes, m3.Vs)}
    for z in (11, 13, 17, 19, 23, 29, 31, 37):
        D = z ** 5
        y = int(round(z ** (10.0 / 3.0)))
        Hz, sprimes, Vs = build_H106(3, z, y, D)
        supports["S%d" % z] = (Hz, z, sprimes, Vs)
    print("Soportes: R1=%d R3=%d S11=%d S13=%d S17=%d S19=%d "
          "S23=%d S29=%d S31=%d S37=%d"
          % tuple(len(supports[t][0]) for t in
                  ("R1", "R3", "S11", "S13", "S17", "S19", "S23",
                   "S29", "S31", "S37")))
    # K0: maquinaria
    D11 = 11 ** 5
    y11 = int(round(11 ** (10.0 / 3.0)))
    m_s11 = Model97(3, 11, y11, D11, "fixed")
    k0 = len(m2.H) == 880 and set(supports["S11"][0]) == set(m_s11.H)
    for z in (11, 13, 17):
        D = z ** 5
        y = int(round(z ** (10.0 / 3.0)))
        Href, _, _ = build_H(3, z, y, D, MU, SPF)
        k0 &= (list(Href) == supports["S%d" % z][0])
    allok &= k0
    print("K0 maquinaria (R2 |H|=880; S11 == Model97; build_H106 "
          "== build_H en z=11,13,17): %s"
          % ("PASS" if k0 else "FAIL"))
    print()
    # ---------------- Tarea 1: medida en 12 puntos ----------------
    plan = [(999982, ["S11", "S13", "S17", "S19", "S23"], "gen"),
            (1048576, ["S11", "S13", "S17", "S19", "S23"], "2^e")]
    pts = {}
    print("Doce puntos (medida completa): B_Delta, momentos de "
          "conteo y cotas:")
    for (N, tags, fam) in plan:
        for tag in tags:
            H, z, sprimes, Vs = supports[tag]
            t0 = time.time()
            m = measure107(H, N, z)
            key = "%s/%s" % (fam, tag)
            pts[key] = m
            print("     [%s] N=%d |Hcal|=%d B_Delta=%.6e "
                  "B_Delta/T=%.3e\n"
                  "        M_w=%.6e M_phi=%.6e crude/T=%.4e "
                  "Dhat^C=%.6e (%.1fs)"
                  % (key, N, m["nHcal"], m["BD"],
                     m["BD"] / T_of_N(N), m["Mw"], m["Mphi"],
                     m["crude"] / T_of_N(N), m["DhatC"],
                     time.time() - t0))
    for tag in ("R1", "R3"):
        H, z, sprimes, Vs = supports[tag]
        m = measure107(H, 999982, z)
        pts["gen/" + tag] = m
        print("     [gen/%s] N=999982 |Hcal|=%d B_Delta=%.6e "
              "B_Delta/T=%.3e M_w=%.6e crude/T=%.4e"
              % (tag, m["nHcal"], m["BD"],
                 m["BD"] / T_of_N(999982), m["Mw"],
                 m["crude"] / T_of_N(999982)))
    print()
    # K1: cota puntual (106.6) y B_Delta <= ln^2 N * M_w
    k1 = True
    print("Cota puntual (106.6): max |Delta_r|/(phi(r) omega(r) "
          "ln N) y max |Delta_r|/(omega(r) ln N) [la segunda es "
          "la cota impresa en el doc 97, NO usada en la prueba]:")
    for key in ("gen/R1", "gen/S11", "gen/S13", "gen/S17",
                "gen/S23", "2^e/S23"):
        m = pts[key]
        print("     [%s] ratio_phi=%.4f ratio_sin_phi=%.4f"
              % (key, m["ratio1"], m["ratio2"]))
        if key in ("gen/R1", "gen/S11", "gen/S13"):
            k1 &= m["ratio1"] <= 1.0 + 1e-9
    for key, m in pts.items():
        if m["BD"] > m["crude"] * (1 + 1e-9):
            k1 = False
            print("     K1 FALLO B_Delta > ln^2 N * M_w [%s]" % key)
    allok &= k1
    print("K1 cota puntual (106.6) en R1/S11/S13 y B_Delta <= "
          "ln^2 N * M_w en los 12 puntos: %s"
          % ("PASS" if k1 else "FAIL"))
    # K2: M_w por ramas (A100.1) contra enumeracion directa
    k2 = True
    print("M_w por ramas contra enumeracion directa:")
    for tag in ("R1", "S11", "S13", "S17", "S19", "S23"):
        H, z, sprimes, Vs = supports[tag]
        for N in (999982, 1048576):
            key = ("gen/" if N == 999982 else "2^e/") + tag
            if key not in pts:
                continue
            mb = branch_Mw(sprimes, Vs, N)
            md = pts[key]["Mw"]
            rel = abs(mb - md) / md if md > 0 else abs(mb - md)
            if rel > 1e-9:
                k2 = False
            print("     [%s] ramas=%.10e directa=%.10e rel=%.2e"
                  % (key, mb, md, rel))
    allok &= k2
    print("K2 M_w por ramas == directa: %s"
          % ("PASS" if k2 else "FAIL"))
    # K3: B_Delta/T en ventana (doc 106: max 2.16e-4)
    mx = max(m["BD"] / T_of_N(m["N"]) for m in pts.values())
    k3 = mx <= 2.5e-4
    allok &= k3
    print("K3 max B_Delta/T(N) en los 12 puntos = %.3e "
          "(doc 106: 2.16e-4) -> C_Delta = 1e-3 cubre la ventana "
          "con factor >= %.1f: %s"
          % (mx, 1e-3 / mx, "PASS" if k3 else "FAIL"))
    print()
    # ---------------- cadena del marco y N_0 ----------------
    print("Cadena [T]: B_Delta <= chain(N) = N (ln N)^4/"
          "(4 (ln 2)^2) = %.6f * N (ln N)^4"
          % (1.0 / (4.0 * LN2 * LN2)))
    print("    chain(N)/T(N) = %.8f * (ln N)^8 * N^{-(1-2c)}, "
          "1-2c = %.6f"
          % (G_GRANDE / (4.0 * LN2 * LN2), 1.0 - 2.0 * C_CHEN))
    k6 = True
    for N in (10 ** 6, 10 ** 9, 10 ** 12, 10 ** 15):
        ident = chain_bound(N) / T_of_N(N)
        if abs(ident - ratio_chain(N)) / ratio_chain(N) > 1e-9:
            k6 = False
    prev = None
    for e in range(4, 17):
        v = ratio_chain(10.0 ** e)
        if prev is not None and v > prev:
            k6 = False
        prev = v
    print("    ratio_chain: N=1e6: %.4e | 1e9: %.4e | 1e12: %.4e "
          "| 1e15: %.4e"
          % (ratio_chain(1e6), ratio_chain(1e9),
             ratio_chain(1e12), ratio_chain(1e15)))
    for target in (1.0, 1e-2, 1e-3):
        N0 = find_N0(target)
        res = ratio_chain(N0)
        if abs(res - target) / target > 0.02:
            k6 = False
        print("    N_0(C_Delta=%.0e) = %.4e  (ratio en N_0: %.3e)"
              % (target, N0, res))
    allok &= k6
    print("K6 identidad de la cadena + N_0 por biseccion: %s"
          % ("PASS" if k6 else "FAIL"))
    print()
    # ---------------- Tarea 2: niveles de M_Omega ----------------
    print("Niveles de M_Omega (conteo, sin psi): masa de peso "
          "sum_d Omega_d^C phi(d)^2 por nivel (tipo, k=omega(d));")
    print("cotas clasicas por nivel en el marco N_f = z^10 "
          "contra T(N_f). LS: W*_l (N_f + Q_l^2) N_f ln N_f;")
    print("BDH: W*_l Q_l (Q_l N_f ln N_f + N_f^2/(ln N_f)^2) "
          "[A=2, pagando phi(d) <= Q_l].")
    levdata = {}
    k4 = True
    k5 = True
    for z in (11, 13, 17, 19, 23, 29):
        H, _zz, sprimes, Vs = supports["S%d" % z]
        c = Counting(H, sprimes, Vs, z, 999982)
        Om, Uu, num, den, num2 = c.direct()
        stats = c.omega_stats()
        Nf = float(z) ** 10
        DA = z ** 5
        lnNf = math.log(Nf)
        TNf = T_of_N(Nf)
        levels = {}
        band = [0.0, 0.0, 0.0]
        tot = 0.0
        for d, om in Om.items():
            _, ph = tau_phi(d)
            key = level_key(d, z)
            L = levels.setdefault(key, {"mass": 0.0, "wstar": 0.0,
                                        "Q": 0, "cnt": 0})
            w = om * ph * ph
            L["mass"] += w
            L["wstar"] = max(L["wstar"], om * ph)
            L["Q"] = max(L["Q"], d)
            L["cnt"] += 1
            tot += w
            if d <= DA ** (1.0 / 3.0):
                band[0] += w
            elif d <= DA ** (2.0 / 3.0):
                band[1] += w
            else:
                band[2] += w
        if abs(sum(L["mass"] for L in levels.values()) - tot) \
                / tot > 1e-12:
            k4 = False
        # W* global desde niveles contra omega_stats del doc 106
        wstar_g = max(L["wstar"] for L in levels.values())
        if abs(wstar_g - stats["wstar"]) / stats["wstar"] > 1e-9:
            k5 = False
        for L in levels.values():
            if L["Q"] > DA:
                k5 = False
        rows = []
        totB = 0.0
        totLS = 0.0
        for key in sorted(levels):
            L = levels[key]
            LS = L["wstar"] * (Nf + L["Q"] ** 2) * Nf * lnNf
            BDH = L["wstar"] * L["Q"] * (
                L["Q"] * Nf * lnNf + Nf * Nf / lnNf ** 2)
            B = min(LS, BDH)
            totB += B
            totLS += LS
            rows.append((key, L, LS / TNf, BDH / TNf, B / TNf))
        levdata[z] = (levels, tot, rows, totB / TNf, band)
        print("  z=%d: |D_C|=%d masa total=%.6e W*=%.4f "
              "(doc 106 z=11,13,17: 9.52/15.00/31.00) "
              "N_f=%.3e T(N_f)=%.6e"
              % (z, len(Om), tot, wstar_g, Nf, TNf))
        print("     bandas de tamano: d<=D_A^{1/3}: %.3f | "
              "media: %.3f | d>D_A^{2/3}: %.3f"
              % (band[0] / tot, band[1] / tot, band[2] / tot))
        for key, L, eLS, eBDH, eB in rows:
            print("     nivel %s k=%d: |D_l|=%5d Q_l=%7d "
                  "W*_l=%8.4f masa=%.4f | exceso LS=%.3e "
                  "BDH=%.3e min=%.3e"
                  % (key[0], key[1], L["cnt"], L["Q"], L["wstar"],
                     L["mass"] / tot, eLS, eBDH, eB))
        print("     EXCESO TOTAL (suma de minimos) = %.4e ; "
              "solo-LS = %.4e"
              % (totB / TNf, totLS / TNf))
        if z in (11, 13, 17):
            ref = {11: 9.52, 13: 15.00, 17: 31.00}[z]
            if abs(wstar_g - ref) > 0.05:
                k5 = False
                print("     K5 DISCREPA W* z=%d: %.4f vs %.2f"
                      % (z, wstar_g, ref))
        if totB / TNf <= 1.0:
            k5 = False
            print("     K5 INESPERADO: la cota clasica total "
                  "alcanzaria T en z=%d" % z)
    allok &= k4
    allok &= k5
    print("K4 masas por nivel suman el total (rel < 1e-12): %s"
          % ("PASS" if k4 else "FAIL"))
    print("K5 W* por niveles == doc 106; Q_l <= D_A; la cota "
          "clasica total NO alcanza T (exceso > 1): %s"
          % ("PASS" if k5 else "FAIL"))
    print()
    # calibracion con psi: masa real de M_Omega por nivel
    print("Calibracion con psi (N=999982): fraccion de Dhat^C "
          "(= M_Omega medido) por nivel:")
    for tag in ("S11", "S13", "S17", "S19", "S23"):
        m = pts["gen/" + tag]
        fr = sorted(((k, v / m["DhatC"])
                     for k, v in m["lev"].items()),
                    key=lambda kv: -kv[1])
        print("     [%s] Dhat^C=%.6e | %s"
              % (tag, m["DhatC"], " ".join(
                  "%s%d:%.3f" % (k[0], k[1], v) for k, v in fr)))
    print()
    print("RESUMEN: K0=%s K1=%s K2=%s K3=%s K4=%s K5=%s K6=%s"
          % tuple("PASS" if v else "FAIL"
                  for v in (k0, k1, k2, k3, k4, k5, k6)))
    print("PASS" if allok else "FAIL (ver lineas marcadas)")
    print("Tiempo total: %.1f s" % (time.time() - t00))
    return 0 if allok else 1


if __name__ == "__main__":
    sys.exit(main())
