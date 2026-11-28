#!/usr/bin/env python3
# 114_desp_sign_banda.py  (Doc 114)
# (DESP-SIGN) por apareamiento y la banda degenerada por conteo.
#
# El doc 113 redujo (FL-DESP+) a la FL clasica sobre el modelo mas
#   (DESP-SIGN):  E(r) = S(r) - Shat(r) <= 0  para todo r>1 visible,
# mas la banda u_r < u_0 [P]. Aqui:
#   Add(r) := V^{(r)} \ Vhat^{(r)},  Rem(r) := Vhat^{(r)} \ V^{(r)},
#   E(r) = sum_{Add} a - sum_{Rem} a   [identidad].
# Apareamiento: phi(x) = x / q_min(x) (quitar el primo menor).
#   Lema 1 [T]: si x no esta en Vhat y omega(x) es par, x/q_min
#     tampoco (el primo menor ocupa la posicion omega(x), par: nunca
#     lleva la condicion que falla; las demas no cambian).
#   Lema 2 [T]: V^+(L) es cerrado por divisores (quitar un primo:
#     los primos mayores no se mueven; los menores suben una posicion
#     y la condicion vieja en i-1 domina a la nueva con cociente
#     p_i^3/(p * p_{i-1}^2) < 1). Luego x en V^{(r)} => x/q_min en
#     V^{(r)}.
#   => phi mapea Add^+ en Add^- y vale la identidad de flujo
#     E_Add(r) = -sum_{y en Add^-} |a(y)| (1 - R_A(y)),
#     R_A(y) = sum_{q: yq en Add^+} 1/(q-1) <= sum_{q<min(y)} 1/(q-1).
#   Mertens: sum_{5<=p<=67} 1/(p-1) = 0,99917 < 1 (con 71: 1,01346):
#   R_A <= 1 es [T] si z_A <= 71; el residuo es min(y) >= 71.
# Lado Rem: el modelo hereda al quitar el primo menor [T]; el fallo
# fusionado en un primo > q_min persiste [T]; el residuo son los
# fallos todos en primos <= q_min (vacio en ventana).
# Banda: S(r) <= 1 [flujo analogo sobre V^{(r)}, cerrado por
# divisores] => delta <= 1 - Pi; t^+ <= (1-Pi^2)/phi; la sub-banda
# Pi >= c0 cierra [T] contra sum Pi^2/phi <= Q_c; el resto es (BAND).
#
# Chequeos:
#   K0  cruce con el doc 113: |Add|=6509, |Rem|=1049 totales;
#       max E/Pi y E <= 0 en las 4113 instancias.
#   K1  E(r) = sum_Add a - sum_Rem a exacto; caracterizacion:
#       Add = A1 (evaporacion de paridad) / A2 (holgura) / mixto,
#       con el lema de holgura A2 verificado (primo de r bajo el
#       primo que falla); Rem = R1 (flip de paridad) / R2 (primo
#       de r), con R1 siempre en posicion propia par [T].
#   K2  Lema 2 [M]: d en V, p|d => d/p en V (0 violaciones);
#       V^{(r)} cerrado en d' (muestra).
#   K3  apareamiento Add: phi(x) en Add para todo x en Add^+
#       (1640); identidad de flujo exacta por instancia;
#       max R_A; cotas de Mertens con fracciones exactas.
#   K4  apareamiento Rem: psi(x) en Rem para todo x en Rem^-
#       (39; ninguno con omega=1); todo Rem^- tiene un fallo
#       fusionado en un primo > q_min (mecanismo [T] suficiente);
#       identidad de flujo; max R_R.
#   K5  S(r) <= 1 en las 4113 instancias [M]; delta <= 1 - Pi.
#   K6  banda u_r<u0 (u0=1,5 y 2): max_L sum_banda (1-Pi^2)/phi /Q_c
#       [M]; split Pi>=1/2 [T] vs Pi<1/2 [M]; en u_r<1, delta>0
#       => S(r)=1 exacto [M]; masa real de banda y conteos.
#   K7  estres z=53 (marco A, muestra): E<=0, aterrizajes,
#       max R_A / R_R.
# Uso: python3 114_desp_sign_banda.py | tee salida_doc114_local.txt
import importlib.util
import math
import sys
import time
from collections import Counter
from fractions import Fraction

spec = importlib.util.spec_from_file_location(
    "doc113", "113_fl_desplazado_positivo.py")
d113 = importlib.util.module_from_spec(spec)
spec.loader.exec_module(d113)

in_Vhat = d113.in_Vhat
Level111 = d113.Level111
FRAMES = d113.FRAMES
frame_levels = d113.d112.d111.frame_levels
primes_upto = d113.d112.d111.primes_upto


def lados(lv, r):
    """Clasifica cada d' coprimo con r: Add / Rem / Both / Neither,
    con (pf, mu, phi)."""
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


def model_failures(pf, L, r):
    """Posiciones impares m donde el modelo falla:
    r * q_m^2 * P_m >= L. Lista de (m, q_m)."""
    out = []
    prod = 1
    for i, p in enumerate(pf):
        m = i + 1
        prod *= p
        if m % 2 == 1 and r * p * p * prod >= L:
            out.append((m, p))
    return out


def merged_failures(pf, rpf, L):
    """Fallos de la condicion fusionada de r*d': lista de
    (primo, origen 'x'/'r', posicion propia)."""
    merged = sorted(list(pf) + list(rpf), reverse=True)
    own_x = {p: i + 1 for i, p in enumerate(pf)}
    own_r = {p: i + 1 for i, p in enumerate(rpf)}
    out = []
    prod = 1
    for i, p in enumerate(merged):
        j = i + 1
        prod *= p
        if j % 2 == 1 and p * p * prod >= L:
            if p in own_x and p not in own_r:
                out.append((p, "x", own_x[p]))
            elif p in own_r and p not in own_x:
                out.append((p, "r", own_r[p]))
    return out


_INST = None


def instancias():
    """Itera (tag, z, a, r, ld) sobre los 240 niveles x r visible>1,
    con ld = lados(lv, r) calculado una sola vez."""
    global _INST
    if _INST is None:
        _INST = []
        for (tag, z, lvs) in FRAMES:
            for a in lvs:
                for r in a.vis:
                    if r > 1:
                        _INST.append((tag, z, a, r, lados(a.lv, r)))
    return _INST


# ------------------------------------------------------------------ K0/K1
def check_K0K1():
    print("K0 cruce doc 113 y K1 descomposicion/caracterizacion:")
    ok = True
    tot_add = tot_rem = 0
    err_E = 0.0
    max_E_rel = -1e300
    nE_pos = 0
    add_class = Counter()
    rem_class = Counter()
    a2_sin_primo_bajo = 0
    r1_paridad = [0, 0]
    sin_fallo = 0
    for (tag, z, a, r, ld) in instancias():
        lv = a.lv
        L = a.L
        pfmap = {rr: pf for (rr, pf, mu, phi) in lv.enum}
        rpf = pfmap[r]
        Add = {d: v for d, v in ld.items() if v[0] == "Add"}
        Rem = {d: v for d, v in ld.items() if v[0] == "Rem"}
        tot_add += len(Add)
        tot_rem += len(Rem)
        E = sum(v[2] / v[3] for v in Add.values()) \
            - sum(v[2] / v[3] for v in Rem.values())
        err_E = max(err_E, abs(E - lv.delta[r]
                               + (d113.model_data(lv, r)[0] - lv.Pi[r])))
        # E contra S - Shat directo:
        Shat = d113.model_data(lv, r)[0]
        err_E = max(err_E, abs(E - (lv.S[r] - Shat)))
        max_E_rel = max(max_E_rel, E / lv.Pi[r])
        if E > 1e-12:
            nE_pos += 1
        for (dp, (lado_, dpf, dmu, dphi)) in Add.items():
            fails = model_failures(dpf, L, r)
            if not fails:
                sin_fallo += 1
                continue
            tipos = set()
            for (m, q) in fails:
                j = m + sum(1 for p in rpf if p > q)
                if j % 2 == 0:
                    tipos.add("A1")
                else:
                    tipos.add("A2")
                    if not any(p < q for p in rpf):
                        a2_sin_primo_bajo += 1
            add_class["+".join(sorted(tipos))] += 1
        for (dp, (lado_, dpf, dmu, dphi)) in Rem.items():
            fails = merged_failures(dpf, rpf, L)
            if not fails:
                sin_fallo += 1
                continue
            tipos = set()
            for (p, orig, mown) in fails:
                if orig == "x":
                    tipos.add("R1")
                    r1_paridad[0] += 1
                    if mown % 2 == 0:
                        r1_paridad[1] += 1
                else:
                    tipos.add("R2")
            rem_class["+".join(sorted(tipos))] += 1
    print("  |Add| total = %d (doc 113: 6509) | |Rem| total = %d "
          "(doc 113: 1049)" % (tot_add, tot_rem))
    print("  E = sum_Add a - sum_Rem a: err max = %.2e | max E/Pi = "
          "%+.2e | instancias con E>0: %d" % (err_E, max_E_rel, nE_pos))
    print("  Add por mecanismo: %s | sin fallo del modelo: %d"
          % (dict(add_class), sin_fallo))
    print("  Rem por mecanismo: %s | fallos R1 en posicion propia "
          "par: %d/%d [T]" % (dict(rem_class), r1_paridad[1],
                              r1_paridad[0]))
    print("  fallos A2 sin primo de r bajo el primo que falla: %d "
          "[T: la holgura exige uno]" % a2_sin_primo_bajo)
    ok &= (tot_add, tot_rem) == (6509, 1049) and err_E < 1e-9
    ok &= nE_pos == 0 and sin_fallo == 0 and a2_sin_primo_bajo == 0
    ok &= r1_paridad[0] == r1_paridad[1]
    print("  => %s" % ("PASS" if ok else "FAIL"))
    return ok


# ------------------------------------------------------------------ K2
def check_K2():
    print("K2 Lema 2 [M]: V^+(L) cerrado por divisores:")
    viol = 0
    tot = 0
    for (tag, z, lvs) in FRAMES:
        for a in lvs:
            lv = a.lv
            for (d, pf, mu, phi) in lv.V:
                tot += 1
                for p in pf:
                    if d // p not in lv.Vset:
                        viol += 1
    print("  d en V, p|d, d/p fuera de V: %d violaciones sobre "
          "%d elementos (todos sus primos)" % (viol, tot))
    ok = viol == 0
    print("  => la prueba del doc (desplazamiento de posiciones) "
          "concuerda con la enumeracion. %s" % ("PASS" if ok
                                                else "FAIL"))
    return ok


# ------------------------------------------------------------------ K3
def check_K3():
    print("K3 apareamiento Add: phi(x) = x/q_min(x):")
    ok = True
    n_par = 0
    mal = 0
    max_RA = 0.0
    max_RA_info = None
    err_flujo = 0.0
    for (tag, z, a, r, ld) in instancias():
        Add = {d: v for d, v in ld.items() if v[0] == "Add"}
        pre = {}
        for (dp, (lado_, dpf, dmu, dphi)) in Add.items():
            if dmu == 1 and dp > 1:
                n_par += 1
                yy = dp // dpf[-1]
                if yy not in Add or Add[yy][2] != -1:
                    mal += 1
                pre[yy] = pre.get(yy, 0.0) + 1.0 / (dpf[-1] - 1)
        E_add = sum(v[2] / v[3] for v in Add.values())
        flujo = 0.0
        for (dp, (lado_, dpf, dmu, dphi)) in Add.items():
            if dmu == -1:
                R = pre.get(dp, 0.0)
                flujo -= abs(dmu / dphi) * (1.0 - R)
                if R > max_RA:
                    max_RA = R
                    max_RA_info = (tag, z, a.L, r, dp)
        err_flujo = max(err_flujo, abs(E_add - flujo))
    s67 = sum(Fraction(1, p - 1) for p in primes_upto(67) if p >= 5)
    s71 = s67 + Fraction(1, 70)
    print("  socios comprobados: %d | fuera de Add^- : %d | "
          "identidad de flujo E_Add = -sum |a(y)|(1-R_A(y)): "
          "err max = %.2e" % (n_par, mal, err_flujo))
    print("  max R_A(y) = %.4f en %s" % (max_RA, max_RA_info))
    print("  Mertens exacto: sum_{5<=p<=67} 1/(p-1) = %.5f < 1 | "
          "con p=71: %.5f > 1" % (float(s67), float(s71)))
    ok &= mal == 0 and err_flujo < 1e-9 and max_RA <= 1.0
    ok &= s67 < 1 < s71
    print("  => (DESP-SIGN) en el lado Add queda reducido a "
          "R_A(y)<=1; probado si z_A<=71 (Mertens) y si min(y)<=67. "
          "%s" % ("PASS" if ok else "FAIL"))
    return ok


# ------------------------------------------------------------------ K4
def check_K4():
    print("K4 apareamiento Rem: psi(x) = x/q_min(x):")
    ok = True
    n_par = 0
    mal = 0
    omega1 = 0
    sin_fallo_alto = 0
    max_RR = 0.0
    err_flujo = 0.0
    for (tag, z, a, r, ld) in instancias():
        lv = a.lv
        pfmap = {rr: pf for (rr, pf, mu, phi) in lv.enum}
        rpf = pfmap[r]
        Rem = {d: v for d, v in ld.items() if v[0] == "Rem"}
        pre = {}
        for (dp, (lado_, dpf, dmu, dphi)) in Rem.items():
            if dmu == -1 and dp > 1:
                n_par += 1
                if len(dpf) == 1:
                    omega1 += 1
                yy = dp // dpf[-1]
                if yy not in Rem or Rem[yy][2] != 1:
                    mal += 1
                else:
                    pre[yy] = pre.get(yy, 0.0) + 1.0 / (dpf[-1] - 1)
                # persistencia: ?hay un fallo fusionado en un
                # primo > q_min(x)? (mecanismo [T] suficiente)
                qmin = dpf[-1]
                fails = merged_failures(dpf, rpf, a.L)
                if not any(p > qmin for (p, o, m) in fails):
                    sin_fallo_alto += 1
        E_rem = -sum(v[2] / v[3] for v in Rem.values())
        flujo = 0.0
        for (dp, (lado_, dpf, dmu, dphi)) in Rem.items():
            if dmu == 1:
                R = pre.get(dp, 0.0)
                flujo -= (dmu / dphi) * (1.0 - R)
                max_RR = max(max_RR, R)
        err_flujo = max(err_flujo, abs(E_rem - flujo))
    print("  socios comprobados: %d | fuera de Rem^+ : %d | con "
          "omega=1: %d | sin fallo fusionado en primo > q_min: %d"
          % (n_par, mal, omega1, sin_fallo_alto))
    print("  identidad de flujo E_Rem = -sum a(y)(1-R_R(y)): "
          "err max = %.2e | max R_R(y) = %.4f" % (err_flujo, max_RR))
    ok &= mal == 0 and omega1 == 0 and sin_fallo_alto == 0
    ok &= err_flujo < 1e-9 and max_RR <= 1.0
    print("  => el lado Rem queda reducido a la persistencia en "
          "primos <= q_min (vacia en ventana) y a R_R<=1. %s"
          % ("PASS" if ok else "FAIL"))
    return ok


# ------------------------------------------------------------------ K5
def check_K5():
    print("K5 S(r) <= 1 y delta(r) <= 1 - Pi(r):")
    mx = -1e300
    viol = 0
    for (tag, z, a, r, ld) in instancias():
        lv = a.lv
        mx = max(mx, lv.S[r])
        if lv.delta[r] > 1.0 - lv.Pi[r] + 1e-12:
            viol += 1
    print("  max S(r) = %.6f sobre las 4113 instancias | "
          "violaciones de delta <= 1-Pi: %d" % (mx, viol))
    ok = mx <= 1.0 + 1e-12 and viol == 0
    print("  => [M] en ventana; [T] si z_A<=71 por el flujo de K3 "
          "sobre V^{(r)} (cerrado por divisores, Lema 2). %s"
          % ("PASS" if ok else "FAIL"))
    return ok


# ------------------------------------------------------------------ K6
def check_K6():
    print("K6 banda u_r < u0: cota cruda (1-Pi^2)/phi por nivel:")
    ok = True
    for u0 in (1.5, 2.0):
        mx_cota = 0.0
        mx_menor = 0.0
        mx_real = 0.0
        mx_n = 0
        for (tag, z, lvs) in FRAMES:
            lnz = math.log(z)
            for a in lvs:
                lv = a.lv
                phimap = {rr: phi for (rr, pf, mu, phi) in lv.enum}
                cota = 0.0
                menor = 0.0
                real = 0.0
                n = 0
                for r in a.vis:
                    if r <= 1:
                        continue
                    u_r = math.log(a.L / r) / lnz if r < a.L else -1.0
                    if not (0 <= u_r < u0):
                        continue
                    n += 1
                    Pi = lv.Pi[r]
                    if lv.delta[r] > 0:
                        c = (1.0 - Pi * Pi) / phimap[r]
                        cota += c
                        if Pi < 0.5:
                            menor += c
                    if a.t[r] > 0:
                        real += a.t[r]
                mx_cota = max(mx_cota, cota / a.qc)
                mx_menor = max(mx_menor, menor / a.qc)
                mx_real = max(mx_real, real / a.qc)
                mx_n = max(mx_n, n)
        print("  u0=%.1f: max_L sum_banda cota/Q_c = %.4f (de la "
              "cual Pi<1/2: %.4f) | max_L masa real t^+/Q_c = %.4f "
              "| max r en banda por nivel = %d"
              % (u0, mx_cota, mx_menor, mx_real, mx_n))
    # u_r < 1: delta>0 => S(r) = 1 exacto (degenerado)
    n_banda1 = 0
    n_no_deg = 0
    for (tag, z, a, r, ld) in instancias():
        lv = a.lv
        u_r = math.log(a.L / r) / math.log(z) if r < a.L else -1.0
        if 0 <= u_r < 1.0 and lv.delta[r] > 0:
            n_banda1 += 1
            if abs(lv.S[r] - 1.0) > 1e-9:
                n_no_deg += 1
    print("  u_r<1: instancias con delta>0: %d | con S(r) != 1 "
          "(no degeneradas): %d" % (n_banda1, n_no_deg))
    ok &= n_no_deg == 0
    print("  => en ventana la banda suma O(Q_c) con la cota cruda "
          "[M]; la sub-banda Pi>=1/2 cierra [T] (<= 4 Q_c por "
          "sum Pi^2/phi <= Q_c); el residuo de marco es (BAND) "
          "sobre Pi<1/2. %s" % ("PASS" if ok else "FAIL"))
    return ok


# ------------------------------------------------------------------ K7
def check_K7():
    print("K7 estres z=53 (marco A, muestra de niveles/conductores):")
    z = 53
    D, y = z ** 5, int(z ** (10.0 / 3.0))
    lvs110 = frame_levels(z, 3, (), D, y, nsamp=10)
    n_inst = 0
    peor_E = -1e300
    max_RA = 0.0
    max_RR = 0.0
    mal_add = 0
    mal_rem = 0
    for lv110 in lvs110:
        a = Level111(z, 3, (), lv110.L)
        lv = a.lv
        vis = sorted(r for r in a.vis if r > 1)
        if len(vis) > 120:
            idx = sorted({round(i * (len(vis) - 1) / 119.0)
                          for i in range(120)})
            vis = [vis[i] for i in idx]
        for r in vis:
            n_inst += 1
            ld = lados(lv, r)
            Shat = d113.model_data(lv, r)[0]
            peor_E = max(peor_E, lv.S[r] - Shat)
            Add = {d: v for d, v in ld.items() if v[0] == "Add"}
            Rem = {d: v for d, v in ld.items() if v[0] == "Rem"}
            preA = {}
            for (dp, (lado_, dpf, dmu, dphi)) in Add.items():
                if dmu == 1 and dp > 1:
                    yy = dp // dpf[-1]
                    if yy not in Add:
                        mal_add += 1
                    preA[yy] = preA.get(yy, 0.0) + 1.0 / (dpf[-1] - 1)
            for yy, R in preA.items():
                if yy in Add:
                    max_RA = max(max_RA, R)
            preR = {}
            for (dp, (lado_, dpf, dmu, dphi)) in Rem.items():
                if dmu == -1 and dp > 1:
                    yy = dp // dpf[-1]
                    if yy not in Rem:
                        mal_rem += 1
                    else:
                        preR[yy] = preR.get(yy, 0.0) \
                            + 1.0 / (dpf[-1] - 1)
            for yy, R in preR.items():
                max_RR = max(max_RR, R)
    print("  instancias: %d | peor E = %+.2e | socios fuera: "
          "Add %d, Rem %d | max R_A = %.4f | max R_R = %.4f"
          % (n_inst, peor_E, mal_add, mal_rem, max_RA, max_RR))
    ok = peor_E <= 1e-9 and mal_add == 0 and mal_rem == 0
    ok &= max_RA <= 1.0 and max_RR <= 1.0
    print("  => el apareamiento y los flujos se sostienen fuera de "
          "los marcos del doc 113. %s" % ("PASS" if ok else "FAIL"))
    return ok


# ------------------------------------------------------------------ main
def main():
    t0 = time.time()
    print("== DOC 114: (DESP-SIGN) por apareamiento + banda ==")
    print("Solo stdlib. Definiciones de los docs 109-113 por import.")
    print()
    allok = True
    allok &= check_K0K1()
    print()
    allok &= check_K2()
    print()
    allok &= check_K3()
    print()
    allok &= check_K4()
    print()
    allok &= check_K5()
    print()
    allok &= check_K6()
    print()
    allok &= check_K7()
    print()
    print("Tiempo total: %.1f s" % (time.time() - t0))
    print("PASS" if allok else "FAIL")
    return 0 if allok else 1


if __name__ == "__main__":
    sys.exit(main())
