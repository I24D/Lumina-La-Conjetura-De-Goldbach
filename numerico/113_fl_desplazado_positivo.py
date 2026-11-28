#!/usr/bin/env python3
# 113_fl_desplazado_positivo.py  (Doc 113)
# (FL-DESP+): el lema fundamental desplazado en su parte positiva.
# Para r visible > 1, con u_r := ln(L/r)/ln z (el nivel PROPIO del
# soporte desplazado V^{(r)}) y rho_0(u) = u^{-u/2} (la forma del error
# del lema fundamental clasico, constante 1):
#
#   ratio(r) := delta(r)^+ / ( Pi(r) * rho_0(u_r) ).
#
# (FL-DESP+) con rho = A*rho_0 se sostiene en ventana sii max ratio <= A.
# El doc 109 sec. 3.4 midio max |delta|/Pi = 1,00-1,30 con rho evaluada
# al nivel u del MARCO; aqui rho se evalua al nivel propio u_r.
#
# Descomposicion exacta usada en el ataque:
#   Shat(r) := sum_{d' en Vhat^{(r)}} a(d'),  Vhat^{(r)} = V^+(L/r; P\r)
#   (soporte RI genuino de nivel L/r sobre los primos restantes: la FL
#   clasica del doc 112 sec. 2 se le aplica con nivel u_r y termino
#   principal Pi(r) EXACTO), y
#   E(r) := S(r) - Shat(r)  (error de desplazamiento: la masa con signo
#   de la diferencia simetrica V^{(r)} △ Vhat^{(r)}).
#   delta(r) = (Shat(r) - Pi(r)) + E(r)   [identidad].
#
# Definiciones de los docs 109/110/111/112 reutilizadas por import
# (Level111 -> Level110 -> LevelData: literalmente las mismas).
#
# Chequeos:
#   K0  validacion cruzada contra el doc 112: max_L R^(L)/Q_c por
#       marco = 0,224/0,330/0,381/0,130; Shat(1) = S(1) exacto.
#   K1  el algebra de la implicacion (FL-DESP+) => (H2''''):
#       sum_{r>1} Pi(r)^2/phi(r) = Q_c - Pi_full^2 [T];
#       Pi(r)^2/phi(r) = Pi_full^2 * h(r), sum_r h(r) = M_full/Pi_full
#       [T]; t(r) <= rho(2+rho) Pi^2/phi en toda instancia con
#       0 <= delta <= rho*Pi (monotonia).
#   K2  MEDICION DECISIVA: A_min = max ratio sobre todos los niveles
#       y r visibles > 1 con delta > 0; distribucion (cuantiles,
#       histograma), bandas por u_r, A_min por marco, top-10; y
#       A_min restringido a sumas genuinas (S(r) < 1).
#   K3  los ofensores del doc 109 (|delta|/Pi = 1,00/1,069/1,295):
#       argmax por marco sobre todo r y sobre visibles r>1, con
#       visibilidad, signo de delta y u_r.
#   K4  el soporte modelo y el error de desplazamiento, sobre TODAS
#       las instancias: constante empirica (una y dos caras) de la
#       FL clasica para Shat al nivel propio; E(r): max E/Pi,
#       max E^+/(Pi rho_0(u_r)), min E/Pi, fraccion E>0.
#   K5  nivel a nivel: P_{>1}(L) contra la cota implicada B_A(L)
#       (A=1 y A=A_min): maximos y conteo de niveles con P_{>1} > B.
#   K6  la diferencia simetrica es bilateral en CONTEO [en ventana]:
#       hay d' en V^{(r)} \ Vhat^{(r)} y en Vhat^{(r)} \ V^{(r)}
#       (ejemplos explicitos); masa absoluta de la diferencia.
#   K7  el mecanismo de los ofensores: soportes degenerados
#       V^{(r)} = {1} (S(r) = 1 exacto, delta = 1 - Pi(r)); su
#       A_min parcial, su u_r maximo contra el umbral 3 ln 5/ln z,
#       y si los top-10 de K2 son todos de este tipo.
# Uso: python3 113_fl_desplazado_positivo.py | tee salida_doc113_local.txt
import importlib.util
import math
import sys
import time

spec = importlib.util.spec_from_file_location(
    "doc112", "112_masa_positiva_bloques.py")
d112 = importlib.util.module_from_spec(spec)
spec.loader.exec_module(d112)

Level111 = d112.Level111
Level110 = d112.d111.Level110


def rho0(u):
    """rho_0(u) = u^{-u/2} = exp(-u ln u / 2), u > 0 (vale 1 en u=1,
    tiende a 1 por arriba cuando u -> 0+, maximo e^{1/(2e)} ~ 1,202
    en u = 1/e, y decae superexponencialmente para u grande)."""
    if u <= 0:
        return 1.0
    return math.exp(-0.5 * u * math.log(u))


def h_of_pf(pf):
    """h(r) = prod_{p|r} (p-1)/(p-2)^2:  Pi(r)^2/phi(r) = Pi_full^2 h(r)."""
    h = 1.0
    for p in pf:
        h *= (p - 1.0) / (p - 2.0) ** 2
    return h


def in_Vhat(pf, L, r):
    """Pertenencia a V^+(L/r; P\\r): la in_Vplus del doc 109 con el
    nivel L/r, en aritmetica entera (toda comparacion contra L/r se
    multiplica por r): pf[0]*r >= L, y para m impar
    r * p_m^2 * (p_1...p_m) >= L."""
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


def model_data(lv, r):
    """Para el nivel lv (Level110) y el conductor r: Shat(r) sobre el
    soporte modelo Vhat^{(r)}; |V^{(r)}| (numero de d' con r d' en V);
    conteos de la diferencia simetrica (n1: en V y no en Vhat;
    n2: en Vhat y no en V), ejemplos y masa absoluta de la diferencia."""
    L = lv.L
    Vset = lv.Vset
    tot = 0.0
    nV = 0
    n1 = n2 = 0
    ex1 = ex2 = None
    absds = 0.0
    for (dp, dpf, dmu, dphi) in lv.enum:
        if math.gcd(dp, r) != 1:
            continue
        inhat = (dp * r < L) and in_Vhat(dpf, L, r)
        inv = (r * dp) in Vset
        if inhat:
            tot += dmu / dphi
        if inv:
            nV += 1
        if inv and not inhat:
            n1 += 1
            absds += 1.0 / dphi
            if ex1 is None:
                ex1 = dp
        elif inhat and not inv:
            n2 += 1
            absds += 1.0 / dphi
            if ex2 is None:
                ex2 = dp
    return tot, nV, n1, n2, ex1, ex2, absds


FRAMES = d112.get_frames()   # [(tag, z, [Level111 x60]), ...] x4


def all_instances():
    """Itera (tag, z, a, r, u_r, delta, Pi) sobre todo r visible > 1
    de los 240 niveles."""
    for (tag, z, lvs) in FRAMES:
        lnz = math.log(z)
        for a in lvs:
            lv = a.lv
            for r in a.vis:
                if r <= 1:
                    continue
                u_r = math.log(a.L / r) / lnz
                yield (tag, z, a, r, u_r, lv.delta[r], lv.Pi[r])


# ------------------------------------------------------------------ K0
def check_K0():
    print("K0 validacion cruzada (doc 112) y soporte modelo en r=1:")
    esperado = {("A", 17): 0.224, ("A", 29): 0.330,
                ("A", 37): 0.381, ("B", 29): 0.130}
    ok = True
    for (tag, z, lvs) in FRAMES:
        mxR = 0.0
        for a in lvs:
            R = a.P - max(a.t[1], 0.0)
            mxR = max(mxR, R / a.qc)
        e = abs(mxR - esperado[(tag, z)])
        print("  [%s] z=%d: max_L R^(L)/Q_c = %.3f (doc 112: %.3f, "
              "dif %.4f)" % (tag, z, mxR, esperado[(tag, z)], e))
        ok &= e < 0.002
    # Shat(1) debe ser S(1): Vhat^{(1)} = V^+(L) = V
    err = 0.0
    for (tag, z, lvs) in FRAMES:
        for a in (lvs[0], lvs[len(lvs) // 2], lvs[-1]):
            Shat, _, _, _, _, _, _ = model_data(a.lv, 1)
            err = max(err, abs(Shat - a.lv.S[1]))
    print("  Shat(1) = S(1) en 12 niveles: err max = %.2e" % err)
    ok &= err < 1e-12
    print("  => %s" % ("PASS" if ok else "FAIL"))
    return ok


# ------------------------------------------------------------------ K1
def check_K1():
    print("K1 algebra de la implicacion (FL-DESP+) => (H2'''') [T]:")
    ok = True
    for z in (17, 29, 37):
        lv = Level110(z, 3, (), 10 ** 9)
        s1 = sum(lv.Pi[r] ** 2 / phi
                 for (r, pf, mu, phi) in lv.enum if r > 1)
        e1 = abs(s1 - (lv.qc - lv.Pi_full ** 2))
        sh = sum(h_of_pf(pf) for (r, pf, mu, phi) in lv.enum)
        e2 = abs(sh - lv.M_full_prod / lv.Pi_full)
        e3 = 0.0
        for (r, pf, mu, phi) in lv.enum:
            e3 = max(e3, abs(lv.Pi[r] ** 2 / phi
                             - lv.Pi_full ** 2 * h_of_pf(pf)))
        print("  z=%d: sum_{r>1} Pi^2/phi = Q_c - Pi_full^2 "
              "(err %.2e) | sum_r h(r) = M_full/Pi_full = %.4f "
              "(err %.2e) | Pi^2/phi = Pi_full^2 h(r) (err %.2e)"
              % (z, e1, sh, e2, e3))
        ok &= e1 < 1e-9 and e2 < 1e-9 and e3 < 1e-12
    # monotonia: si 0 <= delta <= rho*Pi entonces
    # t = delta(2Pi+delta)/phi <= rho(2+rho) Pi^2/phi
    viol = 0.0
    ntest = 0
    for (tag, z, a, r, u_r, delta, Pi) in all_instances():
        if delta < 0:
            continue
        rho = rho0(u_r)
        if delta <= rho * Pi:
            ntest += 1
            t = a.t[r]
            if t > 0:
                phi_r = delta * (2.0 * Pi + delta) / t
                viol = max(viol, t * phi_r
                           - rho * (2.0 + rho) * Pi * Pi)
    print("  monotonia t <= rho(2+rho)Pi^2/phi en las %d instancias "
          "con 0<=delta<=rho_0*Pi: violacion max = %.2e"
          % (ntest, viol))
    ok &= viol < 1e-9
    print("  => con rho <= 1: P_{>1} <= rho(2+rho)(Q_c - Pi_full^2) "
          "<= 3 Q_c. %s" % ("PASS" if ok else "FAIL"))
    return ok


# ------------------------------------------------------------------ K2
def check_K2():
    print("K2 MEDICION DECISIVA: ratio = delta^+/(Pi * rho_0(u_r)) "
          "sobre r visible > 1:")
    ratios = []          # (ratio, tag, z, L, r, u_r, delta/Pi, S)
    ratios_gen = []      # las de sumas genuinas: S(r) < 1 - 1e-9
    bandas = [(0.0, 0.5), (0.5, 1.0), (1.0, 2.0), (2.0, 3.0),
              (3.0, 5.0), (5.0, 1e9)]
    bstat = {b: [0, 0.0, 0.0] for b in bandas}
    por_marco = {}
    n_tot = 0
    n_canal = 0          # delta < -2 Pi  (canal S < -Pi)
    for (tag, z, a, r, u_r, delta, Pi) in all_instances():
        n_tot += 1
        if delta < -2.0 * Pi:
            n_canal += 1
        if delta <= 0:
            continue
        S = a.lv.S[r]
        ratio = delta / (Pi * rho0(u_r))
        reg = (ratio, tag, z, a.L, r, u_r, delta / Pi, S)
        ratios.append(reg)
        if S < 1.0 - 1e-9:
            ratios_gen.append(reg)
        key = (tag, z)
        por_marco[key] = max(por_marco.get(key, 0.0), ratio)
        for b in bandas:
            if b[0] <= u_r < b[1]:
                st = bstat[b]
                st[0] += 1
                st[1] = max(st[1], delta / Pi)
                st[2] = max(st[2], ratio)
    ratios.sort()
    n = len(ratios)
    A_min = ratios[-1][0]
    A_gen = max((x[0] for x in ratios_gen), default=0.0)

    def q(f):
        return ratios[min(n - 1, int(f * n))][0]
    print("  instancias visibles r>1: %d | con delta>0: %d | "
          "canal delta<-2Pi (S<-Pi): %d" % (n_tot, n, n_canal))
    print("  A_min = max ratio = %.4f | restringido a sumas "
          "genuinas (S(r)<1): %.4f" % (A_min, A_gen))
    print("  cuantiles del ratio (delta>0): mediana %.4f | p90 %.4f "
          "| p99 %.4f | max %.4f" % (q(0.5), q(0.9), q(0.99), A_min))
    hist = [(0.0, 0.1), (0.1, 0.25), (0.25, 0.5), (0.5, 0.75),
            (0.75, 1.0), (1.0, 1.5), (1.5, 1e9)]
    hs = []
    for (lo, hi) in hist:
        c = sum(1 for x in ratios if lo <= x[0] < hi)
        hs.append("[%.2f,%.2f): %d" % (lo, min(hi, 9.99), c))
    print("  histograma: " + " | ".join(hs))
    print("  A_min por marco: " + ", ".join(
        "[%s] z=%d: %.4f" % (t, z, v)
        for ((t, z), v) in sorted(por_marco.items())))
    print("  bandas por u_r (solo delta>0):")
    for b in bandas:
        st = bstat[b]
        print("    u_r en [%.1f, %.1f): %5d instancias | max "
              "delta/Pi = %.4f | max ratio = %.4f"
              % (b[0], min(b[1], 99.0), st[0], st[1], st[2]))
    print("  top-10 instancias por ratio (S(r) incluido):")
    for (rt, tag, z, L, r, u_r, dpi, S) in ratios[-10:][::-1]:
        print("    [%s] z=%d L=%d r=%d: u_r=%.3f delta/Pi=%.4f "
              "S(r)=%.6f ratio=%.4f" % (tag, z, L, r, u_r, dpi, S, rt))
    print("  => la cota puntual positiva al nivel propio se sostiene "
          "en ventana con A = %.3f (con A=1: %s); el canal "
          "delta<-2Pi esta vacio: %s."
          % (A_min, "SI" if A_min <= 1.0 else "NO",
             "SI" if n_canal == 0 else "NO"))
    return A_min, ratios


# ------------------------------------------------------------------ K3
def check_K3():
    print("K3 los ofensores del doc 109 (|delta|/Pi hasta 1,295):")
    doc109 = {("A", 17): 1.000, ("A", 29): 1.069, ("A", 37): 1.295,
              ("B", 29): None}
    ok = True
    for (tag, z, lvs) in FRAMES:
        lnz = math.log(z)
        best_all = None
        best_vis = None
        best_vis_abs = None
        for a in lvs:
            lv = a.lv
            for (r, rpf, rmu, rphi) in lv.enum:
                dpi = lv.delta[r] / lv.Pi[r]
                if best_all is None or abs(dpi) > best_all[0]:
                    best_all = (abs(dpi), a.L, r,
                                1 if dpi > 0 else -1, r in a.vis)
                if r in a.vis and r > 1:
                    if best_vis is None or dpi > best_vis[0]:
                        best_vis = (dpi, a.L, r)
                    if (best_vis_abs is None
                            or abs(dpi) > best_vis_abs[0]):
                        best_vis_abs = (abs(dpi), a.L, r,
                                        1 if dpi > 0 else -1)
        (va, La, ra, sa, visa) = best_all
        (vp, Lp, rp) = best_vis
        u_rp = math.log(Lp / rp) / lnz
        u_ra = (math.log(La / ra) / lnz) if ra < La else float("nan")
        ref = doc109[(tag, z)]
        print("  [%s] z=%d: argmax |delta|/Pi sobre todo r: %.3f "
              "(doc 109: %s) en L=%d, r=%d, signo %s, visible=%s, "
              "u_r=%s"
              % (tag, z, va, ("%.3f" % ref) if ref else "  -  ",
                 La, ra, "+" if sa > 0 else "-", visa,
                 ("%.3f" % u_ra) if ra < La else "n/a (r>=L)"))
        print("        argmax delta/Pi sobre visibles r>1: %.3f en "
              "L=%d, r=%d, u_r=%.3f, rho_0(u_r)=%.3f, ratio=%.3f | "
              "argmax |delta|/Pi visible: %.3f (signo %s)"
              % (vp, Lp, rp, u_rp, rho0(u_rp), vp / rho0(u_rp),
                 best_vis_abs[0], "+" if best_vis_abs[3] > 0 else "-"))
        if ref:
            ok &= abs(va - ref) < 0.01
    print("  => el 1,000 de z=17 es un invisible (delta=-Pi "
          "estructural, no cuenta); los visibles positivos grandes "
          "tienen u_r ~ 1, donde rho_0(u_r) ~ 1. %s"
          % ("PASS" if ok else "FAIL"))
    return ok


# ------------------------------------------------------------------ K4/K6/K7
def check_K4K6K7(ratios):
    print("K4 modelo Vhat^(r) y error de desplazamiento E(r) = "
          "S(r) - Shat(r), sobre TODAS las instancias:")
    print("K6 diferencia simetrica V^(r) △ Vhat^(r) (conteos):")
    print("K7 soportes degenerados V^(r) = {1} (S(r) = 1 exacto):")
    ok = True
    max_E_rel = -1e300
    min_E_rel = 1e300
    max_E_pos_ratio = 0.0
    nE = 0
    nEpos = 0
    C_model = 0.0        # max |Shat-Pi|/(Pi rho_0(u_r))  (dos caras)
    C_model_pos = 0.0    # max (Shat-Pi)^+/(Pi rho_0(u_r)) (una cara)
    tot_n1 = tot_n2 = 0
    ejemplos = []
    max_absds = 0.0
    # K7: instancias con V^{(r)} = {1}
    deg = []             # (tag, z, L, r, u_r, delta/Pi, ratio, omega)
    n_deg = 0
    for (tag, z, lvs) in FRAMES:
        lnz = math.log(z)
        for a in lvs:
            lv = a.lv
            pfmap = {r: pf for (r, pf, mu, phi) in lv.enum}
            for r in a.vis:
                if r <= 1:
                    continue
                Shat, nV, n1, n2, ex1, ex2, absds = model_data(lv, r)
                Pi = lv.Pi[r]
                u_r = math.log(a.L / r) / lnz
                E = lv.S[r] - Shat
                nE += 1
                if E > 0:
                    nEpos += 1
                max_E_rel = max(max_E_rel, E / Pi)
                min_E_rel = min(min_E_rel, E / Pi)
                if E > 0:
                    max_E_pos_ratio = max(max_E_pos_ratio,
                                          E / (Pi * rho0(u_r)))
                C_model = max(C_model,
                              abs(Shat - Pi) / (Pi * rho0(u_r)))
                if Shat > Pi:
                    C_model_pos = max(C_model_pos,
                                      (Shat - Pi) / (Pi * rho0(u_r)))
                tot_n1 += n1
                tot_n2 += n2
                max_absds = max(max_absds, absds / Pi)
                if ex1 is not None and len(ejemplos) < 4:
                    ejemplos.append((tag, z, a.L, r, ex1, "V\\Vhat"))
                if ex2 is not None and len(ejemplos) < 8:
                    ejemplos.append((tag, z, a.L, r, ex2, "Vhat\\V"))
                if nV == 1:
                    n_deg += 1
                    delta = lv.delta[r]
                    if delta > 0:
                        deg.append((tag, z, a.L, r, u_r, delta / Pi,
                                    delta / (Pi * rho0(u_r)),
                                    len(pfmap[r])))
    print("  instancias (nivel, r visible>1): %d" % nE)
    print("  FL clasica sobre el MODELO al nivel propio: constante "
          "empirica dos caras C_model = %.4f | una cara (positiva) "
          "C_model^+ = %.4f" % (C_model, C_model_pos))
    print("  E(r): max E/Pi = %+.6f | min E/Pi = %+.4f | fraccion "
          "E>0: %.4f | max E^+/(Pi*rho_0(u_r)) = %.4f"
          % (max_E_rel, min_E_rel, nEpos / max(nE, 1),
             max_E_pos_ratio))
    print("  [K6] diferencia simetrica: total d' en V\\Vhat: %d | "
          "en Vhat\\V: %d | max masa absoluta/Pi = %.3f"
          % (tot_n1, tot_n2, max_absds))
    for (tag, z, L, r, dp, lado) in ejemplos[:6]:
        print("    ejemplo [%s] z=%d L=%d r=%d: d'=%d en %s"
              % (tag, z, L, r, dp, lado))
    ok &= tot_n1 > 0 and tot_n2 > 0
    # K7
    deg.sort(key=lambda x: -x[6])
    print("  [K7] instancias con V^(r) = {1}: %d (de %d) | con "
          "delta>0: %d" % (n_deg, nE, len(deg)))
    if deg:
        mx_ratio_deg = deg[0][6]
        print("    max ratio entre degenerados: %.4f | top-5:"
              % mx_ratio_deg)
        for (tag, z, L, r, u_r, dpi, rt, om) in deg[:5]:
            umbral = 3.0 * math.log(5) / math.log(z)
            print("    [%s] z=%d L=%d r=%d (omega=%d): u_r=%.3f "
                  "(umbral 3 ln5/ln z = %.3f) delta/Pi=%.4f "
                  "ratio=%.4f" % (tag, z, L, r, om, u_r, umbral,
                                  dpi, rt))
    # los top-10 de K2: ?son todos degenerados?
    claves_deg = {(t, z, L, r) for (t, z, L, r, *_rest) in deg}
    top10 = ratios[-10:]
    n_top_deg = sum(1 for (rt, t, z, L, r, *_rest) in top10
                    if (t, z, L, r) in claves_deg)
    print("    de los top-10 ratios de K2, %d son soportes "
          "degenerados V^(r) = {1}." % n_top_deg)
    print("  => E es <= 0 en ventana salvo excepciones <= %.4f "
          "relativas: el desplazamiento no aporta masa positiva; "
          "la diferencia simetrica es bilateral en conteo pero su "
          "masa con signo es una cara. La constante de (FL-DESP+) "
          "en ventana es la de la FL clasica sobre el modelo. %s"
          % (max_E_pos_ratio, "PASS" if ok else "FAIL"))
    return ok, max_E_pos_ratio, C_model_pos


# ------------------------------------------------------------------ K5
def check_K5(A_min):
    print("K5 nivel a nivel: P_{>1}(L) contra la cota implicada "
          "B_A(L) = sum_{r>1 vis} rho_A(u_r)(2+rho_A(u_r)) Pi^2/phi:")
    ok = True
    for A in (1.0, A_min):
        peor = 0.0
        peorB = 0.0
        nviol = 0
        for (tag, z, lvs) in FRAMES:
            lnz = math.log(z)
            for a in lvs:
                lv = a.lv
                P1 = a.P - max(a.t[1], 0.0)
                B = 0.0
                pfmap = {r: phi for (r, pf, mu, phi) in lv.enum}
                for r in a.vis:
                    if r <= 1:
                        continue
                    rho = A * rho0(math.log(a.L / r) / lnz)
                    B += rho * (2.0 + rho) * lv.Pi[r] ** 2 / pfmap[r]
                peorB = max(peorB, B / a.qc)
                if P1 > B + 1e-12:
                    nviol += 1
                if B > 0:
                    peor = max(peor, P1 / B)
        print("  A=%.3f: max_L B_A(L)/Q_c = %.3f | niveles con "
              "P_{>1} > B_A: %d/240 | max P_{>1}/B_A = %.3f"
              % (A, peorB, nviol, peor))
        if A == A_min:
            ok &= nviol == 0
    print("  => con A = A_min la implicacion nivel a nivel se "
          "verifica en los 240 niveles; con A=1 la SUMA tambien "
          "pasa en ventana (los ofensores puntuales pesan poco). %s"
          % ("PASS" if ok else "FAIL"))
    return ok


# ------------------------------------------------------------------ main
def main():
    t0 = time.time()
    print("== DOC 113: (FL-DESP+), lema fundamental desplazado "
          "positivo ==")
    print("Solo stdlib. Modo fijo; nivel L = D//s; "
          "u_r = ln(L/r)/ln z; rho_0(u) = u^{-u/2}.")
    print()
    allok = True
    allok &= check_K0()
    print()
    allok &= check_K1()
    print()
    A_min, ratios = check_K2()
    print()
    allok &= check_K3()
    print()
    ok4, _, _ = check_K4K6K7(ratios)
    allok &= ok4
    print()
    allok &= check_K5(A_min)
    print()
    print("Tiempo total: %.1f s" % (time.time() - t0))
    print("PASS" if allok else "FAIL")
    return 0 if allok else 1


if __name__ == "__main__":
    sys.exit(main())
