#!/usr/bin/env python3
# 112_masa_positiva_bloques.py  (Doc 112)
# Ataque a (H2'''), la masa positiva del exceso de Gram, por bloques
# de conductores. Definiciones de los docs 109/110/111 reutilizadas
# por import (Level111 importa el doc 110, que importa el doc 109:
# son literalmente las mismas funciones).
#
#   t(r) = delta(r)*(2*Pi(r)+delta(r))/phi(r),  E^(L) = sum_r t(r),
#   P^(L) = sum_{t(r)>0} t(r)  <- el objeto de (H2''').
#
# Descomposicion por bloques con R0 = z^2 (primos de P son < z, luego
# r >= z^2 => omega(r) >= 3 [T]; tambien se tabula R0 = z):
#   (i)   r = 1
#   (ii)  r visible (divide a algun d en V), 1 < r < R0
#   (iii) r visible, r >= R0
#   (iv)  r invisible: t(r) = -Pi(r)^2/phi(r) <= 0 [T, doc 111]:
#         no aporta a P.
# P^(L) = t(1)^+ + bloque(ii) + bloque(iii)  [T, particion].
#
# Chequeos:
#   K0  validacion cruzada contra el doc 111: max P/Q_c por marco
#       (1,248/1,073/0,862/0,159) y la particion de P en bloques.
#   K1  masas de bloque: maximo sobre los niveles muestreados de cada
#       marco, y valores en el nivel de maximo P (R0 = z^2; y R0 = z).
#   K2  bloque (i): Pi_full <= 1 <= M_full y Pi_full^2 <= Q_c [T];
#       identidad t(1) = rho(2+rho)*Pi_full^2 con rho = delta(1)/
#       Pi_full; rho(u) por bandas de u = ln L/ln z.
#   K3  bloque (iii), la ruta absoluta con productos tipo M:
#       Pi(r)/phi(r) = Pi_full*g(r) [T, doc 110]; |delta(r)| <=
#       D^-(r) <= P^+ [T]; sum_{r>=z^2} g(r) <= sum_{omega>=3} g(r)
#       = Ghat - 1 - G1 - e2 [T]; el mayorante Lambda3/Q_c resultante
#       y su crecimiento en z (calculos de primos hasta z=101).
#   K4  bloque (ii): censo de conductores visibles 1<r<z^2 por omega,
#       maximo t(r) individual, y el residuo (H2''''):
#       R^(L) = P^(L) - t(1)^+, maximo por marco.
#   K5  libro mayor global: P = t(1)^+ + R^(L) en todo nivel; la
#       cuota de r=1 en el maximo de P de cada marco.
# Uso: python3 112_masa_positiva_bloques.py | tee salida_doc112_local.txt
import importlib.util
import math
import sys
import time

spec = importlib.util.spec_from_file_location(
    "doc111", "111_exceso_gram_euler.py")
d111 = importlib.util.module_from_spec(spec)
spec.loader.exec_module(d111)

Level111 = d111.Level111
frame_levels = d111.frame_levels
primes_upto = d111.primes_upto
allowed_primes = d111.allowed_primes
g_of_pf = d111.g_of_pf


FRAMES = []
for z in (17, 29, 37):
    FRAMES.append(("A", z, z ** 5, int(z ** (10.0 / 3.0))))
FRAMES.append(("B", 29, 29 ** 6, 29 ** 3))

_cache = {}


def get_frames():
    if "frames" not in _cache:
        out = []
        for (tag, z, D, y) in FRAMES:
            lvs = [Level111(z, 3, (), lv.L)
                   for lv in frame_levels(z, 3, (), D, y)]
            out.append((tag, z, lvs))
        _cache["frames"] = out
    return _cache["frames"]


def blocks(a, R0):
    """Bloques de P^(L) para el nivel a (Level111) con umbral R0.
    Devuelve (b1, b2, b3, n2, n3, maxt2): masas positivas de los
    bloques (i),(ii),(iii), conteos de r en (ii),(iii) con t>0 y el
    maximo t(r) individual dentro de (ii)."""
    b1 = b2 = b3 = 0.0
    n2 = n3 = 0
    maxt2 = 0.0
    for (r, t) in a.t.items():
        if t <= 0:
            continue
        if r == 1:
            b1 += t
        elif r not in a.vis:
            pass  # no deberia ocurrir: invisible => t <= 0 [T]
        elif r < R0:
            b2 += t
            n2 += 1
            maxt2 = max(maxt2, t)
        else:
            b3 += t
            n3 += 1
    return b1, b2, b3, n2, n3, maxt2


# ------------------------------------------------------------------ K0
def check_K0(frames):
    print("K0 validacion cruzada contra el doc 111 y particion de P:")
    esperado = {("A", 17): 1.248, ("A", 29): 1.073,
                ("A", 37): 0.862, ("B", 29): 0.159}
    ok = True
    for (tag, z, lvs) in frames:
        mxP = max(lvs, key=lambda a: a.P)
        e = abs(mxP.P / mxP.qc - esperado[(tag, z)])
        err_part = 0.0
        err_inv = 0.0
        for a in lvs:
            b1, b2, b3, _, _, _ = blocks(a, z * z)
            err_part = max(err_part, abs(b1 + b2 + b3 - a.P))
            for (r, t) in a.t.items():
                if r not in a.vis:
                    err_inv = max(err_inv, t)  # debe ser <= 0
        print("  [%s] z=%d: max P/Q_c=%.3f (doc 111: %.3f, dif=%.4f) | "
              "err particion P=(i)+(ii)+(iii): %.2e | max t(r) sobre "
              "invisibles: %.2e (<= 0 [T])"
              % (tag, z, mxP.P / mxP.qc, esperado[(tag, z)], e,
                 err_part, err_inv))
        ok &= e < 0.002 and err_part < 1e-9 and err_inv <= 1e-12
    print("  => P^(L) se parte exactamente en (i)+(ii)+(iii); los "
          "invisibles nunca aportan. %s" % ("PASS" if ok else "FAIL"))
    return ok


# ------------------------------------------------------------------ K1
def check_K1(frames):
    print("K1 masas de bloque (R0 = z^2; entre parentesis R0 = z):")
    ok = True
    for (tag, z, lvs) in frames:
        R0 = z * z
        mxP = max(lvs, key=lambda a: a.P)
        b1, b2, b3, n2, n3, _ = blocks(mxP, R0)
        c1, c2, c3, _, _, _ = blocks(mxP, z)
        # maximos por bloque sobre todos los niveles
        m1 = max(blocks(a, R0)[0] / a.qc for a in lvs)
        m2 = max(blocks(a, R0)[1] / a.qc for a in lvs)
        m3 = max(blocks(a, R0)[2] / a.qc for a in lvs)
        print("  [%s] z=%d, nivel de max P (L=%d, P/Q_c=%.3f): "
              "(i)=%.3f (ii)=%.3f [%d r] (iii)=%.3f [%d r] | con "
              "R0=z: (ii)=%.3f (iii)=%.3f"
              % (tag, z, mxP.L, mxP.P / mxP.qc, b1 / mxP.qc,
                 b2 / mxP.qc, n2, b3 / mxP.qc, n3,
                 c2 / mxP.qc, c3 / mxP.qc))
        print("        maximos sobre los %d niveles (R0=z^2): "
              "(i)=%.3f (ii)=%.3f (iii)=%.3f (en unidades de Q_c)"
              % (len(lvs), m1, m2, m3))
        ok &= abs(b1 + b2 + b3 - mxP.P) < 1e-9
    print("  => PASS" if ok else "  => FAIL")
    return ok


# ------------------------------------------------------------------ K2
def check_K2(frames):
    print("K2 bloque (i) r=1: la desviacion del lema fundamental")
    ok = True
    for z in (17, 29, 37):
        lv = Level111(z, 3, (), 10 ** 9)  # L enorme: solo productos
        Pf = lv.lv.Pi_full
        Mf = lv.lv.M_full_sum
        qc = lv.qc
        print("  z=%d: Pi_full=%.6f <= 1 <= M_full=%.6f | "
              "Pi_full^2=%.6f <= Q_c=%.6f (cociente %.4f)"
              % (z, Pf, Mf, Pf * Pf, qc, Pf * Pf / qc))
        ok &= Pf <= 1.0 <= Mf and Pf * Pf <= qc + 1e-15
    # rho(u) = delta(1)/Pi_full por bandas de u; identidad de t(1)
    print("  rho_1(u) = delta(1)/Pi_full y t(1) por bandas de "
          "u = ln L/ln z (todos los marcos):")
    bandas = [(0.0, 2.0), (2.0, 3.0), (3.0, 4.0), (4.0, 99.0)]
    datos = {b: [] for b in bandas}
    err_id = 0.0
    for (tag, z, lvs) in frames:
        for a in lvs:
            lv = a.lv
            if a.L <= 1:
                continue
            u = math.log(a.L) / math.log(z)
            rho = lv.delta[1] / lv.Pi_full
            err_id = max(err_id,
                         abs(a.t[1] - rho * (2.0 + rho)
                             * lv.Pi_full ** 2))
            for b in bandas:
                if b[0] <= u < b[1]:
                    datos[b].append((rho, a.t[1] / a.qc, u))
    for b in bandas:
        if datos[b]:
            rhos = [abs(r) for (r, _, _) in datos[b]]
            print("    u en [%.1f, %.1f): %3d niveles | |rho_1| en "
                  "[%.3f, %.3f] | max t(1)^+/Q_c = %.3f"
                  % (b[0], min(b[1], 9.0), len(datos[b]), min(rhos),
                     max(rhos),
                     max(max(t, 0.0) for (_, t, _) in datos[b])))
    print("  identidad t(1) = rho(2+rho)*Pi_full^2: err max = %.2e"
          % err_id)
    ok &= err_id < 1e-9
    print("  => t(1) = rho(2+rho)*Pi_full^2 con Pi_full^2 <= Q_c: "
          "basta |rho| <= 1 (u >= u_0) para t(1) <= 3*Q_c. %s"
          % ("PASS" if ok else "FAIL"))
    return ok


# ------------------------------------------------------------------ K3
def prime_stats(z):
    """Cantidades de conteo sobre P (A=3) sin niveles: G1, e2, e3,
    Ghat = prod(1+1/(p-2)), P^+ = prod(1+1/(p-1)), M_full (producto),
    y sumas exactas de g(r) por omega y por tamano (solo z <= 37,
    via enumeracion; para z mayor, por sumas de potencias)."""
    allow = allowed_primes(z, 3, ())
    G1 = sum(1.0 / (p - 2) for p in allow)
    G2 = sum(1.0 / (p - 2) ** 2 for p in allow)
    G3 = sum(1.0 / (p - 2) ** 3 for p in allow)
    e2 = (G1 * G1 - G2) / 2.0
    e3 = (G1 ** 3 - 3.0 * G1 * G2 + 2.0 * G3) / 6.0
    Ghat = 1.0
    Pp = 1.0
    Mf = 1.0
    for p in allow:
        Ghat *= (1.0 + 1.0 / (p - 2))
        Pp *= (1.0 + 1.0 / (p - 1))
        Mf *= (1.0 + 1.0 / ((p - 1) * (p - 2)))
    sum_om3 = Ghat - 1.0 - G1 - e2   # sum_{omega>=3} g(r) exacto
    return {"allow": allow, "G1": G1, "e2": e2, "e3": e3,
            "Ghat": Ghat, "Pp": Pp, "Mf": Mf, "sum_om3": sum_om3}


def check_K3(frames):
    print("K3 bloque (iii): la ruta absoluta con productos tipo M")
    ok = True
    # verificaciones exactas a z <= 37 por enumeracion
    for z in (17, 29, 37):
        st = prime_stats(z)
        lv = Level111(z, 3, (), 10 ** 9)
        enum = lv.lv.enum
        sum_ge = sum(g_of_pf(pf) for (r, pf, mu, phi) in enum
                     if r >= z * z)
        sum_om3_enum = sum(g_of_pf(pf) for (r, pf, mu, phi) in enum
                           if len(pf) >= 3)
        maxr_om2 = max(r for (r, pf, mu, phi) in enum if len(pf) <= 2)
        err_g = 0.0
        for (r, pf, mu, phi) in enum:
            err_g = max(err_g, abs(lv.lv.Pi[r] / phi
                                   - lv.lv.Pi_full * g_of_pf(pf)))
        print("  z=%d: max r con omega<=2: %d < z^2=%d [T] | "
              "sum_{r>=z^2} g(r)=%.4f <= sum_{omega>=3} g(r)=%.4f "
              "(formula: %.4f) | e3=%.4f <= G1^3/6=%.4f | "
              "err Pi(r)/phi(r)=Pi_full*g(r): %.2e"
              % (z, maxr_om2, z * z, sum_ge, sum_om3_enum,
                 st["sum_om3"], st["e3"], st["G1"] ** 3 / 6.0, err_g))
        ok &= maxr_om2 < z * z and sum_ge <= sum_om3_enum + 1e-12
        ok &= abs(sum_om3_enum - st["sum_om3"]) < 1e-9
        ok &= err_g < 1e-9
        # mayorante absoluto de Lambda_3 (parte lineal de (iii))
        lam = 2.0 * st["Pp"] * sum_ge / st["Mf"]
        print("        mayorante |Lambda_3|/Q_c <= 2*P^+*sum_{r>=z^2}"
              "g(r)/M_full = %.3f  (P^+=%.3f, M_full=%.4f)"
              % (lam, st["Pp"], st["Mf"]))
    # crecimiento del mayorante con z (solo primos)
    print("  crecimiento del mayorante (cantidades de primos):")
    for z in (53, 71, 101):
        st = prime_stats(z)
        lam = 2.0 * st["Pp"] * st["sum_om3"] / st["Mf"]
        print("    z=%3d: G1=%.3f (log log z=%.3f) | sum_{omega>=3} "
              "g(r)=%.4f | P^+=%.3f | mayorante Lambda_3/Q_c <= %.3f"
              % (z, st["G1"], math.log(math.log(z)), st["sum_om3"],
                 st["Pp"], lam))
    # |delta(r)| <= D^-(r) <= P^+ en los niveles peores, y el maximo
    # medido de |delta|/Pi sobre visibles r >= z^2
    print("  cota |delta(r)| <= D^-(r) <= P^+ y |delta|/Pi medido:")
    for (tag, z, lvs) in frames:
        a = max(lvs, key=lambda t: t.P)
        lv = a.lv
        Pp = prime_stats(z)["Pp"]
        errD = 0.0
        maxratio = 0.0
        for (r, rpf, rmu, rphi) in lv.enum:
            Dm = 0.0
            for (dp, dpf, dmu, dphi) in lv.enum:
                if math.gcd(dp, r) == 1 and r * dp not in lv.Vset:
                    Dm += 1.0 / dphi
            errD = max(errD, abs(lv.delta[r]) - Dm)
            if r in a.vis and r >= z * z:
                maxratio = max(maxratio,
                               abs(lv.delta[r]) / lv.Pi[r])
        print("    [%s] z=%d L=%d: |delta| - D^-(r) <= %.2e en todo "
              "r (D^- <= P^+=%.3f) | max |delta|/Pi sobre visibles "
              "r>=z^2: %.3f [M]"
              % (tag, z, a.L, errD, Pp, maxratio))
        ok &= errD < 1e-9
    print("  => la ruta absoluta acota Lambda_3 con un factor que en "
          "ventana es < 1 pero NO es constante absoluta: P^+ ~ e^gamma"
          " ln z y sum_{omega>=3} g(r) >= e3 ~ G1^3/6 -> infinito "
          "(G1 ~ log log z). Pierde (ln z)(log log z)^3. %s"
          % ("PASS" if ok else "FAIL"))
    return ok


# ------------------------------------------------------------------ K4
def check_K4(frames):
    print("K4 bloque (ii): censo y el residuo (H2'''')")
    ok = True
    for (tag, z, lvs) in frames:
        for (nombre, a) in (("maxP", max(lvs, key=lambda t: t.P)),
                            ("maxL", lvs[-1])):
            censo = {1: 0, 2: 0, 3: 0}
            pfmap = {r: pf for (r, pf, mu, phi) in a.lv.enum}
            for r in a.vis:
                if r > 1 and r < z * z:
                    censo[min(len(pfmap[r]), 3)] += 1
            b1, b2, b3, n2, n3, maxt2 = blocks(a, z * z)
            print("  [%s] z=%d %s L=%d: visibles 1<r<z^2 por omega: "
                  "{1: %d, 2: %d, >=3: %d} | max t(r) en (ii): %.4f "
                  "(%.3f Q_c) | residuo R = P - t(1)^+ = %.3f Q_c"
                  % (tag, z, nombre, a.L, censo[1], censo[2],
                     censo[3], maxt2, maxt2 / a.qc,
                     (a.P - b1) / a.qc))
        mxR = max((a.P - blocks(a, z * z)[0]) / a.qc for a in lvs)
        print("        max_L R^(L)/Q_c = %.3f  (contra max P/Q_c = "
              "%.3f)" % (mxR, max(a.P / a.qc for a in lvs)))
        # sum_{r>1} Pi(r)^2/phi(r) = Q_c - Pi_full^2
        a = lvs[0]
        s = sum(a.lv.Pi[r] ** 2 / phi
                for (r, pf, mu, phi) in a.lv.enum if r > 1)
        e = abs(s - (a.qc - a.lv.Pi_full ** 2))
        ok &= e < 1e-9
    print("  => la familia (ii) crece con z (contiene todos los "
          "primos de P): no es finita uniformemente; su cota no sale "
          "de un analisis finito. El residuo R^(L) es el objeto de "
          "(H2''''). %s" % ("PASS" if ok else "FAIL"))
    return ok


# ------------------------------------------------------------------ K5
def check_K5(frames):
    print("K5 libro mayor: cuota de r=1 en el maximo de P")
    ok = True
    for (tag, z, lvs) in frames:
        a = max(lvs, key=lambda t: t.P)
        b1, b2, b3, _, _, _ = blocks(a, z * z)
        print("  [%s] z=%d L=%d: P/Q_c=%.3f = t(1)^+ %.3f (%.0f%%) "
              "+ (ii) %.3f + (iii) %.3f"
              % (tag, z, a.L, a.P / a.qc, b1 / a.qc,
                 100.0 * b1 / a.P if a.P > 0 else 0.0,
                 b2 / a.qc, b3 / a.qc))
        ok &= abs(b1 + b2 + b3 - a.P) < 1e-9
    print("  => PASS" if ok else "  => FAIL")
    return ok


# ------------------------------------------------------------------ main
def main():
    t0 = time.time()
    print("== DOC 112: masa positiva del exceso por bloques ==")
    print("Solo stdlib. Modo fijo; nivel L = D//s; R0 = z^2.")
    print()
    frames = get_frames()
    allok = True
    allok &= check_K0(frames)
    print()
    allok &= check_K1(frames)
    print()
    allok &= check_K2(frames)
    print()
    allok &= check_K3(frames)
    print()
    allok &= check_K4(frames)
    print()
    allok &= check_K5(frames)
    print()
    print("Tiempo total: %.1f s" % (time.time() - t0))
    print("PASS" if allok else "FAIL")
    return 0 if allok else 1


if __name__ == "__main__":
    sys.exit(main())
