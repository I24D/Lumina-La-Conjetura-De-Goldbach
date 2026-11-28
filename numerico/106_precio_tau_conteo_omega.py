#!/usr/bin/env python3
# 106_precio_tau_conteo_omega.py  (Doc 106)
# Trabajo fijado por el doc 105: (E105.1) afinar la pendiente de
# la pieza C de Bhat5 con soportes S19/S23; (A105.1) conteo exacto
# de los promedios de tau sobre Hcal y de Omega_d por reindexacion
# sobre V+(D/s) (estilo A100.1, doc 101); reevaluacion de L105;
# cociente B_Delta/T(N). Solo stdlib.
#
# Definiciones (docs 100-105, sin cambios):
#   H = soporte (cierre por divisores de m = s*d, A = 3 fijo);
#   H_N = {r en H : (r, N) = 1}; Hcal = {r en H_N : r > R0},
#   R0 = (ln N)^kappa, kappa = 1; piezas P (primo) / C (compuesto);
#   F_d = phi(d) E_d; M_r = sum_{d|r} mu(r/d) F_d;
#   S_r = M_r + Delta_r (correccion imprimitiva, docs 97-100);
#   B_M = sum_{r en Hcal} |S_r|^2/phi(r);
#   W_{de} = sum_{r en Hcal, d|r, e|r} mu(r/d)mu(r/e)/phi(r)
#          = mu(L/d)mu(L/e) U_L, L = mcm(d, e)   [T, doc 105]
#   Omega_d^X = sum_{r en Hcal_X, d|r} (tau(r) - 1)/phi(r),
#   U_d^X     = sum_{r en Hcal_X, d|r} 1/phi(r),
#   Dhat^X = sum_d phi(d)^2 E_d^2 Omega_d^X >= B_M^{X,(M)}  [T]
#   Bhat5^X = (sqrt(Dhat^X) + sqrt(B_Delta^X))^2 >= B_M^{X,(S)}
#   precio tau^C = Dhat^C/Diag^C (Diag con pesos U_d^C).
#   T(N) (doc 100 seccion 5.3, heredado por L105 del doc 105):
#   el puente exige sqrt(G_grande * B_M) << N/(A ln^2 N) con
#   G_grande = 0.0225 y A = N^c, c = 0.039884; luego
#   T(N) := N^{2-2c}/(0.0225 * ln^4 N), 2-2c = 1.920232.
#
# Conteo A105.1 [T, reindexacion]: H se parte sin solapamiento
# (doc 101, lema 2) en la rama lisa V0 := V+(D/s_min) y las ramas
# s*V+(D/s) (s sprime). Para r = s*h (s no divide a h): tau(r) =
# 2 tau(h), 1/phi(r) = 1/((s-1) phi(h)). Con
#   Phi(V)  = sum_{h en V} 1/phi(h),
#   Thet(V) = sum_{h en V} (tau(h)-1)/phi(h),
#   Thet2(V) = sum_{h en V} (2 tau(h)-1)/phi(h) = 2 Thet + Phi,
#   Thet2(V; d) = sum_{h en V, d|h} (2 tau(h)-1)/phi(h),
#   Phi(V; d)   = sum_{h en V, d|h} 1/phi(h),
# y Thet_C(V0), Phi_C(V0) las sumas lisas filtradas a Hcal_C:
#   sum_{r en Hcal_C} (tau(r)-1)/phi(r)
#       = Thet_C(V0) + sum_{s no | N} (Thet2(Vs[s]) - 1)/(s-1),
#   sum_{r en Hcal_C} 1/phi(r)
#       = Phi_C(V0) + sum_{s no | N} (Phi(Vs[s]) - 1)/(s-1),
# (el -1 retira el termino h = 1, que da r = s primo, pieza P;
#  en la rama s, r = s*h con h >= 2 ya cumple r > R0 y la parte
#  lisa de (r, N) = 1 es vacua en los dos N del modelo: los
#  primos de h son < z <= 37 < 79 <= factores impares de N; se
#  verifica por asercion en K3.)
#   tau_bar^C := cociente de ambas sumas (media de tau(r)-1 sobre
#   Hcal_C con peso 1/phi(r), SIN psi: prediccion de conteo).
#   Omega_d^C por casos (d en D, el conjunto de divisores >= 2
#   de Hcal):
#   - d liso (todo primo < z):
#       Omega_d = Om0[d] + sum_{s no|N} Thet2(Vs[s]; d)/(s-1),
#       Om0[d] = sum_{h en V0 cap Hcal_C, d|h} (tau(h)-1)/phi(h);
#   - d = s0*d0 (s0 sprime, d0 liso), s0 no | N:
#       Omega_d = Thet2(Vs[s0]; d0)/(s0-1) si d0 >= 2,
#       Omega_d = (Thet2(Vs[s0]) - 1)/(s0-1) si d0 = 1;
#     y 0 si s0 | N. Analogas para U_d^C con Phi en vez de Thet.
#
# Chequeos:
#   K0  Maquinaria doc 105 (R2 |H|=880; build_H z=11 == Model97)
#       y build_H106 (pf_trial) == build_H en z = 11, 13, 17.
#   K1  Reproduccion del doc 105 en los seis puntos: B_M^(S)
#       (rel < 2e-3), Bhat5/B_M (rel < 2e-3), Diag^C/B_M^{C,(S)}
#       (rel < 3e-2).
#   K2  Formula de Omega_d^C y U_d^C (conteo) contra enumeracion
#       directa desde H, para todo d en D, en R1, R3 y S11..S23
#       (y S29), en ambos N: conjuntos de claves iguales y error
#       relativo maximo < 1e-9.
#   K3  Formula de tau_bar^C contra enumeracion directa (rel
#       < 1e-12) en todos los soportes de conteo; identidad
#       sum_d Omega_d^C = sum_{r en Hcal_C} (tau(r)-1)^2/phi(r);
#       asercion de vacuidad de (h, N) > 1 en las ramas.
#   K4  Mayorante valido en los puntos nuevos: Dhat^C >=
#       B_M^{C,(M)}, Bhat5 >= B_M^(S), pieza P exacta.
#   K5  E105.1 (reporte): pendientes en seis y en diez puntos;
#       certificacion intervalar de la pieza C bajo 0.1.
#   K6  B_Delta/T(N) y cotas elementales (reporte).
# Uso: python3 106_precio_tau_conteo_omega.py
import importlib.util
import math
import sys
import time

spec = importlib.util.spec_from_file_location(
    "s105", "105_forma_diagonal_segundos_momentos.py")
s105 = importlib.util.module_from_spec(spec)
spec.loader.exec_module(s105)
s102 = s105.s104.s103.s102

Model97 = s105.Model97
build_H = s105.build_H
SegPsi = s105.SegPsi
M_r_E = s105.M_r_E
ols = s105.ols
phi_of = s105.phi_of
measure105 = s105.measure105
MU, SPF, SMALL_PRIMES = s105.MU, s105.SPF, s105.SMALL_PRIMES
pf_trial = s102.pf_trial
divisors_from_pf = s102.divisors_from_pf
in_Vplus = s102.in_Vplus
primes_upto = s102.primes_upto

C_CHEN = 0.039884
G_GRANDE = 0.0225


def T_of_N(N):
    return N ** (2 - 2 * C_CHEN) / (G_GRANDE * math.log(N) ** 4)


def build_H106(A, z, y, D):
    """build_H del doc 102 con pf_trial en vez de pf_desc(spf):
    identica construccion (pares (s,d) -> m -> cierre por
    divisores), sin tope de spf. Se valida contra build_H en K0."""
    allpr = primes_upto(max(y, 100))
    sprimes = [p for p in allpr
               if p >= z and p < y and math.gcd(p, A) == 1]
    allow = [p for p in allpr if 2 < p < z and math.gcd(p, A) == 1]
    divs = [1]
    for p in allow:
        divs += [d * p for d in list(divs)]
    Vs = {}
    Mset = set()
    for s in sprimes:
        Ds = D // s
        Vs_s = []
        for d in divs:
            if d >= Ds:
                continue
            if d == 1 or in_Vplus(d, Ds, pf_trial(d, SMALL_PRIMES)):
                Vs_s.append(d)
                Mset.add(s * d)
        Vs[s] = sorted(Vs_s)
    H = {1}
    for m in Mset:
        for h in divisors_from_pf(pf_trial(m, SMALL_PRIMES)):
            H.add(h)
    return sorted(H), sprimes, Vs


def tau_phi(n):
    """(tau(n), phi(n)) para n libre de cuadrados (todo H)."""
    if n == 1:
        return 1, 1
    pf = pf_trial(n, SMALL_PRIMES)
    ph = n
    for p in pf:
        ph = ph // p * (p - 1)
    return 2 ** len(pf), ph


def is_prime106(n):
    if n < 2:
        return False
    pf = pf_trial(n, SMALL_PRIMES)
    return len(pf) == 1 and pf[0] == n


# ------------------------------------------------------- conteo
class Counting:
    """Agregados de nivel y formulas A105.1 para un (soporte, N).
    No usa psi: solo la estructura (H, sprimes, Vs)."""

    def __init__(self, H, sprimes, Vs, z, N, kappa=1.0):
        self.z = z
        self.N = N
        self.R0 = math.log(N) ** kappa
        self.sprimes = sprimes
        self.Vs = Vs
        self.smin = sprimes[0]
        self.V0 = Vs[self.smin]
        self.removed_s = {s for s in sprimes if N % s == 0}
        self.ok_s = [s for s in sprimes if s not in self.removed_s]
        HN = [r for r in H if math.gcd(r, N) == 1]
        self.HN = HN
        self.Hcal = [r for r in HN if r > self.R0]
        self.Hcal_C = [r for r in self.Hcal if not is_prime106(r)]
        # agregados de nivel (totales)
        self.Phi = {}
        self.Thet = {}
        self.Thet2 = {}
        for s in sprimes:
            ph = 0.0
            th = 0.0
            for h in Vs[s]:
                tau, phi = tau_phi(h)
                ph += 1.0 / phi
                th += (tau - 1) / phi
            self.Phi[s] = ph
            self.Thet[s] = th
            self.Thet2[s] = 2 * th + ph
        # agregados por multiplos, por nivel (ramas)
        self.Om_branch = {}
        self.U_branch = {}
        for s in sprimes:
            om = {}
            uu = {}
            for h in Vs[s]:
                tau, phi = tau_phi(h)
                w_om = (2 * tau - 1) / phi
                w_u = 1.0 / phi
                for d in divisors_from_pf(pf_trial(h, SMALL_PRIMES)):
                    if d < 2:
                        continue
                    om[d] = om.get(d, 0.0) + w_om
                    uu[d] = uu.get(d, 0.0) + w_u
            self.Om_branch[s] = om
            self.U_branch[s] = uu
        # rama lisa filtrada a Hcal_C (sumas por multiplos y totales)
        self.Om0 = {}
        self.U0 = {}
        self.Thet_C_V0 = 0.0
        self.Phi_C_V0 = 0.0
        for h in self.V0:
            if h <= self.R0 or math.gcd(h, N) != 1 \
                    or is_prime106(h) or h == 1:
                continue
            tau, phi = tau_phi(h)
            self.Thet_C_V0 += (tau - 1) / phi
            self.Phi_C_V0 += 1.0 / phi
            for d in divisors_from_pf(pf_trial(h, SMALL_PRIMES)):
                if d < 2:
                    continue
                self.Om0[d] = self.Om0.get(d, 0.0) + (tau - 1) / phi
                self.U0[d] = self.U0.get(d, 0.0) + 1.0 / phi

    # --- sumas totales sobre Hcal_C (formula) ---
    def totals_formula(self):
        num = self.Thet_C_V0
        den = self.Phi_C_V0
        for s in self.ok_s:
            num += (self.Thet2[s] - 1.0) / (s - 1)
            den += (self.Phi[s] - 1.0) / (s - 1)
        return num, den

    # --- Omega_d / U_d por la formula de casos ---
    def omega_formula(self, d):
        pf = pf_trial(d, SMALL_PRIMES)
        big = [p for p in pf if p >= self.z]
        if len(big) >= 2:
            return 0.0, 0.0
        if len(big) == 1:
            s0 = big[0]
            d0 = d // s0
            if s0 in self.removed_s:
                return 0.0, 0.0
            if d0 == 1:
                om = (self.Thet2[s0] - 1.0) / (s0 - 1)
                uu = (self.Phi[s0] - 1.0) / (s0 - 1)
            else:
                om = self.Om_branch[s0].get(d0, 0.0) / (s0 - 1)
                uu = self.U_branch[s0].get(d0, 0.0) / (s0 - 1)
            return om, uu
        om = self.Om0.get(d, 0.0)
        uu = self.U0.get(d, 0.0)
        for s in self.ok_s:
            om += self.Om_branch[s].get(d, 0.0) / (s - 1)
            uu += self.U_branch[s].get(d, 0.0) / (s - 1)
        return om, uu

    # --- enumeracion directa (verdad de referencia) ---
    def direct(self):
        if getattr(self, "_dir", None) is not None:
            return self._dir
        Om = {}
        Uu = {}
        num = 0.0
        den = 0.0
        num2 = 0.0
        for r in self.Hcal_C:
            tau, phi = tau_phi(r)
            num += (tau - 1) / phi
            den += 1.0 / phi
            num2 += (tau - 1) ** 2 / phi
            for d in divisors_from_pf(pf_trial(r, SMALL_PRIMES)):
                if d < 2:
                    continue
                Om[d] = Om.get(d, 0.0) + (tau - 1) / phi
                Uu[d] = Uu.get(d, 0.0) + 1.0 / phi
        self._dir = (Om, Uu, num, den, num2)
        return self._dir

    def omega_stats(self):
        """rho_C = max Omega_d/U_d; W* = max Omega_d * phi(d);
        sumas elementales para B_Delta (sobre Hcal entero)."""
        Om, Uu, num, den, num2 = self.direct()
        rho = 0.0
        wstar = 0.0
        for d, u in Uu.items():
            if u > 0:
                rho = max(rho, Om[d] / u)
            _, phid = tau_phi(d)
            wstar = max(wstar, Om[d] * phid)
        s_w1 = 0.0
        s_w2 = 0.0
        for r in self.Hcal:
            pf = pf_trial(r, SMALL_PRIMES)
            om_r = len(pf)
            _, phi = tau_phi(r)
            s_w1 += phi * om_r ** 2
            s_w2 += om_r ** 2 * r * r / phi
        return {"rho": rho, "wstar": wstar, "S_w1": s_w1,
                "S_w2": s_w2, "num": num, "den": den,
                "num2": num2}


def fmt_slope(name, t):
    return ("     %-34s pendiente=%+.4f +- %.4f (SE), IC95%% "
            "+-%.4f, LOO [%.4f, %.4f]"
            % (name, t[0], t[2], t[3], t[4], t[5]))


def main():
    allok = True
    t00 = time.time()
    print("== DOC 106: precio tau con S19/S23 + conteo A105.1 ==")
    print("Solo stdlib. MEDICION, no prueba. Goldbach sigue "
          "abierta.")
    print()
    # ---------- soportes ----------
    m1 = Model97(3, 11, 60, 20000, "fixed")
    m2 = Model97(3, 11, 3000, 200000, "fixed")
    m3 = Model97(3, 3, 349, 200000, "fixed")
    supports = {"R1": (m1.H, 11, None, None),
                "R3": (m3.H, 3, None, None)}
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
    # K0: maquinaria + build_H106 == build_H en z = 11, 13, 17
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
          "== build_H en z=11,13,17): %s" % ("PASS" if k0 else "FAIL"))
    print()
    # ---------- (a) E105.1: medida completa, diez puntos ----------
    plan = [(999982, ["S11", "S13", "S17", "S19", "S23"], "gen"),
            (1048576, ["S11", "S13", "S17", "S19", "S23"], "2^e")]
    pts = {}
    print("Diez puntos (medida completa, doc 105):")
    for (N, tags, fam) in plan:
        for tag in tags:
            H, z, sprimes, Vs = supports[tag]
            t0 = time.time()
            m = measure105(H, N)
            key = "%s/%s" % (fam, tag)
            pts[key] = m
            print("     [%s] N=%d |H_N|=%d |Hcal|=%d "
                  "B_M^(S)=%.6e\n"
                  "        Bhat5/B_M=%.3f (C x%.3f)  "
                  "Dhat^C/B_M^{C,(M)}=%.3f  Diag^C=%.3f  "
                  "tau^C=%.3f  share_C=%.3f  B_Delta=%.4e  "
                  "B_Delta/T=%.3e (%.1fs)"
                  % (key, N, m["nHN"], m["nHcal"], m["BMs"],
                     m["Bhat"] / m["BMs"], m["Bhat_C"] / m["BMs_C"],
                     m["Dhat_C"] / m["BMm_C"],
                     m["Diag_C"] / m["BMm_C"],
                     m["Dhat_C"] / m["Diag_C"],
                     m["BMs_C"] / m["BMs"],
                     m["BD_P"] + m["BD_C"],
                     (m["BD_P"] + m["BD_C"]) / T_of_N(N),
                     time.time() - t0))
    print("Calibracion R1/R3 (N=999982):")
    for tag in ("R1", "R3"):
        H, z, _, _ = supports[tag]
        m = measure105(H, 999982)
        pts["gen/" + tag] = m
        print("     [gen/%s] |H_N|=%d B_M^(S)=%.6e Bhat5/B_M=%.3f "
              "B_Delta=%.4e B_Delta/T=%.3e"
              % (tag, m["nHN"], m["BMs"], m["Bhat"] / m["BMs"],
                 m["BD_P"] + m["BD_C"],
                 (m["BD_P"] + m["BD_C"]) / T_of_N(999982)))
    print()
    # K1: reproduccion del doc 105 en los seis puntos originales
    doc105_BM = {"gen/S11": 3.3658e9, "gen/S13": 7.6433e9,
                 "gen/S17": 2.7875e10, "2^e/S11": 3.7105e9,
                 "2^e/S13": 8.3332e9, "2^e/S17": 2.9761e10}
    doc105_q = {"gen/S11": 3.190, "gen/S13": 4.115,
                "gen/S17": 5.499, "2^e/S11": 3.600,
                "2^e/S13": 4.470, "2^e/S17": 5.619}
    doc105_DiagC = {"gen/S11": 1.495, "gen/S13": 1.488,
                    "gen/S17": 1.533, "2^e/S11": 1.674,
                    "2^e/S13": 1.616, "2^e/S17": 1.576}
    k1 = True
    for tag, ref in doc105_BM.items():
        if abs(pts[tag]["BMs"] - ref) / ref > 2e-3:
            k1 = False
            print("     K1 DISCREPA B_M^(S) [%s]" % tag)
    for tag, ref in doc105_q.items():
        got = pts[tag]["Bhat"] / pts[tag]["BMs"]
        if abs(got - ref) / ref > 2e-3:
            k1 = False
            print("     K1 DISCREPA Bhat5/B_M [%s]: %.3f vs %.3f"
                  % (tag, got, ref))
    for tag, ref in doc105_DiagC.items():
        got = pts[tag]["Diag_C"] / pts[tag]["BMs_C"]
        if abs(got - ref) / ref > 3e-2:
            k1 = False
            print("     K1 DISCREPA Diag^C [%s]" % tag)
    allok &= k1
    print("K1 reproduccion doc 105 (seis puntos): %s"
          % ("PASS" if k1 else "FAIL (ver arriba)"))
    # K4: mayorante en los puntos nuevos
    k4 = True
    for key, m in pts.items():
        if m["Dhat_C"] < m["BMm_C"] * (1 - 1e-12) \
                or m["Bhat"] < m["BMs"] * (1 - 1e-12) \
                or abs(m["Dhat_P"] - m["BMm_P"]) / m["BMm_P"] \
                > 1e-9:
            k4 = False
            print("     K4 FALLO [%s]" % key)
    allok &= k4
    print("K4 mayorante Dhat^C >= B_M^{C,(M)}, Bhat5 >= B_M^(S), "
          "pieza P exacta, en los 12 puntos: %s"
          % ("PASS" if k4 else "FAIL"))
    print()
    # ---------- pendientes: seis vs diez ----------
    def slopes(keys):
        xs = [math.log(pts[k]["nHN"]) for k in keys]
        out = {}
        series = {
            "agregado Bhat5/B_M^(S)":
                [math.log(pts[k]["Bhat"] / pts[k]["BMs"])
                 for k in keys],
            "pieza C, forma S (Bhat5^C)":
                [math.log(pts[k]["Bhat_C"] / pts[k]["BMs_C"])
                 for k in keys],
            "pieza C, forma M (Dhat^C)":
                [math.log(pts[k]["Dhat_C"] / pts[k]["BMm_C"])
                 for k in keys],
            "precio tau^C (Dhat^C/Diag^C)":
                [math.log(pts[k]["Dhat_C"] / pts[k]["Diag_C"])
                 for k in keys],
            "factor Diag^C/B_M^{C,(M)}":
                [math.log(pts[k]["Diag_C"] / pts[k]["BMm_C"])
                 for k in keys],
        }
        for name, ys in series.items():
            out[name] = ols(xs, ys)
        return out
    six_keys = ["gen/S11", "gen/S13", "gen/S17",
                "2^e/S11", "2^e/S13", "2^e/S17"]
    ten_keys = ["gen/S11", "gen/S13", "gen/S17", "gen/S19",
                "gen/S23", "2^e/S11", "2^e/S13", "2^e/S17",
                "2^e/S19", "2^e/S23"]
    s6 = slopes(six_keys)
    s10 = slopes(ten_keys)
    print("Pendientes log-log contra |H_N| — SEIS puntos (doc 105):")
    for name, t in s6.items():
        print(fmt_slope(name, t))
    print("Pendientes log-log contra |H_N| — DIEZ puntos (E105.1):")
    for name, t in s10.items():
        print(fmt_slope(name, t))
    for forma in ("pieza C, forma S (Bhat5^C)",
                  "pieza C, forma M (Dhat^C)"):
        t = s10[forma]
        cert = t[5] < 0.1 and (t[0] + t[3]) < 0.1
        print("     CERTIFICACION pieza C (%s, diez puntos): "
              "LOO_max=%.4f, b+IC95=%.4f -> %s"
              % (forma, t[5], t[0] + t[3],
                 "CERTIFICADA < 0.1" if cert else "NO certificada"))
    print()
    # ---------- (b) A105.1: conteo ----------
    print("Conteo A105.1 (sin psi):")
    cnt = {}
    for tag in ("R1", "R3", "S11", "S13", "S17", "S19", "S23",
                "S29", "S31", "S37"):
        H, z, sprimes, Vs = supports[tag]
        if sprimes is None:
            # R1/R3: reconstruir (s, Vs) via Model97 equivalente
            mm = m1 if tag == "R1" else m3
            sprimes = mm.sprimes
            Vs = mm.Vs
            z = mm.z
        for N in (999982, 1048576):
            if tag in ("S31", "S37") and N == 1048576:
                continue
            t0 = time.time()
            c = Counting(H, sprimes, Vs, z, N)
            cnt[(tag, N)] = c
            num_f, den_f = c.totals_formula()
            Om_d, U_d, num_d, den_d, num2_d = c.direct()
            if den_d > 0:
                tau_f = num_f / den_f
                tau_d = num_d / den_d
                print("     [%s N=%d] |Hcal_C|=%d tau_bar^C "
                      "formula=%.6f enumeracion=%.6f (rel %.2e)  "
                      "A_bar=%.4f (%.1fs)"
                      % (tag, N, len(c.Hcal_C), tau_f, tau_d,
                         abs(tau_f - tau_d) / tau_d,
                         num2_d / num_d, time.time() - t0))
            else:
                print("     [%s N=%d] |Hcal_C|=0 (soporte de solo "
                      "primos; tau_bar^C no definido) (%.1fs)"
                      % (tag, N, time.time() - t0))
    # K3: tau_bar formula vs enumeracion + identidad sum_d Omega
    k3 = True
    for (tag, N), c in cnt.items():
        num_f, den_f = c.totals_formula()
        Om_d, U_d, num_d, den_d, num2_d = c.direct()
        if den_d > 0:
            if abs(num_f / den_f - num_d / den_d) / (num_d / den_d) \
                    > 1e-12:
                k3 = False
                print("     K3 DISCREPA tau_bar [%s N=%d]" % (tag, N))
            s_om = sum(Om_d.values())
            if abs(s_om - num2_d) / num2_d > 1e-9:
                k3 = False
                print("     K3 DISCREPA sum_d Omega [%s N=%d]"
                      % (tag, N))
        elif num_f != 0.0 or num_d != 0.0:
            k3 = False
            print("     K3 DISCREPA vacio [%s N=%d]" % (tag, N))
        # vacuidad: (h, N) = 1 en las ramas supervivientes
        for s in c.ok_s:
            for h in c.Vs[s]:
                if math.gcd(h, N) != 1:
                    k3 = False
                    print("     K3 NO VACUO [%s N=%d s=%d h=%d]"
                          % (tag, N, s, h))
    allok &= k3
    print("K3 tau_bar^C formula == enumeracion (rel < 1e-12); "
          "sum_d Omega_d = sum_r (tau-1)^2/phi; ramas vacuas: %s"
          % ("PASS" if k3 else "FAIL (ver arriba)"))
    # K2: Omega_d / U_d formula contra enumeracion, todo d en D
    k2 = True
    for (tag, N), c in cnt.items():
        if tag in ("S31", "S37"):
            continue
        Om_d, U_d, _, _, _ = c.direct()
        keys_f = set()
        worst = 0.0
        for d in set(Om_d) | set(U_d):
            om_f, u_f = c.omega_formula(d)
            keys_f.add(d)
            for (got, ref) in ((om_f, Om_d.get(d, 0.0)),
                               (u_f, U_d.get(d, 0.0))):
                if ref > 0:
                    worst = max(worst, abs(got - ref) / ref)
                elif got != 0.0:
                    worst = float("inf")
        if keys_f != set(Om_d):
            k2 = False
            print("     K2 CLAVES [%s N=%d]" % (tag, N))
        if worst > 1e-9:
            k2 = False
            print("     K2 DISCREPA [%s N=%d] peor rel %.2e"
                  % (tag, N, worst))
        else:
            print("     K2 [%s N=%d] Omega/U formula == enum., "
                  "|D|=%d, peor rel %.2e"
                  % (tag, N, len(Om_d), worst))
    allok &= k2
    print("K2 formula de Omega_d^C y U_d^C (todo d, R1/R3/S11-"
          "S29, ambos N): %s" % ("PASS" if k2 else "FAIL"))
    print()
    # ajustes del precio de conteo y comparacion con el medido
    print("Ajustes de tau_bar^C (familia gen, N=999982):")
    zs = [11, 13, 17, 19, 23, 29, 31, 37]
    tb = []
    for z in zs:
        c = cnt[("S%d" % z, 999982)]
        num, den = c.totals_formula()
        tb.append(num / den)
    g1 = ols([math.log(math.log(z)) for z in zs],
             [math.log(v) for v in tb])
    g2 = ols([math.log(z) for z in zs], [math.log(v) for v in tb])
    print("     tau_bar^C(z): " + "  ".join(
        "z=%d: %.4f" % (z, v) for (z, v) in zip(zs, tb)))
    print("     ln tau_bar ~ gamma * ln ln z: gamma=%.3f "
          "+- %.3f, LOO [%.3f, %.3f]"
          % (g1[0], g1[3], g1[4], g1[5]))
    print("     ln tau_bar ~ beta * ln z:     beta=%.3f "
          "+- %.3f, LOO [%.3f, %.3f]"
          % (g2[0], g2[3], g2[4], g2[5]))
    xsH = [math.log(len(cnt[("S%d" % z, 999982)].HN)) for z in zs]
    g3 = ols(xsH, [math.log(v) for v in tb])
    print("     pendiente local de ln tau_bar contra ln|H_N|: "
          "%.4f (LOO [%.4f, %.4f])" % (g3[0], g3[4], g3[5]))
    print("Comparacion precio medido (con psi) vs conteo:")
    for key in ten_keys:
        tag = key.split("/")[1]
        N = 999982 if key.startswith("gen") else 1048576
        c = cnt[(tag, N)]
        num, den = c.totals_formula()
        tau_c = num / den
        st = c.omega_stats()
        a_bar = st["num2"] / st["num"]
        price = pts[key]["Dhat_C"] / pts[key]["Diag_C"]
        print("     [%s] precio=%.3f  tau_bar^C=%.3f "
              "(precio/tau_bar=%.3f)  A_bar=%.3f "
              "(precio/A_bar=%.3f)  rho_C=%.3f  W*=%.2f"
              % (key, price, tau_c, price / tau_c, a_bar,
                 price / a_bar, st["rho"], st["wstar"]))
    xs10 = [math.log(pts[k]["nHN"]) for k in ten_keys]
    ys_r = []
    for key in ten_keys:
        tag = key.split("/")[1]
        N = 999982 if key.startswith("gen") else 1048576
        num, den = cnt[(tag, N)].totals_formula()
        price = pts[key]["Dhat_C"] / pts[key]["Diag_C"]
        ys_r.append(math.log(price / (num / den)))
    tr = ols(xs10, ys_r)
    print("     pendiente de ln(precio/tau_bar^C) contra "
          "ln|H_N| (diez puntos): %.4f (LOO [%.4f, %.4f])"
          % (tr[0], tr[4], tr[5]))
    print()
    # ---------- (d) B_Delta / T(N) y cotas elementales ----------
    print("B_Delta frente a T(N) = N^1.920232/(0.0225 ln^4 N):")
    for key in ten_keys + ["gen/R1", "gen/R3"]:
        N = 999982 if key.startswith("gen") else 1048576
        tag = key.split("/")[1]
        m = pts[key]
        bd = m["BD_P"] + m["BD_C"]
        c = cnt[(tag, N)]
        st = c.omega_stats()
        lnN = math.log(N)
        bnd1 = lnN ** 2 * st["S_w1"]
        bnd2 = lnN ** 2 * st["S_w2"]
        print("     [%s] B_Delta=%.4e  B_Delta/T=%.3e  "
              "cota1/T=%.3e  cota2/T=%.3e  (T=%.4e)"
              % (key, bd, bd / T_of_N(N), bnd1 / T_of_N(N),
                 bnd2 / T_of_N(N), T_of_N(N)))
    print()
    tC = s10["pieza C, forma S (Bhat5^C)"]
    tM = s10["pieza C, forma M (Dhat^C)"]
    certS = tC[5] < 0.1 and (tC[0] + tC[3]) < 0.1
    certM = tM[5] < 0.1 and (tM[0] + tM[3]) < 0.1
    print("RESUMEN: E105.1 pieza C diez puntos: forma S "
          "pendiente=%.4f LOO [%.4f, %.4f] certificada=%s; "
          "forma M pendiente=%.4f LOO [%.4f, %.4f] "
          "certificada=%s; gamma(tau_bar)=%.3f"
          % (tC[0], tC[4], tC[5], certS, tM[0], tM[4], tM[5],
             certM, g1[0]))
    print("PASS" if allok else "FAIL (ver lineas marcadas)")
    print("Tiempo total: %.1f s" % (time.time() - t00))
    return 0 if allok else 1


if __name__ == "__main__":
    sys.exit(main())
