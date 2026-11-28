#!/usr/bin/env python3
# 103_criba_restringida_segundo_factor.py  (Doc 103)
# Prueba propuesta en el doc 102: criba grande restringida a H_A
# sobre el segundo factor B_M con S_r corregido (Delta_r conservado).
# Solo stdlib.
#
# Candidata Bhat1 (criba grande de Montgomery-Vaughan restringida):
#   (1/phi(r)) |S_r|^2 <= sum*_{chi mod r} |psi(N,chi)|^2
#       (Cauchy por conductor, constante 1; S_r entra entero, con
#        su Delta_r, via su definicion como suma de caracteres)
#   sum_{r en Hcal} sum*_{chi} |psi(N,chi)|^2
#       <= (N + delta_H^{-1}) * Lambda_2(N) =: Bhat1,
#   Hcal = {r en H_N : r > R0}, R0 = (ln N)^kappa, kappa = 1,
#   delta_H = espaciamiento minimo de {a/r mod 1 : r en Hcal,
#   (a,r)=1}, Lambda_2(N) = sum_{n<=N} Lambda(n)^2 (medido exacto).
#
# Sandwich riguroso para delta_H [T]: H es impar (m = s*d con s, d
# impares), luego no hay enteros consecutivos en Hcal y el exhibido
# 1/r, 1/(r+2) da delta_H <= 2/(r2 (r2+2)), con r2 = max{r en Hcal :
# r+2 en Hcal}. Y fracciones reducidas distintas distan >= 1/Q_H^2,
# Q_H = max Hcal. Luego
#   Bhat_lo := (N + r2(r2+2)/2) * Lambda_2  <=  Bhat1
#            <= (N + Q_H^2) * Lambda_2 =: Bhat_hi.
# Si Bhat_lo/B_M > 10, C1 falla para la candidata verdadera sin
# necesidad de conocer delta_H exactamente. delta_H exacto se
# calcula en R1 y R3 (soportes pequenos) como calibracion.
#
# Criterio C1+C2 del doc 102 (seis puntos: N=999982 y N=2^20;
# soportes S11, S13, S17):
#   C1: Bhat/B_M <= 10 en los seis puntos.
#   C2: pendiente log-log de Bhat/B_M contra |H_N| <= 0,1.
# Se aplica a la version MAS FAVORABLE de la candidata (Bhat_lo).
#
# Variante V103.1 (mayorante diagonal absoluto, [T] por triangulo):
#   |S_r| <= sum_{d|r} phi(d) |E_d| + |Delta_r|
#   Bhat2 := sum_{r en Hcal} (1/phi(r)) *
#            (sum_{d|r} phi(d)|E_d| + |Delta_r|)^2  >=  B_M.
# Su cociente Bhat2/B_M mide el precio de matar TODA la cancelacion
# (cruzados de Mobius y cruzado con Delta). Criterio previo propio:
# el mismo C1+C2 sobre los seis puntos.
#
# Descomposicion primo/compuesto: B_M = B_M^P + B_M^C con S_r
# corregido; Diag del doc 102 (diagonal de |M_r|^2) recomputado;
# cotas por pieza con exhibidos propios (pares de primos para P,
# gemelos compuestos para C).
#
# Chequeos:
#   K0  Maquinaria: Model97 R2 |H|=880; build_H z=11 == Model97.
#   K1  B_M^(M) (sin Delta, definicion del doc 102) reproduce la
#       tabla generica del doc 102 (rel < 2e-3) y Diag^C/B_M^C
#       reproduce 1.48-1.71 (rel < 3e-2).
#   K2  B_M^(S) = B_M^(M) + B_cross + B_Delta exacto por construccion
#       y cambio relativo en la escala 1e-4..1e-2 (reporte).
#   K3  Sandwich de delta en R1/R3: 1/Q_H^2 <= delta_exacto y
#       delta_exacto <= exhibido gemelo (con tolerancia numerica).
#   K4  Bhat2 >= B_M^(S) en todo punto (mayorante valido).
#   K5  C1+C2 para Bhat_lo y Bhat2 en los seis puntos (reporte).
# Uso: python3 103_criba_restringida_segundo_factor.py
import importlib.util
import math
import sys
import time
from fractions import Fraction

spec = importlib.util.spec_from_file_location(
    "s102", "102_experimentos_e101_1_e101_2.py")
s102 = importlib.util.module_from_spec(spec)
spec.loader.exec_module(s102)

Model97 = s102.Model97
build_H = s102.build_H
SegPsi = s102.SegPsi
M_r_E = s102.M_r_E
ols = s102.ols
divisors_of = s102.divisors_of
phi_of = s102.phi_of
pf_desc = s102.pf_desc
is_prime_spf = s102.is_prime_spf
lambda_sieve = s102.lambda_sieve
MU, SPF, SMALL_PRIMES = s102.MU, s102.SPF, s102.SMALL_PRIMES


# ------------------------------------------- Delta_r (doc 100, copia)
def P_e(e, p, N, mu, spf):
    """P_e = sum_{d|e} mu(e/d) A_d,
    A_d = phi(d) * sum_{j: p^j<=N, p^j=N (mod d)} log p.
    (doc 100, sin cambios)"""
    tot = 0.0
    for d in divisors_of(e, spf):
        s = 0.0
        pw = p
        while pw <= N:
            if pw % d == N % d:
                s += math.log(p)
            if pw > N // p:
                break
            pw *= p
        tot += mu[e // d] * phi_of(d, SMALL_PRIMES) * s
    return tot


def delta_r_orth(r, N, mu, spf):
    """Delta_r = sum_{p|r, p^2 no | r} P_{r/p} (doc 100, sin cambios)."""
    tot = 0.0
    for p in pf_desc(r, spf):
        e = r // p
        if e % p == 0:
            continue
        tot += P_e(e, p, N, mu, spf)
    return tot


# ------------------------------------------------------------------
_LAM2 = {}


def lambda2(N):
    if N not in _LAM2:
        lam = lambda_sieve(N)
        _LAM2[N] = sum(x * x for x in lam)
    return _LAM2[N]


def exact_delta(Hcal):
    """Espaciamiento minimo exacto de {a/r mod 1} (Fraction)."""
    fr = set()
    for r in Hcal:
        for a in range(1, r + 1):
            if math.gcd(a, r) == 1:
                fr.add(Fraction(a, r))
    fs = sorted(fr)
    best = Fraction(1, 1)
    for i in range(len(fs)):
        gap = fs[(i + 1) % len(fs)] - fs[i]
        if i == len(fs) - 1:
            gap += 1
        if gap < best:
            best = gap
    return best, len(fs)


def measure103(H, N, kappa=1.0):
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
    Hset = set(HN)
    BMs = BMs_P = BMs_C = 0.0
    BMm = BMm_P = BMm_C = 0.0
    Bcross = BD = 0.0
    Diag_P = Diag_C = 0.0
    Bhat2 = Bhat2_P = Bhat2_C = 0.0
    for r in Hcal:
        Mr = M_r_E(r, E, MU, SPF)
        Dr = delta_r_orth(r, N, MU, SPF)
        Sr = Mr + Dr
        ph = phi[r]
        tS = Sr * Sr / ph
        tM = Mr * Mr / ph
        Bcross += 2.0 * Mr * Dr / ph
        BD += Dr * Dr / ph
        s2 = 0.0
        sabs = 0.0
        for d in divisors_of(r, SPF):
            if d == 1:
                continue
            e2 = E[d] * E[d]
            s2 += phi[d] * phi[d] * e2
            sabs += phi[d] * abs(E[d])
        diag = s2 / ph
        t2 = (sabs + abs(Dr)) ** 2 / ph
        if is_prime_spf(r):
            BMs_P += tS
            BMm_P += tM
            Diag_P += diag
            Bhat2_P += t2
        else:
            BMs_C += tS
            BMm_C += tM
            Diag_C += diag
            Bhat2_C += t2
        BMs += tS
        BMm += tM
        Bhat2 += t2
    Q_H = max(Hcal)
    # exhibido gemelo (r, r+2) para delta_H (H es impar)
    r2 = max(r for r in Hcal if (r + 2) in Hset)
    L2 = lambda2(N)
    Bhat_lo = (N + r2 * (r2 + 2) / 2.0) * L2
    Bhat_hi = (N + Q_H * Q_H) * L2
    # exhibidos por pieza: primos (mejor par consecutivo en Hcal)
    primes_H = sorted(r for r in Hcal if is_prime_spf(r))
    bestP = None
    for i in range(len(primes_H) - 1):
        p, q = primes_H[i], primes_H[i + 1]
        v = q * p / (q - p)          # 1/delta_exhibido
        if bestP is None or v > bestP[0]:
            bestP = (v, p, q)
    Bhat_lo_P = (N + bestP[0]) * L2
    # compuestos: mayor gemelo compuesto
    comp_H = [r for r in Hcal if not is_prime_spf(r)]
    twinsC = [r for r in comp_H if (r + 2) in Hset
              and not is_prime_spf(r + 2)]
    r2C = max(twinsC) if twinsC else -1
    Bhat_lo_C = (N + r2C * (r2C + 2) / 2.0) * L2 if r2C > 0 \
        else float("nan")
    return {
        "N": N, "nHN": len(HN), "nHcal": len(Hcal), "R0": R0,
        "BMs": BMs, "BMs_P": BMs_P, "BMs_C": BMs_C,
        "BMm": BMm, "BMm_P": BMm_P, "BMm_C": BMm_C,
        "Bcross": Bcross, "BD": BD,
        "Diag_P": Diag_P, "Diag_C": Diag_C,
        "Q_H": Q_H, "r2": r2, "Lam2": L2,
        "Bhat_lo": Bhat_lo, "Bhat_hi": Bhat_hi,
        "Bhat2": Bhat2, "Bhat2_P": Bhat2_P, "Bhat2_C": Bhat2_C,
        "Bhat_lo_P": Bhat_lo_P, "Bhat_lo_C": Bhat_lo_C,
        "prime_pair": (bestP[1], bestP[2]),
        "r2C": r2C, "psiN": seg.psiN,
    }


def fmt_pt(tag, m):
    pc = ("pieza C: Bhat_lo^C=%.6e (x%.3e)"
          % (m["Bhat_lo_C"], m["Bhat_lo_C"] / m["BMs_C"])) \
        if m["BMs_C"] > 0 else "pieza C: n/a (B_M^C=0)"
    dc = ("Diag^C/B_M^C(M)=%.3f Diag^C/B_M^C(S)=%.3f"
          % (m["Diag_C"] / m["BMm_C"], m["Diag_C"] / m["BMs_C"])) \
        if m["BMs_C"] > 0 else "Diag^C: n/a"
    return ("     [%s] N=%d |H_N|=%d |Hcal|=%d Q_H=%d r2=%d "
            "Lam2=%.6e\n"
            "        B_M^(S)=%.6e (P %.6e / C %.6e)  B_M^(M)=%.6e  "
            "cambio rel=%+.3e\n"
            "        Bhat_lo=%.6e (x%.3e)  Bhat_hi=%.6e (x%.3e)  "
            "Bhat2=%.6e (x%.3f)\n"
            "        pieza P: Bhat_lo^P=%.6e (x%.3e)  %s  %s"
            % (tag, m["N"], m["nHN"], m["nHcal"], m["Q_H"], m["r2"],
               m["Lam2"], m["BMs"], m["BMs_P"], m["BMs_C"], m["BMm"],
               (m["Bcross"] + m["BD"]) / m["BMm"],
               m["Bhat_lo"], m["Bhat_lo"] / m["BMs"],
               m["Bhat_hi"], m["Bhat_hi"] / m["BMs"],
               m["Bhat2"], m["Bhat2"] / m["BMs"],
               m["Bhat_lo_P"], m["Bhat_lo_P"] / m["BMs_P"],
               pc, dc))


def main():
    allok = True
    t00 = time.time()
    print("== DOC 103: criba grande restringida a H_A sobre B_M ==")
    print("Solo stdlib. MEDICION, no prueba. Goldbach sigue abierta.")
    print()
    m1 = Model97(3, 11, 60, 20000, "fixed")
    m2 = Model97(3, 11, 3000, 200000, "fixed")
    m3 = Model97(3, 3, 349, 200000, "fixed")
    supports = {"R1": m1.H, "R2": m2.H, "R3": m3.H}
    for z in (11, 13, 17):
        D = z ** 5
        y = int(round(z ** (10.0 / 3.0)))
        Hz, _, _ = build_H(3, z, y, D, MU, SPF)
        supports["S%d" % z] = Hz
    print("Soportes: R1=%d R3=%d R2=%d S11=%d S13=%d S17=%d"
          % tuple(len(supports[t])
                  for t in ("R1", "R3", "R2", "S11", "S13", "S17")))
    # K0
    D11 = 11 ** 5
    y11 = int(round(11 ** (10.0 / 3.0)))
    m_s11 = Model97(3, 11, y11, D11, "fixed")
    k0 = len(m2.H) == 880 and set(supports["S11"]) == set(m_s11.H)
    allok &= k0
    print("K0 maquinaria (R2 |H|=880; build_H z=11 == Model97): %s"
          % ("PASS" if k0 else "FAIL"))
    print()
    # puntos de calibracion + delta exacto
    print("Calibracion (N=999982): R1 y R3, con delta_H exacto:")
    calib = {}
    for tag in ("R1", "R3"):
        t0 = time.time()
        m = measure103(supports[tag], 999982)
        calib[tag] = m
        print(fmt_pt("gen/" + tag, m) + " (%.0fs)" % (time.time() - t0))
        R0 = m["R0"]
        Hcal = [r for r in supports[tag]
                if math.gcd(r, 999982) == 1 and r > R0]
        d_ex, nfr = exact_delta(Hcal)
        lo = Fraction(1, m["Q_H"] ** 2)
        ex = Fraction(2, m["r2"] * (m["r2"] + 2))
        k3 = lo <= d_ex and d_ex <= ex
        allok &= k3
        print("        delta exacto = %.6e (1/%.4e), #frac=%d; "
              "sándwich [1/Q_H^2=%.3e, 2/(r2(r2+2))=%.3e] %s"
              % (float(d_ex), 1.0 / float(d_ex), nfr, float(lo),
                 float(ex), "OK" if k3 else "FAIL"))
    print()
    print("Seis puntos del criterio C1+C2 (doc 102):")
    plan = [(999982, ["S11", "S13", "S17"], "gen"),
            (1048576, ["S11", "S13", "S17"], "2^e")]
    six = []
    allpts = {"gen/R1": calib["R1"], "gen/R3": calib["R3"]}
    for (N, tags, fam) in plan:
        for tag in tags:
            t0 = time.time()
            m = measure103(supports[tag], N)
            six.append(("%s/%s" % (fam, tag), m))
            allpts["%s/%s" % (fam, tag)] = m
            print(fmt_pt("%s/%s" % (fam, tag), m)
                  + " (%.0fs)" % (time.time() - t0))
    print()
    # K1: B_M^(M) y Diag^C/B_M^C contra el doc 102 (familia generica)
    doc102_BM = {"gen/R1": 8.244e7, "gen/R3": 1.898e8,
                 "gen/S11": 3.366e9, "gen/R2": None,
                 "gen/S13": 7.644e9, "gen/S17": 2.788e10}
    doc102_Diag = {"gen/R1": 1.705, "gen/S11": 1.495,
                   "gen/S13": 1.488, "gen/S17": 1.533,
                   "2^e/S11": 1.674, "2^e/S13": 1.616,
                   "2^e/S17": 1.576}
    k1 = True
    for tag, ref in doc102_BM.items():
        if ref is None:
            continue
        m = measure103(supports[tag.split("/")[1]], 999982) \
            if tag not in allpts else allpts[tag]
        rel = abs(m["BMm"] - ref) / ref
        if rel > 2e-3:
            k1 = False
            print("     K1 DISCREPA B_M^(M) [%s]: %.6e vs doc 102 "
                  "%.6e (rel %.2e)" % (tag, m["BMm"], ref, rel))
    for tag, ref in doc102_Diag.items():
        if tag not in allpts:
            continue
        m = allpts[tag]
        got = m["Diag_C"] / m["BMm_C"]
        if abs(got - ref) / ref > 3e-2:
            k1 = False
            print("     K1 DISCREPA Diag^C/B_M^C [%s]: %.3f vs "
                  "doc 102 %.3f" % (tag, got, ref))
    allok &= k1
    print("K1 reproduccion del doc 102 (B_M^(M) y Diag^C/B_M^C, "
          "familia generica + Diag en 2^e): %s"
          % ("PASS" if k1 else "FAIL (ver arriba)"))
    # K2: identidad B_M^(S) = B_M^(M) + B_cross + B_Delta
    k2 = True
    for tag, m in allpts.items():
        rec = m["BMm"] + m["Bcross"] + m["BD"]
        if abs(rec - m["BMs"]) / m["BMs"] > 1e-12:
            k2 = False
            print("     K2 FALLO [%s]" % tag)
    allok &= k2
    print("K2 B_M^(S) = B_M^(M) + B_cross + B_Delta exacto: %s"
          % ("PASS" if k2 else "FAIL"))
    # K4: Bhat2 mayorante
    k4 = all(m["Bhat2"] >= m["BMs"] * (1 - 1e-12)
             for m in allpts.values())
    allok &= k4
    print("K4 Bhat2 >= B_M^(S) en todo punto: %s"
          % ("PASS" if k4 else "FAIL"))
    print()
    # C1 + C2 sobre los seis puntos (version favorable: Bhat_lo)
    print("Criterio C1+C2 sobre los seis puntos (Bhat_lo, la "
          "version mas favorable de la candidata):")
    c1 = all(m["Bhat_lo"] / m["BMs"] <= 10 for (_, m) in six)
    xs = [math.log(m["nHN"]) for (_, m) in six]
    ys = [math.log(m["Bhat_lo"] / m["BMs"]) for (_, m) in six]
    b, a, se, ci, l0, l1 = ols(xs, ys)
    print("     C1 (Bhat_lo/B_M <= 10 en los 6): %s "
          "(min=%.3e, max=%.3e)"
          % ("PASS" if c1 else "FAIL",
             min(m["Bhat_lo"] / m["BMs"] for (_, m) in six),
             max(m["Bhat_lo"] / m["BMs"] for (_, m) in six)))
    print("     C2 (pendiente log-log <= 0.1): pendiente=%.4f "
          "+- %.4f (SE), IC95%% +-%.4f, LOO [%.4f, %.4f] -> %s"
          % (b, se, ci, l0, l1, "PASS" if b <= 0.1 else "FAIL"))
    ys2 = [math.log(m["Bhat2"] / m["BMs"]) for (_, m) in six]
    b2, a2, se2, ci2, l02, l12 = ols(xs, ys2)
    c1v = all(m["Bhat2"] / m["BMs"] <= 10 for (_, m) in six)
    print("     Variante V103.1 (Bhat2): C1 %s (min=%.3f, max=%.3f); "
          "C2 pendiente=%.4f +- %.4f, LOO [%.4f, %.4f] -> %s"
          % ("PASS" if c1v else "FAIL",
             min(m["Bhat2"] / m["BMs"] for (_, m) in six),
             max(m["Bhat2"] / m["BMs"] for (_, m) in six),
             b2, se2, l02, l12,
             "PASS" if b2 <= 0.1 else "FAIL"))
    print()
    print("RESUMEN: Bhat1 (MV restringida) C1=%s C2=%s; "
          "V103.1 (Bhat2) C1=%s C2=%s"
          % ("PASS" if c1 else "FAIL", "PASS" if b <= 0.1 else "FAIL",
             "PASS" if c1v else "FAIL",
             "PASS" if b2 <= 0.1 else "FAIL"))
    print("PASS" if allok else "FAIL (ver lineas marcadas)")
    print("Tiempo total: %.1f s" % (time.time() - t00))
    return 0 if allok else 1


if __name__ == "__main__":
    sys.exit(main())
