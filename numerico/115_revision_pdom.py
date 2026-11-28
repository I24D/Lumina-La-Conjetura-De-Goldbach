#!/usr/bin/env python3
# 115_revision_pdom.py  (Doc 115)
# Revision adversarial de P-DOM (proposicion de Lumina Chat en su
# respuesta al doc 112) + dos reparaciones de rigor (conjunto de
# niveles de (H2'')/(H2''') y cola de M_full del doc 110 seccion 3.1).
# Solo stdlib. Importa la maquinaria verificada de los docs 106/109/110
# sin modificarla.
#
# Chequeos:
#   J0  maquinaria: |Hcal| reproduce la tabla del doc 106 (gen).
#   J1  definiciones: Omega_q^P = 1/(q-1) por enumeracion; Hcal_C =
#       {r en Hcal: r no primo}; ningun elemento de H es par ni
#       divisible por 3 (A=3: primos permitidos (p,A)=1, p>2);
#       q en Hcal_P implica (q, A*N) = 1.
#   J2  membresia prima en V+: in_Vplus(p, L) <=> p^3 < L (rejilla).
#   J3  P-DOM finito: razon Omega_q^C/Omega_q^P por enumeracion en
#       z = 11..37, ambos N; minimos por caso (I: q<z; II: q>=z);
#       cotas de construccion LB_I/LB_II verificadas termino a
#       termino (cada r construido esta en Hcal_C); sumas de Mertens
#       S_I(z), S_II(z) contra log(1.9), log(2.5), tambien a
#       z = 10^2, 10^3, 10^4 (sumas puras sobre primos).
#   J4  identidad Theta: Omega_{qm}^C por enumeracion contra la
#       formula corregida (2 tau(h) - 1) y contra la transmitida
#       (2 (tau(h) - 1)); la discrepancia es (1/(q-1)) * sum 1/phi(h).
#   J5  M_full: producto parcial hasta 10^5, cota de cola
#       exp(1/(X-1)), intervalo riguroso; M_full(z) del doc 110.
#   J6  niveles: los niveles medidos en los docs 110-113 son
#       {D//s : s primo admisible muestreado} (subconjunto de L_A).
import importlib.util
import math
import sys
import time


def load(name, path):
    spec = importlib.util.spec_from_file_location(name, path)
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod


t0 = time.time()
d106 = load("d106", "106_precio_tau_conteo_omega.py")
d110 = load("d110", "110_segundo_momento_defectos.py")
build_H106 = d106.build_H106
in_Vplus = d106.in_Vplus
tau_phi = d106.tau_phi
is_prime = d106.is_prime106
primes_upto = d106.primes_upto
pf_trial = d106.pf_trial
SMALL_PRIMES = d106.SMALL_PRIMES

A_CH = 3
OK = True


def note(ok):
    global OK
    OK = OK and ok
    return ok


def divisors_of(n):
    pf = pf_trial(n, SMALL_PRIMES)
    ds = [1]
    for p in pf:
        ds = ds + [d * p for d in list(ds)]
    return ds


_cache = {}


def build_support(z, N):
    if (z, N) in _cache:
        return _cache[(z, N)]
    D = z ** 5
    y = int(round(z ** (10.0 / 3.0)))
    H, sprimes, Vs = build_H106(A_CH, z, y, D)
    R0 = math.log(N)
    HN = [r for r in H if math.gcd(r, N) == 1]
    Hcal = [r for r in HN if r > R0]
    out = (H, sprimes, Vs, Hcal)
    _cache[(z, N)] = out
    return out


def omega_prime_map(Hcal):
    """Omega_q^C solo sobre q primos, por enumeracion sobre Hcal_C."""
    OmC = {}
    for r in Hcal:
        if is_prime(r):
            continue
        tau, phi = tau_phi(r)
        w = (tau - 1) / phi
        for q in pf_trial(r, SMALL_PRIMES):
            OmC[q] = OmC.get(q, 0.0) + w
    return OmC


def omega_full_map(Hcal):
    """Omega_d^C para todo d, por enumeracion sobre Hcal_C."""
    OmC = {}
    for r in Hcal:
        if is_prime(r):
            continue
        tau, phi = tau_phi(r)
        w = (tau - 1) / phi
        for d in divisors_of(r):
            if d >= 2:
                OmC[d] = OmC.get(d, 0.0) + w
    return OmC


# ------------------------------------------------------------------ J0/J1
print("== J0/J1: maquinaria y definiciones ==")
# La tabla del doc 106 seccion 1.1 rotula |Hcal_N| pero cuenta H_N
# (antes del corte r > R0): se verifica |H_N| contra la tabla y se
# reporta |Hcal|; los conjuntos coinciden con Counting del doc 106.
expect = {11: 798, 13: 1621, 17: 5017, 19: 8219, 23: 17955}
j0 = True
for z, exp in expect.items():
    H, sprimes, Vs, Hcal = build_support(z, 999982)
    nHN = sum(1 for r in H if math.gcd(r, 999982) == 1)
    j0 = j0 and (nHN == exp)
    print("  z=%d gen |H_N|=%d (tabla doc 106: %d)  |Hcal|=%d"
          % (z, nHN, exp, len(Hcal)))
print("  J0:", "PASS" if note(j0) else "FAIL")

j1 = True
for z in (11, 17, 23):
    for N in (999982, 1048576):
        H, sprimes, Vs, Hcal = build_support(z, N)
        for r in H:
            if r % 2 == 0 or r % 3 == 0:
                j1 = False
        for q in [r for r in Hcal if is_prime(r)]:
            tau, phi = tau_phi(q)
            if abs((tau - 1) / phi - 1.0 / (q - 1)) > 1e-12:
                j1 = False
            if math.gcd(q, A_CH * N) != 1:
                j1 = False
print("  J1 (Omega^P=1/(q-1); H sin pares ni multiplos de 3; "
      "q en H_P => (q,AN)=1):", "PASS" if note(j1) else "FAIL")

# ------------------------------------------------------------------ J2
print("== J2: primo en V+(L) <=> p^3 < L ==")
j2 = True
for p in primes_upto(200):
    if p < 3:
        continue
    for L in list(range(2, 400)) + [1000, 5000, 10 ** 5, 10 ** 6]:
        if in_Vplus(p, L, [p]) != (p ** 3 < L):
            j2 = False
print("  J2:", "PASS" if note(j2) else "FAIL")

# ------------------------------------------------------------------ J3
print("== J3: P-DOM en ventana finita (enumeracion, sin psi) ==")
print("  %-4s %-8s %5s %8s %8s %8s %8s %8s" %
      ("z", "N", "|Hp|", "min", "minI", "minII", "med", "LBmin"))
for z in (11, 13, 17, 19, 23, 29, 31, 37):
    for N in (999982, 1048576):
        H, sprimes, Vs, Hcal = build_support(z, N)
        HcalC_set = set(r for r in Hcal if not is_prime(r))
        OmC = omega_prime_map(Hcal)
        ratios = {}
        lb = {}
        for q in [r for r in Hcal if is_prime(r)]:
            ratios[q] = OmC.get(q, 0.0) * (q - 1)
            if q < z:
                ssum = 0.0
                hi = z ** 1.9
                for s in primes_upto(int(hi) + 1):
                    if s >= z and s <= hi and math.gcd(s, A_CH * N) == 1:
                        assert q * s in HcalC_set, (q, s, z, N)
                        ssum += 1.0 / (s - 1)
                lb[q] = 3.0 * ssum
            else:
                ssum = 0.0
                lo, hi = z ** 0.2, z ** 0.5
                for p in primes_upto(int(hi) + 1):
                    if p >= lo and p <= hi and (A_CH * N) % p != 0:
                        assert q * p in HcalC_set, (q, p, z, N)
                        ssum += 1.0 / (p - 1)
                lb[q] = 3.0 * ssum
        vals = sorted(ratios.values())
        caseI = [ratios[q] for q in ratios if q < z]
        caseII = [ratios[q] for q in ratios if q >= z]
        med = vals[len(vals) // 2] if vals else float("nan")
        print("  %-4d %-8d %5d %8.3f %8s %8s %8.3f %8.3f" %
              (z, N, len(ratios), vals[0] if vals else -1,
               ("%.3f" % min(caseI)) if caseI else "-",
               ("%.3f" % min(caseII)) if caseII else "-",
               med, min(lb.values()) if lb else -1), flush=True)
        for q in ratios:
            assert lb[q] <= ratios[q] + 1e-9, (q, z, N)
print("  (LB = cota inferior de la construccion de P-DOM, verificada")
print("   termino a termino: cada r construido pertenece a Hcal_C,")
print("   y LB <= razon real en todos los casos: PASS)")
print("  J3b: caracterizacion de los ceros (A=3, primo permitido")
print("       minimo p0=5): razon(q)=0 <=> D/q <= 125 <=> V+(D/q)={1}")
for z in (11, 13, 17, 19, 23, 29, 37):
    N = 999982
    H, sprimes, Vs, Hcal = build_support(z, N)
    D = z ** 5
    OmC = omega_prime_map(Hcal)
    zeros = set(q for q in Hcal if is_prime(q)
                and OmC.get(q, 0.0) == 0.0)
    pred = set(q for q in Hcal if is_prime(q) and D // q <= 125)
    print("    z=%d: %d ceros; predichos por D/q<=125: %d; coinciden: %s"
          % (z, len(zeros), len(pred), zeros == pred))
    assert zeros == pred

print("  Sumas de Mertens de la construccion (sin soporte):")
print("  %-8s %10s %10s %10s %10s" % ("z", "S_I", "3*S_I", "S_II", "3*S_II"))
for z in (11, 17, 23, 29, 37, 100, 1000, 10000):
    si = sum(1.0 / (s - 1) for s in primes_upto(int(z ** 1.9) + 1)
             if s >= z)
    sii = sum(1.0 / (p - 1) for p in primes_upto(int(z ** 0.5) + 1)
              if p >= z ** 0.2)
    print("  %-8d %10.4f %10.4f %10.4f %10.4f" % (z, si, 3 * si, sii,
                                                 3 * sii))
print("  limites: log(1.9)=%.4f, 3log(1.9)=%.4f; log(2.5)=%.4f, "
      "3log(2.5)=%.4f" % (math.log(1.9), 3 * math.log(1.9),
                          math.log(2.5), 3 * math.log(2.5)))


# ------------------------------------------------------------------ J4
print("== J4: identidad Theta sobre modulos conmutados d = q*m ==")
j4a = True   # formula corregida (2 tau - 1)
worst_lumina = 0.0
for z in (17, 23):
    N = 999982
    H, sprimes, Vs, Hcal = build_support(z, N)
    OmC = omega_full_map(Hcal)
    for q in sprimes:
        if N % q == 0:
            continue
        Vq = Vs[q]
        ms = set()
        for h in Vq:
            ms.update(divisors_of(h))
        for m in ms:
            d = q * m
            enum = OmC.get(d, 0.0)
            corr = 0.0
            lum = 0.0
            for h in Vq:
                if h == 1 or h % m != 0:
                    continue
                tau, phi = tau_phi(h)
                corr += (2 * tau - 1) / phi
                lum += 2 * (tau - 1) / phi
            corr /= (q - 1)
            lum /= (q - 1)
            if abs(enum - corr) > 1e-9 * max(1.0, abs(enum)):
                j4a = False
                print("    FALLO corregida: z=%d q=%d m=%d enum=%.6f "
                      "corr=%.6f" % (z, q, m, enum, corr))
            worst_lumina = max(worst_lumina, abs(enum - lum))
print("  formula corregida (2 tau(h) - 1):",
      "PASS" if note(j4a) else "FAIL")
print("  discrepancia maxima de la formula transmitida 2(tau(h)-1): "
      "%.6f (es (1/(q-1)) * sum_{h} 1/phi(h) > 0: subestima Omega)"
      % worst_lumina)

# ------------------------------------------------------------------ J5
print("== J5: cola de M_full ==")
X = 10 ** 5
Mle = 1.0
for p in primes_upto(X):
    if p > 3:
        Mle *= 1.0 + 1.0 / ((p - 1) * (p - 2))
tail = math.exp(1.0 / (X - 1))
print("  prod_{3<p<=10^5} (1 + 1/((p-1)(p-2))) = %.10f" % Mle)
print("  factor de cola <= exp(1/(X-1)) = %.10f" % tail)
print("  M_full(inf) en [%.10f, %.10f]" % (Mle, Mle * tail))
print("  el enunciado del doc 110 'M_full(inf) <= 1,1618' es falso "
      "tal cual: el producto truncado es cota INFERIOR")
for z in (17, 29, 37):
    mz_lt = 1.0
    mz_le = 1.0
    for p in primes_upto(z):
        if p > 3:
            mz_le *= 1.0 + 1.0 / ((p - 1) * (p - 2))
            if p < z:
                mz_lt *= 1.0 + 1.0 / ((p - 1) * (p - 2))
    print("  M_full(z=%d): p<z -> %.4f ; p<=z -> %.4f (doc 110: %s)"
          % (z, mz_lt, mz_le, {17: "1,1405", 29: "1,1514",
                               37: "1,1543"}[z]))

# ------------------------------------------------------------------ J6
print("== J6: los niveles medidos (docs 110-113) ==")
frames = [("A", 17, 17 ** 5, int(17 ** (10.0 / 3.0))),
          ("A", 29, 29 ** 5, int(29 ** (10.0 / 3.0))),
          ("A", 37, 37 ** 5, int(37 ** (10.0 / 3.0))),
          ("B", 29, 29 ** 6, 29 ** 3)]
tot = 0
j6 = True
for tag, z, D, y in frames:
    spr = d110.sample_sprimes(z, y, 3, (), 60)
    levels = sorted({D // s for s in spr})
    tot += len(levels)
    full = [p for p in primes_upto(y - 1)
            if p >= z and math.gcd(p, 3) == 1]
    LA = {D // s for s in full}
    if not set(levels) <= LA:
        j6 = False
    print("  marco %s z=%d: %d niveles medidos, todos en L_A "
          "(|L_A| = %d)" % (tag, z, len(levels), len(LA)))
print("  total niveles medidos: %d" % tot)
print("  J6:", "PASS" if note(j6) else "FAIL")

print()
print("Tiempo total: %.1f s" % (time.time() - t0))
print("PASS" if OK else "FAIL")
