#!/usr/bin/env python3
# 121a_unicidad_rama.py — Verificacion del lema de unicidad de rama (doc 121a)
# Sobre el soporte VERDADERO de la serie (build_H estilo doc 102/106):
#  K1 marcos reales S11..S23 (A=3, D=z^5, y=round(z^(10/3))):
#     (a) max_{r en H_C} #{q en S : q|r} <= 1 y max_{r,m} #{q : qm|r} <= 1
#     (b) alpha~_m = sum_{q en S} Omega_{qm}^C <= Omega_m^C para todo m en D
#     (c) intercambio con el peso VERDADERO (tau(r)-1)/phi(r): identidad
#         exacta (Fraction) para todo m en D
#     (d) el peso de Lumina Chat, tau(r)^{-1}/phi(r), NO da la identidad
#  K2 el juguete 3*5*7 del doc 120b: 105 no esta en H en los marcos reales
#     (3|A: ningun r de H es divisible por 3), y en el marco hipotetico
#     A=1, z=5 construido con la misma maquinaria tampoco (105 necesitaria
#     dos primos >= z); ahi ademas (a) y (b) valen y los conjuntos de r
#     que aportan a Omega_{15} y Omega_{21} son disjuntos.
# Solo stdlib. PASS si todos los chequeos se comportan como afirma el doc.
import importlib.util
import math
import os
import sys
import time
from fractions import Fraction

os.chdir(os.path.dirname(os.path.abspath(__file__)))
spec = importlib.util.spec_from_file_location(
    "s105", "105_forma_diagonal_segundos_momentos.py")
s105 = importlib.util.module_from_spec(spec)
spec.loader.exec_module(s105)
s102 = s105.s104.s103.s102
pf_trial = s102.pf_trial
divisors_from_pf = s102.divisors_from_pf
in_Vplus = s102.in_Vplus
primes_upto = s102.primes_upto
SMALL_PRIMES = s102.SMALL_PRIMES

fails = []


def check(name, ok, detail=""):
    print(("OK    " if ok else "FAIL  ") + name + ("  " + detail if detail else ""))
    if not ok:
        fails.append(name)


def build_H(A, z, y, D):
    """build_H106 del doc 106, verbatim en construccion."""
    allpr = primes_upto(max(y, 100))
    sprimes = [p for p in allpr if p >= z and p < y and math.gcd(p, A) == 1]
    allow = [p for p in allpr if 2 < p < z and math.gcd(p, A) == 1]
    divs = [1]
    for p in allow:
        divs += [d * p for d in list(divs)]
    Mset = set()
    for s in sprimes:
        Ds = D // s
        for d in divs:
            if d >= Ds:
                continue
            if d == 1 or in_Vplus(d, Ds, pf_trial(d, SMALL_PRIMES)):
                Mset.add(s * d)
    H = {1}
    for m in Mset:
        for h in divisors_from_pf(pf_trial(m, SMALL_PRIMES)):
            H.add(h)
    return sorted(H), sprimes


def phi_of_pf(n, pf):
    ph = n
    for p in pf:
        ph = ph // p * (p - 1)
    return ph


def analyse(A, z, y, D, tag):
    H, S = build_H(A, z, y, D)
    Sset = set(S)
    HC = []
    pf_of = {}
    for r in H:
        pf = pf_trial(r, SMALL_PRIMES)
        pf_of[r] = pf
        if len(pf) >= 2:
            HC.append(r)
    # Omega_d^C por enumeracion desde la definicion (doc 105/106)
    Om = {}
    for r in HC:
        pf = pf_of[r]
        w = Fraction(2 ** len(pf) - 1, phi_of_pf(r, pf))
        for d in divisors_from_pf(pf):
            if d >= 2:
                Om[d] = Om.get(d, Fraction(0)) + w
    # (a) conteos
    max_qr = 0
    max_qmr = 0
    for r in HC:
        qs = [p for p in pf_of[r] if p in Sset]
        max_qr = max(max_qr, len(qs))
        for m in divisors_from_pf(pf_of[r]):
            c = sum(1 for q in qs if r % (q * m) == 0)
            max_qmr = max(max_qmr, c)
    check(f"{tag} (a) max #{{q en S : q|r}} <= 1", max_qr <= 1,
          f"max={max_qr} sobre |H_C|={len(HC)}")
    check(f"{tag} (a) max #{{q en S : qm|r}} <= 1", max_qmr <= 1,
          f"max={max_qmr}")
    # (b) dominacion y (c) intercambio exacto, para todo m en D
    bad_dom = []
    bad_id = []
    for m, Om_m in Om.items():
        lhs = Fraction(0)
        rhs = Fraction(0)
        for q in S:
            if q * m in Om:
                lhs += Om[q * m]
        for r in HC:
            if r % m == 0:
                pf = pf_of[r]
                cnt = sum(1 for q in Sset if r % (q * m) == 0)
                rhs += Fraction(2 ** len(pf) - 1, phi_of_pf(r, pf)) * cnt
        if lhs > Om_m:
            bad_dom.append(m)
        if lhs != rhs:
            bad_id.append(m)
    check(f"{tag} (b) alpha~_m <= Omega_m^C para todo m en D",
          not bad_dom, f"|D|={len(Om)}, violaciones={bad_dom[:5]}")
    check(f"{tag} (c) intercambio con peso (tau-1)/phi: identidad exacta",
          not bad_id, f"fallos={bad_id[:5]}")
    return H, S, HC, Om, pf_of


def main():
    t0 = time.time()
    print("== DOC 121a: lema de unicidad de rama sobre el soporte real ==")
    frames = {}
    for z in (11, 13, 17, 19, 23):
        D = z ** 5
        y = int(round(z ** (10.0 / 3.0)))
        frames[z] = analyse(3, z, y, D, f"S{z}")
        print()
    # (d) el peso tau^{-1} de su identidad escrita NO reproduce sum_q Omega
    H, S, HC, Om, pf_of = frames[11]
    m_ex = None
    for m in sorted(Om):
        lhs = sum(Om[q * m] for q in S if q * m in Om)
        rhs_hers = Fraction(0)
        for r in HC:
            if r % m == 0:
                pf = pf_of[r]
                cnt = sum(1 for q in S if r % (q * m) == 0)
                rhs_hers += Fraction(1, (2 ** len(pf)) * phi_of_pf(r, pf)) * cnt
        if lhs != rhs_hers:
            m_ex = (m, lhs, rhs_hers)
            break
    check("S11 (d) su peso tau(r)^{-1}/phi(r) NO da la identidad",
          m_ex is not None,
          ("m=%d: sum_q Omega_{qm}=%.6f vs su intercambio=%.6f"
           % (m_ex[0], float(m_ex[1]), float(m_ex[2]))) if m_ex else "")
    print()
    # K2: el juguete 3*5*7
    H11 = set(frames[11][0])
    check("K2 105 no pertenece a H en S11 (marco real)", 105 not in H11)
    check("K2 ningun r de H en S11 es divisible por 3 (3 | A=3)",
          all(r % 3 != 0 for r in H11))
    Hh, Sh, HCh, Omh, pfh = analyse(1, 5, int(round(5 ** (10.0 / 3.0))),
                                     5 ** 5, "HIP(A=1,z=5)")
    check("K2 105 tampoco pertenece a H en el marco hipotetico A=1, z=5",
          105 not in set(Hh),
          "105 = 3*5*7 exigiria dos primos >= z (5 y 7) en un solo r")
    R15 = {r for r in HCh if r % 15 == 0}
    R21 = {r for r in HCh if r % 21 == 0}
    check("K2 en el marco hipotetico, los r que aportan a Omega_15 y a "
          "Omega_21 son disjuntos", not (R15 & R21),
          f"|R15|={len(R15)}, |R21|={len(R21)}, interseccion={len(R15 & R21)}")
    print()
    print("Tiempo total: %.1f s" % (time.time() - t0))
    print("RESULTADO:", "PASS" if not fails else f"FAIL en {fails}")
    return 0 if not fails else 1


if __name__ == "__main__":
    sys.exit(main())
