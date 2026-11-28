#!/usr/bin/env python3
# 116b_flujo_afinado.py  (Doc 116b, seccion sobre (R-FLOW-71+))
# Residuo del doc 114: R_A(y), R_R(y) <= 1 para min(y) >= 71.
#
# Resultado: el enunciado PUNTUAL es FALSO [T]. Contraejemplos en
# marco verificados con los predicados de la serie (in_Vplus del
# doc 109, in_Vhat del doc 113, convencion de lados() del doc 114):
#   J2: r=5, y=Y primo, L en (max(Y^3, 125*Y*q_max), 5*Y^3]:
#       todo q en P con 5<q<Y es admisible (q en posicion fusionada
#       2, par; el unico primo volteado, 5, ata solo en q < L/(125Y)
#       > Y). R_A(y) = sum_{5<q<Y} 1/(q-1) -> supera 1 (Y ~ 307) y
#       tiende a infinito como log log Y.
#   J3: analogo para R_R con r=rho primo, y=Y1*Y2, el fallo
#       fusionado en la posicion 3 (primo Y2), L ~ rho*Y1*Y2^3.
# Las cotas afiladas que SI valen [T] (derivadas de las condiciones
# fusionadas de r*y*q; J1 las verifica en las 4113 instancias del
# doc 114, 0 violaciones):
#   (*)   si t(q) = #{p|r: p>q} es impar, q^3 < L*s_q/(r*y),
#         s_q = prod_{p|r, p<q} p  (condicion en la posicion de q).
#   (**)  si p|r, p<q y #{p'|r: p'>p} es par (p en posicion par en
#         r*y), entonces q < L/(p^2 * y * r_{>=p})  (primo volteado).
#   (***) en particular q > max(r) => q < L/(max(r)^3 * y).
# El contraejemplo es la UNICA forma de escapar: q en posicion par
# (t(q) par) y ningun primo de r bajo q en posicion par que ate
# dentro del rango (con r primo pequeno, el volteo ata en
# q < L/(r^3 y), fuera del rango q<Y si L ~ r*Y^3).
# J4: medicion extendida z=71 y z=101 [M]: maximos de R_A y R_R
# por cubeta de min(y); en z<=101 la suma de Mertens disponible
# (<= Mertens(97) - 1/4 < 1 en la familia) no alcanza a violar.
# Uso: python3 116b_flujo_afinado.py | tee salida_doc116b_local.txt
import importlib.util
import math
import sys
import time
from fractions import Fraction


def load(name, path):
    spec = importlib.util.spec_from_file_location(name, path)
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod


d109 = load("doc109", "109_ataque_h2_gram.py")
d114 = load("doc114", "114_desp_sign_banda.py")
d113 = d114.d113

primes_upto = d109.primes_upto
allowed_primes = d109.allowed_primes
in_Vplus = d109.in_Vplus
in_Vhat = d113.in_Vhat
sample_sprimes = d109.sample_sprimes if hasattr(d109, "sample_sprimes") \
    else d114.d113.d112.d111.sample_sprimes

PRIMES = primes_upto(100000)
PRIME_SET = set(PRIMES)


def is_prime(n):
    if n < 2:
        return False
    if n in PRIME_SET:
        return True
    i = 2
    while i * i <= n:
        if n % i == 0:
            return False
        i += 1 if i == 2 else 2
    return True


def prod(pf):
    out = 1
    for p in pf:
        out *= p
    return out


def side_of(dpf, r, rpf, L):
    """Convencion de lados() del doc 114, con predicados locales:
    devuelve 'Add'/'Rem'/'Both'/'Neither' para d' con (d', r) = 1."""
    d = prod(dpf)
    if r * d >= L:
        inv = False
        inhat = False
    else:
        merged = tuple(sorted(tuple(dpf) + tuple(rpf), reverse=True))
        inv = in_Vplus(merged, L)
        inhat = in_Vhat(dpf, L, r)
    if inv and not inhat:
        return "Add"
    if inhat and not inv:
        return "Rem"
    if inv and inhat:
        return "Both"
    return "Neither"


def mu_of(pf):
    return -1 if len(pf) % 2 == 1 else 1


def flow_sum(y, ypf, r, rpf, L, primes_z, want_from, want_to):
    """Suma exacta sum_{q: yq en want_to} 1/(q-1) para y en want_from,
    mas la lista de q admisibles. want_from/want_to: ('Add',-1) etc."""
    tot = Fraction(0)
    qs = []
    mn = ypf[-1]
    for q in primes_z:
        if q >= mn:
            break
        if (r * y) % q == 0:
            continue
        qpf = tuple(sorted(ypf + (q,), reverse=True))
        if side_of(qpf, r, rpf, L) == want_to[0] \
                and mu_of(qpf) == want_to[1]:
            tot += Fraction(1, q - 1)
            qs.append(q)
    return tot, qs


# ------------------------------------------------------------------ J0
def check_J0():
    print("J0 sumas de Mertens exactas (fracciones):")
    def mert(n):
        return sum((Fraction(1, p - 1) for p in PRIMES
                    if 5 <= p <= n), Fraction(0))
    m67 = mert(67)
    m71 = mert(71)
    m199 = mert(199)
    m251 = mert(251)
    print("  sum_{5<=p<=67} 1/(p-1) = %s = %.5f  (< 1: %s)"
          % (m67, float(m67), m67 < 1))
    print("  sum_{5<=p<=71} 1/(p-1) = %.5f  (> 1: %s)"
          % (float(m71), m71 > 1))
    print("  Mertens(199) = %.5f | Mertens(251) = %.5f"
          % (float(m199), float(m251)))
    ok = abs(float(m67) - 0.99917) < 5e-6 and m67 < 1 < m71
    print("  => reproduce el umbral del doc 114. %s"
          % ("PASS" if ok else "FAIL"))
    return ok


# ------------------------------------------------------------------ J1
def check_J1():
    print("J1 cotas afiladas (*)/(**)/(***) en las instancias del "
          "doc 114 (0 violaciones esperadas):")
    viol = 0
    npar = 0
    maxRA = 0.0
    for (tag, z, a, r, ld) in d114.instancias():
        lv = a.lv
        primes_z = allowed_primes(z, 3, ())
        rpf = None
        for (d, pf, mu, phi) in lv.enum:
            if d == r:
                rpf = pf
                break
        rmax = rpf[0]
        for (y, (side, ypf, ymu, yphi)) in ld.items():
            if side != "Add" or ymu != -1:
                continue
            tot = 0.0
            for q in primes_z:
                if q >= ypf[-1]:
                    break
                if (r * y) % q == 0:
                    continue
                x = y * q
                if x not in ld:
                    continue
                (s2, xpf, xmu, xphi) = ld[x]
                if s2 != "Add" or xmu != 1:
                    continue
                npar += 1
                tot += 1.0 / (q - 1)
                # (*) t(q) impar => q^3 * r * y < L * s_q
                tq = sum(1 for p in rpf if p > q)
                if tq % 2 == 1:
                    sq = prod(tuple(p for p in rpf if p < q))
                    if not (q ** 3 * r * y < lv.L * sq):
                        viol += 1
                # (**) p|r, p<q, #{p'>p} par => q*p^2*y*r_{>=p} < L
                for p in rpf:
                    if p >= q:
                        continue
                    u = sum(1 for pp in rpf if pp > p)
                    if u % 2 == 0:
                        rge = prod(tuple(pp for pp in rpf if pp >= p))
                        if not (q * p * p * y * rge < lv.L):
                            viol += 1
                # (***) q > max(r) => q * max(r)^3 * y < L
                if q > rmax and not (q * rmax ** 3 * y < lv.L):
                    viol += 1
            maxRA = max(maxRA, tot)
    print("  pares (y,q) admisibles verificados: %d | violaciones: %d"
          " | max R_A en ventana = %.4f (doc 114: 0,6000)"
          % (npar, viol, maxRA))
    ok = viol == 0 and abs(maxRA - 0.6) < 1e-9
    print("  => las tres cotas son necesarias en cada instancia. %s"
          % ("PASS" if ok else "FAIL"))
    return ok


# ------------------------------------------------------------------ J2
def find_level(D, lo, hi, smax):
    """Primer primo s <= smax con D//s en (lo, hi]; None si no hay.
    Busca desde s ~ D/hi hacia arriba (L = D//s decrece con s)."""
    s = max(5, D // (hi + 1) - 100)
    if s % 2 == 0:
        s += 1
    while s <= smax:
        if is_prime(s):
            L = D // s
            if lo < L <= hi:
                return s, L
            if L <= lo:
                break
        s += 2
    return None, None


def check_J2():
    print("J2 contraejemplo R_A en marco (r=5, y=Y primo, z=331):")
    z = 331
    D = z ** 5
    smax = int(z ** (10.0 / 3.0))
    primes_z = allowed_primes(z, 3, ())
    r, rpf = 5, (5,)
    rows = []
    first = None
    for Y in [p for p in PRIMES if 73 <= p <= 307]:
        qmax = max(p for p in PRIMES if p < Y)
        lo = max(Y ** 3, 125 * Y * qmax)
        hi = 5 * Y ** 3
        s, L = find_level(D, lo, hi, smax)
        if s is None:
            rows.append((Y, None, None, None))
            continue
        # y en Add^- ?
        sy = side_of((Y,), r, rpf, L)
        tot, qs = flow_sum(Y, (Y,), r, rpf, L, primes_z,
                           ("Add", -1), ("Add", 1))
        closed = sum((Fraction(1, q - 1) for q in PRIMES
                      if 5 < q < Y), Fraction(0))
        match = (tot == closed)
        rows.append((Y, s, L, tot))
        if sy == "Add" and match and tot > 1 and first is None:
            first = (Y, s, L, tot, qs)
    for (Y, s, L, tot) in rows:
        if Y in (73, 101, 151, 211, 251, 283, 293, 307) or \
                (tot is not None and tot > 1):
            print("  Y=%d: s=%s L=%s R_A(y)=%s"
                  % (Y, s, L, ("%.5f" % float(tot)) if tot is not None
                     else "sin nivel"))
    ok = first is not None
    if ok:
        (Y, s, L, tot, qs) = first
        print("  PRIMER CRUCE: Y=%d, s=%d, L=%d, R_A(y) = %s = %.5f > 1"
              " | #q admisibles=%d (todos los de (5, Y))"
              % (Y, s, L, tot, float(tot), len(qs)))
        # y visible: r | 5*Y en V^+(L)
        vis = in_Vplus((Y, 5), L) and 5 * Y < L
        print("  r=5 visible en el nivel (5*Y en V^+(L)): %s" % vis)
        ok = ok and vis
    print("  => R_A(y) <= 1 es FALSO como enunciado puntual, en "
          "marco. %s" % ("PASS" if ok else "FAIL"))
    return ok


# ------------------------------------------------------------------ J3
def check_J3():
    print("J3 contraejemplo R_R en marco (r=rho, y=Y1*Y2):")
    rho, Y1, Y2 = 4999, 5009, 4099
    assert is_prime(rho) and is_prime(Y1) and is_prime(Y2)
    y = Y1 * Y2
    ypf = (Y1, Y2)
    rpf = (rho,)
    lo = rho * y * 251 ** 3          # q=251 admisible: rho*q^3*y < L
    hi = rho * Y1 * Y2 ** 3          # fallo fusionado pos.3: >= L
    found = None
    for z in range(6050, 6150):
        D = z ** 5
        for s in (5, 7, 11, 13, 17, 19, 23, 29, 31):
            L = D // s
            if lo < L <= hi:
                found = (z, s, L)
                break
        if found:
            break
    if not found:
        print("  sin (z, s) con L en ventana")
        print("  FAIL")
        return False
    (z, s, L) = found
    primes_z = allowed_primes(z, 3, ())
    sy = side_of(ypf, rho, rpf, L)
    tot, qs = flow_sum(y, ypf, rho, rpf, L, primes_z,
                       ("Rem", 1), ("Rem", -1))
    vis = in_Vplus((Y1, rho), L) and Y1 * rho < L
    print("  z=%d s=%d L=%d | y=%d en %s | rho visible: %s"
          % (z, s, L, y, sy, vis))
    print("  R_R(y) = %s = %.5f | #q admisibles=%d, max q=%d"
          % (tot, float(tot), len(qs), max(qs) if qs else -1))
    ok = sy == "Rem" and vis and tot > 1
    print("  => R_R(y) <= 1 tambien es FALSO como enunciado puntual, "
          "en marco. %s" % ("PASS" if ok else "FAIL"))
    return ok


# ------------------------------------------------------------------ J4
def enum_subsets(primes_list, max_prod):
    """Todos los (d, pf_desc) libres de cuadrados con d <= max_prod
    sobre la lista de primos dada (DFS con poda por producto)."""
    out = []
    pls = sorted(primes_list)

    def dfs(i, d, pf):
        if d > 1:
            out.append((d, tuple(sorted(pf, reverse=True))))
        for j in range(i, len(pls)):
            nd = d * pls[j]
            if nd > max_prod:
                break
            dfs(j + 1, nd, pf + [pls[j]])
    dfs(0, 1, [])
    return out


def measure_frame(z, levels_L, r_pfs, ysubs, label):
    """Maximos de R_A / R_R por cubeta de min(y) sobre los (L, r, y)
    dados. ysubs: lista de (y, ypf). Devuelve dict de resultados."""
    primes_z = allowed_primes(z, 3, ())
    best = {}
    counts = {"Add-": 0, "Rem+": 0}
    for L in levels_L:
        for (r, rpf) in r_pfs:
            if r >= L or not in_Vplus(rpf, L):
                continue
            for (y, ypf) in ysubs:
                if math.gcd(y, r) != 1 or r * y >= L:
                    continue
                side = side_of(ypf, r, rpf, L)
                mu = mu_of(ypf)
                if side == "Add" and mu == -1:
                    counts["Add-"] += 1
                    tot, _ = flow_sum(y, ypf, r, rpf, L, primes_z,
                                      ("Add", -1), ("Add", 1))
                    key = ("RA", ypf[-1])
                    if key not in best or tot > best[key][0]:
                        best[key] = (tot, L, r, y)
                elif side == "Rem" and mu == 1:
                    counts["Rem+"] += 1
                    tot, _ = flow_sum(y, ypf, r, rpf, L, primes_z,
                                      ("Rem", 1), ("Rem", -1))
                    key = ("RR", ypf[-1])
                    if key not in best or tot > best[key][0]:
                        best[key] = (tot, L, r, y)
    print("  [%s] instancias Add^-: %d, Rem^+: %d"
          % (label, counts["Add-"], counts["Rem+"]))
    return best


def check_J4():
    print("J4 medicion extendida [M]: z=71 (global, DFS) y z=101 "
          "(min(y) >= 71 completo sobre la muestra):")
    ok = True
    # ---- z=101, min(y) >= 71: subconjuntos de {71,...,97}
    z = 101
    D = z ** 5
    spr = sample_sprimes(z, int(z ** (10.0 / 3.0)), 3, (), 8)
    levels = sorted({D // s for s in spr})
    # niveles dirigidos a la familia r=5, Y=89 y Y=97
    for Y in (89, 97):
        qmax = max(p for p in PRIMES if p < Y)
        s, L = find_level(D, max(Y ** 3, 125 * Y * qmax), 5 * Y ** 3,
                          int(z ** (10.0 / 3.0)))
        if s is not None and L not in levels:
            levels.append(L)
    primes_z = allowed_primes(z, 3, ())
    r_pfs = []
    for p in primes_z:
        r_pfs.append((p, (p,)))
    for i, p in enumerate(primes_z):
        for q in primes_z[i + 1:]:
            r_pfs.append((p * q, tuple(sorted((p, q), reverse=True))))
    big = [p for p in primes_z if p >= 71]
    ysubs = enum_subsets(big, D)
    best101 = measure_frame(z, levels, r_pfs, ysubs,
                            "z=101, min(y)>=71, %d niveles" % len(levels))
    for mn in (71, 73, 79, 83, 89, 97):
        for kind in ("RA", "RR"):
            if (kind, mn) in best101:
                (tot, L, r, y) = best101[(kind, mn)]
                print("    z=101 min(y)=%d: max %s = %.4f "
                      "(L=%d, r=%d, y=%d)"
                      % (mn, kind, float(tot), L, r, y))
    mx101 = max([float(v[0]) for v in best101.values()], default=0.0)
    print("    z=101: maximo global en la muestra = %.4f (< 1: la "
          "Mertens disponible, Mertens(97)-1/4 = %.4f, no alcanza)"
          % (mx101, float(sum((Fraction(1, p - 1) for p in PRIMES
                                if 5 <= p <= 97), Fraction(0)))
             - 0.25))
    # ---- z=71, global con presupuesto DFS
    z = 71
    D = z ** 5
    spr = sample_sprimes(z, int(z ** (10.0 / 3.0)), 3, (), 6)
    levels = sorted({D // s for s in spr})
    primes_z = allowed_primes(z, 3, ())
    r_pfs = []
    for p in primes_z:
        r_pfs.append((p, (p,)))
    for i, p in enumerate(primes_z):
        for q in primes_z[i + 1:]:
            r_pfs.append((p * q, tuple(sorted((p, q), reverse=True))))
    best71 = {}
    cnt71 = 0
    budget = [1200000]
    for L in levels:
        for (r, rpf) in r_pfs:
            if budget[0] <= 0:
                break
            if r >= L or not in_Vplus(rpf, L):
                continue
            bound = L // r
            pls = [p for p in primes_z if r % p != 0]
            ysubs = enum_subsets(pls, bound)
            budget[0] -= len(ysubs)
            for (y, ypf) in ysubs:
                if len(ypf) > 6 or r * y >= L:
                    continue
                side = side_of(ypf, r, rpf, L)
                mu = mu_of(ypf)
                if side == "Add" and mu == -1:
                    cnt71 += 1
                    tot, _ = flow_sum(y, ypf, r, rpf, L, primes_z,
                                      ("Add", -1), ("Add", 1))
                    if "RA" not in best71 or tot > best71["RA"][0]:
                        best71["RA"] = (tot, L, r, y)
                elif side == "Rem" and mu == 1:
                    tot, _ = flow_sum(y, ypf, r, rpf, L, primes_z,
                                      ("Rem", 1), ("Rem", -1))
                    if "RR" not in best71 or tot > best71["RR"][0]:
                        best71["RR"] = (tot, L, r, y)
    print("  [z=71] instancias Add^-: %d | nodos DFS restantes: %d"
          % (cnt71, budget[0]))
    for kind in ("RA", "RR"):
        if kind in best71:
            (tot, L, r, y) = best71[kind]
            print("    z=71: max %s = %.4f (L=%d, r=%d, y=%d)"
                  % (kind, float(tot), L, r, y))
    ok = mx101 < 1.0
    print("  => en z<=101 los maximos medidos quedan < 1 [M]; la "
          "violacion es de escala mayor (J2/J3). %s"
          % ("PASS" if ok else "FAIL"))
    return ok


# ------------------------------------------------------------------ main
def main():
    t0 = time.time()
    print("== DOC 116b: (R-FLOW-71+), cotas afiladas y contraejemplo ==")
    print("Solo stdlib. Predicados de los docs 109/113/114 por import.")
    print()
    allok = True
    allok &= check_J0()
    print()
    allok &= check_J1()
    print()
    allok &= check_J2()
    print()
    allok &= check_J3()
    print()
    allok &= check_J4()
    print()
    print("Tiempo total: %.1f s" % (time.time() - t0))
    print("PASS" if allok else "FAIL")
    return 0 if allok else 1


if __name__ == "__main__":
    sys.exit(main())
