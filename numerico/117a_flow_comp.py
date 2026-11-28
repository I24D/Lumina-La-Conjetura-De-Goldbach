#!/usr/bin/env python3
# 117a_flow_comp.py  (Doc 117a)
# (FLOW-COMP): masa mala vs masa buena del flujo de apareamiento,
# por conductor r.  E_Add(r) = -sum_{y in Add^-} |a(y)|(1 - R_A(y)),
# luego (FLOW-COMP)_Add  <=>  BAD(r) <= GOOD(r) con
#   BAD(r)  = sum |a(y)|(R_A(y)-1)^+ ,  GOOD(r) = sum |a(y)|(1-R_A(y))^+,
# y GOOD - BAD = -E_Add(r) (suma directa, sin calcular R_A).
# K0  sumas de Mertens exactas y margenes uniformes.
# K1  marcos A/B del doc 114 (instancias()): E_Add, E_Rem por r,
#     masa cabeza min(y)<=37 de Add^-; BAD = 0 identico (z <= 53 < 73).
# K2  PRUEBA DECISIVA: marco del contraejemplo del doc 116b
#     (z = 331, r = 5): enumeracion DFS de {y : 5y in V+(L)} en tres
#     niveles de L_A (cruce Y=251, control Y=211, cruce Y=283);
#     BAD exacto (fracciones) sobre min(y) >= 73; lado Rem en el
#     nivel del cruce.
# K3  estres z = 101: niveles de L_A, r visibles con omega(r) <= 2.
# K4  termino malo exacto del contraejemplo Rem del doc 116b (J3).
# Uso: python3 117a_flow_comp.py | tee salida_doc117a_local.txt
import importlib.util
import sys
import time
from fractions import Fraction

sys.setrecursionlimit(1000000)


def load(name, path):
    spec = importlib.util.spec_from_file_location(name, path)
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod


d116b = load("doc116b", "116b_flujo_afinado.py")
d109 = d116b.d109
d114 = d116b.d114
allowed_primes = d109.allowed_primes
in_Vhat = d116b.in_Vhat
flow_sum = d116b.flow_sum
PRIMES = d116b.PRIMES


def mert(n):
    return sum((Fraction(1, p - 1) for p in PRIMES if 5 <= p <= n),
               Fraction(0))


# ------------------------------------------------------------------ K0
def check_K0():
    print("K0 Mertens exacto y margenes:")
    m37, m67, m251 = mert(37), mert(67), mert(251)
    print("  Mertens(37) = %s = %.6f -> margen uniforme min(y)<=37: "
          "1-R_A >= %.6f" % (m37, float(m37), float(1 - m37)))
    print("  Mertens(67) = %.6f (<1 [T] doc 114); Mertens(251) = %.6f"
          % (float(m67), float(m251)))
    ok = m37 < 1 and m67 < 1 and m251 > 1
    print("  => %s" % ("PASS" if ok else "FAIL"))
    return ok


# ------------------------------------------------------------------ K1
def check_K1():
    print("K1 marcos A/B (instancias del doc 114): por (nivel, r):")
    inst = d114.instancias()
    worst_add = -1e300
    min_rem_sum = 1e300
    min_good_add = 1e300
    head_min_ratio = 1e300
    n = 0
    for (tag, z, a, r, ld) in inst:
        eadd = 0.0
        erem = 0.0
        head = 0.0
        tot_minus = 0.0
        for (dp, (side, dpf, dmu, dphi)) in ld.items():
            if side == "Add":
                eadd += dmu / dphi
                if dmu == -1:
                    tot_minus += 1.0 / dphi
                    if dpf and dpf[-1] <= 37:
                        head += 1.0 / dphi
            elif side == "Rem":
                erem += dmu / dphi
        n += 1
        worst_add = max(worst_add, eadd)
        # E_Rem del doc 114 = -sum_{Rem} a(x): (FLOW-COMP)_Rem pide
        # sum_{Rem} a(x) >= 0 (con BAD_R = 0 en estos marcos).
        min_rem_sum = min(min_rem_sum, erem)
        # BAD = 0 en estos marcos (z <= 53: min(y) <= 47 < 73, Mertens)
        min_good_add = min(min_good_add, -eadd)
        if tot_minus > 0:
            head_min_ratio = min(head_min_ratio, head / tot_minus)
    print("  %d instancias (nivel, r): max E_Add = %.3e, max E_Rem"
          "(doc 114) = %.3e" % (n, worst_add, -min_rem_sum))
    print("  min_r GOOD_Add = -E_Add = %.3e ; min_r GOOD_Rem = "
          "sum_{Rem} a = %.3e" % (min_good_add, min_rem_sum))
    print("  BAD = 0 identico en todo el marco (min(y) <= 47 < 73); "
          "fraccion minima de masa Add^- con min(y)<=37: %.4f"
          % head_min_ratio)
    ok = worst_add <= 1e-9 and min_rem_sum >= -1e-9
    print("  => (FLOW-COMP) vale en los %d con BAD=0. %s"
          % (n, "PASS" if ok else "FAIL"))
    return ok


# ------------------------------------------------------- enumeradores DFS
def enum_y_r5(L, primes_desc):
    """Genera (y, ypf, mu, phi) con 5y en V+(L): poda por la condicion
    de posicion impar de in_Vplus(ypf, L) y por 5y < L; la condicion
    propia del 5 (posicion omega(y)+1) se chequea al generar."""
    stack = [(0, 1, (), 1)]
    while stack:
        idx, prod, pf, phi = stack.pop()
        if 5 * prod < L and (len(pf) % 2 == 1 or 125 * prod < L):
            yield prod, pf, (1 if len(pf) % 2 == 0 else -1), phi
        m = len(pf) + 1
        for j in range(idx, len(primes_desc)):
            q = primes_desc[j]
            nprod = prod * q
            if 5 * nprod >= L:
                continue
            if m % 2 == 1 and q * q * nprod >= L:
                continue
            stack.append((j + 1, nprod, pf + (q,), phi * (q - 1)))


def enum_x_vhat_r5(L, primes_desc):
    """Genera (x, xpf, mu, phi) con x en Vhat^{(5)} (in_Vhat) y x < L/5:
    poda por la condicion impar r*p^2*prod < L con r = 5."""
    stack = [(0, 1, (), 1)]
    while stack:
        idx, prod, pf, phi = stack.pop()
        yield prod, pf, (1 if len(pf) % 2 == 0 else -1), phi
        m = len(pf) + 1
        for j in range(idx, len(primes_desc)):
            q = primes_desc[j]
            nprod = prod * q
            if nprod * 5 >= L:
                continue
            if m % 2 == 1 and 5 * q * q * nprod >= L:
                continue
            stack.append((j + 1, nprod, pf + (q,), phi * (q - 1)))


# ------------------------------------------------------------------ K2
def eval_331_add(L, primes_desc, primes_z):
    eadd = 0.0
    head37 = 0.0
    tail73 = 0.0
    bad = Fraction(0)
    max_ra = 0.0
    argmax = None
    ra251 = None
    nbad = 0
    for (y, ypf, mu, phi) in enum_y_r5(L, primes_desc):
        if in_Vhat(ypf, L, 5):
            continue
        eadd += mu / phi
        if mu == -1 and ypf:
            mn = ypf[-1]
            if mn <= 37:
                head37 += 1.0 / phi
            if mn >= 73:
                tail73 += 1.0 / phi
                ra, _ = flow_sum(y, ypf, 5, (5,), L, primes_z,
                                 ("Add", -1), ("Add", 1))
                if y == 251:
                    ra251 = ra
                if ra > max_ra:
                    max_ra, argmax = float(ra), y
                if ra > 1:
                    bad += Fraction(1, phi) * (ra - 1)
                    nbad += 1
    good = float(bad) - eadd
    return dict(eadd=eadd, head37=head37, tail73=tail73, bad=bad,
                good=good, max_ra=max_ra, argmax=argmax, ra251=ra251,
                nbad=nbad)


def eval_331_rem(L, primes_desc, primes_z):
    erem = 0.0
    bad = Fraction(0)
    max_rr = 0.0
    argmax = None
    for (x, xpf, mu, phi) in enum_x_vhat_r5(L, primes_desc):
        merged = xpf + (5,)
        if d109.in_Vplus(merged, L):
            continue
        erem += mu / phi
        if mu == 1 and xpf and xpf[-1] >= 73:
            rr, _ = flow_sum(x, xpf, 5, (5,), L, primes_z,
                             ("Rem", 1), ("Rem", -1))
            if rr > max_rr:
                max_rr, argmax = float(rr), x
            if rr > 1:
                bad += Fraction(1, phi) * (rr - 1)
    good = float(bad) - erem
    return dict(erem=erem, bad=bad, good=good, max_rr=max_rr,
                argmax=argmax)


def check_K2():
    print("K2 PRUEBA DECISIVA z=331, r=5 (marco del contraejemplo "
          "116b):")
    primes_z = allowed_primes(331, 3, ())
    primes_desc = [p for p in primes_z if p != 5][::-1]
    ok = True
    for (L, tag) in ((46948396, "control Y=211"),
                     (79051268, "CRUCE Y=251"),
                     (113296524, "cruce Y=283")):
        t0 = time.time()
        res = eval_331_add(L, primes_desc, primes_z)
        print("  L=%d [%s]: E_Add(5) = %+.6f | BAD = %.3e (%d y malos)"
              " | GOOD = %.6f | BAD/GOOD = %.2e"
              % (L, tag, res["eadd"], float(res["bad"]), res["nbad"],
                 res["good"], float(res["bad"]) / res["good"]))
        print("    max R_A = %.5f en y=%s | R_A(251) = %s | masa "
              "cabeza37 = %.4f | masa cola>=73 = %.3e | %.1f s"
              % (res["max_ra"], res["argmax"],
                 ("%.5f" % float(res["ra251"])) if res["ra251"]
                 is not None else "n/a",
                 res["head37"], res["tail73"], time.time() - t0))
        if res["bad"] > 0:
            slack = float(1 - mert(37)) * res["head37"] / float(
                res["bad"])
            print("    holgura (FC-2): (1-Mertens(37))*cabeza37 / BAD "
                  "= %.2f" % slack)
        ok &= res["eadd"] < 0 and res["good"] > float(res["bad"])
        if tag.startswith("CRUCE"):
            ok &= res["ra251"] is not None \
                and abs(float(res["ra251"]) - 1.00220) < 1e-4
            ok &= res["bad"] > 0
    res = eval_331_rem(79051268, primes_desc, primes_z)
    # E_Rem del doc 114 = -sum_{Rem} a(x) = BAD_R - GOOD_R.
    print("  L=79051268 lado Rem: E_Rem(doc 114) = %+.6f | BAD_R = "
          "%.3e | GOOD_R = %.6f | max R_R(min>=73) = %.5f en y=%s"
          % (-res["erem"], float(res["bad"]),
             res["erem"] + float(res["bad"]),
             res["max_rr"], res["argmax"]))
    ok &= res["erem"] >= 0
    print("  => %s" % ("PASS" if ok else "FAIL"))
    return ok


# ------------------------------------------------------------------ K3
def enum_V(L, primes_desc):
    """Todos los (d, pf, mu, phi) en V+(L), DFS con poda de in_Vplus."""
    out = []
    stack = [(0, 1, (), 1)]
    while stack:
        idx, prod, pf, phi = stack.pop()
        out.append((prod, pf, (1 if len(pf) % 2 == 0 else -1), phi))
        m = len(pf) + 1
        for j in range(idx, len(primes_desc)):
            q = primes_desc[j]
            nprod = prod * q
            if nprod >= L:
                continue
            if m % 2 == 1 and q * q * nprod >= L:
                continue
            stack.append((j + 1, nprod, pf + (q,), phi * (q - 1)))
    return out


def check_K3():
    print("K3 estres z=101, r visibles con omega(r)<=2:")
    z = 101
    D = z ** 5
    primes_z = allowed_primes(z, 3, ())
    primes_desc = primes_z[::-1]
    ok = True
    for s in (100003, 10007, 1013):
        L = D // s
        V = enum_V(L, primes_desc)
        if len(V) > 2500000:
            print("  s=%d L=%d: |V| = %d demasiado grande, se omite"
                  % (s, L, len(V)))
            continue
        cand = [(d, pf, mu, phi) for (d, pf, mu, phi) in V
                if d > 1 and len(pf) <= 2]
        step = max(1, len(cand) // 40)
        cand = cand[::step][:40]
        worst = -1e300
        min_ratio = 1e300
        max_ra = 0.0
        for (r, rpf, rmu, rphi) in cand:
            rset = set(rpf)
            eadd = 0.0
            bad = Fraction(0)
            for (d, pf, mu, phi) in V:
                if d % r != 0:
                    continue
                ypf = tuple(p for p in pf if p not in rset)
                if in_Vhat(ypf, L, r):
                    continue
                ymu = mu * rmu
                yphi = phi // rphi
                eadd += ymu / yphi
                if ymu == -1 and ypf and ypf[-1] >= 73:
                    ra, _ = flow_sum(d // r, ypf, r, rpf, L, primes_z,
                                     ("Add", -1), ("Add", 1))
                    max_ra = max(max_ra, float(ra))
                    if ra > 1:
                        bad += Fraction(1, yphi) * (ra - 1)
            good = float(bad) - eadd
            worst = max(worst, eadd)
            if bad > 0:
                min_ratio = min(min_ratio, good / float(bad))
        print("  s=%d L=%d |V|=%d, %d conductores: max E_Add = %.3e | "
              "max R_A(min>=73) = %.4f | BAD>0 en %s"
              % (s, L, len(V), len(cand), worst, max_ra,
                 "algun r" if min_ratio < 1e300 else "ningun r"))
        ok &= worst <= 1e-9
    print("  => %s" % ("PASS" if ok else "FAIL"))
    return ok


# ------------------------------------------------------------------ K4
def check_K4():
    print("K4 termino malo exacto del contraejemplo Rem (116b J3):")
    y = 20531891  # 5009 * 4099, y en Rem^+, R_R = Mertens(251)
    phi_y = 5008 * 4098
    excess = mert(251) - 1
    term = excess / phi_y
    print("  |a(y)|(R_R-1) = (Mertens(251)-1)/phi(y) = %s = %.3e"
          % (term, float(term)))
    ok = Fraction(0) < term < Fraction(1, 10000000)
    print("  => el unico termino malo conocido en Rem pesa < 1e-7. %s"
          % ("PASS" if ok else "FAIL"))
    return ok


# ------------------------------------------------------------------ main
def main():
    t0 = time.time()
    print("== DOC 117a: (FLOW-COMP), masa mala contra masa buena ==")
    print("Solo stdlib; importa los docs 109/114/116b sin tocarlos.")
    print()
    allok = True
    allok &= check_K0()
    print()
    allok &= check_K1()
    print()
    allok &= check_K2()
    print()
    allok &= check_K3()
    print()
    allok &= check_K4()
    print()
    print("Tiempo total: %.1f s" % (time.time() - t0))
    print("PASS" if allok else "FAIL")
    return 0 if allok else 1


if __name__ == "__main__":
    sys.exit(main())
