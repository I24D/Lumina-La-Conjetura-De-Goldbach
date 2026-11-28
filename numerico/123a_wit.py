#!/usr/bin/env python3
# 123a_wit.py  (Doc 123a)
# (WIT) [P, doc 122b sec. 4]: en la subcapa degenerada visible
# (estratos q0 en {5,7,11}), todo cofactor c con (c,r)=1 y
# cr en V+(L) cumple c en {1, q0, q0+} y N(r) <= 3.
#   K0  validacion: reproduce los conjuntos de c del doc 122b
#       (marcos A z=17..53, B z=29) con la maquinaria real.
#   K1  barrido extendido con evaluador ligero (validado contra
#       K0 a z=53): marcos A z=61,67 y B z=61, niveles de la
#       muestra de la serie; registra c que ocurren y violaciones.
#   K2  barrido puntual dirigido (predicados exactos de la serie,
#       sin enumerar el universo): candidatos r de estratos
#       5/7/11, z en {127,151,191,227,251}, marcos A y B(z=127);
#       busca c fuera de {1,q0,q0+} o N(r)>3 y verifica el
#       contraejemplo predicado a predicado.
#   K3  robustez: para el contraejemplo, barrido de niveles
#       (sprimos s en [z,8000], marco A): cuantos niveles lo
#       reproducen.
# Solo stdlib; importa la maquinaria de los docs 109/118b1/122b.
import importlib.util
import math
import sys
import time


def load(name, path):
    spec = importlib.util.spec_from_file_location(name, path)
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod


m109 = load("m109", "109_ataque_h2_gram.py")
b1 = load("b1", "118b1_degenerada.py")
b122 = load("b122", "122b_slab_v_prima.py")

allowed_primes = m109.allowed_primes
enum_smooth = m109.enum_smooth
in_Vplus = m109.in_Vplus
Pi_of = m109.Pi_of
sub_divisors = m109.sub_divisors
sample_sprimes = m109.sample_sprimes
primes_upto = m109.primes_upto

STRATA = (5, 7, 11)


def q0_of_pf(allow, pf):
    s = set(pf)
    for p in allow:
        if p not in s:
            return p
    return None


def q0plus_of(allow, q0):
    i = allow.index(q0)
    return allow[i + 1] if i + 1 < len(allow) else None


# ---------------------------------------------------------------- K0
EXPECTED = {
    ("A", 17): {5: set(), 7: set(), 11: set()},
    ("A", 29): {5: {1, 5, 7}, 7: set(), 11: set()},
    ("A", 37): {5: {1, 5, 7}, 7: set(), 11: set()},
    ("A", 53): {5: {1, 5, 7}, 7: {1, 7, 11}, 11: set()},
    ("B", 29): {5: {1}, 7: set(), 11: set()},
}


def check_K0():
    ok = True
    for (tag, z, D, y) in b122.FRAMES:
        occ = {q0: set() for q0 in STRATA}
        nmax = 0
        for lv in b1.frame_levels(z, D, y):
            dat = b122.level_data(lv, z)
            for q0 in STRATA:
                occ[q0] |= {c for c, v in dat["Udeg"][q0].items()
                            if v > 0}
            if dat["mult"]:
                nmax = max(nmax, max(dat["mult"].values()))
        exp = EXPECTED[(tag, z)]
        match = all(occ[q0] == exp[q0] for q0 in STRATA)
        ok &= match and nmax <= 3
        print("K0 [%s z=%d] c: %s | Nmax=%d | coincide con 122b: %s"
              % (tag, z, {q0: sorted(occ[q0]) for q0 in STRATA},
                 nmax, match))
    print("K0 => %s" % ("OK" if ok else "FALLO"))
    return ok


# ---------------------------------------------------------------- K1
def cand_cofactors(allow, B):
    """c libres de cuadrados, 1 <= c < B, primos en allow:
    dict c -> tupla creciente de primos."""
    out = {1: ()}

    def rec(start, cur, pf):
        for i in range(start, len(allow)):
            p = allow[i]
            if cur * p >= B:
                break
            out[cur * p] = pf + (p,)
            rec(i + 1, cur * p, pf + (p,))

    rec(0, 1, ())
    return out


def pi_full(allow):
    t = 1.0
    for p in allow:
        t *= (1.0 - 1.0 / (p - 1))
    return t


def pi_r(pifull, pf):
    t = pifull
    for p in pf:
        t /= (1.0 - 1.0 / (p - 1))
    return t


def lean_level(enum, enum_info, allow, pifull, z, L, cand_by_B):
    """Recomputa, sin Level110, los r degenerados visibles en banda
    y sus testigos por enumeracion directa de cofactores.
    Devuelve (res, |V|); res[q0] = {'c': set, 'nmax': int,
    'viol': [(L, r, c o 'N=..')]} para estratos 5/7/11 y 'otros'."""
    V = [(d, pf) for (d, pf, mu, phi) in enum
         if d < L and in_Vplus(pf, L)]
    if len(V) > 2500000:
        return None, len(V)
    vis = set()
    for (d, pf) in V:
        for (r, rmu, rphi) in sub_divisors(pf):
            vis.add(r)
    logz = math.log(z)
    res = {}
    for r in vis:
        pf = enum_info[r][0]
        if pi_r(pifull, pf) >= 0.5:
            continue
        if math.log(L / r) / logz >= 1.5:
            continue
        q0 = q0_of_pf(allow, pf)
        if not (q0 is None or q0 ** 3 * r >= L):
            continue
        if q0 not in STRATA:
            res.setdefault("otros", {"c": set(), "nmax": 0,
                                     "viol": []})
            res["otros"]["nmax"] += 1  # conteo de r, no de testigos
            continue
        B = q0 ** 3
        cand = cand_by_B[B]
        occ = []
        for c, cpf in cand.items():
            if math.gcd(c, r) != 1 or c * r >= L:
                continue
            merged = tuple(sorted(set(pf) | set(cpf), reverse=True))
            if in_Vplus(merged, L):
                occ.append(c)
        e = res.setdefault(q0, {"c": set(), "nmax": 0, "viol": []})
        e["c"] |= set(occ)
        e["nmax"] = max(e["nmax"], len(occ))
        good = {1, q0, q0plus_of(allow, q0)}
        for c in occ:
            if c not in good:
                e["viol"].append((L, r, c))
        if len(occ) > 3:
            e["viol"].append((L, r, "N=%d" % len(occ)))
    return res, len(V)


def lean_frame(z, D, y, cand_by_B, enum, enum_info, allow, pifull):
    spr = sample_sprimes(z, y, 3, ())
    levels = sorted({D // s for s in spr})
    agg = {}
    skipped = 0
    for L in levels:
        res, nv = lean_level(enum, enum_info, allow, pifull, z, L,
                             cand_by_B)
        if res is None:
            skipped += 1
            continue
        for k, e in res.items():
            a = agg.setdefault(k, {"c": set(), "nmax": 0, "viol": []})
            a["c"] |= e["c"]
            a["nmax"] = max(a["nmax"], e["nmax"])
            a["viol"] += e["viol"]
    return agg, len(levels), skipped


def check_K1():
    ok = True
    # validacion del evaluador ligero contra la maquinaria real (z=53)
    z = 53
    allow = allowed_primes(z, 3, ())
    enum = enum_smooth(allow)
    enum_info = {d: (pf, phi) for (d, pf, mu, phi) in enum}
    pifull = pi_full(allow)
    # formula de Pi contra Pi_of en muestra
    import random
    random.seed(12345)
    worst = 0.0
    for (d, pf, mu, phi) in random.sample(enum, 200):
        worst = max(worst, abs(pi_r(pifull, pf)
                               - Pi_of(pf, z, 3, ())) / Pi_of(pf, z, 3, ()))
    ok &= worst < 1e-12
    cand_by_B = {B: cand_cofactors(allow, B) for B in (125, 343, 1331)}
    agg, nl, sk = lean_frame(z, z ** 5, int(z ** (10.0 / 3.0)),
                             cand_by_B, enum, enum_info, allow, pifull)
    got = {q0: agg.get(q0, {"c": set()})["c"] for q0 in STRATA}
    exp = EXPECTED[("A", 53)]
    match = all(got[q0] == exp[q0] for q0 in STRATA)
    ok &= match
    print("K1 validacion: err Pi=%.1e | ligero z=53 c=%s | "
          "coincide con real: %s" % (worst,
                                     {q0: sorted(got[q0])
                                      for q0 in STRATA}, match))
    # barrido extendido
    frames = [("A", 61, 61 ** 5, int(61 ** (10.0 / 3.0))),
              ("A", 67, 67 ** 5, int(67 ** (10.0 / 3.0))),
              ("B", 61, 61 ** 6, 61 ** 3)]
    for (tag, z, D, y) in frames:
        allow = allowed_primes(z, 3, ())
        enum = enum_smooth(allow)
        enum_info = {d: (pf, phi) for (d, pf, mu, phi) in enum}
        pifull = pi_full(allow)
        cand_by_B = {B: cand_cofactors(allow, B)
                     for B in (125, 343, 1331)}
        agg, nl, sk = lean_frame(z, D, y, cand_by_B, enum,
                                 enum_info, allow, pifull)
        for q0 in STRATA:
            a = agg.get(q0)
            if a is None:
                print("K1 [%s z=%d] q0=%d: sin degenerados en banda"
                      % (tag, z, q0))
            else:
                print("K1 [%s z=%d] q0=%d: c=%s Nmax=%d viol=%s"
                      % (tag, z, q0, sorted(a["c"]), a["nmax"],
                         a["viol"][:4]))
        if sk:
            print("K1 [%s z=%d] niveles saltados por |V|: %d"
                  % (tag, z, sk))
    print("K1 => %s" % ("OK" if ok else "FALLO"))
    return ok


# ---------------------------------------------------------------- K2
CAND_R = [
    ((5, 7, 13, 17, 19, 23), 11),
    ((5, 7, 13, 17, 19), 11),
    ((5, 7, 13, 17, 19, 23, 29), 11),
    ((5, 7, 101, 103), 11),
    ((5, 7, 211, 223), 11),
    ((5, 11, 13, 17, 19, 23), 7),
    ((7, 11, 13, 17, 19, 23, 29), 5),
]


def prod(pf):
    t = 1
    for p in pf:
        t *= p
    return t


def witness_set(allow, cand, r, pf_r, L):
    occ = []
    for c, cpf in cand.items():
        if math.gcd(c, r) != 1 or c * r >= L:
            continue
        merged = tuple(sorted(set(pf_r) | set(cpf), reverse=True))
        if in_Vplus(merged, L):
            occ.append(c)
    return sorted(occ)


def nathanson_margins(pf_desc, L):
    """Para m impar: p_m^2 * prod_{j<=m} p_j frente a L."""
    out = []
    pr = 1
    for i, p in enumerate(pf_desc):
        m = i + 1
        pr *= p
        if m % 2 == 1:
            out.append((m, p, p * p * pr, L))
    return out


def verify_instance(z, tag, D, s, L, pf_r):
    """Verificacion completa con los predicados de la serie.
    Devuelve (q0, occ, bad) o None si el candidato no califica."""
    allow = allowed_primes(z, 3, ())
    r = prod(pf_r)
    q0 = q0_of_pf(allow, pf_r)
    if q0 is None or q0 not in STRATA:
        return None
    if not (r < L and q0 ** 3 * r >= L):
        return None
    if math.log(L / r) / math.log(z) >= 1.5:
        return None
    if Pi_of(pf_r, z, 3, ()) >= 0.5:
        return None
    pf_desc = tuple(sorted(pf_r, reverse=True))
    if not in_Vplus(pf_desc, L):
        return None
    cand = cand_cofactors(allow, q0 ** 3)
    occ = witness_set(allow, cand, r, pf_r, L)
    good = {1, q0, q0plus_of(allow, q0)}
    bad = [c for c in occ if c not in good]
    if not bad and len(occ) <= 3:
        return None
    # --- impresion de la verificacion completa ---
    print("K2 CONTRAEJEMPLO z=%d marco %s: D=%d s=%d L=%d"
          % (z, tag, D, s, L))
    print("   r = %d, primos(r) = %s" % (r, pf_r))
    print("   q0(r) = %d (allow[0..]=%s...), q0+ = %s"
          % (q0, allow[:4], q0plus_of(allow, q0)))
    print("   degenerado: q0^3*r = %d >= L = %d : %s"
          % (q0 ** 3 * r, L, q0 ** 3 * r >= L))
    print("   banda: Pi(r) = %.6f < 1/2 : %s ; u_r = %.4f < 3/2 : %s"
          % (Pi_of(pf_r, z, 3, ()),
             Pi_of(pf_r, z, 3, ()) < 0.5,
             math.log(L / r) / math.log(z),
             math.log(L / r) / math.log(z) < 1.5))
    print("   visible: r en V+(L) : %s ; margenes Nathanson de r: %s"
          % (in_Vplus(pf_desc, L),
             [(m, v, "%s<L" % (v < L)) for (m, p, v, _) in
              nathanson_margins(pf_desc, L)]))
    print("   testigos de r: c en %s ; N(r) = %d" % (occ, len(occ)))
    for c in bad:
        cpf = cand[c]
        merged = tuple(sorted(set(pf_r) | set(cpf), reverse=True))
        print("   c = %d (primos %s) FUERA de {1, q0, q0+}: "
              "d = c*r = %d en V+(L): %s"
              % (c, cpf, c * r, in_Vplus(merged, L)))
        print("      margenes Nathanson de d: %s"
              % [(m, p, v, v < L) for (m, p, v, _) in
                 nathanson_margins(merged, L)])
    return (q0, occ, bad)


def check_K2():
    found = {}
    for z in (127, 151, 191, 227, 251):
        ymax = int(z ** (10.0 / 3.0))
        sprA = [p for p in primes_upto(min(ymax, 20000))
                if p >= z and math.gcd(p, 3) == 1]
        for (pf_r, want) in CAND_R:
            if max(pf_r) >= z:
                continue
            allow = allowed_primes(z, 3, ())
            if q0_of_pf(allow, pf_r) != want:
                continue
            if want in found:
                continue
            D = z ** 5
            for s in sprA:
                L = D // s
                got = verify_instance(z, "A", D, s, L, pf_r)
                if got:
                    found[want] = (z, "A", D, s, L, pf_r, got)
                    break
    # marco B a z=127 para el candidato principal de estrato 11
    if 11 not in found:
        z = 127
        D = z ** 6
        sprB = [p for p in primes_upto(30000)
                if p >= z and math.gcd(p, 3) == 1]
        for s in sprB:
            got = verify_instance(z, "B", D, s, D // s,
                                  (5, 7, 13, 17, 19, 23))
            if got:
                found[11] = (z, "B", D, s, D // s,
                             (5, 7, 13, 17, 19, 23), got)
                break
    for q0 in STRATA:
        if q0 in found:
            print("K2 estrato %d: contraejemplo registrado (z=%d, "
                  "marco %s)" % (q0, found[q0][0], found[q0][1]))
        else:
            print("K2 estrato %d: sin contraejemplo en el barrido "
                  "dirigido" % q0)
    print("K2 => %s" % ("CONTRAEJEMPLO" if found else "sin hallazgo"))
    return found


# ---------------------------------------------------------------- K3
def check_K3(found):
    if 11 not in found:
        print("K3 => no aplica (sin contraejemplo de estrato 11)")
        return True
    (z, tag, D0, s0, L0, pf_r, got) = found[11]
    q0, occ0, bad0 = got
    r = prod(pf_r)
    allow = allowed_primes(z, 3, ())
    cand = cand_cofactors(allow, q0 ** 3)
    good = {1, q0, q0plus_of(allow, q0)}
    D = z ** 5
    spr = [p for p in primes_upto(8000)
           if p >= z and math.gcd(p, 3) == 1]
    hits = []
    for s in spr:
        L = D // s
        if not (r < L and q0 ** 3 * r >= L):
            continue
        if math.log(L / r) / math.log(z) >= 1.5:
            continue
        if Pi_of(pf_r, z, 3, ()) >= 0.5:
            continue
        if not in_Vplus(tuple(sorted(pf_r, reverse=True)), L):
            continue
        occ = witness_set(allow, cand, r, pf_r, L)
        bad = [c for c in occ if c not in good]
        if bad or len(occ) > 3:
            hits.append((s, L, occ))
    print("K3 robustez (z=%d, r=%d, marco A): %d niveles con "
          "violacion de (WIT) entre s en [%d, %d]; ejemplos: %s"
          % (z, r, len(hits), z, 8000,
             [(s, occ) for (s, L, occ) in hits[:5]]))
    print("K3 => %s" % ("OK" if hits else "sin reproduccion"))
    return bool(hits)


def main():
    t0 = time.time()
    print("== DOC 123a: (WIT) — intento de prueba y busqueda "
          "adversarial ==")
    ok0 = check_K0()
    ok1 = check_K1()
    found = check_K2()
    ok3 = check_K3(found)
    ok = ok0 and ok1 and bool(found) and ok3
    print()
    print("Tiempo total: %.1f s" % (time.time() - t0))
    print("PASS" if ok else "FAIL")
    return 0 if ok else 1


if __name__ == "__main__":
    sys.exit(main())
