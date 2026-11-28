#!/usr/bin/env python3
# 125b_erem_pos.py  (Doc 125b)
# (1) E_Rem = 0 en la capa alta: reproduccion de los docs 120a/122a
#     con |Rem| por instancia (K0); caza extendida z=53..71 separando
#     L/r <= 343 (teorema, doc 125b sec. 1) de L/r > 343 (residuo
#     (EREM-HI)) (K1); falsacion del Lema 2 (paridad del fallo)
#     fuera de la capa (K2).
# (2) Estrato 5 de (SLAB-V'): aritmetica del cierre via (POS) (K3)
#     y enumeracion extendida con el evaluador del doc 123a (K4).
# Solo stdlib; importa la maquinaria de los docs 109/120a/123a
# sin tocarla.
import importlib.util
import math
import sys
import time


def load(name, path):
    spec = importlib.util.spec_from_file_location(name, path)
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod


d120 = load("d120", "120a_banda_delta_hi.py")
m109 = load("m109", "109_ataque_h2_gram.py")
d123a = load("d123a", "123a_wit.py")

Level110 = d120.Level110
lados = d120.lados
band_instances = d120.band_instances
frame_levels = d120.frame_levels
visible_set = d120.visible_set
allowed_primes = m109.allowed_primes
q0_of_pf = d123a.q0_of_pf

TOL = 1e-12


def rem_info(lv, r):
    """(nRem, E_Rem, lista de (d', dpf)) desde lados() del doc 120a."""
    ld = lados(lv, r)
    elems = [(dp, v[1]) for dp, v in ld.items() if v[0] == "Rem"]
    s = sum(v[2] / v[3] for v in ld.values() if v[0] == "Rem")
    return len(elems), -s, elems


def merged_failure(rpf, dpf, L):
    """Primera posicion impar que falla en la lista fusionada de
    r*d': devuelve (m, j, p_m) con j = # primos de d' en el
    prefijo; None si r*d' esta en V+(L)."""
    ds = set(dpf)
    merged = sorted(set(rpf) | ds, reverse=True)
    prod = 1
    for i, p in enumerate(merged):
        prod *= p
        m = i + 1
        if m % 2 == 1 and p * p * prod >= L:
            j = sum(1 for q in merged[:m] if q in ds)
            return (m, j, p)
    return None


# ------------------------------------------------------------------ K0
def check_K0():
    print("K0 reproduccion docs 120a/122a (marcos estandar):")
    ok = True
    tot = viol = withrem = 0
    cov = res = 0
    deg_tot = deg_rem = deg_bad = 0
    worst = 0.0
    for (tag, z, D, y, cap, near_u) in d120.FRAMES:
        allow = allowed_primes(z, 3, ())
        deg_sample = []
        for L in frame_levels(z, D, y, cap, near_u):
            lv = Level110(z, 3, (), L)
            lnz = math.log(z)
            for (r, pf, phi, ur, Pi, sh, nhat, dif, E, dl) in \
                    band_instances(lv, L, lnz):
                if ur < 1.0:
                    continue
                q0 = q0_of_pf(allow, pf)
                if (nhat == 1) != (q0 ** 3 * r >= L):
                    deg_bad += 1
                if nhat == 1:
                    deg_sample.append((lv, r))
                    continue
                tot += 1
                n, er, _ = rem_info(lv, r)
                worst = max(worst, abs(er))
                if abs(er) > TOL:
                    viol += 1
                if n > 0:
                    withrem += 1
                if L <= 343 * r:
                    cov += 1
                else:
                    res += 1
        for (lv, r) in deg_sample[:25]:
            deg_tot += 1
            n, er, _ = rem_info(lv, r)
            if n > 0:
                deg_rem += 1
    print("  instancias no degeneradas de capa: %d | viol E_Rem!=0: %d"
          " | con Rem!=vacio: %d | max|E_Rem| = %.2e"
          % (tot, viol, withrem, worst))
    print("  split L/r<=343: %d | L/r>343: %d | degeneradas probadas:"
          " %d (con Rem: %d) | cruces nhat/q0 malos: %d"
          % (cov, res, deg_tot, deg_rem, deg_bad))
    ok &= viol == 0 and deg_rem == 0 and deg_bad == 0
    print("  => %s" % ("PASS" if ok else "FAIL"))
    return ok


# ------------------------------------------------------------------ K1
def check_K1():
    print("K1 caza extendida (residuo L/r > 343):")
    ok = True
    for z in (53, 61, 67, 71):
        levels = sorted(set(
            frame_levels(z, z ** 5, int(z ** (10.0 / 3.0)), 3, 2.0))
            | set(frame_levels(z, z ** 5, int(z ** (10.0 / 3.0)),
                               3, 3.0)))
        nres = ncov = vres = vcov = wres = 0
        worst = 0.0
        elems_all = []
        shape_bad = []
        sign_bad = 0
        cov_rem = 0
        for L in levels:
            lv = Level110(z, 3, (), L)
            lnz = math.log(z)
            inst = [(r, pf) for (r, pf, phi, ur, Pi, sh, nhat, dif,
                                 E, dl) in band_instances(lv, L, lnz)
                    if ur >= 1.0 and nhat > 1]
            res_i = [(r, pf) for (r, pf) in inst if L > 343 * r]
            cov_i = [(r, pf) for (r, pf) in inst if L <= 343 * r]
            if len(res_i) > 120:
                res_i = res_i[::len(res_i) // 120 + 1]
            todo = [("res", r) for r, pf in res_i] + \
                   [("cov", r) for r, pf in cov_i[:25]]
            pfmap = dict(inst)
            for (kind, r) in todo:
                n, er, elems = rem_info(lv, r)
                if kind == "res":
                    nres += 1
                    worst = max(worst, abs(er))
                    if er > TOL:
                        sign_bad += 1
                    if n > 0:
                        wres += 1
                        for (dp, dpf) in elems:
                            f = merged_failure(pfmap[r], dpf, L)
                            elems_all.append((z, L, r, dp, f))
                            if not (len(dpf) == 2 and f is not None
                                    and f[1] == 2 and f[0] >= 3):
                                shape_bad.append((z, L, r, dp, f))
                else:
                    ncov += 1
                    if abs(er) > TOL:
                        vcov += 1
                    if n > 0:
                        cov_rem += 1
        print("  z=%d: residuo probadas %d (con Rem %d, E_Rem>0: %d,"
              " max|E_Rem|=%.1e) | cubiertas probadas %d (viol %d,"
              " con Rem %d)"
              % (z, nres, wres, sign_bad, worst, ncov, vcov, cov_rem))
        for b in elems_all[:6]:
            print("    Rem realizado: z=%d L=%d r=%d d'=%d fallo=%s"
                  % b)
        print("    elementos Rem en residuo: %d | d' distintos: %s"
              " | formas no predichas: %s"
              % (len(elems_all),
                 sorted({e[3] for e in elems_all}), shape_bad[:3]))
        ok &= vcov == 0 and cov_rem == 0 and sign_bad == 0 \
            and not shape_bad
    print("  => %s" % ("PASS" if ok else "FAIL"))
    return ok


# ------------------------------------------------------------------ K2
def check_K2():
    print("K2 falsacion del Lema 2 fuera de la capa:")
    ok = True
    levels = []
    for L in frame_levels(37, 37 ** 5, int(37 ** (10.0 / 3.0)))[:3]:
        levels.append((37, L))
    for L in frame_levels(53, 53 ** 5, int(53 ** (10.0 / 3.0)),
                          4, 2.0):
        levels.append((53, L))
    hered_bad = 0
    hered_n = 0
    rem_found = {"banda_baja": 0, "capa": 0, "fuera": 0}
    parity_bad = []
    for (z, L) in levels:
        lv = Level110(z, 3, (), L)
        lnz = math.log(z)
        vis = sorted(visible_set(lv))
        for r in vis:
            hered_n += 1
            if r not in lv.Vset:
                hered_bad += 1
        info = {r: pf for (r, pf, mu, phi) in lv.enum}
        band = {}
        for (r, pf, phi, ur, Pi, sh, nhat, dif, E, dl) in \
                band_instances(lv, L, lnz):
            band[r] = ur
        probe = [(r, "capa" if band[r] >= 1.0 else "banda_baja")
                 for r in sorted(band)]
        otros = [r for r in vis if r not in band and r > 1]
        probe += [(r, "fuera") for r in otros[::29][:50]]
        for (r, kind) in probe:
            n, er, elems = rem_info(lv, r)
            for (dp, dpf) in elems:
                rem_found[kind] += 1
                f = merged_failure(info[r], dpf, L)
                if f is None:
                    parity_bad.append((z, L, r, dp, "sin fallo"))
                else:
                    (m, j, pm) = f
                    if m < 3 or j < 2 or j % 2 == 1:
                        parity_bad.append((z, L, r, dp, f))
    print("  herencia: %d visibles probados, r fuera de V+(L): %d"
          % (hered_n, hered_bad))
    print("  elementos Rem encontrados: %s | violaciones de"
          " paridad (m>=3, j par>=2): %d %s"
          % (rem_found, len(parity_bad), parity_bad[:4]))
    ok &= hered_bad == 0 and not parity_bad
    print("  => %s" % ("PASS" if ok else "FAIL"))
    return ok


# ------------------------------------------------------------------ K3
def primes_upto(n):
    s = [True] * (n + 1)
    s[0] = s[1] = False
    for i in range(2, int(n ** 0.5) + 1):
        if s[i]:
            for j in range(i * i, n + 1, i):
                s[j] = False
    return [i for i, v in enumerate(s) if v]


def check_K3():
    print("K3 aritmetica del cierre del estrato 5 via (POS):")
    ok = True
    b1 = 125 * 7 / 121
    b2 = 125 * 11 / 169
    b3 = 625 / 49
    print("  875/121 = %.6f (< 7.24: %s) | 1375/169 = %.6f"
          " (< 8.14: %s) | 625/49 = %.6f | 125/11 = %.6f"
          % (b1, b1 < 7.24, b2, b2 < 8.14, b3, 125 / 11))
    ok &= b1 < 7.24 and b2 < 8.14 and b3 < 12.76
    pr = primes_upto(100000)
    nxt = {pr[i]: pr[i + 1] for i in range(len(pr) - 1)}
    sup_bad = [p for p in pr if p >= 11 and p in nxt
               and 125.0 * p / nxt[p] ** 2 > b2 + 1e-12]
    eq = [p for p in pr if p >= 11 and p in nxt
          and abs(125.0 * p / nxt[p] ** 2 - b2) < 1e-12]
    mono_bad = [p for p in pr if 13 <= p <= 100000
                and 11 * p * p - 125 * p + 44 <= 0]
    print("  sup 125p/q^2 (pares consecutivos, p>=11): excede"
          " 1375/169 en %s | igualdad solo en p=11: %s |"
          " 11p^2-125p+44>0 para p>=13: %s"
          % (sup_bad[:3], eq == [11], not mono_bad))
    ok &= not sup_bad and eq == [11] and not mono_bad
    # barrido de candidatos c: libres de cuadrados, primos >= 5
    def sqfree_cs(B):
        out = []

        def rec(start, cur):
            if cur >= 5:
                out.append(cur)
            for p in pr:
                if p < start or cur * p >= B:
                    break
                if p >= 5:
                    rec(p + 1, cur * p)
        rec(5, 1)
        return out
    alive = []
    for c in sqfree_cs(125):
        if c in (1, 5, 7):
            continue
        fac = [p for p in pr if p <= c and c % p == 0]
        if 5 in fac:
            dead = c > b3  # k par: c < 625/49; k impar imposible
        else:
            # rama p_k=7: c < 875/121; rama p_k>=11: c < 1375/169
            dead = c > b1 and c > b2
        if not dead:
            alive.append(c)
    print("  candidatos c<125 fuera de {1,5,7} que sobreviven a"
          " las cotas: %s" % alive)
    ok &= not alive
    print("  => %s" % ("PASS" if ok else "FAIL"))
    return ok


# ------------------------------------------------------------------ K4
def check_K4():
    print("K4 enumeracion extendida del estrato 5 (evaluador 123a):")
    ok = True
    frames = [("A", 61, 61 ** 5, int(61 ** (10.0 / 3.0))),
              ("A", 67, 67 ** 5, int(67 ** (10.0 / 3.0))),
              ("B", 61, 61 ** 6, 61 ** 3)]
    for (tag, z, D, y) in frames:
        allow = allowed_primes(z, 3, ())
        enum = m109.enum_smooth(allow)
        enum_info = {d: (pf, phi) for (d, pf, mu, phi) in enum}
        pifull = d123a.pi_full(allow)
        cand = {B: d123a.cand_cofactors(allow, B)
                for B in (125, 343, 1331)}
        agg, nl, sk = d123a.lean_frame(z, D, y, cand, enum,
                                       enum_info, allow, pifull)
        a = agg.get(5)
        if a is None:
            print("  [%s z=%d] estrato 5: sin degenerados en banda"
                  " (niveles %d, saltados %d)" % (tag, z, nl, sk))
            continue
        viol5 = [v for v in a["viol"]]
        good = a["c"] <= {1, 5, 7} and not viol5 and a["nmax"] <= 3
        ok &= good
        print("  [%s z=%d] estrato 5: c=%s Nmax=%d viol=%s"
              " (niveles %d, saltados %d)"
              % (tag, z, sorted(a["c"]), a["nmax"], viol5[:4],
                 nl, sk))
    print("  => %s" % ("PASS" if ok else "FAIL"))
    return ok


def main():
    t0 = time.time()
    ok = True
    for chk in (check_K0, check_K1, check_K2, check_K3, check_K4):
        try:
            ok &= chk()
        except Exception as e:
            ok = False
            print("  EXCEPCION en %s: %r" % (chk.__name__, e))
    print("Tiempo total: %.1f s" % (time.time() - t0))
    print("PASS" if ok else "FAIL")
    return 0 if ok else 1


if __name__ == "__main__":
    sys.exit(main())
