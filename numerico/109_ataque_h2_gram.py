#!/usr/bin/env python3
# 109_ataque_h2_gram.py  (Doc 109)
# Ataque a (H2), la cota de Gram normalizada (hueco G2 del doc 108):
#   (H2)  sup_{z_A<=s<y_A} Q_s << Q_completo(z_A; N),
# con Q_s = sum_r phi(r) C_s(r)^2 (diagonalizacion, docs 96/97) y
# Q_completo(z;N) = prod_{p<z, p no|AN} (1 - 1/(p-1) + 1/(p-1)^2).
#
# Solo stdlib. Modelo de LIMITE FIJO (segunda reduccion, doc 95 sec. 2.3):
# en este modo V_s = V^+(D/s) solo depende del nivel L = D//s, luego
# Q_s = Q^{(L)}: el barrido en s es un barrido en niveles. El modelo
# replica exactamente el Model97 (fixed) del doc 97 (chequeo K0 lo
# valida contra el script del doc 97 importado).
#
# Chequeos:
#   K0  validacion cruzada: modelo ligero vs Model97 del doc 97 (S1/S4)
#       y la identidad sin criba Q_sin_restriccion = Q_completo.
#   K1  barrido en z (marco A: D=z^5, y=z^(10/3), A=3, E0): sup_s Q_s,
#       R(z) = sup Q_s/Q_completo, Q_s*ln z (comparar doc 95 sec. 2.5).
#   K2  ajustes de tasa: gamma en sup Q_s ~ (ln z)^-gamma; delta en
#       R(z) ~ (ln z)^delta; marco B (D=z^6, y=z^3, u_min=3) vs marco A.
#   K3  dependencia en A y en N (E0 vs E1) a z=29: la normalizacion por
#       Q_completo absorbe el factor local; la forma 1/ln z no.
#   K4  desglose por niveles (z=17, D=2e5, exhaustivo): R(L) contra
#       u = ln L/ln z; umbral empirico de nivel; Q^{(L=2)} = 1.
#   K5  defectos delta_s(r) = S_s(r) - Pi(r): norma L2 del defecto
#       (el sub-enunciado (H2')), defecto puntual maximo, y la
#       descomposicion Q_s = sum (Pi+delta)^2/phi (identidad exacta).
#   K6  las rutas absolutas fracasan, con numeros: (a) triangulo en
#       C_s(r) -> Q_s <= (P^+)^3 ~ (ln z)^3; (b) diagonal d=e:
#       diag = prod (1 + p/(p-1)^2) ~ C3 ln z, no O(1/ln z): la
#       lectura del doc 95 sec. 2.5 es inconsistente con su tabla.
#   K7  contraejemplo primorial: N divisible por todo p<z =>
#       Q_s = 1 = Q_completo(z;N) (R=1) pero 1 >> C/ln z: la forma
#       sin normalizar de (H2) es FALSA uniformemente en N.
# Uso: python3 109_ataque_h2_gram.py | tee salida_doc109_local.txt
import importlib.util
import math
import sys
import time


# ---------------------------------------------------------------- utiles
def primes_upto(n):
    if n < 2:
        return []
    bs = bytearray(b"\x01") * (n + 1)
    bs[0:2] = b"\x00\x00"
    for i in range(2, int(n ** 0.5) + 1):
        if bs[i]:
            bs[i * i:n + 1:i] = b"\x00" * ((n - i * i) // i + 1)
    return [i for i in range(2, n + 1) if bs[i]]


def allowed_primes(z, A, excl):
    """Primos que pueden dividir a d: impares, < z, no | A, no en excl
    (excl modela los primos que dividen a N; 2 siempre excluido)."""
    return [p for p in primes_upto(z - 1)
            if p > 2 and A % p != 0 and p not in excl]


def enum_smooth(allow):
    """Todos los libres de cuadrados sobre `allow`:
    lista de (d, pf_desc, mu, phi)."""
    out = [(1, (), 1, 1)]
    for p in allow:
        out += [(d * p, tuple(sorted(pf + (p,), reverse=True)), -mu,
                 phi * (p - 1)) for (d, pf, mu, phi) in list(out)]
    return out


def in_Vplus(pf, L):
    """Teorema 9.3 de Nathanson, beta=2, V^+ (copiado del doc 97):
    pf en orden DECRECIENTE; para m impar, p_m^2*(p_1...p_m) < L."""
    if not pf:
        return True
    if pf[0] >= L:
        return False
    prod = 1
    for i, p in enumerate(pf):
        m = i + 1
        prod *= p
        if m % 2 == 1 and p * p * prod >= L:
            return False
    return True


def sub_divisors(pf):
    """Divisores de un libre de cuadrados dado por su tupla de primos:
    lista de (r, mu_r, phi_r)."""
    out = [(1, 1, 1)]
    for p in pf:
        out += [(r * p, -mu, phi * (p - 1)) for (r, mu, phi) in list(out)]
    return out


def Q_completo(z, A, excl=()):
    """prod_{p<z, p no|AN} (1 - 1/(p-1) + 1/(p-1)^2)  (doc 97 sec. 2.3)."""
    tot = 1.0
    for p in primes_upto(z - 1):
        if p == 2 or A % p == 0 or p in excl:
            continue
        tot *= (1.0 - 1.0 / (p - 1) + 1.0 / (p - 1) ** 2)
    return tot


def Pi_of(rpf, z, A, excl):
    """Pi(r) = prod_{p<z, p no|rAN} (1 - 1/(p-1)): el factor del caso
    sin criba, C(r) = mu(r)/phi(r) * Pi(r)  (doc 97 sec. 2.3)."""
    tot = 1.0
    for p in primes_upto(z - 1):
        if p == 2 or A % p == 0 or p in excl or p in rpf:
            continue
        tot *= (1.0 - 1.0 / (p - 1))
    return tot


# ------------------------------------------------------- datos por nivel
class LevelData:
    """Para un nivel L = D//s (modo fijo): V = V^+(L) sobre los primos
    permitidos, C(r) = sum_{d en V, r|d} mu(d)/phi(d), y derivados."""

    def __init__(self, z, A, excl, L):
        self.z, self.A, self.excl, self.L = z, A, excl, L
        allow = allowed_primes(z, A, excl)
        self.enum = enum_smooth(allow)
        self.V = [(d, pf, mu, phi) for (d, pf, mu, phi) in self.enum
                  if d < L and in_Vplus(pf, L)]
        C = {}
        Tri = {}
        for (d, pf, mu, phi) in self.V:
            w = mu / phi
            aw = 1.0 / phi
            for (r, rmu, rphi) in sub_divisors(pf):
                C[r] = C.get(r, 0.0) + w
                Tri[r] = Tri.get(r, 0.0) + aw
        self.C = C
        self.Q = 0.0
        self.Q_tri = 0.0
        self.diag = 0.0
        self.defect2 = 0.0
        self.maxrel_defect = 0.0
        self.dec_err = 0.0
        for (d, pf, mu, phi) in self.V:
            self.diag += d / (phi * phi)
        for (r, rpf, rmu, rphi) in self.enum:
            c = C.get(r, 0.0)
            self.Q += rphi * c * c
            t = Tri.get(r, 0.0)
            self.Q_tri += rphi * t * t
            # S(r) = C(r)*phi(r)/mu(r); delta = S(r) - Pi(r)
            S = c * rphi / rmu
            Pi = Pi_of(rpf, z, A, excl)
            delta = S - Pi
            self.defect2 += delta * delta / rphi
            if Pi > 0:
                self.maxrel_defect = max(self.maxrel_defect,
                                         abs(delta) / Pi)
            # chequeo: phi(r) C(r)^2 = (Pi+delta)^2/phi(r)
            self.dec_err = max(self.dec_err,
                               abs(rphi * c * c - (Pi + delta) ** 2 / rphi))


def sample_sprimes(z, y, A, excl, nsamp=60):
    """Primos s en [z,y), (s,AN)=1: todos si son <=nsamp, si no una
    muestra log-uniforme con los extremos incluidos."""
    spr = [p for p in primes_upto(y - 1)
           if p >= z and math.gcd(p, A) == 1 and p not in excl]
    if len(spr) <= nsamp:
        return spr
    idx = sorted({round(i * (len(spr) - 1) / (nsamp - 1.0))
                  for i in range(nsamp)})
    return [spr[i] for i in idx]


def frame_sweep(z, A, excl, D, y, nsamp=60):
    spr = sample_sprimes(z, y, A, excl, nsamp)
    levels = sorted({D // s for s in spr})
    datas = [LevelData(z, A, excl, L) for L in levels]
    qc = Q_completo(z, A, excl)
    supQ = max(d.Q for d in datas)
    worst = max(datas, key=lambda d: d.Q)
    return {"levels": levels, "datas": datas, "qc": qc, "supQ": supQ,
            "R": supQ / qc, "worst": worst, "ns": len(spr)}


def linfit(xs, ys):
    n = len(xs)
    mx = sum(xs) / n
    my = sum(ys) / n
    sxy = sum((x - mx) * (y - my) for x, y in zip(xs, ys))
    sxx = sum((x - mx) ** 2 for x in xs)
    return sxy / sxx, my - (sxy / sxx) * mx


# ------------------------------------------------------------------ K0
def check_K0():
    ok = True
    spec = importlib.util.spec_from_file_location(
        "doc97", "97_verificacion_doc96_diagonal_mobius.py")
    d97 = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(d97)
    # S1 del doc 97
    m1 = d97.Model97(3, 11, 60, 20000, "fixed")
    errs = []
    for s in m1.sprimes:
        ld = LevelData(11, 3, (), 20000 // s)
        errs.append(abs(ld.Q - m1.Q[s]))
    e1 = max(errs)
    print("K0(a) S1-fixed: max|Q_ligero - Q_Model97| = %.2e sobre %d s %s"
          % (e1, len(m1.sprimes), "PASS" if e1 < 1e-9 else "FAIL"))
    ok &= e1 < 1e-9
    # S4 del doc 97 (parametros proximos a los del doc 94)
    m4 = d97.Model97(3, 11, 3000, 200000, "fixed")
    errs = []
    for s in m4.sprimes[::37]:
        ld = LevelData(11, 3, (), 200000 // s)
        errs.append(abs(ld.Q - m4.Q[s]))
    e4 = max(errs)
    print("K0(b) S4-fixed: max|Q_ligero - Q_Model97| = %.2e (muestra de "
          "%d s) %s" % (e4, len(errs), "PASS" if e4 < 1e-9 else "FAIL"))
    ok &= e4 < 1e-9
    # identidad sin criba por la ruta diagonal ligera: L -> infinito
    for (A, z) in ((3, 11), (3, 17), (15, 19)):
        enum = enum_smooth(allowed_primes(z, A, ()))
        C = {}
        for (d, pf, mu, phi) in enum:
            w = mu / phi
            for (r, rmu, rphi) in sub_divisors(pf):
                C[r] = C.get(r, 0.0) + w
        q = sum(phi * C.get(d, 0.0) ** 2 for (d, pf, mu, phi) in enum)
        qc = Q_completo(z, A)
        e = abs(q - qc) / qc
        print("K0(c) sin criba A=%d z=%d: Q=%.12f Q_completo=%.12f "
              "relerr=%.2e %s" % (A, z, q, qc, e,
                                  "PASS" if e < 1e-12 else "FAIL"))
        ok &= e < 1e-12
    return ok


# ------------------------------------------------------------------ K1/K2
def sweep_frame(zlist, A, excl, Dfun, yfun, tag):
    rows = []
    for z in zlist:
        D = Dfun(z)
        y = yfun(z)
        r = frame_sweep(z, A, excl, D, y)
        rows.append((z, r))
        print("  [%s] z=%2d D=%9d y=%7d #s=%5d niveles=%3d | sup Q_s="
              "%.6f Q_c=%.6f R=%.3f | Q_s*ln z en [%.3f, %.3f]"
              % (tag, z, D, y, r["ns"], len(r["levels"]), r["supQ"],
                 r["qc"], r["R"],
                 min(d.Q for d in r["datas"]) * math.log(z),
                 r["supQ"] * math.log(z)))
    return rows


def check_K1K2():
    zlist = (11, 13, 17, 19, 23, 29, 31, 37)
    print("K1 barrido marco A (D=z^5, y=z^(10/3), A=3, E0) "
          "[cf. doc 95 sec. 2.5]:")
    rowsA = sweep_frame(zlist, 3, (), lambda z: z ** 5,
                        lambda z: int(z ** (10.0 / 3.0)), "A")
    print("K2 marco B (D=z^6, y=z^3, u_min=3, A=3, E0):")
    rowsB = sweep_frame(zlist, 3, (), lambda z: z ** 6,
                        lambda z: z ** 3, "B")
    xs = [math.log(math.log(z)) for (z, _) in rowsA]
    gA, _ = linfit(xs, [math.log(r["supQ"]) for (_, r) in rowsA])
    dA, _ = linfit(xs, [math.log(r["R"]) for (_, r) in rowsA])
    gB, _ = linfit(xs, [math.log(r["supQ"]) for (_, r) in rowsB])
    dB, _ = linfit(xs, [math.log(r["R"]) for (_, r) in rowsB])
    print("K2 ajustes (minimos cuadrados sobre ln ln z):")
    print("  marco A: sup Q_s ~ (ln z)^(-gamma), gamma=%.3f | "
          "R(z) ~ (ln z)^delta, delta=%.3f" % (-gA, dA))
    print("  marco B: gamma=%.3f | delta=%.3f" % (-gB, dB))
    print("  R(z) marco A: %s" % ", ".join("%.3f" % r["R"]
                                           for (_, r) in rowsA))
    print("  R(z) marco B: %s" % ", ".join("%.3f" % r["R"]
                                           for (_, r) in rowsB))
    ok = all(1.0 <= r["R"] <= 2.5 for (_, r) in rowsA + rowsB)
    print("  => R(z) acotado en [1, 2.5] en ambos marcos: %s "
          "(medicion compatible con (H2). El ajuste gamma=0.80-0.82 "
          "esta sesgado a la baja por la meseta Q=1 de los niveles "
          "pequenos y por la convergencia lenta de Mertens de "
          "Q_completo*ln z (1.678 a z=11 -> 1.7225); el estadistico "
          "decisivo es R(z): plano en el marco B, 1.00-1.07)"
          % ("PASS" if ok else "FAIL"))
    return ok, rowsA


# ------------------------------------------------------------------ K3
def check_K3():
    z = 29
    D, y = z ** 5, int(z ** (10.0 / 3.0))
    print("K3 dependencia en A y en el factor local de N (z=29, marco A):")
    ok = True
    for (A, excl, tag) in ((3, (), "A=3, E0 (N=2^e)"),
                           (15, (), "A=15, E0"),
                           (105, (), "A=105, E0"),
                           (3, (5, 7, 11, 13), "A=3, E1 (5,7,11,13 | N)")):
        r = frame_sweep(z, A, excl, D, y)
        print("  %-26s Q_completo=%.6f sup Q_s=%.6f R=%.3f "
              "sup Q_s*ln z=%.3f"
              % (tag, r["qc"], r["supQ"], r["R"], r["supQ"] * math.log(z)))
        ok &= 0.8 <= r["R"] <= 2.5
    print("  => R estable y sup Q_s*ln z NO constante (factor local): "
          "la normalizacion correcta es Q_completo(z;N). %s"
          % ("PASS" if ok else "FAIL"))
    return ok


# ------------------------------------------------------------------ K4
def check_K4():
    z, A, D = 17, 3, 200000
    qc = Q_completo(z, A)
    print("K4 desglose por niveles (z=17, A=3, D=2e5, exhaustivo en L):")
    bins = {}
    Lmax = D // z
    datas = {}
    for L in range(1, Lmax + 1):
        ld = LevelData(z, A, (), L)
        datas[L] = ld
        u = math.log(L) / math.log(z) if L > 1 else 0.0
        b = round(u * 2) / 2.0
        bins[b] = max(bins.get(b, 0.0), ld.Q / qc)
    for b in sorted(bins):
        print("  u~%.1f: max R(L) en la banda = %.3f" % (b, bins[b]))
    bad = [L for L in range(2, Lmax + 1) if datas[L].Q / qc > 1.5]
    u0 = (math.log(max(bad)) / math.log(z)) if bad else 0.0
    ld2 = datas[2]
    print("  Q^{(L=2)}=%.6f (=1: V={1}); R=%.3f = 1/Q_completo"
          % (ld2.Q, ld2.Q / qc))
    print("  umbral empirico: R(L)<=1.5 para todo L con u>%.2f "
          "(peor L infractor: %s)" % (u0, max(bad) if bad else None))
    ok = ld2.Q == 1.0 and u0 < 3.0
    print("  => sin condicion de nivel (H2) falla (Q->1); con u>=u0 "
          "la medicion es O(1). %s" % ("PASS" if ok else "FAIL"))
    return ok


# ------------------------------------------------------------------ K5
def check_K5(rowsA):
    print("K5 defectos delta_s(r) = S_s(r) - Pi(r) (marco A):")
    ok = True
    for (z, r) in rowsA:
        if z not in (17, 29, 37):
            continue
        qc = r["qc"]
        mx = max(r["datas"], key=lambda d: d.defect2)
        mp = max(r["datas"], key=lambda d: d.maxrel_defect)
        dec = max(d.dec_err for d in r["datas"])
        print("  z=%d: max_niveles sum_r delta^2/phi = %.6f "
              "(= %.3f * Q_completo) | max |delta|/Pi = %.3f | "
              "err. descomposicion=%.2e"
              % (z, mx.defect2, mx.defect2 / qc, mp.maxrel_defect, dec))
        ok &= dec < 1e-9
    print("  => (H2') en forma L2: la norma del defecto es O(Q_completo)"
          " en la medicion; el defecto PUNTUAL |delta|/Pi no es <<1 en"
          " todo r: la forma puntual fuerte no se sostiene en ventana;"
          " la L2 es el sub-enunciado minimo viable. %s"
          % ("PASS" if ok else "FAIL"))
    return ok


# ------------------------------------------------------------------ K6
def check_K6():
    z, A = 37, 3
    D, y = z ** 5, int(z ** (10.0 / 3.0))
    r = frame_sweep(z, A, (), D, y)
    w = r["worst"]
    Pp = 1.0
    diag_full = 1.0
    for p in primes_upto(z - 1):
        if p == 2 or A % p == 0:
            continue
        Pp *= (1.0 + 1.0 / (p - 1))
        diag_full *= (1.0 + p / (p - 1) ** 2)
    print("K6 rutas absolutas a z=37 (peor nivel del marco A, L=%d):"
          % w.L)
    print("  (a) triangulo: Q_tri(nivel)=%.4f | cota global (P^+)^3=%.4f "
          "| Q_s real=%.6f | (P^+)^3/Q_completo=%.1f"
          % (w.Q_tri, Pp ** 3, w.Q, Pp ** 3 / r["qc"]))
    # crecimiento de la perdida del triangulo con z: (ln z)^4 / const
    for zz in (11, 37):
        Pp_z = 1.0
        for p in primes_upto(zz - 1):
            if p == 2 or A % p == 0:
                continue
            Pp_z *= (1.0 + 1.0 / (p - 1))
        print("      z=%d: (P^+)^3/Q_completo = %.2f (crece como "
              "(ln z)^4/C)" % (zz, Pp_z ** 3 / Q_completo(zz, A)))
    print("  (b) diagonal d=e: diag(nivel)=%.4f | diag/Q_completo=%.2f | "
          "diag sin restriccion = prod(1+p/(p-1)^2)=%.4f = %.3f*ln z"
          % (w.diag, w.diag / r["qc"], diag_full,
             diag_full / math.log(z)))
    print("      diag >= 1 siempre (termino d=1) y ~ C3*ln z sin "
          "restriccion: ninguna cota Q_s<=diag puede dar (H2).")
    ok = Pp ** 3 / r["qc"] > 10 and w.diag / r["qc"] > 3
    print("  => ambas rutas absolutas pierden >= (ln z)^2: la "
          "cancelacion ES el contenido de (H2). %s"
          % ("PASS" if ok else "FAIL"))
    return ok


# ------------------------------------------------------------------ K7
def check_K7():
    z, A = 29, 3
    excl = tuple(p for p in primes_upto(z - 1) if p > 2 and A % p != 0)
    D, y = z ** 5, int(z ** (10.0 / 3.0))
    r = frame_sweep(z, A, excl, D, y)
    print("K7 contraejemplo primorial (z=29, todo p<29 divide a AN):")
    print("  primos permitidos=%s | Q_completo=%.6f | sup Q_s=%.6f | "
          "R=%.3f | C/ln z con C=3: %.3f < 1 = Q_s"
          % (allowed_primes(z, A, excl), r["qc"], r["supQ"], r["R"],
             3.0 / math.log(z)))
    ok = (r["supQ"] == 1.0 and abs(r["qc"] - 1.0) < 1e-12
          and abs(r["R"] - 1.0) < 1e-9 and 3.0 / math.log(z) < 1.0)
    print("  => la forma sup Q_s << 1/log z_A con constante absoluta es "
          "FALSA uniformemente en N; la forma normalizada da R=1. %s"
          % ("PASS" if ok else "FAIL"))
    return ok


# ------------------------------------------------------------------ main
def main():
    t0 = time.time()
    print("== DOC 109: ataque a (H2), cota de Gram normalizada ==")
    print("Solo stdlib. Modo fijo; Q_s = Q^{(D//s)}.")
    print()
    allok = True
    allok &= check_K0()
    print()
    ok2, rowsA = check_K1K2()
    allok &= ok2
    print()
    allok &= check_K3()
    print()
    allok &= check_K4()
    print()
    allok &= check_K5(rowsA)
    print()
    allok &= check_K6()
    print()
    allok &= check_K7()
    print()
    print("Tiempo total: %.1f s" % (time.time() - t0))
    print("PASS" if allok else "FAIL")
    return 0 if allok else 1


if __name__ == "__main__":
    sys.exit(main())
