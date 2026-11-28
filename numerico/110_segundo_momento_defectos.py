#!/usr/bin/env python3
# 110_segundo_momento_defectos.py  (Doc 110)
# Expansion exacta del segundo momento de los defectos delta_s(r)
# (sub-enunciado (H2') del doc 109) y su reindexacion cerrada.
#
# Objetos (docs 96/97/109, modo de LIMITE FIJO, nivel L = D//s):
#   C(r) = sum_{d en V^+(L), r|d} mu(d)/phi(d),  Q = sum_r phi(r) C(r)^2
#   S(r) = C(r)*phi(r)/mu(r) = sum_{d' en V^{(r)}} mu(d')/phi(d'),
#     V^{(r)} = {d' : (d',r)=1, r*d' en V^+(L)}  (soporte desplazado)
#   Pi(r) = prod_{p en P, p no|r} (1 - 1/(p-1)) = sum_{(d',r)=1} mu(d')/phi(d')
#   delta(r) = S(r) - Pi(r) = -T(r),
#     T(r) = sum_{d' en T^{(r)}} mu(d')/phi(d'),
#     T^{(r)} = {d' : (d',r)=1, r*d' no en V^+(L)}  (cola del desplazado)
#   Def2 = sum_r delta(r)^2/phi(r)   <- el momento de (H2')
#
# Resultados que el script verifica a precision de maquina:
#   (110.1) Def2 = Q + Q_completo - 2*Pi_full*M_V,
#     Pi_full = prod_{p en P}(1 - 1/(p-1)),
#     M_V = sum_{d en V} g(d)/phi(d),  g(d) = prod_{p|d} 1/(p-2),
#     via el intercambio sum_r mu(r)Pi(r)C(r) = Pi_full * M_V.
#   (110.2) Q_completo = Pi_full * M_full  (Euler, exacto), luego
#     Def2 = (Q - Q_completo) + 2*Pi_full*M_cola,
#     M_cola = sum_{d no en V} g(d)/phi(d) = M_full - M_V >= 0,
#     y 0 <= 2*Pi_full*M_cola <= 2*Q_completo  [conteo puro].
#   Descomposicion de la cola: Def2 = D_cola + X_cola con
#     D_cola = sum_r (1/phi(r)) sum_{d' en T^{(r)}} 1/phi(d')^2 (diagonal),
#     X_cola = el resto (cruzada d' != e').
#
# Chequeos:
#   K0  validacion cruzada contra LevelData del doc 109 (Q y Def2).
#   K1  particion S(r) + T(r) = Pi(r) para todo r (soportes explicitos
#       V^{(r)} y T^{(r)} construidos desde cero); delta = -T.
#   K2  identidad (110.1): intercambio del termino cruzado y Def2.
#       K2b: sum_r Pi(r)^2/phi(r) = Q_completo.
#   K3  identidad (110.2): Q_completo = Pi_full*M_full (Euler y suma),
#       M_cola directo = M_full - M_V, y Def2 por la forma cerrada.
#   K4  descomposicion diagonal/cruzada de la cola en niveles elegidos:
#       D_cola + X_cola = Def2, con tamanos contra Q_completo.
#   K5  conteo: 0 <= 2*Pi_full*M_cola <= 2*Q_completo en TODOS los
#       niveles muestreados; cota absoluta de M_full; cota de conteo
#       para D_cola (pierde ~(ln z)^2: no basta).
#   K6  equivalencia cuantitativa (H2)<->(H2'): en todos los niveles,
#       Def2/Q_c - (R-1) = 2*M_cola/M_full en [0,2], R = Q/Q_c.
#   K7  independencia de datos (para el sec. 4 del doc): delta^{(L)}(r0)
#       varia con el nivel L a Pi(r0) fijo, y delta = -Pi(r0) exacto
#       cuando r0 no divide a ningun d en V; Delta_r (doc 97) no
#       depende de L en absoluto y crece con log N (aritmetica).
# Uso: python3 110_segundo_momento_defectos.py | tee salida_doc110_local.txt
import importlib.util
import math
import sys
import time

spec = importlib.util.spec_from_file_location(
    "doc109", "109_ataque_h2_gram.py")
d109 = importlib.util.module_from_spec(spec)
spec.loader.exec_module(d109)

primes_upto = d109.primes_upto
allowed_primes = d109.allowed_primes
enum_smooth = d109.enum_smooth
in_Vplus = d109.in_Vplus
sub_divisors = d109.sub_divisors
Q_completo = d109.Q_completo
Pi_of = d109.Pi_of
LevelData = d109.LevelData
sample_sprimes = d109.sample_sprimes


# ------------------------------------------------------------ nivel 110
def g_of_pf(pf):
    """g(d) = prod_{p|d} 1/(p-2)  (g(1) = 1)."""
    g = 1.0
    for p in pf:
        g /= (p - 2)
    return g


class Level110:
    """Nivel L: V, C, Q, defectos por la ruta de C (como doc 109) y
    todas las piezas de la forma cerrada (110.1)/(110.2)."""

    def __init__(self, z, A, excl, L):
        self.z, self.A, self.excl, self.L = z, A, excl, L
        allow = allowed_primes(z, A, excl)
        self.allow = allow
        self.enum = enum_smooth(allow)
        self.V = [(d, pf, mu, phi) for (d, pf, mu, phi) in self.enum
                  if d < L and in_Vplus(pf, L)]
        self.Vset = {d for (d, pf, mu, phi) in self.V}
        C = {}
        for (d, pf, mu, phi) in self.V:
            w = mu / phi
            for (r, rmu, rphi) in sub_divisors(pf):
                C[r] = C.get(r, 0.0) + w
        self.C = C
        self.qc = Q_completo(z, A, excl)
        # Pi_full y productos de Euler
        self.Pi_full = 1.0
        self.M_full_prod = 1.0
        for p in allow:
            self.Pi_full *= (1.0 - 1.0 / (p - 1))
            self.M_full_prod *= (1.0 + 1.0 / ((p - 1) * (p - 2)))
        # M_V, M_full por suma, M_cola directo
        self.M_V = 0.0
        self.M_full_sum = 0.0
        self.M_cola_dir = 0.0
        for (d, pf, mu, phi) in self.enum:
            v = g_of_pf(pf) / phi
            self.M_full_sum += v
            if d in self.Vset:
                self.M_V += v
            else:
                self.M_cola_dir += v
        self.M_cola = self.M_full_sum - self.M_V
        # Q, Def2 por la ruta de C; termino cruzado; sum Pi^2/phi
        self.Q = 0.0
        self.Def2 = 0.0
        self.cross = 0.0
        self.Pi2 = 0.0
        self.S = {}
        self.Pi = {}
        self.delta = {}
        for (r, rpf, rmu, rphi) in self.enum:
            c = C.get(r, 0.0)
            self.Q += rphi * c * c
            S = c * rphi / rmu
            Pi = Pi_of(rpf, z, A, excl)
            self.S[r] = S
            self.Pi[r] = Pi
            self.delta[r] = S - Pi
            self.Def2 += (S - Pi) ** 2 / rphi
            self.cross += rmu * Pi * c
            self.Pi2 += Pi * Pi / rphi
        # formas cerradas
        self.forma1 = self.Q + self.qc - 2.0 * self.Pi_full * self.M_V
        self.forma2 = (self.Q - self.qc
                       + 2.0 * self.Pi_full * self.M_cola)

    def tail_decomp(self):
        """Descomposicion de la cola con soportes desplazados
        EXPLICITOS: para cada r, V^{(r)} y T^{(r)} desde cero.
        Devuelve (Def2_cola, D_cola, X_cola, err_particion, err_S)."""
        Def2t = 0.0
        D = 0.0
        X = 0.0
        err_part = 0.0
        err_S = 0.0
        enum = self.enum
        Vset = self.Vset
        for (r, rpf, rmu, rphi) in enum:
            Sdir = 0.0
            Tdir = 0.0
            diagT = 0.0
            for (dp, dpf, dmu, dphi) in enum:
                if math.gcd(dp, r) != 1:
                    continue
                a = dmu / dphi
                if r * dp in Vset:
                    Sdir += a
                else:
                    Tdir += a
                    diagT += 1.0 / (dphi * dphi)
            err_part = max(err_part,
                           abs(Sdir + Tdir - self.Pi[r]))
            err_S = max(err_S, abs(Sdir - self.S[r]))
            Def2t += Tdir * Tdir / rphi
            D += diagT / rphi
            X += (Tdir * Tdir - diagT) / rphi
        return Def2t, D, X, err_part, err_S


def frame_levels(z, A, excl, D, y, nsamp=60):
    spr = sample_sprimes(z, y, A, excl, nsamp)
    levels = sorted({D // s for s in spr})
    return [Level110(z, A, excl, L) for L in levels]


def relerr(a, b):
    return abs(a - b) / max(abs(b), 1e-300)


# ------------------------------------------------------------------ K0
def check_K0():
    print("K0 validacion cruzada contra LevelData del doc 109:")
    ok = True
    for (z, D, y, Lfrac) in ((17, 17 ** 5, int(17 ** (10.0 / 3.0)), 0.5),
                             (29, 29 ** 5, int(29 ** (10.0 / 3.0)), 0.3),
                             (37, 37 ** 5, int(37 ** (10.0 / 3.0)), 0.7)):
        spr = sample_sprimes(z, y, 3, ())
        L = (D // spr[int(Lfrac * (len(spr) - 1))])
        ld = LevelData(z, 3, (), L)
        lv = Level110(z, 3, (), L)
        eQ = abs(lv.Q - ld.Q)
        eD = abs(lv.Def2 - ld.defect2)
        print("  z=%d L=%d: |Q110-Q109|=%.2e |Def2_110-Def2_109|=%.2e %s"
              % (z, L, eQ, eD, "PASS" if eQ < 1e-9 and eD < 1e-9
                 else "FAIL"))
        ok &= eQ < 1e-9 and eD < 1e-9
    return ok


# ------------------------------------------------------- K1/K2/K3/K5/K6
def check_identities():
    print("K1/K2/K3/K5/K6 identidades y conteo, marcos A (z=17,29,37) "
          "y B (z=29):")
    ok = True
    frames = []
    for z in (17, 29, 37):
        frames.append(("A", z, z ** 5, int(z ** (10.0 / 3.0))))
    frames.append(("B", 29, 29 ** 6, 29 ** 3))
    for (tag, z, D, y) in frames:
        levels = frame_levels(z, 3, (), D, y)
        e_part_max = 0.0
        worst_tau = 0.0
        e_forma = 0.0
        for lv in levels:
            # K2b y K3 (por nivel basta z: Pi2 y productos no dependen
            # de L; se verifican en cada nivel igualmente)
            ok &= relerr(lv.Pi2, lv.qc) < 1e-12
            ok &= relerr(lv.Pi_full * lv.M_full_prod, lv.qc) < 1e-12
            ok &= relerr(lv.M_full_sum, lv.M_full_prod) < 1e-12
            # M_cola por diferencia sufre cancelacion catastrofica
            # cuando la cola es minuscula: el chequeo es ABSOLUTO
            # (las sumas son O(1); el acuerdo real es <= 4e-17).
            ok &= abs(lv.M_cola_dir - lv.M_cola) < 1e-11
            # K2: intercambio del cruzado y forma (110.1). Las formas
            # cerradas restan cantidades O(1) (y Def2 puede ser 0
            # exacto si V = enum): tolerancia relativa a la escala
            # Q + Q_c, no a Def2.
            ok &= relerr(lv.cross, lv.Pi_full * lv.M_V) < 1e-9
            esc = lv.Q + lv.qc
            e_forma = max(e_forma, abs(lv.Def2 - lv.forma1) / esc,
                          abs(lv.Def2 - lv.forma2) / esc)
            ok &= abs(lv.Def2 - lv.forma1) < 1e-9 * esc
            # K3: forma (110.2)
            ok &= abs(lv.Def2 - lv.forma2) < 1e-9 * esc
            # K5: 0 <= 2 Pi_full M_cola <= 2 Q_c
            tau = 2.0 * lv.Pi_full * lv.M_cola
            ok &= -1e-15 <= tau <= 2.0 * lv.qc + 1e-15
            worst_tau = max(worst_tau, tau / lv.qc)
            # K6: Def2/Q_c - (R-1) = 2 M_cola/M_full en [0,2]
            equiv = lv.Def2 / lv.qc - (lv.Q / lv.qc - 1.0)
            ok &= abs(equiv - 2.0 * lv.M_cola / lv.M_full_sum) < 1e-9
            ok &= -1e-12 <= equiv <= 2.0 + 1e-12
        # K1: particion con soportes explicitos en 2 niveles por marco
        for lv in (levels[len(levels) // 3], levels[-1]):
            _, _, _, e_part, e_S = lv.tail_decomp()
            e_part_max = max(e_part_max, e_part, e_S)
        ok &= e_part_max < 1e-9
        mx = max(levels, key=lambda t: t.Def2)
        print("  [%s] z=%d: %d niveles | err particion (K1, 2 niveles)"
              "=%.2e | err formas (110.1)/(110.2) rel. a Q+Q_c = %.1e | "
              "max Def2/Q_c=%.3f (L=%d) | max 2Pi_full*M_cola/Q_c"
              "=%.3f | M_full=%.4f Pi_full=%.4f Q_c=%.4f"
              % (tag, z, len(levels), e_part_max, e_forma,
                 mx.Def2 / mx.qc, mx.L,
                 worst_tau, levels[0].M_full_sum, levels[0].Pi_full,
                 levels[0].qc))
    print("  => (110.1) y (110.2) exactas en todos los niveles; el "
          "termino de conteo 2*Pi_full*M_cola esta en [0, 2*Q_c] "
          "siempre. %s" % ("PASS" if ok else "FAIL"))
    return ok


# ------------------------------------------------------------------ K4
def check_K4():
    print("K4 descomposicion diagonal/cruzada de la cola (soportes "
          "desplazados explicitos):")
    ok = True
    cases = []
    for z in (17, 29, 37):
        cases.append(("A", z, z ** 5, int(z ** (10.0 / 3.0))))
    cases.append(("B", 29, 29 ** 6, 29 ** 3))
    for (tag, z, D, y) in cases:
        levels = frame_levels(z, 3, (), D, y)
        picks = {}
        picks["maxDef2"] = max(levels, key=lambda t: t.Def2)
        picks["maxQ"] = max(levels, key=lambda t: t.Q)
        picks["minL"] = levels[0]
        picks["maxL"] = levels[-1]
        seen = set()
        for (name, lv) in picks.items():
            if lv.L in seen:
                continue
            seen.add(lv.L)
            Def2t, Dc, Xc, e_part, e_S = lv.tail_decomp()
            # error ABSOLUTO: en niveles con V = enum completo,
            # Def2 = 0 exacto por la ruta de la cola y ~1e-33 por la
            # ruta de C (redondeo); el error relativo no significa nada.
            e = abs(Def2t - lv.Def2)
            ok &= e < 1e-9 and abs(Dc + Xc - Def2t) < 1e-9 * max(1.0,
                                                                Def2t)
            print("  [%s] z=%d %-7s L=%9d: Def2=%.6f (err abs vs ruta "
                  "C %.1e) | D_cola=%.4f (%.2f Q_c) | X_cola=%+.4f "
                  "(%+.2f Q_c) | (Q-Q_c)/Q_c=%+.3f | 2Pi_full*M_cola/Q_c"
                  "=%.3f"
                  % (tag, z, name, lv.L, Def2t, e, Dc, Dc / lv.qc, Xc,
                     Xc / lv.qc, (lv.Q - lv.qc) / lv.qc,
                     2.0 * lv.Pi_full * lv.M_cola / lv.qc))
    print("  => la diagonal de la cola NO basta: D_cola >> Q_c y la "
          "cruzada cancela; tras reindexar, el contenido no-contable "
          "es (Q-Q_c). %s" % ("PASS" if ok else "FAIL"))
    return ok


# ------------------------------------------------------------------ K5b
def check_K5b():
    print("K5b cotas de conteo puras:")
    # cota absoluta de M_full (A=3: primos impares != 3)
    prod = 1.0
    for p in primes_upto(100000):
        if p <= 3:
            continue
        prod *= (1.0 + 1.0 / ((p - 1) * (p - 2)))
    print("  M_full(z->inf, A=3) <= prod_{3<p<=1e5}(1+1/((p-1)(p-2)))"
          " = %.4f (producto convergente, cota absoluta)" % prod)
    ok = prod < 1.25
    # cota de conteo para D_cola: D <= P^+ * prod(1 + 1/(p-1)^2)
    for z in (17, 29, 37):
        allow = allowed_primes(z, 3, ())
        Pp = 1.0
        sq = 1.0
        for p in allow:
            Pp *= (1.0 + 1.0 / (p - 1))
            sq *= (1.0 + 1.0 / (p - 1) ** 2)
        qc = Q_completo(z, 3)
        print("  z=%d: cota conteo D_cola <= P^+*prod(1+1/(p-1)^2) = "
              "%.3f = %.1f * Q_c  (pierde ~(ln z)^2 contra el objetivo "
              "O(Q_c))" % (z, Pp * sq, Pp * sq / qc))
    print("  => el conteo cierra M_cola por completo y NO cierra "
          "D_cola. %s" % ("PASS" if ok else "FAIL"))
    return ok


# ------------------------------------------------------------------ K7
def check_K7():
    print("K7 independencia de datos delta vs Delta_r (sec. 4):")
    z, A = 17, 3
    D, y = z ** 5, int(z ** (10.0 / 3.0))
    levels = frame_levels(z, A, (), D, y)
    r0 = 35  # 5*7, libre de cuadrados
    Pi_r0 = None
    rows = []
    for lv in levels[::max(1, len(levels) // 6)][:7]:
        for (r, rpf, rmu, rphi) in lv.enum:
            if r == r0:
                Pi_r0 = lv.Pi[r]
                # r0 divide a algun d en V?
                divides = any(d % r0 == 0 for (d, _, _, _) in lv.V)
                rows.append((lv.L, lv.delta[r], divides))
    print("  delta^{(L)}(35) a Pi(35)=%.6f FIJO (no depende de L):"
          % Pi_r0)
    for (L, dv, divides) in rows:
        print("    L=%8d: delta=%+.6f | 35 divide a algun d en V: %s%s"
              % (L, dv, divides,
                 " | delta = -Pi exacto" if (not divides
                     and abs(dv + Pi_r0) < 1e-12) else ""))
    vals = {round(dv, 9) for (_, dv, _) in rows}
    ok = len(vals) >= 2
    # aritmetica de Delta_r para r primo (doc 97, Teorema 97.1):
    # Delta_p = (log p) * floor(log_p N): depende de N, no de L.
    for N in (10 ** 6, 10 ** 12):
        d5 = math.log(5) * int(math.log(N) / math.log(5))
        print("    Delta_5(N=%d) = (log 5)*floor(log_5 N) = %.2f "
              "(crece con N; delta_5 en el marco es O(Pi) ~ %.3f, sin N)"
              % (N, d5, Pi_of((5,), z, A, ())))
    print("  => delta vive en los datos (z, A, excl, L); Delta_r vive "
          "en (r, N, caracteres). Datos disjuntos. %s"
          % ("PASS" if ok else "FAIL"))
    return ok


# ------------------------------------------------------------------ main
def main():
    t0 = time.time()
    print("== DOC 110: segundo momento de los defectos, expandido ==")
    print("Solo stdlib. Modo fijo; nivel L = D//s.")
    print()
    allok = True
    allok &= check_K0()
    print()
    allok &= check_identities()
    print()
    allok &= check_K4()
    print()
    allok &= check_K5b()
    print()
    allok &= check_K7()
    print()
    print("Tiempo total: %.1f s" % (time.time() - t0))
    print("PASS" if allok else "FAIL")
    return 0 if allok else 1


if __name__ == "__main__":
    sys.exit(main())
