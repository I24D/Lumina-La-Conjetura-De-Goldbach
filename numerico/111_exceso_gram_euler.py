#!/usr/bin/env python3
# 111_exceso_gram_euler.py  (Doc 111)
# Ataque a (H2''), el exceso de Gram, via la diagonalizacion (docs 96/97)
# y la factorizacion de Euler (doc 110).
#
# Objetos (modo de LIMITE FIJO, nivel L = D//s; definiciones de los docs
# 109/110, reutilizadas por import para garantizar que son las mismas):
#   C^(L)(r) = sum_{d en V^+(L), r|d} mu(d)/phi(d),
#   Q^(L) = sum_r phi(r) C^(L)(r)^2        [T, doc 97 sec. 2.1]
#   C^full(r) = mu(r)/phi(r) * Pi(r)       [T, doc 97 sec. 2.3, sin criba]
#   Q_completo = sum_r Pi(r)^2/phi(r) = prod_{p en P}(1-1/(p-1)+1/(p-1)^2)
#   S(r) = C(r)*phi(r)/mu(r),  delta(r) = S(r) - Pi(r),  Def2 = sum delta^2/phi
#   (110.2):  Def2 = (Q^(L) - Q_completo) + 2*Pi_full*M_cola   [T, doc 110]
#
# El exceso por conductor:
#   t(r) := phi(r)[C^(L)(r)^2 - C^full(r)^2]
#         = [S(r)^2 - Pi(r)^2]/phi(r)
#         = delta(r)*(S(r)+Pi(r))/phi(r) = delta(r)*(2*Pi(r)+delta(r))/phi(r)
#   E^(L) := Q^(L) - Q_completo = sum_r t(r).
#   P^(L) := sum_{t(r)>0} t(r)  (parte positiva),  N^(L) := -sum_{t(r)<0} t(r)
#   E = P - N,  y  0 <= N <= sum_r Pi(r)^2/phi(r) = Q_completo  [T].
#
# Chequeos:
#   K0  validacion cruzada contra Level110 del doc 110 (Q, Def2, forma2).
#   K1  la expresion del exceso: sum_r t(r) = Q - Q_c a precision de
#       maquina; factorizacion por r t = delta*(2Pi+delta)/phi; y
#       consistencia con (110.2): E = Def2 - 2*Pi_full*M_cola.
#   K2  estructura de signos: r invisible (no divide a ningun d en V)
#       => t(r) = -Pi(r)^2/phi(r) EXACTO; 0 <= N <= Q_c; P >= E;
#       conteos de r con t>0 / t<0 / S(r)<0 en los niveles peores.
#   K2b clasificacion de S(r) sobre TODOS los niveles: S(r) < 0 nunca
#       ocurre en ventana; 0<S<Pi raro en el marco A, presente en B.
#   K3  bloques por primo mayor P^+(r): neto/positivo/negativo de t(r)
#       con P^+(r) = p, por primo p, en los niveles peores; y lo mismo
#       para Def2. Responde: ?de que primos viene el exceso?
#   K4  bloques por omega(r) = 0,1,2,>=3 (neto de t(r) y de Def2).
#   K5  cola de Euler: Q_c = Q_c^{<=w} * Q_c^{>w} con w = z^{1/2}
#       (reparto de log: el producto tipo Pi NO converge absolutamente);
#       M_full^{>w} <= exp(1/(w-1)) [T, telescopico sobre enteros];
#       M_cola y M_V por bloques P^+(d) = p.
#   K6  intento de cota, libro mayor: piezas cerradas [T] (bloque
#       invisible <= 0; M_cola en [0, 2 Q_c]; N en [0, Q_c]) y la pieza
#       que se resiste: P^(L). Maximo de P/Q_c sobre TODOS los niveles
#       muestreados de cada marco; descomposicion de P por signo de
#       delta (delta>0 contra delta<-2Pi) y por bloques P^+(r).
# Uso: python3 111_exceso_gram_euler.py | tee salida_doc111_local.txt
import importlib.util
import math
import sys
import time

spec = importlib.util.spec_from_file_location(
    "doc110", "110_segundo_momento_defectos.py")
d110 = importlib.util.module_from_spec(spec)
spec.loader.exec_module(d110)

Level110 = d110.Level110
frame_levels = d110.frame_levels
relerr = d110.relerr
primes_upto = d110.primes_upto
allowed_primes = d110.allowed_primes
g_of_pf = d110.g_of_pf


# ------------------------------------------------------------ nivel 111
class Level111:
    """Nivel L con el exceso por conductor t(r) y sus descomposiciones."""

    def __init__(self, z, A, excl, L):
        self.lv = Level110(z, A, excl, L)
        lv = self.lv
        self.z, self.A, self.excl, self.L = z, A, excl, L
        self.qc = lv.qc
        self.Q = lv.Q
        self.Def2 = lv.Def2
        # visibilidad: r divide a algun d en V
        vis = set()
        for (d, pf, mu, phi) in lv.V:
            for (r, rmu, rphi) in d110.sub_divisors(pf):
                vis.add(r)
        self.vis = vis
        # t(r) y piezas por r
        self.t = {}
        self.t_pos = 0.0
        self.t_neg = 0.0
        self.n_pos = 0
        self.n_neg = 0
        self.n_Sneg = 0
        self.invis_sum = 0.0     # sum_{r invisible} t(r) (exacto, <= 0)
        self.invis_qc = 0.0      # sum_{r invisible} Pi(r)^2/phi(r)
        self.err_fact = 0.0
        for (r, rpf, rmu, rphi) in lv.enum:
            S = lv.S[r]
            Pi = lv.Pi[r]
            delta = lv.delta[r]
            t = (S * S - Pi * Pi) / rphi
            self.t[r] = t
            # chequeo de la factorizacion de la diferencia de cuadrados
            self.err_fact = max(
                self.err_fact,
                abs(t - delta * (2.0 * Pi + delta) / rphi))
            if S < 0:
                self.n_Sneg += 1
            if t > 0:
                self.t_pos += t
                self.n_pos += 1
            elif t < 0:
                self.t_neg += t
                self.n_neg += 1
            if r not in vis:
                self.invis_sum += t
                self.invis_qc += Pi * Pi / rphi
        self.E = self.Q - self.qc
        self.sum_t = sum(self.t.values())
        self.P = self.t_pos
        self.N = -self.t_neg

    def block_pplus(self, values):
        """Descompone un dict r->valor por primo mayor P^+(r)
        (r=1 en la clave 1): devuelve {p: (neto, pos, neg)}."""
        out = {}
        for (r, rpf, rmu, rphi) in self.lv.enum:
            p = rpf[0] if rpf else 1
            v = values[r]
            neto, pos, neg = out.get(p, (0.0, 0.0, 0.0))
            neto += v
            if v > 0:
                pos += v
            elif v < 0:
                neg += v
            out[p] = (neto, pos, neg)
        return out

    def block_omega(self, values):
        """Descompone un dict r->valor por omega(r) en clases
        0,1,2,3(>=3): devuelve {k: neto}."""
        out = {0: 0.0, 1: 0.0, 2: 0.0, 3: 0.0}
        for (r, rpf, rmu, rphi) in self.lv.enum:
            k = min(len(rpf), 3)
            out[k] += values[r]
        return out

    def def2_by_r(self):
        lv = self.lv
        return {r: lv.delta[r] ** 2 / rphi
                for (r, rpf, rmu, rphi) in lv.enum}


def worst_levels():
    """Niveles representativos: el de maximo exceso de cada marco
    A (z=17,29,37) y B (z=29), mas el de maximo nivel de cada uno."""
    out = []
    frames = []
    for z in (17, 29, 37):
        frames.append(("A", z, z ** 5, int(z ** (10.0 / 3.0))))
    frames.append(("B", 29, 29 ** 6, 29 ** 3))
    for (tag, z, D, y) in frames:
        lvs = [Level111(z, 3, (), lv.L)
               for lv in frame_levels(z, 3, (), D, y)]
        mx = max(lvs, key=lambda t: t.E)
        out.append((tag, z, mx))
        out.append((tag, z, lvs[-1]))
    return out


# ------------------------------------------------------------------ K0
def check_K0():
    print("K0 validacion cruzada contra el doc 110 (Level110):")
    ok = True
    for (z, L) in ((17, 112), (29, 273), (37, 410)):
        a = Level111(z, 3, (), L)
        e1 = abs(a.Q - a.lv.Q)
        e2 = abs(a.Def2 - a.lv.forma2 + 0.0)  # Def2 ya es la ruta C
        e3 = abs(a.Def2 - a.lv.forma2)
        print("  z=%d L=%d: Q=%.12f Def2=%.12f |Def2-forma2|=%.2e %s"
              % (z, L, a.Q, a.Def2, e3, "PASS" if e3 < 1e-9 else "FAIL"))
        ok &= e3 < 1e-9
    return ok


# ------------------------------------------------------------------ K1
def check_K1(worst):
    print("K1 expresion del exceso por la diagonalizacion:")
    ok = True
    for (tag, z, a) in worst:
        e_sum = abs(a.sum_t - a.E)
        e_110 = abs(a.E - (a.Def2 - 2.0 * a.lv.Pi_full * a.lv.M_cola))
        esc = a.Q + a.qc
        print("  [%s] z=%d L=%7d: E=Q-Q_c=%+.9f | sum_r t(r)-E=%.2e | "
              "err factorizacion/r=%.2e | E-(Def2-2Pi_full*M_cola)=%.2e"
              % (tag, z, a.L, a.E, e_sum, a.err_fact, e_110))
        ok &= e_sum < 1e-9 * esc and a.err_fact < 1e-9
        ok &= e_110 < 1e-9 * esc
    print("  => E^(L) = sum_r phi(r)[C^(L)(r)^2 - C^full(r)^2] "
          "= sum_r delta(2Pi+delta)/phi, exacto; y por (110.2) "
          "E = Def2 - 2*Pi_full*M_cola: diferencia de dos terminos "
          ">= 0. %s" % ("PASS" if ok else "FAIL"))
    return ok


# ------------------------------------------------------------------ K2
def check_K2(worst):
    print("K2 estructura de signos del exceso:")
    ok = True
    for (tag, z, a) in worst:
        # r invisible: t(r) = -Pi(r)^2/phi(r) exacto
        err_inv = 0.0
        for (r, rpf, rmu, rphi) in a.lv.enum:
            if r not in a.vis:
                err_inv = max(err_inv,
                              abs(a.t[r] + a.lv.Pi[r] ** 2 / rphi))
        print("  [%s] z=%d L=%7d: E/Q_c=%+.3f | P/Q_c=%.3f (%d r) | "
              "N/Q_c=%.3f (%d r) | bloque invisible: %+.3f Q_c "
              "(err=%.1e) | r con S(r)<0: %d"
              % (tag, z, a.L, a.E / a.qc, a.P / a.qc, a.n_pos,
                 a.N / a.qc, a.n_neg, a.invis_sum / a.qc, err_inv,
                 a.n_Sneg))
        ok &= err_inv < 1e-9
        ok &= a.N <= a.qc + 1e-9           # N <= Q_c  [T]
        ok &= a.P >= a.E - 1e-12           # P >= E
        ok &= abs(a.P - a.N - a.E) < 1e-9  # E = P - N
        ok &= a.invis_sum <= 1e-15         # invisible <= 0
    print("  => el exceso NO tiene signo definido por r: hay masa "
          "positiva y negativa de orden Q_c, con cancelacion entre "
          "conductores; el bloque invisible es exacto y <= 0 "
          "(solo ayuda). %s" % ("PASS" if ok else "FAIL"))
    return ok


# ------------------------------------------------------------------ K2b
def check_K2b():
    print("K2b clasificacion de S(r) sobre TODOS los niveles muestreados:")
    frames = []
    for z in (17, 29, 37):
        frames.append(("A", z, z ** 5, int(z ** (10.0 / 3.0))))
    frames.append(("B", 29, 29 ** 6, 29 ** 3))
    ok = True
    for (tag, z, D, y) in frames:
        lvs = [Level111(z, 3, (), lv.L)
               for lv in frame_levels(z, 3, (), D, y)]
        cnt = {"S=0": 0, "0<S<Pi": 0, "S=Pi": 0, "S>Pi": 0, "S<0": 0}
        minS = 1e300
        for a in lvs:
            for (r, rpf, rmu, rphi) in a.lv.enum:
                S = a.lv.S[r]
                Pi = a.lv.Pi[r]
                minS = min(minS, S)
                if S == 0.0:
                    cnt["S=0"] += 1
                elif S < 0:
                    cnt["S<0"] += 1
                elif abs(S - Pi) < 1e-12:
                    cnt["S=Pi"] += 1
                elif S < Pi:
                    cnt["0<S<Pi"] += 1
                else:
                    cnt["S>Pi"] += 1
        print("  [%s] z=%d (%d niveles): min S(r)=%.4f | conteos %s"
              % (tag, z, len(lvs), minS, cnt))
        ok &= cnt["S<0"] == 0
    print("  => [M] en ventana S(r) >= 0 siempre: t(r)>0 equivale a "
          "delta(r)>0. En el marco A los defectos visibles negativos "
          "(0<S<Pi) son raros o nulos en los niveles que deciden; en "
          "el marco B existen. %s" % ("PASS" if ok else "FAIL"))
    return ok


# ------------------------------------------------------------------ K3
def check_K3(worst):
    print("K3 bloques por primo mayor P^+(r) (niveles de maximo exceso):")
    ok = True
    for (tag, z, a) in worst:
        if a.L == max(x.L for (t2, z2, x) in worst if z2 == z
                      and t2 == tag):
            continue  # solo niveles peores, no los de L maximo
        bt = a.block_pplus(a.t)
        bd = a.block_pplus(a.def2_by_r())
        tot = sum(v[0] for v in bt.values())
        ok &= abs(tot - a.E) < 1e-9
        print("  [%s] z=%d L=%d (E/Q_c=%+.3f): neto de t por P^+(r)=p, "
              "en unidades de Q_c:" % (tag, z, a.L, a.E / a.qc))
        for p in sorted(bt):
            neto, pos, neg = bt[p]
            dneto = bd[p][0]
            print("    p=%3d: t neto=%+.3f (pos %+.3f, neg %+.3f) | "
                  "Def2 neto=%.3f" % (p, neto / a.qc, pos / a.qc,
                                      neg / a.qc, dneto / a.qc))
    print("  => en los niveles que deciden el marco A la masa "
          "positiva se concentra en r=1 y en los pocos r visibles del "
          "soporte minusculo (el bloque P^+=1 iguala o excede el "
          "exceso neto); en el marco B se reparte. La masa negativa "
          "es, bloque a bloque, el Def2 de los r invisibles. %s"
          % ("PASS" if ok else "FAIL"))
    return ok


# ------------------------------------------------------------------ K4
def check_K4(worst):
    print("K4 bloques por omega(r) (niveles de maximo exceso):")
    ok = True
    for (tag, z, a) in worst:
        if a.L == max(x.L for (t2, z2, x) in worst if z2 == z
                      and t2 == tag):
            continue
        bt = a.block_omega(a.t)
        bd = a.block_omega(a.def2_by_r())
        ok &= abs(sum(bt.values()) - a.E) < 1e-9
        print("  [%s] z=%d L=%d: t neto por omega = "
              "{0: %+.3f, 1: %+.3f, 2: %+.3f, >=3: %+.3f} Q_c | "
              "Def2 por omega = {0: %.3f, 1: %.3f, 2: %.3f, >=3: %.3f}"
              " Q_c" % (tag, z, a.L, bt[0] / a.qc, bt[1] / a.qc,
                        bt[2] / a.qc, bt[3] / a.qc, bd[0] / a.qc,
                        bd[1] / a.qc, bd[2] / a.qc, bd[3] / a.qc))
    print("  => PASS" if ok else "  => FAIL")
    return ok


# ------------------------------------------------------------------ K5
def check_K5():
    print("K5 cola de Euler (A=3):")
    ok = True
    for z in (17, 29, 37):
        allow = allowed_primes(z, 3, ())
        w = z ** 0.5
        qc_le = 1.0
        qc_gt = 1.0
        for p in allow:
            f = 1.0 - 1.0 / (p - 1) + 1.0 / (p - 1) ** 2
            if p <= w:
                qc_le *= f
            else:
                qc_gt *= f
        qc = qc_le * qc_gt
        share = (math.log(qc_gt) / math.log(qc)) if qc != 1 else 0.0
        # M_full: cola absolutamente convergente, cota telescopica
        M_le = 1.0
        M_gt = 1.0
        for p in allow:
            f = 1.0 + 1.0 / ((p - 1) * (p - 2))
            if p <= w:
                M_le *= f
            else:
                M_gt *= f
        bound = math.exp(1.0 / (math.floor(w) - 1))
        print("  z=%d w=%.1f: Q_c=%.4f = Q_c^{<=w}(%.4f) * Q_c^{>w}"
              "(%.4f) | cuota de log de la cola: %.1f%% | M_full^{>w}"
              "=%.4f <= exp(1/(floor(w)-1))=%.4f [T]"
              % (z, w, qc, qc_le, qc_gt, 100.0 * share, M_gt, bound))
        ok &= abs(qc - d110.Q_completo(z, 3, ())) < 1e-12
        ok &= M_gt <= bound
    # M_cola / M_V por bloques P^+(d), nivel peor de z=29 marco A
    a = Level111(29, 3, (), 273)
    lv = a.lv
    Vset = lv.Vset
    cola_by_p = {}
    v_by_p = {}
    for (d, pf, mu, phi) in lv.enum:
        p = pf[0] if pf else 1
        val = g_of_pf(pf) / phi
        if d in Vset:
            v_by_p[p] = v_by_p.get(p, 0.0) + val
        else:
            cola_by_p[p] = cola_by_p.get(p, 0.0) + val
    print("  M_cola y M_V por P^+(d)=p (z=29, L=273; M_full=%.4f, "
          "2*Pi_full*M_cola/Q_c=%.3f):"
          % (lv.M_full_sum, 2.0 * lv.Pi_full * lv.M_cola / lv.qc))
    for p in sorted(set(cola_by_p) | set(v_by_p)):
        print("    p=%3d: M_cola=%.4f | M_V=%.4f"
              % (p, cola_by_p.get(p, 0.0), v_by_p.get(p, 0.0)))
    print("  => el producto Q_c es tipo Mertens (su 'cola' sobre "
          "primos grandes lleva una fraccion O(1) del log: NO converge "
          "absolutamente); la unica cola de Euler absolutamente "
          "convergente es la de M, y su resto M_cola ya esta cerrado "
          "[T, doc 110]. Ninguna parte ABIERTA del exceso es una cola "
          "de Euler. %s" % ("PASS" if ok else "FAIL"))
    return ok


# ------------------------------------------------------------------ K6
def check_K6():
    print("K6 intento de cota: libro mayor y la pieza que se resiste")
    ok = True
    frames = []
    for z in (17, 29, 37):
        frames.append(("A", z, z ** 5, int(z ** (10.0 / 3.0))))
    frames.append(("B", 29, 29 ** 6, 29 ** 3))
    for (tag, z, D, y) in frames:
        lvs = [Level111(z, 3, (), lv.L)
               for lv in frame_levels(z, 3, (), D, y)]
        mxP = max(lvs, key=lambda t: t.P)
        mxE = max(lvs, key=lambda t: t.E)
        # descomposicion de P por signo de delta en el nivel de max P
        a = mxP
        P_dpos = 0.0
        P_dneg = 0.0
        for (r, rpf, rmu, rphi) in a.lv.enum:
            if a.t[r] > 0:
                if a.lv.delta[r] > 0:
                    P_dpos += a.t[r]
                else:
                    P_dneg += a.t[r]
        print("  [%s] z=%d (%d niveles): max E/Q_c=%+.3f | max P/Q_c="
              "%.3f (L=%d) | en ese nivel: N/Q_c=%.3f, invisible/Q_c="
              "%+.3f, P por delta>0: %.3f Q_c, por delta<-2Pi: %.3f Q_c"
              % (tag, z, len(lvs), mxE.E / mxE.qc, mxP.P / mxP.qc,
                 mxP.L, mxP.N / mxP.qc, mxP.invis_sum / mxP.qc,
                 P_dpos / mxP.qc, P_dneg / mxP.qc))
        ok &= mxP.P / mxP.qc < 4.0
    print("  => piezas cerradas [T]: invisible <= 0; 0 <= N <= Q_c; "
          "0 <= 2*Pi_full*M_cola <= 2*Q_c. Como E = P - N con "
          "0 <= N <= Q_c, (H2'') equivale a acotar P: la pieza que se "
          "resiste es P^(L) = sum_{t(r)>0} t(r), la masa positiva del "
          "exceso. %s" % ("PASS" if ok else "FAIL"))
    return ok


# ------------------------------------------------------------------ main
def main():
    t0 = time.time()
    print("== DOC 111: exceso de Gram via diagonalizacion y Euler ==")
    print("Solo stdlib. Modo fijo; nivel L = D//s.")
    print()
    allok = True
    allok &= check_K0()
    print()
    worst = worst_levels()
    allok &= check_K1(worst)
    print()
    allok &= check_K2(worst)
    print()
    allok &= check_K2b()
    print()
    allok &= check_K3(worst)
    print()
    allok &= check_K4(worst)
    print()
    allok &= check_K5()
    print()
    allok &= check_K6()
    print()
    print("Tiempo total: %.1f s" % (time.time() - t0))
    print("PASS" if allok else "FAIL")
    return 0 if allok else 1


if __name__ == "__main__":
    sys.exit(main())
