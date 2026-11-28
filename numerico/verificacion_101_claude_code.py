#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Verificación del doc 101 (experimento E100.1 y conteo A100.1).

Procedencia: Claude Code (Anthropic, modelo claude-opus-5-5), 2026-11-28.
Método: reimplementación PROPIA desde las definiciones escritas (doc 86 §1: V+ de
Nathanson, Teorema 9.3, beta=2; doc 94 §0: soporte m = s d y H_A como conjunto de
divisores; doc 100 §5: descomposición BV y objetivo de Chen; doc 101: familia escalada,
regresiones y parte A). No importa ni ejecuta código de Muse ni de Lumina. Los valores
publicados aparecen solo como "afirmación a contrastar".
Estado: PROPUESTA (pendiente de contraste por Muse). Mediciones [M] y comprobaciones de
identidades finitas; las regresiones son ajustes sobre ventanas finitas, no cotas.
No es revisión independiente en el sentido del contrato (README, regla 4): este agente
comparte usuario de Windows con las corridas de Lumina en esta PC.
Solo biblioteca estándar. Uso:  python -B verificacion_101_claude_code.py
Nada de esto demuestra Goldbach, que sigue abierta.
"""
import math
import sys
import time
from fractions import Fraction
from itertools import compress

T0 = time.time()
CHECKS = []
fsum = math.fsum


def check(name, ok, detail=""):
    CHECKS.append((name, bool(ok)))
    print(f"[{'PASS' if ok else 'FAIL'}] {name}" + (f" | {detail}" if detail else ""))


def primes_upto(n):
    s = bytearray([1]) * (n + 1)
    s[0] = s[1] = 0
    for p in range(2, math.isqrt(n) + 1):
        if s[p]:
            s[p * p::p] = bytearray(len(range(p * p, n + 1, p)))
    return [i for i in range(n + 1) if s[i]]


def factor(n):
    f = {}
    p = 2
    while p * p <= n:
        while n % p == 0:
            f[p] = f.get(p, 0) + 1
            n //= p
        p += 1
    if n > 1:
        f[n] = f.get(n, 0) + 1
    return f


# ----------------------------------------------------------------------------
# V+ (doc 86 §1, Nathanson Thm 9.3, beta=2) con nivel L:
#   d = p1...pk libre de cuadrados, p1 > ... > pk, p1 < L, y
#   p_m^2 * (p1...pm) < L para todo m impar.
# Nivel del par (s,d): L_s = floor(D/s) (doc 94, script 94 l.63; doc 101 §0).
# Enumeración recursiva propia por primos decrecientes, con poda.
# ----------------------------------------------------------------------------
def vplus(small_desc, L):
    """Lista de (d, phi(d), primos) con d en V+(L)."""
    out = []
    n = len(small_desc)

    def rec(i0, prod, ph, ell, fs):
        out.append((prod, ph, fs))
        for i in range(i0, n):
            p = small_desc[i]
            if p >= L:
                continue
            np_ = prod * p
            l = ell + 1
            if l % 2 == 1 and p * p * np_ >= L:
                continue
            rec(i + 1, np_, ph * (p - 1), l, fs + (p,))

    if L > 1:
        rec(0, 1, 1, 0, ())
    return out


def squarefree_below(small, L):
    """Número de d libres de cuadrados con primos en `small` y d < L (envolvente sin V+)."""
    cnt = 0
    ps = sorted(small)

    def rec(i0, prod):
        nonlocal cnt
        cnt += 1
        for i in range(i0, len(ps)):
            q = prod * ps[i]
            if q >= L:
                break
            rec(i + 1, q)

    if L > 1:
        rec(0, 1)
    return cnt


class Family:
    """Modelo finito: A=3, N par (se excluyen 2 y los primos de A), sprimes en [z,y)."""

    def __init__(self, z, y, D, A=3, primes=None, real=False):
        self.z, self.y, self.D, self.A = z, y, D, A
        pr = primes if primes is not None else primes_upto(max(y, z) + 1)
        bad = set(factor(2 * A))
        self.small = sorted((p for p in pr if p < z and p not in bad), reverse=True)
        self.sprimes = [p for p in pr if z <= p < y and p not in bad]
        # real=True: nivel real D/s (sin suelo), para medir el efecto de la convención
        self.levels = {s: (Fraction(D, s) if real else D // s) for s in self.sprimes}

    def count(self):
        """Como formula(), pero solo cuenta (sin listas): para ventanas grandes."""
        small = self.small
        n = len(small)

        def rec(L, i0, prod, ph, ell, acc):
            acc[0] += 1
            acc[1] += prod / ph
            for i in range(i0, n):
                p = small[i]
                if p >= L:
                    continue
                np_ = prod * p
                l = ell + 1
                if l % 2 == 1 and p * p * np_ >= L:
                    continue
                rec(L, i + 1, np_, ph * (p - 1), l, acc)

        H = 0
        W = []
        for s in self.sprimes:
            acc = [0, 0.0]
            rec(self.levels[s], 0, 1, 1, 0, acc)
            H += acc[0]
            W.append(s / (s - 1) * acc[1])
        acc = [0, 0.0]
        rec(self.levels[self.sprimes[0]], 0, 1, 1, 0, acc)
        return H + acc[0], fsum(W) + acc[1]

    def formula(self):
        """A100.1: |H| = sum_s |V+(L_s)| + |V+(L_min)|;
        W_1 = sum_s s/(s-1) Phi+(L_s) + Phi+(L_min), Phi+(L) = sum_{d in V+(L)} d/phi(d)."""
        H = 0
        W = []
        for s in self.sprimes:
            V = vplus(self.small, self.levels[s])
            H += len(V)
            W.append(s / (s - 1) * fsum(d / ph for d, ph, _ in V))
        V0 = vplus(self.small, self.levels[self.sprimes[0]])
        H += len(V0)
        W.append(fsum(d / ph for d, ph, _ in V0))
        return H, fsum(W)

    def envelope(self):
        tot = sum(squarefree_below(self.small, self.levels[s]) for s in self.sprimes)
        return tot + squarefree_below(self.small, self.levels[self.sprimes[0]])

    def enumerate_H(self):
        """Definición del doc 94: H_A = {h : h | m para algún m = s d del soporte}.
        Devuelve dict h -> (phi(h), primos de h)."""
        H = {}
        for s in self.sprimes:
            for d, ph, fs in vplus(self.small, self.levels[s]):
                fac = (s,) + fs
                for k in range(1 << len(fac)):
                    hs = tuple(fac[i] for i in range(len(fac)) if k >> i & 1)
                    h = math.prod(hs)
                    if h not in H:
                        H[h] = (math.prod(q - 1 for q in hs), hs)
        return H

    def closure_monotone(self):
        """Comprueba en TODOS los niveles: V+(L) cerrado por divisores y V+(L_s) ⊆ V+(L_min)."""
        top = {d for d, _, _ in vplus(self.small, self.levels[self.sprimes[0]])}
        ok = True
        for s in self.sprimes:
            V = vplus(self.small, self.levels[s])
            S = {d for d, _, _ in V}
            for d, _, fs in V:
                if d not in top or any(d // p not in S for p in fs):
                    ok = False
                    break
            if not ok:
                break
        return ok, len(self.sprimes)


# ----------------------------------------------------------------------------
# Estadística propia: MCO, error estándar, cuantil t de Student por integración
# numérica (Simpson) de la densidad, IC95% y rango leave-one-out.
# ----------------------------------------------------------------------------
def t_cdf(t, df, n=20000):
    c = math.exp(math.lgamma((df + 1) / 2) - math.lgamma(df / 2)) / math.sqrt(df * math.pi)
    f = lambda x: c * (1 + x * x / df) ** (-(df + 1) / 2)
    h = t / n
    s = f(0) + f(t) + 4 * fsum(f((2 * i - 1) * h) for i in range(1, n // 2 + 1)) \
        + 2 * fsum(f(2 * i * h) for i in range(1, n // 2))
    return 0.5 + s * h / 3


def t_quantile(p, df):
    lo, hi = 0.0, 50.0
    for _ in range(80):
        mid = (lo + hi) / 2
        if t_cdf(mid, df) < p:
            lo = mid
        else:
            hi = mid
    return (lo + hi) / 2


def ols(xs, ys):
    n = len(xs)
    mx, my = fsum(xs) / n, fsum(ys) / n
    sxx = fsum((x - mx) ** 2 for x in xs)
    b = fsum((x - mx) * (y - my) for x, y in zip(xs, ys)) / sxx
    a = my - b * mx
    se = math.sqrt(fsum((y - a - b * x) ** 2 for x, y in zip(xs, ys)) / (n - 2) / sxx) if n > 2 else float('nan')
    return a, b, se


def fit(xs, ys):
    a, b, se = ols(xs, ys)
    n = len(xs)
    tq = t_quantile(0.975, n - 2) if n > 2 else float('nan')
    loo = [ols(xs[:i] + xs[i + 1:], ys[:i] + ys[i + 1:])[1] for i in range(n)] if n > 2 else [b]
    return dict(a=a, b=b, se=se, t=tq, ci=tq * se, loo=(min(loo), max(loo)), n=n)


# ----------------------------------------------------------------------------
# Afirmaciones a contrastar (doc 101)
# ----------------------------------------------------------------------------
DOC101_R = {"R1": (54, 68.264353), "R2": (880, 1011.279168), "R3": (68, 69.307280)}
DOC101_FAM = {11: (422, 802, 912.79, 1692), 13: (682, 1628, 1898.02, 4893), 17: (1502, 5033, 6015.88, 17757),
              19: (2090, 8245, 9943.43, 36683), 23: (3687, 17999, 21933.30, 94940), 29: (7381, 47178, 58126.02, 274982)}
DOC101_FIT = {"W": (4.2886, 0.0237, 0.0658, 4.2650, 4.3203), "H": (4.2079, 0.0170, 0.0472, 4.1969, 4.2279)}
DOC101_LOCAL = (0.438, 0.430, 0.452, 0.414, 0.420)
DOC101_PARTA = {  # N: (|H_N|, W_max, S_BV, B_M, delta)
    510510: (424, 424.95, 9.909969e8, 9.909969e8, 0.3483),
    999982: (876, 1006.34, 4.139329e9, 3.777086e9, 0.3420),
    1048576: (880, 1011.28, 4.813446e9, 4.078042e9, 0.3763),
    9699690: (423, 423.90, 1.954892e10, 1.954892e10, 0.2962),
    16777216: (880, 1011.28, 6.737052e10, 5.692229e10, 0.2743),
    67108864: (880, 1011.28, 2.793088e11, 2.468347e11, 0.2624),
    77597520: (423, 423.90, 1.524628e11, 1.524628e11, 0.2557),
    99999998: (672, 736.00, 3.210471e11, 2.946279e11, 0.2594)}
DOC101_BETA = {"S": (0.9904, 0.0716, 0.1752, 0.9214, 1.0286), "B": (0.9904, 0.0588, 0.1439, 0.9327, 1.0215)}
C_STAR = 0.039884   # doc 100 §5.3
G_BIG = 0.0225      # doc 100 §5.3 (f·G_A ≈ 0,05·0,45)


def part_F():
    print("\n=== PARTE F — A100.1: V+, cierre, fórmula de |H| y W_1 = W_max ===")
    print("Definición V+: Nathanson Thm 9.3 (doc 86 §1), beta=2, nivel entero L_s = floor(D/s).")
    models = {"R1": Family(11, 60, 20000), "R2": Family(11, 3000, 200000), "R3": Family(3, 349, 200000)}
    fam = {z: Family(z, round(z ** (10 / 3)), z ** 5) for z in (11, 13, 17, 19, 23, 29)}
    allok = True
    for name, F in list(models.items()) + [(f"z={z}", F) for z, F in fam.items()]:
        ok, nlev = F.closure_monotone()
        allok &= ok
    check("F1 V+(L) cerrado por divisores y V+(L_s) ⊆ V+(L_min) en TODOS los niveles (R1-R3, z=11..29)", allok)
    res = {}
    for name, F in list(models.items()) + [(z, F) for z, F in fam.items()]:
        Hf, Wf = F.formula()
        He = F.enumerate_H()
        We = fsum(h / ph for h, (ph, _) in He.items())
        # Lema 1: W_d = sum_{r in H, d|r} r/phi(r); argmax
        Wd = {}
        for r, (ph, hs) in He.items():
            w = r / ph
            for k in range(1 << len(hs)):
                d = math.prod(hs[i] for i in range(len(hs)) if k >> i & 1)
                Wd[d] = Wd.get(d, 0.0) + w
        dmax = max(Wd, key=Wd.get)
        res[name] = dict(Hf=Hf, He=len(He), Wf=Wf, We=We, dmax=dmax, Wmax=Wd[dmax], ns=len(F.sprimes),
                         env=F.envelope())
        print(f"  {str(name):>4}: #sprimes={len(F.sprimes):5d} |H| fórmula={Hf:6d} enumeración(doc 94)={len(He):6d} "
              f"W_1 fórmula={Wf:.6f} enumeración={We:.6f} argmax W_d = {dmax} env sin V+={res[name]['env']}")
    check("F2 |H|: fórmula A100.1 = enumeración de divisores (definición del doc 94) en los 9 modelos",
          all(r['Hf'] == r['He'] for r in res.values()))
    e = max(abs(r['Wf'] - r['We']) / r['We'] for r in res.values())
    check("F3 W_1: fórmula = enumeración", e < 1e-12, f"max err rel {e:.1e}")
    check("F4 Lema 1: argmax_d W_d = 1 y W_max = W_1 en los 9 modelos",
          all(r['dmax'] == 1 and abs(r['Wmax'] - r['We']) < 1e-9 * r['We'] for r in res.values()))
    dR = max(max(abs(res[k]['Hf'] - v[0]), abs(res[k]['Wf'] - v[1]) * 1e6) for k, v in DOC101_R.items())
    check("F5 R1/R2/R3 del doc 101 §2.3 (|H| y W_1 a 6 decimales)", dR < 1.0)
    okf = all(res[z]['ns'] == v[0] and res[z]['Hf'] == v[1] and abs(res[z]['Wf'] - v[2]) < 0.006 and res[z]['env'] == v[3]
              for z, v in DOC101_FAM.items())
    check("F6 Tabla del doc 101 §2.4 (#sprimes, |H|, W_1, envolvente) en z=11..29", okf)
    # Lema de inclusión (base de la prueba de α = 1/2, §3 del informe): para delta>0,
    # todo d libre de cuadrados con primos < L^delta y d < L^(1-3 delta) está en V+(L).
    # Y el caso de pares: p2 < p1 con p1^3 < L  =>  p1 p2 en V+(L).
    inc_ok = True
    nlev = 0
    for z, F in fam.items():
        ps = sorted(F.small)
        for s in F.sprimes[::7]:
            L = F.levels[s]
            S = {d for d, _, _ in vplus(F.small, L)}
            nlev += 1
            for delta in (0.10, 0.15, 0.20, 0.25):
                cap_p = L ** delta
                cap_d = L ** (1 - 3 * delta)

                def rec(i0, prod):
                    global_ok = True
                    if prod not in S:
                        return False
                    for i in range(i0, len(ps)):
                        if ps[i] >= cap_p or prod * ps[i] >= cap_d:
                            break
                        if not rec(i + 1, prod * ps[i]):
                            return False
                    return global_ok

                inc_ok &= rec(0, 1)
            for i, p1 in enumerate(ps):
                if p1 ** 3 >= L:
                    break
                inc_ok &= all(p1 * p2 in S for p2 in ps[:i])
    check("F7 Lema de inclusión: V+(L) ⊇ {d sin cuadrados, primos < L^δ, d < L^(1-3δ)} (δ=0.10..0.25) y pares p1^3<L",
          inc_ok, f"{nlev} niveles de z=11..29")
    # Variante con nivel real D/s (sin suelo): efecto de la convención de redondeo
    real = {z: Family(z, round(z ** (10 / 3)), z ** 5, real=True).count() for z in fam}
    print("  variante con nivel real D/s (sin suelo): |H| =", {z: v[0] for z, v in real.items()})
    regressions(list(real), [v[1] for v in real.values()], "W_1 nivel real")
    regressions(list(real), [v[0] for v in real.values()], "|H| nivel real")
    return res


def regressions(zs, vals, label):
    xs = [math.log(z) for z in zs]
    ys = [math.log(v) for v in vals]
    f = fit(xs, ys)
    loc = [(ys[i + 1] - ys[i]) / (xs[i + 1] - xs[i]) / 10 for i in range(len(xs) - 1)]
    print(f"  {label}: pendiente={f['b']:.4f} ± {f['se']:.4f} (SE), t(0.975,{f['n']-2})={f['t']:.4f}, "
          f"IC95% ±{f['ci']:.4f}, LOO [{f['loo'][0]:.4f}, {f['loo'][1]:.4f}] ⟹ α={f['b']/10:.4f}, "
          f"IC95% [{(f['b']-f['ci'])/10:.4f}, {(f['b']+f['ci'])/10:.4f}], LOO [{f['loo'][0]/10:.4f}, {f['loo'][1]/10:.4f}]; "
          f"prefactor (N=z^10) = {math.exp(f['a']):.4f}; α locales: {', '.join(f'{x:.3f}' for x in loc)}")
    return f, loc


def part_G(res):
    print("\n=== PARTE G — Regresiones α_W, α_H (ventana z=11..29) y extensión de la ventana ===")
    zs = [11, 13, 17, 19, 23, 29]
    tq4 = t_quantile(0.975, 4)
    tq6 = t_quantile(0.975, 6)
    print(f"  cuantiles t propios: t(0.975,4)={tq4:.6f}, t(0.975,6)={tq6:.6f} (tablas: 2.776445, 2.446912)")
    check("G0 cuantiles t por integración numérica = tablas", abs(tq4 - 2.776445) < 2e-6 and abs(tq6 - 2.446912) < 2e-6)
    fW, locW = regressions(zs, [res[z]['Wf'] for z in zs], "W_1  z=11..29")
    fH, _ = regressions(zs, [res[z]['Hf'] for z in zs], "|H|  z=11..29")
    fE, _ = regressions(zs, [res[z]['env'] for z in zs], "env  z=11..29")
    dev = 0.0
    for f, key in ((fW, "W"), (fH, "H")):
        b, se, ci, l0, l1 = DOC101_FIT[key]
        dev = max(dev, abs(f['b'] - b), abs(f['se'] - se), abs(f['ci'] - ci), abs(f['loo'][0] - l0), abs(f['loo'][1] - l1))
    check("G1 Pendientes, SE, IC95% y LOO de W_1 y |H| reproducen el doc 101 §2.4", dev < 6e-5, f"max desviación {dev:.1e}")
    check("G2 Envolvente sin V+: pendiente 5.2439 (α_env 0.5244)", abs(fE['b'] - 5.2439) < 6e-5)
    check("G3 α locales de W_1 reproducen 0.438, 0.430, 0.452, 0.414, 0.420",
          max(abs(a - b) for a, b in zip(locW, DOC101_LOCAL)) < 6e-4)
    check("G4 Veredicto en la ventana del doc 101: IC95% de α_W entero > 0.42",
          (fW['b'] - fW['ci']) / 10 > 0.42, f"[{(fW['b']-fW['ci'])/10:.4f}, {(fW['b']+fW['ci'])/10:.4f}]")
    # Extensión (cómputo propio): z primos hasta 59, por la fórmula A100.1 (válida por F1 + prueba de cierre)
    ext = [31, 37, 41, 43, 47, 53, 59]
    pr = primes_upto(round(max(ext) ** (10 / 3)) + 2)
    data = {z: (res[z]['Hf'], res[z]['Wf']) for z in zs}
    for z in ext:
        t1 = time.time()
        F = Family(z, round(z ** (10 / 3)), z ** 5, primes=pr)
        Hf, Wf = F.formula()
        data[z] = (Hf, Wf)
        print(f"  extensión z={z}: D={z**5} y={F.y} #sprimes={len(F.sprimes)} |H|={Hf} W_1={Wf:.2f} ({time.time()-t1:.1f} s)")
    allz = zs + ext
    print("  -- ajustes en ventanas ampliadas (MEDICIÓN; puntos deterministas, el 'IC' mide dispersión del ajuste) --")
    out = {}
    for lab, win in (("11..37", [z for z in allz if z <= 37]), ("11..59", allz),
                     ("19..59", [z for z in allz if z >= 19]), ("29..59", [z for z in allz if z >= 29])):
        f1, l1 = regressions(win, [data[z][1] for z in win], f"W_1 z={lab}")
        f2, _ = regressions(win, [data[z][0] for z in win], f"|H| z={lab}")
        out[lab] = (f1, f2, l1)
    f37 = out["11..37"][0]
    check("G5 z=11..37 reproduce el doc 102 (α_W=0.4281, IC95% [0.4243, 0.4318], LOO [0.4264, 0.4291])",
          abs(f37['b'] / 10 - 0.4281) < 6e-5 and abs((f37['b'] - f37['ci']) / 10 - 0.4243) < 6e-5
          and abs((f37['b'] + f37['ci']) / 10 - 0.4318) < 6e-5)
    # Ventana larga (solo conteo): z primos 59..151. Ilustra que la ventana del doc 101 es preasintótica.
    big = [59, 61, 67, 71, 73, 79, 83, 97, 113, 131, 151]
    pr2 = primes_upto(round(max(big) ** (10 / 3)) + 2)
    print("  -- ventana larga (solo conteo, cómputo propio) --")
    for z in big:
        if z in data:
            continue
        t1 = time.time()
        F = Family(z, round(z ** (10 / 3)), z ** 5, primes=pr2)
        data[z] = F.count()
        print(f"  z={z}: primos pequeños={len(F.small)} #sprimes={len(F.sprimes)} |H|={data[z][0]} W_1={data[z][1]:.2f} "
              f"(cota de saturación 2^{len(F.small)}={2**len(F.small)} por nivel; {time.time()-t1:.1f} s)")
    fL, locL = regressions(big, [data[z][1] for z in big], "W_1 z=59..151")
    regressions(big, [data[z][0] for z in big], "|H| z=59..151")
    f79, _ = regressions([z for z in big if z >= 79], [data[z][1] for z in big if z >= 79], "W_1 z=79..151")
    check("G6 En z=79..151 el α_W ajustado supera la medida del doc 101 (0.4289) y el α local medio > 0.44",
          f79['b'] / 10 > 0.4289 and fsum(locL[5:]) / len(locL[5:]) > 0.44,
          f"α_W(79..151)={f79['b']/10:.4f}; locales {', '.join(f'{x:.3f}' for x in locL[5:])}")
    return fW, fH, data, out


def part_C0():
    print("\n=== PARTE G' — Sensibilidad a la parametrización del doc 94 (P-rev 101.3) ===")
    print("  doc 94: D_A = X^(1/2)/(log X)^C0; la familia del doc 101 usa C0 = 0 (D = z^5). Con X = z^10:")
    zs = [11, 13, 17, 19, 23, 29]
    for C0 in (1, 2):
        Hs, Ws = [], []
        for z in zs:
            D = int(z ** 5 / math.log(z ** 10) ** C0)
            F = Family(z, round(z ** (10 / 3)), D)
            h, w = F.formula()
            Hs.append(h)
            Ws.append(w)
        print(f"  C0={C0}: |H| = {Hs}")
        regressions(zs, Ws, f"W_1 C0={C0}")
        regressions(zs, Hs, f"|H| C0={C0}")
    print("  Cambiar N = z^10 por N = A·z^10 (X = N/A como en el doc 94) solo desplaza ln N: la pendiente en ln N no cambia [T].")


# ----------------------------------------------------------------------------
# PARTE H — E100.1 parte A: soporte fijo R2, N par creciente, psi real hasta 1e8
# ----------------------------------------------------------------------------
def part_H():
    print("\n=== PARTE H — doc 101 §3.2-3.3: soporte fijo R2, psi real hasta 10^8 ===")
    LIM = 10 ** 8
    t1 = time.time()
    sieve = bytearray([1]) * (LIM + 1)
    sieve[0] = sieve[1] = 0
    for p in range(2, math.isqrt(LIM) + 1):
        if sieve[p]:
            sieve[p * p::p] = bytearray(len(range(p * p, LIM + 1, p)))
    ppow = []  # potencias p^k, k>=2
    for p in range(2, math.isqrt(LIM) + 1):
        if sieve[p]:
            x = p * p
            while x <= LIM:
                ppow.append((x, math.log(p)))
                x *= p
    print(f"  criba hasta 1e8: {time.time()-t1:.1f} s; potencias de primos (k>=2): {len(ppow)}")

    def psi_prog(N, d, a):
        a %= d
        if a == 0:
            a = d
        tot = fsum(map(math.log, compress(range(a, N + 1, d), sieve[a:N + 1:d])))
        return tot + fsum(w for x, w in ppow if x <= N and x % d == a % d)

    psi999983 = psi_prog(999983, 1, 0)
    check("H0 psi(999983) = 999586.5975 (doc 97)", abs(psi999983 - 999586.5975) < 6e-5, f"{psi999983:.4f}")
    R2 = Family(11, 3000, 200000)
    H = R2.enumerate_H()      # dict r -> (phi, primos)
    rows = {}
    for N, fam in ((510510, "primorial ..17"), (999982, "2·79·6329"), (1048576, "2^20"), (9699690, "primorial ..19"),
                   (16777216, "2^24"), (67108864, "2^26"), (77597520, "2^4·primorial 19"), (99999998, "2·7·23·310559")):
        t1 = time.time()
        HN = {r: v for r, v in H.items() if math.gcd(r, N) == 1}
        psiN = psi_prog(N, 1, 0)
        E = {d: (psi_prog(N, d, N) - psiN / HN[d][0]) if d > 1 else 0.0 for d in HN}
        closed = all(math.prod(hs[i] for i in range(len(hs)) if k >> i & 1) in HN
                     for r, (ph, hs) in HN.items() for k in range(1 << len(hs)))
        SBV = fsum(HN[d][0] * E[d] ** 2 for d in HN)
        R0 = math.log(N)
        BM_terms = []
        nprime_only = all(len(hs) <= 1 for _, (ph, hs) in HN.items())
        for r, (ph, hs) in HN.items():
            if r <= R0:
                continue
            M = 0.0
            for k in range(1 << len(hs)):
                sub = [hs[i] for i in range(len(hs)) if k >> i & 1]
                d = math.prod(sub)
                mu = -1 if (len(hs) - len(sub)) % 2 else 1
                M += mu * HN[d][0] * E[d]
            BM_terms.append(M * M / ph)
        BM = fsum(BM_terms)
        Wd = {}
        for r, (ph, hs) in HN.items():
            for k in range(1 << len(hs)):
                d = math.prod(hs[i] for i in range(len(hs)) if k >> i & 1)
                Wd[d] = Wd.get(d, 0.0) + r / ph
        dmax = max(Wd, key=Wd.get)
        delta = SBV / (len(HN) * N * math.log(N))
        rows[N] = dict(H=len(HN), W=Wd[dmax], dmax=dmax, SBV=SBV, BM=BM, delta=delta, closed=closed,
                       primes_only=nprime_only, psiN=psiN)
        print(f"  N={N:>9} ({fam}): |H_N|={len(HN)} cerrado={closed} W_max={Wd[dmax]:.2f} (argmax d={dmax}) "
              f"S_BV={SBV:.6e} B_M={BM:.6e} δ={delta:.4f} S_BV/(N^1.5 lnN)={SBV/(N**1.5*math.log(N)):.3e} "
              f"psi/N={psiN/N:.4f} ({time.time()-t1:.1f} s)")
    dev = 0.0
    for N, (h, w, s, b, dl) in DOC101_PARTA.items():
        r = rows[N]
        dev = max(dev, abs(r['SBV'] / s - 1), abs(r['BM'] / b - 1))
        if r['H'] != h or abs(r['W'] - w) > 0.006 or abs(r['delta'] - dl) > 6e-5:
            dev = float('inf')
    check("H1 Tabla del doc 101 §3.2 (|H_N|, W_max, S_BV, B_M, δ) en los 8 N", dev < 2e-6, f"max desviación relativa {dev:.1e}")
    check("H2 Lema 1 en los 8 N: argmax_d W_d = 1", all(r['dmax'] == 1 for r in rows.values()))
    prim = [N for N in rows if rows[N]['primes_only']]
    check("H3 Colapso primorial: si H_N solo tiene 1 y primos, B_M = S_BV (identidad M_r = phi(r) E_r)",
          prim == [510510, 9699690, 77597520] and all(abs(rows[N]['BM'] / rows[N]['SBV'] - 1) < 1e-12 for N in prim),
          f"N con H_N primo: {prim}")
    Ns = sorted(rows)
    for key, lab in (("SBV", "S"), ("BM", "B")):
        f = fit([math.log(N) for N in Ns], [math.log(rows[N][key]) for N in Ns])
        b, se, ci, l0, l1 = DOC101_BETA[lab]
        d = max(abs(f['b'] - b), abs(f['se'] - se), abs(f['ci'] - ci), abs(f['loo'][0] - l0), abs(f['loo'][1] - l1))
        print(f"  β({key}) = {f['b']:.4f} ± {f['se']:.4f} (SE), IC95% ±{f['ci']:.4f}, LOO [{f['loo'][0]:.4f}, {f['loo'][1]:.4f}]")
        check(f"H4{lab} β de {key} reproduce el doc 101 §3.2", d < 6e-5, f"max desviación {d:.1e}")
    for lab, sub in (("2^e", [1048576, 16777216, 67108864]), ("primoriales", [510510, 9699690, 77597520])):
        f = fit([math.log(N) for N in sub], [math.log(rows[N]['SBV']) for N in sub])
        print(f"  subfamilia {lab}: β_S_BV = {f['b']:.4f}, LOO [{f['loo'][0]:.4f}, {f['loo'][1]:.4f}]")
    return rows


# ----------------------------------------------------------------------------
# PARTE I — §2.5 y §3.4: extrapolación contra el objetivo de Chen del doc 100
# ----------------------------------------------------------------------------
def part_I(fW, fH, rows):
    print("\n=== PARTE I — doc 101 §2.5 y §3.4: aritmética de la extrapolación (definiciones del doc 100 §5.3) ===")
    lnN = lambda e: e * math.log(10)
    # Objetivo doc 100: B_M << N^2/(A^2 ln^4 N · 0.0225), A = N^c.  Cota BV: W_max · N^{3/2} ln N.
    weq = lambda e: math.exp(e * math.log(10) * (0.5 - 2 * C_STAR)) / (G_BIG * lnN(e) ** 5)
    for e, claim101, claim100 in ((6, 644, None), (12, 853, 0.30), (24, 1051, 1.05e3), (100, None, 7.24e31)):
        print(f"  W_eq(10^{e}) = {weq(e):.4g}   (doc 100 §5.3: {claim100}; doc 101 §2.5: {claim101})")
    check("I1 W_eq del doc 100 §5.3 se reproduce (0.30 en 10^12, 1.05e3 en 10^24, 7.24e31 en 10^100)",
          abs(weq(12) - 0.30) < 0.01 and abs(weq(24) / 1.05e3 - 1) < 0.01 and abs(weq(100) / 7.24e31 - 1) < 0.01)
    check("I2 Los W_eq del doc 101 §2.5 (≈644 en 10^6, ≈853 en 10^12) NO salen de la definición del doc 100",
          abs(weq(12) - 853) > 100 and abs(weq(6) - 644) > 100, f"reales: {weq(6):.3g} y {weq(12):.3g}")
    cW = math.exp(fW['a'])
    aW = fW['b'] / 10
    Wfit = lambda e: cW * math.exp(lnN(e) * aW)
    for e in (12, 24):
        print(f"  W_max ajustado (doc 101: 0.031·N^0.4289) en 10^{e} = {Wfit(e):.4g}; falla de la cota BV = W_max/W_eq = {Wfit(e)/weq(e):.3g}")
    check("I3 'Falla por ~5×' en 10^12 (doc 101 §2.5) NO se sostiene con W_eq del doc 100: la falla es ~10^4",
          Wfit(12) / weq(12) > 1e4, f"{Wfit(12)/weq(12):.3g}×")
    # §3.4 ruta combinada: cota = W_max(N) · δ |H|(N) N ln N, con δ ≈ 0.26 (doc 101)
    cH = math.exp(fH['a'])
    aH = fH['b'] / 10
    delta = 0.26
    pref = cW * cH * delta
    ratio = lambda e: pref * G_BIG * lnN(e) ** 5 * math.exp(lnN(e) * (1 + aW + aH - (2 - 2 * C_STAR)))
    print(f"  exponente total 1+α_W+α_H = {1+aW+aH:.4f} frente a 2−2c = {2-2*C_STAR:.6f}")
    print(f"  cota/objetivo = {pref*G_BIG:.3g} · ln^5 N / N^{(2-2*C_STAR)-(1+aW+aH):.4f}   (doc 101 §3.4: 0,0119 · ln^5 N / N^0,0705)")
    print(f"  doc 101 equivale a dividir por 0.0225 en vez de multiplicar: {pref/G_BIG:.4f}")
    lo, hi = 10.0, 1000.0
    for _ in range(200):
        mid = (lo + hi) / 2
        if ratio(mid) > 1:
            lo = mid
        else:
            hi = mid
    print(f"  cota/objetivo en 10^24 = {ratio(24):.3g} (doc 101: ≈1,2·10^5); baja de 1 hacia N ≈ 10^{lo:.0f} (doc 101: ~10^150)")
    check("I4 Constante de §3.4: 0.0119 = c_W·c_H·δ/0.0225; la correcta (×0.0225) es ~2000 veces menor",
          abs(pref / G_BIG - 0.0119) < 0.0008 and pref * G_BIG < 1e-5, f"correcta {pref*G_BIG:.3g}")
    # Comprobación directa con los números del doc 100 (objetivo en 10^24 = 5.803e40)
    e = 24
    cota = Wfit(e) * delta * cH * math.exp(lnN(e) * aH) * 10 ** e * lnN(e)
    obj = 10 ** (e * (2 - 2 * C_STAR)) / (G_BIG * lnN(e) ** 4)
    print(f"  directo en 10^24: cota={cota:.3e}, objetivo={obj:.4e} (doc 100: 5,803·10^40), cociente={cota/obj:.3g}")
    check("I5 Objetivo de Chen en 10^24 = 5.803e40 (doc 100 §5.3)", abs(obj / 5.803e40 - 1) < 1e-3)


def main():
    try:
        sys.stdout.reconfigure(encoding="utf-8")
    except (AttributeError, ValueError):
        pass
    print(f"Python {sys.version.split()[0]}")
    res = part_F()
    fW, fH, data, out = part_G(res)
    part_C0()
    rows = part_H()
    part_I(fW, fH, rows)
    npass = sum(ok for _, ok in CHECKS)
    print(f"\nRESUMEN: {npass}/{len(CHECKS)} PASS, {len(CHECKS) - npass} FAIL. Tiempo {time.time() - T0:.1f} s")
    return 0 if npass == len(CHECKS) else 1


if __name__ == "__main__":
    sys.exit(main())
