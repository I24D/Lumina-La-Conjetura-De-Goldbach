#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Medición del exponente de conteo de 𝒟_S (rama conmutada de D-NP, docs 115/120b/121).

Procedencia: Claude Code (Anthropic, modelo claude-opus-5-5), 2026-11-28.
Método: reimplementación PROPIA desde las definiciones escritas:
  - V+ de Nathanson (doc 86 §1), beta = 2, nivel L_q = floor(D/q) (doc 94; doc 101 §0);
  - familia del doc 101 §2.4: A = 3, D = z^5, y = round(z^(10/3)), N par (se excluyen 2 y 3);
  - Ω_d^C = Σ_{r∈ℋ_C, d|r} (τ(r)−1)/φ(r)  (doc 106 (106.2), doc 120b §1);
  - Ω_{qm}^C = Θ̃_m(q)/(q−1), Θ̃_m(q) = Σ_{h∈V+(D/q)∖{1}, m|h} (2τ(h)−1)/φ(h)  (doc 115 §7);
  - 𝒟_S = {d = q·m : Ω_d^C > 0} y W_d = Ω_d^C·φ(d)  (doc 121 §2, doc 121b §3).
No importa ni ejecuta código de nadie.
Algoritmo de conteo (propio, exacto): d ∈ V+(L) ⟺ L > T(d), T(d) = máx_{m impar} p_m²·(p1⋯pm)
(T(1) = 1). Luego Σ_q |V+(L_q)| = Σ_{d∈V+(L_max)} #{q ∈ 𝒮 : q ≤ ⌊D/(T(d)+1)⌋}; basta enumerar
V+ de UN nivel y contar primos (criba segmentada). Validado contra enumeración directa.
Estado: PROPUESTA (pendiente de contraste por Muse). Mediciones [M]; los ajustes son sobre ventanas
finitas y NO demuestran un exponente asintótico (para eso, el teorema del informe).
No es revisión independiente según el contrato (README, regla 4): mismo usuario de Windows.
Solo biblioteca estándar. Uso:  python -B verificacion_conteo_ds_claude_code.py
Nada de esto demuestra Goldbach, que sigue abierta.
"""
import math
import sys
import time
from array import array
from itertools import compress

T0 = time.time()
CHECKS = []
fsum = math.fsum
C_STAR = 0.039884
THRESH = 0.5 - 2 * C_STAR          # 0.420232, umbral exigido (doc 100 §5.3, doc 121 §3)


def check(name, ok, detail=""):
    CHECKS.append((name, bool(ok)))
    print(f"[{'PASS' if ok else 'FAIL'}] {name}" + (f" | {detail}" if detail else ""), flush=True)


def small_sieve(n):
    s = bytearray([1]) * (n + 1)
    s[0] = s[1] = 0
    for p in range(2, math.isqrt(n) + 1):
        if s[p]:
            s[p * p::p] = bytearray(len(range(p * p, n + 1, p)))
    return [i for i in compress(range(n + 1), s)]


def prime_counts(queries, lo, hi, seg=1 << 26):
    """Para una lista ORDENADA de enteros x, devuelve #{primos p : lo <= p <= min(x, hi)}
    (0 si x < lo). Criba segmentada propia."""
    base = small_sieve(math.isqrt(hi) + 1)
    out = [0] * len(queries)
    qi = 0
    nq = len(queries)
    while qi < nq and queries[qi] < lo:
        qi += 1
    acc = 0
    a = lo
    while a <= hi and qi < nq:
        b = min(a + seg, hi + 1)
        s = bytearray([1]) * (b - a)
        for p in base:
            if p * p >= b:
                break
            st = max(p * p, (a + p - 1) // p * p)
            if st < b:
                s[st - a::p] = bytearray(len(range(st, b, p)))
        for v in range(a, min(b, 2)):
            s[v - a] = 0
        prev = 0
        while qi < nq and min(queries[qi], hi) < b:
            x = min(queries[qi], hi)
            idx = x - a + 1
            acc += s.count(1, prev, idx)
            prev = idx
            out[qi] = acc
            qi += 1
        acc += s.count(1, prev, b - a)
        a = b
    while qi < nq:
        out[qi] = acc
        qi += 1
    return out


def enumerate_level(small_desc, Lmax):
    """Enumera V+(Lmax) devolviendo, por elemento: T(d), d, phi(d), omega(d)."""
    T_ = array('Q')
    w_ = array('d')
    n = len(small_desc)
    sys.setrecursionlimit(100000)

    def rec(i, prod, ph, ell, T, om):
        T_.append(T)
        # peso de W_d (doc 121): (2τ(h)−1)·h/φ(h) para h > 1;  0 para h = 1
        w_.append(((1 << (om + 1)) - 1) * prod / ph if prod > 1 else 0.0)
        for k in range(i, n):
            p = small_desc[k]
            np_ = prod * p
            l = ell + 1
            Tn = T
            if l & 1:
                t = p * p * np_
                if t >= Lmax:
                    continue
                if t > Tn:
                    Tn = t
            rec(k + 1, np_, ph * (p - 1), l, Tn, om + 1)

    rec(0, 1, 1, 0, 1, 0)
    return T_, w_


def point(z, D=None, y=None):
    """|H|, |𝒟_S|, Σ_{d∈𝒟_S} W_d para la familia del doc 101 en el z dado."""
    D = z ** 5 if D is None else D
    y = round(z ** (10 / 3)) if y is None else y
    small = sorted((p for p in small_sieve(z) if p < z and p not in (2, 3)), reverse=True)
    # primer sprime >= z
    q0 = z
    while any(q0 % p == 0 for p in range(2, math.isqrt(q0) + 1)) or q0 in (2, 3):
        q0 += 1
    Lmax = D // q0
    T_, w_ = enumerate_level(small, Lmax)
    xs = [D // (t + 1) for t in T_]
    order = sorted(range(len(xs)), key=xs.__getitem__)
    qs = [xs[i] for i in order] + []
    extra = sorted([D // 2, D // 126, y - 1])
    cnt = prime_counts(qs, z, y - 1)
    cnt_extra = dict(zip(extra, prime_counts(extra, z, y - 1)))
    sum_pi = sum(cnt)
    sum_w = fsum(w_[i] * c for i, c in zip(order, cnt))
    nS = cnt_extra[y - 1]
    ones = cnt_extra[D // 2] - cnt_extra[D // 126]     # q con V+(L_q) = {1}
    H = sum_pi + len(T_)
    DS = sum_pi - ones
    return dict(z=z, small=len(small), nS=nS, nV=len(T_), H=H, DS=DS, W=sum_w)


# ---------------------------------------------------------------------------
# Validaciones independientes
# ---------------------------------------------------------------------------
def vplus_direct(small_desc, L):
    out = []

    def rec(i0, prod, ph, ell, fs):
        out.append((prod, ph, fs))
        for i in range(i0, len(small_desc)):
            p = small_desc[i]
            if p >= L:
                continue
            np_ = prod * p
            if (ell + 1) % 2 == 1 and p * p * np_ >= L:
                continue
            rec(i + 1, np_, ph * (p - 1), ell + 1, fs + (p,))

    if L > 1:
        rec(0, 1, 1, 0, ())
    return out


def direct_definitions(z):
    """Enumeración literal: ℋ = divisores de q·d (doc 94), ℋ_C, Ω_d^C por (106.2), 𝒟_S, W_d."""
    D = z ** 5
    y = round(z ** (10 / 3))
    pr = small_sieve(y)
    small = sorted((p for p in pr if p < z and p not in (2, 3)), reverse=True)
    sw = [p for p in pr if z <= p < y and p not in (2, 3)]
    Hset = {}
    for q in sw:
        for d, ph, fs in vplus_direct(small, D // q):
            fac = (q,) + fs
            for k in range(1 << len(fac)):
                hs = tuple(fac[i] for i in range(len(fac)) if k >> i & 1)
                Hset[math.prod(hs)] = hs
    Om = {}
    for r, hs in Hset.items():
        if len(hs) < 2:
            continue                      # ℋ_C: compuestos
        tau = 1 << len(hs)
        val = (tau - 1) / math.prod(p - 1 for p in hs)
        for k in range(1, 1 << len(hs)):
            sub = tuple(hs[i] for i in range(len(hs)) if k >> i & 1)
            dd = math.prod(sub)
            Om[dd] = Om.get(dd, 0.0) + val
    DS = {dd: v for dd, v in Om.items() if v > 0 and max(Hset[dd]) >= z}
    W = fsum(v * math.prod(p - 1 for p in Hset[dd]) for dd, v in DS.items())
    return dict(H=len(Hset), DS=len(DS), W=W)


# ---------------------------------------------------------------------------
# Regresión propia
# ---------------------------------------------------------------------------
def t_quantile(p, df):
    c = math.exp(math.lgamma((df + 1) / 2) - math.lgamma(df / 2)) / math.sqrt(df * math.pi)
    f = lambda x: c * (1 + x * x / df) ** (-(df + 1) / 2)

    def cdf(t, n=20000):
        h = t / n
        s = f(0) + f(t) + 4 * fsum(f((2 * i - 1) * h) for i in range(1, n // 2 + 1)) \
            + 2 * fsum(f(2 * i * h) for i in range(1, n // 2))
        return 0.5 + s * h / 3

    lo, hi = 0.0, 50.0
    for _ in range(70):
        mid = (lo + hi) / 2
        if cdf(mid) < p:
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


def report_fit(zs, vals, label):
    xs = [10 * math.log(z) for z in zs]          # ln X con X = z^10 (convención del doc 101)
    ys = [math.log(v) for v in vals]
    a, b, se = ols(xs, ys)
    t = t_quantile(0.975, len(zs) - 2)
    loo = [ols(xs[:i] + xs[i + 1:], ys[:i] + ys[i + 1:])[1] for i in range(len(xs))]
    lo, hi = b - t * se, b + t * se
    verdict = "ENTERO SOBRE" if lo > THRESH else ("ENTERO DEBAJO" if hi < THRESH else "A CABALLO")
    print(f"  {label}: α = {b:.4f}, IC95% [{lo:.4f}, {hi:.4f}] (t={t:.3f}, n={len(zs)}), LOO [{min(loo):.4f}, {max(loo):.4f}] "
          f"→ {verdict} de {THRESH:.6f}", flush=True)
    return b, lo, hi, verdict


def main():
    try:
        sys.stdout.reconfigure(encoding="utf-8")
    except (AttributeError, ValueError):
        pass
    print(f"Python {sys.version.split()[0]}  |  umbral 1/2 − 2c* = {THRESH:.6f}", flush=True)

    print("\n=== V. Validaciones ===", flush=True)
    DOC101_H = {11: 802, 13: 1628, 17: 5033, 19: 8245, 23: 17999, 29: 47178}
    MY101_H = {59: 990997, 79: 3735822, 97: 9270979, 113: 18622272, 131: 35584592, 151: 68050439}
    pts = {z: point(z) for z in list(DOC101_H) + list(MY101_H)}
    check("V1 |H| por el método de umbrales = doc 101 §2.4 (z=11..29) y = mi enumeración del doc 101 (z=59..151)",
          all(pts[z]['H'] == v for z, v in {**DOC101_H, **MY101_H}.items()))
    ok2 = True
    for z in (11, 13, 17, 19):
        dd = direct_definitions(z)
        p = pts.get(z) or point(z)
        good = dd['H'] == p['H'] and dd['DS'] == p['DS'] and abs(dd['W'] / p['W'] - 1) < 1e-12
        ok2 &= good
        print(f"  z={z}: definición literal |ℋ|={dd['H']} |𝒟_S|={dd['DS']} ΣW={dd['W']:.6f} | umbrales |H|={p['H']} "
              f"|𝒟_S|={p['DS']} ΣW={p['W']:.6f}", flush=True)
    check("V2 |𝒟_S| y Σ W_d por umbrales = definición literal (106.2) enumerando ℋ_C (z=11,13,17,19)", ok2)

    print("\n=== M. Medición: |H|, |𝒟_S| y Σ_{𝒟_S} W_d hasta el mayor z que aguanta la PC ===", flush=True)
    grid = [11, 13, 17, 19, 23, 29, 31, 37, 41, 43, 47, 53, 59, 61, 67, 71, 73, 79, 83, 89, 97, 101, 113, 127,
            139, 151, 167, 181, 199, 227, 251, 277, 307, 353, 401, 449, 503, 557, 607, 659, 709]
    for z in grid:
        if z not in pts:
            t1 = time.time()
            pts[z] = point(z)
            pts[z]['t'] = time.time() - t1
        p = pts[z]
        print(f"  z={z:4d} X=z^10={float(z)**10:.2e} primos pequeños={p['small']:3d} #sprimes={p['nS']:>11d} "
              f"|V+(Lmax)|={p['nV']:>9d} |H|={p['H']:>14d} |𝒟_S|={p['DS']:>14d} ΣW={p['W']:.6e}"
              + (f" ({p['t']:.1f} s)" if 't' in p else ""), flush=True)

    print("\n=== E. Exponentes locales (pendiente entre z consecutivos, respecto de ln X, X = z^10) ===", flush=True)
    for key, lab in (("DS", "|𝒟_S|"), ("W", "ΣW_d"), ("H", "|H|")):
        loc = [(math.log(pts[b][key]) - math.log(pts[a][key])) / (10 * (math.log(b) - math.log(a)))
               for a, b in zip(grid, grid[1:])]
        print(f"  {lab}: " + " ".join(f"{z}:{x:.3f}" for z, x in zip(grid[1:], loc)), flush=True)
        if key == "DS":
            locDS = loc
    tail = [x for z, x in zip(grid[1:], locDS) if z >= 151]
    check("E1 Todos los exponentes locales de |𝒟_S| con z ≥ 151 están sobre 0.420232", min(tail) > THRESH,
          f"mín {min(tail):.4f}, máx {max(tail):.4f}")

    print("\n=== F. Ajustes por ventanas (MCO sobre ln X; puntos deterministas: el IC es formal) ===", flush=True)
    wins = {"z=11..37 (ventana del doc 102)": [z for z in grid if z <= 37],
            "z=11..59": [z for z in grid if z <= 59],
            "z=59..151": [z for z in grid if 59 <= z <= 151],
            "z=151..401": [z for z in grid if 151 <= z <= 401],
            "z=307..709 (ventana superior)": [z for z in grid if z >= 307],
            "z=11..709 (global)": grid}
    res = {}
    for name, win in wins.items():
        for key, lab in (("H", "|H|"), ("DS", "|𝒟_S|"), ("W", "ΣW_d")):
            res[(name, key)] = report_fit(win, [pts[z][key] for z in win], f"{name:32s} {lab:6s}")
    b, lo, hi, v = res[("z=11..37 (ventana del doc 102)", "H")]
    check("F1 Ventana del doc 102: α_H = 0.4211, IC95% [0.4185, 0.4238] reproducido (a caballo del umbral)",
          abs(b - 0.4211) < 6e-5 and abs(lo - 0.4185) < 6e-5 and abs(hi - 0.4238) < 6e-5 and v == "A CABALLO")
    b, lo, hi, v = res[("z=307..709 (ventana superior)", "DS")]
    check("F2 Ventana superior z=307..709: el IC95% de α(|𝒟_S|) queda ENTERO SOBRE 0.420232", v == "ENTERO SOBRE",
          f"α={b:.4f} [{lo:.4f}, {hi:.4f}]")
    b2, lo2, hi2, v2 = res[("z=307..709 (ventana superior)", "W")]
    check("F3 Ventana superior: el IC95% de α(Σ W_d) queda ENTERO SOBRE 0.420232", v2 == "ENTERO SOBRE",
          f"α={b2:.4f} [{lo2:.4f}, {hi2:.4f}]")

    print("\n=== C. Conversión de convención: exponente respecto de N con A = N^c* (X = N^(1−c*)) ===", flush=True)
    for name in ("z=11..37 (ventana del doc 102)", "z=307..709 (ventana superior)"):
        b, lo, hi, _ = res[(name, "DS")]
        print(f"  {name}: α_X(|𝒟_S|) = {b:.4f} → α_N = (1−c*)·α_X = {(1-C_STAR)*b:.4f} "
              f"[{(1-C_STAR)*lo:.4f}, {(1-C_STAR)*hi:.4f}]", flush=True)
    print(f"  Límite asintótico [T]: α_X = 1/2 → α_N = {(1-C_STAR)/2:.6f}; ambos > {THRESH:.6f}", flush=True)

    npass = sum(ok for _, ok in CHECKS)
    print(f"\nRESUMEN: {npass}/{len(CHECKS)} PASS, {len(CHECKS) - npass} FAIL. Tiempo {time.time() - T0:.1f} s", flush=True)
    return 0 if npass == len(CHECKS) else 1


if __name__ == "__main__":
    sys.exit(main())
