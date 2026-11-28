#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Verificación de los docs 102 (E101.1, E101.2), 103 (criba grande restringida, B̂₁, B̂₂) y
104 (E103.1, B̂₃).

Procedencia: Claude Code (Anthropic, modelo claude-opus-5-5), 2026-11-28.
Método: reimplementación PROPIA desde las definiciones escritas en los docs 101–104 (soporte H,
H_N, E_d, M_r, Δ_r por (98.3), B_M, W_d, S_BV, BV1*, BV1_full, Diag, Λ₂, δ_ℋ, B̂₁, B̂₂, B̂₃).
No importa ni ejecuta código de nadie. Los valores publicados aparecen solo como afirmación a
contrastar.
Estado: PROPUESTA (pendiente de contraste por Muse). Mediciones [M]; los ajustes son sobre
ventanas finitas. No es revisión independiente según el contrato (README, regla 4).
Solo biblioteca estándar. Uso:  python -B verificacion_102_104_claude_code.py
Nada de esto demuestra Goldbach, que sigue abierta.
"""
import math
import sys
import time
from itertools import compress

T0 = time.time()
CHECKS = []
fsum = math.fsum


def check(name, ok, detail=""):
    CHECKS.append((name, bool(ok)))
    print(f"[{'PASS' if ok else 'FAIL'}] {name}" + (f" | {detail}" if detail else ""), flush=True)


def sieve(n):
    s = bytearray([1]) * (n + 1)
    s[0] = s[1] = 0
    for p in range(2, math.isqrt(n) + 1):
        if s[p]:
            s[p * p::p] = bytearray(len(range(p * p, n + 1, p)))
    return s


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


def mobius(n):
    f = factor(n)
    return 0 if any(a > 1 for a in f.values()) else (-1) ** len(f)


def phi(n):
    out = n
    for p in factor(n):
        out = out // p * (p - 1)
    return out


def divisors(n):
    ds = [1]
    for p, a in factor(n).items():
        ds = [d * p ** k for d in ds for k in range(a + 1)]
    return ds


def vplus(small_desc, L):
    out = []

    def rec(i0, prod, ell, fs):
        out.append((prod, fs))
        for i in range(i0, len(small_desc)):
            p = small_desc[i]
            if p >= L:
                continue
            np_ = prod * p
            if (ell + 1) % 2 == 1 and p * p * np_ >= L:
                continue
            rec(i + 1, np_, ell + 1, fs + (p,))

    if L > 1:
        rec(0, 1, 0, ())
    return out


def support(z, y, D, A=3):
    """H = divisores de q·d (doc 94), con su factorización. Modelo: N par ⇒ se excluyen 2 y A."""
    s = sieve(max(y, z) + 1)
    pr = [i for i in compress(range(len(s)), s)]
    bad = set(factor(2 * A))
    small = sorted((p for p in pr if p < z and p not in bad), reverse=True)
    sw = [p for p in pr if z <= p < y and p not in bad]
    H = {}
    for q in sw:
        for d, fs in vplus(small, D // q):
            fac = (q,) + fs
            for k in range(1 << len(fac)):
                hs = tuple(sorted(fac[i] for i in range(len(fac)) if k >> i & 1))
                H[math.prod(hs)] = hs
    return H


SUPPORTS = {"R1": (11, 60, 20000), "R3": (3, 349, 200000), "R2": (11, 3000, 200000),
            "S11": (11, round(11 ** (10 / 3)), 11 ** 5), "S13": (13, round(13 ** (10 / 3)), 13 ** 5),
            "S17": (17, round(17 ** (10 / 3)), 17 ** 5)}
ZPAR = {"R1": 11, "R3": 3, "R2": 11, "S11": 11, "S13": 13, "S17": 17}


class Psi:
    def __init__(self, NMAX):
        self.s = sieve(NMAX)
        self.pp = []
        for p in range(2, math.isqrt(NMAX) + 1):
            if self.s[p]:
                x = p * p
                while x <= NMAX:
                    self.pp.append((x, math.log(p)))
                    x *= p

    def prog(self, N, d, a):
        a %= d
        a = a or d
        t = fsum(map(math.log, compress(range(a, N + 1, d), self.s[a:N + 1:d])))
        return t + fsum(w for x, w in self.pp if x <= N and x % d == a % d)

    def lambda2(self, N):
        t = fsum(math.log(p) ** 2 for p in compress(range(N + 1), self.s[:N + 1]))
        return t + fsum(w * w for x, w in self.pp if x <= N)


def delta_r(N, r, hs):
    """Δ_r por la forma reducida (98.3) del doc 98, sin caracteres (válida porque (p, r/p) = 1)."""
    out = []
    if len(hs) == 1:
        p = hs[0]
        x, j = p, 0
        while x <= N:
            j += 1
            x *= p
        out.append(j * math.log(p))
    for p in hs:
        e = r // p
        if e == 1:
            continue
        ds = [(d, mobius(e // d) * phi(d)) for d in divisors(e)]
        lp = math.log(p)
        x = p
        while x <= N:
            out.append(lp * sum(c for d, c in ds if (x - N) % d == 0))
            x *= p
    return fsum(out)


def measure(P, H, N, z):
    """Todas las magnitudes de los docs 102–104 para (soporte, N)."""
    HN = {r: hs for r, hs in H.items() if math.gcd(r, N) == 1}
    psiN = P.prog(N, 1, 0)
    ph = {r: math.prod(p - 1 for p in hs) for r, hs in HN.items()}
    E = {d: (P.prog(N, d, N) - psiN / ph[d]) if d > 1 else 0.0 for d in HN}
    R0 = math.log(N)
    divs = {r: [math.prod(sub) for k in range(1 << len(hs)) for sub in [tuple(hs[i] for i in range(len(hs)) if k >> i & 1)]]
            for r, hs in HN.items()}
    W = {d: 0.0 for d in HN}
    for r in HN:
        for d in divs[r]:
            W[d] += r / ph[r]
    Wmax = max(W.values())
    SBV = fsum(ph[d] * E[d] ** 2 for d in HN)
    BV1full = fsum(ph[d] * E[d] ** 2 * W[d] for d in HN)
    out = dict(HN=len(HN), Wmax=Wmax, SBV=SBV, BV1full=BV1full, psiN=psiN)
    acc = {k: [] for k in ("BM", "BMS", "BMP", "BMC", "BV1", "BV1P", "BV1C", "DiagC", "DiagCS", "B2", "B2all", "B3", "B3all",
                           "intra_num", "intra_den", "BMP_S", "BMC_S")}
    nclass = []
    for r, hs in HN.items():
        Mterms = []
        for d in divs[r]:
            Mterms.append((d, mobius(r // d) * ph[d] * E[d]))
        M = fsum(t for _, t in Mterms)
        Dr = delta_r(N, r, hs)
        S = M + Dr
        absB2 = fsum(abs(t) for _, t in Mterms) + abs(Dr)
        # B̂₃ lectura L2: clase = parte lisa en z (primos < z)
        classes = {}
        for d, t in Mterms:
            h0 = math.prod(p for p in factor(d) if p < z) if d > 1 else 1
            classes[h0] = classes.get(h0, 0.0) + t
        absB3 = fsum(abs(v) for v in classes.values()) + abs(Dr)
        acc["intra_num"].append(fsum(abs(v) for v in classes.values()))
        acc["intra_den"].append(fsum(abs(t) for _, t in Mterms))
        nclass.append(len(classes))
        acc["B2all"].append(absB2 ** 2 / ph[r])
        acc["B3all"].append(absB3 ** 2 / ph[r])
        if r <= R0:
            continue
        isP = len(hs) == 1
        bv = (r / ph[r]) * fsum(ph[d] * E[d] ** 2 for d in divs[r])
        acc["BM"].append(M * M / ph[r])
        acc["BMS"].append(S * S / ph[r])
        acc["BV1"].append(bv)
        acc["B2"].append(absB2 ** 2 / ph[r])
        acc["B3"].append(absB3 ** 2 / ph[r])
        if isP:
            acc["BMP"].append(M * M / ph[r])
            acc["BMP_S"].append(S * S / ph[r])
            acc["BV1P"].append(bv)
        else:
            acc["BMC"].append(M * M / ph[r])
            acc["BMC_S"].append(S * S / ph[r])
            acc["BV1C"].append(bv)
            dg = fsum(ph[d] ** 2 * E[d] ** 2 for d in divs[r]) / ph[r]
            acc["DiagC"].append(dg)
    for k, v in acc.items():
        out[k] = fsum(v)
    out["nclass"] = fsum(nclass) / len(nclass)
    out["Q"] = max(HN)
    tw = [r for r in HN if r + 2 in HN]
    out["r2"] = max(tw) if tw else None
    out["Wmax2"] = max(v for d, v in W.items() if d >= 2)
    # Exhibidos por pieza (doc 103 §3): primos consecutivos de ℋ_P (r > ln N) y gemelo compuesto mayor
    Ps = sorted(r for r, hs in HN.items() if len(hs) == 1 and r > R0)
    out["dinvP"] = max(p * q / (q - p) for p, q in zip(Ps, Ps[1:])) if len(Ps) > 1 else 0.0
    Cs = {r for r, hs in HN.items() if len(hs) >= 2 and r > R0}
    twc = [r for r in Cs if r + 2 in Cs]
    out["dinvC"] = max(r * (r + 2) / 2 for r in twc) if twc else 0.0
    out["fracset"] = [r for r in HN if r > R0]
    return out


def ols(xs, ys):
    n = len(xs)
    mx, my = fsum(xs) / n, fsum(ys) / n
    sxx = fsum((x - mx) ** 2 for x in xs)
    b = fsum((x - mx) * (y - my) for x, y in zip(xs, ys)) / sxx
    a = my - b * mx
    se = math.sqrt(fsum((y - a - b * x) ** 2 for x, y in zip(xs, ys)) / (n - 2) / sxx) if n > 2 else float('nan')
    return a, b, se


T975 = {1: 12.706205, 2: 4.302653, 3: 3.182446, 4: 2.776445, 5: 2.570582}


def slope(xs, ys):
    lx = [math.log(x) for x in xs]
    ly = [math.log(y) for y in ys]
    _, b, se = ols(lx, ly)
    loo = [ols(lx[:i] + lx[i + 1:], ly[:i] + ly[i + 1:])[1] for i in range(len(lx))] if len(lx) > 2 else [b]
    return b, se, T975.get(len(lx) - 2, float('nan')) * se, (min(loo), max(loo))


def exact_delta(H):
    fr = sorted((a / r) for r in H if r > 1 for a in range(1, r) if math.gcd(a, r) == 1)
    fr = [0.0] + fr + [1.0]
    return len(fr) - 2, min(b - a for a, b in zip(fr, fr[1:]))


def main():
    try:
        sys.stdout.reconfigure(encoding="utf-8")
    except (AttributeError, ValueError):
        pass
    print(f"Python {sys.version.split()[0]}", flush=True)
    P = Psi(1048576)
    Hs = {k: support(*v) for k, v in SUPPORTS.items()}
    print("Soportes |H|:", {k: len(v) for k, v in Hs.items()}, flush=True)
    check("S0 |H| = R1 54, R3 68, R2 880, S11 802, S13 1628, S17 5033 (doc 101)",
          [len(Hs[k]) for k in ("R1", "R3", "R2", "S11", "S13", "S17")] == [54, 68, 880, 802, 1628, 5033])
    Ns = {"gen": 999982, "2e": 1048576, "col": 720720}
    res = {}
    for fam, N in Ns.items():
        for k, H in Hs.items():
            if fam != "gen" and k in ("R1", "R3"):
                continue
            t1 = time.time()
            res[(fam, k)] = measure(P, H, N, ZPAR[k])
            m = res[(fam, k)]
            slack = m["Wmax"] * m["SBV"] / m["BM"] if m["BM"] else float('nan')
            print(f"  {fam}/{k}: |H_N|={m['HN']} W_max={m['Wmax']:.2f} S_BV={m['SBV']:.4e} B_M(M)={m['BM']:.4e} "
                  f"B_M(S)={m['BMS']:.4e} slack={slack:.1f} S_C={m['BV1']/m['BM']:.3f} S_cut={m['BV1full']/m['BV1']:.3f} "
                  f"S_max={m['Wmax']*m['SBV']/m['BV1full']:.1f} B^C/B={m['BMC']/m['BM']:.3f} "
                  f"excC={(m['BV1C']-m['BMC'])/(m['BV1']-m['BM']) if m['BV1']>m['BM'] else float('nan'):.3f} "
                  f"Diag^C/B^C={(m['DiagC']/m['BMC']) if m['BMC'] else float('nan'):.3f} ({time.time()-t1:.1f} s)", flush=True)

    print("\n=== Doc 102, E101.1 ===", flush=True)
    D102 = {("gen", "R1"): (54, 68.26, 100.0, 3.95, 25.0, 0.932, 1.705), ("gen", "R3"): (67, 68.29, 69.4, 1.01, 67.8, 0.000, None),
            ("gen", "S11"): (798, 907.85, 988.9, 2.20, 449.4, 0.424, 1.495), ("gen", "R2"): (876, 1006.34, 1102.9, 2.32, 476.0, 0.482, 1.480),
            ("gen", "S13"): (1621, 1889.27, 2104.6, 2.63, 801.6, 0.555, 1.488), ("gen", "S17"): (5017, 5995.20, 7041.8, 3.24, 2175.8, 0.689, 1.533),
            ("2e", "S11"): (802, 912.79, 1070.7, 2.20, 486.6, 0.432, 1.674), ("2e", "R2"): (880, 1011.28, 1193.6, 2.33, 511.3, 0.474, 1.663),
            ("2e", "S13"): (1628, 1898.02, 2234.2, 2.68, 832.2, 0.541, 1.616), ("2e", "S17"): (5033, 6015.88, 7218.7, 3.35, 2153.0, 0.693, 1.576)}
    worst = 0.0
    for key, (h, w, sl, sc, smx, bc, dg) in D102.items():
        m = res[key]
        vals = [(m["HN"], h, 0), (m["Wmax"], w, 0.006), (m["Wmax"] * m["SBV"] / m["BM"], sl, 0.06),
                (m["BV1"] / m["BM"], sc, 0.006), (m["Wmax"] * m["SBV"] / m["BV1full"], smx, 0.06), (m["BMC"] / m["BM"], bc, 0.0006)]
        if dg is not None:
            vals.append((m["DiagC"] / m["BMC"], dg, 0.0006))
        for got, want, tol in vals:
            if abs(got - want) > tol:
                worst = max(worst, abs(got - want))
                print(f"    discrepancia {key}: {got} vs {want}", flush=True)
    check("D102.1 Tablas §2.2 del doc 102 (|H_N|, W_max, slack, S_Cauchy, S_max, B^C/B, Diag^C/B^C) en 10 puntos", worst == 0.0)
    gen = ["R1", "R3", "S11", "R2", "S13", "S17"]
    b, se, ci, loo = slope([res[("gen", k)]["HN"] for k in gen], [res[("gen", k)]["Wmax"] * res[("gen", k)]["SBV"] / res[("gen", k)]["BM"] for k in gen])
    print(f"  θ(genérica) = {b:.3f} ± {se:.3f} (SE), IC95% ±{ci:.3f}, LOO [{loo[0]:.3f}, {loo[1]:.3f}]  (doc 102: 0,988 ± 0,053, ±0,147, [0,927, 1,071])", flush=True)
    ok = abs(b - 0.988) < 6e-4 and abs(se - 0.053) < 6e-4 and abs(ci - 0.147) < 6e-4
    e2 = ["S11", "R2", "S13", "S17"]
    b2, se2, ci2, loo2 = slope([res[("2e", k)]["HN"] for k in e2], [res[("2e", k)]["Wmax"] * res[("2e", k)]["SBV"] / res[("2e", k)]["BM"] for k in e2])
    print(f"  θ(2^e) = {b2:.3f} ± {se2:.3f} (SE), IC95% ±{ci2:.3f}, LOO [{loo2[0]:.3f}, {loo2[1]:.3f}]  (doc 102: 1,036 ± 0,005, ±0,020, [1,032, 1,039])", flush=True)
    ok &= abs(b2 - 1.036) < 6e-4 and abs(se2 - 0.005) < 6e-4 and abs(ci2 - 0.020) < 6e-4
    check("D102.2 θ de la holgura contra |H_N| (familias genérica y 2^e)", ok)
    colok = True
    for k in ("R2", "S11", "S13", "S17"):
        m = res[("col", k)]
        colok &= m["BMC"] == 0.0 and abs(m["BM"] / m["SBV"] - 1) < 1e-12 and abs(m["Wmax"] * m["SBV"] / m["BM"] / m["Wmax"] - 1) < 1e-12
        print(f"  colapso N=720720 {k}: |H_N|={m['HN']} W_max={m['Wmax']:.2f} B_M^C={m['BMC']} B_M/S_BV={m['BM']/m['SBV']:.12f}", flush=True)
    check("D102.3 Familia colapso: B_M^C = 0, B_M = S_BV y slack = W_max exactos en R2/S11/S13/S17", colok)

    print("\n=== Doc 103: B̂₁ (criba grande restringida) y B̂₂ ===", flush=True)
    L2 = {fam: P.lambda2(N) for fam, N in Ns.items() if fam != "col"}
    print(f"  Λ₂(999982) = {L2['gen']:.6e} (doc: 1,280280e7); Λ₂(2^20) = {L2['2e']:.6e} (doc: 1,347870e7)", flush=True)
    check("D103.1 Λ₂(N) reproducido", abs(L2['gen'] / 1.280280e7 - 1) < 1e-6 and abs(L2['2e'] / 1.347870e7 - 1) < 1e-6)
    D103 = {("gen", "S11"): (16345, 2933, 2.02e4, 1.02e6), ("gen", "S13"): (37415, 10703, 9.76e4, 2.35e6), ("gen", "S17"): (144445, 81653, 1.53e6, 9.58e6),
            ("2e", "S11"): (16345, 2933, 1.94e4, 9.74e5), ("2e", "S13"): (37415, 10703, 9.44e4, 2.27e6), ("2e", "S17"): (144445, 81653, 1.51e6, 9.45e6)}
    six = list(D103)
    ok = True
    lo_ratio = {}
    for key, (Q, r2, lo, hi) in D103.items():
        m = res[key]
        N = Ns[key[0]]
        blo = (N + m["r2"] * (m["r2"] + 2) / 2) * L2[key[0]] / m["BMS"]
        bhi = (N + m["Q"] ** 2) * L2[key[0]] / m["BMS"]
        lo_ratio[key] = blo
        ok &= m["Q"] == Q and m["r2"] == r2 and abs(blo / lo - 1) < 6e-3 and abs(bhi / hi - 1) < 6e-3
        print(f"  {key}: Q={m['Q']} r2={m['r2']} B̂_lo/B_M={blo:.3e} B̂_hi/B_M={bhi:.3e} (doc: {lo:.2e}, {hi:.2e})", flush=True)
    check("D103.2 Q_ℋ, r₂ y cocientes B̂_lo/B_M, B̂_hi/B_M del doc 103 §2", ok)
    b, se, ci, loo = slope([res[k]["HN"] for k in six], [lo_ratio[k] for k in six])
    print(f"  pendiente B̂_lo/B_M: {b:.4f} ± {se:.4f}, IC95% ±{ci:.4f}, LOO [{loo[0]:.4f}, {loo[1]:.4f}] (doc: 2,3701 ± 0,0322, ±0,0895, [2,3518, 2,3950])", flush=True)
    check("D103.3 C2 de B̂₁: pendiente 2,3701", abs(b - 2.3701) < 6e-4 and abs(se - 0.0322) < 6e-4)
    for k in ("R1", "R3"):
        nf, dl = exact_delta(res[("gen", k)] and {r: 0 for r in Hs[k] if math.gcd(r, 999982) == 1})
        m = res[("gen", k)]
        r = (999982 + 1 / dl) * L2["gen"] / m["BMS"]
        print(f"  δ_ℋ exacto gen/{k}: #fracciones={nf} δ={dl:.4e} δ⁻¹={1/dl:.4e} cociente B̂₁/B_M={r:.3e}", flush=True)
    for col in ("B2", "B2all"):
        rr = {k: res[k][col] / res[k]["BMS"] for k in six}
        b, se, ci, loo = slope([res[k]["HN"] for k in six], [rr[k] for k in six])
        print(f"  B̂₂ ({'r>ln N' if col=='B2' else 'todo r∈ℋ'}): cocientes " + ", ".join(f"{k[0]}/{k[1]}={rr[k]:.3f}" for k in six)
              + f" | R1={res[('gen','R1')][col]/res[('gen','R1')]['BMS']:.3f} R3={res[('gen','R3')][col]/res[('gen','R3')]['BMS']:.3f}"
              + f" | pendiente {b:.4f} ± {se:.4f}, LOO [{loo[0]:.4f}, {loo[1]:.4f}]", flush=True)
        if col == "B2":
            B2ok = (abs(rr[("gen", "S11")] - 1.709) < 6e-3 and abs(rr[("2e", "S17")] - 2.318) < 6e-3 and abs(b - 0.1358) < 6e-4)
    check("D103.4 V103.1: cocientes B̂₂/B_M (1,709…2,318) y pendiente 0,1358 (C1 pasa, C2 falla)", B2ok)

    print("\n=== Doc 104: B̂₃ (lectura L2: clase = parte lisa en z) ===", flush=True)
    D104 = {("gen", "S11"): 1.679, ("gen", "S13"): 1.876, ("gen", "S17"): 2.181, ("2e", "S11"): 1.798, ("2e", "S13"): 2.000, ("2e", "S17"): 2.262}
    for col in ("B3", "B3all"):
        rr = {k: res[k][col] / res[k]["BMS"] for k in six}
        b, se, ci, loo = slope([res[k]["HN"] for k in six], [rr[k] for k in six])
        q23 = {k: res[k]["B2" if col == "B3" else "B2all"] / res[k][col] for k in six}
        b23, se23, _, loo23 = slope([res[k]["HN"] for k in six], [q23[k] for k in six])
        print(f"  B̂₃ ({'r>ln N' if col=='B3' else 'todo r∈ℋ'}): " + ", ".join(f"{k[0]}/{k[1]}={rr[k]:.3f}" for k in six)
              + f" | R1={res[('gen','R1')][col]/res[('gen','R1')]['BMS']:.3f} | pendiente {b:.4f} ± {se:.4f} (IC ±{ci:.4f}), LOO [{loo[0]:.4f}, {loo[1]:.4f}]"
              + f" | pendiente B̂₂/B̂₃ {b23:.4f} ± {se23:.4f}", flush=True)
        if col == "B3":
            B3ok = all(abs(rr[k] - v) < 6e-3 for k, v in D104.items()) and abs(b - 0.1326) < 6e-4 and abs(se - 0.0195) < 6e-4
            dec_ok = abs(b + b23 - slope([res[k]["HN"] for k in six], [res[k]["B2"] / res[k]["BMS"] for k in six])[0]) < 1e-9
    check("D104.1 Cocientes B̂₃/B_M del doc 104 §2 y pendiente 0,1326 ± 0,0195 (C1 pasa, C2 falla)", B3ok)
    check("D104.2 Descomposición exacta de pendientes: pend(B̂₂/B_M) = pend(B̂₂/B̂₃) + pend(B̂₃/B_M)", dec_ok)
    rho = {k: res[k]["intra_num"] / res[k]["intra_den"] for k in six}
    print("  ρ_intra: " + ", ".join(f"{k[0]}/{k[1]}={v:.3f}" for k, v in rho.items()) + "  (doc 104: 0,987–0,992)", flush=True)
    print("  clases por conductor: " + ", ".join(f"{k[0]}/{k[1]}={res[k]['nclass']:.2f}" for k in six), flush=True)
    check("D104.3 ρ_intra ∈ [0,987, 0,992] en los seis puntos", all(0.9865 <= v <= 0.9925 for v in rho.values()))
    l1 = all(all(math.gcd(r, 3) == 1 for r in Hs[k]) for k in Hs)
    check("D104.4 Lectura L1 (nueva(d) = d/(d,3)): todo r ∈ ℋ es coprimo con 3 ⇒ B̂₃ = B̂₂ idénticamente", l1)

    print("\n=== Extras: [P-rev] resolubles por cómputo ===", flush=True)
    # doc 103 §3: desglose primo/compuesto del sándwich
    D103_3 = {("gen", "S11"): (3.25e4, 4.76e4, 1.495), ("gen", "S13"): (5.27e4, 1.76e5, 1.488), ("gen", "S17"): (1.19e5, 2.22e6, 1.533),
              ("2e", "S11"): (3.18e4, 4.50e4, 1.674), ("2e", "S13"): (4.95e4, 1.74e5, 1.615), ("2e", "S17"): (1.19e5, 2.18e6, 1.576)}
    ok = True
    for key, (lp, lc, dgs) in D103_3.items():
        m = res[key]
        N = Ns[key[0]]
        rp = (N + m["dinvP"]) * L2[key[0]] / m["BMP_S"]
        rc = (N + m["dinvC"]) * L2[key[0]] / m["BMC_S"]
        dS = m["DiagC"] / m["BMC_S"]
        ok &= abs(rp / lp - 1) < 6e-3 and abs(rc / lc - 1) < 6e-3 and abs(dS - dgs) < 6e-4
        print(f"  {key}: B̂_lo^P/B_M^P={rp:.3e} B̂_lo^C/B_M^C={rc:.3e} Diag^C/B_M^C(S)={dS:.3f}  (doc 103 §3: {lp:.2e}, {lc:.2e}, {dgs})", flush=True)
    check("D103.5 Desglose primo/compuesto del doc 103 §3 (exhibidos por pieza) y Diag^C/B_M^C con S_r", ok)
    # [P-rev] 102.4: max_{d>=2} W_d
    for key in [("gen", k) for k in gen] + [("2e", k) for k in e2]:
        m = res[key]
        print(f"  max_(d>=2) W_d {key}: {m['Wmax2']:.2f} (W_1 = {m['Wmax']:.2f}, cociente {m['Wmax2']/m['Wmax']:.4f})", flush=True)
    # [P-rev] 103.2: δ_ℋ exacto en S11 y S13 (familia genérica)
    for k in ("S11", "S13"):
        t1 = time.time()
        m = res[("gen", k)]
        nf, dl = exact_delta({r: 0 for r in m["fracset"]})
        r = (999982 + 1 / dl) * L2["gen"] / m["BMS"]
        lo = (999982 + m["r2"] * (m["r2"] + 2) / 2) * L2["gen"] / m["BMS"]
        hi = (999982 + m["Q"] ** 2) * L2["gen"] / m["BMS"]
        print(f"  δ_ℋ exacto gen/{k}: #fracciones={nf} δ⁻¹={1/dl:.4e} ⇒ B̂₁/B_M = {r:.3e} (sándwich [{lo:.3e}, {hi:.3e}]; "
              f"B̂₁/B̂_lo = {r/lo:.3f}) ({time.time()-t1:.1f} s)", flush=True)
    # [P-rev] 104.3: B̂₃ en soportes mayores S19, S23
    big = {"S19": (19, round(19 ** (10 / 3)), 19 ** 5), "S23": (23, round(23 ** (10 / 3)), 23 ** 5)}
    allpts = {key: res[key] for key in six}
    for k, v in big.items():
        Hb = support(*v)
        for fam in ("gen", "2e"):
            t1 = time.time()
            allpts[(fam, k)] = measure(P, Hb, Ns[fam], v[0])
            m = allpts[(fam, k)]
            print(f"  {fam}/{k}: |H_N|={m['HN']} B̂₂/B_M={m['B2']/m['BMS']:.3f} B̂₃/B_M={m['B3']/m['BMS']:.3f} "
                  f"B^C/B={m['BMC']/m['BM']:.3f} ({time.time()-t1:.1f} s)", flush=True)
    for fam in ("gen", "2e"):
        ks = [(fam, k) for k in ("S11", "S13", "S17", "S19", "S23")]
        b, se, ci, loo = slope([allpts[k]["HN"] for k in ks], [allpts[k]["B3"] / allpts[k]["BMS"] for k in ks])
        print(f"  pendiente B̂₃/B_M familia {fam}, S11..S23 (5 puntos): {b:.4f} ± {se:.4f}, IC95% ±{ci:.4f}, LOO [{loo[0]:.4f}, {loo[1]:.4f}]", flush=True)
    ks = [k for k in allpts]
    b, se, ci, loo = slope([allpts[k]["HN"] for k in ks], [allpts[k]["B3"] / allpts[k]["BMS"] for k in ks])
    print(f"  pendiente B̂₃/B_M, diez puntos (S11..S23 × 2 familias): {b:.4f} ± {se:.4f}, LOO [{loo[0]:.4f}, {loo[1]:.4f}]", flush=True)
    check("D104.5 Con soportes mayores (S19, S23) la pendiente de B̂₃/B_M sigue sobre 0,1 (C2 sigue fallando)", b - 0.0 > 0.1,
          f"{b:.4f}")

    npass = sum(ok for _, ok in CHECKS)
    print(f"\nRESUMEN: {npass}/{len(CHECKS)} PASS, {len(CHECKS) - npass} FAIL. Tiempo {time.time() - T0:.1f} s", flush=True)
    return 0 if npass == len(CHECKS) else 1


if __name__ == "__main__":
    sys.exit(main())
