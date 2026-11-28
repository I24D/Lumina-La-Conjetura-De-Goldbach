#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Contraste numérico del doc 107 (L105(a): B_Δ ≤ C_Δ·T(N)) y de sus calibraciones.

Procedencia: Claude Code (Anthropic, modelo claude-opus-5-5), 2026-11-28.
Método: reimplementación PROPIA desde las definiciones de los docs 97, 98, 100, 101, 106 y 107:
soporte ℋ (doc 101, divisores de q·d), H_N, Δ_r por (98.3) sin caracteres y, para el
contraejemplo, también como S_r − M_r con S_r por la ruta sin caracteres (independiente de
(98.3)); B_Δ = Σ_{r∈ℋ} Δ_r²/φ(r); T(N) = N^(2−2c)/(0,0225 ln⁴N); M_ω, M_φ; Ω_d^C por (106.2)
enumerando ℋ_C; niveles (L|S, k = ω(d)); ℳ_Ω real con ψ. No importa ni ejecuta código de nadie.
No revisa la prueba escrita (encargo: recomputar sus calibraciones), aunque comprueba su álgebra.
Estado: PROPUESTA (pendiente de contraste por Muse). [M] salvo donde se indica [T].
No es revisión independiente según el contrato (README, regla 4).
Solo biblioteca estándar. Uso:  python -B verificacion_107_claude_code.py
Nada de esto demuestra Goldbach, que sigue abierta.
"""
import math
import sys
import time
from itertools import compress

T0 = time.time()
CHECKS = []
fsum = math.fsum
C = 0.039884
G = 0.0225


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


def support(z, y, D):
    s = sieve(max(y, z) + 1)
    pr = [i for i in compress(range(len(s)), s)]
    small = sorted((p for p in pr if p < z and p not in (2, 3)), reverse=True)
    sw = [p for p in pr if z <= p < y and p not in (2, 3)]
    H = {}
    for q in sw:
        for d, fs in vplus(small, D // q):
            fac = (q,) + fs
            for k in range(1 << len(fac)):
                hs = tuple(sorted(fac[i] for i in range(len(fac)) if k >> i & 1))
                H[math.prod(hs)] = hs
    return H, sw, small


def fam(z):
    return (z, round(z ** (10 / 3)), z ** 5)


SUP = {"R1": (11, 60, 20000), "R3": (3, 349, 200000), "S11": fam(11), "S13": fam(13), "S17": fam(17),
       "S19": fam(19), "S23": fam(23), "S29": fam(29), "R2": (11, 3000, 200000)}
ZP = {"R1": 11, "R3": 3, "S11": 11, "S13": 13, "S17": 17, "S19": 19, "S23": 23, "S29": 29, "R2": 11}


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
        self.lam = None

    def prog(self, N, d, a):
        a %= d
        a = a or d
        t = fsum(map(math.log, compress(range(a, N + 1, d), self.s[a:N + 1:d])))
        return t + fsum(w for x, w in self.pp if x <= N and x % d == a % d)

    def residues(self, N, r):
        cnt = [0.0] * r
        for p in compress(range(N + 1), self.s[:N + 1]):
            cnt[p % r] += math.log(p)
        for x, w in self.pp:
            if x <= N:
                cnt[x % r] += w
        return cnt


def delta_983(N, r, hs):
    """Δ_r por (98.3), forma sin caracteres (válida porque (p, r/p) = 1)."""
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
        x = p
        while x <= N:
            out.append(math.log(p) * sum(c for d, c in ds if (x - N) % d == 0))
            x *= p
    return fsum(out)


def T_of(N):
    return N ** (2 - 2 * C) / (G * math.log(N) ** 4)


def omega_C(H):
    """Ω_d^C = Σ_{r∈ℋ_C, d|r} (τ(r)−1)/φ(r), por enumeración literal (106.2)."""
    Om = {}
    for r, hs in H.items():
        if len(hs) < 2:
            continue
        val = ((1 << len(hs)) - 1) / math.prod(p - 1 for p in hs)
        for k in range(1, 1 << len(hs)):
            sub = tuple(hs[i] for i in range(len(hs)) if k >> i & 1)
            Om[math.prod(sub)] = Om.get(math.prod(sub), 0.0) + val
    return Om


def main():
    try:
        sys.stdout.reconfigure(encoding="utf-8")
    except (AttributeError, ValueError):
        pass
    print(f"Python {sys.version.split()[0]}", flush=True)
    P = Psi(1048576)
    supp = {k: support(*v) for k, v in SUP.items()}
    check("K0 Maquinaria: |H| R2=880, S11=802, S13=1628, S17=5033, S23=17999 (doc 101)",
          [len(supp[k][0]) for k in ("R2", "S11", "S13", "S17", "S23")] == [880, 802, 1628, 5033, 17999])
    Ns = {"gen": 999982, "2e": 1048576}

    print("\n=== §1.2 Lema (106.6) |Δ_r| ≤ φ(r)ω(r)ln N, y la errata del doc 97 ===", flush=True)
    # Convención: el doc 107 escribe Σ_{r∈ℋ}, pero sus calibraciones se reproducen con
    # ℋ_N ∩ {r > ln N} (el conjunto de B_M). Se calculan las dos y se informa de ambas.
    BD, BDall = {}, {}
    worst106 = (0.0, None)
    worst106all = (0.0, None)
    worst97 = (0.0, None)
    Mw = {}
    Mphi = {}
    Mw_hi = {}
    for famn, N in Ns.items():
        lnN = math.log(N)
        for k in ("R1", "R3", "S11", "S13", "S17", "S23"):
            H = {r: hs for r, hs in supp[k][0].items() if math.gcd(r, N) == 1 and r > 1}
            terms, terms_all = [], []
            for r, hs in H.items():
                Dr = delta_983(N, r, hs)
                ph = math.prod(p - 1 for p in hs)
                terms_all.append(Dr * Dr / ph)
                q1 = abs(Dr) / (ph * len(hs) * lnN)
                q2 = abs(Dr) / (len(hs) * lnN)
                if q1 > worst106all[0]:
                    worst106all = (q1, (famn, k, r, Dr))
                if r > lnN:
                    terms.append(Dr * Dr / ph)
                    if q1 > worst106[0]:
                        worst106 = (q1, (famn, k, r, Dr))
                if q2 > worst97[0]:
                    worst97 = (q2, (famn, k, r, Dr))
            BD[(famn, k)] = fsum(terms)
            BDall[(famn, k)] = fsum(terms_all)
            Mw[(famn, k)] = fsum(math.prod(p - 1 for p in hs) * len(hs) ** 2 for r, hs in H.items())
            Mw_hi[(famn, k)] = fsum(math.prod(p - 1 for p in hs) * len(hs) ** 2 for r, hs in H.items() if r > lnN)
            Mphi[(famn, k)] = fsum(math.prod(p - 1 for p in hs) for r, hs in H.items())
    print(f"  r > ln N: máx |Δ_r|/(φ(r)ω(r)ln N) = {worst106[0]:.4f} en {worst106[1][:3]} (Δ = {worst106[1][3]:.4f}); doc 107: 0,0513 en r = 17", flush=True)
    print(f"  todo r ∈ ℋ_N: máx = {worst106all[0]:.4f} en {worst106all[1][:3]} (Δ = {worst106all[1][3]:.4f}) — la cota (106.6) también vale", flush=True)
    check("K1a Cota (106.6) se cumple con holgura; máximo 0,0513 en r = 17 (Δ_17 = 11,3329) sobre r > ln N, y ≤ 1 sobre todo ℋ_N",
          abs(worst106[0] - 0.0513) < 6e-5 and worst106[1][2] == 17 and abs(worst106[1][3] - 11.3329) < 6e-5 and worst106all[0] <= 1)
    print(f"  máx |Δ_r|/(ω(r)ln N) (la cota impresa del doc 97, sin φ) = {worst97[0]:.1f} en {worst97[1][:3]}", flush=True)
    # Contraejemplo r = 458815, N = 999982, por dos rutas
    r0, N0 = 458815, 999982
    hs0 = tuple(sorted(factor(r0)))
    d1 = delta_983(N0, r0, hs0)
    cnt = P.residues(N0, r0)
    ds = [(d, mobius(r0 // d) * phi(d)) for d in divisors(r0)]
    S = fsum(cnt[a] * sum(c for d, c in ds if (a - N0) % d == 0) for a in range(r0) if cnt[a] and math.gcd(a, r0) == 1)
    psiN = P.prog(N0, 1, 0)
    M = fsum(mobius(r0 // d) * phi(d) * (fsum(cnt[a] for a in range(N0 % d, r0, d)) - psiN / phi(d)) for d in divisors(r0))
    d2 = S - M
    fac = abs(d1) / (len(hs0) * math.log(N0))
    share = d1 * d1 / phi(r0) / BD[("gen", "S23")]
    print(f"  r = {r0} = {'·'.join(map(str, hs0))}, φ = {phi(r0)}, ω = {len(hs0)}: Δ (98.3) = {d1:.6f}; Δ = S_r − M_r (ruta sin caracteres) = {d2:.6f}", flush=True)
    print(f"  factor sobre ω(r)·ln N = {fac:.1f}; Δ²/φ = {d1*d1/phi(r0):.4e} = {100*share:.1f}% de B_Δ(gen/S23)", flush=True)
    check("K1b Contraejemplo a la cota impresa del doc 97: r = 458815, |Δ_r| = 76522,068094 por dos rutas, factor 1846,3, 21% de B_Δ",
          r0 in supp["S23"][0] and abs(d1 - 76522.068094) < 1e-4 and abs(d2 - d1) < 1e-4 and abs(fac - 1846.3) < 0.06
          and abs(share - 0.21) < 0.006)
    hol = min(math.log(Ns[k[0]]) ** 2 * Mw_hi[k] / BD[k] for k in BD)
    hol_all = min(math.log(Ns[k[0]]) ** 2 * Mw[k] / BDall[k] for k in BD)
    print(f"  holgura mínima de B_Δ ≤ (ln N)²M_ω: {hol:.3g} (r > ln N en ambos lados); {hol_all:.3g} (todo ℋ_N)", flush=True)
    check("K1c B_Δ ≤ (ln N)²·M_ω en los 12 puntos, con holgura ≥ 2·10⁴ (convención r > ln N)", hol >= 2e4)

    print("\n=== §1.4 Momentos M_ω, M_φ (familia gen) por ramas y por enumeración ===", flush=True)
    DOCMW = {"S11": (7.14e6, 1.64e6), "S13": (3.95e7, 7.43e6), "S17": (5.64e8, 8.95e7), "S19": (1.73e9, 2.52e8), "S23": (1.07e10, 1.37e9)}
    ok = True
    for k, (mw, mp) in DOCMW.items():
        z, y, D = SUP[k]
        H, sw, small = supp[k]
        N = Ns["gen"]
        # por ramas: ℋ = 𝒱₀ ⊔ ⨆ s·𝒱_s, filtrando (r, N) = 1
        V0 = vplus(small, D // sw[0])
        bw = fsum(math.prod(p - 1 for p in fs) * len(fs) ** 2 for d, fs in V0 if d > 1 and math.gcd(d, N) == 1)
        bp = fsum(math.prod(p - 1 for p in fs) for d, fs in V0 if d > 1 and math.gcd(d, N) == 1)
        for s in sw:
            if N % s == 0:
                continue
            for d, fs in vplus(small, D // s):
                if math.gcd(d, N) != 1:
                    continue
                ph = (s - 1) * math.prod(p - 1 for p in fs)
                bw += ph * (len(fs) + 1) ** 2
                bp += ph
        if k not in ("S19",):
            ew, ep = Mw.get(("gen", k), None), Mphi.get(("gen", k), None)
        else:
            Hn = {r: hs for r, hs in H.items() if math.gcd(r, N) == 1 and r > 1}
            ew = fsum(math.prod(p - 1 for p in hs) * len(hs) ** 2 for hs in Hn.values())
            ep = fsum(math.prod(p - 1 for p in hs) for hs in Hn.values())
        ok &= abs(bw / ew - 1) < 1e-12 and abs(bp / ep - 1) < 1e-12 and abs(bw / mw - 1) < 6e-3 and abs(bp / mp - 1) < 6e-3
        print(f"  {k}: M_ω ramas={bw:.4e} enum={ew:.4e} (doc {mw:.3g}); M_φ ramas={bp:.4e} enum={ep:.4e} (doc {mp:.3g})", flush=True)
    check("K2 M_ω y M_φ por ramas = enumeración (error relativo < 10⁻¹²) y = valores del doc 107 §1.4", ok)

    print("\n=== §1.5 B_Δ/T(N) en los 12 puntos y la cadena (107.1)–(107.2) ===", flush=True)
    for famn, N in Ns.items():
        print(f"  T({N}) = {T_of(N):.4e}", flush=True)
    rat = {k: BD[k] / T_of(Ns[k[0]]) for k in BD}
    for k, v in rat.items():
        print(f"  {k[0]}/{k[1]}: B_Δ(r>ln N) = {BD[k]:.4e}  B_Δ/T = {v:.3e}  | todo ℋ_N: {BDall[k]:.4e} ({BDall[k]/T_of(Ns[k[0]]):.3e})", flush=True)
    kmax = max(rat, key=rat.get)
    check("K3 máx B_Δ/T = 2,159·10⁻⁴ en gen/S23; R3 ≈ 2,3·10⁻⁷ y R1 ≈ 4,4·10⁻⁷; C_Δ = 10⁻³ cubre con factor ≥ 4,6",
          kmax == ("gen", "S23") and abs(rat[kmax] - 2.159e-4) < 6e-7 and abs(rat[("gen", "R3")] / 2.3e-7 - 1) < 0.03
          and abs(rat[("gen", "R1")] / 4.4e-7 - 1) < 0.03 and 1e-3 / rat[kmax] >= 4.6, f"T(999982)={T_of(999982):.4e}")
    k1 = 1 / (4 * math.log(2) ** 2)
    k2 = G * k1
    print(f"  (107.1): 1/(4 ln²2) = {k1:.6f}; (107.2): 0,0225/(4 ln²2) = {k2:.7f}; exponente 1−2c = {1-2*C:.6f}", flush=True)
    f = lambda lN: k2 * lN ** 8 * math.exp(-(1 - 2 * C) * lN)
    # decreciente para ln N > 8/(1−2c)
    print(f"  cociente decreciente para N > e^(8/(1−2c)) = {math.exp(8/(1-2*C)):.0f}", flush=True)
    n0 = {}
    for Cd in (1.0, 1e-2, 1e-3):
        lo, hi = 8 / (1 - 2 * C), 200.0
        for _ in range(200):
            mid = (lo + hi) / 2
            if f(mid) > Cd:
                lo = mid
            else:
                hi = mid
        n0[Cd] = math.exp(hi)
        print(f"  C_Δ = {Cd:g}: N₀ = {n0[Cd]:.3e}", flush=True)
    check("K6 Constantes de (107.1)–(107.2) y N₀ = 3,82·10⁹ / 7,40·10¹² / 2,36·10¹⁴ por bisección propia; "
          "«decreciente para N > 5970» es verdadero (umbral exacto 5964,4)",
          abs(k1 - 0.520342) < 6e-7 and abs(k2 - 0.0117077) < 6e-8 and math.exp(8 / (1 - 2 * C)) <= 5970
          and abs(n0[1.0] / 3.82e9 - 1) < 3e-3 and abs(n0[1e-2] / 7.40e12 - 1) < 3e-3 and abs(n0[1e-3] / 2.36e14 - 1) < 3e-3)
    # la cota (107.1) y la medida en cada punto
    ok = all(BD[k] <= k1 * Ns[k[0]] * math.log(Ns[k[0]]) ** 4 for k in BD)
    check("K6b B_Δ ≤ 0,520342·N(ln N)⁴ en los 12 puntos (la cota [T] de (107.1))", ok)

    print("\n=== §2.1 Masa de peso por niveles (sin ψ) y d > D_A^(2/3) ===", flush=True)
    DOCL = {11: (0.0000, 0.019, 0.381, 0.600, None, None, 0.798), 13: (0.0001, 0.016, 0.329, 0.636, 0.019, None, 0.857),
            17: (0.0003, 0.014, 0.282, 0.652, 0.042, 0.011, 0.916), 19: (0.0003, 0.011, 0.239, 0.644, 0.070, 0.037, 0.931),
            23: (0.0001, 0.007, 0.165, 0.660, 0.100, 0.068, 0.954), 29: (0.0001, 0.004, 0.118, 0.643, 0.124, 0.111, 0.973)}
    best = {}
    for label, filt in (("ℋ (sin filtro N)", None), ("ℋ_N, N = 999982", 999982)):
        dev = 0.0
        for z in DOCL:
            k = f"S{z}"
            H = supp[k][0] if filt is None else {r: hs for r, hs in supp[k][0].items() if math.gcd(r, filt) == 1}
            Om = omega_C(H)
            DA = z ** 5
            tot = fsum(v * phi(d) ** 2 for d, v in Om.items())
            lev = {}
            big = 0.0
            small_ = 0.0
            for d, v in Om.items():
                hs = H[d]
                key = ("S" if hs[-1] >= z else "L", len(hs))
                w = v * phi(d) ** 2
                lev[key] = lev.get(key, 0.0) + w
                if d > DA ** (2 / 3):
                    big += w
                if d <= DA ** (1 / 3):
                    small_ += w
            Lsum = fsum(v for kk, v in lev.items() if kk[0] == "L") / tot
            row = (Lsum,) + tuple((lev.get(("S", j), 0.0) / tot) for j in range(1, 6)) + (big / tot,)
            want = DOCL[z]
            for a, b in zip(row, want):
                if b is not None:
                    dev = max(dev, abs(a - b))
            if filt == 999982 or (filt is None and z in (11, 29)):
                print(f"  [{label}] z={z}: L={row[0]:.4f} S1={row[1]:.3f} S2={row[2]:.3f} S3={row[3]:.3f} S4={row[4]:.3f} S5={row[5]:.3f} "
                      f"d>D_A^(2/3)={row[6]:.3f} d≤D_A^(1/3)={small_/tot:.3f} S2–S5={sum(row[2:6]):.3f} |𝒟_C|={len(Om)}", flush=True)
        best[label] = dev
        print(f"  [{label}] desviación máxima frente a la tabla del doc 107 §2.1: {dev:.4f}", flush=True)
    check("K4 Tabla de masas por nivel del doc 107 §2.1 reproducida (a 3 decimales) en z = 11..29", min(best.values()) < 6e-4,
          f"mejor convención: {min(best, key=best.get)}")

    print("\n=== §2.1 punto 4: calibración con ψ (masa real de ℳ_Ω por nivel, N = 999982) ===", flush=True)
    DOCPSI = {11: 0.906, 13: 0.908, 17: 0.877, 19: 0.846, 23: 0.809}
    ok = True
    N = 999982
    psiN = P.prog(N, 1, 0)
    for z, want in DOCPSI.items():
        H = {r: hs for r, hs in supp[f"S{z}"][0].items() if math.gcd(r, N) == 1}
        Om = omega_C(H)
        lev = {}
        for d, v in Om.items():
            E = P.prog(N, d, N) - psiN / phi(d)
            hs = H[d]
            key = ("S" if hs[-1] >= z else "L", len(hs))
            lev[key] = lev.get(key, 0.0) + v * phi(d) ** 2 * E * E
        tot = fsum(lev.values())
        s23 = (lev.get(("S", 2), 0) + lev.get(("S", 3), 0)) / tot
        s45 = (lev.get(("S", 4), 0) + lev.get(("S", 5), 0)) / tot
        Lm = fsum(v for kk, v in lev.items() if kk[0] == "L") / tot
        ok &= abs(s23 - want) < 6e-4
        print(f"  z={z}: S2+S3 = {s23:.3f} (doc {want}); S3 = {lev.get(('S',3),0)/tot:.3f}; S4+S5 = {s45:.3f}; L = {Lm:.4f}", flush=True)
    check("K4b Calibración con ψ: S2+S3 = 0,906 / 0,908 / 0,877 / 0,846 / 0,809 (z = 11..23)", ok)

    print("\n=== §2.2 / K5: W* y cota de criba grande por nivel (marco N_f = z^10) ===", flush=True)
    DOCW = {11: 9.52, 13: 15.00, 17: 31.00, 19: 31.00, 23: 31.00}
    DOCSUM = {11: 28.0, 29: 189.6}
    DOCLS = {11: 3.39e7, 13: 1.25e8, 17: 5.47e8, 19: 7.76e8, 23: 1.36e9, 29: 2.72e9}
    okW = True
    okLS = True
    for label, filt in (("ℋ", None), ("ℋ_N", 999982)):
        for z in (11, 13, 17, 19, 23, 29):
            H = supp[f"S{z}"][0] if filt is None else {r: hs for r, hs in supp[f"S{z}"][0].items() if math.gcd(r, filt) == 1}
            Om = omega_C(H)
            Wl = {}
            Ql = {}
            for d, v in Om.items():
                hs = H[d]
                key = ("S" if hs[-1] >= z else "L", len(hs))
                Wl[key] = max(Wl.get(key, 0.0), v * phi(d))
                Ql[key] = max(Ql.get(key, 0), d)
            Wst = max(Wl.values())
            Nf = float(z) ** 10
            Lam2 = Nf * math.log(Nf)
            ls = fsum(Wl[kk] * (Nf + Ql[kk] ** 2) * Lam2 for kk in Wl) / T_of(Nf)
            print(f"  [{label}] z={z}: W* = {Wst:.2f}; Σ_ℓ W*_ℓ = {fsum(Wl.values()):.1f}; exceso solo-LS = {ls:.3g} (doc: {DOCLS[z]:.3g})", flush=True)
            if label == "ℋ_N":
                if z in DOCW:
                    okW &= abs(Wst - DOCW[z]) < 6e-3
                if z in DOCSUM:
                    okW &= abs(fsum(Wl.values()) - DOCSUM[z]) < 0.06
                okLS &= abs(ls / DOCLS[z] - 1) < 0.02 and ls > 1
    check("K5a W* = 9,52 / 15,00 / 31,00 (z = 11/13/17–23) y Σ_ℓ W*_ℓ = 28,0 (z = 11) … 189,6 (z = 29)", okW)
    check("K5b Exceso de la cota de criba grande por nivel frente a T(N_f) (con Λ₂ ≈ N ln N): reproduce la columna solo-LS y es > 1", okLS)

    npass = sum(ok for _, ok in CHECKS)
    print(f"\nRESUMEN: {npass}/{len(CHECKS)} PASS, {len(CHECKS) - npass} FAIL. Tiempo {time.time() - T0:.1f} s", flush=True)
    return 0 if npass == len(CHECKS) else 1


if __name__ == "__main__":
    sys.exit(main())
