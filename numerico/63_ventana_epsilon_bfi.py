#!/usr/bin/env python3
"""Doc 63: ventana epsilon en (1/4,2/7] — ¿BFI supera a BV para K=2?
Marco: docs 58 (BFI/Maynard Def.1 Teor.A, Prop 58.1), 61 (K=2, Prop 61.4),
62 (Prop 62.1, Q*~N^{2eps}, s=1/eps, seleccion por criba).
Parte 1: K=2 exacto (Fraction) en z de la ventana a N=1e12 (z=31,52).
Parte 2: parrilla (N,eps) en la ventana: z, Sz, A, delta, Q*, soporte vs BV, s.
Parte 3: Prop 63.1 — BFI nunca supera a BV para K=2 (dos casos, exacto).
Parte 4: seleccion por criba (X=2e6, numpy) en z de la ventana.
Uso: python3 numerico/63_ventana_epsilon_bfi.py
"""
from fractions import Fraction as F
from itertools import combinations
import math
from math import fsum

FAIL = []

def primes_upto(n):
    if n < 2:
        return []
    s = bytearray(b'\x01') * (n + 1)
    s[0:2] = b'\x00\x00'
    for i in range(2, int(n ** 0.5) + 1):
        if s[i]:
            s[i*i:n+1:i] = b'\x00' * ((n - i*i) // i + 1)
    return [i for i in range(2, n + 1) if s[i]]

def phi(n, ps):
    r = n; m = n
    for p in ps:
        if p * p > m:
            break
        if m % p == 0:
            r = r // p * (p - 1)
            while m % p == 0:
                m //= p
    if m > 1:
        r = r // m * (m - 1)
    return r

def z_of_N_eps(E, num, den):
    """z = floor(N^(eps/2)), N=10^E, eps=num/den (exacto, aritmetica entera)."""
    target_exp = E * num
    z = int(round(10 ** (target_exp / (2 * den))))
    z = max(z, 1)
    while z ** (2 * den) > 10 ** target_exp:
        z -= 1
    while (z + 1) ** (2 * den) <= 10 ** target_exp:
        z += 1
    return z

def k2_lam(z, t):
    ps = primes_upto(z)
    lam = {1: F(1)}
    for l in ps:
        lam[l] = F(-t)
    for a, b in combinations(ps, 2):
        lam[a * b] = F(t * t)
    return ps, lam

def W_of_n(n, lam):
    return sum(v for d, v in lam.items() if n % d == 0)

print("== 1. K=2 exacto en z de la ventana (N=1e12: eps=1/4 -> z=31; eps=2/7 -> z=52) ==")
for z, t in [(31, F(1, 2)), (52, F(1, 2)), (31, F(1, 10)), (52, F(1, 10))]:
    ps, lam = k2_lam(z, t)
    Sz = sum(F(1, p) for p in ps)
    S2 = sum(F(1, p * p) for p in ps)
    A = sum(v / F(d) for d, v in lam.items())
    Aform = F(1) - t * Sz + (t * t / 2) * (Sz * Sz - S2)
    es = sorted(set([1] + ps + [a*b for a, b in combinations(ps, 2)]))
    M = {e: sum(v / F(d) for d, v in lam.items() if d % e == 0) for e in es}
    L = sum(v / F(phi(d, ps)) for d, v in lam.items())
    tot = sum(M[e] / F(phi(e, ps)) for e in es)
    ok = (A == Aform) and (L == tot)
    bad = sum(1 for n in range(1, 3001) if W_of_n(n, lam) < 0)
    minWp = min(W_of_n(p, lam) for p in primes_upto(3000))
    line = f"  z={z} t={t}: A={float(A):.4f} ident={'OK' if L == tot else 'FALLO'} W>=0[1,3000] neg={bad} minW(p)={float(minWp):.3f} soporte max d={max(lam)}<=z^2={z*z}"
    print(line)
    if not ok or bad or minWp <= 0 or max(lam) > z * z:
        FAIL.append(f"K2 exacto z={z} t={t}")
    for l in ps[:3]:
        Mform = (-t + t * t * (Sz - F(1, l))) / l
        dl = abs(M[l]) * (l - 1) / A
        print(f"    l={l}: M={float(M[l]):+.4f} form={float(Mform):+.4f} {'OK' if M[l] == Mform else 'FALLO'} delta={float(dl):.3f}")
        if M[l] != Mform:
            FAIL.append(f"K2 z={z} t={t} l={l}")

print("== 2. Parrilla (N,eps) en la ventana: z, Sz, A, delta, Q*, soporte vs BV, s ==")
print("   N      eps    z     Sz    A(t=1/2) d2(t=1/2) A(1/Sz) d2(1/Sz)  Q*/N^(2eps)  soporte<N^(1/2)?  Q*<=N^(4/7)?  s=1/eps")
for E in [12, 24, 40]:
    for num, den in [(1, 4), (13, 50), (27, 100), (2, 7)]:
        eps = num / den
        z = z_of_N_eps(E, num, den)
        ps = primes_upto(z)
        Sz = fsum(1.0 / p for p in ps)
        S2 = fsum(1.0 / (p * p) for p in ps)
        t2 = 0.5
        A2 = 1 - t2 * Sz + (t2 * t2 / 2) * (Sz * Sz - S2)
        M2 = (-t2 + t2 * t2 * (Sz - 0.5)) / 2
        d2 = abs(M2) / A2
        t = 1.0 / Sz
        A = 0.5 * (1 - S2 / (Sz * Sz))
        # delta_2 = |M_2|/A con M_2 = -1/(4 Sz^2) exacto (doc 61): factor 4, no 2
        d = 1.0 / (4 * Sz * Sz) / A
        # Q*: tres primos mayores <= z
        p1, p2, p3 = ps[-1], ps[-2], ps[-3]
        Mz = p2 * p3
        qmin = max(z * z, Mz * Mz)
        Dmax = p1 * p2
        qmax = Dmax * Dmax
        N = 10 ** E
        # chequeos exactos con enteros: Qmax <= N^(4/7) ssi Qmax^7 <= N^4;
        # Qmin > N^(1/2) es SOLO asintotico (Q* = c(z) N^(2eps), c<1): a N finito
        # el borde eps=1/4 queda por debajo; se mide, no se exige (ver parte 3).
        qmin_gt_half = (qmin * qmin > N)
        qmax_le_47 = (qmax ** 7 <= N ** 4)
        soporte_lt_half = (z ** 4 < N)  # z^2 < N^(1/2) ssi z^4 < N
        ratio_q = (math.log(qmin) / (E * math.log(10)), math.log(qmax) / (E * math.log(10)))
        print(f"   10^{E:<3} {num}/{den:<3} {z:<6} {Sz:.3f}  {A2:.4f}    {d2:.3f}      {A:.4f}   {d:.3f}     {ratio_q[0]:.3f}-{ratio_q[1]:.3f} (2eps={2*eps:.3f})  {'SI' if soporte_lt_half else 'NO'}           {'SI' if qmax_le_47 else 'NO'}         {1/eps:.2f}   Qmin>N^(1/2):{'SI' if qmin_gt_half else 'no'}")
        if not qmax_le_47:
            FAIL.append(f"Qmax>N^(4/7) E={E} eps={num}/{den}")
        if not soporte_lt_half:
            FAIL.append(f"soporte>=N^(1/2) E={E} eps={num}/{den}")
        # cross-check fsum Sz vs Fraction en z<=100
        if z <= 100:
            SzF = float(sum(F(1, p) for p in ps))
            if abs(SzF - Sz) > 1e-12:
                FAIL.append(f"Sz fsum vs exacto E={E} eps={num}/{den}")

print("== 3. Prop 63.1: BFI nunca supera a BV para K=2 (dos casos) ==")
print("  Enunciado: para la familia K=2 con soporte d <= N^eps, el error pesado")
print("  sum_d lambda_d E_d esta cubierto por BV siempre que BFI lo estaria;")
print("  la cobertura de BFI es subconjunto de la de BV. Prueba por casos:")
print("  (i) eps<=1/2: soporte <= N^(1/2) -> BV controla la suma directamente")
print("      (|lambda|<=1). (ii) eps>2/7: Q* > N^(4/7) -> la hipotesis BFI")
print("      (well-factorable de nivel <= N^(4/7)) falla. En (2/7,1/2] BFI")
print("      aplicaria a nivel Q* pero la suma solo llega al soporte < N^(1/2).")
print("  Chequeo (i): soporte z^2 <= N^(1/2) ssi z^4 <= N (exacto)")
for E in [12, 24]:
    for num, den in [(1, 4), (2, 7), (1, 2)]:
        z = z_of_N_eps(E, num, den)
        ok = (z ** 4 <= 10 ** E)
        print(f"    E={E} eps={num}/{den}: z={z} z^4<=N: {'OK' if ok else 'FALLO'}")
        if not ok:
            FAIL.append(f"Prop63.1(i) E={E} eps={num}/{den}")
print("  Chequeo (ii): eps>2/7 -> Qmin^7 > N^4, i.e. Q* > N^(4/7) (exacto)")
for E in [12, 24]:
    for num, den in [(1, 2), (3, 5)]:
        z = z_of_N_eps(E, num, den)
        ps = primes_upto(z)
        p1, p2, p3 = ps[-1], ps[-2], ps[-3]
        qmin = max(z * z, (p2 * p3) ** 2)
        ok = (qmin ** 7 > 10 ** (4 * E))
        print(f"    E={E} eps={num}/{den}: z={z} Qmin={qmin} Qmin^7>N^4: {'OK' if ok else 'FALLO'}")
        if not ok:
            FAIL.append(f"Prop63.1(ii) E={E} eps={num}/{den}")
print("  Chequeo ventana: Q* <= N^(4/7) (BFI aplica) pero soporte < N^(1/2) (BV ya cubre)")
for E in [12, 24]:
    for num, den in [(1, 4), (2, 7)]:
        z = z_of_N_eps(E, num, den)
        ps = primes_upto(z)
        p1, p2 = ps[-1], ps[-2]
        qmax = (p1 * p2) ** 2
        bfi_ok = (qmax ** 7 <= 10 ** (4 * E))
        bv_covers = (z ** 4 < 10 ** E)
        print(f"    E={E} eps={num}/{den}: BFI aplica (Qmax^7<=N^4): {'SI' if bfi_ok else 'NO'}; BV ya cubre (z^4<N): {'SI' if bv_covers else 'NO'} -> {'VACUO' if bfi_ok and bv_covers else 'FALLO'}")
        if not (bfi_ok and bv_covers):
            FAIL.append(f"ventana vacua E={E} eps={num}/{den}")
print("  Cruce finito: menor eps con Qmin > N^(1/2) (exacto, Qmin^2 > N), barrido eps=0.25..2/7")
for E in [12, 24, 40]:
    _cross = None
    for _k in range(50, 58):
        _z = z_of_N_eps(E, _k, 200)
        _ps = primes_upto(_z)
        if len(_ps) < 3:
            continue
        _qmin = max(_z * _z, (_ps[-2] * _ps[-3]) ** 2)
        if _qmin * _qmin > 10 ** E:
            _cross = _k / 200
            break
    print(f"    N=10^{E}: eps* = {_cross if _cross else '>2/7 (no cruza en la ventana)'} (Q* = c(z) N^(2eps), c<1: el borde 1/4 queda por debajo a N finito)")
print("  Parametro s: K=2 como peso de criba: s_K2 = log(z^2)/log(z) = 2 identico (barrera de paridad, f(2)=0)")
for z in [31, 52, 1000]:
    s_k2 = math.log(z * z) / math.log(z)
    print(f"    z={z}: s_K2={s_k2:.6f} {'OK' if abs(s_k2 - 2) < 1e-9 else 'FALLO'}")
    if abs(s_k2 - 2) >= 1e-9:
        FAIL.append(f"s_K2 z={z}")

print("== 4. Seleccion por criba (X=2e6) en z de la ventana ==")
import numpy as np
X = 2000000
ps_big = primes_upto(X)
Om = np.zeros(X + 1, dtype=np.int16)
for p in ps_big:
    pw = p
    while pw <= X:
        Om[pw::pw] += 1
        if pw > X // p:
            break
        pw *= p
Om = Om[1:]
for z in [7, 31, 52, 1000, 2680]:
    ps = primes_upto(z)
    rz = np.zeros(X, dtype=np.int8)
    for l in ps:
        rz[l-1::l] += 1
    rz = rz.astype(np.float64)
    for tf in [0.5, 0.1]:
        W = 1.0 - tf * rz + tf * tf * rz * (rz - 1) / 2
        m2 = (Om == 2); m3 = (Om == 3); m4 = (Om >= 4)
        e1 = float(np.mean(W[Om == 1])); e2 = float(np.mean(W[m2])); e3 = float(np.mean(W[m3])); e4 = float(np.mean(W[m4]))
        fr2 = float(np.mean(rz[m2] == 0))
        cont = float(np.mean(W[m3 | m4])) / e2 if e2 > 0 else float('nan')
        print(f"  z={z:<5} t={tf}: frac r=0={float(np.mean(rz == 0)):.3f} E[W]={float(np.mean(W)):.3f} E[W|Om=1]={e1:.3f} E[W|Om=2]={e2:.3f} E[W|Om=3]={e3:.3f} E[W|Om>=4]={e4:.3f} | P2 frac r=0={fr2:.3f} ratio P2/P3={e2/e3:.3f} contam P3+/P2={cont:.3f}")
        if abs(e1 - 1.0) > 0.01:
            FAIL.append(f"seleccion z={z} t={tf}: E[W|primo] != 1")

print()
print("RESULTADO:", "PASS" if not FAIL else f"FALLO: {FAIL}")
