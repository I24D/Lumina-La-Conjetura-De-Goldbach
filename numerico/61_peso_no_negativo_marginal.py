#!/usr/bin/env python3
"""Doc 61: peso no negativo W(n)>=0 con cancelacion marginal (cambio de marco).
Parte 1: identidades exactas con Fraction. Parte 2: densidades rho_l por criba.
Uso: python3 numerico/61_peso_no_negativo_marginal.py
"""
from fractions import Fraction as F
from itertools import combinations

FAIL = []

def primes_upto(n):
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

def margins(lam, es):
    return {e: sum(v / F(d) for d, v in lam.items() if d % e == 0) for e in es}

def check_identity(name, lam, es, ps):
    L = sum(v / F(phi(d, ps)) for d, v in lam.items())
    M = margins(lam, es)
    tot = sum(M[e] / F(phi(e, ps)) for e in es)
    ok = (L == tot)
    print(f"  [{name}] L={float(L):.6f} sum M_e/phi(e)={float(tot):.6f} ident={'OK' if ok else 'FALLO'}")
    if not ok:
        FAIL.append(name + ' identidad')
    return L, M

def rough_lam(z):
    ps = primes_upto(z)
    lam = {1: F(1)}
    for k in range(1, len(ps) + 1):
        for c in combinations(ps, k):
            d = 1
            for p in c:
                d *= p
            lam[d] = F((-1) ** k)
    return ps, lam

def k2_lam(z, t):
    ps = primes_upto(z)
    lam = {1: F(1)}
    for l in ps:
        lam[l] = F(-t)
    for a, b in combinations(ps, 2):
        lam[a * b] = F(t * t)
    return ps, lam

def wt_lam(z, t):
    ps = primes_upto(z)
    lam = {1: F(1)}
    for k in range(1, len(ps) + 1):
        for c in combinations(ps, k):
            d = 1
            for p in c:
                d *= p
            lam[d] = F((-t) ** k)
    return ps, lam

def W_of_n(n, lam):
    return sum(v for d, v in lam.items() if n % d == 0)

print("== 1. Indicador rugoso: M_l == -A/(l-1) exacto ==")
for z in [5, 7, 11]:
    ps, lam = rough_lam(z)
    es = sorted(lam.keys())  # todos los e squarefree z-suaves (M_e=0 fuera)
    A = sum(v / F(d) for d, v in lam.items())
    L, M = check_identity(f"rugoso z={z}", lam, es, ps)
    for l in ps:
        ok = (M[l] == -A / (l - 1))
        print(f"  l={l}: M={float(M[l]):+.6f} -A/(l-1)={float(-A/(l-1)):+.6f} {'OK' if ok else 'FALLO'}")
        if not ok:
            FAIL.append(f"rugoso z={z} l={l}")
    bad = sum(1 for n in range(1, 301) if W_of_n(n, lam) not in (F(0), F(1)))
    print(f"  W(n) en {{0,1}} para n<=300: {'OK' if bad == 0 else f'FALLO {bad}'}")
    if bad:
        FAIL.append(f"rugoso z={z} W 0/1")

print("== 2. Trivial W==1: margenes cero ==")
lam = {1: F(1)}
L, M = check_identity("trivial", lam, [1, 2, 3, 5], primes_upto(5))
print(f"  A=1 M_2={M[2]} M_3={M[3]} {'OK' if M[2] == 0 and M[3] == 0 else 'FALLO'}")
if M[2] != 0 or M[3] != 0:
    FAIL.append("trivial margenes")

print("== 3. Familia W_t (sin truncar): M_l == -tA/(l-t) exacto ==")
for z, t in [(5, F(1, 2)), (7, F(1, 2))]:
    ps, lam = wt_lam(z, t)
    es = sorted(lam.keys())  # todos los e squarefree z-suaves (M_e=0 fuera)
    A = sum(v / F(d) for d, v in lam.items())
    L, M = check_identity(f"W_t z={z} t={t}", lam, es, ps)
    for l in ps:
        ok = (M[l] == -t * A / (l - t))
        print(f"  l={l}: M={float(M[l]):+.6f} -tA/(l-t)={float(-t*A/(l-t)):+.6f} {'OK' if ok else 'FALLO'}")
        if not ok:
            FAIL.append(f"W_t z={z} l={l}")

print("== 4. Familia K=2 (truncada w<=2): formulas, W>=0, trade-off delta~t ==")
for z, t in [(5, F(1, 2)), (7, F(1, 2)), (11, F(1, 2)),
             (7, F(1, 10)), (11, F(1, 10))]:
    ps, lam = k2_lam(z, t)
    Sz = sum(F(1, p) for p in ps)
    S2 = sum(F(1, p * p) for p in ps)
    es = sorted(set([1] + ps + [a*b for a, b in combinations(ps, 2)]))
    A = sum(v / F(d) for d, v in lam.items())
    Aform = F(1) - t * Sz + (t * t / 2) * (Sz * Sz - S2)
    L, M = check_identity(f"K2 z={z} t={t}", lam, es, ps)
    okA = (A == Aform)
    print(f"  A={float(A):.4f} formula={float(Aform):.4f} {'OK' if okA else 'FALLO'}")
    if not okA:
        FAIL.append(f"K2 z={z} t={t} A")
    for l in ps:
        Mform = (-t + t * t * (Sz - F(1, l))) / l
        ok = (M[l] == Mform)
        dl = abs(M[l]) * (l - 1) / A
        print(f"  l={l}: M={float(M[l]):+.4f} formula={float(Mform):+.4f} {'OK' if ok else 'FALLO'} delta={float(dl):.3f}")
        if not ok:
            FAIL.append(f"K2 z={z} t={t} l={l}")
    bad = 0; minW = None; minWp = None
    for n in range(1, 3001):
        w = W_of_n(n, lam)
        if minW is None or w < minW:
            minW = w
        if w < 0:
            bad += 1
    for p in primes_upto(3000):
        w = W_of_n(p, lam)
        if minWp is None or w < minWp:
            minWp = w
    print(f"  W>=0 [1,3000]: negativos={bad} minW={float(minW):.4f} minW(p)={float(minWp):.4f} {'OK' if bad == 0 and minWp > 0 else 'FALLO'}")
    if bad or minWp <= 0:
        FAIL.append(f"K2 z={z} t={t} W>=0")
    print(f"  soporte max d={max(lam)} <= z^2={z*z} {'OK' if max(lam) <= z*z else 'FALLO'}")

print("== 5. K=2 con t=1/Sz: M_l == -1/(l^2 Sz^2) exacto, delta grande ==")
for z in [7, 11, 13]:
    ps, lam = k2_lam(z, F(1))
    Sz = sum(F(1, p) for p in ps)
    t = F(1) / Sz
    ps, lam = k2_lam(z, t)
    A = sum(v / F(d) for d, v in lam.items())
    _, M = check_identity(f"K2 t=1/Sz z={z}", lam, sorted(set([1] + ps + [a*b for a, b in combinations(ps, 2)])), ps)
    for l in ps:
        ok = (M[l] == -F(1) / (l * l * Sz * Sz))
        dl = abs(M[l]) * (l - 1) / A
        print(f"  l={l}: M={float(M[l]):+.4f} -1/(l^2Sz^2)={float(-F(1)/(l*l*Sz*Sz)):+.4f} {'OK' if ok else 'FALLO'} delta={float(dl):.3f}")
        if not ok:
            FAIL.append(f"K2 1/Sz z={z} l={l}")

print("== 6. Densidades rho_l por criba (X=2e6, numpy) ==")
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
w7 = np.zeros(X + 1, dtype=np.int8)
for l in [2, 3, 5, 7]:
    w7[l::l] += 1
n = np.arange(1, X + 1, dtype=np.float64)
inv = 1.0 / n

def rho(W):
    tot = float(np.sum(W * inv))
    out = {}
    for l in [2, 3, 5, 7]:
        out[l] = float(np.sum(W[l-1::l] * inv[l-1::l])) / tot
    return out

Wrough = (w7[1:] == 0).astype(np.float64)
Wom2 = (Om[1:] <= 2).astype(np.float64)
r = rho(Wrough)
print("  rugoso z=7: " + " ".join(f"rho_{l}={r[l]:.4f} (esperado 0)" for l in [2, 3, 5, 7]))
for l in [2, 3, 5, 7]:
    if abs(r[l]) > 1e-9:
        FAIL.append(f"rho rugoso l={l}")
r = rho(Wom2)
print("  Omega<=2:   " + " ".join(f"rho_{l}={r[l]:.4f} (1/l={1/l:.3f})" for l in [2, 3, 5, 7]))
for l in [2, 3, 5, 7]:
    if not r[l] < 1.0 / l - 0.05:
        FAIL.append(f"rho Omega<=2 l={l} no deficiente")
# K=2 t=1/Sz, z=7: rho predicho por Prop 61.2
z = 7
ps7 = primes_upto(z)
Sz = sum(F(1, p) for p in ps7)
t = F(1) / Sz
_, lam = k2_lam(z, t)
A = sum(v / F(d) for d, v in lam.items())
_, M = check_identity("K2 t=1/Sz z=7 (rho)", lam,
                      sorted(set([1] + ps7 + [a*b for a, b in combinations(ps7, 2)])), ps7)
tf = float(t)
w7f = w7[1:].astype(np.float64)
WK2 = 1.0 - tf * w7f + tf * tf * w7f * (w7f - 1) / 2
r = rho(WK2)
print(f"  K2 z=7 t=1/Sz={tf:.4f} (A={float(A):.4f}):")
for l in [2, 3, 5, 7]:
    pred = 1.0 / l + float(M[l]) * (1 - 1.0 / l) / float(A)
    print(f"    rho_{l} medido={r[l]:.4f} predicho={pred:.4f} dif={abs(r[l]-pred):.4f}")
    if abs(r[l] - pred) > 0.10:
        FAIL.append(f"rho K2 l={l}")

print("== 7. Termino principal fisico L_N (solo (d,N)=1) vs L ==")
for z, t in [(7, F(1, 10)), (7, F(1, 2)), (11, F(1, 10))]:
    ps, lam = k2_lam(z, t)
    L = sum(v / F(phi(d, ps)) for d, v in lam.items())
    for S in [{2}, {2, 5}, {2, 3, 5, 7}, {3, 5}]:
        LN = sum(v / F(phi(d, ps)) for d, v in lam.items()
                 if all(d % l != 0 for l in S))
        print(f"  z={z} t={t} N con primos {sorted(S)}: L={float(L):.4f} L_N={float(LN):.4f} correccion={float(LN-L):+.4f}")
    # rugoso: L=0 pero L_N>0 si N par
    _, lamr = rough_lam(z)
    Lr = sum(v / F(phi(d, ps)) for d, v in lamr.items())
    LNr = sum(v / F(phi(d, ps)) for d, v in lamr.items() if d % 2 != 0)
    print(f"  rugoso z={z}: L={float(Lr):.4f} L_N(N par)={float(LNr):.4f}")

print()
print("RESULTADO:", "PASS" if not FAIL else f"FALLO: {FAIL}")
