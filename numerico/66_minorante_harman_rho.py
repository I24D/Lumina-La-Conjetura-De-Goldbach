#!/usr/bin/env python3
"""Doc 66, etapa 1 SL-65.2: minorante rho de Harman (Maynard III, Teor. 1.3, sec. 9).
Mide: funcion de Buchstab omega (Lema 6.5), perdida L_S = (S_{2,3}+S_{2,5})/pi
por las integrales del paper (sec. 9) y por criba exacta a X=2e6,
convergencia Buchstab por pares a X creciente, fraccion con rho<0, masa B,
densidad c_rho y constante del conteo unilateral.
Regiones (sec. 9, integrales de densidad, verificadas contra L4704):
  I_a (S_{2,3}): u en [4/21,2/7], v en [2/7-u/2, min(u,3/7-u)]
                 (p1p2<=z2, p1p2^2>=z3).
  I_b (S_{2,5}): u en [2/7,3/7], v en [4/7-u, min(u,(1-u)/2)]
                 (p1p2>=z3).
Integrando: omega((1-u-v)/v)/(u*v^2). u=log p1/log x, v=log p2/log x.
Uso: python3 numerico/66_minorante_harman_rho.py  -> PASS/FAIL + tabla.
"""
import math
import numpy as np

FAIL = []

def check(name, cond, detail=""):
    print(("PASS " if cond else "FAIL ") + name + (" | " + str(detail) if detail else ""))
    if not cond:
        FAIL.append(name)

# ---------- 1. Buchstab omega (Lema 6.5, leido) ----------
# omega(u) = 1/u en [1,2]; (u*omega(u))' = omega(u-1) para u>=2.
H = 0.001
UGRID = np.arange(1.0, 10.0 + H / 2, H)
F = np.empty_like(UGRID)  # F(u) = u*omega(u)
F[UGRID <= 2.0] = 1.0
OM = np.empty_like(UGRID)
OM[UGRID <= 2.0] = 1.0 / UGRID[UGRID <= 2.0]
idx2 = int(round((2.0 - 1.0) / H))
for i in range(idx2, len(UGRID) - 1):
    j = i - int(round(1.0 / H))
    F[i + 1] = F[i] + H * 0.5 * (OM[j] + OM[j + 1])
    OM[i + 1] = F[i + 1] / UGRID[i + 1]

def omega(u):
    return float(np.interp(u, UGRID, OM))

def omega_arr(t):
    return np.interp(t, UGRID, OM, left=0.0, right=0.0)

print("== 1. Buchstab omega (Lema 6.5) ==")
check("omega(1)=1", abs(omega(1.0) - 1.0) < 1e-9, omega(1.0))
check("omega(2)=1/2", abs(omega(2.0) - 0.5) < 1e-6, omega(2.0))
check("omega(3)=(1+ln2)/3", abs(omega(3.0) - (1 + math.log(2)) / 3) < 1e-4, omega(3.0))
EG = math.exp(-0.5772156649015329)
check("omega(10)~e^-gamma", abs(omega(10.0) - EG) < 0.002, (omega(10.0), EG))

# ---------- 2. Integrales de perdida del paper (sec. 9) ----------
def paper_integral(n, use_omega=True):
    us = np.linspace(1.0 / 7, 3.0 / 7, n)
    vs = np.linspace(1.0 / 7, 3.0 / 7, n)
    U, V = np.meshgrid(us, vs, indexing="ij")
    t = (1.0 - U - V) / V
    w = omega_arr(t) if use_omega else np.ones_like(t)
    dens = np.where(t >= 1.0, w / (U * V * V), 0.0)
    ma = (U >= 4.0 / 21) & (U <= 2.0 / 7) & (V >= 2.0 / 7 - U / 2) & (V <= np.minimum(U, 3.0 / 7 - U))
    mb = (U >= 2.0 / 7) & (U <= 3.0 / 7) & (V >= 4.0 / 7 - U) & (V <= np.minimum(U, (1.0 - U) / 2))
    du = us[1] - us[0]
    dv = vs[1] - vs[0]
    return float(np.sum(dens * ma) * du * dv), float(np.sum(dens * mb) * du * dv)

Ia1, Ib1 = paper_integral(1600, True)
Ia2, Ib2 = paper_integral(3200, True)
Ia_c, Ib_c = paper_integral(1600, False)  # cota cruda omega<=1
LS_INT = Ia2 + Ib2
print("== 2. Integrales de perdida (Maynard III sec. 9) ==")
print("I_a: %.4f (1600) %.4f (3200) ; I_b: %.4f (1600) %.4f (3200)" % (Ia1, Ia2, Ib1, Ib2))
check("I_a converge <1%", abs(Ia2 - Ia1) / max(Ia2, 1e-9) < 0.01, (Ia1, Ia2))
check("I_b converge <1%", abs(Ib2 - Ib1) / max(Ib2, 1e-9) < 0.01, (Ib1, Ib2))
print("L_S integral = I_a+I_b = %.4f ; c_rho = 1-L_S = %.4f" % (LS_INT, 1 - LS_INT))
check("L_S integral < 7/8 (consistente Teor. 1.3(2))", LS_INT < 0.875, round(LS_INT, 4))
check("c_rho integral > 1/8", 1 - LS_INT > 0.125, round(1 - LS_INT, 4))
print("cota cruda (omega<=1): I_a=%.4f I_b=%.4f suma=%.4f ; 1-suma=%.4f" %
      (Ia_c, Ib_c, Ia_c + Ib_c, 1 - Ia_c - Ib_c))
check("cota cruda 1-suma >= 1/8 (Teor. 1.3(2))", 1 - Ia_c - Ib_c >= 0.125, round(1 - Ia_c - Ib_c, 4))

# ---------- 3. Criba exacta a X=2e6 ----------
X = 2_000_000
LX = math.log(X)
print("== 3. Criba exacta X=%d ==" % X)
spf = np.arange(X + 1, dtype=np.int64)
spf[0] = 0
spf[1] = 10 ** 9
lim = int(X ** 0.5) + 1
is_prime = bytearray(b"\x01") * (X + 1)
is_prime[0:2] = b"\x00\x00"
for p in range(2, lim):
    if is_prime[p]:
        is_prime[p * p:X + 1:p] = b"\x00" * ((X - p * p) // p + 1)
        spf[p * p:X + 1:p] = np.minimum(spf[p * p:X + 1:p], p)
piX = int(sum(is_prime))
print("pi(X) = %d" % piX)
check("pi(2e6)=148933", piX == 148933, piX)

def rough_count(Y, z):
    if Y < 1:
        return 0
    return int(np.sum(spf[1:Y + 1] >= z))

# Lema 6.5: #{n<=X : P^-(n)>=z} ~ X*omega(u)/log z, u = ln X/ln z (u acotado).
for z in [7, 101, 502]:
    u = math.log(X) / math.log(z)
    cnt = rough_count(X, z)
    pred = X * omega(u) / math.log(z)
    rel = abs(cnt - pred) / pred
    print("Lema 6.5 z=%d u=%.2f exacto=%d pred=%.0f rel=%.3f" % (z, u, cnt, pred, rel))
    check("Lema 6.5 z=%d rel<0.15" % z, rel < 0.15, round(rel, 3))

# S_{2,3} y S_{2,5} exactos por suma sobre pares, regiones del paper.
z1 = X ** (1.0 / 7)
z2 = X ** (3.0 / 7)
primes = [p for p in range(2, X + 1) if is_prime[p]]
small = [p for p in primes if z1 <= p < z2]
print("z1=%.2f z2=%.2f ; primos en [z1,z2): %d" % (z1, z2, len(small)))
Sa = 0
Sb = 0
for ii, p2 in enumerate(small):
    v = math.log(p2) / LX
    for p1 in small[ii + 1:]:
        u = math.log(p1) / LX
        pr = p1 * p2
        if pr > X:
            continue
        c = rough_count(X // pr, p2)
        if 4.0 / 21 <= u <= 2.0 / 7 and 2.0 / 7 - u / 2 <= v <= min(u, 3.0 / 7 - u):
            Sa += c
        if 2.0 / 7 <= u <= 3.0 / 7 and 4.0 / 7 - u <= v <= min(u, (1.0 - u) / 2):
            Sb += c
LS_SIEVE = (Sa + Sb) / piX
print("S_a=%d (%.4f pi) S_b=%d (%.4f pi) ; L_S criba = %.4f ; L_S integral = %.4f" %
      (Sa, Sa / piX, Sb, Sb / piX, LS_SIEVE, LS_INT))
check("L_S criba en (0, 7/8)", 0 < LS_SIEVE < 0.875, round(LS_SIEVE, 4))
check("L_S criba < L_S integral (efecto escala finita, signo esperado)", LS_SIEVE < LS_INT,
      (round(LS_SIEVE, 4), round(LS_INT, 4)))

# Fraccion con S_{2,3}+S_{2,5}>=1 (=> rho(n)<0): n con P^-(n)=p2 y otro factor
# primo p1 en (p2,z2) con (p1,p2) en las regiones del paper.
neg = np.zeros(X + 1, dtype=bool)
for ii, p2 in enumerate(small):
    v = math.log(p2) / LX
    for p1 in small[ii + 1:]:
        u = math.log(p1) / LX
        pr = p1 * p2
        if pr > X:
            continue
        in_a = (4.0 / 21 <= u <= 2.0 / 7) and (2.0 / 7 - u / 2 <= v <= min(u, 3.0 / 7 - u))
        in_b = (2.0 / 7 <= u <= 3.0 / 7) and (4.0 / 7 - u <= v <= min(u, (1.0 - u) / 2))
        if not (in_a or in_b):
            continue
        idx = np.arange(pr, X + 1, pr)
        m = idx // pr
        keep = spf[m] >= p2
        neg[idx[keep]] = True
frac_neg = float(np.sum(neg[1:])) / X
print("n<=X con S23+S25>=1 (rho<0): %.4f" % frac_neg)
check("fraccion rho<0 > 0.01", frac_neg > 0.01, round(frac_neg, 4))
check("fraccion rho<0 < 0.50", frac_neg < 0.50, round(frac_neg, 4))

# Masa B (Lema 6.8): primos en B_2, d=0.001.
d = 0.001
Bints = [(X ** (3.0 / 7), X ** (3.0 / 7 + 8 * d)),
         (X ** (0.5 - 4 * d), X ** (0.5 + 4 * d)),
         (X ** (4.0 / 7 - 8 * d), X ** (4.0 / 7))]
bprimes = [p for p in primes if any(a <= p <= b for a, b in Bints)]
hasB = np.zeros(X + 1, dtype=bool)
for p in bprimes:
    hasB[p:X + 1:p] = True
fracB = float(np.sum(hasB[1:])) / X
print("primos en B: %d ; n<=X con factor primo en B: %.4f" % (len(bprimes), fracB))
check("masa B O(delta) < 0.05", fracB < 0.05, round(fracB, 4))

# ---------- 4. Convergencia Buchstab por pares a X creciente ----------
def primes_upto(n):
    sieve = bytearray(b"\x01") * (n + 1)
    sieve[0:2] = b"\x00\x00"
    for p in range(2, int(n ** 0.5) + 1):
        if sieve[p]:
            sieve[p * p:n + 1:p] = b"\x00" * ((n - p * p) // p + 1)
    return np.array([p for p in range(2, n + 1) if sieve[p]], dtype=float)

def buchstab_pair_sum(XX):
    LXX = math.log(XX)
    zz1 = XX ** (1.0 / 7)
    zz2 = XX ** (3.0 / 7)
    pr = primes_upto(int(zz2) + 1)
    sm = pr[(pr >= zz1) & (pr < zz2)]
    tot = 0.0
    for ii, p2 in enumerate(sm):
        v = math.log(p2) / LXX
        p1s = sm[ii + 1:]
        if len(p1s) == 0:
            continue
        u = np.log(p1s) / LXX
        Y = XX / (p1s * p2)
        ok = Y >= p2
        if not np.any(ok):
            continue
        up = np.where(ok, np.log(np.maximum(Y, 1.0)) / math.log(p2), 0.0)
        cb = np.where(ok, Y * omega_arr(up) / math.log(p2), 0.0)
        ma = ok & (u >= 4.0 / 21) & (u <= 2.0 / 7) & (v >= 2.0 / 7 - u / 2) & (v <= np.minimum(u, 3.0 / 7 - u))
        mb = ok & (u >= 2.0 / 7) & (u <= 3.0 / 7) & (v >= 4.0 / 7 - u) & (v <= np.minimum(u, (1.0 - u) / 2))
        tot += float(np.sum(cb[ma | mb]))
    return tot / (XX / LXX)

print("== 4. Convergencia Buchstab por pares ==")
LS_big = None
for XX in [2e6, 1e8, 1e10, 1e12]:
    Ls = buchstab_pair_sum(XX)
    print("X=%.0e: L_S Buchstab-pares=%.4f (integral=%.4f)" % (XX, Ls, LS_INT))
    if XX == 1e12:
        LS_big = Ls
check("L_S(1e12) a <10% del integral", abs(LS_big - LS_INT) / LS_INT < 0.10,
      (round(LS_big, 4), round(LS_INT, 4)))
check("L_S crece con X hacia el integral", LS_big > LS_SIEVE, (round(LS_big, 4), round(LS_SIEVE, 4)))

# ---------- 5. Densidad c_rho y constantes ----------
print("== 5. Densidad c_rho y constantes ==")
c_rho = 1 - LS_INT
print("c_rho = 1 - L_S - O(delta) = %.4f - O(delta), delta<1/1000" % c_rho)
print("cota cruda del paper: O(delta) <= %.4f (pues 1-suma_cruda-O(delta)>=1/8)" % (1 - Ia_c - Ib_c - 0.125))
print("c_rho asintotico en [%.3f, %.3f]" % (c_rho - (1 - Ia_c - Ib_c - 0.125), c_rho))
print("defecto 1-c_rho asintotico = %.4f (+O(delta)) ; cota rigurosa <= 7/8" % LS_INT)
print("constante unilateral c_rho en [%.3f, %.3f] (doc 65 decia ~1/64=0.0156)" %
      (c_rho - (1 - Ia_c - Ib_c - 0.125), c_rho))
print("constante rho*rho c_rho^2 (NO VALIDA: rho<0 en compuestos): [%.4f, %.4f]" %
      ((1.0 / 8) ** 2, c_rho ** 2))
check("c_rho asintotico > 1/8", c_rho > 0.125, round(c_rho, 4))
check("c_rho asintotico < 1", c_rho < 1.0, round(c_rho, 4))

print("== RESULTADO: %s ==" % ("PASS" if not FAIL else "FAIL: " + ",".join(FAIL)))
