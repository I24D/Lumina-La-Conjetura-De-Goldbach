#!/usr/bin/env python3
"""Doc 69, etapa 1: lema del momento mixto, via Vaughan/Heath-Brown de S_rho.
Descomposicion S_rho_tilde = S_1 - S_{2,3} - S_{2,5} (doc 67). Mide: normas
L^2 de las piezas (delgadas o gruesas), la barrera de Cauchy por piezas
(triangulo L^2: suma ||T_j||_2 >= ||S_rho||_2), la aritmetica sup x L^1
(sup requerido vs Vinogradov) y la diagonal del cuarto momento mixto.
Uso: python3 numerico/69_momento_mixto_vaughan.py -> PASS/FAIL.
"""
import math
import numpy as np

FAIL = []

def check(name, cond, detail=""):
    print(("PASS " if cond else "FAIL ") + name + (" | " + str(detail) if detail else ""))
    if not cond:
        FAIL.append(name)

X = 2_000_000
LX = math.log(X)
print("== 0. Criba X=%d (doc 67) ==" % X)
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
is_np = np.array(is_prime, dtype=np.int16)
piX = int(sum(is_prime))
check("pi(2e6)=148933", piX == 148933, piX)

z1 = X ** (1.0 / 7)
z2 = X ** (3.0 / 7)
primes = [p for p in range(2, X + 1) if is_prime[p]]
small = [p for p in primes if z1 <= p < z2]
S23 = np.zeros(X + 1, dtype=np.int16)
S25 = np.zeros(X + 1, dtype=np.int16)
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
        keep = spf[idx // pr] >= p2
        idxk = idx[keep]
        if in_a:
            S23[idxk] += 1
        else:
            S25[idxk] += 1
D = S23.astype(np.int32) + S25.astype(np.int32)
print("== 1. Normas L^2 de las piezas (delgadas o gruesas) ==")
s23 = int(np.sum(S23.astype(np.int32) ** 2, dtype=np.int64))
s25 = int(np.sum(S25.astype(np.int32) ** 2, dtype=np.int64))
sd2 = int(np.sum(D ** 2, dtype=np.int64))
sum_rho2 = piX + sd2
print("sum S23^2=%d (%.3f pi) ; sum S25^2=%d (%.3f pi) ; sum D^2=%d ; sum rho~^2=%d (%.3f pi)"
      % (s23, s23 / piX, s25, s25 / piX, sd2, sum_rho2, sum_rho2 / piX))
check("piezas delgadas: sum S23^2, sum S25^2 < 1.5 pi", s23 < 1.5 * piX and s25 < 1.5 * piX,
      (round(s23 / piX, 3), round(s25 / piX, 3)))
check("sum rho~^2 = 1.467 pi (doc 67)", abs(sum_rho2 / piX - 1.467) < 0.01, round(sum_rho2 / piX, 3))
n1 = math.sqrt(piX)
n23 = math.sqrt(s23)
n25 = math.sqrt(s25)
nrho = math.sqrt(sum_rho2)
print("||S1||=%.1f ||S23||=%.1f ||S25||=%.1f ; suma=%.1f >= ||S_rho~||=%.1f"
      % (n1, n23, n25, n1 + n23 + n25, nrho))
check("triangulo L^2: suma ||T_j|| >= ||S_rho~||", n1 + n23 + n25 >= nrho,
      round((n1 + n23 + n25) / nrho, 3))

print("== 2. Barrera de Cauchy por piezas (Prop. 69.1) ==")
N = X
lN = math.log(N)
target = N / lN ** 2
c23 = math.sqrt(piX * s23)
c25 = math.sqrt(piX * s25)
tot = piX + c23 + c25
cauchy67 = math.sqrt(piX * sum_rho2)
print("int|S1|^2<=pi=%.0f ; C23=%.0f ; C25=%.0f ; total por piezas=%.0f" % (piX, c23, c25, tot))
print("target N/log^2 N=%.0f ; Cauchy doc67=%.0f ; cociente total/target=%.1f ; Cauchy67/target=%.1f"
      % (target, cauchy67, tot / target, cauchy67 / target))
check("Cauchy por piezas >= Cauchy doc 67 (no mejora)", tot >= cauchy67, round(tot / cauchy67, 3))
check("total por piezas >> target (>5x)", tot / target > 5, round(tot / target, 1))
check("termino int|S1|^2 solo ya es ~pi = %.1f x target" % (piX / target),
      piX / target > 5, round(piX / target, 1))

print("== 3. Aritmetica sup x L^1 (Prop. 69.2) ==")
# Sup requerido: int|S1| <= (N/log N)^1/2, luego sup = o(N^1/2/log^3/2 N).
# Vinogradov (A fijo, asintotico): sup << N/log^A N; cociente al requerido =
# N^1/2/(log N)^(A-3/2) -> infinito. Se usa A=3 (realista a N accesible;
# A=10 solo vale para N >> 10^39, fuera de regimen: la forma cruda engana).
for NN in [10 ** 12, 10 ** 24]:
    l = math.log(NN)
    req = NN ** 0.5 / l ** 1.5
    vin3 = NN / l ** 3
    ideal = NN ** 0.5
    print("N=1e%d: sup requerido=%.3g ; Vinogradov A=3=%.3g ; cociente=%.3g ; "
          "ideal raiz x L^1 / target=%.1f"
          % (round(math.log10(NN)), req, vin3, vin3 / req, (ideal * (NN / l) ** 0.5) / (NN / l ** 2)))
    check("N=1e%d: Vinogradov A=3 esta >> por encima del sup requerido (>100x)" % round(math.log10(NN)),
          vin3 / req > 100, "%.2g" % (vin3 / req))
    check("N=1e%d: ni la cancelacion de raiz cuadrada ideal basta (>100x target)" % round(math.log10(NN)),
          (ideal * (NN / l) ** 0.5) / (NN / l ** 2) > 100, round((ideal * (NN / l) ** 0.5) / (NN / l ** 2), 1))

print("== 4. Diagonal del cuarto momento mixto (Prop. 69.3) ==")
need4 = N ** 2 / lN ** 4
d23 = piX * s23
d25 = piX * s25
print("diagonal S23: pi*sum S23^2=%.3g ; diagonal S25: %.3g ; necesario o(N^2/log^4)=%.3g"
      % (d23, d25, need4))
# Diagonal = pi * sum S^2 ~ c*N^2/log^2 N frente a o(N^2/log^4 N): sobra
# c*(log N)^2 -> infinito; a X=2e6 el factor es 17x (S23) y 96x (S25).
check("diagonal S23 >> necesario (>10x)", d23 / need4 > 10, round(d23 / need4, 1))
check("diagonal S25 >> necesario (>10x)", d25 / need4 > 10, round(d25 / need4, 1))

print("== RESULTADO: %s ==" % ("PASS" if not FAIL else "FAIL: " + ",".join(FAIL)))
