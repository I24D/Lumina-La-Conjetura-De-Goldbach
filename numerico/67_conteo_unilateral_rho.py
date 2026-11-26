#!/usr/bin/env python3
"""Doc 67, etapa 2 SL-65.2: conteo unilateral con rho de Harman.
Rehace criba/circulo con rho para U(N) = sum_{p<N} rho(N-p).
Usa el proxy sin parte B: rho_tilde(n) = 1_p(n) - S_{2,3}(n) - S_{2,5}(n)
(con las regiones del doc 66, Maynard III sec. 9). El rho verdadero difiere
de rho_tilde solo por la parte B, de masa relativa O(delta) (doc 66: 0.0431).
Mide: estructura exacta de rho_tilde, lema de soporte (defecto cero si
P^-(n) < z1), densidad en clases reducidas (serie singular), conteo
unilateral exacto vs termino principal c_rho*S(N)*N/log^2 N, y la
cuantificacion L^2 de por que Cauchy no controla los arcos menores.
Uso: python3 numerico/67_conteo_unilateral_rho.py  -> PASS/FAIL + tabla.
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
print("== 0. Criba exacta X=%d ==" % X)
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

# ---------- 1. S_{2,3} y S_{2,5} por n (regiones del doc 66) ----------
z1 = X ** (1.0 / 7)
z2 = X ** (3.0 / 7)
primes = [p for p in range(2, X + 1) if is_prime[p]]
small = [p for p in primes if z1 <= p < z2]
S23 = np.zeros(X + 1, dtype=np.int16)
S25 = np.zeros(X + 1, dtype=np.int16)
Sa = Sb = 0
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
            Sa += int(np.sum(keep))
        else:
            S25[idxk] += 1
            Sb += int(np.sum(keep))
D = (S23.astype(np.int16) + S25.astype(np.int16))
rho_tilde = is_np - D
print("z1=%.2f z2=%.2f ; S_a=%d (%.4f pi) S_b=%d (%.4f pi)" %
      (z1, z2, Sa, Sa / piX, Sb, Sb / piX))
LS = (Sa + Sb) / piX
print("L_S criba = %.4f (doc 66: 0.325)" % LS)
check("L_S reproduce doc 66 (0.325 +- 0.005)", abs(LS - 0.325) < 0.005, round(LS, 4))
check("Sa/pi = 0.066 +- 0.005", abs(Sa / piX - 0.066) < 0.005, round(Sa / piX, 4))
check("Sb/pi = 0.259 +- 0.005", abs(Sb / piX - 0.259) < 0.005, round(Sb / piX, 4))

print("== 1. Estructura exacta de rho_tilde ==")
check("rho_tilde(p)=1 en todo primo", bool(np.all(rho_tilde[is_np == 1] == 1)))
check("rho_tilde(n)<=0 en todo compuesto", bool(np.all(rho_tilde[(is_np == 0)] <= 0)))
dmax = int(np.max(D))
print("max D(n) = %d ; n con D>=1: %d" % (dmax, int(np.sum(D >= 1))))
check("D(n) acotado (max <= 6)", dmax <= 6, dmax)
# Lema de soporte: si P^-(n) < z1, el factor pequeno sobrevive en el cofactor
# y rho(cofactor, p2)=0, luego S_{2,3}(n)=S_{2,5}(n)=0.
rough_small = (np.arange(X + 1) > 1) & (spf < z1)
check("lema soporte: S23=S25=0 si P^-(n)<z1",
      bool(np.all(S23[rough_small] == 0) and np.all(S25[rough_small] == 0)),
      int(np.sum(rough_small)))
check("rho_tilde(n)=1_p(n) si P^-(n)<z1",
      bool(np.all(rho_tilde[rough_small] == is_np[rough_small])))

print("== 2. Densidad en clases reducidas (serie singular) ==")
idx = np.arange(1, X + 1)
total_rho = int(np.sum(rho_tilde[1:]))
print("sum rho_tilde = %d ; (1-L_S)*pi = %.0f" % (total_rho, (1 - LS) * piX))
check("sum rho_tilde = (1-L_S)*pi", total_rho == piX - Sa - Sb, total_rho)
for p in [3, 5, 7]:
    sums = np.bincount(idx % p, weights=rho_tilde[1:].astype(float), minlength=p)
    devs = [abs(sums[a] / total_rho - 1.0 / (p - 1)) for a in range(1, p)]
    print("p=%d: clase 0 = %.0f (esperado 1) ; desv. max reducidas = %.3f" %
          (p, sums[0], max(devs)))
    check("p=%d clase 0 = 1 exacto" % p, sums[0] == 1.0, sums[0])
    check("p=%d desv. reducidas < 0.15" % p, max(devs) < 0.15, round(max(devs), 3))
for p in [11, 13]:
    sums = np.bincount(idx % p, weights=rho_tilde[1:].astype(float), minlength=p)
    devs = [abs(sums[a] / total_rho - 1.0 / (p - 1)) for a in range(1, p)]
    print("p=%d (p>z1, sin lema): clase 0 = %.0f ; desv. max reducidas = %.3f" %
          (p, sums[0], max(devs)))

# ---------- 3. Conteo unilateral exacto ----------
print("== 3. Conteo unilateral U(N) = sum_{p<N} rho_tilde(N-p) ==")
C2 = 0.66016181584686957392

def sing(N):
    s = 2.0 * C2
    n = N
    for q in primes:
        if q * q > n:
            break
        if n % q == 0:
            if q > 2:
                s *= (q - 1.0) / (q - 2.0)
            while n % q == 0:
                n //= q
    if n > 2:
        s *= (n - 1.0) / (n - 2.0)
    return s

Dcum = np.cumsum(D.astype(np.int64))
picum = np.cumsum(is_np.astype(np.int64))
primes_np = np.array(primes, dtype=np.int64)
print("N        U(N)    R(N)    S(N)   U/R     1-L_S(N)  1-L_S*kappa  kappa  r/(1-L_S)")
r_prev = None
q_prev = None
for N in [200000, 500000, 1000000, 1500000, 2000000]:
    p_arr = primes_np[primes_np < N]
    U = int(np.sum(rho_tilde[N - p_arr]))
    R = int(np.sum(is_np[N - p_arr]))
    SN = sing(N)
    lN = math.log(N)
    main0 = SN * N / lN ** 2
    LS_N = Dcum[N] / picum[N]
    # kappa(N) = <pr/phi(pr)> pesado por c(par) = rough_count(N//pr, p2):
    # corrige el conteo de primos p = N mod pr (densidad 1/phi(pr), no 1/pr).
    num = den = 0.0
    for ii, p2 in enumerate(small):
        v = math.log(p2) / LX
        for p1 in small[ii + 1:]:
            u = math.log(p1) / LX
            pr = p1 * p2
            if pr > N:
                continue
            in_a = (4.0 / 21 <= u <= 2.0 / 7) and (2.0 / 7 - u / 2 <= v <= min(u, 3.0 / 7 - u))
            in_b = (2.0 / 7 <= u <= 3.0 / 7) and (4.0 / 7 - u <= v <= min(u, (1.0 - u) / 2))
            if not (in_a or in_b):
                continue
            c = int(np.sum(spf[1:N // pr + 1] >= p2))
            num += c * pr / ((p1 - 1) * (p2 - 1))
            den += c
    kappa = num / den if den else 1.0
    r = U / main0
    q = r / (1 - LS_N)
    print("%7d %7d %7d  %.3f  %.4f  %.4f    %.4f      %.4f  %.4f" %
          (N, U, R, SN, U / R, 1 - LS_N, 1 - LS_N * kappa, kappa, q))
    check("N=%d: U<=R (rho_tilde<=1_p)" % N, U <= R, (U, R))
    check("N=%d: U>0" % N, U > 0, U)
    check("N=%d: U/R ~ 1-L_S(N) (10%%)" % N, abs(U / R - (1 - LS_N)) < 0.10,
          (round(U / R, 4), round(1 - LS_N, 4)))
    check("N=%d: U/R ~ 1-L_S*kappa (5%%)" % N, abs(U / R - (1 - LS_N * kappa)) < 0.05,
          (round(U / R, 4), round(1 - LS_N * kappa, 4)))
    if r_prev is not None:
        check("N=%d: r decrece con N (L_S crece)" % N, r < r_prev, (round(r_prev, 4), round(r, 4)))
        check("N=%d: r/(1-L_S) decrece hacia 1" % N, q < q_prev, (round(q_prev, 4), round(q, 4)))
    r_prev = r
    q_prev = q

# ---------- 4. Cuantificacion L^2: por que Cauchy no controla arcos menores ----------
print("== 4. Cuantificacion L^2 (Cauchy vs termino principal) ==")
# rho_tilde^2 = 1_p + D^2 (en primos D=0; en compuestos 1_p=0). Sin cumsum
# int64: el PC falla al asignar 3 temporales de 16 MB seguidos (cap de
# memoria del proceso). Se liberan los arreglos grandes antes de medir.
del idx, spf, Dcum, picum
import gc
gc.collect()
sum_rho2 = piX + int(np.sum(D.astype(np.int32) ** 2, dtype=np.int64))
print("sum rho_tilde^2 = %d ; /pi(X) = %.3f" % (sum_rho2, sum_rho2 / piX))
check("sum rho^2 / pi(X) en (1, 2)", 1.0 < sum_rho2 / piX < 2.0, round(sum_rho2 / piX, 3))
N = 2000000
lN = math.log(N)
cauchy = math.sqrt(piX * sum_rho2)
main_N = (1 - LS) * sing(N) * N / lN ** 2
ratio_c = cauchy / main_N
print("N=2e6: Cauchy=(pi*sum rho^2)^1/2=%.0f ; main=%.0f ; cociente=%.1f ; 2*ln N=%.1f" %
      (cauchy, main_N, ratio_c, 2 * lN))
check("Cauchy excede el termino principal >10x", ratio_c > 10, round(ratio_c, 1))

print("== RESULTADO: %s ==" % ("PASS" if not FAIL else "FAIL: " + ",".join(FAIL)))
