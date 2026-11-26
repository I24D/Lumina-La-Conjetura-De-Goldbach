#!/usr/bin/env python3
"""Doc 76: mapa de rutas restantes al momento mixto (sin correlacion de pares).
Reusa la criba exacta de rho_tilde del doc 67 (X=2e6) y mide:
 (1) Prop 76.3: la norma L^2 esta forzada por los factores locales
     (piso de Cauchy c*N/log N, 6.75x el objetivo a X=2e6 con c=0.465);
 (2) Prop 76.1/76.2: el sup necesitaria Theta<1/2 (refutado: gamma_1
     tiene beta=1/2) y la densidad de ceros no toca la diagonal
     (ceros en linea contados en el archivo del doc 72);
 (3) escala diagonal del cuarto momento mixto (245x) y cotas inferiores
     rigurosas de las piezas no negativas del defecto (doc 69.3);
 (4) escala de (74.4): sum U(N)^2 ~ X^3/log^4 X (via FFT).
Uso: python3 numerico/76_rutas_momento_mixto.py -> PASS/FAIL.
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
print("== 0. Criba exacta X=%d (misma que doc 67) ==" % X)
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
LS = (Sa + Sb) / piX
print("L_S criba = %.4f (doc 66: 0.325)" % LS)
check("L_S reproduce doc 66", abs(LS - 0.325) < 0.005, round(LS, 4))

print("== 1. Prop 76.3: L^2 forzada por los factores locales ==")
sum_rho2 = piX + int(np.sum(D.astype(np.int32) ** 2, dtype=np.int64))
print("sum rho_tilde^2 = %d (doc 67: 218454)" % sum_rho2)
check("sum rho^2 = 218454", sum_rho2 == 218454, sum_rho2)
c_asym = 0.465  # doc 66, c_rho asintotico
target = X / LX ** 2
floor_cauchy = c_asym * X / LX
print("piso Cauchy c*N/log N = %.0f ; objetivo N/log^2 N = %.0f ; cociente = %.2f (= c*log N)"
      % (floor_cauchy, target, floor_cauchy / target))
check("piso Cauchy = 6.75x el objetivo", abs(floor_cauchy / target - 6.75) < 0.05,
      round(floor_cauchy / target, 2))
cauchy = math.sqrt(piX * sum_rho2)
print("Cauchy real = %.0f = %.1fx objetivo ; >= piso %.0f" % (cauchy, cauchy / target, floor_cauchy))
check("Cauchy real >= piso forzado", cauchy >= floor_cauchy, (round(cauchy), round(floor_cauchy)))
check("Cauchy real ~19x objetivo (doc 67)", abs(cauchy / target - 19.0) < 0.5,
      round(cauchy / target, 1))
# Forzado por clase: por Cauchy en cada clase reducida mod l,
# sum_{p=a(l)} rho(p)^2 >= (sum_{p=a(l)} rho(p))^2 / pi(X;l,a).
idx = np.arange(1, X + 1)
for ell in [3, 5, 7]:
    tot = 0.0
    for a in range(1, ell):
        m = (idx % ell == a) & (is_np[1:] == 1)
        cnt = int(np.sum(m))
        s = float(np.sum(rho_tilde[1:][m]))
        tot += s * s / cnt if cnt else 0.0
    print("l=%d: sum_clases (sum rho)^2/pi = %.0f ; pi(X) = %d (falta 1: el primo l en clase 0)" % (ell, tot, piX))
    check("l=%d: piso por clases = pi(X)-1" % ell, abs(tot - (piX - 1)) < 1.0, round(tot))
    check("l=%d: piso por clases >> c^2*N/log N" % ell, tot > c_asym ** 2 * X / LX,
          round(c_asym ** 2 * X / LX))
del idx
import gc
gc.collect()

print("== 2. Prop 76.1/76.2: sup refutado; diagonal fuera de la densidad ==")
try:
    import mpmath as mp
    g1 = mp.zetazero(1)
    print("gamma_1 = %.10f + %.10fi" % (mp.re(g1), mp.im(g1)))
    check("beta(gamma_1) = 1/2 exacto (Theta<1/2 refutado)", abs(float(mp.re(g1)) - 0.5) < 1e-12,
          float(mp.re(g1)))
    check("gamma_1 = 14.1347251417", abs(float(mp.im(g1)) - 14.1347251417) < 1e-6,
          round(float(mp.im(g1)), 10))
except ImportError:
    print("SKIP mpmath no disponible (check manual: gamma_1 = 14.1347251417, beta=1/2)")
for Nbig, lo, hi in [(10 ** 12, 140, 150), (10 ** 24, 400, 420)]:
    gap = math.log(Nbig) ** 1.5
    print("N=1e%d: hueco sup ideal (log N)^(3/2) = %.1f" % (12 if Nbig == 10 ** 12 else 24, gap))
    check("N=1e%d: hueco sup ideal en rango" % (12 if Nbig == 10 ** 12 else 24), lo < gap < hi, round(gap, 1))
# Ceros en linea hasta T=100 (archivo del doc 72): la densidad N(sigma,T)
# con sigma>1/2 no los cuenta; son la diagonal.
try:
    zeros = [float(l.split()[0]) for l in open("72_ceros_zeta.txt") if l.strip()]
except FileNotFoundError:
    zeros = [float(l.split()[0]) for l in open("/home/hatch/workspace/goldbach/72_ceros_zeta.txt") if l.strip()]
cnt100 = sum(1 for g in zeros if g <= 100.0)
NT = (100.0 / (2 * math.pi)) * (math.log(100.0 / (2 * math.pi)) - 1) + 0.875
print("ceros con gamma<=100: %d ; N(T) Riemann-von Mangoldt = %.1f" % (cnt100, NT))
check("ceros en linea hasta 100 ~ N(T)", abs(cnt100 - NT) < 3, (cnt100, round(NT, 1)))

print("== 3. Cuarto momento: escala diagonal y cotas de piezas no negativas ==")
need4 = X ** 2 / LX ** 4
diag_rho = piX ** 2
print("diagonal rho (sum_p rho)^2 = %.3e ; necesario o(N^2/log^4 N) = %.3e ; cociente = %.0f (ESCALA, no cota inf: rho firma)"
      % (diag_rho, need4, diag_rho / need4))
check("diagonal rho ~245x lo necesario (escala)", abs(diag_rho / need4 - 245) < 10,
      round(diag_rho / need4))
s23sq = int(np.sum(S23.astype(np.int32) ** 2, dtype=np.int64))
s25sq = int(np.sum(S25.astype(np.int32) ** 2, dtype=np.int64))
d23 = piX * s23sq
d25 = piX * s25sq
print("diagonal S23 = %.3e (%.1fx) ; diagonal S25 = %.3e (%.1fx) ; necesario = %.3e"
      % (d23, d23 / need4, d25, d25 / need4, need4))
check("diagonal S23 ~17.3x (doc 69.3)", abs(d23 / need4 - 17.3) < 0.6, round(d23 / need4, 1))
check("diagonal S25 ~95.8x (doc 69.3)", abs(d25 / need4 - 95.8) < 1.5, round(d25 / need4, 1))
del S23, S25
gc.collect()

print("== 4. Escala de (74.4): sum U(N)^2 ~ X^3/log^4 X ==")
X0 = 100_000
n = np.arange(0, 2 * X0 + 1)
a = is_np[:2 * X0 + 1].astype(float)
b = rho_tilde[:2 * X0 + 1].astype(float)
L = 1
while L < len(a) + len(b):
    L *= 2
conv = np.fft.irfft(np.fft.rfft(a, L) * np.fft.rfft(b, L), L)[:2 * X0 + 1]
U = conv[2 * np.arange(X0 // 2, X0 + 1)]  # N pares en [X0, 2*X0]
convR = np.fft.irfft(np.fft.rfft(a, L) * np.fft.rfft(a, L), L)[:2 * X0 + 1]
R = convR[2 * np.arange(X0 // 2, X0 + 1)]
lx0 = math.log(X0)
scale = X0 ** 3 / lx0 ** 4
sU = float(np.sum(U ** 2))
sR = float(np.sum(R ** 2))
print("X0=1e5: sum U^2 = %.3e ; sum R^2 = %.3e ; X^3/log^4 X = %.3e ; cocientes %.2f / %.2f"
      % (sU, sR, scale, sU / scale, sR / scale))
check("sum U^2 orden X^3/log^4 X (0.1-10x)", 0.1 < sU / scale < 10, round(sU / scale, 2))
check("sum R^2 orden X^3/log^4 X (0.1-10x)", 0.1 < sR / scale < 10, round(sR / scale, 2))
check("U(N)<=R(N) en todo N del rango", bool(np.all(U <= R + 1e-6)), "")

print("== 5. D1/D2: el split de valores grandes es vacuo; el sup en q no baja ==")
for Nbig in [10 ** 12, 10 ** 24]:
    lN = math.log(Nbig)
    Vstar = Nbig ** 0.5 / lN ** 1.5
    meas = (Nbig / lN) / Vstar ** 2  # Chebyshev: ||S1||_2^2 / V^2
    print("N=1e%d: V* = N^(1/2)/(log N)^(3/2) ; cota de medida del conjunto grande = (log N)^2 = %.1f (>1: vacua)"
          % (12 if Nbig == 10 ** 12 else 24, meas))
    check("N=1e%d: medida de {|S1|>V*} vacua (>1)" % (12 if Nbig == 10 ** 12 else 24), meas > 1,
          round(meas, 1))
    qfloor = Nbig ** 0.8 / Nbig ** 0.5  # sup dependiente de q >= N^(4/5) vs ideal N^(1/2)
    print("N=1e%d: piso del sup en q (N^(4/5)/N^(1/2) = N^(3/10)) = %.0f"
          % (12 if Nbig == 10 ** 12 else 24, qfloor))
    check("N=1e%d: sup en q >= 10^3 x el ideal" % (12 if Nbig == 10 ** 12 else 24), qfloor > 1000,
          round(qfloor))

print("== RESULTADO: %s ==" % ("PASS" if not FAIL else "FAIL: " + ",".join(FAIL)))
