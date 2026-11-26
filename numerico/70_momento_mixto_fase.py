#!/usr/bin/env python3
"""Doc 70, etapa 2: lema del momento mixto, via CONJUNTA con fase.
La identidad exacta U(N) = R(N) - D(N) (rho = 1_p - D) da, en arcos
menores, E_U = E_R - E_D: el error de arcos menores de rho ES el error
de arcos menores de Goldbach menos el termino cruzado del defecto.
Mide: la identidad exacta por criba, la cota de Cauchy de |E_D|,
la identidad de fase en una malla FFT, y la escala de la convolucion
para la gran criba. Uso: python3 numerico/70_momento_mixto_fase.py
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
print("== 0. Criba X=%d (docs 67/69) ==" % X)
is_prime = bytearray(b"\x01") * (X + 1)
is_prime[0:2] = b"\x00\x00"
lim = int(X ** 0.5) + 1
for p in range(2, lim):
    if is_prime[p]:
        is_prime[p * p:X + 1:p] = b"\x00" * ((X - p * p) // p + 1)
spf = np.arange(X + 1, dtype=np.int64)
spf[0] = 0
spf[1] = 10 ** 9
for p in range(2, lim):
    if is_prime[p]:
        spf[p * p:X + 1:p] = np.minimum(spf[p * p:X + 1:p], p)
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
sd2 = int(np.sum(D ** 2, dtype=np.int64))
print("sum D^2=%d (%.3f pi)" % (sd2, sd2 / piX))
check("sum D^2 = 0.467 pi (doc 69)", abs(sd2 / piX - 0.467) < 0.01, round(sd2 / piX, 3))

print("== 1. Identidad exacta U(N) = R(N) - D(N) (Prop. 70.1) ==")
# R(N) = #{p<N: N-p primo}; U(N) = sum_{p<N} rho~(N-p); D(N) = sum_{p<N} D(N-p).
# rho~ = 1_p - D termino a termino => U = R - D exacto (enteros).
is_np = np.array(is_prime, dtype=np.int32)
for N in [200000, 500000, 1000000, 1500000, 2000000]:
    R = 0
    U = 0
    Dn = 0
    for p in primes:
        if p >= N:
            break
        m = N - p
        if is_prime[m]:
            R += 1
        U += int(is_prime[m]) - int(D[m])
        Dn += int(D[m])
    print("N=%d: R=%d U=%d D(N)=%d ; U+D-R=%d ; D/R=%.3f" % (N, R, U, Dn, U + Dn - R, Dn / R))
    check("N=%d: U = R - D exacto" % N, U + Dn - R == 0, "U=%d R=%d D=%d" % (U, R, Dn))
    check("N=%d: D/R en (0.15, 0.6) (defecto de tamano comparable)" % N,
          0.15 < Dn / R < 0.6, round(Dn / R, 3))

print("== 2. Cota de Cauchy de |E_D| (el termino cruzado) ==")
# |E_D| = |int_m S1 S_D e(-N a)| <= ||S1||_2 ||S_D||_2 = sqrt(pi * sum D^2).
N = X
lN = math.log(N)
target = N / lN ** 2
cauchy_ED = math.sqrt(piX * sd2)
print("||S1||_2 ||S_D||_2 = %.0f ; target N/log^2 N = %.0f ; cociente = %.1f"
      % (cauchy_ED, target, cauchy_ED / target))
check("Cauchy de |E_D| >> target (>5x): el cruzado no es pequeno por Cauchy",
      cauchy_ED / target > 5, round(cauchy_ED / target, 1))
for NN in [10 ** 12, 10 ** 24]:
    l = math.log(NN)
    print("N=1e%d: cociente asintotico 0.683*log N = %.1f -> infinito"
          % (round(math.log10(NN)), 0.683 * l))
    check("N=1e%d: cociente Cauchy/target crece con log N" % round(math.log10(NN)),
          0.683 * l > 10, round(0.683 * l, 1))

print("== 3. Identidad de fase en la malla FFT (E_U = E_R - E_D) ==")
# Malla M = 2^21 puntos a_j = j/M. S(a) = sum a_n e(n a), e(x)=exp(2 pi i x).
# np.fft.fft usa signo -; S(a_j) = conj(fft(a))[j]. Arcos mayores: q <= Q0,
# |a - a/q| <= 1/(q Q1). La suma de Riemann sobre la malla verifica la
# identidad algebraica E_U = E_R - E_D exactamente (linealidad); las
# magnitudes son sumas de malla, no las integrales verdaderas (vease doc).
M = 2 ** 21
Q0 = int(LX ** 2)
Q1 = X / LX ** 2
print("M=%d ; Q0=%d ; Q1=%.0f" % (M, Q0, Q1))
a1 = np.zeros(M, dtype=np.float64)
a1[:X] = is_np[:X]
aD = np.zeros(M, dtype=np.float64)
aD[:X] = D[:X]
S1 = np.conj(np.fft.fft(a1))
SD = np.conj(np.fft.fft(aD))
Srho = S1 - SD
major = np.zeros(M, dtype=bool)
for q in range(1, Q0 + 1):
    hw = M / (q * Q1)
    for a in range(q):
        if a == 0 and q > 1:
            continue
        if math.gcd(a, q) != 1 and not (a == 0 and q == 1):
            continue
        c = M * a / q
        lo = max(0, int(math.ceil(c - hw)))
        hi = min(M - 1, int(math.floor(c + hw)))
        if lo <= hi:
            major[lo:hi + 1] = True
minor = ~major
print("fraccion de malla en arcos mayores: %.4f" % (major.mean()))
check("arcos mayores son una fraccion pequena (<10%)", major.mean() < 0.10,
      round(float(major.mean()), 4))
# Verificacion de circulo completo: (1/M) sum_j S1 Srho e(-N a_j) = U(N) exacto.
jj = np.arange(M)
ph = np.exp(-2j * np.pi * N * jj / M)
full = (S1 * Srho * ph).sum() / M
# U(N) a N=X por criba (seccion 1, N=2000000): recomputar rapido.
UN = 0
for p in primes:
    if p >= N:
        break
    m = N - p
    UN += int(is_prime[m]) - int(D[m])
print("circulo completo FFT: %.1f + %.1fi ; U(N) por criba: %d ; error abs: %.3f"
      % (full.real, full.imag, UN, abs(full - UN)))
check("FFT circulo completo = U(N) exacto (tolerancia 1e-3)", abs(full - UN) < 1e-3,
      round(abs(full - UN), 6))
ER = (S1[minor] ** 2 * ph[minor]).sum() / M
ED = (S1[minor] * SD[minor] * ph[minor]).sum() / M
EU = (S1[minor] * Srho[minor] * ph[minor]).sum() / M
print("|E_R|_malla=%.1f ; |E_D|_malla=%.1f ; |E_U|_malla=%.1f ; target=%.0f"
      % (abs(ER), abs(ED), abs(EU), target))
print("identidad malla: |E_U - (E_R - E_D)| = %.3g" % abs(EU - (ER - ED)))
check("identidad de fase en la malla: E_U = E_R - E_D (tol 1e-6 relativa)",
      abs(EU - (ER - ED)) / (abs(ER) + abs(ED) + 1.0) < 1e-6,
      "%.2g" % abs(EU - (ER - ED)))

print("== 4. Escala de la convolucion (gran criba del producto) ==")
# c_n = sum_p rho~(n-p) = U(n) (n par). Gran criba sobre c: sum_{q<=Q} sum_a
# |c(a/q)|^2 << (N+Q^2) sum |c_n|^2. La masa L^2 de c fija la escala.
Mc = 2 ** 22
F1c = np.fft.fft(a1[:X].astype(np.float64), Mc)
Frhoc = np.fft.fft((is_np[:X] - D[:X]).astype(np.float64), Mc)
Uconv = np.fft.ifft(F1c * Frhoc).real
Uconv = np.round(Uconv[:2 * X + 1]).astype(np.int64)
# Verificacion cruzada con la seccion 1 en los 5 N.
ok_conv = True
for N in [200000, 500000, 1000000, 1500000, 2000000]:
    Un = 0
    for p in primes:
        if p >= N:
            break
        Un += int(is_prime[N - p]) - int(D[N - p])
    if int(Uconv[N]) != Un:
        ok_conv = False
    print("N=%d: U por convolucion FFT=%d ; U por criba=%d" % (N, int(Uconv[N]), Un))
check("convolucion FFT = U(N) por criba en los 5 N", ok_conv)
sumU2 = int(np.sum(Uconv.astype(np.int64) ** 2, dtype=np.int64))
scale3 = X ** 3 / LX ** 4
scale2 = X ** 2 / LX ** 4
print("sum U(n)^2 = %.3g ; X^3/log^4 X = %.3g ; cociente = %.3f"
      % (sumU2, scale3, sumU2 / scale3))
print("sum U(n)^2 / (X^2/log^4 X) = %.3g  (>>1: la masa es de orden X^3, no X^2)"
      % (sumU2 / scale2))
check("sum U^2 es de orden X^3/log^4 X (cociente O(1), entre 0.01 y 10)",
      0.01 < sumU2 / scale3 < 10, round(sumU2 / scale3, 3))
check("sum U^2 / (X^2/log^4 X) >> 1000 (la gran criba no baja la escala)",
      sumU2 / scale2 > 1000, "%.2g" % (sumU2 / scale2))

print("== 5. Aritmetica GRH-ideal (Prop. 69.2(b) extendida) ==")
# Aun con sup ideal N^{1/2} (mejor caso concebible bajo GRH): sup x L^1 da
# N/(log N)^{1/2}; cuarto momento via sup^2 x L^2 da N^2 (log N)^{-1} vs
# o(N^2/log^4 N) necesario: sobra (log N)^3.
for NN in [10 ** 12, 10 ** 24]:
    l = math.log(NN)
    supL1 = l ** 1.5
    fourth = l ** 3
    print("N=1e%d: sup-ideal x L^1 / target = (log N)^{3/2} = %.0f ; "
          "cuarto momento ideal / necesario = (log N)^3 = %.0f"
          % (round(math.log10(NN)), supL1, fourth))
    check("N=1e%d: ni el sup ideal salva sup x L^1 (>100x)" % round(math.log10(NN)),
          supL1 > 100, round(supL1, 0))
    check("N=1e%d: ni el sup ideal salva el cuarto momento (>1000x)" % round(math.log10(NN)),
          fourth > 1000, round(fourth, 0))

print("== RESULTADO: %s ==" % ("PASS" if not FAIL else "FAIL: " + ",".join(FAIL)))
