#!/usr/bin/env python3
"""Doc 77: segundo momento sobre N (74.4).
Formula (74.4) con precision, comprueba si se sigue del segundo momento
clasico de R(N) (Montgomery-Vaughan / Lema 4.2 del survey Bhowmik-Grimmelt)
mas un segundo momento para D(N), y mide la escala de los segundos
momentos de los errores de R, U, D por FFT (convolucion) a X = 2e5, 5e5, 1e6.
Reutiliza la construccion de rho_tilde del doc 67 (Maynard III sec. 9).
Uso: python3 numerico/77_segundo_momento_N.py  -> PASS/FAIL + tablas.
"""
import math
import gc
import numpy as np

FAIL = []

def check(name, cond, detail=""):
    print(("PASS " if cond else "FAIL ") + name + (" | " + str(detail) if detail else ""))
    if not cond:
        FAIL.append(name)

C2 = 0.66016181584686957392
SMAX = 2_000_000

print("== 0. Criba hasta %d ==" % SMAX)
is_prime = bytearray(b"\x01") * (SMAX + 1)
is_prime[0:2] = b"\x00\x00"
lim = int(SMAX ** 0.5) + 1
spf = np.arange(SMAX + 1, dtype=np.int64)
spf[0] = 0
spf[1] = 10 ** 9
for p in range(2, lim):
    if is_prime[p]:
        is_prime[p * p:SMAX + 1:p] = b"\x00" * ((SMAX - p * p) // p + 1)
        spf[p * p:SMAX + 1:p] = np.minimum(spf[p * p:SMAX + 1:p], p)
is_np = np.array(is_prime, dtype=np.int16)
piSMAX = int(sum(is_prime))
check("pi(2e6)=148933", piSMAX == 148933, piSMAX)
primes_all = [p for p in range(2, SMAX + 1) if is_prime[p]]

def build_rho(S):
    """rho_tilde a escala S (doc 66/67): 1_p - S_{2,3} - S_{2,5}."""
    LX = math.log(S)
    z1 = S ** (1.0 / 7)
    z2 = S ** (3.0 / 7)
    small = [p for p in primes_all if z1 <= p < z2 and p <= S]
    S23 = np.zeros(S + 1, dtype=np.int16)
    S25 = np.zeros(S + 1, dtype=np.int16)
    Sa = Sb = 0
    for ii, p2 in enumerate(small):
        v = math.log(p2) / LX
        for p1 in small[ii + 1:]:
            u = math.log(p1) / LX
            pr = p1 * p2
            if pr > S:
                continue
            in_a = (4.0 / 21 <= u <= 2.0 / 7) and (2.0 / 7 - u / 2 <= v <= min(u, 3.0 / 7 - u))
            in_b = (2.0 / 7 <= u <= 3.0 / 7) and (4.0 / 7 - u <= v <= min(u, (1.0 - u) / 2))
            if not (in_a or in_b):
                continue
            idx = np.arange(pr, S + 1, pr)
            keep = spf[idx // pr] >= p2
            idxk = idx[keep]
            if in_a:
                S23[idxk] += 1
                Sa += int(np.sum(keep))
            else:
                S25[idxk] += 1
                Sb += int(np.sum(keep))
    D = S23.astype(np.int16) + S25.astype(np.int16)
    rho = is_np[:S + 1] - D
    return rho, D, Sa, Sb

def sing(N):
    s = 2.0 * C2
    n = N
    for q in primes_all:
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

def conv(a, b):
    n = 1
    while n < len(a) + len(b):
        n *= 2
    A = np.fft.rfft(a.astype(np.float64), n)
    B = np.fft.rfft(b.astype(np.float64), n)
    c = np.fft.irfft(A * B, n)
    return np.rint(c).astype(np.int64)

print("== 1. Segundos momentos de los errores (FFT) ==")
print("X       ratio_R   ratio_U   ratio_U(c_rho) ratio_D   bad_R bad_U bad_D  a_sec  LS")
prev = {}
for X in [200000, 500000, 1000000]:
    S = 2 * X
    rho, D, Sa, Sb = build_rho(S)
    piS = int(np.sum(is_np[:S + 1]))
    LS = (Sa + Sb) / piS
    # sumas L^2 (input de Parseval): sum rho^2 = pi + sum D^2
    sumD2 = int(np.sum(D.astype(np.int32) ** 2, dtype=np.int64))
    sum_rho2 = piS + sumD2
    R = conv(is_np[:S + 1], is_np[:S + 1])
    U = conv(is_np[:S + 1], rho)
    del rho, D
    gc.collect()
    Ns = np.arange(X, 2 * X + 1, 2)  # pares en [X, 2X]
    lN = np.log(Ns.astype(np.float64))
    SN = np.array([sing(int(n)) for n in Ns])
    main_R = SN * Ns / lN ** 2
    main_U = (1.0 - LS) * main_R
    main_D = LS * main_R
    main_U_crho = 0.465 * main_R  # con el c_rho asintotico (doc 66)
    eR = R[Ns].astype(np.float64) - main_R
    eU = U[Ns].astype(np.float64) - main_U
    eD = (R[Ns] - U[Ns]).astype(np.float64) - main_D
    # e_U = e_R - e_D exacto (D = R - U, main_D = main_R - main_U)
    check("X=%d: e_U = e_R - e_D (identidad)" % X,
          float(np.max(np.abs(eU - (eR - eD)))) < 1e-6,
          float(np.max(np.abs(eU - (eR - eD)))))
    scale = X ** 3 / (math.log(X) ** 4)
    sR = float(np.sum(eR ** 2))
    sU = float(np.sum(eU ** 2))
    sU_c = float(np.sum((U[Ns].astype(np.float64) - main_U_crho) ** 2))
    sD = float(np.sum(eD ** 2))
    rR, rU, rUc, rD = sR / scale, sU / scale, sU_c / scale, sD / scale
    # conteos de malos: |error| >= (1/2) * termino principal
    badR = int(np.sum(np.abs(eR) >= 0.5 * main_R))
    badU = int(np.sum(np.abs(eU) >= 0.5 * main_U))
    badD = int(np.sum(np.abs(eD) >= 0.5 * main_D))
    # regresion del error de R sobre el termino secundario S(N)*N/log^3 N
    sec = SN * Ns / lN ** 3
    a_sec = float(np.sum(eR * sec) / np.sum(sec ** 2))
    # ajuste de dos terminos (secundario + siguiente orden) y residuo:
    # e ~ a*S(N)N/log^3 N + b*S(N)N/log^4 N ; el residuo es el error
    # fluctuante, la parte que el lema de segundo momento debe controlar.
    sec2 = SN * Ns / lN ** 4
    def fit2(e, s1, s2):
        A11 = float(np.sum(s1 * s1)); A12 = float(np.sum(s1 * s2))
        A22 = float(np.sum(s2 * s2))
        b1 = float(np.sum(e * s1)); b2 = float(np.sum(e * s2))
        det = A11 * A22 - A12 * A12
        a = (b1 * A22 - b2 * A12) / det
        b = (A11 * b2 - A12 * b1) / det
        res = e - a * s1 - b * s2
        return a, b, float(np.sum(res ** 2))
    aR, bR, resR = fit2(eR, sec, sec2)
    aU, bU, resU = fit2(eU, sec, sec2)
    aD, bD, resD = fit2(eD, sec, sec2)
    meanS2 = float(np.mean(SN ** 2))
    pred_sec = (14.0 / 3.0) * meanS2 / (math.log(X) ** 2) * (aR / 2.0) ** 2
    print("%7d  %.4f    %.4f    %.4f         %.4f    %3d   %3d   %3d   %.3f  %.4f" %
          (X, rR, rU, rUc, rD, badR, badU, badD, a_sec, LS))
    print("        residuo tras 2 terminos: rR=%.4f rU=%.4f rD=%.4f ; (aR,bR)=(%.2f,%.2f) (aU,bU)=(%.2f,%.2f) ; mean S^2=%.2f pred_sec=%.4f" %
          (resR / scale, resU / scale, resD / scale, aR, bR, aU, bU, meanS2, pred_sec))
    n_even = len(Ns)
    check("X=%d: ratio_R < 0.25 (escala del termino secundario)" % X, rR < 0.25, round(rR, 4))
    check("X=%d: ratio_U < 0.2" % X, rU < 0.2, round(rU, 4))
    check("X=%d: ratio_D < 0.05" % X, rD < 0.05, round(rD, 4))
    check("X=%d: ratio_U(c_rho) < 1" % X, rUc < 1.0, round(rUc, 4))
    check("X=%d: bad_R = 0" % X, badR == 0, badR)
    check("X=%d: bad_U <= 0.02%% del intervalo" % X, badU <= 0.0002 * n_even, badU)
    check("X=%d: bad_D <= 0.02%% del intervalo" % X, badD <= 0.0002 * n_even, badD)
    check("X=%d: residuo_R < ratio_R/3" % X, resR / scale < rR / 3.0,
          (round(resR / scale, 4), round(rR, 4)))
    check("X=%d: residuo_U < ratio_U/3" % X, resU / scale < rU / 3.0,
          (round(resU / scale, 4), round(rU, 4)))
    check("X=%d: sU <= 2*sR + 2*sD (E_U = E_R - E_D)" % X,
          sU <= 2 * sR + 2 * sD + 1e-6, round(sU / (2 * sR + 2 * sD), 3))
    check("X=%d: sum rho^2 / pi(S) en (1, 2)" % X,
          1.0 < sum_rho2 / piS < 2.0, round(sum_rho2 / piS, 3))
    check("X=%d: sum D^2 / pi(S) en (0.3, 0.6)" % X,
          0.3 < sumD2 / piS < 0.6, round(sumD2 / piS, 3))
    if prev:
        check("X=%d: ratio_R decrece con X" % X, rR < prev['rR'],
              (round(prev['rR'], 4), round(rR, 4)))
        check("X=%d: ratio_U decrece con X" % X, rU < prev['rU'],
              (round(prev['rU'], 4), round(rU, 4)))
        check("X=%d: ratio_D decrece con X" % X, rD < prev['rD'],
              (round(prev['rD'], 4), round(rD, 4)))
        check("X=%d: ratio_U(c_rho) decrece con X" % X, rUc < prev['rUc'],
              (round(prev['rUc'], 4), round(rUc, 4)))
        check("X=%d: residuo_U decrece con X" % X, resU / scale < prev['resU'],
              (round(prev['resU'], 4), round(resU / scale, 4)))
    check("X=%d: a_sec en (2.0, 3.0)" % X, 2.0 < a_sec < 3.0, round(a_sec, 3))
    check("X=%d: aR en (1.5, 2.5) (coef. secundario ~2)" % X, 1.5 < aR < 2.5, round(aR, 3))
    prev = {'rR': rR, 'rU': rU, 'rD': rD, 'rUc': rUc, 'resU': resU / scale}
    del R, U
    gc.collect()

print("== 2. Escala teorica vs medida (X=1e6) ==")
# Cota de Parseval+Vinogradov con A=2: sum |E_U|^2 << X^3/(log X)^5.
# El segundo momento medido del error COMPLETO (E_U + error de arcos
# mayores O(N/(log N)^A)) debe estar muy por debajo de esa cota.
X = 1000000
lX = math.log(X)
bound_A2 = X ** 3 / lX ** 5
bound_A3 = X ** 3 / lX ** 7
print("X=1e6: X^3/log^4 X = %.3e ; X^3/log^5 X (A=2) = %.3e ; X^3/log^7 X (A=3) = %.3e"
      % (X ** 3 / lX ** 4, bound_A2, bound_A3))
print("El segundo momento medido (sec. 1, ratio_U * X^3/log^4 X) debe ser")
print("<< X^3/log^5 X: ratio_U medido ~ %.4f ; residuo_U ~ %.4f ; 1/log X = %.4f" %
      (prev['rU'], prev['resU'], 1.0 / lX))
check("residuo_U(1e6) < 1/log X (el error fluctuante esta muy por debajo)",
      prev['resU'] < 1.0 / lX, (round(prev['resU'], 4), round(1.0 / lX, 4)))

print("== RESULTADO: %s ==" % ("PASS" if not FAIL else "FAIL: " + ",".join(FAIL)))
