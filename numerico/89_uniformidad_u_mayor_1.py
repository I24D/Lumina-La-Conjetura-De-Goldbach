#!/usr/bin/env python3
# 89 - Uniformidad de la incidencia completa I_k con u>1: delimitacion.
# Solo stdlib. Control finito, no evidencia asintotica.
# H0: reproduce constantes de docs 87/88 (f(5), F(4), J(0), Bal(0), c_{1/10}).
# H1: geometria de la criba vs c (A = N^c): u_max, u_min, s_lower; umbral c <= 2/17.
# H2: balance Bal(c) y su cero c* (biseccion).
# H3: conteo y suma reciproca de u con rad(u)|Q (fuerza bruta) vs asintotica.
# H4: masa del termino principal sumada sobre u (U polilog) vs termino u=1.
# H5: cruce N*(c,A') donde el error BV alcanza al termino principal por u.
# H6: C0(c,N) necesario para el ajuste BV (crece con N si c>0; constante si polilog).
import math

EG = math.exp(0.5772156649015329)  # e^gamma
C110 = 0.4909952010              # doc 88 (88.2), se re-verifica en H0

def simpson(f, a, b, n):
    if n % 2:
        n += 1
    h = (b - a) / n
    s = f(a) + f(b)
    for i in range(1, n):
        s += (4 if i % 2 else 2) * f(a + i * h)
    return s * h / 3.0

# ---- F y f por el sistema de retardo, malla fina (trapecio) ----
H = 0.0005
NGRID = 12000  # hasta x = 6.0
xs = [i * H for i in range(NGRID + 1)]
F = [0.0] * (NGRID + 1)
f = [0.0] * (NGRID + 1)
for i, x in enumerate(xs):
    if 0 < x <= 3.0:
        F[i] = 2.0 * EG / x
# f en (2,4]: s f(s) = int_2^s F(t-1) dt, F(t-1) = 2 EG/(t-1)
cum = 0.0
for i in range(1, NGRID + 1):
    x = xs[i]
    if 2.0 < x <= 4.0:
        # trapecio de G(t)=2EG/(t-1) en [x-H, x]
        g0 = 2.0 * EG / (x - H - 1.0)
        g1 = 2.0 * EG / (x - 1.0)
        cum += 0.5 * (g0 + g1) * H
        f[i] = cum / x
# F en (3,4]: u F(u) = 2 EG + int_3^u f(t-1) dt
cum = 0.0
for i in range(1, NGRID + 1):
    x = xs[i]
    if 3.0 < x <= 4.0:
        # f(t-1) interpolada en la malla (t-1 en [2,3], f ya calculada ahi)
        j0 = int(round((x - H - 1.0) / H))
        j1 = int(round((x - 1.0) / H))
        cum += 0.5 * (f[j0] + f[j1]) * H
        F[i] = (2.0 * EG + cum) / x
# f en (4,5]: s f(s) = 4 f(4) + int_4^s F(t-1) dt
f4 = f[int(round(4.0 / H))]
cum = 4.0 * f4
for i in range(1, NGRID + 1):
    x = xs[i]
    if 4.0 < x <= 5.0:
        j0 = int(round((x - H - 1.0) / H))
        j1 = int(round((x - 1.0) / H))
        cum += 0.5 * (F[j0] + F[j1]) * H
        f[i] = cum / x
# F en (4,5]: u F(u) = 4 F(4) + int_4^u f(t-1) dt
F4 = F[int(round(4.0 / H))]
cum = 0.0
for i in range(1, NGRID + 1):
    x = xs[i]
    if 4.0 < x <= 5.0:
        j0 = int(round((x - H - 1.0) / H))
        j1 = int(round((x - 1.0) / H))
        cum += 0.5 * (f[j0] + f[j1]) * H
        F[i] = (4.0 * F4 + cum) / x
# f en (5,6]: s f(s) = 5 f(5) + int_5^s F(t-1) dt
f5_grid = f[int(round(5.0 / H))]
cum = 5.0 * f5_grid
for i in range(1, NGRID + 1):
    x = xs[i]
    if 5.0 < x <= 6.0:
        j0 = int(round((x - H - 1.0) / H))
        j1 = int(round((x - 1.0) / H))
        cum += 0.5 * (F[j0] + F[j1]) * H
        f[i] = cum / x

def interp(arr, x):
    if x <= 0:
        return float('nan')
    pos = x / H
    i = int(pos)
    if i >= NGRID:
        return arr[NGRID]
    t = pos - i
    return arr[i] * (1 - t) + arr[i + 1] * t

def F_eval(u):
    return interp(F, u)

def f_eval(s):
    return interp(f, s)

# ---- Chequeos independientes por Simpson (formas cerradas) ----
def f_closed(s):
    # s en [2,4]
    return 2.0 * EG * math.log(s - 1.0) / s

def F_34(u):
    # u en [3,4]: (1/u)[2EG + int_3^u f(t-1) dt], f(t-1)=2EG ln(t-2)/(t-1)
    val = simpson(lambda t: 2.0 * EG * math.log(t - 2.0) / (t - 1.0), 3.0, u, 2000)
    return (2.0 * EG + val) / u

def f_45(s):
    # s en [4,5]: (1/s)[4 f(4) + int_4^s F(t-1) dt], F en [3,4] por F_34
    val = simpson(lambda t: F_34(t - 1.0), 4.0, s, 2000)
    return (4.0 * f_closed(4.0) + val) / s

F4_simp = F_34(4.0)
f5_simp = f_45(5.0)
c110_simp = simpson(lambda b: math.log(2.0 - 3.0 * b) / (b * (1.0 - b)), 0.1, 1.0 / 3.0, 4000)

def J_of_c(c):
    def u_of_b(b):
        return 10.0 * ((0.5 - c) - b * (1.0 - c)) / (1.0 - c)
    return simpson(lambda b: F_eval(u_of_b(b)) / b, 0.1, 1.0 / 3.0, 2000)

def Bal(c, phiQ_over_Q=1.0):
    s_low = (5.0 - 10.0 * c) / (1.0 - c)
    Jc = J_of_c(c)
    tail = 0.5 * 0.4 * EG * C110 * phiQ_over_Q / (1.0 - c)
    return f_eval(s_low) - 0.5 * Jc - tail

def umax(c):
    return (4.0 - 9.0 * c) / (1.0 - c)

def umin(c):
    return (5.0 / 3.0 - (20.0 / 3.0) * c) / (1.0 - c)

def slower(c):
    return (5.0 - 10.0 * c) / (1.0 - c)

print("H0_reproduce | f5_grid=%.10f f5_simp=%.10f (doc88: 0.9982417256)" % (f5_grid, f5_simp))
print("H0_reproduce | F4_grid=%.10f F4_simp=%.10f (doc87: 1.021641553)" % (F[int(round(4.0 / H))], F4_simp))
print("H0_reproduce | c110_simp=%.10f (doc88: 0.4909952010)" % c110_simp)
J0 = J_of_c(0.0)
print("H0_reproduce | J0=%.10f (doc87/88: 1.5238762503) Bal0=%.6f (doc88: 0.061404)" % (J0, Bal(0.0)))

print("H1_geometria | c umax umin s_lower")
for c in [0.0, 0.02, 0.05, 0.08, 0.10, 2.0 / 17.0]:
    print("H1_geometria | %.4f %.6f %.6f %.6f" % (c, umax(c), umin(c), slower(c)))
print("H1_geometria | umbral c=2/17: umin=%.10f (debe ser 1), umax=%.10f (debe ser 10/3), s_lower=%.10f (debe ser 13/3)"
      % (umin(2.0 / 17.0), umax(2.0 / 17.0), slower(2.0 / 17.0)))

print("H2_balance | c Bal(c) f(s_lower) J(c) tail")
for c in [0.0, 0.02, 0.04, 0.05, 0.06, 0.08, 0.10, 2.0 / 17.0]:
    s_low = slower(c)
    Jc = J_of_c(c)
    tail = 0.5 * 0.4 * EG * C110 / (1.0 - c)
    print("H2_balance | %.4f %.6f %.6f %.6f %.6f" % (c, Bal(c), f_eval(s_low), Jc, tail))
# biseccion del cero de Bal en [0, 2/17]
lo, hi = 0.0, 2.0 / 17.0
if Bal(hi) < 0:
    for _ in range(60):
        mid = 0.5 * (lo + hi)
        if Bal(mid) > 0:
            lo = mid
        else:
            hi = mid
    print("H2_balance | c_estrella=%.6f Bal(c*)=%.2e" % (0.5 * (lo + hi), Bal(0.5 * (lo + hi))))
else:
    print("H2_balance | Bal(2/17)=%.6f > 0: el balance NO se anula antes del umbral geometrico" % Bal(hi))

# ---- H3: conteo y suma reciproca de u con rad(u)|Q ----
def gen_u(primes, U):
    # todos los u <= U con rad(u)|Q (Q = producto de `primes`), por exponentes
    out = []
    def gen(idx, cur):
        if idx == len(primes):
            out.append(cur)
            return
        p = primes[idx]
        v = cur
        while v <= U:
            gen(idx + 1, v)
            if v > U // p:
                break
            v *= p
    gen(0, 1)
    return sorted(set(out))

for Q, primes in [(3, [3]), (15, [3, 5]), (105, [3, 5, 7])]:
    k = len(primes) + 1
    print("H3_conteo_u | Q=%d k=%d U N(U) asintotico sum1/u prod l/(l-1) cola" % (Q, k))
    for U in [100, 1000, 10000, 100000, 1000000]:
        us = gen_u(primes, U)
        n = len(us)
        asym = (math.log(U) ** (k - 1)) / (math.factorial(k - 1) * math.prod(math.log(p) for p in primes))
        s = sum(1.0 / u for u in us)
        prod = math.prod(p / (p - 1.0) for p in primes)
        print("H3_conteo_u | Q=%d U=%d N=%d asym=%.2f sum1/u=%.6f limite=%.6f cola=%.6f"
              % (Q, U, n, asym, s, prod, prod - s))

# ---- H4: masa sumada sobre u polilog vs u=1 ----
lnN = math.log(1e12)
for Q, primes in [(3, [3]), (15, [3, 5])]:
    U = lnN ** 2
    us = gen_u(primes, int(U))
    s = sum(1.0 / u for u in us)
    prod = math.prod(p / (p - 1.0) for p in primes)
    print("H4_masa | Q=%d N=1e12 t=2 U=%.1f sum1/u=%.6f (u=1 aporta 1) limite=%.6f" % (Q, U, s, prod))

# ---- H5: cruce N*(c,A') con (A'-1) ln ln N = c ln N ----
def cruce(c, Aprima):
    # resolver (A'-1) ln L = c L con L = ln N, biseccion en L
    def g(L):
        return (Aprima - 1.0) * math.log(L) - c * L
    lo, hi = 10.0, 1e5
    if g(lo) < 0:
        return None
    for _ in range(200):
        mid = 0.5 * (lo + hi)
        if g(mid) > 0:
            lo = mid
        else:
            hi = mid
    L = 0.5 * (lo + hi)
    return L / math.log(10.0)  # log10 N*

print("H5_cruce | c A' log10(N*)")
for c in [0.02, 0.05, 0.10, 2.0 / 17.0]:
    for Aprima in [5, 10, 20]:
        print("H5_cruce | %.4f %d %.1f" % (c, Aprima, cruce(c, Aprima)))
# razon termino/error en el rango polilog a N=1e12: (log N)^(A'-1-t)
print("H5_cruce | polilog N=1e12 t=2 A'=10: main/error=(log N)^7=%.3e" % (lnN ** 7))

# ---- H6: C0 necesario para el ajuste BV ----
def C0_nec(c, N, B=10.0):
    return B + (c / 2.0) * (math.log(N) / math.log(math.log(N)))

print("H6_C0 | c N C0_nec(B=10)")
for c in [0.05, 2.0 / 17.0]:
    for N in [1e6, 1e12, 1e24]:
        print("H6_C0 | %.4f 1e%d %.2f" % (c, int(math.log10(N)), C0_nec(c, N)))
print("H6_C0 | polilog (k=2,t=2): C0 = B + (k-1+t)/2 = %.1f constante" % (10.0 + 1.5))
print("PASS")
