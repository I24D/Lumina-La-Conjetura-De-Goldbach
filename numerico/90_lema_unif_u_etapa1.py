#!/usr/bin/env python3
# 90 - Lema UNIF-U etapa 1: por que el resto por-u domina.
# Solo stdlib. Control finito, no evidencia asintotica.
# H0: reproduce constantes del doc 89 (f(5), F(4), c_{1/10}, J0, Bal0, c*).
# H1: cociente resto/termino principal R = N^c/(ln N)^(A'-1) en malla + cruce N*.
# H2: A'(N) necesario para que el error BV sea menor que el termino principal.
# H3: cancelacion de raiz cuadrada en promedio: cociente y umbral c < 1/5.
# H4: precio del C0 creciente: exponente efectivo de criba e_eff(c,N) -> 1/2-c.
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

# ---- F y f por el sistema de retardo, malla fina (trapecio) [doc 89] ----
H = 0.0005
NGRID = 12000  # hasta x = 6.0
xs = [i * H for i in range(NGRID + 1)]
F = [0.0] * (NGRID + 1)
f = [0.0] * (NGRID + 1)
for i, x in enumerate(xs):
    if 0 < x <= 3.0:
        F[i] = 2.0 * EG / x
cum = 0.0
for i in range(1, NGRID + 1):
    x = xs[i]
    if 2.0 < x <= 4.0:
        g0 = 2.0 * EG / (x - H - 1.0)
        g1 = 2.0 * EG / (x - 1.0)
        cum += 0.5 * (g0 + g1) * H
        f[i] = cum / x
cum = 0.0
for i in range(1, NGRID + 1):
    x = xs[i]
    if 3.0 < x <= 4.0:
        j0 = int(round((x - H - 1.0) / H))
        j1 = int(round((x - 1.0) / H))
        cum += 0.5 * (f[j0] + f[j1]) * H
        F[i] = (2.0 * EG + cum) / x
f4 = f[int(round(4.0 / H))]
cum = 4.0 * f4
for i in range(1, NGRID + 1):
    x = xs[i]
    if 4.0 < x <= 5.0:
        j0 = int(round((x - H - 1.0) / H))
        j1 = int(round((x - 1.0) / H))
        cum += 0.5 * (F[j0] + F[j1]) * H
        f[i] = cum / x
F4 = F[int(round(4.0 / H))]
cum = 0.0
for i in range(1, NGRID + 1):
    x = xs[i]
    if 4.0 < x <= 5.0:
        j0 = int(round((x - H - 1.0) / H))
        j1 = int(round((x - 1.0) / H))
        cum += 0.5 * (f[j0] + f[j1]) * H
        F[i] = (4.0 * F4 + cum) / x
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

def f_closed(s):
    return 2.0 * EG * math.log(s - 1.0) / s

def F_34(u):
    val = simpson(lambda t: 2.0 * EG * math.log(t - 2.0) / (t - 1.0), 3.0, u, 2000)
    return (2.0 * EG + val) / u

def f_45(s):
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

def c_estrella():
    lo, hi = 0.0, 2.0 / 17.0
    for _ in range(60):
        mid = 0.5 * (lo + hi)
        if Bal(mid) > 0:
            lo = mid
        else:
            hi = mid
    return 0.5 * (lo + hi)

# ---- H0: reproduce doc 89 ----
J0 = J_of_c(0.0)
cs = c_estrella()
print("H0_reproduce | f5_grid=%.10f f5_simp=%.10f (doc88: 0.9982417256)" % (f5_grid, f5_simp))
print("H0_reproduce | F4_grid=%.10f F4_simp=%.10f (doc87: 1.021641553)" % (F[int(round(4.0 / H))], F4_simp))
print("H0_reproduce | c110_simp=%.10f (doc88: 0.4909952010)" % c110_simp)
print("H0_reproduce | J0=%.10f (doc87/88: 1.5238762503) Bal0=%.6f (doc88: 0.061404) c*=%.6f (doc89: 0.039884)"
      % (J0, Bal(0.0), cs))

# ---- H1: cociente resto/termino principal ----
def log10_R(c, log10N, Aprima):
    lnN = log10N * math.log(10.0)
    return c * log10N - (Aprima - 1.0) * math.log10(lnN)

def cruce(c, Aprima):
    # (A'-1) ln L = c L, L = ln N; devuelve log10 N*
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
    return (0.5 * (lo + hi)) / math.log(10.0)

print("H1_cociente | c A' log10R(N=1e6) log10R(1e12) log10R(1e24) log10R(1e76) log10N*")
for c in [0.01, 0.02, 0.039884, 0.05, 0.10, 2.0 / 17.0]:
    for Aprima in [5, 10]:
        print("H1_cociente | %.6f %2d %9.2f %9.2f %9.2f %9.2f %8.1f"
              % (c, Aprima, log10_R(c, 6, Aprima), log10_R(c, 12, Aprima),
                 log10_R(c, 24, Aprima), log10_R(c, 76, Aprima), cruce(c, Aprima)))

# ---- H2: A'(N) necesario ----
def Aprima_nec(c, log10N):
    lnN = log10N * math.log(10.0)
    return 2.0 + c * lnN / math.log(lnN)

print("H2_Aprima_nec | c log10N=12 log10N=24 log10N=76 log10N=205 log10N=1605")
for c in [0.02, 0.05, 0.10, 2.0 / 17.0]:
    print("H2_Aprima_nec | %.6f %8.2f %8.2f %8.2f %8.2f %8.2f"
          % (c, Aprima_nec(c, 12), Aprima_nec(c, 24), Aprima_nec(c, 76),
             Aprima_nec(c, 205), Aprima_nec(c, 1605)))

# ---- H3: cancelacion de raiz cuadrada en promedio ----
def log10_Rsq(c, log10N):
    lnN = log10N * math.log(10.0)
    return (5.0 * c / 4.0 - 0.25) * log10N + 2.0 * math.log10(lnN)

def cruce_sq(c):
    # (5c/4-1/4) L + 2 ln L = 0, L = ln N; devuelve log10 N* (solo si c<1/5)
    def g(L):
        return (5.0 * c / 4.0 - 0.25) * L + 2.0 * math.log(L)
    lo, hi = 2.0, 1e6
    if g(hi) > 0:
        return None  # c >= 1/5: el cociente nunca baja de 1
    for _ in range(200):
        mid = 0.5 * (lo + hi)
        if g(mid) > 0:
            lo = mid
        else:
            hi = mid
    return (0.5 * (lo + hi)) / math.log(10.0)

print("H3_raiz_cuadrada | c log10Rsq(N=1e12) log10Rsq(N=1e100) log10N*_sq")
for c in [0.05, 0.10, 0.15, 0.19, 0.20, 0.21, 0.25]:
    cr = cruce_sq(c)
    print("H3_raiz_cuadrada | %.2f %9.2f %9.2f %s"
          % (c, log10_Rsq(c, 12), log10_Rsq(c, 100), ("%.1f" % cr) if cr is not None else "nunca"))

# ---- H4: precio del C0 creciente ----
def e_eff(c, log10N):
    lnN = log10N * math.log(10.0)
    return 0.5 - c - 10.0 * math.log(lnN) / lnN

def log10_N0(c):
    # (1/2-c) L = 10 ln L, L = ln N; devuelve log10 N0 (D_A = 1)
    def g(L):
        return (0.5 - c) * L - 10.0 * math.log(L)
    lo, hi = 10.0, 1e5
    for _ in range(200):
        mid = 0.5 * (lo + hi)
        if g(mid) < 0:
            lo = mid
        else:
            hi = mid
    return (0.5 * (lo + hi)) / math.log(10.0)

print("H4_C0_precio | c e_eff(1e12) e_eff(1e24) e_eff(1e100) e_eff(1e1000) limite(1/2-c) log10N0")
for c in [0.0, 0.02, 0.05, 0.10, 2.0 / 17.0]:
    print("H4_C0_precio | %.6f %9.4f %9.4f %9.4f %9.4f %9.4f %8.1f"
          % (c, e_eff(c, 12), e_eff(c, 24), e_eff(c, 100), e_eff(c, 1000), 0.5 - c, log10_N0(c)))
print("PASS")
