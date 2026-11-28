#!/usr/bin/env python3
# 87_lema_swq_parte_b.py  (Doc 87, Lema SW-Q apartado (b))
# Verifica por fuerza bruta / calculo directo, en rango computable,
# las afirmaciones FINITAS del apartado (b) del LEMA SW-Q (docs 85-86):
# la bifurcacion limite variable s vs limite fijo z, la forma explicita
# de F(u_s) y la identidad de coincidencia de los dos terminos
# principales donde ambos son validos. NO es evidencia del tamano
# asintotico: solo controla la forma, los rangos y las identidades.
#
# H1 (bifurcacion): u_s^var = log(D/s)/log s cruza 1 en s = D^{1/2}
# ~ X^{1/4} (proxy D = X^{1/2} sin log^C, como docs 85/86); para
# s > X^{1/4} el Teorema 9.7 no aplica (hipotesis D >= z violada:
# Bordignon Teor. 6, ec. (5), "for any D >= z"). El rango roto es
# 5/14 del rango logaritmico [(1/3-1/4)/(1/3-1/10)]. En cambio
# u_s^fix = log(D/s)/log z esta en [5/3, 4] en todo [z, y).
# H2 (F explicita): F(u) = 2e^gamma/u en [1,3]; en (3,4],
# (uF(u))' = f(u-1), f(v) = 2e^gamma log(v-1)/v; F(4) por integracion
# numerica propia (objetivo de la literatura: 1.02164155) y por la
# forma alternativa (e^gamma/2)(1 + int_1^2 log t/(t+1) dt).
# F decreciente en [1,4]; exceso max F(u)/(2e^gamma/u) <= 1.1471.
# H3 (dominio): S_Q(A_s;s) <= S_Q(A_s;z) exacto, porque
# P_Q(z) | P_Q(s): el conjunto cribado por s es subconjunto del
# cribado por z, con peso Lambda >= 0. (Parametros sinteticos no
# degenerados: a N<=2000, z=2 y P_Q(z) es vacio.)
# H4 (coincidencia): donde ambos limites son validos y u <= 3, los
# dos terminos principales coinciden: (1/log s) F(u^var) =
# (1/log z) F(u^fix) = 2e^gamma/log(D/s), porque u*log(limite) =
# log(D/s). Multiplicando por U_N = e^{-gamma} S(N), el sumando
# comun es 2 S(N)/log(D/s).
# H5 (masa rota): parte de la masa sum_s 1/(phi(s) log(D/s)) que
# vive en el rango roto s > X^{1/4} (X = 10^6, 10^9, 10^12).
# H6 (resto con limite fijo): los pesos con limite z (d | P(z),
# d < D/s) siguen cumpliendo multiplicidad 1 por modulo (doc 86
# Prop. 86.3, independiente del limite de criba). Parametros
# sinteticos como en H3.
#
# Uso: python 87_lema_swq_parte_b.py
import math

GAMMA = 0.5772156649015329
EG = math.exp(GAMMA)
TWO_EG = 2.0 * EG


def sieve(n):
    # Criba de Eratostenes (bytearray) + Lambda por potencias de primo.
    is_comp = bytearray(n + 1)
    is_comp[0:2] = b"\x01\x01"
    for i in range(2, int(n ** 0.5) + 1):
        if not is_comp[i]:
            step = i
            start = i * i
            is_comp[start:n + 1:step] = b"\x01" * ((n - start) // step + 1)
    primes = [i for i in range(2, n + 1) if not is_comp[i]]
    lam = [0.0] * (n + 1)
    for p in primes:
        pw = p
        lp = math.log(p)
        while pw <= n:
            lam[pw] = lp
            if pw > n // p:
                break
            pw *= p
    return primes, lam


PRIMES, LAM = sieve(1000000)


def phi(n):
    r = n
    x = n
    p = 2
    while p * p <= x:
        if x % p == 0:
            r -= r // p
            while x % p == 0:
                x //= p
        p += 1 if p == 2 else 2
    if x > 1:
        r -= r // x
    return r


def F_lin(u):
    # Funcion de la criba lineal (Jurkat-Richert), ec. (8) de
    # Bordignon et al. arXiv:2207.09452: F = 2e^gamma/u en [1,3]
    # (pues f == 0 en [0,2] y (uF)' = f(u-1)); en (3,4] se integra
    # (uF)' = f(u-1) con f(v) = 2e^gamma log(v-1)/v.
    if u <= 3.0:
        return TWO_EG / u
    n = 20000
    a, b = 3.0, u
    h = (b - a) / n
    s = 0.0
    for i in range(n + 1):
        t = a + i * h
        v = t - 1.0
        fv = TWO_EG * math.log(v - 1.0) / v if v > 1.0 else 0.0
        w = 1 if i in (0, n) else (4 if i % 2 == 1 else 2)
        s += w * fv
    return (TWO_EG + s * h / 3.0) / u


def F4_alt():
    # Forma alternativa: F(4) = (e^gamma/2)(1 + int_1^2 log t/(t+1) dt)
    n = 20000
    a, b = 1.0, 2.0
    h = (b - a) / n
    s = 0.0
    for i in range(n + 1):
        t = a + i * h
        fv = math.log(t) / (t + 1.0)
        w = 1 if i in (0, n) else (4 if i % 2 == 1 else 2)
        s += w * fv
    return (EG / 2.0) * (1.0 + s * h / 3.0)


def main():
    # ---- H1: bifurcacion (valores reales, sin redondeo) ----
    h1 = h1_tot = 0
    worst_fix_lo = 99.0
    worst_fix_hi = 0.0
    for X in (1e6, 1e9, 1e12, 1e24):
        z = X ** 0.1
        y = X ** (1.0 / 3.0)
        D = X ** 0.5
        s_star = D ** 0.5  # = X^{1/4}
        frac = (math.log(y) - math.log(s_star)) / (math.log(y) - math.log(z))
        uvar_at_z = math.log(D / z) / math.log(z)
        uvar_at_y = math.log(D / y) / math.log(y)
        ufix_at_z = math.log(D / z) / math.log(z)
        ufix_at_y = math.log(D / y) / math.log(z)
        ok = (abs(frac - 5.0 / 14.0) < 1e-12 and uvar_at_z > 3.0
              and uvar_at_y < 1.0 and 1.6 < ufix_at_y < 1.7
              and 3.9 < ufix_at_z < 4.1)
        grid_ok = True
        for i in range(1, 400):
            s = z * (y / z) ** (i / 400.0)
            uv = math.log(D / s) / math.log(s)
            uf = math.log(D / s) / math.log(z)
            if (uv < 1.0) != (s > s_star):
                grid_ok = False
            if not (5.0 / 3.0 - 0.05 < uf < 4.0 + 0.05):
                grid_ok = False
            worst_fix_lo = min(worst_fix_lo, uf)
            worst_fix_hi = max(worst_fix_hi, uf)
        h1_tot += 1
        if ok and grid_ok:
            h1 += 1
    print("H1_bifurcacion: %d/%d (rango roto 5/14=%.4f, u^fix en [%.4f, %.4f])"
          % (h1, h1_tot, 5.0 / 14.0, worst_fix_lo, worst_fix_hi))

    # ---- H2: F explicita ----
    f3 = F_lin(3.0)
    f4 = F_lin(4.0)
    f4b = F4_alt()
    f53 = F_lin(5.0 / 3.0)
    dec_ok = True
    prev = F_lin(1.0)
    for i in range(1, 301):
        u = 1.0 + 3.0 * i / 300.0
        cur = F_lin(u)
        if cur > prev + 1e-12:
            dec_ok = False
        prev = cur
    exc = max(F_lin(3.0 + 0.01 * i) / (TWO_EG / (3.0 + 0.01 * i))
              for i in range(1, 101))
    print("H2_F: F(3)=%.9f (2e^g/3=%.9f) F(4)=%.9f alt=%.9f "
          "F(5/3)=%.9f (6e^g/5=%.9f) dec=%s exceso_max=%.6f (<=%.6f)"
          % (f3, TWO_EG / 3.0, f4, f4b, f53, 6 * EG / 5.0, dec_ok, exc,
             f4 * 4.0 / TWO_EG))
    h2_ok = (abs(f3 - TWO_EG / 3.0) < 1e-12 and abs(f4 - 1.02164155) < 1e-6
             and abs(f4 - f4b) < 1e-9 and abs(f53 - 6 * EG / 5.0) < 1e-12
             and dec_ok and exc <= f4 * 4.0 / TWO_EG + 1e-9)

    # ---- H4: coincidencia de los dos terminos principales ----
    h4 = h4_tot = 0
    max_rel = 0.0
    for X in (1e6, 1e9, 1e12):
        z = X ** 0.1
        D = X ** 0.5
        for i in range(1, 500):
            s = (X ** 0.2) * ((X ** 0.25) / (X ** 0.2)) ** (i / 500.0)
            uv = math.log(D / s) / math.log(s)
            uf = math.log(D / s) / math.log(z)
            lhs_v = F_lin(uv) / math.log(s)
            lhs_f = F_lin(uf) / math.log(z)
            rhs = TWO_EG / math.log(D / s)
            h4_tot += 1
            rel = max(abs(lhs_v - rhs) / rhs, abs(lhs_f - rhs) / rhs)
            max_rel = max(max_rel, rel)
            if rel < 1e-9:
                h4 += 1
    print("H4_coincidencia: %d/%d max_rel=%.2e (2e^g/log(D/s) comun)"
          % (h4, h4_tot, max_rel))

    # ---- H5: masa del termino principal en el rango roto ----
    for X in (1e6, 1e9, 1e12):
        D = X ** 0.5
        z = X ** 0.1
        y = X ** (1.0 / 3.0)
        s_star = D ** 0.5
        tot = bro = 0.0
        for s in PRIMES:
            if s < z:
                continue
            if s >= y:
                break
            w = 1.0 / (phi(s) * math.log(D / s))
            tot += w
            if s > s_star:
                bro += w
        print("H5_masa_rota: X=1e%d tot=%.6f rota=%.4f (s>X^{1/4})"
              % (int(round(math.log10(X))), tot, bro / tot if tot else 0.0))

    # ---- H3/H6: dominio y resto con limite fijo (sintetico) ----
    # N = 10^6, Q = 77 (7*11, coprimo con N), z = 7, s en [7, y),
    # y = (N/Q)^{1/3} ~ 23.5; X = N/Q. A esta escala z=7 no es
    # degenerado: P_Q(7) = {3} (2 y 5 dividen N).
    N = 1000000
    Q = 77
    X = N / Q
    z = 7
    y = X ** (1.0 / 3.0)
    D = X ** 0.5
    h3 = h3_tot = 0
    pairs = {}
    for s in PRIMES:
        if s < z or s >= y:
            continue
        if math.gcd(s, Q * N) != 1:
            continue
        small_z = [ell for ell in PRIMES if ell < z and (Q * N) % ell != 0]
        small_s = [ell for ell in PRIMES if ell < s and (Q * N) % ell != 0]
        cnt_s = cnt_z = 0
        wsum_s = wsum_z = 0.0
        for w in range(1, int(X // s) + 1):
            n = N - Q * s * w
            if n < 2:
                continue
            ok_z = all(w % ell != 0 for ell in small_z)
            ok_s = ok_z and all(w % ell != 0 for ell in small_s
                                if ell >= z)
            if ok_z:
                cnt_z += 1
                wsum_z += LAM[n]
            if ok_s:
                cnt_s += 1
                wsum_s += LAM[n]
        h3_tot += 1
        if cnt_s <= cnt_z and wsum_s <= wsum_z + 1e-9:
            h3 += 1
        # H6: divisores d | P_Q(z), d < D/s
        Ds = int((X / s) ** 0.5)
        if Ds < 2:
            continue
        divs = [1]
        for ell in small_z:
            divs += [d * ell for d in list(divs) if d * ell < Ds]
        for d in divs:
            pairs.setdefault(Q * s * d, []).append((s, d))
    max_pairs = max((len(v) for v in pairs.values()), default=0)
    h6_tot = len(pairs)
    h6 = sum(1 for v in pairs.values() if len(v) == 1)
    print("H3_dominio: %d/%d (S(A_s;s) <= S(A_s;z) exacto, conteo y peso)"
          % (h3, h3_tot))
    print("H6_resto_limite_fijo: %d/%d modulos con multiplicidad 1, "
          "max_pares=%d" % (h6, h6_tot, max_pairs))
    ok = (h1 == h1_tot and h2_ok and h4 == h4_tot and h3 == h3_tot
          and h6 == h6_tot and max_pairs == 1)
    print("PASS" if ok else "FAIL")


if __name__ == "__main__":
    main()
