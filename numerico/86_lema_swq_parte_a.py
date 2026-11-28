#!/usr/bin/env python3
# 86_lema_swq_parte_a.py  (Doc 86, Lema SW-Q apartado (a))
# Verifica por fuerza bruta, en rango computable, las identidades FINITAS
# del apartado (a) del LEMA SW-Q (doc 85, sec. 5) leidas en Nathanson
# cap. 9-10 peso a peso. NO es evidencia de la estimacion asintotica:
# solo controla la forma de los pesos, la multiplicidad de los modulos
# y la identidad exacta del resto.
#
# H1 (Teorema 9.3, beta=2): los pesos lambda^+(d) = mu(d) si d en V+,
# 0 si no; lambda^-(d) analogo con V-. Se comprueba: soporte d<D,
# d libre de cuadrados, d|P(D); |lambda|<=1; lambda(1)=1; y las
# desigualdades de criba (9.8): sum_{d|m} lambda^-(d) <= 0 <=
# sum_{d|m} lambda^+(d) para todo m.
# H2 (Prop. 86.2, coeficiente afilado): cada modulo m=Qsd aparece a lo
# sumo UNA vez en sum_s sum_d lambda_d(s) E(N;Qsd,N): el s es el mayor
# factor primo de m/Q en [z,y). Se comprueba maximo de pares (s,d)
# por modulo == 1 (y <= 5, cota del doc 85).
# H3 (resto exacto, correccion de (85.6)): r_s(d) = |A_{s,d}| -
# |A_s|/phi(d) es EXACTAMENTE E(N;Qsd,N) - E(N;Qs,N)/phi(d) -
# Lambda(N)*(1-1/phi(d)), con |A_{s,d}| = psi(N;Qsd,N)-Lambda(N)
# (84.2, version PLANA sin pre-criba (w,Q)=1, como en el doc 85).
# Se miden tambien los dos terminos de correccion.
# H4 (densidad (9.34), ilustracion finita): R(u,z) =
# [prod_{u<=p<z, p no | QN} (1-1/(p-1))^{-1}] * (log u / log z);
# se reporta el maximo (a escala finita no tiende a 1; es solo
# control de forma, no prueba de (9.34)).
# H5 (pre-criba (w,Q)=1, cierra doc 85 sec. 3(d)): el conteo
# pre-cribado |A^Q_{s,d}| = sum_{e|Q} mu(e) [psi(N;Qsde,N)-Lambda(N)]
# es exacto por Mobius sobre e|Q.
#
# Escala: N par en [100,2000], k=2,3; H2/H3/H5 en submuestra N%10==0.
# D_s = floor((X/s)^{1/2}) (proxy finito sin log^C, como doc 85).
# Uso: python 86_lema_swq_parte_a.py
import math

MAXN = 2000
NMIN = 100
KS = (2, 3)


def sieve(n):
    spf = list(range(n + 1))
    for i in range(2, int(n ** 0.5) + 1):
        if spf[i] == i:
            for j in range(i * i, n + 1, i):
                if spf[j] == j:
                    spf[j] = i
    is_prime = [False] * (n + 1)
    lam = [0.0] * (n + 1)
    mu = [0] * (n + 1)
    mu[1] = 1
    for i in range(2, n + 1):
        if spf[i] == i:
            is_prime[i] = True
            lam[i] = math.log(i)
            mu[i] = -1
    for i in range(2, n + 1):
        p = spf[i]
        j = i // p
        if j % p == 0:
            lam[i] = math.log(p)
            mu[i] = 0
        else:
            lam[i] = 0.0
            mu[i] = -mu[j]
    primes = [i for i in range(2, n + 1) if is_prime[i]]
    return primes, is_prime, lam, mu


PRIMES, IS_PRIME, LAM, MU = sieve(MAXN)


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


def selector_Q(N, k):
    out = 1
    cnt = 0
    for p in PRIMES:
        if p == 2:
            continue
        if N % p == 0:
            continue
        out *= p
        cnt += 1
        if cnt == k - 1:
            break
    return out


def psi_mod(N, m, a):
    start = a % m
    if start == 0:
        start = m
    s = 0.0
    for n in range(start, N + 1, m):
        s += LAM[n]
    return s


def prime_factors_desc(d):
    # factores primos distintos de d, en orden decreciente
    out = []
    x = d
    p = 2
    while p * p <= x:
        if x % p == 0:
            out.append(p)
            while x % p == 0:
                x //= p
        p += 1 if p == 2 else 2
    if x > 1:
        out.append(x)
    return sorted(out, reverse=True)


def in_V(d, D, parity):
    # Teorema 9.3 con beta=2: d = p1...pk, p1>...>pk, y para todo m
    # de paridad dada, p_m < sqrt(D/(p1...pm)). parity=1: m impar (V+);
    # parity=0: m par (V-). d=1 siempre esta.
    if d == 1:
        return True
    if MU[d] == 0:
        return False
    ps = prime_factors_desc(d)
    if ps[0] >= D:
        return False
    prod = 1
    for i, p in enumerate(ps):
        m = i + 1
        prod *= p
        if (m % 2 == 1) == (parity == 1):
            if not (p * p * prod < D):
                return False
    return True


def main():
    # ---- H1: construccion de Teorema 9.3 (beta=2) ----
    h1 = 0
    h1_tot = 0
    for (D, Pset) in [(10, {3, 5, 7}), (30, {3, 5, 7, 11, 13, 17, 19, 23, 29}),
                      (100, {p for p in PRIMES if 2 < p < 100})]:
        lam_p = {}
        lam_m = {}
        for d in range(1, D):
            if MU[d] == 0:
                continue
            ps = prime_factors_desc(d)
            if any(p not in Pset for p in ps):
                continue
            if in_V(d, D, 1):
                lam_p[d] = MU[d]
            if in_V(d, D, 0):
                lam_m[d] = MU[d]
        lam_p[1] = 1
        lam_m[1] = 1
        ok = True
        # soporte, acotacion, lambda(1)=1
        for d, v in list(lam_p.items()) + list(lam_m.items()):
            if abs(v) > 1 or d >= D or MU[d] == 0:
                ok = False
        if lam_p.get(1) != 1 or lam_m.get(1) != 1:
            ok = False
        # desigualdades de criba (9.8) para todo m | P(D), m <= 200:
        # la prueba de Teorema 9.3 supone m divisor de P(D); para un
        # primo p fuera de P, sum_{d|p} lambda^-(d) = 1, y (9.8) no se
        # afirma. Se restringe a m libres de cuadrados con primos en P.
        for m in range(2, 201):
            if MU[m] == 0:
                continue
            if any(p not in Pset for p in prime_factors_desc(m)):
                continue
            sp = sum(lam_p.get(d, 0) for d in range(1, m + 1) if m % d == 0)
            sm = sum(lam_m.get(d, 0) for d in range(1, m + 1) if m % d == 0)
            if not (sm <= 0 <= sp):
                ok = False
                break
        h1_tot += 1
        if ok:
            h1 += 1
    print("H1_Th93: %d/%d (pesos mu o 0, |lam|<=1, (9.8) para m<=200)" % (h1, h1_tot))

    # ---- H4: densidad (9.34), ilustracion sintetica con z=100 ----
    # R(u,z) = [prod_{u<=p<z, p no | QN} (1-1/(p-1))^{-1}] * log u/log z.
    # La afirmacion asintotica es R -> 1; a z=100 se mide la forma.
    max_R = 0.0
    Q0 = selector_Q(2000, 2)
    for u in range(10, 100):
        prod = 1.0
        for p in PRIMES:
            if p < u or p >= 100:
                continue
            if (Q0 * 2000) % p == 0:
                continue
            prod *= 1.0 / (1.0 - 1.0 / (p - 1))
        R = prod * (math.log(u) / math.log(100))
        if R > max_R:
            max_R = R
    print("H4_densidad_934: max_R=%.3f (ilustracion finita, no prueba)" % max_R)

    # ---- H2/H3/H5: switching con Q ----
    h2 = h2_tot = 0
    max_pairs = 0
    worst_pairs = None
    h3 = h3_tot = 0
    max_corr2 = 0.0
    max_corr3 = 0.0
    max_E1 = 0.0
    h5 = h5_tot = 0
    for k in KS:
        for N in range(NMIN, MAXN + 1, 2):
            Q = selector_Q(N, k)
            X = N / Q
            z = max(2, int(X ** 0.1))
            y = max(z + 1, int(X ** (1.0 / 3.0)))
            logX = math.log(X)
            # H4: densidad (9.34) finita. A escala N<=2000, z=2 y el
            # rango es vacio; se ilustra con z=100 sintetico (primes
            # p no | QN, N=2000, k=2) para controlar la forma del
            # cociente de Mertens, no la afirmacion asintotica.
            pass
            if N % 10 != 0 and N not in (128, 256, 512, 1024):
                continue
            # pares (s,d) por modulo, para H2
            pairs = {}
            for s in PRIMES:
                if s < z or s >= y:
                    continue
                if math.gcd(s, Q * N) != 1:
                    continue
                Ds = int((X / s) ** 0.5)
                if Ds < 2:
                    continue
                # d | P_Q(s), d < Ds: divisores libres de cuadrados con
                # primos < s, primos no | QN
                small = [ell for ell in PRIMES if ell < s and (Q * N) % ell != 0]
                divs = [1]
                for ell in small:
                    divs += [d * ell for d in list(divs) if d * ell < Ds]
                for d in divs:
                    m = Q * s * d
                    pairs.setdefault(m, []).append((s, d))
            for m, lst in pairs.items():
                h2_tot += 1
                if len(lst) <= 5:
                    h2 += 1
                if len(lst) > max_pairs:
                    max_pairs = len(lst)
                    worst_pairs = (N, k, m, lst)
            # H3/H5: resto exacto y pre-criba (submuestra s,d)
            for s in PRIMES:
                if s < z or s >= y:
                    continue
                if math.gcd(s, Q * N) != 1:
                    continue
                Ds = int((X / s) ** 0.5)
                if Ds < 2:
                    continue
                small = [ell for ell in PRIMES if ell < s and (Q * N) % ell != 0]
                divs = [1]
                for ell in small:
                    divs += [d * ell for d in list(divs) if d * ell < Ds]
                # |A_s| y psi base (version plana, sin (w,Q)=1)
                As = 0.0
                for w in range(1, int(X // s) + 1):
                    n = N - Q * s * w
                    if n >= 2:
                        As += LAM[n]
                psi_Qs = psi_mod(N, Q * s, N)
                E_Qs = psi_Qs - N / phi(Q * s)
                lamN = LAM[N]
                for d in divs:
                    h3_tot += 1
                    # |A_{s,d}| por fuerza bruta
                    Asd = 0.0
                    for w in range(d, int(X // s) + 1, d):
                        n = N - Q * s * w
                        if n >= 2:
                            Asd += LAM[n]
                    r = Asd - As / phi(d)
                    E1 = psi_mod(N, Q * s * d, N) - N / phi(Q * s * d)
                    rhs = E1 - E_Qs / phi(d) - lamN * (1.0 - 1.0 / phi(d))
                    if abs(r - rhs) <= 1e-9 * max(1.0, abs(r), abs(rhs)):
                        h3 += 1
                    c2 = abs(E_Qs) / phi(d)
                    c3 = abs(lamN) * (1.0 - 1.0 / phi(d))
                    if c2 > max_corr2:
                        max_corr2 = c2
                    if c3 > max_corr3:
                        max_corr3 = c3
                    if abs(E1) > max_E1:
                        max_E1 = abs(E1)
                    # H5: pre-criba (w,Q)=1 por Mobius sobre e|Q
                    if math.gcd(d, Q) == 1:
                        h5_tot += 1
                        AsdQ = 0.0
                        for w in range(d, int(X // s) + 1, d):
                            if math.gcd(w, Q) != 1:
                                continue
                            n = N - Q * s * w
                            if n >= 2:
                                AsdQ += LAM[n]
                        mob = 0.0
                        e = 1
                        # divisores de Q (Q libre de cuadrados)
                        qprimes = prime_factors_desc(Q)
                        divsQ = [1]
                        for pp in qprimes:
                            divsQ += [t * pp for t in list(divsQ)]
                        for e in divsQ:
                            mob += MU[e] * (psi_mod(N, Q * s * d * e, N) - lamN)
                        if abs(AsdQ - mob) <= 1e-9 * max(1.0, abs(AsdQ), abs(mob)):
                            h5 += 1
    print("H2_coef_afilado: %d/%d modulos, max_pares=%d (peor %s)" %
          (h2, h2_tot, max_pairs, worst_pairs))
    print("H3_resto_exacto: %d/%d, max|E1|=%.3f max_corr2=%.3f max_corr3=%.3f" %
          (h3, h3_tot, max_E1, max_corr2, max_corr3))
    print("H5_precriba_Q: %d/%d (Mobius sobre e|Q exacto)" % (h5, h5_tot))
    ok = (h1 == h1_tot and h2 == h2_tot and max_pairs == 1 and
          h3 == h3_tot and h5 == h5_tot)
    print("PASS" if ok else "FAIL")


if __name__ == "__main__":
    main()
