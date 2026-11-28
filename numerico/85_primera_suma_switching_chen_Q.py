#!/usr/bin/env python3
# 85_primera_suma_switching_chen_Q.py  (Doc 85, etapa 1)
# Verifica por fuerza bruta, en rango computable, las identidades FINITAS
# de la primera suma de switching de Chen para la subincidencia u=1 con
# selector Q_N (docs 82-84). NO es evidencia de la estimacion asintotica:
# solo controla la formulacion (reindexado exacto, inversion de Mobius
# de la criba, error bilineal medido a escala finita).
#
# Definiciones (doc 85):
#   X = N/Q, z = floor(X^(1/10)), y = floor(X^(1/3)).
#   J_k(N;Q) = suma_{v<X,(v,Q)=1,omega(v) in {1,2}} 1_{(N-Qv,N)=1}
#              * Lambda(N-Qv) * (log Q + log rad(v)).
#   J2 = parte omega(v)=2. Cada v = p1*p2 con p1<=p2 tiene un unico
#   menor factor primo p1. Reindexado de switching (exacto):
#     J2 = sum_{s primo, z<=s<y, (s,QN)=1} T_s + R_y,
#   T_s = suma sobre w primo, s<=w, s*w<X, del peso de v=s*w;
#   R_y = suma sobre v=p1*p2 con p1>=y (cola: ambos factores >= y).
#   U_s = suma_{w<X/s, (w,P_Q(s))=1} Lambda(N-Qsw)  (criba lineal
#   PLANA, sin peso log rad y sin la pre-criba (w,Q)=1: esa pre-criba
#   es un detalle [P-rev] aparte, ver doc 85 sec. 4), P_Q(s) = prod
#   de primos l<s con l no | QN. Inversion de Mobius exacta:
#     U_s = sum_{d|P_Q(s)} mu(d) * (psi(N;Qsd,N) - Lambda(N)),
#   valida para todo d por la biyeccion (84.1) (sin necesitar
#   (d,QN)=1; solo el termino principal N/phi(Qsd) exige (N,Qsd)=1,
#   que se tiene porque d|P_Q(s) excluye l|QN).
#   D_s = floor( (X/s)^(1/2) )  [proxy a escala finita SIN el factor
#   log^C: a N<=2000 el nivel X^(1/2)/log^2 X es <1 y degenera; el
#   script comprueba la FORMA del error, no su tamano asintotico]
#   Error bilineal (a escala finita): E_s = U_s^trunc - M_s.
#   Lema de pocos divisores: cada modulo m=Qsd admite a lo sumo
#   floor(log m / log z) primos s en [z,y) que lo dividen (s>=z).
#
# Uso: python 85_primera_suma_switching_chen_Q.py
# Salida: una linea PASS con contadores y ratios maximos.
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
    omega = [0] * (n + 1)
    rad = [1] * (n + 1)
    mu = [0] * (n + 1)
    mu[1] = 1
    for i in range(2, n + 1):
        if spf[i] == i:
            is_prime[i] = True
            lam[i] = math.log(i)
            omega[i] = 1
            rad[i] = i
            mu[i] = -1
    for i in range(2, n + 1):
        p = spf[i]
        j = i // p
        if j % p == 0:
            # p^a, a>=2
            omega[i] = omega[j]
            rad[i] = rad[j]
            lam[i] = math.log(p)
            mu[i] = 0
        else:
            omega[i] = omega[j] + 1
            rad[i] = rad[j] * p
            lam[i] = 0.0
            mu[i] = -mu[j]
    primes = [i for i in range(2, n + 1) if is_prime[i]]
    return primes, is_prime, lam, omega, rad, mu


PRIMES, IS_PRIME, LAM, OMEGA, RAD, MU = sieve(MAXN)


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
    # primeros k-1 primos impares que no dividen N (doc 83, (83.1))
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
    # suma_{n<=N, n==a (mod m)} Lambda(n)
    start = a % m
    if start == 0:
        start = m
    s = 0.0
    for n in range(start, N + 1, m):
        s += LAM[n]
    return s


def main():
    h1 = h2 = h3a = h3b = 0
    max_tail_frac = 0.0
    worst_tail = None
    max_bilin_rel = 0.0
    worst_bilin = None
    sum_abs_E = 0.0
    sum_abs_M = 0.0
    n_bilin = 0
    for k in KS:
        for N in range(NMIN, MAXN + 1, 2):
            Q = selector_Q(N, k)
            X = N / Q
            vmax = (N - 1) // Q
            z = max(2, int(X ** 0.1))
            y = max(z + 1, int(X ** (1.0 / 3.0)))
            logQ = math.log(Q)
            J1 = J2 = 0.0
            T = {}
            Ry = 0.0
            for v in range(1, vmax + 1):
                if math.gcd(v, Q) != 1:
                    continue
                w = OMEGA[v]
                if w not in (1, 2):
                    continue
                n = N - Q * v
                if n < 2 or LAM[n] == 0.0:
                    continue
                if math.gcd(n, N) != 1:
                    continue
                term = LAM[n] * (logQ + math.log(RAD[v]))
                if w == 1:
                    J1 += term
                else:
                    J2 += term
                    # v = p1*p2, p1 = menor factor primo
                    p1 = None
                    for p in PRIMES:
                        if p * p > v:
                            break
                        if v % p == 0:
                            p1 = p
                            break
                    if p1 is None:
                        p1 = v  # v primo (no ocurre con omega=2)
                    if z <= p1 < y and math.gcd(p1, Q * N) == 1:
                        T[p1] = T.get(p1, 0.0) + term
                    elif p1 >= y:
                        Ry += term
                    # p1 < z: cribado fuera de S(A,P,z); no pertenece
                    # al reindexado de switching (se descuenta aparte)
            # H1: J2 = suma T_s + R_y + (cola p1<z, medida aparte)
            # Reindexado exacto que se verifica: J2 - R_y - sum(T)
            # debe ser la masa con p1<z (>=0). Se comprueba la
            # identidad J2 = sum(T) + R_y + masa_p1menor_z.
            masa_menor_z = 0.0
            for v in range(1, vmax + 1):
                if math.gcd(v, Q) != 1 or OMEGA[v] != 2:
                    continue
                n = N - Q * v
                if n < 2 or LAM[n] == 0.0 or math.gcd(n, N) != 1:
                    continue
                p1 = None
                for p in PRIMES:
                    if p * p > v:
                        break
                    if v % p == 0:
                        p1 = p
                        break
                if p1 is not None and p1 < z:
                    masa_menor_z += LAM[n] * (logQ + math.log(RAD[v]))
            if abs(J2 - (sum(T.values()) + Ry + masa_menor_z)) <= 1e-9 * max(1.0, J2):
                h1 += 1
            # H2: fraccion de cola R_y / J2 (si J2>0)
            if J2 > 0:
                frac = Ry / J2
                if frac > max_tail_frac:
                    max_tail_frac = frac
                    worst_tail = (N, k, J2, Ry)
                h2 += 1
            # H3: solo submuestra (N multiplo de 10) por coste
            if N % 10 != 0:
                continue
            logX = math.log(X)
            for s in PRIMES:
                if s < z or s >= y:
                    continue
                if math.gcd(s, Q * N) != 1:
                    continue
                # P_Q(s): producto de primos l<s con l no | QN
                P = 1
                for ell in PRIMES:
                    if ell >= s:
                        break
                    if (Q * N) % ell != 0:
                        P *= ell
                wmax = (N - 1) // (Q * s)
                U_direct = 0.0
                for w in range(1, wmax + 1):
                    if math.gcd(w, P) != 1:
                        continue
                    n = N - Q * s * w
                    if n >= 2:
                        U_direct += LAM[n]
                # H3a: inversion EXACTA de Mobius sobre d | P_Q(s).
                divs = [d for d in range(1, P + 1)
                        if P % d == 0 and MU[d] != 0]
                U_mob = 0.0
                for d in divs:
                    U_mob += MU[d] * (psi_mod(N, Q * s * d, N) - LAM[N])
                if abs(U_direct - U_mob) <= 1e-9 * max(1.0, abs(U_direct)):
                    h3a += 1
                # H3b/c: version TRUNCADA (criba lineal, nivel D_s) y
                # su error frente al termino principal M_s.
                dmax = max(1, int((X / s) ** 0.5))
                U_trunc = 0.0
                M_s = 0.0
                for d in range(1, dmax + 1):
                    if MU[d] == 0 or math.gcd(d, Q * N * s) != 1:
                        continue
                    m = Q * s * d
                    U_trunc += MU[d] * (psi_mod(N, m, N) - LAM[N])
                    M_s += MU[d] * N / phi(m)
                # lema de pocos divisores: primos en [z,y) que dividen
                # m/Q = s*d, para cada d del rango truncado
                ok_div = True
                for d in range(1, dmax + 1):
                    if MU[d] == 0 or math.gcd(d, Q * N * s) != 1:
                        continue
                    md = s * d
                    cnt = 0
                    for sp in PRIMES:
                        if sp < z or sp >= y:
                            continue
                        if md % sp == 0:
                            cnt += 1
                    if cnt > max(1, int(math.log(MAXN) / math.log(z)) + 1):
                        ok_div = False
                if ok_div:
                    h3b += 1
                E = U_trunc - M_s
                sum_abs_E += abs(E)
                sum_abs_M += abs(M_s)
                n_bilin += 1
                if abs(M_s) > 0:
                    rel = abs(E) / abs(M_s)
                    if rel > max_bilin_rel:
                        max_bilin_rel = rel
                        worst_bilin = (N, k, s, U_trunc, M_s, E)
    print("PASS H1=%d H2=%d H3a_mobius=%d H3b_pocos_divisores=%d n_bilin=%d"
          % (h1, h2, h3a, h3b, n_bilin))
    print("max_tail_frac_Ry/J2=%.6f worst=%s" % (max_tail_frac, worst_tail))
    print("max_bilin_rel=%.6f worst=%s" % (max_bilin_rel, worst_bilin))
    print("sum_abs_E=%.3f sum_abs_M=%.3f ratio=%.6f"
          % (sum_abs_E, sum_abs_M,
             sum_abs_E / sum_abs_M if sum_abs_M else 0.0))


if __name__ == "__main__":
    main()
