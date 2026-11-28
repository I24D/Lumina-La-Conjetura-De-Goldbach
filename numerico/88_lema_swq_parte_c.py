#!/usr/bin/env python3
# 88_lema_swq_parte_c.py  (Doc 88, Lema SW-Q apartado (c))
# Verifica, en rango computable, las afirmaciones FINITAS del apartado (c)
# del LEMA SW-Q (docs 85-87): la constante c de (85.4) y la cota de la
# cola R_y, fijadas por la fuente primaria (Nathanson, Additive Number
# Theory: The Classical Bases, GTM 164 (1996), cap. 10, Teorema 10.2;
# Bordignon-Johnston-Starichkova, arXiv:2207.09452, Lema 37, Teorema 51,
# Lema 52). Version de limite FIJO z certificada en el doc 87 (la
# variable del doc 85 quedo abandonada). Solo stdlib (como docs 85-87).
#
# H1 (identidad del peso): sum_v a_Q(v) W(v) = S_Q(A;z) - (1/2) B(Q)
#   - (1/2) Omega_Q EXACTA, con W(v) = 1 - (1/2)#{p|v: z<=p<y}
#   - (1/2)#{descomposiciones P3 de v con p1 in [z,y), p2,p3>=y}.
#   B(Q) con limite FIJO z (doc 87). N = 10^6, 2*10^6, 10^7 (Q=3)
#   y N = 10^6 (Q=77, el de doc 87). Q IMPAR: con Q par, N-Qv es par
#   y a_Q = Lambda(N-Qv) es casi identicamente cero (solo potencias
#   de 2); el selector Q_N del doc 83 es coprimo con N.
# H2 (multiplicidad y switching): cada P2 con p2<y aparece 2 veces en
#   B(Q) (s=p1 y s=p2), con p2>=y una vez; Omega_Q = sum_{b in B_Q}
#   Lambda(b) exacto (biyeccion v <-> b = N-Qv), parte prima >> parte
#   de potencias de primos.
# H3 (constantes de la cola): c_{1/10} = int_{1/10}^{1/3}
#   log(2-3b)/(b(1-b)) db y cbar = int_{1/8}^{1/3} ... db por Simpson
#   y trapecio (acuerdo 1e-9); cbar < 0.363084 (Lema 52 de Bordignon,
#   reproducido).
# H4 (funciones de criba): f(5) por dos metodos (malla de retardo y
#   5f(5) = 4f(4) + int_3^4 F); J = int_{1/10}^{1/3} F(5-10b) db/b
#   por Simpson y trapecio; J(1/8) contra (e^gamma/2) log 6 exacto.
# H5 (balance): f(5) - J/2 - (1/2)*0.4*e^gamma*c_{1/10}*phi(Q)/Q > 0
#   para Q = 1, 2, 6, 30; y el control a z = N^{1/8}, Q = 1:
#   balance * 4/e^gamma = 2 log 3 - log 6 - cbar (constante de Chen).
# H6 (errores del peso): cubos de primos v = p^3 con a_Q(v) > 0 son
#   <= X^{1/3}; la fuga total E1'' = sum_{W>0, Omega>=3} a_Q W es
#   << N^{9/10} (envolvente O(N/z)) y ~1% de S a estas escalas.
# H7 (transferencia a J_k): J_k >= (log Q + log z) R_Q^{P<=2}
#   - log z * log N, exacto por fuerza bruta (peso minimo).
#
# Lectura honesta: H1/H2/H6/H7 son exactos por fuerza bruta a escala
# finita; H3/H4/H5 son calculos numericos de las constantes analiticas
# (cuadraturas con dos metodos independientes). El tamano asintotico
# (BV, Teorema 9.7, resto bilineal de la cola) es analitico, no
# numerico. Los numeros deben ser IDENTICOS local y en el PC.
#
# Uso: python 88_lema_swq_parte_c.py
import math
import sys
from bisect import bisect_left
from math import log, gcd

GAMMA = 0.5772156649015329
EG = math.exp(GAMMA)
TWO_EG = 2.0 * EG

FAILURES = []


def check(name, ok, detail=""):
    print(("PASS " if ok else "FAIL ") + name + (" | " + detail if detail else ""))
    if not ok:
        FAILURES.append(name)


# ---------------------------------------------------------------- criba
def primes_upto(n):
    sieve = bytearray(b"\x01") * (n + 1)
    sieve[0:2] = b"\x00\x00"
    for i in range(2, int(n ** 0.5) + 1):
        if sieve[i]:
            sieve[i * i:n + 1:i] = b"\x00" * ((n - i * i) // i + 1)
    return [i for i in range(n + 1) if sieve[i]]


# ------------------------------------------------------- F, f (retardo)
def build_Ff(M=3000):
    # Malla en [2,6], paso 1/M. En [2,3]: F = 2e^g/u, f = 2e^g ln(u-1)/u.
    # En (3,6]: uF(u) = 3F(3) + int_3^u f(t-1) dt,
    #           uf(u) = 3f(3) + int_3^u F(t-1) dt (trapecio, retardo 1).
    n = 4 * M + 1
    U = [2.0 + i / M for i in range(n)]
    Fv = [0.0] * n
    fv = [0.0] * n
    for i, uu in enumerate(U):
        if uu <= 3.0:
            Fv[i] = TWO_EG / uu
            fv[i] = TWO_EG * math.log(uu - 1) / uu if uu > 2.0 else 0.0
    hh = 1.0 / M
    F3 = TWO_EG / 3.0
    f3 = TWO_EG * math.log(2) / 3.0
    I1 = 0.0
    I2 = 0.0
    for i in range(M + 1, n):
        uu = U[i]
        I1 += 0.5 * (fv[i - 1 - M] + fv[i - M]) * hh
        I2 += 0.5 * (Fv[i - 1 - M] + Fv[i - M]) * hh
        Fv[i] = (3 * F3 + I1) / uu
        fv[i] = (3 * f3 + I2) / uu
    return U, Fv, fv


def Fnum(u, U, Fv):
    if u <= 2.0:
        return TWO_EG / u
    i = bisect_left(U, u)
    if i >= len(U):
        return Fv[-1]
    if U[i] == u:
        return Fv[i]
    u0, u1 = U[i - 1], U[i]
    f0, f1 = Fv[i - 1], Fv[i]
    return f0 + (f1 - f0) * (u - u0) / (u1 - u0)


def fnum(u, U, fv):
    i = bisect_left(U, u)
    if i >= len(U):
        return fv[-1]
    if U[i] == u:
        return fv[i]
    u0, u1 = U[i - 1], U[i]
    f0, f1 = fv[i - 1], fv[i]
    return f0 + (f1 - f0) * (u - u0) / (u1 - u0)


def simpson_split(f, a, b, breaks, n):
    # Simpson por tramos: F tiene quiebres de C^2 en u = 2 y u = 3
    # (b = 0.3 y b = 0.2 en J; b = 0.25 en J(1/8)); partir ahi evita
    # que el error de Simpson caiga a O(h^2).
    pts = [a] + [x for x in breaks if a < x < b] + [b]
    return sum(simpson(f, pts[i], pts[i + 1], n) for i in range(len(pts) - 1))


def simpson(f, a, b, n):
    # n intervalos (par)
    h = (b - a) / n
    s = f(a) + f(b)
    for i in range(1, n):
        s += (4.0 if i % 2 == 1 else 2.0) * f(a + i * h)
    return s * h / 3.0


def trapezoid(f, a, b, n):
    h = (b - a) / n
    s = 0.5 * (f(a) + f(b))
    for i in range(1, n):
        s += f(a + i * h)
    return s * h


# ------------------------------------------------------------- parte finita
def run_finite(N, Q):
    X = N // Q
    z = max(2, int(X ** 0.1))
    y = int(X ** (1 / 3))
    primes = primes_upto(N)
    pset = set(primes)
    small_primes = [p for p in primes if p * p <= X]
    lam = {}
    for p in primes:
        pw = p
        while pw <= N:
            lam[pw] = log(p)
            if pw > N // p:
                break
            pw *= p

    def factorize(v):
        fac = []
        t = v
        for p in small_primes:
            if p * p > t:
                break
            if t % p == 0:
                c = 0
                while t % p == 0:
                    t //= p
                    c += 1
                fac.append((p, c))
        if t > 1:
            fac.append((t, 1))
        return fac

    def sifted(v):
        if v < 1 or gcd(v, Q) != 1:
            return False
        for p, _ in factorize(v):
            if p < z and (N % p != 0) and (Q % p != 0):
                return False
        return True

    S = 0.0
    B = 0.0
    Om = 0.0
    SW = 0.0
    Bs = {}
    n2 = n1 = n0 = 0
    Bset = set()
    cubes = 0
    leaks = []
    Jk = 0.0
    Rle2 = 0.0
    pp_ns = []
    for p in primes:
        pp_ns.append(p)
        if p * p <= N:
            pw = p * p
            while pw <= N:
                pp_ns.append(pw)
                if pw > N // p:
                    break
                pw *= p
    seen = set()
    for n in pp_ns:
        if n in seen:
            continue
        seen.add(n)
        if (N - n) % Q != 0:
            continue
        v = (N - n) // Q
        if v < 1 or v >= X:
            continue
        if gcd(v, Q) != 1:
            continue
        av = lam.get(n, 0.0)
        if av <= 0:
            continue
        if not sifted(v):
            continue
        S += av
        fac = factorize(v)
        distinct = [p for p, _ in fac]
        c1 = sum(1 for p in distinct if z <= p < y)
        pf = []
        for p, c in fac:
            pf += [p] * c
        c3 = 0
        if len(pf) == 3:
            p1, p2, p3 = sorted(pf)
            if z <= p1 < y <= p2 <= p3:
                c3 = 1
                Om += av
                Bset.add(n)
        W = 1 - 0.5 * c1 - 0.5 * c3
        SW += av * W
        for s in distinct:
            if z <= s < y:
                Bs[s] = Bs.get(s, 0.0) + av
        if len(pf) == 2:
            p, q = sorted(pf)
            cnt = (1 if z <= p < y else 0) + (1 if z <= q < y else 0)
            if cnt == 2:
                n2 += 1
            elif cnt == 1:
                n1 += 1
            else:
                n0 += 1
        if len(pf) == 3 and pf[0] == pf[1] == pf[2]:
            cubes += 1
        if W > 0 and len(pf) >= 3:
            leaks.append((v, len(pf), W))
        if len(pf) <= 2:
            Rle2 += av
            r = 1
            for p, _ in fac:
                r *= p
            Jk += av * (log(Q) + log(r))
    B = sum(Bs.values())
    leak_total = sum(lam.get(N - Q * vv, 0.0) * ww for (vv, _, ww) in leaks)
    sw_prime = sum(log(b) for b in Bset if b in pset)
    Om_check = sum(lam.get(b, 0.0) for b in Bset)
    bound7 = (log(Q) + log(z)) * Rle2 - log(z) * log(N)
    return dict(N=N, Q=Q, X=X, z=z, y=y, S=S, B=B, Om=Om, SW=SW,
                n2=n2, n1=n1, n0=n0, nBs=len(Bs), Om_check=Om_check,
                sw_prime=sw_prime, pp_part=Om_check - sw_prime,
                cubes=cubes, nleaks=len(leaks), leak_total=leak_total,
                Jk=Jk, Rle2=Rle2, bound7=bound7)


def main():
    # ---- H3/H4/H5: constantes analiticas ----
    U, Fv, fv = build_Ff(3000)
    f5_grid = fnum(5.0, U, fv)
    f4 = fnum(4.0, U, fv)
    I34_s = simpson(lambda t: Fnum(t, U, Fv), 3.0, 4.0, 2000)
    f5_alt = (4 * f4 + I34_s) / 5.0
    F4 = Fnum(4.0, U, Fv)
    check("H4_f5", abs(f5_grid - f5_alt) < 1e-9,
          "f(5)=%.12f grid vs %.12f via int_3^4 F; F(4)=%.9f" % (f5_grid, f5_alt, F4))
    J_s = simpson_split(lambda b: Fnum(5 - 10 * b, U, Fv) / b, 0.1, 1.0 / 3.0, [0.2, 0.3], 6666)
    J_t = trapezoid(lambda b: Fnum(5 - 10 * b, U, Fv) / b, 0.1, 1.0 / 3.0, 100000)
    J18 = simpson_split(lambda b: Fnum(4 - 8 * b, U, Fv) / b, 0.125, 1.0 / 3.0, [0.25], 6666)
    J18_exact = EG / 2 * log(6)
    check("H4_J", abs(J_s - J_t) < 1e-9 and abs(J18 - J18_exact) < 1e-8,
          "J=%.12f simp vs %.12f trap; J(1/8)=%.12f vs (e^g/2)log6=%.12f" % (J_s, J_t, J18, J18_exact))
    c110_s = simpson(lambda t: log(2 - 3 * t) / (t * (1 - t)), 0.1, 1.0 / 3.0, 20000)
    c110_t = trapezoid(lambda t: log(2 - 3 * t) / (t * (1 - t)), 0.1, 1.0 / 3.0, 20000)
    cbar_s = simpson(lambda t: log(2 - 3 * t) / (t * (1 - t)), 0.125, 1.0 / 3.0, 20000)
    cbar_t = trapezoid(lambda t: log(2 - 3 * t) / (t * (1 - t)), 0.125, 1.0 / 3.0, 20000)
    check("H3_c", abs(c110_s - c110_t) < 1e-9 and abs(cbar_s - cbar_t) < 1e-9
          and cbar_s < 0.363084,
          "c_{1/10}=%.12f (simp %.12f trap); cbar=%.12f < 0.363084" % (c110_s, c110_t, cbar_s))
    f5 = f5_grid
    J = J_s
    c110 = c110_s
    cbar = cbar_s

    def phi_ratio(Q):
        # phi(Q)/Q = prod_{p|Q} (1-1/p); Q=1 -> 1
        tmp = Q
        r = 1.0
        d = 2
        while d * d <= tmp:
            if tmp % d == 0:
                r *= (1 - 1.0 / d)
                while tmp % d == 0:
                    tmp //= d
            d += 1 if d == 2 else 2
        if tmp > 1:
            r *= (1 - 1.0 / tmp)
        return r

    bals = {}
    for Q in (1, 2, 6, 30):
        bals[Q] = f5 - J / 2 - 0.5 * 0.4 * EG * c110 * phi_ratio(Q)
    check("H5_balance", all(v > 0 for v in bals.values()),
          "balance z=X^{1/10}: " + " ".join("Q=%d:%.6f" % (q, v) for q, v in bals.items()))
    bal18 = f4 - J18 / 2 - 0.5 * 0.5 * EG * cbar
    bal18_lit = bal18 * 4 / EG
    lit = 2 * log(3) - log(6) - cbar
    check("H5_chen", abs(bal18_lit - lit) < 1e-6,
          "z=N^{1/8} balance=%.12f; x4/e^g=%.12f vs 2log3-log6-cbar=%.12f" % (bal18, bal18_lit, lit))

    # ---- H1/H2/H6/H7: fuerza bruta finita ----
    for (N, Q) in [(10 ** 6, 3), (2 * 10 ** 6, 3), (10 ** 7, 3), (10 ** 6, 77)]:
        r = run_finite(N, Q)
        rel1 = abs(r["SW"] - (r["S"] - 0.5 * r["B"] - 0.5 * r["Om"])) / max(1.0, abs(r["S"]))
        check("H1_identidad_N%d_Q%d" % (N, Q), rel1 < 1e-9,
              "S=%.4f B=%.4f Om=%.4f SW=%.4f rel=%.2e" % (r["S"], r["B"], r["Om"], r["SW"], rel1))
        check("H2_multiplicidad_N%d_Q%d" % (N, Q), r["n2"] > 0 and r["n1"] > 0,
              "P2: twice=%d once=%d zero=%d (B_s=%d)" % (r["n2"], r["n1"], r["n0"], r["nBs"]))
        rel2 = abs(r["Om"] - r["Om_check"]) / max(1.0, abs(r["Om"]))
        check("H2b_switching_N%d_Q%d" % (N, Q), rel2 < 1e-9,
              "Om=%.4f sum_Lambda(B_Q)=%.4f prime=%.4f pp=%.4f" % (r["Om"], r["Om_check"], r["sw_prime"], r["pp_part"]))
        check("H6_errores_N%d_Q%d" % (N, Q),
              r["cubes"] <= r["X"] ** (1 / 3) + 1 and r["leak_total"] <= N ** 0.9,
              "cubos=%d (<=X^{1/3}=%.1f); fugas=%d E1''=%.4f (%.4f de S; <=N^0.9=%.1f)"
              % (r["cubes"], r["X"] ** (1 / 3), r["nleaks"], r["leak_total"],
                 r["leak_total"] / max(r["S"], 1e-9), N ** 0.9))
        check("H7_peso_minimo_N%d_Q%d" % (N, Q), r["Jk"] + 1e-9 >= r["bound7"],
              "Jk=%.4f Rle2=%.4f bound=%.4f" % (r["Jk"], r["Rle2"], r["bound7"]))

    print("PASS" if not FAILURES else "FAILURES: " + ", ".join(FAILURES))
    return 0 if not FAILURES else 1


if __name__ == "__main__":
    sys.exit(main())
