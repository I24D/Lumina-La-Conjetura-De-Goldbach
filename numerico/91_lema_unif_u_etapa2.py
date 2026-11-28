#!/usr/bin/env python3
# 91 - Lema UNIF-U etapa 2: las cinco rutas, medidas una por una.
# Solo stdlib. Reutiliza F/f por retardo (docs 87/88) y V^+ (doc 86).
import math
from bisect import bisect_left
from fractions import Fraction

TWO_EG = 2 * math.exp(0.5772156649015329)
EG = math.exp(0.5772156649015329)
FAILURES = []


def check(name, ok, detail=""):
    print(("PASS " if ok else "FAIL ") + name + (" | " + detail if detail else ""))
    if not ok:
        FAILURES.append(name)


def primes_upto(n):
    bs = bytearray(b"\x01") * (n + 1)
    bs[0:2] = b"\x00\x00"
    for i in range(2, int(n ** 0.5) + 1):
        if bs[i]:
            bs[i * i:n + 1:i] = b"\x00" * ((n - i * i) // i + 1)
    return [i for i in range(2, n + 1) if bs[i]]


def mobius_upto(n):
    mu = [1] * (n + 1)
    pr = primes_upto(n)
    is_c = [False] * (n + 1)
    for p in pr:
        for m in range(p * p, n + 1, p * p):
            is_c[m] = True
    for p in pr:
        for m in range(p, n + 1, p):
            mu[m] *= -1
    for m in range(2, n + 1):
        if is_c[m]:
            mu[m] = 0
    return mu


def build_Ff(M=3000):
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
    return Fv[i - 1] + (Fv[i] - Fv[i - 1]) * (u - u0) / (u1 - u0)


def fnum(u, U, fv):
    i = bisect_left(U, u)
    if i >= len(U):
        return fv[-1]
    if U[i] == u:
        return fv[i]
    u0, u1 = U[i - 1], U[i]
    return fv[i - 1] + (fv[i] - fv[i - 1]) * (u - u0) / (u1 - u0)


def simpson(f, a, b, n):
    if n % 2:
        n += 1
    h = (b - a) / n
    s = f(a) + f(b)
    for i in range(1, n):
        s += (4 if i % 2 else 2) * f(a + i * h)
    return s * h / 3


def in_V(d, D, parity, MU):
    # Doc 86, Teor. 9.3, beta=2: d = p1>...>pk libre de cuadrados;
    # para m de la paridad dada: p_m^2 * (p1...pm) < D.
    if d == 1:
        return True
    if MU[d] == 0:
        return False
    ps = []
    x = d
    p = 2
    while p * p <= x:
        if x % p == 0:
            ps.append(p)
            while x % p == 0:
                x //= p
        p += 1 if p == 2 else 2
    if x > 1:
        ps.append(x)
    ps.sort(reverse=True)
    if ps[0] >= D:
        return False
    prod = 1
    for i, pp in enumerate(ps):
        m = i + 1
        prod *= pp
        if (m % 2 == 1) == (parity == 1):
            if not (pp * pp * prod < D):
                return False
    return True


def main():
    U, Fv, fv = build_Ff()
    f5 = fnum(5.0, U, fv)
    F4 = Fnum(4.0, U, Fv)
    c110 = simpson(lambda b: math.log(2 - 3 * b) / (b * (1 - b)), 0.1, 1 / 3, 2000)
    J0 = simpson(lambda b: Fnum(5 - 10 * b, U, Fv) / b, 0.1, 1 / 3, 2000)
    check("H0-reproduce (doc87/88)", abs(f5 - 0.9982417256) < 1e-6 and abs(F4 - 1.021641553) < 1e-6
          and abs(c110 - 0.4909952010) < 1e-7 and abs(J0 - 1.5238762503) < 1e-6,
          "f5=%.10f F4=%.9f c110=%.10f J0=%.10f" % (f5, F4, c110, J0))

    # ---- H0: criterio absoluto vs relativo ----
    # R_abs(N) = N^c/(ln N)^{A'-2}: cociente cota absoluta / objetivo.
    rows = []
    for c in (0.02, 0.05, 0.10, 2 / 17):
        vals = []
        for lg10N in (12, 24, 100, 560, 3600):
            lnN = lg10N * math.log(10)
            vals.append(c * lg10N - 8 * math.log10(lnN))
        rows.append((c, vals))
    # Doc 90, Prop. 90.1: el cociente primero CAE y cruza 1 solo en N*;
    # a escala computable es <1, a N -> infinito es >1 para todo c>0 fijo.
    ok = all(v[0] < 0 for _, v in rows) and all(v[4] > 0 for _, v in rows) \
        and all(v[3] > v[2] for _, v in rows)
    check("H0 (criterio)", ok,
          "log10 R_abs c=0.02/0.05/0.10/2-17 a N=1e12/1e24/1e100/1e560/1e3600: " +
          "; ".join("c=%.4g: %s" % (c, "/".join("%.2f" % v for v in vals)) for c, vals in rows))

    # ---- H1: holgura de nivel BFI / Maynard II ----
    cs = [0.02, 0.039884, 0.05, 0.10, Fraction(2, 17)]
    # (calculo exacto con Fractions para c racionales; 0.039884 aproximado)
    def head(c, lvl):
        return float(lvl) - (0.5 + c / 2)
    det = "; ".join("c=%s: 4/7:%+.4f 7/12:%+.4f 3/5:%+.4f" % (
        ("2/17" if c == 2 / 17 else "%.4f" % c), head(float(c), Fraction(4, 7)),
        head(float(c), Fraction(7, 12)), head(float(c), Fraction(3, 5))) for c in [0.02, 0.039884, 0.05, 0.10, 2 / 17])
    check("H1 (nivel BFI)", all(head(float(c), Fraction(4, 7)) > 0 for c in [0.02, 0.039884, 0.05, 0.10, 2 / 17]), det)

    # ---- H2: dependencia en s de V^+(D/s) (X=1e12) ----
    X = 1e12
    z = X ** 0.1
    y = X ** (1 / 3)
    D = X ** 0.5
    MU = mobius_upto(20000)
    Pz = [p for p in primes_upto(20) if p < z and p > 2]
    cand = [1]
    for p in Pz:
        cand += [d * p for d in list(cand) if d * p < 20000]
    cand = sorted(set(cand))
    res = {}
    for s in (17, 19, 101, 9973):
        L = D / s
        vs = [d for d in cand if d < L and in_V(d, L, 1, MU)]
        res[s] = vs
    counts = {s: len(v) for s, v in res.items()}
    wit = 15 in res[17] and 15 not in res[9973]
    check("H2 (V^+ depende de s)", wit and counts[17] > 5 * counts[9973],
          "|V^+| s=17/19/101/9973: %s; testigo d=15 en V^+(D/17) y no en V^+(D/9973): %s; max d s=17: %d, s=9973: %d"
          % (counts, wit, max(res[17]), max(res[9973])))

    # ---- H3: coste de los bloques diadicos ----
    # J_bloque(X) = sum_bloques F(u_S) * ln(1 + ln2/ln S), u_S = log(D/(2S))/log z.
    def J_block(X):
        lnX = math.log(X)
        zz = X ** 0.1
        yy = X ** (1 / 3)
        DD = X ** 0.5
        lnz = math.log(zz)
        tot = 0.0
        S = zz
        while S < yy:
            uS = (math.log(DD) - math.log(2 * S)) / lnz
            if uS < 1.0:
                uS = 1.0
            tot += Fnum(uS, U, Fv) * math.log(1 + math.log(2) / math.log(S))
            S *= 2
        return tot
    bal3 = lambda J: f5 - 0.5 * J - 0.5 * 0.4 * EG * c110 * (2 / 3)
    jb = {e: J_block(10.0 ** e) for e in (6, 12, 24, 100)}
    check("H3 (bloques)", abs(jb[100] - J0) < 0.02 and bal3(jb[100]) > 0.05,
          "J_bloque X=1e6/1e12/1e24/1e100: %s (J=%.10f); Bal_bloque(Q=3): %s (doc88: 0.11965)"
          % ("/".join("%.6f" % jb[e] for e in (6, 12, 24, 100)),
             J0, "/".join("%.5f" % bal3(jb[e]) for e in (6, 12, 24, 100))))

    # ---- H4: Maynard III (ventana y error debil) ----
    # Masa (beta, eta) con densidad (1/beta) d beta d eta; beta in [0.1,1/3],
    # eta in [0, 0.5-beta]; log(Asd)/log N = c + (1-c)(beta+eta).
    def mass_frac(c, cond):
        tot = 0.0
        num = 0.0
        NB = 240
        for i in range(NB):
            b0 = 0.1 + (1 / 3 - 0.1) * i / NB
            b1 = 0.1 + (1 / 3 - 0.1) * (i + 1) / NB
            bm = 0.5 * (b0 + b1)
            w = (b1 - b0) / bm
            etamax = 0.5 - bm
            if etamax <= 0:
                continue
            tot += w * etamax
            # fraccion de eta que cumple cond(beta, eta)
            NE = 200
            hit = 0
            for j in range(NE):
                eta = etamax * (j + 0.5) / NE
                if cond(bm, eta):
                    hit += 1
            num += w * etamax * hit / NE
        return num / tot if tot else 0.0
    fr = {}
    for c in (0.02, 0.05, 0.10):
        fr[c] = mass_frac(c, lambda b, e, c=c: c + (1 - c) * (b + e) > 0.501)
    # Ventana Teor. 1.2, delta=1/1000, c=0.05: eta en ventana y nivel total <= 0.501.
    dl = 0.001
    c = 0.05
    def in_win12(b, e):
        lo = (0.1 + 12 * dl - (1 - c) * b) / (1 - c)
        hi = (0.1 - 4 * dl - 0.6 * (1 - c) * b) / (1 - c)
        return lo < e < hi and c + (1 - c) * (b + e) <= 0.501
    cov12 = mass_frac(c, in_win12)
    rw = {}
    for lg10N in (12, 24, 100):
        lnN = lg10N * math.log(10)
        rw[lg10N] = dl * (10.0 ** (0.05 * lg10N)) * lnN ** 2
    check("H4 (Maynard III)", all(fr[c] > 0.01 for c in fr) and cov12 < 0.10,
          "masa sobre x^0.501: c=0.02:%.3f 0.05:%.3f 0.10:%.3f; cobertura Teor.1.2 (c=0.05): %.4f; R_debil=delta*N^c*(lnN)^2 N=1e12/1e24/1e100: %.2g/%.2g/%.2g"
          % (fr[0.02], fr[0.05], fr[0.10], cov12, rw[12], rw[24], rw[100]))

    # ---- H5: Tipo II (Q^2) y segundo momento a a fija ----
    q2 = {c: (10.0 ** (c * 12)) for c in (0.02, 0.05, 0.10, 2 / 17)}
    # --- H5b: segundo momento fijo-en-a (modelo aleatorio, NO teorema) ---
    # Premisa (NO teorema; BDH promedia sobre a): sum_{q<=Q} |E|^2 << N (ln N)^C.
    # Cauchy: sum|E| <= sqrt(Q) sqrt(N (ln N)^C) = N^{3/4+c/4} (ln N)^{C/2};
    # cociente al termino principal (N/A)/log^2 N = N^{5c/4 - 1/4} (ln N)^{C/2+2}.
    # -> 0 sii c < 1/5 (precision del doc 90); a escala computable aun es > 1.
    C = 5
    rrand = {}
    for c in (0.05, 2 / 17):
        vals = []
        for lg10N in (12, 24, 100, 200):
            lnN = lg10N * math.log(10)
            vals.append((5 * c / 4 - 1 / 4) * lg10N + (C / 2 + 2) * math.log10(lnN))
        rrand[c] = vals
    # GRH por q: sum|E| << Q sqrt(N) (ln N)^2 = N^{1+c/2}(ln N)^2;
    # ratio al objetivo (N/A)/log^2 N = N^{3c/2} (ln N)^4.
    grh = {}
    for lg10N in (12, 24, 100):
        lnN = lg10N * math.log(10)
        grh[lg10N] = 1.5 * 0.05 * lg10N + 4 * math.log10(lnN)
    check("H5 (Tipo II / 2o momento)",
          all(v > 0 for v in q2.values())
          and all(v > 0 for v in rrand[0.05][:1] + rrand[2 / 17][:1])
          and rrand[0.05][2] < 0 and rrand[2 / 17][3] < 0
          and all(v > 0 for v in grh.values()),
          "Q^2/N=N^c a N=1e12: %s; log10 R_aleat c=0.05 N=1e12/1e24/1e100/1e200: %.2f/%.2f/%.2f/%.2f; c=2/17: %.2f/%.2f/%.2f/%.2f; log10 ratio GRH N=1e12/1e24/1e100: %.2f/%.2f/%.2f"
          % ("/".join("%.2f" % q2[c] for c in (0.02, 0.05, 0.10, 2 / 17)),
             rrand[0.05][0], rrand[0.05][1], rrand[0.05][2], rrand[0.05][3],
             rrand[2 / 17][0], rrand[2 / 17][1], rrand[2 / 17][2], rrand[2 / 17][3],
             grh[12], grh[24], grh[100]))

    # ---- H6: densidad de ceros / franja VK ----
    w = {}
    wl = {}
    for lg10N in (12, 24, 100):
        lnN = lg10N * math.log(10)
        w[lg10N] = 1.0 / (57.54 * lnN ** (2 / 3) * math.log(lnN) ** (1 / 3))
        wl[lg10N] = 1.0 / (10.5 * 0.5 * lnN + 61.5 * lnN ** (2 / 3) * math.log(lnN) ** (1 / 3))
    check("H6 (franja VK)", all(w[lg] < 0.02 for lg in w),
          "w(N) Ford 1/57.54 N=1e12/1e24/1e100: %.5f/%.5f/%.5f; c/w c=0.02: %.1f/%.1f/%.1f, c=0.05: %.1f/%.1f/%.1f; w_L(1e12)=%.5f; GM 30/13=%.4f vs Huxley 12/5=%.1f"
          % (w[12], w[24], w[100],
             0.02 / w[12], 0.02 / w[24], 0.02 / w[100],
             0.05 / w[12], 0.05 / w[24], 0.05 / w[100], wl[12], 30 / 13, 12 / 5))

    print("PASS" if not FAILURES else "FAILURES: %s" % FAILURES)


if __name__ == "__main__":
    main()
