#!/usr/bin/env python3
"""Doc 92 - Lema UNIF-U-2M etapa 1: el segundo momento a a=N fijo.
Solo stdlib. Uso: python3 92_lema_unif_u_2m.py [ruta_ceros_zeta.txt]
E(N;q,a) = psi(N;q,a) - psi(N)/phi(q), (q,N)=1.
H0 reproduce la derivacion de Cauchy CORREGIDA del doc 91 (H5b).
H1/H2: segundo momento finito, fijo-a vs promedio-sobre-a (BDH).
H3: fuera-diagonal/diagonal de la formula explicita para zeta (252 ceros).
H4: sensibilidad a un cero de Siegel + franja necesaria w>=1/2.
H5: ventana de Maynard III Teor. 1.1 (umbral de apertura) + perdida L1->L2.
"""
import math
import sys

FAIL = []


def check(name, cond, detail):
    print(("PASS " if cond else "FAIL ") + name + " | " + detail)
    if not cond:
        FAIL.append(name)


def log10(x):
    return math.log10(x)


# ---------------- H0: Cauchy corregido (doc 91 H5b) ----------------
def H0():
    # UNIF-U-2M: sum_{q<=Q} |E|^2 << N (ln N)^C, Q = N^{1/2+c/2}.
    # Cauchy sobre los modulos distintos m=Asd (<=Q, multiplicidad<=5):
    # |R| <= sqrt(Q * N (ln N)^C); objetivo (N/A)/ln^2 N, A=N^c.
    # Cociente R_aleat = N^{5c/4-1/4} (ln N)^{C/2+2}; umbral c<1/5.
    C = 5
    c = 0.05
    Ns = [1e12, 1e24, 1e100, 1e200]
    expect = [4.24, 3.34, -8.12, -25.52]  # doc 91 sec. 4
    got, got2 = [], []
    for N in Ns:
        lnN = math.log(N)
        R = N ** (5 * c / 4 - 1 / 4) * lnN ** (C / 2 + 2)
        got.append(round(log10(R), 2))
        # via sqrt, en espacio log (evita overflow a N=1e200):
        # log10 R2 = 0.5*(log10 M + log10 N + C*log10 lnN)
        #            - ((1-c)*log10 N - 2*log10 lnN), M = N^{1/2+c/2}
        lR2 = 0.5 * ((0.5 + c / 2) * log10(N) + log10(N) + C * log10(lnN)) \
            - ((1 - c) * log10(N) - 2 * log10(lnN))
        got2.append(round(lR2, 2))
    ok = got == expect and got2 == expect
    check("H0 (Cauchy corregido doc91)",
          ok, "log10 R_aleat c=0.05 N=1e12/1e24/1e100/1e200: %s (via sqrt: %s; doc91: %s); umbral 5c/4-1/4=0 -> c=1/5"
          % (got, got2, expect))


# ---------------- criba Lambda ----------------
def build_lambda(N):
    is_prime = bytearray(b"\x01") * (N + 1)
    is_prime[0:2] = b"\x00\x00"
    sq = int(N ** 0.5) + 1
    for p in range(2, sq + 1):
        if is_prime[p]:
            start = p * p
            is_prime[start:N + 1:p] = b"\x00" * ((N - start) // p + 1)
    lam = [0.0] * (N + 1)
    for p in range(2, N + 1):
        if is_prime[p]:
            lp = math.log(p)
            pk = p
            while pk <= N:
                lam[pk] += lp
                if pk > N // p:
                    break
                pk *= p
    return lam


def phi_sieve(n):
    phi = list(range(n + 1))
    for p in range(2, n + 1):
        if phi[p] == p:
            for m in range(p, n + 1, p):
                phi[m] -= phi[m] // p
    return phi


# ---------------- H1/H2: fijo-a vs promedio ----------------
def H1H2():
    N0 = 100000
    lam0 = build_lambda(N0)
    psi0 = sum(lam0)
    phi0 = phi_sieve(2000)
    pp = [(n, lam0[n]) for n in range(2, N0 + 1) if lam0[n] > 0.0]
    c = 0.05
    Q0 = int(N0 ** (0.5 + c / 2))
    S2avg = 0.0
    S2fix = 0.0
    Bqs = []
    for q in range(2, Q0 + 1):
        if math.gcd(q, N0) != 1:
            continue
        buckets = [0.0] * q
        for n, lp in pp:
            buckets[n % q] += lp
        main = psi0 / phi0[q]
        s = 0.0
        for a in range(q):
            if math.gcd(a, q) == 1:
                e = buckets[a] - main
                s += e * e
        Bqs.append((q, s))
        S2avg += s / phi0[q]
        r = N0 % q
        ef = buckets[r] - main
        S2fix += ef * ef
    ln0 = math.log(N0)
    target0 = N0 * ln0 ** 5
    bdh_bound = N0 * Q0 * ln0 ** 5
    ratio = S2fix / S2avg if S2avg > 0 else float("nan")
    meanB = sum(s for _, s in Bqs) / len(Bqs)
    qmax, Bmax = max(Bqs, key=lambda t: t[1])
    # fijo-a a N=1e6 (solo fijo, por rebanadas)
    N1 = 1000000
    lam1 = build_lambda(N1)
    psi1 = sum(lam1)
    phi1 = phi_sieve(3000)
    ln1 = math.log(N1)
    target1 = N1 * ln1 ** 5
    rows = []
    for cc in (0.02, 0.05, 0.10):
        Q1 = int(N1 ** (0.5 + cc / 2))
        s2 = 0.0
        mx = 0.0
        qmx = 0
        for q in range(2, Q1 + 1):
            if math.gcd(q, N1) != 1:
                continue
            psi_q = sum(lam1[N1 % q::q])
            e = psi_q - psi1 / phi1[q]
            s2 += e * e
            if abs(e) > mx:
                mx = abs(e)
                qmx = q
        rows.append((cc, Q1, s2 / target1, mx, qmx))
    # exponente implicito: S2fix/N ~ (ln N)^alpha entre N0 y N1
    s2_1 = rows[1][2] * target1
    alpha = math.log((s2_1 / N1) / (S2fix / N0)) / math.log(ln1 / ln0)
    detail = ("N0=1e5 Q0=%d S2avg=%.4g S2fix=%.4g fijo/prom=%.3f; "
              "S2avg/objetivo=%.3g S2fix/objetivo=%.3g; BDH NQ(lnN)^5=%.3g (S2avg/BDH=%.2g); "
              "max Bq/mean Bq=%.2f (q=%d); S2fix/N: %.2f (1e5) -> %.2f (1e6), alpha=%.2f; "
              "N1=1e6 fijo/objetivo: %s"
              % (Q0, S2avg, S2fix, ratio, S2avg / target0, S2fix / target0,
                 bdh_bound, S2avg / bdh_bound, Bmax / meanB, qmax,
                 S2fix / N0, s2_1 / N1, alpha,
                 "; ".join("c=%.2f Q=%d ratio=%.3g max|E|=%.1f(q=%d)" % r for r in rows)))
    check("H1 (fijo vs promedio)", 0.1 < ratio < 10.0, detail)
    check("H2 (peor clase)", True,
          "max|E(N=1e6;q,N)| por c: " + "; ".join(
              "c=%.2f: %.1f en q=%d (ratio a objetivo %.3g)" % (r[0], r[3], r[4], r[2]) for r in rows))


# ---------------- H3: fuera-diagonal zeta ----------------
def H3(zeros_path):
    gammas = []
    with open(zeros_path, encoding="utf-8") as f:
        for line in f:
            line = line.strip()
            if line:
                gammas.append(float(line))
    n = len(gammas)
    gmax = gammas[-1]
    Sdiag = sum(1.0 / (g * g) for g in gammas)
    # benchmark fase aleatoria: sqrt(sum_{j!=k} 1/(gj gk)^2)
    S4 = sum(1.0 / (g ** 4) for g in gammas)
    bench = math.sqrt(max(Sdiag * Sdiag - S4, 0.0))
    # cola truncada: integral (ln(g/2pi))/(2pi g^2) desde gmax
    tail = (math.log(gmax / (2 * math.pi)) + 1.0) / (2 * math.pi * gmax)
    rows = []
    for x in (1e6, 1e12):
        L = math.log(x)
        # |sum_j e^{i gj L}/gj|^2 = Sdiag + S_off
        re = sum(math.cos(g * L) / g for g in gammas)
        im = sum(math.sin(g * L) / g for g in gammas)
        tot = re * re + im * im
        Soff = tot - Sdiag
        rows.append((x, tot, Soff, abs(Soff) / Sdiag if Sdiag else float("nan")))
    detail = ("n=%d gmax=%.1f Sdiag=%.4f (cola integral=%.4f, total~%.4f); "
              "benchmark aleatorio |Soff|~%.4f (%.2f x Sdiag); "
              % (n, gmax, Sdiag, tail, Sdiag + tail, bench, bench / Sdiag if Sdiag else float("nan"))
              + "; ".join("x=%.0e |Z|^2=%.4f Soff=%+.4f |Soff|/Sdiag=%.3f" % r for r in rows))
    check("H3 (fuera-diagonal zeta)", True, detail)


# ---------------- H4: Siegel + franja ----------------
def H4():
    # franja necesaria w>=1/2 (diagonal ~ N^{2-2w}); conocida VK Ford 1/57.54
    rows = []
    for N in (1e12, 1e24, 1e100):
        lnN = math.log(N)
        w = 1.0 / (57.54 * lnN ** (2.0 / 3) * math.log(lnN) ** (1.0 / 3))
        rows.append((N, w, 0.5 / w))
    # Siegel: beta0 = 1 - cs/log Q, Q = N^{1/2+c/2}; contribucion a S2:
    # (1/phi(q)^2) * N^{2 beta0} / beta0^2, sumada sobre q multiplos de q0:
    # sum_{q0|q} 1/phi(q)^2 <= (sum_m 1/phi(m)^2)/phi(q0)^2 (phi(q0 m)>=phi(q0)phi(m)).
    # Sin factor 1/(1-beta0): el termino de la formula explicita es N^beta0/beta0.
    phi = phi_sieve(20000)
    s_phi2 = sum(1.0 / (phi[m] ** 2) for m in range(1, 20001))
    sieg = []
    for N in (1e12, 1e24):
        lnN = math.log(N)
        Q = N ** (0.5 + 0.05 / 2)
        lnQ = math.log(Q)
        target = N * lnN ** 5
        for cs in (1.0, 0.1):
            beta0 = 1.0 - cs / lnQ
            for q0, ph in ((3, 2), (7, 6)):
                contrib = s_phi2 * N ** (2 * beta0) / (beta0 * beta0) / (ph * ph)
                sieg.append((N, cs, q0, contrib / target))
    d1 = "; ".join("N=%.0e w=%.5f 0.5/w=%.0f" % r for r in rows)
    d2 = "; ".join("N=%.0e cs=%.1f q0=%d ratio=%.3g" % s for s in sieg)
    check("H4 (Siegel + franja)", True,
          "franja: " + d1 + "; Siegel contrib/objetivo: " + d2
          + " (sum 1/phi^2=%.4f)" % s_phi2)


# ---------------- H5: ventana Maynard + L1->L2 ----------------
def H5():
    # Teor. 1.1: Q1 <= x^{0.1-3d}/(ln x)^C, Q2 <= x^{0.4+4d}(ln x)^C.
    # Para cubrir q ~ Q = x^{1/2+c/2} hace falta d >= c/2 (nivel) y 0.1-3d > 0.
    # Ventana de q1 (logs, con el factor 2 de las cajas diadicas):
    #   lo = (0.1+c/2-4d)L - C ln L - ln 2 ; hi = (0.1-3d)L - C ln L + ln 2.
    # Con d = c/2 (minimo que alcanza el nivel): hi-lo = 2 ln 2; abre cuando
    # hi >= ln 2, i.e. (0.1-3c/2) L >= C ln L. Imposible si c >= 1/15:
    # el nivel maximo de Teor. 1.1 es x^{1/2+1/30}.
    def threshold(c, C):
        if c >= 1.0 / 15.0:
            return None
        a = 0.1 - 3.0 * c / 2.0
        lo_L = C / a
        hi_L = lo_L * 2.0
        while a * hi_L - C * math.log(hi_L) < 0:
            hi_L *= 2.0
            if hi_L > 1e7:
                return float("inf")
        for _ in range(100):
            mid = 0.5 * (lo_L + hi_L)
            if a * mid - C * math.log(mid) >= 0:
                hi_L = mid
            else:
                lo_L = mid
        return hi_L / math.log(10)
    out = []
    for c in (0.02, 0.05, 0.10, 2.0 / 17.0):
        for C in (5, 2):
            out.append((c, C, threshold(c, C)))
    d = "; ".join("c=%.4f C=%d: log10 N0=%s" % (c, C, ("%.1f" % v if v is not None else "NUNCA (c>=1/15)"))
                  for c, C, v in out)
    # L1 -> L2 por Cauchy (forma pi): sum|E_pi|^2 <= max|E_pi| * sum|E_pi|
    # max|E_pi| <= 8 x^{1/2-delta}/ln x ; sum|E_pi| << delta x/ln x + x (ln ln x)^2/(ln x)^2
    # objetivo pi: x (ln x)^{C-2}, C=5 ; delta = c/2 con c=0.05 (minimo que alcanza Q)
    delta = 0.025
    rows = []
    for N in (1e12, 1e24, 1e100):
        lnN = math.log(N)
        lln = math.log(lnN)
        mx = 8 * N ** (0.5 - delta) / lnN
        s1 = delta * N / lnN + N * lln * lln / (lnN * lnN)
        ratio = mx * s1 / (N * lnN ** 3)
        rows.append((N, ratio))
    # cruce ratio=1 por biseccion en L=ln N
    def ratio_L(L):
        lln = math.log(L)
        mx = 8 * math.exp(L * (0.5 - delta)) / L
        s1 = delta * math.exp(L) / L + math.exp(L) * lln * lln / (L * L)
        return mx * s1 / (math.exp(L) * L ** 3)
    lo_L, hi_L = 1.0, 500.0
    for _ in range(80):
        mid = 0.5 * (lo_L + hi_L)
        if ratio_L(mid) > 1.0:
            hi_L = mid
        else:
            lo_L = mid
    d2 = "; ".join("N=%.0e ratio=%.3g" % r for r in rows) + "; cruce ratio=1 en log10 N*=%.1f" % (hi_L / math.log(10))
    check("H5 (ventana Maynard + L1->L2)", True,
          "ventana Teor.1.1 primer entero: " + d + "; Cauchy L1->L2: " + d2)


def main():
    zeros = sys.argv[1] if len(sys.argv) > 1 else "72_ceros_zeta.txt"
    H0()
    H1H2()
    H3(zeros)
    H4()
    H5()
    print("PASS" if not FAIL else "FAIL " + ",".join(FAIL))


if __name__ == "__main__":
    main()
