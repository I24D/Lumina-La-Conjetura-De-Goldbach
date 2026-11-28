#!/usr/bin/env python3
"""Doc 93 - Recorte de Rankin, incidencia agregada y errata de (92.1).
Solo stdlib. Uso: python3 93_rankin_cola_incidencia_agregada.py
H0 reproduce constantes docs 87-89. H1 cota de Rankin (analitica).
H2 Rankin vs conteo directo (colas exactas). H3 incidencia agregada
vs u a escala finita (N=1e5, 1e6, Q=3). H4 rango sub-polinomico:
encendido geometrico y de balance, C0 creciente, conteo de u<=y.
H5 errata (92.1): contraejemplo q|N; forma q=N de Lumina (no basta);
codigo doc92 con/sin restriccion. H6 umbral de Cauchy corregido c<1/3.
H7 mecanismo de la implicacion GRH: ortogonalidad y desplazamiento.
"""
import math

FAIL = []


def check(name, cond, detail):
    print(("PASS " if cond else "FAIL ") + name + " | " + detail)
    if not cond:
        FAIL.append(name)


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


# ---- F y f por el sistema de retardo (maquinaria del doc 89, sin cambios) ----
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


# ---------------- H0: reproduce docs 87-89 ----------------
def H0():
    J0 = J_of_c(0.0)
    lo, hi = 0.0, 2.0 / 17.0
    for _ in range(60):
        mid = 0.5 * (lo + hi)
        if Bal(mid) > 0:
            lo = mid
        else:
            hi = mid
    cstar = 0.5 * (lo + hi)
    ok = (abs(f5_grid - 0.9982417256) < 1e-6 and abs(F4_simp - 1.021641553) < 1e-6
          and abs(c110_simp - C110) < 1e-9 and abs(J0 - 1.5238762503) < 1e-6
          and abs(Bal(0.0) - 0.061404) < 1e-5 and abs(cstar - 0.039884) < 1e-5)
    check("H0 (reproduce docs 87-89)",
          ok, "f5=%.10f F4=%.10f c110=%.10f J0=%.10f Bal0=%.6f c*=%.6f"
          % (f5_grid, F4_simp, c110_simp, J0, Bal(0.0), cstar))


# ---------------- H1: cota de Rankin, analitica ----------------
def rankin_product(primes, delta):
    # prod_{l|Q} (1 - l^{-(1-delta)})^{-1}  (suma de Euler de u^{delta-1})
    p = 1.0
    for l in primes:
        p *= 1.0 / (1.0 - l ** (-(1.0 - delta)))
    return p


def H1():
    rows = []
    for N in (1e6, 1e12, 1e24, 1e100, 1e200):
        L = math.log(N)
        delta = L ** (-0.5)
        fac = math.exp(-(L ** 0.25))  # y^{-delta} con y = exp(L^{3/4})
        pr3 = rankin_product([3], delta)
        ck3 = 3.0 / 2.0
        rows.append((N, L ** 0.25, fac, pr3, pr3 / ck3, fac * pr3,
                     L ** (-2), L ** (-5)))
    # cruce fac < L^{-A}: L^{1/4} = A ln L, biseccion en L
    def cross(A):
        lo, hi = 10.0, 1e12
        for _ in range(200):
            mid = 0.5 * (lo + hi)
            if mid ** 0.25 - A * math.log(mid) > 0:
                hi = mid
            else:
                lo = mid
        return 0.5 * (lo + hi) / math.log(10.0)
    cr = {A: cross(A) for A in (2, 5, 10)}
    d = "; ".join("N=%.0e L^1/4=%.2f y^-d=%.3g prod=%.4f prod/Ck=%.4f cota=%.3g L^-2=%.3g L^-5=%.3g"
                  % r for r in rows)
    d2 = "; ".join("A=%d: log10 N*=%.0f" % (A, cr[A]) for A in (2, 5, 10))
    # la cota analitica del producto: prod <= Ck * exp(4 delta sum ln l/(l-1))
    okprod = True
    for N in (1e12, 1e100):
        L = math.log(N)
        delta = L ** (-0.5)
        for primes in ([3], [3, 5], [2, 3, 5]):
            ck = math.prod(p / (p - 1.0) for p in primes)
            bound = ck * math.exp(4 * delta * sum(math.log(p) / (p - 1.0) for p in primes))
            if rankin_product(primes, delta) > bound * (1 + 1e-12):
                okprod = False
    check("H1 (Rankin: cola u>y)",
          okprod and rows[-1][5] < 0.02 and cr[2] > 0,
          d + " || cruces fac<L^-A: " + d2)


# ---------------- H2: Rankin vs conteo directo (colas exactas) ----------------
def gen_u(primes, U):
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


def H2():
    parts = []
    ok = True
    for primes in ([3], [3, 5]):
        k = len(primes) + 1
        total = math.prod(p / (p - 1.0) for p in primes)
        for N in (1e12, 1e24):
            L = math.log(N)
            Ly = L ** 0.75
            y = math.exp(Ly)
            delta = L ** (-0.5)
            us = gen_u(primes, int(y))
            partial = sum(1.0 / u for u in us)
            tail_exact = total - partial
            rankin = math.exp(-(L ** 0.25)) * rankin_product(primes, delta)
            # conteo directo: int_y^inf N(t)/t^2 dt con N(t)~(ln t)^{k-1}/((k-1)! prod ln l).
            # Es una COTA SUPERIOR de la cola (la cola exacta es la integral MENOS
            # el termino de frontera N(y)/y; la diferencia se cancela a escala
            # finita y no se reporta como estimacion: la cola exacta es la medida).
            den = math.factorial(k - 1) * math.prod(math.log(p) for p in primes)
            if k == 2:
                direct = (Ly + 1.0) / (y * math.log(primes[0]))
            else:
                direct = (Ly * Ly + 2 * Ly + 2.0) / (y * den)
            if not (0 < tail_exact <= direct <= rankin):
                ok = False
            parts.append("Q=%d N=%.0e y=%.3g cola_exacta=%.3g rankin=%.3g directo=%.3g directo/exacta=%.1f rankin/exacta=%.0f"
                         % (math.prod(primes), N, y, tail_exact, rankin, direct,
                            direct / tail_exact, rankin / tail_exact if tail_exact > 0 else float("nan")))
    check("H2 (Rankin vs conteo directo)", ok, " || ".join(parts))

# ---------------- criba Lambda / spf (reutilizada del doc 92) ----------------
def build_sieve(N):
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
    return is_prime, lam


def omega_rad(Xmax, is_prime):
    # omega[v] = nº de primos distintos; rad[v] = producto de primos distintos
    omega = [0] * (Xmax + 1)
    rad = [1] * (Xmax + 1)
    for p in range(2, Xmax + 1):
        if is_prime[p]:
            for m in range(p, Xmax + 1, p):
                omega[m] += 1
                rad[m] *= p
    return omega, rad


# ---------------- H3: incidencia agregada vs u (escala finita) ----------------
def H3():
    parts = []
    ok = True
    ln3 = math.log(3.0)
    for N in (100000, 1000000):
        is_prime, lam = build_sieve(N)
        Xmax = N // 3
        omega, rad = omega_rad(Xmax, is_prime)
        L = math.log(N)
        y = math.exp(L ** 0.75)
        us = []
        u = 1
        while u <= N // 3:
            us.append(u)
            u *= 3
        rows = []
        tot = tot1 = tot2 = 0.0
        for u in us:
            Xu = N / (3.0 * u)
            Ju = P1 = P2 = PP = 0.0
            for v in range(1, int(Xu)):
                if v % 3 == 0:
                    continue
                w = omega[v]
                if w not in (1, 2):
                    continue
                n = N - 3 * u * v
                if n < 2:
                    break
                if math.gcd(n, N) != 1:
                    continue
                lv = lam[n]
                if lv == 0.0:
                    continue
                c = lv * (ln3 + math.log(rad[v]))
                Ju += c
                if is_prime[v]:
                    P1 += c
                elif w == 2:
                    P2 += c
                else:
                    PP += c
            rows.append((u, Ju, P1, P2, PP))
            tot += Ju
            tot1 += P1
            tot2 += P2
        head = sum(J for u, J, _, _, _ in rows if u <= 9)
        tail_y = sum(J for u, J, _, _, _ in rows if u > y)
        ju1 = rows[0][1]
        # forma 1/u: J_u * u dentro de un factor 5 de J_1 para u <= 27
        shape = [Ju * u / ju1 for u, Ju, _, _, _ in rows if u <= 27 and Ju > 0]
        shape_ok = all(0.2 < s < 5.0 for s in shape) if shape else False
        p2frac = tot2 / tot if tot > 0 else float("nan")
        if not (head / tot > 0.75 and p2frac > 0.5 and shape_ok):
            ok = False
        parts.append("N=%d y=%.0f Jtot=%.1f frac(u=1)=%.3f frac(u<=9)=%.3f frac(u>y)=%.4f "
                     "P1=%.3f P2=%.3f PP=%.3f; J_u*u/J_1(u<=27)=%s"
                     % (N, y, tot, ju1 / tot, head / tot, tail_y / tot if tot else 0.0,
                        tot1 / tot, p2frac, sum(r[4] for r in rows) / tot,
                        ["%.2f" % s for s in shape]))
    check("H3 (incidencia agregada vs u)", ok, " || ".join(parts))


# ---------------- H4: rango sub-polinomico ----------------
def H4():
    # balance con c_eff = (ln N)^(-1/4) (u <= y); biseccion del cero.
    # Bal(c) solo esta definido para c < 2/17 (geometria); la biseccion
    # empieza por encima del encendido geometrico (L > 9683).
    def bal_ceff(L):
        return Bal(L ** (-0.25))
    lo, hi = 20000.0, 1e7
    for _ in range(200):
        mid = 0.5 * (lo + hi)
        if bal_ceff(mid) > 0:
            hi = mid
        else:
            lo = mid
    Lstar = 0.5 * (lo + hi)
    # encendido geometrico: (1/15)L - (17/30)logA - B ln L = 0, logA = L^{3/4}+ln3
    def g_geo(L, B=10.0):
        logA = L ** 0.75 + math.log(3.0)
        return (1.0 / 15.0) * L - (17.0 / 30.0) * logA - B * math.log(L)
    lo2, hi2 = 10.0, 1e7
    for _ in range(200):
        mid = 0.5 * (lo2 + hi2)
        if g_geo(mid) > 0:
            hi2 = mid
        else:
            lo2 = mid
    Lgeo = 0.5 * (lo2 + hi2)
    rows = []
    for L in (math.log(1e12), math.log(1e24), math.log(1e100),
              5000 * math.log(10.0), 180000 * math.log(10.0)):
        ceff = L ** (-0.25)
        C0 = 10.0 + (L ** 0.75) / (2.0 * math.log(L))
        cnt2 = (L ** 0.75) / math.log(3.0)  # #(u<=y) asintotico, k=2
        err_ratio = cnt2 / (L ** 9)  # #(u<=y)/(ln N)^{A'-1}, A'=10
        geo = g_geo(L) > 0
        bal = Bal(ceff) if ceff < 2.0 / 17.0 else float("nan")
        bal3 = Bal(ceff, 2.0 / 3.0) if ceff < 2.0 / 17.0 else float("nan")
        rows.append((L / math.log(10.0), ceff, geo, bal, bal3, C0, cnt2, err_ratio))
    d = ("cero balance c_eff: log10 N0=%.0f (L*=%.0f); encendido geometrico: log10 N0=%.0f (L=%.0f)"
         % (Lstar / math.log(10.0), Lstar, Lgeo / math.log(10.0), Lgeo))
    d += " || " + "; ".join(
        "log10N=%.0f c_eff=%.4f geom=%s Bal=%s Bal(Q=3)=%s C0=%.2f #(u<=y)~%.1f err_ratio=%.3g"
        % (r[0], r[1], r[2],
           ("%.4f" % r[3]) if r[3] == r[3] else "n/a(sin geometria)",
           ("%.4f" % r[4]) if r[4] == r[4] else "n/a",
           r[5], r[6], r[7]) for r in rows)
    check("H4 (rango sub-polinomico u<=y)",
          Lgeo < Lstar and bal_ceff(180000 * math.log(10.0)) > 0
          and all(r[7] < 1e-6 for r in rows),
          d)


# ---------------- H5: errata (92.1) ----------------
def H5():
    # (a) contraejemplo q|N: E(N;2,N) y otros q|N a N=1e5, 1e6
    out = []
    ok = True
    for N in (100000, 1000000):
        _, lam = build_sieve(N)
        psiN = sum(lam)
        lnN = math.log(N)
        target = N * lnN ** 5
        for q in (2, 4, 5, 10, 25):
            if N % q != 0:
                continue
            # psi(N;q,0): multiplos de q (solo potencias de primos que dividen q)
            psi_q = sum(lam[n] for n in range(q, N + 1, q))
            phi_q = q
            p = q
            for pr in (2, 5):
                if p % pr == 0:
                    phi_q = phi_q // pr * (pr - 1)
                    while p % pr == 0:
                        p //= pr
            e = psi_q - psiN / phi_q
            out.append("N=%d q=%d |E|^2=%.3g objetivo=%.3g ratio=%.3g"
                       % (N, q, e * e, target, e * e / target))
            if q == 2 and N == 1000000 and e * e <= target:
                ok = False
    # (b) doc 92 H1 con y sin la restriccion (q,N)=1 (solo clase fija a=N)
    _, lam0 = build_sieve(100000)
    psi0 = sum(lam0)
    phi0 = [0] * 3000
    for i in range(3000):
        phi0[i] = i
    for p in range(2, 3000):
        if phi0[p] == p:
            for m in range(p, 3000, p):
                phi0[m] -= phi0[m] // p
    N0 = 100000
    Q0 = int(N0 ** (0.5 + 0.05 / 2))
    s2_res = s2_nores = 0.0
    for q in range(2, Q0 + 1):
        psi_q = sum(lam0[N0 % q::q])
        e = psi_q - psi0 / phi0[q]
        if math.gcd(q, N0) == 1:
            s2_res += e * e
        s2_nores += e * e
    out.append("H1-doc92 N=1e5: S2fix CON restriccion=%.4g (doc92: 3.119e5) SIN=%.4g ratio=%.3g"
               % (s2_res, s2_nores, s2_nores / s2_res))
    if not (0.99 < s2_res / 3.119e5 < 1.01 and s2_nores / s2_res > 100):
        ok = False
    # (c) forma q=N de Lumina: N0 primo, E = Lambda(N0) - psi(N0)/phi(N0)
    Np = 999983
    _, lamp = build_sieve(Np)
    psip = sum(lamp)
    e = math.log(Np) - psip / (Np - 1)
    lnNp = math.log(Np)
    out.append("q=N0 primo N0=%d: E=%.4f |E|^2=%.1f vs N0(lnN0)^5=%.3g ratio=%.3g (NO rompe la suma)"
               % (Np, e, e * e, Np * lnNp ** 5, e * e / (Np * lnNp ** 5)))
    if not (abs(e) < math.log(Np) and e * e < Np * lnNp ** 5):
        ok = False
    check("H5 (errata 92.1: (q,N)=1 necesaria)", ok, " || ".join(out))


# ---------------- H6: umbral de Cauchy corregido c<1/3 ----------------
def H6():
    C = 5
    rows = []
    for N in (1e12, 1e24, 1e100, 1e200):
        lnN = math.log(N)
        R = N ** (3 * 0.05 / 4 - 0.25) * lnN ** (C / 2 + 2)
        Rold = N ** (5 * 0.05 / 4 - 0.25) * lnN ** (C / 2 + 2)
        rows.append((N, math.log10(R), math.log10(Rold)))
    # cruces: (3c/4-1/4) L + 4.5 ln L = 0 (corregido) vs (5c/4-1/4) L + 4.5 ln L = 0
    def cross(c, coef):
        a = coef * c - 0.25
        if a >= 0:
            return None
        lo, hi = 10.0, 1e6
        for _ in range(200):
            mid = 0.5 * (lo + hi)
            if a * mid + 4.5 * math.log(mid) > 0:
                lo = mid
            else:
                hi = mid
        return 0.5 * (lo + hi) / math.log(10.0)
    cr = []
    for c in (0.05, 0.10, 2.0 / 17.0):
        cr.append((c, cross(c, 0.75), cross(c, 1.25)))
    d = "; ".join("N=%.0e log10 R_corr=%.2f (doc92: %.2f)" % r for r in rows)
    d += " || cruces log10 N* corr/doc92: " + "; ".join(
        "c=%.4f: %.1f / %.1f" % r for r in cr)
    d += " || D_A/Q=N^-c a N=1e12 c=0.05: %.3f; c*=%.6f<1/3; 2/17=%.4f<1/3" % (
        1e12 ** (-0.05), 0.039884, 2.0 / 17.0)
    ok = all(v is not None and v > 0 for _, v, _ in cr) and 0.039884 < 1.0 / 3.0
    check("H6 (Cauchy corregido: umbral c<1/3)", ok, d)


# ---------------- H7: mecanismo de la implicacion GRH ----------------
def H7():
    out = []
    ok = True
    for N in (100000, 1000000):
        _, lam = build_sieve(N)
        psiN = sum(lam)
        # r=3: chi no principal, chi(1)=1, chi(2)=-1
        psi_chi_direct = sum(lam[n] * (1 if n % 3 == 1 else (-1 if n % 3 == 2 else 0))
                             for n in range(2, N + 1) if lam[n] > 0)
        b1 = sum(lam[n] for n in range(1, N + 1, 3) if n % 3 == 1)
        b2 = sum(lam[n] for n in range(2, N + 1, 3) if n % 3 == 2)
        e1 = b1 - psiN / 2.0
        e2 = b2 - psiN / 2.0
        via_e = e1 - e2
        # r=5: chi primitivo con chi(2)=i
        chi5 = {0: 0j, 1: 1 + 0j, 2: 0 + 1j, 3: 0 - 1j, 4: -1 + 0j}
        psi_chi5 = sum(lam[n] * chi5[n % 5] for n in range(2, N + 1) if lam[n] > 0)
        buckets = [sum(lam[n] for n in range(a, N + 1, 5)) for a in range(5)]
        via_e5 = sum(chi5[a] * (buckets[a] - psiN / 4.0) for a in range(1, 5))
        d3 = abs(psi_chi_direct - via_e)
        d5 = abs(psi_chi5 - via_e5)
        out.append("N=%d r=3: directo=%.4f viaE=%.4f dif=%.2g; r=5: dif=%.2g"
                   % (N, psi_chi_direct, via_e, d3, d5))
        if not (d3 < 1e-6 and d5 < 1e-6):
            ok = False
    # desplazamiento: |psi(x;r,a)-psi(N';r,a)| <= (#n en (x,N'] congr a) * ln N'
    N = 100000
    _, lam = build_sieve(N + 100)
    for r, a in ((3, 1), (3, 2), (5, 3)):
        x = N
        n2 = x + 1
        while n2 % 2 != 0 or n2 % r != a % r:
            n2 += 1
        lhs = abs(sum(lam[n] for n in range(x + 1, n2 + 1) if n % r == a % r))
        cnt = len([n for n in range(x + 1, n2 + 1) if n % r == a % r])
        rhs = cnt * math.log(n2)
        out.append("despl r=%d a=%d: |dpsi|=%.3f <= %d*lnN'=%.3f" % (r, a, lhs, cnt, rhs))
        if not (lhs <= rhs + 1e-9):
            ok = False
    check("H7 (mecanismo GRH: ortogonalidad+desplazamiento)", ok, " || ".join(out))


def main():
    H0()
    H1()
    H2()
    H3()
    H4()
    H5()
    H6()
    H7()
    print("PASS" if not FAIL else "FAIL " + ",".join(FAIL))


if __name__ == "__main__":
    main()
