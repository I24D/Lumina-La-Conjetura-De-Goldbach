#!/usr/bin/env python3
"""58_nivel_distribucion_peso_relajado.py -- Doc 58.
Nivel de distribucion del peso relajado c_t (familia A del doc 56:
c0 = atomos psi0 del doc 54, c1 = constante 1).

Definiciones (docs 38, 46):
  E_d(N) = A_d(N) - N/phi(d), A_d(N) = sum_{n<N, n=N (mod d)} Lambda(n).
  Error pesado: E_t(N) = sum_{i,j} c_t(i,j) * E_{q_i q_j}(N).
  Termino principal: M = N * sum c_t/phi(qr)  (M1 = sum c_t/phi(qr)).
  Nivel theta: el mayor theta tal que el error pesado con qr <= N^theta
  es despreciable frente a M.

Cotas que se miden (exactas donde es barato, float+fsum a h grande
con verificacion exacta a h=12,24):
  R_BV(A)  = 2*max|c_t| / ((ln N)^A * M1): cota BV (doc 38) del error
             con qr <= N^(1/2), dividida por M1 (la N se cancela).
  R_triv   = ln(N) * sum_{qr>N^(1/2)} |c_t|/(qr) / M1: cota trivial
             del error con qr > N^(1/2), dividida por M1.
  Si R_triv << 1, la parte por encima del nivel BV la controla la cota
  trivial; si R_triv >> 1, hace falta dispersion/BFI.
  theta_eff(B,N) = 1/2 - B*ln(ln N)/ln(N): nivel efectivo finito-N
  con la perdida (log N)^B de BV.
  maxfac = max_d #{(i,j): q_i q_j = d}: en la malla hay primos
  representativos degenerados (q_i = q_{i+1}), asi que maxfac > 2;
  se mide y se usa en la cota BV (no se supone <= 2).

Estructura (Proposicion 58.1): c_t es suma de 3 terminos de rango 1
con factores soportados en primos de [N^(1/10), N^(1/2)]; el peso
inducido lambda_t(d) = sum_{qr=d} c_t(q,r) NO es well-factorable de
ningun nivel Q <= N (ver doc 58 para la prueba; el script verifica
los hechos de soporte: min/max primo y c_t(esquina)>0).
"""
from fractions import Fraction
import bisect
import math
from collections import Counter

CHEN = Fraction(6706, 10000)  # 0.6706 (doc 52)


def primes_upto(n):
    bs = bytearray(b"\x01") * (n + 1)
    bs[0:2] = b"\x00\x00"
    for i in range(2, int(n ** 0.5) + 1):
        if bs[i]:
            bs[i * i:n + 1:i] = b"\x00" * ((n - i * i) // i + 1)
    return [i for i in range(n + 1) if bs[i]]


PRIMES = primes_upto(10 ** 6 + 5000)


def rep_prime(alpha, E=12):
    tgt = 10 ** (E * float(alpha))
    return PRIMES[bisect.bisect_left(PRIMES, math.ceil(tgt))]


def mesh(v, h):
    lo = Fraction(1, v)
    hi = Fraction(1, 2)
    step = (hi - lo) / (h + 1)
    return [lo + (t + 1) * step for t in range(h)], step


def psi0(a, b):
    return (a * b) ** 2 + (a * b) ** 3


def build_atoms(alpha, psi):
    h = len(alpha)
    c = [[Fraction(0)] * h for _ in range(h)]
    for i in range(h - 1):
        for j in range(h - 1):
            w = psi(alpha[i], alpha[j])
            c[i][j] += w
            c[i][j + 1] -= w
            c[i + 1][j] -= w
            c[i + 1][j + 1] += w
    return c


def l1(M):
    return sum(abs(M[i][j]) for i in range(len(M)) for j in range(len(M)))


def normalize(M):
    m = l1(M)
    assert m > 0
    return [[M[i][j] / m for j in range(len(M))] for i in range(len(M))], m


def build_ct(alpha, t):
    """Familia A del doc 56: c_t = (1-t)*c0/||c0||_1 + t*const/||const||_1."""
    h = len(alpha)
    c0, m0 = normalize(build_atoms(alpha, psi0))
    c1 = [[Fraction(1, h * h)] * h for _ in range(h)]
    ct = [[(1 - t) * c0[i][j] + t * c1[i][j]
           for j in range(h)] for i in range(h)]
    return ct, m0


def phi2(q, r):
    return q * (q - 1) if q == r else (q - 1) * (r - 1)


def rank1_factors(alpha):
    """u1,v1,u2,v2 de la descomposicion c0 = u1(x)v1 + u2(x)v2.
    F_k(i) = alpha_i^k si 0<=i<=h-2, 0 si i=h-1 o i<0;
    u_k(i) = F_k(i)-F_k(i-1). Exacto (Fraction)."""
    h = len(alpha)

    def F(k, i):
        return Fraction(0) if (i < 0 or i > h - 2) else alpha[i] ** k
    u1 = [F(2, i) - F(2, i - 1) for i in range(h)]
    u2 = [F(3, i) - F(3, i - 1) for i in range(h)]
    return u1, u2


def fmt(x):
    return f"{float(x):.4e}"


def analyze(v, h, E, t):
    alpha, step = mesh(v, h)
    N = 10 ** E
    lnN = math.log(N)
    q = [rep_prime(a, E) for a in alpha]
    ct, m0 = build_ct(alpha, t)
    # --- exacto: normas, max, sumas, masas por alpha+beta<=1/2 ---
    M1n = l1(ct)
    # ||c_t||_1 < 1 por cancelacion parcial entre la parte de atomos
    # (con entradas negativas) y la constante positiva (doc 56: la
    # correccion es <=3% a h=96). No es 1 exacto; se registra y se usa
    # como denominador de las fracciones de masa.
    assert Fraction(9, 10) < M1n <= 1, ("||c_t||_1 fuera de rango", M1n)
    T = sum(ct[i][j] for i in range(h) for j in range(h))
    assert T == t, ("sum c_t != t", v, h, t)
    mx = max(abs(ct[i][j]) for i in range(h) for j in range(h))
    half = Fraction(1, 2)
    mass_below = sum(abs(ct[i][j]) for i in range(h) for j in range(h)
                     if alpha[i] + alpha[j] <= half)
    mass_above = M1n - mass_below
    T_below = sum(ct[i][j] for i in range(h) for j in range(h)
                  if alpha[i] + alpha[j] <= half)
    T_above = T - T_below
    # --- rango 1 exacto: c_t = (1-t)/m0*(u1xv1+u2xv2) + t*(1/h)x(1/h) ---
    u1, u2 = rank1_factors(alpha)
    c0 = build_atoms(alpha, psi0)
    for i in range(h):
        for j in range(h):
            assert c0[i][j] == u1[i] * u1[j] + u2[i] * u2[j], (i, j)
            assert ct[i][j] == ((1 - t) * c0[i][j] / m0
                                + t * Fraction(1, h * h)), (i, j)
    # --- factorizaciones del producto qr ---
    cnt = Counter(q[i] * q[j] for i in range(h) for j in range(h))
    maxfac = max(cnt.values())
    # maxfac crece a N pequeno por primos representativos degenerados
    # (muchos alpha comparten un mismo primo); se mide, no se supone.
    assert maxfac <= 200, ("factorizaciones fuera de rango", maxfac)
    # --- soporte: min/max primo; testigo de la obstruccion 58.1 ---
    assert ct[h - 1][h - 1] > 0, "esquina (h-1,h-1) no positiva"
    assert min(q) >= N ** 0.1 - 1 and max(q) <= N ** 0.5 + 1
    # --- M1 y S_abs: exactos a h<=24, float+fsum a h mayor ---
    def pair_sums(exact):
        M1b = M1a = S_b = S_a = 0
        if exact:
            M1b = Fraction(0)
            M1a = Fraction(0)
            S_b = Fraction(0)
            S_a = Fraction(0)
        else:
            M1b = []
            M1a = []
            S_b = []
            S_a = []
        for i in range(h):
            for j in range(h):
                pr = q[i] * q[j]
                below = (pr * pr <= N)
                if exact:
                    if below:
                        M1b += ct[i][j] / phi2(q[i], q[j])
                        S_b += abs(ct[i][j]) / pr
                    else:
                        M1a += ct[i][j] / phi2(q[i], q[j])
                        S_a += abs(ct[i][j]) / pr
                else:
                    if below:
                        M1b.append(float(ct[i][j]) / phi2(q[i], q[j]))
                        S_b.append(float(abs(ct[i][j])) / pr)
                    else:
                        M1a.append(float(ct[i][j]) / phi2(q[i], q[j]))
                        S_a.append(float(abs(ct[i][j])) / pr)
        if exact:
            return M1b, M1a, S_b, S_a
        return (math.fsum(M1b), math.fsum(M1a),
                math.fsum(S_b), math.fsum(S_a))
    exact = (h <= 24)
    M1b, M1a, S_b, S_a = pair_sums(exact)
    M1 = M1b + M1a
    assert M1 > 0, ("M1 <= 0", v, h, t)
    R_triv = lnN * float(S_a) / float(M1)
    # Cota BV: sum_{qr<=D} |c_t| |E| <= maxfac*max|c_t| * sum_{d<=D}|E_d|.
    R_bv5 = maxfac * float(mx) / ((lnN ** 5) * float(M1))
    R_bv10 = maxfac * float(mx) / ((lnN ** 10) * float(M1))
    # --- nivel efectivo finito-N ---
    th = {}
    cov = {}
    for B in (1, 2, 3, 5):
        te = 0.5 - B * math.log(lnN) / lnN
        th[B] = te
        cov[B] = float(sum(abs(ct[i][j]) for i in range(h) for j in range(h)
                           if float(alpha[i] + alpha[j]) <= te) / M1n) \
            if te > 0 else 0.0
    return dict(v=v, h=h, E=E, t=t, M1=M1, M1b=M1b, M1a=M1a,
                S_b=S_b, S_a=S_a, mx=mx, M1n=M1n,
                mass_below=mass_below / M1n,
                mass_above=mass_above / M1n,
                T_below=T_below, T_above=T_above,
                R_triv=R_triv, R_bv5=R_bv5, R_bv10=R_bv10,
                maxfac=maxfac, qmin=min(q), qmax=max(q),
                th=th, cov=cov, exact=exact)


def main():
    print("== cross-check float vs exacto (h=12,24; v=10, E=12, t=0.05) ==")
    for h in (12, 24):
        alpha, _ = mesh(10, h)
        q = [rep_prime(a, 12) for a in alpha]
        ct, _ = build_ct(alpha, Fraction(1, 20))
        M1e = sum((ct[i][j] / phi2(q[i], q[j]))
                  for i in range(h) for j in range(h))
        M1f = math.fsum(float(ct[i][j]) / phi2(q[i], q[j])
                        for i in range(h) for j in range(h))
        assert abs(float(M1e) - M1f) / float(M1e) < 1e-12, (h, M1e, M1f)
        print(f"h={h}: M1 exacto={float(M1e):.6e} float={M1f:.6e} OK")
    print()
    print("v h E t | M1 | M1_abajo/M1 | |c|_abajo | sumc_abajo/t | "
          "R_triv | R_BV(5) | R_BV(10) | max|c_t| | fac | "
          "theta_eff/cobertura(B=1,2,3,5)")
    rows = []
    for (v, hs, E) in [(10, [12, 24, 48, 96], 12), (8, [48], 12),
                       (10, [96], 7)]:
        for h in hs:
            for t in (Fraction(1, 20), Fraction(1, 10)):
                r = analyze(v, h, E, t)
                rows.append(r)
                print(f"v={v} h={h} E={E} t={float(t):.2f}: "
                      f"M1={fmt(r['M1'])} "
                      f"M1ab/M1={float(r['M1b']) / float(r['M1']):.3f} "
                      f"|c|ab={float(r['mass_below']):.3f} "
                      f"sumcab/t={float(r['T_below']) / float(t):.3f} "
                      f"R_triv={r['R_triv']:.4f} "
                      f"R_BV5={r['R_bv5']:.2e} R_BV10={r['R_bv10']:.2e} "
                      f"max|c|={fmt(r['mx'])} fac={r['maxfac']} "
                      f"th={ {B: round(r['th'][B], 3) for B in (1, 2, 3, 5)} } "
                      f"cov={ {B: round(r['cov'][B], 3) for B in (1, 2, 3, 5)} }")
    print()
    print("PASS: nivel medido; rango 1 exacto (3 terminos) en todos los "
          "casos; sum c_t=t exacto; ||c_t||_1<1 registrado; maxfac medido; "
          "esquina positiva; float==exacto a h<=24.")


if __name__ == "__main__":
    main()
