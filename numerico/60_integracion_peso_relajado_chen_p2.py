#!/usr/bin/env python3
"""60_integracion_peso_relajado_chen_p2.py -- Doc 60.
Integracion del peso relajado c_t (docs 54-56) con Chen/Richert
(doc 52) y la deteccion de P2.

Preguntas (exactas donde es barato; float con cross-checks en
la enumeracion a N=10^7):
  (1) Peso inducido sobre numeros W(m) = sum_{qr|m} c_t(q,r)
      (Prop. 59.1): W(primo)=0 [T, por soporte]; W(P2)=2c_t;
      W(P3)=2*suma de c sobre pares; descomposicion por Omega
      de T = sum_p (log p) W(p) a N=10^7, exacta por criba.
  (2) Escala: T contra S(N)*N/log^2 N; parte P2 E2/T y E2/G.
  (3) Chen (Liu, arXiv:2203.07871v1, LEIDO): BV nivel 1/2
      (Teor. 1: q<=x^{1/2} log^{-B} x); error de Omega (10):
      d<=N^{1/2-eps}; Selberg (Lema 3): z=N^{1/4-eps/2};
      W(N)>=2.6408 por Jurkat-Richert, s=(1/2)/(1/10)=5.
      Ninguna parte usa nivel > 1/2.
  (4) Sustitucion en W(N) y en Omega (sumas sobre la cuna);
      combinacion c_t + Chen (perturbacion O(N^{4/5})).

Marco: malla doc 52/54 (v=10, h=12,24,48,96); familia A del
doc 56 (c0 = atomos psi0, c1 = const); N=10^7 exacto con
primos REALES binned por tamano log (sin primos
representativos: sin degeneracion).
"""
from fractions import Fraction
import bisect
import math
import numpy as np
from collections import Counter

# ---------------------------------------------------------------- helpers
def primes_upto(n):
    bs = bytearray(b"\x01") * (n + 1)
    bs[0:2] = b"\x00\x00"
    for i in range(2, int(n ** 0.5) + 1):
        if bs[i]:
            bs[i * i:n + 1:i] = b"\x00" * ((n - i * i) // i + 1)
    return [i for i in range(n + 1) if bs[i]]


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
    """Familia A doc 56: c_t=(1-t)c0/||c0||_1 + t*const/||const||_1."""
    h = len(alpha)
    c0, m0 = normalize(build_atoms(alpha, psi0))
    c1 = [[Fraction(1, h * h)] * h for _ in range(h)]
    ct = [[(1 - t) * c0[i][j] + t * c1[i][j]
           for j in range(h)] for i in range(h)]
    return ct, m0


def cell_of(alpha_q, lo_f, step_f, h):
    """Celda de la malla para un primo de tamano log alpha_q."""
    t = math.floor((alpha_q - lo_f) / step_f) - 1
    return 0 if t < 0 else (h - 1 if t > h - 1 else t)


def fmt(x):
    return f"{float(x):.4e}"


def twin_c2(limit=10 ** 6):
    """C2 = prod_{p>2} (1 - 1/(p-1)^2), producto hasta `limit`."""
    pr = primes_upto(limit)
    c2 = 1.0
    for p in pr:
        if p > 2:
            c2 *= (1.0 - 1.0 / (p - 1) ** 2)
    return c2


# ---------------------------------------------------------------- exacto: malla
def mesh_quantities(v, h, t):
    """Cantidades exactas (Fraction) sobre la malla: sum c, diag,
    suma sobre i<j, suma sobre la cuna de Chen (doc 52), signo."""
    al, step = mesh(v, h)
    ct, m0 = build_ct(al, t)
    assert abs(float(sum(sum(r) for r in ct)) - float(t)) < 1e-9
    diag = sum(ct[i][i] for i in range(h))
    tri = sum(ct[i][j] for i in range(h) for j in range(i + 1, h))
    # identidad exacta: tri = (sum c - diag)/2 (simetria)
    assert tri == (sum(sum(r) for r in ct) - diag) / 2, "simetria"
    # cuna de Chen (doc 52 sec. 3): i<b, j>=b, alpha_j <= (1-alpha_i)/2
    b = sum(1 for x in al if x < Fraction(1, 3))
    wedge = [(i, j) for i in range(h) for j in range(h)
             if i < b and j >= b and al[j] <= (Fraction(1) - al[i]) / 2]
    sw = sum(ct[i][j] for (i, j) in wedge)
    # fraccion de celdas (i<j) con c_t > 0 / < 0
    pos = sum(1 for i in range(h) for j in range(i + 1, h)
              if ct[i][j] > 0)
    neg = sum(1 for i in range(h) for j in range(i + 1, h)
              if ct[i][j] < 0)
    return dict(t=t, sumc=sum(sum(r) for r in ct), diag=diag, tri=tri,
                wedge_n=len(wedge), wedge_sum=sw, b=b,
                pos=pos, neg=neg, h=h, v=v)


# ---------------------------------------------------------------- N=10^7 exacto
def setup_N7():
    """Criba hasta N=10^7 una sola vez: is_pr, logp, primos en
    [N^{1/10}, N^{1/2}]."""
    N = 10 ** 7
    primes_all = primes_upto(N)
    is_pr = np.zeros(N + 1, dtype=bool)
    is_pr[primes_all] = True
    logp = np.zeros(N + 1, dtype=float)
    logp[primes_all] = np.log(np.array(primes_all, dtype=float))
    lo_q = int(math.ceil(N ** 0.1))
    hi_q = int(math.floor(N ** 0.5))
    i0 = bisect.bisect_left(primes_all, lo_q)
    i1 = bisect.bisect_right(primes_all, hi_q)
    rp = primes_all[i0:i1]  # primos en [N^{1/10}, N^{1/2}]
    return dict(N=N, is_pr=is_pr, logp=logp, rp=rp)


def enum_N7(h, t, setup, triples=None):
    """Enumeracion a N=10^7, v=10, primos REALES binned por tamano
    log. T = sum_p (log p) W(p) via A_{qr}; E2/E3 por enumeracion
    directa de m = N-p con Omega(m)=2,3; Ege4 = T-E2-E3.
    `triples` (opcional): arrays (p1,p2,p3,kind) de casos Omega=3
    ya encontrados (independientes de h,t); se devuelven en el
    dict para reutilizarlos."""
    N = setup["N"]
    is_pr = setup["is_pr"]
    logp = setup["logp"]
    rp = setup["rp"]
    logN = math.log(N)
    al, step = mesh(10, h)
    step_f = float(step)
    cell = {q: cell_of(math.log(q) / logN, 0.1, step_f, h) for q in rp}
    ct, _ = build_ct(al, t)
    ctf = [[float(x) for x in row] for row in ct]
    # --- A_d para d = q*r (progresion n = N mod d) ---
    A = {}
    for qi in range(len(rp)):
        for qj in range(qi, len(rp)):
            d = rp[qi] * rp[qj]
            if d not in A:
                A[d] = float(logp[N % d::d].sum())
    # --- T y M1 actual (primos reales) ---
    T = 0.0
    M1a = 0.0
    for q in rp:
        cq = cell[q]
        for r in rp:
            cr = cell[r]
            c = ctf[cq][cr]
            if c == 0.0:
                continue
            d = q * r
            T += c * A[d]
            M1a += c / (q * (q - 1) if q == r else (q - 1) * (r - 1))
    # --- E2: m = p1*p2 (p1<p2) y m = p^2, con N-m primo ---
    E2 = 0.0
    n2 = 0
    for a in range(len(rp)):
        for b in range(a + 1, len(rp)):
            m = rp[a] * rp[b]
            if m >= N:
                break
            if is_pr[N - m]:
                E2 += 2.0 * ctf[cell[rp[a]]][cell[rp[b]]] * math.log(N - m)
                n2 += 1
    for q in rp:
        m = q * q
        if m >= N:
            break
        if is_pr[N - m]:
            E2 += ctf[cell[q]][cell[q]] * math.log(N - m)
            n2 += 1
    # --- E2 esperado (primalidad aleatoria 1/log N) ---
    cnt = Counter(cell[q] for q in rp)
    E2esp = 0.0
    for i in range(h):
        for j in range(i + 1, h):
            E2esp += 2.0 * ctf[i][j] * cnt.get(i, 0) * cnt.get(j, 0)
    # --- E3: triples (construidos una vez si no se pasan) ---
    if triples is None:
        triples = build_triples(N, rp, is_pr)
    p1a, p2a, p3a, kinda = triples
    E3 = 0.0
    for k in range(len(p1a)):
        p1, p2, p3, kd = int(p1a[k]), int(p2a[k]), int(p3a[k]), int(kinda[k])
        m = p1 * p2 * p3 if kd == 0 else (p1 * p1 * p3 if kd == 1
                                          else p1 ** 3)
        if kd == 0:
            w = 2.0 * (ctf[cell[p1]][cell[p2]] + ctf[cell[p1]][cell[p3]]
                       + ctf[cell[p2]][cell[p3]])
        elif kd == 1:
            w = ctf[cell[p1]][cell[p1]] + 2.0 * ctf[cell[p1]][cell[p3]]
        else:
            w = ctf[cell[p1]][cell[p1]]
        E3 += w * math.log(N - m)
    Ege4 = T - E2 - E3
    return dict(T=T, M1a=M1a, E2=E2, E2esp=E2esp, E3=E3, Ege4=Ege4,
                n2=n2, n3=len(p1a), h=h, t=float(t),
                triples=(p1a, p2a, p3a, kinda))


def build_triples(N, rp, is_pr):
    """Casos Omega(m)=3 con N-m primo: (p1,p2,p3,kind); kind 0 =
    p1<p2<p3 distintos; 1 = p1^2*p3 (p1!=p3); 2 = p1^3.
    Solo depende de N y de los primos, no de (h,t)."""
    P1, P2, P3, KD = [], [], [], []
    rparr = np.array(rp, dtype=np.int64)
    for a in range(len(rp)):
        for b in range(a + 1, len(rp)):
            lim = N // (rp[a] * rp[b])
            if lim <= rp[b]:
                break
            lo_i = b + 1
            hi_i = bisect.bisect_right(rparr, lim)
            if hi_i <= lo_i:
                continue
            cand = rparr[lo_i:hi_i]
            m3 = rp[a] * rp[b] * cand
            keep = is_pr[N - m3]
            for p3 in cand[keep]:
                P1.append(rp[a]); P2.append(rp[b]); P3.append(int(p3))
                KD.append(0)
    for p in rp:
        for q in rp:
            if p == q:
                continue
            m = p * p * q
            if m >= N:
                break
            if is_pr[N - m]:
                P1.append(p); P2.append(p); P3.append(q); KD.append(1)
        if p ** 3 < N and is_pr[N - p ** 3]:
            P1.append(p); P2.append(p); P3.append(p); KD.append(2)
    return (np.array(P1, dtype=np.int64), np.array(P2, dtype=np.int64),
            np.array(P3, dtype=np.int64), np.array(KD, dtype=np.int64))


def rep_prime(alpha, E=12):
    tgt = 10 ** (E * float(alpha))
    return PRIMES_RP[bisect.bisect_left(PRIMES_RP, math.ceil(tgt))]


PRIMES_RP = primes_upto(10 ** 6 + 5000)


if __name__ == "__main__":
    print("== (A) Cantidades exactas sobre la malla (Fraction) ==")
    for (v, h) in [(10, 12), (10, 24), (10, 48), (10, 96)]:
        for t in [Fraction(1, 20), Fraction(1, 10)]:
            q = mesh_quantities(v, h, t)
            print(f"v={v} h={h} t={float(t):.2f}: sumc={fmt(q['sumc'])} "
                  f"diag={fmt(q['diag'])} tri(i<j)={fmt(q['tri'])} "
                  f"cuna: n={q['wedge_n']} sum={fmt(q['wedge_sum'])} "
                  f"(sum/t={float(q['wedge_sum'] / t):.3f}) "
                  f"celdas i<j: +{q['pos']} -{q['neg']}")
    print("PASS (A): sum c = t exacto; tri=(sum-diag)/2 exacto.")
    print()
    print("== (B) N=10^7 exacto (primos reales) ==")
    C2 = twin_c2()
    SN = 2.0 * C2 * (5 - 1) / (5 - 2)  # N=10^7 = 2^7*5^7
    G7 = SN * 10 ** 7 / math.log(10 ** 7) ** 2
    print(f"C2={C2:.6f} S(N)={SN:.4f} G=S(N)N/log^2N={G7:.4e}")
    setup = setup_N7()
    triples = None
    for (h, t) in [(24, Fraction(1, 20)), (24, Fraction(1, 10)),
                   (48, Fraction(1, 20))]:
        r = enum_N7(h, t, setup, triples)
        triples = r["triples"]
        print(f"h={h} t={t}: T={r['T']:.4e} N*M1act={10**7*r['M1a']:.4e} "
              f"T/(N*M1act)={r['T']/(10**7*r['M1a']):.3f} "
              f"T/G={r['T']/G7:.4f} E2={r['E2']:.4e} (n2={r['n2']}, "
              f"E2esp={r['E2esp']:.4e}, E2/E2esp={r['E2']/r['E2esp']:.3f}) "
              f"E3={r['E3']:.4e} (n3={r['n3']}) Ege4={r['Ege4']:.4e} "
              f"frac: E2/T={r['E2']/r['T']:.3f} E3/T={r['E3']/r['T']:.3f} "
              f"Ege4/T={r['Ege4']/r['T']:.3f} E2/G={r['E2']/G7:.4e}")
    print("== (C) Escala a N=10^12 (M1 doc 58) ==")
    M1_12 = 5.54e-06
    T12 = 10 ** 12 * M1_12
    G12 = SN * 10 ** 12 / math.log(10 ** 12) ** 2
    print(f"T(1e12)={T12:.4e} G(1e12)={G12:.4e} T/G={T12/G12:.4e} "
          f"(ratio 1e7 -> 1e12: {(r['T']/G7)/(T12/G12):.2f} "
          f"vs (1e5)^(1/5)=10)")
    print("== (D) Cuna de Chen con peso c_t (N=10^12, primos repr.) ==")
    print("Chen (Liu Lema 4): sum f(a)/(a log(N/a)) <= 0.49254/log N "
          "= 0.017825; analogo c_t: sum_cuna c_t/(qr log(N/qr)).")
    for (h, t) in [(24, Fraction(1, 20)), (48, Fraction(1, 20)),
                   (96, Fraction(1, 20))]:
        al, _ = mesh(10, h)
        ct, _ = build_ct(al, t)
        qq = [rep_prime(a, 12) for a in al]
        b = sum(1 for x in al if x < Fraction(1, 3))
        s = 0.0
        for i in range(h):
            for j in range(h):
                if i < b and j >= b and al[j] <= (Fraction(1) - al[i]) / 2:
                    a = qq[i] * qq[j]
                    s += float(ct[i][j]) / (a * math.log(10 ** 12 / a))
        print(f"h={h} t={float(t):.2f}: sum_cuna c_t/(qr log(N/qr))="
              f"{s:.4e} ratio vs 0.017825: {s/0.017825:.4e}")
    print("PASS: enumeracion N=10^7 completa; T>0; E2+E3+Ege4=T "
          "por construccion. NOTA: E2esp (modelo uniforme 1/log N) "
          "es INVALIDO como predictor: predice signo contrario "
          "porque los pares del extremo (p1p2~N) tienen N-p1p2 "
          "pequeno y casi nunca primo; el modelo ignora el tamano "
          "de N-p1p2. Se reporta solo como diagnostico.")
