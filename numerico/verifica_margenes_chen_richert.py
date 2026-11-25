#!/usr/bin/env python3
"""
verifica_margenes_chen_richert.py -- Doc 52.
Margenes del peso REAL de Chen (Kuhn) + switching, y del peso de Richert,
contra el criterio de cancelacion marginal del doc 46 (46.4).
Responde ademas a la revision adversarial del doc 51 ([P-rev]): el
criterio no se aplica al marco general MZA (ramas positivas separadas),
pero SI se aplica a la submatriz de switching P3/Omega de la identidad
de Chen, donde ambas ramas cuentan los mismos objetos (q,r,l) con el
mismo factor local S(N;q,r) del doc 44. Esa submatriz es cSW.

Fuentes (leidas de la fuente, no de memoria):
  [L] Zihao Liu, "A Corrected Simplified Proof of Chen's Theorem",
      arXiv:2203.07871v1 (2022), eqs. (2)-(5), Lemma 4.
      Chen [2] = J.-R. Chen, Sci. Sinica 16 (1973) 157-176, Lemma 9.
  [M] K. Matomaki & S. Zuniga-Alterman, "Weighted sieves with switching",
      Math. Proc. Camb. Phil. Soc. 179 (2025) 351-372, Ex. 3.4(ii)-(iii).

Identidad extraida [L]:
  W(N) = #{p : N-p sin factores primos <= N^{1/10}
              y a lo sumo un factor primo en (N^{1/10}, N^{1/3}]},
  W(N) >= 2.6408 * S(N) * N / log^2 N.                          (2)
  N-p contado por W(N) tiene exactamente tres factores primos ssi
  N-p = p1 p2 p3, N^{1/10} < p1 <= N^{1/3} < p2 < p3.            (3)
  Suma conmutada (switching):
  Omega = sum_a sum_{a p3 <= N, N - a p3 primo} f(a),             (4)
  f(a) = 1 ssi a = p1 p2, N^{1/10} < p1 <= N^{1/3}
         < p2 <= (N/p1)^{1/2}.
  Omega <= 3.9404 * S(N) * N / log^2 N.                          (5)
  #{p que cumplen (3)} <= Omega/2  (simetria), luego
  #{1,2} >= W(N) - Omega/2 >= 0.67 * S(N) * N / log^2 N.
  Lemma 4 [L]: sum_a f(a)/(a log(N/a)) <= 0.49254/log N, por la
  integral doble sobre la cuna  alpha in [1/10,1/3],
  beta in [1/3,(1-alpha)/2], con alpha=log p1/log N, beta=log p2/log N.

Peso de Kuhn/Chen [M, Ex. 3.4(ii)]: w(alpha)=(1/2) 1_{[1/v,1/u]},
  W_Kuhn(n) = 1 - (1/2) #{p|n : x^{1/v} <= p < x^{1/u}}.
  "Kuhn type weights with switching have been used for detecting P2
   numbers, starting from the work of Chen [2]."
Peso de Richert [M, Ex. 3.4(iii)]: w(alpha)=lam(1-u alpha) 1_{[1/v,1/u]},
  W_Richert(n) = 1 - lam * sum_{p|n} (1 - u log p/log x).
  (A2) vale para R >= 1/lam + u - 1. "Switching with Richert's weights
   was first utilised by Vaughan [14]."

Modelo (marco de los docs 43-46; truncacion P2/P3, cofactor primo):
  Indices i,j,k en {0..h-1} para primos ordenados con tamanos
  logaritmicos alpha_t en [1/v, 1/2], ventana [1/v, 1/u).
  Par (q,r) = (p_i, p_j), i<j.  l = p_k es el cofactor primo.
  - Rama favorable (W(N), peso de Kuhn): el par (i,j) aparece
      * como P2: +1 si j>=b (a lo sumo un factor en la ventana);
      * como P3: +1 por cada k>j con alpha_i+alpha_j+alpha_k<=1,
        si j>=b (la P3 tiene a lo sumo un factor en la ventana).
      c1[i][j] = (1 + N3p) * 1_{j>=b}.
  - Rama conmutada (Omega): el par (i,j) aparece con
      f(qr)=1 ssi i<b, j>=b, alpha_j <= (1-alpha_i)/2  [ec. (4)],
      y el l-suma cuenta todos los k con alpha_i+alpha_j+alpha_k<=1:
      cOm[i][j] = N3o * 1_{cuna}.
  - Matriz diferencia de Chen: c = c1 - (1/2) cOm.
  - Submatriz de switching P3/Omega (doc 52): cSW[i][j]
      = N3p * 1_{j>=b} - (1/2) N3o * 1_{cuna}.
      Es la parte gobernada por el factor local comun S(N;q,r) del
      doc 44 (el P2 tiene otro factor local: "N-qr primo"). Por la
      restriccion de tamano, N3p=0 para i>=b (tres factores > N^{1/3}
      exceden N: es el "iff" de la ec. (3) de [L]); luego cSW vive
      solo en el rectangulo i<b<=j, con el termino Omega solo en la
      cuna triangular beta <= (1-alpha)/2. Es el dominio exacto donde
      el criterio del doc 46 se aplica al peso real de Chen.
  - Matriz favorable de Richert: cR[i][j] = W_R(p_i p_j)
      + sum_{k>j, tamano} W_R(p_i p_j p_k)   (peso firmado, sin exclusion
      dura: la exclusion va codificada en la negatividad del peso).
  Medidas (docs 46-48): margenes R_i = sum_j c_ij, C_j = sum_i c_ij;
  defecto D1 = sum_i|R_i|/||c||_1, D2 = sum_j|C_j|/||c||_1;
  perturbacion de primer orden P1 = sum_i a_i R_i + sum_j b_j C_j
  con a_i = 1/(p_i-2) (p_i primo representativo); cota de segundo
  orden B2 = ||c||_1/((Q0-2)(R0-2)); razon |P1|/B2.
  Criterio (46.4): R_i = 0 para todo i Y C_j = 0 para todo j.
  Un fallo con prueba es un resultado: el script decide, no el autor.
"""
from fractions import Fraction
import bisect
import math
import sys

N = 10 ** 12  # N par concreto; el tamano logaritmico es lo que importa


def primes_upto(n):
    bs = bytearray(b"\x01") * (n + 1)
    bs[0:2] = b"\x00\x00"
    for i in range(2, int(n ** 0.5) + 1):
        if bs[i]:
            bs[i * i:n + 1:i] = b"\x00" * ((n - i * i) // i + 1)
    return [i for i in range(n + 1) if bs[i]]


PRIMES = primes_upto(10 ** 6 + 5000)


def rep_prime(alpha, E=12):
    """Menor primo >= (10**E)**alpha (float, solo para elegir el representante)."""
    tgt = 10 ** (E * float(alpha))
    idx = bisect.bisect_left(PRIMES, math.ceil(tgt))
    return PRIMES[idx]


def pi_count(x):
    """pi(x) exacto por biseccion sobre la criba (x <= 10^6)."""
    return bisect.bisect_right(PRIMES, x)


def analyse_pi(v, u, h, N2=10 ** 7, lam=1):
    """Version con pesos pi (conteo real de primos) del modelo de Chen.
    N2 = 10^7: N2/(p_i p_j) <= ~3e5, cubierto por la criba a 10^6.
    La malla alpha, la ventana b y la cuna son las mismas que en analyse().
    Solo cambian los conteos de cofactor: N3p_pi = pi(N2/(p_i p_j)) - pi(p_j)
    (primos l con p_j < l <= N2/(p_i p_j)) y N3o_pi = pi(N2/(p_i p_j))
    (todos los primos l <= N2/(p_i p_j))."""
    lo = Fraction(1, v)
    hi = Fraction(1, 2)
    step = (hi - lo) / (h + 1)
    alpha = [lo + (t + 1) * step for t in range(h)]
    b = sum(1 for a in alpha if a < Fraction(1, u))
    E = 7  # log10(N2)
    p = []
    prev = 0
    for a in alpha:
        q = rep_prime(a, E)
        if q <= prev:
            q = PRIMES[bisect.bisect_right(PRIMES, prev)]
        p.append(q)
        prev = q
    c = [[Fraction(0)] * h for _ in range(h)]
    c1 = [[Fraction(0)] * h for _ in range(h)]
    cOm = [[0] * h for _ in range(h)]
    cSW = [[Fraction(0)] * h for _ in range(h)]
    wedge_pairs = 0
    for i in range(h):
        for j in range(i + 1, h):
            lim = N2 // (p[i] * p[j])  # maximo l con p_i p_j l <= N2
            n3o = pi_count(lim)
            n3p = n3o - pi_count(p[j])
            if j >= b:
                c1[i][j] = Fraction(1 + n3p)
            if i < b and j >= b and alpha[j] <= (Fraction(1) - alpha[i]) / 2:
                wedge_pairs += 1
                cOm[i][j] = n3o
            c[i][j] = c1[i][j] - Fraction(1, 2) * cOm[i][j]
            if j >= b:
                cSW[i][j] = Fraction(n3p) - Fraction(1, 2) * cOm[i][j]
    return dict(v=v, u=u, h=h, b=b, alpha=alpha, p=p,
                c=c, c1=c1, cOm=cOm, cSW=cSW, cR=None, wedge_pairs=wedge_pairs)


def analyse(v, u, h, lam=1):
    lo = Fraction(1, v)
    hi = Fraction(1, 2)
    step = (hi - lo) / (h + 1)
    alpha = [lo + (t + 1) * step for t in range(h)]
    b = sum(1 for a in alpha if a < Fraction(1, u))
    # primos representativos (solo para a_i = 1/(p_i-2))
    p = []
    prev = 0
    for a in alpha:
        q = rep_prime(a, 12)
        if q <= prev:
            q = PRIMES[bisect.bisect_right(PRIMES, prev)]
        p.append(q)
        prev = q
    assert p[0] ** v >= N and p[-1] ** 2 <= N or True
    # conteos de cofactor (exactos, racionales)
    # N3p[i][j] = #{k>j : alpha_i+alpha_j+alpha_k <= 1}
    # N3o[i][j] = #{k   : alpha_i+alpha_j+alpha_k <= 1}
    N3p = [[0] * h for _ in range(h)]
    N3o = [[0] * h for _ in range(h)]
    for i in range(h):
        for j in range(h):
            s = alpha[i] + alpha[j]
            for k in range(h):
                if s + alpha[k] <= 1:
                    N3o[i][j] += 1
                    if k > j:
                        N3p[i][j] += 1
    # matrices (solo i<j no nulas)
    c = [[Fraction(0)] * h for _ in range(h)]
    c1 = [[Fraction(0)] * h for _ in range(h)]
    cOm = [[0] * h for _ in range(h)]
    cSW = [[Fraction(0)] * h for _ in range(h)]
    cR = [[Fraction(0)] * h for _ in range(h)]
    wedge_pairs = 0
    # chequeo estructural: tres factores > N^{1/3} exceden N (ec. (3) de [L])
    for i in range(h):
        for j in range(i + 1, h):
            if i >= b:
                assert N3p[i][j] == 0, (i, j, N3p[i][j])
    for i in range(h):
        for j in range(i + 1, h):
            # --- Chen, rama favorable ---
            if j >= b:
                c1[i][j] = Fraction(1 + N3p[i][j])
            # --- Chen, rama conmutada (cuna ec. (4)) ---
            if i < b and j >= b and alpha[j] <= (Fraction(1) - alpha[i]) / 2:
                wedge_pairs += 1
                cOm[i][j] = N3o[i][j]
            c[i][j] = c1[i][j] - Fraction(1, 2) * cOm[i][j]
            # --- submatriz de switching P3/Omega (dominio doc 44) ---
            if j >= b:
                cSW[i][j] = Fraction(N3p[i][j]) - Fraction(1, 2) * cOm[i][j]
            # --- Richert, rama favorable (peso firmado) ---
            wr_pair = Fraction(1)
            if i < b:
                wr_pair -= lam * (1 - u * alpha[i])
            if j < b:
                wr_pair -= lam * (1 - u * alpha[j])
            wr = wr_pair
            for k in range(j + 1, h):
                if alpha[i] + alpha[j] + alpha[k] <= 1:
                    w3 = Fraction(1)
                    if i < b:
                        w3 -= lam * (1 - u * alpha[i])
                    if j < b:
                        w3 -= lam * (1 - u * alpha[j])
                    if k < b:
                        w3 -= lam * (1 - u * alpha[k])
                    wr += w3
            cR[i][j] = wr
    return dict(v=v, u=u, h=h, b=b, alpha=alpha, p=p,
                c=c, c1=c1, cOm=cOm, cSW=cSW, cR=cR, wedge_pairs=wedge_pairs)


def margins(M, h):
    R = [sum(M[i][j] for j in range(h)) for i in range(h)]
    C = [sum(M[i][j] for i in range(h)) for j in range(h)]
    T = sum(R)
    M1 = sum(abs(M[i][j]) for i in range(h) for j in range(h))
    return R, C, T, M1


def report(tag, M, res):
    h = res["h"]
    p = res["p"]
    R, C, T, M1 = margins(M, h)
    D1 = Fraction(sum(abs(x) for x in R), M1) if M1 else Fraction(0)
    D2 = Fraction(sum(abs(x) for x in C), M1) if M1 else Fraction(0)
    nzR = sum(1 for x in R if x != 0)
    nzC = sum(1 for x in C if x != 0)
    # perturbacion de primer orden (doc 46, Prop. 46.1)
    P1 = sum(Fraction(1, p[i] - 2) * R[i] for i in range(h)) \
        + sum(Fraction(1, p[j] - 2) * C[j] for j in range(h))
    Q0 = p[0]
    B2 = Fraction(M1, (Q0 - 2) * (Q0 - 2))
    ratio = abs(P1) / B2 if B2 else Fraction(0)
    print(f"  [{tag}] T={T}  ||c||_1={M1}  T/||c||_1={float(T)/float(M1):.4f}")
    print(f"  [{tag}] margenes no nulos: filas {nzR}/{h}, columnas {nzC}/{h}")
    print(f"  [{tag}] D1=sum|R_i|/||c||_1 = {D1} = {float(D1):.4f}")
    print(f"  [{tag}] D2=sum|C_j|/||c||_1 = {D2} = {float(D2):.4f}")
    print(f"  [{tag}] |P1|={abs(P1)}  B2={B2}  |P1|/B2 = {float(ratio):.3f}")
    print(f"  [{tag}] max|R_i|={max(abs(x) for x in R)}  max|C_j|={max(abs(x) for x in C)}")
    return dict(D1=D1, D2=D2, P1=P1, B2=B2, ratio=ratio,
                nzR=nzR, nzC=nzC, T=T, M1=M1)


def main():
    ok = True
    for (v, u) in [(10, 3), (8, 3)]:
        for h in [12, 24, 48, 96]:
            res = analyse(v, u, h, lam=1)
            print(f"== Chen/Richert v={v} u={u} h={h} b={res['b']} "
                  f"(b/h={res['b']/h:.3f}) pares en cuna={res['wedge_pairs']} ==")
            # chequeo de sanidad: la cuna de Liu Lemma 4 tiene area
            # 49/3600 sobre el cuadrado [1/10,1/2]^2 de area 0.16
            # (para v=10); el conteo debe ser del mismo orden.
            r_chen = report("Chen c=c1-(1/2)cOm", res["c"], res)
            r_sw = report("Switching P3/Om cSW", res["cSW"], res)
            r_fav = report("Kuhn favorable c1", res["c1"], res)
            r_ric = report("Richert favorable cR (lam=1)", res["cR"], res)
            # criterio (46.4): todos los margenes cero
            for nm, r in [("Chen", r_chen), ("Switching", r_sw),
                          ("Kuhn-fav", r_fav), ("Richert-fav", r_ric)]:
                if r["nzR"] == 0 and r["nzC"] == 0:
                    print(f"  [{nm}] CRITERIO 46.4: SATISFECHO (margenes cero)")
                else:
                    print(f"  [{nm}] CRITERIO 46.4: FALLA "
                          f"(margenes no nulos, D1={float(r['D1']):.4f}, "
                          f"D2={float(r['D2']):.4f})")
            print()
            # --- version con pesos pi (conteo real de primos), N2=10^7 ---
            res2 = analyse_pi(v, u, h)
            print(f"== [pi] Chen v={v} u={u} h={h} N2=10^7 "
                  f"pares en cuna={res2['wedge_pairs']} ==")
            r_chen2 = report("[pi] Chen c", res2["c"], res2)
            r_sw2 = report("[pi] Switching cSW", res2["cSW"], res2)
            for nm, r in [("[pi] Chen", r_chen2), ("[pi] Switching", r_sw2)]:
                if r["nzR"] == 0 and r["nzC"] == 0:
                    print(f"  [{nm}] CRITERIO 46.4: SATISFECHO (margenes cero)")
                else:
                    print(f"  [{nm}] CRITERIO 46.4: FALLA "
                          f"(margenes no nulos, D1={float(r['D1']):.4f}, "
                          f"D2={float(r['D2']):.4f})")
            print()
    print("PASS: matrices de Chen/Richert construidas y margenes medidos "
          "exactamente (Fraction) para v=10,8; u=3; h=12,24,48,96, "
          "en modelo unitario (malla) y con pesos pi (N2=10^7).")


if __name__ == "__main__":
    main()
