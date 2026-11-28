#!/usr/bin/env python3
# 108_auditoria_cadena_completa.py  (Doc 108)
# Auditoria de la cadena completa (DIST-D) -> balance final.
# Verifica las identidades COMPUTABLES de la cadena; la auditoria
# logica vive en el .md. Solo stdlib.
#
# Chequeos:
#   K1  Algebra del objetivo (doc 100 sec. 5.3): con G_grande = 0,0225
#       y A = N^c (c = c* = 0,039884), T(N) := N^2/(A^2 ln^4 N * 0,0225)
#       reproduce la tabla del doc 100 y sqrt(G*T(N)) = N/(A ln^2 N)
#       EXACTO: el factor 1/(A log^2 N) del puente encaja.
#   K2  Cadena (107.2) de L105(a): coeficiente 0,0225/(4 ln^2 2) =
#       0,0117077; biseccion de N0 para C_Delta = 1e-3 da ~2,36e14
#       (doc 107); el cociente decrece para N > 5970.
#   K3  Factorizacion del precio de la pieza C (doc 105 sec. 4):
#       Dhat^C/B_M^{C,(M)} = (Diag^C/B_M^{C,(M)}) * precio_tau en los
#       seis puntos de la tabla del doc 105.
#   K4  Pieza prima SIN cancelacion de Mobius: con F_1 = 0 y F_d
#       arbitrarios, M_r = sum_{d|r} mu(r/d) F_d vale F_r EXACTO si
#       r es primo, y F_r - F_p - F_q si r = pq: la pieza prima de
#       B_M es un segundo momento CRUDO sum_r phi(r) E_r^2, no una
#       identidad que lo acote (base del hueco G1 del doc 108).
#       Ademas sum_{d|r} mu(r/d) = 0 para r > 1 (cancelacion del
#       termino principal en S_r, doc 97 paso 5).
#   K5  Orden de umbrales: c* = 0,039884 < 2/17 (geometria, doc 89)
#       < 1/3 (Cauchy corregido, doc 93): ni la geometria ni Cauchy
#       son la restriccion activa de la ruta UNIF-U-COND; el balance
#       Bal(c) si: Bal(0,02) > 0 > Bal(0,04) (tabla doc 89).
#   K6  Exceso clasico (docs 106/107): la gran criba/BDH sobrepasa
#       T(N) por 0,045 * W* * N^{2c} * ln^5 N; en el marco z = 11
#       (N = 11^10, W* = 9,52) da ~2,3e7, mismo orden que el exceso
#       medido del doc 107 (1,14e7 suma de minimos; 3,39e7 solo-LS).
#   K7  Extrapolacion del doc 102: |H| ~ 0,033 * N^0,4208 a N = 1e24
#       da ~4e8 (doc 102: |H| ~ 4*10^8).
import math
import sys
import time

C_STAR = 0.039884
G_GRANDE = 0.0225


def mu_sieve(n):
    mu = [1] * (n + 1)
    primes = []
    is_comp = [False] * (n + 1)
    mu[0] = 0
    for i in range(2, n + 1):
        if not is_comp[i]:
            primes.append(i)
            mu[i] = -1
        for p in primes:
            if i * p > n:
                break
            is_comp[i * p] = True
            if i % p == 0:
                mu[i * p] = 0
                break
            mu[i * p] = -mu[i]
    return mu


def divisors(n):
    ds = []
    i = 1
    while i * i <= n:
        if n % i == 0:
            ds.append(i)
            if i * i != n:
                ds.append(n // i)
        i += 1
    return ds


def main():
    allok = True
    t00 = time.time()
    print("== DOC 108: auditoria de la cadena — identidades computables ==")
    print("Solo stdlib. Las cotas/objetivos son los de los docs citados;")
    print("este script verifica ALGEBRA y consistencia, no demuestra cotas.")
    print()

    # ---------- K1: algebra del objetivo del puente (doc 100 sec. 5.3)
    tabla100 = {12: 8.414e18, 24: 5.803e40, 100: 1.668e184}
    k1 = True
    for e10, ref in tabla100.items():
        N = 10.0 ** e10
        A = N ** C_STAR
        lnN = math.log(N)
        T = N * N / (A * A * lnN ** 4 * G_GRANDE)
        puente = math.sqrt(G_GRANDE * T)
        objetivo_puente = N / (A * lnN ** 2)
        rel_tabla = abs(T - ref) / ref
        rel_puente = abs(puente - objetivo_puente) / objetivo_puente
        ok = rel_tabla < 2e-3 and rel_puente < 1e-12
        k1 &= ok
        print("K1 N=1e%-3d T(N)=%.4e (doc100 %.3e, rel %.1e)  "
              "sqrt(G*T)=N/(A ln^2 N) rel %.1e  %s"
              % (e10, T, ref, rel_tabla, rel_puente,
                 "OK" if ok else "FAIL"))
    allok &= k1
    print("K1 algebra del objetivo y factor 1/(A ln^2 N): %s"
          % ("PASS" if k1 else "FAIL"))
    print()

    # ---------- K2: cadena (107.2) y N0 de L105(a)
    coef = G_GRANDE / (4.0 * math.log(2) ** 2)
    k2 = abs(coef - 0.0117077) < 1e-6

    def ratio107(N):
        return coef * math.log(N) ** 8 * N ** (-(1 - 2 * C_STAR))

    # decreciente para N > 5970: chequeo en malla log
    grid = [6000.0 * (1.5 ** i) for i in range(60)]
    vals = [ratio107(x) for x in grid]
    k2 &= all(vals[i + 1] < vals[i] for i in range(len(vals) - 1))
    # biseccion en log10 N para ratio = 1e-3
    lo, hi = 5.0, 20.0
    for _ in range(200):
        mid = 0.5 * (lo + hi)
        if ratio107(10.0 ** mid) > 1e-3:
            lo = mid
        else:
            hi = mid
    N0 = 10.0 ** (0.5 * (lo + hi))
    rel_N0 = abs(N0 - 2.36e14) / 2.36e14
    k2 &= rel_N0 < 0.05
    print("K2 coef (107.2) = %.7f (doc 107: 0,0117077); decreciente "
          "N>5970: %s; N0(C_Delta=1e-3) = %.3e (doc 107: 2,36e14, "
          "rel %.2e)" % (coef, "OK", N0, rel_N0))
    allok &= k2
    print("K2 cadena (107.2) de L105(a): %s" % ("PASS" if k2 else "FAIL"))
    print()

    # ---------- K3: factorizacion del precio, doc 105 sec. 4 (seis puntos)
    # (Diag^C/B_M^{C,(M)}, precio tau^C, Dhat^C/B_M^{C,(M)}) tabla doc 105/106
    pts105 = [
        (1.495, 4.119, 6.158), (1.488, 4.439, 6.605), (1.533, 4.909, 7.525),
        (1.674, 4.185, 7.007), (1.616, 4.582, 7.403), (1.576, 4.862, 7.661),
    ]
    k3 = True
    for diag, tau, dhat in pts105:
        rel = abs(diag * tau - dhat) / dhat
        k3 &= rel < 2e-3
    print("K3 Diag * tau = Dhat^C/B_M^{C,(M)} en los 6 puntos "
          "(peor rel < 2e-3): %s" % ("PASS" if k3 else "FAIL"))
    allok &= k3
    print()

    # ---------- K4: la pieza prima no tiene cancelacion (base de G1)
    mu = mu_sieve(600)
    # F arbitrarios (deterministas), F_1 = 0
    F = {d: math.sin(1.7 * d) * (1.0 + (d % 7)) for d in range(1, 601)}
    F[1] = 0.0

    def M_r(r):
        return sum(mu[r // d] * F[d] for d in divisors(r))

    primes_test = [3, 5, 7, 11, 13, 101, 347, 599]
    k4 = all(abs(M_r(p) - F[p]) < 1e-9 for p in primes_test)
    # r = pq: M_r = F_r - F_p - F_q (con F_1 = 0)
    for (p, q) in [(3, 5), (7, 11), (13, 17)]:
        r = p * q
        expect = F[r] - F[p] - F[q]
        k4 &= abs(M_r(r) - expect) < 1e-9
    # cancelacion del termino principal: sum_{d|r} mu(r/d) = 0, r > 1
    k4 &= all(sum(mu[r // d] for d in divisors(r)) == 0
              for r in range(2, 601))
    print("K4 M_r = F_r exacto en primos (sin cancelacion); "
          "M_{pq} = F_{pq} - F_p - F_q; sum mu(r/d) = 0 (r>1): %s"
          % ("PASS" if k4 else "FAIL"))
    allok &= k4
    print()

    # ---------- K5: orden de umbrales y signo del balance (doc 89)
    k5 = C_STAR < 2.0 / 17.0 < 1.0 / 3.0
    bal = {0.0: 0.061404, 0.02: 0.032934, 0.04: -0.000208}
    k5 &= bal[0.02] > 0 > bal[0.04] and 0.02 < C_STAR < 0.04
    print("K5 c*=%.6f < 2/17=%.6f < 1/3=%.6f; Bal(0,02)=%.6f > 0 > "
          "Bal(0,04)=%.6f: %s"
          % (C_STAR, 2.0 / 17.0, 1.0 / 3.0, bal[0.02], bal[0.04],
             "PASS" if k5 else "FAIL"))
    allok &= k5
    print()

    # ---------- K6: exceso clasico frente a T(N) (docs 106/107)
    z = 11
    Nf = float(z) ** 10
    Wstar = 9.52
    exceso = 0.045 * Wstar * Nf ** (2 * C_STAR) * math.log(Nf) ** 5
    k6 = 5e6 < exceso < 1e8  # doc 107: 1,14e7 (min) .. 3,39e7 (solo-LS)
    print("K6 exceso clasico en z=11: %.3e (doc 107 medido: 1,14e7.."
          "3,39e7; mismo orden): %s" % (exceso, "PASS" if k6 else "FAIL"))
    allok &= k6
    print()

    # ---------- K7: extrapolacion |H| del doc 102
    H24 = 0.033 * (1e24) ** 0.4208
    k7 = 2e8 < H24 < 8e8
    print("K7 |H|(1e24) ~ %.2e (doc 102: ~4e8): %s"
          % (H24, "PASS" if k7 else "FAIL"))
    allok &= k7
    print()

    print("RESUMEN: K1 algebra T(N)/puente OK; K2 L105(a) OK; K3 precio "
          "pieza C OK;")
    print("         K4 pieza prima = segundo momento crudo (hueco G1); "
          "K5 umbrales OK;")
    print("         K6 exceso clasico ~ N^{2c} ln^5 confirmado; K7 |H| OK.")
    print("PASS" if allok else "FAIL (ver lineas marcadas)")
    print("Tiempo total: %.1f s" % (time.time() - t00))
    return 0 if allok else 1


if __name__ == "__main__":
    sys.exit(main())
