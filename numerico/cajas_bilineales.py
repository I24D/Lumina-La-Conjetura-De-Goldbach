#!/usr/bin/env python3
"""
cajas_bilineales.py — diagnóstico [N] para matematicas/06_descomposicion_bilineal.md

Procedencia: Claude Code (claude-opus-5-5), 2026-11-23. Implementación INDEPENDIENTE
de numerico/cola_mobius_scan.py (OpenClaw): itera directamente sobre (c, d), con
vectorización en d, sin la reindexación (5.1). Sirve como comprobación cruzada.

Calcula, para D = floor(N^(1/3)) (como en 04/05) o D dado:
  T_N^prim(D) = sum_{c, (c,N)=1} sum_{D<d<=(N-1)/c} mu(d) log d Lambda(N-cd)
y su reparto en cajas diádicas c in [2^i, 2^(i+1)), d in [(D+1)2^j, (D+1)2^(j+1)).
En cada caja separa:
  P = parte con d primo      (mu(d) log d = -log d)
  Q = parte con d compuesto  (squarefree con omega(d) >= 2)
Por tanto B = P + Q. En la fila c = 1, P = -sum_{p} log p Lambda(N-p) es un bloque
diádico de la propia suma de Goldbach R(N).

ESTADO EPISTEMOLÓGICO: evidencia numérica [N]. No demuestra ningún o(N).
"""
from __future__ import annotations

import argparse
import datetime as dt
import hashlib
import math
import sys
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parent))
from goldbach import criba_simple, serie_singular, guardar  # noqa: E402


def mobius_y_lambda(M: int):
    isp = criba_simple(M)
    primos = np.nonzero(isp)[0]
    mu = np.ones(M + 1, dtype=np.int8)
    mu[0] = 0
    lam = np.zeros(M + 1)
    for p in primos.tolist():
        mu[p::p] *= -1
        pp = p * p
        if pp <= M:
            mu[pp::pp] = 0
        lp = math.log(p)
        q = p
        while q <= M:
            lam[q] = lp
            q *= p
    return mu, lam, isp


def analizar(N: int, D: int | None):
    if D is None:
        D = int(round(N ** (1 / 3)))
        while D ** 3 > N:
            D -= 1
        while (D + 1) ** 3 <= N:
            D += 1
    mu, lam, isp = mobius_y_lambda(N)
    logs = np.zeros(N + 1)
    logs[1:] = np.log(np.arange(1, N + 1))
    beta = mu.astype(np.float64) * logs  # mu(d) log d
    cajas = {}  # (i, j) -> [B, P, Q, n_pares]
    T = 0.0
    cmax = (N - 1) // (D + 1)
    for c in range(1, cmax + 1):
        if math.gcd(c, N) != 1:
            continue
        dmax = (N - 1) // c
        d = np.arange(D + 1, dmax + 1, dtype=np.int64)
        v = beta[d] * lam[N - c * d]
        if not v.any():
            continue
        es_p = isp[d]
        i = c.bit_length() - 1
        j = np.floor(np.log2(d / (D + 1))).astype(np.int64)
        for jj in np.unique(j).tolist():
            m = j == jj
            vb = v[m]
            Bv = float(vb.sum())
            Pv = float(vb[es_p[m]].sum())
            acc = cajas.setdefault((i, jj), [0.0, 0.0, 0.0, 0])
            acc[0] += Bv
            acc[1] += Pv
            acc[2] += Bv - Pv
            acc[3] += int(m.sum())
        T += float(v.sum())
    S = serie_singular(N)
    filas = []
    for (i, j), (Bv, Pv, Qv, npares) in sorted(cajas.items()):
        filas.append({"i": i, "j": j, "c_rango": [2 ** i, 2 ** (i + 1) - 1],
                      "d_rango": [(D + 1) * 2 ** j, (D + 1) * 2 ** (j + 1) - 1],
                      "B_sobre_N": Bv / N, "P_sobre_N": Pv / N, "Q_sobre_N": Qv / N,
                      "pares": npares})
    return D, S, T, filas


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--N", type=float, nargs="+", required=True)
    ap.add_argument("--D", type=int, default=None)
    args = ap.parse_args()
    script = Path(__file__)
    salida = {
        "script": script.name,
        "script_sha256": hashlib.sha256(script.read_bytes()).hexdigest(),
        "fecha_utc": dt.datetime.now(dt.timezone.utc).isoformat(timespec="seconds"),
        "estado_epistemologico": "EVIDENCIA NUMÉRICA — no es demostración",
        "convencion_cajas": "c en [2^i, 2^(i+1)), d en [(D+1)2^j, (D+1)2^(j+1))",
        "resultados": [],
    }
    for Nf in args.N:
        N = int(Nf)
        if N % 2:
            N += 1
        D, S, T, filas = analizar(N, args.D)
        maxB = max(filas, key=lambda f: abs(f["B_sobre_N"]))
        sumP = sum(f["P_sobre_N"] for f in filas)
        sumQ = sum(f["Q_sobre_N"] for f in filas)
        fila_c1 = [f for f in filas if f["i"] == 0]
        print(f"\nN = {N:,}  D = {D}  S(N) = {S:.6f}")
        print(f"  T^prim/N = {T / N:.7f}   T/(N S) = {T / (N * S):.7f}")
        print(f"  suma P/N (d primo) = {sumP:.5f}   suma Q/N (d compuesto) = {sumQ:.5f}")
        print(f"  caja con max |B|/N: i={maxB['i']} j={maxB['j']}  B/N={maxB['B_sobre_N']:.5f} "
              f"(P/N={maxB['P_sobre_N']:.5f}, Q/N={maxB['Q_sobre_N']:.5f})")
        print("  fila c = 1 (P = -bloque diádico de Goldbach):")
        for f in fila_c1:
            print(f"    d en {f['d_rango']}:  B/N={f['B_sobre_N']:+.5f}  "
                  f"P/N={f['P_sobre_N']:+.5f}  Q/N={f['Q_sobre_N']:+.5f}")
        salida["resultados"].append({"N": N, "D": D, "S": S, "T_sobre_N": T / N,
                                     "T_sobre_NS": T / (N * S), "suma_P_sobre_N": sumP,
                                     "suma_Q_sobre_N": sumQ, "cajas": filas})
    ruta = guardar("cajas_bilineales", salida)
    print("->", ruta.name)
    return 0


if __name__ == "__main__":
    sys.exit(main())
