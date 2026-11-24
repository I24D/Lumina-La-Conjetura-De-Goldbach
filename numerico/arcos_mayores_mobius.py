#!/usr/bin/env python3
"""
arcos_mayores_mobius.py — prueba [N] de una fórmula explícita para el sesgo de la
correlación Möbius–Goldbach M0(N) = sum_{n<N} mu(n) Lambda(N-n).

Procedencia: Claude Code (claude-opus-5-5), 2026-11-23.

Hipótesis (heurística de arcos mayores / desarrollo de Ramanujan de Lambda):
    Lambda(m) ≈ sum_{q<=Q} mu(q) c_q(m) / phi(q)      (c_q = suma de Ramanujan)
    =>  M0(N) ≈ P_Q(N) := sum_{q<=Q} mu(q)/phi(q) * sum_{n<N} mu(n) c_q(N-n).
El término q = 1 es la función de Mertens M(N-1).

Se mide, para N en (X/2, X]:
  - corr(M0, P_Q) y la fracción de varianza explicada R^2;
  - si P_Q explica el sesgo por clases N mod 3, 5, 7 y la correlación con S(N);
  - la varianza del residuo M0 - P_Q frente al modelo aleatorio V0.

ESTADO EPISTEMOLÓGICO: evidencia numérica [N]. No es demostración.
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
from goldbach import guardar, serie_singular  # noqa: E402
from cajas_bilineales import mobius_y_lambda  # noqa: E402
from mobius_goldbach_aleatoriedad import conv  # noqa: E402


def ramanujan(q: int, m: np.ndarray) -> np.ndarray:
    """c_q(m) = mu(q/g) phi(q)/phi(q/g), con g = gcd(q, m)."""
    g = np.gcd(q, m)
    out = np.empty(m.shape, dtype=np.float64)
    for gv in np.unique(g).tolist():
        t = q // gv
        out[g == gv] = mu_int(t) * phi_int(q) / phi_int(t)
    return out


def mu_int(n: int) -> int:
    r, p, m = 1, 2, n
    while p * p <= m:
        if m % p == 0:
            m //= p
            if m % p == 0:
                return 0
            r = -r
        p += 1
    return -r if m > 1 else r


def phi_int(n: int) -> int:
    r, p, m = n, 2, n
    while p * p <= m:
        if m % p == 0:
            while m % p == 0:
                m //= p
            r -= r // p
        p += 1
    if m > 1:
        r -= r // m
    return r


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--hasta", type=float, default=1e6)
    ap.add_argument("--Q", type=int, nargs="+", default=[1, 2, 6, 30, 60])
    args = ap.parse_args()
    X = int(args.hasta)
    mu, lam, _ = mobius_y_lambda(X)
    muf = mu.astype(np.float64)
    muf[0] = 0.0
    M0 = conv(muf, lam, X)
    V0 = conv(muf ** 2, lam ** 2, X)
    pares = np.arange(X // 2 + (X // 2) % 2, X + 1, 2)
    m_idx = np.arange(X + 1)
    Qmax = max(args.Q)
    acumulado = np.zeros(X + 1)
    filas = []
    muestra = pares[:: max(1, pares.size // 5000)]
    Smuestra = np.array([serie_singular(int(N)) for N in muestra])
    y = M0[pares]
    for q in range(1, Qmax + 1):
        mq = mu_int(q)
        if mq != 0:
            cq = ramanujan(q, m_idx)
            cq[0] = 0.0
            acumulado += (mq / phi_int(q)) * conv(muf, cq, X)
        if q in args.Q:
            p = acumulado[pares]
            r = y - p
            c = float(np.corrcoef(y, p)[0, 1])
            fila = {"Q": q, "corr(M0,P_Q)": c, "R2": c * c,
                    "var(resid)/media(V0)": float(r.var() / V0[pares].mean()),
                    "var(M0)/media(V0)": float(y.var() / V0[pares].mean())}
            z_res = r / np.sqrt(V0[pares])
            for mod in (3, 5, 7):
                fila[f"media z_resid | {mod}|N"] = float(z_res[pares % mod == 0].mean())
                fila[f"media z_resid | {mod}∤N"] = float(z_res[pares % mod != 0].mean())
            zr_m = (M0[muestra] - acumulado[muestra]) / np.sqrt(V0[muestra])
            fila["corr(z_resid, S(N))"] = float(np.corrcoef(zr_m, Smuestra)[0, 1])
            filas.append(fila)
            print(fila)
    z = y / np.sqrt(V0[pares])
    base = {f"media z | {mod}|N": float(z[pares % mod == 0].mean()) for mod in (3, 5, 7)}
    zm = M0[muestra] / np.sqrt(V0[muestra])
    base["corr(z, S(N))"] = float(np.corrcoef(zm, Smuestra)[0, 1])
    print("sin corregir:", base)
    script = Path(__file__)
    datos = {"script": script.name,
             "script_sha256": hashlib.sha256(script.read_bytes()).hexdigest(),
             "fecha_utc": dt.datetime.now(dt.timezone.utc).isoformat(timespec="seconds"),
             "estado_epistemologico": "EVIDENCIA NUMÉRICA — no es demostración",
             "X": X, "ventana": [X // 2, X], "sin_corregir": base, "por_Q": filas}
    print("->", guardar(f"arcos_mayores_mobius_{X}", datos).name)
    return 0


if __name__ == "__main__":
    sys.exit(main())
