#!/usr/bin/env python3
"""
li_liu_check.py — tarea T-LL2 (ver matematicas/02_analisis_Li_Liu_2026.md).

Calcula la correlación de Goldbach con pesos de Li & Liu (2026),
  G(N) = (1/N) * sum_{n+m=N, (nm,N)=1} Λ(n)Λ(m) φ(n/N) φ(m/N),
y la compara con el término principal c_φ·S(N), con c_φ = ∫φ².
Mide error/δ_N, con δ_N = sqrt(S(N)/N)·(log N)^{3/2}. Si la Hipótesis R
fuera cierta (Teorema 1.2), el cociente estaría acotado por una constante.
Una tendencia creciente con N sería evidencia numérica contra R.

ESTADO EPISTEMOLÓGICO: evidencia numérica [N]; no demuestra ni refuta nada.
Peso usado: φ(u) = exp(1 - 1/(1 - t²)) con t = (u - 1/2)/(1/6), es decir,
un bump C^∞ con soporte en [1/3, 2/3], simétrico, φ(1/2) = 1.
"""
from __future__ import annotations

import datetime as dt
import hashlib
import json
import math
import sys
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parent))
from goldbach import criba_simple, serie_singular, guardar  # noqa: E402


def von_mangoldt(M: int) -> np.ndarray:
    lam = np.zeros(M + 1)
    isp = criba_simple(M)
    for p in np.nonzero(isp)[0]:
        p = int(p)
        lp = math.log(p)
        q = p
        while q <= M:
            lam[q] = lp
            q *= p
    return lam


def phi(u: np.ndarray) -> np.ndarray:
    t = (u - 0.5) / (1.0 / 6.0)
    out = np.zeros_like(u)
    m = np.abs(t) < 1
    out[m] = np.exp(1.0 - 1.0 / (1.0 - t[m] ** 2))
    return out


def c_phi() -> float:
    u = np.linspace(1 / 3, 2 / 3, 2_000_001)
    return float(np.trapezoid(phi(u) ** 2, u))


def G(N: int, lam: np.ndarray) -> float:
    n = np.arange(N // 3, 2 * N // 3 + 2, dtype=np.int64)
    n = n[(n >= 1) & (n < N)]
    w = lam[n] * lam[N - n]
    nz = w != 0
    n, w = n[nz], w[nz]
    cop = np.gcd(n, N) == 1
    n, w = n[cop], w[cop]
    return float(np.sum(w * phi(n / N) * phi((N - n) / N)) / N)


def main() -> int:
    M = int(float(sys.argv[1])) if len(sys.argv) > 1 else 10_000_000
    lam = von_mangoldt(M)
    c = c_phi()
    rng = np.random.default_rng(2026)
    filas = []
    for k in range(10, int(math.log2(M)) + 1):
        lo, hi = 1 << k, min(1 << (k + 1), M)
        muestra = sorted(set((rng.integers(lo // 2, hi // 2, 40) * 2).tolist()))
        for N in muestra:
            S = serie_singular(N)
            g = G(N, lam)
            err = g - c * S
            delta = math.sqrt(S / N) * math.log(N) ** 1.5
            filas.append({"N": N, "S": S, "G": g, "cS": c * S, "error": err,
                          "error_rel": err / (c * S), "error_sobre_delta": err / delta})
    # resumen por tramos diádicos
    resumen = []
    for k in range(10, int(math.log2(M)) + 1):
        f = [r for r in filas if (1 << k) <= r["N"] < (1 << (k + 1))]
        if not f:
            continue
        ed = np.array([abs(r["error_sobre_delta"]) for r in f])
        er = np.array([abs(r["error_rel"]) for r in f])
        resumen.append({"tramo": f"[2^{k},2^{k+1})", "n_muestras": len(f),
                        "max_|err|/delta": round(float(ed.max()), 4),
                        "media_|err|/delta": round(float(ed.mean()), 4),
                        "max_|err_rel|": round(float(er.max()), 5)})
    script = Path(__file__)
    datos = {
        "script": script.name,
        "script_sha256": hashlib.sha256(script.read_bytes()).hexdigest(),
        "fecha_utc": dt.datetime.now(dt.timezone.utc).isoformat(timespec="seconds"),
        "estado_epistemologico": "EVIDENCIA NUMÉRICA — no es demostración",
        "M": M, "c_phi": c, "resumen": resumen, "filas": filas,
    }
    ruta = guardar(f"li_liu_check_{M}", datos)
    print(f"c_phi = {c:.6f}")
    for r in resumen:
        print(f"  {r['tramo']:>16}  max|err|/δ = {r['max_|err|/delta']:>8}  "
              f"media = {r['media_|err|/delta']:>8}  max|err_rel| = {r['max_|err_rel|']}")
    print("->", ruta.name)
    return 0


if __name__ == "__main__":
    sys.exit(main())
