#!/usr/bin/env python3
"""
mobius_goldbach_aleatoriedad.py — experimento [N] para matematicas/07.

Procedencia: Claude Code (claude-opus-5-5), 2026-11-23.

Pregunta: la correlación Möbius–Goldbach
    M1(N) = sum_{n<N} mu(n) log(n) Lambda(N-n)     (caja c=1 de 05/06, sin cortar en D)
    M0(N) = sum_{n<N} mu(n) Lambda(N-n)
¿tiene el tamaño que predice el modelo aleatorio ("mu = signos independientes sobre
los libres de cuadrados") o hay cancelación/conspiración adicional?

Modelo aleatorio: Var[M1(N)] = V1(N) = sum mu(n)^2 log^2(n) Lambda(N-n)^2, y
análogamente V0. Se estudia Z(N) = M(N)/sqrt(V(N)) para N par en (X/2, X].
Si mu fuera "aleatoria" frente a Lambda(N-.), Z ~ N(0,1): media 0 y varianza 1.

Todas las convoluciones se hacen a la vez para todo N <= X por FFT (float64).
Se comprueba la exactitud contra sumas directas en una muestra de N.

ESTADO EPISTEMOLÓGICO: evidencia numérica [N]; no demuestra nada asintótico.
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


def conv(f: np.ndarray, g: np.ndarray, M: int) -> np.ndarray:
    tam = 1 << (2 * M + 2).bit_length()
    F = np.fft.rfft(f, tam)
    G = np.fft.rfft(g, tam)
    return np.fft.irfft(F * G, tam)[: M + 1]


def momentos(z: np.ndarray) -> dict:
    m = float(z.mean())
    v = float(z.var())
    s = float(((z - m) ** 3).mean() / v ** 1.5)
    k = float(((z - m) ** 4).mean() / v ** 2)
    return {"n": int(z.size), "media": m, "varianza": v, "asimetria": s,
            "curtosis": k, "frac_|z|>3": float((np.abs(z) > 3).mean()),
            "max_|z|": float(np.abs(z).max())}


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--hasta", type=float, default=1e7)
    args = ap.parse_args()
    X = int(args.hasta)
    mu, lam, _ = mobius_y_lambda(X)
    n = np.arange(X + 1, dtype=np.float64)
    logn = np.zeros(X + 1)
    logn[1:] = np.log(n[1:])
    muf = mu.astype(np.float64)
    f1 = muf * logn
    f0 = muf.copy()
    f0[0] = 0.0
    q1 = muf ** 2 * logn ** 2
    q0 = muf ** 2
    q0[0] = 0.0
    lam2 = lam ** 2

    M1 = conv(f1, lam, X)
    M0 = conv(f0, lam, X)
    V1 = conv(q1, lam2, X)
    V0 = conv(q0, lam2, X)
    R = conv(lam, lam, X)

    # comprobación de exactitud contra sumas directas
    rng = np.random.default_rng(7)
    err = 0.0
    for N in rng.integers(X // 2, X, 12).tolist():
        d = np.arange(1, N)
        m1 = float(np.sum(f1[d] * lam[N - d]))
        err = max(err, abs(m1 - M1[N]) / max(1.0, abs(m1)))
    print(f"error relativo máximo FFT vs directo (12 muestras): {err:.2e}")

    pares = np.arange(X // 2 + (X // 2) % 2, X + 1, 2)
    impares = np.arange(X // 2 + 1 - (X // 2) % 2, X + 1, 2)
    res = {}
    for nombre, M, V in (("M1_mu_log", M1, V1), ("M0_mu", M0, V0)):
        for etiqueta, Ns in (("N_par", pares), ("N_impar", impares)):
            z = M[Ns] / np.sqrt(V[Ns])
            res[f"{nombre}|{etiqueta}"] = momentos(z)
        z = M[pares] / np.sqrt(V[pares])
        # sesgo según N mod 3 y N mod 5, y correlación con S(N) (muestra)
        por_mod = {}
        for q in (3, 5, 7):
            for r in range(q):
                sel = (pares % q) == r
                por_mod[f"N%{q}=={r}"] = float(z[sel].mean())
        res[f"{nombre}|N_par|media_por_clase"] = por_mod
        muestra = pares[:: max(1, pares.size // 5000)]
        S = np.array([serie_singular(int(N)) for N in muestra])
        zz = M[muestra] / np.sqrt(V[muestra])
        res[f"{nombre}|N_par|corr_con_S(N)"] = float(np.corrcoef(zz, S)[0, 1])
        res[f"{nombre}|N_par|max_|M|/N"] = float(np.max(np.abs(M[pares]) / pares))

    # referencia: el error de Hardy–Littlewood en la misma ventana (muestra)
    muestra = pares[:: max(1, pares.size // 5000)]
    S = np.array([serie_singular(int(N)) for N in muestra])
    err_hl = (R[muestra] - muestra * S)
    res["R_menos_NS|muestra"] = {"media_rel": float(np.mean(err_hl / (muestra * S))),
                                 "max_|rel|": float(np.max(np.abs(err_hl) / (muestra * S))),
                                 "rms/sqrt(N log N)": float(np.sqrt(np.mean(err_hl ** 2 / (muestra * np.log(muestra)))))}

    for k, v in res.items():
        print(k, v)
    script = Path(__file__)
    datos = {"script": script.name,
             "script_sha256": hashlib.sha256(script.read_bytes()).hexdigest(),
             "fecha_utc": dt.datetime.now(dt.timezone.utc).isoformat(timespec="seconds"),
             "estado_epistemologico": "EVIDENCIA NUMÉRICA — no es demostración",
             "X": X, "ventana": [X // 2, X], "error_fft": err, "resultados": res}
    print("->", guardar(f"mobius_goldbach_aleatoriedad_{X}", datos).name)
    return 0


if __name__ == "__main__":
    sys.exit(main())
