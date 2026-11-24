#!/usr/bin/env python3
"""Punto de mira unilateral: deficit puntual de la cola primitiva.

Cuantifica, N por N (par), el cociente unilateral

    rho(N) = T_N^prim(D) / (N * G(N)),   D = floor(N^(1/3)),

donde T_N^prim(D) es la cola primitiva de (4.1)/(4.15) del doc 04
(restriccion (c,N)=1), G(N) la serie singular, y se registra tambien

    r(N) = R(N) / (N * G(N)),   gap(N) = 1 - r(N) - rho(N).

gap debe ser pequeno: valida numericamente la cadena R = S - T (4.1) +
S ~ N*G (4.15) en cada N probado.

Pregunta falsable (doc 14 seccion 5 / doc 04 linea 249, hipotesis H_eps):
    sup_{N par} rho(N)  ? esta acotado lejos de 1 y sin deriva al crecer N.
Si sup rho -> 1 en alguna escala, la ruta unilateral esta muerta en esa
escala. Si sup rho <= 1 - delta con delta visible y estable, el objetivo
analitico "T_prim <= (1-eps) N G(N)" queda calibrado con eps realista.

Procedencia: LUMINA/OpenClaw, 2026-11-24. Estado [N]: experimento, no prueba.
Reutiliza sieves/exact_row de cola_mobius_scan.py (validado en runs previos).
"""

from __future__ import annotations

import argparse
import json
import math
import random
import subprocess
import sys
import time
from datetime import datetime, timezone
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
from cola_mobius_scan import exact_row, sieves  # noqa: E402

HERE = Path(__file__).resolve().parent


def git_rev() -> str:
    try:
        return (
            subprocess.check_output(
                ["git", "rev-parse", "--short", "HEAD"], cwd=HERE.parent
            )
            .decode()
            .strip()
        )
    except Exception:
        return "unknown"


def build_batch(seed: int, dense: bool = False):
    """Even N: exhaustive low range + seeded random samples per scale."""
    rng = random.Random(seed)
    batch = []
    batch.extend(range(4, 2002, 2))  # all even N in [4, 2000]
    if dense:
        scales = [
            (2000, 20000, 600),
            (20000, 200000, 400),
            (200000, 1000000, 300),
        ]
    else:
        scales = [
            (2000, 20000, 300),
            (20000, 200000, 150),
            (200000, 1000000, 40),
        ]
    for lo, hi, count in scales:
        picks = set()
        while len(picks) < count:
            n = 2 * rng.randrange(lo // 2 + 1, hi // 2 + 1)
            if lo < n <= hi:
                picks.add(n)
        batch.extend(sorted(picks))
    return batch


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--seed", type=int, default=20261124)
    ap.add_argument("--limit", type=int, default=1000000)
    ap.add_argument("--dense", action="store_true")
    ap.add_argument("--out", default="")
    args = ap.parse_args()

    t0 = time.time()
    spf, mu, lam, primes = sieves(args.limit)
    t_sieve = time.time() - t0

    batch = build_batch(args.seed, args.dense)
    rows = []
    t1 = time.time()
    for N in batch:
        row = exact_row(N, spf, mu, lam, primes)
        rho = row["T_over_main"]
        r_over = row["R_over_main"]
        rows.append(
            {
                "N": N,
                "D": row["D"],
                "rho": rho,
                "r_over_main": r_over,
                "gap": 1.0 - r_over - rho,
                "T_over_N": row["T_over_N"],
            }
        )
    t_scan = time.time() - t1

    # Estadisticas por escala y globales (solo N pares; batch ya es par).
    def stats(sel):
        if not sel:
            return None
        rhos = sorted((r["rho"], r["N"]) for r in sel)
        n = len(sel)
        return {
            "count": n,
            "sup_rho": rhos[-1][0],
            "sup_rho_N": rhos[-1][1],
            "top10": [
                {"N": nn, "rho": rr} for rr, nn in reversed(rhos[-10:])
            ],
            "median_rho": rhos[n // 2][0],
            "min_rho": rhos[0][0],
            "min_rho_N": rhos[0][1],
            "max_gap": max(r["gap"] for r in sel),
            "min_r_over_main": min(r["r_over_main"] for r in sel),
            "min_r_over_main_N": min(
                sel, key=lambda r: r["r_over_main"]
            )["N"],
        }

    bands = {
        "0004-2000": [r for r in rows if r["N"] <= 2000],
        "2000-20000": [r for r in rows if 2000 < r["N"] <= 20000],
        "20000-200000": [r for r in rows if 20000 < r["N"] <= 200000],
        "200000-1000000": [r for r in rows if 200000 < r["N"] <= 1000000],
    }
    out = {
        "schema": "goldbach_lumina.deficit_unilateral.v1",
        "date": datetime.now(timezone.utc).isoformat(),
        "seed": args.seed,
        "dense": bool(args.dense),
        "limit": args.limit,
        "git_rev": git_rev(),
        "agent": "LUMINA/OpenClaw (telegram)",
        "status": "[N] experimento; no es demostracion",
        "identity": "rho(N)=T_prim(D)/(N*G(N)), D=floor(N^(1/3))",
        "falsable": "sup rho lejos de 1 y sin deriva => H_eps calibrada; sup rho -> 1 => ruta muerta",
        "timings_s": {"sieve": round(t_sieve, 2), "scan": round(t_scan, 2)},
        "bands": {k: stats(v) for k, v in bands.items()},
        "global": stats(rows),
    }

    out_path = (
        Path(args.out)
        if args.out
        else HERE
        / "resultados"
        / (
            datetime.now(timezone.utc).strftime("%Y%m%d_%H%M%S")
            + "_deficit_unilateral_"
            + ("dense_" if args.dense else "")
            + str(args.limit)
            + ".json"
        )
    )
    out_path.write_text(
        json.dumps(out, indent=1, ensure_ascii=False), encoding="utf-8"
    )
    print(json.dumps({k: out[k] for k in ("bands", "global", "timings_s")}, indent=1))
    print("rows:", len(rows), "->", str(out_path))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
