#!/usr/bin/env python3
"""Barrido de las familias de serie singular minima (doc 15, s5).

Pregunta pre-registrada (15 s5): el supremo global de rho(N)
(0.4201 en N=30916, run denso 20261124_224653) vive en las familias
delgadas de G(N) casi minima?

Familias pre-registradas:
    2^a (a=2..19, N<=10^6), 2p, 4p (primos impares, N<=10^5),
    2p^2 (p<=223, N<=10^5).

Brazo suplementario POST-OBSERVACIONAL (anadido tras ver que el supremo
global esta en 30916=2^2*59*131, que no esta en las familias estrictas):
todos los N pares en (25000, 40000], exhaustivo.

Veredicto binario en el JSON: si sup(familia) >= 0.4201 para alguna
familia, las familias capturan el peor caso; si no, el supremo global
esta fuera y la reduccion a familias delgadas queda debilitada.

Procedencia: LUMINA/OpenClaw, 2026-11-24. [N]. No es demostracion.
Reutiliza exact_row de cola_mobius_scan.py (D=floor(N^(1/3)) fijo).
"""

from __future__ import annotations

import json
import subprocess
import sys
import time
from datetime import datetime, timezone
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
from cola_mobius_scan import exact_row, sieves  # noqa: E402

HERE = Path(__file__).resolve().parent
GLOBAL_SUP_DENSE = 0.42014005735698656  # run denso 20261124_224653


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


def build_batch(primes):
    fams = {}
    fams["2^a"] = [2 ** a for a in range(2, 20) if 2 ** a <= 10 ** 6]
    fams["2p"] = [2 * p for p in primes if p > 2 and 2 * p <= 10 ** 5]
    fams["4p"] = [4 * p for p in primes if p > 2 and 4 * p <= 10 ** 5]
    fams["2p^2"] = [
        2 * p * p for p in primes if p > 2 and 2 * p * p <= 10 ** 5
    ]
    # suplementario post-observacional: ventana exhaustiva del peor cluster
    fams["ventana(25000,40000]"] = list(range(25002, 40001, 2))
    return fams


def stats(rows):
    if not rows:
        return None
    rhos = sorted((r["rho"], r["N"]) for r in rows)
    n = len(rows)
    return {
        "count": n,
        "sup_rho": rhos[-1][0],
        "sup_rho_N": rhos[-1][1],
        "top5": [{"N": nn, "rho": rr} for rr, nn in reversed(rhos[-5:])],
        "median_rho": rhos[n // 2][0],
        "min_rho": rhos[0][0],
        "min_rho_N": rhos[0][1],
        "max_gap": max(r["gap"] for r in rows),
        "min_r_over_main": min(r["r_over_main"] for r in rows),
        "captura_sup_global": rhos[-1][0] >= GLOBAL_SUP_DENSE,
    }


def main() -> int:
    t0 = time.time()
    spf, mu, lam, primes = sieves(10 ** 6)
    t_sieve = time.time() - t0

    fams = build_batch(primes)
    fam_stats = {}
    total = sum(len(v) for v in fams.values())
    done = 0
    t1 = time.time()
    for name, ns in fams.items():
        rows = []
        for N in ns:
            row = exact_row(N, spf, mu, lam, primes)
            rho = row["T_over_main"]
            r_over = row["R_over_main"]
            rows.append(
                {
                    "N": N,
                    "rho": rho,
                    "r_over_main": r_over,
                    "gap": 1.0 - r_over - rho,
                }
            )
            done += 1
            if done % 1000 == 0:
                el = time.time() - t1
                print(
                    f"progreso {done}/{total} ({el:.0f}s)", flush=True
                )
        fam_stats[name] = stats(rows)
    t_scan = time.time() - t1

    any_capture = any(
        (s or {}).get("captura_sup_global", False)
        for s in fam_stats.values()
        if s
    )
    out = {
        "schema": "goldbach_lumina.deficit_familias.v1",
        "date": datetime.now(timezone.utc).isoformat(),
        "git_rev": git_rev(),
        "agent": "LUMINA/OpenClaw (telegram)",
        "status": "[N] experimento; no es demostracion",
        "pregunta_pre-registrada": (
            "sup rho de las familias de G minima >= 0.4201 "
            "(sup global denso)?"
        ),
        "global_sup_dense": GLOBAL_SUP_DENSE,
        "brazo_post-observacional": (
            "ventana(25000,40000]: anadida tras ver 30916=2^2*59*131; "
            "NO estaba en el pre-registro"
        ),
        "familias": fam_stats,
        "alguna_familia_captura": any_capture,
        "timings_s": {
            "sieve": round(t_sieve, 2),
            "scan": round(t_scan, 2),
        },
    }
    out_path = (
        HERE
        / "resultados"
        / (
            datetime.now(timezone.utc).strftime("%Y%m%d_%H%M%S")
            + "_deficit_familias.json"
        )
    )
    out_path.write_text(
        json.dumps(out, indent=1, ensure_ascii=False), encoding="utf-8"
    )
    print(json.dumps(fam_stats, indent=1, ensure_ascii=False))
    print("captura_global:", any_capture, "->", str(out_path))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
