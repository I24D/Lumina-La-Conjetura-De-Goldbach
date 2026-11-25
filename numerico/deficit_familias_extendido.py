#!/usr/bin/env python3
"""Pre-registered extension of the thin-family rho scan to N <= 10^6.

Hypothesis fixed before execution (2026-11-24 19:1x ET): in the thin
families N=2p and N=4p, rho(N)=T_prim/(N*G(N)) remains below 0.75 and
does not show a rise of more than 0.15 between the first and last bands.

Refutation criterion: either (i) any sampled rho >= 0.75, or (ii) for
one family sup(last band)-sup(first band) > 0.15.  Passing is numerical
evidence only; it is not a proof or an eventual uniform bound.

Design: deterministic stratified sample (seed 20261124), 100 primes per
nonempty family/band, bands (1e5,3e5], (3e5,6e5], (6e5,1e6].
"""

from __future__ import annotations

import json
import random
import subprocess
import sys
import time
from datetime import datetime, timezone
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
from cola_mobius_scan import exact_row, sieves  # noqa: E402

HERE = Path(__file__).resolve().parent
SEED = 20261124
BANDS = [(100_000, 300_000), (300_000, 600_000), (600_000, 1_000_000)]
PER_CELL = 100
RHO_REFUTE = 0.75
DRIFT_REFUTE = 0.15


def git_rev() -> str:
    try:
        return subprocess.check_output(
            ["git", "rev-parse", "--short", "HEAD"], cwd=HERE.parent
        ).decode().strip()
    except Exception:
        return "unknown"


def cell_stats(rows):
    vals = sorted((r["rho"], r["N"]) for r in rows)
    return {
        "count": len(vals),
        "sup_rho": vals[-1][0],
        "sup_rho_N": vals[-1][1],
        "median_rho": vals[len(vals) // 2][0],
        "min_rho": vals[0][0],
        "top5": [{"N": n, "rho": x} for x, n in reversed(vals[-5:])],
    }


def main() -> int:
    t0 = time.time()
    spf, mu, lam, primes = sieves(1_000_000)
    rng = random.Random(SEED)
    results = {}
    all_rows = []
    for multiplier, family in [(2, "2p"), (4, "4p")]:
        family_cells = []
        for lo, hi in BANDS:
            candidates = [p for p in primes if p > 2 and lo < multiplier * p <= hi]
            chosen = sorted(rng.sample(candidates, min(PER_CELL, len(candidates))))
            rows = []
            for p in chosen:
                N = multiplier * p
                raw = exact_row(N, spf, mu, lam, primes)
                row = {"N": N, "p": p, "rho": raw["T_over_main"]}
                rows.append(row)
                all_rows.append({"family": family, **row})
            family_cells.append({"band": [lo, hi], **cell_stats(rows)})
            print(f"{family} {lo}-{hi}: n={len(rows)} sup={family_cells[-1]['sup_rho']:.6f}", flush=True)
        results[family] = family_cells

    rho_max = max(r["rho"] for r in all_rows)
    rho_arg = max(all_rows, key=lambda r: r["rho"])
    drifts = {
        family: cells[-1]["sup_rho"] - cells[0]["sup_rho"]
        for family, cells in results.items()
    }
    refuted = rho_max >= RHO_REFUTE or any(x > DRIFT_REFUTE for x in drifts.values())
    out = {
        "schema": "goldbach_lumina.deficit_familias_extendido.v1",
        "date": datetime.now(timezone.utc).isoformat(),
        "agent": "LUMINA/OpenClaw; 30-minute loop",
        "status": "[N] experiment; not a proof",
        "git_rev": git_rev(),
        "preregistration": {
            "seed": SEED,
            "bands": BANDS,
            "per_family_band": PER_CELL,
            "hypothesis": "rho<0.75 and last-minus-first band sup <=0.15 in 2p and 4p",
            "refutation": "any rho>=0.75 or either family drift>0.15",
        },
        "families": results,
        "global_sup": {"rho": rho_max, **rho_arg},
        "drift_last_minus_first": drifts,
        "refuted": refuted,
        "elapsed_s": round(time.time() - t0, 2),
    }
    out_path = HERE / "resultados" / (
        datetime.now(timezone.utc).strftime("%Y%m%d_%H%M%S")
        + "_deficit_familias_extendido.json"
    )
    out_path.write_text(json.dumps(out, indent=1, ensure_ascii=False), encoding="utf-8")
    print(json.dumps({"global_sup": out["global_sup"], "drifts": drifts, "refuted": refuted, "elapsed_s": out["elapsed_s"]}, indent=1))
    print(out_path)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
