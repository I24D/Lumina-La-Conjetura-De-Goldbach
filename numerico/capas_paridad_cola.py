#!/usr/bin/env python3
"""Decompose the primitive tail by omega(d), exactly on finite samples.

Provenance: LUMINA/OpenClaw, 2026-11-24. Status [N], not a proof.
"""

from __future__ import annotations

import argparse
import json
import math
from datetime import datetime, timezone
from pathlib import Path

import cola_mobius_scan as base


def omega_squarefree(d: int, spf) -> int:
    j = 0
    while d > 1:
        p = spf[d]
        j += 1
        d //= p
        if d % p == 0:
            return -1
    return j


def row(N, spf, mu, lam, primes):
    D = int(N ** (1 / 3))
    while (D + 1) ** 3 <= N:
        D += 1
    while D ** 3 > N:
        D -= 1
    layers = {}
    for n in range(2, N):
        if lam[n] == 0.0:
            continue
        m = N - n
        for d in base.divisors_from_spf(m, spf):
            if d <= D or not mu[d] or math.gcd(m // d, N) != 1:
                continue
            j = omega_squarefree(d, spf)
            assert j >= 1 and mu[d] == (-1) ** j
            layers[j] = layers.get(j, 0.0) + math.log(d) * lam[n]
    signed = sum(((-1) ** j) * v for j, v in layers.items())
    odd = sum(v for j, v in layers.items() if j % 2)
    even = sum(v for j, v in layers.items() if j % 2 == 0)
    main = N * base.singular_series(N, primes)
    return {
        "N": N,
        "D": D,
        "layers": {str(j): v for j, v in sorted(layers.items())},
        "odd_mass": odd,
        "even_mass": even,
        "T": signed,
        "T_over_main": signed / main,
        "even_over_odd": even / odd if odd else None,
        "j1_over_main": layers.get(1, 0.0) / main,
        "j2_over_main": layers.get(2, 0.0) / main,
    }


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--hasta", type=int, default=1_000_000)
    ap.add_argument("--muestras", type=int, default=80)
    args = ap.parse_args()
    spf, mu, lam, primes = base.sieves(args.hasta)
    rows = [row(N, spf, mu, lam, primes)
            for N in base.sample_points(args.hasta, args.muestras)]
    payload = {
        "provenance": "LUMINA/OpenClaw; exact finite diagnostic; not a proof",
        "generated_utc": datetime.now(timezone.utc).isoformat(),
        "limit": args.hasta,
        "samples": len(rows),
        "max_T": max(rows, key=lambda x: x["T_over_main"]),
        "min_T": min(rows, key=lambda x: x["T_over_main"]),
        "max_even_over_odd": max(rows, key=lambda x: x["even_over_odd"] or -1),
        "rows": rows,
    }
    outdir = Path(__file__).resolve().parent / "resultados"
    outdir.mkdir(exist_ok=True)
    stamp = datetime.now(timezone.utc).strftime("%Y%m%d_%H%M%S")
    out = outdir / f"{stamp}_capas_paridad_cola_{args.hasta}.json"
    out.write_text(json.dumps(payload, indent=2) + "\n", encoding="utf-8")
    print(json.dumps({k: v for k, v in payload.items() if k != "rows"}, indent=2))
    print(out)


if __name__ == "__main__":
    main()
