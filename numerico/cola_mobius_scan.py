#!/usr/bin/env python3
"""Exact finite scan of the primitive Mobius-von Mangoldt tail.

Provenance: LUMINA/OpenClaw, 2026-11-23.
Status: [N]. This is an experiment, not a proof of Goldbach.
"""

from __future__ import annotations

import argparse
import json
import math
from datetime import datetime, timezone
from pathlib import Path


def sieves(limit: int):
    spf = list(range(limit + 1))
    mu = [1] * (limit + 1)
    primes = []
    mu[0] = 0
    is_composite = bytearray(limit + 1)
    for n in range(2, limit + 1):
        if not is_composite[n]:
            primes.append(n)
            spf[n] = n
            mu[n] = -1
        for p in primes:
            v = n * p
            if v > limit:
                break
            is_composite[v] = 1
            spf[v] = p
            if n % p == 0:
                mu[v] = 0
                break
            mu[v] = -mu[n]

    lam = [0.0] * (limit + 1)
    for p in primes:
        lp = math.log(p)
        q = p
        while q <= limit:
            lam[q] = lp
            if q > limit // p:
                break
            q *= p
    return spf, mu, lam, primes


def divisors_from_spf(n: int, spf):
    fac = []
    while n > 1:
        p = spf[n]
        e = 0
        while n % p == 0:
            n //= p
            e += 1
        fac.append((p, e))
    divs = [1]
    for p, e in fac:
        old = divs[:]
        q = 1
        for _ in range(e):
            q *= p
            divs += [d * q for d in old]
    return divs


def singular_series(N: int, primes) -> float:
    c2 = 1.0
    for p in primes:
        if p > 200000:
            break
        if p > 2:
            c2 *= 1.0 - 1.0 / ((p - 1) ** 2)
    value = 2.0 * c2
    n = N
    for p in primes:
        if p * p > n:
            break
        if n % p == 0:
            if p > 2:
                value *= (p - 1) / (p - 2)
            while n % p == 0:
                n //= p
    if n > 2:
        value *= (n - 1) / (n - 2)
    return value


def exact_row(N: int, spf, mu, lam, primes):
    D = int(N ** (1.0 / 3.0))
    while (D + 1) ** 3 <= N:
        D += 1
    while D ** 3 > N:
        D -= 1

    tail = 0.0
    for n in range(2, N):
        ln = lam[n]
        if ln == 0.0:
            continue
        m = N - n
        for d in divisors_from_spf(m, spf):
            if d > D and mu[d] and math.gcd(m // d, N) == 1:
                tail += mu[d] * math.log(d) * ln

    R = sum(lam[n] * lam[N - n] for n in range(1, N))
    ss = singular_series(N, primes)
    main = N * ss
    return {
        "N": N,
        "D": D,
        "R": R,
        "singular_series": ss,
        "main": main,
        "T_prim": tail,
        "T_over_N": tail / N,
        "T_over_main": tail / main,
        "R_over_main": R / main,
    }


def exact_U(N: int, mu, lam):
    """Compute U_N(D) literally; O(N log N) finite diagnostic."""
    D = math.isqrt(math.isqrt(N))  # harmless initial value; corrected below
    D = int(N ** (1.0 / 3.0))
    while (D + 1) ** 3 <= N:
        D += 1
    while D ** 3 > N:
        D -= 1
    cmax = (N - 1) // (D + 1)
    partial = [0.0] * (cmax + 1)
    maxima = [0.0] * (cmax + 1)
    for d in range(D + 1, N):
        md = mu[d]
        if md == 0:
            continue
        for c in range(1, min(cmax, (N - 1) // d) + 1):
            if math.gcd(c, N) != 1:
                continue
            partial[c] += md * lam[N - c * d]
            a = abs(partial[c])
            if a > maxima[c]:
                maxima[c] = a
    U = sum(maxima)
    return {
        "N": N,
        "D": D,
        "U": U,
        "U_logN_over_N": U * math.log(N) / N,
        "U_over_N": U / N,
    }


def sample_points(limit: int, count: int):
    pts = {4, 6, 10, 100, 1000, limit if limit % 2 == 0 else limit - 1}
    lo = math.log(10)
    hi = math.log(limit)
    for i in range(count):
        x = int(round(math.exp(lo + (hi - lo) * i / max(1, count - 1))))
        pts.add(x + (x & 1))
    return sorted(n for n in pts if 4 <= n <= limit and n % 2 == 0)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--hasta", type=int, default=1_000_000)
    ap.add_argument("--muestras", type=int, default=80)
    ap.add_argument("--u", action="store_true", help="also compute exact U at powers of ten")
    args = ap.parse_args()
    spf, mu, lam, primes = sieves(args.hasta)
    rows = [exact_row(N, spf, mu, lam, primes) for N in sample_points(args.hasta, args.muestras)]
    stamp = datetime.now(timezone.utc).strftime("%Y%m%d_%H%M%S")
    outdir = Path(__file__).resolve().parent / "resultados"
    outdir.mkdir(exist_ok=True)
    out = outdir / f"{stamp}_cola_mobius_{args.hasta}.json"
    payload = {
        "provenance": "LUMINA/OpenClaw; exact finite computation; not a proof",
        "generated_utc": datetime.now(timezone.utc).isoformat(),
        "limit": args.hasta,
        "samples": len(rows),
        "max_T_over_main": max(rows, key=lambda r: r["T_over_main"]),
        "min_T_over_main": min(rows, key=lambda r: r["T_over_main"]),
        "max_abs_T_over_N": max(rows, key=lambda r: abs(r["T_over_N"])),
        "rows": rows,
    }
    if args.u:
        u_points = []
        p = 100
        while p <= args.hasta:
            u_points.append(p if p % 2 == 0 else p - 1)
            p *= 10
        payload["U_rows"] = [exact_U(N, mu, lam) for N in u_points]
    out.write_text(json.dumps(payload, indent=2) + "\n", encoding="utf-8")
    print(json.dumps({k: v for k, v in payload.items() if k != "rows"}, indent=2))
    print(out)


if __name__ == "__main__":
    main()
