#!/usr/bin/env python3
"""Finite diagnostics for Mobius sums in reduced/non-reduced residue classes.

Provenance: LUMINA/OpenClaw, 2026-11-24. Status [N].
This tests the exact conductor-lowering identity (9.1) and measures finite
Bombieri-Vinogradov-style maxima. It is not a proof of an asymptotic theorem.
"""

from __future__ import annotations

import argparse
import json
import math
import random
from datetime import datetime, timezone
from pathlib import Path


def mobius_sieve(n: int):
    mu = [0] * (n + 1)
    mu[1] = 1
    primes = []
    composite = bytearray(n + 1)
    for i in range(2, n + 1):
        if not composite[i]:
            primes.append(i)
            mu[i] = -1
        for p in primes:
            v = i * p
            if v > n:
                break
            composite[v] = 1
            if i % p == 0:
                mu[v] = 0
                break
            mu[v] = -mu[i]
    return mu


def prime_support(n: int):
    out = []
    p = 2
    while p * p <= n:
        if n % p == 0:
            out.append(p)
            while n % p == 0:
                n //= p
        p += 1
    if n > 1:
        out.append(n)
    return out


def smooth_numbers(primes, limit: int):
    vals = [1]
    for p in primes:
        old = vals[:]
        q = p
        while q <= limit:
            vals.extend(v * q for v in old if v * q <= limit)
            if q > limit // p:
                break
            q *= p
    return sorted(set(vals))


def lhs(mu, y: int, q: int, r: int):
    return sum(mu[n] for n in range(1, y + 1) if n % q == r % q)


def rhs_91(mu, y: int, q: int, r: int):
    g = math.gcd(r, q)
    # If g is not squarefree, every n in the class is divisible by a square.
    if any((g // p) % p == 0 for p in prime_support(g)):
        return 0
    h, b = q // g, r // g
    total = 0
    for k in smooth_numbers(prime_support(g), y // max(1, g)):
        if math.gcd(k, h) != 1:
            continue
        inv = pow(k, -1, h) if h > 1 else 0
        a = (b * inv) % h if h > 1 else 0
        total += sum(
            mu[u]
            for u in range(1, y // (g * k) + 1)
            if h == 1 or u % h == a
        )
    return mu[g] * total


def maxima_for_q(mu, x: int, q: int):
    cur = [0] * q
    mx = [0] * q
    for n in range(1, x + 1):
        a = n % q
        cur[a] += mu[n]
        if abs(cur[a]) > mx[a]:
            mx[a] = abs(cur[a])
    reduced = max((mx[a] for a in range(q) if math.gcd(a, q) == 1), default=0)
    all_classes = max(mx, default=0)
    winner = max(range(q), key=mx.__getitem__)
    return reduced, all_classes, winner, math.gcd(winner, q)


def run(x: int, qmax: int, trials: int, seed: int):
    mu = mobius_sieve(x)
    rng = random.Random(seed)
    identity_checks = 0
    for _ in range(trials):
        q = rng.randint(1, qmax)
        r = rng.randrange(q)
        y = rng.randint(1, x)
        a, b = lhs(mu, y, q, r), rhs_91(mu, y, q, r)
        if a != b:
            raise AssertionError((x, q, r, y, a, b))
        identity_checks += 1

    e_red = 0
    e_all = 0
    nonreduced_winners = 0
    rows = []
    for q in range(1, qmax + 1):
        red, allc, winner, g = maxima_for_q(mu, x, q)
        e_red += red
        e_all += allc
        nonreduced_winners += int(g > 1)
        rows.append({"q": q, "max_reduced": red, "max_all": allc,
                     "winner": winner, "winner_gcd": g})
    return {
        "x": x,
        "qmax": qmax,
        "identity_checks": identity_checks,
        "E_reduced": e_red,
        "E_all": e_all,
        "E_all_over_reduced": e_all / e_red if e_red else None,
        "E_reduced_logx_over_x": e_red * math.log(x) / x,
        "E_all_logx_over_x": e_all * math.log(x) / x,
        "nonreduced_winner_fraction": nonreduced_winners / qmax,
        "rows": rows,
    }


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--hasta", type=int, nargs="+", default=[10_000, 100_000])
    ap.add_argument("--q-factor", type=float, default=0.25,
                    help="qmax=floor(q_factor*sqrt(x))")
    ap.add_argument("--checks", type=int, default=500)
    ap.add_argument("--seed", type=int, default=20261124)
    args = ap.parse_args()
    data = []
    for x in args.hasta:
        qmax = max(1, int(args.q_factor * math.sqrt(x)))
        data.append(run(x, qmax, args.checks, args.seed + x))
    payload = {
        "provenance": "LUMINA/OpenClaw; exact finite diagnostic; not a proof",
        "created_utc": datetime.now(timezone.utc).isoformat(),
        "parameters": vars(args),
        "scales": data,
    }
    outdir = Path(__file__).resolve().parent / "resultados"
    outdir.mkdir(exist_ok=True)
    stamp = datetime.now(timezone.utc).strftime("%Y%m%d_%H%M%S")
    out = outdir / f"{stamp}_mobius_clases_no_reducidas.json"
    out.write_text(json.dumps(payload, indent=2) + "\n", encoding="utf-8")
    print(json.dumps({"result": str(out), "scales": [
        {k: v for k, v in row.items() if k != "rows"} for row in data
    ]}, indent=2))


if __name__ == "__main__":
    main()
