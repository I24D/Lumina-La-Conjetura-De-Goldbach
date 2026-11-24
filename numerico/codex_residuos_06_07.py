#!/usr/bin/env python3
"""Codex, 2026-11-23 ET; model identifier not exposed. Status [P]/[N].

Sources: project notes 04--07. Finite diagnostics, never a Goldbach proof.
Methods: FFT, independent Ramanujan divisor formula, prefix sums by residue,
exact tail identity S-R-Tbad checked against direct (c,d) loops.
PQ residual != a continuous minor-arc integral. A separate discrete Fourier
projection uses explicitly recorded arcs and grid. No models/services changed.
"""
from __future__ import annotations

import argparse
from datetime import datetime, timezone
import gc
import hashlib
import json
import math
from pathlib import Path
import platform
import time

import numpy as np

from cajas_bilineales import mobius_y_lambda
from goldbach import serie_singular


def factor(n):
    out = {}
    p = 2
    while p * p <= n:
        while n % p == 0:
            out[p] = out.get(p, 0) + 1
            n //= p
        p += 1
    if n > 1:
        out[n] = out.get(n, 0) + 1
    return out


def divisors(n):
    ans = [1]
    for p, exponent in factor(n).items():
        old = ans[:]
        for k in range(1, exponent + 1):
            ans.extend(d * p ** k for d in old)
    return ans


def mu_int(n):
    f = factor(n)
    return 0 if any(e > 1 for e in f.values()) else (-1) ** len(f)


def phi(n):
    ans = n
    for p in factor(n):
        ans = ans // p * (p - 1)
    return ans


def cube_root(n):
    d = int(n ** (1 / 3))
    while (d + 1) ** 3 <= n:
        d += 1
    while d ** 3 > n:
        d -= 1
    return d


def coefficients(qmax, selected):
    """PQ=sum_h a_h(Q)*M(N-1;h,N), exact finite Ramanujan identity."""
    a = np.zeros(qmax + 1)
    result = {}
    for q in range(1, qmax + 1):
        if mu_int(q):
            for h in divisors(q):
                a[h] += h * mu_int(h) / phi(q)
        if q in selected:
            result[q] = a.copy()
    return result


def residue_prefix(values, h):
    out = np.empty(values.size, dtype=values.dtype)
    for r in range(h):
        out[r::h] = np.cumsum(values[r::h], dtype=values.dtype)
    return out


def stats(z, block=2048):
    z = np.asarray(z, dtype=float)
    mean = float(z.mean())
    variance = float(z.var())
    center = z - mean
    nblocks = z.size // block
    means = z[:nblocks * block].reshape(nblocks, block).mean(axis=1) if nblocks else np.array([])
    return {"n": int(z.size), "mean": mean, "variance": variance,
            "skewness": float(np.mean(center ** 3) / variance ** 1.5) if variance else None,
            "kurtosis": float(np.mean(center ** 4) / variance ** 2) if variance else None,
            "fraction_positive": float(np.mean(z > 0)),
            "positive_part_mean": float(np.maximum(z, 0).mean()),
            "negative_part_mean": float(np.maximum(-z, 0).mean()),
            "positive_tail_gt_2": float(np.mean(z > 2)),
            "negative_tail_lt_minus2": float(np.mean(z < -2)),
            "positive_tail_gt_3": float(np.mean(z > 3)),
            "negative_tail_lt_minus3": float(np.mean(z < -3)),
            "quantiles": dict(zip(["p001", "p01", "p05", "p50", "p95", "p99", "p999"],
                                  map(float, np.quantile(z, [.001, .01, .05, .5, .95, .99, .999])))),
            "min": float(z.min()), "max": float(z.max()),
            "block_size_in_even_observations": block,
            "block_mean_quantiles": list(map(float, np.quantile(means, [0, .05, .5, .95, 1]))) if nblocks else [],
            "fraction_block_means_positive": float(np.mean(means > 0)) if nblocks else None,
            "inference": "descriptive only; dependent observations, no iid p-values"}


def classes(z, ns):
    return {str(q): {str(r): float(z[ns % q == r].mean()) for r in range(q)}
            for q in (3, 5, 7)}


def major_mask(nfft, x, qmax):
    mask = np.zeros(nfft // 2 + 1, dtype=bool)
    for q in range(1, qmax + 1):
        for a in range(q // 2 + 1):
            if math.gcd(a, q) != 1:
                continue
            center = a * nfft / q
            radius = nfft * qmax / (q * x)
            lo = max(0, math.ceil(center - radius))
            hi = min(nfft // 2, math.floor(center + radius))
            mask[lo:hi + 1] = True
    return mask


def mangoldt_convolution(mu, lam, x, qmax):
    nfft = 1 << (2 * x + 2).bit_length()
    f = np.fft.rfft(mu.astype(float), nfft)
    g = np.fft.rfft(lam, nfft)
    f *= g
    del g
    total = np.fft.irfft(f, nfft)[:x + 1].copy()
    mask = major_mask(nfft, x, qmax)
    f[~mask] = 0
    projected = np.fft.irfft(f, nfft)[:x + 1].copy()
    del f
    f = np.fft.rfft(mu.astype(float) ** 2, nfft)
    g = np.fft.rfft(lam ** 2, nfft)
    f *= g
    del g
    variance = np.fft.irfft(f, nfft)[:x + 1].copy()
    del f
    return total, variance, projected, {"nfft": nfft, "Q": qmax,
        "major_arc_half_width": "Q/(q*X), reduced a/q; symmetric rfft mask",
        "number_positive_frequency_bins": int(mask.sum()),
        "interpretation": "exact discrete grid projection; no certified continuous-integral error bound"}


def tail_bad(n, d0, mu):
    ans = 0.0
    for p in factor(n):
        prime_power = p
        while prime_power < n:
            m = n - prime_power
            for d in divisors(m):
                if d > d0 and math.gcd(m // d, n) > 1:
                    ans += int(mu[d]) * math.log(d) * math.log(p)
            prime_power *= p
    return ans


def tail_identity(n, mu, lam, beta, primes=None):
    d0 = cube_root(n)
    small = -math.fsum(float(beta[d]) * float(lam[n - d:0:-d].sum())
                       for d in range(1, d0 + 1) if mu[d])
    r = float(np.dot(lam[1:n], lam[n - 1:0:-1]))
    bad = tail_bad(n, d0, mu)
    t = small - r - bad
    if primes is None:
        candidates = np.flatnonzero((mu[1:n] == -1) & (lam[1:n] > 0)) + 1
        primes = np.array([int(p) for p in candidates if factor(int(p)) == {int(p): 1}], dtype=np.int64)
    else:
        primes = primes[primes < n]
    above = primes[primes > d0]
    below = primes[primes <= d0]
    r_above = float(np.dot(np.log(above), lam[n - above]))
    r_below = float(np.dot(np.log(below), lam[n - below]))
    r_powers = r - r_above - r_below
    ss = serie_singular(n)
    return {"N": n, "D": d0, "singular_series": ss, "S_small": small,
            "R_mangoldt": r, "R_prime_first_above_D": r_above,
            "R_prime_first_below_D": r_below, "R_first_proper_prime_power": r_powers,
            "T_bad": bad, "T_primitive": t, "T_over_NS": t / (n * ss),
            "T_m_composite": t + r_above,
            "corrected_R3_difference": (t + r_above) - (small - r_below - r_powers - bad)}


def direct_tail(n, d0, mu, lam, q=None):
    ans = 0.0
    for c in range(1, (n - 1) // (d0 + 1) + 1):
        if math.gcd(c, n) != 1:
            continue
        for d in range(d0 + 1, (n - 1) // c + 1):
            value = lam[n - c * d] if q is None else sum(
                mu_int(v) / phi(v) * sum(h * mu_int(v // h) for h in divisors(v)
                                         if (n - c * d) % h == 0)
                for v in range(1, q + 1))
            ans += int(mu[d]) * math.log(d) * value
    return ans


def model_tail_component(n, d0, h, pref):
    cs = np.arange(1, (n - 1) // (d0 + 1) + 1, dtype=np.int64)
    cs = cs[np.gcd(cs, n) == 1]
    cs = cs[np.gcd(cs, h) == 1]
    inverse = np.array([pow(r, -1, h) if math.gcd(r, h) == 1 else 0
                        for r in range(h)], dtype=np.int64)
    rs = (n % h) * inverse[cs % h] % h
    limits = (n - 1) // cs
    hi = limits - ((limits - rs) % h)
    lo = d0 - ((d0 - rs) % h)
    vals = np.zeros(hi.size)
    valid_hi = hi >= 0
    vals[valid_hi] = pref[hi[valid_hi]]
    valid = lo >= 0
    vals[valid] -= pref[lo[valid]]
    return float(vals.sum())


def selftest():
    mu, lam, _ = mobius_y_lambda(300)
    logs = np.zeros(301); logs[1:] = np.log(np.arange(1, 301))
    beta = mu * logs
    for n in range(4, 201, 2):
        row = tail_identity(n, mu, lam, beta)
        assert abs(row['T_primitive'] - direct_tail(n, row['D'], mu, lam)) < 1e-8
        assert abs(row['corrected_R3_difference']) < 1e-8
    coeff = coefficients(30, [1, 6, 30])
    for q, aa in coeff.items():
        pq = np.zeros(301)
        tail_model = {n: 0.0 for n in (32, 64, 128)}
        for h in np.flatnonzero(aa):
            pref = residue_prefix(mu.astype(np.int32), int(h))
            pq += aa[h] * (pref - mu)
            wp = residue_prefix(beta, int(h))
            for n in tail_model:
                tail_model[n] += aa[h] * model_tail_component(n, cube_root(n), int(h), wp)
        for n in range(4, 100):
            expected = math.fsum(int(mu[d]) * sum(mu_int(v) / phi(v) * sum(
                h * mu_int(v // h) for h in divisors(v) if (n - d) % h == 0)
                for v in range(1, q + 1)) for d in range(1, n))
            assert abs(pq[n] - expected) < 1e-8
        for n, value in tail_model.items():
            assert abs(value - direct_tail(n, cube_root(n), mu, lam, q)) < 1e-8
    counter = tail_identity(32, mu, lam, beta)
    assert abs(counter['T_bad']) > .1
    assert abs(counter['T_bad'] - math.log(2) * math.log(36)) < 1e-10
    print(json.dumps({"selftest": True, "R3_counterexample_N32": counter,
        "R3_missing_term_exact": "Tbad=log(2)*log(36), N=32,D=3"}))


def experiment(x, qs, sample_count):
    started = time.perf_counter()
    mu, lam, is_prime = mobius_y_lambda(x)
    prime_indices = np.flatnonzero(is_prime)
    del is_prime
    ns = np.arange(x // 2 + 1, x + 1, dtype=np.int64)
    ns = ns[ns % 2 == 0]
    qs = sorted(set(qs))
    print(json.dumps({"stage": "sieve_done", "X": x}), flush=True)
    m0, v0, fourier_major, grid = mangoldt_convolution(mu, lam, x, max(qs))
    rng = np.random.default_rng(20261123)
    validation_ns = rng.choice(ns, size=min(16, ns.size), replace=False)
    fft_error = max(abs(float(np.dot(mu[1:n], lam[n - 1:0:-1])) - m0[n])
                    for n in validation_ns)
    picked = np.sort(rng.choice(ns, size=min(6000, ns.size), replace=False))
    singular = np.array([serie_singular(int(n)) for n in picked])
    baseline = m0[ns] / np.sqrt(v0[ns])
    result = {"X": x, "window": "X/2 < N <= X, even N", "fft_max_absolute_error": fft_error,
              "raw_standardized_M0": stats(baseline), "raw_classes": classes(baseline, ns),
              "raw_corr_singular_uniform_seeded_sample": float(np.corrcoef(m0[picked] / np.sqrt(v0[picked]), singular)[0, 1]),
              "correlation_sample": {"n": int(picked.size), "seed": 20261123,
                  "method": "uniform without replacement; avoids stride aliasing modulo 5"},
              "discrete_fourier_projection": grid, "ramanujan": []}
    residual_fourier = (m0[ns] - fourier_major[ns]) / np.sqrt(v0[ns])
    result['discrete_fourier_residual'] = stats(residual_fourier)
    result['discrete_fourier_residual_classes'] = classes(residual_fourier, ns)
    del fourier_major, residual_fourier
    logs = np.zeros(x + 1); logs[1:] = np.log(np.arange(1, x + 1))
    beta = mu * logs
    tails_ns = sorted(set(map(int, rng.choice(ns, size=min(sample_count, ns.size), replace=False))) | {x - x % 2})
    rows = []
    for i, n in enumerate(tails_ns):
        rows.append(tail_identity(n, mu, lam, beta, prime_indices))
        if (i + 1) % 8 == 0:
            print(json.dumps({"stage": "full_tail_rows", "X": x, "done": i + 1, "total": len(tails_ns)}), flush=True)
    coeff = coefficients(max(qs), qs)
    all_h = sorted(set(int(h) for a in coeff.values() for h in np.flatnonzero(a)))
    predictions = {q: np.zeros(x + 1) for q in qs}
    tail_predictions = {q: np.zeros(len(rows)) for q in qs}
    for h in all_h:
        pref = residue_prefix(mu.astype(np.int32), h)
        mh = pref - mu
        for q in qs:
            if coeff[q][h]:
                predictions[q] += coeff[q][h] * mh
        del pref, mh
        wp = residue_prefix(beta, h)
        components = np.array([model_tail_component(row['N'], row['D'], h, wp) for row in rows])
        for q in qs:
            if coeff[q][h]:
                tail_predictions[q] += coeff[q][h] * components
        del wp
    y = m0[ns]
    for q in qs:
        p = predictions[q][ns]
        residual = y - p
        z = residual / np.sqrt(v0[ns])
        c = float(np.corrcoef(y, p)[0, 1])
        res = {"Q": q, "standardized_residual": stats(z), "classes": classes(z, ns),
               "corr_M0_PQ": c, "squared_correlation": c * c,
               "variance_fraction_removed": 1 - float(residual.var() / y.var()),
               "var_residual_over_mean_V0": float(residual.var() / v0[ns].mean()),
               "corr_residual_singular_uniform_seeded_sample": float(np.corrcoef(
                   (m0[picked] - predictions[q][picked]) / np.sqrt(v0[picked]), singular)[0, 1]),
               "interpretation": "finite Ramanujan-model residual; not identified with continuous minor arcs"}
        result['ramanujan'].append(res)
        for i, row in enumerate(rows):
            row.setdefault('ramanujan', {})[str(q)] = {
                "model_T_primitive": float(tail_predictions[q][i]),
                "model_T_over_NS": float(tail_predictions[q][i] / (row['N'] * row['singular_series'])),
                "residual_T_over_NS": float((row['T_primitive'] - tail_predictions[q][i]) / (row['N'] * row['singular_series']))}
    result['full_tail_sample'] = rows
    raw_tail = np.array([r['T_over_NS'] for r in rows])
    result['full_tail_summary'] = {"raw_T_over_NS": stats(raw_tail),
        "ramanujan_residual_T_over_NS": {str(q): stats(np.array([
            r['ramanujan'][str(q)]['residual_T_over_NS'] for r in rows])) for q in qs},
        "number_above_one": int(np.sum(raw_tail >= 1)),
        "method": "S-R-Tbad identity; finite diagnostic dependent on computed R, not an independent HB+ proof",
        "scope": "seeded sample plus N=X, not all even N in window"}
    result['elapsed_seconds'] = time.perf_counter() - started
    print(json.dumps({"stage": "completed", "X": x, "elapsed_seconds": result['elapsed_seconds'],
                      "Q30_residual_stats": result['ramanujan'][-1]['standardized_residual'],
                      "max_full_T_over_NS": float(raw_tail.max())}), flush=True)
    return result


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--hasta', type=float, nargs='+', default=[1e5, 1e6, 1e7])
    parser.add_argument('--Q', type=int, nargs='+', default=[1, 6, 30])
    parser.add_argument('--muestras-cola', type=int, default=24)
    parser.add_argument('--selftest', action='store_true')
    args = parser.parse_args()
    if args.selftest:
        selftest()
        return
    if max(args.hasta) > 1e7 or max(args.Q) > 60:
        parser.error('Resource guard: X<=1e7, Q<=60')
    results = []
    for xf in args.hasta:
        results.append(experiment(int(xf), args.Q, args.muestras_cola))
        gc.collect()
    base = Path(__file__).resolve().parent
    stamp = datetime.now(timezone.utc).strftime('%Y%m%d_%H%M%S_%f')
    file = base / 'resultados' / f'{stamp}_codex_residuos_06_07.json'
    payload = {"provenance": {"agent": "Codex", "model": "identifier not exposed",
        "project_date": "2026-11-23", "run_utc": datetime.now(timezone.utc).isoformat(),
        "status": "[N] not a proof; non-independent review", "python": platform.python_version(),
        "numpy": np.__version__, "script_sha256": hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
        "source_script_hashes": {name: hashlib.sha256((base/name).read_bytes()).hexdigest()
                                 for name in ('cajas_bilineales.py', 'goldbach.py')},
        "sources": ['matematicas/04', 'matematicas/06', 'matematicas/07'],
        "no_models_or_services_changed": True}, "results": results}
    with file.open('x', encoding='utf-8') as f:
        json.dump(payload, f, indent=2, ensure_ascii=False, allow_nan=False)
        f.write('\n')
    print(json.dumps({"result_file": str(file)}), flush=True)


if __name__ == '__main__':
    main()
