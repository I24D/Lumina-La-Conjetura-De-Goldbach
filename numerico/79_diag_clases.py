#!/usr/bin/env python3
"""Diagnostico doc 79: sumas de clase de nu_U y de S(m) para d=7,11.
Pregunta: equidistribuye nu_U sobre clases reducidas (1/phi(d) del total)
o sobre todas (1/d)? Y S(m) sobre la clase N mod d?
"""
import numpy as np

for X in (10 ** 5, 10 ** 6):
    # criba rapida de mu hasta X
    spf = np.zeros(X + 1, dtype=np.int64)
    for i in range(2, X + 1):
        if spf[i] == 0:
            spf[i] = i
            if i * i <= X:
                spf[i * i:X + 1:i] = np.where(spf[i * i:X + 1:i] == 0, i,
                                              spf[i * i:X + 1:i])
    primes = np.nonzero(spf == np.arange(X + 1))[0]
    primes = primes[primes >= 2]
    mu = np.ones(X + 1, dtype=np.int64)
    mu[0] = 0
    for p in primes:
        mu[p::p] *= -1
        if p * p <= X:
            mu[p * p::p * p] = 0
    lam = np.zeros(X + 1)
    for i in range(2, X + 1):
        p = spf[i]
        n = i
        while n % p == 0:
            n //= p
        if n == 1:
            lam[i] = np.log(p)
    U = 30
    nu = np.zeros(X + 1)
    nu[1] = 1.0
    for t in range(1, U + 1):
        if mu[t] != 0:
            nu[t::t] -= mu[t]
    # S(m) con V=30 y V=0
    def S_of(V):
        S = np.zeros(X + 1)
        for a in np.nonzero(lam > 0)[0]:
            if a <= V or a > X:
                continue
            bmax = X // a
            if bmax <= U:
                continue
            S[a * (U + 1):X + 1:a] += lam[a] * nu[U + 1:bmax + 1]
        return S
    S30 = S_of(30)
    Sfull = S_of(0)
    print(f"=== X={X} ===")
    print(f"media nu_U={nu[1:].mean():.4f}  media S30={S30[1:].mean():.4f}  "
          f"media Sfull={Sfull[1:].mean():.4f}  media Lam={lam[1:].mean():.4f}")
    idx = np.arange(X + 1)
    for d in (7, 11):
        tot_nu = nu[1:].sum()
        print(f"-- d={d} (phi={d-1}) total nu={tot_nu:.0f}")
        for s in range(d):
            cs = nu[(idx % d) == s].sum()
            print(f"   clase nu s={s}: {cs:.0f}  (T/phi={tot_nu/(d-1):.0f}, T/d={tot_nu/d:.0f})")
        N0 = 999983
        r = N0 % d
        for name, arr in (("S30", S30), ("Sfull", Sfull), ("Lam", lam)):
            tot = arr[1:].sum()
            cl = arr[(idx % d) == r].sum()
            ratio = ((d / X) * cl) / (tot / X) if tot != 0 else float("nan")
            print(f"   {name}: clase={cl:.0f} total={tot:.0f} ratio={ratio:.3f} "
                  f"(W=1 -> 1.000; W=phi/d -> {(d-1)/d:.3f}; Lam d/phi -> {d/(d-1):.3f})")
