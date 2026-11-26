#!/usr/bin/env python3
"""Doc 79 - Peso movil W(d) = phi(d)/d en el termino principal M_B3 de B3.
Verifica: (1) media de S(m)=(Lambda_{>V}*nu_U)(m) y su formula exacta;
(2) rango y media de W(d)=phi(d)/d y caida de sum_{d~K} mu(d) log d W(d);
(3) densidad de clase de S (ratio = W(d)) y de Lambda (ratio = 1/W(d));
(4) tamano del termino principal M_B3 en una caja a N=1e6.
Ilustracion a escala finita: la pequeez asintotica es PNT (citada en el doc).
Uso: python3 79_peso_movil_Wd.py  -> PASS n/n o FAIL con detalle.
"""
import math
import sys
import numpy as np

X = 10 ** 6
N = 10 ** 6
L = np.log(N)
ok = []
fails = []


def check(name, cond, info=""):
    ok.append(bool(cond))
    print(("PASS " if cond else "FAIL ") + name + (" | " + str(info) if info else ""))
    if not cond:
        fails.append(name)


# ---------- criba: spf, mu, Lambda, phi hasta X ----------
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
phi = np.arange(X + 1, dtype=float)
for p in primes:
    phi[p::p] *= (p - 1.0) / p
phi[0] = 0.0
check("mu(6)=1, mu(30)=-1, mu(12)=0, phi(6)=2, phi(30)=8",
      mu[6] == 1 and mu[30] == -1 and mu[12] == 0 and phi[6] == 2 and phi[30] == 8)
check("sum Lambda(1e6) ~ 1e6", abs(lam[1:].sum() / X - 1) < 0.01,
      f"{lam[1:].sum() / X:.4f}")

# ---------- nu_U(b) = sum_{t|b, t>U} mu(t) ----------
def nu_U(U):
    nu = np.zeros(X + 1)
    nu[1] = 1.0
    for t in range(1, U + 1):
        if mu[t] != 0:
            nu[t::t] -= mu[t]
    return nu


def S_of(U, V):
    """S(m) = sum_{ab=m, a>V, b>U} Lambda(a) nu_U(b). V=0 -> convolucion completa."""
    nu = nu_U(U)
    S = np.zeros(X + 1)
    for a in np.nonzero(lam > 0)[0]:
        if a <= V or a > X:
            continue
        bmax = X // a
        if bmax <= U:
            continue
        S[a * (U + 1):X + 1:a] += lam[a] * nu[U + 1:bmax + 1]
    return S, nu


# ---------- (1) media de S: identidad exacta (Lambda*nu_U) = log * mu_{>U} ----------
logcum = np.concatenate([[0.0], np.cumsum(np.log(np.arange(1, X + 1, dtype=float)))])
# G(y) = sum_{n<=y} log n = log(floor(y)!)
for U in (10, 30):
    Sfull, _ = S_of(U, 0)
    mean_full = Sfull[1:].sum() / X
    t = np.arange(U + 1, X + 1)
    exact = (mu[t] * logcum[X // t]).sum() / X
    s1 = (mu[t] / t).sum()
    s2 = (mu[t] * np.log(t) / t).sum()
    T = np.log(X) * s1 - s2
    print(f"U={U}: media conv completa={mean_full:.4f} exacta={exact:.4f} "
          f"asintotica T={T:.4f} (s1={s1:.4f} s2={s2:.4f})")
    check(f"media conv = formula exacta, U={U}", abs(mean_full - exact) < 1e-6,
          f"dif={mean_full - exact:.2e}")
    S, _ = S_of(U, U)
    meanS = S[1:].sum() / X
    print(f"U={U}: media(S, a>V)={meanS:.4f} (parte a<=V: {mean_full - meanS:.4f})")
    check(f"|media(S)| < 1, U={U}", abs(meanS) < 1.0, f"{meanS:.4f}")

# nu(U) = media de nu_U: decrece con U (PNT para mu)
print("nu(U) = media de nu_U hasta 1e6:")
nus = {}
for U in (10, 30, 100, 300, 1000):
    nus[U] = nu_U(U)[1:].mean()
    print(f"  U={U}: nu={nus[U]:.4f}")
check("|nu(1000)| < |nu(10)|", abs(nus[1000]) < abs(nus[10]),
      f"{nus[1000]:.4f} vs {nus[10]:.4f}")

# ---------- (2) W(d) = phi(d)/d ----------
W = np.zeros(X + 1)
W[1:] = phi[1:] / np.arange(1, X + 1, dtype=float)
for K in (1000, 10000, 100000):
    d = np.arange(K + 1, 2 * K + 1)
    w = W[d]
    s_unw = abs((mu[d] * np.log(d)).sum()) / K
    s_w = abs((mu[d] * np.log(d) * w).sum()) / K
    print(f"K={K}: W min={w.min():.4f} max={w.max():.4f} media={w.mean():.4f} "
          f"|sum mu log d|/K={s_unw:.4f} |sum mu log d W|/K={s_w:.4f}")
    check(f"W en (0,1] K={K}", w.min() > 0 and w.max() <= 1.0, f"min={w.min():.4f}")
    check(f"media W -> 6/pi^2, K={K}", abs(w.mean() - 6 / np.pi ** 2) < 0.03,
          f"{w.mean():.4f} vs 0.6079")
    check(f"|sum mu log d W|/K < 0.06, K={K}", s_w < 0.06, f"{s_w:.4f}")
check("W(7)=6/7, W(30)=8/30, W(210)=48/210",
      abs(W[7] - 6 / 7) < 1e-12 and abs(W[30] - 8 / 30) < 1e-12
      and abs(W[210] - 48 / 210) < 1e-12)

# ---------- (3) densidad de clase: W(d) = cola de Mertens restringida / cola total ----------
# S_nu(U,d) = sum_{t>U,(t,d)=1} mu(t)/t = -sum_{t<=U,(t,d)=1} mu(t)/t  (total -> 0)
def Snu(U, d):
    t = np.arange(1, U + 1)
    m = np.array([math.gcd(int(x), d) == 1 for x in t])
    return -float((mu[t][m] / t[m]).sum())


S30, _ = S_of(30, 30)
meanS30 = S30[1:].sum() / X
meanLam = lam[1:].sum() / X
N0 = 999983  # primo, coprimo con 7, 11, 30
m_idx = np.arange(X + 1)
Snu1 = Snu(30, 1)
for d in (7, 11, 30):
    mask = (m_idx % d) == (N0 % d)
    r_lam = ((d / X) * lam[mask].sum()) / meanLam
    r_S = ((d / X) * S30[mask].sum()) / meanS30
    Wth = Snu(30, d) / Snu1
    print(f"d={d}: ratio Lambda={r_lam:.3f} (esperado d/phi={d / phi[d]:.3f}); "
          f"ratio S={r_S:.3f} (teoria W={Wth:.3f}; phi(d)/d={W[d]:.3f})")
    check(f"ratio Lambda = d/phi(d), d={d}", abs(r_lam - d / phi[d]) < 0.05,
          f"{r_lam:.3f} vs {d / phi[d]:.3f}")
    check(f"ratio S = W(d) (cola Mertens), d={d}",
          np.sign(r_S) == np.sign(Wth) and abs(r_S - Wth) < 0.5,
          f"{r_S:.3f} vs {Wth:.3f}")
# concentracion de nu_U en b = 0 mod d (no afecta la clase N mod d, (N,d)=1)
nu30 = nu_U(30)
for d in (7, 11):
    tot = nu30[1:].sum()
    s0 = nu30[(m_idx % d) == 0].sum()
    print(f"nu_U clase b=0 mod {d}: {s0:.0f} de total {tot:.0f} ({s0 / tot:.2f})")
check("nu_U concentrado en b=0 mod 7", nu30[(m_idx % 7) == 0].sum() > 0.5 * nu30[1:].sum())

# ---------- (4) M_B3 en una caja a N=1e6 ----------
# M_B3 = C (phi(N)/N) L sum_{d~K} mu(d) log d S_nu(U,d)  (la S_nu(U,1) se cancela)
U = V = 30
S, nu = S_of(U, V)
C, K = 512, 1024  # CK = 524288 < N, K > sqrt(N)/L
B3 = 0.0
for d in range(K + 1, 2 * K + 1):
    if mu[d] == 0:
        continue
    cvals = np.arange(C + 1, 2 * C + 1)
    mvals = N - cvals * d
    okm = mvals > 0
    B3 += mu[d] * np.log(d) * S[mvals[okm]].sum()
dbox = np.arange(K + 1, 2 * K + 1)
Snu_d = np.array([Snu(U, int(dd)) for dd in dbox])
M_B3 = C * (phi[N] / N) * L * (mu[dbox] * np.log(dbox) * Snu_d).sum()
M_naive = C * (phi[N] / N) * nus[30] * L * (mu[dbox] * np.log(dbox) * W[dbox]).sum()
print(f"caja C={C} K={K}: B3/N={B3 / N:.4f}  M_B3/N={M_B3 / N:.4f}  "
      f"M_ingenuo(phi(d)/d)/N={M_naive / N:.4f}  |M|/|B3|="
      f"{abs(M_B3 / B3) if B3 != 0 else float('nan'):.3f}")
check("|M_B3| < 0.05 N", abs(M_B3) < 0.05 * N, f"{M_B3 / N:.4f}")
check("B3 caja es o(N) (|B3| < 0.1 N)", abs(B3) < 0.1 * N, f"{B3 / N:.4f}")

print(f"\n{sum(ok)}/{len(ok)} PASS")
sys.exit(0 if all(ok) else 1)
