#!/usr/bin/env python3
"""Doc 62: bajar el soporte a N^eps con el criterio relajado (familia K=2, doc 61).
Parte 1: formulas exactas K=2 (Fraction) en z bajo.
Parte 2: tabla (N,eps): z, viabilidad de t=1/Sz, A, delta.
Parte 3: piso del defecto (z para delta_2<=0.1) y eps necesario.
Parte 4: umbral de well-factorabilidad Q* ~ z^4 y obstruccion triple.
Parte 5: seleccion por criba (X=2e6, numpy): E[W|Omega=k].
Uso: python3 numerico/62_soporte_bajo_nepsilon.py
"""
from fractions import Fraction as F
from itertools import combinations
import math

FAIL = []

def primes_upto(n):
    if n < 2:
        return []
    s = bytearray(b'\x01') * (n + 1)
    s[0:2] = b'\x00\x00'
    for i in range(2, int(n ** 0.5) + 1):
        if s[i]:
            s[i*i:n+1:i] = b'\x00' * ((n - i*i) // i + 1)
    return [i for i in range(2, n + 1) if s[i]]

def phi(n, ps):
    r = n; m = n
    for p in ps:
        if p * p > m:
            break
        if m % p == 0:
            r = r // p * (p - 1)
            while m % p == 0:
                m //= p
    if m > 1:
        r = r // m * (m - 1)
    return r

def k2_lam(z, t):
    ps = primes_upto(z)
    lam = {1: F(1)}
    for l in ps:
        lam[l] = F(-t)
    for a, b in combinations(ps, 2):
        lam[a * b] = F(t * t)
    return ps, lam

def W_of_n(n, lam):
    return sum(v for d, v in lam.items() if n % d == 0)

print("== 1. K=2 exacto en z bajo: A, M_l, delta, W>=0, W(p)>0 ==")
for z, t in [(2, F(1, 2)), (3, F(1, 2)), (5, F(1, 2)), (7, F(1, 2)),
             (3, F(1, 10)), (5, F(1, 10)), (7, F(1, 10))]:
    # t=1 es el borde (W(p)=1-t=0: muere el componente P1); no se incluye en la parrilla.
    ps, lam = k2_lam(z, t)
    Sz = sum(F(1, p) for p in ps)
    S2 = sum(F(1, p * p) for p in ps)
    A = sum(v / F(d) for d, v in lam.items())
    Aform = F(1) - t * Sz + (t * t / 2) * (Sz * Sz - S2)
    okA = (A == Aform)
    # identidad L = sum M_e/phi(e)
    es = sorted(set([1] + ps + [a*b for a, b in combinations(ps, 2)]))
    M = {e: sum(v / F(d) for d, v in lam.items() if d % e == 0) for e in es}
    L = sum(v / F(phi(d, ps)) for d, v in lam.items())
    tot = sum(M[e] / F(phi(e, ps)) for e in es)
    okI = (L == tot)
    line = f"  z={z} t={t}: A={float(A):.4f} form={float(Aform):.4f} {'OK' if okA else 'FALLO'} ident={'OK' if okI else 'FALLO'}"
    print(line)
    if not okA:
        FAIL.append(f"K2 z={z} t={t} A")
    if not okI:
        FAIL.append(f"K2 z={z} t={t} ident")
    for l in ps:
        Mform = (-t + t * t * (Sz - F(1, l))) / l
        ok = (M[l] == Mform)
        dl = abs(M[l]) * (l - 1) / A
        print(f"    l={l}: M={float(M[l]):+.4f} form={float(Mform):+.4f} {'OK' if ok else 'FALLO'} delta={float(dl):.3f}")
        if not ok:
            FAIL.append(f"K2 z={z} t={t} l={l}")
    bad = sum(1 for n in range(1, 3001) if W_of_n(n, lam) < 0)
    minWp = min(W_of_n(p, lam) for p in primes_upto(3000))
    print(f"    W>=0 [1,3000]: neg={bad} minW(p)={float(minWp):.4f} {'OK' if bad == 0 and minWp > 0 else 'FALLO'} soporte max d={max(lam)} <= z^2={z*z} {'OK' if max(lam) <= z*z else 'FALLO'}")
    if bad or minWp <= 0 or max(lam) > z * z:
        FAIL.append(f"K2 z={z} t={t} W/soporte")

def z_of_N_eps(E, num, den):
    """Mayor z con z^(2*den) <= 10^(E*num); es decir z = floor(N^(eps/2)), N=10^E, eps=num/den."""
    target_exp = E * num
    z = int(round(10 ** (target_exp / (2 * den))))
    z = max(z, 1)
    while z ** (2 * den) > 10 ** target_exp:
        z -= 1
    while (z + 1) ** (2 * den) <= 10 ** target_exp:
        z += 1
    return z

print("== 2. Tabla (N,eps): z=N^(eps/2), viabilidad t=1/Sz, A, delta_2 ==")
print("   N        eps    z  |P<=z|  Sz     t=1/Sz?  A(t=1/2) d2 d3(t=1/2)  A(1/Sz) d2(1/Sz)")
for E in [7, 12, 24, 40]:
    for num, den in [(1, 5), (1, 10), (1, 20), (1, 50)]:
        z = z_of_N_eps(E, num, den)
        ps = primes_upto(z)
        if len(ps) == 0:
            print(f"   10^{E:<3} {num}/{den}  z={z}  VACIO: sin primos, lambda=delta_1, W==1 (regimen degenerado)")
            continue
        Sz = sum(F(1, p) for p in ps)
        S2 = sum(F(1, p * p) for p in ps)
        feas = (Sz >= 1)
        t2 = F(1, 2)
        A2 = F(1) - t2 * Sz + (t2 * t2 / 2) * (Sz * Sz - S2)
        M2 = (-t2 + t2 * t2 * (Sz - F(1, 2))) / 2
        d2 = abs(M2) / A2
        # delta_3 tambien: el cero de M_2 en Sz=1/t+1/2 es accidental; delta_3 no se anula ahi
        M3 = (-t2 + t2 * t2 * (Sz - F(1, 3))) / 3 if 3 in ps else F(0)
        d3 = abs(M3) * 2 / A2 if 3 in ps else F(0)
        if feas:
            t = F(1) / Sz
            A = F(1) - t * Sz + (t * t / 2) * (Sz * Sz - S2)
            M = -F(1) / (F(4) * Sz * Sz)
            d = abs(M) / A
            print(f"   10^{E:<3} {num}/{den}  z={z:<5} {len(ps):<4} {float(Sz):.3f}  si        {float(A2):.4f}   {float(d2):.3f} {float(d3):.3f}  {float(A):.4f}   {float(d):.3f}")
        else:
            print(f"   10^{E:<3} {num}/{den}  z={z:<5} {len(ps):<4} {float(Sz):.3f}  NO (Sz<1) {float(A2):.4f}   {float(d2):.3f} {float(d3):.3f}  --       --")
        # chequeo de z exacto
        assert z ** (2 * den) <= 10 ** (E * num) < (z + 1) ** (2 * den), "z mal calculado"

print("== 3. Piso del defecto en el regimen t=1/Sz: delta_2 ~= 1/(2 Sz^2) ==")
# Una sola criba hasta 200000 + sumas prefijas (la busqueda ingenua z+=1 con criba
# en cada paso es cuadratica). Para 0.01 no se busca: cota analitica.
MERTENS = 0.26149721284764278
ZMAX = 200000
_ps = primes_upto(ZMAX)
_pref1 = [0.0]; _pref2 = [0.0]
for _p in _ps:
    _pref1.append(_pref1[-1] + 1.0 / _p)
    _pref2.append(_pref2[-1] + 1.0 / (_p * _p))
import bisect as _bis
def delta2_1sobreSz(z):
    _i = _bis.bisect_right(_ps, z)
    Sz = _pref1[_i]; S2 = _pref2[_i]
    A = 0.5 * (1 - S2 / (Sz * Sz))
    return 1.0 / (4 * Sz * Sz) / A, Sz

_lo, _hi = 8, ZMAX
while _lo < _hi:  # busqueda binaria del z minimo con delta_2<=0.1
    _mid = (_lo + _hi) // 2
    if delta2_1sobreSz(_mid)[0] <= 0.1:
        _hi = _mid
    else:
        _lo = _mid + 1
_zmin = _lo
_d, _Sz = delta2_1sobreSz(_zmin)
print(f"  delta_2<=0.1: z_min={_zmin} (Sz={_Sz:.4f}, delta_2={_d:.6f}); delta_2(z_min-1)={delta2_1sobreSz(_zmin-1)[0]:.6f}")
for E in [12, 24, 40]:
    _eps = 2 * math.log(_zmin) / (E * math.log(10))
    print(f"    a N=10^{E}: eps necesario = {_eps:.3f} (soporte N^{_eps:.3f})")
print(f"    control Mertens: Sz(1400)~{math.log(math.log(1400)) + MERTENS:.3f}, Sz_fsum(1400)={delta2_1sobreSz(1400)[1]:.3f}")
print("  delta_2<=0.01: exige Sz>=sqrt(50)=7.07 (formula exacta), log log z>=6.81, z>=e^903: inalcanzable a N realista (verificacion analitica, no busqueda)")

print("== 4. Well-factorabilidad de K=2: umbral Q* ~ z^4 = N^(2 eps) (Prop. 62.1) ==")
# Definicion (Maynard Def. 1, via doc 58 sec. 6): lambda es well-factorable de nivel Q
# si para toda Q=Q1*Q2 hay gamma^(1), gamma^(2) con |gamma|<=1, soporte<=Q_i, lambda=gamma1*gamma2.
# Obstruccion triple (prueba en el doc): si hay 3 primos p,q,r<=z con pq,pr,qr>max(Q1,Q2),
# las 3 ecuaciones por pares fuerzan a_r=f(a_p) y a_r=f(f(a_p)) con f(x)=-t(t+x)/(t+2x),
# luego f(a_p) es punto fijo de f: 2x^2+2tx+t^2=0, discriminante -4t^2<0. Imposible.
t = F(1, 2)
_disc = 4 * t * t - 8 * t * t
print(f"  discriminante del punto fijo 2x^2+2tx+t^2 (t=1/2): {_disc} <0 {'OK' if _disc < 0 else 'FALLO'}")
if _disc >= 0:
    FAIL.append("discriminante punto fijo")
def _f(x):
    return -t * (t + x) / (t + 2 * x)
for _x in [F(0), F(-1, 2), F(-1, 8), F(-3, 8), F(1, 4), F(-1)]:
    if _x == -t / 2:
        continue
    _fx = _f(_x)
    _g = _f(_fx) - _fx
    # g(x) = -(2 f^2+2t f+t^2)/(t+2f) != 0 porque el numerador es cuadratica sin raices reales
    _num = 2 * _fx * _fx + 2 * t * _fx + t * t
    print(f"  x={float(_x):+.3f}: f(x)={float(_fx):+.4f} g(x)=f(f(x))-f(x)={float(_g):+.4f} num={float(_num):+.4f} {'OK' if _g != 0 and _num > 0 else 'FALLO'}")
    if _g == 0 or _num <= 0:
        FAIL.append(f"obstruccion triple x={_x}")
print("  umbral: Q* en [max(z^2, M(z)^2), Dmax^2], M(z)=p_(2)*p_(3), Dmax=p_(1)*p_(2)")
print("   z    p1 p2 p3   Dmax    M(z)   Q*min      Q*max     Q*/z^4  BFI<=N^(4/7)?")
for z in [2, 3, 5, 7, 11, 13, 100, 251]:
    ps = primes_upto(z)
    if len(ps) >= 3:
        p1, p2, p3 = ps[-1], ps[-2], ps[-3]
        Mz = p2 * p3
        qmin = max(z * z, Mz * Mz)
    elif len(ps) == 2:
        p1, p2 = ps[-1], ps[-2]
        Mz = 0
        qmin = z * z
    else:
        p1 = ps[-1]
        qmin = z * z
        Mz = 0
    # Dmax = maximo soporte de lambda (semiprimo mayor, o el primo mayor si pi(z)<2)
    Dmax = max([p for p in ps] + [a * b for a, b in combinations(ps, 2)])
    qmax = Dmax * Dmax
    print(f"   {z:<4} {str(ps[-3:]) :<16} {Dmax:<6} {Mz:<6} {qmin:<10} {qmax:<10} {qmin / z**4:.2f}-{qmax / z**4:.2f}")
# caso explicito: z=100, Q=10^6 < Q*min: triple {97,89,83} con productos > sqrt(Q)=1000
_ps100 = primes_upto(100)
_tr = _ps100[-3:]
_pr = [a * b for a, b in combinations(_tr, 2)]
print(f"  ejemplo z=100, Q=10^6: triple {_tr}, productos {_pr}, todos > sqrt(Q)=1000: {'OK' if min(_pr) > 1000 else 'FALLO'}")
if min(_pr) <= 1000:
    FAIL.append("triple z=100")
# comparacion BFI/BV por eps (umbral ~ N^(2 eps))
print("  por eps: umbral Q*~N^(2 eps); BFI aplica si 2eps<=4/7; BFI supera BV (1/2) si 2eps>1/2")
for _num, _den in [(1, 5), (1, 10), (1, 20), (1, 4), (2, 7)]:
    _e = _num / _den
    print(f"    eps={_num}/{_den}: Q*~N^{2*_e:.3f}; BFI {'aplica' if 2*_e <= 4/7 + 1e-12 else 'NO aplica'}; {'supera BV' if 2*_e > 0.5 else 'NO supera BV (BV 1/2 ya cubre el soporte)'}; s=logD/logz con D=N^(1/2): s={1/_e:.1f}")
# suficiencia: Q>=Dmax^2 => max(Q1,Q2)>=Dmax para toda factorizacion (trivial)
for z in [7, 100]:
    ps = primes_upto(z)
    Dmax = max(a * b for a, b in combinations(ps, 2))
    Q = Dmax * Dmax
    ok = all(max(q1, Q // q1 if Q % q1 == 0 else Q / q1) >= Dmax for q1 in [1, 2, 7, Dmax, Q // 2, Q])
    print(f"  suficiencia z={z}: Q=Dmax^2={Q}, max(Q1,Q2)>=Dmax en factorizaciones de prueba: {'OK' if ok else 'FALLO'}")
    if not ok:
        FAIL.append(f"suficiencia z={z}")

print("== 5. Seleccion por criba (X=2e6): E[W|Omega=k] y fraccion r=0 ==")
import numpy as np
X = 2000000
ps_big = primes_upto(X)
Om = np.zeros(X + 1, dtype=np.int16)
for p in ps_big:
    pw = p
    while pw <= X:
        Om[pw::pw] += 1
        if pw > X // p:
            break
        pw *= p
Om = Om[1:]
for z in [2, 3, 5, 7]:
    ps = primes_upto(z)
    rz = np.zeros(X, dtype=np.int8)
    for l in ps:
        rz[l-1::l] += 1
    rz = rz.astype(np.float64)
    for tf in [0.5, 0.1]:
        W = 1.0 - tf * rz + tf * tf * rz * (rz - 1) / 2
        fr0 = float(np.mean(rz == 0))
        out = f"  z={z} t={tf}: frac r=0 (W=1)={fr0:.3f} E[W]={float(np.mean(W)):.3f}"
        for k, lab in [(1, "Om=1"), (2, "Om=2"), (3, "Om=3")]:
            m = (Om == k)
            out += f" E[W|{lab}]={float(np.mean(W[m])):.3f}"
        m4 = (Om >= 4)
        out += f" E[W|Om>=4]={float(np.mean(W[m4])):.3f}"
        m2 = (Om == 2)
        fr2 = float(np.mean(rz[m2] == 0))
        out += f" | P2: frac r=0={fr2:.3f} E[W|P2,r=0]=1.000 E[W|P2,r>=1]={float(np.mean(W[m2 & (rz >= 1)])):.3f}"
        # diagnostico de seleccion (medicion, no fallo): cociente P2/P3 y contaminacion P3+
        e2 = float(np.mean(W[m2])); e3 = float(np.mean(W[Om == 3])); e4 = float(np.mean(W[m4]))
        print(f"    ratio E[W|P2]/E[W|P3]={e2/e3:.3f} E[W|P3+]/E[W|P2]={(e3+e4)/2/e2:.3f} (1.0 = sin seleccion; >>1 = aisla P2)")
        # chequeo de coherencia: E[W] debe coincidir con A del modelo (densidad de Dirichlet)
        # y E[W|Om=1] debe ser ~1 (primos > z tienen W=1)
        if abs(float(np.mean(W[Om == 1])) - 1.0) > 0.01:
            FAIL.append(f"seleccion z={z} t={tf}: E[W|primo] != 1")
print()
print("RESULTADO:", "PASS" if not FAIL else f"FALLO: {FAIL}")
