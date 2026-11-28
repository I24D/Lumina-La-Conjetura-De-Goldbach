#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Verificación independiente (Claude Code, brazo verificador de la línea
Goldbach) de las afirmaciones [T] y [M] de los docs 116-120a
(UNIF-U / (H2)): 116a/b/c, 117a/b, 118a/b1/b2, 119a/b, 120a.

Reimplementación propia desde las definiciones ESCRITAS de los docs
109, 110, 113, 114 y de cada doc verificado. No importa código de
ningún otro agente. Del script del doc 109 solo se tomó, leyéndola, la
regla de muestreo de niveles (sample_sprimes), que no está escrita en
los docs y es necesaria para reproducir sus conteos exactos.

Definiciones (docs 109/110/113/114):
  P = primos impares < z que no dividen a A = 3 (luego p >= 5), N = 2^e.
  V+(L): d libre de cuadrados sobre P, d < L y, para m impar de la
    factorización decreciente, p_m^2 (p_1...p_m) < L.  Equivale a
    lam(d) < L con lam(d) = max(d, max_{m impar} p_m^2 P_m).
  V^(r)  = {d' : (d',r)=1, r d' en V+(L)}          (soporte desplazado)
  Vhat^(r) = V+(L/r; P \\ r) en aritmética entera  <=>  r lam(d') < L
  S(r), Shat(r): sumas de a(d') = mu(d')/phi(d') sobre esos soportes.
  Pi(r) = prod_{p en P, p no | r} (1 - 1/(p-1));  delta = S - Pi;
  E = S - Shat;  u_r = ln(L/r)/ln z;  Q_c = prod (1 - 1/(p-1) + 1/(p-1)^2).
  Banda: r visible, u_r < 3/2, Pi(r) < 1/2.
Marcos: A z = 17, 29, 37 (D = z^5, y = int(z^(10/3))) y B z = 29
(D = z^6, y = z^3); niveles L = D//s con s primo en [z, y), 60
muestreados por índice (regla del doc 109).

Solo stdlib.  Uso:  python verificacion_116_120_claude_code.py
"""

import math
import sys
import time
from bisect import bisect_left, bisect_right
from fractions import Fraction
from itertools import compress

A = 3
U0 = 1.5
T0 = time.time()
CHECKS = []


def check(cid, desc, ok, detail=""):
    CHECKS.append((cid, desc, bool(ok), detail))
    print("  [%s] %s: %s%s" % ("PASS" if ok else "FAIL", cid, desc,
                               (" | " + detail) if detail else ""))
    sys.stdout.flush()


def head(txt):
    print("\n" + "=" * 78 + "\n" + txt + "\n" + "=" * 78)
    sys.stdout.flush()


# ------------------------------------------------------------- aritmética
def primes_upto(n):
    if n < 2:
        return []
    bs = bytearray(b"\x01") * (n + 1)
    bs[0] = 0
    bs[1] = 0
    for i in range(2, math.isqrt(n) + 1):
        if bs[i]:
            bs[i * i::i] = bytes(len(range(i * i, n + 1, i)))
    return list(compress(range(n + 1), bs))


_SMALLP = primes_upto(20000)


def primes_open(lo, hi, seg=10 ** 7):
    """Lista (por segmentos) de los primos p con lo < p < hi (lo, hi reales)."""
    a = math.floor(lo) + 1
    b = math.ceil(hi) - 1
    if b < a:
        return
    base = primes_upto(math.isqrt(b) + 1)
    start = max(a, 2)
    while start <= b:
        end = min(b, start + seg - 1)
        n = end - start + 1
        bs = bytearray(b"\x01") * n
        for p in base:
            if p * p > end:
                break
            s0 = max(p * p, ((start + p - 1) // p) * p)
            if s0 > end:
                continue
            bs[s0 - start::p] = bytes(len(range(s0 - start, n, p)))
        yield list(compress(range(start, end + 1), bs))
        start = end + 1


def fsum_primes(lo, hi, f):
    tot = []
    for blk in primes_open(lo, hi):
        tot.append(math.fsum(f(p) for p in blk))
    return math.fsum(tot)


def is_prime(n):
    if n < 2:
        return False
    for p in (2, 3, 5, 7, 11, 13, 17, 19, 23, 29, 31, 37):
        if n % p == 0:
            return n == p
    d, s = n - 1, 0
    while d % 2 == 0:
        d //= 2
        s += 1
    for a in (2, 3, 5, 7, 11, 13, 17, 19, 23, 29, 31, 37):
        x = pow(a, d, n)
        if x in (1, n - 1):
            continue
        for _ in range(s - 1):
            x = x * x % n
            if x == n - 1:
                break
        else:
            return False
    return True


def lam_of(pf):
    """lam(d) para d con factorización DECRECIENTE pf."""
    prod = 1
    lam = 1
    for i, p in enumerate(pf):
        prod *= p
        if i % 2 == 0:
            t = p * p * prod
            if t > lam:
                lam = t
    return max(lam, prod)


def prod_of(pf):
    out = 1
    for p in pf:
        out *= p
    return out


def in_V(pf, L):
    return lam_of(pf) < L


def in_Vhat(pf, r, L):
    return r * lam_of(pf) < L


def merged(pf, rpf):
    return tuple(sorted(tuple(pf) + tuple(rpf), reverse=True))


def mertens(x):
    """Sum_{5 <= p <= x} 1/(p-1), fracción exacta."""
    return sum((Fraction(1, p - 1) for p in _SMALLP if 5 <= p <= x),
               Fraction(0))


def sample_sprimes(z, y, nsamp=60):
    """Regla de muestreo del doc 109 (leída en su script, no importada)."""
    spr = [p for p in primes_upto(y - 1) if p >= z and math.gcd(p, A) == 1]
    if len(spr) <= nsamp:
        return spr
    idx = sorted({round(i * (len(spr) - 1) / (nsamp - 1.0))
                  for i in range(nsamp)})
    return [spr[i] for i in idx]


def frame_levels(z, D, y, nsamp=60):
    return sorted({D // s for s in sample_sprimes(z, y, nsamp)})


FRAMES = [("A", 17, 17 ** 5, int(17 ** (10.0 / 3.0))),
          ("A", 29, 29 ** 5, int(29 ** (10.0 / 3.0))),
          ("A", 37, 37 ** 5, int(37 ** (10.0 / 3.0))),
          ("B", 29, 29 ** 6, 29 ** 3)]


# ---------------------------------------------------- base por z (máscaras)
class Base:
    _cache = {}

    @classmethod
    def get(cls, z):
        if z not in cls._cache:
            cls._cache[z] = Base(z)
        return cls._cache[z]

    def __init__(self, z):
        self.z = z
        P = [p for p in primes_upto(z - 1) if p > 2 and A % p != 0]
        self.P = P
        n = len(P)
        self.n = n
        N = 1 << n
        self.N = N
        val = [1] * N
        mu = [1] * N
        phi = [1] * N
        om = [0] * N
        for m in range(1, N):
            low = m & -m
            i = low.bit_length() - 1
            pm = m ^ low
            val[m] = val[pm] * P[i]
            mu[m] = -mu[pm]
            phi[m] = phi[pm] * (P[i] - 1)
            om[m] = om[pm] + 1
        pf = [tuple(P[i] for i in range(n - 1, -1, -1) if (m >> i) & 1)
              for m in range(N)]
        lam = [lam_of(t) for t in pf]
        self.val, self.mu, self.phi, self.om = val, mu, phi, om
        self.pf, self.lam = pf, lam
        self.a = [mu[m] / phi[m] for m in range(N)]
        self.Pi_full = math.prod(1.0 - 1.0 / (p - 1) for p in P)
        self.Qc = math.prod(1.0 - 1.0 / (p - 1) + 1.0 / (p - 1) ** 2
                            for p in P)
        self.M_full = math.prod(1.0 + 1.0 / ((p - 1) * (p - 2)) for p in P)
        fac = [(p - 1) / (p - 2) for p in P]
        Pi = [self.Pi_full] * N
        for m in range(1, N):
            low = m & -m
            i = low.bit_length() - 1
            Pi[m] = Pi[m ^ low] * fac[i]
        self.Pi = Pi
        self.order = sorted(range(N), key=lambda m: lam[m])
        self.lam_sorted = [lam[m] for m in self.order]
        self.mask_of = {val[m]: m for m in range(N)}
        self.idx = {p: i for i, p in enumerate(P)}

    def minp(self, m):
        return self.P[(m & -m).bit_length() - 1]

    def q0(self, m):
        """Menor primo de P que no divide a m (None si no hay)."""
        for i in range(self.n):
            if not (m >> i) & 1:
                return self.P[i]
        return None

    def Pi_direct(self, m):
        return math.prod(1.0 - 1.0 / (p - 1)
                         for i, p in enumerate(self.P) if not (m >> i) & 1)


class Level:
    def __init__(self, base, L):
        self.b = base
        self.L = L
        self.V = [m for m in range(base.N) if base.lam[m] < L]
        self.Vset = set(self.V)
        self.lnz = math.log(base.z)

    def visible(self):
        vis = set()
        for m in self.V:
            s = m
            while True:
                vis.add(s)
                if s == 0:
                    break
                s = (s - 1) & m
        return vis

    def Vr(self, rm):
        return [m ^ rm for m in self.V if m & rm == rm]

    def Vhat(self, rm):
        b = self.b
        r = b.val[rm]
        k = bisect_right(b.lam_sorted, (self.L - 1) // r)
        return [m for m in b.order[:k] if m & rm == 0]

    def inst(self, rm):
        b = self.b
        Vr = self.Vr(rm)
        Vh = self.Vhat(rm)
        a = b.a
        S = math.fsum(a[m] for m in Vr)
        Sh = math.fsum(a[m] for m in Vh)
        sVr, sVh = set(Vr), set(Vh)
        Add = sVr - sVh
        Rem = sVh - sVr
        Pi = b.Pi[rm]
        r = b.val[rm]
        return {"rm": rm, "r": r, "S": S, "Sh": Sh, "Pi": Pi,
                "delta": S - Pi, "E": S - Sh, "Vr": sVr, "Vh": sVh,
                "Add": Add, "Rem": Rem, "phi": b.phi[rm], "om": b.om[rm],
                "ur": math.log(self.L / r) / self.lnz}


def fused_failures(base, xm, rm, L):
    """Fallos de la condición fusionada de r*x: lista de (p, origen, j)."""
    pfm = base.pf[xm | rm]
    out = []
    prod = 1
    for j, p in enumerate(pfm, 1):
        prod *= p
        if j % 2 == 1 and p * p * prod >= L:
            out.append((p, "x" if (xm >> base.idx[p]) & 1 else "r", j))
    return out


def own_position(base, xm, p):
    return base.pf[xm].index(p) + 1


# ===================================================== barrido principal
def main_sweep():
    head("BARRIDO PRINCIPAL: 4 marcos x 60 niveles (docs 113/114/116/117/118/"
         "119/120a)")
    st = {
        "n_inst": 0, "n_add": 0, "n_rem": 0, "Epos": 0, "maxE": -1e9,
        "Eadd_max": -1e9, "Erem_max": -1e9, "idA": 0.0, "idR": 0.0,
        "vis_neq_V": 0, "hat_char_bad": 0,
        # 116a
        "Vhat_minus": 0, "rem_minus": 0, "fail_total": 0, "fail_le_qmin": 0,
        "fail_qmin": 0, "fail_R2_below": 0, "II_bad": 0, "pairsB": 0,
        "ratioB": 0.0, "residual": 0, "primes_in_rem": 0, "psi_ok": 0,
        "psi_bad": 0, "R1_odd": 0, "one_fail": 0,
        # 116b J1 / 117a K1
        "pairs": 0, "pair_partner_bad": 0, "viol1": 0, "viol2": 0,
        "viol3": 0, "maxRA": 0.0, "maxRA_at": None, "BAD": 0.0,
        "maxRR": 0.0,
        # 113 cross
        "deg_real": 0, "Amin": 0.0,
    }
    band = []          # registros de banda
    layer = []         # capa u_r en [1, 3/2) (120a)
    per_level = []     # sumas por nivel
    rho0 = lambda u: u ** (-u / 2.0)
    for (tag, z, D, y) in FRAMES:
        b = Base.get(z)
        levels = frame_levels(z, D, y)
        for L in levels:
            lv = Level(b, L)
            vis = lv.visible()
            if vis != lv.Vset:
                st["vis_neq_V"] += 1
            pl = {"tag": tag, "z": z, "L": L, "crude": 0.0, "crude_half": 0.0,
                  "real": 0.0, "nband_u": 0, "Qc": b.Qc}
            for rm in sorted(lv.V, key=lambda m: b.val[m]):
                if rm == 0:
                    continue
                I = lv.inst(rm)
                r = I["r"]
                st["n_inst"] += 1
                st["n_add"] += len(I["Add"])
                st["n_rem"] += len(I["Rem"])
                E = I["E"]
                if E > 1e-12:
                    st["Epos"] += 1
                st["maxE"] = max(st["maxE"], E / I["Pi"])
                # identidades de flujo (doc 114)
                EA = math.fsum(b.a[m] for m in I["Add"])
                ER = -math.fsum(b.a[m] for m in I["Rem"])
                st["Eadd_max"] = max(st["Eadd_max"], EA)
                st["Erem_max"] = max(st["Erem_max"], ER)
                # caracterización degenerada del modelo (117b §3.2)
                q0 = b.q0(rm)
                deg_hat = (len(I["Vh"]) == 1)
                if q0 is not None and deg_hat != (q0 ** 3 * r >= L):
                    st["hat_char_bad"] += 1
                if I["Vr"] == {0}:
                    st["deg_real"] += 1
                if I["delta"] > 0:
                    st["Amin"] = max(st["Amin"],
                                     I["delta"] / (I["Pi"] * rho0(I["ur"])))
                # ---- flujo Add: apareamiento x -> x/q_min (116b J1)
                RA = {}
                for x in I["Add"]:
                    if x == 0 or b.mu[x] != 1:
                        continue
                    q = b.minp(x)
                    yy = x ^ (1 << b.idx[q])
                    if not (yy in I["Add"] and b.mu[yy] == -1):
                        st["pair_partner_bad"] += 1
                        continue
                    st["pairs"] += 1
                    RA[yy] = RA.get(yy, 0.0) + 1.0 / (q - 1)
                    # cotas afiladas (*), (**), (***)
                    rpf = b.pf[rm]
                    yv = b.val[yy]
                    tq = sum(1 for p in rpf if p > q)
                    sq = prod_of([p for p in rpf if p < q])
                    if tq % 2 == 1 and not (q ** 3 * r * yv < L * sq):
                        st["viol1"] += 1
                    for p in rpf:
                        if p < q and sum(1 for pp in rpf if pp > p) % 2 == 0:
                            rge = prod_of([pp for pp in rpf if pp >= p])
                            if not (q * p * p * yv * rge < L):
                                st["viol2"] += 1
                    rho = max(rpf)
                    if rho < q < b.minp(yy) and not (q * rho ** 3 * yv < L):
                        st["viol3"] += 1
                # E_Add = -sum_{y in Add-} |a(y)| (1 - R_A(y))
                idA2 = EA - (-math.fsum(abs(b.a[yy]) * (1.0 - RA.get(yy, 0.0))
                                        for yy in I["Add"] if b.mu[yy] == -1))
                st["idA"] = max(st["idA"], abs(idA2))
                for yy, R in RA.items():
                    if R > st["maxRA"]:
                        st["maxRA"] = R
                        st["maxRA_at"] = (tag, z, L, r, b.val[yy])
                    st["BAD"] += abs(b.a[yy]) * max(R - 1.0, 0.0)
                # ---- lado Rem (116a, 114 §4)
                RR = {}
                for x in I["Rem"]:
                    if b.om[x] == 1:
                        st["primes_in_rem"] += 1
                    if b.mu[x] == -1:
                        st["rem_minus"] += 1
                        qm = b.minp(x)
                        fl = fused_failures(b, x, rm, L)
                        if len(fl) == 1:
                            st["one_fail"] += 1
                        if all(p <= qm for (p, o, j) in fl):
                            st["residual"] += 1
                        yy = x ^ (1 << b.idx[qm])
                        if yy in I["Rem"] and b.mu[yy] == 1:
                            st["psi_ok"] += 1
                            RR[yy] = RR.get(yy, 0.0) + 1.0 / (qm - 1)
                        else:
                            st["psi_bad"] += 1
                    for (p, o, j) in fused_failures(b, x, rm, L):
                        if o == "x" and own_position(b, x, p) % 2 == 1:
                            st["R1_odd"] += 1
                idR = ER - (-math.fsum(b.a[yy] * (1.0 - RR.get(yy, 0.0))
                                       for yy in I["Rem"] if b.mu[yy] == 1))
                st["idR"] = max(st["idR"], abs(idR))
                for yy, R in RR.items():
                    st["maxRR"] = max(st["maxRR"], R)
                # V-hat menos (omega impar): Lema A, Lema B, condición (II)
                rpf = b.pf[rm]
                for x in I["Vh"]:
                    if x == 0 or b.om[x] % 2 == 0:
                        continue
                    st["Vhat_minus"] += 1
                    qm = b.minp(x)
                    xv = b.val[x]
                    if not (r * xv * qm * qm < L):
                        st["II_bad"] += 1
                    fl = fused_failures(b, x, rm, L)
                    st["fail_total"] += len(fl)
                    for (p, o, j) in fl:
                        if p <= qm:
                            st["fail_le_qmin"] += 1
                        if p == qm:
                            st["fail_qmin"] += 1
                        if o == "r" and p < qm:
                            st["fail_R2_below"] += 1
                    for p in rpf:
                        if p < qm:
                            st["pairsB"] += 1
                            rge = prod_of([pp for pp in rpf if pp >= p])
                            st["ratioB"] = max(st["ratioB"],
                                               p * p * xv * rge / L)
                # ---- banda y capa
                ur = I["ur"]
                Pi = I["Pi"]
                d = I["delta"]
                t = d * (2 * Pi + d) / I["phi"]
                if ur < U0:
                    pl["crude"] += (1 - Pi * Pi) / I["phi"]
                    pl["real"] += max(t, 0.0)
                    pl["nband_u"] += 1
                    if Pi < 0.5:
                        pl["crude_half"] += (1 - Pi * Pi) / I["phi"]
                if ur < U0 and Pi < 0.5:
                    rec = {"tag": tag, "z": z, "L": L, "r": r, "rm": rm,
                           "om": I["om"], "ur": ur, "Pi": Pi, "S": I["S"],
                           "Sh": I["Sh"], "E": E, "delta": d,
                           "phi": I["phi"], "deg": deg_hat, "q0": q0,
                           "nVr": len(I["Vr"]), "EA": EA,
                           "ERsum": -ER, "Qc": b.Qc,
                           "Add_empty": not I["Add"],
                           "Rem_empty": not I["Rem"], "t": t}
                    band.append(rec)
            per_level.append(pl)
        print("  marco %s z=%d: %d niveles (t=%.1fs)"
              % (tag, z, len(levels), time.time() - T0))
    return st, band, per_level


def report_main(st, band, per_level):
    head("Doc 113/114 (control de la maquinaria) y doc 116a (REM-PERS-small)")
    check("K0.1", "instancias (L, r>1 visible) en 240 niveles = 4113",
          st["n_inst"] == 4113, "medido %d" % st["n_inst"])
    check("K0.2", "|Add| = 6509 y |Rem| = 1049 (doc 113 K6 / 114 K0)",
          st["n_add"] == 6509 and st["n_rem"] == 1049,
          "medido %d / %d" % (st["n_add"], st["n_rem"]))
    check("K0.3", "visibles = V+(L) (Lema 2 del doc 114, cierre por divisores)",
          st["vis_neq_V"] == 0, "niveles con diferencia: %d" % st["vis_neq_V"])
    check("K0.4", "E(r) <= 0 en las 4113 (DESP-SIGN en ventana) y por lados",
          st["Epos"] == 0 and st["Eadd_max"] <= 1e-12 and st["Erem_max"] <= 1e-12,
          "E>1e-12: %d | max E/Pi = %.1e | max E_Add = %.1e | max E_Rem = %.1e"
          % (st["Epos"], st["maxE"], st["Eadd_max"], st["Erem_max"]))
    check("K0.5", "identidades de flujo E_Add, E_Rem (doc 114) a precisión de "
          "máquina", st["idA"] < 1e-14 and st["idR"] < 1e-14,
          "err %.1e / %.1e" % (st["idA"], st["idR"]))
    check("K0.6", "A_min = 1,4273 y 1348 degenerados reales V^(r)={1} (doc 113)",
          abs(st["Amin"] - 1.4273) < 5e-5 and st["deg_real"] == 1348,
          "A_min = %.4f, degenerados reales = %d" % (st["Amin"], st["deg_real"]))
    check("116a.J0", "Rem^- = 39 en los marcos", st["rem_minus"] == 39,
          "medido %d" % st["rem_minus"])
    check("116a.J1", "x en Vhat^- examinados = 5568; condición (II) "
          "r x q_min^2 < L sin violaciones",
          st["Vhat_minus"] == 5568 and st["II_bad"] == 0,
          "medido %d, violaciones %d" % (st["Vhat_minus"], st["II_bad"]))
    check("116a.J2", "fallos fusionados en Vhat^- = 39, uno por x en Rem^-, "
          "0 en primos <= q_min (Lema A: 0 en q_min; Lema B: 0 R2 bajo q_min)",
          st["fail_total"] == 39 and st["one_fail"] == 39
          and st["fail_le_qmin"] == 0 and st["fail_qmin"] == 0
          and st["fail_R2_below"] == 0,
          "fallos %d, x con un fallo %d, <=q_min %d, en q_min %d, R2<q_min %d"
          % (st["fail_total"], st["one_fail"], st["fail_le_qmin"],
             st["fail_qmin"], st["fail_R2_below"]))
    check("116a.J3", "pares (x en Vhat^-, p|r, p<q_min) = 2423; razón máxima "
          "del lado fusionado = 0,7581 < 1",
          st["pairsB"] == 2423 and abs(st["ratioB"] - 0.7581) < 5e-5,
          "medido %d pares, razón máx %.4f" % (st["pairsB"], st["ratioB"]))
    check("116a.J4", "familia residual vacía; 0 primos en Rem; R1 nunca en "
          "posición propia impar",
          st["residual"] == 0 and st["primes_in_rem"] == 0 and st["R1_odd"] == 0,
          "residuales %d, primos en Rem %d, R1 en posición impar %d"
          % (st["residual"], st["primes_in_rem"], st["R1_odd"]))
    check("116a.J5", "aterrizaje psi(x) = x/q_min en Rem^+ : 39/39",
          st["psi_ok"] == 39 and st["psi_bad"] == 0,
          "%d bien, %d mal" % (st["psi_ok"], st["psi_bad"]))
    head("Doc 116b J1 (cotas afiladas) y doc 117a K1 (BAD en los marcos)")
    check("116b.J1", "1640 pares (y,q) admisibles; (*), (**), (***) sin "
          "violaciones; max R_A = 0,6000",
          st["pairs"] == 1640 and st["pair_partner_bad"] == 0
          and st["viol1"] == st["viol2"] == st["viol3"] == 0
          and abs(st["maxRA"] - 0.6) < 1e-12,
          "pares %d (socio fuera %d), viol %d/%d/%d, max R_A = %.4f en %s"
          % (st["pairs"], st["pair_partner_bad"], st["viol1"], st["viol2"],
             st["viol3"], st["maxRA"], st["maxRA_at"]))
    check("117a.K1", "BAD = 0 idéntico en las 4113 instancias; max R_R = 0,2500",
          st["BAD"] == 0.0 and abs(st["maxRR"] - 0.25) < 1e-12,
          "BAD total = %.3g, max R_R = %.4f" % (st["BAD"], st["maxRR"]))
    check("117b.P", "Vhat^(r) = {1} <=> q0^3 r >= L en TODAS las 4113 "
          "instancias (Proposición del doc 117b §3.2)", st["hat_char_bad"] == 0,
          "excepciones %d" % st["hat_char_bad"])
    # ---------------- doc 114 K6 / 116c J2
    head("Doc 114 K6 / 116c J2: banda cruda en los 240 niveles")
    mx_crude = max(p["crude"] / p["Qc"] for p in per_level)
    mx_half = {}
    for p in per_level:
        k = (p["tag"], p["z"])
        mx_half[k] = max(mx_half.get(k, 0.0), p["crude_half"] / p["Qc"])
    mx_real = max(p["real"] / p["Qc"] for p in per_level)
    mx_n = max(p["nband_u"] for p in per_level)
    check("114.K6", "u0=1,5: max cruda 0,6991 Q_c; parte Pi<1/2 0,4408; real "
          "0,3806; max #r 77",
          abs(mx_crude - 0.6991) < 5e-5 and abs(max(mx_half.values()) - 0.4408)
          < 5e-5 and abs(mx_real - 0.3806) < 5e-5 and mx_n == 77,
          "%.4f / %.4f / %.4f / %d" % (mx_crude, max(mx_half.values()),
                                       mx_real, mx_n))
    check("116c.J2", "máximo muestreado de la parte Pi<1/2: 0,1526 (z=29) y "
          "0,4408 (z=37); 0 en z=17",
          abs(mx_half[("A", 29)] - 0.1526) < 5e-5
          and abs(mx_half[("A", 37)] - 0.4408) < 5e-5
          and mx_half[("A", 17)] == 0.0,
          "A17 %.4f | A29 %.4f | A37 %.4f | B29 %.4f"
          % (mx_half[("A", 17)], mx_half[("A", 29)], mx_half[("A", 37)],
             mx_half[("B", 29)]))
    return mx_half


# ----------------------------------------------------- banda (117b/118b1)
def band_reports(band):
    head("Doc 117b: reducción, signo unilateral, sumas por marco, degenerados")
    n = len(band)
    Epos = sum(1 for b in band if b["E"] > 1e-12)
    viol = max((max(b["delta"], 0) - max(b["Sh"] - b["Pi"], 0)) for b in band)
    check("117b.K1", "246 instancias de banda; E>0 en 0; delta^+ <= "
          "(Shat-Pi)^+ (violación máx ~0)", n == 246 and Epos == 0
          and viol < 1e-15, "n=%d, E>0: %d, violación máx %.1e" % (n, Epos, viol))
    by = {}
    for b in band:
        k = min(b["om"], 3)
        c = by.setdefault(k, [0, 0, 0, 0.0, 0.0])
        dif = b["Sh"] - b["Pi"]
        if dif > 1e-15:
            c[0] += 1
            c[3] += dif / b["phi"]
        elif dif < -1e-15:
            c[1] += 1
            c[4] += -dif / b["phi"]
        else:
            c[2] += 1
    txt = "; ".join("w=%d: +%d -%d 0:%d masa+ %.4f masa- %.4f"
                    % (k, v[0], v[1], v[2], v[3], v[4]) for k, v in sorted(by.items()))
    ok = (by.get(1, [0])[0] == 56 and by.get(2, [0])[0] == 187
          and by.get(3, [0])[0] == 3
          and all(v[1] == 0 and v[2] == 0 for v in by.values())
          and abs(by[1][3] - 3.5374) < 5e-5 and abs(by[2][3] - 0.3211) < 5e-5)
    check("117b.K3", "Shat - Pi > 0 en 246/246; por omega 56/187/3 con masas "
          "3,5374 / 0,3211 / 0,0001", ok, txt)
    # sumas por marco: máximo sobre niveles
    lev = {}
    for b in band:
        k = (b["tag"], b["z"], b["L"])
        s = lev.setdefault(k, [0.0, 0.0, 0.0, b["Qc"]])
        s[0] += max(b["delta"], 0) / b["phi"]
        s[1] += max(b["Sh"] - b["Pi"], 0) / b["phi"]
        s[2] += abs(b["Sh"] - b["Pi"]) / b["phi"]
    fr = {}
    for (tag, z, L), s in lev.items():
        f = fr.setdefault((tag, z), [0.0, 0.0, 0.0])
        for i in range(3):
            f[i] = max(f[i], s[i] / s[3])
    exp = {("A", 29): (0.0198, 0.1028, 0.1028), ("A", 37): (0.0933, 0.2103, 0.2103),
           ("B", 29): (0.0128, 0.0128, 0.0128)}
    ok = all(all(abs(fr[k][i] - exp[k][i]) < 5e-5 for i in range(3)) for k in exp) \
        and ("A", 17) not in fr
    check("117b.K4", "máx por marco de sum delta^+/phi, (Shat-Pi)^+/phi, "
          "|Shat-Pi|/phi: A29 0,0198/0,1028/0,1028; A37 0,0933/0,2103/0,2103; "
          "B29 0,0128 x3; A17 vacío", ok,
          "; ".join("%s%d %.4f/%.4f/%.4f" % (k[0], k[1], *v) for k, v in sorted(fr.items())))
    deg = [b for b in band if b["deg"]]
    nondeg = [b for b in band if not b["deg"]]
    mx = max((b["Sh"] - b["Pi"]) / b["Pi"] for b in nondeg)
    check("117b.K5", "188/246 degenerados del modelo; máx (Shat-Pi)^+/Pi no "
          "degenerado = 0,692", len(deg) == 188 and abs(mx - 0.692) < 5e-4,
          "degenerados %d, máx %.4f" % (len(deg), mx))
    # ---- 118b1: Teorema A y estratos
    head("Doc 118b1: Teorema A, histogramas, estratos")
    hist = {}
    bad_A = 0
    for b in deg:
        k = (b["tag"], b["z"])
        h = hist.setdefault(k, {})
        h[b["om"]] = h.get(b["om"], 0) + 1
        if b["om"] <= 2 and not (b["L"] < b["q0"] ** 9):
            bad_A += 1
    exp_h = {("A", 29): {1: 1, 2: 9}, ("A", 37): {1: 31, 2: 47}, ("B", 29): {2: 100}}
    check("118b1.J1", "histogramas omega de los 188 degenerados {1:1,2:9}, "
          "{1:31,2:47}, {2:100}; todos con omega<=2 cumplen L < q0^9",
          hist == exp_h and bad_A == 0, "medido %s, violaciones %d"
          % (dict(sorted(hist.items())), bad_A))
    # sumas degeneradas visibles por estrato: máximo por nivel
    lev = {}
    for b in deg:
        k = (b["tag"], b["z"], b["L"], b["q0"])
        lev[k] = lev.get(k, 0.0) + 1.0 / b["phi"]
    mxs = {}
    for (tag, z, L, q0), s in lev.items():
        kk = (tag, z, q0)
        Qc = Base.get(z).Qc
        mxs[kk] = max(mxs.get(kk, 0.0), s / Qc)
    strata = sorted({k[2] for k in mxs})
    ok = (strata == [5] and abs(mxs[("A", 29, 5)] - 0.1993) < 5e-5
          and abs(mxs[("A", 37, 5)] - 0.3554) < 5e-5
          and abs(mxs[("B", 29, 5)] - 0.0252) < 5e-5)
    check("118b1.J2", "suma degenerada visible máx por nivel (1/phi) solo en "
          "el estrato 5: 0,1993 / 0,3554 / 0,0252 Q_c", ok,
          "estratos presentes %s; %s" % (strata, {k: round(v, 4) for k, v in mxs.items()}))
    return deg


# ------------------------------------------------------------ 119b (SLAB-V)
def slab_reports():
    head("Doc 119b: reindexación divisor-cofactor de (SLAB-V), paso a paso")
    out = {}
    viol_phi = viol_c = viol_cov = viol_N = 0
    worst_err = 0.0
    for (tag, z, D, y) in FRAMES:
        b = Base.get(z)
        for L in frame_levels(z, D, y):
            lv = Level(b, L)
            Vs = lv.Vset
            for q0 in (5, 7, 11):
                B = q0 ** 3
                Pmask = sum(1 << b.idx[p] for p in b.P if p < q0)
                true_s = Psum = 0.0
                Ns = []
                # r degenerados visibles en banda, estrato q0
                for rm in lv.V:
                    if rm == 0:
                        continue
                    r = b.val[rm]
                    if b.q0(rm) != q0 or not (q0 ** 3 * r >= L):
                        continue
                    ur = math.log(L / r) / lv.lnz
                    if not (ur < U0 and b.Pi[rm] < 0.5):
                        continue
                    true_s += 1.0 / b.phi[rm]
                    # testigos d = c r en V: cofactores
                    cofs = [m ^ rm for m in lv.V if m & rm == rm]
                    for c in cofs:
                        if not b.val[c] < B:
                            viol_c += 1
                        if b.phi[c] * b.phi[rm] != b.phi[c | rm]:
                            viol_phi += 1
                    if not cofs:
                        viol_cov += 1
                    # N(r) por la definición del doc: #{c < q0^3: (c,r)=1, cr en V}
                    Ndef = sum(1 for c in range(b.N)
                               if c & rm == 0 and b.val[c] < B and (c | rm) in Vs)
                    if Ndef != len(cofs):
                        viol_N += 1
                    Ns.append(Ndef)
                    Psum += Ndef / b.phi[rm]
                # W(q0) = sum_c U_c
                W = 0.0
                Uc = {}
                for dm in lv.V:
                    if b.val[dm] * B < L:
                        continue
                    s = dm
                    while True:
                        c = s
                        if b.val[c] < B and c & Pmask == 0 and b.q0(dm ^ c) == q0:
                            w = b.phi[c] / b.phi[dm]
                            W += w
                            Uc[b.val[c]] = Uc.get(b.val[c], 0.0) + w
                        if s == 0:
                            break
                        s = (s - 1) & dm
                worst_err = max(worst_err, abs(W - math.fsum(Uc.values())))
                k = (tag, z, q0)
                o = out.setdefault(k, {"true": 0.0, "P": 0.0, "W": 0.0,
                                       "ncof": 0, "Ns": [], "U5": 0.0, "U1": 0.0})
                o["true"] = max(o["true"], true_s / b.Qc)
                o["P"] = max(o["P"], Psum / b.Qc)
                o["W"] = max(o["W"], W / b.Qc)
                o["ncof"] = max(o["ncof"], sum(1 for v in Uc.values() if v > 0))
                o["Ns"] += Ns
                o["U5"] = max(o["U5"], Uc.get(5, 0.0) / b.Qc)
                o["U1"] = max(o["U1"], Uc.get(1, 0.0) / b.Qc)
    check("119b.K0", "phi(c)phi(r) = phi(d) en todo par; todo cofactor de "
          "testigo < q0^3; todo r degenerado visible tiene testigo",
          viol_phi == 0 and viol_c == 0 and viol_cov == 0,
          "violaciones %d / %d / %d" % (viol_phi, viol_c, viol_cov))
    check("119b.N", "N(r) del doc = |V^(r)| (el número de testigos es el "
          "tamaño del soporte desplazado REAL)", viol_N == 0,
          "discrepancias %d" % viol_N)
    allN = out[("A", 29, 5)]["Ns"]
    allN37 = out[("A", 37, 5)]["Ns"]
    allNB = out[("B", 29, 5)]["Ns"]
    mean = lambda v: sum(v) / len(v) if v else float("nan")
    print("    N(r): max %d / %d / %d; media %.2f / %.2f / %.2f (A29/A37/B29)"
          % (max(allN), max(allN37), max(allNB), mean(allN), mean(allN37), mean(allNB)))
    t = {k: (round(v["true"], 4), round(v["P"], 4), round(v["W"], 2))
         for k, v in sorted(out.items())}
    print("    (verdadera, P, W)/Q_c por (marco, z, q0): %s" % t)
    ok = (abs(out[("A", 29, 5)]["true"] - 0.1993) < 5e-5
          and abs(out[("A", 29, 5)]["P"] - 0.5979) < 5e-5
          and abs(out[("A", 37, 5)]["P"] - 0.7109) < 5e-5
          and abs(out[("B", 29, 5)]["P"] - 0.0252) < 5e-5
          and abs(out[("A", 29, 5)]["W"] - 16.26) < 5e-3
          and abs(out[("A", 37, 5)]["W"] - 16.09) < 5e-3
          and abs(out[("B", 29, 5)]["W"] - 5.46) < 5e-3
          and max(allN + allN37 + allNB) == 3)
    check("119b.K1", "tabla: verdadera 0,1993/0,3554/0,0252; P 0,5979/0,7109/"
          "0,0252; W 16,26/16,09/5,46; N(r) <= 3", ok, "ver línea anterior")
    check("119b.K2", "intercambio W = sum_c U_c exacto", worst_err < 1e-12,
          "error máx %.1e" % worst_err)
    print("    cofactores que aparecen en W(5): %s; en W(11): %s"
          % ([out[(f[0], f[1], 5)]["ncof"] for f in FRAMES],
             [out[(f[0], f[1], 11)]["ncof"] for f in FRAMES]))
    print("    U_5 (q0=5) máx/Q_c: %s;  parte c=1 de W(5): %s"
          % ([round(out[(f[0], f[1], 5)]["U5"], 3) for f in FRAMES],
             [round(out[(f[0], f[1], 5)]["U1"], 3) for f in FRAMES]))
    preds = []
    for z in (29, 37):
        b = Base.get(z)
        Ac = math.prod((1 - 1.0 / p) if p == 5 else (1 + 1.0 / (p * (p - 1)))
                       for p in b.P)
        preds.append((z, b.Pi_full * Ac * math.log(125) / b.Qc))
    print("    predicción del modelo Pi_full A(5) ln(5^3) / Q_c: %s (doc: 3,52)"
          % ["z=%d: %.2f" % p for p in preds])
    return out


# ------------------------------------------------- 120a / 119a (capa alta)
def stress53_levels(cap=4, near_u=2.0):
    z = 53
    D, y = z ** 5, int(z ** (10.0 / 3.0))
    lv = frame_levels(z, D, y)
    lnz = math.log(z)
    lv.sort(key=lambda L: abs(math.log(L) / lnz - near_u))
    return sorted(lv[:cap])


def band_records_for(z, levels, tag):
    b = Base.get(z)
    recs = []
    for L in levels:
        lv = Level(b, L)
        for rm in lv.V:
            if rm == 0:
                continue
            r = b.val[rm]
            ur = math.log(L / r) / lv.lnz
            if not (ur < U0 and b.Pi[rm] < 0.5):
                continue
            I = lv.inst(rm)
            recs.append({"tag": tag, "z": z, "L": L, "r": r, "rm": rm,
                         "om": I["om"], "ur": ur, "Pi": I["Pi"], "S": I["S"],
                         "Sh": I["Sh"], "E": I["E"], "delta": I["delta"],
                         "phi": I["phi"], "deg": len(I["Vh"]) == 1,
                         "q0": b.q0(rm), "nVr": len(I["Vr"]),
                         "EA": math.fsum(b.a[m] for m in I["Add"]),
                         "ERsum": math.fsum(b.a[m] for m in I["Rem"]),
                         "Qc": b.Qc, "Add_empty": not I["Add"],
                         "Rem_empty": not I["Rem"]})
    return recs


def layer_reports(band):
    head("Doc 120a: (BAND-delta-HI), capa u_r en [1, 3/2)")
    st53 = stress53_levels()
    print("    niveles de estrés z=53 (4 más cercanos a u_L=2): %s" % st53)
    recs = list(band) + band_records_for(53, st53, "S")
    lay = [b for b in recs if 1.0 <= b["ur"] < U0]
    # tabla por marco
    lev = {}
    for b in lay:
        k = (b["tag"], b["z"], b["L"])
        s = lev.setdefault(k, [0.0, 0.0, 0.0, 0.0, b["Qc"]])
        s[0] += max(b["delta"], 0) / b["phi"]
        s[1] += max(b["Sh"] - b["Pi"], 0) / b["phi"]
        s[2] += abs(b["E"]) / b["phi"]
        if 1.4 <= b["ur"] < 1.5:
            s[3] += max(b["delta"], 0) / b["phi"]
    fr = {}
    for (tag, z, L), s in lev.items():
        f = fr.setdefault((tag, z), [0.0] * 4)
        for i in range(4):
            f[i] = max(f[i], s[i] / s[4])
    exp = {("A", 29): (0.0198, 0.1028, 0.0830, 0.0198),
           ("A", 37): (0.0933, 0.2103, 0.1244, 0.0933),
           ("B", 29): (0.0128, 0.0128, 0.0000, 0.0073),
           ("S", 53): (0.0867, 0.2695, 0.1828, 0.0491)}
    ok = all(all(abs(fr[k][i] - exp[k][i]) < 5e-5 for i in range(4)) for k in exp)
    check("120a.T", "tabla de capa: A29 0,0198/0,1028/0,0830/0,0198; A37 "
          "0,0933/0,2103/0,1244/0,0933; B29 0,0128/0,0128/0/0,0073; "
          "estrés53 0,0867/0,2695/0,1828/0,0491", ok,
          "; ".join("%s%d %s" % (k[0], k[1], "/".join("%.4f" % x for x in v))
                    for k, v in sorted(fr.items())))
    nd = [b for b in lay if not b["deg"] and b["Sh"] - b["Pi"] > 0]
    cover = sum(1 for b in nd if abs(b["E"]) >= (b["Sh"] - b["Pi"]))
    ratios = [abs(b["E"]) / (b["Sh"] - b["Pi"]) for b in nd]
    zero_struct = sum(1 for b in nd if b["Add_empty"] and b["Rem_empty"])
    zero_float = sum(1 for b in nd if b["E"] == 0.0)
    nz = sorted(abs(b["E"]) / (b["Sh"] - b["Pi"]) for b in nd
                if not (b["Add_empty"] and b["Rem_empty"]))
    med_up = nz[len(nz) // 2] if nz else float("nan")
    med_true = (nz[(len(nz) - 1) // 2] + nz[len(nz) // 2]) / 2 if nz else float("nan")
    ok = (len(nd) == 67 and cover == 0 and abs(max(ratios) - 0.9930) < 5e-5
          and zero_struct == 39 and len(nz) == 28)
    check("120a.K1", "67 no degeneradas con exceso; cobertura 0/67; máx "
          "|E|/(Shat-Pi)^+ = 0,9930; E = 0 en 39/67 (58,2 %); 28 con E != 0",
          ok, "n=%d, cobertura %d, máx %.4f, E=0 estructural %d (float %d) = "
          "%.1f %%, E!=0: %d" % (len(nd), cover, max(ratios), zero_struct,
                                 zero_float, 100.0 * zero_struct / len(nd), len(nz)))
    # Diagnóstico: E por la ruta de C(r) (S = phi(r)C(r)/mu(r), sumas en
    # orden creciente de máscara), como la escribe el doc 109, y Shat en
    # el mismo orden: ¿deja ruido de coma flotante donde E = 0 estructural?
    noisy = []
    cache = {}
    for b in nd:
        if not (b["Add_empty"] and b["Rem_empty"]):
            continue
        bs = Base.get(b["z"])
        key = (b["z"], b["L"])
        if key not in cache:
            cache[key] = Level(bs, b["L"])
        lv = cache[key]
        rm = b["rm"]
        c = 0.0
        for m in sorted(lv.V):
            if m & rm == rm:
                c += bs.a[m]
        S_C = c * bs.phi[rm] / bs.mu[rm]
        sh = 0.0
        for m in sorted(lv.Vhat(rm)):
            sh += bs.a[m]
        if S_C - sh != 0.0:
            noisy.append((b["tag"], b["z"], b["L"], b["r"], S_C - sh))
    print("    E=0 estructural pero E != 0 en coma flotante por la ruta de C: "
          "%d caso(s) %s" % (len(noisy), noisy[:3]))
    check("120a.K1n", "la discrepancia 39 frente a 40 se explica por ruido de "
          "coma flotante (E ~ 1e-16 contado como E != 0)",
          zero_struct - len(noisy) == 39 and len(nz) + len(noisy) == 28,
          "estructural %d - ruidosos %d = %d; E!=0 %d + %d = %d"
          % (zero_struct, len(noisy), zero_struct - len(noisy), len(nz),
             len(noisy), len(nz) + len(noisy)))
    check("120a.K1b", "mediana del afeitado donde E != 0 = 0,6271 (62,7 %)",
          abs(med_up - 0.6271) < 5e-5 or abs(med_true - 0.6271) < 5e-5,
          "mediana superior (sorted[n//2]) = %.4f, mediana estándar = %.4f"
          % (med_up, med_true))
    # K2: 43 instancias A29, estrés53, A37 no degeneradas; E_Rem = 0
    sub = [b for b in lay if not b["deg"]
           and (b["tag"], b["z"]) in (("A", 29), ("S", 53), ("A", 37))]
    remnz = sum(1 for b in sub if not b["Rem_empty"])
    ra = sorted(abs(b["EA"]) / (b["Sh"] - b["Pi"]) for b in sub
                if b["Sh"] - b["Pi"] > 0)
    check("120a.K2", "43 instancias no degeneradas (A29, estrés, A37): Rem "
          "vacío en todas (E 100 % lado Add); mediana |E_Add|/(Shat-Pi) = 0,5581",
          len(sub) == 43 and remnz == 0 and abs(ra[len(ra) // 2] - 0.5581) < 5e-5,
          "n=%d, con Rem no vacío %d, mediana (sorted[n//2]) %.4f, n_dif>0 = %d"
          % (len(sub), remnz, ra[len(ra) // 2], len(ra)))
    # Afirmación de redacción: degenerado del modelo => E = 0, delta = 1 - Pi
    degl = [b for b in lay if b["deg"]]
    bad = [b for b in degl if not b["Add_empty"]]
    ex = bad[0] if bad else None
    check("120a.D", "¿todo degenerado del modelo (Vhat={1}) de la capa tiene "
          "E = 0 y delta = 1 - Pi? (frase de 120a §1)", not bad,
          "degenerados de capa %d, con E != 0: %d%s"
          % (len(degl), len(bad), (" | ej.: %s%d L=%d r=%d |V^(r)|=%d E=%.4f "
                                   "delta=%.4f 1-Pi=%.4f" % (ex["tag"], ex["z"], ex["L"], ex["r"],
                                                             ex["nVr"], ex["E"], ex["delta"],
                                                             1 - ex["Pi"])) if ex else ""))
    return recs, st53


def profile_119a(recs):
    head("Doc 119a: perfil (Shat-Pi)/Pi en la capa alta y sumas por marco")
    lay = [b for b in recs if b["ur"] >= 1.0]
    mx = max((b["Sh"] - b["Pi"]) / b["Pi"] for b in lay)
    slices = [(1.00, 1.25), (1.25, 1.40), (1.40, 1.45), (1.45, 1.50)]
    prof = {}
    for (lo, hi) in slices:
        for fr in sorted({(b["tag"], b["z"]) for b in lay}):
            sel = [b for b in lay if (b["tag"], b["z"]) == fr and lo <= b["ur"] < hi]
            if not sel:
                continue
            w = sum(b["Pi"] / b["phi"] for b in sel)
            m = sum((b["Sh"] - b["Pi"]) / b["phi"] for b in sel) / w
            prof.setdefault((lo, hi), []).append(m)
    txt = "; ".join("[%.2f,%.2f): %.2f-%.2f" % (k[0], k[1], min(v), max(v))
                    for k, v in sorted(prof.items()))
    check("119a.K1a", "máx (Shat-Pi)/Pi en u_r>=1 = 1,2947; perfil "
          "ponderado Theta(1)", abs(mx - 1.2947) < 5e-5, "máx %.4f | %s" % (mx, txt))
    lev = {}
    for b in recs:
        k = (b["tag"], b["z"], b["L"])
        s = lev.setdefault(k, [0.0] * 6 + [b["Qc"]])
        dp = max(b["Sh"] - b["Pi"], 0) / b["phi"]
        s[0] += dp
        if b["ur"] >= 1.25:
            s[1] += dp
        if b["ur"] >= 1.40:
            s[2] += dp
        if b["om"] == 1:
            s[3] += dp
        if b["om"] == 2:
            s[4] += dp
        s[5] += b["Pi"] / b["phi"]
    fr = {}
    for (tag, z, L), s in lev.items():
        f = fr.setdefault((tag, z), [0.0] * 6)
        for i in range(6):
            f[i] = max(f[i], s[i] / s[6])
    exp = {("A", 29): (0.1028, 0.1028, 0.1028, 0.1028, 0.0128, 0.0965),
           ("A", 37): (0.2103, 0.2103, 0.1182, 0.2103, 0.0572, 0.2695),
           ("B", 29): (0.0128, 0.0128, 0.0073, 0.0000, 0.0128, 0.0124),
           ("S", 53): (0.3172, 0.2695, 0.1209, 0.2695, 0.0477, 0.3911)}
    ok = all(all(abs(fr[k][i] - exp[k][i]) < 5e-5 for i in range(6)) for k in exp)
    check("119a.K1b", "tabla de 119a §2 (banda total, capas >=1,25 y >=1,40, "
          "omega=1, omega=2, mayorante sum Pi/phi) en A29, A37, B29, estrés53", ok,
          "; ".join("%s%d %s" % (k[0], k[1], "/".join("%.4f" % x for x in v))
                    for k, v in sorted(fr.items())))


# ------------------------------------------------------- 118b2 J1 (25 niv)
def nondeg_25():
    head("Doc 118b2 J1: partición degenerado / no degenerado con 25 niveles")
    res = {}
    for (tag, z, D, y) in FRAMES:
        levels = frame_levels(z, D, y, nsamp=25)
        recs = band_records_for(z, levels, tag)
        f = [0.0] * 5
        cnt = [0, 0, 0]
        per = {}
        for b in recs:
            s = per.setdefault(b["L"], [0.0] * 5)
            dp = max(b["Sh"] - b["Pi"], 0) / b["phi"]
            pp = b["Pi"] / b["phi"]
            if b["deg"]:
                s[4] += dp
                cnt[2] += 1
            elif b["ur"] < 1:
                s[0] += dp
                s[1] += pp
                cnt[0] += 1
            else:
                s[2] += dp
                s[3] += pp
                cnt[1] += 1
        Qc = Base.get(z).Qc
        for s in per.values():
            for i in range(5):
                f[i] = max(f[i], s[i] / Qc)
        res[(tag, z)] = (f, cnt)
    exp = {("A", 29): ((0, 0, 0.0530, 0.0965, 0.0035), (0, 1, 1)),
           ("A", 37): ((0, 0, 0.1570, 0.2695, 0.1822), (0, 13, 28)),
           ("B", 29): ((0, 0, 0.0037, 0.0072, 0.0128), (0, 10, 42))}
    ok = all(all(abs(res[k][0][i] - exp[k][0][i]) < 5e-5 for i in range(5))
             and tuple(res[k][1]) == exp[k][1] for k in exp)
    check("118b2.J1", "celda no degenerada u_r<1 vacía; tabla 0,0530/0,0965/"
          "0,0035 (A29), 0,1570/0,2695/0,1822 (A37), 0,0037/0,0072/0,0128 (B29); "
          "instancias 0/1/1, 0/13/28, 0/10/42", ok,
          "; ".join("%s%d %s %s" % (k[0], k[1], "/".join("%.4f" % x for x in v[0]), v[1])
                    for k, v in sorted(res.items())))


# --------------------------------------- niveles representativos L=z^u
def rep_levels():
    head("Docs 116c §3 / 117b K6 / 118a J1 / 118b2 J2: niveles L = floor(z^u)")
    rows = {}
    for z in (29, 37, 53):
        b = Base.get(z)
        for (uname, u) in (("5/3", 5.0 / 3.0), ("2", 2.0), ("3", 3.0), ("4", 4.0)):
            if z == 53 and u > 2.5:
                continue
            L = int(z ** u)
            lv = Level(b, L)
            crude = real = 0.0
            regpi = [0.0, 0.0, 0.0]
            for rm in lv.V:
                if rm == 0:
                    continue
                r = b.val[rm]
                ur = math.log(L / r) / lv.lnz
                Pi = b.Pi[rm]
                if not (ur < U0 and Pi < 0.5):
                    continue
                I = lv.inst(rm)
                d = I["delta"]
                crude += (1 - Pi * Pi) / I["phi"]
                real += max(d * (2 * Pi + d) / I["phi"], 0.0)
                regpi[min(I["om"], 3) - 1] += Pi / I["phi"]
            rows[(z, uname)] = (L, crude / b.Qc, real / b.Qc,
                                [x / b.Qc for x in regpi])
    for k, v in sorted(rows.items()):
        print("    z=%d u=%s L=%d: cruda %.4f, real %.4f, sum Pi/phi por omega "
              "%s" % (k[0], k[1], v[0], v[1], v[2], ["%.4f" % x for x in v[3]]))
    e116 = {(37, "5/3"): (0.2710, 0.1155), (37, "2"): (0.4408, 0.0948),
            (37, "3"): (0.0791, 0.0493), (37, "4"): (0.0006, 0.0000),
            (29, "3"): (0.0191, 0.0191)}
    ok = all(abs(rows[k][1] - v[0]) < 5e-5 and abs(rows[k][2] - v[1]) < 5e-5
             for k, v in e116.items())
    check("116c.T", "tabla de ventana del doc 116c §3 (cruda/real a L=floor(z^u))",
          ok, "")
    e118 = {(37, "5/3"): (0.1733, 0, 0), (37, "2"): (0.2695, 0, 0),
            (53, "5/3"): (0.4609, 0, 0), (53, "2"): (0.1746, 0.0436, 0),
            (29, "5/3"): (0, 0, 0), (29, "2"): (0, 0, 0)}
    ok = all(all(abs(rows[k][3][i] - v[i]) < 5e-5 for i in range(3))
             for k, v in e118.items())
    tot = [sum(rows[k][3]) for k in e118]
    check("118a.J1", "(REG-Pi) medido: sum_banda Pi/phi = 0,17-0,46 Q_c "
          "(tabla de 118a §5 por omega)", ok,
          "totales %s" % ["%.4f" % t for t in tot])
    # J4 de 116c: p = 7, 11 a z=37, L=37^2
    b = Base.get(37)
    L = 37 ** 2
    lv = Level(b, L)
    txt = []
    ok = True
    expJ4 = {7: (0.4875, 0.6500, 0.1625, 0.1271, 0.0308), 11: (0.4514, 0.5833, 0.1320, 0.0796, 0.0137)}
    for p in (7, 11):
        I = lv.inst(1 << b.idx[p])
        Pi, S, d = I["Pi"], I["S"], I["delta"]
        cr = (1 - Pi * Pi) / I["phi"]
        tt = d * (2 * Pi + d) / I["phi"]
        txt.append("p=%d Pi=%.4f S=%.4f delta=%.4f cruda=%.4f t=%.4f factor %.1f"
                   % (p, Pi, S, d, cr, tt, cr / tt))
        ok &= all(abs(a_ - b_) < 5e-5 for a_, b_ in zip((Pi, S, d, cr, tt), expJ4[p]))
    check("116c.J4", "cancelación que el conteo no ve (z=37, L=37^2)", ok, "; ".join(txt))
    # 118b2 J2b / 118a J4: familia x_p en (1/2, 2/3) a z=53, u_L=2
    b = Base.get(53)
    L = 53 ** 2
    lv = Level(b, L)
    fam = [(p, q) for p in b.P for q in b.P
           if 0.5 < math.log(p) / math.log(53) < 2.0 / 3.0 and 5 <= q < p]
    inb = []
    vis = 0
    for (p, q) in fam:
        rm = (1 << b.idx[p]) | (1 << b.idx[q])
        if rm in lv.Vset:
            vis += 1
        ur = math.log(L / (p * q)) / lv.lnz
        if rm in lv.Vset and ur < U0 and b.Pi[rm] < 0.5:
            inb.append((p * q, round(ur, 3), len(lv.Vhat(rm)) == 1))
    g = sum(1.0 / ((p - 2) * (q - 2)) for (p, q) in fam
            if (p * q) in [x[0] for x in inb])
    strict = [(p, q) for (p, q) in fam if 11 / 20 < math.log(p) / math.log(53) < 3 / 5]
    ok = (len(fam) == 5 and vis == 5 and sorted(x[0] for x in inb) == [77, 91, 143]
          and all(x[2] for x in inb) and abs(g - 0.0505) < 5e-5 and not strict
          and abs(g / b.M_full - 0.0436) < 5e-5)
    check("118a.J4", "familia a z=53: 5 elementos, 5/5 visibles, en banda 77, "
          "91, 143 (degenerados, u_r en [0,750, 0,906]); masa g = 0,0505 = "
          "0,0436 M_full; ventana estricta (11/20, 3/5) vacía", ok,
          "en banda %s, g=%.4f, g/M_full=%.4f" % (inb, g, g / b.M_full))
    return rows


# ---------------------------------------------- 117b K2: doble conteo y W_L
def double_count():
    head("Doc 117b K2: identidad de doble conteo y W_L(385)")
    out = []
    ok = True
    for (z, L) in ((29, 29 ** 2), (29, 29 ** 3), (37, 410)):
        b = Base.get(z)
        lv = Level(b, L)
        lhs = 0.0
        signed = 0.0
        for rm in range(b.N):
            Vh = lv.Vhat(rm)
            Sh = math.fsum(b.a[m] for m in Vh)
            lhs += (b.Pi[rm] - Sh) / b.phi[rm]
            signed += (Sh - b.Pi[rm]) / b.phi[rm]
        rhs = 0.0
        for dm in range(b.N):
            lam = b.lam[dm]
            W = math.fsum(1.0 / b.phi[rm] for rm in range(b.N)
                          if rm & dm == 0 and b.val[rm] * lam >= L)
            rhs += b.a[dm] * W
        out.append("z=%d L=%d: |lhs-rhs|=%.1e, sum_r (Shat-Pi)/phi = %+.4f Q_c"
                   % (z, L, abs(lhs - rhs), signed / b.Qc))
        ok &= abs(lhs - rhs) < 1e-12
        if (z, L) == (29, 841):
            ok &= abs(signed / b.Qc - 1.2722) < 5e-5
        if (z, L) == (37, 410):
            ok &= abs(signed / b.Qc - 1.5653) < 5e-5
    b = Base.get(29)
    dm = b.mask_of[385]
    Ws = []
    for L in (29 ** 2, 29 ** 3):
        Ws.append(math.fsum(1.0 / b.phi[rm] for rm in range(b.N)
                            if rm & dm == 0 and b.val[rm] * b.lam[dm] >= L))
    ok &= b.lam[dm] == 9625 and abs(Ws[0] - 1.2702) < 5e-5 and abs(Ws[1] - 0.2702) < 5e-5
    check("117b.K2", "identidad exacta; signo +1,2722 Q_c y +1,5653 Q_c; "
          "lam(385)=9625, W_L(385) = 1,2702 (L=29^2) y 0,2702 (L=29^3)", ok,
          "; ".join(out) + "; W = %.4f / %.4f" % tuple(Ws))


# ---------------------------------------------- estrés z=53 (doc 114 K7)
def stress_114():
    head("Doc 114 K7 / 116a (estrés z=53, 10 niveles, <=120 conductores)")
    z = 53
    b = Base.get(z)
    D, y = z ** 5, int(z ** (10.0 / 3.0))
    n_inst = remm = vhm = pairs = 0
    ratio = 0.0
    fails_le = 0
    psi_ok = 0
    for L in frame_levels(z, D, y, nsamp=10):
        lv = Level(b, L)
        vis = sorted((m for m in lv.V if m != 0), key=lambda m: b.val[m])
        if len(vis) > 120:
            idx = sorted({round(i * (len(vis) - 1) / 119.0) for i in range(120)})
            vis = [vis[i] for i in idx]
        for rm in vis:
            n_inst += 1
            r = b.val[rm]
            I = lv.inst(rm)
            rpf = b.pf[rm]
            for x in I["Rem"]:
                if b.mu[x] == -1:
                    remm += 1
                    qm = b.minp(x)
                    if (x ^ (1 << b.idx[qm])) in I["Rem"]:
                        psi_ok += 1
            for x in I["Vh"]:
                if x == 0 or b.om[x] % 2 == 0:
                    continue
                vhm += 1
                qm = b.minp(x)
                for (p, o, j) in fused_failures(b, x, rm, L):
                    if p <= qm:
                        fails_le += 1
                for p in rpf:
                    if p < qm:
                        pairs += 1
                        rge = prod_of([pp for pp in rpf if pp >= p])
                        ratio = max(ratio, p * p * b.val[x] * rge / L)
    check("116a.J0s", "estrés: 185 instancias, Rem^- = 36, Vhat^- = 498, "
          "221 pares, razón máx 0,4916, 0 fallos <= q_min, psi 36/36",
          n_inst == 185 and remm == 36 and vhm == 498 and pairs == 221
          and abs(ratio - 0.4916) < 5e-5 and fails_le == 0 and psi_ok == 36,
          "%d / %d / %d / %d / %.4f / %d / %d" % (n_inst, remm, vhm, pairs,
                                                  ratio, fails_le, psi_ok))


# ======================================================= z = 331 (117a K2)
def dfs_union_r5(z, L):
    """Elementos d' (primos de P\\{5}) de Vhat^(5) U V^(5) al nivel L,
    por DFS sobre factorizaciones decrecientes (la unión es cerrada por
    divisores). Devuelve lista de (pf, d, mu, phi, inV, inhat)."""
    P = [p for p in primes_upto(z - 1) if p > 5]
    Pd = sorted(P, reverse=True)
    out = []
    stack = [((), 1, 1, 1, 0, 0)]   # pf, d, mu, phi, lam_own, next idx
    while stack:
        pf, d, mu, ph, lo, i = stack.pop()
        m = len(pf)
        inV = (5 * d < L) and lo < L and (m % 2 == 1 or 125 * d < L)
        inhat = 5 * max(d, lo) < L
        out.append((pf, d, mu, ph, inV, inhat))
        for k in range(i, len(Pd)):
            p = Pd[k]
            nd = d * p
            if 5 * nd >= L:
                continue
            nlo = lo
            if (m + 1) % 2 == 1:
                t = p * p * nd
                if t > nlo:
                    nlo = t
            if nlo >= L:
                continue
            ninV = (nlo < L) and ((m + 1) % 2 == 1 or 125 * nd < L)
            ninhat = 5 * max(nd, nlo) < L
            if ninV or ninhat:
                stack.append((pf + (p,), nd, -mu, ph * (p - 1), nlo, k + 1))
    return out


def z331():
    head("Doc 117a K2 / 116b J2: el marco del contraejemplo (z=331, r=5)")
    z = 331
    D = z ** 5
    yA = int(z ** (10.0 / 3.0))
    M37 = mertens(37)
    res = {}
    for (L, s) in ((46948396, 84629), (79051268, 50261), (113296524, 35069)):
        okL = (D // s == L and is_prime(s) and z <= s < yA)
        U = dfs_union_r5(z, L)
        Vr = {}
        Vh = {}
        info = {}
        for (pf, d, mu, ph, inV, inhat) in U:
            info[d] = (pf, mu, ph)
            if inV:
                Vr[d] = True
            if inhat:
                Vh[d] = True
        Add = [d for d in Vr if d not in Vh]
        Rem = [d for d in Vh if d not in Vr]
        a = lambda d: info[d][1] / info[d][2]
        EA = math.fsum(a(d) for d in Add)
        ER = -math.fsum(a(d) for d in Rem)
        S5 = math.fsum(a(d) for d in Vr)
        Sh5 = math.fsum(a(d) for d in Vh)
        AddS = set(Add)
        RemS = set(Rem)
        RA = {}
        RAx = {}
        for d in Add:
            pf, mu, ph = info[d]
            if mu != 1 or d == 1:
                continue
            q = pf[-1]
            yy = d // q
            assert yy in AddS and info[yy][1] == -1
            RA[yy] = RA.get(yy, 0.0) + 1.0 / (q - 1)
            if info[yy][0][-1] >= 73:
                RAx[yy] = RAx.get(yy, Fraction(0)) + Fraction(1, q - 1)
        GOOD = BAD = HEAD = TAIL = 0.0
        bad_y = []
        for yy in Add:
            if info[yy][1] != -1:
                continue
            w = abs(a(yy))
            R = RA.get(yy, 0.0)
            GOOD += w * max(1 - R, 0)
            if info[yy][0][-1] >= 73 and RAx.get(yy, 0) > 1:
                BAD += w * float(RAx[yy] - 1)
                bad_y.append(yy)
            elif R > 1:
                BAD += w * (R - 1)
                bad_y.append(yy)
            if info[yy][0][-1] <= 37:
                HEAD += w
            if info[yy][0][-1] >= 73:
                TAIL += w
        # Rem
        RR = {}
        for d in Rem:
            pf, mu, ph = info[d]
            if mu != -1:
                continue
            q = pf[-1]
            yy = d // q
            assert yy in RemS and info[yy][1] == 1
            RR[yy] = RR.get(yy, 0.0) + 1.0 / (q - 1)
        GOODR = BADR = 0.0
        for yy in Rem:
            if info[yy][1] != 1:
                continue
            R = RR.get(yy, 0.0)
            GOODR += a(yy) * max(1 - R, 0)
            BADR += a(yy) * max(R - 1, 0)
        # flujo sobre el soporte completo V^(5) (para S(5) <= 1)
        RV = {}
        VrS = set(Vr)
        for d in Vr:
            pf, mu, ph = info[d]
            if mu != 1 or d == 1:
                continue
            q = pf[-1]
            yy = d // q
            RV[yy] = RV.get(yy, 0.0) + 1.0 / (q - 1)
        maxRV = max(RV.values())
        res[L] = {"okL": okL, "n": len(U), "EA": EA, "ER": ER, "GOOD": GOOD,
                  "BAD": BAD, "HEAD": HEAD, "TAIL": TAIL,
                  "bad": sorted(bad_y), "maxRA": max(RA.values()),
                  "argRA": max(RA, key=RA.get), "GOODR": GOODR, "BADR": BADR,
                  "slack": float(1 - M37) * HEAD / BAD if BAD else float("inf"),
                  "S5": S5, "Sh5": Sh5, "maxRV": maxRV,
                  "nRV1": sum(1 for v in RV.values() if v > 1),
                  "idA": abs((GOOD - BAD) + EA), "RA251": RAx.get(251)}
        R = res[L]
        print("    L=%d (s=%d en L_A: %s): |union|=%d, E_Add=%.6f, BAD=%.4e, "
              "GOOD=%.6f, GOOD/BAD=%.2f, HEAD37=%.4f, cola>=73=%.4f, holgura "
              "FC-2=%.2f, malos=%d (%s..%s), max R_A=%.5f en %d; E_Rem=%.6f "
              "BAD_R=%.1e GOOD_R=%.2e; S(5)=%.6f Shat(5)=%.6f max R_V=%.5f "
              "(#R_V>1: %d)"
              % (L, s, okL, R["n"], EA, BAD, GOOD, GOOD / BAD if BAD else 0,
                 HEAD, TAIL, R["slack"], len(R["bad"]), R["bad"][:1], R["bad"][-1:],
                 R["maxRA"], R["argRA"], ER, BADR, GOODR, S5, Sh5, maxRV, R["nRV1"]))
    R = res[79051268]
    exp_bad = [p for p in primes_upto(331) if 251 <= p <= 317]
    check("117a.K2", "nivel decisivo L=79051268: E_Add=-0,029314, BAD=1,086e-3, "
          "GOOD=0,030399 (margen x28), malos = los 13 primos 251..317, "
          "max R_A=1,04522 en 317, HEAD37=0,0321, cola=0,0462, holgura 4,13; "
          "E_Rem=-0,000556, BAD_R=0",
          R["okL"] and abs(R["EA"] + 0.029314) < 5e-7 and abs(R["BAD"] - 1.086e-3) < 5e-7
          and abs(R["GOOD"] - 0.030399) < 5e-7 and R["bad"] == exp_bad
          and abs(R["maxRA"] - 1.04522) < 5e-6 and R["argRA"] == 317
          and abs(R["HEAD"] - 0.0321) < 5e-5 and abs(R["TAIL"] - 0.0462) < 5e-5
          and abs(R["slack"] - 4.13) < 5e-3 and abs(R["ER"] + 0.000556) < 5e-7
          and R["BADR"] == 0 and R["idA"] < 1e-12,
          "margen GOOD/BAD = %.2f; R_A(251) exacto = %s"
          % (R["GOOD"] / R["BAD"], R["RA251"]))
    R1, R3 = res[46948396], res[113296524]
    check("117a.K2b", "niveles de control: E_Add -0,036949 / -0,022764, BAD "
          "1,086e-3 / 7,30e-4, holguras 4,98 / 4,89",
          abs(R1["EA"] + 0.036949) < 5e-7 and abs(R3["EA"] + 0.022764) < 5e-7
          and abs(R1["BAD"] - 1.086e-3) < 5e-7 and abs(R3["BAD"] - 7.30e-4) < 5e-7
          and abs(R1["slack"] - 4.98) < 5e-3 and abs(R3["slack"] - 4.89) < 5e-3,
          "holguras %.2f / %.2f" % (R1["slack"], R3["slack"]))
    rhs_fc2 = float(1 - M37) * R["HEAD"]
    check("117a.S", "fuerza lógica: el lado derecho de (FC-2) es MENOR que el "
          "de (FLOW-COMP) (FC-2 exige más: es la condición suficiente, más fuerte)",
          rhs_fc2 < R["GOOD"], "(1-M37)*HEAD37 = %.5f  <  GOOD = %.5f (%.0f %%)"
          % (rhs_fc2, R["GOOD"], 100 * rhs_fc2 / R["GOOD"]))
    check("114.S1", "S(5) <= 1 sobrevive aquí solo por compensación: hay y con "
          "R_V(y) > 1 en el soporte completo V^(5) (la ruta puntual de S<=1 "
          "también muere)", R["maxRV"] > 1 and R["S5"] <= 1,
          "S(5) = %.6f, max R_V = %.5f, #y con R_V>1 = %d" % (R["S5"], R["maxRV"], R["nRV1"]))
    return res


# ------------------------------------------- 116b J2 / J3 y construcciones
def counterexamples():
    head("Doc 116b J2/J3: contraejemplos con predicados propios")
    z = 331
    D = z ** 5
    yA = int(z ** (10.0 / 3.0))
    P = [p for p in primes_upto(z - 1) if p > 3]
    rows = [(73, 2042687, 1945083, 0.76346), (101, 771269, 5151504, 0.82415),
            (151, 230807, 17214364, 0.90856), (211, 84629, 46948396, 0.97144),
            (251, 50261, 79051268, 1.00220), (283, 35069, 113296524, 1.02855),
            (307, 27479, 144590262, 1.03552)]
    ok = True
    txt = []
    for (Y, s, L, RAexp) in rows:
        okrow = (D // s == L and is_prime(s) and z <= s < yA)
        okrow &= in_V((5,), L) and in_V((Y, 5), L) and not in_Vhat((Y,), 5, L)
        RA = Fraction(0)
        nq = 0
        for q in P:
            if 5 < q < Y and in_V(merged((Y, q), (5,)), L) and not in_Vhat((Y, q), 5, L):
                RA += Fraction(1, q - 1)
                nq += 1
        closed = sum((Fraction(1, q - 1) for q in P if 5 < q < Y), Fraction(0))
        okrow &= (RA == closed) and abs(float(RA) - RAexp) < 5e-6
        ok &= okrow
        txt.append("Y=%d R_A=%.5f (%d q)%s" % (Y, float(RA), nq, "" if okrow else " FAIL"))
        if Y == 251:
            frac_doc = Fraction(22633668515439323928857875, 22584054598420241438576064)
            ok &= (RA == frac_doc) and nq == 50
    check("116b.J2", "tabla del contraejemplo R_A en marco (z=331, r=5, y=Y "
          "primo): niveles en L_A, R_A exacto = forma cerrada; R_A(251) = la "
          "fracción impresa, 50 primos", ok, "; ".join(txt))
    first = None
    for Y in P:
        if sum((Fraction(1, q - 1) for q in P if 5 < q < Y), Fraction(0)) > 1:
            first = Y
            break
    check("116b.J2b", "primer cruce de Mertens(Y-1) - 1/4 > 1 en Y = 251",
          first == 251, "primer Y = %s" % first)
    # J3: la instancia R_R del doc
    z3 = 6052
    D3 = z3 ** 5
    L3 = D3 // 5
    uL = math.log(L3) / math.log(z3)
    inLA = (5 >= z3)
    rho, Y1, Y2 = 4999, 5009, 4099
    y = Y1 * Y2
    v = (L3 == 1623770295721736806 and in_Vhat((Y1, Y2), rho, L3)
         and not in_V(merged((Y1, Y2), (rho,)), L3))
    RR = Fraction(0)
    n = 0
    for q in primes_upto(Y2 - 1):
        if q < 5:
            continue
        if in_Vhat((Y1, Y2, q), rho, L3) and not in_V(merged((Y1, Y2, q), (rho,)), L3):
            RR += Fraction(1, q - 1)
            n += 1
    check("116b.J3a", "la instancia R_R del doc es correcta como aritmética "
          "(y en Rem^+, R_R = Mertens(251) = 1,25620, 52 primos)",
          v and RR == mertens(251) and n == 52,
          "R_R=%.5f, %d primos" % (float(RR), n))
    check("116b.J3b", "¿está L = floor(z^5/5) en L_A? (requiere s primo en "
          "[z, z^(10/3)) y u_L <= 4)", inLA and uL <= 4,
          "s = 5 < z = %d; u_L = %.3f > 4: el nivel NO está en L_A" % (z3, uL))
    # Instancia propia legítima de R_R > 1 en L_A
    found = None
    zz = 4001
    Dz = zz ** 5
    yz = int(zz ** (10.0 / 3.0))
    Pz = [p for p in primes_upto(zz - 1) if p > 3]
    for Y2c in [p for p in Pz if 650 <= p <= 760]:
        for rhoc in [p for p in Pz if Y2c < p <= Y2c + 150]:
            for Y1c in [p for p in Pz if rhoc < p <= rhoc + 150]:
                yy = Y1c * Y2c
                lo = max(rhoc * Y1c ** 3, rhoc * yy * 71 ** 3)
                hi = rhoc * Y1c * Y2c ** 3
                if not (lo < hi <= zz ** 4):
                    continue
                s0 = max(zz, -(-Dz // hi))
                for s in range(s0, min(yz, Dz // lo + 2)):
                    L = Dz // s
                    if L <= lo:
                        break
                    if L <= hi and is_prime(s):
                        found = (Y2c, rhoc, Y1c, s, L)
                        break
                if found:
                    break
            if found:
                break
        if found:
            break
    ok = False
    detail = "sin hallazgo"
    if found:
        Y2c, rhoc, Y1c, s, L = found
        ok = (Dz // s == L and zz <= s < yz and is_prime(s)
              and in_V((rhoc,), L)
              and in_Vhat((Y1c, Y2c), rhoc, L)
              and not in_V(merged((Y1c, Y2c), (rhoc,)), L))
        RRc = Fraction(0)
        n = 0
        for q in Pz:
            if q < Y2c and q not in (rhoc, Y1c) and \
                    in_Vhat((Y1c, Y2c, q), rhoc, L) and \
                    not in_V(merged((Y1c, Y2c, q), (rhoc,)), L):
                RRc += Fraction(1, q - 1)
                n += 1
        ok &= RRc > 1
        detail = ("z=%d, s=%d, L=%d (u_L=%.3f), rho=%d, y=%d*%d: R_R(y) = %.5f "
                  "con %d primos" % (zz, s, L, math.log(L) / math.log(zz), rhoc,
                                     Y1c, Y2c, float(RRc), n))
    check("116b.J3c", "instancia PROPIA de R_R(y) > 1 con nivel en L_A "
          "(el enunciado del doc sobrevive con otro ejemplo)", ok, detail)
    # y COMPUESTO malo (contra 'solo la configuración libre y primo')
    # condiciones: 5 Y1^3 >= L (modelo falla en m=1) y L > max(Y1^3,
    # Y1 Y2 Y3^3 (fusionado m=3), 125 y q_max (fusionado m=5, el primo 5))
    # => Y1^2 > max(Y2 Y3^3 / 5, 25 Y2 Y3 q_max)  => Y1 > 28513
    Y3, Y2b = 251, 257
    found = None
    for Y1c in primes_upto(33000):
        if Y1c * Y1c * 5 <= Y2b * Y3 ** 3 or Y1c * Y1c <= 25 * Y2b * Y3 * 241:
            continue
        zz = Y1c + 2
        Dz = zz ** 5
        yz = int(zz ** (10.0 / 3.0))
        yy = Y1c * Y2b * Y3
        lo = max(Y1c ** 3, Y1c * Y2b * Y3 ** 3, 125 * yy * 241)
        hi = 5 * Y1c ** 3
        if not lo < hi:
            continue
        s0 = max(zz, -(-Dz // hi))
        for s in range(s0, min(yz, Dz // lo + 2)):
            L = Dz // s
            if L <= lo:
                break
            if L <= hi and is_prime(s):
                found = (Y1c, zz, s, L)
                break
        if found:
            break
    ok = False
    detail = "sin hallazgo"
    if found:
        Y1c, zz, s, L = found
        Pz = [p for p in primes_upto(zz - 1) if p > 3]
        pf = (Y1c, Y2b, Y3)
        ok = (zz <= s < int(zz ** (10.0 / 3.0)) and in_V((5,), L)
              and in_V(merged(pf, (5,)), L) and not in_Vhat(pf, 5, L))
        RAc = Fraction(0)
        n = 0
        for q in Pz:
            if 5 < q < Y3 and in_V(merged(pf + (q,), (5,)), L) and \
                    not in_Vhat(pf + (q,), 5, L):
                RAc += Fraction(1, q - 1)
                n += 1
        ok &= RAc > 1
        detail = ("z=%d, s=%d, L=%d (u_L=%.3f), r=5, y=%d*%d*%d (omega 3): "
                  "R_A(y) = %.5f con %d primos" % (zz, s, L, math.log(L) / math.log(zz),
                                                   Y1c, Y2b, Y3, float(RAc), n))
    check("116b.C", "y COMPUESTO con R_A(y) > 1 en L_A (la 'configuración libre' "
          "no exige y primo)", ok, detail)


# ---------------------------------------------- Mertens a gran escala
def mertens_numerics():
    head("Sumas de Mertens a gran escala (116c J3, 118a J2/J3, 118b2 J0)")
    inv1 = lambda p: 1.0 / (p - 1)
    inv2 = lambda p: 1.0 / (p - 2)
    # 116c J3a
    Z = 1e8
    s1 = fsum_primes(Z ** 0.5, Z ** (2 / 3), inv1)
    Z = 1e12
    s2 = fsum_primes(Z ** (1 / 6), Z ** (5 / 9), inv1)
    s3 = fsum_primes(1e8 ** 0.5, 1e8 ** (2 / 3), lambda p: 0.0 if p in (5, 7) else 1.0 / ((p - 1) * 24))
    check("116c.J3a", "B1(2): 0,2866 a Z=10^8; B1(5/3): 1,1917 a Z=10^12; "
          "cofactor 35: 0,0119", abs(s1 - 0.2866) < 5e-5 and abs(s2 - 1.1917) < 5e-5
          and abs(s3 - 0.0119) < 5e-5,
          "%.4f (lím %.4f) / %.4f (lím %.4f) / %.4f (pred %.4f)"
          % (s1, math.log(4 / 3), s2, math.log(10 / 3), s3, math.log(4 / 3) / 24))
    # 116c J3c: K y K M_inf
    vals = []
    for z in (37, 53, 101, 1009, 10007):
        P = [p for p in primes_upto(z - 1) if p > 3]
        Pf = math.prod(1 - 1 / (p - 1) for p in P)
        Qc = math.prod(1 - 1 / (p - 1) + 1 / (p - 1) ** 2 for p in P)
        vals.append((z, Pf * math.log(z), Qc * math.log(z)))
    expK = [1.4668, 1.4600, 1.4690, 1.4790, 1.4809]
    expQ = [1.6931, 1.6894, 1.7035, 1.7181, 1.7205]
    ok = all(abs(v[1] - e) < 5e-5 for v, e in zip(vals, expK)) and \
        all(abs(v[2] - e) < 5e-5 for v, e in zip(vals, expQ))
    C2 = 0.6601618158468696
    K = 3 * math.exp(-0.5772156649015329) * C2 / 0.75
    check("116c.J3c", "Pi_full ln z -> K y Q_c ln z -> K M_inf (z = 37 ... 10007)",
          ok, "; ".join("z=%d: %.4f, %.4f" % v for v in vals)
          + " | K teórico = 3 e^-gamma C2/(3/4) = %.5f" % K)
    # M_inf
    X = 10 ** 7
    lm = 0.0
    for blk in primes_open(4, X):
        lm += math.fsum(math.log1p(1.0 / ((p - 1) * (p - 2))) for p in blk)
    Mx = math.exp(lm)
    hi = Mx * math.exp(1.0 / (X - 2))
    check("116c.Minf", "M_inf en [1,16181601, 1,16182763] (doc 115)",
          1.16181601 <= Mx and hi <= 1.16182763 + 1e-12,
          "M(10^7) = %.8f, cota superior con cola %.8f; K M_inf = %.5f"
          % (Mx, hi, K * Mx))
    # 118a J2: exterior e interior
    ext = []
    inn = []
    for e in (4, 6, 8, 10):
        z = 10.0 ** e
        ext.append(fsum_primes(z ** 0.55, z ** 0.6, inv2))
        inn.append(fsum_primes(4.5, z ** 0.5 + 0.5, inv2))
    check("118a.J2", "exterior 0,0850 -> 0,0870 = log(12/11); interior "
          "1,2262 -> 2,1324 (z = 10^4 ... 10^10)",
          abs(ext[0] - 0.0850) < 5e-5 and abs(ext[-1] - 0.0870) < 5e-5
          and abs(inn[0] - 1.2262) < 5e-5 and abs(inn[-1] - 2.1324) < 5e-5,
          "ext %s | int %s | log(12/11) = %.4f"
          % (["%.4f" % x for x in ext], ["%.4f" % x for x in inn], math.log(12 / 11)))
    # 118a J3: rebanada omega = 1 con peso 1/(p-2)
    r = []
    for e in (10, 12):
        z = 10.0 ** e
        r.append((fsum_primes(z ** (1 / 6), z ** (5 / 9), inv2),
                  fsum_primes(z ** 0.5, z ** (2 / 3), inv2)))
    check("118a.J3", "rebanada omega=1 a z=10^12: 1,1935 (B1=1,2040) y 0,2876 "
          "(log 4/3=0,2877)", abs(r[1][0] - 1.1935) < 5e-5 and abs(r[1][1] - 0.2876) < 5e-5,
          "z=10^10: %.4f, %.4f | z=10^12: %.4f, %.4f" % (r[0] + r[1]))
    check("118a.C", "constante de la rebanada omega=1: log(10/3)/M_inf = 1,036",
          abs(math.log(10 / 3) / Mx - 1.036) < 5e-4,
          "%.4f" % (math.log(10 / 3) / Mx))
    # 118b2 J0
    a6 = fsum_primes(1e6 ** 0.55, 1e6 ** 0.6, inv2)
    a8 = fsum_primes(1e8 ** 0.55, 1e8 ** 0.6, inv2)
    fr = []
    for e in (4, 6, 8):
        z = 10.0 ** e
        qs = [q for blk in primes_open(4.5, z ** 0.6) for q in blk]
        pref = [0.0]
        for q in qs:
            pref.append(pref[-1] + 1.0 / (q - 2))
        num = den = 0.0
        for p in qs:
            if not (z ** 0.55 < p < z ** 0.6):
                continue
            k_all = bisect_left(qs, p)
            k_res = bisect_right(qs, z / p)
            num += pref[k_res] / (p - 2)
            den += pref[k_all] / (p - 2)
        fr.append(num / den)
    # variantes de definición para localizar la discrepancia de J0(b)
    var = []
    for e in (4, 6, 8):
        z = 10.0 ** e
        qs = [q for blk in primes_open(4.5, z ** 0.6) for q in blk]
        pref = [0.0]
        for q in qs:
            pref.append(pref[-1] + 1.0 / (q - 2))
        ps = [p for p in qs if z ** 0.55 < p < z ** 0.6]
        den = sum(pref[bisect_left(qs, p)] / (p - 2) for p in ps)
        v1 = sum(pref[bisect_right(qs, z ** 0.4)] / (p - 2) for p in ps) / den
        v2 = sum(pref[bisect_right(qs, z ** 0.45)] / (p - 2) for p in ps) / den
        var.append((v1, v2))
    print("    J0(b) variantes: corte fijo x_q <= 0,40: %s; x_q <= 0,45: %s"
          % (["%.3f" % v[0] for v in var], ["%.3f" % v[1] for v in var]))
    alt1 = fsum_primes(1e6 ** 0.55, 1e6 ** 0.65, inv2)
    alt2 = fsum_primes(1e6 ** 0.45, 1e6 ** 0.60, inv2)
    print("    J0(c) los valores impresos 0,166 y 0,282 son las masas de los "
          "intervalos (0,55; 0,65) -> %.4f y (0,45; 0,60) -> %.4f a z=10^6, "
          "no las de la subcaja impresa" % (alt1, alt2))
    box1 = fsum_primes(1e6 ** 0.56, 1e6 ** 0.65, inv2)
    box2 = fsum_primes(1e6 ** 0.50, 1e6 ** 0.55, inv2)
    check("118b2.J0ab", "exterior 0,0869/0,0868 (z=10^6/10^8); fracción con "
          "u_r>=1: 0,718/0,781/0,808",
          abs(a6 - 0.0869) < 5e-5 and abs(a8 - 0.0868) < 5e-5
          and all(abs(f - e) < 5e-4 for f, e in zip(fr, (0.718, 0.781, 0.808))),
          "%.4f / %.4f | %s" % (a6, a8, ["%.3f" % f for f in fr]))
    check("118b2.J0c", "masas de la caja omega=3 a z=10^6: 0,166 y 0,282 "
          "(impresas)", abs(box1 - 0.166) < 5e-4 and abs(box2 - 0.282) < 5e-4,
          "medido %.4f (lím log(0,65/0,56) = %.4f) y %.4f (lím log(0,55/0,50) = %.4f)"
          % (box1, math.log(0.65 / 0.56), box2, math.log(0.55 / 0.5)))
    return Mx


# ------------------------------------------------- regiones (exacto)
def regions():
    head("Regiones de Mertens: verificación exacta de cajas (116c §2.4, 118b2 §1.3)")
    F = Fraction

    def box_ok(box, ustar, u0, vis_only=True):
        """box: lista de (a_i, b_i) decrecientes; condiciones en todo punto."""
        for i in range(len(box) - 1):
            if not box[i + 1][1] <= box[i][0]:
                return False
        if not box[0][1] <= 1:
            return False
        if not sum(a for a, b in box) >= ustar - u0:
            return False
        if not sum(b for a, b in box) <= ustar:
            return False
        for j in range(len(box)):
            if j % 2 == 0:
                if not sum(b for a, b in box[:j + 1]) + 2 * box[j][1] <= ustar:
                    return False
        return all(a < b for a, b in box)

    boxes = {F(5, 2): [(F(70, 100), F(80, 100)), (F(35, 100), F(45, 100))],
             F(3): [(F(85, 100), F(95, 100)), (F(65, 100), F(75, 100))],
             F(7, 2): [(F(90, 100), F(95, 100)), (F(80, 100), F(85, 100)),
                       (F(50, 100), F(55, 100))],
             F(4): [(F(95, 100), F(97, 100)), (F(88, 100), F(90, 100)),
                    (F(68, 100), F(69, 100))]}
    txt = []
    ok = True
    for u, bx in boxes.items():
        good = box_ok(bx, u, F(3, 2))
        mass = math.prod(math.log(float(b) / float(a)) for a, b in bx)
        ok &= good and mass > 0
        txt.append("u*=%s: %s, masa %.2e" % (u, "válida" if good else "NO", mass))
    check("116c.R", "cajas propias de medida positiva para u* = 5/2, 3, 7/2, 4 "
          "(omega 2, 2, 3, 3): la obstrucción cubre [5/3, 4]", ok, "; ".join(txt))
    # 118b2: hermana omega=3 (u0* = 1) y escalera k=2 (u0* = 1/2), u_L = 2
    u = F(2)
    sis = [(F(56, 100), F(65, 100)), (F(50, 100), F(55, 100))]
    c1 = sis[0][1] * 3 < u          # m=1
    c3 = sis[0][1] + sis[1][1] < u  # m=3 con x_q -> 0
    low = sum(a for a, b in sis)    # u_r < 2 - low
    st = []
    eta = F(1, 200)
    centers = [F(2, 3), F(2, 3), F(2, 9), F(2, 9)]
    lad = []
    for i, c in enumerate(centers):
        hi = c * (1 - eta * (1 + 2 * (i % 2)))
        lo = c * (1 - eta * (2 + 2 * (i % 2)))
        lad.append((lo, hi))
    okl = all(lad[i + 1][1] <= lad[i][0] for i in range(3))
    for j in (0, 2):
        okl &= sum(b for a, b in lad[:j + 1]) + 2 * lad[j][1] < u
    xq = F(1, 100)
    okl &= sum(b for a, b in lad) + 3 * xq < u
    ur_hi = u - sum(a for a, b in lad)
    ok = c1 and c3 and (u - low) < F(1) and okl and ur_hi < F(1, 2)
    check("118b2.R", "hermana omega=3: visible y con u_r < 0,94 (corta u0* > "
          "0,94); escalera k=2: Nathanson estricto en m=1,3,5 y u_r < 1/2",
          ok, "hermana u_r < %.3f; escalera u_r < %.4f" % (float(u - low), float(ur_hi)))


# ------------------------------------------- 118b1 J3 y la frontera q0+
def strata_frontier():
    head("Doc 118b1 J3: E(z), kappa0, q0+ y su dependencia en z")
    out = []
    for z in (17, 29, 37, 53):
        P_lt = [p for p in primes_upto(z - 1) if p > 3]
        P_le = [p for p in primes_upto(z) if p > 3]
        E_lt = math.prod(p / (p - 1) for p in P_lt)
        E_le = math.prod(p / (p - 1) for p in P_le)
        k_lt = math.prod(1 - 1 / (p - 1) ** 2 for p in P_lt)
        k_le = math.prod(1 - 1 / (p - 1) ** 2 for p in P_le)
        out.append((z, E_lt / math.log(z), E_le / math.log(z), k_lt, k_le,
                    2 * E_le / math.log(z) ** 3))
    for o in out:
        print("    z=%d: E/ln z = %.4f (p<z) / %.4f (p<=z); kappa0 = %.4f (p<z) "
              "/ %.4f (p<=z); cola 2E/(ln z)^3 = %.3f" % o)
    ok = all(abs(o[2] - e) < 5e-5 for o, e in zip(out, (0.6517, 0.6267, 0.6207, 0.6169)))
    ok &= abs(out[0][4] - 0.8926) < 5e-5 and abs(out[-1][4] - 0.8834) < 5e-5
    ok &= all(abs(o[5] - e) < 5e-4 for o, e in zip(out, (0.162, 0.111, 0.095, 0.078)))
    check("118b1.J3", "E/ln z = 0,6517/0,6267/0,6207/0,6169, kappa0 0,8926 -> "
          "0,8834, colas 0,162/0,111/0,095/0,078: se reproducen con p <= z "
          "(no con P = primos < z)", ok, "")

    # J4: cota Euler del estrato 5 y suma degenerada sin visibilidad
    eul = []
    for z in (29, 37):
        b = Base.get(z)
        eul.append(math.prod(p / (p - 1) for p in b.P) * 0.8 / b.Qc)
    tot = {}
    for (tag, z, D, y) in FRAMES[1:]:
        b = Base.get(z)
        best_pos = best_all = 0.0
        for L in frame_levels(z, D, y):
            lnz = math.log(z)
            s_pos = s_all = 0.0
            for rm in range(1, b.N):
                r = b.val[rm]
                q0 = b.q0(rm)
                if q0 is None or q0 ** 3 * r < L or b.Pi[rm] >= 0.5:
                    continue
                if math.log(L / r) / lnz >= U0:
                    continue
                s_all += 1.0 / b.phi[rm]
                if r < L:
                    s_pos += 1.0 / b.phi[rm]
            best_pos = max(best_pos, s_pos / b.Qc)
            best_all = max(best_all, s_all / b.Qc)
        tot[(tag, z)] = (best_pos, best_all)
    print("    J4: cota Euler estrato 5 / Q_c = %s; suma degenerada en banda sin "
          "visibilidad (r < L / todo el enum): %s"
          % (["%.2f" % e for e in eul],
             {k: tuple(round(x, 3) for x in v) for k, v in tot.items()}))
    okJ4 = abs(eul[0] - 3.25) < 5e-3 and abs(eul[1] - 3.72) < 5e-3
    vals = [tot[("A", 29)], tot[("A", 37)], tot[("B", 29)]]
    okJ4 &= any(all(abs(v[i] - e) < 5e-4 for v, e in zip(vals, (0.716, 1.405, 0.025)))
                for i in (0, 1))
    check("118b1.J4", "cota Euler del estrato 5 = 3,25 / 3,72 Q_c (P = primos < z); "
          "suma degenerada sin visibilidad 0,716 / 1,405 / 0,025 Q_c", okJ4, "")

    def q0dag(z):
        lz3 = math.log(z) ** 3
        prod = 1
        for q in [p for p in primes_upto(10 ** 4) if p > 3]:
            if prod >= lz3:
                return q
            prod *= q
        return None
    tab = [(z, q0dag(z)) for z in (17, 29, 1443, 1444, 10 ** 6, 10 ** 12, 10 ** 30, 10 ** 100)]
    check("118b1.F", "q0+(z) = min{q : P(q) >= (ln z)^3} es 13 solo para "
          "26,3 < z <= 1443,6 y crece con z", tab[1][1] == 13 and tab[2][1] == 13
          and tab[3][1] == 17 and tab[-1][1] > 13,
          "; ".join("z=%s: %s" % (("%.0e" % z) if z > 10 ** 5 else z, q) for z, q in tab))
    check("118b1.A", "aritmética del texto: '5^9 > 2*10^7' (118b1 §2)",
          5 ** 9 > 2 * 10 ** 7, "5^9 = %d" % 5 ** 9)


# --------------------------------------------- 119a: exponente y FL
def fl_check():
    head("Doc 119a: exponente de la condición m=1 y término principal FL")
    b = Base.get(37)
    L = 37 ** 2
    lv = Level(b, L)
    rm = 1 << b.idx[7]
    Vh = sorted(b.val[m] for m in lv.Vhat(rm))
    prim_cube = [q for q in b.P if q != 7 and 7 * q ** 3 < L]
    prim_sq = [q for q in b.P if q != 7 and 7 * q ** 2 < L]
    check("119a.E", "la condición m=1 de Vhat es r q^3 < L (no r q^2 < L): "
          "a z=37, L=37^2, r=7, Vhat = {1, 5}",
          Vh == [1, 5] and prim_cube == [5] and prim_sq != [5],
          "Vhat^(7) = %s; primos con 7q^3<L: %s; con 7q^2<L serían %s"
          % (Vh, prim_cube, prim_sq))
    eg = math.exp(0.5772156649015329)
    txt = []
    ok = True
    for z in (10 ** 3, 10 ** 4, 10 ** 5):
        for u in (1.25, 1.45):
            Dl = z ** u
            w = Dl ** (1 / 3)
            Pw = [p for p in primes_upto(int(w) + 1) if p > 3 and p ** 3 < Dl]
            Pd = sorted(Pw, reverse=True)
            # Shat = sum_{d en V+(D)} mu/phi sobre los primos de P (el
            # conductor r se toma > D^(1/3): su exclusión no cambia V+(D))
            tot = 0.0
            cnt = 0

            def rec(d, a_, lo, i, m):
                nonlocal tot, cnt
                tot += a_
                cnt += 1
                for k in range(i, len(Pd)):
                    p = Pd[k]
                    nd = d * p
                    if nd >= Dl:
                        continue
                    nlo = lo
                    if (m + 1) % 2 == 1:
                        t = p * p * nd
                        if t >= Dl:
                            continue
                        nlo = max(lo, t)
                    rec(nd, -a_ / (p - 1), nlo, k + 1, m + 1)
            sys.setrecursionlimit(10000)
            rec(1, 1.0, 0, 0, 0)
            Pz = [p for p in primes_upto(z - 1) if p > 3]
            # conductor r = menor primo de P por encima de max(D^(1/3), z^(1/2))
            r = next(p for p in Pz if p > max(w, z ** 0.5))
            Pi = math.prod(1 - 1 / (p - 1) for p in Pz if p != r)
            ratio = tot / Pi
            txt.append("z=10^%d u=%.2f: Shat/Pi = %.3f (|V+|=%d; F(u)=2e^g/u=%.3f)"
                       % (round(math.log10(z)), u, ratio, cnt, 2 * eg / u))
            ok &= ratio > 1.5
    check("119a.FL", "Shat/Pi > 1 de forma estable (defecto del modelo = "
          "término principal, hacia F(u) = 2e^gamma/u)", ok, "; ".join(txt))


# ------------------------------------------------------- aritmética suelta
def arithmetic():
    head("Aritmética de los textos")
    M67 = mertens(67)
    M71 = mertens(71)
    M37 = mertens(37)
    check("116b.J0", "Mertens(67) = 160107799/160240080 = 0,99917 < 1; "
          "Mertens(71) = 1,01346 > 1",
          M67 == Fraction(160107799, 160240080) and M71 > 1
          and abs(float(M71) - 1.01346) < 5e-6,
          "%s = %.5f; %.5f" % (M67, float(M67), float(M71)))
    check("117a.K0", "Mertens(37) = 15899/18480 = 0,860335; 1 - Mertens(37) = "
          "2581/18480 = 0,139665", M37 == Fraction(15899, 18480)
          and 1 - M37 == Fraction(2581, 18480), "%s" % M37)
    # (R-FLOW) es [T] para min(y) <= 71: q < 71 => q <= 67
    check("114.R", "R_A(y) <= Mertens(q < min(y)) < 1 también para min(y) = 71 "
          "(q <= 67): el umbral real es min(y) >= 73", M67 < 1, "")
    phiy = (5009 - 1) * (4099 - 1)
    term = (mertens(251) - 1) / phiy
    check("117a.K4", "término malo Rem del J3 = (Mertens(251)-1)/phi(20531891) "
          "= 1,248e-8", abs(float(term) - 1.248e-8) < 5e-12, "%.4e" % float(term))
    Pi17 = math.prod(1 - 1 / (p - 1) for p in (5, 7, 11, 13))
    Pi19 = Pi17 * (1 - 1 / 16)
    check("116c.J1", "Pi_full(17) = 0,515625 > 1/2 > Pi_full(19) = 0,4834",
          abs(Pi17 - 0.515625) < 1e-12 and abs(Pi19 - 0.4834) < 5e-5, "%.6f / %.4f" % (Pi17, Pi19))
    check("116c.B1", "B1(5/3) = log(10/3) = 1,2040; B1(2) = log(4/3) = 0,2877",
          abs(math.log(10 / 3) - 1.2040) < 5e-5 and abs(math.log(4 / 3) - 0.2877) < 5e-5, "")
    check("119b.A", "aritmética del texto: '5*7*11*13 = 15015' (119b §1)",
          5 * 7 * 11 * 13 == 15015, "5*7*11*13 = %d (15015 = 3*5*7*11*13)" % (5 * 7 * 11 * 13))
    # N(r) acotado: número de cofactores admisibles por estrato
    for q0 in (5, 7, 11):
        Pq = [p for p in primes_upto(q0 ** 3) if p >= q0]
        cnt = 0

        def rec(prod, i):
            nonlocal cnt
            cnt += 1
            for k in range(i, len(Pq)):
                if prod * Pq[k] >= q0 ** 3:
                    break
                rec(prod * Pq[k], k + 1)
        rec(1, 0)
        print("    estrato q0=%d: cofactores c < q0^3 libres de cuadrados con "
              "primos >= q0: %d (cota absoluta de N(r))" % (q0, cnt))


def bad_prime_decay():
    head("Doc 117a §1b: masa mala de la familia prima libre (r0 = 5) sin PNT")
    rows = []
    for e in (12, 15, 18, 21, 24):
        L = 10 ** e
        top = round(L ** (1.0 / 3.0))
        while top ** 3 >= L:
            top -= 1
        lo = (L / 5) ** (1.0 / 3.0)
        M = 0.0          # sum_{5 <= p < Y} 1/(p-1)
        bad = 0.0
        for blk in primes_open(4, top + 1):
            for Y in blk:
                if Y >= lo:
                    RA = M - 0.25          # q en (5, Y): todos admisibles
                    if RA > 1:
                        bad += (RA - 1) / (Y - 1)
                M += 1.0 / (Y - 1)
        ll = math.log(math.log(L)) / math.log(L)
        # predicción de Mertens elemental: (ln ln Y - 1,7153) ln(ln b/ln a),
        # con sum_{5<=p<=x} 1/(p-1) = ln ln x - 0,46535 + o(1)
        ymid = (lo * top) ** 0.5
        pred = (math.log(math.log(ymid)) - 1.71535) * math.log(math.log(top) / math.log(lo))
        rows.append((e, bad, pred, ll, bad / ll))
    for r_ in rows:
        print("    L=10^%d: masa mala máxima de la familia prima = %.4e "
              "(predicción de Mertens %.4e); (log log L)/log L = %.4e; "
              "cociente %.3f" % r_)
    tmax = math.exp(1 + 1.71535)
    print("    la predicción (ln t - 1,7153)/t, t = ln L^(1/3), tiene su máximo "
          "en t = e^2,715 = %.1f (L ~ 10^%.0f) y decae después como "
          "log log L / log L" % (tmax, 3 * tmax / math.log(10)))
    ok = all(abs(r_[1] - r_[2]) / r_[1] < 0.01 for r_ in rows)
    check("117a.PNT", "la masa mala de la familia prima sigue la predicción de "
          "Mertens ELEMENTAL al 1 %: masa = ln5 (ln ln L - 2,81)/ln L (1+o(1)) "
          "-> 0, sin PNT", ok, "cocientes masa/(loglogL/logL) %s (-> ln 5 = 1,609)"
          % ["%.3f" % r_[4] for r_ in rows])


def main():
    sys.stdout.reconfigure(encoding="utf-8")
    print("Verificación docs 116-120a (Claude Code). Python %s" % sys.version.split()[0])
    arithmetic()
    st, band, per_level = main_sweep()
    report_main(st, band, per_level)
    band_reports(band)
    stress_114()
    recs, st53 = layer_reports(band)
    profile_119a(recs)
    nondeg_25()
    rep_levels()
    double_count()
    slab_reports()
    strata_frontier()
    fl_check()
    regions()
    z331()
    counterexamples()
    bad_prime_decay()
    mertens_numerics()
    head("RESUMEN")
    npass = sum(1 for c in CHECKS if c[2])
    print("  Cada chequeo contrasta una afirmación del texto (o una de este "
          "informe). FAIL = la afirmación del texto NO se sostiene tal como "
          "está escrita: es un hallazgo, interpretado en el informe.")
    for c in CHECKS:
        if not c[2]:
            print("  FAIL %s: %s | %s" % (c[0], c[1], c[3]))
    print("  %d/%d afirmaciones confirmadas; %d hallazgos | %.1f s"
          % (npass, len(CHECKS), len(CHECKS) - npass, time.time() - T0))


if __name__ == "__main__":
    main()
