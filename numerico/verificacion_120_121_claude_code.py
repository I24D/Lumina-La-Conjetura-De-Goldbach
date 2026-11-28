#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Verificación de las prioridades de la cola [P-rev] 116–121: el peso de D-NP del doc 121/121b
(W_d = Ω_d^C·φ(d)) y los veredictos (A2), (A6) de la auditoría DIST-D (doc 120b), con la
recuperación de (A6) del doc 121a.

Procedencia: Claude Code (Anthropic, modelo claude-opus-5-5), 2026-11-28.
Método: reimplementación PROPIA desde las definiciones de los docs 101, 106, 115 §7, 120b §1 y
121b (caracteres por construcción propia vía CRT y raíces primitivas; Λ por criba; V+ de
Nathanson con nivel floor(D/q)). No importa ni ejecuta código de nadie.
Estado: PROPUESTA (pendiente de contraste por Muse). [T] para el álgebra, [M] para lo numérico.
No es revisión independiente según el contrato (README, regla 4).
Solo biblioteca estándar. Uso:  python -B verificacion_120_121_claude_code.py
Nada de esto demuestra Goldbach, que sigue abierta.
"""
import cmath
import math
import sys
import time
from itertools import compress, product

T0 = time.time()
CHECKS = []
fsum = math.fsum


def check(name, ok, detail=""):
    CHECKS.append((name, bool(ok)))
    print(f"[{'PASS' if ok else 'FAIL'}] {name}" + (f" | {detail}" if detail else ""), flush=True)


def sieve(n):
    s = bytearray([1]) * (n + 1)
    s[0] = s[1] = 0
    for p in range(2, math.isqrt(n) + 1):
        if s[p]:
            s[p * p::p] = bytearray(len(range(p * p, n + 1, p)))
    return s


def factor(n):
    f = {}
    p = 2
    while p * p <= n:
        while n % p == 0:
            f[p] = f.get(p, 0) + 1
            n //= p
        p += 1
    if n > 1:
        f[n] = f.get(n, 0) + 1
    return f


def phi(n):
    out = n
    for p in factor(n):
        out = out // p * (p - 1)
    return out


LAMMAX = 1000
_S = sieve(LAMMAX)
LAM = [0.0] * (LAMMAX + 1)
for p in compress(range(LAMMAX + 1), _S):
    x = p
    while x <= LAMMAX:
        LAM[x] = math.log(p)
        x *= p


class Chars:
    """Todos los caracteres de (Z/nZ)*, n impar, por CRT con raíces primitivas (construcción propia)."""

    def __init__(self, n):
        self.n = n
        comps = []
        for p, a in sorted(factor(n).items()):
            m = p ** a
            ph = m // p * (p - 1)
            qs = list(factor(ph))
            g = next(g for g in range(2, m) if math.gcd(g, m) == 1 and all(pow(g, ph // q, m) != 1 for q in qs))
            ind = {}
            x = 1
            for t in range(ph):
                ind[x] = t
                x = x * g % m
            comps.append((m, ph, ind))
        self.comps = comps
        self.labels = list(product(*[range(c[1]) for c in comps])) if comps else [()]

    def val(self, k, x):
        x %= self.n
        if math.gcd(x, self.n) != 1:
            return 0j
        ang = 0.0
        for (m, ph, ind), ki in zip(self.comps, k):
            ang += ki * ind[x % m] / ph
        return cmath.exp(2j * math.pi * ang)


def S_mq(N, m, q, chi_q, kq):
    return sum(LAM[n] * chi_q.val(kq, n) for n in range(2, N + 1) if (n - N) % m == 0 and n % q)


def Psi(N, f):
    return sum(LAM[n] * f(n) for n in range(2, N + 1))


def part_121b():
    print("\n=== Doc 121b: reducción CRT de D-NP ===", flush=True)
    cases = [(200, 3, 7), (200, 21, 13), (300, 7, 11), (100, 11, 7), (300, 13, 17)]
    ok_conj, err_noconj = True, {}
    lam = {}
    for N, m, q in cases:
        Cm, Cq = Chars(m), Chars(q)
        e_c = e_n = 0.0
        lhs_tot = rhs_tot = 0.0
        for kq in Cq.labels:
            if all(x == 0 for x in kq):
                continue
            S = S_mq(N, m, q, Cq, kq)
            Ps = [Psi(N, lambda n, kp=kp: Cm.val(kp, n) * Cq.val(kq, n)) for kp in Cm.labels]
            f_c = sum(Cm.val(kp, N).conjugate() * P for kp, P in zip(Cm.labels, Ps)) / phi(m)
            f_n = sum(Cm.val(kp, N) * P for kp, P in zip(Cm.labels, Ps)) / phi(m)
            e_c = max(e_c, abs(S - f_c))
            e_n = max(e_n, abs(S - f_n))
            lhs_tot += abs(S) ** 2
            rhs_tot += sum(abs(P) ** 2 for P in Ps) / phi(m)
        ok_conj &= e_c < 1e-9
        err_noconj[(N, m, q)] = e_n
        lam[(N, m, q)] = rhs_tot / lhs_tot
        print(f"  (N,m,q)=({N},{m},{q}): |S − fórmula con ψ̄(N)| máx = {e_c:.1e}; sin conjugar = {e_n:.2e}; "
              f"CS: Σ|S|² ≤ (1/φ(m))ΣΣ|Ψ|² {'✓' if lhs_tot <= rhs_tot * (1 + 1e-12) else '✗'}; λ = {lam[(N, m, q)]:.2f}", flush=True)
    check("W1 Ortogonalidad S_{m,q}(χ) = (1/φ(m))Σ_ψ ψ̄(N)Ψ_N(ψχ) exacta en los 5 casos; la variante sin conjugar falla en (200,21,13)",
          ok_conj and err_noconj[(200, 21, 13)] > 1.0, f"sin conjugar en (200,21,13): {err_noconj[(200, 21, 13)]:.2f}")
    check("W2 Pérdida de Cauchy–Schwarz λ = 2,40 (m=3), 4,06 (m=7), 11,64 (m=21), 14,52 (m=13)",
          abs(lam[(200, 3, 7)] - 2.40) < 6e-3 and abs(lam[(300, 7, 11)] - 4.06) < 6e-3 and abs(lam[(200, 21, 13)] - 11.64) < 6e-3
          and abs(lam[(300, 13, 17)] - 14.52) < 6e-3, f"m=11: {lam[(100, 11, 7)]:.2f}")
    # Conteo CRT [T]: #{η mod qm : η_q ≠ χ0} = φ(m)(q−2)
    okc = all(len(Chars(m).labels) * (len(Chars(q).labels) - 1) == phi(m) * (q - 2) for _, m, q in cases)
    check("W3 Conteo CRT: φ(m)(q−2) pares (ψ,χ≠χ0) = #{η mod qm : η_q ≠ χ0}", okc)
    # (m,N)>1: (200, 35, 11)
    N, m, q = 200, 35, 11
    Cq = Chars(q)
    kq = next(k for k in Cq.labels if any(k))
    S = S_mq(N, m, q, Cq, kq)
    pred = math.log(5) * Cq.val(kq, 3)
    check("W4 (m,N) = 5 en (200,35,11): S = log5·χ(3) ≠ 0 mientras ψ̄(200) = 0 para todo ψ mod 35 (la fórmula daría 0)",
          abs(S - pred) < 1e-12 and abs(S) > 1, f"S = {S:.4f}")


def vplus(small_desc, L):
    out = []

    def rec(i0, prod, ell, fs):
        out.append((prod, fs))
        for i in range(i0, len(small_desc)):
            p = small_desc[i]
            if p >= L:
                continue
            np_ = prod * p
            if (ell + 1) % 2 == 1 and p * p * np_ >= L:
                continue
            rec(i + 1, np_, ell + 1, fs + (p,))

    if L > 1:
        rec(0, 1, 0, ())
    return out


def family(z):
    y = round(z ** (10 / 3))
    D = z ** 5
    s = sieve(y + 1)
    pr = [i for i in compress(range(len(s)), s)]
    small = sorted((p for p in pr if p < z and p not in (2, 3)), reverse=True)
    sw = [p for p in pr if z <= p < y and p not in (2, 3)]
    H = {}
    branch = {}
    for q in sw:
        V = vplus(small, D // q)
        branch[q] = V
        for d, fs in V:
            fac = (q,) + fs
            for k in range(1 << len(fac)):
                hs = tuple(sorted(fac[i] for i in range(len(fac)) if k >> i & 1))
                H[math.prod(hs)] = hs
    return H, sw, small, branch


def omega_C(H):
    Om = {}
    for r, hs in H.items():
        if len(hs) < 2:
            continue
        val = ((1 << len(hs)) - 1) / math.prod(p - 1 for p in hs)
        for k in range(1, 1 << len(hs)):
            sub = tuple(hs[i] for i in range(len(hs)) if k >> i & 1)
            d = math.prod(sub)
            Om[d] = Om.get(d, 0.0) + val
    return Om


def theta(branch_q, m):
    """Θ̃_m(q) = Σ_{h∈V+(D/q)∖{1}, m|h} (2τ(h)−1)/φ(h)."""
    return fsum(((1 << (len(fs) + 1)) - 1) / math.prod(p - 1 for p in fs) for h, fs in branch_q if h > 1 and h % m == 0)


def part_weights():
    print("\n=== Doc 115 §7 y doc 121b §3: identidad puntual y el peso W_d = Ω_d^C·φ(d) ===", flush=True)
    worst = 0.0
    fams = {}
    for z in (11, 13, 17):
        H, sw, small, branch = family(z)
        Om = omega_C(H)
        fams[z] = (H, sw, small, branch, Om)
        for q in sw:
            smooth_m = {m for m, _ in branch[q]}
            for m in smooth_m:
                th = theta(branch[q], m)
                d = q * m
                om = Om.get(d, 0.0)
                worst = max(worst, abs(th / (q - 1) - om) / max(om, 1e-300) if om else abs(th))
    check("P1 Identidad puntual Ω_{qm}^C = Θ̃_m(q)/(q−1) para todo (q, m) en S11, S13, S17 (doc 115 §7)", worst < 1e-12,
          f"máx err rel {worst:.1e}")
    # Término a término en el marco z = 17, N = 200 (doc 121b K3)
    H, sw, small, branch, Om = fams[17]
    N = 200
    ratios = []
    holds = True
    for q in (17, 19, 23):
        Cq = Chars(q)
        for m in (7, 11, 13):
            th = theta(branch[q], m)
            if th <= 0:
                continue
            Cm = Chars(m)
            d = q * m
            sumS = 0.0
            sumPsi = 0.0
            for kq in Cq.labels:
                if not any(kq):
                    continue
                sumS += abs(S_mq(N, m, q, Cq, kq)) ** 2
                for kp in Cm.labels:
                    sumPsi += abs(Psi(N, lambda n, kp=kp, kq=kq: Cm.val(kp, n) * Cq.val(kq, n))) ** 2
            if sumS <= 1e-12:
                continue
            term = th * phi(m) ** 2 * sumS
            corr = Om[d] * phi(d) * sumPsi
            wrong = Om[d] / phi(d) * sumPsi
            lam_ = (sumPsi / phi(m)) / sumS
            holds &= term <= corr * (1 + 1e-12)
            ratios.append(term / wrong)
            print(f"  (q,m)=({q},{m}): término real / forma Ω/φ = {term/wrong:.0f}; φ(d)²/λ = {phi(d)**2/lam_:.0f}; "
                  f"término ≤ forma corregida Ω·φ·ΣΨ²: {'✓' if term <= corr*(1+1e-12) else '✗'}", flush=True)
    check("P2 El término real de D_NP cumple la cota con peso Ω_d^C·φ(d) en todas las parejas, y supera la forma con Ω_d^C/φ(d) "
          "por 1622×–4768× (= φ(d)²/λ)", holds and len(ratios) == 9 and abs(min(ratios) - 1622) < 1 and abs(max(ratios) - 4768) < 1,
          f"{len(ratios)} parejas, rango [{min(ratios):.0f}, {max(ratios):.0f}]")
    return fams


def part_120b(fams):
    print("\n=== Doc 120b: (A1), (A2), (A7) y (A6) ===", flush=True)
    cases = [(100, 3, 7), (100, 4, 7), (200, 5, 7), (200, 6, 11), (500, 3, 7), (500, 9, 13), (500, 4, 11), (300, 7, 13)]
    ok1 = ok2c = ok2w = ok7 = True
    for N, m, q in cases:
        psiN = sum(LAM[n] for n in range(2, N + 1))
        A = sum(LAM[n] for n in range(2, N + 1) if (n - N) % m == 0)
        L = sum(LAM[n] for n in range(2, N + 1) if (n - N) % m == 0 and n % q == 0)
        B = {a: sum(LAM[n] for n in range(2, N + 1) if (n - N) % m == 0 and n % q == a) for a in range(q)}
        ok1 &= abs(sum(B[a] for a in range(1, q)) - (A - L)) < 1e-9
        Em = A - psiN / phi(m)
        Eqm = sum(LAM[n] for n in range(2, N + 1) if (n - N) % (q * m) == 0) - psiN / phi(q * m)
        Delta = B[N % q] - (A - L) / (q - 1)
        corr = Delta + (Em - L) / (q - 1)
        wrong = Delta + Em - L / (q - 1)
        ok2c &= abs(corr - Eqm) < 1e-9
        ok2w &= abs((wrong - Eqm) - (q - 2) / (q - 1) * Em) < 1e-9
        # (A7) Parseval: Σ_χ |S(χ)|² = (q−1) Σ_{a≠0} B(a)²
        Cq = Chars(q)
        par = sum(abs(sum(B[a] * Cq.val(k, a) for a in range(1, q))) ** 2 for k in Cq.labels)
        ok7 &= abs(par - (q - 1) * sum(B[a] ** 2 for a in range(1, q))) < 1e-7 * max(1, par)
        print(f"  (N,m,q)=({N},{m},{q}): E(N;qm,N) = {Eqm:.6f}; (A2) corregida = {corr:.6f}; (A2) escrita − real = {wrong-Eqm:.6f} "
              f"= ((q−2)/(q−1))E_m = {(q-2)/(q-1)*Em:.6f}", flush=True)
    check("A1 Σ_{a∈(ℤ/q)*} B_{m,q}(a) = A_m − L_{m,q} en los 8 casos", ok1)
    check("A2 (A2) escrita (Δ + E_m − L/(q−1)) es falsa y su error es exactamente ((q−2)/(q−1))·E_m; "
          "la corregida E(N;qm,N) = Δ_{q,m} + (E_m − L_{m,q})/(q−1) es exacta [T]", ok2c and ok2w)
    check("A7 Parseval Σ_χ|S(χ)|² = (q−1)Σ_{a≠0}B(a)² en los 8 casos", ok7)
    # (A6): α̃_m = Σ_q Ω_{qm}^C ≤ Ω_m^C (recuperada en 121a por unicidad de rama); α_m de la versión auditada lleva (q−1)
    for z in (11, 13, 17):
        H, sw, small, branch, Om = fams[z]
        uniq = all(sum(1 for p in hs if p >= z) <= 1 for hs in H.values())
        viol_t = viol_w = 0
        nm = 0
        ratio_w = 0.0
        smooth = {d for d, hs in H.items() if d > 1 and hs[-1] < z}
        for m in smooth:
            if Om.get(m, 0.0) == 0.0:
                continue
            nm += 1
            at = fsum(Om.get(q * m, 0.0) for q in sw)
            aw = fsum(Om.get(q * m, 0.0) * (q - 1) ** 2 for q in sw)   # α_m := Σ_q Θ̃_m(q)(q−1) = Σ_q (q−1)² Ω_{qm}^C
            viol_t += at > Om[m] * (1 + 1e-12)
            viol_w += aw > Om[m]
            ratio_w = max(ratio_w, aw / Om[m])
        print(f"  S{z}: unicidad de rama (≤ 1 primo ≥ z por r) = {uniq}; m lisos con Ω_m^C>0: {nm}; "
              f"violaciones de α̃_m ≤ Ω_m^C: {viol_t}; violaciones con el α_m auditado: {viol_w} (cociente máx {ratio_w:.3g})", flush=True)
        if z == 11:
            ok6 = uniq and viol_t == 0 and viol_w > 0
        else:
            ok6 &= uniq and viol_t == 0 and viol_w > 0
    check("A6 Con la unicidad de rama [T], α̃_m = Σ_q Ω_{qm}^C ≤ Ω_m^C sin violaciones (S11–S17); "
          "el α_m auditado (factor (q−1) invertido) viola la cota", ok6)
    # Contraejemplo abstracto del doc 120b: ℛ = {3·5·7}, m = 3
    tau, ph = 8, 2 * 4 * 6
    check("A6b Contraejemplo de 120b (ℛ = {105}, m = 3): Σ_{q∈{5,7}} Ω_{qm} = 0,2917 > Ω_m = 0,1458 — necesita dos primos de "
          "conmutación en un mismo r, imposible en ℋ", abs(2 * (tau - 1) / ph - 0.2917) < 6e-5 and abs((tau - 1) / ph - 0.1458) < 6e-5)
    # Doc 121a (d): en S11, m = 5, Σ_q Ω_{qm}^C = 1,091744; con peso τ(r)⁻¹/φ(r) el intercambio da 0,072937
    H, sw, small, branch, Om = fams[11]
    at = fsum(Om.get(q * 5, 0.0) for q in sw)
    alt = fsum(1 / ((1 << len(hs)) * math.prod(p - 1 for p in hs)) for r, hs in H.items()
               if len(hs) >= 2 and r % 5 == 0 and hs[-1] >= 11)
    print(f"  S11, m = 5: Σ_q Ω_(5q)^C = {at:.6f} (doc 121a: 1,091744); intercambio con τ(r)⁻¹/φ(r) = {alt:.6f} (doc 121a: 0,072937)", flush=True)
    check("A6c Doc 121a (d): Σ_q Ω_{5q}^C = 1,091744 en S11 y el intercambio con τ⁻¹/φ da 0,072937 (factor ~15)",
          abs(at - 1.091744) < 6e-7 and abs(alt - 0.072937) < 6e-7, f"factor {at/alt:.1f}")


def main():
    try:
        sys.stdout.reconfigure(encoding="utf-8")
    except (AttributeError, ValueError):
        pass
    print(f"Python {sys.version.split()[0]}", flush=True)
    part_121b()
    fams = part_weights()
    part_120b(fams)
    npass = sum(ok for _, ok in CHECKS)
    print(f"\nRESUMEN: {npass}/{len(CHECKS)} PASS, {len(CHECKS) - npass} FAIL. Tiempo {time.time() - T0:.1f} s", flush=True)
    return 0 if npass == len(CHECKS) else 1


if __name__ == "__main__":
    sys.exit(main())
