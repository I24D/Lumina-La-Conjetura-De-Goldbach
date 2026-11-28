#!/usr/bin/env python3
# 121b_dnp_crt.py — verificacion adversarial de la reduccion CRT de D_NP
# propuesta por Lumina Chat (respuesta al doc 120).
# K1 ortogonalidad mod m; K2 Cauchy-Schwarz + conteo CRT + perdida;
# K3 pesos: forma final Omega/phi(d) vs la que la derivacion da;
# K4 colapso de la formula si (m, N) > 1.
# Solo stdlib (+ maquinaria de la serie 109/110 para V+ y Theta~).
import cmath
import importlib.util
import math
import sys

spec = importlib.util.spec_from_file_location(
    "d110", "110_segundo_momento_defectos.py")
d110 = importlib.util.module_from_spec(spec)
spec.loader.exec_module(d110)
Level110 = d110.Level110

NMAX = 300


def sieve_lam(nmax):
    lam = [0.0] * (nmax + 1)
    is_p = [True] * (nmax + 1)
    is_p[0] = is_p[1] = False
    for i in range(2, nmax + 1):
        if is_p[i]:
            for j in range(2 * i, nmax + 1, i):
                is_p[j] = False
            pk = i
            while pk <= nmax:
                lam[pk] = math.log(i)
                if pk > nmax // i:
                    break
                pk *= i
    return lam


LAM = sieve_lam(NMAX)


def prim_root(p):
    # menor raiz primitiva mod p (p primo)
    n = p - 1
    facs = set()
    t = n
    d = 2
    while d * d <= t:
        while t % d == 0:
            facs.add(d)
            t //= d
        d += 1
    if t > 1:
        facs.add(t)
    g = 2
    while any(pow(g, n // f, p) == 1 for f in facs):
        g += 1
    return g


def chars_prime(p):
    """Lista de caracteres mod p: cada uno es dict a->valor (a coprimo),
    con 0 si p | a. Indice 0 = principal."""
    g = prim_root(p)
    log = {}
    v = 1
    for e in range(p - 1):
        log[v] = e
        v = v * g % p
    out = []
    for k in range(p - 1):
        w = cmath.exp(2j * math.pi * k / (p - 1))
        ch = {}
        for a in range(1, p):
            ch[a] = w ** log[a]
        out.append(ch)
    return out


def chars_mod(m, facs):
    """Caracteres mod m libre de cuadrados, via CRT de sus primos.
    facs = lista de primos distintos. Devuelve lista de funciones."""
    per = [chars_prime(p) for p in facs]
    idx_ranges = [range(len(c)) for c in per]
    combos = [()]
    for r in idx_ranges:
        combos = [c + (i,) for c in combos for i in r]
    chars = []
    for combo in combos:
        def make(combo=combo):
            def f(n):
                val = 1.0 + 0j
                for p, chs, k in zip(facs, per, combo):
                    a = n % p
                    if a == 0:
                        return 0.0 + 0j
                    val *= chs[k][a]
                return val
            return f
        chars.append(make())
    return chars


def factorize_small(m):
    fs = []
    t = m
    d = 2
    while d * d <= t:
        if t % d == 0:
            fs.append(d)
            while t % d == 0:
                t //= d
        d += 1
    if t > 1:
        fs.append(t)
    return fs


def S_direct(N, m, q, chi):
    """S_{m,q}(chi) = sum_{n<=N, n=N (mod m), q∤n} Lambda(n) chi(n)."""
    tot = 0.0 + 0j
    for n in range(2, N + 1):
        if LAM[n] == 0.0 or n % m != N % m or n % q == 0:
            continue
        tot += LAM[n] * chi[n % q]
    return tot


def psiN_of_pair(N, psi, chi, q):
    """Psi_N(psi*chi) = sum_{n<=N} Lambda(n) psi(n) chi(n)."""
    tot = 0.0 + 0j
    for n in range(2, N + 1):
        if LAM[n] == 0.0 or n % q == 0:
            continue
        v = psi(n)
        if v == 0:
            continue
        tot += LAM[n] * v * chi[n % q]
    return tot


def check_K1():
    print("K1 ortogonalidad: S = (1/phi(m)) sum_psi psibar(N) Psi_N(psi chi)")
    ok = True
    cases = [(200, 3, 7), (200, 21, 13), (300, 7, 11), (100, 11, 7),
             (300, 13, 17)]
    for (N, m, q) in cases:
        assert math.gcd(m, N) == 1 and math.gcd(m, q) == 1
        psis = chars_mod(m, factorize_small(m))
        chis = chars_prime(q)
        phim = len(psis)
        err = err_wrong = 0.0
        for k in range(1, q - 1):
            chi = chis[k]
            rhs = sum(psis[j](N).conjugate() * psiN_of_pair(N, psis[j], chi, q)
                      for j in range(phim)) / phim
            rhs_w = sum(psis[j](N) * psiN_of_pair(N, psis[j], chi, q)
                        for j in range(phim)) / phim
            lhs = S_direct(N, m, q, chi)
            err = max(err, abs(lhs - rhs))
            err_wrong = max(err_wrong, abs(lhs - rhs_w))
        print("  N=%d m=%d q=%d: err=%.2e | con psi(N) sin conjugar: %.2e"
              % (N, m, q, err, err_wrong))
        ok &= err < 1e-8
    print("  => %s" % ("PASS" if ok else "FAIL"))
    return ok


def check_K2():
    print("K2 Cauchy-Schwarz + conteo CRT + factor de perdida lambda:")
    ok = True
    for (N, m, q) in [(200, 3, 7), (200, 21, 13), (300, 7, 11),
                      (300, 13, 17)]:
        psis = chars_mod(m, factorize_small(m))
        chis = chars_prime(q)
        phim = len(psis)
        lhs = 0.0
        dbl = 0.0
        per_ok = True
        for k in range(1, q - 1):
            chi = chis[k]
            s = S_direct(N, m, q, chi)
            lhs += abs(s) ** 2
            inner = sum(abs(psiN_of_pair(N, psis[j], chi, q)) ** 2
                        for j in range(phim))
            dbl += inner
            per_ok &= abs(s) ** 2 <= inner / phim + 1e-9
        rhs = dbl / phim
        lam = rhs / lhs if lhs > 0 else float("inf")
        cnt_pairs = phim * (q - 2)
        cnt_eta = phim * (q - 1) - phim  # phi(qm) - #{eta_q = chi0}
        print("  N=%d m=%d q=%d: sum|S|^2=%.4f <= (1/phi)doble=%.4f | "
              "lambda=%.3f | pares=%d eta=%d | CS por chi: %s"
              % (N, m, q, lhs, rhs, lam, cnt_pairs, cnt_eta,
                 "OK" if per_ok else "FALLA"))
        ok &= lhs <= rhs + 1e-9 and per_ok and cnt_pairs == cnt_eta
    print("  => %s" % ("PASS" if ok else "FAIL"))
    return ok


def theta_tilde(lv, m):
    """Theta~_m(q) del doc 115 §7 sobre el nivel lv = V+(D_A/q)."""
    tot = 0.0
    for (h, pf, mu, phi) in lv.V:
        if h != 1 and h % m == 0:
            tot += (2 * (2 ** len(pf)) - 1) / phi
    return tot


def check_K3():
    print("K3 pesos: termino D_NP vs forma final suya (Omega/phi(d)) "
          "y la derivada (Omega*phi(d)):")
    ok = True
    z, A = 17, 3
    D_A = z ** 5
    N = 200
    for q in (17, 19, 23):
        lv = Level110(z, A, (), D_A // q)
        chis = chars_prime(q)
        for m in (7, 11, 13, 77):
            if math.gcd(m, N) != 1 or math.gcd(m, q) != 1:
                continue
            th = theta_tilde(lv, m)
            if th <= 0:
                continue
            psis = chars_mod(m, factorize_small(m))
            phim = len(psis)
            sumS = 0.0
            dbl = 0.0
            for k in range(1, q - 1):
                chi = chis[k]
                sumS += abs(S_direct(N, m, q, chi)) ** 2
                dbl += sum(abs(psiN_of_pair(N, psis[j], chi, q)) ** 2
                           for j in range(phim))
            Om = th / (q - 1)              # Omega_{qm}^C, doc 115 §7 [T]
            phid = (q - 1) * phim          # phi(qm)
            term = th * phim ** 2 * sumS   # termino de D_NP
            hers = Om / phid * dbl         # su forma final
            corr = Om * phid * dbl         # lo que la derivacion da
            lam = (dbl / phim) / sumS if sumS > 0 else float("inf")
            print("  q=%d m=%d: Theta~=%.4f term=%.4f | suya=%.4f "
                  "(term/suya=%.1f) | derivada=%.4f (term<=deriv: %s) | "
                  "lambda=%.2f phi(d)^2=%d%s"
                  % (q, m, th, term, hers, term / hers if hers > 0 else 0.0,
                     corr, term <= corr + 1e-9, lam, phid ** 2,
                     "  [vacio: S=0 en esta ventana]" if sumS == 0 else ""))
            ok &= term <= corr + 1e-9
            if sumS > 0:
                # fila discriminante: su forma exige lambda >= phi(d)^2
                ok &= term > hers and lam < phid ** 2
    print("  => %s" % ("PASS" if ok else "FAIL"))
    return ok


def check_K4():
    print("K4 la hipotesis (m, N) = 1 es necesaria: (N, m, q) = "
          "(200, 35, 11), (m, N) = 5:")
    N, m, q = 200, 35, 11
    psis = chars_mod(m, factorize_small(m))
    chis = chars_prime(q)
    phim = len(psis)
    maxlhs = 0.0
    maxrhs = 0.0
    for k in range(1, q - 1):
        chi = chis[k]
        maxlhs = max(maxlhs, abs(S_direct(N, m, q, chi)))
        rhs = sum(psis[j](N).conjugate() * psiN_of_pair(N, psis[j], chi, q)
                  for j in range(phim)) / phim
        maxrhs = max(maxrhs, abs(rhs))
    print("  max_chi |S| = %.4f (log 5 = %.4f) | max_chi |formula| = %.2e "
          "(psi(200) = 0 para todo psi mod 35)" % (maxlhs, math.log(5),
                                                   maxrhs))
    ok = maxlhs > 1.0 and maxrhs < 1e-12
    print("  => %s" % ("PASS" if ok else "FAIL"))
    return ok


def main():
    t0 = __import__("time").time()
    print("== DOC 121b: reduccion CRT de D_NP ==")
    ok = True
    ok &= check_K1()
    print()
    ok &= check_K2()
    print()
    ok &= check_K3()
    print()
    ok &= check_K4()
    print()
    print("Tiempo total: %.1f s" % (__import__("time").time() - t0))
    print("PASS" if ok else "FAIL")
    return 0 if ok else 1


if __name__ == "__main__":
    sys.exit(main())
