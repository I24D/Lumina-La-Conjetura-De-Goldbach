# 120b_auditoria_distd.py — Verificacion numerica de la auditoria DIST-D (doc 120b)
# Audita (A1), (A2) [tal como esta + corregida], (A3), (A4) termino a termino,
# (A6) [contraejemplo de conteo] y (A7), con primos reales, N <= 500.
# Solo stdlib. PASS si cada chequeo se comporta como afirma el doc 120b.
import math
from fractions import Fraction

def sieve_lam(N):
    lam = [0.0]*(N+1)
    is_p = [True]*(N+1)
    is_p[0] = is_p[1] = False
    for i in range(2, int(N**0.5)+1):
        if is_p[i]:
            for j in range(i*i, N+1, i):
                is_p[j] = False
    for p in range(2, N+1):
        if is_p[p]:
            pk = p
            while pk <= N:
                lam[pk] = math.log(p)
                if pk > N//p: break
                pk *= p
    return lam

def phi(n):
    r, x, p = n, n, 2
    while p*p <= x:
        if x % p == 0:
            while x % p == 0: x //= p
            r -= r//p
        p += 1
    if x > 1: r -= r//x
    return r

def tau(n):
    c, x, p = 1, n, 2
    while p*p <= x:
        e = 0
        while x % p == 0: x //= p; e += 1
        c *= (e+1); p += 1
    if x > 1: c *= 2
    return c

def psi_class(lam, N, d, a):
    # psi(N; d, a) = sum_{n<=N, n = a mod d} Lambda(n); a taken mod d
    return sum(lam[n] for n in range(1, N+1) if n % d == a % d)

def dlog_table(q):
    # generador de (Z/qZ)^* y tabla de logaritmos discretos
    units = [a for a in range(1, q) ]
    for g in units:
        seen, x = set(), 1
        for _ in range(q-1):
            seen.add(x); x = x*g % q
        if len(seen) == q-1:
            gen = g; break
    dl = {}
    x = 1
    for j in range(q-1):
        dl[x] = j; x = x*gen % q
    return gen, dl

def chi(k, a, q, dl):
    # caracter chi_k(a) = exp(2 pi i k dlog(a)/(q-1)); a unidad mod q
    return complex(math.cos(2*math.pi*k*dl[a % q]/(q-1)),
                   math.sin(2*math.pi*k*dl[a % q]/(q-1)))

fails = []
def check(name, cond, detail=""):
    print(("OK  " if cond else "FAIL") + f"  {name}  {detail}")
    if not cond: fails.append(name)

CASES = [(100, 3, 7), (100, 4, 7), (200, 5, 7), (200, 6, 11),
         (500, 3, 7), (500, 9, 13), (500, 4, 11), (300, 7, 13)]

for (N, m, q) in CASES:
    assert math.gcd(m, q) == 1 and N % q != 0
    lam = sieve_lam(N)
    psiN = sum(lam)
    A_m = psi_class(lam, N, m, N)
    # L_{m,q}: potencias q^j <= N con q^j = N (mod m)
    L = 0.0; pk = q
    while pk <= N:
        if pk % m == N % m: L += math.log(q)
        if pk > N//q: break
        pk *= q
    # B_{m,q}(a): masa de la clase a mod q dentro de n = N (mod m)
    B = [0.0]*q
    for n in range(1, N+1):
        if n % m == N % m: B[n % q] += lam[n]
    sumB = sum(B[a] for a in range(1, q))
    # (A1)
    check(f"(A1) N={N},m={m},q={q}", abs(sumB - (A_m - L)) < 1e-9,
          f"sumB={sumB:.6f} A_m-L={A_m-L:.6f}")
    # E de la serie y variantes
    E_qm = psi_class(lam, N, q*m, N) - psiN/phi(q*m)
    E_m_series = A_m - psiN/phi(m)
    E_m_hers   = A_m - psiN*phi(m)
    D118 = B[N % q] - (A_m - L)/(q-1)
    D117 = psi_class(lam, N, q*m, N) - A_m - L/(q-1)
    # (A2) tal como esta (118), con E_m de la serie y con la suya
    r1 = D118 + E_m_series - L/(q-1)
    r2 = D118 + E_m_hers   - L/(q-1)
    check(f"(A2)-serie falla N={N},m={m},q={q}", abs(r1 - E_qm) > 1e-6,
          f"disc={r1-E_qm:+.6f} (predicha (q-2)/(q-1)*E_m={(q-2)/(q-1)*E_m_series:+.6f})")
    check(f"(A2)-suya falla N={N},m={m},q={q}", abs(r2 - E_qm) > 1e-6,
          f"disc={r2-E_qm:+.6f}")
    # (A2) corregida: E = D118 + (E_m^serie - L)/(q-1)
    rc = D118 + (E_m_series - L)/(q-1)
    check(f"(A2)-corregida N={N},m={m},q={q}", abs(rc - E_qm) < 1e-9,
          f"disc={rc-E_qm:+.2e}")
    # las dos Delta difieren
    check(f"Delta117 != Delta118 N={N},m={m},q={q}", abs(D117 - D118) > 1e-6,
          f"difieren en {D117-D118:+.6f}")
    # (A3) y (A7) con caracteres mod q
    gen, dl = dlog_table(q)
    S = []
    for k in range(q-1):
        s = 0j
        for n in range(1, N+1):
            if n % m == N % m and n % q != 0 and lam[n] != 0.0:
                s += lam[n]*chi(k, n % q, q, dl)
        S.append(s)
    # (A3): Delta = 1/(q-1) sum_{k!=0} conj(chi_k(N)) S_k
    rhs = sum(chi(k, N % q, q, dl).conjugate()*S[k] for k in range(1, q-1))/(q-1)
    check(f"(A3) N={N},m={m},q={q}", abs(rhs - D118) < 1e-8,
          f"disc={abs(rhs-D118):.2e}")
    sumS2 = sum(abs(S[k])**2 for k in range(1, q-1))
    check(f"(A3)-CS N={N},m={m},q={q}", (q-1)*abs(D118)**2 <= sumS2 + 1e-9,
          f"(q-1)|D|^2={(q-1)*abs(D118)**2:.6f} <= sum|S|^2={sumS2:.6f}")
    # (A7) Parseval
    parseval = (q-1)*sum(B[a]**2 for a in range(1, q)) - (A_m - L)**2
    check(f"(A7) N={N},m={m},q={q}", abs(sumS2 - parseval) < 1e-7,
          f"sum|S|^2={sumS2:.6f} Parseval={parseval:.6f}")
    check(f"(A7) S_0 = A_m - L N={N},m={m},q={q}",
          abs(S[0] - (A_m - L)) < 1e-9, f"S0={S[0].real:.6f}")
    # (A4) termino a termino con peso simbolico w = Theta~ = 0.7
    w = 0.7
    F_m = phi(m)*E_m_series
    lhs = (q-1)*phi(m)**2*w*E_qm**2
    np_corr = phi(m)**2*w*sumS2
    mean_corr = w*(F_m - phi(m)*L)**2/(q-1)
    mean_hers = w*(q-1)*(F_m - phi(m)*L)**2
    check(f"(A4)-reduc-corregida N={N},m={m},q={q}",
          lhs <= 2*np_corr + 2*mean_corr + 1e-9,
          f"lhs={lhs:.4f} <= {2*np_corr + 2*mean_corr:.4f}")
    ratio = mean_hers/mean_corr if mean_corr > 0 else float('inf')
    check(f"(A4)-D_MEAN suya = (q-1)^2 x corregida N={N},m={m},q={q}",
          abs(ratio - (q-1)**2) < 1e-6, f"ratio={ratio:.1f} (q-1)^2={(q-1)**2}")

# (A6) contraejemplo de conteo para Omega segun (106.2):
# Omega_d = sum_{r en R, d|r} (tau(r)-1)/phi(r) sobre un juguete R = {r0}
m0, q1, q2 = 3, 5, 7
r0 = m0*q1*q2
w0 = Fraction(tau(r0)-1, phi(r0))
Om = lambda d: w0 if r0 % d == 0 else Fraction(0)
sumq = Om(m0*q1) + Om(m0*q2)
check("(A6) sum_q Omega_{qm} > Omega_m (conteo 106.2)", sumq > Om(m0),
      f"sum_q={float(sumq):.4f} > Omega_m={float(Om(m0)):.4f}")
# y el factor: doc 115 da Omega_{qm} = Theta~/(q-1); su alpha usa Theta~(q-1) = (q-1)^2 Omega
th = Fraction(13, 10); qq = 7
Om_from_th = th/(qq-1)
check("(A6) Theta~(q-1) != Omega_{qm} (doc 115)", th*(qq-1) != Om_from_th,
      f"Theta~(q-1)={float(th*(qq-1)):.3f} vs Omega={float(Om_from_th):.3f}; "
      f"lo correcto es Theta~/(q-1)")

print()
print("RESULTADO:", "PASS" if not fails else f"FAIL en {fails}")
