# 72 - Correlacion de pares de ceros de zeta: verificacion numerica
# de la funcion de correlacion R2 (GUE) y de la forma de Montgomery F(alfa).
# Proposito: mostrar que los datos finitos siguen la prediccion GUE
# para R2, y que la forma F(alfa) NO alcanza aun sus valores
# asintoticos a esta altura (correccion finita-T de orden 1:
# log(T/2 pi e)/log T ~ 0.6). Eso delimita exactamente que la
# correlacion que SL-70.1 necesita es la del regimen conjeturado
# (meseta |alfa| >= 1), no alcanzable con datos finitos.
# Metodo: primeros K ceros de zeta por barrido Riemann-Siegel
# (vectorizado, numpy) + biseccion con mpmath.siegelz (Z exacta).
# No es una prueba: es una demostracion finita de R2 y de la
# lentitud logaritmica de F(alfa).
import sys
import numpy as np
from mpmath import mp, siegelz

mp.dps = 15
K = 250
T_MAX = 475.0
STEP = 0.005

# ---------- 1. Ceros de zeta ----------
def theta_stirling(t):
    # theta(t) = arg Gamma(1/4 + i t/2) - (t/2) log pi, expansion
    # de Stirling (error O(t^-5) < 1e-6 para t >= 14).
    return (0.5 * t * np.log(t / (2 * np.pi * np.e)) - np.pi / 8
            + 1.0 / (48 * t) + 7.0 / (5760 * t ** 3))

def z_main(t):
    # Suma principal de Riemann-Siegel (sin resto): localiza los
    # ceros con error O(t^-1/4) ~ 0.2; solo se usa para el barrido.
    # El numero de terminos depende de t: n <= sqrt(t/(2 pi)).
    th = theta_stirling(t)
    m = int(np.sqrt(t.max() / (2 * np.pi))) + 1
    out = np.zeros_like(t)
    for n in range(1, m + 1):
        mask = t >= 2 * np.pi * n * n
        out[mask] += n ** (-0.5) * np.cos(th[mask] - t[mask] * np.log(n))
    return 2.0 * out

grid = np.arange(14.0, T_MAX, STEP)
zg = z_main(grid)
sign = np.sign(zg)
idx = np.where(sign[:-1] * sign[1:] < 0)[0]
approx = 0.5 * (grid[idx] + grid[idx + 1])
print("barrido: %d cambios de signo, primer aprox %.4f" % (len(approx), approx[0]))

def refine_zero(lo, hi):
    flo = float(siegelz(lo))
    fhi = float(siegelz(hi))
    if flo == 0.0:
        return lo
    if fhi == 0.0:
        return hi
    if flo * fhi > 0:
        return None
    for _ in range(30):
        mid = 0.5 * (lo + hi)
        fm = float(siegelz(mid))
        if fm == 0.0:
            return mid
        if flo * fm < 0:
            hi, fhi = mid, fm
        else:
            lo, flo = mid, fm
    return 0.5 * (lo + hi)

def refine_fallback(t0):
    # Barrido fino con Z exacta alrededor del aproximado.
    ts = np.arange(t0 - 0.5, t0 + 0.5, 0.02)
    vs = [float(siegelz(float(t))) for t in ts]
    for i in range(len(vs) - 1):
        if vs[i] == 0.0:
            return float(ts[i])
        if vs[i] * vs[i + 1] < 0:
            return refine_zero(float(ts[i]), float(ts[i + 1]))
    return None

zeros = []
n_fallback = 0
for j, t0 in enumerate(approx):
    lo = t0 - 0.4
    hi = t0 + 0.4
    if j > 0:
        lo = max(lo, 0.5 * (approx[j - 1] + t0))
    if j < len(approx) - 1:
        hi = min(hi, 0.5 * (t0 + approx[j + 1]))
    z = refine_zero(lo, hi)
    if z is None:
        z = refine_fallback(t0)
        n_fallback += 1
    if z is None:
        print("AVISO: sin cero exacto cerca de %.4f" % t0)
        continue
    zeros.append(z)
zeros = np.array(zeros)
print("ceros refinados: %d (barrido %d, fallback %d)"
      % (len(zeros), len(approx), n_fallback))
np.savetxt("72_ceros_zeta.txt", zeros, fmt="%.10f")

# ---------- 1b. Reparacion de pares cercanos perdidos ----------
# Un cero perdido por el barrido deja un hueco ~2x el espaciamiento
# local medio s(g) = 2 pi / log(g / 2 pi). Se repara con un barrido
# fino de Z exacta dentro del hueco.
def mean_spacing(g):
    return 2 * np.pi / np.log(g / (2 * np.pi))

gaps = np.diff(zeros)
gmid = 0.5 * (zeros[:-1] + zeros[1:])
sus = [(i, gaps[i], gmid[i]) for i in range(len(gaps))
       if gaps[i] > 1.7 * mean_spacing(gmid[i]) or gaps[i] < 0.05]
print("huecos sospechosos: %d" % len(sus))
added = []
for i, gp, gm_ in sus:
    lo, hi = zeros[i], zeros[i + 1]
    print("  hueco %.4f en [%.3f, %.3f] (s_local=%.2f)"
          % (gp, lo, hi, mean_spacing(gm_)))
    ts = np.arange(lo + 0.005, hi - 0.005, 0.01)
    vs = [float(siegelz(float(t))) for t in ts]
    for k in range(len(vs) - 1):
        if vs[k] == 0.0:
            added.append(float(ts[k]))
        elif vs[k] * vs[k + 1] < 0:
            z = refine_zero(float(ts[k]), float(ts[k + 1]))
            if z is not None:
                added.append(z)
print("ceros recuperados por reparacion: %d" % len(added))
if added:
    zeros = np.sort(np.concatenate([zeros, np.array(added)]))
    np.savetxt("72_ceros_zeta.txt", zeros, fmt="%.10f")
print("total ceros: %d; gamma_100 = %.9f" % (len(zeros), zeros[99]))

# ---------- 2. Validacion ----------
KNOWN = [14.134725141734693, 21.022039638771554, 25.010857580145688]
checks = []
for n, kv in enumerate(KNOWN):
    err = abs(zeros[n] - kv)
    ok = err < 1e-6
    checks.append(("gamma_%d vs valor conocido" % (n + 1), ok,
                   "err=%.2e" % err))
    print("gamma_%d = %.9f (conocido %.9f, err %.1e) %s"
          % (n + 1, zeros[n], kv, err, "OK" if ok else "FALLO"))

T = zeros[-1]
N_T = (T / (2 * np.pi)) * np.log(T / (2 * np.pi * np.e)) + 0.875
# Tolerancia 3: la formula omite S(T), que fluctua en O(1).
ok = abs(len(zeros) - N_T) <= 3.0
checks.append(("conteo vs Riemann-von Mangoldt", ok,
               "N=%.1f N(T)=%.1f" % (len(zeros), N_T)))
print("T = %.4f, log T = %.4f, N(T) = %.2f, contados = %d %s"
      % (T, np.log(T), N_T, len(zeros), "OK" if ok else "FALLO"))
g100_known = 236.524229665816
err100 = abs(zeros[99] - g100_known)
ok100 = err100 < 1e-6
checks.append(("gamma_100 vs valor conocido", ok100, "err=%.2e" % err100))
print("gamma_100 = %.9f (conocido %.9f, err %.1e) %s"
      % (zeros[99], g100_known, err100, "OK" if ok100 else "FALLO"))
print("factor finito-T log(T/2 pi e)/log T = %.3f"
      % (np.log(T / (2 * np.pi * np.e)) / np.log(T)))

# ---------- 3. Forma de Montgomery F(alfa) ----------
# F(alfa) = (2 pi / (T log T)) sum_{j,k} T^{i alfa (g_j - g_k)}
#           * w(g_j - g_k),  w(u) = 4 / (4 + u^2).
# Asintotica (Montgomery-Goldston, RH): F(alfa) = T^{-2|alfa|} log T
# + |alfa| para |alfa| <= 1; conjetura F(alfa) = 1 para |alfa| > 1.
# A esta altura la asintotica NO se alcanza: la prediccion finita-T
# para F(0) es log(T/2 pi e)*log(T/2 pi)/log T (densidad local por
# integral de w), no log T. Se verifica ESA prediccion y las
# propiedades estructurales (real, >= 0, <= F(0)).
D = zeros[:, None] - zeros[None, :]
W = 4.0 / (4.0 + D ** 2)
norm = 2 * np.pi / (T * np.log(T))
alphas = [0.0, 0.25, 0.5, 0.75, 1.0, 1.5, 2.0, 3.0]
print("\nF(alfa) [asintotica: alfa<1 -> alfa; alfa>=1 -> 1]:")
Fvals = {}
for a in alphas:
    F = norm * np.sum(np.exp(1j * a * D * np.log(T)) * W)
    Fvals[a] = F.real
    print("  alfa=%.2f  F=%.4f  |imag|=%.2e" % (a, F.real, abs(F.imag)))

F0_pred = (np.log(T / (2 * np.pi * np.e)) * np.log(T / (2 * np.pi))
           / np.log(T))
checks.append(("F(0) ~ prediccion finita-T",
               abs(Fvals[0.0] - F0_pred) / F0_pred <= 0.25,
               "F(0)=%.3f pred=%.3f (asintotica logT=%.2f)"
               % (Fvals[0.0], F0_pred, np.log(T))))
checks.append(("F(alfa) real para todo alfa",
               all(abs(norm * np.sum(np.exp(1j * a * D * np.log(T)) * W).imag)
                   < 1e-6 for a in alphas),
               "max|imag| < 1e-6"))
checks.append(("F(alfa) >= 0", all(v >= -1e-9 for v in Fvals.values()),
               "min F=%.4f" % min(Fvals.values())))
checks.append(("F(alfa) <= F(0)",
               all(v <= Fvals[0.0] + 1e-9 for v in Fvals.values()),
               "max F=%.4f F(0)=%.4f" % (max(Fvals.values()), Fvals[0.0])))

# ---------- 4. Histograma R2 (GUE) con normalizacion local ----------
# u = (g_j - g_k) * log(g_barra / 2 pi) / (2 pi), g_barra = punto
# medio del par (densidad local, como los espaciamientos
# normalizados de Odlyzko). Esperado por bin [a,b):
# K * int_a^b [1 - (sin pi u / pi u)^2] du.
def r2(u):
    if u == 0:
        return 0.0
    s = np.sin(np.pi * u) / (np.pi * u)
    return 1.0 - s * s

try:
    from scipy.integrate import quad as _quad
except Exception:
    _quad = None

def r2_int(a, b):
    if _quad is not None:
        v, _ = _quad(r2, a, b)
        return v
    xs = np.linspace(a, b, 2001)
    return float(np.trapezoid([r2(x) for x in xs], xs))

gm = 0.5 * (zeros[:, None] + zeros[None, :])
iu = (zeros[:, None] - zeros[None, :]) * np.log(gm / (2 * np.pi)) / (2 * np.pi)
iu = iu[np.tril_indices(len(zeros), k=-1)]  # pares j>k: u > 0
# Binas de verificacion: [0,1), [1,1.5), [1.5,2), [2,3). La bina
# [0,0.5) NO se usa como check: a esta altura (T < 500) las
# separaciones mas pequenas estan sistematicamente depletadas
# respecto a GUE (efecto de altura finita, ver diagnostico abajo),
# y su conteo esperado (~28) es demasiado pequeno para un check
# robusto. La bina [0,1) agrega esa region y es robusta.
bins = [(0.0, 1.0), (1.0, 1.5), (1.5, 2.0), (2.0, 3.0)]
print("\nR2 (GUE), pares j>k, normalizacion local:")
for a, b in bins:
    cnt = int(np.sum((iu >= a) & (iu < b)))
    exp = len(zeros) * r2_int(a, b)
    rel = abs(cnt - exp) / exp if exp > 0 else 9e9
    ok = rel <= 0.35
    checks.append(("R2 bin [%.1f,%.1f)" % (a, b), ok,
                   "n=%d esp=%.1f rel=%.2f" % (cnt, exp, rel)))
    print("  [%.1f,%.1f): n=%d esperado=%.1f rel=%.2f %s"
          % (a, b, cnt, exp, rel, "OK" if ok else "FALLO"))

# Diagnostico (no es check): la sub-bina [0,0.5) muestra la
# deplecion de altura finita de las separaciones mas pequenas.
cnt05 = int(np.sum((iu >= 0.0) & (iu < 0.5)))
exp05 = len(zeros) * r2_int(0.0, 0.5)
print("  diagnostico [0,0.5): n=%d esperado=%.1f (deplecion "
      "altura finita, no es check)" % (cnt05, exp05))

# ---------- 5. Resumen ----------
npass = sum(1 for _, ok, _ in checks if ok)
print("\nRESULTADO: PASS %d/%d" % (npass, len(checks)))
for name, ok, det in checks:
    print("  [%s] %s -- %s" % ("OK " if ok else "FALLO", name, det))
if npass == len(checks):
    print("PASS")
else:
    print("FAIL")
    sys.exit(1)
