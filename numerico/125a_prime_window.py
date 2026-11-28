# Doc 125a — (PRIME-WINDOW): lema de ventana [T], asintotica de
# Mertens y medicion contra Q_c en los marcos extendidos del 123b.
# K0  reproduce los cinco valores K4 del doc 123b (M_near exacto).
# K1  lema de ventana: M_near(p0) = S(z) - S(L^(1/3)) exacto y
#     constante en p0 sobre la capa; conjuntos de primos de capa
#     (criterio propio) == maquinaria 123b; estadistica de DeltaE.
# K2  por nivel: PW, PW[1.35,1.5), SM (maquinaria), Q_c, cocientes
#     PW/Q_c, PW/SM y prediccion de Mertens log(3/u_L)*log[(u_L/3)/
#     (u_L-3/2)]; chequeo cruzado Q_c propio vs lv.qc.
# K3  Q_c(z)*log z -> c_Q; cociente predicho (0.1167/c_Q)*log z
#     contra el maximo medido de PW/Q_c.
# Solo stdlib; importa la maquinaria de los docs 110/120a/123b.
import importlib.util
import math
import sys
import time

t0 = time.time()
spec = importlib.util.spec_from_file_location(
    "d123b", "123b_outmass.py")
d123 = importlib.util.module_from_spec(spec)
spec.loader.exec_module(d123)
frame_levels = d123.frame_levels
Level110 = d123.Level110
level_sums = d123.level_sums
layer_records = d123.layer_records
TOL = 1e-9


def primes_upto(n):
    if n < 2:
        return []
    sieve = bytearray([1]) * (n + 1)
    sieve[0:2] = b"\x00\x00"
    for i in range(2, int(n ** 0.5) + 1):
        if sieve[i]:
            sieve[i * i::i] = b"\x00" * len(sieve[i * i::i])
    return [i for i, v in enumerate(sieve) if v]


def allowed(z):
    """Primos permitidos del marco: impares < z, != 3 (A=3)."""
    return [p for p in primes_upto(z - 1) if p > 2 and p != 3]


def S_sum(x, primes):
    """S(x) = sum_{q permitido, q < x} 1/(q-1); `primes` = todos
    los primos impares != 3 hasta al menos x (sin tope z)."""
    return math.fsum(1.0 / (q - 1) for q in primes if q < x)


def Qc(z):
    tot = 1.0
    for p in allowed(z):
        tot *= (1.0 - 1.0 / (p - 1) + 1.0 / (p - 1) ** 2)
    return tot


def Pi_prime(p0, z):
    tot = 1.0
    for p in allowed(z):
        if p != p0:
            tot *= (1.0 - 1.0 / (p - 1))
    return tot


def Qset(p0, L, z):
    return [q for q in allowed(z)
            if q > p0 and q ** 3 >= L and p0 * q < L]


def M_near(p0, L, z):
    return math.fsum(1.0 / (q - 1) for q in Qset(p0, L, z))


def layer_primes_mine(L, z):
    lnz = math.log(z)
    out = []
    for p0 in allowed(z):
        if p0 ** 3 >= L:
            continue
        ur = math.log(L / p0) / lnz
        if not (1.0 <= ur < 1.5):
            continue
        if Pi_prime(p0, z) >= 0.5:
            continue
        # no degenerado: existe primo permitido q != p0, q^3 < L/p0
        if not any(q != p0 and q ** 3 * p0 < L for q in allowed(z)):
            continue
        out.append(p0)
    return out


BIG_PRIMES = [p for p in primes_upto(100000) if p > 2 and p != 3]


def E_of(x):
    return S_sum(x, BIG_PRIMES) - math.log(math.log(x))


# ------------------------------------------------------------------ K0
def check_K0():
    print("K0 valores K4 del doc 123b (M_near exacto):")
    casos = [(37, 1286, 7, 0.4159), (53, 2585, 7, 0.3309),
             (61, 3450, 11, 0.3674), (67, 4024, 11, 0.3840),
             (71, 4421, 11, 0.3992)]
    ok = True
    for (z, L, p0, ref) in casos:
        v = M_near(p0, L, z)
        good = abs(v - ref) < 1e-4
        ok &= good
        print("  z=%d L=%d p0=%d: M_near = %.6f (doc 123b: %.4f) %s"
              % (z, L, p0, v, ref, "OK" if good else "DISCREPA"))
    print("  => %s" % ("PASS" if ok else "FAIL"))
    return ok


# ------------------------------------------------------- marcos K1/K2
def frames_ext():
    for z in (37, 53, 61, 67, 71):
        levels = frame_levels(z, z ** 5, int(z ** (10.0 / 3.0)),
                              4, 2.0)
        levels = sorted(set(levels) | {z ** 5 // 5})
        yield z, levels


def check_K1():
    print("K1 lema de ventana sobre las instancias primas de capa:")
    ok = True
    max_id = 0.0
    max_spread = 0.0
    dEs = []
    for z, levels in frames_ext():
        lnz = math.log(z)
        qc_mine = Qc(z)
        for L in levels:
            lv = Level110(z, 3, (), L)
            ok &= abs(lv.qc - qc_mine) / qc_mine < 1e-12
            mach = sorted(r for (r, pf, phi, ur, Pi, dif, E, dl, d)
                          in layer_records(lv, L, lnz) if len(pf) == 1)
            mine = layer_primes_mine(L, z)
            if mach != mine:
                ok = False
                print("  z=%d L=%d: conjuntos difieren: maq=%s "
                      "propio=%s" % (z, L, mach, mine))
            ref = S_sum(z, BIG_PRIMES) - S_sum(L ** (1.0 / 3.0),
                                               BIG_PRIMES)
            vals = []
            for p0 in mach:
                v = M_near(p0, L, z)
                vals.append(v)
                max_id = max(max_id, abs(v - ref))
                dEs.append(v - (math.log(math.log(z))
                                - math.log(math.log(L ** (1.0 / 3.0)))))
            if vals:
                max_spread = max(max_spread, max(vals) - min(vals))
        print("  z=%d: E(z) = S(z)-loglog z = %.4f | Q_c = %.4f | "
              "Q_c*log z = %.4f"
              % (z, E_of(z), qc_mine, qc_mine * lnz))
    print("  identidad |M_near - (S(z)-S(L^1/3))| max = %.2e | "
          "dispersion en p0 max = %.2e" % (max_id, max_spread))
    print("  DeltaE = M_near - log(3/u_L): min = %.4f, max = %.4f"
          % (min(dEs), max(dEs)))
    ok &= max_id < 1e-9 and max_spread < 1e-9
    print("  => %s" % ("PASS" if ok else "FAIL"))
    return ok


def sweep_K2():
    print("K2 PW por nivel contra Q_c y SM (marcos extendidos):")
    C123 = {37: 0.1956, 53: 0.3152, 61: 0.1955, 67: 0.2089,
            71: 0.2598}
    ok = True
    resumen = {}
    for z, levels in frames_ext():
        lnz = math.log(z)
        best = -1.0
        bestL = None
        bestFrac = None
        poolPW = poolSM = 0.0
        for L in levels:
            lv = Level110(z, 3, (), L)
            s = level_sums(lv, L, lnz)
            qc = lv.qc
            uL = math.log(L) / lnz
            pw = pw135 = 0.0
            for (r, pf, phi, ur, Pi, dif, E, dl, d) in \
                    layer_records(lv, L, lnz):
                if len(pf) != 1:
                    continue
                term = M_near(r, L, z) / (r - 1)
                pw += term
                if 1.35 <= ur < 1.5:
                    pw135 += term
            poolPW += pw
            poolSM += s["SM"]
            frac = pw / s["SM"] if s["SM"] > 0 else float("nan")
            if pw / qc > best:
                best, bestL, bestFrac = pw / qc, L, frac
            pred = float("nan")
            if 1.5 < uL < 3.0:
                pred = (math.log(3.0 / uL)
                        * math.log((uL / 3.0) / (uL - 1.5)))
            print("  z=%d L=%d (u_L=%.3f): PW = %.4f | PW/Q_c = "
                  "%.4f | SM/Q_c = %.4f | PW/SM = %.3f | PW135/PW "
                  "= %.3f | Mertens pred = %.4f"
                  % (z, L, uL, pw, pw / qc, s["SM"] / qc, frac,
                     pw135 / pw if pw > 0 else float("nan"), pred))
        resumen[z] = (best, bestL, bestFrac,
                      poolPW / poolSM if poolSM > 0 else float("nan"))
        print("  z=%d: max PW/Q_c = %.4f en L=%s (C(z) doc 123b = "
              "%.4f; PW = %.1f%% del maximo) | PW/SM en ese nivel "
              "= %.3f | pooled PW/SM = %.3f"
              % (z, best, bestL, C123[z], 100.0 * best / C123[z],
                 bestFrac, resumen[z][3]))
    print("  => PASS (medicion; sin criterio de corte)")
    return ok, resumen


def check_K3(resumen):
    print("K3 Q_c(z)*log z y cociente predicho vs medido:")
    for z in (37, 53, 61, 67, 71, 101, 151, 211, 307, 503, 1009):
        print("  z=%4d: Q_c = %.5f | Q_c*log z = %.4f"
              % (z, Qc(z), Qc(z) * math.log(z)))
    cQ = Qc(1009) * math.log(1009)
    print("  con c_Q ~ %.3f (z=1009): predicho max PW/Q_c ~ "
          "(0.1167/c_Q)*log z:" % cQ)
    for z in (37, 53, 61, 67, 71):
        print("  z=%d: predicho = %.4f | medido = %.4f"
              % (z, 0.1167 / cQ * math.log(z), resumen[z][0]))
    print("  => PASS (medicion; sin criterio de corte)")
    return True


def main():
    ok = True
    ok &= check_K0()
    print()
    ok &= check_K1()
    print()
    ok2, resumen = sweep_K2()
    ok &= ok2
    print()
    check_K3(resumen)
    print()
    print("Tiempo total: %.1f s" % (time.time() - t0))
    print("PASS" if ok else "FAIL")
    return 0 if ok else 1


if __name__ == "__main__":
    sys.exit(main())
