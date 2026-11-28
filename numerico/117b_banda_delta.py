#!/usr/bin/env python3
# 117b_banda_delta.py  (Doc 117b)
# (BAND-delta): reduccion a (DESP-SIGN) + (BAND-MOD); ataque a (BAND-MOD)
# por doble conteo y medicion decisiva en marcos A (z=17,29,37) y B.
# Solo stdlib; importa la maquinaria del doc 110 sin tocarla y reproduce
# in_Vhat del doc 113 verbatim (aritmetica entera).
import importlib.util
import math
import sys
import time

spec = importlib.util.spec_from_file_location(
    "d110", "110_segundo_momento_defectos.py")
d110 = importlib.util.module_from_spec(spec)
spec.loader.exec_module(d110)

Level110 = d110.Level110
sub_divisors = d110.sub_divisors
sample_sprimes = d110.sample_sprimes


def in_Vhat(pf, L, r):
    """Doc 113, verbatim: pertenencia a V^+(L/r; P\\r) en enteros."""
    if not pf:
        return True
    if pf[0] * r >= L:
        return False
    prod = 1
    for i, p in enumerate(pf):
        m = i + 1
        prod *= p
        if m % 2 == 1 and r * p * p * prod >= L:
            return False
    return True


def lam_of_pf(pf, d):
    """lambda(d') = max(d', max_{m impar} p_m^2 P_m): d' en Vhat^(r)
    sii (d',r)=1 y r*lambda(d') < L."""
    lam = d
    prod = 1
    for i, p in enumerate(pf):
        m = i + 1
        prod *= p
        if m % 2 == 1:
            lam = max(lam, p * p * prod)
    return lam


def visible_set(lv):
    vis = set()
    for (d, pf, mu, phi) in lv.V:
        for (r, rmu, rphi) in sub_divisors(pf):
            vis.add(r)
    return vis


def Shat_of(lv, r):
    """Shat(r) = sum_{d' en Vhat^(r)} mu(d')/phi(d'); tambien |Vhat|."""
    L = lv.L
    tot = 0.0
    cnt = 0
    for (dp, dpf, dmu, dphi) in lv.enum:
        if math.gcd(dp, r) != 1:
            continue
        if (dp * r < L) and in_Vhat(dpf, L, r):
            tot += dmu / dphi
            cnt += 1
    return tot, cnt


def That_of(lv, r):
    """T^hat(r) = sum_{d' no en Vhat^(r), (d',r)=1} mu(d')/phi(d')."""
    L = lv.L
    tot = 0.0
    for (dp, dpf, dmu, dphi) in lv.enum:
        if math.gcd(dp, r) != 1:
            continue
        if not ((dp * r < L) and in_Vhat(dpf, L, r)):
            tot += dmu / dphi
    return tot


FRAMES = [("A", 17, 17 ** 5, int(17 ** (10.0 / 3.0))),
          ("A", 29, 29 ** 5, int(29 ** (10.0 / 3.0))),
          ("A", 37, 37 ** 5, int(37 ** (10.0 / 3.0))),
          ("B", 29, 29 ** 6, 29 ** 3)]
U0 = 1.5


# ------------------------------------------------------------------ K0
def check_K0():
    print("K0 identidades base: Shat(1)=S(1); Shat-Pi = -T^hat; "
          "delta = (Shat-Pi) + E:")
    ok = True
    for (z, L) in ((29, 29 ** 2), (37, int(37 ** (5.0 / 3.0)))):
        lv = Level110(z, 3, (), L)
        sh1, _ = Shat_of(lv, 1)
        e1 = abs(sh1 - lv.S[1])
        e2 = 0.0
        e3 = 0.0
        for (r, rpf, rmu, rphi) in lv.enum[:80]:
            sh, _ = Shat_of(lv, r)
            th = That_of(lv, r)
            e2 = max(e2, abs((sh - lv.Pi[r]) + th))
            E = lv.S[r] - sh
            e3 = max(e3, abs(lv.delta[r] - ((sh - lv.Pi[r]) + E)))
        print("  z=%d L=%d: |Shat(1)-S(1)|=%.2e |Shat-Pi+T^hat|=%.2e "
              "|desc|=%.2e" % (z, L, e1, e2, e3))
        ok &= e1 < 1e-12 and e2 < 1e-12 and e3 < 1e-12
    print("  => %s" % ("PASS" if ok else "FAIL"))
    return ok


# ------------------------------------------------------------------ K2
def check_K2():
    print("K2 doble conteo de la suma con signo (niveles sueltos):")
    ok = True
    for (z, L) in ((29, 29 ** 2), (37, int(37 ** (5.0 / 3.0)))):
        lv = Level110(z, 3, (), L)
        # LHS: sum_r T^hat(r)/phi(r) sobre todo el enum
        lhs = 0.0
        for (r, rpf, rmu, rphi) in lv.enum:
            lhs += That_of(lv, r) / rphi
        # RHS: sum_{d'} mu(d')/phi(d') * W(d'),
        # W(d') = sum_{r en enum, (r,d')=1, r*lam(d') >= L} 1/phi(r)
        rhs = 0.0
        for (dp, dpf, dmu, dphi) in lv.enum:
            lam = lam_of_pf(dpf, dp)
            W = 0.0
            for (r, rpf, rmu, rphi) in lv.enum:
                if math.gcd(r, dp) == 1 and r * lam >= L:
                    W += 1.0 / rphi
            rhs += (dmu / dphi) * W
        err = abs(lhs - rhs)
        print("  z=%d L=%d: LHS=%.6e RHS=%.6e err=%.2e | LHS/Q_c=%.4f"
              % (z, L, lhs, rhs, err, lhs / lv.qc))
        ok &= err < 1e-9 * max(1.0, abs(lhs))
    # Dependencia de umbral: W(d') cambia con L (en el doc 110 la suma
    # interna analoga es sobre r|d fijo, independiente de L).
    dp, dpf = 385, (11, 7, 5)
    lam = lam_of_pf(dpf, dp)
    Ws = []
    for (z, L) in ((29, 29 ** 2), (29, 29 ** 3)):
        lv = Level110(z, 3, (), L)
        W = sum(1.0 / rphi for (r, rpf, rmu, rphi) in lv.enum
                if math.gcd(r, dp) == 1 and r * lam >= L)
        Ws.append(W)
    print("  d'=385, lambda=%d: W(L=29^2)=%.4f, W(L=29^3)=%.4f "
          "(distintos: la suma interna depende del umbral L/lambda)"
          % (lam, Ws[0], Ws[1]))
    ok &= abs(Ws[0] - Ws[1]) > 1e-6
    print("  => %s" % ("PASS" if ok else "FAIL"))
    return ok


# ------------------------------------------------------- K1/K3/K4/K5
def sweep():
    print("K1/K3/K4/K5 barrido de banda (u_r<1.5, Pi<1/2, r visible>1):")
    ok = True
    # por marco: maximos sobre niveles de las tres razones / Q_c
    frame_max = {}
    # clasificacion de signo global
    sign_by_om = {}     # om -> [n+, n0, n-, masa+, masa-]
    sign_w1_ur = {}     # bin u_r -> [n+, n-, masa+, masa-]
    deg_mass = {}       # z -> max masa degenerada / Q_c
    deg_cnt = 0
    relmax = 0.0        # max (Shat-Pi)+/Pi fuera de degenerados
    rel_frac_num = 0.0  # masa (Shat-Pi)+ con (Shat-Pi) <= Pi
    rel_frac_den = 0.0
    worst_red = 0.0     # violacion de delta+ <= (Shat-Pi)+
    negE = 0
    tot_inst = 0
    for (tag, z, D, y) in FRAMES:
        spr = sample_sprimes(z, y, 3, (), 60)
        levels = sorted({D // s for s in spr})
        fmax = [0.0, 0.0, 0.0]
        for L in levels:
            lv = Level110(z, 3, (), L)
            info = {r: (pf, phi) for (r, pf, mu, phi) in lv.enum}
            lnz = math.log(z)
            s_d = s_m = s_a = 0.0
            for r in visible_set(lv):
                if r == 1:
                    continue
                pf, phi = info[r]
                ur = math.log(L / r) / lnz
                if ur >= U0:
                    continue
                Pi = lv.Pi[r]
                if Pi >= 0.5:
                    continue
                sh, nhat = Shat_of(lv, r)
                dif = sh - Pi
                E = lv.S[r] - sh
                if E > 1e-12:
                    negE += 1
                tot_inst += 1
                worst_red = max(worst_red,
                                lv.delta[r] - dif if lv.delta[r] > 0
                                and dif < lv.delta[r] else 0.0)
                if lv.delta[r] > 0:
                    s_d += lv.delta[r] / phi
                if dif > 0:
                    s_m += dif / phi
                s_a += abs(dif) / phi
                # signo por omega
                om = min(len(pf), 4)
                row = sign_by_om.setdefault(om, [0, 0, 0, 0.0, 0.0])
                if dif > 1e-15:
                    row[0] += 1
                    row[3] += dif / phi
                elif dif < -1e-15:
                    row[1] += 1
                    row[4] += -dif / phi
                else:
                    row[2] += 1
                if len(pf) == 1:
                    b = 0 if ur < 0.5 else (1 if ur < 1.0 else 2)
                    rw = sign_w1_ur.setdefault(b, [0, 0, 0.0, 0.0])
                    if dif > 1e-15:
                        rw[0] += 1
                        rw[2] += dif / phi
                    elif dif < -1e-15:
                        rw[1] += 1
                        rw[3] += -dif / phi
                if nhat == 1:
                    deg_cnt += 1
                    dm = (1.0 - Pi) / phi / lv.qc
                    if dm > deg_mass.get(z, 0.0):
                        deg_mass[z] = dm
                    tot_d = (1.0 - Pi) / phi
                    s_m -= 0.0  # ya contado en s_m via dif=1-Pi
                else:
                    if dif > 0:
                        relmax = max(relmax, dif / Pi)
                        rel_frac_den += dif / phi
                        if dif <= Pi:
                            rel_frac_num += dif / phi
            fmax[0] = max(fmax[0], s_d / lv.qc)
            fmax[1] = max(fmax[1], s_m / lv.qc)
            fmax[2] = max(fmax[2], s_a / lv.qc)
        frame_max[(tag, z)] = fmax
        print("  [%s] z=%d (%d niveles): max_L sum_banda/Q_c: "
              "delta+ = %.4f | (Shat-Pi)+ = %.4f | |Shat-Pi| = %.4f"
              % (tag, z, len(levels), fmax[0], fmax[1], fmax[2]))
    print("  instancias de banda: %d | E>0 (contra DESP-SIGN): %d | "
          "violacion reduccion delta+<=(Shat-Pi)+: %.2e"
          % (tot_inst, negE, worst_red))
    ok &= negE == 0 and worst_red < 1e-12
    print("  signo de Shat-Pi por omega(r) [n+, n-, n0, masa+, masa-]:")
    for om in sorted(sign_by_om):
        rw = sign_by_om[om]
        print("    omega=%s: %d, %d, %d, %.4f, %.4f"
              % (om, rw[0], rw[1], rw[2], rw[3], rw[4]))
    print("  signo para omega=1 por u_r [n+, n-, masa+, masa-]:")
    for b in sorted(sign_w1_ur):
        rw = sign_w1_ur[b]
        print("    u_r en %s: %d, %d, %.4f, %.4f"
              % (("[0,0.5)", "[0.5,1)", "[1,1.5)")[b],
                 rw[0], rw[1], rw[2], rw[3]))
    print("  degenerados Vhat={1}: %d instancias; max (1-Pi)/phi/Q_c "
          "por z: %s" % (deg_cnt,
                         {z: round(v, 4) for z, v in deg_mass.items()}))
    print("  (Shat-Pi)+ <= Pi(r) concentra %.3f de la masa (Shat-Pi)+ "
          "no degenerada; max (Shat-Pi)+/Pi = %.3f"
          % (rel_frac_num / max(rel_frac_den, 1e-300), relmax))
    print("  => %s" % ("PASS" if ok else "FAIL"))
    return ok, frame_max


def main():
    t0 = time.time()
    print("== DOC 117b: (BAND-delta) via (DESP-SIGN) + (BAND-MOD) ==")
    ok = True
    ok &= check_K0()
    print()
    ok &= check_K2()
    print()
    ok2, _ = sweep()
    ok &= ok2
    print()
    ok &= check_K6()
    print()
    print("Tiempo total: %.1f s" % (time.time() - t0))
    print("PASS" if ok else "FAIL")
    return 0 if ok else 1


# ------------------------------------------------------------------ K6
def check_K6():
    print("K6 descomposicion degenerado/no-degenerado y suma Pi/phi "
          "(niveles representativos):")
    ok = True
    for z in (29, 37, 53):
        allow = [p for p in d110.primes_upto(z - 1)
                 if p != 2 and 3 % p != 0]
        for ustar in (5.0 / 3.0, 2.0):
            L = int(z ** ustar)
            lv = Level110(z, 3, (), L)
            info = {r: (pf, phi) for (r, pf, mu, phi) in lv.enum}
            lnz = math.log(z)
            deg = nondeg = pisum = 0.0
            nd = nn = 0
            relmax = relmax1 = 0.0
            char_ok = True
            for r in visible_set(lv):
                if r == 1:
                    continue
                pf, phi = info[r]
                ur = math.log(L / r) / lnz
                if ur >= U0:
                    continue
                Pi = lv.Pi[r]
                if Pi >= 0.5:
                    continue
                sh, nhat = Shat_of(lv, r)
                dif = sh - Pi
                pisum += Pi / phi
                q0 = next(p for p in allow if p not in pf)
                char_ok &= (nhat == 1) == (q0 ** 3 * r >= L)
                if nhat == 1:
                    deg += (1.0 - Pi) / phi
                    nd += 1
                else:
                    nn += 1
                    if dif > 0:
                        nondeg += dif / phi
                        relmax = max(relmax, dif / Pi)
                        if ur >= 1.0:
                            relmax1 = max(relmax1, dif / Pi)
            print("  z=%d u_L=%.3f L=%d: deg/Q_c=%.4f (n=%d) | "
                  "nondeg+/Q_c=%.4f (n=%d) | sum Pi/phi/Q_c=%.4f | "
                  "max (Shat-Pi)/Pi nondeg=%.3f (u_r>=1: %.3f) | "
                  "carac. deg %s"
                  % (z, math.log(L) / lnz, L, deg / lv.qc, nd,
                     nondeg / lv.qc, nn, pisum / lv.qc, relmax,
                     relmax1, "OK" if char_ok else "FALLA"))
            ok &= char_ok
    print("  => %s" % ("PASS" if ok else "FAIL"))
    return ok


if __name__ == "__main__":
    sys.exit(main())
