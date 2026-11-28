#!/usr/bin/env python3
# 123b_outmass.py  (Doc 123b)
# (OUT-MASS): estructura de Out(r) = complemento de V^(r) u Vhat^(r),
# ruta con signo sum_capa (M-P) = sum delta - sum E_Rem, medicion
# extendida (z hasta 71) y el intento de cota puntual cerrada.
#   K0  identidades [T] verificadas: delta = E_Rem + M - P por
#       instancia; Vhat^(r) = fibra de V^+(L/r) (predicados
#       identicos); sum_{r en enum} delta(r)/phi(r) = 0 (via
#       sum S/phi = 1 y sum Pi/phi = 1); forma cerrada de la masa
#       lejana M_far; cota puntual M(r) <= (B_r - Pi(r))/2.
#   K1  marcos estandar del doc 120a (continuidad con 122a):
#       por nivel SM, SP, Snet = SM - SP contra Q_c; anclas
#       max SM/Q_c = 0.0533 / 0.1956 / 0.0059 / 0.3152.
#   K2  marcos consistentes (cap=4, near_u=2.0) z = 29..71 y
#       niveles grandes L = z^5//5 (z = 53, 61, 67): tendencia
#       de la constante en z y fraccion de cancelacion.
#   K3  concentracion: instancias top, bandas de u_r, omega(r),
#       reparto lejos/cerca de M.
# Solo stdlib; importa la maquinaria de los docs 110/120a.
import importlib.util
import math
import sys
import time

spec = importlib.util.spec_from_file_location(
    "d120a", "120a_banda_delta_hi.py")
d120 = importlib.util.module_from_spec(spec)
spec.loader.exec_module(d120)
spec110 = importlib.util.spec_from_file_location(
    "d110", "110_segundo_momento_defectos.py")
d110 = importlib.util.module_from_spec(spec110)
spec110.loader.exec_module(d110)

Level110 = d120.Level110
band_instances = d120.band_instances
frame_levels = d120.frame_levels
lados = d120.lados
FRAMES = d120.FRAMES
in_Vplus = d110.in_Vplus
in_Vhat = d120.in_Vhat
TOL = 1e-9


def in_Vplus_at(pf, L, r):
    """in_Vplus(pf, L/r) leido sobre Q (comparaciones exactas):
    el mundo modelo en r es la fibra de V^+ al nivel reducido."""
    from fractions import Fraction
    Lh = Fraction(L, r)
    if not pf:
        return True
    if pf[0] >= Lh:
        return False
    prod = 1
    for i, p in enumerate(pf):
        m = i + 1
        prod *= p
        if m % 2 == 1 and p * p * prod >= Lh:
            return False
    return True


def decompose2(lv, r, rpf, Pi_r):
    """Descomposicion de Out(r): M/P totales, lejos (r*d' >= L) y
    cerca (r*d' < L); E_Rem, E_Add; sumas T_<, U_< sobre (d',r)=1
    con r*d' < L; B_r = prod_{p permitido, p no|r} (1 + 1/(p-1))."""
    L = lv.L
    ld = lados(lv, r)
    M = P = Mf = Mn = Pf = Pn = 0.0
    E_Add = sum_Rem = 0.0
    Tlt = Ult = 0.0
    n_out_m = n_out_p = n_near_m = 0
    for dp, (cls, dpf, dmu, dphi) in ld.items():
        a = dmu / dphi
        if dp * r < L:
            Tlt += a
            Ult += 1.0 / dphi
        if cls == "Add":
            E_Add += a
        elif cls == "Rem":
            sum_Rem += a
        elif cls == "Neither":
            far = dp * r >= L
            if dmu == -1:
                M += 1.0 / dphi
                n_out_m += 1
                if far:
                    Mf += 1.0 / dphi
                else:
                    Mn += 1.0 / dphi
                    n_near_m += 1
            else:
                P += 1.0 / dphi
                n_out_p += 1
                if far:
                    Pf += 1.0 / dphi
                else:
                    Pn += 1.0 / dphi
    Br = 1.0
    rset = set(rpf)
    for p in lv.allow:
        if p not in rset:
            Br *= (1.0 + 1.0 / (p - 1))
    return dict(M=M, P=P, M_far=Mf, M_near=Mn, P_far=Pf, P_near=Pn,
                E_Rem=-sum_Rem, E_Add=E_Add, T_lt=Tlt, U_lt=Ult,
                B_r=Br, Pi_r=Pi_r, n_out_m=n_out_m, n_out_p=n_out_p,
                n_near_m=n_near_m)


def layer_records(lv, L, lnz):
    """Instancias de capa alta no degenerada (u_r >= 1) con desc."""
    for inst in band_instances(lv, L, lnz):
        (r, pf, phi, ur, Pi, sh, nhat, dif, E, dl) = inst
        if ur < 1.0 or nhat == 1:
            continue
        d = decompose2(lv, r, pf, Pi)
        yield (r, pf, phi, ur, Pi, dif, E, dl, d)


# ------------------------------------------------------------------ K0
def check_K0():
    print("K0 identidades [T] sobre la capa y sobre enum completo:")
    ok = True
    frames = [("A", 29, 29 ** 5, int(29 ** (10.0 / 3.0)), None, None),
              ("A-estres", 53, 53 ** 5, int(53 ** (10.0 / 3.0)), 4, 2.0)]
    for (tag, z, D, y, cap, near_u) in frames:
        lnz = math.log(z)
        for L in frame_levels(z, D, y, cap, near_u):
            lv = Level110(z, 3, (), L)
            # suma con signo sobre TODOS los conductores del enum
            sS = math.fsum(lv.S[r] / rphi
                           for (r, rpf, rmu, rphi) in lv.enum)
            sP = math.fsum(lv.Pi[r] / rphi
                           for (r, rpf, rmu, rphi) in lv.enum)
            sd = math.fsum(lv.delta[r] / rphi
                           for (r, rpf, rmu, rphi) in lv.enum)
            e_id = e_far = e_tot = 0.0
            mism = 0
            viol = 0
            n = 0
            for (r, pf, phi, ur, Pi, dif, E, dl, d) in \
                    layer_records(lv, L, lnz):
                n += 1
                e_id = max(e_id,
                           abs(dl - (d["E_Rem"] + d["M"] - d["P"])))
                # forma cerrada de la masa lejana:
                # unsigned_far = B_r - U_< ; signed_far = Pi - T_<
                Mf_form = ((d["B_r"] - d["U_lt"])
                           - (Pi - d["T_lt"])) / 2.0
                e_far = max(e_far, abs(Mf_form - d["M_far"]))
                # masa total mu=-1 coprima: (B_r - Pi)/2 >= M(r)
                tot_m = (d["B_r"] - Pi) / 2.0
                if d["M"] > tot_m + 1e-9:
                    viol += 1
                # identidad de predicados Vhat = V^+(L/r):
                # instancias de capa (<=6) + muestra de visibles
                sample = []
                if n <= 6:
                    sample.append(r)
                if n == 1:
                    vis = sorted(d120.visible_set(lv))[:400]
                    step = max(1, len(vis) // 10)
                    sample += [v for v in vis[::step]
                               if v != r and v > 1][:10]
                for rr in sample:
                    for (dp, dpf, dmu, dphi) in lv.enum:
                        if math.gcd(dp, rr) != 1:
                            continue
                        a1 = (dp * rr < L) and in_Vhat(dpf, L, rr)
                        a2 = in_Vplus_at(dpf, L, rr)
                        if a1 != a2:
                            mism += 1
            print("  [%s] z=%d L=%d: sum S/phi = %.12f | sum Pi/phi"
                  " = %.12f | sum delta/phi = %.3e | capa n=%d"
                  % (tag, z, L, sS, sP, sd, n))
            print("      |delta-(E_Rem+M-P)| max = %.2e | M_far "
                  "forma cerrada err = %.2e | Vhat=V^+(L/r) "
                  "discrepancias = %d | cota puntual violada = %d"
                  % (e_id, e_far, mism, viol))
            ok &= abs(sS - 1.0) < 1e-7 and abs(sP - 1.0) < 1e-7
            ok &= abs(sd) < 1e-7 and e_id < TOL and e_far < 1e-7
            ok &= mism == 0 and viol == 0
    print("  => %s" % ("PASS" if ok else "FAIL"))
    return ok


def level_sums(lv, L, lnz):
    """Agregados de capa por nivel (unidades /phi, sin Q_c)."""
    SM = SP = Snet = Sd = SEr = Sone = 0.0
    recs = []
    for (r, pf, phi, ur, Pi, dif, E, dl, d) in layer_records(lv, L, lnz):
        SM += d["M"] / phi
        SP += d["P"] / phi
        Snet += (d["M"] - d["P"]) / phi
        Sd += dl / phi
        SEr += d["E_Rem"] / phi
        Sone += 1.0 / phi
        recs.append((r, pf, phi, ur, d))
    return dict(SM=SM, SP=SP, Snet=Snet, Sd=Sd, SEr=SEr, Sone=Sone,
                recs=recs, e_split=abs(Snet - (Sd - SEr)))


# ------------------------------------------------------------------ K1
ANCHORS = {("A", 29): 0.0533, ("A", 37): 0.1956, ("B", 29): 0.0059,
           ("A-estres", 53): 0.3152}


def sweep_K1():
    print("K1 marcos estandar (doc 120a): max sobre niveles de "
          "SM/Q_c, SP/Q_c, Snet/Q_c; cancelacion pooled SP/SM:")
    ok = True
    for (tag, z, D, y, cap, near_u) in FRAMES:
        lnz = math.log(z)
        bestM = bestP = bestN = -1.0
        pM = pP = pN = 0.0
        e_split = 0.0
        sone_best = 0.0
        for L in frame_levels(z, D, y, cap, near_u):
            lv = Level110(z, 3, (), L)
            s = level_sums(lv, L, lnz)
            qc = lv.qc
            bestM = max(bestM, s["SM"] / qc)
            bestP = max(bestP, s["SP"] / qc)
            bestN = max(bestN, s["Snet"] / qc)
            if s["SM"] / qc >= bestM - 1e-12:
                sone_best = s["Sone"] / qc
            pM += s["SM"] / qc
            pP += s["SP"] / qc
            pN += s["Snet"] / qc
            e_split = max(e_split, s["e_split"])
        print("  [%s] z=%d: max SM = %.4f | max SP = %.4f | max "
              "Snet = %.4f | pooled SM = %.4f SP = %.4f Snet = %.4f"
              " | SP/SM pooled = %.3f | sum 1/phi/Q_c en nivel de "
              "max SM = %.4f | err split = %.2e"
              % (tag, z, bestM, bestP, bestN, pM, pP, pN,
                 pP / max(pM, 1e-300), sone_best, e_split))
        key = (tag, z)
        if key in ANCHORS:
            good = abs(bestM - ANCHORS[key]) < 2e-3
            ok &= good
            print("      ancla doc 122a: %.4f -> %s"
                  % (ANCHORS[key], "OK" if good else "DISCREPA"))
        ok &= e_split < 1e-8
    print("  => %s" % ("PASS" if ok else "FAIL"))
    return ok


# ------------------------------------------------------------------ K2
def sweep_K2():
    print("K2 marcos consistentes extendidos (D = z^5, y = z^(10/3),"
          " cap=4, near_u=2.0) y niveles grandes L = z^5//5:")
    ok = True
    trend = {}
    allrecs = []
    for z in (29, 37, 53, 61, 67, 71):
        lnz = math.log(z)
        levels = frame_levels(z, z ** 5, int(z ** (10.0 / 3.0)),
                              4, 2.0)
        levels = sorted(set(levels) | {z ** 5 // 5})
        bestM = bestP = bestN = best1 = -1.0
        bestL = None
        pM = pP = 0.0
        farM = totM = 0.0
        for L in levels:
            lv = Level110(z, 3, (), L)
            s = level_sums(lv, L, lnz)
            qc = lv.qc
            bestM = max(bestM, s["SM"] / qc)
            bestP = max(bestP, s["SP"] / qc)
            bestN = max(bestN, s["Snet"] / qc)
            if s["SM"] / qc >= bestM:
                bestL = L
            best1 = max(best1, s["Sone"] / qc)
            pM += s["SM"] / qc
            pP += s["SP"] / qc
            ok &= s["e_split"] < 1e-8
            for (r, pf, phi, ur, d) in s["recs"]:
                farM += d["M_far"] / phi / qc
                totM += d["M"] / phi / qc
                allrecs.append((z, L, r, pf, phi, ur, d, qc))
        trend[z] = bestM
        print("  z=%d: max SM = %.4f (en L=%s) | max SP = %.4f | "
              "max Snet = %.4f | pooled SP/SM = %.3f | M_far/M "
              "pooled = %.3f | max sum 1/phi/Q_c = %.4f"
              % (z, bestM, bestL, bestP, bestN,
                 pP / max(pM, 1e-300),
                 farM / max(totM, 1e-300), best1))
    print("  tendencia de la constante C(z) = max SM/Q_c "
          "(z=29: capa vacia en los niveles capados; el marco "
          "estandar da 0.0533 — se excluye de las pendientes):")
    zs = [z for z in sorted(trend) if trend[z] > 0]
    for i, z in enumerate(zs):
        c = trend[z]
        extra = ""
        if i > 0:
            z0 = zs[i - 1]
            slope = (math.log(c / trend[z0])
                     / math.log(math.log(z) / math.log(z0)))
            extra = " | pendiente en log-log(z) vs z prev: %.2f" \
                    % slope
        print("    z=%3d: C = %.4f | C/log z = %.4f | C/(log z)^2 "
              "= %.4f%s" % (z, c, c / math.log(z),
                            c / math.log(z) ** 2, extra))
    print("  => %s" % ("PASS" if ok else "FAIL"))
    return ok, allrecs


# ------------------------------------------------------------------ K3
def sweep_K3(allrecs):
    print("K3 concentracion de M (instancias de los marcos K2):")
    tot = sum(d["M"] / phi / qc
              for (z, L, r, pf, phi, ur, d, qc) in allrecs)
    ranked = sorted(allrecs,
                    key=lambda t: -(t[6]["M"] / t[4] / t[7]))
    top = ranked[:12]
    share5 = sum(t[6]["M"] / t[4] / t[7] for t in ranked[:5]) \
        / max(tot, 1e-300)
    share20 = sum(t[6]["M"] / t[4] / t[7] for t in ranked[:20]) \
        / max(tot, 1e-300)
    print("  instancias pooled: %d | SM pooled = %.4f | top-5 "
          "cargan %.1f%% | top-20 cargan %.1f%%"
          % (len(allrecs), tot, 100 * share5, 100 * share20))
    for (z, L, r, pf, phi, ur, d, qc) in top:
        print("    z=%d L=%d r=%d (omega=%d, u_r=%.3f): M/phi/Q_c "
              "= %.4f | M = %.4f | cerca/total = %.3f | n_aus- = %d"
              " (cerca %d)"
              % (z, L, r, len(pf), ur, d["M"] / phi / qc, d["M"],
                 d["M_near"] / max(d["M"], 1e-300), d["n_out_m"],
                 d["n_near_m"]))
    bins = {}
    for (z, L, r, pf, phi, ur, d, qc) in allrecs:
        b = round(ur, 1)
        bins[b] = bins.get(b, 0.0) + d["M"] / phi / qc
    print("  SM por banda de u_r (redondeada a 0.1): "
          + " | ".join("%.1f: %.4f" % (b, bins[b])
                       for b in sorted(bins)))
    bom = {}
    for (z, L, r, pf, phi, ur, d, qc) in allrecs:
        bom[len(pf)] = bom.get(len(pf), 0.0) + d["M"] / phi / qc
    print("  SM por omega(r): "
          + " | ".join("%d: %.4f" % (k, bom[k]) for k in sorted(bom)))
    return True



# ------------------------------------------------------------------ K4
def check_K4():
    print("K4 estructura fina en las instancias top (r primo): "
          "reparto de M por omega(d') y criterio exacto de los "
          "ausentes cercanos con d' primo:")
    ok = True
    tops = [(37, 1286, 7), (53, 2585, 7), (61, 3450, 11),
            (67, 4024, 11), (71, 4421, 11)]
    for (z, L, r) in tops:
        lv = Level110(z, 3, (), L)
        ld = lados(lv, r)
        Mn1 = Mn3 = Mf1 = Mf3 = 0.0
        near_primes = set()
        for dp, (cls, dpf, dmu, dphi) in ld.items():
            if cls != "Neither" or dmu != -1:
                continue
            if dp * r < L:
                if len(dpf) == 1:
                    Mn1 += 1.0 / dphi
                    near_primes.add(dp)
                else:
                    Mn3 += 1.0 / dphi
            else:
                if len(dpf) == 1:
                    Mf1 += 1.0 / dphi
                else:
                    Mf3 += 1.0 / dphi
        # criterio [T]: q primo permitido, q != r, q^3*r >= L
        # (falla V^+(L/r)), max(r,q)^3 >= L (falla V^+(L)), q*r < L
        crit = {q for q in lv.allow
                if q != r and q ** 3 * r >= L
                and max(r, q) ** 3 >= L and q * r < L}
        match = (crit == near_primes)
        ok &= match
        print("  z=%d L=%d r=%d: M_cerca(omega=1) = %.4f | "
              "M_cerca(omega>=3) = %.4f | M_lejos(omega=1) = %.4f "
              "| M_lejos(omega>=3) = %.4f | conjunto de primos "
              "cercanos ausentes = %d, criterio exacto: %s"
              % (z, L, r, Mn1, Mn3, Mf1, Mf3, len(near_primes),
                 "OK" if match else "DISCREPA"))
    print("  => %s" % ("PASS" if ok else "FAIL"))
    return ok

def main():
    t0 = time.time()
    print("== DOC 123b: (OUT-MASS) estructura, signo, medicion "
          "extendida, cota ==")
    ok = True
    ok &= check_K0()
    print()
    ok &= sweep_K1()
    print()
    ok2, allrecs = sweep_K2()
    ok &= ok2
    print()
    sweep_K3(allrecs)
    print()
    ok &= check_K4()
    print()
    print("Tiempo total: %.1f s" % (time.time() - t0))
    print("PASS" if ok else "FAIL")
    return 0 if ok else 1


if __name__ == "__main__":
    sys.exit(main())

