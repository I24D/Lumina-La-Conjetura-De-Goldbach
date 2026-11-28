#!/usr/bin/env python3
# 122a_banda_delta_hi_prima.py  (Doc 122a)
# (BAND-delta-HI'): enunciado exacto con el deficit del flujo Add
# (doc 114), medicion ampliada de la desigualdad sumada sobre la
# capa alta no degenerada, y el lema minimo: la identidad
# delta = E_Rem + M - P (masas por signo de Out = complemento de
# V u Vhat) y el residuo (OUT-MASS).
#   K0  identidades sobre la capa: delta = (Shat-Pi) + E;
#       E = E_Add + E_Rem; deficit de flujo = -E_Add (doc 114);
#       delta = E_Rem + M - P; provisiones (i) E_Rem = 0,
#       (ii) E_Add <= 0; BAD = 0 y R_A <= 1 en la capa.
#   K1  por marco y u0* en {1.0, 1.2, 1.4}: maximos sobre niveles
#       de exceso X, deficit D, no cubierto U, masa M (todo /Q_c),
#       y cociente D/X por nivel.
#   K2  subtotal E = 0: conteo, clasificacion (|Add|, |Rem|) y
#       fraccion de U que cargan.
#   K3  lemas puntuales: (P1) exceso <= deficit, (P2) exceso <=
#       k*deficit, y delta^+ <= M bajo (i).
#   K4  (FC-2) en la capa: BAD(r) y HEAD37(r) por conductor.
# Solo stdlib; importa la maquinaria de los docs 110/117b/120a.
import importlib.util
import math
import sys
import time

spec = importlib.util.spec_from_file_location(
    "d120a", "120a_banda_delta_hi.py")
d120 = importlib.util.module_from_spec(spec)
spec.loader.exec_module(d120)

Level110 = d120.Level110
band_instances = d120.band_instances
frame_levels = d120.frame_levels
lados = d120.lados
FRAMES = d120.FRAMES
U0S = (1.0, 1.2, 1.4)
TOL = 1e-12


def decompose(lv, r, primes):
    """Descomposicion completa de un conductor de capa (doc 114/122a)."""
    ld = lados(lv, r)
    E_Add = 0.0
    sum_Rem = 0.0
    M = 0.0
    P = 0.0
    nAdd = nRem = 0
    addminus = []
    for dp, (cls, dpf, dmu, dphi) in ld.items():
        a = dmu / dphi
        if cls == "Add":
            E_Add += a
            nAdd += 1
            if dmu == -1:
                addminus.append((dp, dpf, dphi))
        elif cls == "Rem":
            sum_Rem += a
            nRem += 1
        elif cls == "Neither":
            if dmu == -1:
                M += 1.0 / dphi
            else:
                P += 1.0 / dphi
    flow = GOOD = BAD = 0.0
    maxRA = 0.0
    HEAD37 = 0.0
    nAddPlus = 0
    for dp, (cls, dpf, dmu, dphi) in ld.items():
        if cls == "Add" and dmu == 1:
            nAddPlus += 1
    for (y, ypf, yphi) in addminus:
        ay = 1.0 / yphi
        yset = set(ypf)
        RA = 0.0
        for q in primes:
            if q in yset:
                continue
            v = ld.get(y * q)
            if v is not None and v[0] == "Add" and v[2] == 1:
                RA += 1.0 / (q - 1)
        maxRA = max(maxRA, RA)
        flow += ay * (1.0 - RA)
        if RA <= 1.0:
            GOOD += ay * (1.0 - RA)
        else:
            BAD += ay * (RA - 1.0)
        if ypf and min(ypf) <= 37:
            HEAD37 += ay
    return dict(E_Add=E_Add, E_Rem=-sum_Rem, M=M, P=P, nAdd=nAdd,
                nRem=nRem, flow=flow, GOOD=GOOD, BAD=BAD,
                maxRA=maxRA, HEAD37=HEAD37, nAddminus=len(addminus),
                nAddPlus=nAddPlus)


def layer_instances(lv, L, lnz):
    """Instancias de capa alta no degenerada (u_r >= 1) con desc."""
    primes = sorted({p for (_, pf, _, _) in lv.enum for p in pf})
    for inst in band_instances(lv, L, lnz):
        (r, pf, phi, ur, Pi, sh, nhat, dif, E, dl) = inst
        if ur < 1.0 or nhat == 1:
            continue
        d = decompose(lv, r, primes)
        yield (r, pf, phi, ur, Pi, dif, E, dl, d)


# ------------------------------------------------------------------ K0
def check_K0():
    print("K0 identidades y provisiones sobre la capa (z=29 y "
          "estres z=53):")
    ok = True
    e1 = e2 = e3 = e4 = 0.0
    n = n_rem = n_addpos = n_bad = 0
    maxRA = 0.0
    frames = [("A", 29, 29 ** 5, int(29 ** (10.0 / 3.0)), None, None),
              ("A-estres", 53, 53 ** 5, int(53 ** (10.0 / 3.0)), 4, 2.0)]
    for (tag, z, D, y, cap, near_u) in frames:
        lnz = math.log(z)
        for L in frame_levels(z, D, y, cap, near_u):
            lv = Level110(z, 3, (), L)
            for (r, pf, phi, ur, Pi, dif, E, dl, d) in \
                    layer_instances(lv, L, lnz):
                n += 1
                e1 = max(e1, abs(dl - (dif + E)))
                e2 = max(e2, abs(E - (d["E_Add"] + d["E_Rem"])))
                e3 = max(e3, abs(d["flow"] + d["E_Add"]))
                e4 = max(e4, abs(dl - (d["E_Rem"] + d["M"] - d["P"])))
                if abs(d["E_Rem"]) > 1e-9:
                    n_rem += 1
                if d["E_Add"] > 1e-9:
                    n_addpos += 1
                if d["BAD"] > 1e-12:
                    n_bad += 1
                maxRA = max(maxRA, d["maxRA"])
    print("  instancias de capa descompuestas: %d" % n)
    print("  |delta-((Shat-Pi)+E)| = %.2e | |E-(E_Add+E_Rem)| = %.2e"
          % (e1, e2))
    print("  |flujo - (-E_Add)| = %.2e | |delta-(E_Rem+M-P)| = %.2e"
          % (e3, e4))
    print("  provision (i) violada (|E_Rem|>1e-9): %d | provision "
          "(ii) violada (E_Add>1e-9): %d | BAD>0: %d | max R_A = %.4f"
          % (n_rem, n_addpos, n_bad, maxRA))
    ok &= max(e1, e2, e3, e4) < 1e-9 and n_rem == 0 and n_bad == 0
    print("  => %s" % ("PASS" if ok else "FAIL"))
    return ok


# ------------------------------------------------------------- K1–K4
def sweep():
    print("K1 medicion sumada por marco y u0* (max sobre niveles, "
          "/ Q_c):")
    ok = True
    # pooled (u0* = 1): estadisticas por instancia
    pool = dict(n=0, n_exc=0, n_E0=0, n_E0x=0, U=0.0, U_E0=0.0,
                M=0.0, dl=0.0,
                viol_P1=0, viol_P2=0, viol_M=0, n_rem0=0,
                E0_empty=0, E0_flow=0, BADmax=0.0, n_BAD=0,
                HEAD37min=None, e_flow=0.0, e_id=0.0, maxRA=0.0,
                nAddPlus=0, nAddMinus=0)
    for (tag, z, D, y, cap, near_u) in FRAMES:
        levels = frame_levels(z, D, y, cap, near_u)
        lnz = math.log(z)
        # por u0*: mejor nivel segun U y rango de cocientes D/X
        best = {u: dict(U=-1.0, X=0.0, D=0.0, M=0.0, rat=None,
                        rmin=None, rmax=None, maxX=0.0, maxD=0.0,
                        maxM=0.0)
                for u in U0S}
        for L in levels:
            lv = Level110(z, 3, (), L)
            qc = lv.qc
            acc = {u: [0.0, 0.0, 0.0, 0.0] for u in U0S}  # X,D,U,M
            for (r, pf, phi, ur, Pi, dif, E, dl, d) in \
                    layer_instances(lv, L, lnz):
                exc = dif if dif > 0 else 0.0
                defi = -d["E_Add"]
                unc = max(exc - abs(d["E_Add"]), 0.0)
                vals = (exc / phi, defi / phi, unc / phi,
                        d["M"] / phi)
                for u in U0S:
                    if ur >= u:
                        a = acc[u]
                        for i in range(4):
                            a[i] += vals[i]
                # pooled solo u0* = 1 (ur >= 1 ya garantizado)
                pool["n"] += 1
                pool["e_flow"] = max(pool["e_flow"],
                                     abs(d["flow"] + d["E_Add"]))
                pool["e_id"] = max(
                    pool["e_id"],
                    abs(dl - (d["E_Rem"] + d["M"] - d["P"])))
                pool["maxRA"] = max(pool["maxRA"], d["maxRA"])
                pool["nAddPlus"] += d["nAddPlus"]
                pool["nAddMinus"] += d["nAddminus"]
                if E == 0.0:
                    pool["n_E0x"] += 1
                pool["U"] += unc / phi / qc
                pool["M"] += d["M"] / phi / qc
                pool["dl"] += (dl / phi / qc if dl > 0 else 0.0)
                if abs(d["E_Rem"]) <= 1e-9:
                    pool["n_rem0"] += 1
                    if dl > d["M"] + 1e-9:
                        pool["viol_M"] += 1
                if d["BAD"] > pool["BADmax"]:
                    pool["BADmax"] = d["BAD"]
                if d["BAD"] > 1e-12:
                    pool["n_BAD"] += 1
                if d["nAddminus"] > 0:
                    h = d["HEAD37"]
                    pool["HEAD37min"] = h if pool["HEAD37min"] is \
                        None else min(pool["HEAD37min"], h)
                if exc > 0:
                    pool["n_exc"] += 1
                    if defi < exc - 1e-12:
                        pool["viol_P1"] += 1
                    if defi <= 0.0:
                        pool["viol_P2"] += 1
                    if abs(E) <= 1e-12:
                        pool["n_E0"] += 1
                        pool["U_E0"] += unc / phi / qc
                        if d["nAdd"] == 0 and d["nRem"] == 0:
                            pool["E0_empty"] += 1
                        else:
                            pool["E0_flow"] += 1
            for u in U0S:
                X, Dd, U, M = (v / qc for v in acc[u])
                b = best[u]
                b["maxX"] = max(b["maxX"], X)
                b["maxD"] = max(b["maxD"], Dd)
                b["maxM"] = max(b["maxM"], M)
                if X > 0:
                    rat = Dd / X
                    b["rmin"] = rat if b["rmin"] is None else \
                        min(b["rmin"], rat)
                    b["rmax"] = rat if b["rmax"] is None else \
                        max(b["rmax"], rat)
                if U > b["U"]:
                    b.update(U=U, X=X, D=Dd, M=M,
                             rat=(Dd / X if X > 0 else None))
        for u in U0S:
            b = best[u]
            print("  [%s] z=%d u0*=%.1f: en el nivel de max U: "
                  "U = %.4f | X = %.4f | D = %.4f | M = %.4f | "
                  "D/X = %s || D/X por nivel en [%.3f, %.3f] | "
                  "maxX = %.4f | maxM = %.4f"
                  % (tag, z, u, b["U"], b["X"], b["D"], b["M"],
                     ("%.3f" % b["rat"]) if b["rat"] is not None
                     else "s/n",
                     b["rmin"] if b["rmin"] is not None else -1,
                     b["rmax"] if b["rmax"] is not None else -1,
                     b["maxX"], b["maxM"]))
    print("K2 subtotal E = 0 (pooled, u0*=1, masas /phi/Q_c sumadas "
          "sobre niveles):")
    print("  instancias de capa: %d | con exceso > 0: %d | con "
          "E = 0 y exceso > 0: %d (E == 0.0 exacto: %d)"
          % (pool["n"], pool["n_exc"], pool["n_E0"],
             pool["n_E0x"]))
    print("  identidades pooled: |flujo - (-E_Add)| max = %.2e | "
          "|delta - (E_Rem+M-P)| max = %.2e | max R_A = %.6f"
          % (pool["e_flow"], pool["e_id"], pool["maxRA"]))
    print("  estructura de Add en la capa: elementos Add^+ "
          "totales = %d | elementos Add^- totales = %d"
          % (pool["nAddPlus"], pool["nAddMinus"]))
    print("  U total = %.4f | U de instancias E = 0 = %.4f (%.1f%%)"
          % (pool["U"], pool["U_E0"],
             100.0 * pool["U_E0"] / max(pool["U"], 1e-300)))
    print("  E = 0 con |Add| = |Rem| = 0 (V = Vhat): %d | E = 0 "
          "con flujo no vacio: %d" % (pool["E0_empty"],
                                      pool["E0_flow"]))
    print("K3 lemas puntuales (pooled): exceso > deficit en %d/%d "
          "instancias con exceso (P1 falso); deficit = 0 con "
          "exceso > 0 en %d (P2 falso para todo k)"
          % (pool["viol_P1"], pool["n_exc"], pool["viol_P2"]))
    print("  delta^+ <= M bajo (i): violaciones %d de %d instancias "
          "con E_Rem = 0" % (pool["viol_M"], pool["n_rem0"]))
    print("  M pooled = %.4f contra delta^+ pooled = %.4f "
          "(cociente %.2f)" % (pool["M"], pool["dl"],
                               pool["M"] / max(pool["dl"], 1e-300)))
    print("K4 (FC-2) en la capa: BAD max = %.3e | conductores con "
          "BAD > 0: %d | HEAD37 min (Add^- no vacio) = %s"
          % (pool["BADmax"], pool["n_BAD"],
             ("%.4f" % pool["HEAD37min"]) if pool["HEAD37min"]
             is not None else "s/n"))
    ok &= pool["viol_M"] == 0 and pool["n_BAD"] == 0
    print("  => %s" % ("PASS" if ok else "FAIL"))
    return ok


def main():
    t0 = time.time()
    print("== DOC 122a: (BAND-delta-HI') exacta, medicion sumada, "
          "(OUT-MASS) ==")
    ok = True
    ok &= check_K0()
    print()
    ok &= sweep()
    print()
    print("Tiempo total: %.1f s" % (time.time() - t0))
    print("PASS" if ok else "FAIL")
    return 0 if ok else 1


if __name__ == "__main__":
    sys.exit(main())
