#!/usr/bin/env python3
# 116a_rem_pers.py  (Doc 116, seccion 116a)
# Residuo (REM-PERS-small) del doc 114: la familia es VACIA [T].
#
# Afirmacion verificada aqui por enumeracion (la prueba va en el .md):
# para x en Vhat^{(r)} con mu(x) = -1 (omega(x) impar), NINGUN fallo
# fusionado de r*x cae en un primo <= q_min(x):
#   - p = q_min(x) (unico primo propio <= q_min): su posicion propia
#     es omega(x), IMPAR; por el teorema R1 del doc 114 sec. 1.2 [T]
#     todo fallo R1 cae en posicion propia par. No puede fallar.
#   - p | r, p = q_min(x): imposible, (x, r) = 1.
#   - p | r, p < q_min(x) (R2 bajo): la condicion del modelo en la
#     posicion propia omega(x) (impar, si lleva condicion) da
#     r * q_min^2 * x < L; un fallo fusionado en p daria
#     L <= p^2 * x * r / rad_{<p}(r) <= p^2 * x * r, luego
#     q_min^2 < p^2, contra p < q_min. No puede fallar.
# Luego TODO fallo de todo x en Rem^- esta en un primo > q_min(x),
# la persistencia [T] del doc 114 aplica sin excepciones, psi(x) en
# Rem^+ siempre, y (REM-PERS-small) es vacuo: contribucion 0.
# Corolario: ningun primo x esta en Rem (Lema C del .md).
#
# Chequeos:
#   J0  cruce con el doc 114: Rem^- en marcos = 39 (K4); instancias
#       totales 4113 (marcos) + 185 (estres z=53, K7).
#   J1  hipotesis de los lemas en ventana: r en V^+(L) siempre; la
#       condicion del modelo en la posicion omega(x) se cumple:
#       r * q_min^2 * x < L para todo x en Vhat^- (Rem o Both).
#   J2  nucleo: sobre todo x en Vhat^- (marcos A/B + estres), ningun
#       fallo fusionado en primo <= q_min(x); en particular ningun
#       R2 con p | r, p < q_min, y ningun fallo en q_min mismo.
#   J3  Lema B directo: para todo x en Vhat^- y todo p | r con
#       p < q_min(x): p^2 * x * r / rad_{<p}(r) < L (p no puede
#       fallar); se reporta la razon maxima contra L.
#   J4  familia residual vacia: x en Rem^- con todos sus fallos en
#       primos <= q_min: 0. Todo Rem^- tiene >= 1 fallo (todos
#       > q_min). Ningun primo en Rem (cualquier signo de mu).
#   J5  aterrizaje incondicional: psi(x) = x/q_min(x) en Rem^+
#       para todo x en Rem^- (marcos 39/39; estres: conteo libre).
# Uso: python3 116a_rem_pers.py | tee salida_doc116a_local.txt
import importlib.util
import math
import sys
import time

spec = importlib.util.spec_from_file_location(
    "doc114", "114_desp_sign_banda.py")
d114 = importlib.util.module_from_spec(spec)
spec.loader.exec_module(d114)

lados = d114.lados
merged_failures = d114.merged_failures
Level111 = d114.Level111
frame_levels = d114.frame_levels


def rad_below(rpf, p):
    """rad_{<p}(r) = producto de los primos de r menores que p."""
    out = 1
    for q in rpf:
        if q < p:
            out *= q
    return out


def nuevo_acc():
    return dict(n_inst=0, r_fuera_V=0, VhatMinus=0, bad_II=0,
                n_fails=0, fail_le_qmin=0, fail_at_qmin=0,
                R2_below=0, pairsB=0, badB=0, maxratioB=0.0,
                RemMinus=0, Rem_sin_fallos=0, residual=0,
                prime_in_Rem=0, partner_ok=0, partner_bad=0,
                ejemplos=[])


def analiza(tag, z, a, r, ld, acc):
    lv = a.lv
    L = a.L
    acc["n_inst"] += 1
    if r not in lv.Vset:
        acc["r_fuera_V"] += 1
    pfmap = {rr: pf for (rr, pf, mu, phi) in lv.enum}
    rpf = pfmap[r]
    for (dp, (lado_, dpf, dmu, dphi)) in ld.items():
        if lado_ == "Rem" and len(dpf) == 1:
            acc["prime_in_Rem"] += 1
        if dmu != -1 or dp <= 1 or lado_ not in ("Rem", "Both"):
            continue
        acc["VhatMinus"] += 1
        qmin = dpf[-1]
        # J1: condicion del modelo en la posicion propia omega(x)
        if not (r * qmin * qmin * dp < L):
            acc["bad_II"] += 1
        # J3: Lema B directo sobre todo p | r bajo q_min
        for p in rpf:
            if p < qmin:
                lhs = p * p * dp * (r // rad_below(rpf, p))
                acc["pairsB"] += 1
                if not (lhs < L):
                    acc["badB"] += 1
                acc["maxratioB"] = max(acc["maxratioB"], lhs / L)
        # J2: fallos fusionados, todos deben ser > qmin
        fails = merged_failures(dpf, rpf, L)
        for (p, o, m) in fails:
            acc["n_fails"] += 1
            if p <= qmin:
                acc["fail_le_qmin"] += 1
                if len(acc["ejemplos"]) < 5:
                    acc["ejemplos"].append((tag, z, L, r, dp, p, o))
            if p == qmin:
                acc["fail_at_qmin"] += 1
            if o == "r" and p < qmin:
                acc["R2_below"] += 1
        if lado_ == "Rem":
            acc["RemMinus"] += 1
            if not fails:
                acc["Rem_sin_fallos"] += 1
            if fails and all(p <= qmin for (p, o, m) in fails):
                acc["residual"] += 1
            yy = dp // qmin
            if yy in ld and ld[yy][0] == "Rem" and ld[yy][2] == 1:
                acc["partner_ok"] += 1
            else:
                acc["partner_bad"] += 1


def instancias_estres():
    """Replica la construccion del K7 del doc 114 (z=53, marco A)."""
    z = 53
    D, y = z ** 5, int(z ** (10.0 / 3.0))
    out = []
    for lv110 in frame_levels(z, 3, (), D, y, nsamp=10):
        a = Level111(z, 3, (), lv110.L)
        vis = sorted(r for r in a.vis if r > 1)
        if len(vis) > 120:
            idx = sorted({round(i * (len(vis) - 1) / 119.0)
                          for i in range(120)})
            vis = [vis[i] for i in idx]
        for r in vis:
            out.append(("S53", z, a, r, lados(a.lv, r)))
    return out


def reporta(nombre, acc):
    print("  [%s] instancias=%d | x en Vhat^-: %d | fallos "
          "fusionados totales: %d" % (nombre, acc["n_inst"],
                                      acc["VhatMinus"], acc["n_fails"]))
    print("       r fuera de V^+(L): %d | violaciones de (II) "
          "r*qmin^2*x<L: %d" % (acc["r_fuera_V"], acc["bad_II"]))
    print("       fallos en primo <= qmin: %d (en qmin: %d, R2 bajo: "
          "%d) | ejemplos: %s" % (acc["fail_le_qmin"],
                                  acc["fail_at_qmin"], acc["R2_below"],
                                  acc["ejemplos"]))
    print("       Lema B: pares (x, p|r, p<qmin) = %d | violaciones: "
          "%d | max [p^2*x*r/rad]/L = %.4f" % (acc["pairsB"],
                                               acc["badB"],
                                               acc["maxratioB"]))
    print("       Rem^-: %d | sin fallos: %d | residuales (todos los "
          "fallos <= qmin): %d | primos en Rem: %d"
          % (acc["RemMinus"], acc["Rem_sin_fallos"], acc["residual"],
             acc["prime_in_Rem"]))
    print("       psi(x) en Rem^+: %d | fuera: %d"
          % (acc["partner_ok"], acc["partner_bad"]))


def main():
    t0 = time.time()
    print("== DOC 116a: residuo (REM-PERS-small) — familia vacia ==")
    print("Solo stdlib. Maquinaria de los docs 113/114 por import.")
    print()
    ok = True

    accF = nuevo_acc()
    for (tag, z, a, r, ld) in d114.instancias():
        analiza(tag, z, a, r, ld, accF)
    print("J0-J5 marcos A/B (instancias del doc 114):")
    reporta("marcos", accF)
    okF = (accF["n_inst"] == 4113 and accF["RemMinus"] == 39
           and accF["r_fuera_V"] == 0 and accF["bad_II"] == 0
           and accF["fail_le_qmin"] == 0 and accF["fail_at_qmin"] == 0
           and accF["R2_below"] == 0 and accF["badB"] == 0
           and accF["residual"] == 0 and accF["Rem_sin_fallos"] == 0
           and accF["prime_in_Rem"] == 0 and accF["partner_bad"] == 0
           and accF["partner_ok"] == accF["RemMinus"])
    print("  => marcos: %s" % ("PASS" if okF else "FAIL"))
    ok &= okF
    print()

    accS = nuevo_acc()
    for (tag, z, a, r, ld) in instancias_estres():
        analiza(tag, z, a, r, ld, accS)
    print("J0-J5 estres z=53 (construccion del K7 del doc 114):")
    reporta("estres", accS)
    okS = (accS["n_inst"] == 185 and accS["r_fuera_V"] == 0
           and accS["bad_II"] == 0 and accS["fail_le_qmin"] == 0
           and accS["fail_at_qmin"] == 0 and accS["R2_below"] == 0
           and accS["badB"] == 0 and accS["residual"] == 0
           and accS["Rem_sin_fallos"] == 0
           and accS["prime_in_Rem"] == 0 and accS["partner_bad"] == 0
           and accS["partner_ok"] == accS["RemMinus"])
    print("  => estres: %s (Rem^- en estres: %d; el doc 114 conto "
          "36 socios)" % ("PASS" if okS else "FAIL", accS["RemMinus"]))
    ok &= okS
    print()
    print("Tiempo total: %.1f s" % (time.time() - t0))
    print("PASS" if ok else "FAIL")
    return 0 if ok else 1


if __name__ == "__main__":
    sys.exit(main())
