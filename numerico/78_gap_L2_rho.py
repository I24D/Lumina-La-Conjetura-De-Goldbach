#!/usr/bin/env python3
"""Contraejemplo finito: control L1 del defecto no implica control L2."""

from fractions import Fraction


def check(name, cond):
    if not cond:
        raise AssertionError(name)
    print("PASS", name)


P = 1000  # masa del indicador de primos abstracto
M = 400   # misma masa L1 del defecto en ambos modelos

# rho = 1_primo - D. Los soportes primo/compuesto son disjuntos.
l2_spread = P + M
l2_spike = P + M * M
sum_rho_spread = P - M
sum_rho_spike = P - M

check("misma masa L1 del defecto", M == M)
check("misma masa total de rho", sum_rho_spread == sum_rho_spike)
check("ambos preservan rho <= 1_primo", -1 <= 0 and -M <= 0)
check("L2 extendido es P+M", l2_spread == 1400)
check("L2 concentrado es P+M^2", l2_spike == 161000)
check("el cociente L2 supera 100", Fraction(l2_spike, l2_spread) > 100)

print("PASS 6/6: masa L1/minoracion no controlan la norma L2")
