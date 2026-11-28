#!/usr/bin/env python3
"""Independent audit of Doc 123: POS arithmetic and joint pooled concentration."""
import contextlib
import importlib.util
import io
import math
from fractions import Fraction

# Exact POS bounds for q0=5, separated by parity and the two lowest primes.
bounds = {
    "odd_if_5_not_in_c": Fraction(125, 49),
    "even_if_p_k_5": Fraction(625, 49),
    "even_if_p_k_7": Fraction(875, 121),
    "even_if_p_k_ge_11": Fraction(1375, 169),
}
assert bounds["odd_if_5_not_in_c"] < 3
assert bounds["even_if_p_k_5"] < 13
assert bounds["even_if_p_k_7"] < 8
assert bounds["even_if_p_k_ge_11"] < 9

spec = importlib.util.spec_from_file_location("m123b", "123b_outmass.py")
m = importlib.util.module_from_spec(spec)
spec.loader.exec_module(m)

# Recreate the exact pooled sample used by Doc 123b K2.
with contextlib.redirect_stdout(io.StringIO()):
    ok, recs = m.sweep_K2()
assert ok

def weight(t):
    z, L, r, pf, phi, ur, d, qc = t
    return d["M"] / phi / qc

total = math.fsum(weight(t) for t in recs)
prime = math.fsum(weight(t) for t in recs if len(t[3]) == 1)
band = math.fsum(weight(t) for t in recs if 1.35 <= t[5] < 1.5)
joint = math.fsum(weight(t) for t in recs
                  if len(t[3]) == 1 and 1.35 <= t[5] < 1.5)

# The five dominant prime-conductor examples: no nearby omega(d') >= 3.
tops = [(37, 1286, 7), (53, 2585, 7), (61, 3450, 11),
        (67, 4024, 11), (71, 4421, 11)]
near_composite = []
for z, L, r in tops:
    lv = m.Level110(z, 3, (), L)
    ld = m.lados(lv, r)
    mass = math.fsum(1.0 / dphi for dp, (cls, pf, mu, dphi) in ld.items()
                     if cls == "Neither" and mu == -1
                     and dp * r < L and len(pf) >= 3)
    near_composite.append((z, L, r, mass))
assert all(x[3] < 1e-12 for x in near_composite)

print("POS q0=5 exact bounds:")
for name, value in bounds.items():
    print("  %s = %s = %.9f" % (name, value, float(value)))
print("POS consequence: c in {1,5,7}; no c=11 residual.")
print("Pooled K2 mass: total=%.12f" % total)
print("  omega(r)=1: %.12f (%.6f%%)" % (prime, 100*prime/total))
print("  1.35<=u_r<1.5: %.12f (%.6f%%)" % (band, 100*band/total))
print("  intersection: %.12f (%.6f%%)" % (joint, 100*joint/total))
print("Near composite mass at five K4 instances:")
for row in near_composite:
    print("  z=%d L=%d r=%d: %.3e" % row)
print("PASS")
