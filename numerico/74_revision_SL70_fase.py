"""Verificador finito del diagnostico del doc 74.

PRE-REGISTRO (antes de ejecutar):
H: la traslacion ciclica preserva espectro de potencia y autocorrelacion,
   pero cambia un coeficiente puntual.
Refutacion: diferencia >1e-12 en potencia/autocorrelacion, o ausencia de
cambio en el coeficiente objetivo.
No es evidencia para Goldbach ni para Montgomery.
"""
import cmath

M = 17
shift = 5
a = [0j] * M
b = [0j] * M
a[0] = 1
b[shift] = 1

def dft(x):
    return [sum(x[n] * cmath.exp(-2j * cmath.pi * j * n / M)
                for n in range(M)) for j in range(M)]

def cyclic_autocorr(x):
    return [sum(x[n].conjugate() * x[(n + h) % M] for n in range(M))
            for h in range(M)]

A, B = dft(a), dft(b)
power_error = max(abs(abs(A[j])**2 - abs(B[j])**2) for j in range(M))
auto_error = max(abs(x-y) for x, y in zip(cyclic_autocorr(a), cyclic_autocorr(b)))
point_changed = (a[0] != b[0])
assert power_error < 1e-12
assert auto_error < 1e-12
assert point_changed
print(f"PASS M={M} shift={shift} power_error={power_error:.3g} auto_error={auto_error:.3g} a0={a[0].real:g} b0={b[0].real:g}")
