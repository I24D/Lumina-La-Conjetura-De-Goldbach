#!/usr/bin/env python3
# 75 - SL-65.1: agudeza del umbral 1/10 de Maynard III (Prop. 5.1 / sec. 11).
# Verifica con Fractions exactos: la desigualdad ligante, el hueco en la
# esquina, el termino T2 no ligante, el regimen peor caso, y la equivalencia
# con la forma alternativa Q^3 R^4 <= x^(8/5). Sin numerico pesado: exponentes.
from fractions import Fraction as F

ok = 0
def check(name, cond):
    global ok
    assert cond, "FAIL: " + name
    ok += 1
    print("PASS", name)

# --- 1. La desigualdad ligante: Q^(15/2) R^10 / x^4, con QR = x^(1/2+d) ---
# exponente del cociente = (15/2) th_Q + 10 th_R - 4; con th_Q = 1/2+d-th_R:
# = -1/4 + (15/2) d + (5/2) th_R <= 0  <=>  th_R <= 1/10 - 3d.
d = F(1, 1000)
th_R_max = F(1, 10) - 3 * d
check("umbral 1/10-3d", th_R_max == F(97, 1000))
# verificacion simbolica del algebra:
# (15/2)(1/2+d-th_R)+10 th_R-4 = -1/4+(15/2)d+(5/2)th_R ; en th_R=1/10-3d da 0
lhs = F(-1, 4) + F(15, 2) * d + F(5, 2) * th_R_max
check("exponente 0 en el umbral", lhs == 0)

# --- 2. La esquina: th_R en [0.100,0.105], th_Q en [0.405,0.410] ---
# hueco gamma = (15/2) th_Q + 10 th_R - 4 > 0 en toda la esquina (potencia x^gamma)
esquina = {
    "fondo-izq (0.405,0.100)": (F(405, 1000), F(100, 1000)),
    "fondo (0.41,0.100)": (F(410, 1000), F(100, 1000)),
    "centro (0.4075,0.1025)": (F(4075, 10000), F(1025, 10000)),
    "tope (0.41,0.105)": (F(410, 1000), F(105, 1000)),
    "tope-izq (0.405,0.105)": (F(405, 1000), F(105, 1000)),
}
for name, (tq, tr) in esquina.items():
    g = F(15, 2) * tq + 10 * tr - 4
    print(f"  gamma {name} = {float(g):.4f}  (hueco x^{g})")
    check("hueco>0 " + name, g > 0)
# valores exactos esperados
check("gamma fondo-izq 3/80", F(15, 2) * F(405, 1000) + 10 * F(100, 1000) - 4 == F(3, 80))
check("gamma centro 13/160", F(15, 2) * F(4075, 10000) + 10 * F(1025, 10000) - 4 == F(13, 160))

# --- 3. th_R maximo permitido en funcion de th_Q ---
def th_R_max_at(tq):
    return (4 - F(15, 2) * tq) / 10
check("th_Rmax(0.4075)=0.094375", th_R_max_at(F(4075, 10000)) == F(94375, 1000000))
check("th_Rmax(0.41)=0.0925", th_R_max_at(F(410, 1000)) == F(925, 10000))
check("th_Rmax(0.405)=0.09625", th_R_max_at(F(405, 1000)) == F(9625, 100000))
# hueco en exponente th_R en el centro: 0.1025-0.094375 = 0.008125
check("hueco th_R centro 0.008125",
      F(1025, 10000) - th_R_max_at(F(4075, 10000)) == F(8125, 1000000))

# --- 4. Termino T2 no ligante: Q^(13/2) R^(15/2) / (N x^3) < 1 en la esquina ---
# exponente = (13/2) th_Q + (15/2) th_R - th_N - 3 ; th_N en [th_Q+2d, 1/2-3d]
tq, tr = F(4075, 10000), F(1025, 10000)
dd = tq + tr - F(1, 2)  # delta = 0.0125
for nm, tn in [("Nmax", F(1, 2) - 3 * dd), ("Nmin", tq + 2 * dd)]:
    e = F(13, 2) * tq + F(15, 2) * tr - tn - 3
    print(f"  exp T2 centro {nm} = {float(e):.5f}")
    check("T2<1 centro " + nm, e < 0)

# --- 5. Peor caso R0=d1=E1=1, R1''~R: regimen del Lema 11.8 ---
# R1'' >= (N/Q)^(2/3) con Nmax: (1/2-3d-th_Q)*(2/3) < th_R
dd = tq + tr - F(1, 2)
lhs_exp = (F(1, 2) - 3 * dd - tq) * F(2, 3)
check("R1''~R en regimen 11.8", lhs_exp < tr)
print(f"  (N/Q)^(2/3) exponente = {float(lhs_exp):.4f} < th_R = {float(tr)}")

# --- 6. Forma alternativa: Q^3 R^4 <= x^(8/5) <=> R <= x^(1/10-3d) ---
# Q^3 R^4 = (QR)^3 R = x^(3/2+3d) R <= x^(8/5) <=> R <= x^(1/10-3d)
check("Q^3R^4 equiv 1/10-3d", F(8, 5) - (F(3, 2) + 3 * d) == F(1, 10) - 3 * d)
# R^5 Q^6 <= x^3 <=> Q <= x^(1/2-5d), implicado por R >= x^(6d)
check("R^5Q^6 es R>=x^6d", (3 - (F(5, 2) + 5 * d)) == F(1, 2) - 5 * d)

# --- 7. Ganancia hipotetica: mejora x^-gamma en la cota S6 ---
# th_R max = 0.094375 + gamma/10 en th_Q=0.4075
for objetivo, g_need in [(F(100, 1000), F(5625, 100000)),
                         (F(1025, 10000), F(8125, 100000)),
                         (F(105, 1000), F(10625, 100000))]:
    g = (F(15, 2) * F(4075, 10000) + 10 * objetivo - 4)
    check(f"gamma para th_R={float(objetivo)}", g == g_need)
    print(f"  gamma necesario para th_R={float(objetivo)}: {float(g):.4f}")

print(f"\nPASS {ok}/{ok} - 75_umbral_un_decimo_agudo.py")
