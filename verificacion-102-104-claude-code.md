# Verificación de los docs 102, 103 y 104 (holgura de Cauchy, criba grande restringida, B̂₃)

**Agente:** Claude Code (Anthropic, modelo `claude-opus-5-5`), brazo verificador de la línea Goldbach.
**Fecha:** 2026-11-28.
**Encargo:** de Muse, por orden de Dal (encargo 2 de 4). Verificar E101.1 y E101.2 (doc 102), el fracaso de la criba grande restringida (doc 103) y la pendiente 0,1326 de B̂₃ (doc 104).
**Método:** reimplementación propia desde las definiciones de los docs 101–104, sin importar código de nadie. Incluye: soporte H y H_N; E_d por ψ real; M_r; Δ_r por (98.3) sin caracteres; B_M; W_d; S_BV; BV1* y BV1_full; Diag; Λ₂; δ_ℋ; B̂₁; B̂₂; B̂₃.
Script: `numerico/verificacion_102_104_claude_code.py`. Salida: `numerico/resultados/20261128_verificacion_102_104_claude_code.txt`. **14/14 PASS**, 8 s.
**Estado:** PROPUESTA hasta el contraste de Muse.
**Etiquetas:** [T] teorema, [M] medición, [P] pendiente.

**Goldbach sigue abierta.**

---

## §1. Resumen

| Punto | Veredicto | Evidencia |
|---|---|---|
| **Doc 102** | | |
| E101.1: tablas de §2.2 (\|H_N\|, W_max, S_BV, B_M, slack, S_Cauchy, S_cut, S_max, B^C/B, exceso de compuestos, Diag^C/B^C) en 10 puntos | **Confirmado** a todos los dígitos impresos [M] | D102.1 |
| θ de la holgura: 0,988 ± 0,053 (genérica) y 1,036 ± 0,005 (2^e), con IC y LOO | **Confirmado** [M] | D102.2 |
| Familia colapso N = 720720: B_M^C = 0, B_M = S_BV, slack = W_max | **Confirmado**, y es identidad [T] | D102.3 |
| Factorización slack = S_Cauchy·S_cut·S_max y Diag^P = B_M^P | **Confirmadas** [T] (telescópica; un solo término en r primo) | |
| E101.2: α_W(11..37) = 0,4281 [0,4243; 0,4318] | **Confirmado** el número [M] | G5 de mi verificación del 101 |
| E101.2: "α(z) se estabiliza en ≈ 0,42–0,43" | **Refutado**: α → 1/2 [T]; en z = 307..709 los exponentes locales valen 0,455–0,472 [M] | Teorema 101-CC; `verificacion-conteo-ds` |
| §4: margen de la ruta directa, "≈ 1,2·10⁵ en 10²⁴" y "3 órdenes" | **Refutado en las cifras** (hereda la errata del doc 101 §3.4). Con la aritmética corregida el margen sería mayor. Es extrapolación [M] | §2.3 |
| **Doc 103** | | |
| Derivación de B̂₁ (Cauchy por conductor y criba grande con espaciamiento) | **Confirmada** [T] | §3.1 |
| Sándwich de δ_ℋ | **Confirmado** [T] | §3.1 |
| Λ₂(N), Q_ℋ, r₂, B̂_lo/B_M y B̂_hi/B_M en 6 puntos; pendiente 2,3701 ± 0,0322 | **Confirmados** [M] | D103.1 a D103.3 |
| δ_ℋ exacto en R1 y R3 | **Confirmado** (δ y cocientes). El nº de fracciones difiere en 32 por convención | §3.2 |
| Desglose primo/compuesto (§3) | **Confirmado** [M] | D103.5 |
| B̂₂ (V103.1): 1,709–2,318; pendiente 0,1358 ± 0,0201 | **Confirmado** [M] | D103.4 |
| [P-rev] 103.2: δ_ℋ exacto en soportes grandes, "no computable" | **Resuelto en S11 y S13**: B̂₁ real está 9,3× y 5,4× por encima del extremo favorable | §3.3 |
| **Doc 104** | | |
| Lecturas L1 y L2 de "parte nueva" | **Confirmado**: bajo L1, B̂₃ = B̂₂ [T]. El [P-rev] 104.2 se resuelve contra el doc 94 íntegro | §4.1 |
| Estructura de clases {h₀, h₀s} con signos opuestos | **Confirmada** [T] | §4.1 |
| B̂₃/B_M: 1,679–2,262; pendiente 0,1326 ± 0,0195, IC ±0,0542, LOO [0,1168; 0,1441] | **Confirmado** [M] | D104.1 |
| Descomposición 0,1358 = 0,0032 + 0,1326; ρ_intra = 0,987–0,992 | **Confirmada** [T por linealidad de MCO; M] | D104.2, D104.3 |
| [P-rev] 104.3: B̂₃ en S19 y S23 | **Hecho**: la pendiente de 10 puntos es 0,1168, LOO [0,1074; 0,1214], y C2 sigue fallando. Matiz: en la familia 2^e sola queda a caballo | §4.2 |
| [P-rev] 102.4: max_{d≥2} W_d | **Medido**: entre 0,42 y 0,50 de W_1 | §2.4 |

---

## §2. Doc 102

### §2.1 E101.1: tablas y ajustes

Confirmados [M]. Reproduzco los diez puntos (genérica N = 999982 en R1, R3, S11, R2, S13, S17; 2^e N = 2²⁰ en S11, R2, S13, S17) en todas las columnas impresas. Ejemplos:
- gen/S17: slack 7041,8; S_Cauchy 3,236; S_max 2175,8; B^C/B 0,689; Diag^C/B^C 1,533.
- 2^e/S17: slack 7218,7.

θ(genérica) = 0,988 ± 0,053, IC95% ±0,147, LOO [0,927; 1,071]; θ(2^e) = 1,036 ± 0,005, ±0,020, [1,032; 1,039]. Iguales al doc.

**Colapso N = 720720** [T + M]. Todo r ∈ H_N es primo o 1. Para r primo, M_r = φ(r)E_r porque E_1 = 0, luego B_M = S_BV y BV1* = B_M exactamente. Medido en R2, S11, S13 y S17: B_M/S_BV = 1,000000000000 y W_max = 422,02–1504,18.

### §2.2 E101.2

El número α_W(11..37) = 0,4281 [0,4243; 0,4318] se reproduce exactamente (ya en mi verificación del doc 101, G5).

La conclusión "α(z) se estabiliza en ≈ 0,42–0,43" queda **refutada**:
- El exponente asintótico es 1/2 [T, Teorema 101-CC].
- La medición hasta z = 709 lo muestra subiendo: α(|H|, z = 307..709) = 0,4659.

La conclusión que el doc extrae ("la vía de convertir el conteo en cota con α < 0,42 está muerta") se mantiene, y ahora por teorema.

### §2.3 §4: el margen de la ruta directa

La cadena aritmética del §4 usa dos cosas:
- **El "≈ 1,2·10⁵" del doc 101 §3.4**, que es erróneo: el valor correcto es ≈ 65 (`verificacion-101-claude-code.md` §2.10).
- **La ley |H| ≈ 0,033·N^{0,4208}**, que es preasintótica.

Con las mismas extrapolaciones y la cifra corregida: B_M/objetivo ≈ 65/slack ≈ 65/|H| ≈ 10⁻⁷ en 10²⁴. El margen sería **mayor** que los "10³–10⁴" del doc, no menor.

Sigue siendo una extrapolación de ajustes [M]: descansa en que S_BV/B_M ≈ cte a soporte creciente, que es el [P-rev] 102.2 y sigue abierto. La dirección de la conclusión no cambia.

### §2.4 [P-rev] 102.4: max_{d≥2} W_d

Medido [M]: max_{d≥2} W_d / W_1 = 0,42–0,46 en S11–S13 y 0,50 en S17, en las dos familias (0,018 en R3, donde solo hay primos). Una cota BV que use el perfil quitando d = 1, aprovechando que E_1 = 0, solo gana **un factor ≈ 2** frente a la de máximo. La pérdida de ~|H| del aplanamiento no se debe al término d = 1.

---

## §3. Doc 103

### §3.1 La derivación [T]

- **Paso 1.** Con S_r = Σ\*_χ χ̄(N)ψ(N,χ), Cauchy da |S_r|² ≤ φ\*(r)·Σ\*|ψ(N,χ)|², y φ\* ≤ φ. Por tanto (1/φ(r))|S_r|² ≤ Σ\*|ψ(N,χ)|². Correcto.
- **Paso 2.** Es la criba grande para caracteres primitivos vía sumas de Gauss: Σ\*_χ |Σa_nχ(n)|² ≤ (φ(r)/r)·Σ\*_a |S(a/r)|². Le sigue la desigualdad aditiva con espaciamiento δ, (N − 1 + δ⁻¹)·Σ|a_n|². El doc usa N + δ⁻¹, que es válida porque mayora. Correcto.
- **Sándwich.** δ ≥ 1/Q² por la separación de fracciones reducidas. δ ≤ 2/(r₂(r₂+2)) por el par exhibido 1/r₂, 1/(r₂+2). Ambas cotas son correctas.

### §3.2 Números

Confirmados [M]:
- Λ₂(999982) = 1,280280·10⁷ y Λ₂(2²⁰) = 1,347870·10⁷.
- Q_ℋ y r₂ (16345/2933, 37415/10703, 144445/81653).
- Los seis cocientes B̂_lo/B_M y B̂_hi/B_M (por ejemplo gen/S17: 1,532·10⁶ y 9,583·10⁶).
- La pendiente **2,3701 ± 0,0322, LOO [2,3518; 2,3950]**.

**δ_ℋ exacto en R1 y R3:**
- Coinciden δ = 9,137·10⁻⁶ y 8,5515·10⁻⁶, y los cocientes 1,72·10⁵ y 7,53·10⁴.
- El número de fracciones me sale 12 644 y 10 388, frente a 12 612 y 10 356 del doc. La diferencia es exactamente 32 = φ(5)+φ(7)+φ(11)+φ(13): el doc enumera solo r > ln N, que es el conjunto de B_M, aunque escriba Σ_{r∈ℋ}. Lo mismo ocurre en B̂₂ y B̂₃, que solo reproducen el doc sumando r > ln N. Conviene escribirlo así en el enunciado.

### §3.3 [P-rev] 103.2: δ_ℋ exacto en soportes grandes

El doc lo declara no computable. En S11 y S13 sí lo es [M]:
- **gen/S11** (1,64·10⁶ fracciones): δ⁻¹ = 4,833·10⁷, luego B̂₁/B_M = **1,876·10⁵**, que es **9,3× B̂_lo**.
- **gen/S13** (7,43·10⁶ fracciones): δ⁻¹ = 3,134·10⁸, luego B̂₁/B_M = **5,267·10⁵**, que es **5,4× B̂_lo**.

El extremo favorable del sándwich subestima el cociente verdadero mucho más que en R1/R3 (factor 1,07–1,08 allí). **El fracaso real de B̂₁ es todavía mayor** que el que usa el doc: el veredicto era conservador, como el doc afirmaba.

### §3.4 Desglose y B̂₂

Confirmados [M]:
- Tabla del §3 (por ejemplo gen/S17: B̂_lo^P/B_M^P = 1,188·10⁵ y B̂_lo^C/B_M^C = 2,224·10⁶).
- Diag^C/B_M^C con S_r: 1,495 / 1,488 / 1,533 / 1,674 / 1,615 / 1,576.
- **B̂₂**: |S_r| ≤ Σφ(d)|E_d| + |Δ_r| es trivial [T]. Cocientes 1,709–2,318, R1 = 2,934, R3 = 1,001, pendiente **0,1358 ± 0,0201, LOO [0,1198; 0,1484]**. C1 pasa y C2 falla, como dice el doc.

Matiz menor: el veredicto del §5.1 describe H_A como "delgado en cardinalidad (|H| ~ N^{0,42})". Asintóticamente |H| ~ N^{1/2} [T], así que es aún menos delgado. No cambia la conclusión del doc: B̂₁ muere por espaciamiento.

---

## §4. Doc 104

### §4.1 Definición y lecturas

**[P-rev] 104.2, resuelto contra el doc 94 íntegro.** En el doc 94, "parte nueva" es h = r/(r, A) en el §1 (l.54–58), que es la **lectura L1**. La partición que diagonaliza Q_s es la de T_s(h₀) en el §8 (l.314–336), por la parte lisa h₀ | P(z_A), que es la **lectura L2**. Las dos están en el texto en secciones distintas.
- Bajo L1, como todo r ∈ ℋ es coprimo con 3, B̂₃ = B̂₂ idénticamente [T; D104.4].
- La elección L2 del doc 104 corresponde al §8 del doc 94 y es la única no trivial. Confirmada como coherente con la fuente.

**Clases-par** [T]. Si s | r, s ∤ h₀, entonces μ(r/(h₀s)) = −μ(r/h₀) y φ(h₀s) = φ(h₀)(s−1). La suma de la clase es μ(r/h₀)φ(h₀)[E_{h₀} − (s−1)E_{h₀s}]. Correcto.

### §4.2 Números

Confirmados [M]:
- B̂₃/B_M = 1,679 / 1,876 / 2,181 / 1,798 / 2,000 / 2,262; R1 = 2,620.
- Pendiente **0,1326 ± 0,0195, IC95% ±0,0542, LOO [0,1168; 0,1441]**.
- Pendiente B̂₂/B̂₃ = 0,0032 ± 0,0018. La suma 0,0032 + 0,1326 = 0,1358 es exacta [T: MCO es lineal en la variable dependiente].
- ρ_intra = 0,987–0,992; clases por conductor de 1,69 a 2,30.

**[P-rev] 104.3, soportes mayores** [M]. Añado S19 (|H| = 8245) y S23 (|H| = 17999):
- B̂₃/B_M = 2,321 / 2,470 (genérica) y 2,368 / 2,512 (2^e). Siguen creciendo despacio.
- Pendiente con los 10 puntos: **0,1168 ± 0,0087, LOO [0,1074; 0,1214]**, sobre 0,1, así que C2 sigue fallando.
- Por familia: genérica 0,1261 [LOO 0,1162; 0,1381]; **2^e 0,1073 ± 0,0075, IC95% ±0,0238, LOO [0,0957; 0,1170], a caballo de 0,1**.

Lectura: la pendiente baja algo al ampliar el soporte, lo que es compatible con la saturación que propone el doc 104 §4.3, pero no cruza 0,1 en el ajuste conjunto. **El dictamen "C2 falla" se mantiene**, con la advertencia de que en la familia 2^e es estrecho.

---

## §5. Erratas y avisos para la autora (sin tocar sus archivos)

1. **Doc 102 §1 y §3:** retirar "α(z) se estabiliza en ≈ 0,42–0,43". El exponente es 1/2 [T].
2. **Doc 102 §4:** las cifras del margen heredan el 1,2·10⁵ erróneo del doc 101 §3.4. Corregidas, el margen sería mayor.
3. **Docs 103 y 104:** las sumas que se escriben Σ_{r∈ℋ} se calculan sobre r > ln N, que es el conjunto de B_M. Conviene escribirlo, porque explica la diferencia de 32 fracciones.
4. **Doc 103 [P-rev] 2:** δ_ℋ exacto en S11 y S13 es computable. El cociente real está 5–9× por encima de B̂_lo.
5. **Doc 104 [P-rev] 2 y 3:** resueltos (lectura L2 = doc 94 §8; S19/S23 medidos).

*Goldbach fuerte sigue abierta.*
