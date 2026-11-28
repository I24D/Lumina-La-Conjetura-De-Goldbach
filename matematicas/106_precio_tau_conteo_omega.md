# Documento 106 — Precio τ con soportes S19/S23 y conteo A105.1 de Ω_d

**Veredicto en una frase.** La pendiente de la pieza compuesta de B̂₅ **no** queda
certificada bajo 0,1 con los soportes grandes — al contrario: con diez puntos
sube a **0,0985** y su motor, el precio τ, queda con pendiente **0,1046**,
certificada *por encima* de 0,1 por LOO en la ventana ampliada —; pero el
conteo A105.1 muestra que ese precio es un efecto de conteo puro en su forma
correcta: no el promedio τ̄ (ley (ln z)^0,582) sino la razón sesgada
Ā = Σ(τ−1)²/φ / Σ(τ−1)/φ (ley (ln z)^1,142), que lo predice sin ψ con un
factor residual 1,13 que no crece; con ello la salvedad del doc 105 se cierra
en forma reformulada (lo que L105(b) necesita es un precio polilogarítmico, y
eso ahora es una ley de conteo exacta), el momento que queda por acotar queda
escrito explícitamente, y B_Δ resulta la pata elemental del lema.

Cadena: doc 100 (Δ_r y segundo factor) → doc 101 (conteo A100.1; BV muerta)
→ doc 102 (la pérdida es el aplanamiento) → doc 103 (criba restringida falla;
B̂₂ pasa C1) → doc 104 (B̂₃; línea de mayorantes absolutos agotada) → doc 105
(forma diagonal B̂₅; Lema L105 con salvedad de pendiente) → **este documento
ejecuta el trabajo que el doc 105 §6 dejó fijado**: experimento E105.1
(soportes S19/S23), conteo A105.1 (fórmulas para τ̄ y Ω_d), reevaluación del
Lema L105 y verificación de B_Δ/T(N).

Script: `106_precio_tau_conteo_omega.py` (solo stdlib; importa la maquinaria
verificada de los docs 97–105 sin modificarla y añade `build_H106`, variante
de `build_H` con factorización por división de prueba, necesaria porque el
SPF del doc 105 llega a 1,42·10⁶ < max ℋ(S29/S31/S37)). Salida:
`salida_doc106_local.txt`. Tiempo total: 14,9 s. Chequeos K0–K6: **PASS**
(K0: maquinaria — R2 |ℋ| = 880, S11 idéntico a Model97, `build_H106` idéntico
a `build_H` en z = 11, 13, 17; K1: reproducción exacta de los seis puntos y
las pendientes del doc 105; K4: mayorantes y pieza P exacta en los 12
puntos; K2/K3: §2; K5/K6: §1 y §4).

**Disciplina.** [T] = teorema demostrado aquí o en docs 85–105; [M] =
medición reproducible del script (no es prueba); [P-rev] = afirmación tomada
de la literatura o de un documento no releído por completo, pendiente de
revisión. Nada de lo medido se presenta como demostrado. No se usa GRH.
Goldbach sigue abierta. Este documento no modifica los docs 85–105.

---

## 1. Experimento E105.1: la pendiente de la pieza C con S19/S23

### 1.1. Los diez puntos [M]

Misma medida completa del doc 105 (con ψ vía Model97), soportes
z = 11, 13, 17, **19, 23** (A = 3), N = 999982 (gen) y N = 2²⁰ (2^e):

| punto | \|ℋ_N\| | B̂₅/B_M | pieza C (B̂₅^C/B_M^C) | D̂^C/B_M^{C,(M)} | Diag^C | precio τ^C | share_C |
|---|---|---|---|---|---|---|---|
| gen/S11 | 798 | 3,190 | 6,163 | 6,158 | 1,495 | 4,119 | 0,424 |
| gen/S13 | 1621 | 4,115 | 6,614 | 6,605 | 1,488 | 4,439 | 0,555 |
| gen/S17 | 5017 | 5,499 | 7,532 | 7,525 | 1,533 | 4,909 | 0,689 |
| gen/S19 | 8219 | 6,186 | 8,077 | 8,071 | 1,545 | 5,224 | 0,733 |
| gen/S23 | 17955 | 7,126 | 8,872 | 8,866 | 1,531 | 5,790 | 0,778 |
| 2^e/S11 | 802 | 3,600 | 7,014 | 7,007 | 1,674 | 4,185 | 0,432 |
| 2^e/S13 | 1628 | 4,470 | 7,411 | 7,403 | 1,616 | 4,582 | 0,541 |
| 2^e/S17 | 5033 | 5,619 | 7,667 | 7,661 | 1,576 | 4,862 | 0,693 |
| 2^e/S19 | 8245 | 6,351 | 8,285 | 8,279 | 1,568 | 5,280 | 0,734 |
| 2^e/S23 | 17999 | 7,354 | 9,118 | 9,113 | 1,552 | 5,872 | 0,783 |

Los seis primeros valores de cada columna reproducen el doc 105 (K1). C1
sobrevive a la extensión: B̂₅/B_M ∈ [3,190, 7,354] ≤ 10 en los diez puntos.

### 1.2. Pendientes: seis puntos contra diez puntos [M]

Regresión log-log contra ln|ℋ_N| (protocolo de los docs 102–105: SE, IC95%
con t de Student, LOO):

| cantidad | 6 puntos (doc 105) | 10 puntos (E105.1) |
|---|---|---|
| agregado B̂₅/B_M | 0,2658 ± 0,0311; LOO [0,235, 0,282] | 0,2423 ± 0,0136; LOO [0,225, 0,250] |
| **pieza C, forma S** | 0,0784; LOO [0,0517, 0,1028] | **0,0985 ± 0,0133; IC95% ± 0,0307; LOO [0,0883, 0,1121]** |
| **pieza C, forma M (D̂^C)** | 0,0785; LOO [0,0518, 0,1028] | **0,0986 ± 0,0133; IC95% ± 0,0308; LOO [0,0885, 0,1123]** |
| precio τ^C | 0,0870 ± 0,0093; LOO [0,0807, 0,0923] | **0,1046 ± 0,0056; IC95% ± 0,0130; LOO [0,1005, 0,1070]** |
| factor Diag^C | −0,0085; LOO [−0,0288, 0,0159] | −0,0060 ± 0,0103; LOO [−0,0152, 0,0056] |

### 1.3. Respuesta a la pregunta (a): NO

**¿La pendiente de la pieza C queda certificada intervalarmente bajo 0,1
con los soportes grandes? No.** Y la dirección del movimiento es la
contraria a la que el doc 105 dejaba abierta como esperanza:

1. La estimación puntual *sube*: 0,0785 (seis puntos) → **0,0985/0,0986**
   (diez puntos), pegada al umbral. Los intervalos no certifican por ningún
   criterio: LOO máximo 0,1121/0,1123 > 0,1; b + IC95% = 0,1292/0,1294 > 0,1.
2. El motor está peor que la pieza: la pendiente del **precio τ** con diez
   puntos es **0,1046**, con LOO **[0,1005, 0,1070] entero por encima de
   0,1** — en el sentido LOO usado en los docs 104–105, la pendiente del
   precio queda certificada *sobre* el umbral (su IC95% [0,0916, 0,1176]
   aún roza el umbral por abajo; se declara con esa precisión). El factor
   Diag^C sigue plano y su intervalo se estrecha en torno a 0 (−0,0060,
   LOO ⊂ [−0,0152, 0,0056]): toda la pendiente de la pieza es el precio,
   como ya decía el doc 105, pero el precio crece *más* de lo que la
   ventana S11–S17 sugería.
3. Las pendientes locales del precio entre soportes gen consecutivos no
   desaceleran en la ventana: 0,106 (11→13), 0,089 (13→17), 0,126 (17→19),
   0,132 (19→23). El 0,0785 del doc 105 era un artefacto de ventana corta:
   los puntos S19/S23 quedan *sobre* la extrapolación de la tendencia
   S11–S17.

**Consecuencia para el criterio C2-por-piezas:** tal como fue formulado en
el doc 103 (pendiente contra |ℋ_N| a N ≈ 10⁶ fijo ≤ 0,1), B̂₅ **no lo
pasa**, y más soportes del mismo tipo no lo rescatan: a N fijo la familia
está agotada en z = 37 y el conteo de §2 muestra que la ley subyacente no
produce pendientes decrecientes en esta ventana. La §3 explica por qué este
fracaso del criterio no es un fracaso del Lema L105: el criterio medía el
proxy equivocado.

---

## 2. Conteo A105.1: fórmulas exactas para τ̄^C y Ω_d^C

### 2.1. Derivación [T]

Se trabaja con la descomposición por ramas del doc 101 (verificada
exhaustivamente allí, J1): con 𝒱₀ = V⁺(D/s_min) la rama del primo cambiante
mínimo y 𝒱_s = V⁺(D/s) la de cada s ∈ 𝒮 (primos cambiantes coprimos con A),

ℋ = 𝒱₀ ⊔ ⨆_{s∈𝒮} s·𝒱_s,   ℋ_C = (𝒱₀∖{1}) ⊔ ⨆_{s∈𝒮} (s·𝒱_s ∖ {s}),

disjunta. Para un nivel L (conjunto de primos distintos, cerrado por
divisores, con τ(h) constante en h si y solo si |L| = ω(h)), se definen los
conteos de nivel, computables sin enumerar ℋ (recursión sobre los primos
permitidos, coste ≪ |𝒱|):

C(L) = |{h : supp(h) = L}|,   Φ(L) = Σ_{h : supp(h)=L} 1/φ(h) = ∏_{p∈L} 1/(p−1).

(Φ(L) es producto porque 1/φ es multiplicativa y sobre un nivel fijo el
conjunto de h es un producto cartesiano de potencias p^j ≥ 1 acotadas por
D/s; con cotas distintas por primo se usa la forma general
Φ(L) = ∏_{p∈L} (Σ_j p^{−j}·p/(p−1)) evaluada por la misma recursión — el
script implementa la recursión general, no el producto simplificado.)
Como τ(h) = 2^{|L|} sobre el nivel L, todo promedio de potencias de τ−1
sobre una rama es una suma sobre niveles. En particular, con
𝒩(𝒱) = Σ_{h∈𝒱, h compuesto} (τ(h)−1)/φ(h) y 𝒟(𝒱) = Σ_{h∈𝒱, h compuesto} 1/φ(h):

**τ̄^C = [𝒩(𝒱₀) + Σ_{s∈𝒮′} 𝒩(s𝒱_s)] / [𝒟(𝒱₀) + Σ_{s∈𝒮′} 𝒟(s𝒱_s)]**   (106.1)

donde 𝒮′ excluye los s | N (doc 96/100: en ℋ_N esos s desaparecen) y, para
la rama s𝒱_s, τ(sh) = 2τ(h) si s ∤ h y τ(h) si s | h, con 1/φ(sh) análogo:
la recursión trata s como un primo más con su cota. Esto da τ̄^C **sin ψ y
sin enumerar ℋ**.

Para Ω_d^C = Σ_{r∈ℋ_C, d|r} (τ(r)−1)/φ(r) y U_d^C = Σ_{r∈ℋ_C, d|r} 1/φ(r),
se reindexa por niveles: d | r con supp(r) = L ⊇ supp(d) se escribe
r = d·u con u recorrido por los cocientes del nivel; agrupando por
(d₀, M) con d₀ = parte suave (primos < z) de d y M = supp(d)∖suave:

**Ω_d^C = Σ_{L ⊇ M} (2^{|L|} − 1)·Q(d₀, L),   U_d^C = Σ_{L ⊇ M} Q(d₀, L)**   (106.2)

con Q(d₀, L) = Σ_{u : d₀u ∈ ℋ_C, supp(d₀u) = L} 1/φ(d₀u), suma sobre los
cocientes compatibles con las ramas (u recorre, por rama, los elementos de
𝒱 o s𝒱_s divisibles por d₀/s según el caso, con (u, d₀) arbitrario — la
recursión del script enumera los cocientes por niveles, nunca ℋ entero).
Si d tiene un primo ≥ z que no es un s válido, o d ∤ ningún r ∈ ℋ_C,
Ω_d^C = U_d^C = 0. Las dos identidades de cierre:

**Σ_d Ω_d^C = Σ_{r∈ℋ_C} (τ(r)−1)²/φ(r),   Σ_d U_d^C = Σ_{r∈ℋ_C} (τ(r)−1)/φ(r)**   (106.3)

(intercambio de sumas: cada r aporta a sus τ(r) divisores). Se define la
razón sesgada de conteo

**Ā := Σ_d Ω_d^C / Σ_d U_d^C = Σ_r (τ(r)−1)²/φ(r) / Σ_r (τ(r)−1)/φ(r)**   (106.4)

— el promedio de (τ(r)−1) bajo la medida (τ(r)−1)/φ(r) sobre ℋ_C, en vez de
bajo 1/φ(r) (que da τ̄^C).

### 2.2. Verificación contra la enumeración [M]

- **K3:** τ̄^C por (106.1) contra enumeración directa de ℋ_C: error
  relativo ≤ 3,5·10⁻¹⁴ en R1 y S11–S37 (ambos N donde aplica), incluida la
  identidad Σ_d Ω_d^C de (106.3). Valores (familia gen): τ̄^C = 3,3078
  (z=11), 3,4770 (13), 3,6500 (17), 3,7622 (19), 3,8912 (23), 4,0348 (29),
  4,1044 (31), 4,2136 (37).
- **K2:** Ω_d^C y U_d^C por (106.2) contra enumeración, **para todo d**,
  en R1, R3 y S11–S29, ambos N: peor error relativo ≤ 1,6·10⁻¹⁴
  (|𝒟| hasta 47 177 en S29). R3 (solo primos): 𝒟 vacío, fórmulas vacuamente
  exactas.
- El coste es el de la medida sin ψ: S37 completo en 2,2 s. **El precio τ
  queda así predecible sin ψ en soportes donde la medida completa es
  inviable** (S29/S31/S37 no tienen medida con ψ en este doc: se declara).

### 2.3. Las leyes de crecimiento [M]

Ajustes sobre la familia gen (z = 11…37 para τ̄ y Ā; diez puntos para el
precio):

- τ̄^C ~ (ln z)^γ con **γ = 0,582 ± 0,027**, LOO [0,567, 0,591]. (Ajuste
  alternativo z^β: β = 0,196 ± 0,017; en esta ventana ambas leyes ajustan,
  pero la forma (ln z)^γ es la que tiene mecanismo de conteo: productos
  sobre primos < z de factores (1 + c/(p−1)), tipo Mertens.) La pendiente
  local de τ̄ contra ln|ℋ_N| es 0,0464 (LOO [0,0439, 0,0476]), consistente
  con γ/(ln z · dln|ℋ|/dln z) ≈ 0,58/(3,0·4,2) ≈ 0,046.
- Ā ~ (ln z)^{γ_Ā} con **γ_Ā = 1,142 ± 0,020**, LOO [1,116, 1,158]
  (Ā = 3,651, 4,020, 4,424, 4,691, 5,012, 5,396, 5,566, 5,865 para
  z = 11…37). Numéricamente γ_Ā ≈ 2γ(τ̄) dentro del error de ajuste
  (2·0,582 = 1,164 frente a 1,142 ± 0,020): se registra como observación,
  no como identidad.
- Comparación con el precio medido (con ψ), diez puntos:
  precio/τ̄^C crece de 1,245 a 1,508 (pendiente de ln(precio/τ̄) =
  **0,0532**, LOO [0,0487, 0,0563], entero sobre 0);
  **precio/Ā ∈ [1,098, 1,169] con pendiente de ln(precio/Ā) = 0,0042 ±
  0,0060, LOO [−0,0007, 0,0082]** — plano. Predicción sin ψ:
  precio ≈ 1,13·Ā reproduce los diez precios medidos con error ≤ 3,5%.

### 2.4. Respuesta a la pregunta (b)

**El crecimiento del precio τ es un efecto de conteo puro, pero la
estadística de conteo correcta no es el promedio τ̄: es la razón sesgada
Ā.** El mecanismo es visible en las definiciones: el precio es
Σ_d F_d²Ω_d / Σ_d F_d²U_d, y Ω_d/U_d es ya, para cada d, un promedio de
(τ(r)−1) sobre los múltiplos de d — los conductores con τ grande pesan en
Ω_d en proporción a su propio (τ−1). Ese sesgo por tamaño es exactamente
Ā = E[(τ−1)²]/E[τ−1] bajo la medida 1/φ, y explica la diferencia de
exponentes: τ̄ ~ (ln z)^0,582 (la restricción 𝒱⁺ rebaja el promedio —
el sobre ingenuo de todos los enteros daría exponente 1) mientras
Ā ~ (ln z)^1,142 conserva el crecimiento del segundo momento. Lo que no es
conteo de primer orden — la correlación de F_d² con Ω_d/U_d — cuesta un
factor 1,13 ± 0,03 **que no crece** en la ventana (pendiente 0,004 ± 0,006).
En síntesis: precio τ = 1,13 · Ā(z) · (1 + o(1)) con Ā dado por la fórmula
exacta (106.4) y ley (ln z)^1,142: **polilogarítmico en la variable del
marco**, no una potencia de |ℋ| con exponente fijo. La pendiente 0,1046 de
§1 es la sombra de esa ley polilogarítmica proyectada sobre ln|ℋ_N| a
N fijo: d ln precio/d ln|ℋ| ≈ γ_Ā/(ln z · d ln|ℋ|/d ln z), que vale ~0,10
en la ventana y decae como 1/ln z.

---

## 3. Reevaluación del Lema L105

### 3.1. La salvedad de la pendiente: cerrada, en forma reformulada

El doc 105 dejó L105 con la salvedad de que la pendiente de la pieza C
(0,0785) no estaba certificada bajo 0,1. Tras E105.1 + A105.1 el estado es:

- **Como criterio numérico C2-por-piezas, la salvedad se cierra en
  negativo:** la pendiente no baja de 0,1 con soportes mayores; sube a
  0,0985 y su motor queda LOO-certificado sobre 0,1 (§1.3). Quien exija
  C2 tal cual debe declarar B̂₅ no alcanzable *por ese criterio*.
- **Como condición sobre el Lema L105, la salvedad se disuelve:** L105(b)
  pide D̂^C ≤ C_Ω·T(N) con C_Ω **absoluta**. El criterio C2 era un filtro
  de alcanzabilidad (doc 103 §2): detecta si el precio crece como potencia
  de |ℋ| — lo que rompería la cota al escalar el soporte — o más despacio.
  El conteo de §2 decide esa pregunta en la forma que importa: el precio
  es 1,13·Ā ≍ (ln z)^1,142, **polilogarítmico en z = N^{1/10} del marco**,
  i.e. ≪ (ln N)^{1,15}: absorbible por el ln⁴N del denominador de T(N) y
  por el margen medido del segundo factor (doc 102: B_M quedaría 10³–10⁴
  veces bajo el objetivo a N = 10²⁴ bajo las tasas medidas; B̂₅ cuesta
  ≤ 9,2 en ventana). La pendiente contra ln|ℋ| a N = 10⁶ fijo mezcla dos
  escalas (z crece, N no) y no es la cantidad que L105(b) necesita. La
  condición que sustituye a la salvedad es explícita y de conteo:
  **ρ_C y W* polilogarítmicos** (§3.2), con leyes medidas
  (ln z)^1,142 para el promedio sesgado y valores ρ_C = max_d Ω_d^C/U_d^C
  = 7, 15, 31, 31, 31 (z = 11…23) [M].

### 3.2. El momento exacto que queda por demostrar

Sustituyendo Ω_d^C por su fórmula de conteo (106.2), D̂^C ≤ C_Ω·T(N) es la
desigualdad concreta

**ℳ_Ω(N) := Σ_{d∈𝒟} Ω_d^C · φ(d)² · E(N; d, N)² ≤ C_Ω · T(N)**   (106.5)

con 𝒟 = {d ≥ 2 : d | r para algún r ∈ ℋ} (el conjunto de divisores del
soporte, |𝒟| medido: 577 en S11 hasta 47 177 en S29 [M]), Ω_d^C dado por
(106.2) — una función de conteo explícita, sin ψ — y
T(N) = N^{1,920232}/(0,0225·ln⁴N) (doc 100 §5.3, heredado por L105).
Reducciones disponibles [T]:

- Ω_d^C/U_d^C es un promedio ponderado de (τ(r)−1) sobre los múltiplos de
  d en ℋ_C, luego Ω_d^C ≤ ρ_C·U_d^C con ρ_C ≤ max_{r∈ℋ_C} (τ(r)−1), y
  D̂^C ≤ ρ_C · Σ_d U_d^C F_d² (la misma forma con los pesos U_d^C de
  (106.2)).
- Parseval por conductor: φ(d)²E(N;d,N)² ≤ φ(d)·Σ*_{χ mod d} |S(χ)|², así
  ℳ_Ω ≤ W* · Σ_{d∈𝒟} Σ*_{χ mod d} |S(χ)|² con
  **W* := max_d Ω_d^C·φ(d)** — medido: 9,52 (z=11), 15,00 (13), 31,00
  (17, 19, 23) [M]; W* es una cantidad de conteo (Ω_d^C es (106.2) y
  φ(d) aritmética elemental).

**¿Lo acotan la gran criba clásica o BDH? No en la fuerza necesaria, y la
cuenta es explícita.** La gran criba sobre el bloque de caracteres da
Σ_{d≤D_A} Σ*_χ |S(χ)|² ≤ (N + D_A²)·Σ_n Λ(n)² ≈ 2N² ln N a D_A = N^{1/2},
luego ℳ_Ω ≤ 2W*·N² ln N: sobrepasa el objetivo C_Ω·T(N) por el factor
0,045·W*·**N^{2c}**·ln⁵N con 2c = 0,079768 — exactamente el ahorro N^{2c}
que define el marco (docs 89/100). La forma clásica de BDH [P-rev: cita
estándar, no reprobada aquí] acota la varianza total
Σ_{q≤Q} Σ_a |E(x;q,a)|² ≪ Qx log x + x²/(log x)^A para
Q ≤ x/(log x)^B: a nivel Q = N^{1/2} el término dominante es
N²/(log N)^A — un ahorro **polilogarítmico** sobre N², no el ahorro
**polinómico** N^{2c} que (106.5) exige. Y ninguna de las dos afirmaciones
clásicas cubre (106.5) literalmente: pesos Ω_d^C φ(d) no constantes,
soporte 𝒟 delgado y de forma especial (|𝒟| ~ N^{0,42}, divisores del
soporte cribado, no todos los módulos ≤ D_A), y una sola clase
a ≡ N (mod d) en vez de la suma sobre clases.

**Lo que queda por demostrar, exactamente:** una cota de varianza de tipo
BV/BDH para el momento ponderado (106.5) con ahorro polinómico N^{2c}
sobre la escala trivial N² a nivel D_A = N^{1/2} — o una vía que explote
lo que los enunciados clásicos no usan: la delgadez y la estructura de
conteo de 𝒟 y de Ω_d^C (ahora explícitos por (106.2)). La parte de conteo
de L105(b) está cerrada por §2 [T]+[M]; la parte abierta es puramente de
distribución de ψ(N;·) en la clase única sobre los módulos de 𝒟. (Nota de
no-confusión: la ruta BV del doc 100 era de *primer* momento, vía W_max;
la de L105 es de *segundo* momento; BV no se aplica a (106.5) ni BDH
clásico basta, por la cuenta anterior.)

### 3.3. Estado del Lema L105 tras el doc 106

| pata | estado |
|---|---|
| B_M ≤ B̂₅ | [T] (doc 105) |
| pesos Ω_d^C, U_d^C | fórmulas exactas de conteo (106.2) [T], verificadas ≤ 1,6·10⁻¹⁴ [M] |
| precio τ | = 1,13·Ā, Ā exacto (106.4), ley (ln z)^1,142 [M]; polilog en el marco |
| D̂^C ≤ C_Ω·T(N) | reducida a (106.5); gran criba corta por N^{2c}ln⁵; BDH clásico corto en su término N²; **abierta** |
| B_Δ ≤ C_Δ·T(N) | §4: numéricamente ≤ 2,2·10⁻⁴·T; prueba elemental asemblable [T] |

---

## 4. B_Δ ≤ C_Δ·T(N): la pata elemental

### 4.1. Medición [M]

B_Δ/T(N) en los doce puntos (T(999982) = 4,0526·10⁸; T(2²⁰) = 4,3787·10⁸):
de 2,3·10⁻⁷ (R3) y 4,4·10⁻⁷ (R1) a **2,16·10⁻⁴** (gen/S23, el máximo).
El cociente **no es plano** a N fijo — crece con el soporte, aproximadamente
como |ℋ| (×43,5 de S11 a S23 gen frente a ×22,5 de |ℋ|) — porque T(N) no
depende del soporte y B_Δ sí; pero está **acotado con margen enorme**:
C_Δ = 10⁻³ cubre toda la ventana con factor ≥ 4,6, y B_Δ ~ 10⁻⁶·B_M
(doc 105) con B_M a su vez 10³–10⁴ bajo el objetivo en la extrapolación
del doc 102.

### 4.2. La prueba esperada es elemental [T]

Por la forma en caracteres de Δ_r (docs 97/98, verificada por dos rutas en
el doc 100), Δ_r = Σ*_χ χ̄(N)(ψ(N,χ) − ψ(N,χ*)) y
ψ(N,χ) − ψ(N,χ*) = −Σ_{p|r, p∤cond χ*} Σ_{j≥1, p^j≤N} χ*(p)^j log p, de
donde, acotando cada primo por Σ_j log p ≤ ln N y sumando sobre ≤ φ(r)
caracteres:

**|Δ_r| ≤ φ(r)·ω(r)·ln N**   (106.6)

y por tanto B_Δ ≤ ln²N · Σ_{r∈ℋ} φ(r)ω(r)². El factor de conteo se acota
con el material del doc 101: en el marco ℋ ⊂ [1, D_A] por construcción de
los pesos, |ℋ| tiene la fórmula exacta de A100.1 (|ℋ| ~ N^{0,4208}), y
ω(r) ≤ ln N/ln 2, luego Σ_{r∈ℋ} φ(r)ω(r)² ≤ ω_max²·|ℋ|·D_A ≪
polilog·N^{0,9208} y

B_Δ ≪ polilog(N)·N^{0,9208} ≪ T(N) = N^{1,920232}/(0,0225 ln⁴N),

con margen ~N^{1,0}/polilog. **Ningún ingrediente de distribución**: ni
BV, ni BDH, ni criba grande — solo la definición de Δ_r como suma sobre
potencias de primos y conteo del soporte. Se declara con precisión el
matiz numérico: las cotas crudas (106.6) evaluadas en la ventana dan
cota/T entre 3,4 (gen/S11) y 9,5·10³ (gen/S23) — sobrepasan T(N) *en la
ventana* porque la ventana evalúa T en el N del modelo (10⁶) mientras el
soporte S23 está escalado para el N del marco (z^{10} ≈ 4·10¹³), y porque
(106.6) desperdicia la cancelación de Δ_r (la cota cruda es ~10⁷–10⁸
veces el B_Δ medido). El exceso es un artefacto de ventana, no una
obstrucción: en el marco, la misma cota cierra con un N entero de margen.

**Respuesta a (d):** B_Δ/T(N) está acotado en ventana (≤ 2,2·10⁻⁴, no
plano a N fijo por el artefacto de escala) y su prueba esperada es la
elemental: cota puntual (106.6) + conteo A100.1. Es la pata fácil de L105.

---

## 5. Qué NO se ha demostrado

- L105(b) sigue **abierto**: (106.5) no está probado; §3.2 muestra que ni
  la gran criba ni el BDH clásico lo dan en la fuerza necesaria. Lo nuevo
  es la *forma exacta* del hueco (ahorro N^{2c} en un momento ponderado
  explícito sobre un soporte delgado con pesos de conteo cerrados), no su
  cierre.
- El criterio C2-por-piezas del doc 103 queda **no satisfecho** por B̂₅ en
  su formulación original (§1.3); §3.1 argumenta que el criterio medía un
  proxy a N fijo y que la cantidad relevante para L105(b) — precio
  polilogarítmico — sí está establecida como ley de conteo medida. Esta
  reformulación es una *interpretación* del filtro, declarada como tal:
  quien mantenga C2 literal debe registrar B̂₅ como no alcanzable por C2
  aunque L105(b) no lo requiera.
- Las leyes (ln z)^0,582 y (ln z)^1,142 son ajustes [M] en z ∈ [11, 37],
  no asintóticas probadas; la predicción precio ≈ 1,13·Ā está verificada
  solo en los diez puntos con ψ disponibles.
- Nada de esto demuestra Goldbach: aun con L105 completo faltaría el
  puente del doc 94/96 (problema independiente, reaperturado en el doc 96)
  y el resto de la cadena de Chen en el marco Q-twisted.

## 6. Archivos

- `106_precio_tau_conteo_omega.md` — este documento.
- `106_precio_tau_conteo_omega.py` — script (K0–K6 PASS, 14,9 s).
- `salida_doc106_local.txt` — salida completa de la corrida local.
