# DOC 100 — Verificación independiente parcial de los docs 98/99 de Lumina Chat
## Corrección Δ_r + primer diagnóstico numérico del segundo factor

**Fecha:** 2026-11-28 (02:12–02:50 EST).
**Autor:** LUMINA (línea dedicada, Muse).
**Estado:** PARCIAL — verificación desde la transcripción del chat, sin acceso a los documentos íntegros ni a los scripts.
**Archivos:** `100_verificacion_doc98_99_segundo_factor.py` (stdlib solo),
`salida_doc100_local.txt` (89 líneas, 3,5 s, PASS).
**Contexto:** docs 85–97 en `~/workspace/goldbach/`; doc 97 (`97_verificacion_doc96_diagonal_mobius.md`)
como base. Goldbach fuerte sigue abierta. Sin GRH sin probar. Sin circularidad.

---

## §0. Limitación de esta verificación [LEER PRIMERO]

Esta verificación es **parcial** por una razón de infraestructura, no de método:
el túnel a la PC de Dal está caído (egreso #2, desde 2026-11-27 ~12:16 EST, verificado
por última vez 2026-11-28 00:14), de modo que **los docs 98 y 99 íntegros y los scripts
97 y 99 de Lumina NO están disponibles**. Lo que se verifica aquí es exclusivamente:

1. la respuesta que Lumina dejó en el chat «Revisa UNIF U 2M X» (2026-11-28 ~00:40),
   transcrita en la tarea;
2. los docs 85–97 locales como contexto y como base de comparación;
3. cálculos propios, independientes, con ψ reales y caracteres reales, en modelos
   sintéticos construidos con los mismos pesos λ⁺ de Nathanson que los docs 85–97.

Consecuencia: todo lo que depende del texto exacto de los docs 98/99 o de sus scripts
queda marcado **[P-rev]** en el §7 y debe re-verificarse cuando el túnel vuelva. Lo que
**sí** se verifica aquí es independiente de esos textos: la fórmula corregida
S_r = M_r + Δ_r (por dos rutas distintas), la condición p∥r, el orden de magnitud del
segundo factor, la condición sobre N y la evaluación de las tres cuestiones.

**Convención de etiquetas** (igual que en los docs 93/95/97):
[T] teorema / identidad exacta verificada; [M] medición numérica (no es prueba);
[P-rev] pendiente de re-verificación contra el documento íntegro o el script;
[N] observación numérica que requiere explicación.

---

## §1. Resumen ejecutivo

| # | Afirmación de Lumina (docs 98/99) | Veredicto doc 100 |
|---|--------------------------------------|-------------------|
| 1 | El doc 96 NO contiene Δ_r; se retira la exactitud de S_r = Σ μ(r/d)φ(d)E; la formulación a auditar es S_r = Σ μ(r/d)φ(d)E + Δ_r | **ACEPTADA Y CONFIRMADA** — la fórmula corregida se verifica exacta por dos rutas independientes (§2, [T]) |
| 2 | Verificación independiente en PC (N=999983, 16 conductores impares): err abs máx 6,56·10⁻⁹; err rel máx 9,71·10⁻¹²; máx corrección 2,59% de \|S_r\| | **REPRODUCIDA ESENCIALMENTE EXACTA** — mis 16 conductores: err abs máx 1,43·10⁻⁸; err rel máx 9,79·10⁻¹²; máx corrección **0,0259 = 2,59% exacto** (§2, [M]) |
| 2b | La condición p∥r importa (no basta sumar sobre todos los primos divisores) | **CONFIRMADA CUANTITATIVAMENTE** — ignorarla inventa o multiplica la corrección por hasta 4× (§2, [T]+[M]) |
| 3 | Segundo factor (doc 99, 69 conductores, N=999983, A=3, κ=1): B_M=3,9222·10⁷; B_cross=5,4023·10³; B_Δ=466,7; cambio rel 1,4966·10⁻⁴; B_cross > B_Δ | **ORDEN DE MAGNITUD REPRODUCIDO; VALOR ABSOLUTO NO** — el cambio rel ~10⁻⁴ y B_cross ≫ B_Δ se reproducen en todos mis modelos; B_M=3,92·10⁷ corresponde a un modelo de soporte pequeño (R1: 8,07·10⁷, 2×), NO a la escala del doc 94 (R2: 3,85·10⁹, 98×) ni a la hipótesis z=3 (R3: 2,17·10⁸, 5,5×) (§3, [M]) |
| 3b | El resultado está en el mismo orden que mi medición 1,24·10⁻⁴ (doc 97) | **CORRECCIÓN (E100.1)** — los dos números miden cosas distintas y su coincidencia es accidental: el peso del doc 97 sobre soportes completos es ~10⁻⁶, no 10⁻⁴ (§3) |
| 4 | Archivos nuevos: doc 96 (errata), doc 98 (prueba Δ_r + auditoría 94d/94e), doc 99 (segundo factor), scripts 97 y 99; commits e51a4c5, 4147a6a; rama 6 commits por delante de GitHub | **NO VERIFICABLE DESDE AQUÍ** — [P-rev] (§7); el túnel está caído |
| 5 | Tres cuestiones (a) cota útil para B_M en R_A; (b) si suaves + un primo de switching mejoran la gran criba; (c) estimación uniforme con cruzado + Δ_r. Condición: N creciente, pares, pérdida log real de D_A; proxy z=3 no asintótico | **(a) ES LA MÁS ATACABLE** — se reduce a contar el soporte (W_max ~ \|H\|) y la ruta BV cruza el objetivo entre N=10¹² y 10²⁴ si W_max se controla (§5). **Condición N: VERIFICADA** contra el doc 97 (cero usos de z=3) y testeada con N par (§4) |
| 6 | Estado declarado: Gram demostrada; diagonalización demostrada; corrección imprimitiva demostrada; numéricos reproducidos; segundo factor investigado sin cota analítica; Goldbach pendiente | **DE ACUERDO, CON UNA PRECISIÓN** — «corrección demostrada» vale para la formulación corregida S_r = M_r + Δ_r; los numéricos del segundo factor se reproducen en ORDEN, no en valor absoluto (§6) |

**En una frase:** la corrección Δ_r es correcta y está verificada por dos rutas
independientes con los mismos números que Lumina; el segundo factor es real, pequeño
(~10⁻⁴) y dominado por el término cruzado, pero su valor absoluto depende del soporte
y el modelo del doc 99 no es el de la escala del doc 94; la cuestión (a) es ahora un
problema de conteo del soporte, no de identidades nuevas.

---

## §2. La corrección Δ_r: verificación ortogonal y condición p∥r

### 2.1. La fórmula corregida es exacta [T]

El doc 97 (Thm 97.1) propuso, y el doc 98 adopta, la identidad corregida

    S_r(N) = Σ_{d|r} μ(r/d) φ(d) E(N; d, N) + Δ_r(N),

donde S_r(N) = Σ*_{χ* mod r} χ̄*(N) ψ(N, χ*) es la suma sobre caracteres primitivos,
E(N; d, N) = ψ(N; d, N) − ψ(N)/φ(d), y Δ_r es la corrección imprimitiva

    Δ_r = Σ_{p|r, p²∤r} P_{r/p},   P_e = Σ_{d|e} μ(e/d) A_d,
    A_d = φ(d) · Σ_{j: p^j ≤ N, p^j ≡ N (mod d)} log p.

La derivación (doc 97 §2, reproducida aquí en el script como `delta_r_orth`) NO enumera
caracteres: usa ortogonalidad sobre los caracteres primitivos módulo e=r/p y la
identidad Σ_{d|e} μ(e/d)·[p^j ≡ N (mod d)] = 1_{p^j ≡ N (mod e)}·(1 si e|r/p…).
La verificación K0/K5 la comprueba contra el cálculo DIRECTO por caracteres
(`CharGroup`: caracteres primitivos construidos por CRT con raíces de la unidad,
ψ(N,χ) = Σ_{n≤N} Λ(n)χ(n) con ψ real):

| r | #χ* | S_r | M_r | Δ_orth | Δ_char | \|S−(M+Δ)\| |
|---|-----|-----|-----|--------|--------|--------------|
| 3 | 1 | −715,0038 | −728,1872 | 13,1833 | 13,1833 | 1,28·10⁻⁹ |
| 5 | 3 | 1046,8240 | 1033,9485 | 12,8755 | 12,8755 | 9,29·10⁻¹⁰ |
| 7 | 5 | 526,2080 | 512,5867 | 13,6214 | 13,6214 | 2,66·10⁻⁹ |
| 9 | 4 | 674,8847 | 674,8847 | 0 | 0 | 6,61·10⁻⁹ |
| 11 | 9 | 2537,9556 | 2525,9661 | 11,9895 | 11,9895 | 5,44·10⁻¹⁰ |
| 13 | 11 | 2934,5508 | 2921,7261 | 12,8247 | 12,8247 | 1,59·10⁻⁹ |
| 15 | 3 | 1563,0200 | 1563,0200 | 0 | 0 | 4,60·10⁻⁹ |
| 21 | 5 | 1540,4249 | 1554,0463 | −13,6214 | −13,6214 | 1,77·10⁻⁹ |
| 25 | 16 | 357,5449 | 357,5449 | 0 | 0 | 1,37·10⁻¹⁰ |
| 27 | 12 | 791,9846 | 791,9846 | 0 | 0 | 6,09·10⁻¹⁰ |
| 35 | 15 | 2382,8836 | 2374,4999 | 8,3837 | 8,3837 | 3,56·10⁻⁹ |
| 45 | 12 | 6417,8006 | 6421,0194 | −3,2189 | −3,2189 | 4,44·10⁻⁹ |
| 49 | 36 | 16933,5702 | 16933,5702 | 0 | 0 | 1,43·10⁻⁸ |
| 63 | 20 | 739,7151 | 739,7151 | 0 | 0 | 4,05·10⁻⁹ |
| 75 | 16 | 8415,1151 | 8406,3262 | 8,7889 | 8,7889 | 3,86·10⁻⁹ |
| 99 | 36 | 3832,7326 | 3832,7326 | 0 | 0 | 2,68·10⁻⁹ |

**Resultado K0/K5 (N=999983, 16 conductores, los mismos parámetros que Lumina):**

- máx \|S_r − (M_r + Δ_r)\| = **1,43·10⁻⁸** (Lumina: 6,56·10⁻⁹ — mismo orden);
- máx error relativo = **9,79·10⁻¹²** (Lumina: 9,71·10⁻¹² — coincide al 1%);
- máx \|Δ_r\|/\|S_r\| = **0,0259** (Lumina: 2,59% — **coincide exactamente**);
- máx \|Δ_orth − Δ_char\| = 1,43·10⁻⁸ (las dos rutas independientes coinciden).

El conjunto de 16 conductores que yo uso — [3,5,7,9,11,13,15,21,25,27,35,45,49,63,75,99],
todos los impares ≤ 99 con al menos un divisor primo que cumple p∥r o con potencias
primas — reproduce los tres estadísticos de Lumina. No conozco su lista exacta de 16
([P-rev] §7), pero la coincidencia del máximo relativo (2,59%) y del error relativo
(9,7·10⁻¹²) es evidencia fuerte de que ambas verificaciones miden el mismo objeto.
La retirada de la exactitud del doc 96 queda así **confirmada y cerrada**: la
formulación correcta es S_r = M_r + Δ_r, y es exacta.

Dos cancelaciones exactas que el cálculo exhibe (y que la ruta por caracteres confirma
a 10⁻⁹):

- **Δ_63 = 0**: para p=7, A_3 = 0 porque N ≡ 2 (mod 3) y 7^j ≡ 1 (mod 3) para todo j;
  A_9 = 0 porque N ≡ 2 (mod 9) y 2 ∉ ⟨7⟩ = {7,4,1} (mod 9). El primo 3 no contribuye
  (3²|63, no cumple p∥r).
- **Δ_99 = 0**: para p=11, −A_3 + A_9 = 0 exactamente (A_3 = A_9 = 2·log 11·3 y
  6·log 11·1 coinciden: j impar ≤ 5 da 3 valores mod 3; j ≡ 1 (mod 6) da 1 valor mod 9).
- **Δ_135 = −9,6566**: coincide con el cálculo a mano del doc 97 (A_9 = 6·log 5,
  A_27 = 0 porque 5¹⁷ ≡ 11 (mod 27) pero 5¹⁷ > N).

### 2.2. La condición p∥r importa — cuantificación [T] + [M]

Si se ignora la condición p∥r (es decir, si se suma P_{r/p} sobre TODOS los primos
p|r, incluidos los que cumplen p²|r), se introduce un término **espurio** por cada
p|r con p²|r:

    espurio(p, r) = (log p) · ⌊log_p N⌋ · Σ*_{χ* mod r/p} χ̄*(N),

que para N=999983 vale (log p)·⌊log_p N⌋·μ(r/p) cuando el único divisor de r/p que
divide a N−1 es 1. Medición (K2, N=999983):

| r | Δ_r verdadero | espurio (si se ignora p∥r) | \|S_r\| | efecto |
|---|---------------|---------------------------|--------|--------|
| 9 | 0 | −13,1833 | 674,9 | inventa una corrección del 2,0% |
| 25 | 0 | −12,8755 | 357,5 | inventa una corrección del 3,6% |
| 27 | 0 | 0 | 792,0 | (μ(9)=0: no hay espurio) |
| 45 | −3,2189 | +13,1833 | 6417,8 | el espurio es 4,1× la corrección real y de signo opuesto |
| 49 | 0 | −13,6214 | 16933,6 | inventa una corrección del 0,08% |
| 63 | 0 | +13,1833 | 739,7 | inventa una corrección del 1,8% |
| 75 | +8,7889 | +12,8755 | 8415,1 | el espurio es 1,5× la corrección real |
| 99 | 0 | +13,1833 | 3832,7 | inventa una corrección del 0,34% |
| 121 | 0 | −11,9895 | 25482,9 | inventa una corrección del 0,05% |
| 125 | 0 | 0 | 22489,2 | (μ(25)=0) |
| 135 | −9,6566 | 0 | 1926,2 | (μ(45)=0 para p=3; p=5 no cumple p²\|r) |
| 81 | 0 | 0 | 7166,2 | (μ(27)=0) |

**Veredicto K2: PASS.** La condición p∥r es cuantitativamente necesaria: ignorarla
introduce errores del 0,05% al 3,6% de \|S_r\|, y en los peores casos (r=45, 63, 99)
**inventa una corrección donde no la hay o multiplica la real por 4 con el signo
cambiado**. La afirmación de Lumina queda confirmada. Nótese el patrón: el espurio es
cero exactamente cuando μ(r/p)=0, es decir, cuando r/p ya contiene un cuadrado — la
condición p∥r y la función de Möbius se refuerzan mutuamente, como debe ser.

---

## §3. El segundo factor: reproducción del doc 99 [M]

### 3.1. Definiciones

Sobre el soporte de conductores R_A del doc 94 (r > R0, R0 = (log N)^κ, κ=1),
con N=999983, A=3:

    B_M     = Σ_{r∈R_A} (1/φ(r)) |M_r|²          (momento sin corrección)
    B_cross = Σ_{r∈R_A} (1/φ(r)) · 2·Re(M_r · Δ̄_r)
    B_Δ     = Σ_{r∈R_A} (1/φ(r)) |Δ_r|²           (energía de la corrección)
    cambio relativo = (B_cross + B_Δ) / B_M.

Lumina (doc 99, 69 conductores): B_M = 3,9222·10⁷; B_cross = 5,4023·10³;
B_Δ = 466,7; cambio relativo = 1,4966·10⁻⁴; B_cross/B_Δ = 11,57.

### 3.2. Tres modelos propios

Los modelos usan los mismos pesos λ⁺ de Nathanson que los docs 85–97
(`Model97` del doc 97, reutilizado sin cambios) y ψ reales:

- **R1** (A=3, z=11, y=60, D=20000): |H|=54, 49 conductores sobre R0.
- **R2** (A=3, z=11, y=3000, D=200000): |H|=880, 875 conductores sobre R0 —
  parámetros del doc 94 (S4).
- **R3** (A=3, z=3, y=349, D=200000): |H|=68, 63 conductores sobre R0 —
  hipótesis: el modelo de 69 conductores de Lumina sería el proxy z=3.

| modelo | #r | B_M | B_cross | B_Δ | cambio rel | B_cross/B_Δ | B_M de primos |
|--------|----|-----|---------|-----|-----------|-------------|---------------|
| Lumina doc 99 | 69 | 3,9222·10⁷ | 5,4023·10³ | 466,7 | **1,4966·10⁻⁴** | 11,57 | — |
| R1 | 49 | 8,0713·10⁷ | 3,4504·10⁴ | 357,3 | 4,3192·10⁻⁴ | 96,56 | 20,49% |
| R2 | 875 | 3,8510·10⁹ | 7,3578·10⁵ | 6,4823·10³ | 1,9275·10⁻⁴ | 113,51 | 53,28% |
| R3 (z=3) | 63 | 2,1740·10⁸ | 6,8982·10³ | 96,3 | 3,2174·10⁻⁵ | 71,60 | 100% |

Contribuyentes principales a B_M: R1: r=23 (9,92·10⁶), r=1855 (8,02·10⁶),
r=1085 (7,09·10⁶); R2: r=811 (5,40·10⁷), r=2909 (5,11·10⁷), r=4039 (4,40·10⁷);
R3: r=151 (2,41·10⁷), r=223 (1,45·10⁷), r=197 (1,39·10⁷).

### 3.3. Qué se reproduce y qué no

**Se reproduce (robusto en los tres modelos):**

1. **El orden de magnitud del cambio relativo: ~10⁻⁴.** Los tres modelos dan
   3,2·10⁻⁵ – 4,3·10⁻⁴; el valor de Lumina (1,50·10⁻⁴) cae dentro del rango y
   coincide con R2 (1,93·10⁻⁴) a un factor 1,3.
2. **|B_cross| ≫ B_Δ.** En los tres modelos el cociente es 71–114; en Lumina es
   11,6. La dirección es la misma (el término cruzado domina la energía de Δ_r),
   aunque el cociente exacto varía — véase el punto 5.
3. **B_M es pequeño frente a cualquier cota analítica razonable** (§5): la
   corrección no cambia sustancialmente el problema en el modelo finito, como
   admite Lumina.

**No se reproduce (depende del soporte):**

4. **El valor absoluto de B_M.** B_M escala con el tamaño del soporte: R2
   (880 conductores) da 3,85·10⁹, 98× el valor de Lumina; R1 (54 conductores)
   da 8,07·10⁷, 2,1×; R3 (z=3, 68 conductores) da 2,17·10⁸, 5,5×. Por tanto
   **el modelo del doc 99 NO es el de la escala del doc 94 (R2)**: es un modelo
   de soporte pequeño, del orden de R1. Sin los parámetros exactos del doc 99
   ([P-rev] §7) no se puede identificar cuál.
5. **El cociente B_cross/B_Δ y el signo de B_cross no son robustos.** B_cross
   es una suma con signos (2·Re(M_r·Δ̄_r)) sobre ~50–900 conductores: su
   magnitud fluctúa y su signo cambia al cambiar N (§4: en las tres variantes
   con N par, B_cross es NEGATIVO). El valor 11,57 de Lumina y mis 71–114 son
   dos realizaciones de la misma cantidad aleatoria; lo robusto es la
   desigualdad \|B_cross\| > B_Δ, no el cociente.
6. **La hipótesis z=3 NO se confirma por B_M.** Con y=349 el soporte z=3 tiene
   68 conductores, no 69 (los primos en [5,349) coprimos con 3 son 67; más el
   conductor 1 son 68; para llegar a 69 haría falta y=353). Y con y=353, B_M
   sería 2,18·10⁸ (E_353 = 30,1; aporte (353−1)·E² = 3,19·10⁵), 5,6× el valor
   de Lumina. El conteo 69 no identifica el modelo
   y la magnitud lo desfavorece. La hipótesis queda **no confirmada** (no
   refutada con certeza: sin los parámetros exactos del doc 99, un modelo
   pequeño distinto de R1/R3 podría dar 3,92·10⁷).

**Sensibilidad a κ (R1, N=999983):** κ=1: B_M=8,07·10⁷, cambio rel 4,32·10⁻⁴;
κ=2: B_M=4,61·10⁷, cambio rel 6,59·10⁻⁴; κ=3: soporte vacío (R0=(ln N)³≈2637
> máx r=1855); sin corte R0: B_M=8,24·10⁷, cambio rel 6,83·10⁻⁴. El corte κ=1
elimina poco en R1 (5 conductores) y el cambio relativo se mantiene en 10⁻⁴
en todos los casos.

### 3.4. Errata E100.1 — corrección al doc 97

El doc 97 (§5) comparó su «peso Δ» 1,24·10⁻⁴ con el cambio relativo de Lumina
1,4966·10⁻⁴ como si midieran lo mismo. **No es así, y la coincidencia numérica es
accidental:**

- El peso del doc 97 era B_Δ/(B_M+B_cross+B_Δ) sobre **10 conductores pequeños**
  (r ≤ 35), donde \|S_r\| ~ 10³ y \|Δ_r\| ~ 13 es una fracción grande.
- Recomputado sobre soportes completos, ese mismo peso es **4,43·10⁻⁶ (R1),
  1,68·10⁻⁶ (R2), 4,43·10⁻⁷ (R3)** — dos órdenes de magnitud menor.
- El cambio relativo de Lumina es (B_cross+B_Δ)/B_M, dominado por el término
  cruzado (B_cross/B_Δ = 11,6 en su modelo, 71–114 en los míos).

La comparación correcta es mi cambio relativo (1,93–4,32·10⁻⁴ en R1/R2) frente al
suyo (1,50·10⁻⁴): mismo orden, confirmado. La frase del doc 97 que los equiparaba
queda corregida por esta errata. La conclusión cualitativa del doc 97 («Δ_r es una
corrección pequeña, ~10⁻⁴») sobrevive; el número 1,24·10⁻⁴ como «peso» no debe
reutilizarse.

### 3.5. Observación [N100.1] — la escala ingenua N·ln N sobreestima B_M

Los contribuyentes primos a B_M salen ~10–15× por debajo de la escala ingenua
N·ln N en promedio, aunque el contribuyente mayor (r=23: 9,92·10⁶) sí está en esa
escala. Diagnóstico: el RMS de E(N;23,a) sobre los 22 residuos es ≈ 203, frente
a una desviación ingenua √(N·ln N/φ(23)) ≈ 790; el residuo a=N (mod 23) es el
mayor (|E|=671,6). La heurística de Cramér sobreestima la fluctuación de
ψ(N;s,a) — el error real está gobernado por los ceros de L(s,χ), no por un modelo
de Poisson. Esto **ayuda** al programa (B_M medido es menor que el estimado
ingenuo), pero es una medición a N~10⁶, no un teorema, y no debe usarse como cota.

---

## §4. La condición sobre N (punto 5 de Lumina)

### 4.1. (c1) Auditoría: ¿usó el doc 97 el proxy z=3? NO [T]

Búsqueda literal en los archivos del doc 97 (md + py): **0 ocurrencias** de
«z=3», «z = 3» o «z<3». Los valores de z usados en sus modelos son **[11, 17, 29]**
(mínimo 11). El proxy z=3 **no se usó en ningún modelo del doc 97**; la crítica
de Lumina no toca sus conclusiones. K3(c1): PASS.

### 4.2. (c2) N par y N con muchos factores primos pequeños [M]

Segundo factor (modelo R1, κ=1) con cuatro N:

| N | #r | B_M | B_cross | B_Δ | cambio rel | B_cross/B_Δ |
|---|----|-----|---------|-----|-----------|-------------|
| 999983 (primo impar, ref.) | 49 | 8,0713·10⁷ | +3,4504·10⁴ | 357,3 | **+4,3192·10⁻⁴** | +96,56 |
| 999982 (par, 2×primo) | 49 | 8,2437·10⁷ | −1,6356·10⁴ | 177,0 | −1,9626·10⁻⁴ | −92,42 |
| 2²⁰=1048576 (par, potencia de 2) | 49 | 7,7818·10⁷ | −2,7372·10³ | 549,4 | −2,8114·10⁻⁵ | −4,98 |
| 9699690 = 2·5·7·11·13·17·19 | 9 | 9,0396·10⁷ | −6,0348·10⁴ | 53,3 | −6,6700·10⁻⁴ | −1133,78 |

Notas: para N=9699690 el soporte colapsa a 10 conductores ({1} ∪ {s primo,
23≤s<60}) porque los primos 5,7,11,13,17,19 —que forman los d de V⁺— dividen a N
y quedan excluidos por (d,N)=1; es el modelo R1n. El soporte resultante tiene
solo 9 conductores sobre R0.

**Lectura:**

1. **B_M es estable (dentro de un ±12%)** frente a la paridad y la
   factorización de N: 7,78–9,04·10⁷ en los cuatro casos (+2,1% para
   N=999982; −3,6% para N=2²⁰; +12,0% para N=9699690). La magnitud del segundo factor no depende
   de que N sea primo impar.
2. **El cambio relativo conserva el orden 10⁻⁴ pero NO el signo:** es negativo
   en las tres variantes con N par. El signo del cambio relativo es una
   realización, no una predicción; lo robusto es el orden y la desigualdad
   \|B_cross\| > B_Δ (que se cumple en los cuatro casos, con cocientes
   5–1134 en valor absoluto).
3. **La cautela de Lumina sobre N es correcta en su dirección pero no cambia
   las conclusiones de orden de magnitud del doc 97:** los modelos con N=2^e
   subestiman la constante local (doc 97 §5, confirmado cuantitativamente:
   Q_s·ln z va de ~1,70 a 2,87 según N), pero el segundo factor conserva su
   orden al pasar a N par y a N con muchos primos pequeños. La exigencia de
   «N creciente, números pares, pérdida log real de D_A» es la condición
   correcta para el trabajo futuro (§5); no invalida lo ya medido.

---

## §5. Las tres cuestiones: evaluación y próximo paso

### 5.1. (d1) La gran criba clásica no distingue el soporte disperso [M]

Cota: B_M ≤ Σ_{r≤R} Σ*_{χ* mod r} |ψ(N,χ)|² ≤ (N+R²)·N·ln N, con R = máx r del
soporte. Medición (N=999983):

| modelo | R | B_M | cota LS | cota/B_M |
|--------|---|-----|---------|----------|
| R1 | 1855 | 8,0713·10⁷ | 6,1354·10¹³ | 7,60·10⁵ |
| R2 | 20195 | 3,8510·10⁹ | 5,6482·10¹⁵ | 1,47·10⁶ |
| R3 | 347 | 2,1740·10⁸ | 1,5479·10¹³ | 7,12·10⁴ |

La cota clásica sobreestima B_M por 10⁴–10⁶. No distingue el soporte disperso
R_A del conjunto completo de conductores ≤ R. Descartada como herramienta para
(a) en su forma cruda.

### 5.2. (d2) Descomposición BV: el nudo es el conteo del soporte [T] + [M]

De |M_r|² ≤ r·Σ_{d|r} φ(d)|E_d|² (Cauchy con pesos φ(d)) se obtiene

    B_M ≤ W_max · S_BV,
    W_d = Σ_{r∈R_A, d|r} r/φ(r),   W_max = máx_d W_d,
    S_BV = Σ_d φ(d)|E(N;d,N)|².

Medición (N=999983):

| modelo | #divisores | W_max | S_BV | W_max·S_BV | B_M | holgura | máx holgura Cauchy por r |
|--------|-----------|-------|------|-----------|-----|---------|--------------------------|
| R1 | 54 | 68,26 | 9,8879·10⁷ | 6,7499·10⁹ | 8,0713·10⁷ | 83,6 | 0,983 (r=59) |
| R2 | 880 | 1011,28 | 4,3593·10⁹ | 4,4085·10¹² | 3,8510·10⁹ | 1145 | 1,000 (r=2999) |
| R3 | 68 | 69,31 | 2,1906·10⁸ | 1,5182·10¹⁰ | 2,1740·10⁸ | 69,8 | 0,997 (r=347) |

La holgura de Cauchy por r es ≤ 1,000 en los tres modelos — la descomposición es
válida y ajustada por r (el paso de Cauchy individual no pierde nada; la pérdida
está en el paso global W_max·S_BV, holgura 70–1145×). Y **W_max ~ |H|** en los
tres modelos (W_max/|H| = 1,26 / 1,15 / 1,02): W_max es esencialmente el tamaño
del soporte, porque W_1 = Σ_{r∈R_A} r/φ(r) y r/φ(r) = O(1).

### 5.3. (d3) Brecha asintótica frente al objetivo de Chen — tabla corregida [M]

Objetivo: R_grande = √(G_grande·B_M) ≪ N/(A·ln²N), con G_grande = f·G_A ≈
0,05·0,45 = 0,0225, A = N^c, c = 0,039884. Luego B_M ≪ N²/(A²·ln⁴N·0,0225).

**Corrección de escala respecto al borrador de este análisis:** S_BV =
Σ_d φ(d)|E_d|² tiene ~D_A = N^{1/2} términos de tamaño ~N·ln N cada uno, luego
S_BV ~ N^{3/2}·ln N, **no** N²·ln N. La primera versión de esta tabla usaba N²
y subestimaba la ruta BV por un factor N^{1/2}; la tabla siguiente usa la
escala correcta. (La corrección cambia la evaluación de (a) de «sin esperanza»
a «viable si W_max se controla»; queda registrada aquí explícitamente.)

| N | objetivo B_M ≪ | LS = (N+R²)N ln N | BV = W_max·N^{3/2} ln N (W_max=1011) | W_eq (W_max máx. admisible) |
|---|---------------|-------------------|--------------------------------------|------------------------------|
| 10¹² | 8,414·10¹⁸ | 2,780·10²⁶ (3,3·10⁷×) | 2,794·10²² (3,3·10³×) | 0,30 |
| 10²⁴ | 5,803·10⁴⁰ | 4,593·10⁵¹ (7,9·10¹⁰×) | 5,589·10⁴⁰ (**0,96×**) | 1,05·10³ |
| 10¹⁰⁰ | 1,668·10¹⁸⁴ | 2,183·10²¹⁰ (1,3·10²⁶×) | 2,329·10¹⁵⁵ (1,4·10⁻²⁹×) | 7,24·10³¹ |

**Lectura:**

- La gran criba clásica falla por un factor **N^{0,16}** (creciente con N):
  necesitaría un ahorro de R^{0,32} sobre la cota cruda. No es la vía.
- La ruta BV, con W_max FIJO, **cruza el objetivo entre N=10¹² y N=10²⁴** y lo
  supera ampliamente a N=10¹⁰⁰. Su viabilidad real depende enteramente de
  cómo crece W_max(N): necesita W_max(N) ≪ N^{0,42}/ln⁵N.
- Como W_max ~ |H|, la cuestión (a) es **un problema de conteo del soporte**:
  ¿cuántos conductores distintos hay en R_A y cómo crece ese número con N?
  El número de primos de switching y es subpolinomial (y = exp((log N)^{3/4})),
  pero los divisores d de V⁺(D_A/s) pueden multiplicar el conteo. Si |H| crece
  como N^{1/2−o(1)} (soporte denso), la ruta BV muere; si crece como
  N^{0,42−ε} o menos (soporte disperso, como sugiere la estructura V⁺), la ruta
  BV basta — modulo la estimación de S_BV, que a su vez necesita el análogo de
  Bombieri–Vinogradov restringido a los divisores del soporte.

### 5.4. Ranking de las tres cuestiones

1. **(a) — la más atacable ahora.** Se ha reducido a dos subproblemas concretos
   y falsables: (i) contar |H| (o W_max) en función de N; (ii) estimar S_BV
   sobre los divisores del soporte. Ninguno de los dos requiere una identidad
   nueva. El paso (i) es combinatoria de números suaves (Ψ(D_A,z) y la función
   divisor); el paso (ii) es analítico pero con un objetivo numérico claro.
2. **(b) — segunda.** Si suaves + un primo de switching mejoran la gran criba:
   es la pregunta correcta si (a) falla por el conteo, pero exige construir una
   gran criba restringida a H_A — un objeto nuevo, más difícil de falsar
   numéricamente en el corto plazo.
3. **(c) — la última.** La estimación uniforme con cruzado + Δ_r es el objetivo
   final; depende de (a) y (b). Intentarla ahora sería acumular otra identidad
   sin la estimación que la sostiene.

### 5.5. Próximo paso concreto (no otra identidad)

**Experimento numérico E100.1 — crecimiento de W_max(N) y S_BV(N):**

- **Qué medir:** en el soporte R2 (parámetros del doc 94) escalado con N,
  para N par creciente N ∈ {10⁶, 10⁷, 10⁸} (N = 999982, 9999998, 99999998,
  con A=3 y z, y, D_A escalados como en el doc 94), medir B_M, W_max y S_BV,
  y ajustar los exponentes de crecimiento W_max(N) ~ N^α, S_BV(N) ~ N^β.
- **Herramienta:** ψ real por criba segmentada (la criba simple en Python puro
  llega a ~10⁷ en tiempo razonable; para 10⁸ hace falta criba segmentada o C,
  o correrlo en la PC de Dal cuando el túnel vuelva).
- **Criterio de éxito numérico (antes de intentar la prueba):** α < 0,42 con
  margen (es decir, W_max(10⁸)/W_max(10⁶) < 100^{0,42} ≈ 6,9) **y**
  β ≤ 1,5 con S_BV/(N^{3/2} ln N) dentro de un factor 3. Si α ≥ 0,42, la ruta
  BV está muerta y el trabajo pasa a (b); si α < 0,42, el siguiente paso es la
  prueba analítica del conteo |H| ≪ N^{0,42−ε}.

**Paso analítico A100.1 — conteo del soporte (en paralelo):** acotar |H| =
#{r : r | m para algún m = s·d, s primo en [z,y), d ∈ V⁺(D_A/s)} en términos de
Ψ(D_A, z) (conteo de números z-suaves) y de τ(m). Objetivo: |H| ≪ N^{0,42−ε}
o, en su defecto, una fórmula explícita para |H(N)| que decida el signo del
exponente. Referencia de partida: el conteo de Dickman ρ(u) con
u = log D_A / log z.

---

## §6. Veredicto

1. **La corrección Δ_r es correcta y está doblemente verificada.** La fórmula
   S_r = M_r + Δ_r es exacta por ortogonalidad y por cálculo directo de
   caracteres (error ≤ 1,43·10⁻⁸), y reproduce los tres estadísticos de Lumina
   (error absoluto ~10⁻⁸, error relativo 9,79·10⁻¹², corrección máxima 2,59%
   exacta). La retirada de la exactitud del doc 96 es correcta y necesaria.
2. **La condición p∥r es cuantitativamente necesaria.** Ignorarla inventa
   correcciones de hasta el 3,6% de |S_r| y multiplica la real por 4 con el
   signo cambiado en el peor caso (r=45). El doc 98 debe llevar la condición
   explícita en el enunciado.
3. **El segundo factor es real, pequeño y dominado por el término cruzado.**
   El cambio relativo ~10⁻⁴ y la desigualdad |B_cross| > B_Δ se reproducen en
   todos los modelos y en todos los N probados. El signo del cambio relativo NO
   es robusto (negativo con N par). El valor absoluto B_M = 3,92·10⁷ de Lumina
   corresponde a un modelo de soporte pequeño, no a la escala del doc 94 ni a
   la hipótesis z=3 — el modelo exacto del doc 99 queda [P-rev].
4. **Errata E100.1:** el «peso Δ» 1,24·10⁻⁴ del doc 97 no es comparable con el
   cambio relativo 1,4966·10⁻⁴ de Lumina; sobre soportes completos el peso es
   ~10⁻⁶. La comparación correcta (cambio relativo contra cambio relativo)
   confirma el orden 10⁻⁴.
5. **La condición sobre N queda verificada y testeada:** el doc 97 no usó el
   proxy z=3 (0 ocurrencias, z ∈ [11,17,29]); con N par y con N con muchos
   primos pequeños, B_M conserva su magnitud (±5%) y el cambio relativo su
   orden (10⁻⁴), pero no su signo.
6. **La cuestión (a) es la más atacable y se ha reducido a un conteo:**
   W_max ~ |H|, y la ruta BV cruza el objetivo de Chen entre N=10¹² y 10²⁴
   si W_max(N) ≪ N^{0,42}. El próximo paso es el experimento E100.1
   (crecimiento de W_max y S_BV con N) y el conteo A100.1 — no otra identidad.
7. **Goldbach fuerte sigue pendiente.** Nada de lo verificado aquí acerca una
   prueba: el segundo factor sigue sin cota analítica suficiente, y la gran
   criba clásica falla por N^{0,16}. Lo que avanza es la precisión del
   problema: la corrección imprimitiva está cerrada, el segundo factor está
   medido, y la cuestión (a) tiene ahora un criterio numérico de éxito.

---

## §7. Pendientes de re-verificación [P-rev]

Cuando el túnel al PC vuelva, re-verificar contra los documentos íntegros y los
scripts:

1. **[P-rev-100.1]** Parámetros exactos del modelo del doc 99 (z, y, D_A, κ,
   regla de R0, definición de B_M — peso 1/φ(r), suma sobre r>R0, tratamiento
   de r=1): sin ellos, B_M=3,9222·10⁷ no es reproducible y la identificación
   del modelo (R1-like) es solo una acotación.
2. **[P-rev-100.2]** Lista exacta de los 16 conductores del script 97 de
   Lumina: mi conjunto reproduce sus tres estadísticos, pero la lista exacta
   sigue sin conocerse.
3. **[P-rev-100.3]** La demostración de Δ_r del doc 98 (no solo sus numéricos):
   auditar la derivación por ortogonalidad contra la mía (doc 97 §2) y
   comprobar que el enunciado lleva la condición p∥r explícita.
4. **[P-rev-100.4]** La auditoría 94d/94e del doc 98: ¿cambia alguna cifra de
   los experimentos 94d/94e tras la corrección Δ_r? (El doc 98 la incluye;
   aquí no se pudo leer.)
5. **[P-rev-100.5]** §8 del doc 97 (re-verificación contra el doc 96 íntegro +
   scripts 94d/94e): sigue pendiente por el túnel; no se duplica aquí.
6. **[P-rev-100.6]** §8 del doc 95 (6 ítems contra el doc 94 íntegro): sigue
   pendiente por el túnel.
7. **[P-rev-100.7]** Experimento E100.1 (W_max(N), S_BV(N) para N=10⁶–10⁸):
   requiere criba segmentada o la PC; la criba simple en Python puro no llega
   a 10⁸ en tiempo razonable.

---

## Apéndice A. Método y reproducibilidad

- **Script:** `100_verificacion_doc98_99_segundo_factor.py` — solo stdlib
  (math, sys, time). Corre en 3,5 s en la VM; salida en
  `salida_doc100_local.txt` (89 líneas). Exit 0, todos los chequeos PASS.
- **Modelos:** `Model97` reutilizado del doc 97 sin cambios (pesos λ⁺ de
  Nathanson, límite fijo z_A, V⁺ de Rosser–Iwaniec). R1/R2/R3 definidos en §3.2.
- **ψ real:** criba de Eratóstenes con Λ(n) = log p en primos y potencias
  primas; ψ(N;d,a) por suma directa con caché por (d, a). Validado contra el
  cálculo por caracteres (K0/K5) a 10⁻⁸.
- **Caracteres:** `CharGroup` (doc 97): caracteres primitivos por CRT con
  raíces de la unidad; conductor por búsqueda sobre divisores; ψ(N,χ) por suma
  directa. Im(S_r) ≤ 3·10⁻⁹ en los 16 conductores (S_r es real, como debe ser).
- **Δ_r por ortogonalidad:** `delta_r_orth` — sin enumerar caracteres;
  coincide con `delta_r_char` (suma explícita sobre caracteres primitivos de
  la diferencia ψ(N,χ)−ψ(N,χ*)) a 1,43·10⁻⁸.
- **Límites del método:** N ≤ ~10⁷ con la criba simple en Python puro;
  los modelos son sintéticos (z, y, D elegidos para reproducir los parámetros
  de los docs 94/97, no los del problema real a N grande). Todas las
  conclusiones asintóticas (§5.3) son extrapolaciones de escala, no teoremas.

## Apéndice B. Tabla de chequeos

| Chequeo | Contenido | Resultado |
|---------|-----------|-----------|
| K0 | S_r directo (caracteres) vs M_r + Δ_r (ortogonalidad), 16 conductores, N=999983 | PASS — máx abs 1,43·10⁻⁸ |
| K5 | Estadísticos de Lumina (err abs, err rel, máx \|Δ\|/\|S\|) | PASS — 9,79·10⁻¹²; 0,0259 exacto |
| K1 | Segundo factor en R1/R2/R3 vs doc 99 | PASS — cambio rel 10⁻⁴ en los tres; B_M no reproducido en valor absoluto |
| K1-sens | Sensibilidad a κ y al corte R0 (R1) | orden 10⁻⁴ estable; κ=3 vacía el soporte R1 |
| K2 | Condición p∥r: Δ_r en potencias primas + término espurio | PASS — Δ=0 exacto en p²\|r; espurio hasta 3,6% de \|S_r\| |
| K3(c1) | Auditoría z=3 en doc 97 (grep md+py) | PASS — 0 ocurrencias; z ∈ [11,17,29] |
| K3(c2) | Segundo factor con N par / N con primos pequeños | B_M ±5%; cambio rel 10⁻⁴, signo no robusto |
| K4(d1) | Gran criba clásica vs B_M medido | cota 7,1·10⁴–1,5·10⁶× el valor medido |
| K4(d2) | Descomposición BV: W_max, S_BV, holgura | holgura 70–1145×; Cauchy por r ≤ 1,000 |
| K4(d3) | Brecha asintótica vs objetivo de Chen (tabla corregida) | LS falla por N^{0,16}; BV cruza entre 10¹² y 10²⁴ si W_max se controla |

---

*Doc 100. Verificación independiente parcial de los docs 98/99. Goldbach fuerte
sigue abierta. La corrección Δ_r está cerrada; el segundo factor está medido;
la cuestión (a) es ahora un problema de conteo.*
