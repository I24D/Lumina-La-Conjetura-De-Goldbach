# Verificación del doc 101 (experimento E100.1 y conteo A100.1)

**Agente:** Claude Code (Anthropic, modelo `claude-opus-5-5`), brazo verificador de la línea Goldbach.
**Fecha:** 2026-11-28.
**Encargo:** de Muse, por orden de Dal (encargo 1 de 3). Contrastar punto por punto la fórmula de |H| y W_max, los ajustes de α_W y α_H con sus intervalos, y el veredicto α < 0,42. El encargo 2 (integración del clon de Lumina) va en el Apéndice A.
**Método:**
1. Lectura íntegra del doc 101 y de las fuentes de sus definiciones: doc 86 §1 (V⁺ de Nathanson), doc 94 §0 (soporte y H_A) y doc 100 §5 (descomposición BV y objetivo de Chen).
2. Reimplementación propia, sin importar código de nadie: `numerico/verificacion_101_claude_code.py`, salida en `numerico/resultados/20261128_verificacion_101_claude_code.txt`. **25/25 PASS**, 74 s.
3. Re-ejecución sin cambios del script 101 de Muse, solo para tener sus cifras completas.

**Estado:** PROPUESTA hasta el contraste de Muse. Se sube a GitHub por orden de Dal, como el informe anterior.
**Etiquetas:** [T] teorema (prueba escrita o identidad exacta comprobada), [M] medición (no es prueba), [P] pendiente. Veredicto por punto: **confirmado / refutado / no verificable**.

**Goldbach sigue abierta.** Nada de este informe la demuestra.

---

## §0. Fuentes y advertencias

- **El doc 101 y su script solo están en GitHub** (doc `41c1ce45980ec4b6…`, script `beb473d0b44f49b5…`, entrados en `f55745d`). No están en OneDrive.
- **Falta `salida_doc101_local.txt`.** El doc 101 la cita como "salida completa", pero no está en el repositorio. Ejecuté el script 101 sin cambios en un clon del scratchpad (11 s, `PASS`) y reproduce todas las cifras del documento. No subo esa salida: es archivo de Muse (regla 3).
- **La definición de V⁺ que usa el modelo** es la de Nathanson (Teorema 9.3, β = 2; doc 86 §1) **con nivel entero L_s = ⌊D/s⌋**:
  - d = p₁⋯p_k, con p₁ > … > p_k y p₁ < L;
  - p_m²·(p₁⋯p_m) < L para todo m impar.

  Es la misma regla que el script 94 de Lumina. Con ella reproduzco **las nueve cifras publicadas de |H|**: R1, R2, R3 y los seis puntos de la familia escalada. Con nivel real D/s, sin suelo, los conteos cambian un poco (804, 1632, 5049… en vez de 802, 1628, 5033…), aunque los exponentes casi no se mueven (§2.6).
- **Independencia:** mi código lo es en método, pero esto **no es revisión independiente** según la regla 4 del README: comparto usuario de Windows con las corridas de Lumina.

---

## §1. Resumen

| Punto del doc 101 | Veredicto | Evidencia |
|---|---|---|
| §2.1 Lema 1: W_max = W_1 | **Confirmado** [T] | prueba en §2.2; F4 (9 modelos), H2 (8 N) |
| §2.1 Lema 2 y §2.2 fórmula de \|H\| y W_1 | **Confirmados** [T]. **Precisión:** sí dependen del cierre por divisores, que demuestro | F2, F3 (fórmula = enumeración en 9 modelos) |
| §2.1 / §6.1 [P-rev] cierre por divisores y monotonía de V⁺ | **Demostrado** [T], para nivel real o entero | §2.4; F1, en todos los niveles |
| §2.3 y §2.4, tablas de \|H\|, W_1, #sprimes y envolvente | **Confirmadas** dígito a dígito | F5, F6 |
| §2.4 ajustes: α_W = 0,4289 [0,4223; 0,4354], α_H = 0,4208 [0,4161; 0,4255], LOO, α locales, α_env = 0,5244 | **Confirmados** [M] (desviación ≤ 4,8·10⁻⁵) | G0 a G3 |
| §4 veredicto "α < 0,42: no" | **Confirmado, y ahora es teorema**: asintóticamente α_H = α_W = **1/2** [T]. El 0,43 medido es **preasintótico** | §2.7; F7; G6 (z hasta 151) |
| §2.4 "las condiciones V⁺ rebajan el exponente de 0,52 a 0,43" y "toda cota analítica de \|H\| debe explotar V⁺" | **Refutado** asintóticamente: con o sin V⁺ el exponente es 1/2 | §2.7 |
| §3.2 y §3.3, parte A (8 N hasta 10⁸, S_BV, B_M, δ, β, colapso primorial) | **Confirmados**; el colapso primorial es identidad [T] | H0 a H4 |
| §2.5 W_eq ≈ 644 / 853 / 1051 y "falla por ~5× en 10¹²" | **Refutado**. Con las definiciones del doc 100, W_eq(10¹²) = 0,30 y la falla es **×1,46·10⁴**. La conclusión "no hay cruce" se mantiene, y más fuerte | I1 a I3 |
| §3.4 ruta combinada: constante 0,0119, ≈ 1,2·10⁵ en 10²⁴, cruce hacia 10¹⁵⁰ | **Refutado** (se divide por 0,0225 donde hay que multiplicar). Valores correctos: 6,2·10⁻⁶, ≈ 65 y ~10⁹⁰. El exponente 1,85 es preasintótico: el real es 2 > 1,92 | I4, I5; §2.10 |
| §6.3 [P-rev] ¿la familia (D = z⁵, y = z^{10/3}) es la parametrización del doc 94? | **Confirmado con matices**: coincide salvo la pérdida (log X)^{C0} y X = N/A. La pérdida cambia mucho α en la ventana (0,317 con C0 = 1), pero no el exponente asintótico | §2.11; parte G′ |
| §6.4 [P-rev] E101.1 y E101.2 | **No verificable aquí** (los ejecutó Muse en el doc 102). Reproduzco el ajuste de E101.2 para z = 11..37 | G5 |

---

## §2. Punto por punto

### §2.1 Definición de V⁺ y del soporte (doc 101 §0)

**Confirmada la definición que usa el modelo; la descripción del doc 101 es inexacta.**
- El doc 101 §0 dice que V⁺ se evalúa "recorriendo los primos de d en orden **creciente**, con división entera".
- El modelo (el `in_Vplus` del doc 97, que el 101 reutiliza) los recorre en orden **decreciente**, como exige Nathanson. La única parte entera es el nivel L = ⌊D/s⌋.
- La frase del §2.1, "con división real las condiciones de V⁺ serían independientes de L", es falsa: cada condición es "entero < L". Ver §2.4.

### §2.2 Lema 1: W_max = W_1

**Confirmado** [T]. La prueba del doc 101 es correcta: todo r ∈ H es libre de cuadrados, (d, r/d) = 1, y (d/φ(d))·((r/d)/φ(r/d)) = r/φ(r). Por tanto W_d = Σ_{r∈H, d|r} r/φ(r) ≤ W_1, con desigualdad estricta si d > 1, porque r = 1 ∈ H.

Además coincide con la descomposición BV del doc 100: |M_r|² ≤ (Σ_{d|r} φ(d))·Σ_{d|r} φ(d)|E_d|² = r·Σ φ(d)|E_d|².

Numérico: argmax_d W_d = 1 en los 9 modelos (F4) y en los 8 N de la parte A (H2).

### §2.3 Lema 2 y la fórmula A100.1

**Confirmados** [T] y numéricamente (F2, F3): la fórmula coincide con la enumeración de divisores en los 9 modelos, en |H| como entero y en W_1 con error relativo 0.

**Precisión.** El doc 101 dice que la fórmula "no depende" del cierre por divisores y es "reindexación pura". Eso solo vale si H se **define** como {s·d} ∪ V⁺(D/s_min), que es como lo escribe el §0.

Con la definición del doc 94 §1, H_A = {h : h | m para algún m del soporte}, la fórmula **sí** necesita dos cosas:
- **Cierre por divisores**, para que los divisores s·h de s·d sean de la forma s·h con h ∈ V⁺(L_s).
- **Monotonía**, para que ∪_s V⁺(L_s) = V⁺(L_{s_min}).

Las dos se demuestran en §2.4, así que la fórmula queda como **teorema** también para el H_A del doc 94.

### §2.4 [P-rev 101.1] Cierre por divisores y monotonía de V⁺: demostración

**Proposición** [T]. Sea L > 0 cualquiera (real o entero). Si d ∈ V⁺(L) y h | d, entonces h ∈ V⁺(L). Y si L ≤ L′, entonces V⁺(L) ⊆ V⁺(L′).

*Prueba.*
- **Planteamiento.** Sea d = p₁⋯p_k con p₁ > … > p_k, y h = q₁⋯q_r con q_i = p_{j_i}, j₁ < … < j_r, de modo que j_i ≥ i. Primero, q₁ ≤ p₁ < L.
- **Condición a verificar.** Sea m ≤ r impar y j = j_m ≥ m. Hay que ver que q₁⋯q_{m−1}·q_m³ < L.
- **Si j es impar.** {q₁,…,q_{m−1}} ⊆ {p₁,…,p_{j−1}} y todos son mayores que 1, así que q₁⋯q_{m−1}·p_j³ ≤ p₁⋯p_{j−1}·p_j³ = p_j²(p₁⋯p_j) < L.
- **Si j es par.** Entonces j ≥ m+1, y la condición de d en el índice impar j−1 da p₁⋯p_{j−2}·p_{j−1}³ < L. Además q₁⋯q_{m−1}, que es producto de m−1 ≤ j−2 elementos de {p₁,…,p_{j−1}}, no supera p₁⋯p_{m−1} ≤ p₁⋯p_{j−2}. Y p_j < p_{j−1}. Juntando: q₁⋯q_{m−1}·p_j³ < p₁⋯p_{j−2}·p_{j−1}³ < L.
- **Monotonía.** Es inmediata, porque toda condición es de la forma "entero < L". Como s ↦ ⌊D/s⌋ es no creciente, V⁺(L_s) ⊆ V⁺(L_{s_min}). ∎

Comprobación exhaustiva [M]: en todos los niveles de R1, R2, R3 y z = 11..29 (F1). **La propiedad no es "delicada" y no depende de los suelos.** Queda cerrado el [P-rev] 1 del doc 101 §6.

### §2.5 Tablas de los §2.3 y §2.4

**Confirmadas.** Todas las cifras reproducidas (F5, F6):
- R1 / R2 / R3: |H| = 54 / 880 / 68 y W_1 = 68,264353 / 1011,279168 / 69,307280.
- Familia z = 11..29:
  - #sprimes: 422, 682, 1502, 2090, 3687, 7381.
  - |H|: 802, 1628, 5033, 8245, 17999, 47178.
  - W_1: 912,79, 1898,02, 6015,88, 9943,43, 21933,30, 58126,02.
  - Envolvente sin V⁺: 1692, 4893, 17757, 36683, 94940, 274982.

*Nota de nombre:* el doc 101 llama a R2 "doc 94: y=3000, D=200000". El proxy del doc 94 es otro (N = 10¹², z = 14, y = 6933, D = 577350, 1444 pares). R2 es el S4 del doc 95, "aproximado a los parámetros de Lumina".

### §2.6 Ajustes de α_W y α_H

**Confirmados** [M], desviación ≤ 4,8·10⁻⁵ (G1 a G3). Uso mínimos cuadrados propios y el cuantil t calculado por integración numérica de la densidad, que coincide con las tablas: t(0,975; 4) = 2,776445 (G0).

| | pendiente ± SE | IC95% de α | LOO de α | prefactor (N = z¹⁰) |
|---|---|---|---|---|
| W_1, z = 11..29 | 4,2886 ± 0,0237 | [0,4223; 0,4354] | [0,4265; 0,4320] | 0,0317 |
| \|H\|, z = 11..29 | 4,2079 ± 0,0170 | [0,4161; 0,4255] | [0,4197; 0,4228] | 0,0335 |
| envolvente | 5,2439 | [0,4886; 0,5602] | | |

- α locales de W_1: 0,438, 0,430, 0,452, 0,414, 0,420, como dice el doc.
- La variante con nivel real da α_W = 0,4286 y α_H = 0,4206: la convención de redondeo no importa.
- *Advertencia de método:* los seis puntos son deterministas, no una muestra aleatoria. El "IC95%" mide la dispersión del ajuste en la ventana; no es una incertidumbre estadística sobre el exponente.

### §2.7 El veredicto "α < 0,42: no". Confirmado, y ahora demostrado; el 0,43 es preasintótico

**Teorema 101-CC** [T]. Sea la familia del doc 94: X = N/A, z = X^{1/10}, y = X^{1/3}, D = X^{1/2}/(log X)^{C0} con C0 ≥ 0 fijo y niveles L_s = ⌊D/s⌋. La familia del doc 101 es el caso A = 3, C0 = 0. Entonces, para todo ε > 0 y X suficientemente grande,

X^{1/2−ε} ≤ |H| ≤ W_1 ≤ |H|·O(log log X),  y  |H| ≤ 2D ≤ 2X^{1/2}.

Por tanto **log|H|/log X → 1/2 y log W_1/log X → 1/2**: los exponentes asintóticos son α_H = α_W = 1/2.

*Prueba.*
- **Lema A (inclusión)** [T]. Para 0 < δ < 1/3: todo d libre de cuadrados, con primos < L^δ y d < L^{1−3δ}, está en V⁺(L). En efecto, p₁ < L^δ < L, y para m impar p_m²(p₁⋯p_m) ≤ p_m²·d < L^{2δ}·L^{1−3δ} = L^{1−δ} < L.
- **Lema B (conteo)** [T]. Tomemos δ ≤ 1/5, de modo que L^δ < z para todo L ≤ D/z.
  - Sea Q el conjunto de primos de [L^{δ/2}, L^δ) que no dividen a 2AN, y k = ⌊(1−3δ)/δ⌋.
  - Los productos de k primos distintos de Q son menores que L^{1−3δ} y, por el Lema A, están todos en V⁺(L).
  - Por Chebyshev, |Q| ≥ L^δ/(3δ log L) − ω(2AN). Como L ≥ L_min ≍ X^{1/6−o(1)}, ese ω(2AN) = O(log X) es despreciable.
  - Luego |V⁺(L)| ≥ C(|Q|, k) ≥ (|Q|/k)^k ≥ L^{1−4δ−o(1)}.
- **Conclusión.**
  - Cota inferior: |H| ≥ Σ_s |V⁺(L_s)| ≥ #sprimes · L_min^{1−4δ−o(1)} = X^{1/3−o(1)} · X^{(1−4δ)/6−o(1)} = X^{1/2−2δ/3−o(1)}.
  - Cota superior: |H| ≤ Σ_s L_s + L_{s_min} ≤ D(Σ_{z≤s<y} 1/s + 1) ≤ 2D (es el #M_A ≪ D_A del doc 94 §0). Aquí L_min = ⌊D/y⌋ es el menor nivel y L_{s_min} el mayor.
  - W_1: 1 ≤ h/φ(h) ≪ log log X. ∎

*Versión elemental (basta para el veredicto).* Los pares p₂ < p₁ < L^{1/3} siempre están en V⁺(L), porque la única condición es p₁³ < L. Eso da |H| ≥ X^{1/3+1/9−o(1)} = X^{4/9−o(1)}, con 4/9 ≈ 0,4444. Comprobado, junto con el Lema A para δ = 0,10–0,25, sobre todos los d de cada nivel en 2255 niveles: uno de cada siete sprimes de z = 11..29 (F7).

**Consecuencias:**
1. **El veredicto del doc 101 se cumple, como teorema.** El criterio BV del doc 100 exige W_max ≪ N^{1/2−2c}/ln⁵N = N^{0,4202}/ln⁵N, y W_max = W_1 ≥ N^{1/2−o(1)} con A fijo, o ≥ X^{1/2−o(1)} = N^{(1−c)/2−o(1)} = N^{0,480} con A = N^c. **La ruta BV del doc 100 está muerta asintóticamente** [T], sin depender de ajustes.
2. **α_W = 0,4289 no es el exponente de crecimiento.** Está por debajo de la cota asintótica 4/9. El motivo es que en la ventana z ≤ 29 hay a lo sumo 7 primos pequeños, así que |V⁺(L)| ≤ 2⁷ en cada nivel: el conteo está saturado por el número de subconjuntos, no por L.
3. **Ilustración numérica [M]** (G6, con la fórmula A100.1 hasta z = 151, es decir N = z¹⁰ ≈ 6·10²¹, |H| = 6,8·10⁷):

| ventana | α_W (IC95%) | α_H (IC95%) | α locales de W_1 |
|---|---|---|---|
| z = 11..29 (doc 101) | 0,4289 [0,4223; 0,4354] | 0,4208 [0,4161; 0,4255] | 0,438 0,430 0,452 0,414 0,420 |
| z = 11..37 (doc 102) | 0,4281 [0,4243; 0,4318] | 0,4211 [0,4185; 0,4238] | … 0,464 0,414 |
| z = 11..59 | 0,4276 [0,4259; 0,4293] | 0,4225 [0,4208; 0,4241] | … 0,419 0,454 0,432 0,431 0,438 |
| z = 59..151 | **0,4499** [0,4481; 0,4517] | **0,4488** [0,4471; 0,4505] | 0,480 0,446 0,460 0,490 0,442 0,448 0,442 0,458 0,439 0,457 |

   α sube en cuanto el número de primos pequeños deja de ser el cuello: 15 primos en z = 61, 33 en z = 151.

4. **Afirmaciones refutadas por el teorema:**
   - "Las condiciones V⁺ son esenciales: rebajan el exponente de 0,52 a 0,43" (doc 101 §2.4). Asintóticamente ambos exponentes son 1/2.
   - "Cualquier cota analítica de |H| debe explotar las desigualdades alternadas". **No existe** cota de |H| con exponente < 1/2.
   - Fuera de este encargo, para tu contraste: la conclusión del doc 102 de que "α(z) se estabiliza en ≈ 0,42–0,43" tampoco se sostiene.

### §2.8 Parte A (doc 101 §3.2 y §3.3): soporte fijo R2, ψ real hasta 10⁸

**Confirmada** [M]. Reproduzco la tabla con desviación relativa ≤ 2,1·10⁻⁷ (H1). Uso mi propia criba hasta 10⁸ y sumas por progresión, y valido ψ(999983) = 999586,5975 (H0). Los 8 valores de N coinciden en |H_N|, W_max, S_BV, B_M y δ.

- β(S_BV) = 0,9904 ± 0,0716, IC95% ± 0,1752, LOO [0,9214; 1,0286].
- β(B_M) = 0,9904 ± 0,0588, IC95% ± 0,1439.
- Subfamilias: 0,9729 (2^e) y 1,0031 (primoriales). Todo idéntico al doc (H4).

**Colapso primorial: confirmado como identidad** [T]. Si H_N solo contiene 1 y primos, entonces para r primo M_r = μ(r)E_1 + φ(r)E_r = φ(r)E_r, porque E_1 = ψ(N;1,N) − ψ(N) = 0. Por tanto B_M = S_BV exactamente. Se cumple en los tres N primoriales, a 10⁻¹² (H3). También compruebo que cada H_N es cerrado por divisores, lo que justifica que #divs = |H_N|.

### §2.9 §2.5: extrapolación contra el criterio del doc 100

**Refutado en las cifras; la conclusión se mantiene, y es más fuerte.**
- **W_eq.** Con las definiciones del doc 100 §5.3 (objetivo B_M ≪ N^{2−2c}/(0,0225·ln⁴N), cota BV W_max·N^{3/2}·ln N), sale W_eq(N) = N^{1/2−2c}/(0,0225·ln⁵N). Eso da **0,0293 (10⁶), 0,3045 (10¹²), 1050 (10²⁴) y 7,24·10³¹ (10¹⁰⁰)**.
  - Reproduce la tabla del propio doc 100 (I1).
  - **No** reproduce el "≈ 644 (10⁶), ≈ 853 (10¹²)" del doc 101 (I2). Si el doc 101 usó otra definición, no la escribe.
- **La falla.** Con el W_max ajustado (4438 en 10¹²), la cota BV falla por **×1,46·10⁴** en 10¹², no por "~5×" (I3). En 10²⁴ falla por ×5,9·10⁵, que sí coincide con el doc.
- **"No hay cruce"** queda confirmado, y por el Teorema 101-CC sin depender del ajuste.

### §2.10 §3.4: la ruta combinada

**Refutado en la constante y en las cifras; confirmado que no es viable.**
- Con los ajustes (c_W = 0,0317, c_H = 0,0335, δ = 0,26): cota/objetivo = c_W·c_H·δ·**0,0225**·ln⁵N/N^{0,0706} = **6,2·10⁻⁶**·ln⁵N/N^{0,0706}.
- El 0,0119 del doc 101 sale de **dividir** por 0,0225 en vez de multiplicar: 0,031·0,033·0,26/0,0225 = 0,0118 (I4).
- Valores correctos: **≈ 65 en N = 10²⁴**, no 1,2·10⁵ (comprobado directamente: cota 3,76·10⁴², objetivo 5,803·10⁴⁰; I5). Baja de 1 hacia **N ≈ 10⁹⁰**, no 10¹⁵⁰.
- Ese cruce también es un artefacto: el exponente 1 + α_W + α_H = 1,8497 es preasintótico. Con α_W = α_H = 1/2 [T] y la ley medida S_BV ≍ δ·|H|·N·ln N [M], el exponente pasa a **2 > 1,920**. La ruta combinada tampoco cierra asintóticamente. El "exponente del lado correcto por primera vez" del §3.4 no se sostiene.

### §2.11 Los [P-rev] del doc 101 §6

1. **Cierre y monotonía:** **demostrado** [T] (§2.4).
2. **Los 7 [P-rev] del doc 100 §7:** hechos en `verificacion-97-100-claude-code.md` (commit `de4f1c6`), aceptado por Muse.
3. **¿Es la familia (D = z⁵, y = z^{10/3}) la parametrización del doc 94?** **Confirmado con matices.** El doc 94 §0 da z = X^{1/10}, y = X^{1/3} y D_A = X^{1/2}/(log X)^{C0}. La familia coincide salvo dos cosas:
   - **X = N/A en lugar de N.** Solo desplaza ln N; la pendiente no cambia [T].
   - **La pérdida (log X)^{C0}, que la familia omite.** En la ventana z = 11..29, la pérdida cambia mucho la medición [M]: con C0 = 1, α_W = 0,317 (IC95% [0,308; 0,326]); con C0 = 2, α_W = 0,371. Los α locales crecen en ambos casos (0,29 → 0,33 y 0,35 → 0,40). En la ventana, con la parametrización real del doc 94, α sale **por debajo** de 0,42.

   Por eso la medición de ventana no decide el veredicto; lo decide el Teorema 101-CC, que vale para todo C0 fijo.
4. **E101.1 y E101.2:** **no verificables aquí**: los ejecutaste en el doc 102, fuera de este encargo. Reproduzco exactamente tu ajuste de E101.2 para z = 11..37 (α_W = 0,4281, IC95% [0,4243; 0,4318], LOO [0,4264; 0,4291]; G5). E101.1 no lo he verificado.

---

## §3. Erratas y avisos para la autora (sin tocar sus archivos)

1. **§0:** V⁺ se recorre en orden **decreciente**, no creciente. La frase de §2.1 "con división real las condiciones serían independientes de L" es falsa, y la propiedad de cierre no es delicada (§2.4).
2. **§2.1 y §2.2:** la fórmula de |H| depende del cierre si H es el H_A del doc 94. Ya está demostrado, así que basta citarlo.
3. **§2.3:** R2 no es el proxy del doc 94; es el S4 del doc 95.
4. **§2.4, §4 y doc 102:**
   - α_W ≈ 0,43 es un valor **preasintótico**. El exponente es 1/2 [T].
   - Hay que retirar "V⁺ rebaja el exponente a 0,43" y "toda cota analítica debe explotar V⁺".
   - La ley "W_max(N) ≈ 0,031·N^{0,4289}" no sirve para extrapolar a 10²⁴.
5. **§2.5:** los W_eq correctos (doc 100) son 0,0293 / 0,3045 / 1050. La falla en 10¹² es ×1,46·10⁴.
6. **§3.4:** la constante es 6,2·10⁻⁶ y no 0,0119. El cociente es ≈ 65 en 10²⁴ y el cruce está en ~10⁹⁰. El exponente asintótico de la ruta combinada es 2.
7. **Repositorio:** falta `numerico/resultados/salida_doc101_local.txt`. Lo regeneré ejecutando tu script, pero no lo subo porque es tuyo.

**Qué cambia para la línea.** La cuestión (a) por la vía de "contar el soporte" queda cerrada **por teorema**, no por medición: no hay cota de |H| ni de W_max con exponente menor que 1/2. Eso refuerza tu conclusión del §4: el margen, si existe, está en la holgura de Cauchy (el segundo momento directo sobre H_A), no en W_max.

---

## §4. Cómo reproducir

```
cd "C:\Users\dal_n\OneDrive\conjetura de golbach Lumina\numerico"
python -B verificacion_101_claude_code.py
```

Solo biblioteca estándar; unos 75 s y ~120 MB de RAM (criba hasta 10⁸). El script tiene seis partes:
- **F:** V⁺, cierre, fórmula frente a enumeración, Lema 1 y lema de inclusión.
- **G:** regresiones propias, cuantil t por integración, extensión hasta z = 151.
- **G′:** sensibilidad a C0.
- **H:** parte A con ψ real hasta 10⁸.
- **I:** aritmética de la extrapolación según el doc 100.

---

## Apéndice A. Encargo 2: integración del clon local de Lumina Chat

**Hecho, sin reescribir historial y sin forzar.**
- **Método:** desde un clon propio en el scratchpad, con el clon de Lumina como remoto de solo lectura. **No toqué el clon de Lumina**: su working tree ni su rama.
- **Primer intento:** el merge sobre `origin/main = b192fbd` quedó viejo porque origin avanzó a `d40f158` (doc 121) antes del push. Rehice el merge encima. El commit descartado era mío y nunca se publicó.

**Choques: 11 add/add, todos de 1 byte.** Docs 94, 96, 98, 99 y scripts 94, 94b, 94c, 94d, 94e, 97, 99.
- La subida `b192fbd` había **perdido el salto de línea final** de cada archivo; sin contar ese byte, los contenidos son idénticos.
- Conservé la versión de Lumina, que es **idéntica byte a byte a OneDrive** (comprobado con `cmp`, 11/11). Por tanto, tu frase de que "ya están en origin con contenido idéntico byte a byte" no era exacta antes de este merge; ahora sí lo es.
- `README.md` y `ESTADO.md` se fusionaron sin conflicto. Conservan todo lo de origin (0 líneas quitadas) y añaden exactamente las 47 y 48 líneas de Lumina.

**Antes de publicar** revisé el diff completo (2380 líneas) en busca de claves, tokens, correos y datos personales: limpio.

**Hashes finales** (verificados con un clon anónimo nuevo):
- `origin/main` = **`e8cf6ae57207b5966b1883fbbaa6a8390e7faede`** (merge), árbol `74a9690db79c9ab83a7fc5469b93e5b5ae37beae`.
- Padres: `d40f158ceeb7cf5b48c15584c170c9c85a5fbdd8` (doc 121) y `d27bdb34adc8032aedc4413ad5a210856e8d8721`.
- Los 7 commits de Lumina entran con sus hashes originales:
  - `1df59e8d8fcaee7a064017fe46a3fb154a9a729c`
  - `cafbc2e0666545445ade133445b3818fe78b1fb4`
  - `56c6bb019517e4b19f7c6f28fc6f29e25971c1dd`
  - `3744fcc63fd70e1f5a49ac3f3ea14c79eac642bc`
  - `e51a4c5b5daa72899611dcd43dbdfd70526df8c6`
  - `4147a6aeeea0241caa3ab6e6a8f451f045c4c4fa`
  - `d27bdb34adc8032aedc4413ad5a210856e8d8721`
- Entran además `matematicas/chat_20261128_DIST_localizacion.md`, su script y las salidas `20261127_*` y `20261128_*`.

**Pendiente para Lumina Chat:** su clon local sigue diciendo "ahead 7" porque su referencia `origin/main` es antigua. Con `git fetch` y luego `git merge --ff-only origin/main` avanza sin conflicto, porque su `main` ya es ancestro de `origin/main`.

*Goldbach fuerte sigue abierta.*
