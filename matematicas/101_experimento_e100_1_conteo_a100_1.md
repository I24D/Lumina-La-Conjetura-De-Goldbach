# DOC 101 — Experimento E100.1 y conteo A100.1: crecimiento de W_max y S_BV, fórmula exacta de conteo y veredicto sobre la ruta BV para la cuestión (a)

**Serie:** verificación y continuación de los docs 94–100 (Gram de Lumina, diagonalización, puente de Möbius, segundo factor).
**Fecha:** 2026-11-28. **Autora:** Lumina Muse.
**Archivos:** `101_experimento_e100_1_conteo_a100_1.py` (solo stdlib; corre en local, **PASS** en 18 s), `salida_doc101_local.txt` (salida completa).
**Estado:** PARCIAL en el mismo sentido que los docs 95/97/100: los textos completos de los docs 98/99 y los scripts 94d/94e/97/99 de Lumina siguen sin estar disponibles en local (túnel caído); este doc trabaja sobre las definiciones del modelo finito tal como quedaron fijadas en los docs 94–97 y 100, que sí están completos en local. Los puntos [P-rev] heredados se listan en §6.
**Disciplina:** se distingue en todo momento MEDICIÓN de PRUEBA. Nada de lo aquí medido demuestra Goldbach; la conjetura fuerte sigue abierta. Sin circularidad, sin GRH sin probar.

**Etiquetas:** [T] = teorema (demostrado en el texto o en un doc anterior citado), [N] = medición numérica reproducible con el script adjunto, [P-rev] = pendiente de revisión contra el texto completo.

---

## 0. Alcance

El doc 100 dejó la cuestión (a) — *¿existe una cota útil para B_M restringido a los conductores admisibles R_A?* — reducida a un problema de conteo, y propuso dos tareas:

- **E100.1 (experimento):** medir el crecimiento de W_max(N) y de S_BV(N) con N par creciente en el modelo finito de los docs 94–100 (A = 3, soporte del doc 94, pesos de Nathanson), ajustar el exponente α en W_max(N) ~ N^α, y contrastarlo con el criterio de la ruta BV del doc 100: la ruta cruza el objetivo de Chen entre N = 10¹² y 10²⁴ **si** α < 0,42.
- **A100.1 (conteo):** contar |H_A| y W_max desde la definición del soporte, con fórmula, no solo por enumeración; verificar la fórmula contra la enumeración y extrapolar su tasa.

Este doc ejecuta ambas. El resultado principal es negativo para la ruta BV tal como fue formulada, y es un resultado de **medición** con intervalos estrechos, no una impresión: **α_W = 0,4289 con IC95% [0,4223, 0,4354] — el intervalo completo queda por encima de 0,42**.

Recordatorio de notación (docs 94–100). El soporte es H = {m = s·d} ∪ V⁺(D/s_min), donde s recorre los «sprimes» (primos en [z, y) coprimos con A·N... en el modelo: primos p ∈ [z, y) con p ∤ A), d recorre V⁺(D/s) ⊆ divisores libres de cuadrados de P_A(z) (primos < z, ∤ A), y V⁺(L) es el conjunto de d que superan las desigualdades alternadas de Nathanson al nivel L (lado-q × lado-p comparado con L, recorriendo los primos de d en orden creciente, con división entera). W_d = Σ_{r∈H, d|r} (d/φ(d))·((r/d)/φ(r/d)), W_max = max_d W_d, S_BV = Σ_d φ(d)|E_d|² sobre los divisores d de elementos de H coprimos con N, B_M = Σ_{r > (ln N)^κ} (1/φ(r))|M_r|² con M_r = Σ_{d|r} μ(r/d)φ(d)E_d (κ = 1, como en el doc 100), E_d = error de ψ en progresiones (ψ real).

---

## 1. Resumen ejecutivo

| Pregunta | Respuesta | Tipo |
|---|---|---|
| ¿La fórmula de conteo A100.1 reproduce la enumeración? | **Sí, exactamente**: \|H\| exacto y W_1 = W_max con error relativo ≤ 5,6·10⁻¹⁶ en R1/R2/R3, y contra enumeración explícita independiente en la familia escalada (z = 13, 17, 19, 23) | [T] + [N] |
| ¿W_max = W_1 siempre? | **Sí**: el máximo de la descomposición BV se alcanza en d = 1, por una identidad exacta término a término | [T] |
| ¿α < 0,42 (W_max ~ N^α, familia escalada de los docs)? | **No**: α_W = 0,4289, IC95% [0,4223, 0,4354], LOO [0,4265, 0,4320]. Para \|H\|: α_H = 0,4208, IC95% [0,4161, 0,4255] (a caballo de 0,42) | [N] |
| ¿Cómo crece S_BV a soporte fijo con N par? | S_BV ~ δ·#divs·N·ln N con δ ∈ [0,256, 0,376] decreciente suave; exponente medido β = 0,9904 ± 0,1752 (IC95%), subfamilias puras 0,973 (2^e) y 1,003 (primoriales). La ley N^{3/2} ln N del doc 100 es una **envolvente de crecimiento del soporte**, no el crecimiento a soporte fijo: S_BV/(N^{3/2} ln N) decae de 0,32 a 0,012 entre N = 10⁶ y 10⁸ | [N] |
| ¿La ruta BV es viable para (a)? | **No, en la forma formulada en el doc 100.** Con W_max medido, la cota BV no cierra ni en N = 10¹² (falla por ~5×) y en 10²⁴ falla por ~6·10⁵, con brecha creciente. La ruta combinada (ambos factores medidos) tiene exponente favorable (1,850 < 1,920) pero pierde por constantes y por ln⁵N hasta N ~ 10¹⁵⁰ | [N] |
| ¿Dónde queda el margen? | En la **holgura de Cauchy** (doc 100: B_M real es 70–1145× menor que W_max·S_BV). El siguiente paso es atacar el segundo momento directamente (cuestión (b)), no vía W_max | §4 |

---

## 2. A100.1 — El conteo: fórmula exacta

### 2.1 Dos lemas

**Lema 1 (W_max = W_1) [T].** Todo r ∈ H es libre de cuadrados. Para d | r, los conjuntos de primos de d y de r/d forman una partición de los primos de r, luego

(d/φ(d)) · ((r/d)/φ(r/d)) = ∏_{p|d} (1−1/p)⁻¹ · ∏_{p|r/d} (1−1/p)⁻¹ = r/φ(r).

Por tanto W_d = Σ_{r∈H, d|r} r/φ(r) ≤ Σ_{r∈H} r/φ(r) = W_1, con desigualdad estricta para d > 1 (el término r = 1 ∈ H, que vale 1, no es divisible por d > 1, y todos los términos son positivos). ∎

*Verificación [N]:* en los tres modelos del doc 94 y en los 8 valores de N de la parte A (§3), el argmax de la descomposición BV es d = 1 sin excepción (chequeos J2, J6), y W_1 calculado por la fórmula de §2.2 coincide con el W_max de la descomposición BV con error relativo ≤ 5,6·10⁻¹⁶.

**Lema 2 (reindexación) [T].** Cada m = s·d ∈ H con s sprime tiene factorización única: s es el único primo ≥ z de m (los primos de d son < z). Además s·d ∉ V⁺(D/s_min) (sus elementos tienen todos los primos < z). Luego el conteo sobre H se reindexa sin solapamiento sobre los pares (s, d). ∎

**Propiedad estructural (verificada, demostración [P-rev]).** V⁺(L) es cerrado por divisores y monótono en el nivel: si d ∈ V⁺(L) y h | d, entonces h ∈ V⁺(L); y V⁺(D/s) ⊆ V⁺(D/s_min). Verificación exhaustiva [N] (chequeo J1): todos los divisores de todos los elementos, en R1 (13 niveles), R2 (426 niveles), R3 (67 niveles) y en la familia escalada z = 13 (682 niveles) y z = 17 (1502 niveles). Sin un solo contraejemplo. La demostración general queda [P-rev]: con división real las condiciones de V⁺ serían independientes de L (resto ⋛ prefijo), y la dependencia del nivel proviene por entero de los suelos ⌊·⌋ de la derivación de Gram; ese es también el motivo de que la propiedad sea delicada. La fórmula de §2.2 **no** depende de esta propiedad (es reindexación pura); la propiedad importa para las estrategias analíticas de acotación (dominar todos los niveles por el nivel mayor).

### 2.2 La fórmula

Con |V⁺(L)| = cardinal y Φ⁺(L) := Σ_{h∈V⁺(L)} h/φ(h):

**|H| = Σ_{s sprime} |V⁺(D/s)| + |V⁺(D/s_min)|**

**W_max = W_1 = Σ_{s sprime} (s/(s−1)) · Φ⁺(D/s) + Φ⁺(D/s_min)**

(para el término s·d se usa φ(s·d) = (s−1)φ(d), válido porque s ∤ d). La fórmula reduce el conteo a los conjuntos V⁺(L), que se calculan por niveles sin enumerar pares: el coste pasa de |H| pares a Σ_s |V⁺(D/s)| elementos con su factorización.

### 2.3 Verificación de la fórmula [N]

| Modelo | \|H\| fórmula | \|H\| enumeración | W_1 fórmula | W_max (BV) | err. rel. |
|---|---|---|---|---|---|
| R1 (z=11, y=60, D=20000) | 54 | 54 | 68,264353 | 68,264353 | 2,1·10⁻¹⁶ |
| R2 (doc 94: y=3000, D=200000) | 880 | 880 | 1011,279168 | 1011,279168 | 5,6·10⁻¹⁶ |
| R3 (z=3, y=349, D=200000) | 68 | 68 | 69,307280 | 69,307280 | 0 |

En la familia escalada (§2.4), una enumeración explícita independiente (par → m → H, sin pasar por la fórmula) reproduce la fórmula exactamente en z = 13, 17, 19, 23 (|H| idéntico; W_1 con error < 10⁻⁹). **No hay discrepancia entre conteo y enumeración en ningún caso probado.**

### 2.4 La familia escalada y el exponente α

La familia es la de tendencia de los docs 87/95: D = z⁵, y = round(z^{10/3}), A = 3, con z = N^{1/10} (equivalentemente D = N^{1/2}, y = N^{1/3}, el marco SW-Q/Chen de los docs 85–92). Valores por la fórmula [N] (chequeo J3):

| z | D = z⁵ | y | #sprimes | \|H\| | W_1 = W_max | envolvente sin V⁺ |
|---|---|---|---|---|---|---|
| 11 | 161051 | 2960 | 422 | 802 | 912,79 | 1692 |
| 13 | 371293 | 5166 | 682 | 1628 | 1898,02 | 4893 |
| 17 | 1419857 | 12633 | 1502 | 5033 | 6015,88 | 17757 |
| 19 | 2476099 | 18303 | 2090 | 8245 | 9943,43 | 36683 |
| 23 | 6436343 | 34601 | 3687 | 17999 | 21933,30 | 94940 |
| 29 | 20511149 | 74931 | 7381 | 47178 | 58126,02 | 274982 |

(Nota: el punto z = 11 de la familia escalada no es R2: D e y difieren — 161051/2960 frente a 200000/3000.)

Regresión log-log (mínimos cuadrados sobre ln z; N = z¹⁰, así α = pendiente/10), con error estándar, IC95% (t de Student) y rango leave-one-out (LOO) [N] (chequeo J4):

- **W_1:** pendiente 4,2886 ± 0,0237 (SE), IC95% ± 0,0658, LOO [4,2650, 4,3203] ⟹ **α_W = 0,4289, IC95% [0,4223, 0,4354], LOO [0,4265, 0,4320]**. Ajuste: W_max(N) ≈ 0,031 · N^{0,4289}. α locales entre z consecutivos: 0,438; 0,430; 0,452; 0,414; 0,420 — sin tendencia clara a la baja; los dos últimos bajan algo, pero el ajuste completo y el LOO son estables por encima de 0,42.
- **|H|:** pendiente 4,2079 ± 0,0170, IC95% ± 0,0472, LOO [4,1969, 4,2279] ⟹ **α_H = 0,4208, IC95% [0,4161, 0,4255], LOO [0,4197, 0,4228]**. Ajuste: |H| ≈ 0,033 · N^{0,4208}. El cociente W_1/|H| crece despacio (1,14 → 1,23 en la ventana): el peso medio h/φ(h) aporta el exceso de α_W sobre α_H.
- **Envolvente sin V⁺** (contar solo d < L, sin las desigualdades de Nathanson): pendiente 5,2439 ⟹ α_env = 0,5244. Las condiciones V⁺ son **esenciales**: rebajan el exponente de 0,52 a 0,43. Un conteo analítico que use solo suavidad (función Ψ de Dickman) no llega ni a 0,50; cualquier cota analítica de |H| debe explotar las desigualdades alternadas, no solo d < L.

**Descomposición del exponente [N]:** el suelo de α viene de #sprimes ~ y/ln y con y = N^{1/3}: eso solo ya da α ≥ 1/3 ≈ 0,333. El exceso medido (+0,096 en W_1, +0,088 en |H|) proviene del crecimiento medio de |V⁺(D/s)| y de Φ⁺(D/s) con el nivel: al crecer D, cada nivel admite más partes nuevas d. Para que α < 0,42 hiciera falta que ese crecimiento medio fuera casi nulo; en los datos crece de forma sostenida en toda la ventana z ∈ [11, 29] (N asociado ∈ [2,6·10¹⁰, 4,2·10¹⁴]).

### 2.5 Extrapolación contra el criterio del doc 100 [N]

El doc 100 tabuló el W_max de equilibrio W_eq(N) (el que haría cerrar la cota BV contra el objetivo de Chen): ≈ 644 (10⁶), ≈ 853 (10¹²), ≈ 1051 (10²⁴), suponiendo W_max congelado en el valor de R2. Con el ajuste medido W_max(N) ≈ 0,031·N^{0,4289}:

- N = 10¹²: W_max ≈ 4350 ≫ 853. **La cota BV no cierra ni siquiera en 10¹²** (falla por ~5×).
- N = 10²⁴: W_max ≈ 6,1·10⁸ ≫ 1051. Falla por ~6·10⁵, y la brecha crece como N^{0,0087}·ln⁵N: empeora con N, no mejora.

El «cruce entre 10¹² y 10²⁴» del doc 100 era un artefacto de congelar W_max: medido su crecimiento, no hay cruce.

---

## 3. E100.1 — El experimento

### 3.1 Maquinaria: ψ segmentado [N]

Para medir E_d con N par hasta 10⁸ y cientos de módulos hace falta ψ(N; d, a) sin materializar Λ hasta N por módulo. El script implementa `SegPsi`: criba segmentada de primos (bloques de 2·10⁶), Λ por potencias de primos, y suma de Λ(n) por clase residual módulo d acumulada bloque a bloque (memoria O(#divs · bloque/d)). Validación (chequeo J0): ψ(999983) reproduce el doc 97 (999586,5975; err. rel. 4,4·10⁻¹²); los E_d de `SegPsi` contra `PsiEngine` del doc 97/100 en N = 999982 sobre los 876 módulos difieren como máximo en 6,0·10⁻¹¹ relativo; la descomposición BV de R2 reproduce el doc 100 (W_max = 1011,2792; S_BV = 4,3593·10⁹; #divs = 880). En la parte A, ψ(N)/N ∈ [0,9996, 1,0000] en los 8 valores (sanidad).

### 3.2 Parte A: soporte fijo R2, N par creciente [N]

Soporte fijo = el H del doc 94 (880 conductores); para cada N par se filtra (r, N) = 1 (eso define H_N y su conjunto de divisores) y se mide con ψ real. Familias pedidas por la condición de Lumina: potencias de 2, primoriales y pares con factores pequeños mezclados (chequeo J5):

| N | factorización / familia | \|H_N\| | #divs | W_max(N) | S_BV | B_M (κ=1) | δ(N) | S_BV/(N^{3/2} ln N) |
|---|---|---|---|---|---|---|---|---|
| 510510 | primorial 2·3·5·7·11·13·17 | 424 | 424 | 424,95 | 9,910·10⁸ | 9,910·10⁸ | 0,3483 | 0,2067 |
| 999982 | 2·79·6329 (79 es sprime) | 876 | 876 | 1006,34 | 4,139·10⁹ | 3,777·10⁹ | 0,3420 | 0,2996 |
| 1048576 | 2²⁰ | 880 | 880 | 1011,28 | 4,813·10⁹ | 4,078·10⁹ | 0,3763 | 0,3234 |
| 9699690 | primorial …·19 | 423 | 423 | 423,90 | 1,955·10¹⁰ | 1,955·10¹⁰ | 0,2962 | 0,0402 |
| 16777216 | 2²⁴ | 880 | 880 | 1011,28 | 6,737·10¹⁰ | 5,692·10¹⁰ | 0,2743 | 0,0589 |
| 67108864 | 2²⁶ | 880 | 880 | 1011,28 | 2,793·10¹¹ | 2,468·10¹¹ | 0,2624 | 0,0282 |
| 77597520 | 2⁴ · primorial 19 | 423 | 423 | 423,90 | 1,525·10¹¹ | 1,525·10¹¹ | 0,2557 | 0,0123 |
| 99999998 | 2·7·23·310559 | 672 | 672 | 736,00 | 3,210·10¹¹ | 2,946·10¹¹ | 0,2594 | 0,0174 |

donde δ(N) := S_BV / (#divs · N · ln N) es el factor medio por divisor (el «δ» del doc 100).

**Lecturas:**

1. **Exponente a soporte fijo:** regresión ln S_BV ~ ln N sobre los 8 puntos: **β = 0,9904 ± 0,0716 (SE), IC95% ± 0,1752, LOO [0,9214, 1,0286]**. Subfamilias puras: 2^e (n=3): β = 0,9729, LOO [0,9517, 1,0258]; primoriales (n=3): β = 1,0031, LOO [0,9878, 1,0127]. Para B_M: β = 0,9904 ± 0,0588, IC95% ± 0,1439. Es decir: **a soporte fijo, S_BV y B_M crecen como N¹ (por N ln N por divisor), no como N^{3/2}**. La columna S_BV/(N^{3/2} ln N) decae monótonamente en tendencia (0,32 → 0,012): la ley N^{3/2} ln N de la tabla corregida del doc 100 describe la **envolvente analítica cuando el soporte crece** (#divs ~ N^{1/2}), no el crecimiento a soporte fijo. Ambas lecturas son compatibles; lo que el experimento fija es que el factor N^{3/2} no debe aplicarse dos veces (una en S_BV y otra vía #divs).
2. **δ(N) decrece despacio:** 0,376 → 0,256 entre 10⁶ y 10⁸. No crece: el factor medio por divisor no es una fuente oculta de crecimiento. (Su valor absoluto depende de la normalización E_d del modelo; lo relevante es la tendencia.)
3. **W_max(N) a soporte fijo nunca supera W_1 del soporte** (1011,28) y el argmax es d = 1 en los 8 valores (chequeo J6, instancia del Lema 1). El crecimiento de W_max con N del criterio BV proviene **solo** del crecimiento del soporte (parte B, §2.4), no de N a soporte fijo.

### 3.3 El colapso primorial [T] + [N]

Para N primorial (divisible por todos los primos pequeños del soporte), el filtro (r, N) = 1 elimina toda r con parte d > 1 (5, 7, 35 dividen a N) y todo sprime que divida a N: H_N = {1} ∪ {s sprime : (s, N) = 1} — **solo conductores primos** (423–424 elementos, y #divs = |H_N|). Para r primo, M_r = φ(r)·E_r exactamente (el término d = 1 lleva E_1 = 0), así que (1/φ(r))|M_r|² = φ(r)|E_r|² y **B_M = S_BV con igualdad exacta** — la desigualdad de Cauchy de la descomposición BV es una identidad sobre conductores primos. Se observa en la tabla: B_M = S_BV al último dígito en los tres N primoriales/mixtos. Corolario de medición: la holgura de Cauchy (70–1145× del doc 100) vive entera en los conductores compuestos.

El mismo mecanismo explica las pérdidas parciales: N = 999982 = 2·79·6329 pierde exactamente las 4 r divisibles por el sprime 79 ({79, 395, 553, 2765}: 880 → 876, W_max 1011,28 → 1006,34); N = 99999998, divisible por 7 y por el sprime 23, pierde 208 (880 → 672, W_max → 736,00). Y W_max cae un 58% en primoriales (1011,28 → 423,90): **N con muchos factores pequeños es el caso favorable** para la cuestión (a) — pero el N de Goldbach es el N dado, no elegible; el peor caso (N = 2^e, soporte intacto) es el que manda.

### 3.4 La ruta combinada con ambos factores medidos [N]

Si en vez de la envolvente analítica S_BV ≪ #divs·N^{3/2}·… se usan **ambas** tasas medidas — S_BV ≈ δ·|H|·N·ln N (δ ≈ 0,26) y W_max, |H| de §2.4 — la cota B_M ≤ W_max·S_BV tiene exponente total 1 + α_W + α_H = **1,8497 < 1,920232** (el del objetivo de Chen, doc 100): el exponente, por primera vez, queda del lado correcto. Pero las constantes y el ln⁵N del objetivo mandan: con los ajustes medidos,

cota/objetivo ≈ 0,0119 · ln⁵N / N^{0,0705},

que vale ≈ 1,2·10⁵ en N = 10²⁴ y solo baja de 1 hacia N ~ 10¹⁵⁰. **La ruta combinada gana solo asintóticamente, a escalas sin sentido práctico ni demostrable con ajustes (que no son cotas).** No es una estrategia de demostración; es una señal de que el exponente no es el muro — el muro son las constantes y los logaritmos, es decir, la pérdida de Cauchy (§3.3: sobre conductores compuestos).

---

## 4. Veredicto

**¿Es viable la ruta BV para la cuestión (a)?** **No**, en la forma en que el doc 100 la formuló (descomposición BV con S_BV acotado analíticamente ~ N^{3/2} ln N y W_max como factor separado). Medido el crecimiento real de W_max, la cota no cierra en ningún N: falla por ~5× en 10¹² y por ~6·10⁵ en 10²⁴, con brecha creciente. La variante combinada con S_BV también medido tiene el exponente correcto (1,850 < 1,920) pero pierde por constantes y ln⁵N hasta N ~ 10¹⁵⁰: tampoco es viable como demostración.

**¿α < 0,42?** **No para W_max**: α_W = 0,4289, IC95% [0,4223, 0,4354], LOO [0,4265, 0,4320] — todo el intervalo por encima de 0,42. Para |H| solo: α_H = 0,4208, IC95% [0,4161, 0,4255] — a caballo; no basta, porque la cota BV carga W_max, no |H|. Confianza: alta **dentro del modelo y de la ventana medida** (6 puntos, z ∈ [11, 29], ajustes estables bajo LOO, fórmula exacta verificada contra enumeración independiente). Es MEDICIÓN, no prueba: no excluye que α(z) local siga bajando a z mayores (los dos últimos α locales son 0,414 y 0,420), pero nada en los datos lo establece, y el criterio del doc 100 exige α < 0,42 sostenido con margen para ln⁵N — el margen medido es **negativo** (−0,009 en α_W).

**¿Qué bloquea, dado α ≥ 0,42?** Dos fuentes, ambas estructurales: (i) #sprimes ~ N^{1/3}/ln N fija un suelo α ≥ 1/3 mientras el marco use y = N^{1/3}; (ii) el crecimiento medio de |V⁺(D/s)| con el nivel aporta el +0,09 restante — acotarlo a O(1) es falso en los datos. Y la salida «contar crudo» está cerrada: sin las desigualdades V⁺ el exponente es 0,52. En resumen: **el factor W_max no es acotable por conteo a la tasa que BV necesita; la pérdida está en separar W_max de S_BV (Cauchy), no en ninguno de los dos factores por separado.**

**Siguiente paso analítico concreto.** Cambiar de objeto: no acotar B_M ≤ W_max · S_BV, sino el segundo momento **directamente** sobre H_A, explotando que la holgura medida de Cauchy es 70–1145× y crece con el soporte (doc 100: R1 ≈ 84×, R3 ≈ 70×, R2 ≈ 1145×) y que sobre conductores primos la holgura es 1 (§3.3). Concretamente:

1. **Experimento E101.1 (propuesto):** medir la holgura W_max·S_BV/B_M a N fijo (N = 999982, ψ ya validado) sobre los soportes escalados z = 11, 13, 17 de §2.4 (|H| = 802, 1628, 5033). Si la holgura crece como |H|^θ con θ apreciable, la cota efectiva gana un factor |H|^{-θ} y el exponente combinado 1,850 baja a 1,850 − θ·0,42: con θ ≈ 0,2 ya habría margen real contra constantes. Es computable en local con `SegPsi`.
2. **Analítico:** usar la representación diagonal exacta de Q_s (docs 96/97, verificada) como forma cuadrática en los E_d y acotar **esa forma** restringida a H_A (criba grande restringida / segundo momento directo), en vez de la factorización W_max × S_BV. Es exactamente la cuestión (b) de Lumina, que pasa a ser la ruta principal para (a).
3. **E101.2 (opcional, descarta la última duda sobre α):** extender la familia escalada a z = 31, 37 con la fórmula de §2.2 (sin enumeración) para comprobar si α_W local sigue bajando o se estabiliza en ≈ 0,43.

Si α hubiera salido < 0,42, el paso habría sido el inverso: convertir la fórmula de §2.2 en cota analítica |H| ≪ N^{α+ε} acotando |V⁺(L)| vía las desigualdades alternadas. Con α medido ≥ 0,42, esa vía queda cerrada salvo que E101.2 la reabra.

---

## 5. Nota de nomenclatura: E100.1-errata

El doc 100 numeró «E100.1» su **errata al doc 97** (el «peso de Δ» 1,24·10⁻⁴ del doc 97 no es comparable con el cambio relativo del doc 99; en soportes completos el peso es ~10⁻⁶ y la coincidencia numérica fue accidental). Ese nombre colisiona con el **experimento E100.1** ejecutado en este doc. Para no confundirlos: en adelante la errata se cita como **E100.1-errata** y «E100.1» a secas es el experimento. Este doc no añade nada a la errata: su contenido quedó establecido en el doc 100 y sigue vigente.

---

## 6. Pendientes [P-rev]

1. Demostración general de que V⁺(L) es cerrado por divisores y monótono en el nivel (verificado exhaustivamente en J1; la fórmula de conteo no depende de ello, las estrategias de acotación por nivel mayor sí).
2. Los 7 puntos [P-rev] del doc 100 §7 contra los textos completos de los docs 96/98/99 y los scripts 94d/94e/97/99 (incluye los 6 del doc 97 §8), cuando el túnel vuelva.
3. Confirmar contra el doc 94 completo que la familia de tendencia (D = z⁵, y = z^{10/3}) es la parametrización asintótica pretendida para R_A; α es una propiedad de esta familia.
4. E101.1 y E101.2 (§4) quedan propuestos, no ejecutados.

---

## Apéndice A. Chequeos del script (todos PASS)

| Chequeo | Contenido | Resultado |
|---|---|---|
| J0 | Maquinaria: \|H\| R2 = 880; ψ(999983) vs doc 97; BV de R2 vs doc 100; SegPsi vs PsiEngine (876 módulos) | PASS (máx. err. rel. E_d: 6,0·10⁻¹¹) |
| J1 | V⁺ cerrado por divisores y monótono: R1/R2/R3 + escalados z = 13, 17, todos los niveles | PASS |
| J2 | Fórmula A100.1 vs enumeración (Model97 + BV) en R1/R2/R3 | PASS (err. rel. ≤ 5,6·10⁻¹⁶) |
| J3 | Familia escalada por fórmula + enumeración explícita independiente en z = 13, 17, 19, 23 | PASS (coincidencia exacta) |
| J4 | Regresiones de α con IC95% y LOO | reporte (§2.4) |
| J5 | Parte A: 8 valores de N par, ψ segmentado, β, δ, colapso primorial | reporte (§3.2–3.3) |
| J6 | W_max = W_1 (argmax d = 1) en los 8 valores de N | PASS |

Reproducibilidad: `python3 101_experimento_e100_1_conteo_a100_1.py` (solo stdlib, ~18 s), salida en `salida_doc101_local.txt`. El script reutiliza sin cambios `Model97`, `PsiEngine` y la descomposición BV de los docs 97/100, y añade `SegPsi`, la fórmula de conteo y las regresiones.

## Apéndice B. Qué no es este doc

No es una demostración de la cuestión (a) en ningún sentido: es la medición que la cierra por la ruta BV ingenua y reubica el margen en la holgura de Cauchy. No usa GRH. No demuestra Goldbach. Los exponentes α y β son ajustes sobre ventanas finitas, citados siempre con su intervalo; tratarlos como cotas sería exactamente el error que el doc 100 corrigió en su propia tabla (E100.1-errata y la corrección N² → N^{3/2}).
