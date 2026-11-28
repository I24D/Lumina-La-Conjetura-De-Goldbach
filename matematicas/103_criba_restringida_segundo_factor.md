# DOC 103 — Criba grande restringida a H_A sobre B_M: la candidata de Montgomery–Vaughan muere por el espaciamiento δ_ℋ⁻¹; el mayorante diagonal absoluto pasa C1 y roza C2

**Serie:** verificación y continuación de los docs 94–102 (Gram de Lumina, diagonalización, puente de Möbius, segundo factor, conteo A100.1, holgura de Cauchy).
**Fecha:** 2026-11-28. **Autora:** Lumina Muse.
**Archivos:** `103_criba_restringida_segundo_factor.py` (solo stdlib; corre en local, **PASS** en 2,1 s), `salida_doc103_local.txt` (salida completa).
**Estado:** PARCIAL en el mismo sentido que los docs 95/97/100/101/102: los textos completos de los docs 98/99 y los scripts 94d/94e/97/99 de Lumina siguen sin estar disponibles en local (túnel caído); este doc trabaja sobre las definiciones del modelo finito tal como quedaron fijadas en los docs 94–97, 100, 101 y 102, que sí están completos en local. Los puntos [P-rev] heredados se listan en §6.
**Disciplina:** se distingue en todo momento MEDICIÓN de PRUEBA. Nada de lo aquí medido demuestra Goldbach; la conjetura fuerte sigue abierta. Sin circularidad, sin GRH sin probar.

**Etiquetas:** [T] = teorema (demostrado en el texto o en un doc anterior citado), [N] = medición numérica reproducible con el script adjunto, [P-rev] = pendiente de revisión contra el texto completo.

---

## 0. Alcance

El doc 102 fijó el siguiente paso analítico: formular y tratar de probar la cota de criba grande **restringida** a H_A sobre el segundo factor, y dejó un criterio numérico previo de alcanzabilidad: una candidata B̂ está en la vía correcta si y solo si, en los seis puntos (N = 999982 y N = 2²⁰; soportes S11, S13, S17),

- **(C1)** B̂/B_M ≤ 10 en los seis puntos, y
- **(C2)** la pendiente log-log de B̂/B_M contra |H_N| es ≤ 0,1.

Este doc formula la candidata pedida — la desigualdad de la criba grande en forma de Montgomery–Vaughan, restringida al soporte H_A, aplicada a la forma diagonal (en caracteres) de B_M con Δ_r conservado —, la evalúa en los seis puntos con un sándwich riguroso para su constante de espaciamiento, la descompone en primos y compuestos, e identifica el término que la rompe. Como la candidata **no pasa** C1+C2, se propone y se pre-evalúa la variante siguiente (mayorante diagonal absoluto, V103.1) y se formula la refinación que le sigue (V103.2) con su propio criterio.

**Notación (docs 94–102).** H = soporte (cierre por divisores de los m = s·d); H_N = {r ∈ H : (r, N) = 1}; ℋ = {r ∈ H_N : r > R0} con R0 = (ln N)^κ, κ = 1; E_d = error de ψ en progresiones (ψ real, `SegPsi`); M_r = Σ_{d|r} μ(r/d)φ(d)E_d; S_r = M_r + Δ_r(N) con Δ_r la corrección imprimitiva de los docs 97–100 (identidad corregida [T]); el objeto a acotar es el momento corregido

B_M = Σ_{r∈ℋ} (1/φ(r))·|S_r|².

En los seis puntos, B_M difiere del B_M del doc 102 (definido con M_r solo) en un cambio relativo (B_cross + B_Δ)/B_M ∈ [−5,0·10⁻⁵, +8,4·10⁻⁵] [N, chequeo K2]: la corrección Δ_r está conservada en todo lo que sigue y no mueve ninguna cifra a la precisión relevante.

---

## 1. La candidata B̂₁: criba grande de Montgomery–Vaughan restringida a ℋ

### 1.1 Derivación, con todas las constantes [T]

**Paso 1 (Cauchy por conductor, constante 1).** S_r se define como suma sobre caracteres primitivos, S_r = Σ\*_{χ mod r} χ̄(N)ψ(N, χ). Cauchy con pesos 1 da |S_r|² ≤ φ\*(r)·Σ\*_{χ mod r} |ψ(N, χ)|², donde φ\*(r) = #{χ primitivos mod r} ≤ φ(r). Luego

(1/φ(r))·|S_r|² ≤ Σ\*_{χ mod r} |ψ(N, χ)|²,

con constante exactamente 1. Δ_r entra entero dentro de S_r: no se separa M_r de Δ_r (separarlos pagaría el término cruzado, doc 100: |B_cross| ≫ B_Δ).

**Paso 2 (Montgomery–Vaughan restringida).** La criba grande en su forma de espaciamiento, aplicada solo a los módulos de ℋ (no a todos los módulos ≤ Q):

Σ_{r∈ℋ} Σ\*_{χ mod r} |ψ(N, χ)|² ≤ (N + δ_ℋ⁻¹)·Λ₂(N) =: **B̂₁**,

donde Λ₂(N) = Σ_{n≤N} Λ(n)² y δ_ℋ es el espaciamiento mínimo del conjunto de fracciones ℱ_ℋ = {a/r mod 1 : r ∈ ℋ, (a, r) = 1}. La constante que multiplica a (N + δ_ℋ⁻¹) es 1. Λ₂(N) se mide exactamente en el modelo: Λ₂(999982) = 1,280280·10⁷, Λ₂(2²⁰) = 1,347870·10⁷ (es decir Λ₂ ≈ N·ln N con constante medida 12,80 y 12,86 frente a ln N = 13,82 y 13,86) [N].

**La candidata completa:** B_M ≤ B̂₁ = (N + δ_ℋ⁻¹)·Λ₂(N). Todo el contenido de la restricción a H_A está en δ_ℋ: si las fracciones de ℋ estuvieran mejor espaciadas que las de Farey completo, δ_ℋ⁻¹ ≪ Q² y la restricción ganaría. La pregunta es cuánto vale δ_ℋ.

### 1.2 El sándwich de δ_ℋ [T]

No hace falta conocer δ_ℋ exactamente para dictaminar C1+C2; basta acotarlo por ambos lados con exhibidos computables:

- **Cota inferior de δ_ℋ:** dos fracciones reducidas distintas a/r, a'/r' con r, r' ≤ Q_ℋ := max ℋ distan |ar' − a'r|/(rr') ≥ 1/Q_ℋ². Luego δ_ℋ ≥ 1/Q_ℋ² y B̂₁ ≤ (N + Q_ℋ²)·Λ₂(N) =: B̂_hi.
- **Cota superior de δ_ℋ (exhibido gemelo):** todo elemento de H es impar (m = s·d con s, d impares), así que ℋ no contiene enteros consecutivos; pero si r, r+2 ∈ ℋ, las fracciones 1/r y 1/(r+2) ∈ ℱ_ℋ distan exactamente 2/(r(r+2)). Con r₂ := max{r ∈ ℋ : r+2 ∈ ℋ}: δ_ℋ ≤ 2/(r₂(r₂+2)), luego δ_ℋ⁻¹ ≥ r₂(r₂+2)/2 y

**B̂_lo := (N + r₂(r₂+2)/2)·Λ₂(N) ≤ B̂₁ ≤ (N + Q_ℋ²)·Λ₂(N) =: B̂_hi.**

Consecuencia lógica: si B̂_lo/B_M > 10 en un punto, C1 falla en ese punto **para la candidata verdadera**, sea cual sea δ_ℋ. Y C2 se evalúa sobre B̂_lo, la versión más favorable de la familia. El veredicto de §2 es, por construcción, conservador a favor de la candidata.

**Calibración con δ_ℋ exacto (R1, R3; enumeración completa de fracciones) [N, chequeo K3]:**

| Punto | #fracciones | δ_ℋ exacto | δ_ℋ⁻¹ | sándwich [1/Q_ℋ², 2/(r₂(r₂+2))] | B̂(δ exacto)/B_M | B̂_lo/B_M |
|---|---|---|---|---|---|---|
| gen/R1 | 12 612 | 9,137·10⁻⁶ | 1,094·10⁵ | [2,906·10⁻⁷, 4,287·10⁻⁵] | 1,72·10⁵ | 1,59·10⁵ |
| gen/R3 | 10 356 | 8,551·10⁻⁶ | 1,169·10⁵ | [8,305·10⁻⁶, 2,055·10⁻⁵] | 7,53·10⁴ | 7,07·10⁴ |

Dos lecturas: (i) el sándwich contiene al δ exacto en ambos casos y el extremo favorable B̂_lo subestima el cociente verdadero solo por un factor 1,07–1,08 — el veredicto con B̂_lo no es un artefacto del sándwich; (ii) en R3 (soporte de solo primos) δ_ℋ está a un 3% del suelo absoluto 1/Q_ℋ²: el conjunto de fracciones de un soporte «delgado» puede ser casi tan denso como Farey completo a su escala.

---

## 2. Evaluación en los seis puntos: C1 y C2 fallan [N]

| Punto | \|H_N\| | Q_ℋ | r₂ | δ_ℋ⁻¹ ≥ | B_M | B̂_lo/B_M | B̂_hi/B_M |
|---|---|---|---|---|---|---|---|
| gen/S11 | 798 | 16 345 | 2 933 | 4,30·10⁶ | 3,3658·10⁹ | 2,02·10⁴ | 1,02·10⁶ |
| gen/S13 | 1 621 | 37 415 | 10 703 | 5,73·10⁷ | 7,6433·10⁹ | 9,76·10⁴ | 2,35·10⁶ |
| gen/S17 | 5 017 | 144 445 | 81 653 | 3,33·10⁹ | 2,7875·10¹⁰ | 1,53·10⁶ | 9,58·10⁶ |
| 2^e/S11 | 802 | 16 345 | 2 933 | 4,30·10⁶ | 3,7105·10⁹ | 1,94·10⁴ | 9,74·10⁵ |
| 2^e/S13 | 1 628 | 37 415 | 10 703 | 5,73·10⁷ | 8,3332·10⁹ | 9,44·10⁴ | 2,27·10⁶ |
| 2^e/S17 | 5 033 | 144 445 | 81 653 | 3,33·10⁹ | 2,9761·10¹⁰ | 1,51·10⁶ | 9,45·10⁶ |

- **C1: FAIL en los seis puntos, por ≥ 4 órdenes de magnitud.** Incluso el extremo favorable del sándwich da B̂_lo/B_M ∈ [1,94·10⁴, 1,53·10⁶], contra el umbral 10. La candidata verdadera está entre B̂_lo y B̂_hi, es decir, entre 10⁴ y 10⁷ veces el momento que pretende acotar.
- **C2: FAIL.** Pendiente log-log de B̂_lo/B_M contra |H_N| (seis puntos, OLS): **2,3701 ± 0,0322 (SE), IC95% ± 0,0895, LOO [2,3518, 2,3950]** — contra el umbral 0,1. El cociente crece como |H_N|^2,37: la candidata se degrada con el soporte más rápido que cualquier candidata anterior de la serie (el producto BV del doc 102 tenía pendiente ≈ 1,0).

**Diagnóstico en una frase:** la restricción a H_A no compra nada en la forma de Montgomery–Vaughan, porque lo que la criba grande cobra no es el número de módulos sino su espaciamiento de Farey, y H_A — delgado en cardinalidad (|H| ~ N^0,42, doc 101) — es **denso en Farey en la parte alta de su rango**: contiene pares gemelos (r, r+2) hasta r₂ = 81 653 en S17, y eso fija δ_ℋ⁻¹ ≳ r₂²/2 ~ Q_ℋ²/6, del mismo orden que el Q² de la criba sin restringir.

---

## 3. Descomposición primo/compuesto [N]

Aplicando la misma construcción por pieza (para primos, el exhibido es el mejor par de primos consecutivos en ℋ, δ⁻¹ ≥ pp'/(p'−p); para compuestos, el mayor gemelo compuesto):

| Punto | B̂_lo^P/B_M^P | B̂_lo^C/B_M^C | cociente C/P | Diag^C/B_M^C (M) | Diag^C/B_M^C (S) |
|---|---|---|---|---|---|
| gen/S11 | 3,25·10⁴ | 4,76·10⁴ | 1,5 | 1,495 | 1,495 |
| gen/S13 | 5,27·10⁴ | 1,76·10⁵ | 3,3 | 1,488 | 1,488 |
| gen/S17 | 1,19·10⁵ | 2,22·10⁶ | 18,7 | 1,533 | 1,533 |
| 2^e/S11 | 3,18·10⁴ | 4,50·10⁴ | 1,4 | 1,674 | 1,674 |
| 2^e/S13 | 4,95·10⁴ | 1,74·10⁵ | 3,5 | 1,616 | 1,615 |
| 2^e/S17 | 1,19·10⁵ | 2,18·10⁶ | 18,4 | 1,576 | 1,576 |

**Lecturas:**

1. **El fracaso no se concentra en los compuestos: es uniforme.** La pieza prima — donde el doc 102 cerró el problema con una identidad exacta y donde ninguna criba debería perder nada — recibe de B̂₁ un cociente de 3·10⁴–1,2·10⁵. La candidata no distingue la estructura que los docs 101/102 aislaron; pierde por espaciamiento en todas partes. La pieza compuesta sí degrada más rápido (el cociente C/P crece de ~1,4 a ~19 con el soporte, porque los gemelos compuestos llegan más arriba que los pares de primos útiles), pero es una diferencia de grado sobre un fracaso común de 4–6 órdenes.
2. **La cancelación Diag^C/B_M^C ≈ 1,5 se mantiene intacta con el momento corregido.** Recomputada sobre S_r = M_r + Δ_r, vale lo mismo que sobre M_r a las tres cifras impresas (1,488–1,674 según el punto; chequeo K1 contra el doc 102: PASS). La corrección Δ_r no toca la cancelación de los cruzados de Möbius. Pero B̂₁ **no la ve**: la candidata salta de S_r a caracteres y nunca pasa por la base E donde la cancelación vive (Diag es una forma en los E_d; la criba de caracteres es ciega a ella). Cualquier cota que aspire a C1+C2 debe formularse en la base donde Diag^C/B_M^C ≈ 1,5 es visible, no en la base de caracteres.

---

## 4. El término que rompe, y la variante siguiente

### 4.1 El término que rompe: δ_ℋ⁻¹

En los seis puntos, el término de espaciamiento domina la candidata: δ_ℋ⁻¹ ≥ r₂(r₂+2)/2 vale 4,30·10⁶ (S11), 5,73·10⁷ (S13) y 3,33·10⁹ (S17), frente a N ≈ 10⁶ — es decir, el término «Q²» de la criba supera al término «N» por factores 4,3, 57 y 3 300, y crece con el soporte porque r₂ crece con el rango de ℋ (pendiente medida del cociente total: 2,37). El término Λ₂(N) ≈ N ln N es inocente (constante medida ~12,8, sin crecimiento anómalo); el paso de Cauchy por conductor es constante 1. **Toda la pérdida está en δ_ℋ⁻¹, y es estructural: es el precio de tratar los caracteres de ℋ como si sus fracciones pudieran estar tan mal espaciadas como las de Farey — y en ℋ, arriba, efectivamente lo están (gemelos).** Ninguna restricción adicional del soporte que conserve los conductores grandes de H_A puede arreglar δ_ℋ, porque los gemelos (r, r+2) son conductores legítimos del soporte, no ruido eliminable.

### 4.2 Variante V103.1 — mayorante diagonal absoluto: pasa C1, roza C2 [N]

La variante que el diagnóstico sugiere es abandonar la base de caracteres y mayorar en la base E, pagando deliberadamente toda la cancelación para medir su precio. Por la desigualdad triangular en la suma de Möbius y en Δ_r [T]:

|S_r| ≤ Σ_{d|r} φ(d)|E_d| + |Δ_r| ⟹ **B̂₂ := Σ_{r∈ℋ} (1/φ(r))·(Σ_{d|r} φ(d)|E_d| + |Δ_r|)² ≥ B_M** (chequeo K4: PASS en todo punto).

Evaluación con el mismo criterio C1+C2 (es el criterio previo propio de la variante):

| Punto | B̂₂/B_M | | Punto | B̂₂/B_M |
|---|---|---|---|---|
| gen/S11 | 1,709 | | 2^e/S11 | 1,826 |
| gen/S13 | 1,920 | | 2^e/S13 | 2,053 |
| gen/S17 | 2,231 | | 2^e/S17 | 2,318 |
| (gen/R1) | 2,934 | | (gen/R3) | 1,001 |

- **C1: PASS** — B̂₂/B_M ∈ [1,71, 2,32] en los seis puntos, contra el umbral 10. Matar toda la cancelación cuesta solo un factor ~2: la cancelación total (cruzados de Möbius + cruzado con Δ) vale menos que el factor 1,5 de Diag sugería por separado, porque B̂₂ también pierde el signo dentro de cada |E_d| individual… y aun así queda en ~2.
- **C2: FAIL, estrecho.** Pendiente log-log: **0,1358 ± 0,0201 (SE), LOO [0,1198, 0,1484]** — el intervalo LOO entero queda sobre el umbral 0,1, pero por poco (la pendiente es ~7× menor que la de BV y ~17× menor que la de B̂₁). El precio de la cancelación crece despacio con el soporte, como |H_N|^0,14: al crecer el soporte, los conductores acumulan más divisores y el triángulo por conductor pierde un poco más.

**Estado de V103.1:** es un mayorante válido y numéricamente excelente como *forma*, pero (i) no pasa C2, y (ii) su cierre analítico exigiría acotar una suma de **primeros momentos absolutos** Σ_d φ(d)|E_d|·u_d sobre los divisores del soporte — un objeto para el que la criba grande clásica tampoco sirve (cobra el mismo δ⁻¹ si se pasa por caracteres; y |E_d| no tiene la cancelación entre residuos que BDH explota). V103.1 es el banco de pruebas, no la cota.

### 4.3 Variante V103.2 — triángulo solo entre clases de parte nueva (propuesta, con criterio)

La pendiente 0,14 de V103.1 es la señal a atacar: hay que dejar de pagar el triángulo **dentro** de los grupos donde la serie ya sabe diagonalizar con signo. Propuesta concreta:

**B̂₃ := Σ_{r∈ℋ} (1/φ(r))·( Σ_{h} | Σ_{d|r, nueva(d)=h} μ(r/d)φ(d)E_d | + |Δ_r| )²**,

donde nueva(d) es la parte nueva h = d/(d, A)-reducida de los docs 94/96 (la misma partición que diagonaliza el lado de coeficientes: Q_s = Σ φ(r)C_s(r)²). El triángulo se aplica solo *entre* clases h; dentro de cada clase la suma conserva los signos de Möbius y los E_d con signo. B̂₃ ≤ B̂₂ siempre, y B̂₃ ≥ B_M [T, mismo argumento].

**Criterio numérico previo (antes de intentar cualquier cierre analítico):** el mismo C1+C2 del doc 102 sobre los seis puntos, aplicado a B̂₃. **Hipótesis falsable que lo motiva:** si la pendiente 0,14 de V103.1 proviene de la cancelación entre divisores de clases distintas (que crece con el número de clases por conductor), B̂₃ la aplana por debajo de 0,1 conservando C1; si proviene de cancelación *dentro* de cada clase, B̂₃ heredará la pendiente y la vía del mayorante absoluto queda cerrada del todo. **Experimento E103.1 (propuesto):** computar B̂₃ en los seis puntos y dictaminar. Es computable en local con la maquinaria existente (la partición por parte nueva ya está en `Model97` vía Hsmooth/T).

---

## 5. Veredicto

1. **La criba grande de Montgomery–Vaughan restringida a H_A (B̂₁) está muerta, y el motivo es nuevo respecto a los docs 100–102:** no es el conteo del soporte ni la holgura de Cauchy — es el **espaciamiento de Farey** δ_ℋ. H_A es delgado en cardinalidad pero denso en Farey en su parte alta (gemelos (r, r+2) hasta 81 653), así que δ_ℋ⁻¹ ~ Q_ℋ²/2 y restringir los módulos no restringe la constante. C1 falla por ≥ 4 órdenes en los seis puntos (extremo favorable del sándwich: 1,9·10⁴–1,5·10⁶); C2 falla con pendiente 2,370 ± 0,090.
2. **El fracaso es uniforme por piezas:** la pieza prima recibe 3·10⁴–1,2·10⁵ pese a estar cerrada por identidad; la compuesta degrada más rápido (cociente C/P hasta ~19) pero sobre el mismo fracaso. La cancelación Diag^C/B_M^C ≈ 1,5 **se mantiene idéntica** con Δ_r conservado (1,488–1,674, = valores del doc 102), y es invisible para cualquier candidata en base de caracteres.
3. **El mayorante diagonal absoluto B̂₂ (V103.1) es la primera forma de esta línea que pasa C1** (1,71–2,32 en los seis puntos): matar toda la cancelación cuesta ~2×, no 10⁴×. Pero **falla C2 por poco** (pendiente 0,136, LOO [0,120, 0,148]) y su cierre analítico pediría primeros momentos absolutos de E_d, igualmente fuera del alcance de la criba clásica.
4. **La siguiente variante es B̂₃ (V103.2):** triángulo solo entre clases de parte nueva, signos conservados dentro de cada clase vía la diagonalización conocida. Criterio previo: C1+C2 en los seis puntos (experimento E103.1). La hipótesis es falsable y el experimento es barato.
5. **No se formula lema a probar en este doc:** la condición del encargo era formularlo solo si la candidata pasaba C1+C2. B̂₁ no pasa; V103.1 pasa C1 pero no C2; el enunciado exacto queda condicionado al resultado de E103.1 sobre B̂₃. **Goldbach sigue abierta.** Lo que este doc cierra: la vía «criba grande restringida en base de caracteres» completa (clásica en el doc 100, restringida aquí). Lo que abre: la familia de mayorantes diagonales con signo parcial (V103.2), con el precio de la cancelación ya medido (~2× y pendiente 0,14 a eliminar).

---

## 6. Pendientes [P-rev]

1. Los heredados del doc 102 §5: textos completos de los docs 96/98/99 y scripts 94d/94e/97/99; demostración general de cierre por divisores de V⁺; la extrapolación de la holgura a N = 10²⁴; max_{d≥2} W_d sin medir por separado.
2. δ_ℋ exacto en los soportes grandes (S11–S17) no es computable por enumeración de fracciones (Σ_{r∈ℋ} φ(r) ~ 10⁷–10⁹); el sándwich de §1.2 es el sustituto riguroso y la calibración de R1/R3 fija su holgura en un factor ≤ 1,08 sobre el cociente. Un algoritmo de espaciamiento sin enumeración total (tipo Farey por intervalos) queda como mejora técnica posible, no necesaria para el veredicto.
3. La forma MV usada es la clásica (constante 1, término N + δ⁻¹). Variantes con pesos de Selberg *en el dual de caracteres* no se han evaluado: el diagnóstico de §4.1 (δ_ℋ es propiedad del soporte, no de los pesos sobre n) sugiere que tampoco escapan del espaciamiento, pero no está medido; queda subsumido por la propuesta V103.2, que cambia de base en vez de ponderar la misma.
4. E103.1 (computar B̂₃ y dictaminar C1+C2) queda propuesto, no ejecutado.
5. El cierre analítico de V103.1/V103.2 (primeros momentos absolutos de E_d sobre 𝒟, o su versión por clases) no tiene aún herramienta identificada en la serie; es el análogo en la base E del problema de espaciamiento en la base de caracteres.

---

## Apéndice A. Chequeos del script (todos PASS)

| Chequeo | Contenido | Resultado |
|---|---|---|
| K0 | Maquinaria: Model97 R2 \|H\| = 880; build_H z = 11 == Model97 | PASS |
| K1 | B_M^(M) reproduce la tabla genérica del doc 102 (R1 8,244·10⁷; R3 1,898·10⁸; S11 3,366·10⁹; S13 7,644·10⁹; S17 2,788·10¹⁰; rel < 2·10⁻³) y Diag^C/B_M^C reproduce 1,488–1,674 en genérica y 2^e (rel < 3·10⁻²) | PASS |
| K2 | B_M^(S) = B_M^(M) + B_cross + B_Δ exacto por construcción; cambio relativo en los seis puntos ∈ [−5,0·10⁻⁵, +8,4·10⁻⁵] | PASS |
| K3 | Sándwich de δ_ℋ en R1/R3: 1/Q_ℋ² ≤ δ exacto ≤ 2/(r₂(r₂+2)) | PASS |
| K4 | B̂₂ ≥ B_M^(S) en todo punto (mayorante válido) | PASS |
| K5 | C1+C2 para B̂_lo: C1 FAIL (mín 1,94·10⁴), C2 FAIL (2,370); para B̂₂: C1 PASS (1,71–2,32), C2 FAIL (0,136) | reporte (§2, §4.2) |

Reproducibilidad: `python3 103_criba_restringida_segundo_factor.py` (solo stdlib, ~2,1 s), salida en `salida_doc103_local.txt`. El script reutiliza sin cambios `Model97`, `build_H`, `SegPsi`, `M_r_E` y `ols` de los docs 97/100/101/102 (importados del script 102), copia `P_e`/`delta_r_orth` del doc 100, y añade la medición del momento corregido, el sándwich de espaciamiento con δ exacto en soportes pequeños, los exhibidos por pieza y el mayorante B̂₂.

## Apéndice B. Qué no es este doc

No es una demostración de la cuestión (a) ni de Goldbach: es la ejecución de la prueba numérica previa que el doc 102 exigió antes de intentar la demostración de la criba restringida, y su resultado es negativo para la forma en base de caracteres y mixto (C1 sí, C2 no) para el mayorante diagonal absoluto. No usa GRH. Los cocientes y pendientes son mediciones sobre ventanas finitas, citados con su intervalo cuando procede; el sándwich de δ_ℋ, en cambio, es una implicación rigurosa [T]: el fracaso de C1 para B̂₁ no depende de la ventana, solo de que ℋ contenga los gemelos exhibidos.
