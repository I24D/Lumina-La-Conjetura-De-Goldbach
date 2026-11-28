# DOC 105 — La candidata diagonal de segundos momentos B̂₅ pasa (i) y C1, y pasa C2 en la lectura por piezas (la honesta tras el doc 104); el precio lo pone el conteo de divisores, no la cancelación

**Serie:** verificación y continuación de los docs 94–104 (Gram de Lumina, diagonalización, puente de Möbius, segundo factor, conteo A100.1, holgura de Cauchy, criba restringida, mayorantes absolutos).
**Fecha:** 2026-11-28. **Autora:** Lumina Muse.
**Archivos:** `105_forma_diagonal_segundos_momentos.py` (solo stdlib; corre en local, **PASS** en 1,6 s), `salida_doc105_local.txt` (salida completa).
**Estado:** PARCIAL en el mismo sentido que los docs 95/97/100–104: los textos completos de los docs 98/99 y los scripts 94d/94e/97/99 de Lumina siguen sin estar disponibles en local (túnel caído); este doc ejecuta la redirección fijada en el doc 104 §5.5 sobre las definiciones del modelo finito de los docs 94–97 y 100–104, que sí están completos en local. Los puntos [P-rev] heredados y los nuevos se listan en §6.
**Disciplina:** se distingue en todo momento MEDICIÓN de PRUEBA. Nada de lo aquí medido demuestra Goldbach; la conjetura fuerte sigue abierta. Sin circularidad, sin GRH sin probar.

**Etiquetas:** [T] = teorema (demostrado en el texto o en un doc anterior citado), [N] = medición numérica reproducible con el script adjunto, [P-rev] = pendiente de revisión contra el texto completo.

---

## 0. Alcance

El doc 104 cerró la línea de mayorantes absolutos en la base E (B̂₂ y B̂₃ pasan C1 y fallan C2, y su cierre exigiría primeros momentos absolutos ponderados de E_d, sin herramienta clásica) y fijó la redirección: trabajar B_M en su **forma diagonal** con **solo segundos momentos**, con criterio previo doble — (i) expresabilidad solo en segundos momentos de E_d sobre 𝒟 (sin primeros momentos absolutos, sin sumas de caracteres sobre ℋ) y (ii) C1+C2 en los seis puntos — antes de formular lema alguno. Este doc construye esa candidata (B̂₅, §1), verifica el criterio (i) término a término, decide con argumento — antes de los números, §2 — qué lectura de C2 es la honesta tras el hallazgo de mezcla del doc 104, evalúa (§3), descompone el precio (§4) y dictamina con el enunciado exacto del lema a probar (§5).

**Notación (docs 94–104).** H = soporte (cierre por divisores de los m = s·d); H_N = {r ∈ H : (r, N) = 1}; ℋ = {r ∈ H_N : r > R0} con R0 = (ln N)^κ, κ = 1; 𝒟 = {d ≥ 2 : d | r para algún r ∈ ℋ}; E_d = error de ψ en progresiones (ψ real, `SegPsi`); F_d := φ(d)E_d (F_1 = 0); M_r = Σ_{d|r} μ(r/d)F_d; S_r = M_r + Δ_r(N) con Δ_r la corrección imprimitiva de los docs 97–100 (identidad corregida [T]); B_M = Σ_{r∈ℋ} (1/φ(r))·|S_r|², con piezas P (r primo) y C (r compuesto) y formas (M) (solo M_r) y (S) (con Δ_r). Los seis puntos del criterio (doc 102): N = 999982 y N = 2²⁰, soportes S11, S13, S17. Criterio: **(C1)** B̂/B_M ≤ 10 en los seis puntos; **(C2)** pendiente log-log de B̂/B_M contra |H_N| ≤ 0,1 (lectura a decidir en §2).

---

## 1. La candidata diagonal B̂₅ y el criterio (i)

### 1.1 La forma diagonal de B_M^(M) [T]

Expandiendo |M_r|² y sumando sobre ℋ, B_M^(M) es una forma cuadrática en el vector F = (F_d)_{d∈𝒟}:

**B_M^(M) = Σ_{d,e ∈ 𝒟} F_d F_e W_{de}**, con **W_{de} := Σ_{r∈ℋ, d|r, e|r} μ(r/d)μ(r/e)/φ(r)**.

**Factorización [T].** Todo elemento de H es libre de cuadrados (m = s·d con d libre de cuadrados y s primo distinto). Para d|r, e|r y L := mcm(d, e): r/d = (r/L)·(L/d) con factores coprimos (r/L es coprimo con L), luego μ(r/d) = μ(r/L)μ(L/d), análogamente μ(r/e), y μ(r/L)² = 1. Por tanto el sumando no depende de r salvo por 1/φ(r):

**W_{de} = μ(L/d)·μ(L/e)·U_L**, donde **U_L := Σ_{r∈ℋ, L|r} 1/φ(r) ≥ 0**,

y en particular **|W_{de}| = U_{mcm(d,e)}** (0 si ningún r ∈ ℋ es múltiplo del mcm). La diagonal es W_{dd} = U_d, y Diag = Σ_d F_d²U_d (la diagonal de los docs 102–104, aquí recuperada como la diagonal de W). Las tres identidades — forma cuadrática, factorización con signo, y Diag por suma en d = Diag por suma en r — se verifican numéricamente a precisión de máquina en R1 y R3 (chequeo K2: PASS, peor error relativo 0,00 en la factorización).

### 1.2 El mayorante por filas: AM–GM sobre W, sin tocar los signos de F [T]

Para cada par (d, e), |F_dF_eW_{de}| ≤ |W_{de}|·(F_d² + F_e²)/2. Sumando sobre todos los pares y usando la simetría |W_{de}| = |W_{ed}| (visible en la factorización), cada F_d² recoge el coeficiente Σ_e |W_{de}|:

**B_M^(M) ≤ D̂ := Σ_{d∈𝒟} φ(d)²E_d²·Ω_d**, con **Ω_d := Σ_{e∈𝒟} |W_{de}|**.

**Forma cerrada de Ω_d [T, conteo].** Σ_{e∈𝒟} U_{mcm(d,e)} = Σ_{r∈ℋ} (1/φ(r))·#{e ∈ 𝒟 : mcm(d, e) | r}. Si d ∤ r el conjunto es vacío; si d | r, mcm(d, e) | r ⟺ e | r, y #{e ≥ 2 : e | r} = τ(r) − 1 (τ = función divisor). Luego

**Ω_d = Σ_{r∈ℋ, d|r} (τ(r) − 1)/φ(r)**,

verificado contra las sumas por fila directas de la matriz W en R1 y R3 (K2: PASS). El mayorante no paga (Σ_d |F_d|)² como B̂₂, ni sumas de caracteres como B̂₁: paga, por cada d, el peso U_d multiplicado por un factor de conteo de divisores — el «precio τ» de §4.

**Por piezas [T].** W = W^P + W^C según r primo o compuesto, y la misma construcción da D̂^X ≥ B_M^{X,(M)} para X ∈ {P, C} con Ω_d^X := Σ_{r∈ℋ_X, d|r} (τ(r)−1)/φ(r). La pieza prima es **exacta**: un r primo solo tiene un divisor ≥ 2 (él mismo), W^P es diagonal en 𝒟, Ω_d^P = [d primo ∈ ℋ]/φ(d) y **D̂^P = Diag^P = B_M^{P,(M)}** (verificado a rel < 10⁻⁹ en todo punto, K4: PASS) — la candidata respeta el cierre por identidad de los primos, como B̂₃ y a diferencia de B̂₁.

**Δ_r, por Minkowski por pieza [T].** Con B_Δ^X := Σ_{r∈ℋ_X} Δ_r²/φ(r) (segundo momento de la corrección explícita):

B_M^{X,(S)} = Σ_{r∈ℋ_X} |M_r + Δ_r|²/φ(r) ≤ (√B_M^{X,(M)} + √B_Δ^X)² ≤ (√D̂^X + √B_Δ^X)² =: **B̂₅^X**,

y **B̂₅ := B̂₅^P + B̂₅^C ≥ B_M** (K3: PASS en todo punto, incluidos R1 y R3). Medido: B_Δ/B_M^(S) ~ 10⁻⁶ en los seis puntos [N], así que el plegado de Δ_r cuesta menos que la tercera cifra de cada cociente — pero está en la desigualdad, no despreciado (la exigencia del doc 102 §126.3: conservar Δ_r).

### 1.3 Verificación formal del criterio (i), término a término

| Término de B̂₅ | Naturaleza | ¿Segundo momento? |
|---|---|---|
| D̂^X = Σ_d φ(d)²E_d²Ω_d^X | combinación positiva de E_d² (d ∈ 𝒟) con pesos Ω_d^X | **Sí** — segundos momentos de E_d, con pesos aritméticos explícitos que dependen solo del soporte (sumas de (τ(r)−1)/φ(r) sobre r ∈ ℋ_X): cantidades de conteo, controlables por la maquinaria de A100.1 (doc 101) |
| B_Δ^X = Σ_{r∈ℋ_X} Δ_r²/φ(r) | segundo momento sobre el soporte de la corrección imprimitiva Δ_r | **Sí, con declaración**: no es un momento de E_d, es el segundo momento de una **función explícita** (suma sobre potencias de primos p^j ≤ N con p | r, docs 97–100). Es el único ingrediente no-E de la candidata; es elemental y clásico (no requiere BV ni criba), y su tamaño medido es ~10⁻⁶ del momento total [N]. Se declara, no se esconde |
| (√D̂^X + √B_Δ^X)² | composición de los dos anteriores (Minkowski) | No introduce ningún objeto nuevo |

**Ningún término usa** primeros momentos absolutos Σ φ(d)|E_d|·u_d (lo que mató el cierre de B̂₂/B̂₃, doc 104 §5.4b), ni sumas de caracteres sobre ℋ (lo que mató a B̂₁ por δ_ℋ, doc 103), ni triángulo dentro de la base E. El criterio (i) del doc 104 queda **satisfecho formalmente**. Lo que la candidata sí necesita — y no se regala — es acotar D̂^C: un segundo momento ponderado de E sobre 𝒟. Eso es exactamente la clase de objeto que la tecnología tipo BV/BDH controla (§5.3), a diferencia de los primeros momentos absolutos.

---

## 2. Qué lectura de C2 es la honesta (argumento, antes de los números)

El doc 104 demostró que la pendiente agregada de B̂₂/B̂₃ es un **efecto mezcla**: el cociente agregado es el promedio ponderado por B_M de los cocientes por pieza (1,000 en P y ~2,7 en C para B̂₃), y los pesos se desplazan (share_C 0,42 → 0,69) al crecer el soporte. Tres hechos estructurales, todos [T] o verificados por punto (K5):

1. El cociente agregado de cualquier candidata separable por piezas es el promedio ponderado de sus cocientes por pieza con pesos B_M^{X}/B_M; en particular **cociente agregado ≤ max(cociente^P, cociente^C) en cada punto** — el agregado no puede divergir si las piezas no divergen; su «pendiente» en una ventana finita mide composición, no degradación por conductor.
2. La función de C2 en la serie es detectar **divergencia** del precio con el soporte (B̂₁: pendiente 2,37 por δ_ℋ⁻¹, divergencia real; BV: ~1,0). Una candidata con precios por pieza constantes q^P ≠ q^C para siempre «fallaría» la lectura agregada literal en cuanto la composición se moviese — es decir, la lectura literal rechaza candidatas uniformemente dentro de un factor constante del objetivo, lo que contradice la función del criterio.
3. El objetivo final de Chen se aplica sobre el total, pero el total se controla por (precio de la peor pieza) × (masa total): si ambas piezas tienen precio acotado con pendiente ≤ 0,1 y el nivel agregado cumple C1, el balance sobre el total hereda el precio de la peor pieza, no una divergencia. El riesgo real sería que share_C → 1 **y** el precio^C creciera; lo primero es inofensivo si lo segundo no ocurre, y lo segundo es exactamente C2 sobre la pieza C.

**Decisión:** la lectura honesta de C2 es **por piezas** (pendientes de B̂^P/B_M^P y B̂^C/B_M^C ≤ 0,1), **junto con** C1 en agregado y la cota puntual del punto 1. La lectura agregada literal se reporta igualmente, completa, en §3 — y falla — pero queda reclasificada como lo que el doc 104 demostró que es: un artefacto de composición. Como corroboración independiente se computa además el **contrafactual de composición fija**: la pendiente del cociente agregado con share_C congelado en su media de los seis puntos (0,556); si el argumento de mezcla es correcto, esa pendiente debe coincidir con el crecimiento intra-pieza y quedar ≤ 0,1. Este párrafo se escribió con el criterio fijado por el análisis del doc 104, no a la vista de los resultados de B̂₅.

---

## 3. Evaluación en los seis puntos [N]

| Punto | \|H_N\| | B_M^(S) | B̂₅/B_M^(S) | B̂₅^P/B_M^{P,(S)} | B̂₅^C/B_M^{C,(S)} | D̂^C/B_M^{C,(M)} | Diag^C/B_M^{C,(M)} | precio τ^C | share_C |
|---|---|---|---|---|---|---|---|---|---|
| gen/S11 | 798 | 3,3658·10⁹ | 3,190 | 1,001 | 6,163 | 6,158 | 1,495 | 4,119 | 0,424 |
| gen/S13 | 1 621 | 7,6433·10⁹ | 4,115 | 1,000 | 6,614 | 6,605 | 1,488 | 4,439 | 0,555 |
| gen/S17 | 5 017 | 2,7875·10¹⁰ | 5,499 | 1,000 | 7,532 | 7,525 | 1,533 | 4,909 | 0,689 |
| 2^e/S11 | 802 | 3,7105·10⁹ | 3,600 | 1,001 | 7,014 | 7,007 | 1,674 | 4,185 | 0,432 |
| 2^e/S13 | 1 628 | 8,3332·10⁹ | 4,470 | 1,000 | 7,411 | 7,403 | 1,616 | 4,582 | 0,541 |
| 2^e/S17 | 5 033 | 2,9761·10¹⁰ | 5,619 | 1,000 | 7,667 | 7,661 | 1,576 | 4,862 | 0,693 |

(B_M^(S) y Diag^C/B_M^C reproducen los docs 103/104, chequeo K1: PASS. Calibración fuera de criterio: gen/R1 B̂₅/B_M = 7,919 — también ≤ 10, con share_C = 0,932 —; gen/R3: 1,001, soporte de solo primos.)

- **C1: PASS** — B̂₅/B_M^(S) ∈ [3,190, 5,619] en los seis puntos, contra el umbral 10. Es el primer pase de C1 de una candidata **expresable en segundos momentos**: B̂₂/B̂₃ pagaban menos nivel (1,71–2,32) pero sin ruta de cierre; B̂₅ paga un factor ~2–2,5 más de nivel a cambio de que cada uno de sus términos sea un segundo momento (§4 descompone ese sobreprecio: es el precio τ).
- **C2 agregado (lectura literal): FAIL** — pendiente **0,2658 ± 0,0311 (SE), IC95% ± 0,0864, LOO [0,2349, 0,2823]**, entero sobre 0,1. Se reporta sin maquillar: bajo la redacción literal del doc 102, B̂₅ no pasa C2. La §2 argumenta por qué esa lectura dejó de ser la honesta tras el doc 104, y el contrafactual de abajo lo corrobora numéricamente.
- **C2 por piezas (lectura defendida en §2): PASS en estimación puntual, con margen delgado en C** —
  - Pieza P: pendiente **0,0000** en forma (M) (exacta [T], K4) y **−0,0001** en forma (S): precio 1,000–1,001 en los seis puntos.
  - Pieza C: pendiente **0,0785 ± 0,0330 (SE), IC95% ± 0,0916, LOO [0,0518, 0,1028]** en forma (M); **0,0784** en forma (S). La estimación puntual queda bajo el umbral 0,1, pero el intervalo LOO **roza el umbral por arriba en su extremo** (0,1028) y el IC95% lo contiene: con seis puntos, la pendiente de la pieza C no se puede certificar < 0,1 con la misma fuerza con que el doc 104 certificó la de B̂₃ > 0,1 (LOO entero sobre el umbral). Dictamen preciso: **pasa por estimación puntual; la certificación intervalar queda abierta** y depende del término identificado en §4.
- **Contrafactual de composición fija: pendiente 0,0703 ± 0,0295, LOO [0,0465, 0,0921] — entero bajo 0,1.** Con share_C congelado en 0,556, el cociente agregado de B̂₅ crece con pendiente 0,070: el crecimiento intra-pieza es el de la pieza C (0,0785) atenuado por la pieza P exacta, y el resto de la pendiente agregada literal (0,2658 − ~0,07 ≈ 0,20) es el desplazamiento de composición share_C 0,42 → 0,69 actuando sobre la brecha de precios q^C − q^P ≈ 6,5. La descomposición es coherente con el mecanismo del doc 104 y con el punto 1 de §2 (el cociente agregado ≤ max(q^P, q^C) = q^C en los seis puntos, K5: PASS).

---

## 4. Dónde se concentra el precio: el conteo de divisores, no la cancelación [N]

El precio de la pieza C se factoriza exactamente (aditividad de pendientes por la linealidad de OLS, como en el doc 104 §3):

D̂^C/B_M^{C,(M)} = (Diag^C/B_M^{C,(M)}) × (D̂^C/Diag^C) = **(1,49–1,67) × (4,12–4,91)**,
pendientes: **0,0785 = (−0,0085) + 0,0870**.

- **El factor Diag (la cancelación de Möbius) es plano y barato:** 1,49–1,67 con pendiente −0,0085 ± 0,0275, LOO [−0,0288, 0,0159] — compatible con cero, como ya midieron los docs 102–104. La candidata diagonal **captura** esa cancelación en vez de pagarla dos veces: D̂ parte de Diag, no de la suma de valores absolutos.
- **El precio τ (D̂^C/Diag^C = promedio ponderado de Ω_d^C/U_d^C) es el precio dominante y el único que crece:** 4,12–4,91 con pendiente ≈ 0,087. Es un precio de **conteo**, no de análisis: Ω_d/U_d promedia (τ(r) − 1) sobre los r ∈ ℋ_C divisibles por d, y los conductores acumulan divisores al crecer el soporte (doc 104: divisores por conductor 2,38 → 3,58 entre S11 y S17). AM–GM por filas paga, por cada par (d, e), la posibilidad de alineación total de signos entre F_d y F_e a través de todos los conductores comunes; el factor (τ(r) − 1) es el número de socios e que cada r ofrece a cada uno de sus divisores.
- **Consecuencia para el cierre:** el término que hay que vigilar — y el que decidirá si la pendiente de la pieza C se queda bajo 0,1 fuera de la ventana — no es analítico sino aritmético: el crecimiento de los promedios de τ(r) − 1 sobre ℋ. Por construcción del soporte, ω(r) ≤ 1 + #{primos permitidos < z} (r | s·d con un único primo ≥ z y d formado por primos < z), así que τ(r) está acotado por una cantidad de conteo sobre V⁺ — el ingrediente 2 del lema (§5.3) es una cota de tipo A100.1 sobre esos promedios, no una conjetura.

---

## 5. Veredicto

1. **La candidata diagonal B̂₅ satisface el criterio (i) formalmente [T]:** cada término es un segundo momento — de E_d con pesos explícitos del soporte (D̂^X), o de la corrección explícita Δ_r (B_Δ^X, declarado en §1.3 como el único ingrediente no-E, elemental y de tamaño medido ~10⁻⁶). Nada de primeros momentos absolutos, nada de caracteres sobre ℋ, nada de triángulo en la base E. Es la primera candidata de la cuestión (b) que pasa el filtro de expresabilidad que el doc 104 impuso.
2. **C1 pasa [N]:** B̂₅/B_M ∈ [3,190, 5,619] en los seis puntos (umbral 10), y 7,919 en R1 fuera de criterio. El nivel es ~2–2,5× peor que el de B̂₂/B̂₃ — ese es el coste medido de la expresabilidad en segundos momentos, y lo paga entero el precio τ (§4). El margen del doc 102 (B_M quedaría 10³–10⁴× bajo el objetivo de Chen a las tasas medidas) absorbe un factor ≤ 5,7 sin cambiar de orden.
3. **C2 pasa en la lectura por piezas y falla en la literal agregada, y ambas cosas quedan dichas:** agregado 0,2658 (LOO entero sobre 0,1); por piezas P = 0,0000 exacto y C = 0,0785 (LOO [0,0518, 0,1028], rozando el umbral); contrafactual de composición fija 0,0703 (LOO entero bajo 0,1). La decisión de §2 — la lectura por piezas es la honesta, porque el agregado es un promedio ponderado acotado por la peor pieza (K5) y su pendiente mide composición (doc 104), no degradación — se sostiene con el contrafactual. **La salvedad que queda:** la pendiente de la pieza C no está certificada intervalarmente bajo 0,1; su motor es el precio τ, un objeto de conteo (§4), y afinarla pide soportes mayores (S19/S23), no otra candidata.
4. **Se formula el lema a probar** (la condición del encargo y del doc 104: formularlo solo si la candidata pasa (i) + C1 + C2 en la versión defendida):

   **Lema L105.** Sea el soporte H_A construido como en la serie (A fijo; s primo en [z, y); d ∈ V⁺(D/s); H = cierre por divisores), N par, ℋ = {r ∈ H_N : r > (ln N)^κ}, 𝒟 su conjunto de divisores ≥ 2, y T(N) el objetivo de Chen para el segundo factor (la cota sobre B_M que el puente de Cauchy de los docs 94/96 requiere para que el resto primitivo sea o(N/(A log²N))).
   (a) **[TEOREMA — demostrado en §1 de este doc]** B_M ≤ Σ_{X∈{P,C}} (√D̂^X + √B_Δ^X)², con D̂^X y B_Δ^X como en §1.2; y D̂^P = B_M^{P,(M)} exactamente.
   (b) **[A PROBAR]** Existen constantes absolutas C_Ω y C_Δ, independientes de N y del soporte dentro de la familia, tales que
   **D̂^C = Σ_{d∈𝒟} φ(d)²E_d²Ω_d^C ≤ C_Ω·T(N)** y **B_Δ^P + B_Δ^C ≤ C_Δ·T(N)**.
   Conclusión del lema: B_M ≤ (√C_Ω + √C_Δ + o(1))²·T(N) — el segundo factor queda acotado a la escala que el balance de Chen exige, con el precio medido en este doc (≤ 5,7 en la ventana) como candidato a las constantes.
5. **Ingredientes analíticos que (b) requiere** (la lista pedida; todos clásicos o de conteo, ninguno es un primer momento absoluto):
   1. **Segundo momento ponderado de E sobre 𝒟, tipo Barban–Davenport–Halberstam / criba grande en forma dual:** Σ_{d∈𝒟} w_d φ(d)²E(N; d, N)² con w_d = Ω_d^C. Es la suma sobre módulos de errores al cuadrado en una clase por módulo — la clase de objeto que BV/BDH controlan — y **no** la forma de espaciamiento de Farey que murió en el doc 103 (δ_ℋ era el precio de sumar |ψ(N, χ)|² sobre caracteres; aquí no hay caracteres). El trabajo real: BDH suma sobre todas las clases (a, d); extraer la clase diagonal a ≡ N con pesos φ(d)² y w_d del soporte es la adaptación a hacer, al nivel de distribución de 𝒟.
   2. **Cotas de conteo para los pesos Ω_d^C (tipo A100.1, doc 101):** promedios de τ(r) − 1 y sumas Σ_{r∈ℋ_C, d|r} 1/φ(r) sobre el soporte; en concreto, que el precio τ efectivo, Ω_d^C/U_d^C promediado con los pesos φ(d)²E_d²U_d^C, quede acotado por una constante (medido: 4,12–4,91, pendiente 0,087 en la ventana). Por construcción ω(r) ≤ 1 + #{primos permitidos < z}, así que es conteo sobre V⁺, no análisis.
   3. **Segundo momento de Δ_r:** elemental — Δ_r es una suma explícita sobre potencias de primos (docs 97–100); su momento sobre ℋ es un ejercicio de conteo (medido: ~10⁻⁶ del momento total).
   4. **Nada para la pieza P:** cerrada por identidad [T].
   La composición (share_C) **no** es un ingrediente: por el punto 1 de §2, el total hereda el precio de la peor pieza vía promedio ponderado; no hay que controlar cómo se reparte B_M entre piezas, solo los precios por pieza.
6. **Lo que este doc no hace:** no demuestra (b) — formula su enunciado exacto y verifica que la candidata que lo alimenta pasa los filtros numéricos y de expresabilidad de la serie. Si el ingrediente 1 o el 2 fallara al intentarlo, el término que rompe ya está identificado (el precio τ y su pendiente 0,087, §4), y la variante siguiente sería un AM–GM **ponderado** por filas (pesos u_{de} optimizados por bloque de mcm, explotando que |W_{de}| = U_{mcm(d,e)} depende solo del mcm) con el mismo criterio previo (i) + C1 + C2 por piezas. **Goldbach sigue abierta.**

---

## 6. Pendientes [P-rev]

1. Los heredados del doc 104 §6: textos completos de los docs 96/98/99 y scripts 94d/94e/97/99; demostración general de cierre por divisores de V⁺; la extrapolación de la holgura a N = 10²⁴; max_{d≥2} W_d sin medir por separado; la lectura L2 de nueva(d) contra el doc 94 íntegro.
2. **La pendiente de la pieza C de B̂₅ (0,0785, LOO hasta 0,1028) no está certificada bajo 0,1.** Afinarla pide medir B̂₅ en soportes S19/S23 (coste de `SegPsi` sobre conjuntos de divisores mayores) o una cota de conteo directa del precio τ (ingrediente 2 de §5.5 hecho número). Es el chequeo que decide si el pase de C2 por piezas es robusto o un accidente de ventana.
3. El ingrediente 1 de §5.5 (BDH ponderado con la clase diagonal a ≡ N y pesos del soporte) está identificado como tecnología, no verificado contra la literatura en detalle: la adaptación exacta (nivel de distribución alcanzado sobre 𝒟, dependencia de los pesos Ω en z) queda como primer paso del intento de demostración de L105(b), no como afirmación de este doc.
4. La candidata usa AM–GM simétrico por filas; la variante ponderada de §5.6 (pesos por bloque de mcm) no está medida. Si el precio τ resultara ser el obstáculo del ingrediente 2, esa variante es la siguiente en la cola, con su criterio previo ya fijado.

---

## Apéndice A. Chequeos del script (todos PASS)

| Chequeo | Contenido | Resultado |
|---|---|---|
| K0 | Maquinaria: Model97 R2 \|H\| = 880; build_H z = 11 == Model97 | PASS |
| K1 | B_M^(S) reproduce los docs 103/104 en los seis puntos (rel < 2·10⁻³); Diag^C/B_M^C(S) reproduce 1,488–1,674 (rel < 3·10⁻²) | PASS |
| K2 | En R1 y R3: forma cuadrática Σ F_dF_eW_{de} = B_M^(M) (rel < 10⁻⁹); W_{de} = μ(L/d)μ(L/e)U_L con peor error relativo 0,00; sumas por fila Σ_e\|W_{de}\| = Ω_d (fórmula τ); Diag por suma en d = Diag por suma en r | PASS |
| K3 | Mayorante válido: D̂^X ≥ B_M^{X,(M)}, B̂₅^X ≥ B_M^{X,(S)}, B̂₅ ≥ B_M^(S) en todo punto (seis + R1 + R3) | PASS |
| K4 | Pieza prima exacta: D̂^P = B_M^{P,(M)} (rel < 10⁻⁹) en todo punto | PASS |
| K5 | Cociente agregado = promedio ponderado de los cocientes por pieza y ≤ max(piezas) en los seis puntos | PASS |
| K6 | C1+C2 para B̂₅: C1 PASS (3,190–5,619 ≤ 10); C2 agregado FAIL (0,2658); C2 por piezas PASS en estimación puntual (P 0,0000; C 0,0785, LOO [0,0518, 0,1028]) | reporte (§3) |

Reproducibilidad: `python3 105_forma_diagonal_segundos_momentos.py` (solo stdlib, ~1,6 s), salida en `salida_doc105_local.txt`. El script reutiliza sin cambios la maquinaria de los docs 102–104 (importada del script 104: `Model97`, `build_H`, `SegPsi`, `M_r_E`, `ols`, `delta_r_orth`) y añade la construcción de los pesos Ω^X por la fórmula τ, el mayorante D̂ por piezas, el plegado de Δ_r por Minkowski, la verificación directa de la matriz W en R1/R3 y el contrafactual de composición fija.

## Apéndice B. Qué no es este doc

No es una demostración de la cuestión (a), de la cuestión (b) ni de Goldbach: es la construcción y evaluación de la primera candidata de la ruta diagonal que satisface el criterio de expresabilidad del doc 104, con su lema (L105) formulado y sus ingredientes identificados. La desigualdad B_M ≤ B̂₅ sí es un teorema [T] dentro del modelo (AM–GM + Minkowski sobre identidades exactas, verificado en K2–K4); lo que queda abierto es acotar B̂₅ a la escala T(N), que es L105(b). No usa GRH. Los cocientes, pendientes y precios son mediciones sobre ventanas finitas, citados con su intervalo cuando procede; la factorización de W, la forma cerrada de Ω_d y la exactitud de la pieza prima son implicaciones rigurosas [T].
