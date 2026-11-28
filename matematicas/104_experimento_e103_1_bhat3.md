# DOC 104 — Experimento E103.1: B̂₃ (triángulo solo entre clases de parte nueva) pasa C1 y hereda la pendiente de B̂₂; la línea de mayorantes absolutos en la base E queda agotada

**Serie:** verificación y continuación de los docs 94–103 (Gram de Lumina, diagonalización, puente de Möbius, segundo factor, conteo A100.1, holgura de Cauchy, criba restringida).
**Fecha:** 2026-11-28. **Autora:** Lumina Muse.
**Archivos:** `104_experimento_e103_1_bhat3.py` (solo stdlib; corre en local, **PASS** en 1,6 s), `salida_doc104_local.txt` (salida completa).
**Estado:** PARCIAL en el mismo sentido que los docs 95/97/100–103: los textos completos de los docs 98/99 y los scripts 94d/94e/97/99 de Lumina siguen sin estar disponibles en local (túnel caído); este doc ejecuta el experimento E103.1 propuesto en el doc 103 sobre las definiciones del modelo finito fijadas en los docs 94–97 y 100–103, que sí están completos en local. Los puntos [P-rev] heredados y los nuevos se listan en §6.
**Disciplina:** se distingue en todo momento MEDICIÓN de PRUEBA. Nada de lo aquí medido demuestra Goldbach; la conjetura fuerte sigue abierta. Sin circularidad, sin GRH sin probar.

**Etiquetas:** [T] = teorema (demostrado en el texto o en un doc anterior citado), [N] = medición numérica reproducible con el script adjunto, [P-rev] = pendiente de revisión contra el texto completo.

---

## 0. Alcance

El doc 103 cerró la criba grande restringida en base de caracteres (B̂₁ muere por el espaciamiento δ_ℋ) y pre-evaluó la variante V103.1: el mayorante diagonal absoluto B̂₂ (triángulo en toda la suma de Möbius y en Δ_r) **pasa C1** (B̂₂/B_M ∈ [1,71, 2,32] en los seis puntos) pero **falla C2 por poco** (pendiente log-log 0,1358 ± 0,0201, LOO [0,1198, 0,1484], umbral 0,1). Su §4.3 propuso la variante V103.2 — B̂₃, triángulo solo *entre* clases de parte nueva, signos conservados *dentro* de cada clase — y dejó el experimento **E103.1** con una hipótesis falsable sobre el origen de la pendiente 0,14. Este doc define B̂₃ con precisión (declarando la ambigüedad de la definición y eligiendo la lectura fiel), lo evalúa contra C1+C2 en los seis puntos, somete la hipótesis a prueba con una descomposición exacta de la pendiente, localiza dónde se concentra el exceso, y dictamina sobre la línea de mayorantes.

**Notación (docs 94–103).** H = soporte (cierre por divisores de los m = s·d); H_N = {r ∈ H : (r, N) = 1}; ℋ = {r ∈ H_N : r > R0} con R0 = (ln N)^κ, κ = 1; E_d = error de ψ en progresiones (ψ real, `SegPsi`); M_r = Σ_{d|r} μ(r/d)φ(d)E_d; S_r = M_r + Δ_r(N) con Δ_r la corrección imprimitiva de los docs 97–100 (identidad corregida [T]); el objeto a acotar es B_M = Σ_{r∈ℋ} (1/φ(r))·|S_r|². Los seis puntos del criterio (doc 102): N = 999982 y N = 2²⁰, soportes S11, S13, S17. Criterio: **(C1)** B̂/B_M ≤ 10 en los seis puntos; **(C2)** pendiente log-log de B̂/B_M contra |H_N| ≤ 0,1.

---

## 1. Definición de B̂₃ y la ambigüedad de «parte nueva»

### 1.1 La fórmula del doc 103

El doc 103 §4.3 define:

**B̂₃ := Σ_{r∈ℋ} (1/φ(r))·( Σ_{h} | Σ_{d|r, nueva(d)=h} μ(r/d)φ(d)E_d | + |Δ_r| )²**,

donde nueva(d) es «la parte nueva h = d/(d, A)-reducida de los docs 94/96 (la misma partición que diagonaliza el lado de coeficientes: Q_s = Σ φ(r)C_s(r)²)». El triángulo se aplica solo entre clases h; dentro de cada clase la suma conserva los signos de Möbius y los E_d con signo; Δ_r entra triangulado por separado, como en B̂₂. Que B̂₃ ≥ B_M (triángulo entre clases sobre una partición de la suma M_r) y B̂₃ ≤ B̂₂ (refina el triángulo total) vale para **cualquier** partición [T]; ambas desigualdades se verifican por conductor en el script (chequeo K2: PASS en los ocho puntos medidos, incluidos R1 y R3).

### 1.2 Dos lecturas de nueva(d); se adopta L2 y se declara la ambigüedad

La glosa «d/(d, A)-reducida» admite dos lecturas, y el experimento depende de cuál se tome:

- **Lectura literal (L1):** nueva(d) = d/(d, A) con A = 3, el módulo fijo del modelo. En este modelo todo elemento de H es coprimo con 3 (los m = s·d se construyen con s primo ≥ z coprimo con A y d formado por primos 2 < p < z coprimos con A; el primo 3 nunca aparece). Luego para todo d | r con r ∈ ℋ, (d, 3) = 1 y **nueva(d) = d: cada clase es unitaria y B̂₃ = B̂₂ idénticamente** [T; verificado computacionalmente en los ocho puntos, chequeo K3: PASS]. Bajo L1, E103.1 está decidido de forma trivial (C1 pasa, C2 falla con la pendiente de B̂₂) y la hipótesis del doc 103 no es falsable. Además, bajo L1 sería falsa la afirmación del propio doc 103 de que la partición «ya está en `Model97` vía Hsmooth/T»: nada en `Model97` computa d/(d, 3).
- **Lectura operativa (L2, ADOPTADA):** la partición que efectivamente diagonaliza el lado de coeficientes en la maquinaria existente (`Model97`, docs 94/96/97) es la descomposición por la **parte lisa en z**: cada módulo se escribe con su parte h0 ∈ `Hsmooth` (primos < z) y a lo sumo un primo de conmutación s ≥ z; los coeficientes se agrupan por h0 (T_s(h0) = Σ_{d∈V_s, h0|d} μ(d)/φ(d), Q_s = Σ_{h0} φ(h0)T_s(h0)²). Trasladada a la suma M_r, la clase de un divisor d | r es su parte lisa:

  **nueva_z(d) := ∏_{p^e ‖ d, p < z} p^e**, con z el parámetro de construcción del soporte (11, 13, 17 en los seis puntos; 11 en R1; 3 en R3).

  Como r | m y cada m tiene un único primo ≥ z (su s, con multiplicidad 1) y r es libre de cuadrados, cada clase contiene a lo sumo dos divisores: **{h0, h0·s}** con s el primo ≥ z de r (si h0·s | r). Los dos términos de una clase tienen signos de Möbius **opuestos** [T]: si s | r, s ∤ h0, entonces μ(r/(h0·s)) = −μ(r/h0) y φ(h0·s) = φ(h0)(s−1), luego la suma de la clase es μ(r/h0)·φ(h0)·[E_{h0} − (s−1)·E_{h0·s}]. Es exactamente el par cuya cancelación la diagonalización de coeficientes explota en el lado K_A. Esta es la lectura que el doc 103 indica al remitir a `Hsmooth`/`T`, y la única que hace el experimento falsable. **Toda la §2–§5 usa L2.** Si el texto completo del doc 94 (hoy inaccesible, [P-rev]) definiera nueva(d) de otro modo, B̂₃ debe recomputarse bajo esa definición; la dependencia queda registrada en §6.

**Estructura medida de las clases (L2) [N]:** en los seis puntos, el 99,8–99,9% de los r ∈ ℋ tienen alguna clase de tamaño 2; clases por conductor: 1,69 (S11) → 2,30 (S17); divisores por conductor: 2,38 → 3,58; prácticamente toda la masa de Σ_h|suma de clase| está en clases-par (99,5–100%). En R3 (soporte de solo primos, z = 3) cada r tiene una sola clase efectiva y B̂₃ = B̂₂ (cociente 1,001, como el doc 103 obtuvo para B̂₂): consistencia en el caso degenerado.

---

## 2. Evaluación en los seis puntos: C1 pasa, C2 falla [N]

| Punto | \|H_N\| | B_M | B̂₂/B_M (doc 103) | B̂₃/B_M | B̂₂/B̂₃ |
|---|---|---|---|---|---|
| gen/S11 | 798 | 3,3658·10⁹ | 1,709 | 1,679 | 1,018 |
| gen/S13 | 1 621 | 7,6433·10⁹ | 1,920 | 1,876 | 1,023 |
| gen/S17 | 5 017 | 2,7875·10¹⁰ | 2,231 | 2,181 | 1,023 |
| 2^e/S11 | 802 | 3,7105·10⁹ | 1,826 | 1,798 | 1,016 |
| 2^e/S13 | 1 628 | 8,3332·10⁹ | 2,053 | 2,000 | 1,026 |
| 2^e/S17 | 5 033 | 2,9761·10¹⁰ | 2,318 | 2,262 | 1,025 |

(B_M y B̂₂/B_M reproducen el doc 103 con rel < 2·10⁻³ y < 10⁻² respectivamente, chequeo K1: PASS; calibración: gen/R1 B̂₃/B_M = 2,620 frente a B̂₂/B_M = 2,934; gen/R3: 1,001 = B̂₂.)

- **C1: PASS** — B̂₃/B_M ∈ [1,679, 2,262] en los seis puntos, contra el umbral 10. Conservar los signos dentro de las clases mejora el nivel de B̂₂ solo por un factor 1,016–1,026.
- **C2: FAIL** — Pendiente log-log de B̂₃/B_M contra |H_N| (seis puntos, OLS): **0,1326 ± 0,0195 (SE), IC95% ± 0,0542, LOO [0,1168, 0,1441]**. El intervalo LOO entero queda sobre el umbral 0,1, igual que con B̂₂ (0,1358, LOO [0,1198, 0,1484]). B̂₃ reduce la pendiente en **0,0032**, un 2,4% de la pendiente original: la refinación diseñada para aplanarla la deja intacta a todos los efectos.

**Consecuencia según la regla de la serie:** B̂₃ no pasa C1+C2, luego **no se formula lema a probar** para V103.2 (la condición del encargo y del doc 103 era formularlo solo en caso de pase).

---

## 3. La hipótesis falsable del doc 103, sometida a prueba [N, con descomposición exacta]

La hipótesis del doc 103 §4.3, en su forma predictiva: *si la pendiente 0,14 de B̂₂ proviene de la cancelación entre divisores de clases distintas (que crece con el número de clases por conductor), B̂₃ la aplana por debajo de 0,1 conservando C1; si proviene de cancelación dentro de cada clase, B̂₃ heredará la pendiente y la vía del mayorante absoluto queda cerrada del todo.*

**Descomposición exacta.** Como OLS es lineal en la variable dependiente y log(B̂₂/B_M) = log(B̂₂/B̂₃) + log(B̂₃/B_M) punto a punto, las pendientes se suman sin error (verificado a 10⁻⁹, chequeo K4: PASS):

pendiente(B̂₂/B_M) = pendiente(B̂₂/B̂₃) + pendiente(B̂₃/B_M)
**0,1358 = 0,0032 + 0,1326**,

donde B̂₂/B̂₃ es el precio del triángulo **dentro** de clase (lo único que B̂₃ deja de pagar respecto a B̂₂) y B̂₃/B_M es el precio del triángulo **entre** clases más la separación de Δ_r (lo que B̂₃ sigue pagando). La pendiente del primer sumando es 0,0032 ± 0,0018 (SE), IC95% ± 0,0050, LOO [0,0018, 0,0043]: despreciable y casi plana.

**Por qué dentro de clase no hay casi nada que recuperar [T + N].** La suma de una clase-par es μ(r/h0)φ(h0)[E_{h0} − (s−1)E_{h0·s}] (§1.2): el término en d = h0·s domina al de d = h0 por el factor de peso (s−1) ≥ z−1 ≥ 10, y E_{h0}, E_{h0·s} son errores en módulos distintos, sin relación sistemática que produzca cancelación. El cociente medido ρ_intra := Σ_r Σ_h |suma de clase| / Σ_r Σ_d |término| vale **0,987–0,992** en los seis puntos (0,933 en R1, donde B̂₂/B̂₃ llega a 1,120): conservar los signos dentro de clase retiene ~99% de la suma absoluta. La estructura de conmutación fuerza clases desequilibradas — φ(h0·s)/φ(h0) = s−1 — y una clase desequilibrada no puede cancelarse a sí misma.

**Dictamen sobre la hipótesis.** En su forma predictiva, **refutada**: B̂₃ no aplana la pendiente por debajo de 0,1 (0,1326, LOO entero sobre 0,1); la reducción conseguida (0,0032) es un orden de magnitud menor que la necesaria (≥ 0,036). El mecanismo de la primera rama tampoco es el medido: la pendiente no residía en la cancelación que el firmado intra-clase recupera (esa vale 0,0032 de pendiente y un factor ~1,02 de nivel). Y la segunda rama acierta en el resultado — B̂₃ **hereda** la pendiente — pero no en el mecanismo que enunciaba: no hay cancelación intra-clase apreciable que heredar (ρ_intra ≈ 0,99); lo heredado es el precio del triángulo *entre* las sumas de clase, que B̂₃ sigue pagando entero. La dicotomía del doc 103 era incompleta; la descomposición de arriba la sustituye. Por la regla de decisión del propio doc 103 («si B̂₃ hereda la pendiente, la vía del mayorante absoluto queda cerrada del todo»), la herencia medida activa el cierre, que §5 fundamenta.

---

## 4. Dónde se concentra el exceso: piezas planas y efecto mezcla [N]

Descomponiendo B̂₃ por piezas (P = r primo, C = r compuesto), como en el doc 103 §3:

| Punto | share_C(B_M) | B̂₃^P/B_M^P | B̂₃^C/B_M^C |
|---|---|---|---|
| gen/S11 | 0,424 | 1,000 | 2,600 |
| gen/S13 | 0,555 | 1,000 | 2,579 |
| gen/S17 | 0,689 | 1,000 | 2,714 |
| 2^e/S11 | 0,432 | 1,000 | 2,846 |
| 2^e/S13 | 0,541 | 1,000 | 2,847 |
| 2^e/S17 | 0,693 | 1,000 | 2,822 |

- **La pieza prima es exacta:** B̂₃^P/B_M^P = 1,000 en los seis puntos (pendiente −0,0001): para r primo, M_r = φ(r)E_r tiene un único término no nulo (E_1 = 0) y no hay triángulo que pague, salvo la separación de Δ_r, invisible a tres cifras. Como en los docs 102/103, el problema de los primos está cerrado por identidad y B̂₃ lo respeta.
- **Todo el exceso está en los compuestos, a precio por pieza constante:** B̂₃^C/B_M^C ∈ [2,58, 2,85] con pendiente **0,0104 ± 0,0268 (SE), LOO [−0,0059, 0,0299]** — compatible con cero. Por conductor compuesto, B̂₃ paga un factor ~2,6–2,85 que **no crece** con el soporte en la ventana medida.
- **La pendiente total 0,1326 es un efecto de mezcla (composición), no de crecimiento por conductor:** al crecer |H_N|, la fracción de B_M que cargan los compuestos sube de 0,42 a 0,69 en ambas familias (la masa del segundo factor se desplaza hacia conductores compuestos, coherente con el doc 102: la holgura vive en los compuestos), y el cociente total — promedio ponderado de 1,000 (P) y ~2,7 (C) — sube arrastrado por los pesos. El cociente B̂₃/B_M está por tanto acotado por el precio de la pieza C (~2,9) en esta familia de soportes: su «crecimiento» es de saturación, no de divergencia — a diferencia del fracaso de B̂₁ (pendiente 2,370, divergente por δ_ℋ⁻¹). Esto no salva C2 (el criterio se aplica tal cual sobre los seis puntos y falla), pero cambia el diagnóstico: B̂₃ no se degrada por conductor; se degrada porque el objeto mismo se vuelve más compuesto.

---

## 5. Veredicto

1. **B̂₃ (V103.2, lectura L2) pasa C1 y falla C2:** B̂₃/B_M ∈ [1,679, 2,262] (umbral 10); pendiente 0,1326 ± 0,0195, LOO [0,1168, 0,1441] (umbral 0,1). No se formula lema a probar. Bajo la lectura literal L1, B̂₃ = B̂₂ [T] y el resultado es el del doc 103: mismo dictamen.
2. **La hipótesis falsable del doc 103 queda refutada en su forma predictiva** (aplanamiento < 0,1): la descomposición exacta 0,1358 = 0,0032 + 0,1326 localiza el 97,6% de la pendiente en el factor que B̂₃ sigue pagando (triángulo entre sumas de clase), y solo el 2,4% en el que deja de pagar. Dentro de clase no hay cancelación explotable (ρ_intra ≈ 0,99) porque las clases-par están forzadas a ser desequilibradas por el factor (s−1) ≥ 10 [T]. La dicotomía del doc 103 se corrige: herencia de la pendiente, sí; mecanismo intra-clase, no.
3. **El exceso se concentra en los conductores compuestos** (precio por pieza 2,58–2,85, plano; primos exactos a 1,000), y la pendiente total es el efecto mezcla del desplazamiento de B_M hacia los compuestos (share_C 0,42 → 0,69). El fracaso de C2 para B̂₂/B̂₃ es de tipo saturación (cociente acotado por ~2,9 en esta familia), no de divergencia como el de B̂₁.
4. **La línea de mayorantes absolutos en la base E queda AGOTADA como ruta a un lema demostrable**, por dos motivos independientes, cada uno suficiente:
   (a) *Numérico:* las dos únicas formas de la familia con estructura (B̂₂, triángulo total; B̂₃, el firmado más fino compatible con un mayorante absoluto sobre una partición) fallan C2 con pendientes 0,136 y 0,133 — la refinación diseñada recuperó el 2,4% de la pendiente. La libertad restante de la familia es la partición: particiones más gruesas pagan más triángulo entre bloques, y el límite de firmar más fino es B_M mismo, donde el «mayorante» deja de mayorar nada útil. No queda grado de libertad con recorrido medido.
   (b) *Cierre analítico (doc 103 §4.2 y §6.5, intacto tras E103.1):* cerrar B̂₂ o B̂₃ exige acotar **primeros momentos absolutos** ponderados Σ_d φ(d)|E_d|·u_d sobre el conjunto de divisores 𝒟. Ninguna herramienta clásica los controla: BV/BDH y la criba grande controlan segundos momentos (y en base de caracteres mueren aquí por δ_ℋ, doc 103); |E_d| no tiene cancelación entre residuos que explotar. Incluso un mayorante hipotético con pendiente 0 carecería hoy de ruta de demostración. (Esto responde también a la pregunta del encargo sobre los momentos necesarios: los de B̂₃ serían primeros momentos absolutos de E_d — no clásicos.)
   Precisión de honestidad: «agotada» significa que ningún miembro de la familia pasa C1+C2 **y** tiene ruta de cierre; no significa que B̂₃/B_M diverja (punto 3).
5. **Siguiente paso (redirección, no otra variante de la misma familia):** abandonar los mayorantes y trabajar sobre B_M en su propia forma diagonal (docs 96/97/102): B_M = Diag + términos cruzados con signo en la base E, conservando Δ_r, donde **todos** los objetos son segundos momentos de E — la clase que la tecnología tipo BV/BDH sí controla — y donde la diagonal ya tiene el tamaño correcto (Diag^C/B_M^C ≈ 1,5, docs 102/103). El problema abierto vuelve a ser el cruzado con signo, ahora con tres rutas cerradas con evidencia: factorización BV (muere por W_max, docs 101/102), criba en base de caracteres (muere por δ_ℋ, doc 103), mayorantes absolutos (mueren por C2 + cierre, docs 103/104). **Criterio previo para cualquier candidata futura, a fijar antes de intentar demostración alguna:** (i) expresable solo con segundos momentos de E_d sobre 𝒟 (sin primeros momentos absolutos, sin sumas de caracteres sobre ℋ); (ii) C1+C2 en los seis puntos. Lo que resucitaría la línea de mayorantes: una cota para primeros momentos absolutos ponderados de E_d sobre conjuntos de divisores, o una partición de {d | r} en clases equilibradas — esta última imposible en la estructura de conmutación actual, porque φ(h0·s)/φ(h0) = s−1 ≥ z−1 ≥ 10 fuerza el desequilibrio [T]. **Goldbach sigue abierta.**

---

## 6. Pendientes [P-rev]

1. Los heredados del doc 103 §6: textos completos de los docs 96/98/99 y scripts 94d/94e/97/99; demostración general de cierre por divisores de V⁺; la extrapolación de la holgura a N = 10²⁴; max_{d≥2} W_d sin medir por separado; δ_ℋ exacto en soportes grandes (sustituido por el sándwich); variantes MV con pesos de Selberg en el dual de caracteres (subsumidas por el cambio de base, no medidas).
2. **La lectura L2 de nueva(d) es una elección interpretativa declarada (§1.2).** Si el texto completo del doc 94 define la parte nueva de otro modo (p. ej. la literal L1, bajo la cual B̂₃ = B̂₂ [T]), el contenido de E103.1 se reduce al dictamen del doc 103 sobre B̂₂. La verificación contra el doc 94 íntegro queda pendiente de que el túnel vuelva.
3. La lectura de «saturación» del punto 4.3 (cociente acotado por el precio de la pieza C) es una extrapolación desde precios por pieza planos en la ventana medida (seis puntos, |H_N| ≤ 5 033), no un teorema: si en soportes mayores el precio por pieza C empezara a crecer, la pendiente total empeoraría por los dos canales. Medir B̂₃ en soportes S19/S23 cuando el coste lo permita queda como comprobación posible, no necesaria para el dictamen C1+C2.
4. El cierre analítico de la redirección de §5.5 (cruzado con signo en la base E vía segundos momentos) no tiene aún candidata formulada: fijar su forma y su criterio previo es el siguiente trabajo de la serie, no parte de este doc.

---

## Apéndice A. Chequeos del script (todos PASS)

| Chequeo | Contenido | Resultado |
|---|---|---|
| K0 | Maquinaria: Model97 R2 \|H\| = 880; build_H z = 11 == Model97 | PASS |
| K1 | B_M^(S) reproduce el doc 103 en los seis puntos (rel < 2·10⁻³) y B̂₂/B_M también (rel < 10⁻²): 1,709 / 1,920 / 2,231 / 1,826 / 2,053 / 2,318 | PASS |
| K2 | Por conductor: \|M_r\| ≤ Σ_h \|suma de clase\| ≤ Σ_d \|término\|; en agregado B_M ≤ B̂₃ ≤ B̂₂ en los ocho puntos (seis + R1 + R3) | PASS |
| K3 | Degeneración L1: (d, 3) = 1 para todo d \| r, r ∈ ℋ, en los ocho puntos (bajo L1, nueva(d) = d y B̂₃ = B̂₂) | PASS |
| K4 | Aditividad de pendientes: 0,135844 = 0,003209 + 0,132634 (tol 10⁻⁹) | PASS |
| K5 | C1+C2 para B̂₃: C1 PASS (1,679–2,262 ≤ 10), C2 FAIL (0,1326 > 0,1) | reporte (§2) |

Reproducibilidad: `python3 104_experimento_e103_1_bhat3.py` (solo stdlib, ~1,6 s), salida en `salida_doc104_local.txt`. El script reutiliza sin cambios la maquinaria de los docs 102/103 (importada del script 103: `Model97`, `build_H`, `SegPsi`, `M_r_E`, `ols`, `delta_r_orth`) y añade la partición por parte lisa en z, el mayorante B̂₃ con descomposición por piezas, los diagnósticos de clase (ρ_intra, clases/divisores por conductor, masa en pares) y la descomposición exacta de la pendiente.

## Apéndice B. Qué no es este doc

No es una demostración de la cuestión (a) ni de Goldbach: es la ejecución del experimento previo E103.1 que el doc 103 exigió antes de intentar cualquier cierre analítico de V103.2, y su resultado es negativo para el criterio de la serie (C1 sí, C2 no), con la hipótesis que lo motivaba refutada en su forma predictiva. No usa GRH. Los cocientes, pendientes y diagnósticos son mediciones sobre ventanas finitas, citados con su intervalo cuando procede; en cambio, la degeneración de la lectura L1 (B̂₃ = B̂₂), las desigualdades B_M ≤ B̂₃ ≤ B̂₂ y la estructura de la suma de clase μ(r/h0)φ(h0)[E_{h0} − (s−1)E_{h0·s}] son implicaciones rigurosas [T] dentro del modelo.
