# DOC 102 — Experimentos E101.1 y E101.2: la holgura de Cauchy vive en el aplanamiento W_d → W_max y en los conductores compuestos; α(z) se estabiliza sobre 0,42

**Serie:** verificación y continuación de los docs 94–101 (Gram de Lumina, diagonalización, puente de Möbius, segundo factor, conteo A100.1).
**Fecha:** 2026-11-28. **Autora:** Lumina Muse.
**Archivos:** `102_experimentos_e101_1_e101_2.py` (solo stdlib; corre en local, **PASS** en 2,6 s), `salida_doc102_local.txt` (salida completa).
**Estado:** PARCIAL en el mismo sentido que los docs 95/97/100/101: los textos completos de los docs 98/99 y los scripts 94d/94e/97/99 de Lumina siguen sin estar disponibles en local (túnel caído); este doc trabaja sobre las definiciones del modelo finito tal como quedaron fijadas en los docs 94–97, 100 y 101, que sí están completos en local. Los puntos [P-rev] heredados se listan en §5.
**Disciplina:** se distingue en todo momento MEDICIÓN de PRUEBA. Nada de lo aquí medido demuestra Goldbach; la conjetura fuerte sigue abierta. Sin circularidad, sin GRH sin probar.

**Etiquetas:** [T] = teorema (demostrado en el texto o en un doc anterior citado), [N] = medición numérica reproducible con el script adjunto, [P-rev] = pendiente de revisión contra el texto completo.

---

## 0. Alcance

El doc 101 cerró la ruta BV para la cuestión (a) tal como la formuló el doc 100 (α_W = 0,4289 > 0,42) y reubicó el margen en la **holgura de Cauchy**: el segundo momento real B_M es 70–1145× menor que el producto W_max·S_BV, y sobre conductores primos la desigualdad es una identidad exacta. El doc 101 propuso dos experimentos, que este doc ejecuta:

- **E101.1 (holgura de Cauchy):** a N fijo, con soportes escalados, medir la holgura slack = W_max·S_BV/B_M, su crecimiento con |H|, su descomposición por etapas y su atribución a conductores primos frente a compuestos.
- **E101.2 (α en z mayor):** repetir la medición de α_W y α_H con z = 31 y z = 37 (el doc 101 llegó a z = 29) y decidir si α(z) tiende a cruzar 0,42 o se estabiliza por encima.

**Notación (docs 94–101).** Soporte H (cierre por divisores de los m = s·d, s sprime en [z, y), d ∈ V⁺(D/s)); H_N = {r ∈ H : (r, N) = 1}; E_d = error de ψ en progresiones (ψ real); M_r = Σ_{d|r} μ(r/d)φ(d)E_d; B_M = Σ_{r∈H_N, r>(ln N)^κ} |M_r|²/φ(r) con κ = 1; W_d = Σ_{r∈H_N, d|r} r/φ(r); W_max = max_d W_d = W_1 [T, doc 101]; S_BV = Σ_{d|r, r∈H_N} φ(d)|E_d|².

**La cadena que se descompone [T].** Para cada r, Cauchy con los pesos que dan el Lema 1 del doc 101 ((d/φ(d))·((r/d)/φ(r/d)) = r/φ(r)) da |M_r|² ≤ r·Σ_{d|r} φ(d)|E_d|². Sumando sobre el conjunto de B_M y luego sobre todo H_N:

B_M ≤ BV1* ≤ BV1_full ≤ W_max · S_BV,

donde BV1* = Σ_{r en el conjunto de B_M} (r/φ(r))·Σ_{d|r} φ(d)|E_d|² y BV1_full = Σ_d φ(d)|E_d|²·W_d. La holgura total se factoriza exactamente:

slack = S_Cauchy · S_cut · S_max, con S_Cauchy = BV1*/B_M (Cauchy por conductor), S_cut = BV1_full/BV1* (inflación por el corte κ), S_max = W_max·S_BV/BV1_full (aplanamiento W_d → W_max).

**Descomposición primo/compuesto.** B_M = B_M^P + B_M^C según r primo o compuesto (ídem BV1*, S_BV según d). Y la suma diagonal sin términos cruzados Diag_r = (1/φ(r))·Σ_{d|r} φ(d)²|E_d|², que cumple Diag^P = B_M^P exactamente [T: para r primo, M_r = φ(r)E_r tiene un solo término no nulo].

---

## 1. Resumen ejecutivo

| Pregunta | Respuesta | Tipo |
|---|---|---|
| ¿Cómo crece la holgura con \|H\|? | **Linealmente**: slack ~ \|H_N\|^θ con θ = 0,988 ± 0,147 (IC95%, familia genérica) y θ = 1,036 ± 0,020 (familia 2^e). Sin saturación en la ventana | [N] |
| ¿Qué etapa pierde? | **El aplanamiento W_d → W_max** (S_max = 25–2176, θ_max ≈ 0,82–0,91). El corte κ no pierde nada (S_cut = 1,000–1,018) y el Cauchy por conductor pierde poco (S_Cauchy = 1,0–3,9) | [N] |
| ¿La holgura se concentra en compuestos? | **Sí, cuantificado**: el 99,9–100% del exceso de la etapa Cauchy proviene de conductores compuestos en todos los puntos con compuestos; sobre soportes de solo primos el exceso es ~0 y B_M = S_BV exacto (dos familias independientes) | [N] + [T] |
| ¿α(z) cruza 0,42 al crecer z? | **No hay evidencia de cruce**: α_W(z=11..37) = 0,4281, IC95% [0,4243, 0,4318] (más estrecho que el del doc 101 y entero sobre 0,42); las pendientes locales oscilan 0,412–0,464 sin tendencia. α se estabiliza en ≈ 0,42–0,43 | [N] |
| ¿La ruta diagonal directa tiene margen real? | **Sí, y es grande**: la holgura medida (≈ \|H\|, sin saturar) supera en ~3 órdenes de magnitud el factor ~10⁵–10⁶ que le falta a la ruta BV en N = 10²⁴ (extrapolación de ajustes, §4) | [N] |

---

## 2. E101.1 — La holgura de Cauchy, medida y descompuesta

### 2.1 Diseño

Tres familias de N par (condición de Lumina, doc 99/101): **genérico** N = 999982 = 2·79·6329 (el N validado de los docs 100/101); **2^e** N = 1048576 = 2²⁰ (soporte intacto: todo r es impar, (r, N) = 1 siempre — el peor caso); **colapso** N = 720720 = 2⁴·3²·5·7·11·13 (muchos factores primos pequeños: el filtro (r, N) = 1 elimina toda parte d > 1 y deja solo conductores primos — control exacto). Soportes: R1 (54), R3 (68), R2 (880) de los docs 94/100 y los escalados S11 (802), S13 (1628), S17 (5033) de la familia D = z⁵, y = z^{10/3} del doc 101. ψ real vía `SegPsi` (validado en K0 contra el doc 101: W_max = 1006,34, S_BV = 4,139·10⁹, B_M = 3,777·10⁹ en R2/N = 999982, reproducidos; y contra `PsiEngine`, máx. dif. rel. 1,35·10⁻¹³ en la muestra).

### 2.2 Tabla de mediciones [N]

Familia genérica, N = 999982:

| Soporte | \|H_N\| | W_max | S_BV | B_M | slack | S_Cauchy | S_cut | S_max | B_M^C/B_M | exc. Cauchy de compuestos | Diag^C/B_M^C |
|---|---|---|---|---|---|---|---|---|---|---|---|
| R1 | 54 | 68,26 | 1,208·10⁸ | 8,244·10⁷ | 100,0 | 3,95 | 1,010 | 25,0 | 0,932 | 0,999 | 1,705 |
| R3 | 67 | 68,29 | 1,929·10⁸ | 1,898·10⁸ | 69,4 | 1,01 | 1,018 | 67,8 | 0,000 | 0,000 | — |
| S11 | 798 | 907,85 | 3,666·10⁹ | 3,366·10⁹ | 988,9 | 2,20 | 1,000 | 449,4 | 0,424 | 0,999 | 1,495 |
| R2 | 876 | 1006,34 | 4,139·10⁹ | 3,777·10⁹ | 1102,9 | 2,32 | 1,000 | 476,0 | 0,482 | 0,999 | 1,480 |
| S13 | 1621 | 1889,27 | 8,515·10⁹ | 7,644·10⁹ | 2104,6 | 2,63 | 1,000 | 801,6 | 0,555 | 1,000 | 1,488 |
| S17 | 5017 | 5995,20 | 3,274·10¹⁰ | 2,788·10¹⁰ | 7041,8 | 3,24 | 1,000 | 2175,8 | 0,689 | 1,000 | 1,533 |

Familia 2^e, N = 1048576 (soporte intacto):

| Soporte | \|H_N\| | W_max | slack | S_Cauchy | S_max | B_M^C/B_M | Diag^C/B_M^C |
|---|---|---|---|---|---|---|---|
| S11 | 802 | 912,79 | 1070,7 | 2,20 | 486,6 | 0,432 | 1,674 |
| R2 | 880 | 1011,28 | 1193,6 | 2,33 | 511,3 | 0,474 | 1,663 |
| S13 | 1628 | 1898,02 | 2234,2 | 2,68 | 832,2 | 0,541 | 1,616 |
| S17 | 5033 | 6015,88 | 7218,7 | 3,35 | 2153,0 | 0,693 | 1,576 |

(S_cut = 1,000 en los cuatro puntos; fracción del exceso Cauchy de compuestos: 0,999–1,000.)

Familia colapso, N = 720720 (control): en los cuatro soportes (R2, S11, S13, S17), **B_M^C = 0 y B_M = S_BV al último dígito** (chequeo K1), S_Cauchy = 1,00, y slack = W_max con cociente 1,000000 (W_max = 422–1504 según el soporte filtrado). La predicción del doc 101 §3.3 se reproduce exactamente en una familia de N distinta de los primoriales originales.

### 2.3 Lecturas

1. **La holgura crece como |H|¹·⁰.** Regresión log-log de slack contra |H_N|: familia genérica θ = 0,988 ± 0,053 (SE), IC95% ± 0,147, LOO [0,927, 1,071]; familia 2^e θ = 1,036 ± 0,005, IC95% ± 0,020, LOO [1,032, 1,039]. La holgura es, a primer orden, **el tamaño del soporte**: slack/|H_N| ∈ [1,0, 1,9] en toda la ventana, sin tendencia. Equivalentemente slack/W_max = S_BV/B_M ∈ [1,02, 1,20] en todos los puntos no triviales: el producto BV pierde W_max y un factor adicional de solo 2–20%.
2. **La pérdida dominante es el aplanamiento, no Cauchy.** S_max aporta θ ≈ 0,91 (genérica) / 0,82 (2^e) del crecimiento y factores de 25× a 2176×. S_Cauchy queda entre 1,0 y 3,9: en la familia genérica su crecimiento no es significativo (θ = 0,085 ± 0,344 de IC95%, LOO incluye valores negativos); en 2^e crece despacio (θ = 0,220 ± 0,072). S_cut ≈ 1,000: el corte κ = 1 es gratis. Traducción: la cota intermedia BV1_full = Σ_d φ(d)|E_d|²W_d ya solo sobreestima B_M por un factor S_Cauchy·S_cut ≈ 1,0–3,4; **el factor ~|H| entero se pierde en el último paso, al sustituir el perfil W_d por su máximo**. (Observación: el término d = 1 del perfil es el mayor, W_1, pero multiplica a φ(1)|E_1|² = 0 — el aplanamiento cobra el máximo exactamente donde el momento no tiene masa.)
3. **La atribución a compuestos es casi total.** La fracción del exceso de la etapa Cauchy (BV1* − B_M) atribuible a conductores compuestos es 0,999–1,000 en los seis puntos genéricos y los cuatro de 2^e con compuestos. El exceso primo es el término trivial Σ_{r primo} |E_r|² (de bound_r − término = (r − φ(r))|E_r|² = |E_r|²), despreciable frente al compuesto. Sobre soportes sin compuestos el exceso desaparece: R3 (z = 3 ⇒ d = 1 forzoso, soporte de solo primos) tiene B_M^C = 0 y S_Cauchy = 1,01 en N genérico, y la familia colapso da S_Cauchy = 1,00 exacto. **La predicción del doc 101 — «la holgura vive entera en los conductores compuestos» — queda confirmada y cuantificada en dos familias independientes** (primoriales del doc 101 y colapso 720720 aquí, más el caso R3).
4. **La cancelación de Möbius sobre compuestos es real pero modesta.** Diag^C/B_M^C ∈ [1,48, 1,71]: la suma diagonal (sin términos cruzados) sobreestima el momento compuesto en ~50–70%; los términos cruzados de M_r cancelan entre un tercio y dos quintos de la masa diagonal. No es la fuente de la holgura grande (que es el aplanamiento), pero es un factor ~1,5 que cualquier cota directa debe capturar o pagar.
5. **Los compuestos cargan cada vez más del momento verdadero.** B_M^C/B_M crece con el soporte: 0,42 (S11) → 0,69 (S17) en ambas familias. El problema de la cuestión (a) es, asintóticamente, un problema sobre conductores compuestos; los primos están cerrados por la identidad exacta.

---

## 3. E101.2 — α(z) con z = 31, 37: no cruza

Fórmula de conteo A100.1 [T, doc 101], con la variante ordenada de `vplus_from_table` validada contra la fórmula original en z = 29 (|H| exacto, W_1 con error relativo < 10⁻⁹, chequeo K4) y reproducción exacta de α_W/α_H del doc 101 en z = 11..29:

| z | D = z⁵ | y | #sprimes | \|H\| | W_1 = W_max |
|---|---|---|---|---|---|
| 29 | 20511149 | 74931 | 7381 | 47178 | 58126,02 |
| 31 | 28629151 | 93585 | 9028 | 64082 | 79203,30 |
| 37 | 69343957 | 168787 | 15386 | 132785 | 164856,28 |

Regresiones log-log (α = pendiente/10, N = z¹⁰) [N]:

| Cantidad | Ventana | α | IC95% | LOO |
|---|---|---|---|---|
| α_W | z = 11..29 (doc 101) | 0,4289 | [0,4223, 0,4354] | [0,4265, 0,4320] |
| α_W | **z = 11..37** | **0,4281** | **[0,4243, 0,4318]** | [0,4264, 0,4291] |
| α_W | z = 19..37 (mitad superior) | 0,4224 | [0,4144, 0,4304] | [0,4214, 0,4251] |
| α_H | z = 11..37 | 0,4211 | [0,4185, 0,4238] | [0,4206, 0,4216] |
| α_H | z = 19..37 | 0,4179 | [0,4097, 0,4261] | [0,4170, 0,4211] |

Pendientes locales de W_1 entre z consecutivos: 0,438; 0,430; 0,452; 0,414; 0,420; **0,464 (29→31); 0,414 (31→37)**. Las de |H|: 0,424; 0,421; 0,444; 0,409; 0,416; **0,459; 0,412**.

**Lectura.** Los dos puntos nuevos no muestran descenso: la local 29→31 *sube* a 0,464 y la 31→37 vuelve a 0,414 — oscilación alrededor de ≈ 0,43, sin tendencia monótona en siete locales. El ajuste extendido *estrecha* el intervalo del doc 101 y lo deja entero por encima de 0,42 ([0,4243, 0,4318]); el ajuste de la mitad superior baja el punto estimado a 0,4224 pero su LOO sigue sobre 0,42 y su IC95% abraza el umbral sin regalar el margen que el criterio del doc 100 exige (α < 0,42 sostenido *con margen* para absorber ln⁵N). Para |H| el panorama es el mismo a caballo de 0,42, y de todos modos la cota BV carga W_max, no |H| (doc 101 §4). **Conclusión de E101.2: α(z) se estabiliza en ≈ 0,42–0,43; no hay evidencia de que cruce 0,42 para z grande, y la última duda razonable que el doc 101 dejó abierta queda cerrada en contra de la ruta BV.** La vía «convertir la fórmula de conteo en cota analítica con α < 0,42» está muerta salvo un cambio de marco (otro y, otro D) que ya no sería esta familia.

---

## 4. Veredicto

**¿Tiene la ruta diagonal directa (cuestión (b)) margen real para cerrar la cuestión (a) sin la factorización BV?** **Sí.** La aritmética, con los ajustes medidos y sus caveats:

- La ruta BV (tasas medidas de ambos factores, doc 101 §3.4) sobreestima el objetivo de Chen por ≈ 1,2·10⁵ en N = 10²⁴ (≈ 6·10⁵ con la envolvente analítica de S_BV, doc 101 §2.5).
- La holgura medida crece como |H|^1,0 sin saturación (§2.3.1). Extrapolada con el ajuste |H| ≈ 0,033·N^0,4208 a N = 10²⁴ (|H| ≈ 4·10⁸), la holgura sería ≈ 5·10⁸ — **tres órdenes de magnitud más** que el factor que le falta a BV. Es decir: bajo las tasas medidas, B_M mismo quedaría ~10³–10⁴ veces *por debajo* del objetivo de Chen en 10²⁴. El margen no es estrecho: lo que BV no puede probar no es falso por poco, es holgado por mucho.
- La extrapolación se apoya en que la holgura es una propiedad del soporte, no de N: a soporte fijo, el doc 101 (parte A) midió S_BV/B_M ≈ 1,1–1,18 estable hasta N = 10⁸ y W_max independiente de N; aquí S_BV/B_M ∈ [1,02, 1,20] en dos familias. Es MEDICIÓN extrapolada, no prueba: los ajustes son sobre ventanas finitas y la frase anterior no es una cota.

**Dónde debe atacar la ruta directa (lo que E101.1 deja fijado):**

1. **No aplanar.** El factor ~|H| se pierde al sustituir W_d por W_max. Cualquier cota candidata que pase por «max_d W_d» hereda la pérdida entera. La forma a acotar es la cuadrática en el vector E con su perfil, o mejor, directamente en la base de conductores: B_M = Σ_{r∈H_A, r>ln N} (1/φ(r))·|S_r|² con S_r = Σ_{d|r} μ(r/d)φ(d)E_d + Δ_r(N) (identidad corregida, docs 96–100).
2. **El soporte es R_A = H_A con su conteo exacto** (doc 101: |H| ≈ 0,033·N^0,4208, W_max = W_1): una criba grande *restringida* a este conjunto de conductores — cuya densidad y estructura (s primo de switching × d ∈ V⁺) son conocidas — en vez de la criba grande sobre todos los conductores ≤ Q. El peor caso es el de soporte intacto (N = 2^e); el caso primorial está cerrado por identidad.
3. **Los compuestos son el problema; los primos no.** La cota debe capturar, sobre conductores compuestos, (i) la cancelación de los términos cruzados de Möbius (factor medido Diag^C/B_M^C ≈ 1,5) y (ii) el término Δ_r (≤ 2,59% por conductor, ~10⁻⁴ en el momento, doc 100: conservarlo, no despreciarlo en las desigualdades).

**Siguiente paso analítico concreto.** Formular y tratar de probar la cota de criba grande restringida

Σ_{r∈H_A, r>(ln N)} (1/φ(r)) |S_r|² ≪ N/(ln N)² · (objetivo de Chen del doc 100),

con S_r en la forma corregida con Δ_r, explotando que los conductores de H_A tienen un único primo ≥ z (el de switching) y parte suave en V⁺ — es decir, separar la suma en la parte prima (cerrada, identidad) y la compuesta, y en la compuesta usar la representación diagonal Q_s = Σ φ(r)C_s(r)² (docs 96/97) para el lado de coeficientes en vez de Cauchy sobre E.

**Criterio numérico de alcanzabilidad (antes de intentar la prueba).** Una desigualdad candidata B̂ para B_M está en la vía correcta si y solo si, evaluada sobre los datos medidos de este doc (N = 999982 y N = 2²⁰; soportes S11, S13, S17):

- **(C1)** B̂/B_M ≤ 10 en los seis puntos, y
- **(C2)** la pendiente log-log de B̂/B_M contra |H_N| es ≤ 0,1 (no crece con el soporte).

Referencia de contraste: el producto BV da B̂/B_M = slack ≈ 10³–7·10³ con pendiente ≈ 1,0 (falla C1 y C2); la cota intermedia BV1_full da ≈ 1,0–3,4 con pendiente ≈ 0,09–0,22 (pasa C1, roza C2) pero sigue pagando el perfil W_d, cuyo conteo tiene exponente α_W > 0,42 — por eso BV1 no basta analíticamente aunque numéricamente esté cerca. Una candidata que pase C1+C2 en estos seis puntos y cuya dependencia de N sea la del objetivo es la primera señal de que la cota es alcanzable; una que falle C2 está disfrazando la pérdida BV con otro nombre.

**Goldbach sigue abierta.** Lo que este doc cierra: la duda sobre α(z) (no cruza) y la localización exacta de la holgura (aplanamiento + compuestos). Lo que abre: la cota directa sobre B_M vía criba restringida a H_A, con criterio numérico previo C1+C2.

---

## 5. Pendientes [P-rev]

1. Los heredados del doc 101 §6: demostración general de que V⁺(L) es cerrado por divisores y monótono en el nivel; los 7 puntos contra los textos completos de los docs 96/98/99 y los scripts 94d/94e/97/99; confirmar contra el doc 94 completo que la familia D = z⁵, y = z^{10/3} es la parametrización asintótica pretendida para R_A.
2. La extrapolación de la holgura a N = 10²⁴ (§4) asume que slack es propiedad del soporte con S_BV/B_M ≈ const: verificado a soporte fijo hasta N = 10⁸ (doc 101) y en dos familias a N ≈ 10⁶ (este doc); no verificado con soporte y N creciendo a la vez más allá de la ventana.
3. El criterio C1+C2 es un filtro numérico de alcanzabilidad, no una condición suficiente: pasarlo no demuestra la cota; fallarlo sí descarta la candidata.
4. max_{d≥2} W_d (el perfil efectivo tras anular d = 1 vía E_1 = 0) no se midió por separado; es el dato que fijaría cuánto gana exactamente una cota BV con perfil frente a la de máximo.

---

## Apéndice A. Chequeos del script (todos PASS)

| Chequeo | Contenido | Resultado |
|---|---|---|
| K0 | Maquinaria: Model97 R2 \|H\| = 880; build_H == Model97 (R2 y escalado z = 11); SegPsi reproduce el doc 101 en R2/N = 999982 (W_max, S_BV, B_M); SegPsi vs PsiEngine (muestra, máx. dif. rel. 1,35·10⁻¹³); fórmula de conteo vs tabla del doc 101 (z = 11..29) | PASS |
| K1 | Diag^P = B_M^P en todo punto; colapso primorial en N = 720720: B_M^C = 0 y B_M = S_BV en los 4 soportes | PASS |
| K2 | Cadena B_M ≤ BV1* ≤ BV1_full ≤ W_max·S_BV en los 14 puntos (N, soporte) | PASS |
| K3 | argmax W_d = 1 (W_max = W_1 sobre H_N) en los 14 puntos | PASS |
| K4 | E101.2: α_W/α_H en z = 11..29 reproducen el doc 101 (0,4289/0,4208); variante ordenada de vplus == fórmula original en z = 29 | PASS |
| K5 | E101.1: tablas completas y regresiones θ (genérica θ_slack = 0,988; 2^e θ_slack = 1,036) | reporte (§2) |

Reproducibilidad: `python3 102_experimentos_e101_1_e101_2.py` (solo stdlib, ~2,6 s), salida en `salida_doc102_local.txt`. El script reutiliza sin cambios `Model97`, `PsiEngine`, `SegPsi` y la fórmula de conteo de los docs 97/100/101, y añade `build_H` (validado contra Model97), la medición de la cadena completa con descomposición primo/compuesto, y la extensión de la familia escalada a z = 31, 37.

## Apéndice B. Qué no es este doc

No es una demostración de la cuestión (a) ni de Goldbach: es la medición que localiza la holgura de Cauchy con atribución cuantitativa (≥ 99,9% del exceso de Cauchy en compuestos; el factor ~|H| en el aplanamiento) y que cierra la última duda sobre α(z). No usa GRH. Los exponentes θ y α son ajustes sobre ventanas finitas, citados siempre con su intervalo; la extrapolación a 10²⁴ de §4 es una estimación de margen, no una cota.
