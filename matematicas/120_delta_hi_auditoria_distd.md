# Doc 120 — La cancelación puntual no existe: (BAND-δ-HI′); y la auditoría de la etapa DIST-D de Lumina Chat: la arquitectura sobrevive, tres normalizaciones no

**Fecha:** 2026-11-28. **Líneas:** UNIF-U / (H2) (§1) y auditoría cruzada de (DIST-𝒟) (§2).
**Documentos de prueba (este doc es la síntesis):** `120a_banda_delta_hi.md` (script PASS, K0–K3, 2,2 s), `120b_auditoria_distd.md` (script PASS, 90/90 chequeos, N ≤ 500, primos reales).
**Etiquetas:** [T] = teorema, [M] = medición/verificación numérica, [P] = pendiente.

---

## 0. Veredicto en breve

1. **(BAND-δ-HI): la cancelación puntual no ocurre** (§1, doc 120a). E(r) *afeita* el exceso del modelo (mediana 62,7 % donde no es cero), pero es exactamente 0 en el 58,2 % de las instancias con exceso positivo, y la cobertura completa es 0/67. El enunciado sobrevive reformulado como **(BAND-δ-HI′)** [P]: dominación *en suma* del exceso por el déficit del flujo Add sobre la capa — E es 100 % lado Add ahí (E_Rem ≡ 0 en las 43 instancias descompuestas). El margen ×28 de (FLOW-COMP) **no** transfiere: son demandas independientes sobre E.
2. **Auditoría DIST-D** (§2, doc 120b): la etapa de Lumina Chat **sobrevive corregida** — ℳ ≤ 2·D_NP + 2·D_MEAN^corr con la parte de potencias en o(T(N)) — pero tres piezas tal como las escribió son falsas: la normalización de E_m en (A2), el factor (q−1)² de más en su D_MEAN (A4), y el (q−1) invertido con dominación inexistente en (A6). Los frentes reales de su línea quedan: D_NP y la parte-F con α̃_m, **sin dominación gratuita**.

Inventario de (H2) tras este doc: **(FC-2), (BAND-δ-HI′), (BAND-LOW-NONDEG), (SLAB-V′)**. Goldbach permanece abierto.

---

## 1. (BAND-δ-HI′): E afeita, no cancela (doc 120a)

**La suma de capa en ventana es pequeña** [M]: Σ δ⁺/φ ≤ 0,093·Q_c en marcos A/B y estrés z = 53. Pero la familia obstructora ℱ sigue ausente de la banda en ventana, así que el enunciado asintótico no se decide por medición: ni refutado ni apoyado como cancelación.

**El mecanismo imaginado muere** [M]: cobertura completa de (Ŝ−Π)⁺ por |E|: **0/67**. E(r) = 0 exacto en el 58,2 % de las instancias no degeneradas con exceso positivo; donde E ≠ 0, cubre una mediana del 62,7 % (máx 99,3 %). La imagen correcta no es "E cancela el exceso" sino "E recorta el exceso donde el desplazamiento muerde, y no está en la mayoría de los conductores".

**Estructura** [M/T parcial]: en la capa, E es íntegramente lado Add: E_Rem ≡ 0 en las 43 instancias descompuestas. El residuo mínimo es por tanto una afirmación sobre el flujo Add en suma:

**(BAND-δ-HI′)** [P]: sobre la capa alta de la banda, el déficit del flujo Add domina en suma al exceso del modelo, módulo E_Rem = 0 y el signo de (DESP-SIGN) restringido a la capa.

**Independencia** [T de lectura]: el margen ×28 con que (FLOW-COMP) se sostuvo en su contraejemplo (doc 117a) es signo agregado en un solo conductor (r = 5); (BAND-δ-HI′) pide magnitud contra (Ŝ−Π)⁺ en conductores de banda. Ninguna de las dos demandas implica la otra: la línea del flujo y la de la banda comparten el objeto E pero no se rescatan entre sí.

## 2. Auditoría adversarial de la etapa DIST-D de Lumina Chat (doc 120b)

Objeto: las afirmaciones (A1)–(A7) de su respuesta al doc 118 (con la preparación de la del 117), auditadas contra las definiciones corregidas del doc 115 y por verificación numérica con primos reales (90/90 chequeos, N ≤ 500).

| Afirmación | Veredicto |
|---|---|
| (A1) Σ_a B_{m,q}(a) = A_m − L_{m,q} | **CONFIRMADA** [T][M] |
| (A2) E(N; qm, N) = Δ_{q,m} + E_m − L_{m,q}/(q−1) | **REFUTADA tal como está** |
| (A3) Δ por caracteres + Cauchy–Schwarz | **CONFIRMADA** [T][M] |
| (A4) ℳ ≤ 2·D_NP + 2·D_MEAN | **CONFIRMADA CON ENMIENDA material** |
| (A5) ℒ ≪ N^{4/5}(log N)^4 = o(T(N)) | **CONFIRMADA CON ENMIENDAS** |
| (A6) α_m = Σ_q Θ̃_m(q)(q−1) = Σ_q Ω_{qm}^C ≤ Ω_m^C | **REFUTADA** |
| (A7) Parseval de S_{m,q} | **CONFIRMADA** [T][M] |

**(A2) — el error de normalización, otra vez.** Su E_m resta ψ(N)·φ(m) en vez de ψ(N)/φ(m), que es la normalización de la serie (doc 101). Además a la identidad le falta un 1/(q−1) sobre el bloque (E_m − L_{m,q}) completo. La forma corregida **E(N; qm, N) = Δ_{q,m} + (E_m − L_{m,q})/(q−1) es teorema** y verifica numéricamente a 5,6·10⁻¹⁶. Y conste: sus dos Δ (respuesta al 117 vs respuesta al 118) son funciones distintas; la buena es la de clases B_{m,q} del 118.

**(A4) — la desigualdad vale; su D_MEAN, no.** La reducción ℳ ≤ 2·D_NP + 2·D_MEAN es correcta como desigualdad, pero el D_MEAN que ella escribió es **(q−1)² veces** el término medio verdadero (razón medida exacta en los chequeos): la derivación escrita no es la demostrada. Con D_MEAN^corr (sin el (q−1)² espurio) la reducción queda en pie y es la que debe usarse.

**(A5) — potencias de q: cerradas, con dos declaraciones debidas.** Σ_{m|h} φ(m)² ≤ ζ(2)·h², sin pérdida τ (mejor que lo afirmado). La cota ℒ ≪ N^{4/5}(log N)^4 exige declarar D_A ≍ N^{1/2} y absorber un log log; y su ℒ controla la parte-L de D_MEAN^corr, no la de su D_MEAN original. Con esas enmiendas: **las potencias de q quedan en o(T(N))** — este cierre de su línea es real.

**(A6) — (q−1) invertido y una dominación que no existe.** Contra el doc 115 §7, lo correcto es α̃_m = Σ_q Θ̃_m(q)/(q−1) = Σ_q Ω_{qm}^C. Y la desigualdad Σ_q Ω_{qm}^C ≤ Ω_m^C **no está en el doc 115 y es falsa** para la (106.2) sin un lema de unicidad de rama (contraejemplo numérico en el script). Consecuencia: la parte-F de D_MEAN^corr no tiene dominación gratuita por la masa original; o se prueba el lema de unicidad de rama, o la parte-F es un frente independiente.

**Balance.** La etapa DIST-D sobrevive corregida: ℳ_{Ω,S}^C ≤ 2·D_NP + 2·D_MEAN^corr, parte-L en o(T(N)), y los frentes reales son **D_NP** (dispersión de caracteres, qm | (n₁−n₂), sin valores absolutos término a término) y la **parte-F con α̃_m**. Es el mismo patrón del doc 115: su arquitectura es correcta y sus normalizaciones resbalan — tercera vez que una constante o un factor (q−1) decide si la cadena existe. Por eso esta auditoría va antes de su siguiente etapa, no después.

## 3. Veredicto de cadena actualizado

- (H2): (FC-2) y (BAND-δ-HI′) son demandas independientes sobre E (signo/masa del flujo; magnitud en suma sobre la capa). (BAND-LOW-NONDEG) y (SLAB-V′) sin cambios.
- (DIST-𝒟): la ruta DIST-D queda con dos frentes (D_NP; parte-F con α̃_m) y un cierre parcial verificado (potencias de q). La corrección (A2)/(A4)/(A6) es obligatoria antes de construir encima.
- Con el doc 115: **(DIST-𝒟) + (H2) ⇒ UNIF-U-COND**. Siguen fuera: G0. **Goldbach permanece abierto.**

## 4. Reproducibilidad

Pruebas en `120a_banda_delta_hi.py` (PASS, K0–K3) y `120b_auditoria_distd.py` (PASS, 90/90), salidas `salida_doc120a_local.txt` y `salida_doc120b_local.txt`. Este doc es la síntesis; no introduce cómputo nuevo.
