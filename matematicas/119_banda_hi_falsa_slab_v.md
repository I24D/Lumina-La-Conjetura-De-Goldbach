# Doc 119 — La capa alta del modelo también es falsa: la banda vive o muere con la cancelación de E; (SLAB-V) se reduce a (SLAB-V′)

**Fecha:** 2026-11-28. **Línea:** UNIF-U / (H2), continuación de los docs 117–118.
**Documentos de prueba (este doc es la síntesis):** `119a_banda_hi.md` (script PASS, K0/K1, 0,7 s), `119b_slab_v.md` (script PASS, K0–K3, 0,5 s).
**Etiquetas:** [T] = teorema, [M] = medición en ventana, [P] = pendiente.

---

## 0. Veredicto en breve

El doc 118 mató el mayorante Σ Π/φ y dejó la capa alta de la banda como (BAND-MOD-HI): atacar (Ŝ−Π)⁺ directamente. Este doc lo intenta por todas las rutas y el resultado es más profundo:

1. **(BAND-MOD-HI) es FALSA asintóticamente** (§1, doc 119a) [T módulo el término principal FL del modelo; M en ventana]: el peso real (Ŝ−Π)⁺ **no** mejora al mayorante — a u fijo < 2, el defecto del modelo es el término principal W(u) − 1 > 0, no un error. La banda no puede cerrarse desde el modelo solo, en ninguna capa alta.
2. El residuo que sobrevive es **(BAND-δ-HI)** [P] (§2): la capa alta enunciada directamente sobre δ(r) = (Ŝ−Π) + E(r). Su contenido es la **cancelación** con E(r) ≤ 0 — la dominación puntual δ⁺ ≤ (Ŝ−Π)⁺ del doc 117b es cierta pero inútil ahí. La línea de la banda y la línea del flujo quedan **acopladas**: ambas dependen del tamaño de |E| contra (Ŝ−Π)⁺.
3. **(SLAB-V) se reduce [T] a (SLAB-V′)** (§3, doc 119b): la reindexación divisor–cofactor es exacta, el factor 1/ln z que faltaba es Π_full (la densidad de V⁺(L)), y el residuo es una familia finita de sumas sobre soportes desplazados V^(c) conservando la degeneración del cociente — sin cuantificador de visibilidad.

Inventario de (H2) tras este doc: **(FC-2), (BAND-δ-HI), (BAND-LOW-NONDEG), (SLAB-V′)**. Goldbach permanece abierto.

---

## 1. (BAND-MOD-HI) falsa: el defecto del modelo es término principal (doc 119a)

**El perfil [M].** Sobre u_r ∈ [1, 3/2), el cociente medido (Ŝ−Π)⁺/Π vale 0,5–1,3 (máximo 1,2947, por debajo de la cota A_min = 1,4273 del doc 113): el exceso del modelo es del **orden de Π mismo**, no una fracción pequeña.

**La razón [T, módulo FL].** A u fijo < 2, Ŝ(r)/Π(r) no tiende a 1: el modelo sin desplazar tiene por término principal la función W(u) > 1 (la condición m = 1 de V̂ excluye primos de peso Mertens log(2/u) que el producto Π conserva). El "defecto" Ŝ − Π ≈ (W(u) − 1)·Π es asintótico, no un resto. En consecuencia, sobre la familia ℱ del doc 118a restringida a {u_r ≥ u₀*}, la suma directa Σ (Ŝ−Π)⁺/φ hereda la divergencia log log z para todo u₀* < 29/20 — igual que el mayorante.

**La ruta L² muere igual [T].** Toda desigualdad del tipo (x)⁺ ≤ x²/Π + Π devuelve un término Σ Π/φ sobre la parte complementaria (la ruta muerta del doc 118), y el segundo momento mismo sobre ℱ es ≍ (W−1)²·Π/φ: no hay ganancia cuadrática porque el defecto no es ruido, es señal.

**Nota de ventana [M].** En z ≤ 53 la familia ℱ aún no entra en la banda (Π(pq) > 1/2), y la banda total medida es 0,10–0,32·Q_c: la divergencia log log (log log 53 ≈ 1,38) sigue siendo invisible numéricamente; el peso de ℱ se lee por el perfil, no por la suma.

## 2. (BAND-δ-HI): la cancelación con E es el contenido (doc 119a)

Con el modelo insuficiente en la capa alta, el enunciado que sobrevive es el original sobre el defecto completo:

**(BAND-δ-HI)** [P]: Σ_{r ∈ banda, u_r ≥ u₀*} δ(r)⁺/φ(r) ≤ C·Q_c, con δ(r) = (Ŝ(r) − Π(r)) + E(r).

Como (Ŝ−Π)⁺ ≍ Π sobre las familias obstructoras y E(r) ≤ 0 (si (DESP-SIGN)), el enunciado exige que |E(r)| cancele a (Ŝ−Π)⁺ salvo masa sumable C·Q_c: es una afirmación **cuantitativa** sobre el error de desplazamiento, no solo de signo. (BAND-δ) no queda refutada — la medición de δ en ventana siempre fue pequeña — pero pierde el apoyo de la reducción del doc 117b en la capa alta: a partir de aquí, E se mantiene **dentro** de δ en toda la argumentación de capa alta.

Esto acopla las dos líneas restantes de (H2): (DESP-SIGN)/(FLOW-COMP)/(FC-2) gobiernan el signo y la masa del lado Add de E; (BAND-δ-HI) pide su tamaño contra el exceso del modelo en la banda. Un mismo objeto, E, con dos demandas.

## 3. (SLAB-V) → (SLAB-V′): la reindexación divisor–cofactor (doc 119b)

**Reindexación exacta [T].** Todo testigo d ∈ V⁺(L) de un r degenerado (q₀³r ≥ L, d múltiplo de r, d < L) tiene cofactor c = d/r < q₀³, y vale la identidad 1/φ(r) = φ(c)/φ(d). La suma degenerada del estrato q₀ es la suma de pares P(q₀) = Σ_r N(r)/φ(r), con N(r) el número de testigos. Medición [M]: N(r) ≤ 3, concentrada en c ∈ {1, 5, 7} (marco A) o c = 1 (marco B); P(5)/suma verdadera = 3,00 / 2,00 / 1,00 según marco; P(7) = P(11) = 0 en ventana.

**El factor 1/ln z, localizado [T/M].** Es Π_full, la densidad de V⁺(L): cada suma de ventana U_c a razón fija q₀³ aporta solo el factor ln(q₀³) = O(1), y la predicción del modelo Π_full·A(c)·ln(q₀³) reproduce los U_c medidos dentro de un factor 2–3. La comparación con Q_c ~ Π_full²-escala queda así explicada sin misterio.

**La obstrucción exacta [T].** El mayorante cerrado W que suelta la degeneración del cociente no cierra: W(5) ≈ 16·Q_c, con pérdida 23–217× y soporte de cofactores creciente con z (10 → 17 → 19). Lo que hay que conservar es la degeneración del cociente misma.

**(SLAB-V′)** [P]: Σ_c U_c^deg ≤ C(q₀)·Q_c para q₀ ∈ {5, 7, 11}, donde U_c^deg son sumas de Mertens sobre los soportes desplazados V^(c) con la condición de degeneración del cociente incorporada — familia finita (c < 1331), sin cuantificador de visibilidad. Es la misma especie de objeto que los soportes desplazados de los docs 112–113: no se afirma implicación alguna entre (SLAB-V′) y (FL-DESP), pero viven en el mismo mundo.

## 4. Veredicto de cadena actualizado

- (DESP-SIGN) ⟸ (FLOW-COMP) + (REM-PERS-small) [114–116]; (REM-PERS-small) cerrado [T]; (FLOW-COMP) ⟸ **(FC-2)** [117a], bloqueado por la entrada de fuerza PNT y la cota inferior de HEAD37(r).
- (BAND-δ) ⟸ (DESP-SIGN) + (BAND-MOD) **solo en la capa baja** [117b, enmendado]: capa baja no degenerada = **(BAND-LOW-NONDEG)** [P]; capa baja degenerada = **(SLAB-V′)** [P] (vía (SLAB-V), reducido en §3). En la capa alta la reducción no aplica: el enunciado es **(BAND-δ-HI)** [P], directo sobre δ con E dentro.
- Muertos en esta serie (todos [T], con la salvedad de §1): (REG-Π), (REG-Π-LOW), (BAND-MOD-HI) — y con (BAND-MOD-HI), toda esperanza de cerrar la banda desde el modelo sin desplazar. **No refutados:** (BAND-δ), (DESP-SIGN), (H2).
- Con el doc 115: **(DIST-𝒟) + (H2) ⇒ UNIF-U-COND**.

**Inventario de residuos de (H2) tras el doc 119: (FC-2), (BAND-δ-HI), (BAND-LOW-NONDEG), (SLAB-V′).** El centro de gravedad se ha movido: tres de los cuatro hablan de E o de soportes desplazados — el desplazamiento, que en el doc 113 parecía un error inocente (E ≤ 0 en 4113 instancias), es ahora el objeto central de la línea.

Siguen fuera, sin cambios: (DIST-𝒟) en la línea de Lumina Chat, y la separación de capas G0. **Goldbach permanece abierto.**

## 5. Reproducibilidad

Pruebas y verificaciones en `119a_banda_hi.py` (PASS, K0/K1) y `119b_slab_v.py` (PASS, K0–K3), con salidas `salida_doc119a_local.txt` y `salida_doc119b_local.txt`. Este doc es la síntesis y el veredicto de cadena; no introduce cómputo nuevo.
