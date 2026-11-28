# Doc 118 — El mayorante Π/φ muere a toda profundidad: (REG-Π) y (REG-Π-LOW) falsas; la banda se ataca directa sobre Ŝ

**Fecha:** 2026-11-28. **Línea:** UNIF-U / (H2), continuación de los docs 116–117.
**Documentos de prueba (este doc es la síntesis):** `118a_reg_pi.md` (script PASS, J0–J4, 7,0 s), `118b1_degenerada.md` (script PASS, J0–J4, 0,7 s), `118b2_no_degenerada.md` (script PASS, J0–J2b, 12,3 s).
**Etiquetas:** [T] = teorema, [M] = medición en ventana, [P] = pendiente.

---

## 0. Veredicto en breve

El doc 117 dejó la banda apoyada en **(REG-Π)**: Σ_{r∈banda} Π(r)/φ(r) ≤ C·Q_c. Este doc lo somete a prueba y el mayorante **muere**:

1. **(REG-Π) es FALSA uniformemente [T]** (§1, doc 118a): la familia ω = 2 {pq : x_p ∈ (11/20, 3/5), 5 ≤ q < p} entra entera en la banda a niveles u_L → 2 y su masa diverge como log log z. Solo la rebanada ω = 1 cierra [T], con constante log(10/3)/M_∞ ≈ 1,036.
2. **(REG-Π-LOW) es FALSA para todo corte u₀* ∈ (0, 3/2) [T]** (§2, doc 118b2): el mecanismo de 118a tiene hermanas a toda profundidad. Corrección al doc 117b §3.2: la loseta baja tiene medida de Mertens **positiva**; la medida cero es solo de la subcapa degenerada.
3. **La subcapa degenerada no cierra por conteo elemental** (§3, doc 118b1): los estratos q₀ ≥ 13 cierran [T] por telescopio de Euler, y el **Teorema A [T]** expulsa a los degenerados de ω ≤ 2 en el régimen L ≥ q₀⁹; el residuo es **(SLAB-V)** [P], la contabilidad de Mertens de la rebanada superior de V⁺(L) en los estratos finitos q₀ < 13.

**Ni (BAND-MOD) ni (BAND-δ) quedan refutados.** Lo que muere es la ruta del mayorante: a partir de aquí la banda se ataca **directamente sobre Ŝ − Π**, en dos enunciados: **(BAND-MOD-HI)** [P] (capa alta) y **(BAND-LOW-NONDEG)** [P] (capa baja no degenerada), más (SLAB-V) en la degenerada. (H2) sigue abierta; su inventario actualizado está en §4.

---

## 1. (REG-Π): falsa [T]; la rebanada ω = 1 cierra (doc 118a)

**Obstrucción [T].** A niveles u_L → 2, la familia ℱ = {pq : x_p ∈ (11/20, 3/5), 5 ≤ q < p} (x_p = log p/log z) entra entera en la banda, y su masa g = Σ_ℱ Π(r)/φ(r)·(normalización) cumple g ≥ log(12/11)·(log log z + O(1)) → ∞. El cociente a Q_c crece como (log(12/11)/M_∞)·log log z. La divergencia es solo log log z: por eso la ventana z ≤ 53 la veía plana en 0,17–0,46 — una divergencia log log es invisible en cualquier ventana alcanzable, y solo la construcción teórica la detecta. La parte ω = 2 medida a z = 53 (0,0436·Q_c) es exactamente esta familia (masa g 0,0505 en J1/J4).

**Lo que sí cierra [T].** La rebanada ω = 1, por Mertens directo: su aporte a Σ_banda Π/φ es ≤ (log(10/3)/M_∞)·Q_c ≈ 1,036·Q_c, verificado con sumas a Z = 10¹² contra B₁ con error ≤ 0,0104. (Es la misma constante B₁(5/3) del doc 116c.)

**Consecuencia.** El cierre de la capa alta de (BAND-MOD) que el doc 117b dio módulo (REG-Π) queda sin apoyo: la capa alta hereda el residuo **(BAND-MOD-HI)** [P] — (BAND-MOD) sobre {u_r ≥ u₀*} atacado directamente sobre Ŝ, sin pasar por Σ Π/φ. Las rutas naturales ahí son L² (acotar (Ŝ−Π)⁺ vía (Ŝ−Π)²/Π, que sí vive en el mundo de Q), no L¹.

## 2. (REG-Π-LOW): falsa a toda profundidad [T] (doc 118b2)

La familia ℱ de §1 concentra su masa divergente en u_r ∈ (7/5, 29/20) — capa alta — pero su **mecanismo** no: una hermana ω = 3 y una escalera de ω creciente producen la misma divergencia log log z a **toda** profundidad u_r < u₀*, para cualquier corte u₀* ∈ (0, 3/2). Luego la parte no degenerada de (BAND-LOW) no puede descansar en ningún (REG-Π-LOW): hereda el residuo directo **(BAND-LOW-NONDEG)** [P], en paralelo a (BAND-MOD-HI).

**Corrección al doc 117b §3.2 [T]:** la loseta baja de la banda tiene medida de Mertens positiva; la medida cero corresponde únicamente a la subcapa degenerada (§3).

**Medición [M].** La celda no degenerada con u_r < 1 está **vacía en ventana**, y los 3 miembros en banda de ℱ a z = 53 (ventana ampliada, acuerdo exacto con el J4 de 118a) son degenerados, con u_r ∈ [0,750, 0,906]. Es decir: donde el mayorante diverge, la masa real de Ŝ en ventana es degenerada o nula — otra confirmación de que el problema era el mayorante, no la masa.

## 3. La subcapa degenerada: estratos altos cerrados, Teorema A, y (SLAB-V) (doc 118b1)

La subcapa degenerada {Ŝ(r) = 1} ⟺ q₀³·r ≥ L (q₀ = menor primo de 𝒫 no divisor de r; doc 117b) no cierra por conteo elemental, con la obstrucción localizada con exactitud:

- **Estratos q₀ ≥ q₀† = 13: CERRADOS [T]** por telescopio de Euler, con constante de cola → 0.
- **Teorema A [T]:** si L ≥ q₀⁹, ningún degenerado visible tiene ω(r) ≤ 2. La hipótesis es exactamente aguda: las 188 instancias degeneradas de ventana tienen todas ω ≤ 2 y todas L < q₀⁹ — son artefactos del régimen pequeño; en el régimen asintótico los degenerados tienen ω ≥ 3.
- **La ruta Euler sin ventana pierde (log z)²** en el estrato q₀ = 5: da 3,25–3,72·Q_c y creciendo, contra una suma verdadera ≤ 0,355·Q_c que viene **toda** del estrato 5.
- **Residuo (SLAB-V) [P]:** la contabilidad de Mertens de la rebanada superior de V⁺(L) restringida a los estratos finitos q₀ < 13 — un problema finito en q₀ (5, 7, 11) pero asintótico en L, donde la estructura de V⁺(L), no el producto de Euler desnudo, debe dar la cota.

## 4. Veredicto de cadena actualizado

- (BAND-δ) ⟸ (DESP-SIGN) + (BAND-MOD) [117b, T]. (DESP-SIGN) ⟸ (FLOW-COMP) + (REM-PERS-small) [114–116]; (REM-PERS-small) cerrado [T]; (FLOW-COMP) ⟸ **(FC-2)** [117a, T], bloqueado por la entrada de fuerza PNT y la cota inferior de HEAD37(r).
- (BAND-MOD) se parte en: capa alta → **(BAND-MOD-HI)** [P] (directo sobre Ŝ; rutas L²); capa baja degenerada → **(SLAB-V)** [P] (estratos q₀ ≥ 13 cerrados [T]; Teorema A [T]); capa baja no degenerada → **(BAND-LOW-NONDEG)** [P] (directo sobre Ŝ).
- Muertos en este doc, todos [T]: (REG-Π), (REG-Π-LOW) para todo corte, y con ellos toda ruta de la banda vía el mayorante Σ Π/φ. **No refutados:** (BAND-MOD), (BAND-δ), (DESP-SIGN), (H2).
- Con el doc 115: **(DIST-𝒟) + (H2) ⇒ UNIF-U-COND**.

**Inventario de residuos de (H2) tras el doc 118: (FC-2), (BAND-MOD-HI), (SLAB-V), (BAND-LOW-NONDEG).** Cuatro enunciados, todos directos sobre los objetos reales (el flujo compensado, Ŝ, la rebanada de V⁺), ninguno apoyado ya en mayorantes que la serie haya refutado. El patrón se mantiene por sexta vez: el enunciado que muere era un sustituto; el verdadero es el directo, y la medición nunca vio el problema que el sustituto inventaba.

Siguen fuera, sin cambios: (DIST-𝒟) en la línea de Lumina Chat, y la separación de capas G0. **Goldbach permanece abierto.**

## 5. Reproducibilidad

Pruebas y verificaciones en `118a_reg_pi.py` (PASS, J0–J4), `118b1_degenerada.py` (PASS, J0–J4) y `118b2_no_degenerada.py` (PASS, J0–J2b), con salidas `salida_doc118a_local.txt`, `salida_doc118b1_local.txt` y `salida_doc118b2_local.txt`. Este doc es la síntesis y el veredicto de cadena; no introduce cómputo nuevo.
