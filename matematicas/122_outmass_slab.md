# Doc 122 — La banda se reubica: (BAND-δ-HI′) se reduce a (OUT-MASS) por una identidad nueva; (SLAB-V′) abierta con la obstrucción exacta y un sub-lema de rigidez

**Fecha:** 2026-11-28. **Línea:** UNIF-U / (H2).
**Documentos de prueba (este doc es la síntesis):** `122a_banda_delta_hi_prima.md` (script PASS, K0–K4, 1,3 s), `122b_slab_v_prima.md` (script PASS, K0–K3, 5,4 s).
**Etiquetas:** [T] = teorema, [M] = medición/verificación numérica, [P] = pendiente.

---

## 0. Veredicto en breve

1. **La lectura ingenua de (BAND-δ-HI′) es falsa** (§1, doc 122a): el déficit Add no cubre al exceso ni siquiera en suma, y sobre la capa el apareamiento está vacío — no hay déficit que invoque. Pero la desigualdad sumada vale en ventana, y una **identidad nueva, δ(r) = E_Rem(r) + M(r) − P(r)** [T], reubica el problema: el residuo mínimo de la banda alta es **(OUT-MASS)** [P], la masa μ=−1 de los elementos ausentes de ambos mundos sobre la capa, que en ventana vale con C ≈ 0,32.
2. **(SLAB-V′) queda abierta con la obstrucción exacta identificada** (§2, doc 122b): no existe telescopio para q₀ fijo, ni factorización en el cofactor, ni promedio en L que la cierre; su cierre sería un teorema de Mertens sobre V^{(c)} con degeneración de cociente. La medición extendida (z ≤ 53) la apoya sin excepción y revela una rigidez nueva: solo ocurren cofactores c ∈ {1, q₀, q₀⁺} — sub-lema candidato **(WIT)** [P].

Inventario de (H2) tras este doc: **(FC-2), (OUT-MASS), (BAND-LOW-NONDEG), (SLAB-V′)** — con (WIT) como sub-lema candidato dentro de (SLAB-V′). Goldbach permanece abierto.

---

## 1. (BAND-δ-HI′) → (OUT-MASS): el déficit no era (doc 122a)

**La refutación de la lectura ingenua** [M]. La razón déficit/exceso D/X es < 1 en todo nivel con masa, y 0,000 en los niveles peores de z = 37 y del marco B, con tendencia de erosión en z y en u₀*. La razón es estructural: sobre la capa el apareamiento del doc 114 está **vacío** (Add⁺ = 0 elementos, Add⁻ = 41, R_A ≡ 0): lo que el doc 120a llamó "déficit Add" es masa cruda de Add, no un déficit de apareamiento. Y ningún lema puntual exceso ≤ f(déficit) existe: refutado en 67/67 instancias; 40 instancias con 𝒟_A = 0 cargan el 58,7 % de la masa U, todas con |Add| = |Rem| = 0 — ahí V⁽ʳ⁾ = V̂⁽ʳ⁾ y el exceso es un defecto intrínseco del modelo, no algo que un flujo pueda pagar.

**La identidad que reubica el problema** [T, verificada a 2,4·10⁻¹⁵]:

δ(r) = E_Rem(r) + M(r) − P(r),

donde M(r) y P(r) son las masas μ=−1 y μ=+1 de los elementos **ausentes de ambos mundos** (ni en V⁽ʳ⁾ ni en V̂⁽ʳ⁾). Con E_Rem = 0 en la capa [M, doc 120a], la parte positiva de δ sobre la capa la gobierna M: el residuo mínimo es

**(OUT-MASS)** [P]: Σ_{r ∈ capa} M(r)/φ(r) ≤ C·Q_c.

En ventana vale con C ≈ 0,32 (máx Σ M/φ = 0,3152·Q_c). El precio declarado: al quedarse solo con M se descarta la cancelación M↔P, que en ventana vale un 55,8 % (factor pooled ~2,26) — (OUT-MASS) es suficiente pero no necesaria, y su margen real es más fino que el 0,32 bruto si la cancelación persiste asintóticamente.

**Independencia fina** [T de lectura + M]: (FC-2) y la línea de la banda son independientes también en mecanismo de decaimiento — objetos, mecanismos y cuantificadores disjuntos; en ventana (FC-2) es vacuo (BAD ≡ 0) mientras la masa de la banda es positiva. Ninguna de las dos líneas rescata a la otra a ningún nivel.

## 2. (SLAB-V′): la obstrucción exacta y la rigidez (WIT) (doc 122b)

**Medición extendida** [M, z ≤ 53]: (SLAB-V′) se apoya sin excepción. P(5)/Q_c ≤ 0,783 con crecimiento en z desacelerado (0,598 → 0,711 → 0,783), incompatible con la pérdida (log z)² de los mayorantes refutados. El estrato 7 se enciende a z = 53 (P(7) = 1,76·Q_c; era 0 a z ≤ 37 — su ausencia era artefacto de ventana); el 11 sigue en 0. La masa vive solo en niveles pequeños de ℒ_A.

**La obstrucción, demostrada como imposibilidad de las vías estándar** [T]: (i) no hay telescopio para q₀ fijo — el telescopio del doc 118b1 vive en la dirección q₀ (1/P(q₀†) ≤ (log z)⁻³) y está agotado exactamente en el corte q₀ < 13; (ii) no hay factorización en el cofactor c — las condiciones de Nathanson sobre la lista fusionada lo impiden, y la medición lo confirma: jamás ocurre un c compuesto; (iii) no hay promedio en L disponible — el enunciado es uniforme en L; (iv) el factor 1/φ(P(q₀)) da solo las constantes 1, 1/4, 1/24: no es dominación. El cierre de (SLAB-V′), si existe, es un teorema de Mertens sobre los soportes desplazados V^{(c)} **con la degeneración de cociente conservada** — la misma especie de enunciado que (FL-DESP) en la línea del flujo.

**La rigidez nueva** [M]: los únicos cofactores que ocurren en toda la medición son c ∈ {1, q₀, q₀⁺}, con multiplicidad N(r) ≤ 3. Si esa rigidez es teorema, la familia "finita de sumas" de (SLAB-V′) se reduce a tres formas por estrato. Queda registrada como sub-lema candidato:

**(WIT)** [P]: en la subcapa degenerada visible, el cofactor c = d/r solo toma los valores 1, q₀ y q₀⁺ (el siguiente primo admisible tras q₀).

## 3. Estado de la cadena tras el doc 122

- (H2): (FC-2) [decaimiento de la masa mala del flujo]; (OUT-MASS) [masa μ=−1 de los ausentes de ambos mundos en la capa alta]; (BAND-LOW-NONDEG) [sin cambios]; (SLAB-V′) [abierta: Mertens sobre V^{(c)} con degeneración; sub-lema (WIT)]. (BAND-δ-HI′) queda absorbida por (OUT-MASS) vía la identidad δ = E_Rem + M − P.
- (DIST-𝒟): sin cambios en este doc (D-NP y la parte-F ≡ ℳ_Ω^{liso}; el balance de factibilidad de D-NP está en curso en la línea de Lumina Chat).
- Con el doc 115 y la precisión del doc 108: (DIST-D + DIST-P) + (H2) ⇒ UNIF-U-COND; G0 fuera. **Goldbach permanece abierto.**

## 4. Reproducibilidad

Pruebas en `122a_banda_delta_hi_prima.py` (PASS, K0–K4) y `122b_slab_v_prima.py` (PASS, K0–K3), salidas `salida_doc122a_local.txt` y `salida_doc122b_local.txt`. Este doc es la síntesis; no introduce cómputo nuevo.
