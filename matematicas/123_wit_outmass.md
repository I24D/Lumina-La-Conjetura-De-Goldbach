# Doc 123 — (WIT) refutado con contraejemplos y reemplazado por la familia (POS); (OUT-MASS) reducida: la parte con signo es una identidad y el término dominante es una suma exacta sobre primos

**Fecha:** 2026-11-28. **Línea:** UNIF-U / (H2).
**Documentos de prueba (este doc es la síntesis):** `123a_wit.md` (script PASS, K0–K3, 10,4 s), `123b_outmass.md` (script PASS, K0–K4, 91,1 s).
**Etiquetas:** [T] = teorema, [M] = medición/verificación numérica, [P] = pendiente.

---

## 0. Veredicto en breve

1. **(WIT) está refutado como lema uniforme** (§1, doc 123a): la rigidez de cofactores era un artefacto de la ventana z ≤ 53. Pero el intento de demostración no fue en vano: produce la familia posicional **(POS)** [T] que caracteriza exactamente qué cofactores pueden ocurrir — c es siempre **primo** [T], acotado por una desigualdad posicional exacta — y dos contraejemplos verificados fijan dónde termina la rigidez (estratos 7 y 11).
2. **(OUT-MASS) queda reducida a residuos menores nombrados** (§2, doc 123b): la suma con signo sobre la capa es una **identidad exacta** [T] (Σ δ(r)/φ(r) = 0 sobre la enumeración completa), y el 99,9 % de la masa vive en conductores primos con u_r ∈ [1,35; 3/2), donde el término dominante tiene **fórmula exacta**: una suma Σ 1/(q−1) sobre una ventana explícita de primos. La constante medida no crece hasta z = 71. Los tres caminos de cierre terminan en sumas de Mertens sobre conjuntos delgados: ese es el obstáculo exacto, ahora con nombre propio por ruta.

Inventario de (H2) tras este doc: **(FC-2), (NET-LAYER) + (P-MASS)** [los dos residuos de (OUT-MASS)], **(BAND-LOW-NONDEG), (SLAB-V′)** [con la familia (POS) en lugar de (WIT)]. Goldbach permanece abierto.

---

## 1. (WIT) refutado; la familia (POS) lo reemplaza (doc 123a)

**Lo que el intento de prueba sí demuestra** [T]. De la posición impar de r en la lista fusionada se sigue la familia de desigualdades (POS): c·p_m² < q₀³·∏_{j>m} p_j para cada posición m. Dos consecuencias son teorema: **jamás ocurre un cofactor compuesto**, y en el estrato 5 con k impar, c = 1 forzosamente. La cota de última posición, c < q₀³/p_k², es **exacta, no una pérdida técnica**: deja pasar c = q₀⁺⁺ en el estrato 7 (13 < 13,72) y c ≤ 47 en el estrato 11.

**Los contraejemplos** [M, verificados predicado a predicado con las funciones de la serie]. Estrato 7 a z = 61 (L = 39179, r = 115 = 5·23): testigos {1, 7, 11, 13}, N(r) = 4 — ya hay un cofactor fuera de {1, q₀, q₀⁺} y la multiplicidad supera 3. Estrato 11 a z = 127 (L = 127⁴, r = 5·7·101·103): testigos {1, 11, 13, 17, 19, 23}, N(r) = 6, robusto en 17 niveles. (WIT) vale en z ≤ 53 [M] y para el estrato 5 queda [P] con un único caso residual nombrado.

**Consecuencia para (SLAB-V′).** El enunciado no cambia [P], pero mueren su reducción a tres sumas por estrato y la cota N(r) ≤ 3. El ataque de Mertens debe tratar la familia completa de cofactores — que ahora no es arbitraria: está caracterizada por (POS) (c primo, cota posicional). La familia es más grande que tres formas, pero sigue siendo una familia de primos con control posicional explícito.

## 2. (OUT-MASS): identidad con signo y término dominante exacto (doc 123b)

**Estructura** [T]. V̂⁽ʳ⁾ es la fibra de V⁺(L/r); la condición de ausencia de ambos mundos es Out(r) ⟺ r·λ(d′) ≥ L y Λ(r, d′) ≥ L. Las piezas lejanas de M(r) tienen forma cerrada por productos de Euler (verificada a 3,8·10⁻¹⁶) y cargan solo el 7–18 % de la masa: manda la pieza cercana.

**La identidad con signo** [T]. Sobre la enumeración completa, Σ δ(r)/φ(r) = 0 exactamente (vía Σ S/φ = Σ Π/φ = 1). Por tanto la suma con signo de δ sobre la capa es menos la del complemento, y la ruta con signo abarata el objetivo de (OUT-MASS) por un factor medido ×1,6–3,6 [M]. Lo que queda tras cobrar la identidad: **(NET-LAYER)** (la suma neta con signo sobre la capa, en su forma reducida) y **(P-MASS)** (la masa μ=+1 de los ausentes, que la cancelación no regala) [P].

**El hallazgo central** [T + M]. El 99,9 % de la masa pooled de M vive en conductores **primos** con u_r ∈ [1,35; 3/2), y ahí la pieza cercana tiene fórmula exacta:

M_near = Σ_{q ∈ 𝒬(p₀)} 1/(q−1), con 𝒬(p₀) = {q primo: q > p₀, q³ ≥ L, p₀·q < L},

reproducida sin discrepancias (p. ej. 0,4159 = la suma sobre los 7 primos de 𝒬 en el marco (37, 1286, 7)). La constante global C(z) oscila en [0,196; 0,315] para 37 ≤ z ≤ 71, sin tendencia de crecimiento, con máximos siempre en niveles u_L ≈ 2.

**Los residuos mínimos, por ruta** [P]. Ruta puntual: la cota M(r) ≤ (B_r − Π(r))/2 es [T] pero insuficiente, y deja **(LAYER-COUNT)**. Ruta exacta: deja **(PRIME-WINDOW)** — acotar la suma sobre capas de las sumas de ventana Σ_{q ∈ 𝒬(p₀)} 1/(q−1) — más tres restos menores (R1–R3). Los tres caminos terminan en sumas de Mertens sobre conjuntos delgados: ese es el obstáculo exacto de (OUT-MASS), ya sin niebla alrededor.

## 3. Estado de la cadena tras el doc 123

- (H2): (FC-2) [sin cambios]; (OUT-MASS) reducida a (NET-LAYER) + (P-MASS), con (PRIME-WINDOW) como su forma dominante exacta y (LAYER-COUNT) como residuo de la ruta puntual; (BAND-LOW-NONDEG) [sin cambios]; (SLAB-V′) con la familia de cofactores caracterizada por (POS) [c primo, cota posicional], sin la reducción de (WIT).
- (DIST-𝒟): en su frontera delimitada por el balance de factibilidad (doc de Lumina Chat, 2026-11-28): falta un segundo momento con ahorro X^{2c*} sobre los pesos exactos; la contabilidad de la multiplicidad quedó cerrada con la identidad A_{qr_0}.
- Con el doc 115 y la precisión del doc 108: (DIST-D + DIST-P) + (H2) ⇒ UNIF-U-COND; G0 fuera. **Goldbach permanece abierto.**

## 4. Reproducibilidad

Pruebas en `123a_wit.py` (PASS, K0–K3) y `123b_outmass.py` (PASS, K0–K4), salidas `salida_doc123a_local.txt` y `salida_doc123b_local.txt`. Este doc es la síntesis; no introduce cómputo nuevo.
