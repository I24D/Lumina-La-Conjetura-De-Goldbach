# Doc 121 — Verificación de la respuesta de Lumina Chat al doc 120: el lema de unicidad de rama es teorema y (A6) queda recuperada; la reducción CRT de D-NP vale con el peso corregido Ω_d^C·φ(d)

**Fecha:** 2026-11-28. **Línea:** auditoría cruzada de (DIST-𝒟).
**Documentos de prueba (este doc es la síntesis):** `121a_unicidad_rama.md` (script PASS, 57,7 s, aritmética exacta), `121b_dnp_crt.md` (script PASS, K1–K4).
**Etiquetas:** [T] = teorema, [M] = medición/verificación numérica, [P] = pendiente.

---

## 0. Veredicto en breve

1. **El lema de unicidad de rama es [T] y (A6) queda recuperada con hipótesis explícitas** (§1, doc 121a). El Lema 2 del doc 101 §2.1 dice literalmente lo que Lumina Chat afirmó; sobre el soporte real, α̃_m ≤ Ω_m^C vale. Una precisión: su identidad de intercambio escrita usa el peso equivocado y es falsa como identidad; con el peso verdadero de (106.2) la desigualdad sobrevive intacta. Recuperada ≠ cerrada: la parte-F queda *identificada* con el momento liso ℳ_Ω^{liso}, sin ganancia cuantitativa frente a T(N).
2. **La reducción CRT de D-NP vale en su mecanismo y falla en su forma final escrita** (§2, doc 121b). Ortogonalidad y Cauchy–Schwarz/CRT confirmados; el peso correcto del objetivo es **Ω_d^C·φ(d)**, no Ω_d^C/φ(d) — la cuarta inversión de φ de la línea. El objetivo corregido coincide exactamente con el momento pesado por su W_d = Ω_d^C·φ(d) de la respuesta al doc 119: la reducción CRT y la ruta Baier–Zhao aterrizan en el mismo objeto.

Estado de (DIST-𝒟) tras este doc: dos frentes limpios — D-NP como momento de caracteres pesado por W_d sobre 𝒟_S (puente a Baier–Zhao; quedan sus hipótesis y el exponente de conteo), y la parte-F identificada con ℳ_Ω^{liso} (sin acotar). (H2) sin cambios: (FC-2), (BAND-δ-HI′), (BAND-LOW-NONDEG), (SLAB-V′). Goldbach permanece abierto.

---

## 1. Unicidad de rama: (A6) recuperada (doc 121a)

**El texto.** El Lema 2 del doc 101 §2.1 enuncia la factorización única del soporte tal como Lumina Chat la citó. Y S_A es 𝒮, el conjunto de sprimes de la serie: z ≤ q < y, q ∤ A. De ese lema se deduce [T]: para todo r ∈ ℋ_C y todo m, #{q ∈ 𝒮 : qm | r} ≤ 1 — ningún conductor del soporte real lleva dos primos de conmutación.

**La dominación.** Rehaciendo el intercambio de sumas con las definiciones verdaderas, α̃_m = Σ_q Ω_{qm}^C ≤ Ω_m^C vale sobre el soporte real [T con el lema; M: sin violaciones para todo m ∈ 𝒟 en los marcos S11–S23, y también en un marco hipotético A = 1, z = 5 construido para darle al contraejemplo su mejor oportunidad]. El 3·5·7 del doc 120b queda fuera por dos vías: necesitaría dos primos ≥ z en un mismo r, y en los marcos reales además 3 | A, con Ω_3^C = 0.

**La precisión debida.** La identidad de intercambio tal como ella la escribió usa el peso τ(r)⁻¹/φ(r) y como identidad es falsa (factor ~15 medido en S11, m = 5). El peso verdadero de (106.2) es (τ(r)−1)/φ(r); con él, la desigualdad — que es lo que la cadena usa — sobrevive intacta. La moraleja es la de siempre en esta línea: las desigualdades aguantan; las identidades hay que escribirlas con el peso exacto.

**Alcance exacto.** Con la dominación recuperada, la parte-F de D_MEAN^corr queda identificada con ℳ_Ω^{liso}: el momento Ω original restringido a los conductores lisos m. No es circularidad (es una subfamilia genuina, la rama no conmutada) y no es una cota: nadie ha ganado todavía un exponente frente a T(N) en esa masa. Ella misma lo advirtió; queda registrado como el enunciado correcto del frente F.

## 2. La reducción CRT de D-NP: mecanismo confirmado, peso corregido (doc 121b)

**Ortogonalidad [T][M].** S_{m,q}(χ) = (1/φ(m)) Σ_{ψ mod m} ψ̄(N)·Ψ_N(ψχ) es correcta **con el conjugado** ψ̄(N): el caso discriminante (N, m, q) = (200, 21, 13) refuta la variante sin conjugar. La condición q∤n la imponen los propios caracteres no principales mod q.

**Cauchy–Schwarz y conteo CRT [T][M].** El conteo es exacto: φ(m)(q−2) caracteres en ambos lados de la desigualdad, sin pérdidas ni duplicados por la biyección CRT (ψ, χ) ↔ η. El coste de Cauchy–Schwarz es real y queda medido: factor de pérdida λ entre 2,4 y 14,5 en los casos pequeños — la reducción no es gratis, y ese factor tendrá que figurar en cualquier balance final contra T(N).

**La hipótesis (m, N) = 1 [T].** Está garantizada por construcción: el filtro H_N del doc 101 (y doc 115 §1) restringe la cadena a conductores coprimos con N. Y es necesaria: en (200, 35, 11), con (m, N) > 1, la fórmula daría 0 mientras S = log 5·χ(3) ≠ 0 (chequeo K4). Ningún m de la cadena queda fuera de la reducción por esta vía.

**El peso final [T][M] — la corrección.** La identidad Θ̃_m(q) = (q−1)·Ω_{qm}^C es puntual (doc 115 §7), así que la conversión es término a término; y da Θ̃_m(q)·φ(m) = (q−1)φ(m)·Ω_{qm}^C = **φ(d)·Ω_d^C** con d = qm. La forma que ella escribió, Ω_d^C/φ(d), es falsa: en el marco genuino z = 17 el término verdadero supera al suyo por un factor 1622×–4768×, exactamente φ(d)²/λ. Es la cuarta inversión de φ de la línea (tras las del doc 115 y las de (A2)/(A6)): el patrón ya no es casualidad, es el punto débil sistemático de la colaboración, y por eso cada forma final pasa por auditoría antes de usarse.

**El objetivo corregido.** D_NP ≤ Σ_{d=qm ∈ 𝒟_S} Ω_d^C·φ(d) · Σ_{η mod d, η_q ≠ χ0} |Ψ_N(η)|². Es decir: el segundo momento de caracteres sobre la rama conmutada, pesado por **W_d = Ω_d^C·φ(d)** — exactamente el peso que ella misma introdujo en su ruta Baier–Zhao de la respuesta al doc 119. Las dos rutas son una: la reducción CRT dice *qué* momento hay que acotar; Baier–Zhao es la herramienta candidata para acotarlo por estructura del conjunto. El precio honesto: el objetivo es φ(d)² más pesado que el anunciado, y la pérdida λ de Cauchy–Schwarz va dentro.

## 3. Estado de la cadena tras el doc 121

- (DIST-𝒟), vía DIST-D corregida: ℳ ≤ 2·D_NP + 2·D_MEAN^corr; potencias de q en o(T(N)) [cerrado]; parte-F ≡ ℳ_Ω^{liso} [identificada, sin acotar]; D-NP reducido [T] al momento W_d sobre 𝒟_S [sin acotar; puente Baier–Zhao con hipótesis y exponente de conteo por resolver: 0,420232 exigido frente a 0,4211 medido].
- (H2): (FC-2), (BAND-δ-HI′), (BAND-LOW-NONDEG), (SLAB-V′) — sin cambios en este doc.
- Con el doc 115: (DIST-𝒟) + (H2) ⇒ UNIF-U-COND ⇒ C_k + C_{k+1} > 0 para c < c*. G0 sigue fuera. **Goldbach permanece abierto.**

## 4. Reproducibilidad

Pruebas en `121a_unicidad_rama.py` (PASS, aritmética exacta) y `121b_dnp_crt.py` (PASS, K1–K4), salidas `salida_doc121a_local.txt` y `salida_doc121b_local.txt`. Este doc es la síntesis; no introduce cómputo nuevo. Los scripts acompañan al doc en el envío a Lumina Chat, atendiendo a su nota de reproducibilidad del doc 120.
