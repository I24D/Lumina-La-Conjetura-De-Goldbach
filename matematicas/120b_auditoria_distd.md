# Doc 120b — Auditoría adversarial de la etapa DIST-D de Lumina Chat (resp. a los docs 117–118)

**Fecha:** 2026-11-28. **Línea:** UNIF-U / (DIST-𝒟), verificación cruzada según el protocolo conjunto.
**Objeto:** las afirmaciones (A1)–(A7) de la respuesta de Lumina Chat al doc 118 (descomposición Δ, reducción a D_NP + D_MEAN, cota de potencias ℒ, identidad α_m, Parseval).
**Etiquetas:** [T] = teorema/demostrado aquí, [M] = medición numérica, [P] = pendiente.
**Método:** contraste contra las definiciones corregidas del doc 115 (`115_pdom_g1.md`), álgebra término a término, y verificación numérica con primos reales (script `120b_auditoria_distd.py`).

---

## 0. Veredictos en breve

| Afirmación | Veredicto |
|---|---|
| Normalización de E_m / F_m | **REFUTADA** (la suya); vale la de la serie (§2) |
| (A1) Σ_a B_{m,q}(a) = A_m − L_{m,q} | **CONFIRMADA** [T][M] |
| (A2) E(N; qm, N) = Δ + E_m − L/(q−1) | **REFUTADA** tal como está; **CONFIRMADA CON ENMIENDA**: E = Δ + (E_m − L)/(q−1) (§4) |
| (A3) caracteres + Cauchy–Schwarz | **CONFIRMADA** [T][M] |
| (A4) ℳ ≤ 2 D_NP + 2 D_MEAN | **CONFIRMADA CON ENMIENDA material**: la desigualdad vale, pero su D_MEAN es (q−1)² veces el término medio verdadero (§6) |
| (A5) ℒ ≪ N^{4/5} (log N)^4 | **CONFIRMADA CON ENMIENDAS**: exige D_A ≍ N^{1/2} (no declarado) y un factor log log omitido (§7) |
| (A6) α_m = Σ_q Ω_{qm}^C ≤ Ω_m^C | **REFUTADA**: factor (q−1) invertido contra el doc 115; la desigualdad no está en el doc 115 y es falsa para (106.2) sin un lema extra (§8) |
| (A7) Parseval | **CONFIRMADA** [T][M] |

**Lectura global:** la etapa DIST-D sobrevive como estrategia tras las enmiendas — (A2) corregida + (A3) + (A7) dan una reducción válida, y (A5) sí elimina las potencias de q. Pero la reducción que ella escribió no es la que demostró: hay dos errores de normalización/factor del mismo tipo que el doc 115 ya corrigió (φ y (q−1) que multiplican donde deben dividir), y su análisis de la parte media (§5 de su respuesta) descansa sobre el α_m equivocado.

## 1. Definiciones de referencia (doc 115)

De la serie (docs 101, 105, 106, 115), sin ambigüedad:

- **E(N; d, N) = ψ(N; d, N) − ψ(N)/φ(d)** — el error de ψ en la clase N mód d contra la masa esperada de UNA clase reducida. Así está en el código de la serie (doc 101: `self.E[d] = acc[d] - psiN / phi_of(d, ...)`), y **F_d := φ(d)·E_d** (doc 105), luego F_d = φ(d)ψ(N; d, N) − ψ(N).
- **Ω_d^C = Σ_{r ∈ ℋ_C, d | r} (τ(r) − 1)/φ(r)** (doc 106, (106.2)).
- **Ω_{qm}^C = (1/(q−1))·Θ̃_m(q)**, con Θ̃_m(q) = Σ_{h ∈ V⁺(D_A/q)∖{1}, m | h} (2τ(h) − 1)/φ(h), y **ℳ_{Ω,S}^C = Σ_q (q−1) Σ_m φ(m)² Θ̃_m(q) E(N; qm, N)²** (doc 115 §7, [T], verificado por enumeración).
- En la suma ℳ: q es el primo de conmutación (q ≥ z), m es liso, **(m, q) = 1**, y q ∤ N (doc 115 §1: ningún primo de AN está en ℋ_P). Estas hipótesis se usan abajo y se declaran donde hacen falta.

## 2. La normalización de E_m: el punto de partida — REFUTADA la suya

Ella define (resp. 117 y 118) **E_m = A_m − ψ(N)·φ(m)** y F_m = φ(m)E_m. Contra §1:

- El término principal de A_m = ψ(N; m, N) es ψ(N)/φ(m), no ψ(N)φ(m): su E_m resta **φ(m)² veces** el término correcto. Su E_m no es el E(N; m, N) de la serie salvo si φ(m) = 1.
- Su F_m = φ(m)A_m − ψ(N)φ(m)² tampoco es el F_m de la serie (φ(m)A_m − ψ(N)).

Es el mismo patrón de error que el doc 115 corrigió dos veces en la identidad Θ (un φ que multiplica donde debe dividir), ahora dentro de la descomposición DIST-D. Todo lo que sigue depende de con qué E_m se lea (A2): la auditoría evalúa ambas.

## 3. (A1) — CONFIRMADA [T][M]

Los n ≤ N con n ≡ N (mód m) y Λ(n) ≠ 0 se reparten entre las clases mód q. Los de clase a = 0 son los divisibles por q, y un primo-potencia divisible por q es q^j: su masa dentro de la clase N (mód m) es exactamente L_{m,q} por definición. Luego Σ_{a ∈ (ℤ/qℤ)*} B_{m,q}(a) = A_m − L_{m,q}. ∎ Su razón es correcta. [M: exacto a <10⁻⁹ en los 8 casos del script.]

## 4. (A2) y las dos definiciones de Δ — REFUTADA tal como está; CONFIRMADA CON ENMIENDA

**Las dos Δ no son la misma.** Resp. 117: Δ⁽¹¹⁷⁾ = ψ(N; qm, N) − A_m − L/(q−1). Resp. 118: Δ⁽¹¹⁸⁾ = B_{m,q}(N) − (A_m − L)/(q−1), la desviación de la clase diagonal contra la media de las clases no nulas. Difieren en [−A_m(q−2) − 2L]/(q−1) ≠ 0 genéricamente [M: difieren en los 8 casos]. Solo Δ⁽¹¹⁸⁾ es una fluctuación centrada.

**(A2) con Δ⁽¹¹⁸⁾ y el E_m de la serie:** Δ + E_m − L/(q−1) = E(N; qm, N) + ((q−2)/(q−1))·E_m ≠ E(N; qm, N). [M: la discrepancia numérica coincide con ((q−2)/(q−1))E_m a 6 decimales en los 8 casos.] Con **su** E_m la discrepancia es de orden ψ(N)φ(m) [M: hasta −2924 en N=500]. **(A2) tal como está escrita es falsa con cualquiera de las dos normalizaciones.**

**La identidad correcta [T]:** por CRT ((m,q)=1), ψ(N; qm, N) = B_{m,q}(N mód q), y φ(qm) = (q−1)φ(m); despejando la media de clases de la definición de Δ⁽¹¹⁸⁾:

**E(N; qm, N) = Δ_{q,m} + (E_m − L_{m,q})/(q−1)**, con E_m el de la serie.

Es un teorema por partición de clases — no necesita verificación numérica, aunque la tiene [M: discrepancia ≤ 5,6·10⁻¹⁶ en los 8 casos]. La enmienda es un factor 1/(q−1) sobre TODO el paréntesis (E_m − L), no solo sobre L. Nota: la afirmación de la resp. 117 de que "la identidad" se verificó en 49 casos con error 2,7·10⁻¹⁵ no es reconstruible desde el texto transmitido (con las definiciones transmitidas, la identidad falla); qué se verificó allí exactamente queda [P].

## 5. (A3) — CONFIRMADA [T][M]

Por ortogonalidad de caracteres mód q sobre las clases no nulas: B_{m,q}(a) = (1/(q−1)) Σ_χ χ̄(a) S_{m,q}(χ), donde la restricción q∤n en S es inofensiva (χ(n) = 0 si q | n, para todo χ mód q) y da S_{m,q}(χ₀) = A_m − L_{m,q} (por (A1)). Restando el término χ₀ se obtiene su fórmula para Δ⁽¹¹⁸⁾. ∎ Cauchy–Schwarz: (q−1)|Δ|² ≤ ((q−2)/(q−1))·Σ_{χ≠χ₀}|S|² ≤ Σ_{χ≠χ₀}|S|² — su constante es válida (la fina es (q−2)/(q−1)). [M: identidad a ≤1,2·10⁻¹⁵ y CS holgada en los 8 casos.]

## 6. (A4): la reducción, término a término — CONFIRMADA CON ENMIENDA material

Desde ℳ = Σ_q (q−1) Σ_m φ(m)² Θ̃_m(q) E(N; qm, N)², con la (A2) **corregida**: E = a + b con a = Δ_{q,m} y b = (F_m − φ(m)L_{m,q})/((q−1)φ(m)), F_m el de la serie. (a+b)² ≤ 2a² + 2b²:

- **Parte NP:** 2(q−1)φ(m)²Θ̃·Δ² ≤ 2φ(m)²Θ̃·Σ_{χ≠χ₀}|S_{m,q}(χ)|² por (A3) — el (q−1) de ℳ lo consume exactamente el (q−1) de Cauchy–Schwarz. Su D_NP **cuadra término a término** ✓.
- **Parte media:** 2(q−1)φ(m)²Θ̃·b² = **2·Θ̃_m(q)·(F_m − φ(m)L_{m,q})²/(q−1)**. Su D_MEAN lleva (q−1) **multiplicando**: D_MEAN^{suya} = (q−1)² · D_MEAN^{corr} por término [M: razón medida = (q−1)² exacta].

Conclusión: la desigualdad ℳ ≤ 2D_NP + 2D_MEAN vale (su D_MEAN es mayor que el verdadero), pero **no por la derivación escrita**: no se sigue de su (A2) (falsa), y el peso de la parte media está inflado por (q−1)². La reducción correcta es

**ℳ_{Ω,S}^C ≤ 2·D_NP + 2·D_MEAN^{corr}**,  D_MEAN^{corr} = Σ_{q,m} Θ̃_m(q)·(F_m − φ(m)L_{m,q})²/(q−1).

La enmienda importa: acotar su D_MEAN tal cual sería demostrar más de lo necesario, y su discusión de la parte media (§5 de su respuesta) no controla ninguna de las dos versiones, porque usa el α_m de (A6).

## 7. (A5): la cota de potencias, desigualdad a desigualdad — CONFIRMADA CON ENMIENDAS

- L_{m,q} ≤ log N: **[T]** ✓ (Σ_j log q sobre q^j ≤ N es ≤ log N).
- Σ_{m|h} φ(m)² ≪ h²: **[T] ✓, sin pérdida τ** — contra la sospecha inicial de esta auditoría: φ(m) ≤ m y Σ_{m|h} m² = σ₂(h) ≤ ζ(2)·h², porque Π_{p^a‖h}(1 − p^{−2a−2})/(1 − p^{−2}) ≤ Π_p (1−p^{−2})⁻¹ = π²/6. Constante explícita π²/6.
- La suma en h: lo que la cadena necesita es Σ_{h ∈ V⁺(H)} (2τ(h)−1)/φ(h)·h² con H = D_A/q, que es ≍ Σ_{h≤H} τ(h)·h·(h/φ(h)) ≪ H² (log H)(log log H). Ella cita Σ_{h≤H} hτ(h) ≪ H² log H — correcto en orden, pero **omite el factor h/φ(h) ≪ log log h**: pérdida de un log log, absorbible.
- El cierre: Σ_{q≥z} q⁻³ ~ 1/(2z² log z) con z = N^{1/10} da el factor N^{−1/5}/log N; para llegar a N^{4/5} hace falta **D_A ≍ N^{1/2}**, hipótesis que su cadena no declara. Con D_A ≍ N^{1/2}: ℒ ≪ (log N)²·(log N)·(10/log N)·N^{4/5}·(log log)^{O(1)} — el exponente (log N)⁴ que ella afirma es holgado ✓.

Veredicto: la cota sobre su ℒ es correcta con esas dos enmiendas (hipótesis D_A ≍ N^{1/2} declarada; un log log absorbido). Nota de colocación: su ℒ es exactamente la parte-L de **D_MEAN^{corr}** (vía Σ_m Θ̃_m(q)g(m) = Σ_h (2τ(h)−1)/φ(h) Σ_{m|h} g(m), intercambio exacto [T], y (F − φL)² ≤ 2F² + 2φ²L²) — no la de su D_MEAN, cuya parte-L sería (q−1)² mayor. Es decir: (A5) sí prueba ℒ = o(T(N)) para la reducción **corregida**.

## 8. (A6): α_m contra el doc 115 — REFUTADA

Dos fallos independientes:

1. **Factor invertido [T].** El doc 115 §7 da Ω_{qm}^C = Θ̃_m(q)/(q−1). Luego Θ̃_m(q)·(q−1) = (q−1)²·Ω_{qm}^C ≠ Ω_{qm}^C: su α_m := Σ_q Θ̃_m(q)(q−1) **no** es Σ_q Ω_{qm}^C. La combinación correcta es **α̃_m := Σ_q Θ̃_m(q)/(q−1) = Σ_q Ω_{qm}^C**. Es el tercer error de este tipo en la línea (tras los dos del doc 115): un (q−1) que multiplica donde debe dividir.
2. **La desigualdad no está en el doc 115 y es falsa sin un lema extra [T][M].** El doc 115 no contiene "Σ_q Ω_{qm}^C ≤ Ω_m^C". Y con la definición de conteo (106.2), Σ_q Ω_{qm}^C = Σ_{r ∈ ℋ_C, m|r} (τ(r)−1)/φ(r)·#{q : qm | r}: cada r se cuenta una vez por primo de conmutación distinto que divida a r/m. Contraejemplo [M, script]: sobre R = {r₀}, r₀ = 3·5·7, m = 3: Σ_{q∈{5,7}} Ω_{qm} = 2·(τ(r₀)−1)/φ(r₀) = 0,2917 > Ω_m = 0,1458. La desigualdad solo valdría con un **lema de unicidad de rama** (cada r ∈ ℋ_C tiene a lo sumo un primo de conmutación sobre cada m) que no está demostrado en la serie: queda **[P]**.

Consecuencia: su frase "esta parte queda dominada por una subfamilia de la propia masa original" no está demostrada ni siquiera como dominación; y en la reducción corregida la parte-F de D_MEAN^{corr} es Σ_m α̃_m F_m² (con L dependiente de q, salvo esa salvedad), cuyo control es exactamente el problema abierto, no una dominación gratuita.

## 9. (A7): Parseval — CONFIRMADA [T][M]

Parseval sobre el grupo (ℤ/qℤ)*: Σ_χ |S_{m,q}(χ)|² = (q−1) Σ_{a≠0} B_{m,q}(a)², y el término χ₀ es S(χ₀) = A_m − L_{m,q} (por (A1), consistente con la restricción q∤n de S). Restándolo: su identidad exacta. ∎ [M: igualdad a <10⁻⁷ en los 8 casos, incluidos S(χ₀) y el lado derecho por separado.] La advertencia que ella misma hace — no tomar valores absolutos término a término tras expandir, porque se destruye la sustracción del término principal — es correcta y queda intacta: es la razón de que D_NP sea un problema de dispersión genuino.

## 10. Verificación numérica [M]

Script `120b_auditoria_distd.py` (salida: `salida_doc120b_local.txt`, **PASS, 90/90 chequeos**). Casos: (N, m, q) ∈ {(100,3,7), (100,4,7), (200,5,7), (200,6,11), (500,3,7), (500,9,13), (500,4,11), (300,7,13)}, primos reales vía Λ de criba, caracteres mód q por logaritmo discreto. Cubre: (A1) exacta; (A2) falsa con discrepancia igual a la fórmula predicha ((q−2)/(q−1))E_m (con el E_m de la serie) y discrepancias de orden ψ(N)φ(m) con el suyo; (A2) corregida exacta a ≤5,6·10⁻¹⁶; Δ⁽¹¹⁷⁾ ≠ Δ⁽¹¹⁸⁾; (A3) y su CS; (A7) Parseval completo; (A4) término a término con peso simbólico (la reducción corregida domina y la razón D_MEAN^{suya}/D_MEAN^{corr} = (q−1)² exacta); (A6) contraejemplo de conteo y el factor invertido.

## 11. Consecuencias para la cadena

- La etapa DIST-D, **corregida**, queda: ℳ_{Ω,S}^C ≤ 2·D_NP + 2·D_MEAN^{corr}, con la parte-L de D_MEAN^{corr} en o(T(N)) por (A5). Los dos problemas reales que quedan son los que ella identifica en su §6 — D_NP (dispersión de caracteres, con la sustracción de Parseval intacta) y la parte-F de D_MEAN^{corr} — pero el segundo debe reformularse con α̃_m y sin la dominación de (A6): su control independiente sigue abierto [P].
- Antes de construir encima: usar E_m y F_m de la serie (§1), Δ⁽¹¹⁸⁾, y D_MEAN^{corr}. Con las definiciones tal como ella las transmitió, la reducción no es válida.
- (DIST-𝒟) sigue abierta; esta auditoría no la demuestra ni la refuta. **Goldbach permanece abierto.**

## Reproducibilidad

`120b_auditoria_distd.py` + `salida_doc120b_local.txt` (PASS). Definiciones contrastadas contra `115_revision_pdom.md` §7, `105_forma_diagonal_segundos_momentos.md` §notación y (106.2) de `106_precio_tau_conteo_omega.md`, y el código de `101_experimento_e100_1_conteo_a100_1.py` para E(N; d, N).
