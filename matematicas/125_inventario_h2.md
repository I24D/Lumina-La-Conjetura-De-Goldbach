# Doc 125 — Síntesis: muere la ruta absoluta de la banda alta; la ruta con signo queda como única viva, con el signo de E_Rem rescatado como teorema

**Fecha:** 2026-11-28. **Línea:** UNIF-U / (H2), continuación de los
docs 122a/123b (banda alta), 123a ((POS)) y 124 (auditoría de
Lumina Chat).
**Documentos de prueba (este doc es la síntesis):** `125a_prime_window.md`
(script PASS, K0–K3), `125b_erem_pos.md` (script PASS, K0–K4).
**Etiquetas:** [T] = teorema, [M] = medición/verificación numérica,
[P] = pendiente.

---

## 0. Veredictos en breve

| Pieza | Estado tras el doc 123 | Estado tras el doc 125 |
|---|---|---|
| (PRIME-WINDOW) | [P], forma dominante exacta de (OUT-MASS) | **REFUTADA [T]** en el régimen u_L ≈ 2 (125a) |
| (OUT-MASS) absoluta, Σ M/φ ≤ C·Q_c | [P], reducida a (NET-LAYER) + (P-MASS) con ruta exacta viva | **FALSA [T]** en el régimen u_L ≈ 2 (125a) |
| (LAYER-COUNT) | [P], residuo de la ruta puntual | **MUERTA [T]** como condición suficiente de (OUT-MASS) absoluta |
| E_Rem = 0 uniforme en la capa | provisión (i), [M] en ventana | **FALSA**: contraejemplos exactos; = 0 [T] solo para L/r ≤ 343 (125b) |
| signo de E_Rem | no enunciado | **E_Rem ≤ 0 [T]** para L/r ≤ 1331; residuo (EREM-SIGN) [P] para L/r > 1331 |
| Estrato 5 de (SLAB-V′) | [P] con un caso residual (123a); cierre anunciado en el doc 124 | **CERRADO [T vía (POS)]**: c ∈ {1, 5, 7} (125b) |
| Concentración de la masa de (OUT-MASS) | «99,9 % primo ∩ [1,35; 3/2)» (123b) | **errata**: conjunta 92,6537 %; 99,8993 % es solo la marginal de primos |

Lo que NO cae: **(BAND-δ-HI′) como enunciado** y la **ruta con
signo** — (NET-LAYER) + (P-MASS) — que pasa a ser la única ruta
viva de la banda alta. **Goldbach permanece abierto.**

---

## 1. La muerte de la ruta absoluta (doc 125a)

**El lema de ventana [T].** Sobre la capa, los extremos efectivos
de 𝒬(p₀) no son los de su definición ingenua: como q permitido
exige q < z y en la capa vale L/p₀ ≥ z (por u_{p₀} ≥ 1), el tope
p₀q < L nunca muerde y el tope q > p₀ tampoco (q ≥ L^{1/3} > p₀).
Los extremos que mandan son A = L^{1/3} y B = z, y

$$M_{near}^{(1)}(p_0) \;=\; S(z) - S\!\left(L^{1/3}\right)
\;=\; \log\frac{3}{u_L} \;+\; \Delta E \qquad [T],$$

idéntico para todo p₀ del nivel (verificado a 5,6·10⁻¹⁷,
dispersión exactamente 0 [M]). La fórmula ingenua con
log log(L/p₀) es incorrecta en toda la capa: olvida el tope
q < z, que es el extremo que manda.

**La doble suma [T].** La suma externa sobre los p₀ de la capa
también es Mertens, y como M_near⁽¹⁾ no depende de p₀:

$$PW(L) \;=\; \log\frac{3}{u_L}\cdot
\log\frac{u_L/3}{\,u_L-3/2\,} \;+\; O\!\left(\tfrac{1}{\log z}\right),
\qquad u_L \in (3/2,\,3)\ \text{fijo}.$$

En el régimen de los marcos, u_L ≈ 2: PW → log(3/2)·log(4/3)
= **0,1167**, una constante positiva — es la ventana completa
de primos (L z^{−3/2}, L^{1/3}] la que carga el término
principal, no un sub-rango extremo aislable. Mientras tanto el
presupuesto decae: **Q_c(z) ~ c_Q/log z** con c_Q ≈ 1,72
(Q_c·log z = 1,687–1,718 medido hasta z = 1009 [M]). Luego

$$\frac{PW(L)}{Q_c(z)} \;\sim\; 0{,}068\cdot\log z
\;\longrightarrow\; \infty \qquad [T].$$

Basta Mertens clásico (de hecho basta E(x) = o(1)): el error
sumado sobre la capa es o(1) y no decide nada. La medición
explica por qué la ventana no lo veía: PW ya carga el 80–92 %
del C(z) medido a los niveles máximos [M], pero el crecimiento
~ log z es indetectable como tendencia en z ≤ 71 (log z solo
recorre 3,61 → 4,26); la planitud de C(z) ∈ [0,196; 0,315] es
el comienzo de la recta 0,068·log z, no una contradicción.

**Consecuencia [T].** Como SM ≥ PW término a término,
**(OUT-MASS) en forma absoluta es falsa** en el régimen
u_L ≈ 2. Muere la ruta absoluta del doc 122a (δ⁺ ≤ M con la
masa acotada en valor absoluto) y muere cualquier condición
suficiente para ella — en particular **(LAYER-COUNT)**: si la
suma puntual Σ (B_r − Π(r))/(2φ(r)) cupiera en C·Q_c, daría
(OUT-MASS) absoluta vía la cota puntual [T] del doc 123b; el
blanco es falso, luego la ruta suficiente también.

**Qué sobrevive, exactamente.** (BAND-δ-HI′) es un enunciado
sobre δ, no sobre M: su reducción del doc 122a descartaba la
cancelación M↔P (factor pooled 2,26) y por eso pasaba por la
ruta absoluta; el enunciado mismo no queda tocado. La ruta que
queda es la **con signo**: la identidad exacta
δ(r) = E_Rem(r) + M(r) − P(r) junto con la identidad global
Σ δ(r)/φ(r) = 0 [T, doc 123b] reducen la banda alta a
**(NET-LAYER)** (la suma neta con signo sobre la capa) y
**(P-MASS)** (la masa μ = +1 que la cancelación no regala).
Sin atajo absoluto: la cancelación M↔P es obligatoria.

## 2. E_Rem: teorema parcial, contraejemplos y rescate del signo (doc 125b)

**Teorema parcial [T].** E_Rem(r) = 0 para todo conductor de
capa con **L/r ≤ 343** (y para todo r degenerado, sin
hipótesis). El mecanismo es un lema de paridad nuevo: si
d′ ∈ Rem(r) y m es una posición impar de la lista fusionada
de r·d′ que falla, el número j de primos de d′ en el prefijo
de longitud m es **par** — j = 0 muere por herencia
(Lema 1: r visible ⟹ r ∈ V⁺(L)) y j impar muere contra la
condición de V̂ en la posición j — luego ω(d′) ≥ 2, y el
umbral m = 1 de V⁺(L/r) fuerza L/r > 7³ = 343. Esto cubre la
capa entera si z ≤ 49 y 66 de las 67 instancias de los docs
120a/122a (la 67ª, del estrés z = 53 con L/r > 343, queda
[M] con E_Rem = 0 verificado).

**La forma uniforme es falsa [M, contraejemplos exactos].**
La caza extendida (z ∈ {61, 67, 71}, 70 instancias de
residuo) realiza exactamente la configuración mínima que el
lema de paridad permite: **Rem(r) = {35}** para
r ∈ {11, 13} primos, con fallo (m, j) = (3, 2) y
**E_Rem(r) = −1/24**, en 12 instancias; ningún otro elemento
de Rem aparece. La provisión (i) del doc 120a/122a no puede
promoverse a [T] tal cual.

**El rescate: el signo es teorema [T].** Si L/r ≤ 1331, todo
elemento de Rem tiene ω(d′) = 2 forzado (ω(d′) ≥ 3 exigiría
q₁ ≥ 11, luego L/r > 11³), así μ(d′) = +1 y

$$E_{Rem}(r) \;=\; -\sum_{d'\in Rem(r)} \frac{\mu(d')}{\varphi(d')}
\;\le\; 0 \qquad [T],\quad L/r \le 1331,$$

en particular en todo marco con z ≤ 121. Es lo que la
reducción del doc 122a realmente usa: con la provisión
(i′) E_Rem ≤ 0 en lugar de (i), la identidad
δ = E_Rem + M − P sigue dando δ(r) ≤ M(r) − P(r) ≤ M(r) y
δ(r)⁺ ≤ M(r) puntualmente. La cadena no necesitaba la
anulación, solo el signo — y el signo es [T] en el rango.

**Residuo mínimo (EREM-SIGN) [P].** El signo de E_Rem para
**L/r > 1331** (posible solo con z > 121), donde elementos
de Rem con ω(d′) = 3 (μ = −1) ya son posibles por tamaño;
la configuración mínima de amenaza es d′ = 385 = 5·7·11,
que por sí solo daría E_Rem = +1/240. Nada en la ventana
actual (z ≤ 71) toca ese rango.

## 3. Estrato 5 cerrado y errata del doc 123 incorporada (doc 125b)

**Estrato 5 de (SLAB-V′): CERRADO [T vía (POS)].** La
corrección del doc 124 de Lumina Chat se re-deriva paso a
paso desde la familia (POS) del doc 123a. El punto que el
doc 123a perdió: la lista fusionada tiene primos distintos
en orden estrictamente decreciente, luego p_{k−1} ≥ el
siguiente primo tras p_k — no solo p_{k−1} ≥ p_k. Con eso:
k impar da c = 1 (y 5 | c es imposible); k par con 5 | c da
c < 625/49 < 12,76, luego c = 5; k par con 5 ∤ c y p_k = 7
da p_{k−1} ≥ 11 y c < 125·7/121 = **875/121 < 7,24**, luego
c ∈ {1, 7}; k par con 5 ∤ c y p_k ≥ 11 excluye el candidato
peligroso: c = 11 forzaría p_k = 11, p_{k−1} ≥ 13 y
c < 125·11/169 = **1375/169 < 8,14 < 11**, contradicción —
queda c = 1. Conclusión: **c ∈ {1, 5, 7}** en el estrato 5;
(WIT)-estrato-5 es teorema. La enumeración extendida lo
confirma [M]: conjunto de cofactores exactamente {1, 5, 7},
N máx = 3, sin violaciones. (SLAB-V′) como cota sigue [P];
los estratos 7 y 11 siguen refutados como (WIT) (doc 123a).

**Errata del doc 123 (corrección 2 del doc 124).** La
concentración de la masa pooled de (OUT-MASS) (SM = 4,0720)
se corrige con cifras: marginal de primos ω(r) = 1:
**99,8993 %**; marginal de la rebanada u_r ∈ [1,35; 3/2):
**92,7157 %**; **conjunta primo ∩ rebanada: 92,6537 %**. El
resto se descompone exacto: primos con u_r < 1,35 cargan el
**7,2456 %** — es lo que (R1) absorbe ahora: los ausentes
cercanos con ω(d′) ≥ 3 son un resto vivo sobre esa rebanada,
no un nulo de ventana — y los conductores compuestos el
**0,1007 %**, que es (R3) medido con cifra. (PRIME-WINDOW)
como enunciado queda intacto (cuantifica sobre todo p₀ primo
de la capa, marginal 99,8993 %); lo que baja es su cobertura
«con fórmula exacta»: el régimen donde M_near = M_near⁽¹⁾
íntegro es la conjunta, 92,6537 %. (R2) sin cambios.

## 4. Inventario de (H2) tras el doc 125

**Vivos [P]:**
- **(FC-2)** — sin cambios.
- **Banda alta = (NET-LAYER) + (P-MASS)** — la ruta con
  signo, ahora la única viva: la identidad
  δ = E_Rem + M − P con E_Rem ≤ 0 [T para L/r ≤ 1331] y la
  identidad global con signo del doc 123b. La cancelación
  M↔P ya no es una opción de economía: es el problema.
- **(EREM-SIGN)** — el signo de E_Rem solo para L/r > 1331
  (z > 121); amenaza mínima ω(d′) = 3 nombrada (d′ = 385).
- **(BAND-LOW-NONDEG)** — sin cambios.
- **(SLAB-V′)** — con la familia (POS) caracterizando los
  cofactores (c primo, cota posicional exacta) y el
  **estrato 5 cerrado [T]**: c ∈ {1, 5, 7}.

**Resueltos en este doc:**
- E_Rem = 0 [T] para L/r ≤ 343 y r degenerado; E_Rem ≤ 0
  [T] para L/r ≤ 1331 — la provisión (i) del doc 122a queda
  reemplazada por (i′) con rango explícito.
- (WIT)-estrato-5 [T vía (POS)].

**Muertos [T]:**
- **(PRIME-WINDOW)** como pieza de la ruta exacta hacia
  (OUT-MASS) absoluta: refutada asintóticamente en el
  régimen u_L ≈ 2 — término principal constante 0,1167
  contra presupuesto Q_c ~ c_Q/log z → 0.
- **(OUT-MASS) absoluta** (Σ M/φ ≤ C·Q_c): falsa en ese
  régimen, por SM ≥ PW término a término.
- **(LAYER-COUNT)** como condición suficiente de
  (OUT-MASS) absoluta: el blanco que probaría es falso.
- La forma uniforme «E_Rem = 0 en la capa»: refutada por
  Rem = {35}, r ∈ {11, 13}, E_Rem = −1/24 [M].

## 5. Qué queda en pie de la cadena completa

La cadena del doc 115 no se mueve: **(DIST-𝒟) + (H2) ⇒
UNIF-U-COND**, y con la precisión del doc 108,
(DIST-D + DIST-P) + (H2); G0 sigue fuera. (DIST-𝒟) sigue en
su frontera delimitada por el balance de factibilidad: falta
un segundo momento con ahorro X^{2c*} sobre los pesos
exactos. **(H2) sigue siendo el cuello**, y este doc cambia
su geografía interna: de los frentes de la banda alta, la
mitad absoluta está muerta por teorema — no por falta de
empuje — y la mitad con signo es todo lo que queda. Eso es
un avance de delimitación, no de cierre: (NET-LAYER) y
(P-MASS) deben estimarse explotando la cancelación M↔P que
las rutas muertas descartaban, con (EREM-SIGN) como
provisión de signo pendiente solo en el rango z > 121.
Los demás frentes de (H2) — (FC-2), (BAND-LOW-NONDEG) y
(SLAB-V′) con la familia (POS) y su estrato 5 cerrado —
quedan como estaban, pendientes.

**Goldbach permanece abierto.**

## 6. Reproducibilidad

Pruebas en `125a_prime_window.py` (PASS, K0–K3: reproduce
los cinco valores K4 del doc 123b; verifica el lema de
ventana a 5,6·10⁻¹⁷ con dispersión 0; mide PW/Q_c, PW/SM y
Q_c(z)·log z hasta z = 1009) y `125b_erem_pos.py` (PASS,
K0–K4: reproducción de las 67 instancias con |Rem|; caza
extendida z ∈ {53, 61, 67, 71}; falsación del lema de
paridad fuera de la capa; aritmética del estrato 5;
enumeración extendida del estrato 5), salidas
`salida_doc125a_local.txt` y `salida_doc125b_local.txt`.
Este doc es la síntesis; no introduce cómputo nuevo.
