# Doc 125b — E_Rem = 0: demostración parcial (sub-rango L/r ≤ 343) y el residuo exacto; cierre del estrato 5 vía (POS); errata del doc 123

**Fecha:** 2026-11-28. **Línea:** UNIF-U / (H2), continuación de los
docs 120a/122a (E_Rem), 123a ((POS)) y 123b (concentración).
**Script:** `125b_erem_pos.py` (stdlib; importa la maquinaria de
los docs 109/117b/120a/123a sin tocarla). **Salida:**
`salida_doc125b_local.txt` — chequeos K0–K4.
**Etiquetas:** [T] = teorema, [M] = medición/verificación,
[P] = pendiente.

---

## 0. Veredicto en breve

1. **E_Rem = 0 queda DEMOSTRADO [T] en un sub-rango explícito de
   la capa**: todo conductor de capa (degenerado o no) con
   L/r ≤ 343 tiene Rem(r) = ∅, luego E_Rem(r) = 0. Esto cubre la
   capa entera de todo marco con z ≤ 49 — de las 67 instancias
   de los docs 120a/122a, 66 pasan de [M] a [T] — y el sub-rango
   L/r ∈ [z, 343] para z ≥ 53. El mecanismo es un lema de
   paridad (§1): un elemento de Rem exige que la posición de
   Nathanson que falla en la lista fusionada tenga un número
   **par** ≥ 2 de primos de d′ en su prefijo, lo que fuerza
   ω(d′) ≥ 2 y L/r > 343.
2. **Pero la provisión (i) uniforme es FALSA, y el contraejemplo
   es exactamente la configuración mínima predicha**: en la
   rebanada L/r > 343 la caza extendida (K1) realiza
   Rem(r) = {35} para r ∈ {11, 13} primos, con fallo
   (m, j) = (3, 2) y E_Rem(r) = −1/24, en z = 61, 67, 71
   (12 instancias; ningún otro elemento de Rem en 70
   instancias de residuo probadas). Lo que sobrevive, y es lo
   que la reducción del doc 122a realmente usa, es el signo:
   **E_Rem(r) ≤ 0 [T] en toda la capa con L/r ≤ 1331** (todo
   marco con z ≤ 121), porque ahí todo elemento de Rem tiene
   ω(d′) = 2 forzado, luego μ(d′) = +1 y −Σ_{Rem} a ≤ 0. El
   residuo mínimo que queda es **(EREM-SIGN)** [P]: el signo
   de E_Rem para L/r > 1331 (solo posible con z > 121), donde
   elementos de Rem con ω(d′) = 3 (μ = −1) ya son posibles
   por tamaño.
3. **Estrato 5 de (SLAB-V′): CERRADO [T]** — la corrección del
   doc 124 de Lumina Chat se re-deriva desde (POS) paso a paso
   (§3): c ∈ {1, 5, 7}, y la enumeración extendida lo confirma
   [M, K4]: conjunto de cofactores exactamente {1, 5, 7},
   N máx = 3, sin violaciones en A z = 61, A z = 67 y
   B z = 61. El caso residual que el doc 123a dejó abierto
   muere al usar p_{k−1} ≥ siguiente primo tras p_k, no
   p_{k−1} ≥ p_k.
4. **Errata formal del doc 123** (§4): la concentración
   conjunta primo∩banda [1,35; 3/2) de la masa de (OUT-MASS)
   es **92,6537 %**, no el 99,9 % (que es solo la marginal de
   primos). Consecuencia exacta para (PRIME-WINDOW) y R1–R3.

Goldbach permanece abierto.

## 1. E_Rem = 0: la razón estructural

**Definiciones (serie).** Por conductor r: V⁽ʳ⁾ = {d′ :
(d′, r) = 1, r·d′ ∈ V⁺(L)} (fibra fusionada), V̂⁽ʳ⁾ = {d′ :
(d′, r) = 1, d′ ∈ V⁺(L/r)} (fibra modelo; la equivalencia con
`in_Vhat` es [T], doc 123b §1). Add(r) = V⁽ʳ⁾ ∖ V̂⁽ʳ⁾,
Rem(r) = V̂⁽ʳ⁾ ∖ V⁽ʳ⁾, a(d′) = μ(d′)/φ(d′), y
E_Rem(r) := −Σ_{d′∈Rem(r)} a(d′). En forma entera:
d′ ∈ V̂⁽ʳ⁾ ⟺ r·d′ < L y r·q_i²·Q_i < L para todo i impar
≤ ω(d′) (lista decreciente q_1 > … de d′, Q_i = q_1···q_i);
d′ ∈ Rem(r) ⟺ además r·d′ ∉ V⁺(L), es decir, alguna posición
impar m de la lista fusionada de r·d′ falla:
p_m²·P_m ≥ L (P_m = producto de los m mayores primos de r·d′).

**Lema 1 (herencia) [T].** r visible ⟹ r ∈ V⁺(L). En efecto,
r | d con d ∈ V⁺(L) y ambos libres de cuadrados: los primos
de r son un subconjunto de los de d, luego el i-ésimo primo
de r es ≤ el i-ésimo de d y los productos de prefijos de r
son ≤ los de d, término a término; las condiciones de
Nathanson de d en posiciones impares implican las de r, y
r ≤ d < L. ∎

**Lema 2 (paridad del fallo) [T].** Sean r visible y
d′ ∈ Rem(r). Sea m una posición impar donde la lista
fusionada de r·d′ falla, y sea j el número de primos de d′
entre los m mayores de r·d′ (el prefijo se parte en
P_m = R*·Q_j, con R* el producto de los primos de r del
prefijo y Q_j = q_1···q_j). Entonces **j es par**.
*Prueba.* (i) j = 0: el prefijo es exactamente el prefijo de
longitud m de r, y el fallo dice s_m(r)²·(prefijo) ≥ L,
contra el Lema 1 (m ≤ ω(r), m impar). (ii) j impar ≥ 1: la
condición de V̂ en la posición impar j da
r·q_j²·Q_j < L. Pero p_m ≤ q_j (q_j está en el prefijo y
p_m es su mínimo) y R* | r (producto de primos distintos
de r), luego R* ≤ r, y el fallo daría
L ≤ p_m²·R*·Q_j ≤ q_j²·r·Q_j < L, contradicción. ∎

**Corolario [T].** Todo d′ ∈ Rem(r) cumple ω(d′) ≥ 2
(j par y j ≠ 0 fuerzan j ≥ 2) y la lista fusionada tiene
longitud ≥ m ≥ j + 1 ≥ 3 (m impar, m > j). Además el caso
degenerado es trivial: si V̂⁽ʳ⁾ = {1}, Rem(r) = ∅ porque
1 ∈ V⁽ʳ⁾ por el Lema 1 [T, sin hipótesis sobre L/r].

**Lema 3 (umbral) [T].** Si d′ ∈ V⁺(L/r) tiene ω(d′) ≥ 2,
entonces L/r > 343. En efecto, q_1 necesita un primo
permitido menor q_2 ≥ 5, luego q_1 ≥ 7, y la condición
m = 1 de V⁺(L/r) da q_1³ < L/r, es decir L/r > 7³ = 343. ∎

**Teorema (E_Rem = 0 en el sub-rango bajo de la capa) [T].**
Sea r un conductor de la capa alta (banda, u_r ∈ [1, 3/2),
degenerado o no) con **L/r ≤ 343**. Entonces Rem(r) = ∅ y
E_Rem(r) = 0. *Prueba.* Si r es degenerado, el corolario.
Si no: todo d′ ∈ V⁺(L/r) tiene ω(d′) ≤ 1 por el Lema 3
(contrapositivo), luego ningún d′ puede estar en Rem(r)
por el corolario del Lema 2. ∎

**Alcance exacto.** La capa tiene L/r < z^{3/2}; si
z ≤ 49 entonces z^{3/2} ≤ 343 y el teorema cubre la capa
entera. Las 67 instancias no degeneradas de los docs
120a/122a viven en marcos con z ≤ 37 (A), z = 29 (B) y el
estrés z = 53: las de z ≤ 37 quedan cubiertas [T]; de las
del estrés, las que tengan L/r ≤ 343 también. Para
z ≥ 53 el teorema cubre el sub-rango L/r ∈ [z, 343] de la
capa, y deja fuera la rebanada L/r ∈ (343, z^{3/2}),
u_r ∈ (log 343/log z, 3/2) — para z = 53, u_r ∈ (1,472; 1,5);
para z = 71, u_r ∈ (1,369; 1,5).

**La configuración mínima se realiza [M, K1].** El párrafo
anterior nombraba la amenaza mínima: j = 2, (q₁, q₂) =
(7, 5), d′ = 35, fallo 25·R*·35 ≥ L con R* | r. La caza
extendida la encuentra, con esa forma exacta y ninguna
otra: Rem(r) = {35} para r = 11 (z = 61, L = 3926;
z = 67, L = 4279, 4573, 4911; z = 71, L = 4707, 5035,
5402) y r = 13 (z = 67, L = 4573, 4911; z = 71, L = 4707,
5035, 5402) — 12 instancias, siempre con fallo
(m, j, p_m) = (3, 2, 5) (la lista fusionada de 35r es
(r, 7, 5) y falla en m = 3: 25·35r ≥ L) y siempre
E_Rem(r) = −a(35) = −1/24 < 0. En las 70 instancias de
residuo probadas no aparece ningún otro elemento de Rem.
Por tanto **(EREM-HI) — Rem(r) = ∅ en toda la capa — es
FALSO** [M, contraejemplos verificados], y la provisión
(i) del doc 120a §2 / 122a §1 no puede promoverse a [T]
tal cual. Configuraciones con j ≥ 4 exigen ω(d′) ≥ 4,
luego q_1 ≥ 11 y L/r > 1331: imposibles en la ventana
extendida (z ≤ 71 da L/r < 598), vivas solo para z > 121.

**Lema 4 (el signo sobrevive) [T].** Si L/r ≤ 1331, todo
d′ ∈ Rem(r) tiene ω(d′) = 2. En efecto, ω(d′) ≥ 2 por el
corolario del Lema 2, y ω(d′) ≥ 3 forzaría q_1 ≥ 11 (los
únicos primos permitidos bajo 11 son 5 y 7: no hay tres
distintos con q_1 ≤ 7), luego L/r > q_1³ ≥ 1331 por el
Lema 3. Con ω(d′) = 2, μ(d′) = +1 y a(d′) = 1/φ(d′) > 0,
así que **E_Rem(r) = −Σ_{Rem(r)} a ≤ 0 para todo conductor
de capa con L/r ≤ 1331** — en particular en todo marco
con z ≤ 121, pues la capa tiene L/r < z^{3/2}. ∎

**Consecuencia para la cadena [T].** La identidad del
doc 122a, δ(r) = E_Rem(r) + M(r) − P(r), con la provisión
de signo (i′) E_Rem(r) ≤ 0 en lugar de (i), sigue dando
δ(r) ≤ M(r) − P(r) ≤ M(r) y δ(r)⁺ ≤ M(r) puntualmente:
la reducción (OUT-MASS) ⟹ (BAND-δ-HI′) no necesita
E_Rem = 0, solo su signo, y ese signo es [T] en todo el
rango anterior. Simétricamente, con (ii) E_Add ≤ 0,
δ⁺ ≤ ((Ŝ−Π) − |E_Add|)⁺ se mantiene bajo (i′).

**El residuo mínimo (EREM-SIGN) [P].** Para todo
L ∈ ℒ_A y todo r no degenerado de la capa con
L/r > 1331 (posible solo si z > 121): E_Rem(r) ≤ 0,
equivalentemente Σ_{d′∈Rem(r)} a(d′) ≥ 0. Por el Lema 2
un contraejemplo exige un elemento de Rem con masa
negativa neta; la configuración mínima es ω(d′) = 3
(μ = −1) con j = 2 — p. ej. d′ = 385 = 5·7·11, que por sí
solo daría E_Rem = +1/240 > 0. Nada en la ventana actual
(z ≤ 71) toca este rango: ahí (EREM-SIGN) es [T] por el
Lema 4.

## 2. Medición: reproducción, caza extendida y falsación del Lema 2 [M]

Script §K0–K2. **K0 (reproducción + |Rem|).** Los marcos
estándar reproducen exactamente las 67 instancias no
degeneradas de capa de los docs 120a/122a: E_Rem = 0 y
Rem = ∅ en las 67, con el split que el Teorema predice —
66 con L/r ≤ 343 (ahora [T]) y 1 con L/r > 343 (estrés
z = 53; su E_Rem = 0 queda [M]). De 58 instancias
degeneradas de capa probadas, Rem = ∅ en todas [T por el
corolario, verificado]; el cruce nhat = 1 ⟺ q₀³r ≥ L no
falla nunca.
**K1 (caza extendida en el residuo).** Marcos A con
z ∈ {53, 61, 67, 71}, niveles cerca de u_L ≈ 2 y u_L ≈ 3:
70 instancias de residuo (L/r > 343) probadas — 3, 11, 23
y 33 por z — y 189 instancias del sub-rango cubierto como
control (Rem = ∅ en todas, como manda el Teorema). En el
residuo aparecen exactamente 12 instancias con Rem ≠ ∅,
todas con Rem = {35}, r ∈ {11, 13}, fallo (3, 2, 5) y
E_Rem = −1/24 (detalle en §1); E_Rem > 0 en **ninguna**
instancia probada, en ningún rango.
**K2 (falsación estructural).** En niveles de z = 37 y del
estrés z = 53 se probaron todos los conductores de banda
(cualquier u_r) y una muestra determinista de visibles
fuera de banda: el Lema 1 se verifica en los 56 visibles
probados (ninguno fuera de V⁺(L)) y no se encontró ningún
elemento de Rem en esos niveles — la falsación del Lema 2
es vacua ahí, y su contraste no vacuo es K1: los 12
elementos realizados cumplen (m, j) = (3, 2), par y ≥ 2,
como el lema exige.

## 3. Estrato 5 de (SLAB-V′): re-derivación del cierre del doc 124 [T vía (POS)]

Notación del doc 123a: r degenerado visible en banda de
estrato q₀ = 5, testigo d = cr ∈ V⁺(L), lista fusionada
p_1 > … > p_k, y la familia **(POS_m)**: c·p_m² < q₀³·∏_{j>m}
p_j [T, doc 123a §2]. Los primos de c son ≥ 5, c es libre de
cuadrados, (c, r) = 1 y 5 ∉ r. Las dos posiciones útiles:
(POS_k) si k es impar, (POS_{k−1}) si k es par:
c·p_{k−1}² < 125·p_k. Hecho elemental usado dos veces: la
lista fusionada tiene primos **distintos** en orden
estrictamente decreciente, luego p_{k−1} > p_k y por tanto
p_{k−1} ≥ el siguiente primo tras p_k — en particular
**p_{k−1} ≥ 11 si p_k = 7** y **p_{k−1} ≥ 13 si p_k = 11**.
(El doc 123a usó solo p_{k−1} ≥ p_k en un paso; de ahí su
caso residual.)

**Caso k impar — (POS_k): c·p_k² < 125.** Si 5 | c, el 5
está en d y es su menor primo posible: p_k = 5, y
25c < 125 da c < 5, imposible con 5 | c. Luego 5 ∤ c,
p_k ≥ 7 y 49c < 125 da c < 2,56: **c = 1**.

**Caso k par, 5 | c — (POS_{k−1}) con p_k = 5.**
p_{k−1} ≥ 7, luego c < 625/49 < 12,76. Con 5 | c y c libre
de cuadrados, la única opción < 12,76 es **c = 5**
(5·7 = 35 ya excede).

**Caso k par, 5 ∤ c — (POS_{k−1}): c·p_{k−1}² < 125·p_k,**
con p_k ≥ 7 (ni c ni r contienen 5).
- *p_k = 7:* p_{k−1} ≥ 11, luego
  c < 125·7/121 = **875/121 < 7,24**. Como 5 ∤ c, los
  candidatos son c = 1, 7, o ≥ 11: queda **c ∈ {1, 7}**
  (un c compuesto ≥ 7·11 = 77 también excede).
- *p_k ≥ 11:* el candidato peligroso es c = 11. Si
  c = 11, entonces 11 | d y p_k ≤ 11, luego p_k = 11 y
  p_{k−1} ≥ 13: c < 125·11/169 = **1375/169 < 8,14 < 11**,
  contradicción — **c = 11 excluido**. (Este es exactamente
  el caso residual del doc 123a §2 — k par, 5 ∤ c, menor
  primo de r ≥ 11 — cerrado: con p_{k−1} ≥ p_k = 11 solo
  salía c < 125/11 < 11,37.) Los demás candidatos: c = 7
  forzaría p_k ≤ 7, fuera de este sub-caso; c ≥ 13 muere
  por c < 125·p_k/p_{k−1}² < 125/p_k ≤ 125/11 < 11,37 < 13
  (y de hecho la cota uniforme de este sub-caso es
  1375/169: para primos consecutivos p ≥ 13,
  125p/q² < 125·11/169 porque 11p² − 125p + 44 > 0 —
  verificado en K3 sobre todos los pares consecutivos
  hasta 10⁵). Luego **c = 1** en este sub-caso.

**Conclusión [T vía (POS)]: en el estrato 5, c ∈ {1, 5, 7}**
— (WIT)-estrato-5 es teorema. La enumeración extendida
(marcos A z = 61, 67 y B z = 61, evaluador del doc 123a
validado contra la maquinaria real) lo confirma sin
excepción [M, K4]: en los tres marcos el conjunto de
cofactores del estrato 5 es exactamente {1, 5, 7}, con
N máx = 3 y ninguna violación. Los chequeos aritméticos de
las cotas están en K3: 875/121 = 7,231405 < 7,24;
1375/169 = 8,136095 < 8,14; 625/49 = 12,755102; el supremo
de 125p/q² sobre pares de primos consecutivos con p ≥ 11
(hasta 10⁵) es 1375/169, alcanzado solo en p = 11;
11p² − 125p + 44 > 0 para todo primo p ≥ 13; y ningún
candidato libre de cuadrados c < 125 fuera de {1, 5, 7}
sobrevive a las cotas de las ramas.

## 4. Errata formal del doc 123 (corrección 2 del doc 124)

**Texto que se corrige** (doc 123b §3 y su eco en §5): la
lectura «ω(r) = 1 carga el 99,9 % de la masa y
u_r ∈ [1,35; 3/2) casi toda» como si el régimen dominante
— conductores primos en la rebanada alta de la banda, donde
M_near = M_near⁽¹⁾ íntegro con la fórmula exacta de la §1
del doc 123b — cubriera casi toda la masa de (OUT-MASS).

**Cifras correctas** (doc 124, sobre el pooled del doc
123b: SM = 4,0720): marginal de primos ω(r) = 1:
**99,8993 %**; marginal de la rebanada u_r ∈ [1,35; 3/2):
**92,7157 %**; **conjunta primo ∩ rebanada: 92,6537 %**.
Descomposición exacta del resto: primos con u_r < 1,35:
99,8993 − 92,6537 = **7,2456 %**; conductores compuestos:
**0,1007 %**. (92,6537 + 7,2456 + 0,1007 = 100,0000.)

**Consecuencia exacta para la ruta exacta del doc 123b §4.**
- **(PRIME-WINDOW) no cambia como enunciado**: cuantifica
  sobre *todo* p₀ primo de la capa, y la marginal de primos
  sigue siendo 99,8993 % — la partición
  (PRIME-WINDOW) + (R3) del lado de conductores se mantiene.
- **Lo que se corrige es su cobertura «con fórmula
  exacta»**: el régimen donde las cinco instancias
  dominantes viven y donde R1 se midió nulo es la conjunta,
  el **92,6537 %**, no el 99,9 %.
- **(R1) crece en alcance**: los ausentes cercanos con
  ω(d′) ≥ 3 dejan de ser «nulos en las instancias
  dominantes» a secas; son un resto vivo sobre el 7,2456 %
  de masa en primos con u_r < 1,35 (donde L/r es mayor y el
  argumento de tamaño del doc 123b §1 no los anula), además
  de su papel general cuando L/r supere el triple mínimo
  coprimo.
- **(R2) sin cambios** (la masa lejana tiene forma cerrada;
  su suma sobre la capa es tipo (LAYER-COUNT)).
- **(R3) queda medido con precisión**: 0,1007 % de la masa
  pooled — «pequeño en ventana» con cifra, no con adverbio.

## 5. Veredicto

- **E_Rem: RESUELTO en la ventana, con una corrección.**
  E_Rem(r) = 0 es [T] para L/r ≤ 343 (Rem = ∅; Lema de
  paridad + umbral) y para todo r degenerado; la forma
  uniforme «E_Rem = 0 en la capa» es **FALSA** — refutada
  por Rem = {35} en r ∈ {11, 13}, z ∈ {61, 67, 71}
  [M] — pero la provisión que la cadena necesita es el
  signo, y **E_Rem(r) ≤ 0 es [T] para L/r ≤ 1331** (Lema
  4: ω(d′) = 2 forzado en Rem). Residuo [P]:
  **(EREM-SIGN)**, el signo para L/r > 1331 (z > 121),
  con amenaza mínima ω(d′) = 3 nombrada.
- **Estrato 5: CONFIRMADO** — (WIT)-estrato-5 es [T] vía
  (POS) con la corrección del doc 124 re-derivada paso a
  paso; sin discrepancias; enumeración extendida limpia
  [M, K4]. Estado de (SLAB-V′): los estratos 7 y 11 siguen
  refutados como (WIT) (doc 123a); el 5 queda cerrado como
  conjunto de cofactores {1, 5, 7} — (SLAB-V′) como cota
  sigue [P].
- **Errata del doc 123 incorporada** (§4): concentración
  conjunta 92,6537 %; alcance de (PRIME-WINDOW) y de
  R1–R3 re-delimitado con cifras.

Inventario de (H2) tras esta sección: (FC-2), (OUT-MASS)
con la provisión (i) reemplazada por (i′) E_Rem ≤ 0 —
[T] para L/r ≤ 1331 por el Lema 4, con E_Rem = 0 [T] en
el sub-rango L/r ≤ 343 — y el residuo (EREM-SIGN) [P]
solo para L/r > 1331, (BAND-LOW-NONDEG), (SLAB-V′) con
el estrato 5 cerrado en cofactores.
**Goldbach permanece abierto.**

## 6. Reproducibilidad

`python3 125b_erem_pos.py | tee salida_doc125b_local.txt`
— solo stdlib; importa `120a_banda_delta_hi.py` (maquinaria
de los docs 110/114/117b), `109_ataque_h2_gram.py`
(`in_Vplus`) y `123a_wit.py` (evaluador ligero de
estratos), sin tocarlos. K0: reproducción de la capa de los
docs 120a/122a con |Rem| y E_Rem por instancia. K1: caza
extendida z ∈ {53, 61, 67, 71} separando L/r ≷ 343. K2:
falsación del Lema 2 fuera de la capa. K3: aritmética del
cierre del estrato 5. K4: enumeración extendida del
estrato 5.
