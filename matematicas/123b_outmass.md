# Doc 123b — (OUT-MASS): estructura de los ausentes de ambos mundos, la ruta con signo y el intento de cota

**Fecha:** 2026-11-28. **Línea:** UNIF-U / (H2), continuación del doc 122a.
**Documento hermano:** 123a (línea paralela); esta sección es autocontenida.
**Etiquetas:** [T] = teorema, [M] = medición en ventana, [P] = pendiente.

---

## 0. La pregunta

El doc 122a redujo la banda alta de (H2), vía la identidad
δ(r) = E_Rem(r) + M(r) − P(r) [T], al residuo **(OUT-MASS)** [P]:

$$\sum_{r \in \mathrm{capa}(L)} \frac{M(r)}{\varphi(r)} \;\le\; C\cdot Q_c(L)
\qquad \forall\, L \in \mathcal L_A,$$

donde la capa es la banda alta no degenerada (u_r ≥ u₀* = 1,
V̂⁽ʳ⁾ ≠ {1}) y M(r), P(r) son las masas |a| = 1/φ por signo de
**Out(r)**, el complemento de V⁽ʳ⁾ ∪ V̂⁽ʳ⁾ entre los d′ coprimos
con r. En ventana (doc 122a) vale con C ≈ 0,32 descartando la
cancelación M↔P (55,8 %). Esta sección hace cuatro cosas:
(1) caracteriza Out(r) combinatoriamente y descompone M(r) en
piezas con forma cerrada; (2) explota la ruta con signo
Σ(M−P) = Σδ − ΣE_Rem y mide qué parte de ΣM/φ es cancelación y
qué parte residuo neto; (3) extiende la medición (z y niveles)
para decidir si la constante 0,3152 es estable, logarítmica o
potencial, y dónde se concentra M; (4) intenta la cota [T] y
nombra el lema mínimo que falta.

Convenciones: marco del doc 110 (A = 3, excl = ∅, primos
permitidos = impares < z, ≠ 3); enum = libres de cuadrados sobre
los permitidos, factores en orden decreciente; V⁺(L) por el
Teorema 9.3 de Nathanson (β = 2): d ∈ V⁺(L) sii d < L y para todo
m impar, p_m²·(p_1···p_m) < L sobre la lista decreciente de d.

## 1. Estructura de M(r): quiénes son los ausentes de ambos mundos

**Los dos mundos como fibras de V⁺ [T].** Por definición del
marco, V⁽ʳ⁾ = {d′ : (d′,r) = 1, r·d′ ∈ V⁺(L)}. El mundo modelo
tiene la misma forma a nivel reducido:

$$\hat V^{(r)} \;=\; \{\,d' : (d',r)=1,\ d' \in V^{+}(L/r)\,\}
\qquad [T].$$

En efecto, in_Vhat(pf, L, r) exige q₁·r < L y, para todo m
impar, r·q_m²·Q_m < L sobre la lista decreciente de d′; dividiendo
por r es exactamente in_Vplus(pf, L/r) leído sobre ℚ, y
d′ ∈ V⁺(L/r) ya implica r·d′ < L, así que la condición extra
r·d′ < L de la clasificación es automática. (Verificado sin
discrepancias en K0, sobre instancias de capa y muestra de
visibles.) Con λ(d′) = max(d′, max_{m impar} q_m²Q_m) (doc 113) y
Λ(r,d′) = max(r·d′, max_{m impar} p_m²P_m) sobre la lista
fusionada de r·d′:

$$d' \in \mathrm{Out}(r) \iff (d',r)=1,\quad r\lambda(d') \ge L,
\quad \Lambda(r,d') \ge L \qquad [T].$$

El ausente de ambos mundos es el d′ cuyo propio λ ya desborda el
nivel reducido L/r **y** cuyo producto con r desborda las
condiciones de Nathanson fusionadas al nivel L.

**Reparto lejos/cerca y formas cerradas [T].** Out(r) se parte
en *lejos* (r·d′ ≥ L) y *cerca* (r·d′ < L). Sobre la capa
(u_r ≥ 1) se tiene L/r ≥ z, luego todo primo d′ = q < z cumple
r·q < L: **ningún ausente lejano tiene ω(d′) = 1 en la capa**
[T]. Las masas totales sobre (d′,r) = 1 son productos de Euler:
Σ a(d′) = Π(r) y Σ |a(d′)| = B_r := Π_{p permitido, p∤r}
(1 + 1/(p−1)); por tanto la masa total μ = −1 es (B_r − Π(r))/2
y la μ = +1 es (B_r + Π(r))/2. Con T_<(r) = Σ_{r·d′<L} a(d′) y
U_<(r) = Σ_{r·d′<L} 1/φ(d′), las piezas lejanas tienen forma
cerrada:

$$M_{far} = \tfrac12\bigl[(B_r - U_<) - (\Pi(r) - T_<)\bigr],
\qquad
P_{far} = \tfrac12\bigl[(B_r - U_<) + (\Pi(r) - T_<)\bigr]
\quad [T],$$

verificadas a error ≤ 3,8·10⁻¹⁶ (K0). La pieza cercana no tiene
forma cerrada: es el desacuerdo local entre la fibra de V⁺(L)
y V⁺(L/r). Y es la que manda: la fracción lejana pooled de M es
0,074 (z=37), 0,145 (53), 0,160 (61), 0,165 (67), 0,181 (71)
[M] — la cercana carga el 82–93 %.

**Conductores primos: la pieza cercana es una suma sobre primos
exacta [T].** Sea r = p₀ primo visible; entonces p₀³ < L (la
condición m = 1 de cualquier elemento de V⁺(L) divisible por p₀
la fuerza sobre su mayor primo, y p₀ ≤ ese primo). Para
d′ = q primo permitido, q ≠ p₀: si q < p₀, la lista fusionada
es (p₀, q) y su única condición impar es p₀³ < L — r·q ∈ V⁺(L)
siempre; si q > p₀, r·q ∉ V⁺(L) ⟺ q³ ≥ L, y q ∉ V⁺(L/r) ⟺
q³·p₀ ≥ L, implicada por la anterior. Luego los ausentes
cercanos con ω(d′) = 1 son exactamente

$$\mathcal Q(p_0) = \{\,q \text{ primo permitido} : p_0 < q,\ 
q^3 \ge L,\ p_0 q < L\,\}, \qquad
M_{near}^{(1)}(p_0) = \sum_{q \in \mathcal Q(p_0)} \frac{1}{q-1}
\quad [T].$$

En las cinco instancias dominantes de la ventana extendida
(§3) no existe *ningún* ausente cercano con ω(d′) ≥ 3: ahí
L/r < 455 ≤ todo producto de tres primos permitidos coprimos
con r, así que M_near = M_near⁽¹⁾ íntegro [T por tamaño, en
esas instancias]. Valores [M, K4]: 0,4159 en (z=37, L=1286,
r=7) — exactamente 1/10+1/12+1/16+1/18+1/22+1/28+1/30 —;
0,3309 en (53, 2585, 7); 0,3674 en (61, 3450, 11); 0,3840 en
(67, 4024, 11); 0,3992 en (71, 4421, 11); el conjunto 𝒬
reproducido por el criterio coincide exactamente en las cinco.
La pieza lejana en esas instancias vale 0,0250 / 0,0451 /
0,0683 / 0,0738 / 0,0790, toda con ω(d′) ≥ 3.

**P(r), breve [T].** La asimetría es estructural: un primo
tiene μ = −1, así que P no tiene pieza ω(d′) = 1 — su masa
cercana vive en ω(d′) par (≥ 2) y su pieza lejana tiene la
forma cerrada de arriba. La cancelación M↔P del doc 122a no es
un espejo primo a primo: M cercana es una suma sobre primos;
P cercana es una suma sobre pares y cuádruplos.


## 2. La ruta con signo: Σ(M−P) sobre la capa

Por nivel, la identidad del doc 122a sumada da [T]

$$\sum_{r \in \mathrm{capa}} \frac{M(r)-P(r)}{\varphi(r)}
= \sum_{r \in \mathrm{capa}} \frac{\delta(r)}{\varphi(r)}
- \sum_{r \in \mathrm{capa}} \frac{E_{Rem}(r)}{\varphi(r)},$$

verificada a error ≤ 4,5·10⁻¹⁶ en todos los niveles medidos
(K1/K2). La pieza Σδ admite una identidad global nueva:

$$\sum_{r \in \mathrm{enum}} \frac{\delta(r)}{\varphi(r)} = 0
\qquad [T].$$

Prueba: Σ_r S(r)/φ(r) = Σ_r C(r)/μ(r) = Σ_{d∈V⁺} a(d)·
Σ_{r|d} μ(r) = a(1) = 1, porque la suma interior vale [d = 1];
y Σ_r Π(r)/φ(r) = 1 por el producto de Euler — por primo, la
alternativa p∤r aporta (p−2)/(p−1) y p|r aporta 1/(p−1), que
suman 1. (Verificado en K0: ambas sumas dan 1,000000000000 y
su diferencia ≤ 1,8·10⁻¹⁶ en cada nivel.) Consecuencia [T]:

$$\sum_{r \in \mathrm{capa}} \frac{\delta(r)}{\varphi(r)}
= - \sum_{r \in \mathrm{enum} \setminus \mathrm{capa}}
\frac{\delta(r)}{\varphi(r)}.$$

La suma con signo de la capa es *menos* la del complemento
(r = 1, degenerados, banda baja u_r < 1, u_r ≥ 3/2, Π(r) ≥ 1/2,
no visibles). Ningún objeto nombrado de la línea acota hoy esa
suma del complemento: la pieza con signo queda **reducida, no
cerrada**, al enunciado (NET-LAYER) [P] — Σ_capa δ/φ ≤ C₁·Q_c,
equivalente a acotar el complemento. La otra pieza,
Σ_capa E_Rem/φ, es 0 en ventana [M] (provisión (i) del
doc 122a; no es [T]).

**La partición y lo que compra [M].** Por nivel vale la
identidad exacta SM = Snet + SP (masas /φ sumadas sobre la
capa), así que (OUT-MASS) ⟸ (NET-LAYER) + (P-MASS), con
(P-MASS) [P] el espejo μ = +1: Σ_capa P/φ ≤ C₂·Q_c. En
ventana la cancelación es sustancial pero no decisiva:
SP/SM pooled = 0,629 (z=29 A), 0,494 (37), 0,376 (29 B),
0,757 (53 estrés); en los marcos extendidos, 0,517–0,759.
El máximo sobre niveles del neto vale en cada marco
estándar 0,0198 / 0,0933 / 0,0037 / 0,0867 contra el máximo
del bruto 0,0533 / 0,1956 / 0,0059 / 0,3152: la ruta con
signo abarata el objetivo un factor 1,6–3,6 a escala de los
máximos, a cambio de dos residuos abiertos, uno de ellos
(P-MASS) del mismo tipo que el original. Ganancia real,
cierre no.

## 3. Medición extendida [M]

Marcos consistentes (D = z⁵, y = z^{10/3}, cap = 4 niveles
cerca de u_L = 2, más el nivel grande L = z⁵//5) para
z ∈ {37, 53, 61, 67, 71}, y los marcos estándar del doc 120a
como anclas (reproducidas exactas: 0,0533 / 0,1956 / 0,0059 /
0,3152). Constante C(z) = max sobre niveles de Σ_capa M/φ/Q_c:

| z | C(z) | nivel del máx | C/log z | C/(log z)² | SP/SM pooled | M_far/M pooled |
|---|---|---|---|---|---|---|
| 37 | 0,1956 | L=1529 | 0,0542 | 0,0150 | 0,517 | 0,074 |
| 53 | 0,3152 | L=2585 | 0,0794 | 0,0200 | 0,756 | 0,145 |
| 61 | 0,1955 | L=4214 | 0,0476 | 0,0116 | 0,741 | 0,160 |
| 67 | 0,2089 | L=4911 | 0,0497 | 0,0118 | 0,694 | 0,165 |
| 71 | 0,2598 | L=5402 | 0,0609 | 0,0143 | 0,759 | 0,181 |

**Lectura.** C(z) no crece en la ventana extendida: oscila en
[0,196, 0,315] sin tendencia (C(71) < C(53)); C/log z y
C/(log z)² tampoco muestran tendencia, y las pendientes
log–log entre z consecutivos oscilan (+5,0, −13,7, +2,9,
+15,9): el ruido de nivel domina cualquier ley suave a esta
resolución. Los datos son compatibles con constante acotada
o de crecimiento muy lento (≤ log log), e incompatibles con
un crecimiento potencial visible; la ventana no distingue
constante de log log. Todos los máximos se alcanzan en
niveles pequeños L ≈ z² (u_L ≈ 2); el nivel grande z⁵//5
nunca gana en ningún z.

**Concentración [M, K3].** Sobre 2127 instancias pooled de
los marcos extendidos (SM pooled = 4,0720): el top-5 carga el
16,4 % y el top-20 el 56,3 %. Por ω(r): **ω(r) = 1 carga
4,0679 de 4,0720 — el 99,9 %**; ω = 4 carga 0,0041 y el
resto es despreciable. Por banda de u_r: [1,35, 1,45) carga
3,0458 y [1,45, 3/2) carga 0,7296: la masa vive pegada al
borde superior de la banda, en conductores *primos*, en
niveles u_L ≈ 2 — exactamente el régimen donde la fórmula
exacta de la §1 aplica (las cinco instancias top son
r = 7 y r = 11 con u_r ∈ [1,40, 1,49]).

## 4. Intento de cota [T] y el lema que falta

**Ruta puntual: cota cerrada sí, suficiencia no.** Como
Out(r) ⊆ {d′ : (d′,r) = 1}, la §1 da la cota puntual [T]

$$M(r) \;\le\; \frac{B_r - \Pi(r)}{2} \;\le\; \frac{B_{full}(z)}{2},$$

verificada sin violaciones en todas las instancias de capa
de K0. Por tanto (OUT-MASS) seguiría de

**(LAYER-COUNT) [P]:** Σ_{r∈capa(L)} (B_r − Π(r))/(2φ(r))
≤ C·Q_c(L) para todo L ∈ ℒ_A — a fortiori, de acotar
Σ_capa 1/φ(r).

(LAYER-COUNT) es un enunciado solo sobre la distribución de
conductores de la capa (visibles, u_r ∈ [1, 3/2), Π(r) < 1/2,
pesados por 1/φ(r)): no interviene ninguna fibra. Su calibre
en ventana no crece: max sobre niveles de Σ_capa 1/φ/Q_c =
0,5687 (37), 0,8225 (53), 0,4467 (61), 0,4541 (67), 0,6182
(71) [M]. Con B_full/2 ≈ 1,2–1,3 la ruta puntual daría
(OUT-MASS) con C ≈ 1,0 a escala de ventana — el orden
correcto pero un factor ~3,4 más flojo que el 0,3152 medido —
*si* (LAYER-COUNT) fuera teorema. No lo es: ningún objeto de
la línea (ni (FC-2), ni (FLOW-COMP), ni (SLAB-V′)) implica un
conteo de conductores de capa, y el doble conteo
Σ_capa M(r)/φ(r) = Σ_{d′: μ=−1} |a(d′)|·Σ_{r ∈ capa, d′∈Out(r)}
1/φ(r) traslada el problema a la multiplicidad de cada
ausente sobre la capa, que es la misma dificultad vista desde
el otro lado.

**Ruta exacta: el lema mínimo con nombre.** Por las §§1 y 3,
en ventana (OUT-MASS) lo carga una sola pieza con fórmula
exacta. El residuo mínimo por esta ruta es

**(PRIME-WINDOW) [P]:** para todo L ∈ ℒ_A,
Σ_{p₀ primo ∈ capa(L)} (1/(p₀−1)) · Σ_{q ∈ 𝒬(p₀)} 1/(q−1)
≤ C·Q_c(L), con 𝒬 el conjunto explícito de la §1;

más tres restos a acotar por separado [P]: (R1) los ausentes
cercanos con ω(d′) ≥ 3 — nulos en las instancias dominantes
por tamaño, reales cuando L/r supere el triple mínimo
coprimo; (R2) la masa lejana — con forma cerrada (§1) pero
cuya suma sobre la capa es del mismo tipo que
(LAYER-COUNT); (R3) los conductores compuestos — el 0,1 %
de la masa en ventana, pero «pequeño en ventana» no es una
cota. La suma interna de (PRIME-WINDOW) es una suma de
Mertens sobre primos en la ventana explícita
(max(p₀, L^{1/3}), L/p₀): su tamaño natural es una diferencia
de log log entre extremos de la ventana, y la externa es otra
suma sobre primos en (L·z^{−3/2}, L·z^{−1}]: el lema es una
doble suma de Mertens sobre ventanas explícitas — plausible
a la escala medida, no demostrado, y no implicado por
ningún enunciado actual de (H2).

**Dónde se detiene el rigor.** (OUT-MASS) queda reducido a
residuos estrictamente menores y mejor nombrados —
(LAYER-COUNT) por la ruta puntual, (PRIME-WINDOW) + (R1–R3)
por la ruta exacta, (NET-LAYER) + (P-MASS) por la ruta con
signo — pero ninguno de los tres caminos cierra hoy: los
tres terminan en sumas de Mertens sobre conjuntos delgados
(la capa, las ventanas de primos) para las que la línea aún
no tiene un teorema de conteo.

## 5. Veredicto

**(OUT-MASS): reducida a residuos menores nombrados, con la
parte con signo cerrada como identidad y la estructura del
término dominante exacta; apoyada en ventana extendida sin
crecimiento de la constante. No demostrada; no refutada.**

- Estructura [T]: V̂⁽ʳ⁾ es la fibra de V⁺(L/r); Out(r) se
  caracteriza por r·λ(d′) ≥ L y Λ(r,d′) ≥ L; las piezas
  lejanas tienen forma cerrada por productos de Euler y en
  la capa ningún ausente lejano tiene ω(d′) = 1.
- El término dominante tiene fórmula exacta [T]: en
  conductores primos, M_near es la suma de Mertens
  Σ_{q ∈ 𝒬(p₀)} 1/(q−1) sobre la ventana explícita
  𝒬(p₀) = {q > p₀, q³ ≥ L, p₀q < L} — y en ventana ese
  término *es* el problema: ω(r) = 1 carga el 99,9 % de la
  masa pooled y u_r ∈ [1,35, 3/2) casi toda.
- Parte con signo [T]: Σ_{r ∈ enum} δ(r)/φ(r) = 0
  exactamente; la suma de la capa es menos la del
  complemento. La partición SM = Snet + SP abarata el
  objetivo un factor 1,6–3,6 a escala de los máximos [M],
  a cambio de (NET-LAYER) y (P-MASS), ambos [P].
- Medición extendida [M]: C(z) ∈ [0,196, 0,315] para
  37 ≤ z ≤ 71 sin crecimiento neto; los máximos viven en
  niveles u_L ≈ 2, nunca en el nivel grande.
- Residuos mínimos, por ruta: (LAYER-COUNT) [puntual];
  (PRIME-WINDOW) + (R1) ω ≥ 3 cercano + (R2) suma lejana +
  (R3) compuestos [exacta]; (NET-LAYER) + (P-MASS) [signo].
  Los tres son sumas de Mertens sobre conjuntos delgados:
  ese es el obstáculo exacto que queda.

Inventario de (H2) tras esta sección: (FC-2), (OUT-MASS)
reducida a {(LAYER-COUNT)} / {(PRIME-WINDOW), (R1–R3)} /
{(NET-LAYER), (P-MASS)}, (BAND-LOW-NONDEG), (SLAB-V′).
**Goldbach permanece abierto.**

## 6. Reproducibilidad

Script `123b_outmass.py` (solo stdlib; importa la maquinaria
de los docs 110/117b/120a sin tocarla), salida
`salida_doc123b_local.txt`. PASS completo (K0–K4) en 91,1 s:
K0 identidades [T], K1 anclas del doc 122a reproducidas,
K2 tendencia extendida, K3 concentración, K4 fórmula exacta
en conductores primos.
