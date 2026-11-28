# Doc 125a — (PRIME-WINDOW): la doble suma de Mertens tiene término principal constante y Q_c decae — refutación asintótica en el régimen u_L ≈ 2

**Fecha:** 2026-11-28. **Línea:** UNIF-U / (H2), continuación de los docs 122a/123b.
**Documento hermano:** 125b (línea paralela); esta sección es autocontenida.
**Etiquetas:** [T] = teorema, [M] = medición en ventana, [P] = pendiente.

---

## 0. La pregunta

El doc 123b redujo (OUT-MASS), por la ruta exacta, a **(PRIME-WINDOW)**
[P] más los restos (R1–R3): en conductores primos de la capa, la masa
cercana de los ausentes de ambos mundos con ω(d′) = 1 es la suma exacta
M_near(p₀) = Σ_{q ∈ 𝒬(p₀)} 1/(q−1) [T, doc 123b §1], y esa pieza carga
la inmensa mayoría de la masa medida (ω(r) = 1: 99,9 % pooled; la
conjunta primo∩banda [1,35; 3/2): 92,6537 %, auditoría del doc 124).
Esta sección: (1) enuncia (PRIME-WINDOW) sin ambigüedad; (2) evalúa la
suma de ventana por Mertens sobre primos, con el error que la serie
puede pagar; (3) suma sobre la capa y decide el enunciado; (4) mide el
cociente contra Q_c en los marcos extendidos del doc 123b (z ≤ 71).

Convenciones (docs 110/123b): A = 3, excl = ∅; primos permitidos =
impares < z, ≠ 3; nivel L ∈ ℒ_A con parámetro z;
u_r = log(L/r)/log z; Q_c = Q_c(z) = Π_{p permitido}
(1 − 1/(p−1) + 1/(p−1)²) (doc 109/110: depende solo de z);
capa(L) = banda alta no degenerada: r visible (r ∈ V⁺(L)),
u_r ∈ [1, 3/2), Π(r) < 1/2, V̂⁽ʳ⁾ ≠ {1}.

## 1. Enunciado exacto de (PRIME-WINDOW) y su lugar en (OUT-MASS)

Sea capa₁(L) = {p₀ ∈ capa(L) : p₀ primo}. Para p₀ ∈ capa₁(L),
p₀ ∈ V⁺(L) ⟺ p₀³ < L [T, Nathanson m = 1 sobre el primo mayor], y

$$\mathcal Q(p_0) = \{\,q \text{ permitido} : q > p_0,\ q^3 \ge L,\
p_0 q < L\,\}, \qquad
M_{near}^{(1)}(p_0) = \sum_{q \in \mathcal Q(p_0)} \frac{1}{q-1}
\quad [T].$$

**(PRIME-WINDOW) [P] (doc 123b §4):** para todo L ∈ ℒ_A, con
C absoluta,

$$PW(L) \;:=\; \sum_{p_0 \in \mathrm{capa}_1(L)}
\frac{1}{p_0-1}\, M_{near}^{(1)}(p_0) \;\le\; C\cdot Q_c(z).$$

Descomposición de la suma total de (OUT-MASS), SM(L) :=
Σ_{r ∈ capa} M(r)/φ(r) [T, por las definiciones del doc 123b]:

SM(L) = PW(L) + R1(L) + R2(L) + R3(L), con

- R1(L) = Σ_{p₀ ∈ capa₁} (p₀−1)⁻¹ M_near^{ω≥3}(p₀): ausentes cercanos
  con ω(d′) ≥ 3 en conductores primos (nulos en las instancias
  dominantes de ventana por tamaño [M, doc 123b §1]);
- R2(L) = Σ_{r ∈ capa} M_far(r)/φ(r): la masa lejana, con forma
  cerrada por productos de Euler [T, doc 123b §1]; fracción pooled
  de M: 0,074–0,181 en los marcos extendidos [M];
- R3(L) = Σ_{r ∈ capa, ω(r) ≥ 2} M_near(r)/φ(r): la cercana en
  conductores compuestos; 0,1 % pooled en ventana [M].

En particular PW(L) ≤ SM(L) término a término [T]: PW es una
sub-suma con todos los sumandos ≥ 0. Lo que se decida sobre PW
con C absoluta se hereda como cota inferior de SM.

## 2. La suma de ventana: Mertens sobre primos, con el tope q < z

Sea S(x) := Σ_{q permitido, q < x} 1/(q−1) (x real; el conjunto
permitido sin el tope z: primos impares ≠ 3). Como q permitido
exige q < z, la suma sobre 𝒬(p₀) es una diferencia de S con
extremos truncados:

**Lema de ventana [T].** Para p₀ ∈ capa₁(L) vale p₀ < L^{1/3}
(visibilidad) y L/p₀ ≥ z (porque u_{p₀} ≥ 1). Por tanto los
extremos efectivos de 𝒬(p₀) son A = L^{1/3} y B = z — el tope
p₀q < L nunca muerde (L/p₀ ≥ z > q) y el tope q > p₀ tampoco
(q ≥ L^{1/3} > p₀) — y

$$M_{near}^{(1)}(p_0) \;=\; S(z) \;-\; S\!\left(L^{1/3}\right)
\qquad [T],$$

idéntico para todo p₀ de la capa al mismo nivel (L, z). En
particular la fórmula ingenua log log(L/p₀) − log log(max(p₀,
L^{1/3})) es incorrecta sobre la capa: olvida el tope q < z del
conjunto permitido, que es el extremo que manda (L/p₀ ≥ z
siempre en la capa).

**Mertens.** S(x) = log log x + M′ + E(x) con E(x) → 0, donde
M′ = M − 3/2 + Σ_{q ≥ 5} 1/(q(q−1)): en efecto,
Σ_{p < x} 1/p = log log x + M + O(1/log x) (Mertens II; el error
O(1/log x) se sigue de Σ_{p≤x} (log p)/p = log x + O(1), fuerza
Chebyshev, por sumación parcial), 1/(q−1) = 1/q + 1/(q(q−1))
con Σ 1/(q(q−1)) absolutamente convergente (cola
O(1/(x log x))), y quitar los primos 2 y 3 resta 1 + 1/2.
Basta esta forma clásica — de hecho basta E(x) = o(1) — para
todo lo que sigue; ninguna forma fuerte de Mertens hace falta.
Con ella:

$$M_{near}^{(1)}(p_0) \;=\; \log\frac{\log z}{\log L^{1/3}}
\;+\; \bigl[E(z) - E(L^{1/3})\bigr]
\;=\; \log\frac{3}{u_L} \;+\; \Delta E \qquad [T],$$

con u_L = log L/log z y ΔE = O(1/log L^{1/3}) → 0 si L → ∞.
En ventana ΔE no es despreciable numéricamente (a z = 53,
L = 2585: término principal log(3/1,979) = 0,4167 contra el
valor exacto 0,3309; ΔE = −0,0858, el −21 %): la §4 lo mide.
Lo que decide el veredicto no es ΔE sino el término principal.

## 3. La suma sobre la capa: término principal constante contra Q_c ≍ c/log z

**La suma externa también es Mertens [T].** Los p₀ de capa₁(L)
son los primos permitidos en (L z^{−3/2}, min(L^{1/3}, L z^{−1})]
que además cumplen Π(p₀) < 1/2 y no degeneración. Para
u_L ∈ (3/2, 3) el extremo superior es L^{1/3} y

$$\sum_{p_0 \in \mathrm{capa}_1(L)} \frac{1}{p_0-1}
\;=\; S\!\left(L^{1/3}\right) - S\!\left(L z^{-3/2}\right)
\;+\; \varepsilon_\Pi
\;=\; \log\frac{u_L/3}{\,u_L - 3/2\,} \;+\; \Delta E'
\;+\; \varepsilon_\Pi,$$

donde ε_Π ≥ 0 es la masa de los primos excluidos solo por
Π(p₀) ≥ 1/2 o degeneración: ε_Π = 0 en toda la ventana medida
(§4), y asintóticamente Π(p₀) ≍ c/log z < 1/2 siempre y
L/p₀ > 7³ en la capa, así que ε_Π = 0 para z grande [T por el
tamaño de Π]. Como M_near⁽¹⁾ no depende de p₀ (§2):

$$PW(L) \;=\; \log\frac{3}{u_L}\;\cdot\;
\log\frac{u_L/3}{\,u_L-3/2\,} \;+\; O\!\left(\tfrac{1}{\log z}\right)
\qquad [T],\quad u_L \in (3/2,\,3)\ \text{fijo}.$$

En el régimen de los marcos, u_L ≈ 2: PW → log(3/2)·log(4/3)
= 0,4055 · 0,2877 = **0,1167**, una constante positiva. Los p₀
que hacen grande el factor interno son los de u_{p₀} → 3/2
(extremo L z^{−3/2} de la suma externa, donde log log apenas
se mueve) — pero no son pocos ni extremos: es toda la ventana
de primos (L z^{−3/2}, L^{1/3}], cuya masa de Mertens es la
constante log(4/3) íntegra. No hay sub-rango residual que
aislar: el término principal es la ventana completa.

**El presupuesto decae [T].** log Q_c(z) = Σ_{p perm}
log(1 − 1/(p−1) + 1/(p−1)²) = −Σ 1/(p−1) + O(1)
= −log log z + O(1), porque log(1−x+x²) + x = O(x²) es
absolutamente sumable y Σ_{p<z} 1/(p−1) = log log z + M₁ + o(1).
Luego **Q_c(z) ~ c_Q/log z** con c_Q > 0 (la §4 mide
Q_c(z)·log z → c_Q ≈ 1,7 en ventana).

**Conclusión [T].** En cualquier nivel con u_L ∈ (3/2, 3)
fijo — en particular los niveles u_L ≈ 2 que los marcos de la
serie incluyen para todo z y donde se alcanzan los máximos
medidos de SM/Q_c (doc 123b §3) —

$$\frac{PW(L)}{Q_c(z)} \;\sim\;
\frac{\log(3/u_L)\,\log\!\bigl(\frac{u_L/3}{u_L-3/2}\bigr)}{c_Q}
\;\log z \;\longrightarrow\; \infty,$$

y como SM ≥ PW término a término, lo mismo vale para SM/Q_c.
**(PRIME-WINDOW) con C absoluta contra Q_c es falsa en este
régimen** — no es un enunciado abierto: su término principal
es una constante positiva y su presupuesto tiende a 0 como
1/log z. El error de Mertens sumado sobre la capa es
O(1/log z) = o(1): no rescata ni condena nada; la vía Mertens
basta y sobra para decidir, y decide en contra. El crecimiento
~ log z es indetectable como crecimiento en z ≤ 71 (log z solo
recorre 3,61 → 4,26), lo que explica la planitud medida de
C(z) ∈ [0,196; 0,315] sin contradecirla: es el comienzo de la
recta (0,1167/c_Q)·log z ≈ 0,069·log z.

Alcance exacto de la refutación: muere (PRIME-WINDOW) y, con
ella, (OUT-MASS) tal como está enunciada (C absoluta, escala
Q_c) en el régimen u_L ≈ 2, y la ruta exacta completa del
doc 123b. **No** muere (BAND-δ-HI′): δ⁺ ≤ M era la cota
holgada que descartaba la cancelación M↔P (factor pooled
2,26, doc 122a); la ruta con signo — (NET-LAYER) + (P-MASS),
doc 123b §2 — queda intacta, así como la ruta puntual
(LAYER-COUNT) en tanto enunciado sobre Σ 1/φ (cuyo calibre
medido tampoco crecía: 0,45–0,82; su comparación es contra
Q_c con la cota puntual B_r/2, y su estatus tras este doc es
el de siguiente frente, no el de corolario). Fuera del régimen
u_L ∈ (3/2, 3): si u_L → 3 el factor log(3/u_L) → 0 y el
análisis cambia; los marcos no viven ahí.

## 4. Medición de control [M]

Script `125a_prime_window.py` (solo stdlib; importa la
maquinaria de los docs 110/120a/123b sin tocarla), salida
`salida_doc125a_local.txt`. Marcos extendidos del doc 123b
(D = z⁵, y = z^{10/3}, cap = 4 niveles cerca de u_L = 2, más
L = z⁵//5), z ∈ {37, 53, 61, 67, 71}.

**K0 [M].** Los cinco valores K4 del doc 123b se reproducen
exactos: M_near = 0,415891 / 0,330884 / 0,367356 / 0,384023 /
0,399174 contra 0,4159 / 0,3309 / 0,3674 / 0,3840 / 0,3992.

**K1 [M, verifica T].** El lema de ventana de la §2 vale a
precisión de máquina sobre todas las instancias primas de
capa de los marcos: |M_near − (S(z) − S(L^{1/3}))| ≤ 5,6·10⁻¹⁷
y la dispersión de M_near entre los p₀ de un mismo nivel es
exactamente 0 — M_near es constante sobre la capa, como dice
el lema. Los conjuntos de primos de capa por el criterio
propio coinciden con los de la maquinaria del doc 123b en
todos los niveles, y Q_c propio = lv.qc a error relativo
< 10⁻¹². El error de Mertens: E(z) = S(z) − log log z ∈
[−0,452, −0,446] (estable, camino de M′) y ΔE = M_near −
log(3/u_L) ∈ [−0,0907; +0,0017]: en ventana ΔE ≈ −0,09
recorta el término principal un ~21 %, como anticipó la §2;
es o(1) y no toca el veredicto.

**K2 [M].** PW por nivel (máximos por z; entre paréntesis el
nivel): PW/Q_c = 0,1797 (37, L=1437), 0,2722 (53, L=2585),
0,1641 (61, L=3450), 0,1744 (67, L=4024), 0,2081 (71,
L=5035) — es decir **PW carga el 80–92 % del C(z) del
doc 123b a su nivel máximo**, y PW/SM vale 0,80–0,94 en esos
niveles (pooled por z: 0,81–0,93). La sub-banda [1,35; 3/2)
carga el 55–100 % de PW según el nivel. La predicción de
Mertens por nivel, log(3/u_L)·log[(u_L/3)/(u_L−3/2)]
(0,095–0,141 en los niveles medidos), sobreestima el PW
exacto (0,053–0,116) justo por los ΔE, ΔE′ ≈ −0,09 aún
visibles a esta escala; al nivel máximo de z = 53 el PW
exacto ya es 0,1158, a un 1 % de la constante asintótica
0,1167. El nivel grande L = z⁵//5 da PW = 0 en todo z
(capa₁ vacía), como en el doc 123b.

**K3 [M].** Q_c(z)·log z = 1,693 / 1,689 / 1,687 / 1,697 /
1,695 (z = 37…71) y 1,704 / 1,703 / 1,706 / 1,712 / 1,710 /
1,718 (z = 101…1009): converge a c_Q ≈ 1,72, confirmando
Q_c(z) ~ c_Q/log z. La recta asintótica predicha para el
cociente, (0,1167/c_Q)·log z ≈ 0,068·log z (0,245–0,290 en
la ventana), pasa por el rango medido de max PW/Q_c
(0,164–0,272, con el ruido de nivel discreto ya documentado
en el doc 123b §3): la medición es consistente con la
asintótica de la §3 y con que el crecimiento ~ log z sea
indetectable como tendencia a z ≤ 71.

## 5. Veredicto

**(PRIME-WINDOW): REFUTADA asintóticamente [T] en el régimen
u_L ≈ 2** — su suma tiene la asintótica exacta
log(3/u_L)·log[(u_L/3)/(u_L − 3/2)] + O(1/log z), una
constante positiva (0,1167 en u_L = 2), mientras
Q_c(z) ~ c_Q/log z → 0; el cociente crece como
(0,1167/c_Q)·log z. El obstáculo exacto no es un error de
Mertens ni un sub-rango de p₀: es el término principal íntegro
de la doble suma contra un presupuesto que decae. En ventana
la refutación es invisible como crecimiento (C(z) plano) pero
visible como nivel: PW ya carga la casi totalidad del C(z)
medido a los niveles máximos. (OUT-MASS) con C absoluta
contra Q_c cae con ella en este régimen; (BAND-δ-HI′) y la
ruta con signo (NET-LAYER) + (P-MASS) no quedan tocadas.
**Goldbach permanece abierto.**

Inventario de (H2) tras esta sección: (FC-2), (OUT-MASS)
reducida a {(LAYER-COUNT)} / {(NET-LAYER), (P-MASS)} — la
ruta exacta {(PRIME-WINDOW), (R1–R3)} cerrada por refutación
del término dominante —, (BAND-LOW-NONDEG), (SLAB-V′).

## 6. Reproducibilidad

Script `125a_prime_window.py`, salida
`salida_doc125a_local.txt`: K0 reproduce los cinco valores
K4 del doc 123b; K1 verifica el lema de ventana (identidad
exacta y constancia de M_near en p₀) sobre todas las
instancias primas de capa de los marcos; K2 mide PW/Q_c,
PW/SM y la predicción de Mertens por nivel; K3 tabula
Q_c(z)·log z.
