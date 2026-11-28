# Doc 122a — (BAND-δ-HI′): enunciado exacto con déficit Add, medición ampliada, el lema que falta y su relación con (FC-2)

**Fecha:** 2026-11-28. **Línea:** UNIF-U / (H2), continuación del doc 120a.
**Documento hermano:** ninguno; esta sección es autocontenida.
**Etiquetas:** [T] = teorema, [M] = medición en ventana, [P] = pendiente.

---

## 0. La pregunta

El doc 120a dejó (BAND-δ-HI) reformulada como **(BAND-δ-HI′)**: el
déficit del flujo Add debe dominar *en suma*, sobre la capa alta de la
banda, al exceso del modelo (Ŝ−Π)⁺. Esta sección hace cuatro cosas:
(1) escribe (BAND-δ-HI′) como desigualdad exacta con los objetos del
flujo del doc 114; (2) mide la desigualdad sumada separando exceso,
déficit y su cociente por nivel, con el subtotal de los conductores
con E = 0; (3) formula el lema mínimo suficiente y lo somete a la
ventana; (4) decide si (BAND-δ-HI′) y (FC-2) comparten el lema de
decaimiento subyacente o son independientes también a ese nivel.

Convenciones de la serie: nivel L ∈ ℒ_A con parámetro z; Q_c = Q_c(L)
(doc 110); banda = {r visible, r > 1: u_r < 3/2, Π(r) < 1/2} con
u_r = log(L/r)/log z; capa alta = banda con u_r ≥ u₀*; parte no
degenerada = V̂⁽ʳ⁾ ≠ {1} (la parte degenerada es contabilidad de
(SLAB-V′), doc 120a §1, y queda fuera de esta sección).

## 1. Enunciado exacto de (BAND-δ-HI′) con los objetos del flujo

Objetos del doc 114, por conductor r: V⁽ʳ⁾ y V̂⁽ʳ⁾ los conjuntos
fusionado y modelo sobre el nivel; **Add(r)** = V⁽ʳ⁾ ∖ V̂⁽ʳ⁾,
**Rem(r)** = V̂⁽ʳ⁾ ∖ V⁽ʳ⁾; a(d′) = μ(d′)/φ(d′);
E(r) = S(r) − Ŝ(r) = Σ_{Add} a − Σ_{Rem} a [T], con
E_Add(r) := Σ_{d′∈Add(r)} a(d′) y E_Rem(r) := −Σ_{d′∈Rem(r)} a(d′).
El apareamiento φ: Add⁺ → Add⁻ del doc 114 (Add⁻ = {y ∈ Add: μ(y) =
−1}, Add⁺ = {y ∈ Add: μ(y) = +1, y > 1}) da la **identidad de flujo**
[T, doc 114 §3]:

$$E_{Add}(r) \;=\; -\sum_{y \in Add^-(r)} |a(y)|\,\bigl(1 - R_A(y)\bigr),
\qquad R_A(y) = \sum_{q\,:\,yq \in Add^+} \frac{1}{q-1}.$$

**Déficit Add de r** := 𝒟_A(r) := −E_Add(r) = Σ_{y∈Add⁻(r)} |a(y)|
(1 − R_A(y)) (suma con signo [T]); es la masa que el apareamiento
*deja de ganar* en r: cada socio y ∈ Add⁻ retiene la fracción
(1 − R_A(y)) de su masa tras absorber a sus preimágenes de Add⁺.
Con GOOD(r) = Σ |a(y)|(1−R_A(y))⁺ y BAD(r) = Σ |a(y)|(R_A(y)−1)⁺
(doc 117a), 𝒟_A(r) = GOOD(r) − BAD(r); si BAD(r) = 0 y (DESP-SIGN)
vale en r, 𝒟_A(r) = |E_Add(r)| ≥ 0.

**(BAND-δ-HI′) [P] — forma exacta.** Para todo L ∈ ℒ_A:

$$\sum_{\substack{r \in \mathrm{banda}(L),\; u_r \ge u_0^*\\
\hat V^{(r)} \ne \{1\}}}
\frac{\bigl((\hat S-\Pi)(r)^{+} - \mathcal D_A(r)\bigr)^{+}}{\varphi(r)}
\;\le\; C\cdot Q_c(L),$$

con u₀* el umbral de la FL del modelo (u₀* < 3/2; en ventana se mide
u₀* ∈ {1,0, 1,2, 1,4}), declarado módulo las dos provisiones del
doc 120a §2, ambas [M] en ventana y ninguna [T]: (i) E_Rem(r) = 0
sobre los conductores de banda de la capa; (ii) E_Add(r) ≤ 0 ahí.
Bajo (i)+(ii), 𝒟_A(r) = |E_Add(r)|, el término interno es
((Ŝ−Π)⁺ − |E_Add|)⁺ = δ(r)⁺ [T, álgebra puntual, doc 120a §2], y
(BAND-δ-HI′) es exactamente la parte no degenerada de (BAND-δ-HI).
Sin las provisiones, la forma con 𝒟_A con signo es la que habría
que demostrar; la medición de §2 verifica (i), (ii) y BAD = 0 en
toda la capa medida, así que en ventana las tres formas coinciden.

## 2. Medición ampliada de la desigualdad sumada [M]

Script `122a_banda_delta_hi_prima.py` (PASS, K0–K4, 1,3 s).
Marcos A (z = 17, 29, 37), B (z = 29) y estrés (z = 53), los
niveles del doc 120a. K0 verifica sobre la capa, a precisión de
máquina: las identidades δ = (Ŝ−Π) + E (error 0,0) y
E = E_Add + E_Rem (error ≤ 2,2·10⁻¹⁶); la identidad de flujo
𝒟_A(r) = Σ_{y∈Add⁻} |a(y)|(1−R_A(y)) contra −E_Add(r) (error
0,0 exacto en las 67 instancias pooled — la razón aparece abajo);
BAD(r) = 0 y R_A(y) ≤ 1 en todo lo medido; y las provisiones (i)
E_Rem = 0 y (ii) E_Add ≤ 0 sobre la capa (0 violaciones en las 67
instancias no degeneradas de capa, que son todas las que hay:
toda instancia de capa medida tiene exceso > 0).

**El apareamiento está vacío sobre la capa [M].** En las 67
instancias, Add⁺ tiene **0 elementos en total** y Add⁻ tiene 41:
todo elemento de Add(r) en la capa tiene μ = −1 y ninguna
preimagen yq ∈ Add⁺, así que R_A(y) ≡ 0 y el déficit es la masa
cruda de Add, 𝒟_A(r) = Σ_{y∈Add(r)} |a(y)| = |E_Add(r)|, sin
absorción alguna. La maquinaria de flujo del doc 114 no llega a
actuar sobre la capa en ventana: no hay nada que aparear.

Por nivel y marco se separan tres sumas sobre la capa no
degenerada (u₀* = 1): exceso X = Σ (Ŝ−Π)⁺/φ, déficit
D = Σ 𝒟_A/φ y no cubierto U = Σ ((Ŝ−Π)⁺ − 𝒟_A)⁺/φ = Σ δ⁺/φ,
junto con la masa M de §3. En el nivel de máximo U de cada marco
(unidades de Q_c; z = 17 no tiene instancias de capa):

| marco | X | D | U | M | D/X |
|---|---|---|---|---|---|
| A, z=29 | 0,0530 | 0,0332 | 0,0198 | 0,0533 | 0,627 |
| A, z=37 | 0,0933 | 0,0000 | 0,0933 | 0,1567 | 0,000 |
| B, z=29 | 0,0037 | 0,0000 | 0,0037 | 0,0059 | 0,000 |
| estrés, z=53 | 0,2695 | 0,1828 | 0,0867 | 0,3152 | 0,678 |

**El déficit no cubre al exceso ni siquiera en suma.** El
cociente D/X es < 1 en *todo* nivel con X > 0 (rangos por nivel:
z=29 A: el único nivel con masa, 0,627; z=37: [0,000, 0,878];
z=29 B: 0,000; z=53: [0,678, 0,747]) y en los niveles donde U es
máximo vale 0,627, **0,000**, **0,000** y 0,678: en dos de los
cuatro marcos el nivel peor no recibe ni un átomo de déficit.
(BAND-δ-HI′) vale en ventana no porque el flujo gane, sino
porque U — lo que el flujo no cubre — es él mismo
≤ 0,0933·Q_c: la misma pequeñez de ventana de los docs 118–120a.
La tendencia con z es de **erosión**, no de margen creciente:
X/Q_c en el nivel de máximo U crece 0,0530 → 0,0933 → 0,2695
(z = 29 → 37 → 53), el máximo de X sobre niveles crece
0,0530 → 0,1570 → 0,2695, el de M crece 0,0533 → 0,1956 →
0,3152, y D/X no mejora (en z=37 el nivel peor tiene D = 0).
Con u₀*: en z = 29, 37 y B el nivel de máximo U es el mismo para
u₀* = 1,0, 1,2 y 1,4 — la masa ya vive en u_r ≥ 1,4, como vio el
doc 120a — y en z=53 subir a u₀* = 1,4 baja U a 0,0491 pero deja
X = 0,1209 y *empeora* D/X a 0,594 (rango [0,594, 0,732] contra
[0,678, 0,747]): acercarse al borde u = 3/2 no recupera
cobertura, la pierde.

**El subtotal E = 0.** De las 67 instancias de capa con exceso
positivo, 40 tienen E(r) = 0 (39 con E == 0,0 exacto y una más
con |E| ≤ 10⁻¹²); su déficit es 𝒟_A(r) = 0, así que su exceso
pasa íntegro a U: cargan el **58,7 %** de la masa no cubierta
pooled (0,9020 de 1,5360, en unidades Σ/φ/Q_c sumadas sobre
niveles). Qué pasa ahí exactamente, medido: **las 40 tienen
|Add| = |Rem| = 0** — ninguna tiene flujo no vacío que se cancele.
En ellas V⁽ʳ⁾ = V̂⁽ʳ⁾ como conjuntos: el desplazamiento no existe
y el exceso es el defecto intrínseco S(r) − Π(r) del conjunto
verdadero. Recíprocamente, E(r) ≠ 0 ⟺ Add(r) ≠ ∅ en ventana
(27 instancias, 41 elementos, todos μ = −1). En el primer caso
(BAND-δ-HI′) no pide nada al flujo en esos conductores: pide
que el defecto intrínseco de los conductores donde modelo y
realidad coinciden sea sumable — territorio de la banda misma,
no del apareamiento.

## 3. El lema que falta: ni puntual ni de déficit solo; la identidad de los ausentes

**Candidatos puntuales: falsos, maximalmente.** (P1) (Ŝ−Π)⁺(r) ≤
𝒟_A(r) falla en *todas* las instancias medidas con exceso positivo
y E_Rem = 0 — 67 de 67: en ellas δ(r)⁺ > 0 es exactamente la
diferencia. (P2) (Ŝ−Π)⁺(r) ≤ κ·𝒟_A(r) falla para todo κ en las
40 instancias con 𝒟_A(r) = 0 y exceso > 0: ninguna cota del
exceso por una función del déficit en el mismo conductor puede
ser el lema. **La dominación en suma tampoco es déficit contra
exceso**: §2 muestra D < X en todos los niveles con masa. Lo que
sostiene la desigualdad en ventana no es compensación entre
conductores vía flujo — es que U es pequeño porque *ambas* masas
viven en la escala invisible de la ventana.

**La identidad que reubica el problema [T].** Sea Out(r) el
complemento de V⁽ʳ⁾ ∪ V̂⁽ʳ⁾ en el nivel (los d′ coprimos con r que
no están ni en el conjunto fusionado ni en el modelo), y sean
M(r) = Σ_{d′∈Out(r), μ(d′)=−1} |a(d′)|, P(r) = Σ_{d′∈Out(r),
μ(d′)=+1} |a(d′)| sus masas por signo. Como Π(r) es la suma completa
sobre el nivel y S(r) la suma sobre V⁽ʳ⁾, particionar el
complemento de V⁽ʳ⁾ en Rem(r) ⊎ Out(r) da, término a término:

$$\delta(r) \;=\; E_{Rem}(r) \;+\; M(r) \;-\; P(r) \qquad [T].$$

(K0 la verifica a precisión de máquina sobre la capa.) Bajo la
provisión (i), δ(r) = M(r) − P(r) en la capa y δ(r)⁺ ≤ M(r)
puntualmente. El exceso no cubierto no es, en el fondo, un déficit
de flujo: es la masa μ = −1 de los elementos *ausentes de ambos
mundos*, menos su contraparte μ = +1. El lema mínimo suficiente
queda así desacoplado del apareamiento:

**(OUT-MASS) [P]:** Σ_{r ∈ capa(L)} M(r)/φ(r) ≤ C·Q_c(L) para
todo L ∈ ℒ_A. (OUT-MASS) + (i) ⟹ (BAND-δ-HI′) [T], y es
estrictamente más fuerte que (BAND-δ-HI′) en general: descarta la
cancelación M ↔ P dentro de Out, que la identidad deja abierta.

**¿Vale (OUT-MASS) en ventana? Sí, con C ≈ 0,32 — y perdiendo
un factor ~2,3 por el camino.** El máximo sobre niveles de
Σ_{capa} M/φ es 0,0533 (z=29 A), 0,1956 (z=37), 0,0059 (z=29 B)
y 0,3152 (z=53, estrés), siempre en unidades de Q_c: (OUT-MASS)
vale en ventana con C = 0,32, aunque con la misma tendencia
creciente con z que todo lo demás (0,053 → 0,196 → 0,315).
Pooled sobre niveles, Σ M/φ/Q_c = 3,4733 contra
Σ δ⁺/φ/Q_c = 1,5360: M es 2,26 veces el no cubierto, es decir
P cancela el 55,8 % de M dentro de Out. La compensación real es
mixta y queda medida: la desigualdad sumada se sostiene porque
M es ya pequeño en ventana (≤ 0,32·Q_c por nivel) *y* porque P
se come algo más de la mitad de M; (OUT-MASS) descarta esa
cancelación y aun así pasa en ventana — es un lema suficiente
honesto, no ajustado. Y δ(r)⁺ ≤ M(r) se verifica sin
violaciones en las 67 instancias (todas con E_Rem = 0) [M].
El doble conteo que (OUT-MASS) pediría controlar es explícito:
Σ_r M(r)/φ(r) = Σ_{d′: μ(d′)=−1} |a(d′)| · Σ_{r ∈ capa, (d′,r)=1,
d′∈Out(r)} 1/φ(r) — la multiplicidad de cada ausente μ = −1 sobre
la capa, pesada por 1/φ(r). Ningún objeto de (FC-2) aparece aquí.

## 4. Conexión con (FC-2): independientes también en la estructura de decaimiento

(FC-2) (doc 117a §4) es un enunciado **por conductor y dentro de
Add⁻**: BAD(r) = Σ_{y∈Add⁻} |a(y)|(R_A(y)−1)⁺ ≤ (1 − Mertens(37))·
HEAD37(r) para todo r visible. Su decaimiento es el de la *razón de
apareamiento*: malo ⟹ min(y) ≥ 73 [T], cada término malo decae
como (Mertens(min(y)−1) − 1)·|a(y)| y el total de la familia prima
libre decae como (log log L)/log L *con PNT* en la ventana cúbica
[(L/r₀)^{1/3}, L^{1/3}) [T, doc 117a §1b]. Sus dos bloqueadores
son ese decaimiento PNT y una cota inferior uniforme de HEAD37(r).

(BAND-δ-HI′), tras §3, es un enunciado **sumado sobre una capa
delgada con pesos 1/φ(r)** sobre la masa de Out — elementos que no
están en Add ni en Rem. Los dos enunciados no comparten objeto:
(FC-2) nunca sale de Add⁻(r); la masa de (BAND-δ-HI′) nunca entra
en Add. No comparten mecanismo de decaimiento: el de (FC-2) es
Mertens/PNT sobre R_A − 1 dentro del apareamiento; el que
(OUT-MASS) necesitaría es un conteo de ausentes dobles (fuera de
V y de V̂ a la vez) sobre la capa, para el que ni siquiera hay
todavía una familia obstructora identificada. Y no comparten
cuantificadores: ∀r puntual contra Σ sobre capa.

La ventana lo certifica de la manera más cruda [M]: en todos los
marcos medidos (z ≤ 53 < 73), BAD(r) ≡ 0 idénticamente — de
hecho es [T] en estos marcos, por el umbral del doc 117a §1a:
min(y) < z ≤ 53 ≤ 67 ⟹ R_A(y) ≤ Mertens(67) < 1 — así que
(FC-2) es *vacuo* en ventana, su masa mala no existe aún
(el script lo confirma: BAD máx = 0,0 en la capa, con
HEAD37(r) ≥ 0,1000 donde Add⁻ ≠ ∅), mientras la
masa no cubierta de (BAND-δ-HI′) es estrictamente positiva
(U/Q_c hasta 0,0933 en el nivel peor). Un enunciado vacío y uno con contenido no
pueden estar compartiendo el lema que los haría ciertos. Lo único
común es el ambiente: ambas masas viven de |a(d′)| = 1/φ(d′) y de
sumas tipo Mertens que crecen como log log — la razón genérica por
la que todo en esta línea es pequeño en ventana y nada decae
visiblemente. **Veredicto de la §4: independencia también a nivel
de estructura de decaimiento** — cerrar (FC-2) no acercaría
(BAND-δ-HI′), y cerrar (OUT-MASS) no tocaría (FC-2). (Esto refina
el doc 120a §3: allí la independencia era del *margen* de
(FLOW-COMP); aquí es del *mecanismo*.)

## 5. Veredicto

**(BAND-δ-HI′): apoyada en ventana como desigualdad sumada;
reducida a un residuo menor y mejor nombrado.**

- Enunciado exacto fijado (§1) con el déficit 𝒟_A(r) del flujo
  del doc 114; en ventana las provisiones (i), (ii) y BAD = 0
  valen en toda la capa medida [M, K0], así que la forma medida
  es la forma enunciada.
- La lectura "el déficit domina al exceso" es **falsa en
  ventana** en su forma ingenua: D < X en todos los niveles de
  masa [M] (D/X ≤ 0,878 siempre; 0,000 en los niveles peores de
  z=37 y del marco B), el déficit no cubre ni en suma, y la
  tendencia con z y con u₀* es de erosión, no de margen. Sobre
  la capa, además, el apareamiento está vacío (Add⁺ = 0
  elementos, R_A ≡ 0): el "déficit" medido es masa cruda de
  Add, no flujo parcialmente cancelado. La desigualdad vale
  porque el no cubierto U es pequeño en ventana
  (≤ 0,0933·Q_c), no porque el flujo gane.
- Ningún lema puntual exceso ≤ f(déficit) existe: refutado en
  ventana con contraejemplos en *todas* las instancias con
  exceso (67/67 para f = identidad; 40 instancias con
  𝒟_A = 0 para todo f con f(0) = 0) [M, §3]. El subtotal
  E = 0 (58,7 % de U) no es un fallo del flujo: ahí el
  desplazamiento no actúa (|Add| = |Rem| = 0) y el exceso es
  defecto intrínseco del conjunto verdadero.
- Residuo mínimo que queda, vía la identidad δ = E_Rem + M − P
  [T]: **(OUT-MASS)** [P] — la masa μ = −1 de los ausentes de
  ambos mundos, sumada sobre la capa, ≤ C·Q_c — más la provisión
  (i) E_Rem = 0 en la capa. (OUT-MASS) **vale en ventana con
  C ≈ 0,32** [M] (máx Σ M/φ = 0,3152·Q_c, z=53), perdiendo un
  factor pooled ~2,26 frente a δ⁺ porque descarta la
  cancelación M↔P (P se come el 55,8 % de M). Es un enunciado
  de conteo de ausentes, sin apareamientos: estrictamente menor
  que (BAND-δ-HI′) en estructura, y el primero de la línea de
  la banda que no pasa por E en absoluto.
- (FC-2) y (BAND-δ-HI′)/(OUT-MASS): **independientes** también
  en mecanismo de decaimiento [T de estructura + M de ventana].

Inventario de (H2) tras esta sección: (FC-2), **(OUT-MASS)**
(como forma mínima de (BAND-δ-HI′), módulo E_Rem = 0 en la capa),
(BAND-LOW-NONDEG), (SLAB-V′). **Goldbach permanece abierto.**

## 6. Reproducibilidad

Script `122a_banda_delta_hi_prima.py` (solo stdlib; importa la
maquinaria de los docs 110/114/117b/120a sin tocarla), salida
`salida_doc122a_local.txt`.
