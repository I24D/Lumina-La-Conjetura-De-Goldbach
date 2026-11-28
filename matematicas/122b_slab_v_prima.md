# Doc 122b — (SLAB-V′): enunciado exacto, medición extendida (z ≤ 53) e intento de cierre

**Fecha:** 2026-11-28. **Línea:** UNIF-U / (H2), sección del doc 122
sobre el residuo (SLAB-V′) del doc 119b §4.
**Script:** `122b_slab_v_prima.py` (stdlib; importa la maquinaria de
los docs 110/117b/118b1 sin tocarla). **Salida:**
`salida_doc122b_local.txt` — chequeos K0–K3.
**Etiquetas:** [T] = teorema, [M] = medición en ventana, [P] = pendiente.

---

## 0. Veredicto en breve

**(SLAB-V′) queda ABIERTA [P], con el obstáculo exacto
identificado y la medición extendida a su favor.** El único
telescopio de la línea vive en la dirección $q_0$ (doc 118b1
§3) y está agotado en el corte: para $q_0$ fijo no hay cola
que telescopar, y no existe factorización en el cofactor $c$
ni promedio en $L$ que lo sustituya [T la obstrucción, §3] —
el cierre es un teorema de Mertens sobre $V^{(c)}$ con la
degeneración del cociente conservada, misma especie que
(FL-DESP). La medición de la familia completa ($z = 17$–$53$,
K0–K3 PASS) la apoya sin excepción: $P(5)/Q_c \le 0{,}79$ con
crecimiento en $z$ desacelerado ($0{,}598 \to 0{,}711 \to
0{,}783$), el estrato 7 se enciende a $z=53$ ($P(7) =
1{,}76\,Q_c$), el 11 sigue en 0, la masa vive solo en niveles
pequeños, y los únicos cofactores que ocurren son $c \in
\{1, q_0, q_0^+\}$ con $N(r) \le 3$ — rigidez que queda como
sub-lema candidato (WIT) [P].

## 1. Enunciado exacto de (SLAB-V′) [transcripción del doc 119b, §§1 y 4]

**Marco.** $z$ fijo, $\mathcal P = \{p \text{ primo} : 5 \le p < z\}$
(el conjunto `allow` de la serie), $A = 3$, $\mathcal L_A$ los
niveles $L = \lfloor D_A/s \rfloor$ de la serie. $V^+(L)$ es el
conjunto de $d$ libres de cuadrados con primos en $\mathcal P$,
$d < L$, que cumplen las condiciones de Nathanson de la serie
(el `lv.V` de la maquinaria). $Q_c = \Pi_{\mathrm{full}} \cdot
M_{\mathrm{full}}$ con $\Pi_{\mathrm{full}} = \prod_{p \in
\mathcal P} (1 - \frac{1}{p-1})$ la densidad de $V^+(L)$.

**Los objetos.** $r$ es *visible* si $r \mid d$ para algún
$d \in V^+(L)$; está *en banda* si $u_r := \log(L/r)/\log z < 3/2$
y $\Pi(r) < 1/2$; es *degenerado* si $\widehat V^{(r)} = \{1\}$,
equivalentemente [T, doc 117b §3.2] $q_0(r)^3 \cdot r \ge L$ con
$q_0(r) := \min(\mathcal P \setminus \{\text{primos de } r\})$
— el *estrato* de $r$. Todo $r$ degenerado de estrato $q_0$
contiene todo primo de $\mathcal P$ menor que $q_0$ y $q_0 \nmid r$.
La familia finita es $q_0 \in \{5, 7, 11\}$ (los estratos
$q_0 \ge q_0^\dagger$, $q_0^\dagger = \min\{q : P(q) \ge
(\log z)^3\} = 13$ en ventana, están cerrados [T] por el
telescopio de Euler del doc 118b1 §3).

**Las sumas.** Para $q_0 \in \{5, 7, 11\}$ y cada cofactor $c$,
$$ U_c^{\mathrm{deg}}(q_0) := \sum_{\substack{r' \text{ deg., visible, en banda} \\
\mathrm{estrato}(r') = q_0,\ c r' \in V^+(L)}} \frac{1}{\varphi(r')}, $$
la restricción a cocientes degenerados de la suma de ventana
$U_c$ del doc 119b §1 sobre el soporte desplazado $V^{(c)}$.
Condiciones que el par $(c, r')$ cumple automáticamente [T,
119b §1]: $c < q_0^3$ (de $r' \ge L/q_0^3$ y $cr' < L$),
$(c, r') = 1$ y $\varphi(cr') = \varphi(c)\varphi(r')$ (todo es
libre de cuadrados), $(c, P(q_0)) = 1$ con $P(q_0) =
\prod_{p \in \mathcal P, p < q_0} p$ (los primos $< q_0$ están
todos en $r'$), y $\omega(c) \le 3$ ($5 \cdot 7 \cdot 11 \cdot
13 = 5005 > 1331 \ge q_0^3$ acota los primos $\ge 5$ de $c$).
La familia es finita **en formas**: $q_0 \in \{5,7,11\}$,
$c < q_0^3 \le 1331$ libre de cuadrados con $\omega(c) \le 3$;
los valores primos concretos de $c$ recorren $\mathcal P$.

**(SLAB-V′)** [P, doc 119b §4]: para cada $q_0 \in \{5, 7, 11\}$,
$$ P(q_0) := \sum_{c < q_0^3} U_c^{\mathrm{deg}}(q_0)
\;\le\; C(q_0)\, Q_c \quad \text{uniformemente en } L \in
\mathcal L_A. $$

**Forma de pares equivalente** [T, 119b §1]:
$P(q_0) = \sum_{r \text{ deg., vis., estrato } q_0}
\frac{N(r)}{\varphi(r)}$ con $N(r) := \#\{c < q_0^3 : (c, r) = 1,\
cr \in V^+(L)\} \ge 1$. La única pérdida de la reindexación
respecto de (SLAB-V) es el factor $N(r)$; la ruta que suelta la
degeneración (mayorante $W(q_0)$ de 119b) está refutada como
ruta [M, 119b §2: pierde $23$–$217\times$]. El factor
$1/\ln z$ que el enunciado necesita es $\Pi_{\mathrm{full}}$,
ya presente en cada $U_c$ [M, 119b §2: predicción
$\Pi_{\mathrm{full}} A(c) \ln(q_0^3)$ dentro de factor $2$–$3$].

## 2. Medición [M]: la familia completa, $z = 17$–$53$

El script extiende el K3 del doc 119b: marcos A ($z = 17, 29,
37, 53$; $D = z^5$, niveles $\lfloor D/s \rfloor$) y B ($z = 29$;
$D = z^6$). Por nivel calcula $U_c^{\mathrm{deg}}(q_0)$ para
todo $c$ que ocurre, $P(q_0)$, la suma verdadera, $N(r)$, el
$\omega$ mínimo de los $r$ degenerados y la rebanada
$\sum_{d \in V^+, d \ge L/125} 1/\varphi(d)$; chequeos K0
(identidad $\varphi(c)\varphi(r) = \varphi(d)$, cobertura
$N(r) \ge 1$), K1 (doble contabilidad de $P$: por cofactores y
por $r$ con multiplicidad), K3 (Teorema A: si $L \ge q_0^9$,
ningún degenerado visible del estrato $q_0$ con $\omega \le 2$).

**Resultados** (máximos sobre niveles, en unidades de $Q_c$;
K0/K1 exactos en todos los marcos: 0 violaciones de
$\varphi(c)\varphi(r) = \varphi(d)$, 0 sin cobertura, doble
contabilidad con error $\le 5{,}6 \cdot 10^{-17}$):

| marco | $q_0$ | verdadera | $P(q_0)$ | en $L$ | $N$ máx | $U_c^{\mathrm{deg}}$ que ocurren |
|---|---|---|---|---|---|---|
| A $z=17$ | 5, 7, 11 | 0 | 0 | — | — | ninguno |
| A $z=29$ | 5 | 0,1993 | 0,5979 | 1360 | 3 | $c=1,5,7$: 0,1993 c/u |
| A $z=37$ | 5 | 0,3554 | 0,7109 | 410 | 3 | $c=1$: 0,3554; $c=5$: 0,3554; $c=7$: 0,2133 |
| A $z=53$ | 5 | 0,3917 | 0,7834 | 747 | 3 | $c=1$: 0,3917; $c=5$: 0,3917; $c=7$: 0,2350 |
| A $z=53$ | 7 | 0,5875 | 1,7626 | 1364 | 3 | $c=1,7,11$: 0,5875 c/u |
| B $z=29$ | 5 | 0,0252 | 0,0252 | 24399 | 1 | $c=1$: 0,0252 |
| todos | 11 | 0 | 0 | — | — | ninguno |

Hechos que la medición establece:

1. **Concentración extrema en el cofactor** [M]: en todos los
   marcos, los únicos $c$ con $U_c^{\mathrm{deg}} > 0$ son
   $c \in \{1, q_0, q_0^+\}$, con $q_0^+$ el siguiente primo de
   $\mathcal P$ tras $q_0$ (estrato 5: $\{1,5,7\}$; estrato 7:
   $\{1,7,11\}$). Ningún $c \ge 13$ aporta jamás, pese a que la
   forma lo permite ($c < 1331$, $\omega(c) \le 3$). $N(r) \le 3$
   en todo lo medido, y en los máximos $P/\text{verdadera} \in
   \{1, 2, 3\}$ exacto: cada $r$ degenerado tiene 1, 2 o 3
   testigos y los tres estratos de testigos aparecen completos.
2. **La masa vive en niveles pequeños** [M]: el máximo de
   $P(q_0)$ ocurre siempre en $L \in [410, 24399]$ (el extremo
   inferior de $\mathcal L_A$ en ventana) y $P \to 0$ en el
   nivel mayor de cada marco (p. ej. $z=53$: $P(5)/Q_c =
   0{,}0020$ en $L = 7{,}9 \cdot 10^6$). La degeneración con
   visibilidad es un fenómeno de $L$ pequeño en ventana.
3. **Crecimiento en $z$ desacelerado** [M]: $P(5)/Q_c =
   0{,}598 \to 0{,}711 \to 0{,}783$ para $z = 29 \to 37 \to 53$.
   Si creciera como $\log z$, el segundo incremento debería ser
   $\approx 0{,}166$ (escalando el primero por $\Delta \log z$);
   el medido es $0{,}072$, menos de la mitad. Compatible con
   saturación o crecimiento $\log \log$; incompatible con la
   pérdida $(\log z)^2$ de las rutas mayorantes.
4. **El estrato 7 se enciende a $z=53$** [M]: era 0 a $z \le 37$
   y aparece con $P(7) = 1{,}76\,Q_c$ ($L = 1364$). Su
   ausencia en ventana pequeña era artefacto, no teorema — la
   uniformidad de (SLAB-V′) no puede apoyarse en ella. El
   estrato 11 sigue en 0 en toda la ventana.
5. **Teorema A reverificado** (K3): 0 violaciones bajo
   hipótesis ($L \ge q_0^9$, niveles que existen a $z=53$);
   los degenerados visibles con $\omega \le 2$ solo aparecen
   con $L < q_0^9$ (conteos por marco: 0/10/78/100/185 —
   artefactos de ventana, como manda el doc 118b1).
6. **La rebanada sin degeneración no controla** [M]: la suma
   de rebanada $\sum_{d \in V^+, d \ge L/q_0^3} 1/\varphi(d)$
   vale $0{,}10$–$3{,}43\,Q_c$ (máximos por marco/estrato) —
   es $\Theta(Q_c)$ en ventana pero con constante hasta
   $4\times$ la de $P$ y $> 0$ incluso donde $P = 0$ ($z=17$:
   $1{,}70$–$2{,}48\,Q_c$): la pequeñez de $P$ es trabajo de
   la conjunción rebanada × visibilidad × degeneración, no de
   la rebanada sola.


## 3. Intento de cierre: no hay telescopio para $q_0$ fijo [T la obstrucción]

**Por qué el caso $q_0 \ge 13$ era especial.** El cierre del
doc 118b1 §3 es un telescopio **en la dirección $q_0$**: la
cota del estrato es $E(z)\,\frac{q_0 - 1}{q_0 P(q_0)}$, y su
pequeñez en el corte viene entera de $1/P(q_0^\dagger) \le
(\log z)^{-3}$ — el producto $P(q_0) = \prod_{p < q_0} p$
*crece* con $q_0$ y estrangula la cola. Para $q_0 \in
\{5, 7, 11\}$ fijo, $P(q_0) \in \{1, 5, 35\}$ es una constante:
no hay cola, no hay dirección en la que telescopar, y la misma
cuenta da $\Theta(E(z)) = \Theta(\log z)$ contra $Q_c =
\Theta(\Pi_{\mathrm{full}}) = \Theta(1/\log z)$ — la pérdida
$(\log z)^2$ del doc 118b1 §4, que es un teorema sobre la ruta:
todo mayorante que factorice por el producto de Euler sin
ventana pierde $(\log z)^2$ en los estratos finitos.

**No hay telescopio en $c$.** Un telescopio en el cofactor
exigiría que $\sum_c U_c^{\mathrm{deg}}$ factorice como
producto sobre los primos de $c$. No factoriza: $U_c^{\mathrm
{deg}}$ acopla $c$ y $r'$ por dos condiciones no
multiplicativas — $cr' \in V^+(L)$ (las condiciones de
Nathanson se evalúan sobre la lista ordenada *fusionada* de
primos de $cr'$) y la degeneración de $r'$ con su testigo.
La medición lo confirma desde el otro lado: si hubiera
estructura de producto, ocurrirían $c$ con $\omega(c) = 2, 3$
($35, 55, 5 \cdot 7 \cdot 11, \dots$) con masa comparable;
ocurren solo $c \in \{1, q_0, q_0^+\}$, primos o 1 (§2.1).

**No hay telescopio en $L$.** (SLAB-V′) es uniforme en
$L \in \mathcal L_A$, no una suma sobre niveles: no se puede
promediar en $L$ lo que debe valer nivel a nivel. (Y la
medición §2.2 muestra que el peor nivel es el pequeño, justo
donde un promedio escondería la masa.)

**El factor $1/\varphi(P(q_0))$ no cierra.** El estrato fuerza
$P(q_0) \mid r$, luego $1/\varphi(r) = 1/[\varphi(P(q_0))\,
\varphi(r/P(q_0))]$: ganancias constantes $1, 1/4, 1/24$ para
$q_0 = 5, 7, 11$ frente a una suma de comparación sin el
primo forzado. Pero la comparación no es una dominación: la
visibilidad no desciende por división ($r$ visible no implica
$r/P(q_0)$ visible, el testigo no se divide entre $P(q_0)$ en
general), y una constante $1/24$ no rescata una ruta que
pierde $(\log z)^2$. El factor explica la jerarquía entre
estratos en una contabilidad de Euler; no explica — ni
sustituye — la pequeñez medida.

**Qué haría falta.** La pequeñez de $P(q_0)$ proviene de la
conjunción de tres condiciones — rebanada de razón $q_0^3$,
visibilidad, degeneración — cuya región límite tiene medida
de Mertens cero ($u_r \le 3 \ln q_0 / \ln z \to 0$, doc 117b
§3.2). Cada mayorante que suelta una de las tres está refutado
como ruta con pérdida exacta: ventana+visibilidad [118b1 §4,
$(\log z)^2$], degeneración del cociente [119b §2,
$23$–$217\times$], rebanada sola [§2.6 aquí]. Un cierre de
(SLAB-V′) es, por construcción, un teorema de Mertens sobre
los soportes desplazados $V^{(c)}$ *con* la degeneración del
cociente: la misma especie de objeto que (FL-DESP) en los
docs 112–113, no una identidad de Euler. No se afirma
implicación entre ambos; se afirma que ningún atajo
telescópico existe en esta contabilidad.

## 4. Residuo mínimo y veredicto

El residuo es el propio (SLAB-V′) del §1, con dos precisiones
que la medición añade sin cambiar el enunciado:

- **Forma reducida en ventana** [M]: basta acotar, por
  estrato, las sumas $U_c^{\mathrm{deg}}$ con $c \in
  \{1, q_0, q_0^+\}$ — son las únicas que ocurren a $z \le 53$
  — con constantes medidas $C(5) \approx 0{,}79$,
  $C(7) \approx 1{,}77$ (recién encendido), $C(11) = 0$ en
  ventana. Que el conjunto de testigos sea siempre
  $\{1, q_0, q_0^+\}$ es **[P] como sub-lema** (WIT): si se
  demostrara, (SLAB-V′) quedaría reducida a $\le 3$ sumas por
  estrato, finita en el sentido más fuerte; seguiría sin
  estar cerrada, porque cada una de esas sumas es una suma de
  Mertens sobre $V^{(c)}$ con cociente degenerado.
- **Dónde mirar**: el peor caso en ventana es $L$ pequeño
  (extremo inferior de $\mathcal L_A$), no el régimen puente;
  cualquier intento de cierre debe ser uniforme ahí o
  explicar por qué el régimen de la cadena excluye esos
  niveles.

**VEREDICTO: (SLAB-V′) queda ABIERTA [P]**, con el obstáculo
exacto identificado: (i) el único telescopio de la línea vive
en la dirección $q_0$ y está agotado en el corte $q_0^\dagger$
[T]; (ii) no existe factorización en $c$ ni promedio en $L$
disponible [T la obstrucción, §3]; (iii) su cierre es un
teorema de Mertens con degeneración de cociente sobre
$V^{(c)}$, misma especie que (FL-DESP). La medición extendida
($z \le 53$, familia completa) la apoya sin excepción:
constantes $\le 1{,}77\,Q_c$, crecimiento en $z$
desacelerado, concentración total en $c \in \{1, q_0,
q_0^+\}$. (BAND-MOD) sigue heredando (SLAB-V′) por su subcapa
degenerada. Goldbach no se declara demostrado.

## 5. Reproducibilidad

`python3 122b_slab_v_prima.py | tee salida_doc122b_local.txt`
— stdlib; importa `118b1_degenerada.py` (que encadena
117b→110). K0: identidad $\varphi(c)\varphi(r) = \varphi(d)$
en todo par y cobertura $N(r) \ge 1$. K1: doble contabilidad
de $P(q_0)$ (por cofactores y por $r$ con multiplicidad) y
verdadera $\le P$. K2: tablas por $(q_0, c)$, concentración,
tendencias en $z$ y en $L$. K3: Teorema A bajo hipótesis
$L \ge q_0^9$. **PASS global** (5,4 s).
