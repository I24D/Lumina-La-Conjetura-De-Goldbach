# Doc 119b — (SLAB-V): la rebanada superior de V⁺(L) por reindexación divisor–cofactor

**Fecha:** 2026-11-28. **Línea:** UNIF-U / (H2), sección del doc 119
sobre el residuo (SLAB-V) del doc 118b1 §5.
**Script:** `119b_slab_v.py` (stdlib; importa la maquinaria de los
docs 110/117b/118b1 sin tocarla). **Salida:**
`salida_doc119b_local.txt` — chequeos K0–K3.
**Etiquetas:** [T] = teorema, [M] = medición en ventana, [P] = pendiente.

---

## 0. Veredicto en breve

La reindexación divisor–cofactor pedida por el doc 118b1 §5 se
ejecuta completa y es **exacta** [T, §1]: todo $r$ degenerado
visible del estrato $q_0$ tiene cofactor $c = d/r < q_0^3$ en
cualquiera de sus testigos, $1/\varphi(r) = \varphi(c)/\varphi(d)$,
y la suma degenerada es la suma de pares
$P(q_0) = \sum_r N(r)/\varphi(r)$ con multiplicidad $N(r)$ —
medida $\le 3$ y concentrada en $c \in \{1, 5, 7\}$ [M, §3]. El
factor $1/\ln z$ que (SLAB-V) necesita queda **identificado**: es
$\Pi_{\mathrm{full}}$, la densidad de $V^+(L)$, presente en cada
suma de ventana $U_c$ de razón fija $q_0^3$ (la ventana aporta
solo $\ln(q_0^3) = O(1)$). Pero el mayorante cerrado $W(q_0)$ que
suelta la degeneración del cociente **no cierra** [M, §2]:
$W(5)/Q_c \approx 16$ con soporte de cofactores creciente en $z$,
pérdida $23$–$217\times$ localizada exactamente en los cocientes
no degenerados. **Veredicto: (SLAB-V) reducido [T] a (SLAB-V′)
[P, §4]** — la misma contabilidad con la degeneración del cociente
conservada, $\sum_c U_c^{\mathrm{deg}} \le C(q_0)\,Q_c$ — familia
finita de sumas sobre soportes desplazados $V^{(c)}$, sin
cuantificador de visibilidad. (BAND-MOD) hereda (SLAB-V′).

## 1. La reindexación divisor–cofactor [T]

Sea $q_0 \in \{5, 7, 11\}$ fijo, $B := q_0^3 \le 1331$, y $r$ un
degenerado visible del estrato $q_0$ en banda: $r \mid d$ para algún
$d \in V^+(L)$, $r < L$, $q_0^3 r \ge L$, $r$ contiene todo primo de
$\mathcal P$ menor que $q_0$ y $q_0 \nmid r$.

**Cofactor pequeño.** Para cualquier testigo $d$ de $r$:
$c := d/r < q_0^3 = B$, pues $d < L \le q_0^3 r$. Además $d \ge r \ge
L/B$: todo testigo vive en la **rebanada superior**
$\{d \in V^+(L) : L/B \le d < L\}$.

**Identidad de pesos.** Todo elemento del enum es libre de cuadrados,
luego $(c, r) = 1$ y $\varphi(d) = \varphi(c)\,\varphi(r)$:
$$ \frac{1}{\varphi(r)} = \frac{\varphi(c)}{\varphi(d)}. $$

**Compatibilidad de estrato.** Los primos de $\mathcal P$ menores que
$q_0$ están todos en $r$ y ninguno en $c$: $(c, P(q_0)) = 1$ con
$P(q_0) = \prod_{p \in \mathcal P,\, p < q_0} p$. Agrupando los pares
por el estrato del cociente, defínase el peso mayorante
$$ \Phi^{(q_0)}(d) := \sum_{\substack{c \mid d,\ c < B \\ (c,\, P(q_0)) = 1 \\
\mathrm{estrato}(d/c) = q_0}} \varphi(c),
\qquad
W(q_0) := \sum_{\substack{d \in V^+(L) \\ d \ge L/B}} \frac{\Phi^{(q_0)}(d)}{\varphi(d)}, $$
donde el estrato de $d/c$ es su propio $q_0(d/c)$. Cada $r$
degenerado visible del estrato aporta su término
$1/\varphi(r) = \varphi(c)/\varphi(d)$ a $W(q_0)$ vía cualquiera de
sus testigos, luego
$$ \sum_{r\ \mathrm{deg,\ vis,\ estrato}\ q_0} \frac{1}{\varphi(r)}
\;\le\; W(q_0). $$

**Intercambio de suma.** Agrupando por el cofactor:
$$ W(q_0) = \sum_{\substack{c < B \\ (c,\, P(q_0)) = 1}} U_c,
\qquad
U_c := \sum_{\substack{r' :\, c r' \in V^+(L),\, c r' \ge L/B \\
\mathrm{estrato}(r') = q_0}} \frac{1}{\varphi(r')}, $$
sumas de ventana sobre los soportes desplazados $V^{(c)}$ del doc
109, con ventana multiplicativa de razón exactamente $B = q_0^3$
fija. La lista de cofactores es finita y $z$-libre: $c < 1331$ con
primos en $\mathcal P$ fuerza $\omega(c) \le 3$
($5 \cdot 7 \cdot 11 \cdot 13 = 15015 > 1331$).

**La suma de pares (multiplicidad).** Sin soltar la degeneración
del cociente, la reindexación es una identidad con multiplicidad:
$$ P(q_0) := \sum_{c < B} U_c^{\mathrm{deg}}
= \sum_{r\ \mathrm{deg,\ vis,\ estrato}\ q_0} \frac{N(r)}{\varphi(r)},
\qquad
N(r) := \#\{c < q_0^3 : (c, r) = 1,\ cr \in V^+(L)\}, $$
con $U_c^{\mathrm{deg}}$ la restricción de $U_c$ a cocientes $r'$
degenerados visibles en banda. Todo $r$ degenerado visible tiene
$N(r) \ge 1$ [T, su testigo], luego verdadera $\le P(q_0) \le
W(q_0)$, y la única pérdida de la reindexación es el factor $N(r)$.

## 2. Dónde pierde el mayorante ingenuo $W$: obstrucción exacta [T/M]

La tentación es cerrar (SLAB-V) acotando $W(q_0)$: cada $U_c$ es una
suma de ventana de razón fija $q_0^3$ sobre un soporte desplazado, y
su única densidad $z$-dependiente posible es $\Pi_{\mathrm{full}}$
— la predicción del modelo $\Pi_{\mathrm{full}} \cdot A(c) \cdot
\ln(q_0^3)$, con $A(c) = \prod_{p \in \mathcal P, p \nmid c}
(1 + \frac{1}{p(p-1)}) \prod_{p \in \mathcal P, p \mid c}
(1 - \frac{1}{p})$, reproduce los $U_c$ medidos dentro de un factor
$2$–$3$ (K2: p. ej. $c=5$, $q_0=5$: medido $2{,}33$–$2{,}49\,Q_c$
contra predicho $3{,}52\,Q_c$). **El factor $1/\ln z$ sí está en
cada $U_c$: es $\Pi_{\mathrm{full}}$, la densidad de $V^+(L)$.** El
problema no es cada término, es la suma:

1. **El soporte de cofactores crece con $z$** [M]: el número de $c$
   que ocurren en $W(5)$ es $10/17/19$ a $z = 17/29/37$, y en
   $W(11)$ es $4/16/29$. En el régimen puente todo $c < q_0^3$ con
   primos en $\mathcal P$ acaba ocurriendo, y $\sum_c A(c)$ sobre
   esa lista fija pero numerosa ($\sim$ decenas) ya da
   $W(5)/Q_c \approx 16$ en ventana — sin crecer como $(\log z)^2$,
   pero con constante inaceptable y sin mecanismo que la reduzca.
2. **La pérdida está localizada** [M]: $W(5)/P(5) = 27{,}2$
   ($z=29$ A), $22{,}6$ ($z=37$ A), $217$ (B). Todo el exceso de
   $W$ viene de cocientes $r' = d/c$ **no degenerados**: para los
   estratos $7$ y $11$, donde la suma verdadera es $0$ en todos los
   marcos, $W$ vale hasta $9{,}38\,Q_c$ y $2{,}04\,Q_c$ de puro
   artefacto; y a $z=17$ (banda vacía) $W(5) = 10{,}35\,Q_c$ sin un
   solo $r$ degenerado.

Conclusión de la sección: **soltar la degeneración del cociente es
la obstrucción exacta de la ruta $W$**. La contabilidad correcta
es la suma de pares $P(q_0)$ de §1, que la conserva.

## 3. Medición [M] (K0–K3)

Identidades exactas en todos los marcos (K0/K2): $\varphi(c)
\varphi(r) = \varphi(d)$ en el 100% de los pares (0 violaciones);
todo cofactor de un testigo de un $r$ degenerado es $< q_0^3$
(0 violaciones); cobertura completa (ningún $r$ degenerado visible
sin par); intercambio $W = \sum_c U_c$ con error $\le 3{,}6 \cdot
10^{-15}$. Multiplicidad medida: $N(r) \le 3$, media $1{,}20$
($z=29$ A), $1{,}41$ ($z=37$ A), $1{,}00$ (B).

Sumas por estrato (máximos sobre niveles, en unidades de $Q_c$):

| marco | $q_0$ | verdadera | pares $P$ | $W$ | $P/$verdadera |
|---|---|---|---|---|---|
| A $z=29$ | 5 | 0,1993 | 0,5979 | 16,26 | 3,00 |
| A $z=37$ | 5 | 0,3554 | 0,7109 | 16,09 | 2,00 |
| B $z=29$ | 5 | 0,0252 | 0,0252 | 5,46 | 1,00 |
| todos | 7, 11 | 0 | 0 | $\le 9{,}38$ / $\le 2{,}04$ | — |

**Cofactores que cargan la masa degenerada** (K3, máximos por
nivel): estrato 5, marco A: solo $c \in \{1, 5, 7\}$ — a $z=29$,
$0{,}1993\,Q_c$ cada uno (los tres testigos de cada $r$ son $r$,
$5r$, $7r \in V^+(L)$); a $z=37$: $c=1$ y $c=5$ con $0{,}3554$ y
$c=7$ con $0{,}2133$. Marco B: solo $c=1$ ($N(r) = 1$: el único
testigo es $d = r$). La parte $c=1$ de $W(5)$ es $0{,}332/0{,}355/
0{,}065\,Q_c$ (A $z=29$ / A $z=37$ / B): la rebanada misma,
$\sum_{d \in V^+ \cap \mathrm{reb}} 1/\varphi(d)$, ya es $O(Q_c)$
en ventana. Estratos $\ge 13$ y None: verdadera $0$ en ventana,
como en el doc 118b1.

## 4. Veredicto: (SLAB-V) reducido a (SLAB-V′) [T la reducción]

**(SLAB-V) no queda cerrado [T]**, pero la reindexación lo reduce
exactamente, y el residuo es estrictamente menor:

**(SLAB-V′):** para cada $q_0 \in \{5, 7, 11\}$,
$\;P(q_0) = \sum_{c < q_0^3} U_c^{\mathrm{deg}} \le C(q_0)\, Q_c$
uniformemente en $L \in \mathcal L_A$ — la contabilidad de Mertens
de las sumas de ventana degeneradas sobre los soportes desplazados
$V^{(c)}$, con la degeneración del cociente **conservada**.

En qué es menor que (SLAB-V), punto por punto:

- El cuantificador de visibilidad desaparece [T]: ya no se suma
  sobre $r$ con un testigo existencial; se suma sobre pares
  $(d, c)$ con $d \in V^+(L)$ explícito, y la multiplicidad $N(r)$
  es un conteo de cofactores $c < q_0^3$ con $(c, r) = 1$ y
  $cr \in V^+(L)$ — medido $\le 3$, y la masa degenerada se
  concentra en $c \in \{1, 5, 7\}$ (marco A) o $c = 1$ (marco B).
- El factor $1/\ln z$ queda identificado [T/M]: viene de
  $\Pi_{\mathrm{full}}$, la densidad de $V^+$, presente en cada
  $U_c$ (predicción del modelo dentro de factor $2$–$3$ de lo
  medido); la ventana aporta solo $\ln(q_0^3) = O(1)$ y la lista de
  cofactores es $z$-libre. Esto responde la pregunta de la sección:
  la prueba de (SLAB-V) no necesita producir ningún otro factor —
  necesita exactamente la cota tipo FL superior para las
  $U_c^{\mathrm{deg}}$, una familia finita de sumas sobre
  $V^{(c)}$ con cociente degenerado.
- La ruta alternativa (mayorante $W$ sin degeneración del
  cociente) queda **refutada como ruta** [M, §2]: pierde un factor
  $23$–$217$ en ventana y su soporte de cofactores crece con $z$.
  Es el análogo exacto, dentro de la subcapa, de la pérdida
  $(\log z)^2$ de la ruta Euler del doc 118b1 §4: en ambos casos el
  conteo que ignora la estructura fina (allí la ventana y la
  visibilidad; aquí la degeneración del cociente) no cierra,
  aunque el enunciado sea compatible con toda la medición
  ($P(5) \le 0{,}72\,Q_c$ en ventana).

(BAND-MOD) hereda (SLAB-V′) por su subcapa degenerada, en el lugar
de (SLAB-V). Goldbach no se declara demostrado.

## 5. Reproducibilidad

`python3 119b_slab_v.py | tee salida_doc119b_local.txt` — stdlib;
importa `118b1_degenerada.py` (que encadena 117b→110→109). K0:
identidad $\varphi(c)\varphi(r) = \varphi(d)$ en todo par y
cofactor de testigo $< q_0^3$. K1: cobertura ($N(r) \ge 1$ en todo
$r$ degenerado visible) y cadena verdadera $\le P \le W$ por
estrato contra $Q_c$. K2: intercambio $W = \sum_c U_c$ exacto y
predicción $\Pi_{\mathrm{full}} A(c) \ln(q_0^3)$. K3: desglose de
$P$ por cofactor y parte $c=1$. **PASS global** (0,5 s).
