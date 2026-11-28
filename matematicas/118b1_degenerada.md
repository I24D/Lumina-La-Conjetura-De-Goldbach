# Doc 118b1 — La subcapa degenerada de (BAND-LOW): cierre parcial y obstrucción exacta

**Fecha:** 2026-11-28. **Línea:** UNIF-U / (H2), micro-sección del
ataque a (BAND-LOW) del doc 117b §3.2.
**Script:** `118b1_degenerada.py` (stdlib; importa la maquinaria de los
docs 110/117b sin tocarla). **Salida:** `salida_doc118b1_local.txt` —
chequeos J0–J4, **todos PASS** (0,7 s).
**Etiquetas:** [T] = teorema, [M] = medición en ventana, [P] = pendiente.

---

## 0. Veredicto en breve

La subcapa degenerada es $\{r \text{ visible, en banda} :
\widehat V^{(r)} = \{1\}\} = \{r : q_0^3\, r \ge L\}$, con
$q_0 = \min(\mathcal P \setminus \mathrm{primos}(r))$ [T, doc 117b
§3.2]; ahí $(\widehat S - \Pi)^+ = 1 - \Pi(r) \le 1$, luego su
contribución a (BAND-LOW) es $\le \sum_{\mathrm{deg}} 1/\varphi(r)$.
Resultado de esta micro-sección:

1. **Estratos $q_0 \ge q_0^\dagger$: CERRADOS [T]** (§3) por Euler con
   telescopio, con constante de cola $\to 0$; $q_0^\dagger = 13$ en los
   marcos de la serie.
2. **Teorema A [T, §2]:** si $L \ge q_0^9$, ningún degenerado visible
   del estrato $q_0$ tiene $\omega(r) \le 2$ — en el régimen puente la
   obstrucción prima del doc 116c no vive en esta subcapa. La
   hipótesis es **exactamente aguda en ventana** [M]: las 188
   instancias medidas tienen $\omega \le 2$ y *todas* cumplen
   $L < q_0^9$.
3. **La obstrucción exacta a cerrar los estratos finitos por conteo
   elemental es la ruta, no el enunciado** (§4): la cota Euler sin
   ventana del estrato $q_0 = 5$ es $\Theta(\log z) =
   \Theta((\log z)^2)\,Q_c$ (medido: $3{,}25$–$3{,}72\,Q_c$ y
   creciendo), mientras la suma verdadera visible es $\le 0{,}36\,Q_c$
   [M]. Lo que falta es exactamente **(SLAB-V)** [P, §5]: la
   contabilidad de Mertens de la rebanada superior de $V^+(L)$.

## 1. Agrupación por $q_0$: la suma se factoriza [T]

Todo $r$ degenerado contiene **todos** los primos de $\mathcal P$
menores que $q_0$ (por definición de $q_0$) y $q_0 \nmid r$. Con
$P(q_0) = \prod_{p < q_0,\, p \in \mathcal P} p$ y $r = P(q_0)\,m$:
$(m, P(q_0)) = 1$, $q_0 \nmid m$, todo primo de $m$ es $> q_0$, y la
degeneración $q_0^3 r \ge L$ con $r < L$ da la ventana multiplicativa

$$ \frac{L}{q_0^3\, P(q_0)} \le m < \frac{L}{P(q_0)}, \qquad
\sum_{r\ \mathrm{deg}} \frac{1}{\varphi(r)}
\;\le\; \sum_{q_0 \in \mathcal P} \frac{1}{\varphi(P(q_0))}\,
S(q_0), \quad S(q_0) := \sum_{m\ \mathrm{en\ ventana}} \frac{1}{\varphi(m)} $$

(igualdad salvo visibilidad, que solo quita términos). Además
$P(q_0) \le r < L$ fuerza $\theta(q_0) \lesssim \log L$: $q_0 =
O(\log L)$, hecho que se usa en §2.

## 2. Teorema A [T]: con $L \ge q_0^9$ no hay degenerados visibles con $\omega \le 2$

*Prueba.* Sea $r$ visible, $r \mid d \in V^+(L)$. El mayor primo de $d$
cumple Nathanson en $m=1$: $p_1(d)^3 < L$ (la condición $q^3 < L/r$
del doc 117b §3.2 con $r = 1$); luego todo primo de $r$ es
$< L^{1/3}$ (verificado sin excepción en J1: 0 violaciones de $p_1$).
Si $\omega(r) \le 2$ entonces $r < L^{2/3}$. Si además $L \ge q_0^9$,
la degeneración da $r \ge L/q_0^3 \ge L^{2/3}$: contradicción. ∎

**Agudeza [M].** La hipótesis no sobra a escala de ventana: en los
marcos A/B hay 188 instancias degeneradas visibles, con histogramas
$\omega = \{1{:}\,1,\ 2{:}\,9\}$ ($z=29$ A), $\{1{:}\,31,\ 2{:}\,47\}$
($z=37$), $\{2{:}\,100\}$ ($z=29$ B) — y **las 188 cumplen
$L < q_0^9$** (J1, 0 violaciones de la implicación
$\omega \le 2 \Rightarrow L < q_0^9$). A esas escalas ($L \le 2\cdot
10^7 < 5^9$) los primos y semiprimos degenerados *sí* son visibles;
en el régimen puente, $P(q_0) < L$ fuerza $q_0 = O(\log L)$, luego
$q_0^9 = O((\log L)^9) \ll L$ y la hipótesis es automática: allí la
obstrucción de conductores primos del doc 116c §2 queda fuera de la
subcapa degenerada por teorema, no por medición.

## 3. Estratos de $q_0$ grande: cerrados por Euler [T]

Sin ventana, $S(q_0) \le \prod_{p \in \mathcal P,\, p > q_0}
\frac{p}{p-1} = E(z) / \prod_{p \le q_0} \frac{p}{p-1}$ con
$E(z) := \prod_{p \in \mathcal P} \frac{p}{p-1}$. Componiendo con
$1/\varphi(P(q_0)) = \prod_{p<q_0} \frac{1}{p-1}$ el producto telescopa
a $\frac{1}{\varphi(P(q_0))} S(q_0) \le E(z)\,\frac{q_0-1}{q_0 P(q_0)}$,
con $E(z)\,\Pi_{\mathrm{full}} = \kappa_0(z) :=
\prod_{p \in \mathcal P} (1 - \frac{1}{(p-1)^2})$. Los términos decaen
por factor $\ge q_0$ por paso, luego la cola desde $Q$ es
$\le 2E(z)/P(Q)$. Con $q_0^\dagger := \min\{q : P(q) \ge (\log z)^3\}$
y $E(z) \le 0{,}66 \log z$ (J3: cocientes $0{,}6517/0{,}6267/0{,}6207/
0{,}6169$ a $z = 17/29/37/53$, límite $e^\gamma/3 \approx 0{,}5937$):

$$ \sum_{q_0 \ge q_0^\dagger} \frac{S(q_0)}{\varphi(P(q_0))}
\le \frac{2 E(z)}{(\log z)^3} \le \frac{1{,}32}{(\log z)^2}
= \frac{1{,}32}{K M_\infty \log z}\, Q_c\, (1+o(1)), $$

es decir cociente a $Q_c$ de $O(1/\log z) \to 0$: **cerrado con
holgura**. Valores (J3): $q_0^\dagger = 11$ ($z=17$), $13$
($z = 29, 37, 53$); cotas de cola absolutas $0{,}162/0{,}111/0{,}095/
0{,}078$ (a $z=29$: $0{,}22\,Q_c$); $\kappa_0(z) = 0{,}8926 \searrow
0{,}8834$, $\kappa_0(\infty) \approx 0{,}88$.

## 4. La obstrucción exacta: la ruta Euler pierde $(\log z)^2$ [T/M]

Para el estrato $q_0 = 5$ la cota de §3 sin ventana vale
$E(z)\cdot\frac{4}{5} = \Theta(\log z)$, mientras $Q_c =
\Theta(1/\log z)$: la ruta pierde exactamente $(\log z)^2$ y **ninguna
variante que ignore la ventana y la visibilidad puede cerrar los
estratos finitos $q_0 < q_0^\dagger$**. Cuantificado en ventana (J4):
la cota Euler del estrato 5 es $3{,}25\,Q_c$ ($z=29$) y $3{,}72\,Q_c$
($z=37$), creciendo, contra la suma verdadera visible de
$0{,}199\,Q_c$ y $0{,}355\,Q_c$ (J2). La suma *total* en banda
degenerada (visible + no visible, sobre todo el enum) es
$0{,}716\,Q_c$ ($z=29$ A), $1{,}405\,Q_c$ ($z=37$), $0{,}025\,Q_c$
($z=29$ B) [M, J4]: es $\Theta(Q_c)$ con constante $> 1$ a $z=37$,
así que incluso el enunciado sin visibilidad no tiene constante
pequeña en ventana — la diferencia entre $1{,}405$ y $0{,}355$ es
exactamente el trabajo de la visibilidad, y el Teorema A dice dónde
entra en el régimen puente (expulsa $\omega \le 2$).

## 5. El residuo mínimo: (SLAB-V) [P]

Lo que queda, exactamente, es: para cada $q_0 < q_0^\dagger$ fijo,

**(SLAB-V):** $\displaystyle\sum_{r\ \mathrm{deg,\ visible,\ estrato}\
q_0} \frac{1}{\varphi(r)} \le C(q_0)\, Q_c$ uniformemente en
$L \in \mathcal L_A$.

Por el apareamiento divisor–cofactor ($r \mid d \in V^+(L)$,
$r \ge L/q_0^3 \Rightarrow$ cofactor $c = d/r < q_0^3$), (SLAB-V) es
la contabilidad de la **rebanada superior** $\{d \in V^+(L) :
\log(L/d) \le 3 \log q_0\}$ de $V^+$ con pesos
$\sum_{c \mid d,\, c < q_0^3} \varphi(c)/\varphi(d)$. Es la región de
medida cero del doc 117b §3.2 — el caso en que la dicotomía del doc
116c §4 dice que el conteo cierra — pero exige la contabilidad de
Mertens de la rebanada *con las condiciones de Nathanson del testigo*,
que es lo que la ruta Euler de §3 no ve. Mismo estatus que (REG-Π)
en el doc 117b: ruta identificada, contabilidad no ejecutada.

## 6. Medición [M] (J0/J2)

Caracterización cruzada sin excepción (J0: $\widehat V^{(r)} = \{1\}$
vía conteo de $\widehat S$ ⟺ $q_0^3 r \ge L$, en las 246 instancias de
banda). Sumas degeneradas visibles máximas por nivel, **todas del
estrato $q_0 = 5$**: $0$ ($z=17$, banda vacía), $0{,}1993\,Q_c$
($z=29$ A), $0{,}3554\,Q_c$ ($z=37$), $0{,}0252\,Q_c$ ($z=29$ B).
Ningún estrato $q_0 \ge 7$ aporta en ventana. Compatible con
(SLAB-V) con $C(5) \approx 0{,}36$ a estas escalas, sin decidirla.

## 7. Veredicto

La subcapa degenerada **no queda cerrada [T] en su totalidad**, y la
obstrucción está identificada con precisión: cerrados los estratos
$q_0 \ge q_0^\dagger$ [T, constante $\to 0$] y excluido $\omega \le 2$
en el régimen puente [T, Teorema A]; el residuo es finito
($q_0 \in \{5, 7, 11\}$ en la práctica) y es exactamente (SLAB-V), la
rebanada superior de $V^+$ — ni conteo elemental (pierde
$(\log z)^2$, §4) ni la obstrucción prima del doc 116c (fuera por el
Teorema A). (BAND-LOW) hereda este residuo intacto por su subcapa
degenerada. Goldbach no se declara demostrado.

## 8. Reproducibilidad

`python3 118b1_degenerada.py | tee salida_doc118b1_local.txt` —
stdlib; importa `117b_banda_delta.py` (que encadena 110→109). J0:
caracterización degenerada por conteo de $\widehat S$ contra
$q_0^3 r \ge L$. J1: $p_1^3 < L$ en todo visible y agudeza
$\omega \le 2 \Rightarrow L < q_0^9$. J2: sumas por estrato contra
$Q_c$. J3: $E(z)$, $\kappa_0(z)$, $q_0^\dagger$, cola. J4: total sin
visibilidad y cota Euler del estrato 5. **PASS global.**
