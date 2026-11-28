# Doc 116c — El residuo (BAND): la cota cruda es falsa en el marco asintótico

Sección 116c del doc 116. Objeto único: el residuo **(BAND)** del doc 114 §7,

$$\sum_{r\ \mathrm{visible},\ u_r<u_0,\ \Pi(r)<1/2}
\frac{1-\Pi(r)^2}{\varphi(r)} \;\le\; C\cdot Q_{\mathrm{completo}}
\qquad\text{uniformemente en el marco,}$$

con $u_0 = 3/2$, $u_r = \ln(L/r)/\ln z$, $\Pi(r) = \prod_{p\in\mathcal P,\,p\nmid r}
(1-\frac1{p-1})$, $Q_{\mathrm{completo}} = \sum_r \Pi(r)^2/\varphi(r)$ (doc 110),
y $t(r)^+ \le (1-\Pi(r)^2)/\varphi(r)$ la cota cruda que legitima el doc 114 §6
vía $S(r)\le 1$. Definiciones exactas: docs 109/110 ($\mathcal P$ = primos
impares $<z$, $p\nmid A$, $A=3$; $r$ libre de cuadrados sobre $\mathcal P$;
visible $\iff r\in V^+(L)$, pues $V^+(L)$ es cerrado por divisores, doc 114
Lema 2). Niveles: $\mathcal L_A = \{\lfloor D_A/s\rfloor\}$, $D_A = z^5$,
$s$ primo de intercambio en $[z, z^{10/3}]$; luego
$u_L := \ln L/\ln z$ recorre valores $\to u^*$ para cualquier
$u^*\in[5/3,\,4]$ cuando $z\to\infty$.

## 0. Veredicto en breve

**(BAND), tal como está formulada para la cota cruda, es FALSA en el marco
asintótico [T].** La obstrucción es exacta y explícita: los conductores
primos $r=p$ del slice superior. A niveles con $u_L\to u^*\in(3/2,\,9/4)$
— en particular $u^*=5/3$ y $u^*=2$, ambos en $\mathcal L_A$ — la familia
$\omega=1$ sola aporta

$$\liminf_{z\to\infty}\sum_{\mathrm{banda}}\frac{1-\Pi(r)^2}{\varphi(r)}
\;\ge\; B_1(u^*) \;=\;
\log\frac{\min(1,\,u^*/3)}{u^*-u_0} \;>\;0,
\qquad B_1(5/3)=\log\tfrac{10}{3}=1{,}2040,\quad
B_1(2)=\log\tfrac43=0{,}2877,$$

mientras $Q_{\mathrm{completo}} \sim K M_\infty/\ln z \to 0$ con
$K\approx1{,}48$ y $M_\infty\in[1{,}16181601,\,1{,}16182763]$ (doc 115).
Luego el cociente suma$/Q_c$ crece como $(B_1/KM_\infty)\ln z\to\infty$:
ninguna constante $C$ sirve. Lo mismo vale, vía cofactores fijos y cajas
explícitas, en **todo** nivel $u^*\in[5/3,4]$ y para **todo** $\omega\ge1$
(§2). El conteo cierra exactamente las sub-bandas cuya región límite de
Mertens tiene medida cero (§4): la dicotomía es exacta.

Lo que NO se refuta: la masa real $\sum t(r)^+$ ni (H2). Lo que muere es la
ruta de conteo vía $\delta(r)\le 1-\Pi(r)$: en la familia obstructora
$S(p)$ es del orden de $\Pi(p)$, no de $1$ (en ventana: $S=0{,}650$ frente a
$\Pi=0{,}487$; el factor perdido por término es $(1-\Pi^2)/(S^2-\Pi^2)$,
que $\to\infty$ cuando $\Pi\to0$ con $S/\Pi$ acotado). El residuo mínimo
exacto que queda es de defectos, no de conteo:

**(BAND-δ):** $\displaystyle\sum_{r\in\mathrm{banda}}
\frac{\delta(r)^+}{\varphi(r)} \le C\cdot Q_{\mathrm{completo}}$,
$\delta(r)=S(r)-\Pi(r)$,

que implica la forma que la cadena consume, (BAND-S)
$\sum_{\mathrm{banda}} t(r)^+\le \frac32 C\,Q_c$ [T, módulo $S(r)\le1$ del
doc 114 §6]. (BAND-δ) es (FL-DESP) en media $L^1$ sobre $u_r<u_0$.
El teorema condicional del doc 114 se actualiza sustituyendo (BAND) por
(BAND-δ). Goldbach permanece abierto; (H2) no queda ni probada ni refutada.

## 1. Caracterización de {Π(r) < 1/2} [T]

Sea $\Pi_{\mathrm{full}} = \prod_{p\in\mathcal P}(1-\frac1{p-1})$
($=\Pi(1)$, el producto completo). Como $\Pi(r)$ excluye exactamente los
factores con $p\mid r$ (definición del doc 109/110, verificada a precisión
de máquina en J0: error $\le1{,}7\times10^{-16}$ en $z=17,29,37$):

$$\Pi(r) \;=\; \frac{\Pi_{\mathrm{full}}}
{\prod_{p\mid r}\bigl(1-\frac1{p-1}\bigr)}
\;=\; \Pi_{\mathrm{full}}\prod_{p\mid r}\Bigl(1+\frac1{p-2}\Bigr)
\;\ge\; \Pi_{\mathrm{full}},$$

con igualdad $\iff r=1$. En particular la dirección es la de *división*:
quitar factores $<1$ agranda el producto, y $\Pi(r)$ **crece** con el
soporte de $r$. Por tanto

$$\Pi(r)<\tfrac12 \iff
\prod_{p\mid r}\Bigl(1+\frac1{p-2}\Bigr) < \frac{1}{2\Pi_{\mathrm{full}}},
\qquad\text{y tomando logaritmos,}\qquad
\sum_{p\mid r}\frac1{p-1}\;\le\;
\log\frac{1}{2\Pi_{\mathrm{full}}} + O(1),$$

pues $\frac1{p-1}\le\log(1+\frac1{p-2})\le\frac1{p-2}$: la condición es una
cota sobre la masa armónica del soporte, $\sum_{p\mid r}1/(p-1)$. Es decir:
$\Pi(r)$ es pequeño cuando $r$ contiene **pocos primos pequeños** — los
$p\nmid r$ que quedan en el producto incluyen a los primos chicos y hunden
$\Pi(r)$ hacia $\Pi_{\mathrm{full}}$.

Consecuencias [T]:

- $\Pi_{\mathrm{full}}(z)$ decrece cruzando $1/2$ entre $z=17$ y $z=19$
  (J1: $\Pi_{\mathrm{full}}(17)=0{,}515625>1/2>0{,}4834=
  \Pi_{\mathrm{full}}(19)$). A $z=17$ la banda $\{\Pi<1/2\}$ es **vacía**
  ($\Pi(r)\ge\Pi_{\mathrm{full}}>1/2$ para todo $r$; J2 lo confirma:
  suma $0$ en los cuatro niveles).
- Para un primo: $\Pi(p)=\Pi_{\mathrm{full}}\cdot\frac{p-1}{p-2}$
  (J1, error $\le1{,}1\times10^{-16}$). Como $\frac{p-1}{p-2}\le\frac43$
  ($p\ge5$) y $\Pi_{\mathrm{full}}\to0$ (§2.1), **todo** $p\in\mathcal P$
  cumple $\Pi(p)<1/2$ para $z$ grande: la condición $\Pi<1/2$ no excluye
  a ningún primo del slice superior asintóticamente. En ventana sí
  excluye: a $z=29$ el primo $p=7$ del intervalo de $u_L=2$ tiene
  $\Pi(7)=0{,}5578>1/2$ y queda fuera (J2: banda $\omega=1$ vacía a
  $z=29$, $u_L=2$); a $z=37$ ya entran $p=7,11$
  ($\Pi=0{,}4875$, $0{,}4514$).
- Por Mertens, $\log\Pi_{\mathrm{full}} = -\sum_{p\in\mathcal P}
  \frac1{p-1} + O(1) = -\log\log z + O(1)$, luego
  $\Pi_{\mathrm{full}} \sim K/\ln z$ con $K>0$ (J3c: $\Pi_{\mathrm{full}}
  \cdot\ln z = 1{,}4668,\,1{,}4600,\,1{,}4690,\,1{,}4790,\,1{,}4809$ a
  $z=37,53,101,1009,10007$; $K\approx1{,}48$). El umbral
  $1/(2\Pi_{\mathrm{full}})\sim\frac{\ln z}{2K}\to\infty$: la condición
  sobre el soporte se *debilita* con $z$; asintóticamente admite
  $\sum_{p\mid r}1/(p-1)$ hasta $\log\log z - \log(2K) + O(1)$, es decir,
  casi todo el rango posible.

## 2. La obstrucción [T]: la suma cruda de la banda no es O(Q_c)

### 2.1 El lado derecho tiende a 0

$Q_{\mathrm{completo}} = \Pi_{\mathrm{full}}\cdot M_{\mathrm{full}}$
(identidad de Euler, doc 110 [T]) y $M_{\mathrm{full}}\to M_\infty\in
[1{,}16181601,\,1{,}16182763]$ (doc 115, con cola telescópica rigurosa).
Con $\Pi_{\mathrm{full}}\sim K/\ln z$ (§1):

$$Q_{\mathrm{completo}} \;\sim\; \frac{K M_\infty}{\ln z}\;\longrightarrow\;0,
\qquad KM_\infty \approx 1{,}72$$
(J3c: $Q_c\cdot\ln z = 1{,}6931,\,1{,}6894,\,1{,}7035,\,1{,}7181,\,1{,}7205$
a $z=37,\dots,10007$; y $1{,}4809\times1{,}16182 = 1{,}7205$: consistente).

### 2.2 La familia ω = 1 en niveles u* ∈ (u₀, (3/2)u₀)

Fija $u^*\in(3/2,\,9/4)$ (con $u_0=3/2$: $(u_0,\frac32 u_0)$) y niveles
$L\in\mathcal L_A$ con $u_L\to u^*$ (existen: $\log_z s$ recorre
$[1,10/3]$). Sea $x_p=\ln p/\ln z$ y el intervalo abierto

$$I(u^*) \;=\; \bigl(u^*-u_0,\ \min(1,\,u^*/3)\bigr)\;\ne\;\emptyset .$$

Para todo primo $p\in\mathcal P$ con $x_p\in I(u^*)$, cuando $z\to\infty$:

- **Visible:** $p\in V^+(L)\iff p^3<L\iff x_p<u_L/3$ ✓ (condición de
  Nathanson en $m=1$; es la única para $\omega=1$).
- **Banda:** $u_r = u_L - x_p + o(1) < u_0 \iff x_p > u^*-u_0$ ✓.
- **Π(p) < 1/2:** $\Pi(p)=\Pi_{\mathrm{full}}\frac{p-1}{p-2}\le
  \frac43\Pi_{\mathrm{full}}\to0$ ✓ (§1).
- **Sumando:** $(1-\Pi(p)^2)/(p-1) = (1-o(1))/(p-1)$ uniformemente.

Por Mertens, $\sum_{x_p\in I} \frac1{p-1} \to
\log\frac{\min(1,u^*/3)}{u^*-u_0} = B_1(u^*) > 0$. Por tanto

$$\liminf \sum_{\mathrm{banda}}\frac{1-\Pi^2}{\varphi} \;\ge\; B_1(u^*)>0
\qquad\text{frente a}\qquad C\cdot Q_c \sim \frac{CKM_\infty}{\ln z}\to0,$$

para **toda** constante $C$: (BAND) es falsa. Valores:
$B_1(5/3)=\log\frac{10}{3}=1{,}2040$ (intervalo $x_p\in(1/6,\,5/9)$),
$B_1(2)=\log\frac43=0{,}2877$ ($x_p\in(1/2,\,2/3)$). Verificación numérica
directa de las sumas de primos (J3a): a $Z=10^8$, intervalo de $u_L=2$:
$0{,}2866$ (error $0{,}0011$); a $Z=10^{12}$, intervalo de $u_L=5/3$:
$1{,}1917$ (error $0{,}0123$).

### 2.3 Cofactores fijos: todo ω ≥ 1 está obstruido en esos niveles

Sea $m_0>1$ fijo, libre de cuadrados sobre $\mathcal P$, y $r = p\,m_0$
con $x_p\in I(u^*)$. Para $z$ grande: $p\nmid m_0$ ✓; la visibilidad de
$r$ añade solo condiciones de Nathanson con primos fijos
($p < L/(\text{fijo})$: sin restricción asintótica, pues $x_p<u^*/3<u^*$);
la banda se desplaza en $\log_z m_0\to0$ ✓; $\Pi(r)\to
\Pi_{\mathrm{full}}\prod_{q\mid m_0}\frac{q-1}{q-2}\cdot(1+o(1))\to0$ ✓;
y el sumando es $(1-o(1))/((p-1)\varphi(m_0))$. Luego la subfamilia
$\{p\,m_0\}$ aporta $\to B_1(u^*)/\varphi(m_0) > 0$. J3a2: con $m_0=35$
($\omega$ total $3$): predicho $\log(4/3)/24 = 0{,}0120$, medido
$0{,}0119$ a $Z=10^8$. Ningún $\omega$ se salva por conteo.

### 2.4 Los niveles restantes u* ∈ [9/4, 4]

Ahí $I(u^*)=\emptyset$ para $\omega=1$ ($B_1=0$), pero cajas explícitas
$[c_i-\delta,\,c_i+\delta]$ dentro de la región
$\{\sum x_i > u^*-u_0,\ \sum_{i\le j}x_i + 2x_j < u^*$ ($j$ impar)$\}$,
con $x_i$ acotados lejos de $0$, dan cotas inferiores positivas del
límite de Mertens (J3b; las cotas de caja son muy flojas — solo certifican
$>0$):

| $u^*$ | familia | $B(u^*)\ge$ |
|---|---|---|
| 5/3 | ω=1 (exacta) | 1,2040 |
| 2 | ω=1 (exacta) | 0,2877 |
| 5/2 | ω=2 (caja) | 4,0×10⁻³ |
| 3 | ω=2 (caja) | 2,4×10⁻³ |
| 7/2 | ω=3 (caja) | 2,5×10⁻⁵ |
| 4 | ω=3 (caja) | 1,3×10⁻⁵ |

Conclusión: para **todo** $u^*\in[5/3,4]$,
$\liminf\sum_{\mathrm{banda}}(\text{cruda}) \ge B(u^*) > 0 \gg Q_c$.
La pérdida de la cota cruda no es $(\ln z)^2$ «de la ruta» (doc 114):
es intrínseca al mayorante — la suma misma viola la cota.

## 3. Medición en ventana [M] (J2, J4, J5)

Descomposición de la banda cruda por $\omega(r)$, niveles representativos
$L=\lfloor z^{u^*}\rfloor$, $u_0=1{,}5$ (cocientes a $Q_c$):

| marco | $z$ | $u_L$ | crudo/$Q_c$ | real/$Q_c$ | domina |
|---|---|---|---|---|---|
| A | 17 | cualquiera | 0 | 0 | banda vacía ($\Pi_{\mathrm{full}}>1/2$) |
| A | 29 | 3 | 0,0191 | 0,0191 | ω=2 |
| A | 37 | 5/3 | 0,2710 | 0,1155 | ω=1 (todo) |
| A | 37 | 2 | 0,4408 | 0,0948 | ω=1 (todo) |
| A | 37 | 3 | 0,0791 | 0,0493 | ω=2 |
| A | 37 | 4 | 0,0006 | 0,0000 | ω=3 |
| B | 29 | 3 | 0,0191 | 0,0191 | ω=2 |

Sobre los niveles muestreados del doc 114 ($L=D//s$, 60 por marco), el
máximo de la parte $\Pi<1/2$ es $0{,}1526$ ($z=29$) y $0{,}4408$ ($z=37$):
el $0{,}4408$ **coincide exactamente** con la columna «de la cual
$\Pi(r)<1/2$» del doc 114 §7 — validación cruzada de ambas mediciones.
El total $0{,}6991$ del doc 114 incluye la sub-banda $\Pi\ge1/2$, cerrada
[T] (§4). La ventana no muestra la divergencia de §2 porque a $z\le37$,
$Q_c\approx0{,}47$ y la suma $\omega=1$ aún está por debajo de su límite
($0{,}2067$ frente a $B_1(2)=0{,}2877$ a $z=37$, nivel $u_L=2$): el
cociente crece $\sim\ln z$, y a $z=37$ vale $\approx0{,}44$.

La cancelación que el conteo no ve (J4; $z=37$, $L=37^2$, familia ω=1):

| $p$ | $\Pi(p)$ | $S(p)$ | $\delta(p)$ | crudo | $t(p)^+$ real | factor |
|---|---|---|---|---|---|---|
| 7 | 0,4875 | 0,6500 | +0,1625 | 0,1271 | 0,0308 | 4,1× |
| 11 | 0,4514 | 0,5833 | +0,1320 | 0,0796 | 0,0137 | 5,8× |

$S(p)$ está a distancia $O(\Pi(p))$ de $\Pi(p)$, no cerca de $1$: la cota
$\delta\le1-\Pi$ del doc 114 §6 sustituye $S$ por su máximo trivial y
pierde exactamente el factor $(1-\Pi^2)/(S^2-\Pi^2)$ por término. La masa
real de la banda en ventana sigue $\le0{,}3806\,Q_c$ (doc 114): (H2) no
queda tocada por esta sección; queda tocada la *ruta*.

## 4. Qué cierra el conteo, y el residuo mínimo exacto

**Cerrado por conteo [T]:**
- Sub-banda $\{\Pi(r)\ge c_0\}$ (doc 114): $(1-\Pi^2)/\varphi \le
  (\Pi^2/\varphi)/c_0^2$ y $\sum\Pi^2/\varphi = Q_c$ ⇒ cota $Q_c/c_0^2$.
  Reverificado (J5): a $c_0=1/2$, medido $\le0{,}43\,Q_c$ en los niveles
  de arriba, frente a la cota teórica $4\,Q_c$.
- Sub-bandas de región límite de medida cero: p. ej. la familia
  $\omega=1$ a niveles $u_L\ge\frac32 u_0$ (intervalo $I(u^*)$ vacío),
  y la banda entera a $z\le17$ (§1). El conteo directo — suma sobre la
  región con pesos $1/\varphi$ — tiene por límite la integral de Mertens
  de la región (veces potencias de $\log\log z$ por factores pequeños
  libres); cierra $\iff$ ese límite es $0$. La banda completa tiene
  región de medida positiva a todo $u^*\in[5/3,4]$ (§2): **ninguna
  variante del conteo del mayorante puede cerrar (BAND)**. Esta es la
  obstrucción exacta.

**Residuo mínimo.** Lo que la cadena (docs 113–114) consume de la banda
es la masa real, no el mayorante:

**(BAND-S):** $\sum_{r\in\mathrm{banda}} t(r)^+ \le C\,Q_{\mathrm{completo}}$,
$\quad t(r)=\delta(r)\bigl(2\Pi(r)+\delta(r)\bigr)/\varphi(r)$.

Y sobre la banda, $\Pi(r)<1/2$ y $S(r)\le1$ (doc 114 §6: [T] si
$z_A\le71$, en general módulo (R-FLOW)) dan
$t(r)^+ = \delta^+\,(\Pi+S)/\varphi \le \frac32\,\delta(r)^+/\varphi(r)$
[T]. Luego basta — y es la forma mínima que no pasa por el mayorante
refutado —

**(BAND-δ):** $\displaystyle\sum_{r\in\mathrm{banda}}
\frac{\delta(r)^+}{\varphi(r)} \;\le\; C\cdot Q_{\mathrm{completo}},$

un enunciado de **defectos en media $L^1$** sobre $u_r<u_0$: (FL-DESP)
promediada en la banda, no conteo. Nota sobre la sub-banda $u_r<1$: el
«la cota cruda es exacta» del doc 114 §7 era [M] de ventana; el mismo
argumento de región (ω=2 a $u^*=5/3$: $x_1+x_2>2/3$, $x_1<5/9$, región de
medida positiva) muestra que la suma cruda sobre $\{u_r<1\}$ también
tiene $\liminf>0\gg Q_c$, así que tampoco ahí hay cierre por conteo; si
$S(r)=1$ persiste asintóticamente en esa familia es territorio de
(DESP-SIGN), no de esta sección.

**Cadena actualizada:** (R-FLOW-71+) + (REM-PERS-small) + **(BAND-δ)**
⟹ (DESP-SIGN) + banda ⟹ (FL-DESP⁺) ⟹ (H2⁗) ⟹ (H2), con la constante del
doc 113 (≈6,4) más $\frac32 C_\delta$. (BAND) en su forma cruda sale de
la lista de residuos: está refutada, no abierta. Siguen fuera, sin
cambios: (DIST-𝒟) (tras el doc 115) y la separación de capas G0.
Goldbach permanece abierto.

## 5. Reproducibilidad

`116c_banda_conteo.py` (solo stdlib; importa la maquinaria de los docs
109/110 sin tocarla): chequeos J0–J5, salida en
`salida_doc116c_local.txt` — **PASS**, 0,7 s. J0/J1: identidades de
$\Pi$ a precisión de máquina y umbral $\Pi_{\mathrm{full}}=1/2$ entre
$z=17$ y $19$. J2: descomposición por $\omega$ y máximo muestreado
($0{,}4408$ a $z=37$ = columna $\Pi<1/2$ del doc 114). J3: sumas de
Mertens a $Z=10^8$/$10^{12}$, cofactor fijo $m_0=35$, cajas de §2.4 y
estabilización $\Pi_{\mathrm{full}}\ln z\to K$, $Q_c\ln z\to KM_\infty$.
J4: $S(p)$ frente a la cota cruda. J5: sub-banda $\{\Pi\ge1/2\}\le4Q_c$.
