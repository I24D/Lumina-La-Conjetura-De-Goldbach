# Doc 117b — (BAND-δ) vía (DESP-SIGN) + (BAND-MOD): reducción exacta, y (BAND-MOD) reducido a la capa de nivel bajo

Sección 117b del doc 117. Objeto único: el residuo **(BAND-δ)** del doc
116c §4,

$$\sum_{r\in\mathrm{banda}}\frac{\delta(r)^+}{\varphi(r)}
\;\le\; C\cdot Q_{\mathrm{completo}},
\qquad
\mathrm{banda}=\{r\ \mathrm{visible}:\ u_r<u_0=\tfrac32,\ \Pi(r)<\tfrac12\},$$

con $\delta(r)=S(r)-\Pi(r)=(\widehat S(r)-\Pi(r))+E(r)$ [T, doc 113 §2],
$\widehat S$ la suma del modelo RI sin desplazar
$\widehat V^{(r)}=V^+(L/r;\ \mathcal P\setminus\mathrm{primos}(r))$ y $E(r)$
el error de desplazamiento. Definiciones: docs 109/110/113; niveles
$\mathcal L_A$ (doc 115). (BAND-δ) implica (BAND-S)
$\sum_{\mathrm{banda}}t(r)^+\le\frac32 C\,Q_c$ [T, mód. $S(r)\le1$ del
doc 114 §6; doc 116c §4].

## 0. Veredicto en breve

1. **Reducción exacta [T] (§1):** (BAND-δ) se sigue de **(DESP-SIGN)**
   ($E(r)\le0$) más un enunciado sobre el modelo solo,
   **(BAND-MOD):** $\sum_{\mathrm{banda}}(\widehat S(r)-\Pi(r))^+/\varphi(r)
   \le C\,Q_c$. La reducción es puntual ($\delta^+\le(\widehat S-\Pi)^+$
   si $E\le0$), está verificada a precisión de máquina (K1), y es
   **acíclica**: (DESP-SIGN) descansa sobre (FLOW-COMP) + (REM-PERS-small)
   (docs 116a/116b), que no tocan la banda; (BAND-MOD) no involucra $E$.
2. **La suma con signo no tiene forma cerrada (§2):** la identidad de
   intercambio del doc 110 no se transfiere a $\widehat S$ — al permutar,
   la suma interna es una suma de Mertens **truncada por el umbral**
   $r\lambda(d')\ge L$, no un producto de Euler (demostrado por la
   dependencia explícita en $L$). Y no hace falta: en la banda,
   $\widehat S-\Pi$ es **de un solo signo** [M: 246/246 instancias
   positivas], luego la parte positiva *es* el valor absoluto y no hay
   cancelación que explotar. La suma con signo sobre todo el enum es
   también positiva en los niveles probados ($+1{,}27\,Q_c$, $+1{,}57\,Q_c$).
3. **(BAND-MOD) queda reducido a un residuo mínimo nombrado (§3):**
   sobre $\{u_r\ge u_0^*\}$ cierra [T] por la FL del modelo (docs
   112–113) módulo la suma de región **(REG-Π)**
   $\sum_{\mathrm{banda}}\Pi(r)/\varphi(r)\le C\,Q_c$ — medida en
   $0{,}17$–$0{,}46\,Q_c$ sin crecimiento (K6), con ruta de prueba por
   el método de regiones de Mertens del doc 116c §2 —. El residuo es la
   capa de nivel bajo **(BAND-LOW)**: la misma cota sobre
   $\{r\in\mathrm{banda}:u_r<u_0^*\}$, que contiene la subcapa
   exactamente degenerada $\widehat V^{(r)}=\{1\}$, caracterizada [T]
   en §3.2. Su región límite de Mertens tiene medida cero (capa
   frontera $u_r\to0$): la dicotomía del doc 116c §4 apunta a su favor,
   pero la cota de densidad no está probada. **(BAND-LOW) [P].**
4. **Medición decisiva [M] (§4, K1/K3/K4/K5):** sobre los 240 niveles
   muestreados de los marcos A ($z=17,29,37$) y B, los máximos de
   $\sum_{\mathrm{banda}}(\widehat S-\Pi)^+/\varphi$ son
   $0$, $0{,}1028$, $0{,}2103$, $0{,}0128$ veces $Q_c$ — acotados en
   ventana, sin explosión; la ventana no decide (BAND-MOD) (mismo
   caveat preasintótico del doc 116c §3), pero es compatible con él y
   con (BAND-δ) vía la reducción (máximos de $\sum\delta^+/\varphi$:
   $0{,}0198$ a $0{,}0933\,Q_c$).

Cadena: (DESP-SIGN) + (BAND-LOW) + (REG-Π) ⟹ (BAND-MOD) ⟹ (BAND-δ)
⟹ (BAND-S). Goldbach permanece abierto; (H2) no queda probada.

## 1. La reducción [T]: (BAND-δ) ⟸ (DESP-SIGN) + (BAND-MOD)

**Teorema (reducción).** Supóngase (DESP-SIGN): $E(r)\le0$ para todo
$r$ visible. Entonces, puntualmente en todo $r$,

$$\delta(r)^+ \;\le\; \bigl(\widehat S(r)-\Pi(r)\bigr)^+,$$

y por tanto (BAND-MOD) implica (BAND-δ) con la misma constante.

*Prueba.* $\delta(r)=(\widehat S(r)-\Pi(r))+E(r)\le\widehat S(r)-\Pi(r)$
porque $E(r)\le0$ (doc 113 §2). La función $x\mapsto x^+$ es monótona:
$\delta(r)^+\le(\widehat S(r)-\Pi(r))^+$. Sumando sobre la banda con
peso $1/\varphi(r)$ se obtiene la implicación. ∎

**Orden de dependencias (sin circularidad).** En la cadena vigente,
(DESP-SIGN) se apoya en (FLOW-COMP) (doc 116b, compensación global del
flujo de Add/Rem a $r$ fijo) y en (REM-PERS-small), cerrado [T] (doc
116a). Ninguno de los dos usa la banda ni (BAND-δ): la banda entra en
la ruta del doc 113 *después* de (DESP-SIGN), como sumando separado de
(FL-DESP⁺). Y (BAND-MOD) es un enunciado sobre $\widehat S$ solo: el
modelo sin desplazar, sin $E$ en ninguna parte. La reducción no
alimenta nada que (FLOW-COMP) consuma. ∎

**Verificación (K1) [M]:** en las 246 instancias de banda de los 240
niveles muestreados: $E(r)>0$ en **0** instancias (consistente con
(DESP-SIGN), doc 113/114), y la violación máxima de
$\delta^+\le(\widehat S-\Pi)^+$ es $1{,}1\times10^{-16}$ (cero a
precisión de máquina). Las sumas confirman el orden:
$\max_L\sum\delta^+/\varphi \le \max_L\sum(\widehat S-\Pi)^+/\varphi$
en los cuatro marcos (§4).

## 2. La suma con signo: doble conteo exacto, sin forma cerrada

### 2.1 Por qué el intercambio del doc 110 no se transfiere [T]

En el doc 110, la identidad $\sum_r\mu(r)\Pi(r)C(r)=\Pi_{\mathrm{full}}
M_V$ funciona porque $C(r)=\sum_{d\in V^+,\ r\mid d}\mu(d)/\varphi(d)$
es una convolución sobre el soporte **fijo** $V^+(L)$: al permutar, la
suma interna corre sobre los divisores de un $d$ fijo y factoriza como
$\prod_{p\mid d}1/(p-2)=g(d)$. Para $\widehat S$ el soporte depende de
$r$ dos veces (nivel $L/r$ y primos excluidos), y la permutación da,
con $\lambda(d')=\max\bigl(d',\max_{m\ \mathrm{impar}}p_m^2P_m\bigr)$
— de modo que $d'\in\widehat V^{(r)}\iff(d',r)=1$ y $r\lambda(d')<L$
(aritmética entera del doc 113) —:

$$\sum_{r}\frac{\widehat T(r)}{\varphi(r)}
=\sum_{d'}\frac{\mu(d')}{\varphi(d')}\,W_L(d'),
\qquad
W_L(d')=\!\!\sum_{\substack{r\in\mathrm{enum},\ (r,d')=1\\
r\lambda(d')\ge L}}\frac1{\varphi(r)},
\quad \widehat T(r)=\Pi(r)-\widehat S(r).$$

La identidad es exacta (K2: error $\le3{,}3\times10^{-16}$ en dos
niveles), pero $W_L(d')$ es una **suma de Mertens truncada por el
umbral** $L/\lambda(d')$: no factoriza. Demostración de la dependencia:
para $d'=385$ ($\lambda=9625$), $W_L=1{,}2702$ a $L=29^2$ y $0{,}2702$
a $L=29^3$ (K2) — la suma interna cambia con $L$, mientras que la
análoga del doc 110 es independiente de $L$. Esta es la obstrucción
exacta a una forma cerrada tipo (110.1)/(110.2) para el modelo: el
nivel propio $u_r$ acopla $r$ y $d'$ por tamaño, no por divisibilidad.

### 2.2 El signo en la banda: unilateral [M], y la suma con signo no cancela

Clasificación del signo de $\widehat S(r)-\Pi(r)$ sobre las 246
instancias de banda (K3), por $\omega(r)$ — conteos $(n_+,n_-,n_0)$ y
masas $\sum(\widehat S-\Pi)^\pm/\varphi$ sin normalizar:

| $\omega(r)$ | $n_+$ | $n_-$ | $n_0$ | masa $+$ | masa $-$ |
|---|---|---|---|---|---|
| 1 | 56 | 0 | 0 | 3,5374 | 0 |
| 2 | 187 | 0 | 0 | 0,3211 | 0 |
| 3 | 3 | 0 | 0 | 0,0001 | 0 |

**En ventana, $\widehat S(r)>\Pi(r)$ en el 100% de la banda**: la
parte positiva coincide con el valor absoluto nivel a nivel (K4:
$\max_L\sum|\widehat S-\Pi|/\varphi=\max_L\sum(\widehat S-\Pi)^+/
\varphi$ en los cuatro marcos). No se declara [T]: fuera de la banda
($u_r$ grande) el error de la FL sí oscila; la unilateralidad es un
fenómeno de la capa $u_r<3/2$ (cf. §3.2: en la subcapa degenerada es
[T]). Para $\omega=1$, visibilidad fuerza $u_r>\frac23 u_L\ge10/9$
[T: $p\in V^+(L)\Rightarrow p^3<L$], y las 56 instancias caen en
$u_r\in[1,3/2)$, todas positivas.

La suma con signo tampoco ofrece rescate fuera de la banda: sobre
todo el enum, $\sum_r(\widehat S(r)-\Pi(r))/\varphi(r)=+1{,}2722\,Q_c$
($z=29$, $L=841$) y $+1{,}5653\,Q_c$ ($z=37$, $L=410$) [M, K2] — del
mismo orden que $Q_c$ y **positiva**: el modelo sobrepasa a $\Pi$ en
media $1/\varphi$, no cancela. (BAND-MOD) es por tanto una cuestión de
tamaño del exceso del modelo, no de cancelación.

## 3. (BAND-MOD): cierre sobre $\{u_r\ge u_0^*\}$ y el residuo (BAND-LOW)

### 3.1 La parte alta cierra [T] módulo (REG-Π)

La FL clásica aplicada al modelo (doc 112 §2, doc 113 §2) da
$|\widehat S(r)-\Pi(r)|\le\rho_{\mathrm{FL}}(u_r)\,\Pi(r)$ para
$u_r\ge u_0^*$, con $u_0^*>0$ absoluto el umbral de la FL citada y
$\rho_{\mathrm{FL}}$ su función de error. Sobre
$\{r\in\mathrm{banda}:u_r\ge u_0^*\}$:

$$\sum\frac{(\widehat S-\Pi)^+}{\varphi}
\le\rho_{\mathrm{FL}}(u_0^*)\sum_{\mathrm{banda}}\frac{\Pi(r)}{\varphi(r)},$$

y basta la suma de región

**(REG-Π):** $\displaystyle\sum_{r\in\mathrm{banda}}
\frac{\Pi(r)}{\varphi(r)}\le C\,Q_{\mathrm{completo}}.$

Nota de estructura: $\Pi(r)/\varphi(r)=\Pi_{\mathrm{full}}\prod_{p\mid r}
\frac1{p-2}$ (identidad del doc 116c §1), luego (REG-Π) pide que la
masa de Mertens de la región banda con pesos $\prod1/(p-2)$ sea $O(1)$:
la región $\{u_r<3/2\}\cap\{\mathrm{visible}\}$ es, nivel a nivel, una
unión finita de politopos acotados en el espacio
$x_i=\ln p_i/\ln z$ (el método del doc 116c §2 la computa familia a
familia; la condición $\Pi(r)<1/2$ acota además
$\sum_{p\mid r}1/(p-1)\lesssim\log\log z$, doc 116c §1). Ruta de
prueba identificada, contabilidad de Mertens no ejecutada aquí:
**(REG-Π) [P]**, con medición [M] en K6: cocientes
$\sum_{\mathrm{banda}}\Pi/\varphi\,/\,Q_c$ de $0{,}17$ a $0{,}46$ en
niveles representativos de $z=29,37,53$, sin crecimiento visible.

### 3.2 La capa baja y su subcapa degenerada [T]

El residuo es la capa que la FL no cubre:

**(BAND-LOW) [P]:** $\displaystyle\sum_{r\in\mathrm{banda},\
u_r<u_0^*}\frac{(\widehat S(r)-\Pi(r))^+}{\varphi(r)}
\le C\,Q_{\mathrm{completo}}.$

Dentro de ella, la subcapa exactamente degenerada está caracterizada:

**Proposición [T].** $\widehat V^{(r)}=\{1\}\iff q_0^{\,3}\,r\ge L$,
donde $q_0=\min(\mathcal P\setminus\mathrm{primos}(r))$. En ese caso
$\widehat S(r)=1$ y $\widehat S(r)-\Pi(r)=1-\Pi(r)>0$.

*Prueba.* Un primo $q\in\mathcal P\setminus r$ está en
$\widehat V^{(r)}$ sii $q^3<L/r$ (condición de Nathanson en $m=1$,
aritmética entera del doc 113; $q<L/r$ es más débil). Si ningún primo
entra, ningún $d'$ con $\omega(d')\ge1$ entra: su primo mayor
$p_1\ge q_0$ viola la condición en $m=1$ ($p_1^3\ge q_0^3\ge L/r$).
Recíprocamente, si $q_0^3<L/r$ entonces $q_0\in\widehat V^{(r)}$.
Verificado sin excepción en todas las instancias de banda de K6. ∎

En ventana la capa baja domina la banda: 188 de las 246 instancias
muestreadas son degeneradas (K3/K5), con término máximo
$(1-\Pi)/\varphi$ de $0{,}1028\,Q_c$ ($z=29$) y $0{,}1822\,Q_c$
($z=37$) por instancia (K5); en niveles representativos la suma
degenerada vale $0$–$0{,}52\,Q_c$ (K6, conteos $n\le3$ por nivel:
chequeo puntual, no tendencia). Pero la capa es una **capa frontera**:
$u_r<u_0^*$ fijo con $u_r=\ln(L/r)/\ln z$, y la subcapa degenerada
exige $u_r\le\ln(q_0^3)/\ln z\to0$. En el límite de Mertens de la
región banda (doc 116c §2/§4), $\{u_r<u_0^*\}$ colapsa sobre la
frontera $\sum x_i=u^*$: región de **medida cero**, exactamente el
caso en que la dicotomía del doc 116c §4 dice que el conteo cierra.
Lo que falta para [T] es la cota de densidad de la masa
$(\widehat S-\Pi)^+/\varphi$ cerca de la frontera — eso es
(BAND-LOW), ni más ni menos. Nótese el paralelismo con el doc 114 §7:
allí la degeneración era del soporte *desplazado* ($S(r)=1$,
$u_r<1$); aquí es del *modelo*, y está caracterizada, no solo medida.

### 3.3 La parte no degenerada lleva el factor $\Pi$ [M]

Fuera de la subcapa degenerada, el exceso del modelo es relativo a
$\Pi(r)$, no absoluto (K5): $\max(\widehat S-\Pi)^+/\Pi(r)=0{,}692$
sobre las instancias no degeneradas de la banda, y el 100% de la masa
$(\widehat S-\Pi)^+$ no degenerada cumple $\widehat S-\Pi\le\Pi(r)$.
En niveles representativos (K6): $\max(\widehat S-\Pi)/\Pi=0{,}539$
($z=37$) y $0{,}870$ ($z=53$), siempre en instancias con $u_r\ge1$.
Esto es lo que la cota cruda del doc 114 perdía y el doc 116c
diagnosticó: el exceso real es $O(\Pi(r))$ donde el modelo es genuino,
y $\Theta(1)$ solo en la capa frontera de §3.2.

## 4. Medición decisiva [M] (K1/K3/K4/K5/K6)

Máximos sobre los 60 niveles muestreados de cada marco
($L=\lfloor D_A/s\rfloor$), cocientes a $Q_c$:

| marco | $z$ | $\sum\delta^+/\varphi$ | $\sum(\widehat S-\Pi)^+/\varphi$ | $\sum|\widehat S-\Pi|/\varphi$ |
|---|---|---|---|---|
| A | 17 | 0 | 0 | 0 |
| A | 29 | 0,0198 | 0,1028 | 0,1028 |
| A | 37 | 0,0933 | 0,2103 | 0,2103 |
| B | 29 | 0,0128 | 0,0128 | 0,0128 |

($z=17$: banda vacía, $\Pi_{\mathrm{full}}>1/2$, doc 116c §1.)
Lectura honesta: los cocientes crecen de $z=29$ a $z=37$ y los niveles
representativos a $z=53$ (K6) dan totales degenerado+no degenerado de
$0{,}52$ y $0{,}20\,Q_c$ con conteos mínimos ($n\le5$): la ventana es
**compatible** con (BAND-MOD) y con (BAND-δ) — ningún cociente supera
$\sim0{,}5\,Q_c$ — pero no los decide, por el mismo motivo que en el
doc 116c §3: a $z\le53$, $Q_c\approx0{,}47$ y las sumas aún están por
debajo de sus límites de Mertens. La decisión es el enunciado
(BAND-LOW) de §3.2, no más cómputo en ventana.

## 5. Veredicto

- **(BAND-δ) ⟸ (DESP-SIGN) + (BAND-MOD)** [T, §1], sin circularidad.
- **(BAND-MOD)** reducido a **(BAND-LOW)** [P, §3.2] más **(REG-Π)**
  [P con ruta de prueba identificada y medición favorable, §3.1];
  sobre $\{u_r\ge u_0^*\}$ el cierre es [T] dados (REG-Π) y la FL del
  modelo. La suma con signo no tiene forma cerrada [T, §2.1] y no
  cancela en la banda [M, §2.2]: el exceso del modelo es unilateral
  y de tamaño $O(\Pi(r))$ fuera de la capa frontera [M, §3.3].
- Cadena actualizada: (FLOW-COMP) ⟹ (DESP-SIGN) [con (REM-PERS-small)
  cerrado]; (DESP-SIGN) + (BAND-LOW) + (REG-Π) ⟹ (BAND-δ) ⟹ (BAND-S),
  que es lo que la ruta de los docs 113–114 consume de la banda.
  (H2) sigue abierta; Goldbach no se declara demostrado.

## 6. Reproducibilidad

`117b_banda_delta.py` (solo stdlib; importa la maquinaria del doc 110
y reproduce `in_Vhat` del doc 113): chequeos K0–K6, salida en
`salida_doc117b_local.txt` — **PASS**, 0,6 s. K0: identidades
$\widehat S(1)=S(1)$, $\widehat S-\Pi=-\widehat T$, descomposición de
$\delta$ a precisión de máquina. K1: reducción puntual y $E\le0$ en
la banda. K2: identidad de doble conteo exacta y dependencia de
umbral de $W_L$. K3: signo por $\omega(r)$ y por $u_r$. K4: las tres
sumas por marco. K5: degenerados y factor $\Pi$. K6: niveles
representativos $z=29,37,53$ con caracterización degenerada
verificada y (REG-Π) medido.
