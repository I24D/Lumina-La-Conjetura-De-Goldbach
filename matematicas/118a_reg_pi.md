# Doc 118a — (REG-Π) es falsa: la rebanada ω = 1 es teorema, y la familia ω = 2 con cofactor libre la refuta

**Fecha:** 2026-11-28. **Línea:** UNIF-U / (H2), continuación de los docs
116c y 117b. **Script:** `118a_reg_pi.py` (stdlib; importa la maquinaria
de los docs 110/117b sin tocarla). **Salida:**
`salida_doc118a_local.txt` — chequeos J0–J4, **todos PASS** (7,0 s).
**Etiquetas:** [T] = teorema, [M] = medición, [P] = pendiente.

---

## 0. Veredicto en breve

El doc 117b §3.1 redujo la capa alta de (BAND-MOD) a la FL del modelo
más la suma de región **(REG-Π)**
$\sum_{r\in\mathrm{banda}}\Pi(r)/\varphi(r)\le C\,Q_{\mathrm{completo}}$,
con medición favorable (0,17–0,46 $Q_c$). Este doc la resuelve, y el
resultado es negativo para el enunciado completo:

1. **(REG-Π) es FALSA uniformemente [T]** (§3): a niveles
   $u_L\to2$ de $\mathcal L_A$, la familia $\omega=2$
   $\{pq: x_p\in(11/20,\,3/5),\ 5\le q<p\}$ está entera en la banda
   para $z$ grande y su masa es
   $\ge\log(12/11)\cdot(\log\log z+O(1))\to\infty$ en unidades de
   $M_{\mathrm{full}}$. El cociente a $Q_c$ crece como
   $\frac{\log(12/11)}{M_\infty}\log\log z$: ninguna constante $C$
   sirve. La divergencia es solo $\log\log z$ — por eso la ventana
   $z\le53$ (donde $\log\log z\le1{,}38$) la veía plana.
2. **La rebanada $\omega=1$ SÍ es teorema [T]** (§2), con constante
   explícita: $\limsup \sum_{\omega=1}\Pi/\varphi\,/\,Q_c
   \le \log(10/3)/M_\infty \approx 1{,}036$, por Mertens — es el mismo
   cálculo $B_1$ del doc 116c §2.2 con peso $1/(p-2)$.
3. **Residuo mínimo exacto** (§4): la cadena del doc 117b pierde el
   mayorante (REG-Π) y conserva (BAND-MOD) sobre la capa alta como
   enunciado directo, **(BAND-MOD-HI)** [P]. Ni (BAND-MOD) ni
   (BAND-δ) quedan refutados: el paso que muere es la mayoración
   $(\widehat S-\Pi)^+\le\rho_{\mathrm{FL}}\,\Pi(r)$ sumada contra
   $1/\varphi$.

Goldbach no se declara demostrado; (H2) no queda ni probada ni
refutada por esta sección.

## 1. La identidad y qué pide (REG-Π) en realidad [T]

Por la caracterización del doc 116c §1,
$\Pi(r)=\Pi_{\mathrm{full}}\prod_{p\mid r}\frac{p-1}{p-2}$, luego

$$\frac{\Pi(r)}{\varphi(r)}
\;=\;\Pi_{\mathrm{full}}\prod_{p\mid r}\frac{1}{p-2}
\;=\;\Pi_{\mathrm{full}}\,g(r),
\qquad g(r):=\prod_{p\mid r}\frac1{p-2},$$

y $Q_{\mathrm{completo}}=\Pi_{\mathrm{full}}M_{\mathrm{full}}$ (doc
110), con $M_{\mathrm{full}}(z)\uparrow M_\infty\in
[1{,}16181601,\,1{,}16182763]$ (doc 115). Por tanto

$$\text{(REG-Π)}\quad\Longleftrightarrow\quad
\sum_{r\in\mathrm{banda}} g(r)\;\le\;C\,M_{\mathrm{full}}(z)\;\asymp\;C,$$

una masa de Mertens de la región banda con pesos $g$. Verificación
(J0): identidad a error $\le2{,}8\times10^{-17}$ y
$Q_c=\Pi_{\mathrm{full}}M_{\mathrm{full}}$ a $\le6{,}6\times10^{-16}$
en $z=29,37$.

Dos advertencias de estructura. (i) La condición $\Pi(r)<1/2$ en
forma aditiva es una cota **superior** sobre la masa armónica del
soporte, $\sum_{p\mid r}\log\bigl(1+\frac1{p-2}\bigr)
<\log\frac{1}{2\Pi_{\mathrm{full}}}\sim\log\log z$ (doc 116c §1):
excluye soportes con *mucha* masa armónica, pero la familia de §3
tiene masa armónica $O(1)$ por elemento y pasa entera. (ii) Sin la
región, no hay nada que hacer: $\sum_{\text{todo }r}g(r)
=\prod_{p\in\mathcal P}\bigl(1+\frac1{p-2}\bigr)\asymp\ln z$, es
decir $\sum_{\text{todo}}\Pi/\varphi\,/\,Q_c\asymp\ln z\to\infty$; la
comparación término a término vía $1/\Pi(r)\le1/\Pi_{\mathrm{full}}$
pierde exactamente ese $\ln z$. Toda la ganancia debía venir de la
región — y la región no basta, porque su masa $g$ ya diverge en la
subfamilia de §3.

## 2. La rebanada ω = 1: teorema [T] por Mertens

Sea $u_L=\ln L/\ln z\in[5/3,\,4]$ ($L\in\mathcal L_A$). Para $r=p$
primo: visible $\iff p^3<L\iff x_p<u_L/3$ (doc 116c §2.2), y banda
$\iff u_p=u_L-x_p<3/2\iff x_p>u_L-3/2$ (la condición $\Pi(p)<1/2$
solo puede *quitar* elementos: para la cota superior se omite). La
banda $\omega=1$ es no vacía $\iff u_L<9/4$, y

$$\sum_{\omega=1}\frac{\Pi(p)}{\varphi(p)}
=\Pi_{\mathrm{full}}\!\!\sum_{x_p\in I(u_L)}\frac1{p-2},
\qquad I(u_L)=\bigl(u_L-\tfrac32,\ \tfrac{u_L}3\bigr).$$

Por Mertens ($\sum_{p\le x}\frac1{p-2}=\log\log x+O(1)$, la misma
entrada que usa el doc 116c §2.2), la suma interior tiende a
$B_1(u_L)=\log\frac{u_L/3}{u_L-3/2}$, decreciente en
$[5/3,\,9/4)$ con máximo $B_1(5/3)=\log\frac{10}3=1{,}2040$. Luego

$$\limsup_{z\to\infty}\,
\frac{\sum_{\omega=1,\,\mathrm{banda}}\Pi/\varphi}{Q_c}
\;\le\;\frac{\log(10/3)}{M_\infty}\;\approx\;1{,}036,$$

y (REG-Π) vale sobre $\omega=1$ con cualquier $C>1{,}037$ para $z$
grande. **(REG-Π-ω1) [T].** Verificación (J3): a $z=10^{12}$ las
sumas dan $1{,}1935$ (frente a $1{,}2040$) y $0{,}2876$ (frente a
$\log\frac43=0{,}2877$).

## 3. La refutación: ω = 2 con el cofactor libre [T]

Existen niveles $L\in\mathcal L_A$ con $u_L\to2$ (doc 116c §1).
Fíjense esos niveles y la familia

$$\mathcal F=\bigl\{\,r=pq:\ p,q\in\mathcal P,\ 
x_p\in\bigl(\tfrac{11}{20},\,\tfrac35\bigr),\ 5\le q<p\,\bigr\}.$$

Para $z$ grande, todo $r\in\mathcal F$ está en la banda:

- **Visible:** $\omega(r)=2$; la única condición de Nathanson en
  posición impar es la de $m=1$: $p^3<L$, y
  $x_p<3/5<u_L/3$ (pues $u_L>9/5$); $r<L$ porque
  $x_p+x_q<6/5<u_L$. Luego $r\in V^+(L)$, y $V^+(L)$ cerrado por
  divisores (doc 114, Lema 2) da visibilidad.
- **Banda (tamaño):** $u_r=u_L-x_p-x_q<3/2$ porque
  $x_p+x_q>11/20\ge u_L-3/2$ (pues $u_L\le41/20$).
- **Banda (Π):**
  $\Pi(r)=\Pi_{\mathrm{full}}\frac{p-1}{p-2}\frac{q-1}{q-2}
  \le\frac85\,\Pi_{\mathrm{full}}\to0<1/2$
  ($p\ge7\Rightarrow\frac{p-1}{p-2}\le6/5$;
  $q\ge5\Rightarrow\frac{q-1}{q-2}\le4/3$).

Su masa $g$ factoriza ($p$ es el primo mayor, único en cada $r$):

$$\sum_{r\in\mathcal F} g(r)
=\sum_{x_p\in(11/20,\,3/5)}\frac1{p-2}
\sum_{5\le q<p}\frac1{q-2}
\;\ge\;
\underbrace{\sum_{x_p\in(11/20,\,3/5)}\frac1{p-2}}_{\to\ \log(12/11)}
\cdot
\underbrace{\sum_{5\le q\le z^{1/2}}\frac1{q-2}}_{=\ \log\log z\,+\,O(1)},$$

ambas por Mertens. Por tanto
$\sum_{\mathrm{banda}}g(r)\ge\log(12/11)\,(\log\log z+O(1))\to\infty$
y

$$\frac{\sum_{\mathrm{banda}}\Pi(r)/\varphi(r)}{Q_{\mathrm{completo}}}
=\frac{\sum_{\mathrm{banda}}g(r)}{M_{\mathrm{full}}(z)}
\;\ge\;\frac{\log(12/11)}{M_\infty}\,\log\log z\;(1+o(1))
\;\longrightarrow\;\infty.$$

**Ninguna constante $C$ sirve: (REG-Π) es falsa.** El mecanismo es
la *multiplicidad del cofactor pequeño*: en el doc 116c §2.3 un
cofactor **fijo** $m_0$ aportaba $B_1(u^*)/\varphi(m_0)$, acotado;
aquí el cofactor $q$ recorre el segmento inicial y
$\sum_q\frac1{q-2}\sim\log\log z$. (Las familias $\omega\ge3$
análogas divergen aún más rápido; basta $\omega=2$.) Además
$u_r\to(u_L-3/5,\,u_L-11/20)\to(7/5,\,29/20)\subset(1,\,3/2)$: la
obstrucción vive en la capa alta $\{u_r\ge1\}$ — exactamente donde
el doc 117b §3.1 aplica el mayorante —, no en la capa frontera.

*Por qué la ventana no lo veía [M].* A $z=53$, $\log\log z=1{,}377$
y la constante es pequeña: la parte $\omega=2$ medida a $u_L=2$ es
$0{,}0436\,Q_c$ (J1), y es **exactamente** esta familia: su masa $g$
en banda es $0{,}0505$, idéntica a $0{,}0436\times M_{\mathrm{full}}
(53)$ (J1/J4). A $z=53$ solo 3 de los 5 elementos de la familia con
$x_p\in(1/2,\,2/3)$ están ya en la banda — $\Pi_{\mathrm{full}}(53)
=0{,}3677$ aún excluye a los de factor $\frac{p-1}{p-2}\frac{q-1}{q-2}$
grande (hasta $0{,}588>1/2$) — y la visibilidad es 5/5 (J4). La
entrada completa es asintótica, como en el doc 116c §3.

## 4. Qué sobrevive y el residuo mínimo exacto

- **(REG-Π) queda retirada** como ruta: su rebanada $\omega=1$ es
  [T] (§2) pero no controla la cadena, porque la masa
  $(\widehat S-\Pi)^+$ sobre $\omega\ge2$ es justamente la parte no
  controlada — y §3 muestra que el mayorante $\Pi/\varphi$ diverge
  ya en $\omega=2$.
- **(BAND-MOD) y (BAND-δ) NO quedan refutados.** Lo que muere es la
  mayoración $(\widehat S-\Pi)^+\le\rho_{\mathrm{FL}}(u_0^*)\,
  \Pi(r)$ sumada contra $1/\varphi$: en la familia $\mathcal F$ el
  exceso real del modelo es una cuestión sobre $\widehat S$, no
  sobre $\Pi$ (el doc 117b §3.3 midió
  $(\widehat S-\Pi)^+/\Pi(r)\le0{,}692$ en las instancias no
  degeneradas de la banda; si en $\mathcal F$ el cociente fuese
  $\asymp\rho>0$ con $\rho$ fijo, (BAND-MOD) heredaría la
  divergencia — decidirlo es el contenido del residuo, no de esta
  sección).
- **Residuo mínimo exacto. (BAND-MOD-HI) [P]:**
  $\sum_{r\in\mathrm{banda},\ u_r\ge u_0^*}
  (\widehat S(r)-\Pi(r))^+/\varphi(r)\le C\,Q_{\mathrm{completo}}$,
  a atacar directamente sobre el modelo (signo y tamaño de
  $\widehat S-\Pi$ en la familia de cofactor libre), sin pasar por
  $\sum\Pi/\varphi$.
- Cadena actualizada: (DESP-SIGN) + (BAND-LOW) + (BAND-MOD-HI)
  $\Rightarrow$ (BAND-MOD) $\Rightarrow$ (BAND-δ) $\Rightarrow$
  (BAND-S). Sin circularidad: esta sección no usa (BAND-δ) ni
  (BAND-LOW); solo Mertens, la FL no se usa, y la identidad de §1.

## 5. Medición [M] (J1): la banda por ω, cocientes a $Q_c$

| $z$ | $u_L$ | $\omega=1$ | $\omega=2$ | $\omega\ge3$ | total |
|---|---|---|---|---|---|
| 29 | 5/3 | 0 | 0 | 0 | 0 |
| 29 | 2 | 0 | 0 | 0 | 0 |
| 37 | 5/3 | 0,1733 | 0 | 0 | 0,1733 |
| 37 | 2 | 0,2695 | 0 | 0 | 0,2695 |
| 53 | 5/3 | 0,4609 | 0 | 0 | 0,4609 |
| 53 | 2 | 0,1746 | 0,0436 | 0 | 0,2182 |

(A $z=29$ la banda es vacía en estos niveles: $\Pi_{\mathrm{full}}$
aún grande frente a la visibilidad, cf. doc 117b §4.) Los totales
reproducen el rango 0,17–0,46 del doc 117b y son compatibles con la
refutación: a $z\le53$ el factor $\log\log z$ vale $\le1{,}38$ y la
constante de la familia es $\log(12/11)/M_\infty\approx0{,}075$.

## 6. Reproducibilidad

`python3 118a_reg_pi.py | tee salida_doc118a_local.txt` — stdlib;
importa `117b_banda_delta.py` (que importa la maquinaria del doc
110) sin modificarlos. J0: identidad $\Pi/\varphi=\Pi_{\mathrm{full}}
g$ y $Q_c=\Pi_{\mathrm{full}}M_{\mathrm{full}}$. J1: banda por
$\omega$ en $z=29,37,53$. J2: sumas de Mertens de la familia
obstructora a $z=10^4,\dots,10^{10}$ (exterior
$0{,}0850\to0{,}0870=\log(12/11)$; interior
$1{,}2262\to2{,}1324\approx\log\log z-1$). J3: rebanada $\omega=1$
a $z=10^{10},10^{12}$ contra $B_1$. J4: la familia contra
`visible_set`, $\Pi$ y $u_r$ reales a $z=53$. **PASS global.**
