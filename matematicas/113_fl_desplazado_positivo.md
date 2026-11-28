# Doc 113 — (FL-DESP⁺): el lema fundamental desplazado en su parte positiva

**Misión Goldbach — LUMINA (Muse), 2026-11-28.**
Script: `113_fl_desplazado_positivo.py` (este directorio; importa
`Level111`/`Level110` de los docs 111/110 por la cadena de imports:
literalmente las mismas definiciones de los docs 109–112).
Salida: `salida_doc113_local.txt` — **PASS** (K0–K7).
Solo stdlib; modo fijo (los 240 niveles de los marcos A/B, sin
muestreo).

**Pregunta del doc.** El doc 112 cerró el bloque $r=1$ de la masa
positiva [T] con el lema fundamental (FL) clásico y redujo (H2) a

$$
(H2⁗):\qquad
P_{>1}^{(L)}=\sum_{\substack{r>1\ \mathrm{visible}\\ \delta(r)>0}}t(r)
\;\le\;C\cdot Q_{\mathrm{completo}},\qquad
t(r)=\frac{\delta(r)\,(2\Pi(r)+\delta(r))}{\varphi(r)},
$$

con $\delta(r)=S(r)-\Pi(r)$,
$S(r)=\sum_{d'\in V^{(r)}}a(d')$ sobre el soporte desplazado
$V^{(r)}=\{d':(d',r)=1,\ rd'\in V^{+}(L)\}$ (defs. docs 109–112).
La hipótesis del doc: el nivel correcto para el desplazado no es el
$u$ del marco (donde el doc 109 §3.4 midió
$\max|\delta|/\Pi=1{,}00$–$1{,}30$) sino el **nivel propio**
$u_r:=\ln(L/r)/\ln z$ del soporte $V^{(r)}$, y con
$\rho(u_r)$ la parte positiva puntual podría sostenerse. Este doc
formula (FL-DESP⁺), mide la cota al nivel propio en los 240 niveles,
disecciona los ofensores y lleva el intento de prueba hasta su paso
exacto de ruptura. **Veredicto anticipado: la cota se sostiene en
ventana con constante pequeña $A=1{,}43$ — no con $A=1$ —; la
constante es la de la FL clásica a nivel pequeño, NO un efecto del
desplazamiento: el error de desplazamiento medido es
$E(r)\le0$ en las 4113 instancias (máx $E/\Pi=6{,}4\times10^{-16}$,
cero a precisión de máquina). La prueba no cierra: el residuo
mínimo es el signo del desplazamiento, (DESP-SIGN). (H2) sigue sin
demostrarse; Goldbach sigue abierta.**

Etiquetas: **[T]** teorema/identidad demostrada aquí o cita clásica
con hipótesis verificadas · **[M]** medición del script ·
**[P]** pendiente.

---

## 1. (FL-DESP⁺): formulación precisa y su implicación [T]

Sea $\rho_0(u):=u^{-u/2}=e^{-u\ln u/2}$ la forma del error de la FL
clásica (doc 112 §2: Halberstam–Richert Cap. 3; Tenenbaum Cap. 6,
error $O(u^{-u/2})$, dimensión $\kappa=1$). Nótese que $\rho_0$ no es
monótona en $(0,\infty)$: $\rho_0(1)=1$, $\rho_0\to1^{+}$ cuando
$u\to0^{+}$, máximo $e^{1/(2e)}\approx1{,}202$ en $u=1/e$, y decae
superexponencialmente para $u$ grande.

$$
(\mathrm{FL\!-\!DESP}^{+}):\qquad
\delta(r)\;\le\;\Pi(r)\cdot\rho(u_r)
\quad\text{para todo } r>1 \text{ visible, en todo nivel } L
\text{ del marco,}
$$

con $\rho$ una función de tipo FL ($\rho(u)\to0$ cuando
$u\to\infty$) y $u_r=\ln(L/r)/\ln z$, bajo la condición de nivel del
marco (doc 109 §2: $u=\ln L/\ln z\ge u_{\min}$, $D_A\simeq N^{1/2}$,
$y_A=N^{1/3}$). Solo la cara positiva: los $r$ invisibles tienen
$\delta(r)=-\Pi(r)<0$ (doc 109 §3.1 [T]) y no cuentan; el canal
$\delta<-2\Pi$ (i.e. $S<-\Pi$) está vacío en ventana [M, K2: 0 de
4113], así que la cara positiva es la única que alimenta $t(r)>0$
vía $\delta>0$.

**Implicación (FL-DESP⁺) ⇒ (H2⁗) [T].** Si $0\le\delta(r)\le\rho_r
\Pi(r)$ con $\rho_r:=\rho(u_r)$, la monotonía de
$x\mapsto x(2\Pi+x)$ en $x\ge0$ da
$t(r)\le\rho_r(2+\rho_r)\,\Pi(r)^2/\varphi(r)$ (verificado sin
violación en las 2892 instancias de ventana con
$0\le\delta\le\rho_0\Pi$ [M, K1]). Sumando y usando la identidad del
doc 110 §1,
$\sum_{r>1}\Pi(r)^2/\varphi(r)=Q_{\mathrm{completo}}-\Pi_{\mathrm
{full}}^2$ [T, K1: error $0$], con
$\Pi(r)^2/\varphi(r)=\Pi_{\mathrm{full}}^2\,h(r)$,
$h(r)=\prod_{p\mid r}(p-1)/(p-2)^2$, y
$\sum_r h(r)=M_{\mathrm{full}}/\Pi_{\mathrm{full}}=2{,}2118$ /
$2{,}6422$ / $2{,}8415$ para $z=17/29/37$ [T, K1]:

$$
P_{>1}^{(L)}\;\le\;\bar\rho\,(2+\bar\rho)\cdot
\big(Q_{\mathrm{completo}}-\Pi_{\mathrm{full}}^2\big)
\;\le\;\bar\rho\,(2+\bar\rho)\cdot Q_{\mathrm{completo}},
\qquad \bar\rho:=\sup_u\rho(u).
$$

Con $\rho=A\rho_0$: $\bar\rho=1{,}202\,A$ y
$C=A\cdot1{,}202\,(2+1{,}202A)$; para el $A$ medido en §3,
$C\approx6{,}4$. (Si $\rho\le1$ en el rango relevante, $C\le3$.)
La implicación es álgebra exacta; todo el contenido está en
(FL-DESP⁺) misma.

## 2. La descomposición exacta: modelo RI + error de desplazamiento [T]

El soporte desplazado difiere de un soporte RI genuino en dos cosas
(tarea (3)): **(a)** la exclusión de los primos de $r$ y **(b)** la
condición de Nathanson aplicada al producto $rd'$ (los primos de
$r$ entran en la factorización decreciente fusionada). La pieza (a)
es inofensiva y exacta: defínase el **soporte modelo**

$$
\widehat V^{(r)}:=V^{+}\big(L/r;\ \mathcal P\setminus\mathrm{primos}(r)\big),
\qquad
\widehat S(r):=\sum_{d'\in\widehat V^{(r)}}a(d'),
$$

un soporte de Rosser–Iwaniec genuino de nivel $L/r$ sobre los primos
restantes. Su término principal es **exactamente** $\Pi(r)$:
$\prod_{p\in\mathcal P\setminus r}(1-\frac1{p-1})=\Pi(r)$ por
definición (doc 109 §1). La FL clásica del doc 112 §2 se aplica a
$\widehat S(r)$ con nivel $u_r$ y el mismo $\kappa=1$ (Mertens sobre
$\mathcal P\setminus r$ difiere de Mertens sobre $\mathcal P$ en
$O(1)$ factores locales ya absorbidos en $\Pi(r)$):
$|\widehat S(r)-\Pi(r)|\le\rho_{\mathrm{FL}}(u_r)\Pi(r)$ con la
$\rho_{\mathrm{FL}}$ y la constante absoluta de la FL citada.
Control del código: $\widehat S(1)=S(1)$ con error $0{,}00$ en 12
niveles [M, K0] ($\widehat V^{(1)}=V^{+}(L)=V$).

La pieza (b) es el **error de desplazamiento**

$$
E(r):=S(r)-\widehat S(r),\qquad\text{i.e.}\qquad
\delta(r)=\big(\widehat S(r)-\Pi(r)\big)+E(r)\quad[\mathrm T],
$$

la masa con signo de la diferencia simétrica
$V^{(r)}\mathbin\triangle\widehat V^{(r)}$: en $V^{(r)}$ la
condición sobre $d'$ se lee en la factorización **fusionada** de
$rd'$; en $\widehat V^{(r)}$, en la de $d'$ sola, escalada por $r$
(arimética entera en el script: $p_m^2(p_1\cdots p_m)\ge L/r\iff
r\,p_m^2(p_1\cdots p_m)\ge L$). Por tanto

$$
(\mathrm{FL\!-\!DESP}^{+})\ \Longleftarrow\
\text{FL clásica sobre }\widehat S\ +\ E(r)\le0
\ \text{(o }E(r)\le\Pi(r)\eta(u_r)\text{)}.
$$

## 3. Medición decisiva: la cota al nivel propio [M, K2]

Ratio $\delta(r)^{+}/(\Pi(r)\rho_0(u_r))$ sobre las **4113**
instancias $(L,r)$ con $r>1$ visible (240 niveles); 3085 con
$\delta>0$:

- $A_{\min}=\max\mathrm{ratio}=\mathbf{1{,}4273}$: con $A=1$ la
  cota puntual positiva al nivel propio **NO** se sostiene; con
  constante pequeña $A=1{,}43$, **SÍ**.
- Restringido a **sumas genuinas** ($S(r)<1$):
  $A_{\min}^{\mathrm{gen}}=\mathbf{0{,}8844}\le1$. Toda violación
  de $A=1$ viene de sumas degeneradas (§4).
- Distribución (sobre $\delta>0$): mediana $0{,}4516$, p90
  $0{,}7898$, p99 $1{,}3160$; histograma: 2883 instancias con ratio
  $<1$, 202 en $[1,1{,}5)$, ninguna $\ge1{,}5$.
- Por marco: $z=17$: $0{,}8293$; $z=29$ (A): $1{,}2523$; $z=37$ (A):
  $1{,}4273$; $z=29$ (B): $1{,}3809$. Crece con $z$ en ventana — el
  mecanismo es el umbral finito de §4, no un fenómeno a $u_r$ fijo.
- Bandas por $u_r$ (máx ratio): $[0{,}5,1)$: $1{,}1718$;
  $[1,2)$: $1{,}4273$; $[2,3)$: $0{,}6640$; $[3,5)$: $0{,}0524$.
  Sin instancias con $\delta>0$ para $u_r<0{,}5$ ni $u_r\ge5$.
- Nivel a nivel (K5): la cota sumada implicada
  $B_A(L)=\sum_{r>1}\rho_A(2+\rho_A)\Pi^2/\varphi$ domina a
  $P_{>1}(L)$ en los **240 niveles** tanto con $A=A_{\min}$
  ($\max P_{>1}/B_A=0{,}487$; $\max B_A/Q_c=1{,}360$) como con
  $A=1$ ($\max P_{>1}/B_A=0{,}786$): puntualmente $A=1$ falla,
  pero los ofensores pesan tan poco ($\Pi^2/\varphi$) que la suma
  pasa incluso con $A=1$ en ventana.

**Los ofensores del doc 109, diseccionados [M, K3].**
$\max|\delta|/\Pi$ del doc 109 se reproduce exacto
($1{,}000/1{,}069/1{,}295$):

| marco | ofensor doc 109 | visible | signo $\delta$ | $u_r$ | ratio |
|---|---|---|---|---|---|
| A $z=17$ | $1{,}000$ en $L=112$, $r=5$ | **no** | $-$ ($\delta=-\Pi$) | $1{,}097$ | — (no cuenta) |
| A $z=29$ | $1{,}069$ en $L=12378$, $r=437$ | sí | $+$ | $0{,}993$ | $1{,}065$ |
| A $z=37$ | $1{,}295$ en $L=39557$, $r=899$ | sí | $+$ | $1{,}048$ | $1{,}327$ |
| B $z=29$ | $1{,}069$ en $L=24399$, $r=437$ | sí | $+$ | $1{,}195$ | $1{,}188$ |

El $1{,}000$ de $z=17$ es el invisible estructural del doc 109
§3.4 ($\delta=-\Pi<0$): no alimenta la parte positiva. El máximo
visible positivo de $z=17$ es $0{,}600$ ($r=143$, $u_r=0{,}991$,
ratio $0{,}597$). Los visibles grandes tienen $u_r\approx1$, donde
$\rho_0(u_r)\approx1$: el nivel propio explica su tamaño pero no
los reduce bajo $A=1$ — §4 explica por qué.

## 4. El mecanismo de los ofensores: soportes degenerados [T+M]

**Los 10 ratios máximos de K2 son, todos, soportes degenerados
$V^{(r)}=\{1\}$** [M, K7]: $S(r)=1{,}000000$ exacto,
$\widehat S(r)=1$ exacto, $E(r)=0$ exacto. Hay 1348 instancias
degeneradas de 4113 (1347 con $\delta>0$), y el máximo ratio entre
ellas es $1{,}4273=A_{\min}$: **toda la constante de (FL-DESP⁺) en
ventana la fijan los soportes de un solo elemento.** Para ellos

$$
\delta(r)=1-\Pi(r)\quad[\mathrm T],\qquad
\frac{\delta(r)}{\Pi(r)}=\Pi(r)^{-1}-1,
$$

y (FL-DESP⁺) se reduce a $1-\Pi(r)\le\Pi(r)\rho(u_r)$: una
afirmación sobre la constante de la FL **clásica** a nivel
$u_r\approx1{,}2$–$1{,}4$ sobre una suma RI de un término — de
hecho la constante empírica de la propia FL sobre el modelo es
$C_{\mathrm{model}}^{+}=1{,}4729$ [M, K4], alcanzada también en un
soporte modelo degenerado ($\widehat S=1$: [A] $z=37$, $L=1358$,
$r=11$ primo, $u_r=1{,}334$; ahí $S=0{,}5833$ y $E/\Pi=-0{,}9231$).
Conclusión de la sección: $\rho_0$ con constante $1$ es falsa a
nivel $u\approx1{,}3$ **ya para sumas RI genuinas**; la constante
$1{,}43$–$1{,}47$ es la constante finita de la FL clásica a nivel
pequeño. El desplazamiento no interviene: en los ofensores,
$E=0$.

**Por qué los degenerados viven a $u_r$ pequeño [T].**
$V^{(r)}=\{1\}$ exige que todo $d'>1$ coprimo falle alguna
condición fusionada impar. Para el primo más pequeño disponible
$p_0=\min(\mathcal P\setminus r)$:

- Si $\omega(r)$ es par y $5\nmid r$ ($p_0=5$): $5$ cae en la
  posición fusionada $\omega(r)+1$, impar, y su condición es
  $5^2\cdot(r\cdot5)=125\,r<L$; si $V^{(r)}=\{1\}$ debe fallar:
  $u_r\le\ln125/\ln z=3\ln5/\ln z$ **[T]**. Verificado en las 363
  instancias de esta familia: $\max(u_r-3\ln5/\ln z)=-0{,}0019$
  [M]. (Los top-5 degenerados: $u_r=1{,}331$ contra umbral
  $1{,}337$ en $z=37$; $1{,}431$ contra $1{,}434$ en $z=29$ B.)
- Si $\omega(r)$ es impar y $5\nmid r$: $5$ cae en posición par
  (sin condición) y las posiciones impares de la fusión son las
  del propio $r$, que se cumplen por $r\in V$; luego
  $5\in V^{(r)}$ en cuanto $5r<L$ **[T]**: no puede haber
  degenerados de esta familia — medidos: **0** [M]. (Conteos por
  (paridad de $\omega$, $5\mid r$): (impar, sí) 513, (par, sí)
  471, (par, no) 363, (impar, no) 0.)

Corolario de lectura: a $u_r\ge c>0$ **fijo**, los degenerados de
estas familias desaparecen cuando $z\to\infty$ (el umbral
$3\ln p_0/\ln z\to0$); el crecimiento de $A_{\min}$ con $z$ en
ventana viene de que $\Pi(r)^{-1}-1\sim\ln z/C_{\mathrm{loc}}$
crece sobre degenerados cuyo $u_r$ se arrima al umbral. No es
evidencia de que $A_{\min}\to\infty$ a nivel propio fijo; la
pregunta asintótica correcta es la constante sobre sumas genuinas,
y ahí la medición da $\le0{,}8844$ con el máximo en la banda
$[1,2)$ y decaimiento rápido ($0{,}052$ en $[3,5)$).

## 5. El error de desplazamiento: una sola cara en ventana [M, K4/K6]

Sobre las 4113 instancias, con $\widehat S$ calculado exactamente:

$$
\max_r \frac{E(r)}{\Pi(r)}=6{,}4\times10^{-16},\qquad
\min_r \frac{E(r)}{\Pi(r)}=-0{,}9231,\qquad
\max_r \frac{E(r)^{+}}{\Pi(r)\rho_0(u_r)}=5{,}8\times10^{-15}.
$$

**$E(r)\le0$ en ventana, a precisión de máquina**: la fracción de
instancias con $E>0$ es $2{,}07\%$, todas con
$E/\Pi\le6{,}4\times10^{-16}$ (ruido de punto flotante). Es decir
$S(r)\le\widehat S(r)$ siempre: la condición de Nathanson sobre el
producto nunca deja pasar, en masa con signo, más que el soporte RI
de nivel $L/r$. En consecuencia la constante de (FL-DESP⁺) en
ventana ($1{,}4273$) es $\le$ la de la FL clásica sobre el modelo
($C_{\mathrm{model}}^{+}=1{,}4729$, dos caras $1{,}4729$): el
desplazado hereda la constante del modelo y el desplazamiento solo
resta.

Y sin embargo la diferencia simétrica es **bilateral en conteo**
[M, K6]: 6509 elementos $d'$ en $V^{(r)}\setminus\widehat V^{(r)}$
contra 1049 en $\widehat V^{(r)}\setminus V^{(r)}$ en total, con
masa absoluta hasta $0{,}923\,\Pi(r)$; ejemplos explícitos:
$d'=7\in V^{(5)}\setminus\widehat V^{(5)}$ ([A] $z=17$, $L=344$:
la fusión $(7,5)$ pasa — $7^2\cdot7=343<344$ — pero la condición
propia escalada $5\cdot7^3\ge L$ falla) y
$d'=35\in\widehat V^{(11)}\setminus V^{(11)}$ ([A] $z=17$,
$L=5239$). **Ningún argumento de monotonía o inclusión signa
$E(r)$**: el signo medido es un hecho de la masa con signo
$a(d')=\mu(d')/\varphi(d')$, no de los soportes.

## 6. Intento de prueba: el paso exacto que se rompe

**(a) Lo que sí cierra [T].** La exclusión de primos está
contabilizada exactamente en $\Pi(r)$ (§2) y la FL clásica cierra
$\widehat S(r)$ al nivel propio $u_r$ con la cita del doc 112 §2
— la misma verificación de hipótesis (Mertens $\Rightarrow$
$\kappa=1$; nivel $L/r$; conjunto de cribado
$\mathcal P\setminus r$). Por §2, (FL-DESP⁺) queda reducida a
signar el desplazamiento.

**(b) La recursión de Buchstab sobre $V^{(r)}$: dónde se rompe
[T].** La prueba clásica de la FL (la construcción del mayorante
$\lambda^{+}$ de Halberstam–Richert, o la recursión de Buchstab
sobre la que descansa) hereda sobre sub-soportes de la **misma**
familia al pelar el primo más pequeño $p_0$ de $d'$. Para el
desplazado, la rama divisible sí es hereditaria, por asociatividad
de la fusión:
$d'=p_0d''\in V^{(r)}\iff d''\in V^{(rp_0)}$ [T, es la identidad
$rd'=rp_0d''\in V$]. Pero la **rama complementaria**
$\{d'\in V^{(r)}:p_0\nmid d'\}=\{d':(d',rp_0)=1,\ rd'\in V\}$ no
es ningún $V^{(r')}$: es un objeto **desplazamiento + exclusión
extra**. La familia de un parámetro $\{V^{(r)}\}$ no es cerrada
bajo Buchstab; ése es el paso exacto que se rompe. (En el soporte
RI la rama complementaria es $V^{+}(x;\mathcal P)$ con el primo
excluido del lado del soporte, y la recursión clásica la absorbe
porque las condiciones de Nathanson son intrínsecas a $d'$; en el
desplazado dependen de la fusión con $r$.)

**(c) La familia ampliada no rescata la inducción [T].**
Ampliando a pares (desplazamiento $r$, exclusiones $E\supseteq
\mathrm{primos}(r)$), $V^{(r,E)}=\{d':(d',E)=1,\ rd'\in V\}$, la
familia sí es formalmente cerrada bajo Buchstab en ambas ramas.
Pero la comparación $S(r,E)\le\widehat S(\cdot,E)$ no se induce:
al pelar $p_0$, la rama divisible entra multiplicada por
$a(p_0)=-1/(p_0-1)<0$, que **invierte** la desigualdad en esa
rama, y las ramas se mezclan; controlar la inducción exigiría la
FL **a dos caras** sobre toda la familia ampliada — es decir, el
enunciado completo otra vez. Además la forma $u^{-u/2}$ de la FL
clásica no nace de una inducción desnuda sobre soportes: nace de
la construcción del mayorante $\lambda^{+}$ apoyado en
$V^{+}(D_A)$ con suavizado sobre niveles, y no existe en la
literatura un mayorante desplazado análogo. Declarar probada
(FL-DESP⁺) por esta vía sería circular o falso: no se hace.

**(d) Residuo mínimo.** Lo que falta, exactamente, es

$$
(\mathrm{DESP\!-\!SIGN}):\qquad
E(r)=S(r)-\widehat S(r)\;\le\;0
\qquad\text{para todo } r>1 \text{ visible,}
$$

o su forma debilitada suficiente
$E(r)\le\Pi(r)\,\eta(u_r)$ con $\eta$ de tipo FL. Con
(DESP-SIGN), (FL-DESP⁺) sigue de la FL clásica sobre el modelo
(§2, §6a) con $\rho=\rho_{\mathrm{FL}}$; y con (FL-DESP⁺),
(H2⁗) sigue por §1. (DESP-SIGN) es **falso de probar por
inclusión** (§5: la diferencia es bilateral) y **verdadero en
ventana a precisión de máquina** (§5): es un lema de masa con
signo sobre la diferencia de dos condiciones de Nathanson —
genuinamente el mismo tipo de muro que (DIST-$\mathcal D$) en la
otra línea (doc 107), aunque de naturaleza combinatoria, no
distribucional.

## 7. Qué separa la medición de la prueba

1. **(DESP-SIGN) sin prueba** (§6b–d): uniformidad en $r$ y en el
   marco. La medición cubre $z\le37$ y los marcos A/B completos;
   el signo es exacto en ventana pero el argumento combinatorio
   que lo fuerce (¿una inyección con signo entre
   $V^{(r)}\setminus\widehat V^{(r)}$ y
   $\widehat V^{(r)}\setminus V^{(r)}$ que respete $a(d')$?) no
   existe todavía: la masa absoluta llega a $0{,}92\,\Pi$, así
   que la cancelación es casi ajustada y una inyección grosera
   no basta.
2. **La banda $u_r<u_0$**, donde la FL clásica es vacua (su
   $O(u^{-u/2})$ es asintótico en $u$): contiene a los
   degenerados de §4, para los que $\delta=1-\Pi(r)$ es [T] pero
   $\rho(u_r)\approx1$ no da ganancia. En ventana la **suma**
   sobrevive (K5: pasa con $A=1$ a nivel de $B_A(L)$) porque
   $\Pi(r)^2/\varphi(r)$ es pequeño ahí; convertir eso en prueba
   exige la sumabilidad trivial de la banda dentro del marco
   auditado ($\delta=O(\log z)$ puntual, $\varphi(r)$ grande,
   soporte polilogarítmico) — trabajo estándar **no desarrollado**
   en el marco, [P]. No es un muro conceptual; es deuda de
   escritura con contenido.
3. **La constante de la FL a nivel pequeño**: (FL-DESP⁺) debe
   formularse con la $\rho_{\mathrm{FL}}$ real (constante absoluta
   de Halberstam–Richert/Tenenbaum), no con la forma desnuda
   $\rho_0$ de constante $1$, que es falsa a $u\approx1{,}3$
   incluso para sumas RI genuinas [M: $C_{\mathrm{model}}^{+}=
   1{,}4729$]. Con esa corrección de formulación, la medición es
   compatible con (FL-DESP⁺) en todo lo medido.

## 8. Estado de la cadena tras el doc 113

$$
(\mathrm{H2})\ \xLeftarrow[\text{doc }109]{}\ (\mathrm{H2}')
\ \xLeftrightarrow[\text{doc }110]{}\ (\mathrm{H2}'')
\ \xLeftrightarrow[\text{doc }111]{}\ (\mathrm{H2}''')
\ \xLeftarrow[\text{doc }112]{}\ (H2⁗)
\ \xLeftarrow[\S1\ \mathrm T]{}\ (\mathrm{FL\!-\!DESP}^{+})
\ \xLeftarrow[\S2,\S6a\ \mathrm T]{}\ \text{FL clásica (modelo)}
\ +\ (\mathrm{DESP\!-\!SIGN})\ +\ \text{banda }u_r<u_0\ [\mathrm P].
$$

(H2) **no** queda demostrada: queda reducida, en esta línea, a un
único lema de signo combinatorio (DESP-SIGN) más la deuda de la
banda pequeña. Es un residuo más pequeño y de naturaleza distinta
que el (H2⁗) del doc 112 — el desplazamiento, que era el
sospechoso, resulta inocente en ventana ($E\le0$, y $E=0$ en los
ofensores) — pero sigue siendo un lema nuevo sin prueba. La
cadena condicional del doc 108 no cambia en su conclusión:
(UNIF-U-COND) sigue requiriendo además (DIST-$\mathcal D^{+}$).
**La conjetura fuerte de Goldbach sigue abierta.**

*Verificación: K0 (cruce doc 112: $0{,}224/0{,}330/0{,}381/0{,}130$;
$\widehat S(1)=S(1)$ exacto) · K1 (álgebra de la implicación [T],
errores $\le3\times10^{-17}$; monotonía en 2892 instancias sin
violación) · K2 ($A_{\min}=1{,}4273$; genuinas $0{,}8844$; canal
$\delta<-2\Pi$ vacío) · K3 (ofensores del doc 109 reproducidos y
diseccionados) · K4 ($E\le0$ a precisión de máquina;
$C_{\mathrm{model}}^{+}=1{,}4729$) · K5 (implicación nivel a nivel
en 240/240 con $A_{\min}$ y con $A=1$) · K6 (diferencia simétrica
bilateral: 6509 contra 1049, ejemplos explícitos) · K7 (1348
degenerados; top-10 todos degenerados; umbral $3\ln5/\ln z$
verificado, holgura $-0{,}0019$; familia (impar, $5\nmid r$): 0).
Salida completa en `salida_doc113_local.txt`: **PASS**.*
