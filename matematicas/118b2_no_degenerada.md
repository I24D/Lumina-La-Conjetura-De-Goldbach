# Doc 118b2 — La parte no degenerada de (BAND-LOW) y su reconciliación con la refutación de (REG-Π)

**Fecha:** 2026-11-28. **Línea:** UNIF-U / (H2), continuación de los docs 117b y 118a.
**Script:** `118b2_no_degenerada.py` (stdlib). **Salida:** `salida_doc118b2_local.txt`.
**Etiquetas:** [T] = teorema, [M] = medición en ventana, [P] = pendiente.

---

## 0. Veredicto en breve

El doc 117b partió (BAND-MOD) en la capa alta $\{u_r\ge u_0^*\}$ —
cerrada allí módulo (REG-Π) por la FL del modelo — y la capa baja
(BAND-LOW), $\{u_r<u_0^*\}$. El doc 118a refutó (REG-Π) sobre toda la
banda [T] con la familia $\omega=2$
$\mathcal F=\{pq: x_p\in(11/20,\,3/5),\ 5\le q<p\}$ a niveles
$u_L\to2$. Esta sección responde: ¿vive esa obstrucción también en
la capa baja, y qué queda de la parte **no degenerada** de
(BAND-LOW)?

1. **Reconciliación (§1).** La masa divergente de $\mathcal F$ se
   concentra en $u_r\in(7/5,\,29/20)$ [J0]: $\mathcal F$ es una
   obstrucción de **capa alta** para todo corte $u_0^*\le7/5$ (el
   propio doc 118a la ubica en $\{u_r\ge1\}$). **Pero el mecanismo
   no es intrínsecamente de capa alta:** (REG-Π) restringida a la
   capa baja, **(REG-Π-LOW), es FALSA para todo $u_0^*\in(0,\,3/2)$
   [T]** — familias hermanas (una $\omega=3$ y una escalera de
   $\omega$ creciente) producen la misma divergencia
   $\sim\log\log z$ a **toda** profundidad $u_r<u_0^*$. La ruta del
   mayorante muere también en la capa baja.
2. **Residuo mínimo (§2).** La parte no degenerada de (BAND-LOW)
   queda como enunciado directo sobre $\widehat S$,
   **(BAND-LOW-NONDEG) [P]**, en paralelo exacto al (BAND-MOD-HI)
   del doc 118a. Ni siquiera una cota puntual
   $(\widehat S-\Pi)^+\le C_1\Pi$ en la capa la cerraría: su suma
   contra $1/\varphi$ es la que diverge.
3. **Medición (§3) [M].** En ventana, la banda no degenerada vive
   entera en $u_r\ge1$ (la celda $u_r<1$ no degenerada está vacía
   en los marcos muestreados); la familia $\mathcal F$ aún no
   aparece en su ventana estricta a $z\le53$, y en la ventana
   ampliada del J4 de 118a sus 3 miembros en banda son
   **degenerados** y están en $u_r\in[0{,}750,\,0{,}906]$: su hogar
   asintótico $(7/5,\,29/20)$ y su no-degeneración son fenómenos de
   $z\to\infty$, cuantificados en J0.

(BAND-MOD) no queda refutado; (H2) y Goldbach, intactos.

## 1. Reconciliación: la obstrucción de 118a y la capa baja

### 1.1 El valor de $u_0^*$ en el doc 117b

El doc 117b **no fija $u_0^*$ numéricamente**: es «$u_0^*>0$
absoluto, el umbral de la FL citada» (117b §3.1; la FL de
Halberstam–Richert Cap. 3 / Tenenbaum, $\rho(u)=O(u^{-u/2})$, docs
112 §2 y 113 §2). La arquitectura solo usa que la FL está en vigor
sobre $\{u_r\ge u_0^*\}$ con $\rho_{\mathrm{FL}}(u_0^*)$ constante
absoluta, y que $u_0^*<3/2$ (la capa baja es subcapa propia de la
banda $u_r<3/2$). La reconciliación se hace, por tanto,
**paramétricamente en $u_0^*$**, y la respuesta es uniforme.

### 1.2 Dónde vive $\mathcal F$ [T, J0]

A niveles $u_L\to2$, la masa $g$ de $\mathcal F$ factoriza (doc
118a §3): $\sum_{x_p\in(11/20,3/5)}1/(p-2)\to\log(12/11)$ [J0(a):
$0{,}0869$/$0{,}0868$ a $z=10^6$/$10^8$, frente a $0{,}0870$] por la
suma interna $\sum_{5\le q<p}1/(q-2)\sim\log\log z$. La parte de la
suma interna con $x_q\ge a>0$ fijo es $O(1)$ (vale
$\log\frac{3/5}{a}+O(1)$); **toda la divergencia viene de
$x_q\to0$**, es decir de $u_r=u_L-x_p-x_q\to u_L-x_p\in(7/5,\,
29/20)$. J0(b) lo cuantifica: la fracción de la masa interna con
$u_r\ge1$ (i.e. $x_q\le1-x_p$) crece $0{,}718\to0{,}781\to0{,}808$
para $z=10^4,10^6,10^8$ (tiende a $1$ como $1-O(1/\log\log z)$).
Luego $\mathcal F$ obstruye la capa $\{u_r\ge u_0^*\}$ para todo
$u_0^*\le7/5$, y no toca $\{u_r<u_0^*\}$ si $u_0^*\le7/5$.


### 1.3 (REG-Π-LOW): falsa para todo corte [T]

**(REG-Π-LOW)** := $\sum_{r\in\mathrm{banda},\,u_r<u_0^*}
\Pi(r)/\varphi(r)\le C\,Q_{\mathrm{completo}}$, equivalente a
$\sum g(r)\le C\,M_{\mathrm{full}}$ sobre la capa (identidad del doc
118a §1). Basta refutarla en los niveles $u_L\to2$ de
$\mathcal L_A$. El mecanismo de 118a — *un primo pequeño libre
multiplica una masa de Mertens positiva* — opera a toda
profundidad de la capa:

- **Hermana $\omega=3$** (cubre $u_0^*>2/3$): la subcaja
  $x_1\in(0{,}56,\,0{,}65)$, $x_2\in(0{,}50,\,0{,}55)$ cumple
  $x_2<x_1<2/3=u_L/3$ (visibilidad, condición de Nathanson en
  $m=1$) y $x_1+x_2>1{,}06>u_L-u_0^*$ si $u_0^*>0{,}94$ — y con
  $x_1+x_2>u_L-u_0^*$ en general si $u_0^*>2/3$, ajustando la caja
  —; sus masas de Mertens por coordenada son $>0$ [J0(c): $0{,}166$
  y $0{,}282$ a $z=10^6$ en la subcaja impresa]. Con el tercer
  primo $q\to$ pequeño libre, la familia
  $\{p_1p_2q\}$ está en la banda (Π$(r)\to0$), es visible (la
  condición en $m=3$ es $x_1+x_2+3x_q<u_L$, cierta pues
  $x_1+x_2<1{,}2<2$), tiene $u_r\to u_L-x_1-x_2\in(2/3,\,u_0^*)$ y
  masa $g\ge(\text{masa de la caja})\cdot(\log\log z+O(1))\to\infty$.
- **Escalera** (cubre todo $u_0^*>0$): primos por parejas en
  $x\approx 2/3^j$ ($j=1,\dots,k$), más un primo pequeño libre. Las
  condiciones de Nathanson en $m=2j-1$ son
  $\sum_{i<m}x_i+3x_m<u_L$ y en el centro valen $=u_L$ exactamente
  (telescopio: $2(1-3^{-j})+2\cdot3^{-j}=2$), luego se cumplen
  estrictas acercándose desde dentro; la suma de coordenadas
  fijas es $2(1-3^{-k})$ y $u_r\to2\cdot3^{-k}$, que es $<u_0^*$
  para $k$ grande. Cada coordenada aporta masa de Mertens $>0$
  (intervalo fijo en escala $x$) y el primo pequeño aporta
  $\log\log z$. [J0(d) verifica el caso $k=2$ ($\omega=5$) con
  $u_0^*=1/2$: suma $1{,}776>3/2$, $u_r=0{,}224<1/2$, Nathanson
  $m=1,3$ estrictas; la condición en $m=5$ es la propia cota de
  la suma.] Variando el centro se obtiene divergencia con
  $u_r$ en cualquier subintervalo de $(0,\,u_0^*)$.

**Conclusión [T]:** (REG-Π-LOW) es falsa para todo
$u_0^*\in(0,\,3/2)$. La obstrucción de 118a no es intrínsecamente
de capa alta: solo lo es su familia concreta $\mathcal F$.

*Corrección a 117b §3.2.* La capa $\{u_r<u_0^*\}$ a nivel fijo
$u_L$ es la región $\{u_L-u_0^*<\sum x_i\}$ de la región banda:
una **loseta de medida de Mertens positiva**, no una frontera —
la dicotomía del doc 116c §4 predice entonces divergencia del
conteo, y §1.3 la exhibe. Lo que sí tiene región límite de medida
cero es la **subcapa degenerada**: su condición
$q_0^{\,3}r\ge L$ fuerza $u_r\le 3\ln q_0/\ln z\to0$ cuando $q_0$
es acotado (doc 117b §3.2, proposición). La frase de 117b §3.2
debe leerse restringida a esa subcapa (objeto de 118b1).

## 2. La parte no degenerada de (BAND-LOW)

### 2.1 Por qué la ruta (2a) no existe

La vía (2a) — cota puntual
$(\widehat S(r)-\Pi(r))^+\le C_1\Pi(r)$ para
$\widehat V^{(r)}\ne\{1\}$ en la capa, vía la FL clásica del
modelo, más (REG-Π-LOW) — falla por sus dos patas:

- **La suma diverge aunque la cota puntual valiera.** Sumar la
  cota contra $1/\varphi$ es exactamente (REG-Π-LOW), falsa [T]
  (§1.3). El mayorante $\Pi/\varphi$ es demasiado grande en la
  capa baja, igual que en la alta (doc 118a §4).
- **La cota puntual no está probada en la capa.** La FL citada es
  vacua para $u_r<u_0^*$ por construcción de la capa (doc 113
  §7.2: su $O(u^{-u/2})$ es asintótico en $u$); la cola de
  $\widehat S$ tras excluir $\{1\}$ no hereda ninguna cota de la
  FL ahí. En ventana el cociente medido es $\le0{,}692$ (doc 117b
  §3.3 [M]), siempre en instancias con $u_r\ge1$ (véase §3).

### 2.2 Residuo mínimo exacto (rama 2b)

**(BAND-LOW-NONDEG) [P]:**
$$\sum_{\substack{r\in\mathrm{banda},\ u_r<u_0^*\\
\widehat V^{(r)}\ne\{1\}}}
\frac{(\widehat S(r)-\Pi(r))^+}{\varphi(r)}
\;\le\;C\,Q_{\mathrm{completo}},
\qquad\text{para el corte }u_0^*\text{ de 117b.}$$

Es un enunciado directo sobre $\widehat S$ en soportes no
degenerados de nivel bajo, sin mayorante: el paralelo exacto del
(BAND-MOD-HI) del doc 118a §4 para la capa alta. Dos precisiones:

- Las familias que refutan (REG-Π-LOW) son **eventualmente no
  degeneradas**: sus $u_r$ son constantes $>0$ fijas, mientras la
  degeneración exige $u_r\le3\ln q_0/\ln z\to0$ ($q_0$ el menor
  primo de $\mathcal P$ fuera de $r$, acotado en esas familias).
  Luego (BAND-LOW-NONDEG) hereda la pregunta abierta de fondo: si
  el cociente $(\widehat S-\Pi)^+/\Pi$ sobre esas familias se
  mantiene $\asymp\rho>0$ — y entonces (BAND-LOW-NONDEG)
  heredaría la divergencia $\log\log z$ — o decae con $z$. Es el
  mismo contenido residual que en (BAND-MOD-HI), no uno menor.
- El enunciado es uniforme en el corte: como (REG-Π-LOW) cae
  para todo $u_0^*$, ninguna elección de $u_0^*$ dentro de la
  arquitectura de 117b convierte la parte no degenerada en
  conteo. Lo que sí depende del corte es solo cuánto de la banda
  cae en este residuo y cuánto en (BAND-MOD-HI).

## 3. Medición [M]

Script `118b2_no_degenerada.py` (PASS, J0–J2b, 12,3 s; importa la
maquinaria del doc 110 y reproduce las funciones del doc 117b).

**J1 — marcos A/B** (25 niveles muestreados por marco,
$L=\lfloor D_A/s\rfloor$; máximos sobre niveles de las sumas en
unidades de $Q_c$; la caracterización degenerada
$\widehat V^{(r)}=\{1\}\iff q_0^{\,3}r\ge L$ se reverificó sin
excepción en todas las instancias):

| marco | $z$ | no-deg $u_r<1$: $(\widehat S-\Pi)^+$ | no-deg $u_r<1$: $\Pi$ | no-deg $1\le u_r<3/2$: $(\widehat S-\Pi)^+$ | no-deg $1\le u_r<3/2$: $\Pi$ | deg. | instancias (baja/alta/deg) |
|---|---|---|---|---|---|---|---|
| A | 17 | 0 | 0 | 0 | 0 | 0 | 0/0/0 |
| A | 29 | 0 | 0 | 0,0530 | 0,0965 | 0,0035 | 0/1/1 |
| A | 37 | 0 | 0 | 0,1570 | 0,2695 | 0,1822 | 0/13/28 |
| B | 29 | 0 | 0 | 0,0037 | 0,0072 | 0,0128 | 0/10/42 |

Lectura: **la celda no degenerada con $u_r<1$ está vacía en
ventana** en todos los marcos. Por tanto, para cortes $u_0^*\le1$
la suma de (BAND-LOW-NONDEG) es $0$ en ventana, y para
$u_0^*\in(1,\,3/2)$ coincide con la celda alta de la tabla
(máx. $0{,}157\,Q_c$ a $z=37$). El máximo cociente
$(\widehat S-\Pi)^+/\Pi=0{,}692$ del doc 117b se alcanza en
$u_r\ge1$, coherente con §2.1: en ventana, la capa baja no
degenerada aún no existe — es un fenómeno preasintótico más,
como el de 118a §3.

**J2/J2b — la familia $\mathcal F$ en ventana.** En su ventana
estricta $x_p\in(11/20,\,3/5)$ no hay **ninguna** instancia en
banda a $z\in\{29,37,53\}$, $u_L\in\{5/3,\,2\}$: a $z=37,53$ la
ventana $(z^{11/20},z^{3/5})$ no contiene primos (efecto de
rejilla), y a $z=29$ el único candidato ($r=35$) queda excluido
por $\Pi(r)\ge1/2$. En la ventana ampliada del J4 de 118a,
$x_p\in(1/2,\,2/3)$, a $z=53$, $u_L=2$: 5 elementos, **exactamente
3 en banda** — $r=77$ ($u_r=0{,}906$), $r=91$ ($0{,}864$),
$r=143$ ($0{,}750$) — en acuerdo exacto con el J4 de 118a, y los
tres **degenerados** ($\widehat V^{(r)}=\{1\}$). En ventana la
familia vive pues en $u_r<1$ y dentro de la subcapa de 118b1; su
hogar $(7/5,\,29/20)$ y su no-degeneración llegan solo con
$z\to\infty$ [J0(b)].

## 4. Estado final

- **(REG-Π-LOW) es FALSA para todo $u_0^*\in(0,\,3/2)$ [T]**
  (§1.3): la familia concreta de 118a es de capa alta, pero su
  mecanismo — primo pequeño libre × masa de Mertens positiva —
  diverge a toda profundidad de la capa baja (hermana $\omega=3$;
  escalera de $\omega$ creciente). La corrección a 117b §3.2
  queda registrada: la loseta $\{u_r<u_0^*\}$ tiene medida de
  Mertens positiva; la medida cero es solo de la subcapa
  degenerada.
- **La parte no degenerada de (BAND-LOW) NO se cierra módulo
  (REG-Π-LOW).** Residuo mínimo exacto: **(BAND-LOW-NONDEG) [P]**
  (§2.2), enunciado directo sobre $\widehat S$, en paralelo a
  (BAND-MOD-HI) [P] del doc 118a.
- Inventario actualizado de (BAND-MOD): subcapa degenerada
  (118b1) + (BAND-LOW-NONDEG) [P] + (BAND-MOD-HI) [P]. Ambas rutas
  de mayorante por $\Pi/\varphi$ — (REG-Π) y (REG-Π-LOW) — están
  refutadas [T]; lo que queda en ambas capas es el exceso real
  del modelo, no conteo.
- Cadena sin cambios: (DESP-SIGN) + (BAND-MOD) ⟹ (BAND-δ) ⟹
  (BAND-S); (H2) sigue abierta y Goldbach no se declara
  demostrado.

## 5. Reproducibilidad

`python3 118b2_no_degenerada.py | tee salida_doc118b2_local.txt`
— solo stdlib. J0: aritmética de Mertens de la reconciliación
(ventana de $\mathcal F$, fracción de masa con $u_r\ge1$, caja
$\omega=3$, escalera $\omega=5$). J1: medición en marcos A/B con
la partición degenerado/no-degenerado y el corte de referencia
$u_r=1$. J2/J2b: la familia de 118a en ventana, estricta y
ampliada. **PASS global.**
