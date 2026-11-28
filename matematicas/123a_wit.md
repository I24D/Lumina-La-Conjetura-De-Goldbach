# Doc 123a — (WIT) REFUTADO: las desigualdades de posición son exactas y el contraejemplo vive en su ventana

**Fecha:** 2026-11-28. **Línea:** UNIF-U / (H2), sección del doc 123
sobre el sub-lema (WIT) registrado en el doc 122b §4.
**Script:** `123a_wit.py` (stdlib; importa la maquinaria de los
docs 109/110/117b/118b1/122b sin tocarla). **Salida:**
`salida_doc123a_local.txt` — chequeos K0–K3, PASS (10,4 s).
**Etiquetas:** [T] = teorema, [M] = medición/verificación, [P] = pendiente.

---

## 0. Veredicto en breve

**(WIT) está REFUTADO como lema uniforme.** El intento de
demostración (§2) produce la familia de desigualdades de
posición (POS) [T], cuya forma más fuerte sobre el cofactor es
$c\,p_k^2 < q_0^3$ en la última posición impar — una cota
**exacta**, no una pérdida de técnica: deja pasar
$c = q_0^{++}$ en el estrato 7 ($13 < 343/25 = 13{,}72$) y
$c \in \{17, 19, 23\}$ en el estrato 11 ($c < 1331/25 =
53{,}24$). La búsqueda adversarial (§3) realiza ambas
ventanas con instancias verificadas predicado a predicado con
las funciones de la serie:

- **Estrato 7, $z = 61$** (marco A, $L = 39179$): $r = 115 =
  5 \cdot 23$ admite testigos $c \in \{1, 7, 11, 13\}$;
  $c = 13 = q_0^{++}$ y $N(r) = 4 > 3$.
- **Estrato 11, $z = 127$** (marco A, $s = 127$, $L = 127^4$):
  $r = 5 \cdot 7 \cdot 101 \cdot 103$ admite testigos
  $c \in \{1, 11, 13, 17, 19, 23\}$; $N(r) = 6$.

(WIT) sigue siendo cierto en la ventana medida original
$z \le 53$ [M, doc 122b; K0 aquí] — donde el estrato 11 ni se
enciende y el 7 aún no alcanza su ventana — y para el estrato
5 sigue sin aparecer contraejemplo (§3), con el caso residual
exacto nombrado (§2). (SLAB-V′) no se mueve como enunciado
[P]; lo que muere es el atajo de las tres sumas (§4).

## 1. Definiciones exactas [transcripción de la maquinaria]

**Marco.** $z$ fijo; $\mathcal P$ = `allow` $= \{p \text{ primo}:
5 \le p < z\}$ (el 2 está excluido siempre y el 3 divide a
$A = 3$); niveles $\mathcal L_A = \{\lfloor D/s \rfloor : s
\text{ primo}, z \le s < y, (s, 3) = 1\}$ con marco A
($D = z^5$, $y = z^{10/3}$) o B ($D = z^6$, $y = z^3$).

**El soporte.** $V^+(L)$ = libres de cuadrados $d$ con primos
en $\mathcal P$, $d < L$, que cumplen las **condiciones de
Nathanson** (Teo. 9.3, $\beta = 2$; `in_Vplus` del doc 109):
con los primos de $d$ en orden decreciente
$p_1 > p_2 > \dots > p_k$, para todo $m$ impar,
$1 \le m \le k$:
$$ p_m^2 \cdot (p_1 p_2 \cdots p_m) \;<\; L. $$
(La condición $m = 1$ es $p_1^3 < L$.) La **lista fusionada**
de un par $(c, r)$ es la lista decreciente de primos de
$d = cr$: las condiciones se evalúan sobre ella, no sobre $c$
y $r$ por separado — de ahí la no-multiplicatividad del doc
122b §3.

**Los objetos.** $r$ es *visible* si $r \mid d$ para algún
$d \in V^+(L)$; está *en banda* si $u_r := \log(L/r)/\log z
< 3/2$ y $\Pi(r) < 1/2$, con $\Pi(r) = \prod_{p \in \mathcal
P \setminus r} (1 - \frac{1}{p-1})$ (`Pi_of` del doc 109).
$q_0(r) := \min(\mathcal P \setminus \mathrm{primos}(r))$ es
el *estrato* de $r$: todo $r$ de estrato $q_0$ contiene todo
primo de $\mathcal P$ menor que $q_0$ y $q_0 \nmid r$. $r$ es
*degenerado* si $\widehat V^{(r)} = \{1\}$, equivalentemente
[T, doc 117b §3.2] $q_0(r)^3 \cdot r \ge L$ — donde
$\widehat V^{(r)}$ es el soporte desplazado del *modelo* (un
primo $q \notin r$ entra en él sii $q^3 < L/r$), no el soporte
real de testigos: por eso un $r$ degenerado puede tener
testigos $cr \in V^+(L)$ con $c > 1$ aunque $q_0 \notin
\widehat V^{(r)}$.

**El par y la multiplicidad** [T, doc 119b §1]. Para $r$
degenerado visible en banda y un testigo $d \in V^+(L)$,
$c := d/r < q_0^3$, $(c, r) = 1$, $(c, P(q_0)) = 1$ con
$P(q_0) = \prod_{p \in \mathcal P, p < q_0} p$,
$\varphi(cr) = \varphi(c)\varphi(r)$, y
$N(r) := \#\{c < q_0^3 : (c, r) = 1,\ cr \in V^+(L)\} \ge 1$.
Los primos de $c$ están en $\mathcal P \setminus r$, luego
son $\ge q_0$.

**(WIT)** [P, doc 122b §4]: para todo estrato
$q_0 \in \{5, 7, 11\}$, todo $L \in \mathcal L_A$ y todo $r$
degenerado visible en banda de estrato $q_0$: todo testigo
cumple $c \in \{1, q_0, q_0^+\}$ ($q_0^+$ = siguiente primo
de $\mathcal P$ tras $q_0$) y $N(r) \le 3$.

## 2. Intento de demostración: las desigualdades (POS) y su exactitud

Sea $r$ como en (WIT), $c$ un testigo, $d = cr$, con lista
fusionada $p_1 > \dots > p_k$. Para cada posición impar $m$,
la condición de Nathanson y la degeneración $L \le q_0^3 r$
dan, con $r = d/c$:
$$ p_m^2 \prod_{j \le m} p_j \;<\; L \;\le\; q_0^3\, r
\;=\; q_0^3\, d / c,
\qquad\text{es decir}\qquad
\textbf{(POS}_m\textbf{)}\quad c \cdot p_m^2 \;<\; q_0^3
\prod_{j > m} p_j . $$
Las dos posiciones útiles son las de abajo de la lista
(arriba, $\prod_{j>m}$ es grande y la cota es vacía):

- **Última posición** ($m = k$, si $k$ es impar):
  $c \cdot p_k^2 < q_0^3$, con $p_k$ el menor primo de $d$.
- **Penúltima** ($m = k-1$, si $k$ es par):
  $c \cdot p_{k-1}^2 < q_0^3 \cdot p_k$.

**Lo que (POS) sí demuestra** [T]. En los estratos 7 y 11 el
menor primo de $d$ es $p_k = 5$ (forzado: $5 \in r$). En el
estrato 5, $5 \notin r$:

- *Estrato 5, $k$ impar:* si $5 \mid c$, $p_k = 5$ y
  $25c < 125$ da $c < 5$, imposible; si $5 \nmid c$,
  $p_k \ge 7$ (menor primo de $r$) y $49c < 125$ da $c = 1$.
  Luego $c = 1$: en esta paridad (WIT) vale con holgura.
- *Estrato 5, $k$ par:* (POS$_{k-1}$) con $p_k = 5$ (si
  $5 \mid c$) da $c < 125 \cdot 5/49 < 12{,}8$, luego $c = 5$;
  con $5 \nmid c$ y menor primo de $r$ igual a $7$, da
  $c < 125 \cdot 7/121 < 7{,}3$, luego $c \in \{1, 7\}$. El
  **único caso residual** del estrato 5: $k$ par, $5 \nmid c$,
  menor primo de $r$ $\ge 11$ (es decir $7 \in r$ o $7 \mid
  c$ con $c \ge 77$ ya excluido por tamaño... en concreto
  $p_k = 11$), donde (POS$_{k-1}$) solo da $c < 11{,}4$ y
  $c = 11 = q_0^{++}$ no queda excluido por esta vía.

**El paso que no cierra — y por qué no es técnica** [T la
obstrucción]. En el estrato 7, (POS$_k$) con $p_k = 5$ da
$c < 343/25 = 13{,}72$: el primo $c = 13 = q_0^{++}$ la
satisface con margen. En el estrato 11 da $c < 1331/25 =
53{,}24$: satisfecha por $13, 17, 19, 23, 29, 31, 37, 41,
43, 47$. Ninguna otra posición mejora esto: las posiciones
impares superiores tienen $\prod_{j>m} \ge p_k$ y dan cotas
más débiles, y la condición $m = 1$ controla $p_1$, no $c$.
El mecanismo que (WIT) imaginaba — «un primo $\ell \mid c$
con $\ell \ge q_0^{++}$ produce un divisor visible que
contradice $\widehat V^{(r)} = \{1\}$» — **no existe**: la
degeneración es una propiedad del modelo $\widehat V^{(r)}$
(condición $q^3 < L/r$ sobre $m$ solo), mientras que el
testigo real $cr$ se juzga por la lista fusionada, donde
$\ell$ cae en una posición baja de la lista, entre los
primos pequeños de $r$, y su condición de Nathanson es la
débil $p_k^2 d < L$. El caso $c = q_0^2$ ni siquiera es posible ($d$ es libre de
cuadrados), y un $c$ compuesto paga $c \ge q_0 q_0^+$
($\ge 35$; en el estrato 11, $\ge 11 \cdot 13 = 143$):
(POS$_k$) y (POS$_{k-1}$) lo excluyen en los tres estratos
y en ambas paridades — en el 11, la paridad par da
$c < 1331 \cdot 5/49 = 135{,}8 < 143$. Es decir: la parte
«jamás ocurre un $c$ compuesto» de la medición del doc
122b **sí** es consecuencia de (POS); (WIT) muere
únicamente por cofactores primos. **Conclusión del intento: la
cota $c < q_0^3/p_k^2$ es la verdad del asunto en la última
posición; (WIT) pide más de lo que la estructura da.**

## 3. Búsqueda adversarial [M]: la ventana de (POS) se realiza

**K0 — validación.** Con la maquinaria real (Level110 +
`band_degen` + el `level_data` del doc 122b), los conjuntos
de $c$ de los marcos A ($z = 17, 29, 37, 53$) y B ($z = 29$)
reproducen exactamente la tabla del doc 122b §2
($\{1,5,7\}$, $\{1,7,11\}$, $N \le 3$): (WIT) vale en
$z \le 53$.

**K1 — barrido extendido.** Un evaluador ligero (mismos
predicados; validado contra K0 a $z = 53$ con conjuntos
idénticos y error de $\Pi$ de $4{,}8 \cdot 10^{-16}$) barre
los niveles de la muestra de la serie en los marcos A
$z = 61, 67$ y B $z = 61$. Estrato 5: limpio ($c \in
\{1,5,7\}$, $N \le 3$). Estrato 11: sin degenerados en banda.
**Estrato 7: violaciones en marco A ya a $z = 61$** —
$L = 39179$, $r = 115$: $c = 13$ ocurre y $N(r) = 4$; y a
$z = 67$ — $L = 28942$, $r = 85$: idéntico patrón. (En marco
B a $z = 61$ el estrato 7 sigue limpio.)

**K2 — las instancias, verificadas predicado a predicado**
(con `allowed_primes`, `Pi_of`, `in_Vplus` de la serie, sin
enumerar el universo):

*Estrato 7* ($z = 61$, marco A, $L = 39179 \in \mathcal
L_A$): $r = 115 = 5 \cdot 23$. $q_0(r) = 7$; degenerado:
$7^3 \cdot 115 = 39445 \ge L$; banda: $\Pi(r) = 0{,}495081
< 1/2$, $u_r = 1{,}4184 < 3/2$; visible: $r \in V^+(L)$.
Testigos: $c \in \{1, 7, 11, 13\}$, todos con $cr \in
V^+(L)$ verificado. **$c = 13 = q_0^{++}$ viola (WIT), y
$N(r) = 4 > 3$.** La condición que lo deja pasar es la de la
última posición de $d = 1495$ (primos $23, 13, 5$):
$25 \cdot 1495 = 37375 < 39179$ — la ventana $13 < 13{,}72$
de (POS$_k$), realizada con margen del 4,6 %.

*Estrato 11* ($z = 127$, marco A, $s = 127$, $L = 127^4 =
260144641$): $r = 5 \cdot 7 \cdot 101 \cdot 103 = 364105$.
$q_0(r) = 11$; degenerado: $11^3 \cdot 364105 = 484623755
\ge L$; banda: $\Pi(r) = 0{,}495342 < 1/2$, $u_r = 1{,}3566
< 3/2$; visible: $r \in V^+(L)$ (márgenes de Nathanson
$1092727$ y $3568229$, ambos $< L$). Testigos: $c \in
\{1, 11, 13, 17, 19, 23\}$, $N(r) = 6$. **$c = 17, 19, 23$
violan (WIT)**; p. ej. para $c = 17$, $d = 6189785$ con
primos $103, 101, 17, 7, 5$ y márgenes $1092727$,
$51109939$, $154744625$, todos $< L$.

**K3 — robustez.** La instancia del estrato 11 no es un
filo: al barrer los niveles del marco A con $s$ sprimo en
$[127, 8000]$, **17 niveles** reproducen la violación para
el mismo $r$ (los cinco primeros con el conjunto completo
$\{1, 11, 13, 17, 19, 23\}$). La del estrato 7 aparece en
dos marcos y dos valores de $z$ (K1).

**Estrato 5.** Ni el barrido K1 ($z = 61, 67$, marcos A/B)
ni el dirigido K2 ($z \in \{127, \dots, 251\}$) encuentran
contraejemplo: (WIT) para $q_0 = 5$ queda [P], apoyado por
la medición y por (POS) salvo el caso residual de §2
($k$ par, $5 \nmid c$, menor primo de $r$ $\ge 11$,
$c = 11$ candidato único).

## 4. Consecuencia para (SLAB-V′)

(WIT) muere como lema uniforme; esto es lo que le pasa a
(SLAB-V′) del doc 119b §4 / 122b §1:

1. **El enunciado no cambia.** (SLAB-V′) sigue siendo [P]
   tal cual: $P(q_0) = \sum_{c < q_0^3} U_c^{\mathrm{deg}}
   (q_0) \le C(q_0)\,Q_c$ uniformemente en $L \in \mathcal
   L_A$, $q_0 \in \{5, 7, 11\}$. La refutación de (WIT) no
   toca ni una cota: ninguna medición sugiere que $P(q_0)$
   crezca; lo que crece con $z$ es el *soporte* de
   cofactores que la suma recorre.
2. **La reducción a tres sumas por estrato es falsa.** La
   forma reducida que (WIT) prometía,
   $U_1^{\mathrm{deg}} + U_{q_0}^{\mathrm{deg}} +
   U_{q_0^+}^{\mathrm{deg}}$, no es (SLAB-V′): ya a
   $z = 61$ falta $U_{13}^{\mathrm{deg}}(7)$ y a $z = 127$
   faltan $U_{17}^{\mathrm{deg}}(11)$,
   $U_{19}^{\mathrm{deg}}(11)$, $U_{23}^{\mathrm{deg}}(11)$.
   La familia que hay que acotar es la completa: por
   estrato, todos los $c < q_0^3$ libres de cuadrados con
   $(c, P(q_0)) = 1$ y $\omega(c) \le 3$ — finita en formas
   (la lista no depende de $z$; sus valores recorren
   $\mathcal P$), pero con soporte efectivo creciente en
   $z$, como ya pasaba con el mayorante $W$ del doc 119b §2.
3. **La cota $N(r) \le 3$ también muere** ($N = 4$ a
   $z = 61$; $N = 6$ a $z = 127$): la forma de pares
   $P(q_0) = \sum_r N(r)/\varphi(r)$ no admite el atajo
   $P \le 3 \times$ verdadera. Lo que (POS) sí da [T] es
   una cota del tipo $N(r) \le \#\{c : c < q_0^3/p_k^2\}$-es
   decir, $N(r)$ está controlado por la ventana de la
   última posición, no por 3.
4. **Lo que sobrevive de este doc para el ataque de
   Mertens**: (i) la familia (POS) [T], que acota cada
   cofactor por la geometría de la lista fusionada —
   en particular ningún $c$ compuesto aparece en los
   estratos 5 y 7, y en el 11 todo $c$ que ocurre cumple
   $c < 53{,}24$ en la paridad impar; (ii) el estrato 5
   reducido de facto a $\{1, 5, 7\}$ en todo lo medido,
   con un único caso residual nombrado; (iii) la lección
   estructural: la degeneración del cociente no confina
   los testigos a la vecindad de $q_0$, porque el testigo
   se juzga en la lista fusionada y no en el modelo
   $\widehat V^{(r)}$. El ataque a (SLAB-V′) debe tratar
   las $U_c^{\mathrm{deg}}$ como familia completa desde el
   principio.

**VEREDICTO: (WIT) REFUTADO** — contraejemplos verificados
en los estratos 7 ($z = 61$) y 11 ($z = 127$), con la
obstrucción identificada antes por la vía de la prueba:
la cota exacta de última posición $c < q_0^3/p_k^2$ es
estrictamente más débil que (WIT) en cuanto $q_0 \ge 7$.
Para el estrato 5, (WIT) queda [P] con un caso residual
único nombrado. (SLAB-V′) permanece [P] sin cambios de
enunciado; su reducción a tres sumas por estrato queda
descartada. Goldbach no se declara demostrado.

## 5. Reproducibilidad

`python3 123a_wit.py | tee salida_doc123a_local.txt` —
stdlib; importa `109_ataque_h2_gram.py`,
`118b1_degenerada.py`, `122b_slab_v_prima.py`. K0:
reproducción exacta de los conjuntos de $c$ del doc 122b
con la maquinaria real. K1: evaluador ligero validado
contra K0 ($z = 53$) y barrido A $z = 61, 67$ / B $z = 61$.
K2: barrido puntual dirigido y verificación predicado a
predicado de las instancias. K3: robustez en niveles de la
instancia del estrato 11. **PASS global** (10,4 s).
