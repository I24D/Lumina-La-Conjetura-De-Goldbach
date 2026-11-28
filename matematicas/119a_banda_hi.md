# Doc 119a — (BAND-MOD-HI): la capa alta de (BAND-MOD) por la ruta directa sobre Ŝ

**Fecha:** 2026-11-28. **Línea:** UNIF-U / (H2) → (BAND-δ) → (BAND-MOD), continuación de los docs 117b y 118.
**Script:** `119a_banda_hi.py` (stdlib; importa la maquinaria del doc 117b).
**Salida:** `salida_doc119a_local.txt` — K0/K1, **PASS** (0,7 s).
**Etiquetas:** [T] = teorema, [M] = medición en ventana, [P] = pendiente.

---

## 0. Veredicto en breve

**(BAND-MOD-HI) es FALSA asintóticamente [T, módulo el término principal FL del modelo; M en ventana].**
El peso real sobre la capa alta no es menor que el mayorante muerto: la medición da
$(\widehat S-\Pi)^+/\Pi \approx 0{,}5$–$1{,}3$ (media ponderada por $\Pi/\varphi$, máx 1,2947)
sobre $u_r\in[1,3/2)$ — el defecto del modelo a $u$ fijo es un **término principal**
$W(u)-1>0$, no un error. Por tanto la suma directa sobre la familia $\mathcal F$ del
doc 118 hereda la misma divergencia $\log\log z$ que mató a (REG-Π), y las rutas L²
mueren con ella (§3). El residuo mínimo que sobrevive ya no habla de $\widehat S$:
es **(BAND-δ-HI)** [P], la capa alta de (BAND-δ) directa sobre
$\delta(r)=(\widehat S-\Pi)+E(r)$, cuyo contenido es la cancelación con $E(r)\le0$.

## 1. El enunciado y dónde estaba apoyado

**(BAND-MOD)** pide $\sum_{r\in\mathrm{banda}}(\widehat S(r)-\Pi(r))^+/\varphi(r)\le C\,Q_c$,
con banda $=\{r$ visible $: u_r<u_0=3/2,\ \Pi(r)<1/2\}$, $u_r=\ln(L/r)/\ln z$.
El doc 117b §3.1 cerraba la capa alta $\{u_r\ge u_0^*\}$ así: la FL clásica del modelo da
$|\widehat S(r)-\Pi(r)|\le\rho_{\mathrm{FL}}(u_r)\,\Pi(r)$ para $u_r\ge u_0^*$,
con $u_0^*$ el umbral absoluto de la FL citada ($u_0^*<3/2$: la capa incluye
$u_r\approx1{,}45$, cf. las rebanadas $(7/5,\,29/20)$ del doc 118a), y bastaba
**(REG-Π)**: $\sum_{\mathrm{banda}}\Pi/\varphi\le C\,Q_c$. El factor
$\rho_{\mathrm{FL}}(u_0^*)$ era una constante fija cualquiera: **el 100 % del peso
del cierre lo llevaba (REG-Π)**. El doc 118 probó (REG-Π) falsa [T] — la familia
$\mathcal F=\{pq: x_p\in(11/20,3/5),\ 5\le q<p\}$ a niveles $u_L\to2$ tiene masa
$\sum_{\mathcal F}\Pi/\varphi \big/ Q_c \sim \frac{\log(12/11)}{M_\infty}\log\log z\to\infty$ —
y dejó la capa alta como **(BAND-MOD-HI)** [P]: atacar $\sum(\widehat S-\Pi)^+/\varphi$
directo sobre $\widehat S$, sin pasar por $\sum\Pi/\varphi$.

La pregunta decisiva es cuánto vale el peso real $(\widehat S-\Pi)^+$ sobre esa familia
comparado con $\Pi$: si fuera $o(\Pi)$ o $\ll\Pi$, la suma directa convergería aunque
el mayorante no. La respuesta (§2): **vale $\approx\Pi$**, con cociente acotado entre
$\approx0{,}5$ y $\approx1{,}3$ en toda la capa. La razón estructural (§3): a $u_r$ fijo
$<2$, la condición de Nathanson $m=1$ excluye de $\widehat S$ los primos
$q\in[z^{u_r/2},z)$, cuyos términos $\mu/\varphi=-1/(q-1)$ son negativos y $\Pi$ sí
conserva; el defecto $\widehat S-\Pi$ es el término principal $W(u_r)-1>0$ de la FL
del modelo, no su error. Por eso el doc 117b midió $\widehat S-\Pi$ unilateral en la
banda (246/246 positivo): no es casualidad de ventana, es el signo del término principal.

## 2. La medición decisiva [M] (K1)

Suma real $\sum_{\mathrm{banda}}(\widehat S-\Pi)^+/\varphi$ normalizada por $Q_c$,
máximos sobre niveles (marcos A $z=17,29,37$; B $z=29$; estrés A $z=53$ con niveles
$u_L\approx2$). Identidad de control K0: $\widehat S+\widehat T=\Pi$ a $2{,}2\cdot10^{-16}$.

| marco | banda total | capa $u_r\ge1{,}25$ | capa $u_r\ge1{,}40$ | $\omega=1$ | $\omega=2$ | mayorante $\sum\Pi/\varphi$ |
|---|---|---|---|---|---|---|
| A $z=29$ | 0,1028 | 0,1028 | 0,1028 | 0,1028 | 0,0128 | 0,0965 |
| A $z=37$ | 0,2103 | 0,2103 | 0,1182 | 0,2103 | 0,0572 | 0,2695 |
| B $z=29$ | 0,0128 | 0,0128 | 0,0073 | 0,0000 | 0,0128 | 0,0124 |
| A $z=53$ | 0,3172 | 0,2695 | 0,1209 | 0,2695 | 0,0477 | 0,3911 |

(A $z=17$: banda vacía, $\Pi_{\mathrm{full}}=0{,}5156>1/2$ — doc 116c. Las capas son
sumas coherentes por nivel, luego maximizadas; las rebanadas finas están en la salida.)

**Perfil del cociente** $(\widehat S-\Pi)/\Pi$ (media ponderada por $\Pi/\varphi$;
máx entre paréntesis), por rebanada de $u_r$, todos los marcos con banda no vacía:
$[1{,}00,1{,}25)$: $1{,}07$–$1{,}24$ (máx 1,2947); $[1{,}25,1{,}40)$: $0{,}87$–$1{,}22$;
$[1{,}40,1{,}45)$: $0{,}69$–$1{,}07$; $[1{,}45,1{,}50)$: $0{,}55$–$0{,}70$.
El cociente **no decae** hacia la parte alta ni con $z$: es $\Theta(1)$ con suelo
medido $\approx1/2$ en $[1{,}25,3/2)$.

**Quién domina en ventana:** $\omega=1$ (conductores primos), con la rebanada
$\omega=2$ en $0{,}013$–$0{,}057$ — el mismo orden que la masa $\Pi$ de $\omega=2$
medida por el doc 118 a $z=53$ ($0{,}0436$–$0{,}0505$), cociente $\approx1$ incluido.
**La familia $\mathcal F$ no aparece en ventana** (fila 0,0000 en la salida): a
$z\le53$, $\Pi(pq)=\Pi_{\mathrm{full}}\frac{(p-1)(q-1)}{(p-2)(q-2)}>1/2$ para sus
miembros y la condición $\Pi<1/2$ los excluye de la banda; $\mathcal F$ entra solo
asintóticamente ($\Pi_{\mathrm{full}}\to0$). Su peso real se lee por el perfil:
sus miembros tienen $u_r\in(4/5,\,29/20)$ y el cociente a ese $u_r$ es el medido arriba.

**Validaciones cruzadas:** (i) la parte $\omega=1$ real ($0{,}21$–$0{,}27$) respeta
la cota [T] del doc 118a para su mayorante, $(\log(10/3)/M_\infty)\,Q_c\approx1{,}036\,Q_c$;
(ii) el máx del cociente, 1,2947, es $<A_{\min}=1{,}4273$, la constante puntual del
doc 113 para $(\widehat S-\Pi)^+$ — dos mediciones independientes del mismo objeto;
(iii) los momentos L² en ventana: $\sum(\widehat S-\Pi)^2/\varphi = 0{,}05$–$0{,}11\,Q_c$
y $\sum(\widehat S-\Pi)^2/(\Pi\varphi) = 0{,}11$–$0{,}26\,Q_c$ (su lectura asintótica, §3).

## 3. Análisis: por qué muere también la ruta directa, y la L² con ella

**El defecto es término principal [T, vía FL del modelo].** Para $u_r\in[u_0^*,3/2)$ fijo,
la FL clásica aplicada a la suma del modelo (la misma citada en los docs 112–113) da
$\widehat S(r)=\Pi(r)\,\bigl(W(u_r)+o(1)\bigr)$ con $W(u)>1$ para $u<2$: la condición
de Nathanson en la posición impar $m=1$ exige $r\,q^2<L$ para el primo mayor $q$ de
$d'$, así que $\widehat S$ omite todos los $d'$ con primo mayor en $[z^{u_r/2},z)$,
empezando por los primos solos, de peso total de Mertens $\log(2/u_r)+O(1)>0$ y signo
negativo en $\Pi$; la iteración de Buchstab conserva el signo (defecto unilateral,
246/246 en el doc 117b). En particular $(\widehat S-\Pi)^+\ge c\,\Pi$ sobre cualquier
familia con $u_r\in[u_0^*,3/2)$, con $c=\min(W-1)>0$ — en ventana $c\approx1/2$ (§2).
Esto refuta la esperanza $(\widehat S-\Pi)^+\le\Pi\cdot\varepsilon(u_r)$ con
$\varepsilon\to0$: a $u$ fijo no hay $\varepsilon$ pequeño; $\rho_{\mathrm{FL}}(u_0^*)$
del doc 117b era exactamente este $W(u_0^*)-1=\Theta(1)$.

**La divergencia se hereda [T].** Sea $u_0^*<29/20$ (la capa en juego). Sobre
$\mathcal F$ con la restricción $u_r\ge u_0^*$, i.e. $x_q\le 2-u_0^*-x_p$: fijando
$a>0$ con $x_p\in(11/20,\min(3/5,\,2-u_0^*-a))$ (intervalo no vacío), la masa $\Pi$
restringida es $\ge \Pi_{\mathrm{full}}\cdot\log\!\bigl(\tfrac{\min(3/5,2-u_0^*-a)}{11/20}\bigr)
\cdot\bigl(\log\log z+O(1)\bigr)$ — la suma en $q$ hasta $z^a$ sigue dando
$\log\log z+\log a$ (Mertens) — luego por el párrafo anterior
$\sum_{\mathcal F\cap\{u_r\ge u_0^*\}}(\widehat S-\Pi)^+/\varphi \big/ Q_c
\ge c'\log\log z\to\infty$. **(BAND-MOD-HI) es falsa en la misma escala $\log\log z$
que (REG-Π)**, y con ella la capa alta de (BAND-MOD) en su forma literal. En ventana
esto es invisible ($\log\log 53\approx1{,}38$; sumas medidas $0{,}10$–$0{,}32\,Q_c$),
exactamente como el doc 118 advirtió del mayorante.

**Autopsia de la ruta L².** Las desigualdades $(x)^+\le x^2/\Pi$ (solo en
$\widehat S\ge2\Pi$) y $(x)^+\le x^2/(2\Pi)+\Pi/2$ devuelven siempre un término
$\sum\Pi/\varphi$ sobre la parte complementaria $\{\Pi\le\widehat S<2\Pi\}$ — que es
donde vive la masa (cociente medido $\approx1$, i.e. $\widehat S\approx2\Pi$, justo en
la frontera) — es decir, la ruta muerta. Y el segundo momento mismo no salva: sobre
$\mathcal F$, $(\widehat S-\Pi)^2/(\Pi\varphi)\approx(W-1)^2\,\Pi/\varphi$ diverge
igual; el valor en ventana ($0{,}11$–$0{,}26\,Q_c$) es la misma ilusión $\log\log$.
No hay desigualdad puntual ni L² sobre $\widehat S-\Pi$ que cierre la capa alta,
porque el objeto sumado es genuinamente $\asymp\Pi$ ahí.

## 4. Veredicto y residuo mínimo

- **(BAND-MOD-HI): FALSA asintóticamente** [T módulo el término principal FL del
  modelo; M en ventana con cociente $\Theta(1)$]. Muere con ella la forma literal de
  (BAND-MOD) en su capa alta y el apoyo que el doc 117b daba a (BAND-δ) por esa capa.
- **(BAND-δ) NO queda refutada**: trata de $\delta(r)=(\widehat S-\Pi)+E(r)$, y
  $E(r)\le0$ (DESP-SIGN) puede cancelar el término principal $W-1$ — de hecho debe
  hacerlo, porque la masa real de la banda medida con $\delta$ es $\le0{,}3806\,Q_c$
  (doc 114) mientras la de $\widehat S-\Pi$ ya es comparable al mayorante.
- **Residuo mínimo: (BAND-δ-HI) [P]** — $\sum_{r\in\mathrm{banda},\,u_r\ge u_0^*}
  \delta(r)^+/\varphi(r)\le C\,Q_c$, directo sobre $S$ (no sobre $\widehat S$), cuyo
  contenido exacto es la cancelación $(\widehat S-\Pi)^+ + E(r)$ en la capa alta.
  Es un enunciado sobre el par $(\widehat S,E)$, no sobre el modelo solo: ninguna
  propiedad de $\widehat S$ por separado puede darlo, porque su término principal
  es positivo y divergente en la escala de la banda.
- Inventario tras 119a: la capa alta de la banda ya no es un problema de mayorantes
  ni del modelo: es el mismo tipo de compensación que (FLOW-COMP) en el lado del
  flujo. (BAND-MOD) sobrevive, si acaso, solo en su capa baja ($u_r<u_0^*$), fuera
  del alcance de esta sección.

## 5. Reproducibilidad

`python3 119a_banda_hi.py | tee salida_doc119a_local.txt` — stdlib únicamente;
importa `117b_banda_delta.py` (que encadena el doc 110). K0: identidad
$\widehat S+\widehat T=\Pi$ (error $2{,}2\cdot10^{-16}$). K1: barrido de la banda
en marcos A ($z=17,29,37$), B ($z=29$) y estrés A ($z=53$, niveles $u_L\approx2$):
suma real por rebanadas de $u_r$ y por $\omega(r)$, capas $u_r\ge1{,}00/1{,}25/1{,}40$,
mayorante $\sum\Pi/\varphi$, momentos L² y perfil del cociente $(\widehat S-\Pi)/\Pi$.
**PASS global.**
