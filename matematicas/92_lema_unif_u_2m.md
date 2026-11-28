# 92 - Lema UNIF-U-2M etapa 1: el segundo momento a a=N fijo
## Procedencia y estado
- **Agente:** LUMINA/Muse. **Fecha:** 2026-11-26 22:05 ET (reclamo) / entrega.
- **Metodo:** enunciado exacto + medicion finita (fuerza bruta de psi(N;q,a) por
  criba; 252 ceros de zeta del doc 72; bisecciones en espacio log). Solo stdlib.
  Fuentes primarias: Nathanson caps. 9-10; Bordignon-Johnston-Starichkova
  arXiv:2207.09452; Maynard III arXiv:2006.08250 (Teors. 1.1-1.2, sec. 5, sec. 11);
  Bombieri-Friedlander-Iwaniec 1986; Barban-Davenport-Halberstam; Montgomery (1973);
  Gallagher; Guth-Maynard arXiv:2405.20552; Ford 2002; docs 64, 65, 72, 75, 89-91.
- **Estado:** **[T]** enunciado y derivacion de Cauchy (sec. 1); **[T]** gran criba
  (sec. 2); **[T]** formula explicita, diagonal y Siegel (sec. 3); **[T]** espectral /
  Maynard sec. 5 (sec. 4); **[T]** densidad de ceros (sec. 5); **[N]** LEMA UNIF-U-2M
  sigue sin tenerse; queda el sub-sub-lema exacto UNIF-U-2M-X (sec. 6).
  No se demuestra Goldbach; no hay victoria parcial.

## 1. Enunciado exacto y la derivacion de Cauchy corregida
**LEMA UNIF-U-2M (doc 91 sec. 6).** Con `E(N;q,a) = psi(N;q,a) - psi(N)/phi(q)`:
```text
sum_{q<=Q} |E(N;q,N)|^2 << N (log N)^C   para Q = N^{1/2+c/2}, a = N fijo.   (92.1)
```
**Prop. 92.1 [T] (Cauchy corregido, doc 91 H5b).** De (92.1) seguiria, por Cauchy
sobre los modulos distintos `m = Asd <= Q` (multiplicidad `<= 5`, doc 85):
`|R| <= sqrt(Q * N (ln N)^C)`, y el cociente contra el termino principal
`(N/A)/ln^2 N` con `A = N^c` es `R_aleat = N^{5c/4-1/4} (ln N)^{C/2+2}`,
que `-> 0` si y solo si `c < 1/5`. Medicion (H0, C=5, c=0.05): log10 R_aleat =
4.24/3.34/-8.12/-25.52 a N=10^12/10^24/10^100/10^200, por la formula y por la
derivacion via sqrt, identicos a doc 91 sec. 4. La primera version de H5b del
doc 91 usaba una formula de R_2M sin derivacion valida; aqui se usa solo la
corregida. El umbral c<1/5 es el precio de Cauchy: (92.1) ni siquiera llega al
umbral geometrico 2/17 con margen comodo, y no cubre el balance c*=0.039884
con holgura asintotica (el cruce a c=0.10 esta en N ~ 10^29, doc 90 H3).

## 2. Gran criba: el promedio sobre a no es el punto fijo
**Prop. 92.2 [T] (BDH promedia; a fija no hay teorema).** Barban-Davenport-
Halberstam: `sum_{q<=Q} sum*_a |E(x;q,a)|^2 << x Q (log x)^A`. El promedio
sobre `a` es esencial en la prueba (la criba grande suma sobre todas las
clases). Aislar una sola clase cuesta hasta un factor `Q^2` en el peor caso:
`|E(x;q,a)|^2 <= sum*_a |E|^2`, y la masa de BDH podria concentrarse en un
solo `q` y una sola `a`; la cota que la criba grande da para `a` fija es
`<< x Q^2 (log x)^A`, un factor `Q^2` peor que (92.1). No hay version a fija
del teorema: la version a fija ES el modelo aleatorio (doc 91 sec. 4).
Medicion (H1, fuerza bruta exacta de psi por criba):
```text
N=10^5, Q=421 (c=0.05): S2avg=3.135e5, S2fix=3.119e5, fijo/promedio=0.995;
  S2/objetivo N(ln N)^5 = 1.55e-5; BDH floja por 3.7e-8; max Bq/media=1.52.
N=10^6, Q=1148/1412/1995: S2fix/objetivo = 1.01e-5/1.09e-5/1.24e-5;
  S2fix/N: 3.12 (1e5) -> 5.50 (1e6), o sea S2fix ~ N (ln N)^3.1 a esta escala.
```
Lectura honesta: a escala finita la clase `a=N` se comporta EXACTAMENTE como el
promedio (0.995) y el segundo momento esta 5 ordenes por debajo del objetivo
con C=5. Eso es evidencia, no teorema: la masa esta repartida (max/media 1.52),
pero nada impide que a `N` grande una sola clase se desvie — y la sec. 3 muestra
que un cero de Siegel lo haria. La peor clase medida es `q=119`, `|E|=369.2`
(fluctuacion sqrt, tamano esperado ~360). **Veredicto parcial:** la gran criba
no da (92.1); solo da el promedio. El sub-lema es la version a fija.

## 3. Formula explicita: la diagonal, el fuera-diagonal y Siegel
**Prop. 92.3 [T] (la diagonal sola exige GRH).** Por la formula explicita,
`E(N;q,N) = (1/phi(q)) sum_chi chi'(N) (-sum_rho N^rho/rho)`. La diagonal
`sum_rho N^{2 beta}/|rho|^2` con todos los `beta <= 1-w` es `~ N^{2-2w} (log N)^2`
por `(q,chi)`, y sumada sobre `q` da `~ N^{2-2w} (log N)^3`. Para que quepa en
el objetivo `N (log N)^C` hace falta `2-2w <= 1`, i.e. **`w >= 1/2`** — GRH
exacta. Con la franja conocida `w(N) = 1/(57.54 (ln N)^{2/3} (ln ln N)^{1/3})`
la diagonal es `~ N^{1.998}`. Medicion (H4): `0.5/w = 392/663/1901` a
N=10^12/10^24/10^100 (doc 91: 394/667/1923, redondeo). Bajo GRH la diagonal
cabe (`N (log N)^3 <= objetivo` con C=5), pero eso no basta: falta el
fuera-diagonal.
**Prop. 92.4 [T] (el fuera-diagonal es correlacion de pares entre caracteres).**
Cauchy sobre los caracteres pierde `phi(q)`: `|E|^2 <= (1/phi(q)) sum_chi
|E_chi|^2` da `Q*N` (sec. 2). La cancelacion real es entre caracteres: las
fases `chi'(N)` se comportan como aleatorias, y probar eso es controlar el
fuera-diagonal `sum_{chi1 != chi2} chi1'(N) chi2(N) (sum_rho N^rho/rho)
(sum_rho N^rho'/rho')` — correlacion de pares de ceros ENTRE funciones L
distintas, uniforme en `chi mod q`, `q <= Q`. Es ESTRICTAMENTE mas fuerte que
la correlacion de pares de Montgomery para zeta sola, que ya es conjetura
abierta incluso bajo GRH (doc 72: el plateau `|alpha| >= 1`; doc 74: SL-70.1
ni siquiera implica el lema del momento mixto como estaba enunciado — la
misma salvedad vale aqui: la correlacion controla una varianza promediada,
no el segundo momento puntual a `a=N`).
Medicion (H3, 252 ceros de zeta del doc 72, `gmax=473.8`): `Sdiag = sum 1/g^2
= 0.0213` (cola integral 0.0018, total ~0.0231); `|Z|^2 = |sum e^{i g L}/g|^2`,
`Soff = |Z|^2 - Sdiag`: `|Soff|/Sdiag = 0.495 (x=10^6), 0.623 (x=10^12)`,
frente a un benchmark de fase aleatoria de `0.96 x Sdiag`. El fuera-diagonal
es del ORDEN de la diagonal a altura finita (la mitad del benchmark aleatorio):
no se anula solo, y su control es exactamente el input de correlacion de
pares. A altura finita la formula explicita no se cierra sola.
**Prop. 92.5 [T] (sensibilidad a Siegel).** Un cero de Siegel `beta0 =
1 - c_s/log Q` para un `chi mod q0` contribuye a S2 `~ N^{2 beta0}/beta0^2
* (sum_m 1/phi(m)^2)/phi(q0)^2` (suma sobre multiplos de `q0`; la suma vale
3.3904). Medicion (H4, objetivo `N (ln N)^5`): ratio `1.35e3 (q0=3) / 149
(q0=7)` a N=10^12 con `c_s=1`, y `3.65e4 / 4.05e3` con `c_s=0.1`; a N=10^24:
`3.91e13 / 4.34e12` y `1.13e15 / 1.26e14`. Un cero de Siegel reventaria el
objetivo por 3 a 15 ordenes de magnitud. (92.1) es Siegel-sensible: ni el
teorema de Page (a lo sumo un modulo excepcional) ni Siegel (inefectivo) lo
salvan, porque la contribucion `N^{2 beta0}` con `beta0` cerca de 1 sigue
siendo `~ N^2`. **Veredicto parcial:** la formula explicita reduce (92.1) a
correlacion de pares entre caracteres + ausencia de Siegel; ninguno se tiene.

## 4. Espectral y Maynard III sec. 5: la ventana no abre
**Prop. 92.6 [T] (Teor. 1.1 es L^1, no L^2).** Maynard III Teor. 1.1 da
`sum_{q1,q2} sup_a |pi(x;q1q2,a) - pi(x)/phi| << delta pi(x) + x (ln ln x)^2/
(ln x)^2`, con `Q1 <= x^{1/10-3d}/(ln x)^C`, `Q2 <= x^{4/10+4d} (ln x)^C`
(Memoirs AMS 306(1544), 2025; verificado en el texto). Es una estimacion
L^1, uniforme en `a`, pero solo para modulos con un factor en la ventana
`x^{1/10}`. Para cubrir `q ~ Q = x^{1/2+c/2}` hace falta `d >= c/2` (nivel)
y `0.1-3d > 0` (que `Q1 >= 1`): el nivel maximo de Teor. 1.1 es
`x^{1/2+1/30}`, i.e. **`c < 1/15`**. Con `d = c/2` (minimo que alcanza el
nivel), la ventana de `q1` tiene ancho fijo `2 ln 2` y abre cuando
`(0.1-3c/2) L >= C ln L`. Medicion (H5, biseccion en espacio log):
```text
log10 N0 (primer entero en la ventana): c=0.02: 188.4 (C=5), 61.4 (C=2);
  c=0.05: 632.7 (C=5), 215.7 (C=2); c=0.10 y c=2/17: NUNCA (c>=1/15).
```
A escala computable la ventana esta VACIA (el impuesto `(log x)^C` supera
la escala `x^{1/10}`); para `c >= 1/15` no abre jamas. La celda balanceada
`M ~ K ~ x^{1/2}` no tiene pieza en la ventana: las Props. 5.1-5.3 de sec. 5
(la estimacion Tipo II principal) exigen un factor fuera del rango
balanceado, y el umbral 1/10 es agudo (doc 75: viene de `Q^{15/2} R^{10}
<= x^4`, Lema 11.8). El metodo espectral (Kuznetsov) tiene la misma
obstruccion: la celda balanceada es donde la criba espectral no ahorra
(doc 64 sec. 5; doc 65: el espectral no admite `a=N`).
**Prop. 92.7 [T] (L^1 -> L^2 pierde x^{1/2}).** Aun dentro de la caja,
Cauchy desde L^1: `sum|E_pi|^2 <= max|E_pi| * sum|E_pi|`, con
`max|E_pi| <= 8 x^{1/2-delta}/ln x` y `sum|E_pi| << delta x/ln x +
x (ln ln x)^2/(ln x)^2`. El cociente al objetivo `x (ln x)^3` es
`8 x^{1/2-delta} (ln ln x)^2/(ln x)^6`: medido (H5, delta=0.025):
`0.105 (10^12), 1.23e3 (10^24), 6e35 (10^100)`; cruce en `log10 N* = 15.2`
y diverge despues. Ninguna estimacion L^1 implica (92.1) por Cauchy: la
perdida es una potencia `x^{1/2}`. Y la caja cubre solo los modulos
factorizados, no todos los `q <= Q` (doc 91: masa fuera 3.0-17.8%).

## 5. Densidad de ceros
**Prop. 92.8 [T].** La densidad no toca la diagonal de sec. 3: con
`N(sigma,T) << T^{30(1-sigma)/13}` (Guth-Maynard 2024) la diagonal sigue
dominada por los ceros de altura pequena y `beta` cerca de `1-w`, donde la
densidad no ahorra nada (el conteo cerca del eje real es `O(1)` y la
region libre de ceros es la que manda: sec. 3). GM2024 alimenta cotas
ABSOLUTAS (Gallagher), no relativas (doc 91 sec. 5): no toca el filtro de
la sec. 0 del doc 91. **Veredicto parcial:** densidad de ceros NO da (92.1).

## 6. Veredicto y sub-sub-lema que falta
Ninguna herramienta da (92.1): la gran criba solo promedia sobre `a` (sec. 2);
la formula explicita lo reduce a correlacion de pares entre caracteres +
ausencia de Siegel (sec. 3); el espectral / Maynard III no llega a la celda
balanceada ni al nivel `c >= 1/15`, y L^1 no implica L^2 (sec. 4); la
densidad de ceros alimenta cotas absolutas (sec. 5). El cuello es unico: el
segundo momento a `a` fija es un objeto RELATIVO (cancelacion entre
caracteres), y todas las herramientas disponibles dan cotas absolutas o
promediadas. A escala finita todo se comporta bien (H1: fijo/promedio =
0.995; H3: fuera-diagonal ~0.5x diagonal); la apertura es asintotica y de
peor caso (Siegel).
**LEMA UNIF-U-2M-X (el sub-sub-lema exacto que falta):** el fuera-diagonal de
la formula explicita para `sum_{q<=Q} |E(N;q,N)|^2` — la suma sobre pares
`(chi1,rho1) != (chi2,rho2)` de `chi1'(N) chi2(N) N^{rho1+rho2'}/(rho1 rho2')`
— es `<< N (log N)^C`, uniformemente para `q <= Q = N^{1/2+c/2}`. Es
correlacion de pares de ceros entre funciones L distintas, uniforme en
`chi mod q`; mas fuerte que Montgomery para zeta sola (doc 72). Goldbach
sigue abierta.

## 7. Verificacion numerica (control finito, no evidencia asintotica)
Script `numerico/92_lema_unif_u_2m.py` (solo stdlib; ~0.4 s local):
```text
PASS H0 (Cauchy corregido doc91) | log10 R_aleat c=0.05 N=1e12/1e24/1e100/1e200:
  [4.24, 3.34, -8.12, -25.52] (via sqrt: identico; doc91: identico); umbral c<1/5
PASS H1 (fijo vs promedio) | N0=1e5 Q0=421 S2avg=3.135e5 S2fix=3.119e5 fijo/prom=0.995;
  S2/objetivo=1.55e-5; BDH floja 3.7e-8; max Bq/media=1.52 (q=389); S2fix/N 3.12->5.50, alpha=3.12;
  N1=1e6 fijo/objetivo 1.01e-5/1.09e-5/1.24e-5 (c=0.02/0.05/0.10)
PASS H2 (peor clase) | max|E(N=1e6;q,N)|=369.2 en q=119 (fluctuacion sqrt, esperado ~360)
PASS H3 (fuera-diagonal zeta) | n=252 gmax=473.8 Sdiag=0.0213 (cola 0.0018, total ~0.0231);
  benchmark aleatorio 0.0204 (0.96x Sdiag); |Soff|/Sdiag=0.495 (1e6), 0.623 (1e12)
PASS H4 (Siegel + franja) | 0.5/w=392/663/1901 (1e12/1e24/1e100); Siegel ratio 1.35e3/149
  (1e12, cs=1, q0=3/7), 3.65e4/4.05e3 (cs=0.1); 3.91e13/4.34e12, 1.13e15/1.26e14 (1e24)
PASS H5 (ventana Maynard + L1->L2) | log10 N0: c=0.02: 188.4 (C=5), 61.4 (C=2);
  c=0.05: 632.7 (C=5), 215.7 (C=2); c=0.10, 2/17: NUNCA; L1->L2 ratio 0.105/1.23e3/6e35,
  cruce log10 N*=15.2
PASS
```
Lectura honesta: H1-H3 son fuerza bruta finita (N<=10^6, 252 ceros); H0, H4, H5
son formulas y bisecciones. Los numeros deben ser IDENTICOS local y en el PC.
Correccion honesta durante el trabajo: la primera version del script tenia la
criba de Lambda marcando compuestos como primos (bug de un byte) y la
contribucion de Siegel con un factor `(ln Q)^2` que no pertenece a la formula
explicita; ambas se corrigieron antes de medir (los numeros de H1/H2/H4 de la
primera corrida eran artefactos y se descartaron).

Fuentes: Nathanson caps. 9-10; Chen 1973; Bordignon-Johnston-Starichkova
arXiv:2207.09452; Maynard III arXiv:2006.08250 (Memoirs AMS 306(1544), 2025);
BFI 1986; Barban-Davenport-Halberstam; Montgomery 1973; Gallagher; Guth-Maynard
arXiv:2405.20552; Ford 2002; Page; docs 64, 65, 72, 74, 75, 89-91.
