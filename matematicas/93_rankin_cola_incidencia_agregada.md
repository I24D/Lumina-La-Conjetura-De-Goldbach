# 93 - Recorte de Rankin, incidencia agregada y la errata de (92.1)
## Procedencia y estado
- **Agente:** LUMINA/Muse (linea dedicada Goldbach). **Fecha:** 2026-11-27.
- **Metodo:** veredicto de Lumina Chat sobre los docs 85-92 (chat "Revisa UNIF U 2M X",
  2026-11-27) + calculo directo (integracion numerica propia de F/f por el sistema de
  retardo, como doc 89; fuerza bruta finita por criba). Solo stdlib. Fuentes primarias:
  Nathanson caps. 9-10; Chen 1973; Bordignon-Johnston-Starichkova arXiv:2207.09452;
  Halberstam-Richert (1974); Montgomery (1973); Gallagher; docs 85-92.
- **Nota de procedencia.** Los docs 81 y 83 NO estan en el workspace local (viven en
  OneDrive); sus enunciados se toman de las citas en los docs 85-92 y del veredicto de
  Lumina Chat del 2026-11-27. En particular: la identidad (81)
  `I_k = k*C_k + binom(k+1,2)*C_{k+1}` y el resultado del doc 83 (selector polilogaritmico
  + sub-incidencia u=1 bastan para `C_k + C_{k+1} > 0`). La identificacion exacta de
  las capas `C_k`/`C_{k+1}` con las capas P1/P2 de `J_k` (doc 88) queda **[P-rev]**:
  este doc las trata como las dos capas no negativas de (81) y mide el analogo P1/P2
  certificado por el doc 88, sin afirmar que sean el mismo objeto.
- **Estado:** **[T]** recorte de Rankin (sec. 1); **[T]** incidencia agregada controlada
  (sec. 2); **[T]** errata de (92.1) + necesidad de (q,N)=1 (sec. 3); **[T]** umbral de
  Cauchy corregido c<1/3 (sec. 4); **[T]** implicacion GRH de Lumina, verificada con
  una correccion de exponente (sec. 5); **[N]** separacion de las capas (sec. 6);
  **[N]** UNIF-U-2M-R fuera del track de lemas. No se demuestra Goldbach; no hay
  victoria parcial.

## 1. El recorte de Rankin: la cola u>y es exponencialmente pequena [T]
Marco (doc 89 sec. 1): `N` par grande, `Q=Q_N` selector de (83.1), `u` la parte
`Q`-radical de `N-p` (`rad(u)|Q`), `1 <= u <= N/Q`. La masa del termino principal
sumada sobre `u` es proporcional a `sum_u 1/u`, que converge a
`C_k = prod_{l|Q} l/(l-1)` (doc 89 (89.3)). Sea
```text
y = exp((log N)^{3/4}),   delta = (log N)^{-1/2}.                (93.1)
```
**Prop. 93.1 [T] (cola de Rankin).**
```text
sum_{u>y, rad(u)|Q} 1/u  <=  y^{-delta} * prod_{l|Q} (1 - l^{-(1-delta)})^{-1}
   <=  exp(-(log N)^{1/4}) * C_k * (1 + o(1)).                  (93.2)
```
*Prueba.* Para `u>y`, `1_{u>y} <= (u/y)^delta` (truco de Rankin). Luego
`sum_{u>y} 1/u <= y^{-delta} sum_u u^{delta-1}`, y la suma sobre todos los `u`
con `rad(u)|Q` es el producto de Euler `prod_{l|Q} (1-l^{delta-1})^{-1}`
(`delta<1`). Cada factor: `l^{-(1-delta)} = l^{-1} e^{delta ln l} <=
l^{-1}(1+2 delta ln l)` para `delta ln l` pequeno, luego
`(1-l^{-(1-delta)})^{-1} <= (1-l^{-1})^{-1} (1+4 delta ln l/(l-1))`, y el
producto sobre los `<=k-1` primos de `Q` es `<= C_k exp(4 delta sum ln l/(l-1))
= C_k (1+o(1))` (H1 verifica la desigualdad factor a factor). Y
`y^{-delta} = exp(-delta log y) = exp(-(log N)^{1/4})`. ∎
**Corolario [T].** La cola es menor que cualquier `(log N)^{-A}` fijo:
`exp(-(log N)^{1/4}) < (log N)^{-A}` cuando `(log N)^{1/4} > A ln log N`.
Medido (H1): el cruce esta en `log10 N* = 201332 (A=2)`, `28487863 (A=5)`,
`924283253 (A=10)` — asintotico, no computable, pero real; y el factor de
Rankin ya vale `0.145 (N=10^6)`, `0.101 (10^12)`, `0.00973 (10^200)`, con
`prod/C_k = 1.2077 -> 1.0270` en ese rango.
**Prop. 93.2 [T] (Rankin es suficiente pero no afilado).** La cola exacta es
mucho menor que la cota de Rankin: por conteo directo (89.2) +
sumacion parcial, la cola es `~ (log y)^{k-1}/(y (k-1)! prod log l)` (salvo
el termino de frontera, que se cancela a escala finita y no se usa como
estimacion). Medido (H2): `Q=3, N=10^12`: exacta `8.47e-06`, Rankin `0.171`
(`20243x`); `Q=3, N=10^24`: exacta `1.29e-09`, Rankin `0.107` (`8.3e7x`);
`Q=15, N=10^12`: exacta `4.61e-05`, Rankin `0.235` (`5110x`). El conteo
directo (termino integral, cota superior de la cola) la estima dentro de un
orden de magnitud (`directo/exacta = 6.1-23.6`). El valor de Rankin no es la
precision: es la elementaridad — solo usa el producto de Euler, sin el
conteo reticular (89.2). Para el veredicto de la sec. 6 basta cualquiera de
las dos: la cola es despreciable en masa de termino principal.

## 2. La incidencia agregada: los tres regimenes de u [T]
**Los tres regimenes.** (i) **Polilog** `u <= (log N)^t`: SW-Q se transfiere
limpio (doc 89 Prop. 89.1, `C0` constante) y la cota sumada es positiva y
cuantitativa (doc 89 Prop. 89.2). (ii) **Sub-polinomico** `u <= y`: la
transferencia necesita `C0` creciente, `C0 ~ (log N)^{3/4}/(2 ln ln N)`
(medido H4: `11.82 (10^12)`, `12.53 (10^24)`, `15.43 (10^100)`, `641.43
(log10 N = 180000)`); la geometria (`D_A >= y_A z_A`) enciende en
`log10 N0 = 4205` (H4, biseccion de `(1/15)L - (17/30)log A - B ln L = 0`);
el balance `Bal(c_eff)`, `c_eff = (log N)^{-1/4}`, se hace positivo en
`log10 N0 = 171634` (H4, `L* = 395202`, o sea `(1/c*)^4`, doc 89 H2).
A `log10 N = 5000` el balance aun es `-0.1301` (peor caso `Q=1`) /
`-0.0656 (Q=3)`; a `log10 N = 180000` es `+0.0008 / +0.0615`. (iii)
**Potencia** `u ~ N^c`, `c>0` fijo: el balance muere en `c* = 0.039884`
(doc 89) y el resto BV domina al termino principal por `u` (docs 90-91).
**Prop. 93.3 [T] (el agregado sobre u<=y esta controlado sin control por-u).**
En unidades `R` (conteo sin peso, doc 88 (88.4); `J_k` sigue por el peso
minimo `log A + log z_A ~ log N`, luego es de orden `N`):
```text
sum_{u<=y} R_{Qu}^{P<=2} >= sum_{u<=(ln N)^t} R_{Qu}^{P<=2}
   >= Bal(0) * (N*U_N*prod_{l|Q}(l-1)/(l-2)*10/(phi(Q)*log N))
      * sum_{u<=(ln N)^t} 1/u * (1+o(1)) - O(#u * N/log^{A'} N),   (doc 89 (89.4))
sum_{u<=y} R_{Qu}^{P<=2} <= (N/log N) * const * C_k * (1+o(1))
      + O(#(u<=y) * N/log^{A'} N),
   #(u<=y) ~ (log N)^{3(k-1)/4}/((k-1)! prod log l).
```
El error total es `<< N/log N` para `A' > 1+3(k-1)/4` (medido H4:
`err_ratio = #(u<=y)/(ln N)^9 = 1.17e-12 (10^12)`, `4.12e-47
(log10 N = 180000)`). El punto clave es la no-negatividad: cada
`J_k(N;Qu) = sum_v a_{Qu}(v)(log Qu + log rad v)` con `a >= 0` es `>= 0`
termino a termino, luego la cabeza polilog da la cota inferior del total
y la cola (intermedia + Rankin) solo anade masa no negativa cuya
participacion en el termino principal es despreciable (sec. 1). El control
por-u en el rango intermedio NO existe (el error BV absoluto domina al
termino principal por-u, doc 90 (I2)) y NO hace falta para el agregado.
**Medicion (H3, fuerza bruta, Q=3).** `J_u` para `u=3^a`, `N=10^5, 10^6`:
```text
N=10^5: Jtot=441506.2, frac(u=1)=0.691, frac(u<=9)=0.972, frac(u>y)=0.0005,
  P1=0.438, P2=0.558, PP=0.004; J_u*u/J_1 (u<=27) = 1.00, 0.94, 0.83, 0.81.
N=10^6: Jtot=4920540.3, frac(u=1)=0.682, frac(u<=9)=0.968, frac(u>y)=0.0003,
  P1=0.394, P2=0.605, PP=0.001; J_u*u/J_1 = 1.00, 0.96, 0.90, 0.88.
```
Lectura: la forma es `1/u` (dentro del 20% hasta `u=27`); la masa se
concentra en `u` pequenos (`u=1` lleva el 68-69%, frente al 1/1.5 = 67%
del peso reciproco; `u<=9` el 97%, frente al 96% teorico); la cola medida
mas alla de `y` es `<= 0.05%` de la masa (a escala finita dominan los
`X_u` pequenos, no la cola asintotica: el peso de termino principal de la
cola es `~5e-6`). La capa P2 lleva el 56-61% de `J_k` y la P1 el 39-44%
(vease sec. 6).

## 3. Errata: (92.1) debe imponer (q,N)=1 [T]
**Erratum E93.1.** El enunciado (92.1) del doc 92 y el del doc 91 sec. 6
deben leer
```text
sum_{q<=Q, (q,N)=1} |E(N;q,N)|^2 << N (log N)^C,   Q = N^{1/2+c/2}.   (92.1')
```
**Prop. 93.4 [T] (la restriccion es necesaria).** Si `(q,N) = d > 1`, la
clase `a = N mod q` no es coprima con `q`: `psi(N;q,N)` cuenta solo
potencias de primos que dividen a `d`, luego `psi(N;q,N) = O((log N)^2)`,
mientras el termino que se resta es `psi(N)/phi(q) ~ N/phi(q)`. Asi
`|E(N;q,N)| ~ N/phi(q)` y un solo termino no coprimo ya es
`~ N^2/phi(q)^2 >> N (log N)^C`. En particular `q=2` (`N` par):
`|E(N;2,N)| = psi(N) - O(log N) ~ N`. Medido (H5): `N=10^6, q=2`:
`|E|^2 = 9.99e11` frente al objetivo `N(ln N)^5 = 5.03e11` (ratio `1.99`:
rompe la cota); `N=10^5`: ratio `0.495` (el cruce `N > (ln N)^5` esta en
`N ~ 3.3e5`). Para `q|N` en general el ratio es `~ N/(phi(q)^2 (ln N)^5)`.
Luego (92.1) como esta escrita en el doc 92 es FALSA.
**Sobre la forma q=N de Lumina (correccion a la forma relayada).** El
resumen que me llego decia: "`q=N`, `E=-Lambda(N)`, `|E|^2=log^2 N` si `N`
es primo". Calculo exacto (H5, `N0 = 999983` primo):
`E(N0;N0,N0) = Lambda(N0) - psi(N0)/phi(N0) = log N0 - psi(N0)/(N0-1)
= 12.8159`, o sea `Lambda(N0) - 1`, NO `-Lambda(N0)` (signo cambiado), y
`|E|^2 = 164.2` frente a `N0 (ln N0)^5 = 5.03e11` (ratio `3.26e-10`):
ese termino solo NO rompe la cota de la suma. La forma `q=N` relayada no
carga la prueba; la conclusion de Lumina (la restriccion es necesaria) es
correcta y la prueba es `q|N` con `q` pequeno (en particular `q=2`).
Queda registrado asi, sin atribuirle a Lumina la forma corregida.
**El codigo del doc 92 ya lo hacia bien.** Lectura del script: H1/H2
saltan `q` con `gcd(q,N) != 1`. Verificacion (H5): recomputando H1 a
`N=10^5` CON restriccion da `S2fix = 3.119e5` (doc 92: `3.119e5`, exacto)
y SIN restriccion `2.019e10` (ratio `6.47e4`, dominado por `q=2`). Las
mediciones del doc 92 quedan intactas; solo el enunciado se corrige.

## 4. El umbral de Cauchy corregido: c<1/3, no c<1/5 [T]
**Prop. 93.5 [T] (correccion del doc 92 Prop. 92.1).** En
`R_B(A) = sum_s sum_d lambda_d(s) E(N;Asd,N)`, los modulos distintos son
`m = Asd` con `sd <= D_A`: hay a lo sumo `D_A` de ellos, no `Q`. El doc 92
uso `|R| <= sqrt(Q * N (ln N)^C)` con `Q = N^{1/2+c/2}`; lo correcto es
```text
|R_B(A)| <= sqrt(D_A * sum_{m<=D_A} |E(N;Am,N)|^2) * O(1)
         <= N^{3/4-c/4} (log N)^{C/2} * polylog,
```
usando (92.1') para la suma sobre `m <= D_A` (los modulos `Am` con
`(Am,N)=1` estan dentro de `q <= A*D_A = Q`). El cociente contra el
termino principal `(N/A)/log^2 N` es `N^{3c/4-1/4} (log N)^{C/2+2}`,
que `-> 0` si y solo si **`c < 1/3`** (no `c < 1/5`: el doc 92 conto
`Q/D_A = N^c` modulos de mas; en `|R|` un factor `N^{-c/2}`). Medido
(H6, `C=5, c=0.05`): `log10 R_corr = 3.94/2.74/-10.62/-30.52` a
`N=10^12/10^24/10^100/10^200` (doc 92: `4.24/3.34/-8.12/-25.52`); cruces
`log10 N*` corregido/doc92: `42.1/49.3 (c=0.05)`, `53.8/81.9 (c=0.10)`,
`59.4/104.0 (c=2/17)`; `D_A/Q = N^{-c} = 0.251` a `N=10^12, c=0.05`.
Esto CONFIRMA con derivacion independiente la correccion de Lumina
(`c<1/3`); la formula relayada "`(2*1/2+c)^2<2`" no es la forma que se
usa aqui: la comparacion es de exponentes `3/4-c/4` contra `1-c`.
**Consecuencia.** UNIF-U-2M-R bastaria hasta `c<1/3` via Cauchy: cubre
`c* = 0.039884` y el umbral geometrico `2/17 ~= 0.1176` con holgura. El
problema es el segundo momento mismo, no el paso de Cauchy.

## 5. La implicacion GRH de Lumina, verificada (con una correccion) [T]
**Prop. 93.6 [T].** UNIF-U-2M corregido (92.1') implica GRH para caracteres
primitivos no principales de conductor impar. *Prueba.* Fijense `r` impar
y `a` reducida mod `r`, `chi != chi0` primitivo mod `r`. (1) Para cada `x`
grande, elijase `N` par con `|N-x| <= 2r` y `N = a (mod r)` (TChR, `(2,r)=1`).
(2) Dominacion: `|E(N;r,a)|^2 <= sum_{q<=Q,(q,N)=1} |E(N;q,N)|^2`
(el termino `q=r` esta en la suma: `r <= Q` para `N` grande y `(r,N)=1`
porque `N = a (mod r)`), luego `|E(N;r,a)| << N^{1/2} (log N)^{C/2}`.
(3) Desplazamiento: `|psi(x;r,a)-psi(N;r,a)| <= (nº de enteros en (x,N]
congruentes a `a`) * log N = O(log x)`, y el termino principal se desplaza
`O(1)`; asi `|E(x;r,a)| << x^{1/2} (log x)^{C/2}`. (4) Ortogonalidad:
`psi(x,chi) = sum_a chi(a) psi(x;r,a) = sum_a chi(a) E(x;r,a)` (pues
`sum_a chi(a) = 0` para `chi != chi0`), luego
`|psi(x,chi)| <= phi(r) x^{1/2} (log x)^{C/2} <<_r x^{1/2} (log x)^{C/2}`.
(5) Mellin: `-L'/L(s,chi) = s int_1^inf psi(x,chi) x^{-s-1} dx` converge
absolutamente para `Re(s) > 1/2` con esa cota, luego `-L'/L` es holomorfa
ahi y `L(s,chi)` no tiene ceros con `Re(s) > 1/2`. ∎
**Correccion del exponente.** La cota correcta es
`|psi(x,chi)| <<_r x^{1/2} (log x)^{C/2}`, NO `r*x (log x)^{-C/2}` como
llego en el resumen: con `x^1` no se excluiria ningun cero con
`beta en (1/2,1)`; con `x^{1/2}` la implicacion es valida (para todo
`sigma > 1/2`, `x^{1/2} (log x)^{C/2} = o(x^sigma)`). La conclusion de
Lumina se sostiene con el exponente correcto.
**Verificacion finita del mecanismo (H7).** Ortogonalidad exacta por criba:
`r=3`: `psi(N,chi)` directo vs via `E`, diferencia `2.1e-12 (10^5)`,
`4.4e-12 (10^6)`; `r=5` (chi primitivo, `chi(2)=i`): `7.3e-13`, `8.7e-12`.
Desplazamiento: `|dpsi| <= (#n) ln N'` en los tres casos probados
(`11.513 <= 23.026`, etc.). La dominacion (2) es definicional (un termino
`<=` suma de no negativos).
**Consecuencia.** UNIF-U-2M-R es estrictamente MAS FUERTE que GRH: GRH da
la diagonal (doc 92 Prop. 92.3) pero no el fuera-diagonal. Sale del track
de lemas y queda como problema de investigacion aparte. Su unico residuo
para el track: el fuera-diagonal restringido a `(q,N)=1` (UNIF-U-2M-X con
la restriccion) — sigue sin tenerse, y sigue siendo Siegel-sensible
(doc 92 Prop. 92.5: un cero de Siegel con `(q0,N)=1` seguiria rompiendolo).

## 6. Veredicto: basta la incidencia agregada para la identidad sumada; no separa las capas [N]
**(b1) SI basta para la identidad sumada (81).** La cota inferior de `I_k`
viene de la cabeza polilog (doc 89 Prop. 89.2) mas la no-negatividad de la
cola; la cola es despreciable en masa de termino principal (Rankin, sec. 1;
conteo directo, sec. 2); y el agregado esta acotado entre dos cantidades de
orden `N/log N` (en unidades `R`) con constantes explicitas (Prop. 93.3).
Ninguna identidad exhibida hasta ahora exige una cota inferior por-u para
`u` grandes, y UNIF-U-2M-X NO hace falta para la identidad sumada.
**(b2) NO basta para separar las capas.** La identidad (81) es UNA ecuacion
lineal en DOS incognitas no negativas: `(C_k,C_{k+1}) = (0, I_k/binom(k+1,2))`
siempre es consistente con `I_k > 0`. Los pesos favorecen a la segunda capa
(`k=2`: `3` contra `2`; en general `binom(k+1,2) = k(k+1)/2 > k`). Medido
(H3): la capa P2 lleva el 56-61% de `J_k` y la P1 solo el 39-44%: el
escenario "toda la masa en la segunda capa" es el tipico, no un rincon.
**El problema central sigue abierto:** el mecanismo que separa las capas.
Ni la incidencia agregada, ni Rankin, ni UNIF-U-2M-X aunque se tuviera lo
dan: aun resolviendo el segundo momento, habria que justificar la
separacion (punto de Lumina, que se suscribe). Lo que separaria las capas:
una segunda relacion lineal independiente, o una cota inferior directa
para la capa `C_k` (p. ej. un switching solo-P1), o un argumento de
signo/cancelacion. Nada de eso se intenta aqui.
**Que identidades SI exigirian control individual de u grandes.** (a) La
uniformidad en el SELECTOR `Q <= N^theta` (doc 89 sec. 4): sigue
necesitando UNIF-U / UNIF-U-2M-R — esa es la linea viva, y es la de Lumina
(UNIF-U-COND, agrupacion por conductor primitivo). (b) Afirmaciones sobre
un `u` grande especifico: ninguna exhibida en los docs 81-92.
**Coordinacion.** Mi linea (Rankin + incidencia agregada) queda cerrada en
lo que cubre (secs. 1-2). La linea de Lumina Chat (UNIF-U-COND) es la viva
para la uniformidad en selector. La separacion de capas es el problema
central abierto y no pertenece a ninguna de las dos lineas tal como estan.

## 7. Verificacion numerica (control finito, no evidencia asintotica)
Script `93_rankin_cola_incidencia_agregada.py` (solo stdlib; ~2.3 s local):
```text
PASS H0 (reproduce docs 87-89) | f5=0.9982417384 F4=1.0216415525 c110=0.4909952010 J0=1.5238762554 Bal0=0.061404 c*=0.039884
PASS H1 (Rankin: cola u>y) | N=1e+06 L^1/4=1.93 y^-d=0.145 prod=1.8115 prod/Ck=1.2077 cota=0.263 L^-2=0.00524 L^-5=1.99e-06; N=1e+12 L^1/4=2.29 y^-d=0.101 prod=1.6973 prod/Ck=1.1315 cota=0.171 L^-2=0.00131 L^-5=6.21e-08; N=1e+24 L^1/4=2.73 y^-d=0.0654 prod=1.6298 prod/Ck=1.0865 cota=0.107 L^-2=0.000327 L^-5=1.94e-09; N=1e+100 L^1/4=3.90 y^-d=0.0203 prod=1.5585 prod/Ck=1.0390 cota=0.0317 L^-2=1.89e-05 L^-5=1.54e-12; N=1e+200 L^1/4=4.63 y^-d=0.00973 prod=1.5405 prod/Ck=1.0270 cota=0.015 L^-2=4.72e-06 L^-5=4.83e-14 || cruces fac<L^-A: A=2: log10 N*=201332; A=5: log10 N*=28487863; A=10: log10 N*=924283253
PASS H2 (Rankin vs conteo directo) | Q=3 N=1e+12 y=1.71e+05 cola_exacta=8.47e-06 rankin=0.171 directo=6.93e-05 directo/exacta=8.2 rankin/exacta=20243 || Q=3 N=1e+24 y=6.35e+08 cola_exacta=1.29e-09 rankin=0.107 directo=3.05e-08 directo/exacta=23.6 rankin/exacta=82648511 || Q=15 N=1e+12 y=1.71e+05 cola_exacta=4.61e-05 rankin=0.235 directo=0.000283 directo/exacta=6.1 rankin/exacta=5110 || Q=15 N=1e+24 y=6.35e+08 cola_exacta=2.03e-08 rankin=0.142 directo=2.02e-07 directo/exacta=10.0 rankin/exacta=6997058
PASS H3 (incidencia agregada vs u) | N=100000 y=518 Jtot=441506.2 frac(u=1)=0.691 frac(u<=9)=0.972 frac(u>y)=0.0005 P1=0.438 P2=0.558 PP=0.004; J_u*u/J_1(u<=27)=['1.00', '0.94', '0.83', '0.81'] || N=1000000 y=1295 Jtot=4920540.3 frac(u=1)=0.682 frac(u<=9)=0.968 frac(u>y)=0.0003 P1=0.394 P2=0.605 PP=0.001; J_u*u/J_1(u<=27)=['1.00', '0.96', '0.90', '0.88']
PASS H4 (rango sub-polinomico u<=y) | cero balance c_eff: log10 N0=171634 (L*=395202); encendido geometrico: log10 N0=4205 (L=9683) || log10N=12 c_eff=0.4362 geom=False Bal=n/a(sin geometria) Bal(Q=3)=n/a C0=11.82 #(u<=y)~11.0 err_ratio=1.17e-12; log10N=24 c_eff=0.3668 geom=False Bal=n/a(sin geometria) Bal(Q=3)=n/a C0=12.53 #(u<=y)~18.4 err_ratio=3.84e-15; log10N=100 c_eff=0.2567 geom=False Bal=n/a(sin geometria) Bal(Q=3)=n/a C0=15.43 #(u<=y)~53.8 err_ratio=2.96e-20; log10N=5000 c_eff=0.0965 geom=True Bal=-0.1301 Bal(Q=3)=-0.0656 C0=69.43 #(u<=y)~1011.7 err_ratio=2.85e-34; log10N=180000 c_eff=0.0394 geom=True Bal=0.0008 Bal(Q=3)=0.0615 C0=641.43 #(u<=y)~14868.7 err_ratio=4.12e-47
PASS H5 (errata 92.1: (q,N)=1 necesaria) | N=100000 q=2 |E|^2=1e+10 objetivo=2.02e+10 ratio=0.495 || N=100000 q=4 |E|^2=2.5e+09 objetivo=2.02e+10 ratio=0.124 || N=100000 q=5 |E|^2=6.25e+08 objetivo=2.02e+10 ratio=0.0309 || N=100000 q=10 |E|^2=6.26e+08 objetivo=2.02e+10 ratio=0.0309 || N=100000 q=25 |E|^2=2.49e+07 objetivo=2.02e+10 ratio=0.00123 || N=1000000 q=2 |E|^2=9.99e+11 objetivo=5.03e+11 ratio=1.99 || N=1000000 q=4 |E|^2=2.5e+11 objetivo=5.03e+11 ratio=0.496 || N=1000000 q=5 |E|^2=6.24e+10 objetivo=5.03e+11 ratio=0.124 || N=1000000 q=10 |E|^2=6.24e+10 objetivo=5.03e+11 ratio=0.124 || N=1000000 q=25 |E|^2=2.5e+09 objetivo=5.03e+11 ratio=0.00496 || H1-doc92 N=1e5: S2fix CON restriccion=3.119e+05 (doc92: 3.119e5) SIN=2.019e+10 ratio=6.47e+04 || q=N0 primo N0=999983: E=12.8159 |E|^2=164.2 vs N0(lnN0)^5=5.03e+11 ratio=3.26e-10 (NO rompe la suma)
PASS H6 (Cauchy corregido: umbral c<1/3) | N=1e+12 log10 R_corr=3.94 (doc92: 4.24); N=1e+24 log10 R_corr=2.74 (doc92: 3.34); N=1e+100 log10 R_corr=-10.62 (doc92: -8.12); N=1e+200 log10 R_corr=-30.52 (doc92: -25.52) || cruces log10 N* corr/doc92: c=0.0500: 42.1 / 49.3; c=0.1000: 53.8 / 81.9; c=0.1176: 59.4 / 104.0 || D_A/Q=N^-c a N=1e12 c=0.05: 0.251; c*=0.039884<1/3; 2/17=0.1176<1/3
PASS H7 (mecanismo GRH: ortogonalidad+desplazamiento) | N=100000 r=3: directo=128.5754 viaE=128.5754 dif=2.1e-12; r=5: dif=7.3e-13 || N=1000000 r=3: directo=715.0038 viaE=715.0038 dif=4.4e-12; r=5: dif=8.7e-12 || despl r=3 a=1: |dpsi|=11.513 <= 2*lnN'=23.026 || despl r=3 a=2: |dpsi|=0.000 <= 2*lnN'=23.026 || despl r=5 a=3: |dpsi|=11.513 <= 2*lnN'=23.026
PASS
```
Lectura honesta: H0-H2, H4-H7 son formulas, bisecciones e integracion numerica
de `F/f` por dos metodos (malla de retardo + Simpson con formas cerradas;
acuerdo `~1e-8`); H2/H3/H5/H7 son fuerza bruta finita (`N<=10^6`, 252 ceros no
usados aqui). H3 mide la FORMA (fracciones de masa, caida `1/u`, reparto
P1/P2), no el tamano asintotico. Los numeros deben ser IDENTICOS local y en
el PC. Correcciones honestas durante el trabajo: la primera version del
script tenia la criba de `spf` iterando potencias en vez de multiplos (bug
de un bucle; `omega/rad` salian cero) y la estimacion directa de H2 sin el
termino de frontera se cancelaba y daba negativa para `Q=15` (la forma
correcta es la integral como cota superior; la cola exacta es la medida).
Ambas se corrigieron antes de medir.

## 8. Estado y siguiente paso
- **[T]** Recorte de Rankin (Prop. 93.1-93.2): la cola `u>y` es
  `<= exp(-(log N)^{1/4}) C_k (1+o(1))`, menor que cualquier polilog; el
  conteo directo la estima dentro de un orden de magnitud.
- **[T]** Incidencia agregada controlada (Prop. 93.3): la cabeza polilog da
  la cota inferior, la cola es no negativa y despreciable en masa, el
  agregado esta acotado entre dos cantidades de orden `N/log N` (R) / `N`
  (`J_k`). UNIF-U-2M-X no hace falta para la identidad sumada.
- **[T]** Errata E93.1: (92.1) debe imponer `(q,N)=1`; la necesidad se
  prueba con `q|N` (H5); el codigo del doc 92 ya restringia y sus
  mediciones quedan intactas.
- **[T]** Cauchy corregido (Prop. 93.5): umbral `c<1/3` (doc 92 decia
  `c<1/5`); confirma la correccion de Lumina.
- **[T]** Implicacion GRH verificada (Prop. 93.6) con el exponente correcto
  `x^{1/2}`: UNIF-U-2M-R es mas fuerte que GRH y sale del track de lemas.
- **[N]** Separacion de las capas (sec. 6): una ecuacion, dos incognitas;
  la P2 lleva la mayor parte de la masa. Problema central abierto.
- **[N]** UNIF-U-2M-R: problema de investigacion aparte; su residuo es el
  fuera-diagonal restringido a `(q,N)=1`, Siegel-sensible.
- **[P-rev]** Identificacion exacta `C_k`/`C_{k+1}` con las capas P1/P2
  (doc 81 no disponible en local); re-escritura completa del Lema 32 de
  Bordignon con `Q` (heredada del doc 88).
- **Siguiente paso:** UNIF-U-COND (linea de Lumina Chat, agrupacion por
  conductor primitivo) para la uniformidad en selector; el mecanismo de
  separacion de capas queda como problema central abierto, sin ruta
  asignada. Goldbach sigue abierta.

Fuentes: Nathanson caps. 9-10; Chen 1973; Bordignon-Johnston-Starichkova
arXiv:2207.09452; Halberstam-Richert (1974); Montgomery (1973); Gallagher;
Barban-Davenport-Halberstam; Ford 2002; Maynard III arXiv:2006.08250;
Guth-Maynard arXiv:2405.20552; docs 64, 65, 72, 75, 85-92; veredicto de
Lumina Chat, chat "Revisa UNIF U 2M X", 2026-11-27.
