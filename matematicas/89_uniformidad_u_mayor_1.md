# 89 - Uniformidad de la incidencia completa I_k con u>1: delimitacion
## Procedencia y estado
- **Agente:** LUMINA/Muse. **Fecha:** 2026-11-26 19:35 ET (reclamo) / entrega.
- **Metodo:** lectura de los docs 81-88 y calculo directo (integracion numerica propia de F/f por el sistema de retardo; fuerza bruta finita). Fuentes primarias: Nathanson, *Additive Number Theory: The Classical Bases* (1996), caps. 9-10; Bordignon-Johnston-Starichkova, arXiv:2207.09452; Chen 1973; Halberstam-Richert, *Sieve Methods* (1974).
- **Estado:** **[T]** transferencia polilogaritmica y cota inferior sumada sobre u (sec. 2); **[T]** geometria de la criba con umbral c <= 2/17 y balance con cero c* ~= 0.0399 (sec. 3); **[N]** LEMA UNIF-U, el resto bilineal en el rango potencia, enunciado exacto (sec. 4). No se demuestra Goldbach; no hay victoria parcial.

## 1. La descomposicion en u y que significa "uniformidad"
**Marco (docs 81/82).** `N` par grande, `k>=2` fijo, `Q=Q_N` el selector de (83.1) (producto de `k-1` primos impares distintos que no dividen a `N`, `Q<<_k (log N)^(k-1)`). Para cada primo `p` con `N-p` divisible por `Q`:
```text
N - p = Q u v,   rad(u)|Q,   (v,Q)=1,   omega(v) in {1,2}.     (89.1)
```
`u` es la parte `Q`-radical de `N-p` (solo primos de `Q`, con exponentes arbitrarios); `1 <= u <= N/Q`. La descomposicion es una PARTICION: cada `p` da un unico `u`. La incidencia completa es la suma sobre `u` (y sobre selectores `Q`) de las subincidencias; la multiplicidad exacta es `I_k = k C_k + binom(k+1,2) C_{k+1}` (doc 81).

**Que debe ser uniforme.** Para cada `u`, la maquina de Chen con modulo exterior `A = Qu`: sucesion `a_A(v) = Lambda(N-Av)`, `v < X_A = N/A`, `(v,Q)=1`; parametros `z_A = X_A^(1/10)`, `y_A = X_A^(1/3)`, `D_A = X_A^(1/2)/log^{C0} X_A`; suma de switching `B(A)`, termino principal `M_B(A)`, resto `R_B(A) = sum_s sum_{d<D_A/s} lambda_d(s) E(N;Asd,N)`. "Uniformidad en u" significa: el resto es `<<` el termino principal con constantes independientes de `u` (hasta el rango que se cubra) y el balance sigue positivo. El doc 85 sec. 4 ya lo dijo: `u=1` EVITA el problema para `C_k+C_{k+1}>0` (doc 88 lo certifica), pero lo ESCONDE para `I_k` completa.

**Dos hechos de conteo (H3).** El numero de `u <= U` con `rad(u)|Q` es
```text
N(U) = #{u <= U : rad(u)|Q} ~ (log U)^(k-1) / ((k-1)! prod_{l|Q} log l),   (89.2)
```
polilogaritmico en `U` (medido: `Q=3`: `N(10^6)=13` vs asintotico `12.58`; `Q=105`: `203` vs `127.7` — el termino principal subestima a esta escala porque la correccion de frontera es `O((log U)^(k-2))`, cociente `-> 1` como `1+O(1/log U)`). Y lo decisivo:
```text
sum_{u<=U} 1/u  ->  prod_{l|Q} l/(l-1)    (U -> infinito),                (89.3)
```
CONVERGE (medido: `Q=3`: `1.499999` vs limite `1.5`; `Q=15`: `1.874991` vs `1.875`; `Q=105`: `2.187458` vs `2.1875`; cola `< 10^-4` ya en `U=10^5`). La masa del termino principal se concentra en `u` pequenos: `u=1` ya aporta `1` de un total `<= 2.19` (para `k<=4`).

## 2. Transferencia polilogaritmica [T] y cota sumada
**Prop. 89.1 [T] (el Lema SW-Q se transfiere a `u <= (log N)^t`).** Sea `t` fijo y `A = Qu` con `u <= (log N)^t`. Entonces (a), (b), (c) de los docs 86-88 valen con `Q` reemplazado por `A` y `C0 = C0(A',k,t)`:
- Ajuste BV: `A*D_A = sqrt(AN)/log^{C0} X_A <= N^{1/2}/log^B N` porque `sqrt(A) <= (log N)^{(k-1+t)/2}` y basta `C0 >= B + (k-1+t)/2`, constante.
- Geometria: `D_A >= y_A*z_A` porque `1/2 > 13/30` y `X_A = N/A -> infinito` polinomialmente; `u_s^fix` recorre `[5/3,4]-o(1)`.
- Balance IDENTICO: `phi(A)/A = phi(Q)/Q` y `rad(A)=rad(Q)` (la serie singular local (82.6) no ve los exponentes de `u`); la constante de la cola `c_{1/10}` no cambia (`z_A = X_A^{1/10}`); `f(s_lower) -> f(5)`, `J -> 1.5238762503` uniformemente (`log X_A = log N - O(log log N)`).
- Pre-criba y resto exacto: identicos (doc 86), modulos `<= Q^2*A*D_A`, absorbidos en `C0`.

**Prop. 89.2 [T] (cota inferior sumada).** Sumando las cotas por `u` (disjuntas por la particion (89.1)) sobre `u <= U = (log N)^t`:
```text
sum_{u<=U} R_{Qu}^{P<=2} >= Bal(0) * (N*U_N*prod_{l|Q}(l-1)/(l-2)*10/(phi(Q)*log N))
                           * sum_{u<=U} 1/u * (1+o(1)) - O(#u * N/log^{A'} N),   (89.4)
```
con `Bal(0) = 0.061404` (doc 88) y `#u ~ (t log log N)^{k-1}`. El error total es `<< N/log^{A'-k+1} N`; tomando `A' > k` la cota es positiva y cuantitativa, de orden `N/log N`. Medicion (H4): a `N=10^12`, `t=2`, `sum_{u<=U} 1/u = 1.499314` (`Q=3`) y `1.869605` (`Q=15`): el `99.9%` del limite (89.3). El peso `J_k` se transfiere por el argumento del peso minimo (doc 88 (88.4), `log A + log z_A` en vez de `log Q + log z`).

**Correccion honesta.** NO hay amplificacion `(log log N)^{k-1}`: la suma reciproca (89.3) converge, asi que el total es un factor CONSTANTE `prod l/(l-1)` sobre el termino `u=1` (1.5x para `k=2`, 1.875x para `k=3`). La intuicion ingenua ("sumar sobre `u` multiplica por el conteo") es falsa: el conteo `N(U)` crece, pero el peso `1/u` lo compensa exactamente.

## 3. Donde se rompe al pasar a `u ~ N^c`
Sea `A = N^c` (`c>0` fijo), con `C0` ajustado para que quepa en BV: `C0 >= B + (c/2)(log N/log log N)`.

**(a) Geometria de la criba [T].** Con ese ajuste, `log D_A = (1/2-c) log N` y
```text
u_max(c) = (4-9c)/(1-c),   u_min(c) = (5/3-(20/3)c)/(1-c),   s_lower(c) = (5-10c)/(1-c).   (89.5)
```
La version de limite fijo exige `D_A >= y_A*z_A`, es decir `u_min(c) >= 1`, o sea **`c <= 2/17 ~= 0.1176`** (H1: en el umbral, `u_min=1`, `u_max=10/3`, `s_lower=13/3`, exactos). Para `c > 2/17` la parte alta del rango de switching (`s` cerca de `y_A`) tiene `u_s < 1` y el Teorema 9.7 no aplica ahi (mismo fenomeno que el doc 87 Prop. 87.1, ahora movido por `c`).

**(b) Balance [T].** En unidades de `|A|V_Q(z_A)`:
```text
Bal(c) = f(s_lower(c)) - (1/2) J(c) - (1/2)*0.4*e^gamma*c_{1/10}*(phi(Q)/Q)/(1-c),   (89.6)
J(c) = int_{1/10}^{1/3} F(u(b)) db/b,  u(b) = 10[(1/2-c)-b(1-c)]/(1-c).
```
El factor `1/(1-c)` de la cola sale de `|B_A|/|A| ~ 1/log X_A` contra el peso `Lambda ~ log N` (misma derivacion que doc 88 sec. 3 con `X_A` en vez de `X`). Medicion (H2): `Bal(0)=0.061404`, `Bal(0.02)=0.032934`, `Bal(0.04)=-0.000208`; el cero es
```text
c* = 0.039884,                                                           (89.7)
```
biseccion a `1e-16`. **El balance — no la geometria — es la restriccion activa** para los parametros `(1/10,1/3)`: se agota en `c* ~= 0.04`, muy antes del umbral geometrico `2/17`.

**(c) Tamano del error contra termino principal [T].** Por `u` fijo, termino principal `~ N^{1-c}/log N` (en unidades `R`) contra error BV `~ N/log^{A'} N` para todo `A'` fijo: el cociente es `(log N)^{A'-1}/N^c -> 0`. El cruce `N*` donde se igualan, `(A'-1) ln ln N* = c ln N*`, medido (H5):
```text
c=0.02: log10 N* = 1605 (A'=10);  c=0.05: 560;  c=0.10: 248;  c=2/17: 205.
```
Incluso el menor cruce es `10^{76}` (`c=2/17, A'=5`). La obstruccion es ASINTOTICA, invisible a cualquier escala computable — pero real: para `c>0` fijo y `N -> infinito`, el error BV domina al termino principal por `u`. En el rango polilog, en cambio, el cociente es `(log N)^{A'-1-t} -> infinito` (H5: `1.23e10` a `N=10^12`, `t=2`, `A'=10`).

**(d) La constante `C0` [T].** El ajuste BV exige `C0 >= B + (c/2)(log N/log log N)` (H6): para `c=0.05`, `C0 = 10.13 -> 10.21 -> 10.34` al pasar `N=10^6 -> 10^12 -> 10^24`; para `c=2/17`, `10.31 -> 10.81`. No hay constante `C0` uniforme en `N` para el rango potencia. En el rango polilog, `C0 = B + (k-1+t)/2` es constante (Prop. 89.1).

## 4. El lema que falta: LEMA UNIF-U [N]
**LEMA UNIF-U.** *Hipotesis.* `N` par grande; `A = Qu` con `u <= N^c` (`c>0` fijo), `rad(u)|Q`; parametros de Chen re-optimizados (`z_A = X_A^alpha`, `y_A = X_A^beta`, `D_A = X_A^{1/2}/(log X_A)^{C0}`). *Conclusion.* `B(A) <= M_B(A)(1+o(1)) + O(N/log^{A'} N)` con balance positivo; en la forma que realmente se necesita,
```text
sum_s sum_d lambda_d(s) * E(N;Asd,N) << (N/A)/log^2 N * o(1),          (89.8)
```
es decir, el resto bilineal es menor que el termino principal, no solo `<< N/log^{A'} N`.

**Por que no lo da ningun teorema conocido.** BV da `<< N/log^{A'} N`, independiente de `A`: para `A = N^c` eso es mayor que el termino principal `N^{1-c}/log N` por el factor `N^c/(log N)^{A'-1} -> infinito`. GRH en el peor caso por modulo da `sum_{m<=A*D_A} |E| << (A*D_A)*N^{1/2}*polylog = N*sqrt(A)*polylog`, PEOR que BV. Barban-Davenport-Halberstam es media cuadratica sobre `q` y sobre `a`, no fija `a=N`. Lo que lo daria es cancelacion de raiz cuadrada EN PROMEDIO sobre los modulos `A*m`, `sum_{m<=A*D_A} |E(N;Am,N)| << sqrt(N*A*D_A)*polylog = N^{3/4}A^{1/4}*polylog`, que bastaria hasta `c<1` — pero no es un teorema: es esencialmente la conjetura de Montgomery para residuo fijo. Tampoco lo repara re-optimizar `(alpha,beta)`: el fallo esta en el TAMANO del error, no en la geometria (sec. 3(a)-(b) solo mueven `c*`).

**Alcance honesto.** La suma sobre `u` (sec. 2) NO necesita este lema: la masa esta en `u` pequenos y la cola `u > N^{c*}` es despreciable en masa (`sum_{u>N^{c*}} 1/u ~ (log N)^{k-2}/N^{c*}`). El lema hace falta para la uniformidad en el SELECTOR `Q <= N^theta` (doc 81) y para cualquier afirmacion sobre un `u` grande especifico. La extension de la suma hasta `u <= N^{c*}` sigue el mismo argumento que la Prop. 89.2 con `Bal(c*)` en vez de `Bal(0)` (el `o(1)` es uniforme en `u <= N^{c*}` porque `log z_A, log D_A >= const*log N`).

## 5. Verificacion numerica (control finito, no evidencia asintotica)
Script `numerico/89_uniformidad_u_mayor_1.py` (solo stdlib; ~1.2 s local):
```text
PASS H0 (reproduce docs 87/88): f5_grid=0.9982417384 f5_simp=0.9982417245 (doc88: 0.9982417256);
  F4_grid=1.0216415485 F4_simp=1.0216415525 (doc87: 1.021641553); c110_simp=0.4909952010 (exacto);
  J0=1.5238762554 (doc: 1.5238762503); Bal0=0.061404 (doc88: 0.061404)
PASS H1 (geometria): umin(2/17)=1.0000000000, umax(2/17)=3.3333333333, s_lower(2/17)=4.3333333333
PASS H2 (balance): c*=0.039884, Bal(c*)=-2.78e-17; Bal(0.02)=0.032934, Bal(0.04)=-0.000208
PASS H3 (conteo/suma de u): Q=3: sum1/u -> 1.5 (cola 1e-6 en U=10^6); Q=15: -> 1.875 (cola 9e-6);
  Q=105: -> 2.1875 (cola 4.2e-5); N(U) con asintotico (89.2), cociente -> 1 como 1+O(1/log U)
PASS H4 (masa): N=10^12, t=2: sum1/u = 1.499314 (Q=3), 1.869605 (Q=15): 99.9% del limite
PASS H5 (cruce): log10 N* = 1605 (c=0.02,A'=10), 560 (0.05,10), 248 (0.10,10), 205 (2/17,10); 76 (2/17,5)
PASS H6 (C0): c=0.05: 10.13/10.21/10.34 (N=10^6/10^12/10^24); polilog: 11.5 constante
PASS
```
Lectura honesta: H0-H2 son formulas e integracion numerica de `F/f` por dos metodos independientes (malla de retardo + Simpson con formas cerradas; acuerdo `~1e-8`); H3-H4 son fuerza bruta finita; H5-H6 resuelven numericamente las ecuaciones de cruce y ajuste. El script controla formas, rangos y constantes, no el tamano asintotico (que es BV + criba lineal, analitico). Los numeros deben ser IDENTICOS local y en el PC.

## 6. Estado y siguiente paso
- **[T]** Transferencia polilog (Prop. 89.1) y cota inferior sumada sobre `u` (Prop. 89.2): la incidencia completa para `Q` fijo tiene cota `>> N/log N`; la masa esta en `u` pequenos, sin amplificacion.
- **[T]** Geometria: umbral `c <= 2/17`. **[T]** Balance: cero `c* = 0.039884`, restriccion activa para `(1/10,1/3)`.
- **[N]** LEMA UNIF-U (89.8): el resto bilineal en el rango potencia. Obstruccion: tamano del error BV contra termino principal; no la geometria.
- **Siguiente paso natural (no se empieza aqui):** (i) re-optimizar `(alpha,beta)` para empujar `c*` mas alla de `0.0399` (el balance con `(1/10,1/3)` se agota pronto); o (ii) atacar el LEMA UNIF-U por el lado del resto bilineal (estimacion Tipo II / BFI para los modulos `A*m`). Goldbach sigue abierta.
