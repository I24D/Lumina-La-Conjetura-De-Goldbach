# 90 - Lema UNIF-U etapa 1: por que el resto por-u domina
## Procedencia y estado
- **Agente:** LUMINA/Muse. **Fecha:** 2026-11-26 20:55 ET (reclamo, reutilizado del intento anterior) / entrega.
- **Metodo:** calculo directo (integracion numerica propia de F/f por el sistema de retardo, como doc 89; bisecciones). Fuentes: docs 85-89; Nathanson caps. 9-10; Bordignon-Johnston-Starichkova arXiv:2207.09452.
- **Estado:** **[T]** el cociente resto/termino principal medido en malla (sec. 1); **[T]** los dos ingredientes que faltan, enunciados exactos (sec. 2); **[T]** umbral de la cancelacion de raiz cuadrada c<1/5, precision del doc 89 (sec. 3); **[N]** LEMA UNIF-U sigue sin tenerse; el mapa de rutas queda pendiente de medir (sec. 4). No se demuestra Goldbach; no hay victoria parcial.

## 1. El cociente resto/termino principal
**Marco (doc 89 sec. 3(c)).** Por `u` fijo con `A = Qu = N^c`, termino principal `~ N^{1-c}/log N` (en unidades `R`) contra error BV `~ N/log^{A'} N` para todo `A'` fijo. El cociente es
```text
R(c,N;A') = resto/principal = N^c / (ln N)^{A'-1}.                  (90.1)
```
**Prop. 90.1 [T] (el cociente primero CAE).** `d/d(ln N) [ln R] = c - (A'-1)/ln N`: el cociente DECRECE mientras `ln N < (A'-1)/c` (el ahorro logaritmico gana al principio) y solo despues crece. Minimo en `ln N = (A'-1)/c` (para `c=0.02, A'=10`: `N ~ 10^195`). A toda escala computable el resto no solo es despreciable, sino que se hace MAS pequeno al crecer `N`. Medicion (H1, `log10 R`):
```text
c        A'   N=1e6   N=1e12   N=1e24   N=1e76   log10 N* (cruce R=1)
0.01     10   -10.20  -12.85   -15.44   -19.43   3517.6
0.02     10   -10.14  -12.73   -15.20   -18.67   1605.5
0.039884 10   -10.02  -12.49   -14.72   -17.16    727.5
0.05     10    -9.96  -12.37   -14.48   -16.39    559.9
0.10     10    -9.66  -11.77   -13.28   -12.59    248.1
2/17      5    -3.86   -4.35    -4.15    -0.03     76.3
```
Los cruces reproducen el doc 89 H5 (1605, 560, 248, 205, 76). El menor cruce es `10^76` (`c=2/17, A'=5`): la dominacion es ASINTOTICA, invisible a cualquier escala computable — pero real: para `c>0` fijo y `N -> infinito`, `R -> infinito` para todo `A'` fijo.

## 2. Los dos ingredientes que faltan
El lema (89.8) pide `sum_s sum_d lambda_d(s) E(N;Asd,N) << (N/A)/log^2 N * o(1)`: el resto menor que el termino principal. BV da `<< N/log^{A'} N`, ABSOLUTO (independiente de `A`). Faltan dos cosas distintas, y hay que nombrarlas por separado:

**(I1) Nivel [T].** Con `C0` FIJO, los modulos llegan a `A*D_A = N^{1/2+c/2}/polylog`: hace falta nivel de distribucion `theta = 1/2 + c/2`, mas alla de BV (`1/2`). Con `C0` CRECIENTE (el ajuste del doc 89, `C0 = 10 + (c/2)(ln N/ln ln N)`), los modulos caben en BV — pero el exponente efectivo de criba cae de `(1-c)/2` a `1/2-c` (H4): `e_eff(c,N) = 1/2 - c - 10*ln ln N/ln N -> 1/2-c`. Ese es el precio: el nivel y el balance son LA MISMA moneda — comprar aplicabilidad de BV cuesta exactamente el balance que se agota en `c* = 0.039884` (doc 89). Medicion (H4): `D_A > 1` solo a partir de `log10 N0 = 39.1 (c=0), 41.2 (0.02), 44.7 (0.05), 51.9 (0.10), 55.0 (2/17)`: el ajuste empuja el encendido de la criba de `~10^39` a `~10^55`.

**(I2) Pequenez relativa [T].** Aun con los modulos dentro de BV, el error es absoluto. Para que `N/log^{A'} N <= (N/A)/log^2 N` hace falta
```text
A'(N) >= 2 + c * ln N / ln ln N,                                     (90.2)
```
es decir, un ahorro logaritmico que CRECE con `N`. Medicion (H2): para `c=0.05`, `A'_nec = 2.42 (N=10^12) -> 5.83 (10^205) -> 24.49 (10^1605)`; para `c=0.10`, `46.99` a `10^1605`. BV da todo `A'` FIJO, nunca uno creciente. **Enunciado exacto del ingrediente:** hace falta una cota con el factor `1/A` — `sum_{m<=A*D_A} |E(N;Am,N)| << (N/A)/log^2 N * o(1)` — no una cota absoluta mas fuerte.

## 3. Cancelacion de raiz cuadrada en promedio y su umbral
**Prop. 90.2 [T] (umbral c<1/5).** La cancelacion de raiz cuadrada en promedio sobre los modulos `A*m` daria `sum |E| << sqrt(N*A*D_A)*polylog = N^{3/4}A^{1/4}*polylog`. Su cociente contra el termino principal `(N/A)/log^2 N` es
```text
R_sq(c,N) = N^{5c/4-1/4} * (ln N)^2,                                (90.3)
```
que `-> 0` si y solo si `c < 1/5` (en `c=1/5` el factor `(ln N)^2` lo mata: nunca cruza). Medicion (H3, `log10 R_sq`): `c=0.05`: `0.63 (N=1e12) -> -14.03 (1e100)`, cruce `log10 N* = 17.0`; `c=0.10`: cruce `29.3`; `c=0.15`: `70.8`; `c=0.19`: `488.1`; `c>=0.20`: nunca.
**Precision del doc 89.** El doc 89 sec. 4 dijo que la cancelacion de raiz cuadrada "bastaria hasta c<1": eso vale para el objetivo DEBIL en forma BV (error `<< N/log^{A'} N`, donde `N^{3/4+c/4} << N` si `c<1`). Para el objetivo REAL (89.8) — error menor que el termino principal — el umbral es **`c < 1/5`**, no `c<1`. La cancelacion de raiz cuadrada, aunque se tuviera, solo cubre `c<1/5`: ni siquiera llega al umbral geometrico `2/17 ~= 0.1176` con margen comodo (lo cubre, pero el cruce a `c=0.10` esta en `N ~ 10^29`).

## 4. Mapa de rutas (etapa 2: una por una)
- **BV con pesos bien factorizables:** los `lambda_d(s)` son `mu(d)` o 0 (doc 86), no triplemente bien factorizables (doc 18, C7 para `mu(d) log d`); si la suma sobre `s` factoriza de otro modo es otra pregunta. PENDIENTE DE MEDIR.
- **BFI (nivel 4/7):** exige pesos bien factorizables y clase fija `a=N`; el doc 59 midio cobertura de solo 5.4% de la masa M1 para el peso relajado — para SW-Q la cobertura es PENDIENTE DE MEDIR.
- **Maynard III (clases uniformes):** umbral 1/10 agudo (doc 75); la ventana no cubre modulos `A=N^c` con `c>0` fijo en el rango que hace falta (conjetura preliminar: no cubre). PENDIENTE DE MEDIR.
- **Dispersion Tipo II:** la gran criba capa el Tipo II en 1/2 (doc 59); la dispersion que haria falta es para modulos `A*m` con `A=N^c`, fuera de lo conocido. PENDIENTE DE MEDIR.
- **Densidad de ceros:** el doc 76 cerro la ruta de densidad para el momento mixto (el sup necesitaria `Theta<1/2`, refutado); para UNIF-U la densidad daria a lo sumo nivel `1/2+o(1)` y no la pequenez relativa (I2). PENDIENTE DE MEDIR.

## 5. Verificacion numerica (control finito, no evidencia asintotica)
Script `numerico/90_lema_unif_u_etapa1.py` (solo stdlib; ~1.3 s local):
```text
PASS H0 (reproduce doc 89): f5_grid=0.9982417384 f5_simp=0.9982417245 (doc88: 0.9982417256);
  F4_grid=1.0216415485 F4_simp=1.0216415525 (doc87: 1.021641553); c110_simp=0.4909952010 (exacto);
  J0=1.5238762554 (doc: 1.5238762503); Bal0=0.061404 (doc88: 0.061404); c*=0.039884 (doc89)
PASS H1 (cociente): tabla sec. 1; cruces 3517.6 / 1605.5 / 727.5 / 559.9 / 248.1 / 204.5 (A'=10);
  76.3 (c=2/17, A'=5); minimo del cociente en ln N = (A'-1)/c
PASS H2 (A' necesario): c=0.05: 2.42 -> 24.49 (N=1e12 -> 1e1605); c=0.10: 46.99 a 1e1605; crece sin cota
PASS H3 (raiz cuadrada): umbral c<1/5; cruces 17.0 / 29.3 / 70.8 / 488.1 (c=0.05/0.10/0.15/0.19);
  c>=0.20: nunca
PASS H4 (precio C0): e_eff -> 1/2-c; log10 N0 = 39.1 / 41.2 / 44.7 / 51.9 / 55.0
  (c=0 / 0.02 / 0.05 / 0.10 / 2/17)
PASS
```
Lectura honesta: H0 es integracion numerica de `F/f` por dos metodos independientes (malla de retardo + Simpson con formas cerradas; acuerdo `~1e-8`); H1-H4 son formulas y bisecciones. El script controla formas, rangos y constantes, no el tamano asintotico (que es BV + criba lineal, analitico). Los numeros deben ser IDENTICOS local y en el PC.

## 6. Estado y siguiente paso
- **[T]** El cociente resto/principal (90.1) medido en malla: cae hasta `ln N = (A'-1)/c` y cruza 1 solo en `N* >= 10^76`.
- **[T]** Los dos ingredientes: (I1) nivel `theta = 1/2+c/2`, o `C0` creciente que cuesta `c/2` en el exponente y empuja el encendido a `N0 ~ 10^41-10^55`; (I2) pequenez relativa: hace falta `A'(N) >= 2 + c ln N/ln ln N`, es decir, una cota con el factor `1/A`.
- **[T]** La cancelacion de raiz cuadrada en promedio bastaria solo hasta `c<1/5` (precision del doc 89).
- **[N]** LEMA UNIF-U sigue sin tenerse. **Siguiente paso (etapa 2, no se empieza aqui):** examinar las cinco rutas de la sec. 4 una por una y decidir para cada una si da el lema, si no lo da, o si queda abierta con el sub-lema exacto que falta. Goldbach sigue abierta.
