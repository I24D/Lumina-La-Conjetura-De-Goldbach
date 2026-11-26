# 61 - Peso no negativo W(n) >= 0 con cancelacion marginal: la obstruccion de equidistribucion

## Procedencia y estado

- **Agente:** LUMINA/Muse.
- **Fecha:** 2026-11-25 ~20:05 ET. Reclamo de ESTADO.md 19:56 ET (este paso; el bucle no lo duplica).
- **Metodo:** marco del lado de los numeros, analogo de (46.1)-(46.6); construccion explicita K=2; verificacion exacta con aritmetica racional (`Fraction`) + medicion de densidades por criba (X=2*10^6, numpy).
- **Estado:** **[T]** para las Proposiciones 61.1-61.4 y las mediciones. No hay cota inferior de capas, no hay deteccion de P2, no hay demostracion de Goldbach.
- **Nota de numeracion:** existen dos docs 60 (integracion y layer cake del bucle); este es el 61, siguiente numero libre.

## 1. El marco del lado de los numeros (N61.1)-(N61.4)

Sea `W(n) = sum_{d|n} lambda_d`, con `lambda` soportada en `d` squarefree, `|lambda_d| <= 1`. El termino principal (normalizacion `phi`) y los margenes del lado de los numeros son

```text
L(lambda) := sum_d lambda_d/phi(d),                          (N61.1)
M_e(lambda) := sum_{d: e|d} lambda_d/d,                       (N61.2)
A := M_1(lambda) = sum_d lambda_d/d.                          (N61.3)
```

**Proposicion 61.1 [T].** Se tiene la identidad exacta

```text
L(lambda) = sum_e M_e(lambda)/phi(e),                        (N61.4a)
```

suma sobre `e` squarefree (incluido `e=1`, `phi(1)=1`). El termino `e=1` es `A`; el primer orden es `sum_l M_l/(l-1)` sobre primos `l`; el resto es `omega(e) >= 2`.

*Prueba.* Para `d` squarefree, `1/phi(d) = (1/d) prod_{l|d} (1 + 1/(l-1)) = (1/d) sum_{e|d} 1/phi(e)`; intercambiar las sumas da (N61.4a). QED.

El criterio del lado de los numeros, analogo de (46.4), es

```text
M_l = 0 para todo primo l.                                    (N61.4)
```

A diferencia de (46.4), aqui los margenes cero **no** fuerzan `A = 0`: los margenes son sumas de divisibilidad, no sumas de filas (`sum_l M_l = sum_d omega(d) lambda_d/d`, distinto de `A` en general). Por eso el marco de los numeros admite, en principio, termino principal de orden 1 con cancelacion exacta; la obstruccion esta en otro lado (sec. 3).

## 2. Que significa un margen: equidistribucion (Prop. 61.2) [T]

Para `A > 0`, sea `rho_l` el cociente de densidades de Dirichlet de `W` en multiplos de `l`:

```text
rho_l := lim_{x->inf} [sum_{n<=x, l|n} W(n)/n] / [sum_{n<=x} W(n)/n].
```

**Proposicion 61.2 [T].** `M_l = A (l rho_l - 1)/(l-1)`.

*Prueba.* El coeficiente de `log x` en `sum_{n<=x, l|n} W(n)/n` es `(1/l)A + M_l(1-1/l)`: si `l|d`, la suma interior sobre `m` es `log(x/d)` entera; si `l∤d`, solo `m` multiplo de `l` contribuye, `(1/l)log(x/d)`. El de `sum_{n<=x} W(n)/n` es `A`. Dividiendo se obtiene `rho_l = 1/l + M_l(1-1/l)/A`, que es la formula. QED.

**Corolario.** `M_l = 0` si y solo si `rho_l = 1/l`: equidistribucion de `W` en los multiplos de `l`. Defino la fuerza de seleccion `delta_l := 1 - l rho_l`; entonces

```text
M_l = -delta_l * A/(l-1).                                     (N61.5)
```

El margen es exactamente el defecto de equidistribucion, pesado por el termino principal.

## 3. La obstruccion (Prop. 61.3) [T]

**Proposicion 61.3 [T].** Si `W(n) = 0` cuando `l|n` (apoyo rugoso en `l`) y `A > 0`, entonces `rho_l = 0` y `M_l = -A/(l-1)` exactamente (defecto maximo). En general, seleccion `delta_l` cuesta defecto `delta_l * A/(l-1)` por (N61.5).

*Prueba.* Si `W` se anula en multiplos de `l`, el numerador de `rho_l` es cero; sustituir en la Prop. 61.2 da el valor. QED.

**Consecuencia.** La cancelacion marginal exacta (N61.4) es incompatible con cualquier seleccion en `l`: un peso con `M_l = 0` es equidistribuido mod `l`, es decir, ciego a la criba en `l`. El unico ejemplo exacto es el peso trivial `W == const` (que no criba). Si `W` selecciona `Omega(n) <= 2` (primos y P2), entonces `rho_l -> 0` (medido: `rho_2 = 0.251` frente a `1/2` a `X=2*10^6`), luego `delta_l -> 1` y `M_l` es casi `-A/(l-1)`: la seleccion de P2 exige el defecto maximo.

**Advertencia sobre `L(lambda)`.** El termino principal fisico de `sum_p W(N-p)` es `L_N(lambda) := sum_{(d,N)=1} lambda_d/phi(d)`, no `L(lambda)`: la correccion por `(d,N)>1` es `O(1)` y depende de `N`. Para el indicador rugoso, `L = 0` exactamente (el factor `l=2` anula el producto) pero `L_N(N par) = 0.3125` a `z=7` (medido, sec. 5). La identidad (N61.4a) es exacta; su interpretacion fisica exige la restriccion `(d,N)=1`.

## 4. Construccion: la familia K=2 (Prop. 61.4) [T]

Para `z` fijo y `t` en `[0,1]`, sea `lambda_d = mu(d) t^{omega(d)}` para `d` squarefree con primos `<= z` y `omega(d) <= 2`; `0` en otro caso. Soporte `d <= z^2`.

**Proposicion 61.4 [T].** Con `r = omega_z(n) = #{l|n: l<=z}`, se tiene `W(n) = 1 - t r + t^2 C(r,2)`, y:

- `W(n) >= 0` para todo `n` y todo `t` en `[0,1]`: `W` depende solo de `r` entero; `W(r) = (1 - t r/2)^2 + t^2 r(r-2)/4 >= 0` para `r >= 2`, y `W(0)=1`, `W(1)=1-t >= 0`.
- `W(p) = 1-t > 0` para `p <= z` (si `t<1`); `W(p) = 1` para `p > z`.
- `A = 1 - t S_z + (t^2/2)(S_z^2 - S_z^{(2)})`, `M_l = (-t + t^2(S_z - 1/l))/l`, con `S_z = sum_{l<=z} 1/l`, `S_z^{(2)} = sum_{l<=z} 1/l^2`.
- Con `t = 1/S_z`: `A = (1/2)(1 - S_z^{(2)}/S_z^2) = Theta(1)` y `M_l = -1/(l^2 S_z^2)` exacto.
- Nivel: soporte `d <= z^2`; con `z = N^{1/10}`, `d <= N^{1/5}`; Bombieri-Vinogradov da nivel `1/2` ([T-lit], estandar: `|lambda_d| <= 1`).

**Pero:** `W(n) = 1` cuando `omega_z(n) = 0`, incluidos P3, P4, ... con factores grandes. `W` es ciega a los factores grandes: **no detecta P2**. La construccion cumple (a) positivo en primos, (b) margenes pequenos, (c) termino principal de orden 1, pero no la seleccion.

**El trade-off medido (tabla, `delta_l = |M_l|(l-1)/A`):**

| z | t | A | delta_2 | delta_3 | delta_5 |
|---|---|---|---|---|---|
| 7 | 1/10 | 0.887 | 0.053 | 0.069 | 0.081 |
| 7 | 1/2 | 0.532 | 0.311 | 0.362 | 0.385 |
| 7 | 1/S_z | 0.348 | 0.520 | 0.462 | 0.333 |
| 11 | 1/10 | 0.879 | 0.053 | 0.069 | 0.081 |
| 11 | 1/S_z | 0.366 | 0.425 | 0.378 | 0.272 |

El defecto es `delta ~ t` para `t` pequeno: la seleccion ES el defecto. La formula `M_l = -1/(l^2 S_z^2)` no significa defecto pequeno: `S_z = O(1)` (`log log z`) para `z` realista, asi que `delta_2 = 0.52` a `z=7`.

**Familia `W_t` sin truncar (contraste).** `lambda_d = mu(d) t^{omega(d)}` para todo `d` `z`-suave: `W(n) = (1-t)^{omega_z(n)} >= 0`; `M_l = -tA/(l-t)`, `delta_l = t(l-1)/(l-t) ~ t`. Mismo trade-off lineal; sin nivel limpio (soporte hasta el primorial).

## 5. Verificacion (script, PASS en el PC)

`python numerico/61_peso_no_negativo_marginal.py` — PASS en el PC (numeros identicos a la corrida local; salida completa en `salida_doc61.txt`).

- Exacto (`Fraction`): identidad `L = sum M_e/phi(e)` en todos los casos; `M_l = -A/(l-1)` exacto para el rugoso `z=5,7,11`; `M_l = -tA/(l-t)` exacto para `W_t`; formulas de `A` y `M_l` de K=2 exactas; `W >= 0` y `W(p) > 0` en `[1,3000]`; soporte `max d <= z^2`.
- Densidades (criba `X=2*10^6`): rugoso `z=7`: `rho_l = 0.0000` exacto; `Omega<=2`: `rho_2=0.251`, `rho_3=0.166`, `rho_5=0.099`, `rho_7=0.070` (todos ~mitad de `1/l`); K=2 `t=1/S_7`: medido vs predicho por la Prop. 61.2 con diferencia `<= 0.040` (error `O(1/log X)`).
- `L_N` fisico: K=2 `z=7, t=1/10`: `L=0.82`, `L_N=0.91` (N par), `0.93` ({2,5}), `1.00` ({2,3,5,7}); rugoso: `L=0`, `L_N(N par)=0.3125`. `L_N` se mantiene `Theta(1)` y positivo.

## 6. Veredicto

**El cambio de marco SI produce un peso no negativo con cancelacion marginal (relajada) y termino principal de orden 1 (familia K=2): las tres condiciones pueden coexistir. Pero la obstruccion demuestra que la coexistencia se compra renunciando a la seleccion: cancelacion es equidistribucion, seleccion es defecto, con el trade-off lineal exacto `delta*A/(l-1)`.** La version exacta del criterio es imposible para un peso que seleccione P2; la version relajada existe pero es debilmente selectora por construccion. El siguiente paso no es mas ingenieria de margenes sino el mecanismo de seleccion (switching, C10, frentes del bucle). Goldbach sigue abierta.

## 7. Que sigue

1. Bajar el soporte a `N^eps` (doc 60 #1) con el criterio relajado: arregla escala y well-factorabilidad; el cero-en-primos de (60.1) sigue exigiendo el componente P1 que la familia K=2 ya aporta (`W(p) = 1-t > 0`).
2. Frentes del bucle: C10/W(d), uniones, Pintz-Althofer, Motohashi (doc 57 sec. 4), `b_{P,l}` (doc 53).
3. Dispersion uniforme en `a` (doc 59 #1) sigue abierta como problema independiente.

## Verificacion reproducible

```text
python numerico/61_peso_no_negativo_marginal.py
```

PASS en el PC: identidades exactas de malla (rugoso, trivial, `W_t`, K=2 en `z=5,7,11,13` y `t=1/2, 1/10, 1/S_z`); `W>=0` y `W(p)>0` en `[1,3000]`; densidades `rho_l` por criba a `X=2*10^6`; `L_N` fisico para varios conjuntos de primos de `N`.

## Alcance honesto

- Lo que este doc SI dice: el marco (N61.1)-(N61.4) y la identidad exacta; la interpretacion `M_l = A(l rho_l - 1)/(l-1)`; la obstruccion seleccion-defecto con el defecto maximo `-A/(l-1)` para pesos rugosos; la familia K=2 con `W >= 0`, `W(p) > 0`, termino principal `Theta(1)`, nivel `1/2`, y el trade-off lineal medido; la advertencia de que `L` no es el termino fisico sin la restriccion `(d,N)=1`.
- Lo que NO dice: no hay cota inferior de capas, no hay deteccion de P2, no hay demostracion de Goldbach.

## Memoria compartida

Propuesta pendiente: el MCP canonico `lumina_shared_memory` / `lumina_propose` no esta expuesto. No se afirma persistencia remota.
