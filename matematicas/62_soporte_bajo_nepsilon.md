# 62 - Bajar el soporte a N^eps con el criterio relajado: escala y well-factorabilidad si, seleccion no

## Procedencia y estado

- **Agente:** LUMINA/Muse.
- **Fecha:** 2026-11-25 ~20:20 ET. Reclamo de ESTADO.md 20:09 ET (este paso; el bucle no lo duplica).
- **Metodo:** familia K=2 del doc 61 con `z = N^(eps/2)` (soporte `d <= N^eps`); medicion exacta con aritmetica racional (`Fraction`) donde es barato y `math.fsum` con una sola criba para el piso del defecto; criba `X=2*10^6` (numpy) para la seleccion; prueba de la obstruccion triple para la well-factorabilidad.
- **Estado:** **[T]** para las mediciones exactas, el piso del defecto, la Proposicion 62.1 (umbral de well-factorabilidad) y la seleccion medida; **[T-lit]** para BV, BFI (Maynard Def. 1, Teor. A, via doc 58 sec. 6) y Jurkat-Richert `f(s)>0` para `s>2` (via doc 60 sec. 1). No hay cota inferior de capas, no hay deteccion de P2, no hay demostracion de Goldbach.
- **Nota de numeracion:** existen dos docs 60 (integracion y layer cake del bucle); este es el 62, siguiente numero libre tras el 61.

## 1. La construccion: K=2 con z = N^(eps/2)

Se toma la familia K=2 del doc 61 (Prop. 61.4), `lambda_d = mu(d) t^{omega(d)}` para `d` squarefree con primos `<= z` y `omega(d) <= 2`, con el soporte bajado a `d <= N^eps`, es decir `z = N^(eps/2)`. El componente P1 que faltaba en el doc 60 ya lo aporta K=2: `W(p) = 1-t > 0` para `t<1`. Parrilla medida: `eps` en `{1/5, 1/10, 1/20, 1/50}` (el `1/5` es la linea base del doc 61) y `N` en `{10^7, 10^12, 10^24, 10^40}`. El caso principal es `eps = 1/10`: soporte `N^(1/10)` (la mitad del exponente del doc 61, por debajo del nivel BV con margen) y `z = N^(1/20)` todavia `>= 2` para `N >= 10^7`. El borde `t=1` se excluye: ahi `W(p) = 0` y muere el componente P1.

**Degeneracion a N finito (medido, parte 2 del script).** A `N=10^12`, `eps<=1/20` da `z=1`: sin primos, `lambda = delta_1`, `W == 1` (regimen vacio). `eps=1/10` da `z=3`: solo los primos `{2,3}`. El regimen `t=1/S_z` exige `S_z >= 1`, es decir `z >= 5`, es decir `eps >= 0.116` a `N=10^12`: el `eps=1/10` "bajo" es INVIABLE para `t=1/S_z` a `N=10^12`. A `N=10^7` el vacio empieza en `eps<=1/20` y a `N=10^24` en `eps<=1/50`. El regimen `N^eps` con `eps` pequeno es vacio o degenerado a N accesible: ese es el primer precio, antes de medir nada mas.

## 2. Termino principal, defecto y trade-off (medido exacto)

Formulas del doc 61 (Prop. 61.4), verificadas exactas (`Fraction`) en `z=2,3,5,7` y `t=1/2,1/10`: `A = 1 - t S_z + (t^2/2)(S_z^2 - S_z^{(2)})`, `M_l = (-t + t^2(S_z - 1/l))/l`, identidad `L = sum M_e/phi(e)`, `W>=0` en `[1,3000]`, `W(p)>0`, soporte `max d <= z^2`. Tabla (N,eps), exacto salvo `S_z` por `fsum` en `z` grande:

```text
N      eps    z    |P<=z|  Sz     t=1/Sz?  A(t=1/2) d2   d3    A(1/Sz) d2(1/Sz)
10^7   1/5    5    3      1.033  si        0.5667   0.324 0.382  0.3122   0.750
10^7   1/10   2    1      0.500  NO        0.7500   0.333 0.000  --       --
10^12  1/5    15   6      1.344  si        0.4993   0.289 0.330  0.3794   0.365
10^12  1/10   3    2      0.833  NO        0.6250   0.333 0.400  --       --
10^24  1/5    251  54     1.984  si        0.4436   0.146 0.131  0.4426   0.144
10^24  1/10   15   6      1.344  si        0.4993   0.289 0.330  0.3794   0.365
10^40  1/5    10000 1229  2.483  si        0.4726   0.004 0.053  0.4633   0.088
10^40  1/10   100  25     1.803  si        0.4486   0.194 0.197  0.4307   0.179
10^40  1/20   10   4      1.176  si        0.5321   0.311 0.362  0.3477   0.520
```

Lecturas:

- **El defecto lo fija `t`, no `z`.** A `t=1/2`, `delta_2 ~ 0.29-0.33` en todos los casos no degenerados; a `t=1/10`, `delta_2 ~ 0.053`. Bajar `z` no reduce `delta`: el trade-off del doc 61 (`delta ~ t`, `M_l = -delta_l*A/(l-1)`) se mantiene intacto en el regimen de soporte bajo.
- **El `delta_2 = 0.004` a `z=10000, t=1/2` es accidental.** `M_2 = 0` cuando `S_z = 1/t + 1/2 = 2.5`, y `S_z(10000) = 2.483` cae ahi por casualidad; `delta_3 = 0.053` en ese mismo caso muestra que no hay mejora uniforme: cada `l` tiene su propio cruce por cero (`S_z = 1/t + 1/l`) y no se anulan a la vez.
- **El termino principal fisico es `Theta(1)` a cualquier soporte.** Medido exacto: `L_N(N par)` es `1.00, 0.75, 0.66, 0.60` para `z=2,3,5,7` a `t=1/2` (y `0.95, 0.91` a `t=1/10`); `L_N((N,30)=1) = 1.00` para `z<=5`. La escala del doc 60 (O2, `N^(4/5)`) era un problema del peso de pares y NO se aplica a K=2: `A = Theta(1)` siempre. Bajar el soporte no "arregla" la escala de K=2 porque nunca estuvo rota.
- **El defecto en el regimen `t=1/S_z` es PEOR a `z` pequeno.** A `z=5` (`N=10^7, eps=1/5`), `delta_2 = 0.750`: `S_z` apenas supera 1 y la correccion `S_z^{(2)}/S_z^2` encoge `A` a `0.31`. El regimen `t=1/S_z` solo mejora con `z` grande (sec. 3).

## 3. Piso del defecto: defecto pequeno y soporte bajo son incompatibles a N finito

En el regimen `t=1/S_z`, `M_l = -1/(l^2 S_z^2)` exacto (doc 61) y `A = (1/2)(1 - S_z^{(2)}/S_z^2)`, luego `delta_2 = 1/(2 S_z^2 (1 - S_z^{(2)}/S_z^2)) ~= 1/(2 S_z^2)`. Busqueda binaria exacta (una sola criba hasta 200000 + sumas prefijas):

```text
delta_2 <= 0.1: z_min = 2777 (S_z = 2.3350, delta_2 = 0.099997; delta_2(2776) = 0.100031)
  eps necesario: 0.574 a N=10^12, 0.287 a N=10^24, 0.172 a N=10^40
delta_2 <= 0.01: exige S_z >= sqrt(50) = 7.07, log log z >= 6.81, z >= e^903: inalcanzable
```

(Mi primera estimacion a mano, `z~1400`, era optimista: la correccion `S_z^{(2)}/S_z^2` pesa a `S_z` pequeno; el valor exacto es `2777`.) Conclusion: **defecto pequeno y soporte bajo son incompatibles a N realista.** Incluso a `N=10^40`, `delta_2 <= 0.1` exige `eps >= 0.172`, y a `N=10^12` exige `eps >= 0.574`, que ya no es soporte bajo. El `delta -> 0` asintotico (`t=1/S_z`, `N -> infinito`, `eps` fijo) existe, pero el cruce esta astronomicamente lejos: `delta_2 <= 0.01` exige `z >= e^903`. A `N <= 10^12` el defecto es `>= 0.29` en todo caso no degenerado de soporte bajo.

## 4. Well-factorabilidad: umbral Q* ~ z^4 = N^(2 eps) (Proposicion 62.1) [T]

**Proposicion 62.1 [T].** Sea `lambda` la familia K=2 con corte `z` (doc 61, Prop. 61.4), `t` en `(0,1]`. Sea `p_(1) > p_(2) > p_(3)` los tres primos mayores `<= z`, `D_max = p_(1) p_(2)` el maximo soporte y `M(z) = p_(2) p_(3)`. Entonces `lambda` es well-factorable (Maynard Def. 1, via doc 58 sec. 6) de nivel `Q` si y solo si `Q >= Q*` con

```text
max(z^2, M(z)^2) <= Q* <= D_max^2,
```

es decir `Q* = Theta(z^4) = Theta(N^(2 eps))`. En particular existe un umbral finito, a diferencia del peso de pares del doc 58 (Prop. 58.1: no well-factorable a ningun nivel `Q <= N`).

*Prueba.* (Necesidad, condicion prima.) Si `p > max(Q1,Q2)` para un primo `p <= z`, entonces `(gamma^(1)*gamma^(2))(p) = 0 != lambda_p = -t`; luego `max(Q1,Q2) >= z` para la factorizacion simetrica, es decir `Q >= z^2`.
(Necesidad, obstruccion triple.) Si hay tres primos `p,q,r <= z` con `pq, pr, qr > max(Q1,Q2)`, entonces `gamma^(1)` y `gamma^(2)` se anulan en los semiprimos y las ecuaciones por pares dan, con `a_x = gamma^(1)(x)`, `b_x = gamma^(2)(x)` y `a_1 = b_1 = 1` (pues `|gamma| <= 1` y `lambda_1 = 1`): `a_p + b_p = -t` y `a_p b_q + a_q b_p = t^2`. Sustituyendo, `a_q = f(a_p)` con `f(x) = -t(t+x)/(t+2x)`; las tres ecuaciones fuerzan `f(f(a_p)) = f(a_p)`, es decir `f(a_p)` es punto fijo de `f`: `2x^2 + 2tx + t^2 = 0`, discriminante `-4t^2 < 0`. Imposible. (El caso `a_p = -t/2` tampoco resuelve la ecuacion por pares.) Luego no hay descomposicion si existen tres primos con productos por pares `> max(Q1,Q2)`; para la factorizacion simetrica eso es `Q < M(z)^2`.
(Suficiencia.) Si `Q >= D_max^2`, para toda factorizacion `Q = Q1 Q2` se tiene `max(Q1,Q2) >= sqrt(Q) >= D_max`, y la descomposicion trivial (`gamma = lambda` del lado mayor, `delta_1` del otro) funciona, con `|lambda| <= 1`. QED.

Medido (intervalo del umbral; `Q*/z^4` se estrecha hacia 1 con `z`):

```text
z    Q*min      Q*max     Q*/z^4
2    4          4         0.25
3    9          36        0.11-0.44
5    36         225       0.06-0.36
7    225        1225      0.09-0.51
100  5.46e7     7.45e7    0.55-0.75
251  3.32e9     3.66e9    0.84-0.92
```

Ejemplo explicito verificado: `z=100, Q=10^6`: triple `{83,89,97}`, productos `[7387, 8051, 8633]`, todos `> sqrt(Q)=1000`: la obstruccion se aplica. La verificacion numerica del argumento (`f`, `g(x)=f(f(x))-f(x) != 0`, numerador `2f^2+2tf+t^2 > 0` en 6 puntos, discriminante `-1 < 0` a `t=1/2`) esta en el script.

**Pero el arreglo es vacuo para el nivel de la criba.** BFI (Maynard Teor. A, via doc 58) da error controlado hasta `Q` si `lambda` es well-factorable de nivel `Q <= N^(4/7)`: aqui `Q* ~ N^(2 eps)`, asi que BFI aplica cuando `2 eps <= 4/7` y superaria a BV (`1/2`) solo cuando `2 eps > 1/2`, es decir `eps > 1/4`. Para todo soporte bajo (`eps <= 1/4`), BV ya cubre el soporte (`d <= N^eps < N^(1/2)`): el error `sum_d lambda_d E_d` esta controlado por BV sin necesidad de BFI. La well-factorabilidad se arregla en sentido literal (umbral finito `N^(2 eps)` frente a ninguno), pero no compra nivel: BFI solo ayudaria con `eps > 1/4`, que ya no es soporte bajo. El valor de `s` con nivel BV es `s = log D/log z = 1/eps` (`5, 10, 20` para `eps=1/5,1/10,1/20`), todo `> 2`, luego `f(s) > 0` [T-lit, Jurkat-Richert via doc 60]: la criba es "eficiente" precisamente porque casi no criba (sec. 5).

## 5. Seleccion: lo que paga el soporte bajo (medido por criba, X=2*10^6)

`W(n) = 1 - t r + t^2 r(r-1)/2` con `r = omega_z(n)`. Medido por criba (numpy):

```text
z  t    frac r=0  E[W]  E[W|Om=1] E[W|Om=2] E[W|Om=3] E[W|Om>=4]  P2: frac r=0  ratio P2/P3
2  0.5  0.500     0.750  1.000     0.904     0.789     0.624        0.807        1.145
2  0.1  0.500     0.950  1.000     0.981     0.958     0.925        0.807        1.024
3  0.5  0.333     0.625  1.000     0.837     0.661     0.456        0.675        1.268
3  0.1  0.333     0.918  1.000     0.967     0.930     0.878        0.675        1.040
5  0.5  0.267     0.567  1.000     0.796     0.588     0.389        0.591        1.353
5  0.1  0.267     0.900  1.000     0.959     0.913     0.852        0.591        1.051
7  0.5  0.229     0.532  1.000     0.765     0.538     0.355        0.530        1.421
7  0.1  0.229     0.887  1.000     0.953     0.901     0.834        0.530        1.058
```

Lecturas:

- **`E[W|Om=1] = 1.000` siempre, y los P2 con ambos factores `> z` tambien reciben `W = 1` exacto.** El peso es ciego a la distincion que importa para Goldbach: un P2 con factores grandes y un primo reciben el mismo peso. Los P2 con un factor `<= z` reciben `1-t` (`0.5` a `t=1/2`): el peso castiga justo a los P2 que ve.
- **El cociente P2/P3 es `1.42` en el mejor caso (`z=7, t=1/2`) y tiende a `1.02` al bajar `z` o `t`.** La contaminacion P3+ es el `58-96%` de la senal P2. Para aislar P2 haria falta un cociente `>> 1` (o `E[W|Om>=3] ~ 0`); aqui `E[W|Om=3] = 0.54-0.96` y `E[W|Om>=4] = 0.36-0.93`. El peso ordena debilmente `Om=1 > Om=2 > Om=3 > Om>=4`, pero no aisla P2.
- **Al bajar `z`, la seleccion se aplana y la fraccion con `W=1` crece** (`0.229` a `z=7` -> `0.500` a `z=2`): la criba cuenta mas numeros con peso completo porque ve menos primos. El defecto `delta` no mide la seleccion: `delta` lo fija `t` (sec. 2) mientras que la seleccion la fija `z` (cuantos primos ve la criba). Bajar `z` aplana el ordenamiento en `Omega` sin reducir `delta`: se pierde seleccion sin ganar cancelacion.

## 6. Veredicto

**El soporte bajo con el criterio relajado NO da un peso usable para P2.** Por partes:

- **Escala: SI, pero vacuo.** El conteo pesado esta a escala completa `Theta(N/log N)` con `L_N = Theta(1)` (`0.60-1.00`) a cualquier soporte. Pero ya estaba a escala completa con el soporte del doc 61: bajar `eps` no arregla la escala de K=2 porque nunca estuvo rota (la O2 del doc 60 era un problema del peso de pares, no de K=2).
- **Well-factorabilidad: SI, pero vacuo.** El umbral finito `Q* ~ N^(2 eps)` existe y esta demostrado (Prop. 62.1), frente a la obstruccion total de la Prop. 58.1 para el peso de pares. Pero no compra nivel: BV ya cubre el soporte pequeno, y BFI solo superaria a BV con `eps > 1/4`, que ya no es soporte bajo.
- **Seleccion: NO, y el precio es total.** A N finito la criba degenera a 1-2 primos (`z=3` a `N=10^12, eps=1/10`; vacio para `eps<=1/20`), el regimen `t=1/S_z` es inviable (`S_z<1`), el defecto se queda en `O(t) ~ 0.3`, y el peso no distingue un P2 de un primo ni de un P3 (cociente P2/P3 `<= 1.42`, contaminacion `>= 58%`). El componente P1 `W(p) = 1-t > 0` esta presente (la pieza que faltaba en el doc 60), pero no ayuda a detectar P2: primos y P2 (factores `> z`) reciben el mismo peso `1`.
- **En una frase:** bajar el soporte arregla dos problemas tecnicos que K=2 no tenia (escala y well-factorabilidad) y cobra el unico que importa (seleccion). El peso K=2 de soporte bajo es un peso de criba valido para numeros rugosos, no un detector de P2. La linea de los docs 54-61 queda cerrada como ruta a P2 por este camino: el siguiente paso no es mas ingenieria de margenes ni de soporte, sino el mecanismo de seleccion. Goldbach sigue abierta.

## 7. Que sigue

1. **Mecanismo de seleccion** (switching, C10/peso movil `W(d)`, frentes del bucle): lo que falta no es un peso con mejores margenes, sino uno que seleccione `Omega <= 2` sin pagar el defecto maximo de la Prop. 61.3.
2. **La ventana `eps` en `(1/4, 2/7]`**, donde BFI si superaria a BV (`Q* ~ N^(2 eps)`, `s = 1/eps` en `3.5-4`): no es soporte bajo y el defecto ahi sigue siendo `O(0.2)`; queda anotada como ventana sin explorar, no como promesa.
3. **Frentes del bucle:** C10/W(d), uniones `C_k+C_{k+1}`, variante Lema 4 Pintz-Althofer, Motohashi (doc 57 sec. 4), `b_{P,l}` (doc 53). Dispersion uniforme en `a` (doc 59 #1) sigue abierta como problema independiente.

## Verificacion reproducible

```text
python numerico/62_soporte_bajo_nepsilon.py
```

PASS local (43 chequeos OK): formulas exactas K=2 en `z=2,3,5,7` (`Fraction`); identidad `L = sum M_e/phi(e)`; `W>=0`, `W(p)>0`, soporte `<= z^2`; tabla (N,eps) con `z` exacto por aritmetica entera; piso del defecto por busqueda binaria con una sola criba; obstruccion triple (discriminante `-1`, `g(x) != 0` en 6 puntos, ejemplo `z=100, Q=10^6`); suficiencia `Q >= D_max^2`; seleccion por criba a `X=2*10^6`. La corrida en el PC debe dar numeros identicos (salida completa en `salida_doc62.txt`).

## Alcance honesto

- Lo que este doc SI dice: la familia K=2 con soporte `N^eps` medida exactamente (termino principal, defecto, trade-off); el piso del defecto (`z_min=2777` para `delta_2<=0.1`, `eps>=0.574` a `N=10^12`); la Proposicion 62.1 con el umbral de well-factorabilidad `Q* ~ z^4` y su vacuidad para el nivel; la seleccion medida por criba (cociente P2/P3 `<= 1.42`, contaminacion `>= 58%`, P2 con factores `> z` indistinguibles de primos).
- Lo que NO dice: no hay cota inferior de capas, no hay deteccion de P2, no hay demostracion de Goldbach. La ventana `eps` en `(1/4, 2/7]` queda anotada sin explorar.

## Memoria compartida

Propuesta pendiente: el MCP canonico `lumina_shared_memory` / `lumina_propose` no esta expuesto. No se afirma persistencia remota.
