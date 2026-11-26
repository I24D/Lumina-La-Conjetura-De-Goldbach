# 63 - Ventana eps en (1/4, 2/7]: BFI aplica pero es vacuo; Q* no es un nivel

## Procedencia y estado

- **Agente:** LUMINA/Muse.
- **Fecha:** 2026-11-25 ~20:40 ET. Reclamo de ESTADO.md 20:26 ET (este paso; el bucle no lo duplica).
- **Metodo:** medicion exacta (`Fraction`) en z de la ventana + `fsum` en z grande; criba `X=2*10^6` (numpy) para la seleccion; lectura de BFI/Maynard via doc 58 sec. 6 ([T-lit]).
- **Estado:** **[T]** para la Proposicion 63.1, las mediciones, el cruce finito y la seleccion; **[T-lit]** para BFI (Maynard Teor. A), BV y Jurkat-Richert `f(s)>0` para `s>2` (via docs 58/60/62). No hay cota inferior de capas, no hay nivel `> 1/2`, no hay demostracion de Goldbach.
- **Nota de numeracion:** existen dos docs 60 y dos docs 62; este es el 63, siguiente numero libre.
- **Correccion durante el trabajo:** la primera version del script calculaba `delta_2` en el regimen `t=1/S_z` con un factor 2 de mas; se detecto al contrastar con el piso del doc 62 (`z_min=2777`, `delta_2=0.100`) y se corrigio antes de medir la parrilla. La formula correcta es `delta_2 = |M_2|/A = 1/(4 S_z^2 A)` con `M_2 = -1/(4 S_z^2)` exacto (doc 61).

## 1. Lo que la ventana prometia (y la correccion de categoria)

El doc 62 sec. 4 y sec. 7.2 anoto la ventana `eps` en `(1/4, 2/7]` como el lugar donde BFI superaria a BV: `Q* ~ N^(2 eps)` (Prop. 62.1), BFI aplica si `2 eps <= 4/7` y "superaria a BV (1/2)" si `2 eps > 1/2`, es decir `eps > 1/4`. Esta doc mide esa ventana y encuentra que la premisa confunde dos objetos distintos:

- **`Q*` es el umbral de well-factorabilidad** (Prop. 62.1): el nivel a partir del cual `lambda` admite una factorizacion `gamma^(1)*gamma^(2)`. No es un nivel de distribucion.
- **El teorema BFI** (Maynard Teor. A, via doc 58 sec. 6) controla la suma `sum_{q<=Q} lambda(q) E_q` hasta `Q`, **pero la suma solo tiene terminos en el soporte** `d <= N^eps`. El nivel efectivo esta acotado por el soporte, no por `Q*`.

Medido en la ventana (parrilla sec. 3): el soporte `d_max = z^2 = N^eps < N^(1/2)` siempre (`z^4 < N` exacto en toda la parrilla), y `Q* <= N^(4/7)` (`Qmax^7 <= N^4` exacto) — BFI aplica, pero la suma que controla es la misma que BV ya controla. Vacuidad en ambos sentidos.

## 2. Proposicion 63.1 [T]: BFI nunca supera a BV para K=2

**Proposicion 63.1 [T].** Para la familia K=2 con soporte `d <= N^eps`, la cobertura de BFI es subconjunto de la de BV: el error pesado `sum_d lambda_d E_d` esta cubierto por BV siempre que BFI lo estaria.

*Prueba.* La suma corre sobre el soporte `d <= N^eps`. (i) Si `eps <= 1/2`, el soporte esta dentro del nivel BV (`N^(1/2)`): BV controla la suma directamente (`|lambda_d| <= 1`). (ii) Si `eps > 2/7`, entonces `Q* > N^(4/7)` (asintotico: `Q* ~ N^(2 eps)` con `2 eps > 4/7`; verificado exacto a `eps=1/2, 3/5`: `Qmin^7 > N^4`), luego la hipotesis BFI (well-factorable de nivel `<= N^(4/7)`) falla. En el intervalo intermedio `(2/7, 1/2]`, BFI aplicaria a nivel `Q*`, pero la suma solo llega al soporte `< N^(1/2)`, ya cubierto por BV. En todos los casos, nada mas alla de BV. QED.

Chequeos exactos (script, parte 3): caso (i) `z^4 <= N` en `eps` en `{1/4, 2/7, 1/2}` (`N=10^12, 10^24`); caso (ii) `Qmin^7 > N^4` en `eps` en `{1/2, 3/5}`; ventana: BFI aplica (`Qmax^7 <= N^4`) Y BV ya cubre (`z^4 < N`) en los 4 casos (`eps=1/4, 2/7`; `N=10^12, 10^24`) -> VACUO.

**Corolario.** No existe `eps` con BFI `>` BV para K=2. La ventana `(1/4, 2/7]` esta vacia de contenido para el nivel. La frase del doc 62 ("BFI si superaria a BV") queda corregida: `Q*` no es un nivel.

## 3. La parrilla en la ventana (medido)

`N` en `{10^12, 10^24, 10^40}`, `eps` en `{1/4, 13/50, 27/100, 2/7}`; `z = floor(N^(eps/2))` exacto; `S_z`, `S_z^{(2)}` por `fsum` (verificado contra `Fraction` en `z<=100`); `Q*` en `[max(z^2, M(z)^2), Dmax^2]` (Prop. 62.1).

```text
N      eps    z      Sz    A(t=1/2) d2(t=1/2) A(1/Sz) d2(1/Sz)  log(Q*)/log(N)  soporte<N^(1/2)  Q*<=N^(4/7)  s=1/eps
10^12  1/4    31     1.566  0.4678    0.250      0.4090   0.249     0.471-0.492      SI              SI           4.00
10^12  13/50  36     1.566  0.4678    0.250      0.4090   0.249     0.471-0.492      SI              SI           3.85
10^12  27/100 41     1.617  0.4624    0.239      0.4145   0.231     0.510-0.530      SI              SI           3.70
10^12  2/7    51     1.662  0.4583    0.229      0.4188   0.216     0.541-0.551      SI              SI           3.50
10^24  1/4    1000   2.198  0.4484    0.084      0.4532   0.114     0.499-0.500      SI              SI           4.00
10^24  13/50  1318   2.238  0.4506    0.073      0.4549   0.110     0.519-0.519      SI              SI           3.85
10^24  27/100 1737   2.275  0.4529    0.062      0.4563   0.106     0.539-0.540      SI              SI           3.70
10^24  2/7    2682   2.329  0.4570    0.047      0.4583   0.101     0.571-0.571      SI              SI           3.50
10^40  1/4    100000 2.705  0.5056    0.051      0.4691   0.073     0.500-0.500      SI              SI           4.00
10^40  13/50  158489 2.744  0.5127    0.060      0.4700   0.071     0.520-0.520      SI              SI           3.85
10^40  27/100 251188 2.782  0.5199    0.068      0.4708   0.069     0.540-0.540      SI              SI           3.70
10^40  2/7    517947 2.839  0.5314    0.080      0.4719   0.066     0.571-0.571      SI              SI           3.50
```

Lecturas:

- **`Q*/N^(2 eps)` tiende a 1 desde abajo.** El cociente logaritmico `log(Q*)/log(N)` es `0.471-0.551` a `N=10^12` y `0.500-0.571` a `N=10^40`: `Q* = c(z) N^(2 eps)` con `c(z) -> 1` cuando `z` crece (los tres primos mayores se acercan a `z`). A `N` finito el factor `c<1` pesa.
- **Cruce finito `eps*(N)`:** el menor `eps` con `Qmin > N^(1/2)` (exacto, barrido `0.25..2/7`) es `0.27` a `N=10^12` y `0.255` a `N=10^24, 10^40`. En el borde `eps=1/4`, `Q* < N^(1/2)` a todo `N` medido: la afirmacion "`Q* > N^(1/2)` en la ventana" es solo asintotica, y en la mitad baja de la ventana no se cumple ni siquiera como umbral.
- **Defecto en la ventana:** `delta_2(t=1/2) = 0.23-0.25` a `N=10^12` (0.05-0.08 a `N=10^40`, pero es el cero accidental de `M_2` en `S_z=2.5`, doc 62 sec. 2: `delta_3` no se anula ahi); `delta_2(1/S_z) = 0.216-0.249` a `N=10^12`, `0.101-0.114` a `N=10^24`, `0.066-0.073` a `N=10^40`. Sigue siendo `O(0.2)` a `N` realista, como ya decia el doc 62 sec. 7.2. `A = Theta(1)` en toda la ventana (`0.41-0.53`).

## 4. El parametro s: dos marcos, ninguno da capas

Hay que distinguir dos `s` que el doc 62 sec. 4 no separaba:

- **Marco A (criba lineal estandar sobre `A_N`, nivel BV `1/2`, corte `z = N^(eps/2)`):** `s = log(N^(1/2))/log(z) = 1/eps`, en `[3.5, 4)` en la ventana, `f(s) > 0` [T-lit, Jurkat-Richert via docs 60/62]. Pero ya era `> 2` para todo `eps < 1/2` (`eps=1/5` da `s=5`, doc 62 sec. 4): cuenta numeros rugosos, no P2. Nada nuevo en la ventana.
- **Marco B (K=2 como peso de criba, `D` = soporte = `N^eps`, `z = N^(eps/2)`):** `s_K2 = log(z^2)/log(z) = 2.000000` exacto (medido en `z=31, 52, 1000`) — identicamente 2 en todo `eps`, barrera de paridad, `f(2) = 0`. La ventana no mueve `s_K2`.

Las capas (doc 38) necesitan nivel `> 1/2` con corte `> N^(1/4)` (`s < 2` a nivel `1/2`, `f(s)=0`). La ventana no toca ni el nivel de `A_N` (`1/2`, BV) ni el corte de las capas. **El nivel efectivo es `1/2`, sin cambio.**

## 5. Seleccion en la ventana (criba X=2*10^6)

`z=7` es el control del doc 62; `31, 52` son la ventana a `N=10^12`; `1000, 2680` son la ventana a `N=10^24`. `W(n) = 1 - t r + t^2 C(r,2)`, `r = omega_z(n)`.

```text
z     t    frac r=0  E[W]  E[W|Om=1] E[W|Om=2] E[W|Om=3] E[W|Om>=4]  P2 frac r=0  ratio P2/P3  contam P3+/P2
7     0.5  0.229     0.532  1.000     0.765     0.538     0.355        0.530        1.421        0.547
7     0.1  0.229     0.887  1.000     0.953     0.901     0.834        0.530        1.058        0.899
31    0.5  0.153     0.468  1.000     0.675     0.414     0.323        0.350        1.631        0.525
31    0.1  0.153     0.853  1.000     0.935     0.867     0.788        0.350        1.079        0.872
52    0.5  0.139     0.458  1.000     0.651     0.386     0.328        0.301        1.684        0.535
52    0.1  0.139     0.845  1.000     0.930     0.858     0.778        0.301        1.084        0.866
1000  0.5  0.076     0.448  0.999     0.496     0.289     0.423        0.009        1.716        0.760
1000  0.1  0.076     0.802  1.000     0.898     0.809     0.726        0.009        1.110        0.841
2680  0.5  0.074     0.452  0.999     0.467     0.288     0.446        0.000        1.624        0.838
2680  0.1  0.074     0.792  1.000     0.888     0.796     0.715        0.000        1.115        0.837
```

Lecturas:

- **El cociente P2/P3 mejora un poco con `z` (1.42 -> 1.63-1.72 a `t=1/2`) pero se queda en `O(1)`.** No aisla P2.
- **La ceguera `r=0` se encoge con `z`** (P2 con ambos factores `> z`: `0.53 -> 0.35 -> 0.30 -> 0.009 -> 0.000`), **pero la contaminacion P3+ CRECE** (`0.55 -> 0.84`): `W(r) = 1 - t r + t^2 C(r,2)` crece cuadraticamente para `r >= 3` a `t=1/2` (`W(5)=1.0`, `W(6)=1.75`), asi que los numeros con `Omega>=4` reciben peso `> 1` cuando `z` es grande. La ventana cambia una ceguera por otra: no mejora la seleccion, la redistribuye.
- `E[W|Om=1] = 1.000` (`0.999` en `z>=1000`: los primos `<= z` reciben `1-t`).

## 6. Veredicto

**NO: la ventana `eps` en `(1/4, 2/7]` NO da un nivel `> 1/2` usable para capas.** Tres razones independientes:

1. **Prop. 63.1:** BFI nunca supera a BV para K=2. El soporte capa la suma; `Q*` no es un nivel. En la ventana BFI aplica pero es vacuo; fuera de ella (`eps > 2/7`) ni siquiera aplica.
2. **Parametro `s`:** `s_K2 = 2` identico (barrera de paridad, `f(2)=0`) y `s_std = 1/eps > 2` ya lo era para todo `eps < 1/2` sin seleccionar P2.
3. **Seleccion:** no mejora en la ventana (cociente P2/P3 `<= 1.72`, contaminacion `>= 52%`, defecto `O(0.2)` a `N=10^12`).

Lo que SI queda: la correccion de categoria (`Q*` no es un nivel) con prueba; el cruce finito `eps*(N)` medido (`0.27` a `N=10^12`, `0.255` a `N>=10^24`); la parrilla completa en la ventana; y la seleccion medida, que muestra que subir `z` redistribuye la contaminacion en vez de eliminarla. La linea de los docs 54-62 queda cerrada tambien por este camino: ni el peso, ni el soporte, ni la ventana BFI suben el nivel. Goldbach sigue abierta.

## 7. Que sigue

1. **Mecanismo de seleccion** (switching, C10/peso movil `W(d)`, frentes del bucle): lo que falta no es nivel, sino un peso que seleccione `Omega <= 2` sin pagar el defecto maximo (Prop. 61.3).
2. **Frentes del bucle:** C10/W(d), uniones `C_k+C_{k+1}`, variante Lema 4 Pintz-Althofer, Motohashi (doc 57 sec. 4), `b_{P,l}` (doc 53), atomizacion de Kuhn (doc 62 del bucle).
3. **Dispersion uniforme en `a`** (doc 59 #1) sigue abierta como problema independiente.

## Verificacion reproducible

```text
python numerico/63_ventana_epsilon_bfi.py
```

PASS local y PASS en el PC (numeros identicos; `salida_doc63.txt`). Chequeos: K=2 exacto (`Fraction`) en `z=31, 52` (A, M_l, identidad `L = sum M_e/phi(e)`, `W>=0`, `W(p)>0`, soporte `<= z^2`); parrilla (N,eps) con `z` exacto por aritmetica entera; Prop. 63.1 (enteros exactos: `z^4<=N`, `Qmin^7>N^4`, `Qmax^7<=N^4`); `s_K2 = 2` exacto; seleccion por criba a `X=2*10^6`.

## Alcance honesto

- Lo que este doc SI dice: la correccion de categoria sobre `Q*` (Prop. 63.1) con prueba y chequeos exactos; la parrilla completa en la ventana (termino principal, defecto, `Q*`, soporte, `s`); el cruce finito `eps*(N)`; los dos marcos de `s` distinguidos; la seleccion medida en la ventana.
- Lo que NO dice: no hay cota inferior de capas, no hay nivel `> 1/2`, no hay deteccion de P2, no hay demostracion de Goldbach. La ventana `(1/4, 2/7]` queda cerrada como ruta al nivel.

## Memoria compartida

Propuesta pendiente: el MCP canonico `lumina_shared_memory` / `lumina_propose` no esta expuesto. No se afirma persistencia remota.
