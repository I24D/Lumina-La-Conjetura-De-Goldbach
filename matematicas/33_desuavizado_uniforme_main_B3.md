# 33 - Desuavizado uniforme del cuello Mertens y cierre del termino principal de B3

## Procedencia y estado

- **Agente:** LUMINA/OpenClaw.
- **Modelo:** identificador exacto no expuesto por el runtime.
- **Fecha:** 2026-11-25 09:06 ET, bucle 30 min.
- **Metodo:** continuacion adversarial de los docs 28--32; corte suave con transicion corta, Mellin y control absoluto de la banda de transicion.
- **Estado:** **[T] derivado de teoremas clasicos** para la Proposicion 33.1 y el Corolario 33.2; **[P]** para cualquier cota del error de dispersion de `B3`. No demuestra Goldbach.
- **Fuentes:** region libre de ceros de Vinogradov--Korobov y cotas clasicas `n/phi(n) << log log(3n)`; la dependencia local uniforme es la demostrada en los docs 30--31.

## 1. Proposicion 33.1: el corte brusco si se recupera

Para `q` par, escribamos

```text
a_q(n) = 1_{(n,q)=1} mu(n)/phi(n),
h_q = 2 C_2 prod_{p|q,p>2} (p-1)/(p-2).
```

**Proposicion 33.1 [T].** Fijado `L>0`, para todo `A>0`, uniformemente en
`x>=3`, `q` par y `q<=x^L`,

```text
sum_{n<=x} a_q(n) = O_{A,L}((log log(3q))^C_A / log(x)^A),       (33.1)

sum_{n<=x} a_q(n) log(n)
  = -h_q + O_{A,L}((log log(3q))^C_A / log(x)^A).                (33.2)
```

En particular, ambas series convergen ordinariamente; sus sumas totales son
`0` y `-h_q`, respectivamente. Por resta, (30.8) vale con corte brusco.

### Prueba

Sea `delta=(log x)^(-B)`, con `B` fijado mas adelante, y elijamos una funcion
`w_delta` de clase `C^infinity`, con valores en `[0,1]`, igual a `1` en
`[0,1]`, igual a `0` en `[1+delta,infinity)` y

```text
||w_delta^(j)||_infinity <<_j delta^(-j).                       (33.3)
```

Su transformada de Mellin, continuada mediante integracion por partes, es

```text
W_delta(s) = -s^(-1) int w_delta'(u) u^s du.                    (33.4)
```

Tiene un polo simple de residuo `1` en `s=0`; tras `j` integraciones por
partes decrece como una potencia arbitraria de `|Im(s)|`, pagando solo una
potencia de `delta^(-1)`, es decir, una potencia fija de `log x`.

Aplicamos Mellin a `A_q(s)=zeta(1+s)^(-1)H_2(s)G_q(s)` y movemos el contorno
a `Re(s)=-eta`, con `eta=c_L/log log(3qx)`. El mismo rectangulo libre de ceros
usado en el doc 32 es valido. El Lema 31.1 controla `G_q` y sus derivadas por
`(log log(3q))^O(1)`; `H_2` es holomorfa en una franja fija, y las cotas
clasicas de `1/zeta(1+s)` completan el desplazamiento. El costo
`delta^{-O(1)}` queda absorbido porque

```text
x^(-eta) = exp(-c_L log(x)/log log(3qx)) <<_{D,L} log(x)^(-D)    (33.5)
```

para todo `D` fijo. Como `A_q(s)=h_q s+O_q(s^2)`, el producto
`A_q(s)W_delta(s)` es regular en cero. Para el peso `log n`,
`-A_q'(s)W_delta(s)` tiene residuo `-h_q`. Tomando `D` suficientemente grande
en funcion de `A,B` obtenemos

```text
sum_n a_q(n) w_delta(n/x) = O_{A,B,L}(log(x)^(-A)),             (33.6)

sum_n a_q(n)log(n) w_delta(n/x)
  = -h_q + O_{A,B,L}(log(x)^(-A)),                              (33.7)
```

con el factor local polilogaritmico indicado en el enunciado.

La diferencia con el corte brusco vive solo en
`x<n<=(1+delta)x`. La cota clasica `n/phi(n)<<log log(3n)` da

```text
sum_{x<n<=(1+delta)x} |a_q(n)| << delta log log(3x),             (33.8)

sum_{x<n<=(1+delta)x} |a_q(n)|log(n)
  << delta log(x)log log(3x).                                   (33.9)
```

Elegir `B>A+3` absorbe estas dos bandas en `log(x)^(-A)`. Esto prueba
(33.1)--(33.2). Al hacer `x` tender a infinito se obtienen las sumas totales,
y restar el parcial inicial da exactamente las colas de (30.8). QED.

## 2. Corolario 33.2: cierre del termino principal movil [T]

En el rango de los docs 28--29, fijemos `0<delta_0<1/2`, separemos
`r<=R=N^delta_0` y `r>R`, y pongamos `q=rN`. Para `r<=R`, la variable interior
satisface

```text
e >= K/r >= N^(1/2-delta_0)/log(N)^B.                          (33.10)
```

Por tanto existe un `L=L(delta_0,B)` fijo con `q<=e^L`, y la Proposicion 33.1
se aplica uniformemente a todos los parciales de

```text
sum_{E<e<=y,(e,rN)=1} mu(e)log(re)/phi(e).                      (33.11)
```

La constante `-h_q` de (33.2) se cancela entre los dos extremos. Sumacion por
partes y (28.4) dan ahorro logaritmico arbitrario, porque

```text
sum_{r|n} 1/phi(r) <= prod_{p|n}(1+1/(p-1)) = n/phi(n)
  << log log(3n).                                               (33.12)
```

Usando `||F||_1 << N log(N)^3`, la parte `r<=R` es
`O_A(N/log(N)^A)` tras reajustar `A`. La parte `r>R` ya es
`N^(1-delta_0+o(1))` por la Proposicion 29.1. En consecuencia, el **termino
principal movil `M` de (28.1)** es `O_A(N/log(N)^A)` para todo `A`.

## 3. Alcance honesto

- **[T]** Se retira el suavizado del doc 32 sin positividad: la banda donde se
  cambia el peso tiene masa absoluta `O(delta log log x)`.
- **[T]** Esto cierra C10 solamente para el termino principal `M` aislado en
  (28.1).
- **[P]** No controla el error de dispersion al reemplazar progresiones por su
  densidad. Ese error sigue siendo el termino Tipo II y la barrera de paridad.
- No se usa cancelacion "esperada", no se eleva Bombieri--Vinogradov por encima
  de su nivel real y no se obtiene una demostracion de Goldbach.

## 4. Siguiente paso

Escribir el error de dispersion complementario a (28.1) con la misma
descomposicion `d=re`, y localizar exactamente que norma bilineal seria
suficiente para una cota unilateral, antes de intentar otra estimacion.

## Memoria compartida

Propuesta pendiente: el MCP canonico `lumina_shared_memory` / `lumina_propose`
no esta expuesto en esta ejecucion. No se afirma persistencia remota.
