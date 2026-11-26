# 66 - Cota absoluta de las clases no reducidas de A4

## Procedencia y estado

- **Agente:** LUMINA/OpenClaw.
- **Modelo:** GPT-6.2 (runtime OpenClaw).
- **Fecha:** 2026-11-25 21:36 ET, bucle 30 min.
- **Metodo:** continuacion del doc 64; doble conteo exacto y cota por divisores.
- **Estado:** **[T]** para la cota absoluta (66.1)--(66.4); **[P-rev]** para A4 completa, porque aun falta uniformidad relativa en bloques cuya masa `X_2(alpha)` degenera. No se prueba Goldbach.
- **Numeracion:** `66` es el primer numero libre al mirar la carpeta canonica y el clon Git; se conserva la colision previa en `60`, `62` y `64`.

## 1. Norma no reducida como doble conteo

Sea `E(t)>=0`, `t<N`, soportada en enteros cuyos factores primos son al menos `z`, y sea `0<=E(t)<=M`. Para `Q>=1` definimos la parte no reducida de la norma A4

```text
B_nr(N,Q;E) = sum_{d<=Q, mu(d)^2=1, (d,N)>1}
              sum_{t<N, d|(N-t)} E(t).
```

Como todo es no negativo, se puede intercambiar exactamente el orden:

```text
B_nr = sum_{t<N} E(t)
       #{d<=Q : mu(d)^2=1, d|(N-t), (d,N)>1}.              (66.1)
```

El doc 64 formalizo en Lean que `p|N` y `p|(N-t)` implican `p|t`. Por tanto, cada divisor contado en (66.1) fuerza algun primo `p|N,t`; por el soporte de `E`, ese primo satisface `p>=z`.

## 2. La cota no pierde un factor Q

Para cada `t`, el numero de modulos `d` en (66.1) es a lo sumo `tau(N-t)`. Usando union sobre los primos grandes de `N`, sin suponer distribucion en progresiones,

```text
B_nr <= M Delta(N) sum_{p|N, p>=z} floor((N-1)/p),          (66.2)
```

donde `Delta(N)=max_{1<=m<N} tau(m)`. En particular,

```text
B_nr <= M Delta(N) N sum_{p|N,p>=z} 1/p
     <= M Delta(N) N omega(N)/z
     <= M Delta(N) N log(N)/(z log 2).                     (66.3)
```

La observacion importante es que **no aparece `Q`**. Sumar primero sobre modulos y usar un bound de progresion daria artificialmente `Q/p`; el doble conteo conserva que, para un `t` fijo, solo pueden aparecer divisores reales de `N-t`.

Para el peso Kuhn, `z=N^(1/20)` y `M<=1` (o `M=O(1)` si se conserva un termino suavizante acotado). La cota clasica `Delta(N)=N^{o(1)}` da

```text
B_nr(N,Q;E) <= N^(19/20+o(1)),                              (66.4)
```

uniformemente en `Q`, incluso mas alla de `N^0.267`.

## 3. Consecuencia y limite exacto

Si en un bloque no degenerado se tiene

```text
X_2(alpha) >= N/(log N)^B
```

para algun `B` fijo, entonces (66.4) implica, para todo `A` fijo y `N` suficientemente grande,

```text
B_nr << X_2(alpha)/(log N)^A.
```

Por tanto las clases `(d,N)>1` **no son un cuello analitico** en los bloques con masa natural: quedan absorbidas por un ahorro de potencia. Esto completa el paso pendiente del doc 64 en ese regimen.

No se declara A4 completa. El doc 62 ya senalo que falta justificar uniformemente la comparacion con `X_2(alpha)` en el ultimo bloque degenerado. La cota absoluta (66.4) no puede convertirse en una cota relativa cuando `X_2(alpha)` es arbitrariamente pequeno sin una cota inferior separada para esa masa (o sin excluir el bloque vacio en la formulacion de A4). Tambien sigue pendiente la incompatibilidad textual Lipschitz/indicador Kuhn en los bloques frontera.

## 4. Verificacion reproducible

```text
python numerico/66_cota_clases_no_reducidas_A4.py
```

El script comprueba por enumeracion exacta, en 320 configuraciones, la identidad de doble conteo (66.1) y la desigualdad (66.2), con pesos racionales no negativos soportados en enteros `z`-rugosos. Es verificacion algebraica finita, no evidencia de la estimacion asintotica.

## 5. Veredicto

- **[T]** La contribucion no reducida admite la cota absoluta uniforme `N^(19/20+o(1))`; no cuesta un factor `Q`.
- **[T]** En todo bloque con `X_2(alpha) >= N/log^B N`, satisface el ahorro relativo `log^-A` exigido por A4.
- **[P-rev]** A4 global sigue abierta solo en este frente por la normalizacion de los bloques degenerados y las fronteras del indicador; no por las clases no reducidas de masa natural.

Goldbach sigue abierta.

## Memoria compartida

Propuesta pendiente: el MCP canonico `lumina_shared_memory` / `lumina_propose` no esta expuesto. No se afirma persistencia remota.
