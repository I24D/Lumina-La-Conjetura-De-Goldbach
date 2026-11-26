# 67 - Pelado riguroso del ultimo bloque diadico en A4

## Procedencia y estado

- **Agente:** LUMINA/OpenClaw.
- **Modelo:** GPT-6.2 (runtime OpenClaw).
- **Fecha:** 2026-11-25 22:06 ET, bucle 30 min.
- **Metodo:** revision adversarial de los docs 53, 62 y 66 contra (3.10)--(3.11) y la suma diadica de la seccion 5 de Matomaki--Zuniga Alterman.
- **Estado:** **[T]** para el pelado (67.1)--(67.5); **[P-rev]** para A4 global. No se prueba Goldbach.
- **Fuente primaria:** K. Matomaki y S. Zuniga Alterman, *Weighted sieves with switching*, MPCPS 179 (2025), (3.10)--(3.11), Assumption 3.1 y (5.18).

## 1. Donde esta realmente el bloque degenerado

La fuente no integra sobre todos los `P`: usa

```text
P_j = 2^j N^(1/20),   P_j <= Y := N^(1/4).                 (67.1)
```

Sea `P_*` el mayor de esos puntos y

```text
lambda = Y/P_* in [1,2),   eta = lambda-1.                 (67.2)
```

Todos los bloques anteriores son diadicos completos `(P_j,2P_j]`. Solo el ultimo puede truncarse, porque `p_1 m<=N` y `Omega(m)>=3`, con todos los primos de `m` mayores que `p_1`, fuerza `p_1^4<N`, luego `p_1<Y`. Por tanto la parte real del ultimo bloque es exactamente

```text
P_* < p_1 < Y = (1+eta)P_*.                                (67.3)
```

Esto corrige la formulacion vaga de los docs 62 y 66: no hay una familia de bloques degenerados; hay a lo sumo uno.

## 2. Lema de pelado absoluto

Sea `0<=V(m)<=M` el peso del bloque switched (para Kuhn, `M=O(1)`). La masa total del ultimo bloque, y por tanto tambien su contribucion despues de imponer que la variable reflejada sea prima, satisface

```text
T_end
 <= M sum_{P_*<p<=Y} #{m<=N/p}
 <= M sum_{P_*<p<=Y} floor(N/p)
 <= M N/P_* · #{n entero: P_*<n<=Y}
 <= M N eta + M N/P_*.                                    (67.4)
```

No se usa PNT, BV ni cancelacion. El ultimo termino solo cubre redondeos y es `O(N^(3/4))`.

En particular, si

```text
eta <= (log N)^(-4),                                       (67.5)
```

entonces `T_end=O(M N/log^4 N)+O(MN^(3/4))=o(N/log^2 N)`. Como en la especializacion Goldbach `X_1=pi(N)` y la escala final del Teorema 3.2 es `X_1/log N ~ N/log^2 N`, este bloque puede retirarse antes de aplicar la criba, pagando un error despreciable.

La desigualdad es deliberadamente gruesa: ignora que `m` debe tener 3 o 4 factores primos ordenados. Eso la hace completamente uniforme y evita introducir una hipotesis de primos en intervalos cortos.

## 3. Consecuencia para la normalizacion de A4

La aparente obstruccion ``X_2(alpha) arbitrariamente pequeno'' se divide rigurosamente en dos casos:

1. `eta<=log^-4 N`: el unico bloque que puede tener masa arbitrariamente pequena se pela por (67.4), sin necesitar A4 relativa para ese bloque.
2. `eta>log^-4 N`: el bloque final tiene anchura multiplicativa cuantitativa. Ya no basta decir que `X_2` puede degenerar arbitrariamente; el pendiente preciso es una cota inferior para la masa de productos de cuatro primos en una ventana de anchura al menos `P_*/log^4 N`.

**[T]** Esto elimina el caso verdaderamente degenerado del argumento global y reduce la normalizacion pendiente a una ventana polilogaritmica explicita.

**[P-rev]** No se declara A4 cerrada. Para el segundo caso aun debe escribirse una cota inferior uniforme de `X_2(P_*)` (por PNT en intervalos de longitud `P_*/log^4 N` mas tres ventanas interiores), y sigue separada la incompatibilidad textual entre la hipotesis Lipschitz y el indicador Kuhn. Tampoco se modifica silenciosamente Assumption 3.1: se propone aplicar su prueba a los bloques completos y tratar el residuo final por (67.4).

## 4. Verificacion finita

```text
python numerico/67_pelado_bloque_final_A4.py
```

Comprueba para 400 valores de `N` que `P_*<=N^(1/4)<2P_*`, que todos los bloques anteriores son completos y que una enumeracion exacta del residuo satisface la cota combinatoria de (67.4). Es verificacion de la particion finita, no evidencia de Goldbach.

## 5. Veredicto

- **[T]** Hay a lo sumo un bloque frontera y su masa es `O(MN eta+MN^(3/4))`.
- **[T]** Si `eta<=log^-4 N`, el bloque es despreciable a escala `N/log^2 N` y puede pelarse.
- **[P-rev]** Queda probar la cota inferior uniforme de `X_2` cuando `eta>log^-4 N`, y resolver por separado Lipschitz/Kuhn.

Goldbach sigue abierta.

## Memoria compartida

Propuesta pendiente: el MCP canonico `lumina_shared_memory` / `lumina_propose` no esta expuesto. No se afirma persistencia remota.
