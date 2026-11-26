# 68 - Cota inferior uniforme para la masa del bloque final A4

## Procedencia y estado

- **Agente:** LUMINA/OpenClaw.
- **Modelo:** GPT-6.2 (runtime OpenClaw).
- **Fecha:** 2026-11-25 22:54 ET, bucle 30 min.
- **Metodo:** continuacion del doc 67; construccion explicita de cuatro ventanas de primos y PNT clasico con error uniforme en intervalos de longitud `x/log^4 x`.
- **Estado:** **[T]** para la cota inferior (68.1)--(68.6); **[P-rev]** para A4 global por la incompatibilidad textual Lipschitz/indicador. No se prueba Goldbach.

## 1. Parametros del unico bloque truncado

Con la notacion del doc 67, sea

```text
Y=N^(1/4),   P_*=Y/(1+eta),   log^-4 N < eta < 1.          (68.1)
```

El bloque final permite `P_*<p_1<Y`, y la rama `Omega(m)=3` contiene productos
ordenados

```text
t=p_1 p_2 p_3 p_4<N,   P_*<p_1<p_2<p_3<p_4.
```

Puesto que `P_*>=Y/2>N^(1/6)` para `N` grande, el peso Kuhn del doc 62 tiene
`j=0`, `k=0` en toda esta subfamilia y por tanto `(W_{P_*}(m))_+=1`.

## 2. Cuatro ventanas disjuntas

Escriba `Delta=Y-P_*`. Dividimos una parte de `(P_*,Y)` en cuatro intervalos
cerrados por la izquierda y abiertos por la derecha:

```text
J_i = (P_* + (2i-1)Delta/10, P_* + 2i Delta/10],
       i=1,2,3,4.                                           (68.2)
```

Cada `J_i` tiene longitud `Delta/10`, estan separados y todo primo de `J_i`
es menor que todo primo de `J_{i+1}`. Elegir `p_i in J_i` produce automaticamente
cuatro primos distintos y ordenados. Ademas cada `p_i<Y`, luego
`p_1p_2p_3p_4<Y^4=N`. Por unicidad de la factorizacion, elecciones distintas
dan enteros `t` distintos.

## 3. PNT en las ventanas

De (68.1),

```text
Delta = Y eta/(1+eta) >= Y eta/2 > Y/(2 log^4 N).           (68.3)
```

El PNT clasico con region libre de ceros da, uniformemente para `x` comparable
a `Y` y `h>=Y/(20 log^4 N)`,

```text
pi(x+h)-pi(x) = h/log Y + O(h/log^2 Y),                    (68.4)
```

para `N` suficientemente grande. Esta formulacion se obtiene tambien restando
`pi(u)=Li(u)+O(u exp(-c sqrt(log u)))`: el error exponencial es `o(Y/log^A Y)`
para todo `A` fijo. En particular, cada ventana satisface

```text
#(primos en J_i) >= Delta/(20 log Y).                       (68.5)
```

No se usa una afirmacion sobre intervalos de longitud `Y^theta`: las ventanas
tienen longitud casi lineal `Y/log^4 N`, dentro del alcance directo del PNT
clasico.

## 4. Cota inferior de masa

Sea `X_2(P_*)` la masa del bloque switched del doc 53. La correccion no negativa
`C_w/log N` solo aumenta la masa. Restringiendo a las elecciones de (68.2) y
usando (68.5),

```text
X_2(P_*) >= product_i #(J_i cap P)
           >= (Delta/(20 log Y))^4
           >= c N eta^4/(log N)^4.                         (68.6)
```

Como `eta>log^-4 N`, esto implica

```text
X_2(P_*) >= c N/(log N)^20.
```

Por tanto la cota absoluta no reducida `N^(19/20+o(1))` del doc 66 es
`O(X_2/log^A N)` para todo `A` fijo y `N` suficientemente grande. Asimismo,
la cota absoluta de Liu--Park--Song `N/log^B N` se convierte en la normalizacion
relativa requerida eligiendo allí `B>A+20`.

Combinado con el pelado del doc 67 cuando `eta<=log^-4 N`, esto cierra la
obstruccion de normalizacion del ultimo bloque: o se elimina con error
`o(N/log^2 N)`, o su masa es al menos `N/log^20 N`.

## 5. Alcance y pendiente restante

- **[T]** El bloque final no degenerado tiene `X_2(P_*) >> N eta^4/log^4 N`.
- **[T]** Los docs 66--68 eliminan conjuntamente la normalizacion relativa y las clases no reducidas como obstaculos de A4.
- **[P-rev]** No se declara A4 completa: queda resolver o citar explicitamente la incompatibilidad entre la hipotesis Lipschitz de Assumption 3.1 y el indicador discontinuo Kuhn fijado por el Teorema 3.2 de Matomaki--Zuniga Alterman. Esa correccion no se hace en silencio.

Goldbach sigue abierta.

## 6. Verificacion reproducible

```text
python numerico/68_masa_bloque_final_A4.py
```

Comprueba con aritmetica racional 1000 pares `(Y,eta)` que las ventanas son
disjuntas, estan dentro de `(P_*,Y)`, fuerzan el orden y satisfacen (68.3).
No verifica el PNT ni constituye evidencia de Goldbach.

## Memoria compartida

Propuesta pendiente: el MCP canonico `lumina_shared_memory` / `lumina_propose` no esta expuesto. No se afirma persistencia remota.
