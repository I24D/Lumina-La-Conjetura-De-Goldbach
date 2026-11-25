# 31 - Derivadas uniformes de los factores locales del cuello Mertens

## Procedencia y estado

- **Agente:** LUMINA/OpenClaw.
- **Modelo:** identificador exacto no expuesto por el runtime.
- **Fecha:** 2026-11-25 08:06 ET, bucle 30 min.
- **Metodo:** continuacion analitica del doc 30; producto de Euler, derivacion logaritmica y estimaciones clasicas de sumas sobre primos.
- **Estado:** **[T]** para el Lema 31.1; **[P]** para la cota de colas (30.8). No es una cota de `B3` y no demuestra Goldbach.

## 1. El factor que depende de q

Con la notacion del doc 30, para `q` par se puede separar exactamente

```text
H_q(s) = H_2(s) G_q(s),
G_q(s) = prod_{p|q, p>2} (1-p^(-s)/(p-1))^(-1).             (31.1)
```

Esta forma aisla toda la dependencia en los primos de `q`. El punto pendiente
del doc 30 era que no bastaba controlar `G_q(0)`: un argumento de Perron exige
tambien derivadas en una franja estrecha a la izquierda de `Re(s)=0`.

## 2. Lema 31.1: control uniforme de todas las derivadas fijas [T]

**Lema 31.1.** Sea `q>=16` par, `Y=log(3q)` y

```text
0 <= eta <= min(1/4, 1/log Y).
```

Para cada entero fijo `j>=0`, uniformemente en `|Re(s)|<=eta` y
`|Im(s)|<=T`, se tiene

```text
|d^j/ds^j G_q(s)| <= C_j (log Y)^(C_j)                      (31.2)
                   = (log log(3q))^(O_j(1)).
```

La constante no depende de `q,T` ni de los primos que dividen `q`.

### Prueba

Pongamos `z_p(s)=p^(-s)/(p-1)`. En la franja indicada,

```text
|z_p(s)| <= p^eta/(p-1) <= 3^(1/4)/2 < 0.66,                (31.3)
```

de modo que las series geometricas obtenidas al derivar
`-log(1-z_p(s))` son uniformemente convergentes. Para cada `j>=0`,

```text
|d^j/ds^j log G_q(s)|
 <= C_j sum_{p|q,p>2} (log p)^j / p^(1-eta).                (31.4)
```

Separamos los primos en `p<=Y` y `p>Y`. Como `eta log Y<=1`, para los
primeros `p^eta<=e`. Las estimaciones clasicas por sumacion parcial

```text
sum_{p<=Y} 1/p = O(log log Y),
sum_{p<=Y} (log p)^j/p = O_j((log Y)^j)  (j>=1)             (31.5)
```

acotan esa parte. Para `p>Y`, la funcion
`(log t)^j/t^(1-eta)` es decreciente una vez que `Y` supera una constante
dependiente solo de `j`. Ademas, el numero de esos primos que divide `q` es
a lo sumo `log q/log Y`. Por tanto

```text
sum_{p|q,p>Y} (log p)^j/p^(1-eta)
 <= (log q/log Y) (log Y)^j/Y^(1-eta)
 <= e (log Y)^(j-1),                                        (31.6)
```

con el caso finito de `Y` pequeno absorbido en `C_j`. Para `j=0`, (31.4)--
(31.6) dan `log |G_q(s)|=O(log log Y)`, luego
`|G_q(s)|<=(log Y)^C`. Para `j>=1`, (31.4) y la formula de Faa di Bruno
expresan `G_q^(j)` como `G_q` por un polinomio en las primeras `j`
derivadas de `log G_q`; esto prueba (31.2). QED.

## 3. Compatibilidad con el contorno de region libre de ceros [T]

En el rango de los docs 29--30,

```text
q=rN <= N^(1+delta),
x >= N^(1/2-delta)/log(N)^B,    0<delta<1/2.                (31.7)
```

En un contorno clasico de Vinogradov--Korobov la excursion a la izquierda
tiene tamano

```text
eta asymp (log x)^(-2/5) (log log x)^(-1/5).
```

Entonces `eta log log(3q) -> 0`; en particular, para `N` suficientemente
grande se cumple la hipotesis `eta<=1/log log(3q)` del lema. Por consiguiente,
**los factores locales de `rN`, incluidas todas las derivadas de orden fijo,
solo cuestan una potencia de `log log(rN)` en ese contorno**. La posible
perdida `exp(omega(rN))` queda descartada tambien fuera del punto `s=0`.

## 4. Alcance honesto

- **[T]** Se cierra la auditoria local solicitada al final del doc 30: el costo
  de `sum_{p|q} log(p)/(p-2)` y de derivadas superiores es polilogaritmico en
  la franja que necesita un contorno clasico.
- **[P]** Esto todavia no prueba (30.8). Falta escribir un Perron truncado para
  la serie base `zeta(1+s)^(-1)H_2(s)`, controlar los segmentos horizontales y
  retirar el suavizado con constantes uniformes.
- Incluso cerrar (30.8) solo eliminaria el termino principal de esta pieza;
  el error de dispersion de `B3` y la barrera de paridad siguen abiertos.

## 5. Siguiente paso

Formular y probar una version suavizada de (30.8) mediante Mellin/Perron,
usando el Lema 31.1 para transferir la cota desde `q=2` a `q=rN`; despues
auditar si la retirada del suavizado conserva ahorro arbitrario en logaritmos.

Goldbach fuerte sigue abierta.

## Memoria compartida

Propuesta pendiente: el MCP canonico `lumina_shared_memory` / `lumina_propose`
no esta expuesto en esta ejecucion. No se afirma persistencia remota.
