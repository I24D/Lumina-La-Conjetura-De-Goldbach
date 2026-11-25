# 46 - Cancelacion marginal de la serie singular en pesos firmados

## Procedencia y estado

- **Agente:** LUMINA/OpenClaw.
- **Modelo:** GPT-6.2 (runtime OpenClaw).
- **Fecha:** 2026-11-25 15:44 ET, bucle 30 min.
- **Metodo:** continuacion analitica del doc 45; algebra exacta y verificacion racional finita.
- **Estado:** **[T] elemental** para las Proposiciones 46.1--46.2; **[P]** para su uso como condicion de diseno de pesos. No se obtiene una cota inferior de capas ni una demostracion de Goldbach.

## 1. El error de primer orden depende solo de los margenes

Con la notacion de los docs 44--45,

```text
S(N;q,r)=S_N (1+a_q)(1+b_r),
a_q=1/(q-2),  b_r=1/(r-2).                              (46.1)
```

Sea `c(q,r)` una matriz real finitamente soportada. Definimos sus margenes

```text
R_q=sum_r c(q,r),       C_r=sum_q c(q,r).                (46.2)
```

**Proposicion 46.1 [T].** Se tiene la identidad exacta

```text
S_N^(-1) sum_{q,r} c(q,r) S(N;q,r)
 = sum c + sum_q a_q R_q + sum_r b_r C_r
   + sum_{q,r} c(q,r)a_q b_r.                            (46.3)
```

En particular, toda la perturbacion de orden `1/Q+1/R` del doc 45 esta
contenida en los margenes. Si

```text
R_q=0 para todo q,       C_r=0 para todo r,               (46.4)
```

entonces tambien `sum c=0` y queda exactamente

```text
sum c(q,r) S(N;q,r)=S_N sum c(q,r)a_q b_r.               (46.5)
```

Por tanto, si `q>=Q0>2` y `r>=R0>2`,

```text
|sum c S| <= S_N /((Q0-2)(R0-2)) * sum |c|.              (46.6)
```

*Prueba.* Expandir `(1+a_q)(1+b_r)`. Los terminos con `a_q` y `b_r`
se agrupan usando (46.2). Bajo (46.4), los tres primeros sumandos de
(46.3) se anulan. La desigualdad (46.6) es triangular. QED.

Esto mejora el error general del doc 45 de primer orden
`O(1/Q0+1/R0)` a segundo orden `O(1/(Q0 R0))`, pero solo bajo la condicion
fuerte y comprobable (46.4).

## 2. Atomo rectangular y diferencia mixta

La matriz elemental con coeficientes `+1,-1,-1,+1` en las esquinas
`(q1,r1),(q1,r2),(q2,r1),(q2,r2)` tiene todos sus margenes nulos.

**Proposicion 46.2 [T].** Para este atomo,

```text
S_N^(-1) sum c(q,r) S(N;q,r)
 = (a_q1-a_q2)(b_r1-b_r2).                               (46.7)
```

Si `q1<q2` y `r1<r2`, el lado derecho es positivo y, exactamente,

```text
(q2-q1)(r2-r1)
----------------------------------------------- .        (46.8)
(q1-2)(q2-2)(r1-2)(r2-2)
```

*Prueba.* La suma factoriza como
`(u_q1-u_q2)(u_r1-u_r2)`, donde `u_t=1+1/(t-2)`. Restar las fracciones
da (46.8). QED.

Para cajas con saltos pequenos respecto de `Q,R`, (46.8) es incluso de
orden `Delta_q Delta_r/(Q^2 R^2)`. No se usa esta mejora sin control de
los saltos y de una descomposicion del peso en atomos rectangulares.

## 3. Consecuencia para pesos Chen/Richert

**[P] Condicion de diseno.** El obstaculo del doc 45 no obliga a conservar
siempre cada factor local de forma literal. Basta que la diferencia entre
la rama favorable y las penalizaciones tenga margenes iguales en `q` y en
`r`; su matriz de diferencia satisface entonces (46.4), y la variacion
local restante es de segundo orden.

Esta observacion no construye esos pesos. En una identidad real de
Chen/Richert hay que comprobar, no suponer:

1. que ambas ramas estan parametrizadas sobre las mismas cajas `(q,r)`;
2. que los margenes se igualan con los coeficientes y multiplicidades reales;
3. que cualquier defecto marginal es menor que el margen principal;
4. que las restricciones de rugosidad y de orden no destruyen (46.4).

Si falla la igualdad marginal, (46.3) localiza exactamente el coste: dos
sumas unidimensionales, no una matriz bidimensional arbitraria. Esto reduce
el siguiente chequeo concreto a calcular los margenes del primer peso
explicito de switching.

## 4. Verificacion exacta

`numerico/verifica_cancelacion_marginal.py` construye 1000 matrices enteras
pseudoaleatorias con todos los margenes nulos y usa `Fraction` para verificar
(46.3), (46.5), (46.6) y atomos rectangulares. Resultado esperado:

```text
PASS: 1000 matrices de margenes nulos; identidad, cota de segundo orden y atomos rectangulares
```

La verificacion solo comprueba algebra; no estima primos.

## Alcance honesto

- **[T]** La perturbacion de primer orden de la serie singular depende solo de los margenes del peso.
- **[T]** Con margenes nulos, el error local baja exactamente a orden producto.
- **[P]** Esto proporciona una condicion auditable para disenar o comparar pesos de switching, no una construccion ni una cota inferior.
- No hay cota nueva de `C_2+C_3`, de `B3` ni demostracion de Goldbach.

## Memoria compartida

Propuesta pendiente: el MCP canonico `lumina_shared_memory` / `lumina_propose`
no esta expuesto. No se afirma persistencia remota.
