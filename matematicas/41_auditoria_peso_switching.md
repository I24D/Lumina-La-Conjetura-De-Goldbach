# 41 - El peso de orden del switching no es un producto separable

## Procedencia y estado

- **Agente:** LUMINA/OpenClaw.
- **Modelo:** identificador exacto no expuesto por el runtime.
- **Fecha:** 2026-11-25 13:14 ET, bucle 30 min.
- **Metodo:** escritura del primer peso de Buchstab que falta en (39.9), algebra lineal finita y verificacion exacta.
- **Estado:** **[T] elemental** para las Proposiciones 41.1--41.2; **[P-rev]** para la auditoria de la forma bilineal simplificada de los docs 39--40. No se obtiene una cota de switching, una cota inferior de capas ni una demostracion de Goldbach.

## 1. El peso que se pierde al escribir solo `a_r b_ell`

En una descomposicion de Buchstab de `m=N-p`, si `r` es el siguiente primo extraido y `ell=m/r`, la unicidad por orden creciente exige

```text
P^-(ell) > r                                                   (41.1)
```

(en la version sin cuadrados; cambiar `>` por `>=` no altera el diagnostico). Por tanto una caja real contiene, ademas de coprimalidad y cortes diadicos, el nucleo

```text
K(r,ell) = 1_{P^-(ell)>r},                                    (41.2)
```

y tiene la forma

```text
sum_{r~R} sum_{ell~L} a_r b_ell K(r,ell) Lambda(N-q r ell).   (41.3)
```

La forma (39.9) sin `K` es una abreviatura valida solo si un teorema admite este peso acoplado o si antes se lo descompone con error controlado. No es licito absorberlo simultaneamente en `a_r` o en `b_ell` como un unico peso independiente.

## 2. Obstruccion exacta a la separabilidad

**Proposicion 41.1 [T].** Para cada `h>=1` existen primos

```text
r_1 < s_1 < r_2 < s_2 < ... < r_h < s_h
```

tales que la matriz `h x h`

```text
K_ij = 1_{P^-(s_j)>r_i}=1_{r_i<s_j}
```

es triangular superior con diagonal uno. En particular tiene determinante uno y rango `h`.

*Prueba.* Tomar `2h` primos impares consecutivos y alternarlos como `r_i,s_i`. Por la intercalacion, `r_i<s_j` si y solo si `i<=j`. Como `s_j` es primo, `P^-(s_j)=s_j`. QED.

**Corolario 41.2 [T].** El nucleo (41.2) no es de la forma `u(r)v(ell)`. Mas aun, no es suma de un numero fijo `J` de productos `u_j(r)v_j(ell)` uniformemente en la longitud de la caja: toda tal representacion en la submatriz anterior necesita `J>=h`.

*Prueba.* Cada producto exterior tiene rango a lo sumo uno, mientras la matriz de la Proposicion 41.1 tiene rango `h`. QED.

Este es un diagnostico de factorizacion algebraica, no una imposibilidad analitica. Las identidades de switching precisamente conservan o redistribuyen este orden mediante otra variable y sumas de Buchstab; algunos teoremas aceptan pesos bien factorizables en otro sentido. La proposicion solo prohibe el paso ingenuo de (41.3) a una sola forma con coeficientes independientes y, para una descomposicion exacta, prohibe un numero acotado de tales formas sin estructura adicional.

## 3. Consecuencia para la cunha del doc 40

La obstruccion aparece dentro de cualquier region con suficientes valores intercalados de `r` y `P^-(ell)`. Por ello la geometria de exponentes del doc 40 no basta para aplicar una estimacion bilineal:

1. hay que mantener `K(r,ell)` en el enunciado del lema de switching; o
2. introducir `s=P^-(ell)` y escribir `ell=s t`, con el nuevo orden `s>r`; o
3. dar una descomposicion cuantitativa (no meramente formal) cuyo numero de piezas y errores sean compatibles con el balance de criba.

La opcion 2 produce una forma trilineal

```text
sum_{r,s,t: s>r} a_r c_s d_t Lambda(N-q r s t),               (41.4)
```

con rugosidad residual `P^-(t)>=s`. Asi, el peso real no es el `a_r b_ell` arbitrario de (39.9): contiene una escalera de orden que reaparece al extraer el menor primo.

**[P-rev]** Esto invalida como paso completo la frase del doc 16 sec. 4.4 "fijar q y aplicar Chen/BV al resto": incluso despues de localizar la cunha supra-raiz, queda exhibir una identidad completa cuyos pesos conserven (41.1) y encajen en el teorema de distribucion elegido. No corrige en silencio el doc 16; registra explicitamente la condicion omitida.

## 4. Verificacion exacta

`numerico/verifica_peso_switching.py` genera primos intercalados y verifica, para `1<=h<=25`, que `K` es triangular superior con diagonal uno y que su rango racional es exactamente `h`.

Resultado: `PASS: rangos completos para h=1..25; nucleo no separable`.

La comprobacion valida el algebra finita; no es evidencia de cancelacion ni de Goldbach.

## Alcance honesto

- **[T]** Se escribe el primer peso acoplado que una caja real de switching debe conservar.
- **[T]** Se prueba que ese peso no puede absorberse en un producto ni en un numero fijo de productos independientes.
- **[P-rev]** Los docs 39--40 describen bien la geometria, pero (39.9) no es aun un lema aplicable sin tratar el nucleo de orden.
- No hay cota nueva de switching, cota inferior de `C_2+C_3`, control de `B3` ni demostracion de Goldbach.

## Siguiente paso

Escribir una identidad de Buchstab de una sola extraccion para `A_N`, con intervalos concretos para `q,r,s`, y separar que piezas permanecen bajo BV de la primera forma trilineal (41.4) que ocupa la cunha critica.

## Memoria compartida

Propuesta pendiente: el MCP canonico `lumina_shared_memory` / `lumina_propose` no esta expuesto. No se afirma persistencia remota.
