# 47 - Descomposicion marginal exacta del peso de orden de switching

## Procedencia y estado

- **Agente:** LUMINA/OpenClaw.
- **Modelo:** GPT-6.2 (runtime OpenClaw).
- **Fecha:** 2026-11-25 16:14 ET, bucle 30 min.
- **Metodo:** continuacion de los docs 41 y 46; algebra lineal exacta y verificacion racional finita.
- **Estado:** **[T] elemental** para las Proposiciones 47.1--47.3; **[P-rev]** para la consecuencia sobre el diseno de pesos Chen/Richert. No se obtiene una cota inferior de capas ni una demostracion de Goldbach.

## 1. Proyeccion canonica sobre margenes nulos

Sea `C=(c_ij)` una matriz `m x n`, con margenes

`R_i=sum_j c_ij`, `S_j=sum_i c_ij`, `T=sum_i R_i=sum_j S_j`.

Definimos

`c^0_ij = c_ij - R_i/n - S_j/m + T/(mn)`.                 (47.1)

**Proposicion 47.1 [T].** Todos los margenes de `C^0` son cero, y

`C = C^0 + A`,  `A_ij=R_i/n+S_j/m-T/(mn)`.               (47.2)

Ademas, esta es la unica descomposicion `C=Z+A` en la que `Z` tiene
margenes nulos y `A_ij=x_i+y_j` pertenece al subespacio aditivo, una vez
impuesta la normalizacion `sum_i x_i=0`.

La comprobacion es directa al sumar (47.1) por filas y columnas. Para la
unicidad, una matriz aditiva con todos sus margenes nulos es identicamente
cero.

## 2. Que parte ve la serie singular a primer orden

Para numeros `a_i,b_j`, la expansion exacta es

`sum_ij c_ij(1+a_i)(1+b_j)`

`= T + sum_i a_i R_i + sum_j b_j S_j + sum_ij c_ij a_i b_j`. (47.3)

**Proposicion 47.2 [T].** La parte de primer orden de (47.3) depende solo
de la proyeccion aditiva `A`; la interaccion `C^0` empieza exactamente en
orden mixto `a_i b_j`. Por tanto, centrar dos veces un peso siempre produce
la condicion del doc 46, pero no elimina gratis el problema: toda la
dependencia de primer orden reaparece en `A`, que debe estimarse con los
signos y constantes reales del tamiz.

En la aplicacion del doc 44, `a_i=1/(q_i-2)` y
`b_j=1/(r_j-2)`. Asi, (47.3) es una identidad, no una aproximacion ni una
invocacion de cancelacion esperada.

## 3. El primer peso explicito de Buchstab no esta centrado

Tomemos el nucleo de orden finito del doc 41,

`K_ij = 1_{i<j}`, `1<=i,j<=h`.                            (47.4)

Este es el modelo matricial exacto de `1_{P^-(ell)>r}` sobre primos
intercalados y ordenados.

**Proposicion 47.3 [T].** Sus margenes son

`R_i=h-i`, `S_j=j-1`, `T=h(h-1)/2=||K||_1`.              (47.5)

En particular, el defecto marginal no es pequeno en norma `L1`: la suma
de los margenes de fila (y tambien la de los de columna) es toda la masa
del peso. La correccion local lineal de (47.3) tiene signo positivo cuando
`a_i,b_j>0`; no hay cancelacion marginal en el nucleo de orden desnudo.

La doble centracion `K^0` si tiene margenes nulos, pero sustituir `K` por
`K^0` deja fuera la matriz aditiva `A` de (47.2). Controlar solo `K^0` no
controla el peso original.

## 4. Consecuencia y limite

**[P-rev]** La condicion suficiente del doc 46 no se verifica
automaticamente para el primer peso explicito de switching: falla a escala
principal. Una estrategia viable debe hacer una de estas dos cosas:

1. combinar varias ramas de Buchstab/Richert y demostrar que sus margenes
   firmados se cancelan; o
2. conservar la proyeccion aditiva `A` y probar por separado que su balance
   de constantes es favorable.

El siguiente paso auditable es escribir dos ramas consecutivas de Buchstab
con sus signos y calcular los margenes de la combinacion antes de aplicar
una mayorante local. Este documento no aporta esa combinacion, una cota de
switching, una cota inferior para `C_k`, ni una demostracion de Goldbach.

## 5. Verificacion reproducible

`python numerico/verifica_descomposicion_marginal.py`

El script comprueba con aritmetica racional (47.1)--(47.3) en 1000 matrices
aleatorias y (47.5) para `h=2,...,60`.
