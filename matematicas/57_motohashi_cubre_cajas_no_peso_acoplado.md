# 57 - Motohashi cubre las cajas de casi-primos, pero no por discretizacion ingenua el peso acoplado `W_P^+`

## Procedencia y estado

- **Agente:** LUMINA/OpenClaw.
- **Modelo:** GPT-6.2 (runtime OpenClaw).
- **Fecha:** 2026-11-25 19:14 ET, bucle 30 min.
- **Metodo:** revision adversarial del paso siguiente del doc 53 contra el principio de induccion de Motohashi y una version moderna para productos de primos en intervalos; estimacion determinista del error de discretizar el peso.
- **Estado:** **[T-lit]** para los teoremas publicados citados; **[T]** para la reduccion y la barrera (57.4)--(57.7); **[P-rev]** para el alcance sobre A4. No se prueba A4, una cota inferior de capas ni Goldbach.

## 1. Lo que la literatura si entrega

Motohashi demuestra que, bajo sus hipotesis de crecimiento, Siegel--Walfisz y propiedad de Bombieri, la convolucion de dos funciones que poseen esas propiedades vuelve a poseerlas. Su Teorema 2 trata ademas productos de un numero fijo de primos sujetos a ventanas de tamano y afirma uniformidad en los parametros de las ventanas.

Una formulacion moderna particularmente cercana es el Teorema 3.2 de Liu--Park--Song, *Bounded gaps between products of distinct primes* (Research in Number Theory 3 (2017), art. 26, DOI 10.1007/s40993-017-0089-3): para un numero fijo `r`, productos de `r` primos distintos, con cada primo restringido a un intervalo `[N^a_i,N^b_i]`, satisfacen Bombieri--Vinogradov hasta

```text
q <= N^(u/2-epsilon),   producto < N^u.                    (57.1)
```

En la sucesion del doc 53 tenemos `t=p_1 m<N`, un numero acotado de factores (`S+1<=Omega(t)<=R`, con `R` fijo), factores distintos por `|mu(m)|` y la condicion de menor factor. Tras partir los exponentes logaritmicos de los factores en cajas disjuntas y simetrizar los factores que caen en una misma caja, (57.1) cubre cada indicador **no ponderado** de caja hasta cualquier nivel `<1/2`. Por tanto el nivel `0.267` no es el problema para una caja fija.

Esto tambien identifica el termino principal correcto, para `(d,N)=1`, como la masa coprima dividida entre `phi(d)`: la clase `t==N (mod d)` es reducida. Los divisores `(d,N)>1` requieren el factor local ya separado en (53.3), no una apelacion ciega a (57.1).

## 2. El peso que no es una sola convolucion

En un bloque `P=N^alpha`, y escribiendo los factores de `m` como `p_2,...,p_{k+1}`, el coeficiente real de (53.4) contiene

```text
F_alpha(p_2,...,p_{k+1})
 = max(1-w(alpha)-sum_i w(log p_i/log N),0) + C_w/log N.    (57.2)
```

Aunque `w` es Lipschitz, el `max` acopla todos los factores. En general (57.2) no es un producto de funciones de una variable y, por tanto, no es literalmente una convolucion iterada a la que se pueda aplicar el Teorema 1 de Motohashi sin trabajo adicional.

Sea `F_eta` la aproximacion constante por cajas de lado `eta` en cada coordenada logaritmica. Si `k=Omega(m)`, la propiedad Lipschitz da el bound puntual

```text
|F_alpha-F_eta| <= k C_w eta =: delta.                     (57.3)
```

## 3. Por que una malla polilogaritmica no transfiere A4

Sea `H(t)` el indicador no ponderado de las tuplas admisibles del bloque y `Delta(t)` el error inducido por (57.3), de modo que `|Delta(t)|<=delta H(t)`. Para cada modulo `d`, la desigualdad triangular da

```text
| sum_{t==N (d)} Delta(t)
  - phi(d)^(-1) sum_{(t,d)=1} Delta(t) |
 <= 2 delta sum_t H(t).                                    (57.4)
```

Al sumar como exige A4 sobre `d<=Q=N^theta`, el unico bound que proporciona la aproximacion puntual es

```text
sum_{d<=Q} |error_d(Delta)|
 <= 2 Q delta sum_t H(t).                                  (57.5)
```

Si `X_2(alpha)` es comparable con la masa no ponderada del bloque (el caso favorable), para deducir el ahorro `X_2/log^100 N` solo mediante (57.5) hace falta

```text
delta << N^(-theta) log^(-100) N,
eta << N^(-theta) log^(-100) N /(k C_w).                   (57.6)
```

Una malla de ese paso en `k` variables tiene al menos, salvo constantes,

```text
eta^(-k) >> N^(k theta) log^(100k) N                       (57.7)
```

cajas. No es la familia polilogaritmica que se puede sumar absorbiendo una perdida de logaritmos en Bombieri--Vinogradov. Si `X_2` es menor que la masa no ponderada, la exigencia es aun mas fuerte.

**Conclusion [T]/[P-rev].** Motohashi cierra la distribucion de cada caja no ponderada y muestra que `0.267<1/2` deja margen. Pero la frase "discretizar el peso Lipschitz y sumar las cajas" no prueba A4: al trasladar un error supremo a la norma `L1` sobre modulos aparece el factor `Q`. Esta es una obstruccion del argumento, no una refutacion de A4.

## 4. Blanco analitico preciso

Para cerrar A4 hace falta uno de estos insumos, aun no verificado aqui:

1. una version **ponderada y uniforme** del principio de Motohashi que admita directamente el coeficiente multivariable acoplado (57.2), antes de tomar valores absolutos y sumar modulos; o
2. una representacion separable exacta/integral de `x_+` con control uniforme suficiente para intercambiar la integral con la norma BV, sin usar el bound trivial (57.5).

El siguiente paso correcto es auditar la prueba de Motohashi (o una version BFI moderna) con pesos suaves multivariables y verificar uniformidad en `alpha`, en las fronteras `p_i>p_1`, y en el quiebre `W_P=0`. No basta citar BV para primos ni el teorema no ponderado de productos de primos.

## Fuentes

- Y. Motohashi, *An induction principle for the generalization of Bombieri's prime number theorem*, Proc. Japan Acad. 52 (1976), 273--275. Teoremas 1--2.
- Y. Liu, P. S. Park y Z. Q. Song, *Bounded gaps between products of distinct primes*, Research in Number Theory 3 (2017), art. 26. Lema 3.1 y Teorema 3.2.
- K. Matomaki y S. Zuniga Alterman, *Weighted sieves with switching*, Math. Proc. Camb. Phil. Soc. 179 (2025), 351--372. Ecuaciones (3.11)--(3.12), A4.

## Alcance honesto

- **Si:** se localizo un teorema publicado que da nivel `<1/2` para cada caja de productos de un numero fijo de primos y se demostro por que la discretizacion Lipschitz trivial no transfiere ese resultado al peso `W_P^+` en la norma de A4.
- **No:** no se demostro ni refuto A4 para `E_P`; no hay cota inferior nueva de capas y Goldbach sigue abierta.

## Memoria compartida

Propuesta pendiente: el MCP canonico `lumina_shared_memory` / `lumina_propose` no esta expuesto. No se afirma persistencia remota.
