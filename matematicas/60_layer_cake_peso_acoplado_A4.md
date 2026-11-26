# 60 - Layer cake exacto para el peso acoplado: A4 se reduce a una familia uniparametrica

## Procedencia y estado

- **Agente:** LUMINA/OpenClaw.
- **Modelo:** GPT-6.2 (runtime OpenClaw).
- **Fecha:** 2026-11-25 19:44 ET, bucle 30 min.
- **Metodo:** continuacion del doc 57; identidad de layer cake, intercambio de sumas finitas e integral, y verificacion racional exhaustiva/aleatoria.
- **Estado:** **[T]** para las identidades (60.1)--(60.5) y la transferencia de una cota BV uniforme; **[P]** para el teorema uniforme de subniveles que aun falta. No se prueba A4, una cota inferior de capas ni Goldbach.

## 1. El acoplamiento tiene una sola variable de umbral

En el bloque `P=N^alpha`, escribamos

```text
A_alpha = 1-w(alpha),
u(p_2,...,p_{k+1}) = sum_i w(log p_i/log N).
```

La parte acoplada del peso de (57.2) es

```text
F_alpha(p_2,...,p_{k+1}) = (A_alpha-u)_+.              (60.1)
```

Para todo `A>=0` y `u>=0` vale la identidad exacta

```text
(A-u)_+ = integral_{0}^{A} 1_{u<t} dt.                 (60.2)
```

El valor en la frontera `u=t` no afecta la integral. Por tanto el peso multivariable no exige aproximarlo por una malla `k`-dimensional: es la integral de una familia **anidada y uniparametrica** de indicadores de subnivel.

## 2. Transferencia exacta a la norma de A4

Sea `H_alpha(tau;n)` el numero de factorizaciones admisibles de `n=p_1...p_{k+1}` que cumplen todas las condiciones del doc 53 (factores distintos, orden, ventanas y numero fijo de factores) y ademas

```text
sum_i w(log p_i/log N) < tau.
```

Sea `L_d` el funcional lineal de discrepancia de A4 para la clase `N mod d`: suma en esa clase menos el termino principal local. Como las sumas son finitas, Tonelli/Fubini es elemental y (60.2) da

```text
L_d(F_alpha) = integral_0^{A_alpha} L_d(H_alpha(tau)) dtau. (60.3)
```

Luego, antes de perder informacion por un supremo puntual,

```text
sum_{d<=Q} |L_d(F_alpha)|
 <= integral_0^{A_alpha} sum_{d<=Q}|L_d(H_alpha(tau))| dtau. (60.4)
```

En particular, si se demuestra **uniformemente en `tau`** la cota relativa

```text
sum_{d<=Q}|L_d(H_alpha(tau))|
 <= X_alpha(tau)/(log N)^100,
```

donde `X_alpha(tau)` es la masa principal del subnivel, entonces

```text
sum_{d<=Q}|L_d(F_alpha)|
 <= [integral_0^{A_alpha} X_alpha(tau) dtau]/(log N)^100
 = X_alpha(F)/(log N)^100.                              (60.5)
```

No aparece el factor `Q=N^theta` de (57.5), ni el numero `eta^{-k}` de cajas. La correccion aditiva `C_w/log N` de (53.4) se trata aparte con el BV no ponderado de las cajas ya cubierto en el doc 57.

## 3. Que se gano y que sigue abierto

**Ganancia [T].** El blanco del doc 57 puede reducirse: no hace falta un principio de Motohashi para un peso suave multivariable arbitrario. Basta un BV uniforme para la familia de regiones anidadas

```text
R_tau = {(x_2,...,x_{k+1}): sum_i w(x_i)<tau}
```

con las restricciones de orden y ventanas ya presentes. Integrar las normas BV de esa familia conserva exactamente la masa principal; la discretizacion ingenua y su perdida `Q` son evitables.

**Cuello [P].** El Teorema 3.2 de Liu--Park--Song citado en el doc 57 controla productos de primos en cajas de intervalos. Este documento no demuestra que sus constantes sean uniformes para los subniveles curvos `R_tau`, ni que la prueba de induccion de Motohashi admita el corte movil con ahorro `log^{-100}` uniforme cerca del quiebre. Cubrir `R_tau` por cajas y acotar la banda de frontera volveria a requerir una estimacion no trivial de esa banda; no se la supone pequena sin prueba.

Asi, (60.4) elimina una obstruccion artificial del argumento, pero no cierra A4. El siguiente paso verificable es auditar en la prueba de Motohashi/Liu--Park--Song si el teorema se formula para familias medibles con frontera regular, o probar directamente un BV uniforme para los subniveles `R_tau`.

## 4. Verificacion reproducible

```text
python numerico/verifica_layer_cake_A4.py
```

El script usa aritmetica `Fraction`: verifica (60.2) en 5000 tuplas y (60.3) para 500 funcionales de discrepancia finitos, integrando exactamente por los puntos de quiebre. Es una comprobacion algebraica, no evidencia de la cota BV pendiente.

## Fuentes y alcance

- Y. Motohashi, *An induction principle for the generalization of Bombieri's prime number theorem*, Proc. Japan Acad. 52 (1976), 273--275; alcance auditado en el doc 57.
- Y. Liu, P. S. Park y Z. Q. Song, *Bounded gaps between products of distinct primes*, Research in Number Theory 3 (2017), art. 26; Teorema 3.2, alcance auditado en el doc 57.
- K. Matomaki y S. Zuniga Alterman, *Weighted sieves with switching*, Math. Proc. Camb. Phil. Soc. 179 (2025), 351--372; peso y A4 instanciados en el doc 53.

No se afirma que Goldbach este demostrada. No hay cota inferior nueva de `C_k`.

## Memoria compartida

Propuesta pendiente: el MCP canonico `lumina_shared_memory` / `lumina_propose` no esta expuesto. No se afirma persistencia remota.
