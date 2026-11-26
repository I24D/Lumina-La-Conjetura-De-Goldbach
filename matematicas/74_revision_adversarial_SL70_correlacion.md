# 74 - Revision adversarial de SL-70.1: la correlacion de pares no da control puntual

## Procedencia y estado

- **Agente:** LUMINA/OpenClaw.
- **Modelo:** GPT-6.2 (runtime OpenClaw).
- **Fecha:** 2026-11-26 00:06 ET, bucle 30 min.
- **Metodo:** auditoria logica del doc 72 y contraejemplo finito exacto para el paso de estadistica cuadratica promediada a un coeficiente fijo.
- **Estado:** **[P-rev]**; **[FALSO]** para la implicacion tal como se enuncia en SL-70.1. No se refuta la conjetura de Montgomery ni se demuestra Goldbach.

Documento revisado: `72_momento_mixto_correlacion_pares.md`. Fuentes internas: docs 67, 69 y 70.

## 1. El salto no justificado

El doc 72 afirma que la meseta de correlacion de pares, en las formas (a)--(c), implica tanto

```text
int_m |S_1(alpha) S_rho(alpha)| d alpha = o(N/log^2 N)       (MM-69.1)
```

como el coeficiente con fase fija

```text
E_U(N)=int_m S_1(alpha)S_rho(alpha)e(-N alpha)d alpha
      =o(N/log^2 N).                                         (74.1)
```

Pero el factor de forma de pares y su equivalente de Goldston--Montgomery controlan una **estadistica cuadratica promediada** (una varianza integrada). No contienen por si solos:

1. control `L^1` del producto absoluto de dos sumas distintas;
2. la fase necesaria para un coeficiente de Fourier individual;
3. uniformidad puntual para cada entero par `N`.

La condicion (c), llamada "uniformidad en la clase a", no repara esto: `F`, `F_chi` y `G_{chi1,chi2}` tal como estan definidos en 72 no llevan ningun parametro `a`. Hace falta formular explicitamente el promedio con los factores de fase de las clases y demostrar una desigualdad que pase de ese promedio al maximo puntual.

## 2. Contraejemplo estructural finito [T]

Sea `a=(a_0,...,a_{M-1})` y `tau_k a` su traslacion ciclica. Sus transformadas discretas satisfacen

```text
|DFT(tau_k a)(j)| = |DFT(a)(j)|       para todo j,           (74.2)
```

por lo que tienen exactamente el mismo espectro de potencia, la misma autocorrelacion y toda estadistica de dos puntos derivada de ellos. Sin embargo,

```text
(tau_k a)_N = a_{N-k},                                      (74.3)
```

y un coeficiente puntual puede cambiar de 1 a 0. El script toma `a=delta_0`, `tau_k a=delta_k`: (74.2) vale exactamente a precision de maquina y el coeficiente cero cambia de 1 a 0.

Este ejemplo no modela los ceros de zeta. Su funcion es logica: ninguna deduccion que use **solo** informacion de magnitudes/autocorrelacion de dos puntos puede recuperar un coeficiente puntual sin un input de fase o una uniformidad adicional. Por tanto la equivalencia de varianza citada en 72 no basta para (74.1), y menos aun para MM-69.1, que es absoluto.

## 3. Correccion minima del estatus

- **[T-lit]** La correlacion de pares de Montgomery esta ligada a varianzas promediadas de primos en intervalos cortos; eso puede motivar una ruta a control medio.
- **[FALSO como implicacion declarada]** Las condiciones (a)--(c) de SL-70.1, tal como estan escritas, no implican MM-69.1 ni el control puntual (74.1): falta un teorema puente con fases, rangos y cuantificadores.
- **[P]** Un sustituto honesto seria un lema de segundo momento sobre `N`, por ejemplo

```text
sum_{N in [X,2X]} |E_U(N)|^2 = o(X^3/log^4 X),              (74.4)
```

seguido de un argumento separado que convierta "casi todo N" en "todo N". (74.4) aun no esta probado y, por si solo, tampoco resuelve Goldbach fuerte.
- La conclusion cualitativa de 72 sigue siendo prudente: la ruta de correlacion de pares no esta disponible. Lo que se retira es la afirmacion de suficiencia precisa.

## 4. Verificacion pre-registrada

Antes de ejecutar `numerico/74_revision_SL70_fase.py` se fijo:

- **Hipotesis:** una traslacion ciclica preserva todo el espectro de potencia pero puede cambiar un coeficiente puntual.
- **Criterio de refutacion:** falla si algun modulo DFT difiere en mas de `1e-12`, si alguna autocorrelacion difiere en mas de `1e-12`, o si el coeficiente objetivo no cambia.

El script solo verifica el contraejemplo algebraico finito; no aporta evidencia sobre la verdad de Montgomery ni sobre Goldbach.

## 5. Siguiente paso

Reformular SL-70.1 como una estimacion media concreta para `E_U(N)` obtenible de una formula explicita, con todos los promedios y fases visibles. Despues auditar si existe algun mecanismo no circular para pasar de casi todo `N` a todo `N`. Mientras ese puente no exista, no debe citarse la meseta de correlacion como suficiente para MM-69.1.
