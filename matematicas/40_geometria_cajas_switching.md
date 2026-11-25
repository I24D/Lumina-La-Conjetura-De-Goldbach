# 40 - Geometria exacta de las cajas de switching fuera de Bombieri--Vinogradov

## Procedencia y estado

- **Agente:** LUMINA/OpenClaw.
- **Modelo:** identificador exacto no expuesto por el runtime.
- **Fecha:** 2026-11-25 12:44 ET, bucle 30 min.
- **Metodo:** parametrizacion por exponentes de la forma bilineal aislada en el doc 39 y verificacion racional exhaustiva en una malla finita.
- **Estado:** **[T] elemental** para las Proposiciones 40.1--40.3; **[P-rev]** para la consecuencia metodologica. No se obtiene una cota inferior de capas, una estimacion de switching ni una demostracion de Goldbach.

## 1. Parametros de una caja

Fijemos una caja de la forma (39.9)

```text
q ~ N^delta,   r ~ N^rho,   ell ~ N^lambda,
q r ell ~ N.
```

Ignorando solamente factores `N^o(1)`, se tiene

```text
delta + rho + lambda = 1.                         (40.1)
```

Podemos ordenar las dos variables de switching como `r <= ell`, es decir,
`rho <= lambda`. El control directo de progresiones por Bombieri--Vinogradov
al fijar `q` y `r` requiere

```text
q r <= N^(1/2-o(1)),  o equivalentemente rho <= 1/2-delta.  (40.2)
```

Por tanto la region que el switching debe tratar, pero que BV ya no ve, es

```text
1/2-delta < rho <= (1-delta)/2,
lambda = 1-delta-rho.                              (40.3)
```

## 2. La cunha critica tiene anchura exacta delta/2

**Proposicion 40.1 [T].** Para todo `0<delta<1/2`, las cajas ordenadas
`r<=ell` que quedan fuera de BV forman el intervalo (40.3), de anchura
exacta

```text
(1-delta)/2 - (1/2-delta) = delta/2.               (40.4)
```

En particular, la region ausente no es un borde logaritmico cuando
`q=N^delta`: tiene anchura de potencia positiva. Al crecer el factor fijado
`q`, la cunha critica se ensancha linealmente en `delta`.

*Prueba.* La cota superior de `rho` es exactamente `rho<=lambda` junto con
(40.1). La cota inferior estricta es la negacion de (40.2). Restar los
extremos da (40.4). QED.

Ejemplos representativos:

```text
delta=1/12:  5/12 < rho <= 11/24,   ancho 1/24;
delta=1/8:   3/8  < rho <= 7/16,    ancho 1/16;
delta=1/6:   1/3  < rho <= 5/12,    ancho 1/12;
delta=1/4:   1/4  < rho <= 3/8,     ancho 1/8.
```

## 3. Ambas agrupaciones de modulo quedan fuera de BV

**Proposicion 40.2 [T].** En el interior de (40.3),

```text
delta+rho > 1/2,
delta+lambda >= delta+rho > 1/2.                   (40.5)
```

Por ello ni agrupar `qr` como modulo ni intercambiar `r` y `ell` para usar
`q ell` como modulo devuelve la caja al rango de BV.

*Prueba.* La primera desigualdad es la definicion de la region. Como
`rho<=lambda`, la segunda es inmediata. QED.

Esto corrige una posible salida falsa del doc 39: el ordenamiento de las dos
variables no recupera nivel; la segunda agrupacion solo produce un modulo
mayor.

## 4. Las variables son sub-raiz, pero ambos modulos compuestos son supra-raiz

**Proposicion 40.3 [T].** En el interior de la cunha,

```text
rho < 1/2,   lambda < 1/2,
rho+lambda = 1-delta > 1/2.                        (40.6)
```

Ademas `lambda` recorre

```text
(1-delta)/2 <= lambda < 1/2.                       (40.7)
```

Asi, `r` y `ell` son individualmente sub-raiz, pero su producto tiene
longitud `N^(1-delta)` y el producto de cualquiera de ellas con `q` supera
la raiz. Esta es precisamente una caja bilineal genuina: no es una
progresion individual restante que BV pueda absorber.

*Prueba.* De `rho>1/2-delta` y (40.1) se obtiene `lambda<1/2`; luego
`rho<=lambda<1/2`. La igualdad de la suma es (40.1), y (40.7) resulta de
los dos extremos de (40.3). QED.

## 5. Consecuencia para el siguiente insumo

**[P-rev]** Un intento de adaptar Chen/Timofeev a la sucesion reflejada debe
enunciar una estimacion uniforme para toda la cunha (40.3), no solo para el
borde balanceado `rho=lambda`. Tampoco basta aplicar BV dos veces despues de
intercambiar variables, por la Proposicion 40.2.

El siguiente blanco acotado es escribir el peso real de una identidad de
switching sobre estas cajas y comprobar si es bien factorizable; este
documento no afirma esa propiedad ni aporta cancelacion en
`sum a_r b_ell Lambda(N-q r ell)`.

## 6. Verificacion finita

`numerico/verifica_geometria_switching.py` recorre exactamente, con
`fractions.Fraction`, una malla de `delta,rho` de denominador 240. Verifica
las equivalencias (40.1)--(40.7), los cuatro ejemplos y que intercambiar las
variables no devuelve ninguna caja critica al rango BV.

Resultado: `PASS: 3540 cajas racionales; 4 ejemplos; 0 rescates por intercambio`.

La verificacion solo detecta errores algebraicos; no es evidencia de una
cota analitica ni de Goldbach.

## Alcance honesto

- **[T]** Se delimita exactamente la region de exponentes fuera de BV.
- **[T]** Se prueba que el intercambio de las variables no recupera nivel.
- **[P-rev]** Se fija el dominio minimo que una estimacion de switching debe cubrir.
- No se ha probado ninguna estimacion nueva para esa region; Goldbach sigue abierta.

## Memoria compartida

Propuesta pendiente: el MCP canonico `lumina_shared_memory` / `lumina_propose`
no esta expuesto. No se afirma persistencia remota.
