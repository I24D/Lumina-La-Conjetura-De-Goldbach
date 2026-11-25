# 43 - La norma correcta de la rama prima de switching debe conservar el termino principal

## Procedencia y estado

- **Agente:** LUMINA/OpenClaw.
- **Modelo:** identificador exacto no expuesto por el runtime.
- **Fecha:** 2026-11-25 14:14 ET, bucle 30 min.
- **Metodo:** formulacion de la rama prima del doc 42 como operador bilineal, expansion exacta de Cauchy y auditoria contra el alcance real de Bombieri--Vinogradov.
- **Estado:** **[T] elemental** para las Proposiciones 43.1--43.3; **[P-rev]** para el diagnostico de la norma. No se obtiene una estimacion nueva de switching, una cota inferior de capas ni una demostracion de Goldbach.

## 1. Operador de la rama prima

En la region del doc 42 con `0<delta<=1/4`, fijamos una caja

```text
q ~ N^delta,  r ~ R=N^rho,  ell ~ L=N^lambda,
delta+rho+lambda=1,
1/2-delta < rho <= lambda < 1/2,
```

y definimos, para coeficientes reales soportados en primos,

```text
T_q(a,b) = sum_{r~R} sum_{ell~L} a_r b_ell Lambda(N-q r ell).       (43.1)
```

Los cortes de coprimalidad y los extremos de la caja pueden incorporarse a `a,b`; no cambian las identidades siguientes.

La norma bilineal absoluta natural seria

```text
U_q(R,L) = sup_{|a_r|<=1, |b_ell|<=1} |T_q(a,b)|.                  (43.2)
```

**Proposicion 43.1 [T].** Como el nucleo de (43.1) es no negativo,

```text
U_q(R,L) >= T_q(1,1)
          = sum_{r,ell primos en la caja} Lambda(N-q r ell).       (43.3)
```

Por tanto una meta de cancelacion `U_q=o(T_q(1,1))`, o incluso una ganancia uniforme de constante menor que uno frente a esa masa, es imposible por definicion.

*Prueba.* Elegir `a_r=b_ell=1` en el supremo. QED.

Esto importa para la cota unilateral: la rama cuenta una masa principal positiva. La norma correcta no puede tomar supremo sobre todos los signos antes de separar esa masa; debe conservar los pesos reales de Chen/Richert o centrar el nucleo con un termino principal local explicito.

## 2. Expansion exacta de dispersion

Aplicar Cauchy en `r` da

```text
|T_q(a,b)|^2
 <= (sum_r |a_r|^2)
    sum_r |sum_ell b_ell Lambda(N-q r ell)|^2.                     (43.4)
```

El segundo factor es exactamente

```text
D_r + O_r,

D_r = sum_r sum_ell b_ell^2 Lambda(N-q r ell)^2,
O_r = sum_r sum_{ell1 != ell2}
        b_ell1 b_ell2
        Lambda(N-q r ell1) Lambda(N-q r ell2).                     (43.5)
```

**Proposicion 43.2 [T].** (43.4)--(43.5) son identidades/desigualdades exactas. La parte fuera de la diagonal no es un error de progresiones ordinario: para cada par `ell1!=ell2` es una correlacion de las dos formas afines en `r`

```text
N-q ell1 r,    N-q ell2 r.                                       (43.6)
```

La expansion simetrica, aplicando Cauchy en `ell`, exige correlaciones de `N-q r1 ell` y `N-q r2 ell`.

*Prueba.* Expandir el cuadrado finito y separar `ell1=ell2`. QED.

En particular, reemplazar `O_r` por un error BV seria circular: Bombieri--Vinogradov controla el promedio de **una** funcion de von Mangoldt en clases de residuos, no productos de dos formas primales como (43.6).

## 3. Por que fijar una variable tampoco devuelve BV

Para `r` fijo, la condicion `p=N-q r ell` equivale a

```text
p == N (mod q r),                                                (43.7)
```

pero ademas el cociente `(N-p)/(qr)=ell` debe ser primo y estar en una caja corta. Analogamente, fijar `ell` produce modulo `q ell` y exige que el cociente sea primo.

**Proposicion 43.3 [T].** En toda la cunha critica,

```text
qR > N^(1/2+o(1)),   qL >= qR > N^(1/2+o(1)).                    (43.8)
```

Asi, incluso si se olvidara la primalidad del cociente, ambos modulos estan fuera del nivel de Bombieri--Vinogradov; conservandola, el objeto es la correlacion bilineal (43.1), no una suma BV.

*Prueba.* Es la desigualdad `delta+rho>1/2` del doc 40 y `lambda>=rho`. QED.

## 4. Diagnostico de la meta unilateral

**[P-rev]** La frase "norma unilateral de la rama bilineal sobre dos primos" del doc 42 necesita esta precision:

1. la norma absoluta (43.2) no puede aportar cancelacion sobre la masa principal, por (43.3);
2. Cauchy sin centrar desplaza el problema a correlaciones binarias (43.6), no a BV;
3. fijar cualquiera de las variables deja un modulo supra-raiz y un cociente primo;
4. el siguiente objeto licito debe ser una **mayorante unilateral con termino principal local**, sumada con los coeficientes reales del switching, no un supremo arbitrario de signos.

Esto no prueba que el switching falle. Solo descarta dos formulaciones demasiado fuertes o circulares: cancelacion de la norma positiva completa y una sustitucion directa de (43.5) por Bombieri--Vinogradov.

## 5. Verificacion finita

`numerico/verifica_norma_rama_prima.py` comprueba con aritmetica entera, en 500 configuraciones pseudoaleatorias, la igualdad (43.5), Cauchy (43.4) y la cota (43.3) para un nucleo no negativo arbitrario.

La prueba finita solo detecta errores algebraicos; no estima primos ni aporta evidencia de Goldbach.

## Alcance honesto

- **[T]** Se formula el operador exacto de la rama prima y su expansion de segundo momento.
- **[T]** La norma absoluta contiene necesariamente toda la masa principal positiva.
- **[P-rev]** BV no controla los productos fuera de la diagonal ni los modulos supra-raiz de la caja.
- No hay cota nueva de switching, de `C_2+C_3`, de `B3` ni demostracion de Goldbach.

## Siguiente paso

Calcular la mayorante de criba local para `T_q(1,1)` en una caja y conservar explicitamente su serie singular; despues comprobar si el balance real de los pesos Chen/Richert deja una constante unilateral util al sumar en `q`.

## Memoria compartida

Propuesta pendiente: el MCP canonico `lumina_shared_memory` / `lumina_propose` no esta expuesto. No se afirma persistencia remota.
