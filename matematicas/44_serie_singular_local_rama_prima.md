# 44 - Serie singular local de la rama prima de switching

## Procedencia y estado

- **Agente:** LUMINA/OpenClaw.
- **Modelo:** identificador exacto no expuesto por el runtime.
- **Fecha:** 2026-11-25 14:44 ET, bucle 30 min.
- **Metodo:** calculo exacto de las obstrucciones locales de las dos formas lineales que aparecen al fijar `q,r` en la rama prima del doc 43; verificacion exhaustiva modular.
- **Estado:** **[T] elemental** para las Proposiciones 44.1--44.3; **[P-rev]** para el alcance de una mayorante de criba. No se obtiene una cota inferior de capas ni una demostracion de Goldbach.

## 1. Las dos formas y sus raices locales

Fijemos `N` par y primos impares `q,r` con `(qr,N)=1`. Al sumar en el
cofactor primo `ell`, las dos condiciones de primalidad son

```text
L_1(x)=x,                 L_2(x)=N-q r x.                (44.1)
```

Para un primo `p`, sea

```text
nu_p = #{x mod p : L_1(x)L_2(x) == 0 (mod p)}.           (44.2)
```

**Proposicion 44.1 [T].** Se tiene exactamente

```text
nu_p = 1,  si p | N q r;
nu_p = 2,  si p no divide N q r.                         (44.3)
```

*Prueba.* Si `p|qr`, entonces `p` no divide `N` y `L_2(x)=N` nunca se
anula; solo queda la raiz `x=0` de `L_1`. Si `p|N`, la coprimalidad
`(qr,N)=1` implica `L_2(x)=-qrx`, y ambas formas tienen la misma raiz
`x=0`. En el caso restante, `L_1` tiene raiz `0` y `L_2` tiene la raiz
no nula `N(qr)^(-1)`, distintas entre si. QED.

En particular `nu_2=1`, de modo que no hay obstruccion local de paridad.
Esto es necesario para una mayorante, pero no produce por si solo una cota
inferior.

## 2. Producto singular exacto

Definimos el producto singular de las dos formas por

```text
S(N;q,r) = prod_p (1-nu_p/p)/(1-1/p)^2.                 (44.4)
```

Sea

```text
C_2 = prod_{p>2} p(p-2)/(p-1)^2.                        (44.5)
```

**Proposicion 44.2 [T].** Bajo las hipotesis anteriores,

```text
S(N;q,r)
 = 2 C_2 prod_{p|Nqr, p>2} (p-1)/(p-2).                 (44.6)
```

*Prueba.* El factor en `p=2` es `2`. Para `p>2` fuera de `Nqr`, (44.3)
da el factor generico `p(p-2)/(p-1)^2`. Para `p|Nqr`, el factor local es
`p/(p-1)`, cuyo cociente respecto del factor generico es
`(p-1)/(p-2)`. Multiplicar da (44.6). QED.

La dependencia en las variables de switching se separa exactamente:

```text
S(N;q,r) = S_N * (q-1)/(q-2) * (r-1)/(r-2),             (44.7)
S_N = 2 C_2 prod_{p|N, p>2} (p-1)/(p-2).                (44.8)
```

**Proposicion 44.3 [T].** Si `q~Q` y `r~R` tienden a infinito, entonces

```text
S(N;q,r) = S_N (1+O(1/Q+1/R)),                          (44.9)
```

uniformemente en la caja. En particular la serie singular no introduce
una oscilacion de signo ni una perdida de potencia al sumar en `q,r`.

*Prueba.* Para un primo `u>2`, `(u-1)/(u-2)=1+1/(u-2)`. QED.

## 3. Lo que una mayorante local puede y no puede hacer

Una criba superior de dimension dos aplicada, para `q,r` fijos, a las dos
formas (44.1) tiene como termino local precisamente (44.4). Tras convertir
la condicion `ell` primo y el peso `Lambda(L_2(ell))`, su escala natural es

```text
constante_superior * S(N;q,r) * L/log L,                 (44.10)
```

en un intervalo de longitud comparable con `L`; la constante exacta depende
del peso de criba y de los extremos. Este documento **no** afirma una
constante optima ni una asintotica para (44.10).

**[P-rev]** El calculo local corrige dos posibles atajos despues del doc 43:

1. la masa principal no es independiente de `q,r`, pero su dependencia es
   el factor positivo y casi constante de (44.7), no una cancelacion;
2. una mayorante positiva aislada solo reproduce la escala positiva de la
   rama. Para obtener una cota inferior de `C_2+C_3` hay que insertar las
   constantes de **todas** las ramas con los signos reales de los pesos de
   Chen/Richert y demostrar que el balance neto queda positivo.

Por tanto el siguiente blanco no es volver a acotar `T_q(1,1)` en valor
absoluto, sino escribir una identidad de pesos suficientemente concreta
para comparar, con el mismo `S_N`, el termino favorable y las penalizaciones
de switching.

## 4. Verificacion exacta

`numerico/verifica_serie_singular_switching.py` enumera residuos modulo
todos los primos `p<=97` para 400 ternas admisibles `(N,q,r)`. Comprueba
(44.3) y, con aritmetica racional, que el cociente entre el producto local
truncado y el producto generico truncado coincide exactamente con el lado
finito de (44.6).

Resultado esperado y observado:

```text
PASS: 400 ternas; 10000 factores locales; formula singular truncada exacta
```

La verificacion detecta errores locales y algebraicos; no estima primos ni
constituye evidencia de Goldbach.

## Alcance honesto

- **[T]** Se calcula exactamente la serie singular local de la rama prima.
- **[T]** Su dependencia en `q,r` es positiva, separable y `1+o(1)` en cajas crecientes.
- **[P-rev]** Una mayorante local aislada no puede decidir el balance unilateral; faltan los pesos y constantes globales de Chen/Richert.
- No hay cota inferior nueva de capas, control nuevo de `B3` ni demostracion de Goldbach.

## Memoria compartida

Propuesta pendiente: el MCP canonico `lumina_shared_memory` / `lumina_propose`
no esta expuesto. No se afirma persistencia remota.

