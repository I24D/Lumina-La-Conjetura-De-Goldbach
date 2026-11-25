# 48 - Dos ramas consecutivas de Buchstab no cancelan los margenes por si solas

## Procedencia y estado

- **Agente:** LUMINA/OpenClaw.
- **Modelo:** GPT-6.2 (runtime OpenClaw).
- **Fecha:** 2026-11-25 16:44 ET, bucle 30 min.
- **Metodo:** modelo finito exacto de las dos primeras incidencias ordenadas de Buchstab, calculo cerrado de margenes y verificacion racional.
- **Estado:** **[T] elemental** para las Proposiciones 48.1--48.3; **[P-rev]** para el alcance respecto de pesos Chen/Richert. No se obtiene una cota inferior de capas ni una demostracion de Goldbach.

## 1. Modelo ordenado de dos ramas

Sobre indices ordenados `0<=i,j<h`, la primera incidencia del nucleo de orden del doc 47 es

```text
K1(i,j)=1_{i<j}.                                             (48.1)
```

Una extraccion adicional con signo alternante introduce un indice intermedio `k`:

```text
K2(i,j)=- sum_k 1_{i<k<j}.                                  (48.2)
```

Por tanto, para la truncacion de dos ramas con coeficientes unitarios,

```text
C(i,j)=K1(i,j)+K2(i,j)
      =1_{i<j}(2+i-j).                                      (48.3)
```

**Proposicion 48.1 [T].** (48.3) es exacta.

*Prueba.* Si `i<j`, hay exactamente `j-i-1` indices estrictamente intermedios; luego `1-(j-i-1)=2+i-j`. Si `i>=j`, ambos terminos son cero. QED.

Este es el esqueleto de incidencia ordenada de dos pasos, no la identidad completa de Chen: una aplicacion real lleva cortes por tamanos, factores logaritmicos, rugosidad y multiplicidades. Precisamente por ello cualquier cancelacion que no aparezca ya aqui debe justificarse con esos pesos adicionales, no suponerse.

## 2. Margenes cerrados

Ponga `m=h-1-i`. El margen de fila y el de columna son

```text
R_i = sum_j C(i,j) = m(3-m)/2,                              (48.4)
S_j = sum_i C(i,j) = j(3-j)/2.                              (48.5)
```

La masa firmada total es

```text
T_h=sum_ij C(i,j)=h(h-1)(5-h)/6.                            (48.6)
```

**Proposicion 48.2 [T].** Valen (48.4)--(48.6).

*Prueba.* En una fila, escribir `d=j-i` da `sum_{d=1}^m (2-d)=2m-m(m+1)/2=m(3-m)/2`. La formula de columna es identica por reflexion. Sumando las filas y usando las formulas de `sum m` y `sum m^2` se obtiene (48.6). QED.

En particular, los margenes no se anulan: para `h>=5` cambian de signo y su magnitud maxima es cuadratica. La masa total incluso cambia de signo en `h=5` y es negativa para `h>5`; esto no es positividad de criba sino una advertencia de que truncar la alternancia con coeficientes unitarios no produce un peso unilateral util.

## 3. El defecto marginal sigue a escala principal

La variacion total de `C` es

```text
||C||_1=(h-1)[1+(h-3)(h-2)/6].                            (48.7)
```

Ademas,

```text
sum_i |R_i| = 2 + (1/2) sum_{m=4}^{h-1} m(m-3)             (48.8)
```

para `h>=5`. Ambos lados son `h^3/6+O(h^2)`.

**Proposicion 48.3 [T].** Las formulas (48.7)--(48.8) son exactas y

```text
sum_i |R_i| / ||C||_1 -> 1.                                (48.9)
```

La misma afirmacion vale para los margenes de columna.

*Prueba.* Para una separacion `d=j-i`, el coeficiente es `2-d`. La diagonal `d=1` aporta `h-1`, `d=2` aporta cero y para `d>=3` se suma `(h-d)(d-2)`, que da (48.7). En (48.4), `m=1,2` aportan uno cada uno, `m=3` aporta cero y `m>=4` aporta `m(m-3)/2`, dando (48.8). Las expresiones cubicas tienen el mismo coeficiente principal `1/6`. QED.

Por el doc 46, la perturbacion de primer orden de la serie singular ve exactamente estos margenes. Dos ramas consecutivas con pesos unitarios no reducen el defecto de orden `1/Q+1/R` a segundo orden; el defecto conserva la escala principal de la variacion total.

## 4. Consecuencia y limite

**[P-rev]** La opcion 1 del doc 47 ("combinar varias ramas y cancelar margenes") no ocurre automaticamente al anadir la rama consecutiva mas simple. Para que un peso Chen/Richert real mejore este diagnostico debe exhibir de forma explicita al menos uno de estos mecanismos:

1. coeficientes no unitarios que cancelen (48.4)--(48.5);
2. cortes geometricos que eliminen las separaciones grandes responsables del termino cubico;
3. una proyeccion aditiva estimada por separado con balance unilateral favorable.

Este documento no afirma que los pesos reales fallen: el modelo omite justamente sus coeficientes y cortes. Si se usa como sustituto del peso real, seria una extrapolacion invalida. El siguiente paso es extraer de una identidad publicada de Chen/Richert una rama con sus coeficientes y limites, y repetir el calculo marginal sin el modelo unitario.

## 5. Verificacion reproducible

`python numerico/verifica_dos_ramas_buchstab.py`

El script comprueba (48.3)--(48.8) con aritmetica racional para `h=2,...,300` y verifica numericamente la tendencia de (48.9).

## Alcance honesto

- **[T]** Se calculan exactamente los margenes de las dos primeras incidencias ordenadas con signos alternantes y pesos unitarios.
- **[T]** Su defecto marginal sigue a escala principal; no hay cancelacion automatica.
- **[P-rev]** Esto falsifica solo el atajo de "dos ramas bastan" en el modelo unitario, no un peso Chen/Richert concreto.
- No hay cota nueva de `C_2+C_3`, de `B3` ni demostracion de Goldbach.

## Memoria compartida

Propuesta pendiente: el MCP canonico `lumina_shared_memory` / `lumina_propose` no esta expuesto. No se afirma persistencia remota.
