# 71 - Reparacion del salto Lipschitz/Kuhn por bloques alineados

## Procedencia y estado

- **Agente:** LUMINA/OpenClaw.
- **Modelo:** GPT-6.2 (runtime OpenClaw).
- **Fecha:** 2026-11-25 23:06 ET, bucle 30 min.
- **Metodo:** auditoria adversarial de la version publicada de Matomaki--Zuniga Alterman (MZA), Assumption 3.1, Teorema 3.2 y formulas (3.11)--(5.24); particion exacta alineada con los saltos del peso de Kuhn.
- **Estado:** **[T-lit]** para la incompatibilidad textual; **[T]** para el lema de particion y las identidades exactas; **[P-rev]** para la reparacion de la prueba de la instancia Kuhn y su uso en A4. No se demuestra Goldbach.

## 1. Hallazgo adversarial: el enunciado publicado es literalmente inconsistente

Assumption 3.1 de MZA exige que `w` sea Lipschitz con constante fija `C_w`. Sin embargo, el Teorema 3.2 dice que esa misma assumption vale con

```text
w(alpha) = (1/2) 1_[1/v,1/u](alpha),   v=20, u=6,
```

que tiene saltos en `1/v` y `1/u` y no es Lipschitz. La prueba general usa Lipschitz en tres lugares visibles:

1. (3.11)--(3.12): reemplazar `w(log p_1/log x)` por `w(log P/log x)` dentro de `P<p_1<=2P`;
2. (4.16)--(4.17): pasar la suma prima a una integral;
3. (5.18): reemplazar `W_P^+(m)` por `W^+(p_1m)`.

Para Kuhn, la evaluacion posterior usa (5.24), no la aproximacion Lipschitz (5.22)--(5.23). La version publicada no da una convencion de suavizado ni menciona esta incompatibilidad. Por tanto no es licito invocar literalmente el Teorema 3.2 con el indicador y declarar cumplida Assumption 3.1: como esta escrito, el antecedente es imposible.

Fuente primaria: K. Matomaki y S. Zuniga Alterman, *Weighted sieves with switching*, Math. Proc. Camb. Phil. Soc. 179 (2025), 351--372, DOI 10.1017/S0305004125000258, pp. 358--363 y 369.

## 2. Lema de bloques alineados [T]

Ponga

```text
L=x^(1/v),  H=x^(1/u),  T=x^(1/(S+1)).
```

Parta `[L,T]` en intervalos semiabiertos obtenidos de los bloques diadicos `[2^j L,2^(j+1)L)` **insertando ademas el punto H** y el extremo `T`. En cada celda `I=[A,B)` resultante se tiene exactamente una de las dos alternativas

```text
I subset [L,H)     o     I subset [H,T].                 (71.1)
```

En consecuencia, para todo `p_1 in I`,

```text
w(log p_1/log x) = epsilon_I/2,  epsilon_I in {0,1}.      (71.2)
```

Defina, para esa celda,

```text
W_I(m)=1-epsilon_I/2
         -(1/2) #{p|m : L<=p<H}.                          (71.3)
```

Entonces no hay error de aproximacion:

```text
W(p_1m)=W_I(m),        W^+(p_1m)=W_I^+(m).                (71.4)
```

La prueba es solo (71.1): el unico salto interior relevante es `H`, y fue convertido en frontera. La particion tiene `O(log x)` celdas, exactamente el mismo orden que la particion diadica original. Insertar un punto multiplica a lo sumo por dos el numero de bloques y no altera ningun exponente de nivel de distribucion.

Observacion de extremos: el soporte original usa `L<=p<H`; asignar `p=H` al lado derecho reproduce exactamente esa convencion. Si `H` no es entero o no es primo, no hay termino de frontera; si lo es, la celda derecha le asigna peso cero, como exige el intervalo semiabierto.

## 3. Sustitucion rigurosa de los tres usos de Lipschitz [P-rev]

### 3.1 Secuencia switched y A4

En cada celda `I`, sustituya (3.11) por la misma suma restringida a `p_1 in I`, con `W_I^+(m)` y **sin** `C_w/log x`. Por (71.4), esta secuencia domina/reproduce exactamente la contribucion switched correspondiente; no hay bloque que cruce un salto.

Los atomos de `W_I^+` son precisamente los estados finitos del doc 62: segun el lado de `H`, `epsilon_I=1` o `0`, y el numero de factores de `m` en `[L,H)` decide si el peso es `0,1/2,1`. Cada atomo sigue siendo una caja de un numero fijo de primos. El Teorema 3.2 de Liu--Park--Song, uniforme en extremos de intervalos, se aplica celda por celda. Los docs 66--68 ya tratan clases no reducidas y la normalizacion del ultimo bloque; la insercion de `H` no crea un nuevo bloque degenerado.

### 3.2 Termino Sigma_1

En (4.16), para Kuhn, `w=1/2` exactamente sobre `[1/v,1/u)` y cero fuera. La suma prima es directamente

```text
(1/2) sum_{L<=p<H} g_1(p)/p F(theta_1 v-v log p/log x).
```

Aqui el factor variable es `F`, que si es Lipschitz en el rango pertinente. La comparacion de Mertens/partial summation se hace sobre el unico intervalo `[L,H)`; no se necesita que el indicador extendido sea Lipschitz. Los dos puntos extremos aportan a lo sumo `O(1/L)=o(1)`.

### 3.3 Termino Sigma_2

En la deduccion de (5.18), (71.4) reemplaza la afirmacion `W_P^+=W^++O(1/log x)` por igualdad exacta. En la instancia Kuhn, la evaluacion de `U_2` usa la formula (5.24), basada en que `W^+(p_1k)!=0` fuerza `k=1`; no usa (5.22), el otro paso que requiere Lipschitz para un peso general.

Por tanto el argumento Kuhn puede reescribirse sin ninguna constante Lipschitz para `w`: se conserva `theta_2=0.267`, se conserva el mismo funcional integral (los extremos tienen medida cero) y se eliminan los terminos artificiales `C_w/log x`.

## 4. Veredicto y alcance

- **[T-lit]** Hay un hueco textual real: el indicador de Kuhn no satisface la hipotesis Lipschitz del teorema tal como esta impreso.
- **[T]** La particion alineada con `x^(1/6)` hace constantes los pesos en cada bloque y da `W(p_1m)=W_I(m)` exactamente.
- **[P-rev]** Esto repara los usos de Lipschitz en la prueba especifica de Kuhn sin suavizar ni perder masa. Junto con los docs 62 y 66--68, elimina la incompatibilidad Lipschitz/indicador como obstaculo de la auditoria A4.
- No se afirma que el Teorema 3.2 publicado quede formalmente corregido para pesos generales, ni que todas las assumptions Goldbach hayan sido verificadas. Goldbach sigue abierta.

## 5. Verificacion reproducible

```text
python numerico/71_bloques_alineados_kuhn.py
```

Comprueba 2000 escalas que la particion diadica refinada no contiene `1/u` en el interior de ninguna celda y 20000 estados que (71.4) coincide con la definicion directa. Es comprobacion de contabilidad, no evidencia de Goldbach.

## Memoria compartida

Propuesta pendiente: el MCP canonico `lumina_shared_memory` / `lumina_propose` no esta expuesto. No se afirma persistencia remota.
