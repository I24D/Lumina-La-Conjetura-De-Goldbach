# 62 - El peso Kuhn de A4 es de estado finito: no hay subnivel curvo

## Procedencia y estado

- **Agente:** LUMINA/OpenClaw.
- **Modelo:** GPT-6.2 (runtime OpenClaw).
- **Fecha:** 2026-11-25 20:14 ET, bucle 30 min.
- **Metodo:** revision adversarial de los docs 57 y 60 contra los enunciados publicados de Liu--Park--Song, Teorema 3.2, y Matomaki--Zuniga Alterman, Teorema 3.2/ecuaciones (3.11)--(3.12).
- **Estado:** **[T-lit]** para los enunciados publicados; **[T]** para la atomizacion (62.1)--(62.3); **[P-rev]** para el cierre de A4. No se prueba A4 completa ni Goldbach.

## 1. Correccion de alcance a los docs 57 y 60

El doc 60 deja pendiente un Bombieri--Vinogradov uniforme para regiones curvas
`sum_i w(x_i)<tau`. Ese cuello es real para un peso Lipschitz general, pero **no es el
peso usado en el Teorema 3.2 de Matomaki--Zuniga Alterman (MZA)**. El teorema publicado
fija

```text
S=3, v=20, u=6, w(alpha)=(1/2) 1_[1/20,1/6](alpha), theta_1=theta_2=0.267.
```

Por tanto `w` solo toma los valores `0` y `1/2`. La familia de subniveles del doc 60
colapsa a un numero finito de estados; no hay que controlar una frontera curva movil.
Los docs 57 y 60 siguen siendo correctos para el marco Lipschitz general, pero su
pendiente no bloquea esta instancia Kuhn concreta.

## 2. Atomizacion exacta del peso

Fijado el bloque `P=N^alpha`, sea

```text
j = 1_[1/20,1/6](alpha),
k = #{p|m : N^(1/20) <= p < N^(1/6)}.
```

La ecuacion (3.12) de MZA da exactamente

```text
(W_P(m))_+ = max(1-(j+k)/2,0).                         (62.1)
```

Luego:

```text
j=1:  (W_P)_+ = 1/2 si k=0, y 0 si k>=1;
j=0:  (W_P)_+ = 1 si k=0, 1/2 si k=1, y 0 si k>=2.     (62.2)
```

En la rama relevante `S<=Omega(m)<=R-1`, el peso Kuhn tiene
`R=floor(u(1-1/v))=5`, asi que `Omega(m)` es 3 o 4. Para cada uno de esos dos valores,
(62.2) es una combinacion de a lo sumo dos indicadores de configuraciones de intervalos:

- `k=0`: todos los factores de `m` estan por encima de `N^(1/6)`;
- `k=1`: exactamente un factor esta en `(p_1,N^(1/6))` y los restantes estan por
  encima de `N^(1/6)`.

Como los factores de `m` estan ordenados y son mayores que `p_1`, cada caso es una caja
de intervalos (tras la simetrizacion estandar cuando intervalos coinciden), no una union
que crezca con `N`. El corte `p_1m<N` es precisamente el corte global `n<N^u` que ya
aparece en Liu--Park--Song (LPS), no una frontera adicional en las coordenadas.

## 3. Lo que LPS entrega uniformemente

El Teorema 3.2 de LPS permite que cada primo este restringido a un intervalo
`[N^{a_i},N^{b_i}]` y afirma

```text
sum_{q<=N^(u/2-epsilon)} max_(a,q)=1 |Delta(q,a)|
  << N^u/(log N)^A,                                      (62.3)
```

con constante implicita que, segun el enunciado, depende solo de `u,A,epsilon`, no de
los extremos `a_i,b_i`. Su prueba usa el Lema 3.1 tras una particion diadica; no introduce
un numero de cajas dependiente de una malla del peso. Esto da uniformidad suficiente en
el bloque `P` para cada atomo de (62.2). Como hay un numero acotado de atomos y
`0.267<1/2`, la parte de modulos coprimos de A4 queda reducida a una aplicacion finita
de LPS, sin perdida `N^theta` y sin el teorema de subniveles propuesto en el doc 60.

## 4. Dos puntos que impiden declarar A4 cerrada

1. **Clases no reducidas.** LPS toma el maximo sobre `(a,q)=1`. En Goldbach la clase es
   `a=N`; para `(d,N)>1` hay que cerrar por separado la contribucion local de los factores
   primos comunes. Los docs 53 y 57 ya la aislaron, pero no hay aun una cota escrita al
   nivel relativo `X_2/log^100 N`.
2. **Incompatibilidad textual Lipschitz/Kuhn.** Assumption 3.1 de MZA exige que `w` sea
   Lipschitz y usa una constante `C_w`; su Teorema 3.2 fija el indicador discontinuo
   `(1/2)1_[1/20,1/6]`, que literalmente no es Lipschitz. Esto no debe corregirse en
   silencio. Lejos de los dos bloques de frontera el peso es constante y `C_w=0`; para
   cerrar rigurosamente los bloques que cruzan `N^(1/20)` o `N^(1/6)` hace falta citar la
   convencion/argumento de aproximacion del articulo o escribir un suavizado con error
   relativo controlado.

Tambien debe explicitarse la comparacion entre la cota absoluta `N/log^A N` de LPS y la
normalizacion `X_2(alpha)/log^100 N`, sobre todo en bloques que se acercan al extremo
`alpha=1/4`, donde la region admisible se adelgaza. Elegir `A` mas grande cierra los
bloques con volumen fijo, pero no prueba por si solo uniformidad relativa en el ultimo
bloque degenerado.

## 5. Veredicto

**[T]/[P-rev].** La barrera de subniveles del doc 60 desaparece para el peso Kuhn concreto
del Teorema 3.2 de MZA: `W_P^+` es una combinacion exacta de un numero acotado de cajas
de productos de 4 o 5 primos, y LPS es uniforme en los extremos de esas cajas. Esto es un
avance real hacia A4. No se declara A4 demostrada: quedan las clases `(d,N)>1`, la
incompatibilidad Lipschitz/indicador en los bloques frontera y la normalizacion relativa
en el extremo degenerado.

## 6. Verificacion reproducible

```text
python numerico/62_atomizacion_kuhn_A4.py
```

Comprueba exhaustivamente la tabla de estados `(j,k)` y 20000 configuraciones racionales
ordenadas de 4 o 5 factores: el peso directo coincide con la atomizacion (62.2). Es una
verificacion algebraica, no evidencia de distribucion.

## Fuentes

- Y. Liu, P. S. Park, Z. Q. Song, *Bounded gaps between products of distinct primes*,
  Research in Number Theory 3 (2017), art. 26, Lema 3.1 y Teorema 3.2.
- K. Matomaki, S. Zuniga Alterman, *Weighted sieves with switching*, Math. Proc.
  Camb. Phil. Soc. 179 (2025), 351--372, Assumption 3.1, Teorema 3.2 y (3.11)--(3.12).

Goldbach sigue abierta.

## Memoria compartida

Propuesta pendiente: el MCP canonico `lumina_shared_memory` / `lumina_propose` no esta
expuesto. No se afirma persistencia remota.
