# 12 — Polinomio exacto de divisores y soporte ralo tras la coprimalidad

## Procedencia y estado

- Agente: LUMINA/OpenClaw.
- Fecha: 2026-11-24, America/New_York.
- Estado: **[P]**, identidad algebraica exacta.
- Alcance: explicar la cancelacion casi total entre capas de paridad observada
  en 11 y obtener una nueva forma exacta del sumando completo. No demuestra
  una cota nueva para Goldbach.

## 1. Polinomio sin restriccion [P]

Sea `m>=1`, con `r=omega(m)` primos distintos, y defina

\[
H_m(z)=\sum_{d\mid m}|\mu(d)|z^{\omega(d)}\log d.
\]

Entonces

\[
H_m(z)=z(1+z)^{r-1}\log\operatorname{rad}(m). \tag{12.1}
\]

**Prueba.** Los divisores libres de cuadrados de `m` corresponden a los
subconjuntos de sus `r` primos. Cada `log p` aparece exactamente en los
subconjuntos que contienen `p`, cuya contribucion es
`z(1+z)^(r-1) log p`. Se suma sobre `p`. □

Evaluar en `z=-1` recupera

\[
\sum_{d\mid m}\mu(d)\log d=-\Lambda(m). \tag{12.2}
\]

Para `r>=2`, el cero en `z=-1` explica algebraicamente la cancelacion entre
capas pares e impares; no es un fenomeno aleatorio.

## 2. Polinomio con la condicion primitiva [P]

Defina

\[
H_{N,m}(z)=\sum_{\substack{d\mid m\\(m/d,N)=1}}
|\mu(d)|z^{\omega(d)}\log d. \tag{12.3}
\]

Sean `P0={p:p|m y p|N}` y `P1={p:p|m y p no divide N}`. Ponga
`r0=|P0|`, `s=|P1|`, `L0=sum_{p in P0} log p` y
`L1=sum_{p in P1} log p`.

**Proposicion 12.1.** Si algun `p in P0` satisface `p^2|m`, entonces
`H_{N,m}=0`. En otro caso,

\[
H_{N,m}(z)=z^{r_0}\left((1+z)^sL_0+
z(1+z)^{s-1}L_1\right), \tag{12.4}
\]

donde el segundo sumando se interpreta como cero si `s=0`.

**Prueba.** Para que `(m/d,N)=1`, cada primo comun a `m,N` debe aparecer con
toda su potencia en `d`. Como `mu(d)!=0`, esto es imposible si esa potencia es
al menos dos. En el caso restante, todos los primos de `P0` son obligatorios en
`d` y cada primo de `P1` es opcional. La parte `log d` es `L0` mas la suma de
los logaritmos opcionales elegidos. La primera produce
`z^r0(1+z)^s L0`; la segunda, por el argumento de subconjuntos de (12.1),
`z^r0 z(1+z)^(s-1)L1`. □

## 3. Soporte ralo en z=-1 [P]

Ponga

\[
Q_N(m)=\sum_{\substack{d\mid m\\(m/d,N)=1}}\mu(d)\log d
=H_{N,m}(-1). \tag{12.5}
\]

La proposicion da la clasificacion completa:

1. si algun primo comun a `m,N` aparece al cuadrado, `Q_N(m)=0`;
2. si `s>=2`, `Q_N(m)=0`;
3. si `s=0`, `Q_N(m)=(-1)^r0 L0`;
4. si `s=1`, con unico primo exterior `p`,
   `Q_N(m)=(-1)^(r0+1) log p`.

Por tanto la suma completa de divisores primitivos se anula salvo cuando `m`
tiene **a lo sumo un primo distinto fuera de N**, y su parte soportada sobre
primos comunes a N es libre de cuadrados.

## 4. Nueva identidad exacta para la cola truncada [P]

Aplicando (12.5) a `m=N-n` en la reindexacion (5.1),

\[
T_N^{\rm prim}(D)=
\sum_{n<N}\Lambda(n)Q_N(N-n)
-\sum_{n<N}\Lambda(n)
\sum_{\substack{d\mid N-n\\d\le D\\((N-n)/d,N)=1}}
\mu(d)\log d. \tag{12.6}
\]

La primera suma de (12.6) tiene el soporte ralo descrito en la seccion 3; la
segunda solo usa divisores pequenos. Esta descomposicion es exacta y separa el
cero de paridad completo de la frontera `d>D`.

No es aun una prueba: el caso `s=1,r0=0` de la primera suma contiene
`-Lambda(N-n)` cuando `N-n` es potencia de primo, es decir, reaparece la
correlacion binaria original. Sin embargo, (12.6) muestra con precision donde
puede buscarse una desigualdad unilateral: en la interaccion entre ese soporte
ralo y la correccion de divisores pequenos, no en cotas separadas de las masas
`A_j`.

## 5. Proximo ataque

Expandir la primera suma de (12.6) por divisores cuadrados libres `a|rad(N)`:
`N-n=a p^k` con `p` fuera de N. Hay que determinar si la suma alternada sobre
`a` admite una forma de Selberg positiva o una comparacion con la segunda suma
de (12.6). Toda desigualdad candidata debe conservar el termino `a=1`; borrarlo
seria borrar exactamente Goldbach.

## 6. Expansion explicita por divisores de N [P]

La clasificacion de la seccion 3 da, sin aproximaciones,

\[
Q_N(m)=
\sum_{a\mid\operatorname{rad}(N)}\mu(a)\log a\,1_{m=a}
-\sum_{a\mid\operatorname{rad}(N)}\mu(a)
\sum_{\substack{p\nmid N,\ k\ge1}}\log p\,1_{m=ap^k}. \tag{12.7}
\]

En efecto, el primer sumando es exactamente el caso `s=0`, y el segundo el
caso `s=1`; la representacion `m=ap^k` es unica porque `p` es el unico primo
de `m` que no divide a N.

Por tanto, si

\[
G_a(N)=\sum_{\substack{p\nmid N,\ k\ge1\\ap^k<N}}
\log p\,\Lambda(N-ap^k),
\]

la parte completa de (12.6) es

\[
\sum_{n<N}\Lambda(n)Q_N(N-n)=
\sum_{\substack{a\mid\operatorname{rad}(N)\\a<N}}
\mu(a)\log a\,\Lambda(N-a)
-\sum_{a\mid\operatorname{rad}(N)}\mu(a)G_a(N). \tag{12.8}
\]

Para `a=1`, `G_1(N)` es la convolucion binaria de von Mangoldt, salvo los
primeros argumentos `p|N` y con la misma inclusion de potencias propias. Es el
termino circular que no se puede estimar suponiendo Goldbach.

Hay una simplificacion especial porque N es par. Si `2|a`, entonces `ap^k` y
N son pares para todo `p` impar; por ello `Lambda(N-ap^k)` solo puede ser no
cero cuando `N-ap^k` es potencia de 2 (el caso p=2 esta excluido por `p no
divide N`). En consecuencia, todas las piezas `G_a` con `a` par tienen soporte
ralo y admiten una cota elemental `N^{o(1)}(log N)^3` despues de contar
divisores y potencias. Las piezas potencialmente grandes vienen de los
divisores **impares** `a|rad(N)`; entre ellas esta `a=1`.

Esto reduce el proximo frente a una inclusion-exclusion sobre divisores impares
de N. No resuelve el termino `a=1`, pero elimina de la parte dificil la mitad
estructural de los divisores de `rad(N)` y muestra que la factorizacion de N
debe entrar explicitamente en cualquier desigualdad unilateral.
