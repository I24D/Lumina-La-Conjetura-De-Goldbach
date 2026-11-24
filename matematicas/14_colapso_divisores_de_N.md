# 14 — Colapso exacto de la expansión por divisores de N

## Procedencia y estado

- Agente: LUMINA/OpenClaw.
- Fecha: 2026-11-24, America/New_York.
- Estado: **[P]** para las identidades elementales; **[N]** para su contraste.
- Alcance: corrige y simplifica el frente abierto al final de 12. No demuestra
  Goldbach: muestra que la inclusión-exclusión propuesta allí vuelve a aislar
  la correlación binaria original.

## 1. Lema de soporte

Sea `a|rad(N)` y recuerde

\[
G_a(N)=\sum_{\substack{p\nmid N,\ k\ge1\\ap^k<N}}
\log p\,\Lambda(N-ap^k).
\]

Si un primo `ell|a`, entonces `ell|N-ap^k`. Como `Lambda(m)` solo es no nula
cuando `m` es potencia de un único primo, se deduce:

1. si `omega(a)>=2`, entonces `G_a(N)=0` exactamente;
2. si `a=ell` es primo, todo término no nulo satisface
   `N-ell p^k=ell^j`, es decir
   `N/ell=p^k+ell^(j-1)`;
3. por tanto `G_ell(N)=O((log N)^3)` de manera elemental: hay
   `O(log N)` exponentes `j`, y para cada uno queda determinado el posible
   valor de `p^k`; cada peso es a lo sumo `(log N)^2`.

El mismo argumento muestra que `Lambda(N-a)=0` cuando `omega(a)>=2`.

## 2. Fórmula colapsada

Aplicando el lema a (12.8), y usando `mu(ell)=-1`, queda la identidad exacta

\[
\sum_{n<N}\Lambda(n)Q_N(N-n)
=-G_1(N)+\sum_{\ell\mid N\atop \ell\ \mathrm{primo}}
\left(G_\ell(N)-\log\ell\,\Lambda(N-\ell)\right). \tag{14.1}
\]

En particular,

\[
\sum_{n<N}\Lambda(n)Q_N(N-n)=-G_1(N)
+O\!\left(\omega(N)(\log N)^3\right). \tag{14.2}
\]

Así, no queda una inclusión-exclusión grande sobre los divisores impares de
`N`: todos los divisores con al menos dos primos se anulan. El único término
grande es `a=1`, precisamente la convolución de Goldbach (con las exclusiones
de potencias cuyos primos dividen a `N`).

Además, los términos que separan `G_1(N)` de

\[
R(N)=\sum_{m<N}\Lambda(m)\Lambda(N-m)
\]

tienen base prima `ell|N`. Si `m=ell^k`, entonces un sumando no nulo
obliga de nuevo a `N-ell^k=ell^j`. Por el mismo conteo,

\[
G_1(N)=R(N)+O\!\left(\omega(N)(\log N)^3\right),
\]

y (14.2) es, salvo correcciones poli-logarítmicas, `-R(N)`. La expansión
no crea un sustituto más fácil para la correlación: la recupera casi intacta.

## 3. Consecuencia para el programa de RH/positividad

La búsqueda de una desigualdad de Selberg sobre todos los divisores impares
de `rad(N)` era demasiado amplia. La graduación Z/2 de ese conjunto colapsa
antes de producir una nueva desigualdad: no hay masas grandes de grados altos
que comparar, sino un término grande `-G_1` y correcciones poli-logarítmicas.

Esto constituye un **descarte útil** del paso 2 propuesto en 13 §6. Un puente
inspirado en Alpoge tendría que actuar directamente sobre `G_1`, mediante una
estructura aditiva/espectral capaz de conservar el coeficiente individual
`N`; la combinatoria de divisores de `N` no puede aportar por sí sola la
positividad faltante.

## 4. Verificación reproducible [N]

`numerico/colapso_divisores_N.py` evalúa literalmente ambos lados de (14.1),
además de comprobar que todo soporte prohibido con `omega(a)>=2` es cero.
La ejecución registrada usa todos los pares hasta 500 y 300 pares aleatorios
hasta 5000, con semilla fija 20261124.

## 5. Nuevo frente preciso

La próxima ruta no será intentar acotar separadamente las piezas `G_a`. Hay
que buscar una representación espectral o una identidad de cuadrados para
`G_1(N)` en la que:

- el coeficiente individual `N` permanezca visible;
- el término principal sea positivo;
- el error minorante sea estrictamente menor que ese término principal;
- no se presuponga la propia positividad de `G_1(N)`.

La candidata natural es estudiar kernels de Fejér/Riesz aplicados a la suma
exponencial de von Mangoldt, pero cualquier desigualdad deberá ser unilateral
y puntual en `N`, no solo una media sobre `N`.
