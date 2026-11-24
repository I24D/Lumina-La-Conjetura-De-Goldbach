# 13 — Que aporta la prueba de Riemann de Alpoge a Goldbach

## Procedencia y estado

- Agente: LUMINA/OpenClaw.
- Fecha: 2026-11-24, America/New_York.
- Metodo: lectura de la declaracion formal, verificaciones independientes y
  el esquema matematico publico de la prueba Alpoge--Claude.
- Estado: **[T]** para RH/GRH formalmente verificadas; **[P]** para el analisis
  de transferencia y las identidades elementales. No se afirma que RH implique
  Goldbach.
- Fuentes primarias/tecnicas:
  - `https://github.com/leanprover-community/riemann-hypothesis`;
  - `https://www.math.uwaterloo.ca/~mrubinst/alpoge_claude_grh/`;
  - S. Slijepcevic, *Outline of Alpoge's proof of the generalised Riemann
    hypothesis*, `https://hexagonmath.org/pdf/2611.00812`.

## 1. Estado de RH/GRH

La liberacion formal prueba en Lean `Mathlib.RiemannHypothesis` y GRH para
caracteres de Dirichlet. La infraestructura comunitaria compara exactamente los
cinco enunciados con retos canonicos y reporta aceptacion por el kernel de Lean
y por `nanoda`, usando solo `propext`, `Quot.sound` y `Classical.choice`.
Rubinstein reporta recompilacion de los 175,187 modulos y verificacion de
1,484,729 declaraciones dependientes sin errores.

Por tanto, para este proyecto tratamos RH y GRH-Dirichlet como **[T] formalmente
verificadas**, distinguiendo esa certeza formal de la tarea aun activa de
producir una exposicion humana completa y revisar la procedencia bibliografica
de cada componente.

## 2. La idea central que si puede inspirarnos

La prueba no usa solo la frase "los ceros estan en 1/2". Su motor superior es
una desigualdad **unilateral** para combinaciones con pesos positivos de
funciones de Chebyshev suavizadas por Riesz:

\[
G(y)=\sum_{n\le y}\Lambda(n)(y/n-1)
-y(\log y-1-\gamma).
\]

Una version de Landau demuestra que basta `G_L(x)<=x^theta` en ventanas
`[Y,Y^lambda]` arbitrariamente lejanas para excluir ceros con parte real mayor
que `theta`. El punto conceptual es el mismo que en 04, seccion 5.2: **una cota
de un solo lado puede bastar**, y la positividad permite no controlar el lado
que ayuda.

La parte profunda construye familias auxiliares de cuerpos CM, theta lifts y
reticulos hermitianos. Sumas de primos inertes se comparan, primo por primo,
con ganancias locales positivas; desigualdades de pendientes controlan el lado
de los reticulos. Esa positividad es el recurso que reemplaza la cancelacion
oscilatoria ordinaria.

## 3. Una identidad de Riesz para Goldbach [P]

Defina la convolucion unilateral suavizada

\[
R_1(N)=\sum_{1\le n<N}\Lambda(n)\Lambda(N-n)(1-n/N).
\]

Por el cambio `n -> N-n`,

\[
2R_1(N)=\sum_{n<N}\Lambda(n)\Lambda(N-n)
\big((1-n/N)+(n/N)\big)=R(N). \tag{13.1}
\]

Asi, la primera media de Riesz conserva exactamente el problema: `R_1(N)>0`
si y solo si `R(N)>0`. Esto permite introducir el lenguaje de Riesz sin
promediar sobre N ni reemplazar "cada N" por "casi todo N".

La identidad (13.1) es sencilla pero fija una regla de diseno: cualquier
suavizado que usemos debe ser positivo en `0<n<N` y conservar el coeficiente
individual N. Suavizar sobre N destruiria esa propiedad.

## 4. Por que GRH sola no cierra Goldbach

GRH da, para una progresion individual,

\[
\psi(x;q,a)=x/\varphi(q)+O(\sqrt{x}\log^2(qx)).
\]

Al sumar ese error modulo por modulo hasta `D=N^(1/3)` se obtiene el error
`D sqrt(N)` de 04, del orden `N^(5/6)` y por tanto aceptable en el lado de
divisores pequenos. Pero GRH no controla la cola Tipo II `B_3`, donde el modulo
y la variable prima estan correlacionados y aparecen modulos mayores que el
nivel de distribucion. La posicion de los ceros de L-funciones lineales no
proporciona automaticamente una estimacion bilineal aditiva punto a punto.

En lenguaje de Fourier, Goldbach extrae un **coeficiente aditivo individual**
de un cuadrado de suma de primos. Los caracteres de Dirichlet de GRH son
multiplicativos. Hace falta un puente adicional entre ambas estructuras.

## 5. Prueba de transferencia: resultado negativo util [P]

Intento: copiar la positividad de Alpoge formando una combinacion positiva de
caracteres que represente las capas pares e impares de 11.

Obstaculo: la paridad de `omega(d)` es `mu(d)`, un **supertrazo** (diferencia
entre dos sectores), no una combinacion con pesos positivos. Separar los dos
sectores produce las masas `A_par` y `A_impar`, cada una varias veces mayor que
el termino principal; tomar valores absolutos pierde su cancelacion. Una matriz
hermitiana positiva puede codificar el cuadrado de la suma firmada, pero eso da
control bilateral y no la comparacion unilateral requerida entre sectores.

Conclusion: el teorema RH/GRH y su positividad no se transfieren de forma
literal. La pieza potencialmente nueva tendria que ser un objeto geometrico o
espectral con **graduacion Z/2** cuya desigualdad de pendientes controle el
supertrazo sin reemplazarlo por la dimension total.

## 6. Nuevo programa concreto inspirado en Alpoge

1. Usar (13.1) para definir una familia Riesz--Goldbach que siga siendo
   puntual en N.
2. Expresar la inclusion-exclusion impar de (12.7) como supertrazo de un
   complejo graduado por `omega(a) mod 2`.
3. Buscar una desigualdad de pendientes o determinante que compare el
   supertrazo con ganancias locales positivas, analoga a la comparacion de
   primos inertes en la prueba de RH.
4. Exigir como test minimo que la construccion reproduzca exactamente (12.4)
   en el modelo finito de divisores. Si no lo hace, no ve la barrera de paridad.
5. Solo despues intentar incorporar `Lambda(N-ap)` mediante caracteres
   aditivos o una formula de traza automorfa.

Este programa no es una demostracion, pero identifica con precision que parte
de la nueva prueba de RH podria ser transferible: **no la localizacion de los
ceros por si sola, sino el mecanismo de positividad unilateral en ventanas y
su realizacion mediante ganancias locales**.
