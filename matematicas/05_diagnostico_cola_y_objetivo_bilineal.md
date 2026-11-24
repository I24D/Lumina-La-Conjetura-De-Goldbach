# 05 — Diagnostico exacto de la cola y nuevo objetivo bilineal

## Procedencia y estado

- Agente: LUMINA/OpenClaw.
- Fecha: 2026-11-23, America/New_York.
- Metodo: reindexacion algebraica exacta y calculo finito reproducible.
- Estado: **[P]** para las identidades y reducciones escritas; **[N]** para los datos. No es una demostracion de Goldbach.
- Dependencia: `04_reduccion_mobius_goldbach.md`.

## 1. Reindexacion exacta que abarata el calculo [P]

Para `D=floor(N^(1/3))`, la cola primitiva de 04 puede escribirse exactamente como

\[
T_N^{\rm prim}(D)=
\sum_{1\le n<N}\Lambda(n)
\sum_{\substack{d\mid N-n\\ d>D\\ ((N-n)/d,N)=1}}
\mu(d)\log d. \tag{5.1}
\]

**Prueba.** En la definicion de 04 se hace `n=N-cd`. Entonces `d|(N-n)` y
`c=(N-n)/d`; la condicion `(c,N)=1` se convierte en la ultima condicion de
(5.1). La transformacion es biyectiva y conserva el peso. □

Esta forma permite cribar `mu` y `Lambda` una vez y enumerar solamente los
divisores de `N-n` cuando `Lambda(n)` no se anula. El programa
`numerico/cola_mobius_scan.py` implementa (5.1). Una segunda implementacion
literal con los bucles `(c,d)` coincidio para todos los pares `N<=300`.

## 2. Primer diagnostico de la hipotesis HmuLambda [N]

Se calculo literalmente

\[
U_N(D)=\sum_{c\le (N-1)/(D+1),(c,N)=1}
\sup_{D\le t\le (N-1)/c}
\left|\sum_{D<d\le t}\mu(d)\Lambda(N-cd)\right|.
\]

| N | D | U/N | U log(N)/N |
|---:|---:|---:|---:|
| 100 | 4 | 0.48935 | 2.25353 |
| 1,000 | 10 | 0.46857 | 3.23674 |
| 10,000 | 21 | 0.37438 | 3.44820 |
| 100,000 | 46 | 0.28448 | 3.27516 |
| 1,000,000 | 100 | 0.22938 | 3.16907 |

HmuLambda de 04 exige `U log(N)/N -> 0`. En este rango la cantidad permanece
cerca de 3.2 desde `10^3` hasta `10^6`. Esto **no refuta** un limite
asintotico, pero si es evidencia adversa fuerte contra usar HmuLambda como
objetivo principal. El supremo separado para cada `c` parece pagar un sesgo
positivo comparable con `N/log N`, aunque la cola firmada global puede cancelar.

Conclusion operativa: HmuLambda sigue siendo suficiente, pero queda
**despriorizada** como blanco posiblemente falso o demasiado fuerte.

## 3. La cola firmada si muestra cancelacion, pero no monotona [N]

Se evaluo (5.1) en 84 pares espaciados logaritmicamente hasta `10^6`. En la
muestra, el maximo positivo de `T/(N*S(N))` fue `0.29728` en `N=129994`; el
minimo fue `-0.85763` en `N=2940`. En `N=10^6`,

\[
T_N^{\rm prim}/N=0.0642495,\qquad
T_N^{\rm prim}/(N\mathfrak S(N))=0.0364964.
\]

Los valores negativos grandes no perjudican Goldbach. Los datos apoyan estudiar
directamente el lado positivo de la suma firmada, no la suma de supremos
absolutos. No establecen ningun `o(N)`.

## 4. Objetivo bilineal diadico que conserva la cancelacion [P]

Para intervalos diadicos `C<c<=2C`, `K<d<=2K`, defina

\[
B_N(C,K)=\sum_{\substack{C<c\le2C\\(c,N)=1}}
\sum_{\substack{K<d\le2K\\cd<N}}
\mu(d)\log d\,\Lambda(N-cd). \tag{5.2}
\]

La region `d>D, cd<N` se cubre por `O((log N)^2)` cajas diadicas. Por tanto,
para cualquier `A>3`, la estimacion uniforme

\[
B_N(C,K)\ll_A \frac{N}{(\log N)^A} \tag{HB}
\]

en todas las cajas admisibles implica

\[
T_N^{\rm prim}(D)\ll_A \frac{N}{(\log N)^{A-2}}=o(N).
\]

La prueba es particionar ambas variables en potencias de dos y aplicar la
desigualdad triangular **despues** de sumar cada caja. A diferencia de
HmuLambda, (HB) conserva la cancelacion conjunta entre valores de `c` y `d`.
No se afirma que (HB) este demostrada; es el nuevo blanco tecnico.

Para positividad basta una version unilateral agregada:

\[
\sum_{C,K} B_N(C,K)\le (1-\kappa)N\mathfrak S(N). \tag{HB+}
\]

## 5. Siguiente paso matematico

Descomponer `Lambda(N-cd)` en (5.2) mediante una identidad bilineal exacta
(Vaughan o Heath-Brown), manteniendo visibles los rangos `C,K`. Cada termino
debe clasificarse como Tipo I, Tipo II o diagonal, indicando con precision cual
estimacion publicada controla cada rango. No se aceptara "cancelacion esperada"
ni una aplicacion de Bombieri-Vinogradov fuera de su nivel.

## 6. Reproducibilidad y limites

- Script: `python numerico/cola_mobius_scan.py --hasta 1000000 --muestras 40 --u`.
- Resultado: `numerico/resultados/20261124_023525_cola_mobius_1000000.json`.
- Validaciones: `py_compile` PASS; doble bucle PASS para los 149 pares `4<=N<=300`.
- La serie singular se trunca en primos hasta 200000 solo en cocientes
  diagnosticos; `T` y `U` son exactos salvo redondeo de logaritmos IEEE-754.
- Ningun resultado numerico de esta nota prueba un enunciado asintotico.
