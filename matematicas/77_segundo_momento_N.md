# 77 - Segundo momento sobre N (74.4): SE TIENE; el puente es el problema

Fecha: 2026-11-26. Script: `numerico/77_segundo_momento_N.py` — PASS
58/58 local (FFT, X = 2x10^5, 5x10^5, 10^6). No hay cota inferior de
capas, no hay deteccion de P_2, no hay demostracion de Goldbach.

## 0. (74.4), formulado con precision

Escala fija X. Sean S_1(alfa) = sum_{p<2X} e(p*alfa),
S_rho(alfa) = sum_{n<2X} rho(n) e(n*alfa) con rho el minorante de
Harman a escala fija 2X (doc 66; el proxy rho_tilde del doc 67 en la
medicion), y M, m los arcos mayores/menores a escala fija:
M = union sobre q <= Q_0 = (log X)^B de [a/q - 1/(qQ_1),
a/q + 1/(qQ_1)], Q_1 = X/(log X)^B, B fijo grande; m = [0,1] \ M.
Para N par en [X, 2X]:

  E_U(N) := int_m S_1(alfa) S_rho(alfa) e(-N alfa) d alfa.

(Con S_1, S_rho y m fijos a escala X, la integral sobre [0,1]
reproduce U(N) = sum_{p<N} rho(N-p) exactamente para N <= 2X,
pues p < N y N-p < N fuerzan p, N-p < 2X. La variante del doc 70
con arcos dependientes de N difiere solo en factores logaritmicos;
el enunciado (74.4) es robusto a esa eleccion.)

**Enunciado (74.4) (doc 74, sec. 3).**

  sum_{N en [X,2X], N par} |E_U(N)|^2 = o(X^3/log^4 X).

La escala X^3/log^4 X es la escala natural: sum_{N} (termino
principal)^2 ~ c * X^3/log^4 X (medido en el doc 76 por FFT:
sum U^2 = 4.49x, sum R^2 = 6.65x X^3/log^4 X). La pregunta del
doc 76 era si (74.4) se sigue del segundo momento clasico de R(N)
mas un segundo momento para D(N). Respuesta: SI, y por la via
directa es aun mas simple.

## 1. Prop. 77.1 [T]: (74.4) SE TIENE, por Parseval + Vinogradov

**Prop. 77.1.** Para todo A fijo, con B = B(A) grande,

  sum_{N en [X,2X]} |E_U(N)|^2 <<_A X^3/(log X)^{2A+1},

en particular (74.4). Lo mismo vale para E_R y para E_D por
separado.

*Prueba.* Sea f(alfa) = 1_m(alfa) S_1(alfa) S_rho(alfa). Entonces
E_U(N) = f^(N) (coeficiente de Fourier), y por Parseval/Bessel,

  sum_{N en [X,2X]} |E_U(N)|^2 <= sum_{N en Z} |f^(N)|^2
      = int_m |S_1|^2 |S_rho|^2 d alfa
      <= (sup_m |S_1|^2) * int_0^1 |S_rho|^2 d alfa
      = (sup_m |S_1|^2) * sum_n rho(n)^2.

Los dos factores: (i) sup_m |S_1| <<_A X/(log X)^A (Vinogradov;
Prop. 4.1 del survey, para la suma pesada por Lambda, mas
sumacion parcial para pasar a la suma sin peso — es la cota que
el doc 69 Prop. 69.2 usa como [T-lit]); (ii) sum rho(n)^2 <<
X/log X: para el proxy, sum rho_tilde^2 = pi(2X) + sum D^2 =
1.467*pi(2X) (docs 67/70, medido; en el script de este doc:
1.404-1.467 pi(S) segun la escala). Luego la suma es
<< X^2/(log X)^{2A} * X/log X = X^3/(log X)^{2A+1}. Con A = 2 ya
es X^3/(log X)^5 = o(X^3/log^4 X). QED.

La via por descomposicion (la que el doc 76 pedia comprobar) es
la misma prueba aplicada a cada pieza, usando E_U = E_R - E_D
(doc 70, Prop. 70.1): sum |E_R|^2 <= sup^2 * pi(2X) y
sum |E_D|^2 <= sup^2 * sum D^2, con sum D^2 = 0.467*pi(2X)
(doc 70; en el script: 0.404-0.467 pi(S)). El defecto es delgado
en L^2, asi que su segundo momento tiene la misma forma que el
clasico. Luego sum |E_U|^2 <= 2 sum |E_R|^2 + 2 sum |E_D|^2
tambien da (74.4).

*El analogo clasico (leido).* Es exactamente el Lema 4.2 del
survey Bhowmik-Grimmelt: para la suma pesada S(alfa) =
sum Lambda(n) e(alfa n), sum_{N<=X} |r_m(N)|^2 <= int_m |S|^4
<= sup_m |S|^2 * int |S|^2 << (X^2/R + X^{8/5}) X (log X)^5,
y con R = (log X)^A da sum |r_m|^2 << X^3 (log X)^{-A+5}. De
ahi el Teorema 4.3 (Chudakov; van der Corput; Estermann):
|E(X)| <<_A X (log X)^{-A} para todo A. Es decir: (74.4) es la
maquinaria clasica del casi-todo-N aplicada a U(N). NO es un
teorema nuevo: es aplicar el argumento de 1937-1975 al conteo
unilateral. Su valor aqui es cerrar la pregunta: el segundo
momento ESTA disponible; lo que falta es solo el puente.

*Version condicional (leida, survey sec. 3, (3.8)).* Por formula
explicita, sum_{N<=X} |E_2(N)|^2 << X^{3/2+2Theta+epsilon},
donde E_2(N) = r_2(N) - N S(N) es el error pesado y Theta el
supremo de las partes reales de los ceros. Bajo GRH (Theta =
1/2): << X^{5/2+epsilon}, un ahorro de potencia en el segundo
momento. Solo bate la escala natural X^3 cuando Theta < 3/4.
Bajo cuasi-RH Theta < 11/12 (Li-Liu) no da ahorro en el segundo
momento; el resultado E_+ << N^{1/2} (log N)^{100} de Li-Liu
usa el cuarto momento, no el segundo.

## 2. Medicion: el segundo momento esta dominado por el termino secundario conocido

Conteos exactos R(N), U(N), D(N) = R(N) - U(N) por convolucion
FFT (indicador de primos contra si mismo y contra rho_tilde, a
escala fija 2X), terminos principales M_R = S(N) N/log^2 N,
M_U = (1-L_S) M_R, M_D = L_S M_R con L_S medido por escala.
Errores e_R = R - M_R, e_U = U - M_U, e_D = D - M_D
(e_U = e_R - e_D exacto, verificado a 1e-12). Cocientes respecto
a la escala X^3/log^4 X:

| X       | sum e_R^2 | sum e_U^2 | sum e_U^2 (c_rho) | sum e_D^2 | malos R/U/D | a_sec | L_S   |
|---------|-----------|-----------|-------------------|-----------|-------------|-------|-------|
| 2x10^5  | 0.168     | 0.119     | 0.767             | 0.015     | 0 / 2 / 10  | 2.40  | 0.285 |
| 5x10^5  | 0.145     | 0.099     | 0.637             | 0.014     | 0 / 1 / 2   | 2.38  | 0.308 |
| 10^6    | 0.130     | 0.086     | 0.551             | 0.012     | 0 / 0 / 0   | 2.36  | 0.325 |

("malos": N con |error| >= (1/2) termino principal; "a_sec":
coeficiente del ajuste e_R ~ a_sec * S(N) N/log^3 N.)

Lectura:
- Los cocientes decrecen con X en los tres conteos, y estan muy
  por debajo de 1: (74.4) se verifica en escala a X accesible.
- El cociente esta DOMINADO por el termino secundario conocido:
  el ajuste de dos terminos e_R ~ a*S(N)N/log^3 N +
  b*S(N)N/log^4 N da (aR, bR) = (2.25, 1.93), (2.07, 4.16),
  (2.05, 4.44) a X = 2x10^5, 5x10^5, 10^6: aR -> 2.05, bR ->
  4.44, los coeficientes de la expansion de li (1, 2, 6... en la
  forma integral; el 4.44 medido es el coeficiente efectivo del
  termino 1/log^4). La prediccion del cociente desde el termino
  secundario, (14/3)*media(S^2)/log^2 X * (aR/2)^2 con
  media(S^2) = 4.60 medida, es 0.182, 0.134, 0.118 frente a los
  medidos 0.168, 0.145, 0.130: el termino secundario explica el
  segundo momento dentro del ~10%.
- Tras quitar los dos terminos sistematicos, el residuo
  fluctuante es minusculo: sum residuo_R^2 / escala = 0.0006,
  0.0003, 0.0002; residuo_U = 0.0084, 0.0072, 0.0052; residuo_D
  = 0.0079, 0.0070, 0.0051. (El ajuste es en la misma muestra,
  asi que el residuo es una cota inferior del error fluctuante,
  no una estimacion insesgada; aun asi muestra que el grueso del
  segundo momento es sistematico conocido, no error desconocido.)
- Con el c_rho asintotico 0.465 en vez del L_S finito, el cociente
  de U es mayor (0.77 -> 0.55) y decrece: es la deriva lenta de
  L_S(N) -> 0.535 (doc 67 sec. 2), no un error nuevo.
- Los malos tienden a 0 (2 -> 1 -> 0 en U; 10 -> 2 -> 0 en D):
  a X = 10^6 ningun N par en [10^6, 2x10^6] tiene error mayor que
  la mitad del termino principal en ninguno de los tres conteos.

Consecuencia (Chebyshev, con la Prop. 77.1): el numero de N
pares en [X, 2X] con |E_U(N)| >= (c_rho/2) S(N) N/log^2 N es
<<_A X/(log X)^{2A-3} para todo A; para esos N, U(N) > 0 (el
termino principal domina) y por tanto R(N) > 0. Es el teorema de
casi-todo-N para el conteo unilateral, con la misma calidad
logaritmica que el Teorema 4.3 del survey para R(N).

## 3. El puente casi-todo -> todo: cuantificacion y auditoria

**Cuantificacion del hueco.** De (74.4) (forma fuerte de la
Prop. 77.1) se obtiene un conjunto excepcional de tamano
<<_A X/(log X)^A para todo A. El mejor conocido incondicional
es mucho mas fuerte: Pintz 2018, E(X) << X^{0.72} (delta = 0.28,
survey sec. 6); Zhao 2025, E(X) << X^{0.70} (GB-L05,
arXiv:2511.05631). Bajo GRH como teorema: E(X) << X^{1/2+eps}
(Hardy-Littlewood; Goldston 1992 reemplaza el eps por un
logaritmo). El hueco que queda es, por tanto: de X^{0.70} (o de
X/(log X)^A, si solo se usa el segundo momento) a 0. No es un
hueco cuantitativo pequeno: es cualitativo, y ningun mecanismo
conocido lo cierra. Auditoria:

(i) **No hay continuidad en N.** U(N) es un conteo entero; saber
U(N') > 0 para N' proximo a N no dice nada de U(N): los pares
de primos que suman N' y los que suman N son disjuntos. No
existe ninguna relacion de transferencia entre U(N) y U(N').
Un puente necesitaria una desigualdad que ligue el error en N
con errores en N' vecinos, y no existe ninguna.

(ii) **La barrera de paridad bloquea la criba.** La criba (con
switching, Chen) da N = p + P_2 para todo N, no N = p + p: no
distingue un primo de un producto de dos primos. El unico metodo
que podria dar todo-N es el metodo del circulo, y sus arcos
menores son exactamente la obstruccion (Prop. 76.5).

(iii) **GRH no cierra el puente.** Bajo GRH el conjunto
excepcional es X^{1/2+eps} y el segundo momento es X^{5/2+eps}
(sec. 1), pero el error puntual de arcos menores sigue sin
control: el sup ideal N^{1/2} deja (log N)^{3/2} (doc 69,
Prop. 69.2b). Ni siquiera con GRH como teorema hay puente.

(iv) **Intervalos cortos siguen siendo casi-todo.** Grimmelt
2022: todos salvo O_{A,eps}(H (log X)^{-A}) pares en [X-H, X]
son Goldbach cuando X^{7/120+eps} < H <= X (via Li-Liu 2026,
sec. 1). Es casi-todo en intervalos cortos, no todo: el mismo
hueco, mas pequeno.

(v) **El puente es donde viven los ceros.** Bhowmik-Ruzsa,
completando un argumento de Granville: un ahorro de potencia en
el error PROMEDIO implica una franja libre de ceros fija para
zeta (via Li-Liu 2026, sec. 1). Es decir, mejorar el promedio
mas alla de cierto punto ES informacion dura sobre los ceros;
el paso de promedio a puntual no es un tecnicismo, es el
problema. En la direccion contraria, el ahorro de potencia del
conjunto excepcional (Montgomery-Vaughan 1975, Pintz 2018) viene
de estimaciones de densidad de ceros (Gallagher), no del
segundo momento: el segundo momento solo da el ahorro
logaritmico (Teorema 4.3); la potencia exige input de ceros.

**Prop. 77.2 [T].** No existe ningun puente conocido de
casi-todo-N a todo-N para el conteo unilateral (ni para R(N)):
cualquier puente necesitaria el control puntual
|E_U(N)| = o(N/log^2 N) para todo N, que por la Prop. 76.5 es
equivalente, modulo la asintotica del defecto, a los arcos
menores de Goldbach. El puente ES el problema, no un paso
hacia el.

## 4. Veredicto y que sigue

**Veredicto: (74.4) SE TIENE.** Se sigue del segundo momento
clasico de R(N) (Lema 4.2 del survey: Parseval + Vinogradov)
mas el segundo momento de D(N) (misma prueba; el defecto es
delgado en L^2, sum D^2 = 0.467 pi), y por la via directa es
aun mas simple (Prop. 77.1). La medicion confirma la escala:
los segundos momentos de los errores son o(X^3/log^4 X),
dominados por el termino secundario conocido (aR -> 2.05,
bR -> 4.44), con residuo fluctuante ~0.005 y malos -> 0.
Honestidad: esto es maquinaria clasica (Chudakov, van der
Corput, Estermann; Montgomery-Vaughan), no un descubrimiento
nuevo — es aplicar el argumento del casi-todo-N al conteo
unilateral. Lo que cierra es la pregunta: el segundo momento
esta disponible; lo que falta es SOLO el puente, y el puente
es el problema entero (Prop. 77.2, Prop. 76.5).

Que sigue: el frente vivo para Goldbach sigue siendo el de
seleccion del bucle (switching, C10/W(d); Motohashi doc 57
sec. 4; b_{P,l} doc 53; uniones C_k+C_{k+1}; variante Lema 4
Pintz-Althofer). Una idea nueva para el puente tendria que
producir control puntual de arcos menores sin pasar por el
sup, el cuarto momento ni la correlacion de pares — las tres
vias estan cerradas (docs 67, 69, 70, 72, 76). Goldbach sigue
abierta.

Referencias (leidas, no de memoria):
- G. Bhowmik y L. Grimmelt, "The exceptional set of the Goldbach
  problem", arXiv:2607.27282v2; Anal. Math. 52(3) (2026),
  767-808. Leido: sec. 1 (GRH: E << X^{1/2+eps}; MV 1975;
  Pintz delta = 0.28), sec. 3 ((3.8): sum |E_2|^2 <<
  X^{3/2+2Theta+eps}), sec. 4 (Prop. 4.1 Vinogradov-Vaughan;
  Lema 4.2; Teorema 4.3 Chudakov/van der Corput/Estermann),
  sec. 5 (Teorema 5.1 Montgomery-Vaughan; Prop. 5.2 Gallagher).
- Li-Liu, "Fourth moments of prime exponential sums and
  exceptional sets", arXiv:2610.26756v1, sec. 1 (leido en el
  texto): E_+ << N^{1/2} (log N)^{100} bajo cuasi-RH
  Theta < 11/12; Grimmelt 2022 (intervalos cortos); Bhowmik-
  Ruzsa/Granville (ahorro en media => franja libre de ceros);
  clasico E_+ << N^{2Theta-1/2+eps} para Theta < 3/4.
- Zhao, arXiv:2511.05631 (via GB-L05): E(X) << X^{7/10}.
- Docs 66, 67, 69, 70, 72, 74, 75, 76 de esta carpeta.

## Alcance honesto

- Lo que este doc SI dice: (74.4) formulado con precision y
  PROBADO por el argumento clasico (Prop. 77.1); la medicion
  de la escala (el segundo momento esta dominado por el termino
  secundario conocido; residuo fluctuante ~0.005; malos -> 0);
  la cuantificacion del puente (de X/(log X)^A, o X^{0.70} con
  lo mejor conocido, a 0) y la auditoria de que ningun puente
  conocido lo cierra (Prop. 77.2).
- Lo que NO dice: no hay control puntual de arcos menores, no
  hay cota inferior de capas, no hay deteccion de P_2, no hay
  demostracion de Goldbach. (74.4) es casi-todo-N, no todo-N.
