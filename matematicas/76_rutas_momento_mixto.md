# 76 - Mapa de rutas restantes al momento mixto (sin correlacion de pares)

Fecha: 2026-11-26. Script: `numerico/76_rutas_momento_mixto.py` — PASS
27/27 local (X = 2x10^6; FFT X0 = 10^5). No hay prueba del lema del
momento mixto, no hay cota inferior de capas, no hay demostracion de
Goldbach.

## 0. El lema y lo cerrado

El lema que falta (docs 67, 69) es MM-69.1: int_m |S_1 S_rho| =
o(N/log^2 N), con S_1(alfa) = sum_{p<N} e(p*alfa), S_rho(alfa) =
sum_n rho(n) e(n*alfa) (rho el minorante de Harman, doc 66, diadico a
escala N) y m los arcos menores. La version con fase,
E_U = int_m S_1 S_rho e(-N alfa) = o(N/log^2 N), es lo que el metodo
del circulo necesita para el conteo unilateral U(N) = sum_{p<N}
rho(N-p) (doc 70). El objetivo de este doc es mapear TODAS las rutas
restantes al lema que NO pasen por correlacion de pares, dar el
veredicto de cada una con su razon estructural medida, y elegir la
mas prometedora con plan concreto.

Cerrado hasta aqui (recap, con la referencia exacta):
- L^2 / Cauchy (doc 67, ruta (ii)): (pi(N)*sum rho^2)^{1/2} = 180,375
  = 19.0x el objetivo N/log^2 N = 9,501 a X = 2x10^6; un log por
  encima, asintoticamente ~2*log N. Las normas L^2 estan fijadas por
  la diagonal.
- sup x L^1 (doc 69, Prop. 69.2): las cotas disponibles estan N^{1/2}
  por encima de lo requerido, y aun el sup ideal N^{1/2} deja
  (log N)^{3/2} (145x a N = 10^12).
- Cuarto momento (doc 67(iii), doc 69, Prop. 69.3): Li-Liu da
  X^{5/2} frente a X^2/log^4 necesario; la diagonal de las piezas no
  negativas del defecto es 17.3x y 95.8x lo necesario.
- Descomposicion de Vaughan/Heath-Brown por piezas (doc 69, Props.
  69.1-69.3): descomponer solo redistribuye la masa L^2, nunca la
  reduce; las tres estimaciones por piezas fallan.
- Via conjunta con fase (doc 70, Prop. 70.1): circularidad,
  E_U = E_R - E_D, donde E_R es el error de arcos menores de Goldbach
  mismo; E_D es de tamano comparable (D(N)/R(N) = 0.19-0.35, Cauchy
  0.683*N/log N, un log por encima).
- Gran criba del producto (doc 70, sec. 4a): reproduce Cauchy
  (sum U^2 = 0.943*X^3/log^4 X).
- GRH (doc 70, sec. 4b): no rescata ninguna ruta (ni el sup ideal:
  145x; cuarto momento ideal (log N)^3 por encima).
- Correlacion de pares (doc 72): SL-70.1 NO SE TIENE — la meseta de
  Montgomery (|alfa| >= 1) en version uniforme Dirichlet es conjetura
  abierta incluso bajo GRH; su forma Dirichlet equivale a Hooley suave
  (abierta, con extension falsa a q pequeno); la uniformidad en a = N
  es otra capa abierta.
- Correccion adversarial (doc 74, bucle): [FALSO] la implicacion tal
  como se enuncia en SL-70.1. La correlacion de pares controla una
  estadistica cuadratica promediada (una varianza), no el producto
  absoluto MM-69.1 ni la fase de un coeficiente puntual E_U(N): una
  traslacion ciclica preserva todo el espectro de potencia y cambia un
  coeficiente puntual. El sustituto honesto es el segundo momento
  sobre N, sum_{N en [X,2X]} |E_U(N)|^2 = o(X^3/log^4 X) (74.4),
  aun no probado, mas un puente de casi-todo-N a todo-N que no existe.
  La conclusion cualitativa del doc 72 sigue en pie (la ruta no esta
  disponible); lo que se retira es la afirmacion de suficiencia
  precisa.
- SL-65.1 cerrado (doc 75): el umbral 1/10 de Maynard III es agudo
  para el metodo de la sec. 11 (desigualdad ligante
  Q^(15/2)R^10 <= x^4 del Lema 11.8; el hueco es una potencia
  x^gamma con gamma >= 0.0375, no un logaritmo).

## 1. Ruta A: densidad de ceros directa como input

Pregunta: ¿que estimacion N(sigma,T) bastaria para el momento mixto?
¿Existe en la literatura o es tambien abierta?

**Prop. 76.1 [T]. La via sup necesitaria Theta < 1/2, y eso es
REFUTADO, no abierto.** Sea Theta = sup{beta: zeta(beta + i*gamma) =
0}. La cota por formula explicita es sup_m |S_1| << N^Theta
(log N)^C. La ruta sup x L^1 exige sup_m |S_rho| =
o(N^{1/2}/(log N)^{3/2}) (doc 69, Prop. 69.2), es decir Theta < 1/2.
Pero existen ceros con beta = 1/2: el primero es gamma_1 = 1/2 +
14.1347251417i (script, mpmath: parte real 0.5 exacta, imaginaria
14.1347251417). Luego Theta >= 1/2 siempre, y ninguna estimacion de
densidad — que solo acota ceros con beta >= sigma > 1/2 — puede dar
Theta < 1/2. Incluso con GRH (Theta = 1/2, lo mejor posible) el sup
ideal N^{1/2} deja el hueco (log N)^{3/2}: 145.2x a N = 10^12 y
410.8x a N = 10^24 (script). La densidad de ceros no puede rescatar
la ruta sup ni en su mejor caso.

**Prop. 76.2 [T]. La densidad no toca la diagonal.** N(sigma,T)
cuenta ceros con beta >= sigma > 1/2. La diagonal del momento mixto
— la que fija la escala N/log N de Cauchy (docs 67/69) — viene de
los ceros EN la linea (beta = 1/2), cuyo conteo N(T) ~
(T/2*pi)*log T lo fija la formula de Riemann-von Mangoldt y ninguna
estimacion de densidad lo mueve. Medido (script, archivo de ceros
del doc 72): 29 ceros con gamma <= 100, N(T) = 29.0 exacto. El mejor
input disponible, Guth-Maynard 2024 (L. Guth y J. Maynard, "New
large value estimates for Dirichlet polynomials", arXiv:2405.20552;
Ann. of Math. 203(2) (2026), 623-675, Teor. 1.2):
N(sigma,T) << T^{30(1-sigma)/13+o(1)} (mejora el 12/5 de Huxley;
da primos en intervalos cortos de longitud x^{17/30+o(1)}), es
irrelevante para la diagonal: acota ceros fuera de la linea, y la
obstruccion esta en la linea. Incluso la hipotesis de densidad
(A(sigma) = 2, lo mas fuerte posible sin RH) dejaria la diagonal
intacta, porque la diagonal son los ceros en linea, no los de
beta > 1/2.

La sub-via por formula explicita (escribir S_1 y S_rho por ceros y
hacer el termino cruzado directo) se reduce a correlacion de pares
(doc 72, sec. 1) — y por el doc 74, ni siquiera esa implicacion esta
establecida como se enuncio: la formula explicita da una estadistica
promediada, no el coeficiente puntual.

**Veredicto A: CERRADA.** Refutada en la via sup (Theta < 1/2 es
falso); reducida a SL-70.1 en la via explicita, y SL-70.1 no implica
el lema como esta enunciado (doc 74).

## 2. Ruta B: otro minorante con mejor cancelacion de Fourier

Pregunta: ¿que propiedad exacta necesitaria S_rho? ¿Existe un
minorante que bata la escala de Cauchy?

**Prop. 76.3 [T]. La norma L^2 esta forzada por los factores
locales.** Sea rho un minorante con los factores locales correctos
(la serie singular es la estandar, doc 67 Prop. 67.1: la masa de rho
en cada clase reducida mod l es c/(l-1) de la masa total) y densidad
c (sum_{p<=N} rho(p) ~ c*N/log N, con rho(p) >= 0 en primos).
Entonces sum_{n<=N} rho(n)^2 >= c^2*N/log N, luego la estimacion de
Cauchy int|S_1 S_rho| <= (pi(N)*sum rho^2)^{1/2} es >= c*N/log N.
Prueba: por Cauchy-Schwarz en cada clase reducida mod l,
sum_{p=a(l)} rho(p)^2 >= (sum_{p=a(l)} rho(p))^2 / pi(N;l,a)
~= c^2*N/((l-1)*log N); sumando las l-1 clases reducidas se obtiene
c^2*N/log N. Medido (script, rho_tilde, X = 2x10^6): con c = 0.465
(doc 66, c_rho asintotico), el piso es c*N/log N = 64,100 =
6.75x el objetivo N/log^2 N = 9,501 (el cociente es c*log N);
asintoticamente c*log N -> infinito. El Cauchy real de rho_tilde es
180,375 = 19.0x, por encima del piso. La verificacion por clases
(l = 3, 5, 7) reproduce el mecanismo: la suma sobre clases de
(sum rho)^2/pi es pi(X)-1 (falta 1: el primo l en clase 0),
muy por encima de c^2*N/log N = 29,806. Conclusion: NINGUN
minorante con la serie singular correcta bate la escala de Cauchy.
La norma L^2 esta forzada por los factores locales, no por la
construccion de rho. (Ojo: esto acota la ESTIMACION de Cauchy, no
la integral misma; la integral podria ser menor, pero ningun
minorante ayuda a la ruta de Cauchy.)

sup: cualquier minorante de criba tiene S_rho con la misma
estructura Tipo I/II que S_1 (doc 69, sec. 2), luego sup <<
N/(log N)^A (Vinogradov) y el ideal N^{1/2} falla (Prop. 69.2b).
Minorantes de caracteres o pesos cuadrados: el argumento de clases
de la Prop. 76.3 solo usa los promedios de rho(p) por clase, no la
construccion — la L^2 sigue forzada; el sup sigue siendo
Vinogradov/Burgess en el mejor caso. Sin ganancia.

Cuarto momento: la diagonal es la escala natural
(sum_p rho(p))^2 ~ c^2*N^2/log^2 N, 246x lo necesario a X = 2x10^6
(script). Es ESCALA, no cota inferior: para rho la cancelacion es
posible en principio porque rho firma en compuestos. Cota inferior
rigurosa solo para las piezas no negativas del defecto (doc 69.3,
verificado: diagonal S_{2,3} = 17.3x, S_{2,5} = 95.8x lo necesario).
Mejor cota disponible: Li-Liu X^{5/2} frente a X^2/log^4 necesario.
Sin mecanismo.

**Veredicto B: CERRADA.** Cauchy forzado por los factores locales
(Prop. 76.3); sup ideal insuficiente (Prop. 69.2b); cuarto momento
sin mecanismo.

## 3. Ruta C: formula explicita dual

La formula explicita dual — escribir S_1(alfa) por ceros de zeta y
S_rho(alfa) por ceros de L(s,chi) y hacer el termino cruzado
directo — es SL-70.1(b): los factores de forma cruzados
G_{chi1,chi2} (doc 72, sec. 1-2). No es una ruta distinta: ES la
correlacion de pares. Y por el doc 74, aun concediendo la meseta,
las condiciones (a)-(c) de SL-70.1 NO implican MM-69.1 ni
E_U(N) = o(N/log^2 N) como estan enunciadas: el factor de forma
controla una varianza promediada, no el producto absoluto ni la fase
de un coeficiente puntual (el contraejemplo finito de la traslacion
ciclica lo muestra: mismo espectro de potencia, coeficiente puntual
distinto). El sustituto honesto es el segundo momento sobre N (74.4)
mas un puente (sec. 4, D5).

**Veredicto C: CERRADA.** No es distinta de SL-70.1; y SL-70.1 ni
siquiera implica el lema como esta enunciado (doc 74).

## 4. Ruta D: atajos que no pasan por el momento mixto

**D1. Valores grandes / Halasz-Montgomery.** Partir m en
{|S_1| > V} y {|S_1| <= V}. En la parte pequena, sup x L^1 necesita
V = o(N^{1/2}/(log N)^{3/2}); en ese V, la cota de Chebyshev para
la medida del conjunto grande es ||S_1||_2^2/V^2 = (log N)^2 > 1:
VACUA (script: 763.5 a N = 10^12, 3053.9 a N = 10^24). Y acotar la
parte grande por Cauchy necesita el cuarto momento (por Holder con
la medida): Li-Liu X^{5/2} frente a X^2/log^4 necesario. El split
no gana nada: se reduce al cuarto momento. CERRADA.

**D2. Sup dependiente de q / arcos mayores extendidos.** Prop. 76.4
[T]: el factor L^1 int_m |S_rho| es independiente del corte. Los
arcos mayores tienen medida ~(log N)^{2B}/N; quitarlos no reduce la
masa L^1 de los menores, que sigue <= (N/log N)^{1/2} por Cauchy.
Y el sup dependiente de q (Vaughan, "The Hardy-Littlewood Method",
2a ed., cap. 3: << (N/q^{1/2} + N^{4/5} + N^{1/2}q^{1/2})(log N)^C)
tiene piso N^{4/5} — el termino Tipo II no baja con q —, es decir
N^{3/10} por encima del ideal N^{1/2} (script: 3,981x a N = 10^12,
1.6x10^7x a N = 10^24). Ningun corte ayuda a la ruta sup x L^1.
CERRADA.

**D3. Holder (p,q).** Por interpolacion, ||S_rho||_q con q < 2 se
controla por ||S_rho||_1 — la misma incognita L^1. No es una ruta
distinta. CERRADA.

**D4. Transferencia / U^2 (Green-Tao).** Da casi-todo-N (conjunto
excepcional E(X) << X^{0.70}; 0.50 bajo GRH, GB-L05), no todo-N; el
salto de casi-todo a todo ES el problema de arcos menores. CERRADA
como ruta a Goldbach fuerte.

**D5. Segundo momento sobre N (74.4) + puente.** (74.4):
sum_{N en [X,2X]} |E_U(N)|^2 = o(X^3/log^4 X). Escala verificada
(script, FFT a X0 = 10^5): sum U(N)^2 = 4.49x X^3/log^4 X y
sum R(N)^2 = 6.65x X^3/log^4 X — (74.4) pide un ahorro genuino
sobre la escala natural, no es vacua ni trivialmente falsa; es la
forma en segundo momento del teorema clasico de casi-todo-N
(Montgomery-Vaughan 1975, via el survey de Bhowmik-Grimmelt). Pero
(i) no se sigue de la correlacion de pares como se enuncio (doc 74:
falta el lema puente con fases), y (ii) aun probada da casi-todo-N;
el puente a todo-N no tiene mecanismo conocido (barrera de paridad
de la criba; no hay continuidad en N). CERRADA como ruta a Goldbach
fuerte; es el enunciado vivo mas cercano al momento mixto.

## 5. Teorema de cierre (Prop. 76.5)

Los tres errores de arcos menores forman una sola incognita.
E_U = E_R - E_D (doc 70, Prop. 70.1), con E_R el error de arcos
menores de Goldbach mismo y E_D el del conteo del defecto. Modulo
la asintotica del defecto D(N) = L_S*S(N)*N/log^2 N + o(N/log^2 N)
(la heuristica medida en el doc 70, sec. 2: D(N)/R(N) = 0.19-0.35
creciendo hacia L_S), se tiene: E_U pequeno + E_D pequeno si y
solo si la asintotica de Goldbach R(N) = S(N)*N/log^2 N +
o(N/log^2 N) (usando c_rho = 1 - L_S, docs 66-67). Luego la version
con fase del momento mixto es EQUIVALENTE, modulo la asintotica del
defecto, a los arcos menores de Goldbach. No hay un "lema mas facil"
en la linea rho: el momento mixto ES el problema, no un peldaño.
(Refuerzo de la Prop. 70.1: la circularidad no es solo que una via
falle, es que el lema mismo es equivalente al problema; la
asintotica del defecto, a su vez, necesita E_D pequeno, asi que el
sistema {E_U, E_D, E_R} con E_U = E_R - E_D no tiene ninguna
incognita mas facil que las otras dos.)

Consecuencia: la linea rho queda cerrada como ruta a una cota
inferior de Goldbach. Lo que queda de ella es un logro real — un
minorante con el termino principal correcto (Prop. 67.1) y los arcos
mayores controlados (nivel 1/2+delta) — pero no puede dar una cota
inferior sin resolver los arcos menores, que son los de Goldbach
mismo.

## 6. Eleccion: la ruta mas prometedora

Todas las rutas al momento mixto por el metodo del circulo estan
cerradas (secs. 1-4). La eleccion es el **segundo momento sobre N
(74.4)**: es el unico enunciado restante sobre el momento mixto que
es (a) preciso, (b) probablemente verdadero (es la forma en segundo
momento del teorema clasico de casi-todo-N), (c) verificable o
refutable desde resultados conocidos, y (d) su resolucion delimita
exactamente lo que queda (el puente casi-todo -> todo).

Plan concreto (proximo doc 77):
1. Formular (74.4) para E_U con la normalizacion exacta (conteo
   unilateral, termino principal c_rho*S(N)*N/log^2 N).
2. Comprobar si (74.4) se sigue del segundo momento clasico de R(N)
   (Montgomery-Vaughan 1975) mas un segundo momento para D(N):
   calcular o acotar sum |E_D(N)|^2 por Parseval y la diagonal.
3. Cuantificar el hueco del puente (casi-todo -> todo) y auditar los
   puentes conocidos (criba/paridad, conjunto excepcional,
   continuidad en N — ninguno funciona).
4. Medir la escala del segundo momento a X = 10^5-2x10^6 por FFT
   (el script de este doc ya mide sum U^2 y sum R^2).

Lo que NO se hace: no se duplica el frente de seleccion del bucle
(C10/W(d), Motohashi doc 57 sec. 4, b_{P,l}, uniones C_k+C_{k+1},
variante Lema 4 Pintz-Althofer) — ese es el frente vivo para
Goldbach por el lado de la criba, y sigue siendo del bucle.

Especificacion para una idea genuinamente nueva (si alguien la
busca): debe producir cancelacion conjunta en
int_m S_1 S_rho e(-N alfa) invisible a toda norma L^p de S_1 y S_rho
por separado, invisible a las correlaciones de pares de ceros, y no
obtenible descomponiendo S_rho — es decir, debe usar la estructura
conjunta rho = 1_p - D con D >= 0 delgado de un modo que la integral
con fase vea pero el momento absoluto no, sin reducirse a E_R - E_D
(doc 70). Ningun fenomeno asi se conoce.

## Referencias (leidas, no de memoria)

- Docs 66, 67, 69, 70, 72, 74, 75 de esta carpeta; GB-L05
  (GRH como teorema para arcos mayores; exponente excepcional
  0.70 -> 0.50).
- L. Guth y J. Maynard, "New large value estimates for Dirichlet
  polynomials", arXiv:2405.20552 (2024); Ann. of Math. 203(2)
  (2026), 623-675, Teor. 1.2: N(sigma,T) << T^{30(1-sigma)/13+o(1)}.
- H. L. Montgomery, "The pair correlation of zeros of the zeta
  function", Proc. Sympos. Pure Math. XXIV (1973), 181-193;
  D. A. Goldston y H. L. Montgomery (1987) — via doc 72.
- G. Bhowmik y L. Grimmelt, "The exceptional set of the Goldbach
  problem" (survey; Montgomery-Vaughan 1975, Pintz 2018).
- E. C. Titchmarsh / R. C. Vaughan, "The Hardy-Littlewood Method",
  2a ed., caps. 2-3: cota de Vinogradov, cota dependiente de q.
- Li-Liu (cuarto momento) — via docs 16/67.
- J. Maynard, arXiv:2006.08250v1, sec. 9 (rho) y Teor. 1.3.

## Alcance honesto

- Lo que este doc SI dice: el mapa completo de rutas al momento
  mixto que no pasan por correlacion de pares, cada una con su
  veredicto y su razon estructural medida; el teorema de cierre
  (el momento mixto es equivalente a los arcos menores de Goldbach);
  la eleccion del segundo momento sobre N con plan concreto.
- Lo que NO dice: no hay prueba del momento mixto, no hay cota
  inferior de capas, no hay deteccion de P_2, no hay demostracion
  de Goldbach.
