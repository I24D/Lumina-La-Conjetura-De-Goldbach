# 79 - C10 / peso movil W(d): el termino principal exacto M_B3

Fecha: 2026-11-26. Script: `numerico/79_peso_movil_Wd.py` — PASS 26/26 local
(X = 10^6). No hay prueba del lema del momento mixto, no hay cota inferior
de capas, no hay demostracion de Goldbach.

## 0. La pregunta (C10, docs 04/05/06/14, GB-L18)

En una caja diadica (C, K) de R2 (K > sqrt(N)/L^B, CK ~= N), el termino
Tipo II de Vaughan es

    B3 = sum_{c~C} alpha_c sum_{d~K, cd<N} beta_d S(N-cd),
    alpha_c = 1_{(c,N)=1}, beta_d = mu(d) log d,
    S(m) = sum_{ab=m, a>V, b>U} Lambda(a) nu_U(b),  nu_U(b) = sum_{t|b, t>U} mu(t),
    U = V = N^{1/8}/L^{B3}.

GB-L18 (C10) establecio que el termino principal de B3 por caja es
<< N/L^A (via sum_{d~K} mu(d) log d << K/L^A), asi que (1-kappa) no puede
venir de la Prop. 3.1 de Lichtman bajo ninguna lectura. Este doc formula
el termino principal REAL de B3, con su peso movil W(d) explicito, acota
su variacion y decide si esa variacion puede dar (1-kappa).

## 1. La densidad global de S: identidad exacta

**Lema 79.0 (exacto).** Como convolucion de Dirichlet,
nu_U = mu_{>U} * 1 y Lambda * 1 = log, luego

    Lambda * nu_U = log * mu_{>U}   (exacto).

Por tanto, con S_nu(U,1) := sum_{t>U} mu(t)/t,

    (1/x) sum_{m<=x} S(m) = S_nu(U,1) * log x + O((log x)/L^A)      (79.1)

para todo A (PNT para mu: sum_{t<=y} mu(t)/t << 1/(log y)^A,
sum_{t<=y} mu(t) log t/t = -1 + O(1/(log y)^A)).

Verificacion (script, X = 10^6): la media de la convolucion completa
coincide con la formula exacta sum_{t>U} mu(t) * log(floor(X/t)!) / X
hasta 4e-15 (U = 10: -0.9437; U = 30: 0.2835). La media de S con a > V
es -0.7904 (U = 10) y 0.1972 (U = 30); la parte a <= V es -0.1533 y
0.0863. La media de nu_U decae con U (PNT): nu(U) = -0.0905 (U=10),
0.0305 (30), -0.0311 (100), 0.0105 (300), -0.0044 (1000): 20 veces
menor de U=10 a U=1000.

**Honestidad de escala.** A N = 10^6 el regimen asintotico NO esta
alcanzado (en la teoria U = N^{1/8}/L^{B3} ~= 5). Los numeros verifican
las identidades exactas y los mecanismos, no la pequeez asintotica
(nu << 1/L^A, sum mu(d) log d << K/L^A), que es PNT citada, ilustrada
por la caida de nu(U) y de las sumas diadicas (secc. 2).

## 2. La densidad de clase y W(d) (Prop. 79.1)

Sea (N, d) = 1 y T_d la suma de S sobre la clase m = N (mod d) en el
intervalo I_d (longitud C*d, m = N - cd):

    T_d = sum_{c~C, (c,N)=1} S(N - cd)
        = sum_{m en I_d, m = N (mod d), (m,N)=1} S(m).

Derivacion de la densidad de clase. En S(m) = sum_{ab=m} Lambda(a) nu_U(b),
con m = N (mod d) y (N,d)=1: si d | b entonces ab = 0 (mod d), imposible,
asi que solo (a,d) = (b,d) = 1 contribuyen. Para (a,d) = 1, la suma de
clase de nu_U es

    sum_{b<=y, b = s (mod d)} nu_U(b) = (y/phi(d)) * S_nu(U,d) + E,
    S_nu(U,d) := sum_{t>U, (t,d)=1} mu(t)/t,

donde E es el error de redondeo sum_t mu(t) O(1). La suma armonica de
Lambda sobre clases reducidas aporta (phi(d)/d) log x. El factor phi(d)/d
de la suma de Lambda y el 1/phi(d) de la suma de nu_U se cancelan, y
queda la densidad de clase

    delta_d := (d/x) sum_{m<=x, m = N (mod d)} S(m) = S_nu(U,d) * log x,
    delta   := (1/x) sum_{m<=x} S(m)            = S_nu(U,1) * log x.

**Prop. 79.1 (definicion exacta de W(d)).**

    W(d) := delta_d / delta = S_nu(U,d) / S_nu(U,1)
          = [sum_{t>U, (t,d)=1} mu(t)/t] / [sum_{t>U} mu(t)/t].

W(d) es el cociente de la cola de Mertens restringida a (t,d)=1 entre
la cola total. Depende de d (por la restriccion (t,d)=1) y de U.

**Correccion de una conjetura ingenua.** W(d) NO es phi(d)/d.
phi(d)/d es solo el limite heuristico (U -> infinito): el cociente de
las series de Dirichlet en s -> 1+ es d/phi(d) (factor de Euler), pero
la convergencia es lenta (como 1/(log U)^A) y a U finito el cociente
esta dominado por las colas parciales, que oscilan. Medido a U = 30,
X = 10^6 (script secc. 3, teoria calculada exacta desde la criba, sin
ajuste):

    d=7:  ratio S medido 0.370, teoria W = 0.219, phi(d)/d = 0.857
    d=11: ratio S medido -0.404, teoria W = -0.492, phi(d)/d = 0.909
    d=30: ratio S medido -16.661, teoria W = -16.403, phi(d)/d = 0.267

La teoria reproduce signo y orden en los tres casos; d = 30 al 2%.
phi(d)/d falla en los tres, incluido el SIGNO en d = 11 y d = 30.
La referencia Lambda si cumple d/phi(d) en los tres (1.167, 1.103,
3.757 medido contra 1.167, 1.100, 3.750): el peso movil es propio de
nu_U, no de Lambda.

**Mecanismo del movimiento.** nu_U se concentra en b = 0 (mod d):
a d = 7, la clase 0 lleva el 81% del total de nu_U; a d = 11, el 145%
(las reducidas son negativas). La concentracion viene de sumas cortas
sum_{U/d < t <= U} mu(t)/t = O(1) a U finito (exacto a d = 7, U = 30:
la clase 0 predicha 24729, medida 24762). Esa clase no afecta a
T_d ((N,d)=1 implica b != 0 (mod d)), por eso la formula de W(d)
sigue valiendo; lo que mueve W(d) son las colas restringidas
S_nu(U,d), no la concentracion en cero.

## 3. El termino principal exacto (Prop. 79.2)

De la secc. 2, con |I_d| ~= C*d y el factor phi(N)/N de la restriccion
(c,N)=1 (equivalente a (m,N)=1 dado (d,N)=1):

    T_d ~= C * (phi(N)/N) * S_nu(U,d) * log N.

**Prop. 79.2 (M_B3 exacto).**

    M_B3(caja) = C * (phi(N)/N) * L * sum_{d~K} mu(d) log d * S_nu(U,d)
                 + error,

con L = log N. La cola total S_nu(U,1) se CANCELA exactamente: el
peso que entra en M_B3 es S_nu(U,d), no el cociente W(d). El factor
phi(N)/N es la aproximacion de equidistribucion; las correcciones por
las clases p | N son de orden inferior en el cociente (entran por los
e | d con e | N en la suma de Mobius de la secc. 4).

Verificacion en caja (N = 10^6, C = 512, K = 1024, CK = 524288 < N;
script secc. 4): B3/N = -0.0178, M_B3/N = -0.0073 (|M|/|B3| = 0.41).
El termino ingenuo con phi(d)/d da +0.0001: signo erroneo y 70 veces
menor. La formula reproduce signo y orden; a esta escala M_B3 NO es
despreciable frente a B3 (nu(30) = 0.03 aun no es pequeno) - la
pequeez es asintotica. Tampoco se pretende que M_B3 aproxime a B3:
B3 = M_B3 + error, y el error (la fluctuacion Tipo II) no esta
controlado; es justamente la dificultad.

## 4. Variacion de W(d) (Prop. 79.3)

**Prop. 79.3.** (a) Heuristicamente (U -> infinito), W(d) -> d/phi(d)
en [1, O(log log N)]: variacion log log N, acotada. (b) A U finito,
W(d) es un cociente de dos colas << 1/(log U)^A: puede ser negativo
y |W| >> 1 (medido -16.7 en d = 30); en los U donde S_nu(U,1) ~= 0,
|W(d)| es arbitrariamente grande. La variacion del cociente NO esta
acotada. (c) Pero esa variacion es un artefacto del cociente: lo que
entra en M_B3 es S_nu(U,d), y por inversion de Mobius de la
restriccion (t,d)=1,

    S_nu(U,d) = sum_{e|d, mu(e)!=0} (1/e) * S_nu(U/e, e),

de donde, partiendo e <= sqrt(U) y e > sqrt(U),

    |S_nu(U,d)| <<_A (log log d) / (log U)^A     para todo A,

para U >= (log N)^C. De ahi, por valores absolutos,

    |M_B3| << N * (log log N) / (log N)^{A-2},

cota cruda; el tamano real es menor por la cancelacion de mu en
sum_{d~K} mu(d) log d * S_nu(U,d) (medido: |sum|/K = 0.0230, 0.0184,
0.0000 para K = 1000, 10000, 100000, frente a |sum mu log d|/K =
0.0263, 0.0471, 0.0053; y la media de phi(d)/d -> 6/pi^2 = 0.6079 en
los tres K).

## 5. Veredicto sobre (1-kappa) (Prop. 79.4)

**Prop. 79.4.** La variacion de W(d), por salvaje que sea, NO puede
dar (1-kappa): se cancela exactamente de M_B3. El termino principal
es pequeno por las mismas dos razones de C10, ahora en forma afilada:
la media de S (cola de Mertens S_nu(U,1) << 1/L^A) y la suma con mu
(sum mu(d) log d * S_nu(U,d) << K/L^A por PNT). W(d) no amplifica
nada: es un cociente, no un factor. Por tanto (1-kappa) NO puede venir
del termino principal de B3 bajo ninguna lectura - ni de la Prop. 3.1
de Lichtman (C10) ni del termino principal real con W(d). Debe venir
del termino de error (la fluctuacion de B3 alrededor de M_B3), es
decir, del error Tipo II mismo. Consistente con doc 06 secc. 8.2 (la
cota esta exactamente en el umbral 1) y con la Prop. 76.5 (el momento
mixto equivale a los arcos menores).

**Que sigue.** El frente vivo sigue siendo la seleccion del bucle
(switching, Motohashi doc 57 secc. 4, b_{P,ell}, uniones, Pintz-
Althofer). C10 queda cerrado aqui como "no sale del termino
principal". Goldbach sigue abierta.
