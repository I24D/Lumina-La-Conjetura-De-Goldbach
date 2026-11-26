# 75 - SL-65.1: ¿es agudo el umbral 1/10 de Maynard III?

Fecha: 2026-11-26. Siguiente paso del doc 65 (SL-65.1). Script:
`numerico/75_umbral_un_decimo_agudo.py` — PASS 21/21 local (Fractions
exactos). No hay prueba del lema del momento mixto, no hay cota inferior
de capas, no hay demostracion de Goldbach.

## 1. El umbral, con total precision

Marco (docs 64-65): Maynard, "Primes in arithmetic progressions to large
moduli III: Uniform residue classes", arXiv:2006.08250v1 (2020); Mem. Amer.
Math. Soc. 306(1544), 2025. Leido en esta sesion: sec. 11 completa
(Lemas 11.1-11.9 y prueba de la Prop. 5.1), mas el enunciado de la Prop. 5.1.

La Prop. 5.1 (estimacion Tipo II uniforme) exige, con QR = x^(1/2+delta):

  R <= x^(1/10-3delta) / (log x)^C,   Q x^(2delta)(log x)^C <= N <= x^(1/2-3delta)/(log x)^C.

La esquina DUA (doc 65) necesita el factor pequeno con
theta_R = log R / log x en [0.100, 0.105] y theta_Q en [0.405, 0.410].
Pregunta SL-65.1: ¿es el 1/10 una restriccion real del metodo (aguda) o un
artefacto de una estimacion gruesa en la sec. 11?

## 2. De donde sale el 1/10: la cadena exacta

El 1/10 sale de UN SOLO sitio: la cota final del Lema 11.8 (terminos
fuera de la diagonal), ec. (11.47):

  sum_{q~Q} sum_{r0~R0} S6 << (log x)^O(1) * E1^2 N^4 Q^(11/2) R1''^2 R2'^2 R^10 / x^4.

El objetivo (prueba del Lema 11.9) es
sum sum S6 << (E1 R1'' R2' N^2)^2 / ((log x)^A Q^2).
El cociente cota/objetivo es, cancelandose E1, R1'', R2', N:

  Q^(15/2) R^10 / x^4 << 1.                                    (*)

Con QR = x^(1/2+delta), escribir R = x^(1/2+delta)/Q da
Q^(15/2) R^10 / x^4 = x^(1+10delta) / Q^(5/2); pedir << 1 es
Q >> x^(2/5+4delta), es decir R << x^(1/10-3delta). Es la ec. (11.48).
En forma alternativa (enunciado de la Prop. 5.1): Q^3 R^4 <= x^(8/5),
que con QR = x^(1/2+delta) es exactamente R <= x^(1/10-3delta); la otra
condicion alternativa, R^5 Q^6 <= x^3, equivale a Q <= x^(1/2-5delta) y
es consecuencia de R > x^(6delta) (ec. (11.34)), no es ligante.

En exponentes (theta_Q, theta_R): (*) es (15/2) theta_Q + 10 theta_R <= 4.
El 1/10 es 1/2 - 2/5: el 1/2 es el nivel QR = x^(1/2+delta) y el 2/5 es
el exponente que sale de Q^(5/2) >> x^(1+10delta).

De donde salen los exponentes de (11.47): Q^(11/2) = Q^(3/2) (suma sobre
q0 tras la cota de Weil para sumas de Kloosterman) por Q^4 (H''^4, cuatro
variables duales tras completar sumas); R^10 viene de las potencias de
R2' en las sumas de divisores y de R1''^3 en (11.44). El "4" es el numero
de variables duales (l1, l1', l2, l3): una por cada suma-n tras los DOS
Cauchy-Schwarz (Lemas 11.4 y 11.6). El "1/2" de Weil es la cancelacion de
raiz cuadrada para sumas de Kloosterman.

## 3. El paso de dominacion no pierde nada en el peor caso

Prop. 75.1 [T]. La simplificacion (11.46) -> (11.47) no afloja el umbral:
en el rango diadico peor (R0 = d1 = E1 = 1, R1'' ~ R) la cota (11.47)
coincide con el termino T1 de (11.44) y el termino T2 de (11.45) no es
ligante en la esquina.

*Prueba.* (11.46) es (log x)^O(1) [T1 + T2] + x^o(1) [E1^2 T1/R + E1 T2/R]
con T1 = R1''^3 R2'^2 Q^(11/2) N^4 R^9/x^4 y
T2 = R1''^2 R2'^2 N^3 Q^(9/2) R^(15/2)/x^3. El paso a (11.47) multiplica
por E1^2 y se queda con el primer termino. El cociente (11.47)/(E1^2 T1)
es R/R1''; con R0 = d1 = E1 = 1 se tiene E1 R1'' ~ R/(R0 d1) = R, luego
R1'' ~ R y el cociente es 1: sin perdida. Ese rango es admisible en el
Lema 11.8: R1'' >= (N/QR0)^(2/3) se cumple (script: en el centro de la
esquina, con N maxima, (N/Q)^(2/3) tiene exponente 0.0417 < 0.1025 =
theta_R). El cociente (11.47)/(E1^2 T2) es NQR^(5/2)/x >= R^(1/2) >= 1
(usando NQR^2 >= x), luego (11.47) domina a T2. Y T2 por si solo no es
ligante en la esquina: T2/objetivo = Q^(13/2) R^(15/2)/(N x^3), cuyo
exponente en el centro es -0.0525 con N maxima y -0.0100 con N minima
(script): < 1 en todo el rango de N. QED.

Corolario: el umbral no es un artefacto del paso (11.46) -> (11.47).

## 4. El hueco en la esquina es una potencia, no un logaritmo

Prop. 75.2 [T]. En TODA la esquina la desigualdad ligante falla por una
potencia de x: el hueco gamma = (15/2) theta_Q + 10 theta_R - 4 es

  fondo-izq (0.405, 0.100): gamma = 3/80   = 0.0375
  fondo     (0.41,  0.100): gamma = 3/40   = 0.0750
  centro    (0.4075,0.1025): gamma = 13/160 = 0.08125
  tope      (0.41,  0.105): gamma = 1/8    = 0.1250

(script, Fractions exactos). El maximo theta_R permitido a theta_Q =
0.4075 es 0.094375, frente al 0.1025 que pide la esquina: hueco 0.008125
en exponente. Ningun refinamiento de nivel logaritmico en la sec. 11
(promedios de divisores, promedios de mcd, constantes) puede cerrar un
hueco x^gamma con gamma >= 0.0375.

## 5. Por que los exponentes son forzados (no hay holgura estructural)

Prop. 75.3 [T]. Dentro del metodo de la sec. 11 (dos Cauchy-Schwarz +
completar sumas + Weil, factor por factor), los exponentes (15/2, 10, 4)
son forzados. Tres razones independientes:

(a) Weil es agudo para modulos primos. Las sumas de Kloosterman que
aparecen son S(k, l; r) con r primo en nuestra aplicacion (la esquina
usa r primo; en la Prop. 5.1 r es libre de cuadrados con tau(qr) <=
(log x)^B). La cota de Weil |S| <= 2 sqrt(r) es la mejor posible en
general: existen a, b con |S(a,b;p)| proximo a 2 sqrt(p). El metodo
necesita la cota uniforme sobre los argumentos (k, l varian), asi que el
peor caso gobierna. Promediar sobre l no ayuda: la suma sobre l se hace
con valores absolutos dentro (ec. (11.41)), y la estimacion trivial
(numero de l) x (Weil) es la que vale; una cancelacion conjunta en l
seria una segunda completacion de sumas, es decir, volver a la suma
original.

(b) El "4" es el numero de variables duales y lo fijan los dos Cauchy.
Cada Cauchy duplica las variables-n; tras el segundo hay 4 sumas-n que
completar (l1, l1', l2, l3), cada una de longitud H'' ~ QR/N, y cada una
aporta un factor x^-1 al denominador (de ahi x^4) y un factor Q R al
numerador (de ahi Q^4 R^4 en H''^4). Quitar un Cauchy exigiria eliminar
los coeficientes alpha' sin cuadrar, lo que el metodo no sabe hacer: el
Lema 11.4 lo dice explicitamente ("We cannot simultaneously eliminate
the alpha_n2 coefficients because the modulus would increase too much").

(c) Los pasos gruesos de la sec. 11 solo pierden logs. Los promedios de
divisores (Lema 10.9, sumas de tau(k0') y de d2') y los promedios de mcd
(d2'^(1/2) d3'^(1/2) <= d2' + d3') pierden a lo sumo (log x)^O(1) y
x^o(1); las condiciones de los Lemas 11.1-11.2 ((11.4)-(11.7)) son de
nivel log. Nada de eso mueve un exponente.

## 6. Veredicto: el umbral ES agudo para el metodo de la sec. 11

**El 1/10 es agudo.** No es un artefacto de una estimacion gruesa: es la
solucion exacta de la desigualdad ligante Q^(15/2) R^10 <= x^4 del Lema
11.8, con exponentes forzados por la arquitectura del metodo (dos
Cauchy, cuatro variables duales, Weil agudo para modulos primos), y el
hueco en la esquina es una potencia x^gamma con gamma >= 0.0375, no un
logaritmo. El paso de dominacion (11.46) -> (11.47) no esconde holgura
(Prop. 75.1) y el segundo termino no es ligante (Prop. 75.2, sec. 3).

Las unicas salidas estructurales son:

(i) Una cota CONJUNTA (Deligne/Weil II) para el producto de las cinco
sumas de Kloosterman en vez de acotarlas factor por factor. Eso es el
metodo espectral (Maynard I, BFI), que introduce dependencia en la
clase residual a y rompe la uniformidad completa en a — justo lo que
Maynard III existe para conseguir (Introduccion: el metodo espectral
"necessarily introduces a dependence on the residue class"). No se puede
usar aqui sin perder la uniformidad.

(ii) Una arquitectura de dispersion nueva con menos pasos de Cauchy o
menos variables duales. No se conoce ninguna.

## 7. Que ganaria la linea si el umbral se moviera

Si la cota S6 mejorara por un factor x^-gamma (gamma > 0), el umbral
seria theta_R <= 0.094375 + gamma/10 a theta_Q = 0.4075. Para cubrir la
esquina haria falta gamma >= 0.0563 (fondo, theta_R = 0.100), >= 0.0813
(centro) o >= 0.1063 (tope) (script). Pero aun con el umbral movido:

- La Prop. 5.1 extendida a la esquina daria el error Tipo II de la
  esquina con BUEN error ((log x)^-A, el de la Prop. 5.1 misma) — el
  error debil O(delta) es del Teor. 1.1 (paso a primos via Prop. 5.4 y
  criba), no de la Prop. 5.1. Asi que el lema DUA-esquina (doc 65) SI se
  cerraria.
- Pero la esquina es el 0.004% de M1 (Prop. 65.4): aun probado el lema,
  NO daria cota inferior de capa. La ganancia seria un rango Tipo II
  uniforme mas amplio (un teorema en si mismo) y el cierre del lema
  DUA-esquina, no un avance hacia Goldbach.

SL-65.1 queda CERRADO con veredicto: agudo. La linea de dispersion
uniforme (docs 59, 64, 65) queda delimitada en el umbral 1/10 de
Maynard III como su frontera exacta para el metodo conocido.

## 8. Referencias (leidas en esta sesion, salvo indicacion)

- J. Maynard, "Primes in arithmetic progressions to large moduli III:
  Uniform residue classes", arXiv:2006.08250v1 (2020); Mem. Amer. Math.
  Soc. 306(1544), 2025. [Leido: sec. 11 completa (Lemas 11.1-11.9,
  prueba de la Prop. 5.1), enunciado de la Prop. 5.1 y su forma
  alternativa; Introduccion y sec. 5 leidas en el doc 65.]
- Docs 58, 59, 63, 64, 65, 72 de esta carpeta; script
  numerico/75_umbral_un_decimo_agudo.py (PASS 21/21 local).

## Alcance honesto

- Lo que este doc SI dice: el origen exacto del umbral 1/10 (Lema 11.8,
  ecs. (11.47)-(11.48)); la verificacion de que el paso de dominacion no
  esconde holgura en el peor caso; el hueco en la esquina medido como
  potencia; las tres razones por las que los exponentes son forzados; el
  veredicto de agudeza y lo que ganaria la linea si el umbral se moviera.
- Lo que NO dice: no hay prueba ni refutacion del lema DUA-esquina por
  otra via; no hay cota inferior de capas, no hay deteccion de P2, no hay
  demostracion de Goldbach. SL-65.1 queda cerrado.
