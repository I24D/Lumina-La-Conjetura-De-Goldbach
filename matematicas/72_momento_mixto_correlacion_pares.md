# 72 - Sub-lema SL-70.1: la correlacion de pares de ceros suficiente para el momento mixto

Fecha: 2026-11-26. Siguiente paso del doc 70 (SL-70.1). Script:
`numerico/72_momento_mixto_correlacion_pares.py` — PASS 13/13 local
(252 ceros de zeta; R2 verifica GUE en 4 binas). No hay prueba del
lema del momento mixto, no hay cota inferior de capas, no hay
demostracion de Goldbach.

## 1. Por que la correlacion de pares

El lema MM-69.1 (doc 67 sec. 3, doc 69) es
int_m |S_1 S_rho| = o(N/log^2 N), y su version con fase es
E_U = int_m S_1 S_rho e(-N alfa) = o(N/log^2 N).
El doc 69 cerro las tres estimaciones por piezas (Props. 69.1-69.3)
y el doc 70 cerro la via conjunta con fase: Prop. 70.1
(circularidad), E_U = E_R - E_D, donde E_R es el error de arcos
menores de Goldbach mismo. La unica via que queda en pie es la
formula explicita: en arcos menores, S_1(alfa) y S_rho(alfa) se
escriben como sumas sobre ceros de zeta (para S_1) y sobre ceros
de L(s,chi) (para las piezas bilineales de rho, doc 66), y el
momento mixto se convierte en una suma doble sobre pares de ceros.
La diagonal (gamma = gamma') da por si sola la escala N/log N de
Cauchy (docs 67/69): la correlacion de pares no puede ayudar a la
ruta de Cauchy, solo a la ruta con fase, y esa ruta es circular
(doc 70) salvo que la correlacion controle E_R (y E_D) por la
formula explicita. El sub-lema SL-70.1 es exactamente ese input.

La autoridad para la equivalencia es Goldston-Montgomery (1987):
la correlacion de pares de ceros (PCC) es EQUIVALENTE a la
asintotica de varianza de primos en intervalos cortos. La varianza
es el control de arcos menores que E_R necesita (Prop. 70.1).
Para el momento mixto se necesita la version mixta: ceros de zeta
contra ceros de L(s,chi), incluidos los terminos cruzados entre
caracteres distintos.

## 2. SL-70.1: enunciado preciso

**SL-70.1 (forma F).** Sea F(alfa;T) el factor de forma de
Montgomery para zeta,
  F(alfa) = ((T/2 pi) log T)^{-1} sum_{0<gamma,gamma'<=T}
            T^{i alfa (gamma-gamma')} w(gamma-gamma'),
  w(u) = 4/(4+u^2),
F_chi(alfa;T) el analogo para L(s,chi) (chi primitivo mod q), y
G_{chi1,chi2}(alfa;T) el factor de forma cruzado entre ceros de
L(s,chi1) y L(s,chi2). Entonces MM-69.1 (y E_U) se sigue de:
  (a) F(alfa) = 1 + o(1) uniformemente para 1 <= alfa <= A,
      para todo A >= 1 fijo (la meseta);
  (b) F_chi(alfa) = 1 + o(1) uniformemente para 1 <= alfa <= A y
      uniformemente en chi primitivo mod q, q primo <= Q_0 =
      (log N)^B (B fijo, el de los arcos mayores del doc 70),
      con G_{chi1,chi2}(alfa) = o(1) en el mismo rango y
      uniformidad (ceros de L distintas no correlacionados);
  (c) todo ello uniformemente en la clase a, en particular para
      a = N, que varia con N (la misma obstruccion de
      uniformidad en a del doc 64).
Observacion: la version con A fijo es ya la conjetura de
Montgomery; la aplicacion pediria ademas uniformidad en alfa
hasta ~log N / log log N (el rango de alfa crece con N porque
Q_0 = (log N)^B, asi theta -> 1 y el nucleo de Goldston-Montgomery
se extiende a alfa grande). La version con A fijo es la version
MINIMA que falta; la completa es mas fuerte aun.

**SL-70.1 (forma varianza, equivalente).** Por Goldston-Montgomery
(1987) para zeta y por Liu (2026) para el caso Dirichlet, (a) es
la asintotica de varianza de psi en intervalos cortos,
  V~(X,h) ~ h X log(X/h),
y (b) es la conjetura suave de Hooley para la varianza de primos
en progresiones,
  G_W(x,q) ~ ||W||_2^2 x log q
uniformemente para x^epsilon <= q <= x, en los rangos que los
arcos menores exigen. La version mixta (zeta contra L) es la que
el momento mixto necesita; la forma F y la forma varianza son el
mismo input en dos lenguajes.

## 3. Lo que existe en la literatura

Solo resultados leidos; nada de memoria.

- **Montgomery (1973), "The pair correlation of zeros of the zeta
  function", Proc. Sympos. Pure Math. XXIV, 181-193.** Bajo RH,
  F(alfa) = T^{-2|alfa|} log T (1+o(1)) + |alfa| + o(1)
  uniformemente para |alfa| <= 1. Es decir: la parte |alfa| < 1
  esta DEMOSTRADA (bajo RH). La meseta |alfa| >= 1 queda como
  conjetura: F(alfa) = 1 + o(1). El metodo no la alcanza.
- **Goldston-Montgomery (1987), "Pair correlation of zeros and
  primes in short intervals", Prog. Math. 70, 183-203.** La PCC
  (meseta) es equivalente a la asintotica de varianza de primos
  en intervalos cortos. Sin la meseta no hay varianza en el rango
  completo; con |alfa| < 1 solo se obtiene el rango restringido.
- **Baluyot-Goldston-Suriajaya-Turnage-Butterbaugh (2024),
  "An unconditional Montgomery theorem for pair correlation of
  zeros of the Riemann zeta-function", Acta Arith. 214, 357-376
  (arXiv:2306.04799), Teor. 1.** La parte |alfa| <= 1 es ahora
  INCONDICIONAL: F(alfa) = T^{-2 alfa}(log T + O(1)) + alfa +
  O(sqrt(log T)) uniformemente para 0 <= alfa <= 1, sin RH.
  La meseta sigue siendo conjetura. El avance incondicional no
  toca el regimen que SL-70.1 necesita.
- **Chan (2000), "On the pair correlation of zeros of the Riemann
  zeta-function", Proc. London Math. Soc.** Bajo GRH,
  F(alfa) >= 3/2 - |alfa| - epsilon para 1 <= |alfa| <= 3/2.
  Es una COTA INFERIOR: la direccion equivocada. SL-70.1 necesita
  la asintotica (una cota superior), y Chan no la da.
- **Ozlük (1982, tesis; J. Number Theory 59 (1996), 319-351) y
  Yildirim (1991, Manuscripta Math. 72, 325-334).** El analogo
  Dirichlet de |alfa| < 1 bajo GRH (Ozlük: >= 11/12 de los ceros
  simples). La meseta Dirichlet tambien es conjetura.
- **Kandhil-Languasco-Moree (2026), Math. Ann. 394, 43.**
  Correlacion de pares de ceros de L de Dirichlet, incluidos los
  terminos cruzados G_{chi1,chi2}; Conjetura 1 para q = 1 o q
  primo, 0 < eta <= 1. El marco cruzado que (b) necesita existe
  como conjetura, no como teorema.
- **Liu (2026), "Pair correlation of L-zeros and Hooley's
  conjecture", preprint nov. 2026 (Hexagon).** Factor de forma de
  ceros de L(s,chi), promediando sobre caracteres primitivos mod
  q primo. Incondicionalmente |alfa| < 1 coincide con CUE. Bajo
  GRH + primos gemelos de Hardy-Littlewood, |alfa| < 2. Para
  |alfa| >= 1, la prediccion CUE es EQUIVALENTE a la conjetura
  suave de Hooley para q >= x^epsilon. Es la equivalencia exacta
  que usa la forma varianza de SL-70.1.
- **Fiorilli-Martin (2023), "Disproving Hooley's conjecture",
  J. Eur. Math. Soc. 25, 4791-4812, Teor. 1.1.** La conjetura de
  Hooley en su forma original (G(x;q) << x log q para todo
  q -> infinito, x >= q) es FALSA: el bound falla para
  q ~= log log x, y G(x;q) puede ser tan grande como
  x (log q + log log log x)^2 / 4. Ademas, Hooley en un rango
  q ~= x^delta implicaria la franja libre de ceros Re(s) > 1-delta
  para todas las L mod q (mas fuerte que GRH). Hayani (2025)
  extiende la refutacion a la version suavizada para q < x^epsilon.
  Consecuencia: la meseta Dirichlet es equivalente a una conjetura
  abierta (Hooley suave, q >= x^epsilon) cuya extension natural a
  q pequeno es FALSA. La verdad es delicada: ningun argumento que
  funcione tambien para q pequeno puede dar la meseta.
- **Leung (2026).** Bajo las conjeturas de k-tuplas de
  Hardy-Littlewood, psi(x;q,a) es gaussiana con varianza
  x log q / phi(q)^2 para x^{1-epsilon} <= q <= x/log x, cuando
  a VARIA. Es un promedio sobre a, no uniformidad en a fijo:
  no da (c) (a = N fijo, variando con N).

## 4. Veredicto

**SL-70.1 NO SE TIENE.** Lo que existe es la parte |alfa| < 1
(Montgomery 1973 bajo RH; Baluyot et al. 2024 incondicional) y la
conjetura de la meseta |alfa| >= 1, abierta incluso bajo GRH para
zeta. Para L de Dirichlet, la meseta es equivalente (Liu 2026) a
la conjetura suave de Hooley en q >= x^epsilon, abierta, y cuya
extension a q pequeno es falsa (Fiorilli-Martin 2023). La
uniformidad en a (c) es una capa abierta adicional: lo mejor
disponible promedia sobre a (Leung 2026), no fija a = N.

Por tanto SL-70.1 no es un lema: es (una version uniforme de) la
conjetura de correlacion de pares de Montgomery. Pedirla como
input del momento mixto es pedir la conjetura. No hay refutacion
incondicional (la correlacion es compatible con todo lo conocido),
pero tampoco hay ninguna ruta conocida que la pruebe en el rango
que necesitamos: el regimen |alfa| < 1 no basta (es el regimen
equivocado: gobierna intervalos cortos, no los arcos menores con
q hasta N), y la meseta es exactamente lo que falta. La linea rho
queda asi delimitada con precision: su cuello de botella es la
conjetura de Montgomery, en version uniforme Dirichlet.

## 5. Verificacion numerica (demostracion finita, no prueba)

Script `numerico/72_momento_mixto_correlacion_pares.py`, PASS 13/13.
Primeros 252 ceros de zeta (barrido Riemann-Siegel vectorizado +
biseccion con Z exacta de mpmath; 4 ceros cercanos recuperados por
reparacion de huecos). Validacion: gamma_1, gamma_2, gamma_3 con
error <= 3e-10 contra valores conocidos; gamma_100 =
236.524229666 (error 1e-11); conteo 252 contra N(T) = 251.47.

- **R2 (GUE) con normalizacion local** (espaciamientos de Odlyzko):
  [0,1): 140 vs 138.2 (rel 0.01); [1,1.5): 113 vs 122.4 (0.08);
  [1.5,2): 121 vs 123.6 (0.02); [2,3): 251 vs 249.9 (0.00).
  La funcion de correlacion 1-(sin pi u/pi u)^2 se verifica a
  esta altura en u >= 0.5.
- **F(alfa):** F(0) = 1.97 contra la prediccion finita-T
  log(T/2 pi e)*log(T/2 pi)/log T = 2.33 (no log T = 6.16: el
  factor finito-T log(T/2 pi e)/log T = 0.539 domina). Para
  alfa > 0: F = 0.31, 0.43, 0.58, 0.53, 0.51, 0.58, 0.53 en
  alfa = 0.25, 0.5, 0.75, 1, 1.5, 2, 3: todos ~0.5, sin rastro
  de la transicion asintotica (alfa para alfa < 1, 1 para
  alfa >= 1). La forma F converge logaritmicamente lenta: a
  T < 500 el regimen asintotico no se ve.
- **Diagnostico (no es check):** la sub-bina [0,0.5) da 13 contra
  28.5 esperados: las separaciones mas pequenas estan depletadas
  respecto a GUE a esta altura (efecto de altura finita, uniforme
  en bandas de altura; las binas finas muestran la deplecion
  concentrada en u < 0.3). Es la region donde la convergencia a
  GUE es mas lenta.

Leccion numerica: los datos finitos confirman GUE localmente pero
NO pueden dirimir la meseta: la correlacion que SL-70.1 necesita
es un enunciado limite, y a alturas accesibles su forma F ni
siquiera muestra la transicion. Esa es exactamente la razon por
la que es una conjetura abierta y no un lema disponible.

## 6. Veredicto y que sigue

**Veredicto: la correlacion de pares suficiente para el momento
mixto es la meseta de Montgomery (|alfa| >= 1) en version
uniforme Dirichlet, y NO SE TIENE.** Existe |alfa| < 1
(incondicional desde 2024); la meseta es conjetura abierta incluso
bajo GRH; su forma Dirichlet equivale a Hooley suave (abierta,
con extension falsa a q pequeno); la uniformidad en a = N es otra
capa abierta. SL-70.1 queda formulado con precision y cerrado
como ruta: no es un paso tecnico pendiente, es la conjetura misma.

Que sigue: los frentes del bucle sin tocar (seleccion: switching,
C10/W(d); Motohashi doc 57 sec. 4; b_{P,l} doc 53; uniones
C_k+C_{k+1}; variante Lema 4 Pintz-Althofer; SL-65.1: agudeza del
umbral 1/10); o una idea genuinamente nueva para el momento mixto
que no pase por la correlacion de pares. Ninguna de esas vias se
duplica aqui. Goldbach sigue abierta.

Referencias (leidas, no de memoria):
- Docs 66 (rho), 67 (conteo unilateral, Prop. 67.1), 69 (etapa 1:
  Props. 69.1-69.3), 70 (etapa 2: Prop. 70.1, circularidad),
  64-65 (dispersion, DUA-esquina).
- H. L. Montgomery, "The pair correlation of zeros of the zeta
  function", Proc. Sympos. Pure Math. XXIV (1973), 181-193.
- D. A. Goldston y H. L. Montgomery, "Pair correlation of zeros
  and primes in short intervals", Prog. Math. 70 (1987), 183-203.
- S. A. C. Baluyot, D. A. Goldston, A. I. Suriajaya y
  C. L. Turnage-Butterbaugh, "An unconditional Montgomery theorem
  for pair correlation of zeros of the Riemann zeta-function",
  Acta Arith. 214 (2024), 357-376 (arXiv:2306.04799), Teor. 1.
- Chan, "On the pair correlation of zeros of the Riemann
  zeta-function", Proc. London Math. Soc. (2000).
- A. E. Ozluk, tesis (1982); J. Number Theory 59 (1996), 319-351.
- Yildirim, Manuscripta Math. 72 (1991), 325-334.
- Kandhil-Languasco-Moree, Math. Ann. 394 (2026), 43.
- T. Z. Liu, "Pair correlation of L-zeros and Hooley's
  conjecture", preprint nov. 2026 (Hexagon).
- D. Fiorilli y G. Martin, "Disproving Hooley's conjecture",
  J. Eur. Math. Soc. 25 (2023), 4791-4812, Teor. 1.1.
- Hayani (2025), analogo suavizado de Fiorilli-Martin.
- Leung (2026), momentos superiores y ley gaussiana en a.
