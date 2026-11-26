# 64. Dispersion uniforme en `a` para pesos u*v con factores >= N^{1/10}: busqueda, obstruccion y veredicto

Fecha: 2026-11-25. Estado: Goldbach fuerte sigue abierto. Este doc cierra el problema
independiente doc 59 #1 (dispersion uniforme en `a`) como ruta: no existe en la
literatura en la forma que necesitamos, y se demuestra por que ningun metodo conocido
la puede dar. Lo que queda delimitado es el lema faltante exacto (sec. 6).

## 1. El enunciado que haria falta (DUA)

Marco: docs 54-59, 63. Peso relajado c_t, familia A del doc 56 (t = 0.05), malla
alpha en [1/v, 1/2] con v = 10, primos representativos q_i >= N^{alpha_i},
d = q_i q_j, M1 = sum_{i,j} c_t(q_i,q_j)/phi(q_i q_j), M = N*M1. Slice (doc 59):
1/2 < alpha+beta <= 1/2+delta, delta = 1/100. Error por clase:
E_{qr}(a) = psi(N; qr, a) - N/phi(qr), con psi(x; q, a) = sum_{n<=x, n=a (mod q)} Lambda(n).

Definicion (DUA). Existen A > 0 y delta > 0 tales que, uniformemente en N (par,
N -> infinito) y en la clase residual a (mod qr) con (a, qr) = 1:

  | sum_{qr en el slice} c_t(q,r) E_{qr}(a) | << M / (log N)^A,

donde la suma corre sobre los pares (q,r) del peso con 1/2 < alpha+beta <= 1/2+delta,
y M = N*M1 es el termino principal del peso en el slice.

Notas. (i) La uniformidad es sobre a: el mismo a = N aparece para todos los modulos
qr, y N varia con N; no es un a fijo. (ii) DUA es mas fuerte que Bombieri-Vinogradov
en dos sentidos: nivel 1/2+delta (BV: 1/2) y peso bilineal con factores >= N^{1/10}
(no todo q). (iii) El doc 59 cerro la dispersion estandar: la gran criba capa el
Tipo II en 1/2, Maynard cubre 0%, BFI 29/56 cubre 5.4% — y todos para a fija.

## 2. Busqueda en la literatura: que existe y que no

2.1 Bombieri-Vinogradov (1965). Nivel 1/2, suma sobre q, completamente uniforme en a
(max_a dentro de la suma). No llega a 1/2+delta. No sirve para el slice.

2.2 Linnik / dispersion de Fouvry. La dispersion de Linnik-Fouvry es el metodo que da
BV y sus variantes; Fouvry (1984-85) y Fouvry-Iwaniec dan niveles > 1/2 para
sucesiones con coeficientes, pero para a FIJA (la constante depende de a). La propia
introduccion de Maynard III resena estos trabajos como no uniformes en a: su
abstract dice que sus estimaciones son "completely uniform with respect to the
residue classes considered, unlike previous works on primes in arithmetic
progressions to large moduli". No hay version uniforme en a de la dispersion de
Fouvry en la literatura.

2.3 Bombieri-Friedlander-Iwaniec (1986-89), Teor. A, y Maynard I (arXiv:2006.06572).
Nivel hasta 4/7 (BFI, pesos bien factorizables) y 11/21 (Maynard I, q1q2 con
q1 <= x^{1/21}), pero la constante depende de a: a fija, a << x^eps. Maynard I,
Introduccion: el metodo Deshouillers-Iwaniec/Kuznetsov "requires us to only
consider the situation when a is fixed (or a small power of x)". Nuestro a = N no
es fijo ni pequeno: excluido.

2.4 Zhang / Polymath 8a (arXiv:1402.0811). Nivel 1/2+7/300 con uniformidad completa
en a, pero SOLO para modulos LISOS: p | q => p <= x^delta. Nuestros modulos son
productos de dos primos >= N^{1/10}: anti-lisos. Excluido por medicion (sec. 3):
0 modulos lisos.

2.5 Maynard III (arXiv:2006.08250v1; Mem. Amer. Math. Soc. 306(1544), 2025). El
unico articulo que da uniformidad completa en a mas alla de 1/2. Tres teoremas:

- Teor. 1.1: Q1 <= x^{1/10-3delta}/(log x)^C, Q2 <= x^{4/10+4delta}(log x)^C;
  sum_{Q1<=q1<=2Q1} sum_{Q2<=q2<=2Q2} sup_{(a,q1q2)=1} |pi(x;q1q2,a) - pi(x)/phi(q1q2)|
  <<_C delta*pi(x) + x(log log x)^2/(log x)^2. Error debil: solo O(delta) sobre la
  cota trivial (palabras del autor: "we are only winning a factor O(delta) over the
  trivial bound"). Ventana del factor pequeno: ~x^{1/10}. Excluido por medicion:
  ventana vacia (sec. 3).
- Teor. 1.2: 0 < delta < 1/1000, Q1Q2Q3 = x^{1/2+delta},
  x^{40delta} < Q2 < x^{1/20-7delta},
  x^{1/10+12delta}/Q2 < Q3 < x^{1/10-4delta}/Q2^{3/5}; buen error x/(log x)^A,
  pero uniformidad "casi": sup_b sup_{a = b (mod q1q2)} — la clase a es fija
  mod q1q2. Ademas exige TRES factores con Q3 > x^{1/10+12delta}/Q2 > 1: el split
  degenerado q3 = 1 esta excluido por la cota inferior de Q3, y nuestros modulos
  tienen Omega(d) = 2 exactamente. Excluido por partida doble.
- Teor. 1.3: equidistribucion uniforme de un MINORANTE rho (sum rho >= pi(x)/8),
  buen error, modulos q1q2 con Q1 en [x^{2/5+5delta}, x^{3/7}]. Es el unico
  teorema uniforme con buen error que toca nuestro rango. Pero: (a) es para rho,
  no para Lambda — no se aplica a nuestro peso c_t, que multiplica E_{qr} (psi,
  i.e. Lambda); (b) su ventana cubre solo el 4.8% de la masa M1 del slice (sec. 3).
  Apertura real pero despreciable.

2.6 Friedlander-Granville, "Limitations to the equidistribution of primes III",
Compositio Math. 81 (1992), 19-32. Prueba que la formula pi(x;q,a) ~ pi(x)/phi(q)
NO vale uniformemente para a fija en el rango q ~ x/log^N x. Es una limitacion
real, pero su rango esta muy por encima de x^{1/2+delta}: no es una refutacion de
DUA, y no se conoce un analogo en nuestro rango.

Conclusion de la busqueda: ningun teorema publicado da DUA. El unico candidato con
uniformidad completa mas alla de 1/2 y buen error (Maynard III Teor. 1.3) es para
un minorante y cubre una esquina del 4.8% del slice.

## 3. Medicion de cobertura (script numerico/64_dispersion_uniforme_a.py, PASS en la VM y en el PC, numeros identicos)

Familia A, t = 0.05, v = 10, slice (1/2, 1/2+1/100]. "Ventana" = pares con
max(alpha,beta) en [2/5+5delta, 3/7], delta = 1/1000 (Teor. 1.3 de Maynard III).

  h=12, E=12: min_alpha = 0.1308; M1slice/M1 = 0.43%; M1slice_win/M1slice = 0.0000 (0/9 pares); M1win/M1 = 0.09%.
  h=24, E=12: min_alpha = 0.1160; M1slice/M1 = 0.16%; M1slice_win/M1slice = 0.0000 (0/18 pares); M1win/M1 = 0.03%.
  h=96, E=12: min_alpha = 0.1041; M1slice/M1 = 0.08%; M1slice_win/M1slice = 0.0478 (2/219 pares); M1win/M1 = 0.04%.
  h=96, E=7:  min_alpha = 0.1041; M1slice/M1 = 1.00%; M1slice_win/M1slice = 0.0472 (2/219 pares); M1win/M1 = 0.91%.

- Maynard III Teor. 1.1: 0 primos representativos <= N^{1/10-3delta} (delta = 1/1000
  y 1/100); la ventana cierra en 0.097 < min_alpha = 0.1041. El script chequea la
  ventana SIN el divisor (log N)^C del teorema, i.e. un superconjunto de la ventana
  real: vacia ahi => vacia en el teorema.
- Zhang/Polymath (modulos lisos): 0 primos representativos <= N^{1/1000} (y
  <= N^{1/100}); min(q)^1000 > N exacto.
- Omega(d) = 2 para todo par (q_i, q_j): el split en tres factores del Teor. 1.2
  no se puede alimentar.
- Celda balanceada de Vaughan: 2^19 = 524288 <= 10^6 < 2^20 = 1048576, dentro de
  (10^4, 10^8]: la celda M ~ K ~ x^{1/2} existe en el Tipo II de Vaughan
  (U = V = x^{1/3}).

Lectura: la esquina del Teor. 1.3 es delgada pero real: con el slice (1/2, 1/2+1/100]
y el soporte min alpha >= 1/10, la ventana nominal [x^{0.405}, x^{0.4286}] se corta
a la efectiva [x^{0.405}, x^{0.41}] (el otro factor queda en [x^{0.1}, x^{0.105}]);
a h = 12 y 24 la malla es demasiado gruesa para resolverla (0 pares), a h = 96 se
resuelve con 2 pares y 4.78% de la masa M1 del slice, estable entre E = 7 y E = 12
(4.72% vs 4.78%). Pero es el 4.8% de un slice que es el 0.08% de M1: el 0.004%
de M1. Insuficiente para una cota inferior de capa por dos ordenes de magnitud,
antes siquiera de cambiar Lambda por el minorante.

## 4. Obstruccion: por que ningun metodo conocido puede dar DUA

Prop. 64.1 (los dos ingredientes que llegan mas alla de 1/2 no dan uniformidad en
a). La gran criba — el ingrediente de la dispersion — da BV hasta el nivel 1/2
con uniformidad completa en a, pero su termino (N + Q^2) capa el metodo en 1/2:
con Q^2 > N el error de dispersion supera al termino principal. El ingrediente
que rompe 1/2 es el metodo espectral (Deshouillers-Iwaniec/Kuznetsov: sumas de
Kloosterman), y ese llega mas alla de 1/2 solo para a fija (BFI Teor. A; Maynard I,
Introduccion, cita en sec. 2.3). Luego ningun metodo conocido produce uniformidad
en a a nivel > 1/2 para sumas bilineales: no es un reordenamiento pendiente, es un
limite de los dos ingredientes disponibles.

Prop. 64.2 (la celda balanceada es exactamente el caso imposible). Maynard III,
sec. 5 (Props. 5.1-5.2): las convoluciones desbalanceadas (un factor en
[x^{2/5+eps}, x^{1/2-eps}]) se tratan con uniformidad completa; las balanceadas
M ~ K ~ x^{1/2} NO: "we are not able to produce an estimate which is completely
uniform for such terms". Nuestro slice (alpha+beta en (1/2, 1/2+1/100], alpha,
beta >= 1/10) contiene la diagonal balanceada alpha = beta ~ 1/4 y la celda de
Vaughan M ~ K ~ x^{1/2} (sec. 3). Es decir: el corazon del slice es precisamente
el caso que el unico articulo uniforme declara fuera de alcance. Esto no es una
refutacion de DUA (el articulo no prueba imposibilidad, solo que su metodo no
llega), pero es la evidencia mas fuerte disponible de que la uniformidad completa
en la celda balanceada exige una idea nueva.

Prop. 64.3 (especializacion a un solo par: DUA general implicaria una cota
individual nueva). Si DUA valiera para TODO peso bilineal alpha_u beta_v con
u,v >= N^{1/10} en el slice, tomese alpha = delta_{u0}, beta = delta_{v0} con
u0 v0 = q ~ N^{1/2+delta/2}, u0, v0 >= N^{1/10} (posible: u0 ~ N^{1/4},
v0 ~ N^{1/4+delta/2}). Entonces DUA daria
|psi(N;q,N) - N/phi(q)| << N/(phi(q)(log N)^A) uniformemente en a, para q hasta
N^{1/2+delta/2}. Esa cota individual uniforme en a a nivel > 1/2 NO se conoce: la
region libre de ceros clasica la da solo para q <= exp(c (log N)^{3/5}); la gran
criba la da en media sobre a, no para a individual. Una cota asi exigiria una
estimacion de densidad de ceros nueva (del tipo N(sigma, T) con sigma < 1 en el
rango necesario). Luego DUA general es AL MENOS tan dificil como un teorema nuevo
de densidad de ceros. (Nota: esto es para DUA general; para nuestro c_t
especifico la especializacion no aplica directamente, porque c_t es un peso suave
fijo, no un delta. Para c_t, la obstruccion es la Prop. 64.2 + la medicion de la
sec. 3.)

Prop. 64.4 (refutacion: no disponible). Una refutacion de DUA exigiria exhibir a
y un peso con error >> M/(log N)^A. La unica via conocida seria un cero de Siegel
(caracter excepcional: el error seria ~ N^beta/phi(q) con beta cercano a 1). No se
conoce ningun cero de Siegel, y su no-existencia tampoco esta probada. No hay
refutacion incondicional disponible. Friedlander-Granville (sec. 2.6) muestra que
la uniformidad PUEDE fallar, pero en un rango muy superior al nuestro.

## 5. Veredicto

DUA NO EXISTE en la literatura en la forma que necesitamos, y NO SE PUEDE PROBAR
con los metodos conocidos:

- La gran criba capa en 1/2 (Prop. 64.1).
- El metodo espectral (BFI/Maynard I-II) llega mas alla de 1/2 solo para a fija.
- El unico resultado uniforme mas alla de 1/2 con buen error (Maynard III
  Teor. 1.3) es para un minorante y cubre el 4.8% de la masa M1 del slice
  (0.004% de M1).
- La celda balanceada que necesitamos es exactamente la que Maynard III sec. 5
  declara fuera de alcance con uniformidad completa (Prop. 64.2).
- DUA general implicaria una cota individual uniforme en a a nivel > 1/2, que
  exigiria densidad de ceros nueva (Prop. 64.3).
- No hay refutacion incondicional: la unica via seria un cero de Siegel
  (Prop. 64.4).

En una frase: la dispersion uniforme en a que necesitamos no existe, no se puede
probar con lo conocido, y su prueba exigiria o bien una idea nueva para la celda
balanceada, o bien una estimacion de densidad de ceros nueva; su refutacion
exigiria un cero de Siegel.

## 6. Que sigue: el lema faltante, delimitado

Lema faltante (DUA-esquina). Una estimacion de dispersion uniforme en a para pesos
bilineales con factores en la esquina (alpha, beta) en [1/10, 1/10+1/200] x
[2/5+5delta, 2/5+5delta+1/100] (y su simetrica), con error << M/(log N)^A — o bien
una version del Teor. 1.3 de Maynard III para Lambda (no para el minorante rho) en
esa esquina, con el mismo error. Si ese lema existiera, cubriria el 4.8% del
slice; haria falta ademas extenderlo fuera de la esquina, lo cual es el problema
completo. La esquina sola no basta para una cota inferior de capa (es el 0.004%
de M1; el doc 58 ya mostro que el slice entero es el 0.08% de M1 y no basta).

Frentes que quedan (sin duplicar el bucle): (a) mecanismo de seleccion (switching,
C10/W(d)); (b) Motohashi doc 57 sec. 4; (c) uniones C_k + C_{k+1}; (d) variante
del Lema 4 de Pintz-Althofer con sup en sigma > 5/7. La dispersion uniforme en a
queda CERRADA como ruta: no es el camino.

## Referencias (leidas en esta sesion, salvo indicacion)

- J. Maynard, "Primes in arithmetic progressions to large moduli III: Uniform
  residue classes", arXiv:2006.08250v1 (2020); Mem. Amer. Math. Soc. 306(1544),
  2025. Teors. 1.1, 1.2, 1.3; sec. 5 (Props. 5.1-5.2). [Leido: version arXiv v1,
  introduccion y enunciados; Teors. 1.1 y 1.2 verificados textualmente contra la
  copia de Oxford (ora.ox.ac.uk) el 2026-11-25.]
- J. Maynard, "Primes in arithmetic progressions to large moduli I: Fixed residue
  classes", arXiv:2006.06572. Introduccion (cita sobre a fija). [Leido.]
- J. Maynard, "Primes in arithmetic progressions to large moduli II:
  Well-factorable estimates", arXiv:2006.07088v1. [Citado en Maynard III.]
- Y. Zhang, "Bounded gaps between primes", Ann. of Math. 179 (2014); D. H. J.
  Polymath, "Variants of the Selberg sieve, and bounded intervals containing many
  primes", Res. Math. Sci. 1 (2014), art. 12 (arXiv:1402.0811). Nivel 1/2+7/300,
  modulos lisos. [Leido: enunciados y nivel.]
- E. Bombieri, J. Friedlander, H. Iwaniec, "Primes in arithmetic progressions to
  large moduli", Acta Math. 156 (1986); ibid. II, Math. Ann. 277 (1987); ibid.
  III, Compositio Math. 68 (1989). Teor. A (a fija). [Citado en Maynard I/III.]
- E. Fouvry, "Autour du theoreme de Bombieri-Vinogradov", Acta Math. 152 (1984);
  "Repartition des suites dans les progressions arithmetiques", Ann. Inst.
  Fourier 35 (1985). Dispersion, a fija. [Citado en Maynard III, Introd.]
- J. Friedlander, A. Granville, "Limitations to the equidistribution of primes
  III", Compositio Math. 81 (1992), 19-32. [Leido: abstract y teorema.]
- Docs 58, 59, 60, 63 de esta carpeta; script numerico/64_dispersion_uniforme_a.py
  (PASS en la VM y en el PC, numeros identicos).
