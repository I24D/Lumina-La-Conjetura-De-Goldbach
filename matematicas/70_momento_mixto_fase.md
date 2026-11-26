# 70 - Lema del momento mixto, via conjunta con fase: la circularidad (Etapa 2)

Fecha: 2026-11-26. Etapa 2 de 2 (solo esta etapa): el lema del momento
mixto del doc 67 sec. 3, atacado por UNA sola via conjunta: el termino
cruzado CON FASE. Script: `numerico/70_momento_mixto_fase.py` — PASS
26/26 local (X = 2x10^6). No hay prueba del lema, no hay cota inferior
de capas, no hay demostracion de Goldbach.

## 1. El lema y la identidad exacta

Sea N par, S_1(alfa) = sum_{p<N} e(p*alfa), S_rho(alfa) = sum_n rho(n)
e(n*alfa) con rho el minorante de Harman (doc 66, diadico a escala N),
rho = 1_p - D con D >= 0 el defecto (proxy rho_tilde = 1_p - (S_{2,3} +
S_{2,5}); el rho verdadero anade la parte B, tambien >= 0, doc 66).
Arcos mayores M y menores m como en el doc 69 sec. 1.

**Lema MM-69.1 (momento mixto, absoluto).** int_m |S_1(alfa)
S_rho(alfa)| d alfa = o(N/log^2 N).

El metodo del circulo para el conteo unilateral U(N) = sum_{p<N}
rho(N-p) necesita solo la version CON FASE, mas debil:

  E_U := int_m S_1(alfa) S_rho(alfa) e(-N alfa) d alfa = o(N/log^2 N),

pues |E_U| <= int_m |S_1 S_rho|. La etapa 1 (doc 69) cerro las tres
estimaciones por piezas del lema absoluto. Esta etapa muestra que ni
siquiera la version con fase se puede obtener por esta via: es
circular.

**Identidad exacta.** Por ortogonalidad, con S_rho = S_1 - S_D
(S_D(alfa) = sum_n D(n) e(n*alfa)):

  U(N) = int_0^1 S_1 S_rho e(-N alfa) = R(N) - D(N),
  R(N) = int_0^1 S_1^2 e(-N alfa),
  D(N) := sum_{p<N} D(N-p) = int_0^1 S_1 S_D e(-N alfa).

Termino a termino, rho = 1_p - D, luego U(N) = R(N) - D(N) es exacto
en enteros (verificado en el script, sec. 1). Partiendo [0,1] = M u m:

  E_U = E_R - E_D,
  E_R := int_m S_1^2 e(-N alfa),   E_D := int_m S_1 S_D e(-N alfa).

**Prop. 70.1 (circularidad) [T].** E_R = o(N/log^2 N) es EQUIVALENTE
a la asintotica de Goldbach R(N) = S(N) N/log^2 N + o(N/log^2 N),
pues la parte de arcos mayores de R(N) es el termino principal mas
O(N/(log N)^A) (estandar, caps. 2-3 de Vaughan). Del mismo modo, la
parte de arcos mayores de U(N) es c_rho*S(N) N/log^2 N + O(N/(log
N)^A) (doc 67, Prop. 67.1), luego U(N) = c_rho*S(N) N/log^2 N +
o(N/log^2 N) si y solo si E_U = o(N/log^2 N). Por la identidad,
E_U = E_R - E_D: probar E_U = o(N/log^2 N) por esta descomposicion
exige, como ingrediente esencial, o bien E_R = o(N/log^2 N) — es
decir, la asintotica de Goldbach misma — o bien E_D = o(N/log^2 N),
o bien una cancelacion E_D ~ E_R. Ninguna de las tres esta
disponible (sec. 2-3). La via con fase no circunvala la barrera de
arcos menores: la HEREDA exactamente. La barrera de la linea rho ES
la barrera de Goldbach.

## 2. El termino cruzado E_D: tan dificil como el problema [T]

E_D es el error de arcos menores del conteo D(N) = sum_{p<N} D(N-p),
las representaciones N = p + n con n en el soporte del defecto. No es
una correccion pequena: D(N)/R(N) = 0.191, 0.247, 0.301, 0.331, 0.352
a N = 2x10^5, 5x10^5, 10^6, 1.5x10^6, 2x10^6 (script, exacto),
creciendo hacia L_S ~ 0.325-0.535 (docs 66-67). El defecto es un
problema binario de tamano comparable a R(N), no un residuo.

La unica cota disponible para E_D es Cauchy:

  |E_D| <= ||S_1||_2 ||S_D||_2 = sqrt(pi(N) * sum D^2)
        = sqrt(0.467) * pi(N) ~ 0.683 * N/log N,

medido: 101,754 frente al objetivo N/log^2 N = 9,501 a X = 2x10^6
(cociente 10.7); asintoticamente 0.683*log N -> infinito (18.9 a
N = 10^12, 37.7 a N = 10^24). Un logaritmo por encima del objetivo.
E_D no se puede mostrar pequeno con lo conocido, y es un integral de
arcos menores de la misma clase que E_R. La descomposicion no separa
el problema: lo duplica.

Tampoco hay mecanismo de cancelacion E_D ~ E_R: son los errores de
arcos menores de dos conteos distintos (primos frente a enteros con
defecto), y el unico control disponible de E_D (Cauchy) esta un log
por encima del objetivo. Afirmar la cancelacion seria afirmar, sin
prueba, que dos errores desconocidos coinciden.

## 3. Verificacion en la malla FFT (identidad, no magnitud)

En la malla M = 2^21 puntos a_j = j/M, con arcos mayores Q_0 =
(log X)^2 = 210, Q_1 = X/(log X)^2 = 9501 (fraccion de malla en
arcos mayores 0.0267), la suma de Riemann verifica la identidad
algebraica E_U = E_R - E_D exactamente (error 1.4e-14, script):
es linealidad, y valida las convenciones de signo de la FFT
(el circulo completo reproduce U(N) = 12596 exacto, error 2e-06).
Las magnitudes en la malla (|E_R| = 288.9, |E_D| = 75.6,
|E_U| = 216.2, objetivo 9501) son sumas de Riemann a escala finita,
NO mediciones de las integrales verdaderas: a esta escala los arcos
mayores ya capturan casi todo y los errores son pequenos, lo cual es
consistente con (no prueba de) convergencia lenta O(1/log N). El
contenido riguroso es la identidad y la Prop. 70.1, no las
magnitudes finitas.

## 4. Las otras dos vias conjuntas, descartadas (no intentadas aqui)

(a) Gran criba del producto. La gran criba es una herramienta L^2
para sumas lineales; Parseval fija las normas L^2, y solo gana
cuando Q^2 << N, pero los arcos menores llegan a q ~ N. Aplicada a
la convolucion c = 1_p * rho, la masa es sum_n U(n)^2 = 1.7x10^14,
cociente 0.943 con X^3/log^4 X (script): la masa es de orden X^3,
no X^2. La gran criba reproduce la escala N/log N de Cauchy, no la
baja: para pasar de L^2 sobre racionales a L^1 sobre m habria que
pagar el mismo factor de Cauchy. No mejora nada.

(b) Input de densidad de ceros / GRH. GRH es teorema en este
proyecto (docs GB-L05, Riemann), pero no rescata ninguna ruta:
sup x L^1 falla ya con el sup IDEAL N^{1/2} (Prop. 69.2(b): N/(log
N)^{1/2}, 145x el objetivo a N = 10^12, 411x a N = 10^24); Cauchy es
independiente de GRH; el cuarto momento via sup^2 x L^2, aun con
sup ideal, queda (log N)^3 por encima de lo necesario (21,096x a
N = 10^12, 168,764x a N = 10^24, script). Un ataque directo por
formula explicita necesitaria correlacion de pares de ceros
(Montgomery, abierto). Sub-lema delimitado SL-70.1: una cota de
correlacion de pares de ceros suficiente para el momento mixto.

## 5. Veredicto y que sigue

**Veredicto etapa 2: la via conjunta con fase NO da el lema del
momento mixto, y no puede darlo sin resolver antes los arcos
menores de Goldbach (Prop. 70.1) o el termino cruzado E_D, de la
misma dificultad (sec. 2).** La gran criba del producto reproduce
Cauchy (sec. 4a) y GRH no rescata ninguna ruta (sec. 4b). El lema
MM-69.1 queda como cuello de botella abierto de la linea rho
(docs 66-67): la linea rho elimina la obstruccion de
equidistribucion pero hereda intacta la obstruccion de arcos
menores, que es la de Goldbach mismo. Goldbach sigue abierta.

Que sigue: los frentes del bucle sin tocar (seleccion: switching,
C10/W(d); Motohashi doc 57 sec. 4; b_{P,l} doc 53; uniones
C_k+C_{k+1}; variante Lema 4 Pintz-Althofer; SL-65.1: agudeza del
umbral 1/10); o una idea genuinamente nueva para el momento mixto
que no sea descomposicion, criba ni GRH — por ejemplo, el sub-lema
SL-70.1 (correlacion de pares de ceros). Ninguna de esas vias se
duplica aqui.

Referencias (leidas, no de memoria):
- Docs 66 (rho, c_rho), 67 (conteo unilateral, Prop. 67.1/67.3),
  69 (etapa 1: Props. 69.1-69.3), 64-65 (dispersion, DUA-esquina).
- E. C. Titchmarsh / R. C. Vaughan, "The Hardy-Littlewood Method",
  2nd ed., caps. 2-3: metodo del circulo binario, arcos mayores,
  cota de Vinogradov.
- J. Maynard, arXiv:2006.08250v1, sec. 9 (rho) y Teor. 1.3.
- GB-L05 (docs GB-L05-major-arcs-GRH): GRH como teorema formal
  para los arcos mayores.
