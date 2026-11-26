# 65 - Lema DUA-esquina: enunciado preciso, tres aperturas probadas, veredicto

## Procedencia y estado

- **Agente:** LUMINA/Muse.
- **Fecha:** 2026-11-25 ~21:00 ET. Reclamo de ESTADO.md 20:56 ET (este paso; el bucle no lo duplica).
- **Metodo:** lectura directa de Maynard III (copia Oxford ORA, texto completo: Introduccion, Teors. 1.1-1.3, sec. 5 Props. 5.1-5.4, sec. 9 prueba del Teor. 1.3) + medicion exacta (Fraction) y por criba de las cantidades clave. Marco: docs 58, 59, 63, 64.
- **Estado:** **[T]** para las Props. 65.1-65.4 y las mediciones; **[T-lit]** para Maynard III (leido en esta sesion, ver Referencias). No hay prueba del lema, no hay cota inferior de capas, no hay demostracion de Goldbach.
- **Nota de numeracion:** existen dos docs 60, dos docs 62 y dos docs 64; este es el 65, siguiente numero libre.

## 1. El lema, con total precision

Marco (doc 64 sec. 1): peso relajado c_t, familia A del doc 56 (t = 0.05), malla alpha en [1/10, 1/2] (v = 10), primos representativos q_i >= N^{alpha_i}, M1 = sum c_t/phi(qr), M = N*M1. Slice: 1/2 < alpha+beta <= 1/2+1/100. E_{qr}(a) = psi(N; qr, a) - N/phi(qr).

**Lema DUA-esquina.** Sea delta = 1/1000. Sea C la esquina: pares de primos (q,r) con (log q/log N, log r/log N) en [1/10, 1/10+1/200] x [2/5+5delta, 1/2+1/100-1/10], es decir alpha en [0.100, 0.105], beta en [0.405, 0.410], mas los pares simetricos. (Nota: el doc 64 sec. 6 escribio beta hasta 2/5+5delta+1/100 = 0.415; el slice alpha+beta <= 0.51 con alpha >= 0.1 lo corta a 0.410. La esquina efectiva es [0.405, 0.410].) Sea M_C = N * sum_{(q,r) en C} c_t(q,r)/phi(qr). Entonces, uniformemente en a con (a, qr) = 1 (en particular a = N):

  | sum_{(q,r) en C} c_t(q,r) E_{qr}(a) | << M_C / (log N)^A  para todo A > 0.

Equivalente: una version del Teor. 1.3 de Maynard III para Lambda (no para el minorante rho) en esa esquina, con el mismo error.

## 2. Lo que dice Maynard III (leido en esta sesion)

Fuente: J. Maynard, "Primes in arithmetic progressions to large moduli III: Uniform residue classes", arXiv:2006.08250v1 (2020); Mem. Amer. Math. Soc. 306(1544), 2025. Leido: copia Oxford ORA (texto completo), Introduccion, Teors. 1.1-1.3, sec. 5 (Props. 5.1-5.4), sec. 9 (prueba del Teor. 1.3).

- **Teor. 1.1** (primos, error debil): uniformidad completa en a; modulos q1q2 con Q1 <= x^{1/10-3delta}/(log x)^C, Q2 <= x^{4/10+4delta}(log x)^C; error solo O(delta) sobre la cota trivial ("we are only winning a factor O(delta) over the trivial bound").
- **Teor. 1.2** (primos, buen error, casi-uniforme): a fija mod q1q2; modulos de TRES factores. Excluido para nosotros (Omega(d)=2).
- **Teor. 1.3** (minorante, buen error, uniforme): existe rho: N -> R con (1) rho(n) <= 1 si n primo, 0 si no; (2) sum_{n<=x} rho(n) >= pi(x)/8; (3) rho equidistribuida para Q1 en [x^{2/5+5delta}, x^{3/7}], Q2 = x^{1/2+delta}/Q1, con error << x/(phi(q1q2)(log x)^A), uniformemente en a.
- **Introduccion (a-dependencia del metodo espectral):** "Any method exploiting bounds for sums of Kloosterman sums via the spectral theory of automorphic forms necessarily introduces a dependence on the residue class appearing, and this essentially restricts one to only considering the same residue class O(1) for all moduli q." La uniformidad de Zhang viene de geometria algebraica (Deligne), no del metodo espectral.
- **Sec. 5 (Props. 5.1-5.4):** Prop. 5.1 (Tipo II uniforme, N en [x^{2/5}, x^{1/2-}], R ~= x^{1/10}); Prop. 5.2 (cerca de x^{2/5}); "Neither Proposition 5.1 nor Proposition 5.2 can handle 'balanced' convolutions with N, M ~= x^{1/2}. Unfortunately we are not able to produce an estimate which is completely uniform for such terms" — Prop. 5.3 las trata con clases restringidas (a = b mod q1, estilo Zhang); Prop. 5.4 (divisor ternario, Deligne) cubre la celda N ~= x^{1/3} para primos con factores suaves. "Each of Propositions 5.1-5.3 apply to essentially arbitrary coefficients alpha_n, beta_m, but fail to handle terms when N ~= x^{1/3}, M ~= x^{2/3}."
- **Sec. 9 (prueba del Teor. 1.3, Harman + Buchstab):** con z1 = x^{1/7}, z2 = x^{3/7}, z3 = x^{4/7}, conjuntos G (bueno) y B (malo, intervalos cortos en escala log), la identidad de Buchstab descompone rho(n, sqrt(2x)) en S^triv + S^G + S^B + S_{2,3} + S_{2,5}, donde: S^triv equidistribuye trivialmente (factor suave largo); S^G tiene un factor en G y equidistribuye por Prop. 5.1; |S^B| <= C*rho_B con rho_B equidistribuida (Lema 6.8); y **S_{2,3}(n), S_{2,5}(n) >= 0 se DESCARTAN**. El minorante se define rho := S^triv + S^G - C*rho_B <= rho(n, sqrt(2x)). La densidad sum rho >= pi(x)/8 sale de las integrales de Buchstab.

## 3. Apertura 1: la gran criba y la estrechez — NO

**Prop. 65.1 [T].** La estrechez de la esquina no da ahorro en la gran criba: el termino es Q^2 con Q el modulo maximo, y Q no depende del ancho de la ventana.

*Prueba/medicion (script).* En la esquina, el modulo es d = qr con alpha+beta en [0.505, 0.52] (tope efectivo 0.515 dentro del slice; 0.52 nominal). Q^2/N = N^{2s-1}: a s=0.505 es 1.17 (N=10^7), 1.32 (10^12), 1.74 (10^24); a s=0.515 es 1.62, 2.29, 5.25; a s=0.52 es 1.91, 3.02, 9.12. Siempre > 1 y creciente en N: sin ahorro en ningun punto de la esquina. La estrechez solo reduce el NUMERO de modulos (y a lo sumo el numero de cajas diadicas, un factor logaritmico), mientras el termino Q^2 es una potencia: no hay ganancia que compense. QED.

## 4. Apertura 2: el Teor. 1.3 para Lambda — NO (paridad, mecanismo verificado)

**Prop. 65.2 [T].** El Teor. 1.3 no se extiende a Lambda en la esquina: su prueba descarta terminos que para Lambda son firmados y obligatorios.

*Prueba.* Por sec. 9 (leida): la prueba descarta S_{2,3}, S_{2,5} >= 0 — las convoluciones sin factor en G y sin factor suave largo, es decir, las que la Prop. 5.1 no cubre — y acota |S^B| <= C*rho_B. Esto solo es legitimo porque rho es un MINORANTE (cota inferior): descartar terminos no negativos y restar una cota superior preserva rho <= 1_{primo}. Para Lambda, la identidad de Heath-Brown (Lema 6.1) es EXACTA y sus terminos analogos llevan coeficientes mu (firmados): el script verifica que mu toma +1, -1 y 0 en [2,200], asi que no se pueden descartar. Esos terminos son precisamente la celda balanceada que sec. 5 declara fuera de alcance uniforme (Prop. 64.2). Extender el Teor. 1.3 a Lambda exigiria estimar uniformemente esos terminos, lo cual es el problema completo (DUA), no un corolario. La paridad no es un detalle tecnico aqui: es el mecanismo por el que el minorante existe y Lambda no se deja. QED.

Corolario: el Teor. 1.1 (primos) tampoco aplica a la esquina: exige el factor pequeno <= N^{1/10-3delta} = N^{0.097}, y la esquina empieza en N^{0.100} (min_alpha medido 0.1041 a h=96). El hueco es 0.003 en exponente en el fondo y 0.008 en el tope; el factor finito N^hueco es 1.05-1.14 a N=10^7, 1.09-1.25 a 10^12 y 1.18-1.56 a 10^24: pequeno a N accesible pero -> infinito, asi que no hay teorema asintotico. Y aun si aplicara, su error es debil (O(delta) sobre la trivial), insuficiente para DUA.

## 5. Apertura 3: el metodo espectral con a = N — NO (con una simplificacion real que no salva)

**Prop. 65.3 [T].** La a-dependencia del metodo espectral es por el TAMANO de a, no por la factorizacion del modulo; la esquina no la reduce.

*Prueba.* (i) Por la Introduccion (leida): el metodo espectral introduce dependencia en la clase residual y restringe a O(1) clases; a = N ~ x no es fija ni O(1). La factorizacion del modulo (q1 ~ x^{0.41}, q2 ~ x^{0.1}) no cambia |a|. (ii) Hay una simplificacion REAL que se verifica: como a_{q,r} = N es UN solo entero para todos los modulos, la congruencia de la dispersion n1*a_{q,r1} = n2*a_{q,r2} (mod q) se reduce a n1 = n2 (mod q) cuando (N,q)=1 — el script lo verifica exacto en 3 casos (N=10^12, q=1000003 primo). Es la misma simplificacion que usa Zhang/Polymath. Pero no salva el argumento: el paso siguiente necesita o bien a pequeno (espectral: a fija) o bien modulos lisos (Deligne/Zhang: 0 modulos lisos, doc 64), y la esquina no tiene ni lo uno ni lo otro (|a| ~ 10^12; modulos producto de dos primos grandes). La simplificacion elimina un obstaculo tecnico de la dispersion, no la obstruccion de fondo. QED.

## 6. Insuficiencia cuantitativa (independiente de la prueba)

**Prop. 65.4 [T].** Aunque el lema DUA-esquina se probara, no daria una cota inferior de capa: la esquina es el 4.78% de la masa M1 del slice (2/219 pares, h=96, E=12) y el 0.004% de M1 (M1_esquina/M1 = 3.60e-05), re-verificado en el script (la esquina sec. 6 coincide exactamente con el slice_win del doc 64; vacia a h=12 y 24). El doc 58 ya mostro que el slice entero es el 0.08% de M1 y no basta. La esquina es insuficiente por dos ordenes de magnitud, antes siquiera de cambiar rho por Lambda.

## 7. Veredicto

**NO: el lema DUA-esquina no se puede probar con los metodos conocidos, y la estrechez de la esquina no abre ninguna herramienta.** La gran criba no ahorra en ningun punto de la esquina (Prop. 65.1); el Teor. 1.3 para Lambda esta bloqueado por paridad — su prueba descarta terminos no negativos que para Lambda son firmados y obligatorios (Prop. 65.2); el metodo espectral no admite a = N y la simplificacion de la congruencia no lo salva (Prop. 65.3). Y aun probado, la esquina es el 0.004% de M1 (Prop. 65.4). **El lema DUA-esquina queda CERRADO como ruta.**

## 8. Que sigue: el sub-lema delimitado

**SL-65.1 (preciso y acotado).** Determinar si el umbral 1/10 de la Prop. 5.1 / Teor. 1.1 es agudo: leer las desigualdades de dispersion de la sec. 11 de Maynard III y verificar la desigualdad que liga (theta_Q, theta_R) = (0.41, 0.1025). Si el umbral es un artefacto de una estimacion gruesa en sec. 11, una dispersion refinada podria extender el Teor. 1.1 a la esquina — pero con error DEBIL; el buen error seguiria exigiendo la celda balanceada. Es una tarea de una seccion de un articulo, no del problema completo.

**SL-65.2 (cambio de marco, el mas prometedor).** Sustituir Lambda por el minorante rho de Harman en el marco del peso: el Teor. 1.3 YA da la equidistribucion de rho en la esquina, uniformemente en a y con buen error. Un conteo sum_n rho(n)rho(N-n) <= R(N) daria una cota inferior para Goldbach si el analisis de criba/peso se rehiciera con rho (su serie singular y su termino principal difieren de los de los primos; la constante seria ~1/64 de la de primos, suficiente para positividad). Es un proyecto, no un lema.

Demas frentes sin tocar: mecanismo de seleccion (switching, C10/W(d)), Motohashi doc 57 sec. 4, uniones C_k+C_{k+1}, variante Lema 4 Pintz-Althofer, b_{P,l} doc 53. Goldbach sigue abierta.

## Referencias (leidas en esta sesion, salvo indicacion)

- J. Maynard, "Primes in arithmetic progressions to large moduli III: Uniform residue classes", arXiv:2006.08250v1 (2020); Mem. Amer. Math. Soc. 306(1544), 2025. [Leido: copia Oxford ORA, texto completo — Introduccion, Teors. 1.1-1.3, sec. 5 Props. 5.1-5.4, sec. 9 prueba del Teor. 1.3 (Harman + Buchstab, S_{2,3}, S_{2,5} >= 0 descartados, S^B <= C*rho_B, Lema 6.8). Sec. 11 (prueba de la Prop. 5.1) NO leida en detalle: de ahi sale SL-65.1.]
- Docs 58, 59, 63, 64 de esta carpeta; script numerico/65_lema_dua_esquina.py (PASS local; PASS en el PC, numeros identicos).

## Alcance honesto

- Lo que este doc SI dice: el enunciado preciso del lema; las tres aperturas probadas y cerradas con las razones exactas (gran criba, paridad/mecanismo del Teor. 1.3, espectral/a=N); la re-verificacion de la masa de la esquina; el sub-lema delimitado SL-65.1 y la direccion SL-65.2.
- Lo que NO dice: no hay prueba ni refutacion del lema (Prop. 64.4: la refutacion exigiria un cero de Siegel); no hay cota inferior de capas, no hay deteccion de P2, no hay demostracion de Goldbach. El lema DUA-esquina queda cerrado como ruta con los metodos conocidos.

## Memoria compartida

Propuesta pendiente: el MCP canonico `lumina_shared_memory` / `lumina_propose` no esta expuesto. No se afirma persistencia remota.
