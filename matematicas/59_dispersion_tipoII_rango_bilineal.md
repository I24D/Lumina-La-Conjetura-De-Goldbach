# 59 - Dispersion / Tipo II para el rango bilineal: la obstruccion del factor pequeno

## Procedencia y estado

- **Agente:** LUMINA/Muse.
- **Fecha:** 2026-11-25 ~19:30 ET. Reclamo de ESTADO.md 19:23 ET (este paso; el bucle no lo duplica).
- **Metodo:** auditoria de las herramientas publicadas de dispersion (fuentes leidas, no de memoria) + medicion exacta con aritmetica racional (`Fraction`) de la cobertura de cada herramienta sobre el peso relajado.
- **Estado:** **[T]** para la reduccion por reordenamiento (Prop. 59.1), la aritmetica de la esquina BFI, la vacuidad de la ventana de Maynard y las mediciones; **[T-lit]** para los teoremas publicados citados; **[P-rev]** para la obstruccion de la clase fija. No hay cota inferior de capas, no hay nivel > 1/2, no hay demostracion de Goldbach.

## 1. El objetivo preciso: que tendria que dar la dispersion

El peso relajado `c_t` (familia A del doc 56) tiene nivel 1/2 asintotico (doc 58): el error pesado con `qr <= N^(1/2)` es despreciable (BV, `R_BV(10) ~ 1e-9`) y la cola por encima se controla con la cota trivial (`R_triv = 0.094` a `N=10^12`), pero eso no es un teorema de nivel. Por el doc 38, Prop. 38.2, con nivel 1/2 el tamiz lineal da `f(s) = 0` para las capas en cuestion: el nivel 1/2 es exactamente el nivel que no alcanza. Una estimacion de dispersion / Tipo II tendria que dar, para algun `delta > 0`, que el error pesado sobre el slice

```text
S(delta) = {(q,r): N^(1/2) < qr <= N^(1/2+delta)}
```

es `<< M/(log N)^A` para todo `A > 0` (Def. 58.1, `M = N*M1`), es decir, un ahorro de potencia de logaritmo sobre la cota trivial, no solo una fraccion constante de `M`. La pregunta es si alguna herramienta publicada lo da para este peso, cuyo soporte es `q,r` primos en `[N^(1/10), N^(1/2)]` (doc 58, sec. 2: `c_t` es suma de 3 terminos de rango 1 con factores soportados en ese intervalo).

## 2. El tope de la gran criba: por que el Tipo II se detiene en 1/2

**Desigualdad de la gran criba (Davenport-Halberstam) [T-lit].** En la forma aditiva,

```text
sum_{q<=Q} sum*_{(a,q)=1} |sum_{M<n<=M+N} a_n e(an/q)|^2
    << (Q^2 + N) sum_{M<n<=M+N} |a_n|^2,
```

ver Kowalski, "The large sieve and L-functions over finite fields" (notas, Torino), forma trigonometrica y forma multiplicativa; y Zhao, "Large Sieves and Exponential Sums" (tesis, dir. Iwaniec), Cor. 1. El termino `Q^2` es el tope: la gran criba solo ahorra cuando `Q^2 << N`, es decir, `Q << N^(1/2)`.

Es el mecanismo que limita la estimacion Tipo II / bilineal en la prueba de Bombieri-Vinogradov: la suma bilineal `sum_{m~M} sum_{n~N} a_m b_n 1_{mn = a (q)}` se controla con la gran criba aplicada a la sucesion `a*b`, y el termino `Q^2` exige `Q <= x^(1/2)` (salvo logaritmos). La dispersion de Linnik tiene el mismo tope: es la herramienta que DA el nivel 1/2, no la que lo supera. Superarlo exige factorizar el modulo con un factor pequeno (sec. 3): ese es el contenido de BFI y de Maynard, y es exactamente lo que este peso no tiene en el grueso de su masa.

## 3. Las herramientas mas alla de 1/2 y sus condiciones exactas

Tres herramientas publicadas superan 1/2, y las tres exigen que el modulo tenga un factor pequeno (o un peso bien factorizable). Se citan como se leyeron.

**Maynard, Teorema 1.1 y Corolario 1.2 [T-lit].** J. Maynard, "Primes in arithmetic progressions to large moduli I: Fixed residue classes", Teor. 1.1 y Cor. 1.2 (leido): para `a` en Z FIJO y `0 < delta < 1/42`, los modulos `q <= x^(1/2+delta)` que tienen un divisor en el intervalo

```text
[x^(2 delta), x^(1/10 - 7 delta/5)]
```

satisfacen `sum |pi(x;q,a) - pi(x)/phi(q)| << x/(log x)^A` con valores absolutos. Cor. 1.4: `q1 <= x^(1/21)`, `q2 <= x^(10/21-epsilon)`, nivel `11/21`. La propia glosa del articulo: cuando `delta` es pequeno, "la mayoria" de los modulos de tamano `x^(1/2+delta)` tienen un divisor en ese rango. La sec. 15 ("Fouvry-style estimates near x^(1/7)") es la version estilo Fouvry, con un factor cerca de `x^(1/7)`.

**BFI, exponente 29/56 para factorizacion fija [T-lit].** Bombieri-Friedlander-Iwaniec [4], segun se resume en "On primes in arithmetic progressions and bounded gaps between many primes" (ScienceDirect, Introduccion): exponente `29/56` para coeficientes soportados en modulos de la forma `q1 q2` con `q1 <= x^th1`, `q2 <= x^th2`, `th1 < 1/3`, `th2 < 1/5` y `5 th1 + 2 th2 < 2`. La esquina optima es `(th1, th2) = (9/28, 11/56)`: `5*(9/28)+2*(11/56) = 2` (frontera, verificada exacta) y `9/28+11/56 = 29/56`. Es el resultado BFI que NO exige well-factorabilidad completa, solo una factorizacion fija con los factores acotados.

**BFI, Teorema B (well-factorable, 4/7) [T-lit].** BFI [4, Teor. 10], segun Maynard (loc. cit., Teor. B): si `lambda` es well-factorable de nivel `Q < x^(4/7-epsilon)`, entonces `sum lambda_q (pi(x;q,a) - pi(x)/phi(q)) << x/(log x)^A`. Bloqueado para este peso por la Prop. 58.1: `lambda_t` no es well-factorable de ningun nivel `Q <= N`, porque sus dos factores primos son ambos `>= N^(1/10)` y no hay factor pequeno con el que partir el nivel.

Las tres herramientas comparten la misma anatomia: la ganancia viene de partir el modulo `d = d1 d2` con `d1` pequeno. Nuestro modulo es `d = qr` con `q, r` primos `>= N^(1/10)`: sus unicos divisores son `1, q, r, qr`. No hay factor pequeno salvo que `min(q,r)` caiga en la ventana de cada teorema. Esa es la obstruccion que se mide en la sec. 4.

## 4. Cobertura medida de cada herramienta sobre el peso

Caso de referencia: `v=10, h=96, N=10^12, t=0.05` (M1 = 5.54e-06, igual que el doc 58). Slice BFI: `1/2 < alpha+beta <= 29/56`. Slice Maynard: `1/2 < alpha+beta <= 1/2+1/100` (delta = 1/100 < 1/42). Region L_BFI: `min(alpha,beta) <= 11/56` y `max(alpha,beta) <= 9/28`. Ventana Maynard estricta: divisor en `[N^(1/50), N^(0.086)]`. Ventana generosa estilo Fouvry: divisor `<= N^(1/7)`. `R_unc` = cota trivial del error de la parte NO cubierta del slice, dividida por M (sec. 1: para un nivel hace falta `<< 1/(log N)^A`, no una constante).

| v | h | N | t | M1_sliceBFI/M1 | M1_sliceBFI en L / M1_sliceBFI | R_unc_BFI | M1_sliceMay/M1 | M1_sliceMay en ventana 1/7 / M1_sliceMay | R_unc_May | M1_L/M1 | sumc_L/t |
|---|---|---|---|---|---|---|---|---|---|---|---|
| 10 | 12 | 1e12 | 0.05 | 0.0043 | 0.195 | 0.095 | 0.0043 | 0.305 | 0.083 | 0.987 | 0.917 |
| 10 | 24 | 1e12 | 0.05 | 0.0016 | 0.104 | 0.039 | 0.0016 | 0.291 | 0.031 | 0.987 | 0.712 |
| 10 | 48 | 1e12 | 0.05 | 0.0012 | 0.079 | 0.030 | 0.0007 | 0.333 | 0.012 | 0.987 | 0.663 |
| 10 | 96 | 1e12 | 0.05 | 0.0012 | 0.054 | 0.030 | 0.0008 | 0.339 | 0.014 | 0.990 | 0.645 |
| 10 | 96 | 1e12 | 0.10 | 0.0011 | 0.056 | 0.029 | 0.0008 | 0.324 | 0.014 | 0.990 | 0.414 |
| 8 | 48 | 1e12 | 0.05 | 0.0032 | 0.054 | 0.083 | 0.0013 | 0.300 | 0.024 | 0.978 | 0.614 |
| 10 | 96 | 1e7 | 0.05 | 0.0160 | 0.050 | 0.226 | 0.0100 | 0.336 | 0.102 | 0.895 | 0.645 |

Lecturas (caso de referencia `h=96, t=0.05, N=10^12`):

- **Ventana de Maynard: cobertura EXACTA 0.** La ventana es `[x^(2 delta), x^(1/10-7 delta/5)]` con extremo superior `< x^(1/10)` para todo `delta > 0`, y todo factor primo del peso es `>= N^(1/10)`; los divisores de `d = qr` son `1, q, r, qr`, ninguno en la ventana. El script lo verifica: 0 pares en la ventana en todos los casos. Los modulos de Maynard "casi todos" tienen un divisor pequeno; los nuestros son semiprimos con ambos factores grandes, sistematicamente fuera.
- **BFI 29/56: cubre solo el 5.4% de la masa M1 del slice.** El slice es apenas el 0.12% de M1, pero el nivel se decide en el slice: la parte no cubierta tiene cota trivial `R_unc_BFI = 0.030`, una fraccion constante de M (3.0%), no `<< 1/(log N)^A`. La cobertura cae con `h` (0.195 a h=12, 0.054 a h=96): al refinar la malla, el slice queda dominado por pares con ambos factores grandes.
- **Ventana generosa estilo Fouvry (divisor <= N^(1/7)): cubre el 34% de la masa M1 del slice de Maynard.** La parte no cubierta tiene cota trivial `R_unc_May = 0.014` (1.4% de M), tambien constante.
- **El grueso de M1 SI esta en L (99%), pero eso no da el nivel.** La masa del termino principal vive en modulos pequenos (doc 58: 99.7% de M1 esta bajo `N^(1/2)`), y esos modulos pequenos tienen un factor pequeno, asi que BFI los cubre. El nivel, en cambio, es una propiedad de la cola (el slice): ahi la cobertura es del 5.4%. Que el grueso este cubierto es necesario pero no suficiente.
- **A N pequeno es peor:** a `N=10^7`, `R_unc_BFI = 0.226` (23% de M). La cota trivial decrece con N (0.226 a 0.030 de `10^7` a `10^12`) pero en todo N accesible es una fraccion constante de M, no un ahorro de potencia de logaritmo.

## 5. La reduccion por reordenamiento (nueva para este peso) y por que vuelve al tope 1/2

**Proposicion 59.1 [T].** Sea `lambda_t(d) = sum_{qr=d} c_t(q,r)` y `W(n) = sum_{d | (N-n)} lambda_t(d)`. Entonces, con `E_d(N) = A_d(N) - N/phi(d)` y `A_d(N) = sum_{n<N, n=N (mod d)} Lambda(n)`,

```text
sum_{q,r} c_t(q,r) E_{qr}(N) = sum_{n<N} Lambda(n) W(n) - N * M1.
```

*Prueba.* `sum_{q,r} c E_{qr} = sum_{q,r} c * sum_{n=N (qr)} Lambda(n) - N sum c/phi(qr)`. El primer termino es `sum_{n<N} Lambda(n) sum_{q,r: qr | (N-n)} c(q,r)`, por reagrupamiento finito (`n = N (mod qr)` si y solo si `qr | (N-n)`), y la suma interior es `sum_{d | (N-n)} lambda_t(d) = W(n)`. El segundo es `N*M1` por definicion. QED. Verificada exactamente en 7 muestras (productos de 2 y 3 primos de la malla, cuadrados, `2q`, `1`, primo no representativo): `sum_{qr|m} c = sum_{d|m} lambda_t(d)` en todos los casos, y `sum_d lambda_t(d) = t`, `sum_d |lambda_t(d)| <= ||c_t||_1`.

La reduccion no escapa del tope: para estimar `sum_n Lambda(n) W(n)` hay que aplicar la identidad de Vaughan a `Lambda(n)`, y el termino Tipo II es `sum_d lambda_t(d) sum_{mk = N (d)} a_m b_k`, una forma bilineal sobre la progresion `mk = N (mod d)` que la gran criba controla solo para `d <= N^(1/2)` (el termino `Q^2` de la sec. 2). `W(n)` es una suma de divisores sobre `N-n`, no una funcion suave de `n`: no hay estimacion Tipo II publicada que la controle mas alla de 1/2 sin factorizar el modulo. La reduccion identifica el objeto exacto que falta estimar, pero no lo estima: la dispersion para este peso es tan dificil como el problema original, porque el peso vive del lado de los modulos y su reordenamiento lo convierte en un peso divisor desplazado.

## 6. Segunda obstruccion: la clase residual es a = N, no una clase fija [P-rev]

Todos los resultados mas alla de 1/2 de la sec. 3 se enuncian para una clase residual `a` FIJA. Maynard (Introduccion, leida): las estimaciones de Deshouillers-Iwaniec "require us to only consider the situation when a is fixed (or a small power of x)". Fiorilli (Canad. J. Math. 64 (2012), Introduccion, leida): "By fixing a, one can go up to Q = x^(1/2)(log log x)^B"; sus Teoremas 1.1-1.3 dicen "Let a != 0" (fijo). Nuestra clase es `a = N`, que crece con `x` (aqui `x = N`): no es fija ni una potencia pequena de `x`. Por tanto, incluso la cobertura parcial medida en la sec. 4 (el 5.4% del slice por BFI) es CONDICIONAL a una version uniforme en `a` de esos teoremas, que no existe en la literatura. La obstruccion del factor pequeno (sec. 4) y la de la clase fija son independientes: habria que resolver las dos.

## 7. Detectar P_2 pese al nivel 1/2: en principio si, por Chen; la integracion sigue abierta

El nivel 1/2 no bloquea por si solo la deteccion de P_2: el teorema de Chen (1973) detecta numeros con a lo sumo 2 factores primos en `A_N = {N-p}` con nivel 1/2 (BV) mas el principio de switching (doc 52: identidad extraida de Liu, arXiv:2203.07871v1, ecs. (2)-(5), siguiendo a Chen 1973, Lema 9; constante 0.67). Lo que bloquea la deteccion AQUI es que `c_t` es una matriz `(q,r)` del lado de la serie singular, no un peso de criba que cuente numeros `{1,2}` (doc 58, sec. 7, razon 3; doc 56, sec. 5). La dispersion no es el paso que falta para P_2: el paso que falta es integrar `c_t` con Chen/Richert y la deteccion de P_2 (doc 58, sec. 8, #2).

## 8. Veredicto

**La dispersion NO sube el nivel efectivo para el peso completo.** La gran criba capa el Tipo II en 1/2 (termino `Q^2`); la ganancia exige un factor pequeno del modulo, y el peso no lo tiene en el grueso: la ventana de Maynard cubre exactamente el 0% (su extremo superior `< x^(1/10)` es el punto donde empieza nuestro soporte), el resultado BFI de factorizacion fija (nivel 29/56) cubre solo el 5.4% de la masa M1 del slice, y la ventana generosa estilo Fouvry cubre el 34% del slice de Maynard. La parte no cubierta del slice tiene cota trivial del 3.0% de M (1.4% en el slice de Maynard), una fraccion constante, no un ahorro de potencia de logaritmo: no hay nivel > 1/2 en el sentido de la Def. 58.1. Encima, la cobertura parcial es condicional a la clase fija (sec. 6): `a = N` no es fija. La reduccion por reordenamiento (Prop. 59.1) muestra que la dispersion para este peso equivale a un Tipo II con peso divisor desplazado, sin estimacion publicada. **P_2 pese a nivel 1/2:** posible en principio por Chen, pero la integracion con `c_t` sigue abierta. Goldbach sigue abierta.

## 9. Que sigue

1. **Una dispersion uniforme en `a` para pesos `lambda = u*v` con `u,v` soportados en primos `>= N^(1/10)` (sin factor pequeno).** No existe en la literatura: habria que probarla (o demostrar que es imposible, lo que cerraria la via). Es el blanco preciso que deja este doc.
2. **Bajar el soporte de los factores a `N^epsilon`.** La Prop. 58.1 sugiere que un peso cuyos factores bajaran hasta `N^epsilon` si seria well-factorable; eso exige cambiar la construccion del doc 54 (otro potencial, otra malla), no solo relajar margenes.
3. **Integrar con Chen/Richert y la deteccion de P_2** (doc 58, sec. 8, #2): el peso con termino principal de orden 1 (doc 56) y nivel 1/2 mas switching es la ruta de Chen; falta el acoplamiento.

## 10. Verificacion reproducible

```text
python numerico/59_dispersion_tipoII_rango_bilineal.py
```

PASS en el PC: esquina BFI exacta (`5*(9/28)+2*(11/56)=2`, `9/28+11/56=29/56`); ventana Maynard `< x^(1/10)` exacta y vacia en la malla (0 pares); reordenamiento exacto en 7 muestras (`sum lambda_t = t`); `sum c_t = t` exacto; `||c_t||_1` en rango; M1 = 5.54e-06 a `h=96` (igual que el doc 58); float==exacto a `h<=24`.

## Memoria compartida

Propuesta pendiente: el MCP canonico `lumina_shared_memory` / `lumina_propose` no esta expuesto. No se afirma persistencia remota.
