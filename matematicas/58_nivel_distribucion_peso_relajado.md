# 58 - Nivel de distribucion del peso relajado: 1/2 asintotico, insuficiente para capas

## Procedencia y estado

- **Agente:** LUMINA/Muse.
- **Fecha:** 2026-11-25 ~19:15 ET. Reclamo de ESTADO.md 19:09 ET (este paso; el bucle no lo duplica).
- **Metodo:** definicion precisa del nivel para el peso; cotas BV y trivial
  medidas exactamente (Fraction donde es barato; float con `math.fsum`
  a h grande, verificado contra exacto a h=12,24); lectura de fuente
  para la definicion de well-factorable y el teorema BFI.
- **Estado:** **[T]** para la estructura de rango 1, la Proposicion 58.1
  (no well-factorable) y las mediciones; **[T-lit]** para BV (doc 38),
  BFI (Maynard Def. 1, Teorema A; Fiorilli Teor. 1.3) y Motohashi /
  Liu-Park-Song (via doc 57). No hay cota inferior de capas ni
  demostracion de Goldbach.

## 1. Que significa "nivel" para este peso

Para la sucesion reflejada `A_N = {N-p}` (doc 38) sea
`E_d(N) = A_d(N) - N/phi(d)`, con `A_d(N) = sum_{n<N, n=N (mod d)} Lambda(n)`.
El peso relajado `c_t` (familia A del doc 56) induce el peso unidimensional
sobre el modulo producto `lambda_t(d) = sum_{qr=d} c_t(q,r)` (pares ordenados)
y el error pesado

```text
E_t(N) = sum_{i,j} c_t(i,j) * E_{q_i q_j}(N).
```

El termino principal del conteo pesado es `M = N * M1`,
`M1 = sum c_t/phi(qr)` (con `phi(q^2)=q(q-1)` en la diagonal).

**Definicion 58.1.** El nivel (asintotico) `theta(c_t)` es el mayor `theta`
tal que el error pesado con `qr <= N^theta` es despreciable frente a `M`,
es decir `|E_t| << M/(log N)^A` para todo `A>0`. El nivel es una
afirmacion asintotica (N -> infinito); el nivel efectivo a N finito
se trata aparte en la sec. 5.

Dos masas distintas, no confundir: la masa L1 `||c_t||_1` vive sobre
todo en `qr > N^(1/2)` (solo ~5% esta abajo), pero la masa del termino
principal `M1 = sum c_t/phi(qr)` vive casi toda abajo (99.7%), porque
`1/phi(qr)` pesa mucho mas los modulos pequenos. El nivel se pregunta
sobre `M1`, no sobre `||c_t||_1`.

## 2. Teorema 58.2 [T-lit/BV]: nivel >= 1/2

Por el doc 38, Prop. 38.1 (Bombieri-Vinogradov, forma maxima sobre
clases reducidas), para cada `A>0` hay `B=B(A)` tal que con
`D = N^(1/2)/log(N)^B`,

```text
sum_{d<=D, (d,N)=1} |E_d(N)| <<_A N/log(N)^A.                (38.2)
```

En la malla todos los `q_i, r_j >= N^(1/10) >= 7` son coprimos con
`N = 10^E`, asi que no hay excepciones no coprimas. El error pesado
con `qr <= D` cumple

```text
|sum_{qr<=D} c_t E_{qr}| <= maxfac * max|c_t| * sum_{d<=D} |E_d|,
```

donde `maxfac = max_d #{(i,j): q_i q_j = d}` (en la malla hay primos
representativos degenerados, `q_i = q_{i+1}`; se mide, no se supone:
`maxfac = 8` a `N=10^12, h=96`; `60` a `N=10^7, h=96`, artefacto
de la malla gruesa a N pequeno). Definiendo

```text
R_BV(A) = maxfac * max|c_t| / ((ln N)^A * M1),
```

`R_BV(A)` es la cota BV del error con `qr <= D`, dividida por `M1`
(la N se cancela). Medido: `R_BV(5) ~ 0.02`, `R_BV(10) ~ 1e-9`
(v=10, h=96, N=10^12). El error por debajo del nivel BV es
despreciable frente a `M`. Luego `theta(c_t) >= 1/2`.

Estructura que hace bien planteada la pregunta (verificada exactamente):
`c_t` es suma de 3 terminos de rango 1,

```text
c_t = (1-t)/||c0||_1 * (u1(x)v1 + u2(x)v2) + t * (1/h)(x)(1/h),
```

con `u_k(i) = F_k(i)-F_k(i-1)`, `F_k(i) = alpha_i^k` (0 fuera de
`[0,h-2]`), porque la diferencia mixta discreta de una funcion
separable es separable. Los tres factores estan soportados en primos
de `[N^(1/10), N^(1/2)]`. En particular `lambda_t` es suma de 3
convoluciones de sucesiones 1D soportadas en ese intervalo.

## 3. Tabla de medidas (v=10 salvo control; familia A)

`M1ab/M1` = fraccion de `M1` con `qr <= N^(1/2)` (test exacto
`(q_i q_j)^2 <= N`); `|c|ab` = fraccion de `||c_t||_1` abajo;
`sumcab/t` = fraccion de `sum c_t` abajo; `R_triv` = cota trivial
del error con `qr > N^(1/2)` dividida por `M1`, es decir
`ln(N) * sum_{qr>N^(1/2)} |c_t|/(qr) / M1`; `max|c|` = mayor entrada
de `c_t` (esta en el borde: la entrada de esquina es `psi` crudo,
no una diferencia mixta, ver sec. 4).

| v | h | N | t | M1 | M1ab/M1 | \|c\|ab | sumcab/t | R_triv | R_BV(5) | R_BV(10) | max\|c\| | fac |
|---|---|---|---|---|---|---|---|---|---|---|---|---|
| 10 | 12 | 1e12 | 0.05 | 2.89e-06 | 0.993 | 0.052 | 1.034 | 0.263 | 0.010 | 6.4e-10 | 0.238 | 2 |
| 10 | 24 | 1e12 | 0.05 | 3.83e-06 | 0.995 | 0.049 | 0.978 | 0.159 | 0.0077 | 4.8e-10 | 0.238 | 2 |
| 10 | 48 | 1e12 | 0.05 | 4.55e-06 | 0.997 | 0.048 | 0.951 | 0.119 | 0.0065 | 4.0e-10 | 0.238 | 2 |
| 10 | 96 | 1e12 | 0.05 | 5.54e-06 | 0.997 | 0.047 | 0.938 | 0.094 | 0.021 | 1.3e-09 | 0.238 | 8 |
| 10 | 96 | 1e12 | 0.10 | 6.69e-06 | 0.997 | 0.059 | 0.591 | 0.090 | 0.017 | 1.0e-09 | 0.225 | 8 |
| 8 | 48 | 1e12 | 0.05 | 2.19e-06 | 0.993 | 0.043 | 0.857 | 0.238 | 0.014 | 8.4e-10 | 0.238 | 2 |
| 10 | 96 | 1e7 | 0.05 | 1.34e-04 | 0.946 | 0.047 | 0.938 | 1.523 | 0.098 | 9.0e-08 | 0.238 | 60 |
| 10 | 96 | 1e7 | 0.10 | 1.70e-04 | 0.946 | 0.059 | 0.591 | 1.342 | 0.073 | 6.7e-08 | 0.238 | 60 |

Lecturas:

- **El termino principal vive abajo.** `M1ab/M1 = 0.997` a `N=10^12`
  (0.946 a `N=10^7`): casi toda la masa de `M1` esta en `qr <= N^(1/2)`,
  es decir, dentro del nivel BV. La masa L1, en cambio, esta arriba
  (`|c|ab ~ 0.05`): el peso en L1 lo domina la parte de atomos en
  `alpha,beta` grandes, que casi no aporta a `M1` porque `1/phi(qr)`
  es minusculo ahi. `sum c_t` abajo es el 94% de `t` a `h=96`.
- **La cota trivial por encima de 1/2.** A `N=10^12`, `R_triv = 0.094`
  (t=0.05, h=96): el error con `qr > N^(1/2)`, acotado trivialmente
  por `|E_d| <= (N/d) log N`, es el 9.4% de `M`. Decrece con `h`
  (0.26 -> 0.094 de h=12 a h=96) y con `N`. A `N=10^7`, en cambio,
  `R_triv = 1.52 > 1`: la cota trivial NO controla la cola. El control
  de la cola es asintotico y con margen fino (factor ~10 a `N=10^12`),
  no un teorema de nivel.
- **Error total.** A `N=10^12, h=96, t=0.05`: error/M `<= R_triv +
  R_BV(5) ~ 0.115` (con A=10: `~0.094`). El termino principal
  sobrevive: el nivel no mata el peso. A `N=10^7` el error total
  supera `M` por la cola (`R_triv > 1`).

## 4. Dos hechos de soporte (exactos)

- **La mayor entrada esta en el borde.** `max|c_t| = 0.2375` en todos
  los casos (no depende de `h` ni de `N`): es la entrada de esquina
  `(h-1,h-1)`, que recibe un solo atomo con coeficiente `psi` crudo
  (`~0.075/||c0||_1`), no una diferencia mixta `O(s^2)`. Las entradas
  interiores son `O(s^2)`. La cota BV usa `maxfac * max|c_t|`, que
  sigue siendo `O(1)` (1.9 a `h=96, N=10^12`).
- **`||c_t||_1 < 1`.** Medido `0.98-1.00` (no 1 exacto): la parte
  constante positiva cancela parcialmente las entradas negativas de
  los atomos (doc 56 ya lo advertia: la correccion es `<=3%`).
  Las fracciones de masa de la sec. 3 usan el `||c_t||_1` real.

## 5. Nivel efectivo a N finito: el nivel 1/2 no llega al soporte

La perdida `(log N)^B` de BV da el nivel efectivo
`theta_eff(B,N) = 1/2 - B*ln(ln N)/ln(N)`. A `N=10^12`
(`ln N = 27.6`): `theta_eff = 0.38 (B=1), 0.26 (B=2), 0.14 (B=3),
-0.10 (B=5)`. El soporte del peso empieza en `alpha+beta = 0.2`
(modulo minimo `N^(1/5)`). Cobertura de masa L1: `B=1: 1.5%`,
`B=2: 0.2%`, `B>=3: 0`. A `N=10^7`: `theta_eff = 0.33 (B=1),
0.16 (B=2)`, cobertura `0.7%` y `0`. Es decir: **a la escala del
modelo, el nivel BV no alcanza el soporte del peso en ningun sentido
util**; el nivel 1/2 es una afirmacion asintotica (`N -> infinito`),
no una cobertura a `N=10^12`. Esto no contradice la sec. 3: alli el
corte es `N^(1/2)` exacto (asintotico), y la cola se controla por la
cota trivial, no por BV.

## 6. Proposicion 58.1 [T]: lambda_t NO es well-factorable

**Definicion (leida).** Maynard, "Primes in arithmetic progressions
to large moduli II: Well-factorable estimates", arXiv:2006.07088v1,
Def. 1 (Memoirs Amer. Math. Soc. 306 (2025), no. 1543): una sucesion
`lambda_q` es well-factorable de nivel `Q` si, para toda factorizacion
`Q = Q1 Q2` con `Q1,Q2 >= 1`, existen `gamma^(1), gamma^(2)` con
`|gamma| <= 1`, soportadas en `q_i <= Q_i`, tales que
`lambda = gamma^(1) * gamma^(2)` (convolucion de Dirichlet).

**Teorema A (BFI, leido).** Maynard, ibid., Teor. A (Bombieri,
Friedlander, Iwaniec, "Primes in arithmetic progressions to large
moduli", Acta Math. 156 (1986), Teor. 10; ver tambien Fiorilli,
"On a theorem of Bombieri, Friedlander and Iwaniec", Canad. J. Math.
64 (2012), Teor. 1.3): si `lambda` es well-factorable de nivel
`Q <= x^(4/7-epsilon)`, entonces
`sum_{q<=Q,(q,a)=1} lambda(q)(pi(x;q,a)-pi(x)/phi(q)) << x/(log x)^A`.
Maynard extiende a `x^(3/5-epsilon)` para pesos triplemente
well-factorables. En ambos casos la hipotesis es la well-factorabilidad.

**Proposicion 58.1 [T].** `lambda_t` no es well-factorable de ningun
nivel `Q <= N`.

*Prueba.* `lambda_t` esta soportada en semiprimos `d = qr` (y
cuadrados `q^2`) con `q,r` primos en `[N^(1/10), N^(1/2)]`
(hecho de soporte verificado por el script: menor primo `>= N^(1/10)`,
mayor `<= N^(1/2)`, y `c_t(h-1,h-1) > 0`, asi que el semiprimo
`d* = q_max * r_max ~ N` esta en el soporte con `lambda_t(d*) != 0`).
Sea `Q <= N` y tomese una factorizacion `Q = Q1 Q2` con
`1 <= Q1 < N^(1/10)` (posible: si `Q >= N^(1/10)/2` tomese
`Q1 = N^(1/10)/2`; si no, `Q1 = Q/2`). Para `d = qr` en el soporte,
los divisores son `1, q, r, qr`, todos `> Q1` salvo `1`; luego
`(gamma^(1)*gamma^(2))(d) = gamma^(1)(1) * gamma^(2)(d)`. Como
`lambda_t(d) != 0` para el `d*` de arriba, `gamma^(1)(1) != 0` y
`gamma^(2)(d) = lambda_t(d)/gamma^(1)(1)` para todo `d` en el soporte,
asi que `gamma^(2)` esta soportada en todo el soporte, hasta `~N`:
`Q2 >= N`. Pero `Q2 = Q/Q1 < N` si `Q <= N` y `Q1 >= 1` (y si
`Q1 = Q/2`, `Q2 = 2 < N` directamente). Contradiccion. QED.

El argumento vale pieza por pieza para cada uno de los 3 terminos
de rango 1 de la sec. 2 (cada pieza tiene los mismos soportes), asi
que tampoco hay factorizacion por linealidad. **Consecuencia: BFI
(4/7) y Maynard (3/5) no se aplican a este peso.** La obstruccion es
estructural: los dos factores primos son ambos grandes
(`>= N^(1/10)`); no hay factor pequeno con el que partir el nivel.
Un peso cuyos factores bajaran hasta `N^epsilon` si seria
well-factorable; este no.

## 6.1. La ruta Motohashi / Liu-Park-Song (doc 57 del bucle)

El doc 57 del bucle (19:14 ET) audito el principio de induccion de
Motohashi y el Teor. 3.2 de Liu-Park-Song: para un numero fijo `r` de
primos distintos en intervalos, BV vale hasta `q <= N^(u/2-epsilon)`
(producto `< N^u`), es decir, nivel `< 1/2` para cada indicador no
ponderado de caja. Tambien demostro que la discretizacion ingenua del
peso Lipschitz acoplado `W_P^+` no transfiere a A4 (al pasar el error
supremo a la norma L1 sobre modulos aparece un factor `Q`).

Para `c_t` la situacion es distinta en un sentido favorable: `c_t`
es EXACTAMENTE suma de 3 terminos separables (sec. 2), sin
discretizacion ni aproximacion Lipschitz. La maquinaria de cajas de
Motohashi / Liu-Park-Song se aplica directamente a cada pieza de
rango 1 (cada pieza es una convolucion de dos sucesiones 1D sobre
primos en intervalos, `r = 2`): es una segunda ruta independiente
al nivel 1/2, y la obstruccion del doc 57 (el factor `Q` de
discretizar un peso acoplado) no se aplica a `c_t`. Pero tampoco
sube el nivel: Liu-Park-Song da `q <= N^(u/2-epsilon)` con `u < 1`,
es decir, nivel `< 1/2` (y `u = 1` en nuestro caso, producto hasta
`N`, da `1/2`). La ruta Motohashi confirma el nivel 1/2 y nada mas.

Fuentes (leidas por el doc 57, [T-lit] via doc 57): Y. Motohashi,
*An induction principle for the generalization of Bombieri's prime
number theorem*, Proc. Japan Acad. 52 (1976), 273-275, Teors. 1-2;
Y. Liu, P. S. Park y Z. Q. Song, *Bounded gaps between products of
distinct primes*, Research in Number Theory 3 (2017), art. 26,
Lema 3.1 y Teor. 3.2.

## 7. Veredicto

**Nivel: 1/2 asintotico.** El error pesado por debajo de `N^(1/2)`
es despreciable (`R_BV(10) ~ 1e-9`); la cola por encima se controla
con la cota trivial a `N` grande (`R_triv = 0.094` a `N=10^12`),
pero no a `N` pequeno (`R_triv = 1.52` a `N=10^7`). El termino
principal `sum c_t = t` sobrevive al analisis de error a `N=10^12`
(error total `<= ~12%` de `M`).

**NO suficiente para una cota inferior de capas.** Tres razones:

1. El nivel es 1/2, y el doc 38, Prop. 38.2, demuestra que con nivel
   1/2 el tamiz lineal da `f(s) = 0` para las capas en cuestion
   (`s < 2` cuando `z > N^(1/4)`). El nivel 1/2 es exactamente el
   nivel que no alcanza. Subirlo a 4/7 exigiria BFI, bloqueado por
   la Prop. 58.1; la ruta Motohashi / Liu-Park-Song (sec. 6.1) da
   nivel `< 1/2` para las cajas, no mas.
2. El control de la cola no es un teorema de nivel: es la cota
   trivial, con margen fino y que falla a `N=10^7`. No hay estimacion
   de dispersion publicada para el rango bilineal
   `[N^(1/10), N^(1/2)]^2` de este peso.
3. El peso es una matriz `(q,r)` del lado de la serie singular, no
   un peso de criba que detecte `P_2` (doc 56, sec. 5). El nivel,
   aunque fuera mayor, no bastaria por si solo.

Lo que este doc SI dice: el nivel del peso relajado es 1/2
(asintotico), el error pesado esta controlado a `N=10^12` con margen
`~10x`, y la ruta BFI esta cerrada con prueba (Prop. 58.1). Lo que
NO dice: no hay cota inferior de capas, no hay nivel `> 1/2`, no hay
demostracion de Goldbach.

## 8. Que sigue

1. **Estimacion de dispersion / Tipo II para el rango bilineal
   `[N^(1/10), N^(1/2)]^2`, o una version ponderada y uniforme del
   principio de Motohashi** (el blanco preciso del doc 57, sec. 4:
   admitir directamente el coeficiente multivariable acoplado, o
   una representacion separable con control uniforme). Es el paso
   que falta para un nivel genuino `> 1/2` de este peso. No es una
   cita: habria que probarla (o localizar una version publicada del
   principio de induccion de Bombieri para este rango, como pide el
   doc 53). La estructura de rango 1 (sec. 2) es el punto de
   partida: cada pieza es una convolucion `u*v` con `u,v`
   soportadas en primos de `[N^(1/10), N^(1/2)]`.
2. Integrar con Chen/Richert y la deteccion de `P_2` (doc 56,
   "Que sigue" #2): un peso `(q,r)` con `sum c` de orden 1 no basta
   si no cuenta numeros `{1,2}`.
3. El Lema 56.1 (ningun peso mejor en este marco) y la Prop. 58.1
   (BFI no aplica) cierran dos vias; cualquier mejora exige cambiar
   el marco (otro factor local, otra identidad) o probar dispersion.

## 9. Verificacion reproducible

```text
python numerico/58_nivel_distribucion_peso_relajado.py
```

PASS en el PC: cross-check float vs exacto a h=12,24 (M1);
rango 1 exacto (3 terminos) en todos los casos; `sum c_t = t`
exacto; `||c_t||_1 < 1` registrado; `maxfac` medido; esquina
positiva; tabla completa (v=10,8; h=12,24,48,96; N=10^12,10^7;
t=0.05,0.10).

## Memoria compartida

Propuesta pendiente: el MCP canonico `lumina_shared_memory` /
`lumina_propose` no esta expuesto. No se afirma persistencia remota.
