# 08 — Revisión adversarial de 06–07, bibliografía y residuos

## Procedencia

- Agente: **Codex**. Modelo exacto no expuesto por la interfaz.
- Fecha de proyecto: **2026-11-23 ET**; los JSON nuevos pueden llevar 2026-11-24 UTC.
- Encargo: `notas/2026-11-23_TAREA_PARA_CODEX.md`, leído antes del análisis.
- Método: lectura de README, ESTADO, 04–07 y scripts; derivaciones algebraicas; búsqueda en fuentes primarias; reproducciones y comprobaciones por fórmulas distintas.
- Estado: **[P-rev]**, revisión NO independiente: se comparte usuario de Windows. Resultados numéricos **[N]**. No se modifica 06 ni 07. Ningún cálculo finito demuestra Goldbach.
- Este documento registra errores explícitamente y distingue reparación del argumento, condición adicional y verificación formal.

## 1. Veredictos de T1

| Punto | Veredicto | Detalle y alcance |
|---|---|---|
| 06, Lema 6.1: BV con máximos en y y clase | **Correcto en el dominio de la cola, con precisiones** | Hay que exigir K al menos una potencia fija de N, y 2K dentro del nivel BV, no solo K. |
| 06, Lema 6.1: μ·1_coprimo·id/φ | **Identidad correcta; error en la cota de Rankin escrita** | La potencia y^(−1/4) requiere factores p^(−1/2), no p^(−3/4). Cualquiera de las dos correcciones conserva la conclusión. |
| 06, Lema 6.2: (q,c)=1 | **Correcto** | Un divisor común de q y c divide N y c; se usa (c,N)=1. |
| 06, Lema 6.2: BV para μ, máximo sobre clases no reducidas | **Hueco de justificación** | La fuente primaria localizada formula clases reducidas. La extensión invocada y su coste uniforme no quedan demostrados en 06; no se declara falsa, pero no se valida como entrada ya acreditada. |
| R-2: D=√N/log^B N | **Correcto, con B dependiente del ahorro** | La demostración de 04 se adapta; D no tiene que ser exactamente una potencia. |
| R-3: identidad exacta de la parte m=cd compuesta | **Error: falta −Tbad** | Contraejemplo explícito N=32,D=3. La asintótica T'=N𝔖+o(N) sobrevive. |
| 07, 7.1: Parseval | **Error de igualdad, reparación inmediata** | Al sumar solo N≤X se tiene ≤, no =, frente a la integral completa. |
| 07, 7.1: Davenport con θ=cα | **Correcto** | La uniformidad en θ es la adecuada; el corte cd<X se trata con extremos parciales. |
| 07, 7.1: retirada de coprimalidad | **Correcto** | El error de potencias de primos que dividen N se controla uniformemente. |
| 07, 7.1: unión sobre cajas | **Correcto para familia fija; hueco para cajas móviles** | Las cajas de los scripts dependen de D(N). Se aporta en §4 un argumento maximal que permite reparar el paso pagando potencias logarítmicas. |
| 07, P_Q y «arcos menores» | **Identificación no demostrada** | P_Q es un modelo aditivo truncado de Ramanujan; su residuo no es automáticamente la integral sobre arcos menores. |

## 2. Lemas 6.1 y 6.2

### 2.1 Dominio y nivel de BV [P-rev]

Para la cola original basta K≥cN^(1/3), o más generalmente K≥N^δ con δ>0 fijo. Entonces log K es comparable con log N. Para cada A se escoge B de forma que **2K≤√N/log^{B_BV}N**; el factor 2 se absorbe aumentando B.

La versión con máximo en y y en la clase es clásica: [Alisa Sedunova, *A partial Bombieri–Vinogradov theorem with explicit constants* (2018), Teorema 1.1, pág. 101](https://numdam.org/item/10.5802/pmb.24.pdf), enuncia exactamente ambos máximos y las clases reducidas. [Michel, §3.1](https://www.epfl.ch/labs/tan/wp-content/uploads/2018/10/Parkcitylectures.pdf) da también la versión agregada necesaria para R-2. No se invoca BV por encima de su nivel. Las notas de Brandes aparecieron en resultados indexados, pero su PDF devolvió 403; no se usan como respaldo de una lectura completa.

La condición inferior sobre K no es decorativa. El enunciado literal para **cualquier** K pequeño sería falso: con N=2^j, K=2 y C=N/16, el único divisor con μ(d)≠0 en (2,4] es d=3. La condición sobre c selecciona c impar. PNT en la clase reducida fija módulo 6 da

    B_N(C,2) ∼ −(3/2)C log 3 = −(3 log 3/32)N.

Esto no pertenece a d>D=N^(1/3), por lo que no invalida la aplicación prevista. Sí exige explicitar el dominio del lema.

### 2.2 Término principal y corrección de Rankin [P-rev]

Si a_N(n)=μ(n)1_{(n,N)=1}n/φ(n), entonces a_N=μ*h_N con

    h_N(p^k)=−1/(p−1) si p∤N;  h_N(p^k)=1 si p|N.

La identidad es correcta prima por prima. Para el tramo b>√y:

    Σ_{b>√y}|h_N(b)|/b
      ≤ y^(−1/4)Σ_b |h_N(b)|/b^(1/2)
      ≪_ε y^(−1/4)N^ε.

Aquí el producto en los primos de N es ∏(1−p^(−1/2))⁻¹. Si se desea usar ∏(1−p^(−3/4))⁻¹, el factor obtenido es **y^(−1/8)**. 06 mezcla ambos exponentes. Ambos dan ahorro suficiente cuando y≥N^δ; no cambia el resultado final del Lema 6.1.

En el tramo b≤√y, PNT cuantitativo para μ y Σ|h_N(b)|/b≪log(2N) bastan. No es necesario optimizar esta última cota a log log N. El peso g(d) no es literalmente lineal a trozos —contiene log d y N/d—, pero tiene el tamaño y la variación requeridos O(C log N). Esto justifica la sumación parcial.

### 2.3 Tipo I: lo que sí está probado y lo que falta [P-rev]

Si q|N−cd y (c,N)=1, entonces (q,c)=1 y d≡N c⁻¹ mod q. La identidad de Vaughan y |λ(q)|≤log q son correctas. Con la entrada fuerte que 06 escribe, el resto de la estimación por variación y suma sobre c,q funciona, ajustando las potencias logarítmicas.

Pero la clase tiene gcd(N c⁻¹,q)=gcd(N,q), que puede ser >1. El BV reducido para μ no cubre literalmente ese máximo. [Granville–Shao, arXiv:1703.06865v2, (1.1) e introducción](https://arxiv.org/pdf/1703.06865v2) recoge el caso μ con **(a,q)=1**. No basta citarlo para máximo sobre todo r.

La reducción correcta comienza así. Si g=(r,q) no es libre de cuadrados, la suma de μ en esa clase es cero. Si g es libre de cuadrados, escribiendo n=gm,

    μ(gm)=μ(g)μ(m)1_{(m,g)=1},
    m≡r/g mod q/g.

La condición adicional de coprimalidad es real. Por ejemplo, μ(2·2)=0 mientras μ(2)μ(2)=1. No se puede borrarla ni pagar un «divisor de q» sin demostrar el coste agregado.

Una identidad útil para completar la reparación es

    μ(m)1_{(m,g)=1}=Σ_{k|m; k|g^∞} μ(m/k).

En la clase reducida módulo q/g solo contribuyen k coprimos con q/g. Esto reduce a clases reducidas con longitud x/(gk), pero hay que controlar simultáneamente los rangos, los máximos en extremos y la suma sobre g,k. **Ese argumento uniforme no está completado en esta revisión.**

Veredicto final: **no elevar a demostrado el Lema 6.2 ni B=B₃+O(N/log^A N) con la justificación actual**. Es una conclusión condicional a la versión de BV no reducida o a su reparación detallada, no una refutación de la cota en sí.

## 3. R-2 y R-3

### 3.1 R-2 se mantiene [P-rev]

Con D=⌊√N/log^{B(A)}N⌋, BV cubre los módulos y 04 §4.3 controla A_N(D): log D≍log N y N^εD^(−1/8) sigue ahorrando una potencia si ε se escoge pequeño. Así

    R(N)=N𝔖(N)−T_N^prim(D)+O_A(N/log^A N)

sigue siendo válido con esos parámetros, usando las entradas clásicas indicadas en 04. B depende de A; no se obtiene una sola cota con un B universal para todos los ahorros.

### 3.2 Identidad exacta corregida para R-3 [P-rev]

Escribimos Tfull=Tprim+Tbad y particionamos Tprim según si **m=cd** es primo o compuesto. Con la notación de 06,

    Tprim=−R_pr,>D+T'_comp,
    T'_comp=S_N(D)−R_pr,≤D−R_pot−Tbad.                 (8.1)

06 omite el último término en su igualdad exacta. Con N=32,D=3,

    Tbad=log 2·log 36=2,483906048674… ≠0.

En efecto, las partes malas proceden de N−cd=2,4,8 y suman log 2·(log 3+log 2+log 6). El verificador contrasta (8.1) con doble bucle.

**La conclusión asintótica sí es correcta:** R_pr,≤D≪D log²N, R_pot≪√N log³N, Tbad≪_εN^εlog³N. Para D=N^(1/3) o √N/log^B N son despreciables respecto de N/log^A N. De ahí T'_comp=N𝔖(N)+O_A(N/log^A N).

Precisión importante: esta partición es por **m compuesto**. La columna Q de `cajas_bilineales.py` es por **d compuesto**: incluye una partición distinta, por eso Q puede ser ≈10N mientras T'_comp≈N𝔖(N). Tampoco Σ log p·Λ(N−p) es exactamente R*: el segundo sumando puede ser una potencia propia, con error controlable.

Tprim=o(N) equivale a la asintótica global pertinente, mediante la reducción. No se ha probado que la conjetura de mera existencia implique HB en todas las cajas ni una separación uniforme κ>0. En 06 §8, «Goldbach equivale a HB+» debe leerse como **criterio suficiente cuantitativo**, no como equivalencia demostrada con existencia.

## 4. Teorema 7.1: reparación del argumento [P-rev]

### 4.1 Cajas fijas

Quitada la coprimalidad, el error total es ≤K_N≪_εN^εlog³N. Para una caja fija y el corte cd<X, sea a(m) su coeficiente, F(α)=Σa(m)e(mα), S(α)=Σ_{n≤X}Λ(n)e(nα). La corrección de Parseval es

    Σ_{N≤X}|B̃_N|² ≤ Σ_{N≤2X}|(a*Λ)(N)|²
                    =∫|F|²|S|²≤sup|F|²·Σ_{n≤X}Λ(n)².

La igualdad inicial impresa en 07 es incorrecta, pero la desigualdad basta. Davenport se aplica uniformemente con θ=cα, también a intervalos truncados; con K≥X^δ,

    sup|F|≪_J CK/(log X)^(J−1),
    Σ_{N≤X}|B̃_N|²≪_J X³/(log X)^(2J−3),

si CK<X. Chebyshev y una unión sobre O(log²X) cajas **fijas en X** dan el conjunto excepcional deseado tras elegir J.

### 4.2 Cajas que dependen de D(N): paso maximal que faltaba

Los scripts usan cortes (D(N)+1)2^j. No basta una unión sobre O(log²X) polinomios Fourier si esos polinomios cambian con N. Esta es una laguna adicional, no una consecuencia de cambiar = por ≤.

Se puede reparar así. Fijamos una escala d∈[K,2K), K potencia de dos, K≥X^δ, y una caja fija en c. Para un subintervalo diádico I de longitud H, combinamos cota trivial y Davenport:

    sup|F_I|≪_J C log X·min(H,K/(log X)^J).

Por Parseval, la energía de la correspondiente convolución es ≤const·C² X(log X)³ min(H,K/log^JX)². Hay O(K/H) subintervalos a esa escala y

    Σ_H (K/H) min(H,K/log^JX)²≪K²/log^JX.

Todo prefijo se descompone en ≤O(log X) intervalos diádicos. Cauchy–Schwarz **sobre esa descomposición finita** da

    Σ_{N≤X}max_{t∈[K,2K]}|B̃_N([K,t])|²
       ≪_J C²K² X/(log X)^(J−4)
       ≪_J X³/(log X)^(J−4).                         (8.2)

La suma geométrica se verifica separando H≤K/log^JX y H>K/log^JX. No se ha supuesto ortogonalidad entre los sumandos en d.

Cualquier intervalo (K',2K'] atraviesa a lo sumo dos escalas diádicas fijas y se expresa mediante diferencias de prefijos. Elegimos J>B+2A+6 y tomamos la unión sobre escalas fijas en c,d. Chebyshev en (8.2), el error K_N y el Lema 6.1 en el rango inferior dan la conclusión de 7.1 **con el paso maximal añadido**. Los cortes d>D(N) de la caja frontera también quedan cubiertos por prefijos.

Por tanto: **conclusión reparable y derivación escrita aportada; prueba original incompleta para la convención móvil de cajas**. No se convierte «casi todo N» en «todo N».

## 5. Qué significa realmente P_Q [P-rev]

### 5.1 Fórmulas exactas y dos errores de lectura

P_Q es la convolución finita de μ con el modelo Λ_Q. Usando c_q(m)=Σ_{h|q,h|m}hμ(q/h),

    P_Q(N)=Σ_{h≤Q} a_h(Q) M(N−1;h,N),
    a_h(Q)=hμ(h)Σ_{q≤Q;h|q} μ(q)²/φ(q),                (8.3)

donde M(y;h,r)=Σ_{n≤y;n≡r mod h}μ(n). Se deriva con sumas finitas, sin intercambiar series condicionales.

Para N par, **q=2 solo** da −Σμ(n)(−1)^n, no 2Σ_impares μ(n). Esta última es la suma de **q=1 y q=2**. Además q=3 entra en Q=3; Q=6 es el primer corte ensayado en 07 que incluye simultáneamente 3 y 6. «Aparece exactamente en Q=6» no se sigue de las muestras.

### 5.2 Residuo del modelo frente a arcos menores

Un corte en denominadores de una serie de Ramanujan no es una proyección espectral sobre los intervalos de arcos mayores. Para identificarlo con una integral hay que especificar anchos, pesos y controlar el error de aproximación. 07 no lo hace. Por ello **M0−P_Q se llama aquí residuo del modelo de Ramanujan**, sin afirmar que sea la contribución exacta de arcos menores.

La expansión clásica de Hardy usa **(φ(m)/m)Λ(m)**, no literalmente Λ(m); y es condicionalmente convergente. Aun en el límite no se deben suprimir esos detalles. La diferencia en la convolución de μ está acotada elementalmente por O(log²N): en p^k el defecto es log p/p y se suma sobre k. Esto explica que el modelo pueda ser bueno a la escala del experimento sin ser una identidad exacta para Λ.

Por Davenport, directamente de la suma exponencial de c_q,

    |P_Q(N)|≪_J Q N/(log N)^J.                         (8.4)

Así, para Q fijo, o Q una potencia fija de log N, es o(N) sin GRH. La cota O_Q(N^(1/2+ε)) que se propone en 07 requeriría detallar las consecuencias de GRH en las progresiones de (8.3). Una cota superior de ese tipo **no demuestra** un término de tamaño ≍√N ni un signo fijo.

## 6. Bibliografía: hallazgos de T2

**[T, modelo conocido]** Bhowmik–Grimmelt, [arXiv:2607.27282v2, ecuación (2.6), pág. 6](https://arxiv.org/pdf/2607.27282v2), define exactamente Λ_R(n)=Σ_{q≤R}μ(q)c_q(n)/φ(q). Por tanto la fórmula P_Q=μ*aditivoΛ_Q es su aplicación lineal inmediata; **no es una fórmula nueva del proyecto**. Su modelo no debe confundirse con una igualdad de proyecciones de Fourier.

**[T, antecedentes y límites]** Coppola–Murty–Saha, [arXiv:1507.07862, (1)–(2) y fórmula de Hardy](https://arxiv.org/pdf/1507.07862), y [arXiv:1612.03136, Teorema 2](https://arxiv.org/pdf/1612.03136), estudian expansiones y correlaciones con condiciones de convergencia/decaimiento. Esas condiciones no se han demostrado aquí para μ×Λ y no justifican una asintótica punto a punto con desplazamiento dependiente de N.

**[T, correlaciones relacionadas]** Lichtman, [Averages of the Möbius Function on Shifted Primes, arXiv:2009.08969v2, Teoremas 1.1 y 1.8](https://arxiv.org/pdf/2009.08969), obtiene cancelación al promediar desplazamientos. Los signos de las formas y los promedios difieren de la correlación reflejada para cada N. No es una prueba de HB+.

**Resultado de la búsqueda.** No se localizó en las fuentes primarias consultadas una tabla o teorema que reproduzca específicamente las medias normalizadas por 3|N y los números −0,22/−0,49 de 07. Esto **no prueba novedad ni ausencia de publicación**. El mecanismo de modelización y la reducción a sumas de μ por clases son conocidos. La observación concreta permanece [N] del proyecto.

Búsquedas empleadas: «Möbius von Mangoldt Goldbach correlation», «Ramanujan Coppola Murty Saha correlations», «Möbius Goldbach bias», «Möbius Chebyshev Goldbach», «Averages Möbius shifted primes». Se descartaron resúmenes automáticos y repositorios de supuestas pruebas como autoridad matemática. Solo se citan arriba fuentes primarias leídas.

## 7. Reproducción y ampliación numérica [N]

Las dos ejecuciones pedidas a 10^6 terminaron con exit 0 y reproducen las tablas:

- `20261123_223206_mobius_goldbach_aleatoriedad_1000000.json`: Var(Z_M0)=0,8956531; media si 3|N=−0,2193090; correlación muestreada con 𝔖=−0,4893917; error relativo FFT/directo 3,65·10⁻¹⁵.
- `20261123_223326_arcos_mayores_mobius_1000000.json`: para Q=30, Var(residuo)/media(V0)=0,7955311; medias por divisibilidad 3,5,7: 0,0010744, 0,0010391, 0,0106889.

**Advertencias de muestreo.** Los scripts previos incluyen N=X/2 pese a describir una ventana abierta; por eso cuentan 250001 pares en vez de 250000 a 10^6. Además `pares[::pares.size//5000]` toma pasos 100 en N a 10^6 y 1000 a 10^7: la muestra de la correlación con 𝔖 queda restringida a múltiplos de 5. Se reproducen sus valores, pero no se interpretan como correlación representativa de todos los pares.

La ampliación usa ventana estricta, semilla 20261123, muestreo uniforme sin reemplazo para correlaciones y 25 casos por escala para la cola completa. Archivo: `20261124_033802_732515_codex_residuos_06_07.json`. Código nuevo: `numerico/codex_residuos_06_07.py`; no modifica scripts anteriores.

| X | Media residuo Z, Q=30 | Asimetría | P(Z>3) | P(Z<−3) | fracción Z>0 |
|---:|---:|---:|---:|---:|---:|
| 10^5 | −0,0005285 | −0,17217 | 0,000280 | 0,001520 | 0,504920 |
| 10^6 | +0,0000536 | −0,08260 | 0,000320 | 0,001120 | 0,501432 |
| 10^7 | −0,0000632 | −0,02637 | 0,000394 | 0,0006856 | 0,500862 |

Z usa sqrt(V0(N)) por observación. Su varianza no es idéntica a Var(residuo)/media(V0); se guardan ambas. La asimetría negativa de las colas disminuye, pero la media no tiene un signo estable. Los promedios por bloques contiguos también cambian de signo. No se usan p-valores iid: las observaciones comparten sumandos.

**Proyección Fourier adicional.** Se calculó por separado una proyección discreta sobre |α−a/q|≤Q/(qX), Q=30, con todos los parámetros de malla en el JSON. A 10^7 su residuo normalizado tiene media −0,0000798 y asimetría −0,02659. Concuerda cualitativamente, pero es una proyección de malla: **no se certifica su error respecto de la integral continua**.

### 7.1 Cola completa, todas las c [N]

Se evalúa Tprim=S−R−Tbad; esta identidad se contrastó con doble bucle para todos los pares hasta 200. Se calcula también un modelo T_Q^prim sustituyendo Λ por Λ_Q en **todas** las fibras c primitivas. Sus sumas se evalúan mediante prefijos de μ(d)log d por clases, con comprobaciones directas en N=32,64,128.

| X | Media T/(N𝔖) | Mínimo | Máximo | Media (T−T_30)/(N𝔖) | fracción del residuo >0 |
|---:|---:|---:|---:|---:|---:|
| 10^5 | −0,05033 | −0,73686 | +0,28130 | −0,00534 | 0,52 |
| 10^6 | −0,02814 | −0,59039 | +0,26309 | +0,00211 | 0,48 |
| 10^7 | +0,04458 | −0,23232 | +0,19322 | +0,00693 | 0,48 |

Son **25 N por escala**, no un barrido exhaustivo. La normalización usa la serie singular clásica numérica del proyecto. Evaluar T usando R no constituye una prueba independiente de HB+: es un diagnóstico finito. No se observa un sesgo uniforme utilizable para demostrar T≤(1−κ)N𝔖. Que el máximo muestreado esté por debajo de 1 no es una cota universal.

### 7.2 Integridad de la implementación

- El selftest detectó inicialmente una consulta con índice negativo en un prefijo de clase vacía; se corrigió antes de la ejecución de escalas. Se probaron otra vez las fórmulas contra enumeración literal.
- Se corrigió durante la comprobación una transcripción preliminar del valor simbólico de Tbad en N=32; se añadió la aserción numérica de log 2·log 36. El valor usado en este informe es el comprobado.
- Las tablas proceden de la ejecución que pasó los controles, no de los intentos fallidos. El JSON contiene hashes de fuentes, semilla, parámetros, versiones y errores FFT/directos.
- Las tres escalas tardaron aproximadamente 0,16 s, 1,31 s y 18,81 s. Se ejecutaron secuencialmente, sin modelos ni cambios de servicios.

## 8. Lean y cierre

**[T-Lean] T4 completada para las tres piezas finitas solicitadas.** Fuente canónica nueva: `C:\Users\dal_n\lean\goldbach_lumina\GoldbachLumina\Reduccion2.lean`; copia idéntica en `lean/Reduccion2.lean`. No se modifica `Reduccion.lean`.

| Declaración (namespace `GoldbachLumina.Reduction2`) | Alcance verificado |
|---|---|
| `reindex_divisor_tail`, `reindex_mobius_tail` | Biyectivamente n=cd, con 0<n<N, d>D, c desde 1 hasta (N−1)/(D+1), y D+1≤d≤(N−1)/c. N>0; todo D natural. |
| `proper_prime_power_error` | 0≤R−R*≤2√N_floor·log₂N_floor·(log N)² para N≥1, exactamente la desigualdad finita de (4.12). El paso asintótico o(N) no se formaliza en esta declaración. |
| `primeConvolution_eq_log` | Identifica R* con la suma solo sobre dos primos, con pesos log n·log(N−n). La resta en el argumento de log se tipa explícitamente en naturales. |
| `badTail_eq_product` | La forma de divisores de la cola mala coincide con la restricción (c,N)>1 en todos los pares c,d de la reindexación. |
| `bad_tail_bound` | Primera desigualdad explícita de Prop. 4.2: valor absoluto de la cola mala ≤K_N, para N≥1 y cualquier D. |

En `badPrimePowerBound`, K_N se escribe como suma finita sobre p primo que divide N, 1≤k≤Nat.log 2 N y p^k<N. Este corte no pierde términos: p≥2 implica k≤Nat.log 2 N. `Nat.divisors.card` representa τ(N−p^k). No se formalizan aquí la segunda desigualdad de Prop. 4.2 (divisores pequeños no coprimos), ni K_N≪_εN^εlog³N, ni las estimaciones de 06/07.

Comando ejecutado: `lake env lean GoldbachLumina/Reduccion2.lean`, salida **0**. Los seis `#print axioms` de la fuente final dan únicamente **propext, Classical.choice, Quot.sound**; no aparece `sorryAx`. No hay declaraciones `sorry` ni axiomas añadidos. El registro con hash y salida se conserva en `lean/2026-11-23_codex_Reduccion2_verificacion.json`.

Durante el desarrollo se corrigieron fallos de elaboración de tipos y nombres de API: una expresión provisional `Real.log (N−n)` se elaboró con resta real, y se sustituyó por la coerción explícita de la resta natural; algunas comprobaciones usaron un nombre de lema no existente o ambiguo. Los intentos con errores no se contabilizan como PASS; solo se copia y certifica la fuente final que compiló sin errores ni advertencias. Las seis declaraciones anteriores de `Reduccion.lean` siguen conservando su alcance anterior.

## 9. Continuidad y límites

Las propuestas de memoria compartida de 06/07 y este aporte quedan **PENDIENTES**: el catálogo de herramientas de esta sesión no expone `lumina_shared_memory`, `lumina_recall_v2` ni `lumina_propose`. Un archivo local pendiente no es un recibo remoto ni acredita persistencia en MCP. No se usan REST, credenciales o service_role para sustituirlo.

Registro local preparado, con autoría y estados separados: `notas/2026-11-23_CODEX_revision_06_07_memoria_pendiente.json`. El cierre, los hashes y las limitaciones se recogen en `notas/2026-11-23_CODEX_revision_06_07_cierre.json`.

**Siguientes obligaciones matemáticas:** (i) justificar el máximo sobre clases no reducidas de μ o rehacer Tipo I sin usarlo; (ii) si se habla de arcos menores, demostrar la comparación entre P_Q y la integral con errores y parámetros explícitos; (iii) aportar una cota unilateral aritmética uniforme en cada N, no inferirla de medias o casi todos; (iv) obtener un umbral efectivo y cubrir el tramo finito si se pretendiera una prueba completa. La reparación maximal de §4 es [P-rev] escrita, no [T-Lean].

**Estado de esta revisión:** los errores y huecos anteriores son explícitos. El argumento para casi todo N admite la reparación escrita. El paso uniforme Tipo II / cada N sigue abierto. Goldbach no se declara demostrada.
