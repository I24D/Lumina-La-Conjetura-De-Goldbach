# 06 — Descomposición bilineal de las cajas B_N(C,K): qué está controlado y qué no

## Procedencia y estado

- **Agente:** Claude Code. **Modelo:** `claude-opus-5-5`.
- **Fecha:** 2026-11-23 (America/New_York), 2026-11-24 UTC.
- **Encargo:** tarea (HB) coordinada por LUMINA en nombre de Dal. Descomponer Λ(N−cd) dentro de B_N(C,K), clasificar los términos, dar referencias por rango, documentar los rangos no controlados y no usar Bombieri–Vinogradov más allá de su nivel.
- **Método:** lectura de README, ESTADO y documentos 01–05; derivación escrita; diagnóstico numérico propio (`numerico/cajas_bilineales.py`), independiente del de OpenClaw.
- **Estado global:** **[P]**, sin revisión independiente y sin formalizar.
  - Los teoremas publicados que se citan son **[T]**. Las referencias marcadas *(verificar)* son de memoria del modelo y hay que contrastarlas antes de citarlas formalmente.
  - Todas las constantes que dependen de Siegel–Walfisz son **inefectivas**.
- **Archivos:** este documento es nuevo. **No** se modificaron 04 ni 05; las observaciones sobre ellos están en §2.

---

## Resumen

1. **Lo demostrado [P].**
   - (HB) se cumple en las cajas con **K ≤ √N/(log N)^B** (Lema 6.1, Bombieri–Vinogradov aplicado directamente a la variable c).
   - En las cajas con **K > √N/(log N)^B**, la identidad de Vaughan separa B_N en partes degeneradas y de Tipo I, todas ≪ N/(log N)^A (Lema 6.2, BV para μ en la variable d), y una parte de **Tipo II, B₃**, que contiene todo el contenido de la caja.
2. **Lo que falta.** La parte B₃ de Tipo II en **todas** las cajas con K > √N/(log N)^B, es decir, c < √N (log N)^{B+1}. Ningún resultado publicado que conozca la controla.
   - **Si C ≥ N^δ:** bastaría una distribución de nivel K > √N en la clase **fija** N mod d, con pesos μ(d) log d; algo del tipo Elliott–Halberstam. **[C: EH]** lo da.
   - **Si C < N^δ:** es el núcleo de paridad, y ni siquiera EH alcanza (criba asintótica de Bombieri).
3. **Observación estructural (Lema 6.0, [T] elemental).** En la caja con c = 1, la parte con d primo es **−Σ_{p∈(K,2K]} log p·Λ(N−p): un bloque diádico de la propia suma de Goldbach R(N)**. Por tanto (HB) en esas cajas tiene la misma dificultad que la asintótica de Goldbach en bloques diádicos. Numéricamente, la cola es la diferencia de dos cantidades de tamaño ≈ 11·N que se cancelan hasta ≈ 0,09·N (N = 10⁷, §9).
4. **Mejor cota demostrable.**
   - Por caja: |B_N(C,K)| ≤ (3,9+o(1))·CK·log N, por Brun–Titchmarsh.
   - Para la cola entera: **−(3+o(1))·N𝔖(N) ≤ T_N^prim ≤ (1+o(1))·N𝔖(N)**.
   - Goldbach para N grande necesita T ≤ (1−κ)·N𝔖(N). La cota superior conocida está **exactamente en el umbral**: la distancia que falta es la barrera de paridad en su forma más limpia.
5. **Lean 4 + Mathlib instalado y funcionando** (§11). La identidad base Λ = −μ·log ya compila contra Mathlib.

---

## 1. Notación

N par y grande; L = log N; A > 0 arbitrario; B, B₁, … son constantes que dependen de A.

    α_c = 1_{(c,N)=1},   β_d = μ(d) log d,
    B_N(C,K) = Σ_{C<c≤2C} α_c Σ_{K<d≤2K, cd<N} β_d Λ(N−cd)          (05, (5.2))

**Convención de cajas.** Para cubrir exactamente la región d > D de 04 se toman K = (D+1)·2^j/2, es decir, d ∈ ((D+1)2^{j−1}, (D+1)2^j], y C = 2^{i−1}. Con la definición literal de (5.2), la caja que contiene D incluiría valores d ≤ D. Es una **precisión** sobre 05, no un error: basta truncar esa caja en d > D. Hay O(L²) cajas.

**Cajas triviales (R0).** Si CK ≤ N/L^{A+3}, entonces |B_N(C,K)| ≤ CK·L² ≤ N/L^{A+1} **[T, elemental]**. Solo importan las cajas con N/L^{A+3} < CK < N, que están a distancia logarítmica de la "diagonal" CK ≍ N.

---

## 2. Revisión de 04 y 05: informes explícitos

**R-1 [confirmación].** La implicación "(HB) en todas las cajas ⇒ T_N^prim = o(N)" de 05 §4 es correcta: hay O(L²) cajas y se aplica la desigualdad triangular *después* de sumar cada caja. Solo falta la precisión de bordes de §1.

**R-2 [mejora; 04 no contiene ningún error].** La Prop. 4.6 de 04 vale con la misma demostración para **D = ⌊√N/(log N)^{B(A)}⌋**, no solo para D = ⌊N^{1/3}⌋:
- En §4.3 b) de 04, BV se usa con módulos d ≤ D, y su nivel admite D ≤ √N/(log N)^{B}.
- En §4.3 a), la evaluación de A_N(D) solo necesita D → ∞ como potencia de N; con D = √N/L^B el error N^ε D^{−1/8} sigue siendo potencia-pequeño.

Revisé los pasos de 04 §§2–4.3 y su estructura me parece correcta. Etiqueta: **[P-rev]**, con revisión no independiente porque los agentes comparten usuario de Windows. Consecuencia: **la cola pendiente se reduce a d > √N/L^B, c < √N·L^{B+1}**. Es lo mismo que demuestra el Lema 6.1 caja por caja.

**R-3 [observación crítica sobre el alcance de (HB); no es un error].** Ver el Lema 6.0. La cola T_N^prim(D) contiene, con c = 1 y d = p primo, la cantidad −Σ_{D<p<N} log p·Λ(N−p). Es *toda* la contribución de los primos m = p > D a R(N), porque en Λ(p) = −Σ_{d|p} μ(d) log d solo contribuye d = p. Por tanto:

    T_N^prim(D) = −R_{>D}^{pr}(N) + T'_N(D),   R_{>D}^{pr}(N) := Σ_{D<p<N} log p·Λ(N−p),

donde T' es la parte de la cola con m = cd compuesto. La identidad exacta (4.1), junto con R = R_{>D}^{pr} + R_{≤D}^{pr} + (potencias propias), da

    T'_N(D) = S_N(D) − R_{≤D}^{pr}(N) − R^{pot}(N) = N𝔖(N) + O_A(N/L^A)     [P-rev, de (4.15)].

Es decir: **la parte compuesta de la cola es incondicionalmente asintótica a N𝔖(N)**, y "T = o(N)" significa literalmente "R_{>D}^{pr}(N) = N𝔖(N) + o(N)", que es la asintótica de Hardy–Littlewood. 05 no presenta (HB) como algo más fácil, así que no hay error. Pero conviene que el equipo sepa que en las cajas con C pequeño **(HB) no es un objetivo más débil que Goldbach: es Goldbach por bloques diádicos**.

**R-4 [N, comprobación cruzada].** Mi implementación (bucle en (c,d), sin la reindexación (5.1)) da **T_N^prim/N = 0,0642495 en N = 10⁶**, idéntico al valor de 05 §3. Las dos implementaciones son independientes en código, aunque no en autoría humana.

---

## 3. Mapa de regiones y qué estimación controla cada una

Coordenadas: γ = log C/L y κ = log K/L. En la región de 04, κ ∈ [1/3, 1] y γ + κ ≤ 1. Solo importan las cajas con γ + κ ≥ 1 − O(log L/L).

| Región | Rango | Término | Tipo | Estimación publicada que lo controla | Estado |
|---|---|---|---|---|---|
| R0 | CK ≤ N/L^{A+3} | toda la caja | — | cota trivial | **[T]** elemental |
| R1 | K ≤ √N/L^{B} (equivale a C ≥ √N·L^{B−1}) | toda la caja, sin descomponer Λ | Tipo I en c: Λ en progresiones mod d | **Bombieri 1965 / A. I. Vinogradov 1965**, con máx. en y y en la clase (Davenport, *Multiplicative Number Theory*, cap. 28); PNT para μ con coprimalidad | **[P]**, Lema 6.1 |
| R2 | K > √N/L^{B} | B₀: Λ_{≤V}, n = N−cd ≤ V | degenerado ("diagonal" cd ≈ N) | cota trivial + τ(m) ≪ m^ε | **[P]**, §6 |
| R2 | ídem | términos no coprimos y potencias de primos que dividen a N | degenerado | 04, Prop. 4.2 | **[P-rev]** |
| R2 | ídem | B₁, B₂: Tipo I de Vaughan, q ≤ UV ≤ K^{1/2}/L^{B₂} | **Tipo I** (variable lisa b) | **Vaughan 1977** (identidad) + **BV para μ** en la variable d, nivel K^{1/2} (gran criba + Siegel–Walfisz para μ; ver Iwaniec–Kowalski, *Analytic Number Theory*, Thm. 17.4 para formas bilineales *(verificar numeración)*) | **[P]**, Lema 6.2 |
| R2-II-a | N^δ ≤ C < √N·L^{B+1} | B₃: Λ_{>V} * (1*μ_{>U}) | **Tipo II** | **Ninguna conocida.** Haría falta nivel K > √N en la clase fija N mod d, con pesos μ(d) log d | **ABIERTO**; [C: EH] |
| R2-II-b | C < N^δ | B₃ | **Tipo II**, núcleo de paridad | **Ninguna.** Ni siquiera EH (Bombieri 1976) | **ABIERTO** |

Elección de parámetros en R2: U = V = N^{1/8}/L^{B₃}, de modo que UV ≤ K^{1/2}/L^{B₂} para todo K ≥ √N/L^B.

---

## 4. Lema 6.0: la parte prima de la caja con c = 1 [T, elemental]

**Lema 6.0.** Para cada K,

    B_N((1/2,1], K) = −Σ_{K<p≤2K, p<N} log p·Λ(N−p) + Σ_{K<d≤2K, d<N, ω(d)≥2} μ(d) log d·Λ(N−d).

**Prueba.** Con C = 1/2, la única c posible es c = 1, y (1, N) = 1. Si μ(d) ≠ 0, d es libre de cuadrados: o d = p y μ(p) log p = −log p, o ω(d) ≥ 2. □

**Lectura.** El primer término es un bloque diádico de R*(N) (Goldbach con pesos). Por Hardy–Littlewood, su tamaño esperado es ≍ 𝔖(N)·K. (HB) en esa caja exige que el segundo término lo compense hasta O(N/L^A). Probarlo exige separar la contribución de los d con ω(d) impar de la de los d con ω(d) par: **es el problema de paridad**. En §9 se ve numéricamente: en la caja c = 1, d ∈ [3,5·10⁶; 7,1·10⁶] con N = 10⁷, las dos partes valen −0,622·N y +0,618·N.

---

## 5. Lema 6.1: región R1, con BV a su nivel real [P]

**Lema 6.1.** Para cada A > 0 existe B = B(A) tal que, para N par suficientemente grande y toda caja con K ≤ √N/L^B,

    B_N(C,K) ≪_A N/L^A     (constante inefectiva).

**Prueba.**
1. **Quitar la coprimalidad.** Si Λ(m) ≠ 0 con m = N − cd y (m, N) = 1, entonces (c, N) = 1 automáticamente: un primo común a c y N divide a cd = N − m, y por tanto a m. Los m con (m, N) > 1 son potencias de primos que dividen a N. Su aporte total, sumando todas las cajas, es ≤ K_N ≪_ε N^ε L³ (04, Prop. 4.2). Lo mismo pasa con (d, N) > 1. Queda, salvo O(N^{1/2+ε}),

        B = Σ_{K<d≤2K, (d,N)=1} μ(d) log d·[ψ(y₂(d); d, N) − ψ(y₁(d); d, N)],

    con y₁(d) = max(0, N−1−2Cd) e y₂(d) = N−1−Cd (truncando en cd < N).
2. **BV a su nivel real.** Con la forma de Davenport, cap. 28:

        Σ_{q ≤ Q} max_{y≤N} max_{(a,q)=1} |ψ(y;q,a) − y/φ(q)| ≪_{A'} N/L^{A'},   Q = √N/L^{B(A')},

    y el máximo en la clase permite a = N mod d aunque N esté fijo. Como 2K ≤ Q, el error total es ≪ L·N/L^{A'}.
3. **Término principal.** Queda Σ_d μ(d) 1_{(d,N)=1} (d/φ(d))·g(d), con g(d) = (log d)(y₂ − y₁)/d. Es una función lineal a trozos en d, de variación ≪ C·L en (K, 2K]. Por sumación parcial basta

        Σ_{d≤y, (d,N)=1} μ(d) d/φ(d) ≪_J y/(log y)^J,   uniformemente en N, para y ≥ N^{1/3}.

    Se escribe μ·1_{(·,N)=1}·id/φ = μ * h_N, donde:
    - para p ∤ N, h_N(p^k) = −1/(p−1);
    - para p | N, h_N(p^k) = 1.

    Se parte en b ≤ √y y b > √y, como en 04 §4.3 a):
    - en el primer rango, M(x) ≪_J x/(log x)^J y Σ|h_N(b)|/b ≪ log log N;
    - en el segundo, Σ_{b>√y}|h_N(b)|/b ≤ y^{−1/4} ∏_{p|N}(1−p^{−3/4})^{−1}·O(1) = y^{−1/4} N^{o(1)}.

    Así el término principal es ≪ CK·L/L^{J−1} ≤ N/L^{J−2}. □

**Observación.** El Lema 6.1 **no** usa la identidad de Vaughan: c ya es una variable lisa larga. Coincide con R-2: la región R1 forma parte del lado "S" de la hipérbola.

---

## 6. Descomposición de Vaughan en R2 y Lema 6.2 (Tipo I) [P]

**Identidad de Vaughan (1977) [T].** Para U, V ≥ 1 y n ≥ 1,

    Λ(n) = Λ_{≤V}(n) + Σ_{ab=n, a≤U} μ(a) log b − Σ_{aeb=n, a≤U, e≤V} μ(a)Λ(e)
           + Σ_{ab=n, a>V} Λ(a)·ν_U(b),      ν_U(b) = Σ_{t|b, t>U} μ(t).

Se tiene ν_U(1) = 0, ν_U(b) = 0 para 1 < b ≤ U, y |ν_U(b)| ≤ τ(b).

Aplicándola a n = N − cd dentro de (5.2), **con C y K fijos**:

    B_N(C,K) = B₀ + B₁ − B₂ + B₃,
    B₀ = Σ_{c,d} α_c β_d Λ_{≤V}(N−cd),
    B₁ = Σ_{a≤U} μ(a) Σ_{c~C} α_c Σ_{d~K, cd<N, a | N−cd} β_d log((N−cd)/a),
    B₂ = Σ_{q≤UV} λ(q) Σ_{c~C} α_c Σ_{d~K, cd<N, q | N−cd} β_d,   λ(q) = Σ_{ae=q, a≤U, e≤V} μ(a)Λ(e),
    B₃ = Σ_{c~C} α_c Σ_{d~K, cd<N} β_d Σ_{ab=N−cd, a>V, b>U} Λ(a) ν_U(b).

Se tiene |λ(q)| ≤ log q.

**B₀ (degenerado).** Λ_{≤V}(N−cd) ≠ 0 obliga a N−V ≤ cd < N. Por tanto |B₀| ≤ L² Σ_{N−V≤m<N} τ(m) ≪_ε V·N^ε, que es ≪ N^{1/8+ε}. **[P]**

**Lema 6.2 (Tipo I).** Supongamos K ≥ √N/L^B y UV ≤ K^{1/2}/L^{B₂(A)}. Entonces |B₁| + |B₂| ≪_A N/L^A.

**Prueba.**
- Si q | N − cd y (c, N) = 1, entonces (q, c) = 1, porque un factor común dividiría a N. La condición equivale a **d ≡ N·c̄ (mod q)**.
- Para c y q fijos, el peso en d (log d, o log d·log((N−cd)/q)) es monótono a trozos y tiene variación ≪ L². Por sumación parcial,

        |suma interior| ≪ L²·max_{y≤2K} max_{r mod q} |Σ_{K<d≤y, d≡r (q)} μ(d)|.

- **BV para μ [T],** con la misma demostración que para Λ (gran criba, Siegel–Walfisz para μ e identidad de Vaughan para μ):

        Σ_{q≤K^{1/2}/L^{B₂}} max_r max_y |Σ_{d≤y, d≡r(q)} μ(d)| ≪_{A'} K/L^{A'}.

  Las clases no reducidas, (r, q) = g > 1, se reducen escribiendo d = g d′ con (d′, g) = 1. Eso da clases de módulo ≤ q con una condición de coprimalidad, que se absorbe con un divisor de q **[P, detalle rutinario]**. A diferencia de Λ, μ **no** tiene término principal que restar.
- Sumando sobre c ~ C y q ≤ UV: |B₁| + |B₂| ≪ C·L³·K/L^{A'} ≤ N/L^{A'−3}. □

**Consecuencia [P].** En R2, **B_N(C,K) = B₃ + O_A(N/L^A)**: toda la caja está en el término de Tipo II.

**¿Mejora algo la identidad de Heath-Brown (1982)?** *(Heath-Brown, Canad. J. Math. 34 (1982), 1365–1377. La fecha correcta es 1982, no 1978.)* No:
- produce además términos con **dos** variables lisas (Tipo I₂). Aquí equivalen a correlaciones Σ_d μ(d)·τ((N−cd)/q), que con la técnica de hipérbola exigen BV para μ con módulos de hasta √(qN) > K^{1/2}. Es decir, necesitan **más** que BV;
- los términos de Tipo II que deja siguen siendo de la misma naturaleza que B₃.

Vaughan es la elección correcta aquí.

---

## 7. Por qué B₃ (Tipo II) no está controlado en R2

### 7.1 La única variable lisa es c
En B₃, las variables a (Λ), b (ν_U) y d (μ) son ásperas. **Solo c es lisa.** Hay dos posibilidades:
- **Sumar en c.** Convierte B₃ en la distribución de la sucesión de Tipo II (Λ_{>V} * ν_U)(m) en la progresión **m ≡ N (mod d)**, con d ~ K > √N/L^B y pesos con signo μ(d) log d. Eso es **nivel de distribución mayor que 1/2 en una clase fija**.
- **Sumar en cualquier otra variable.** Si se fijan dos de las variables ásperas, la tercera queda ligada a la cuarta por una forma lineal, y reaparece una **correlación binaria**: por ejemplo Σ_d μ(d) log d·Λ((N−cd)/b) con d ≡ N c̄ (mod b). Es del mismo tipo que el problema original, y de paridad.

### 7.2 Qué dan los resultados publicados más allá de 1/2, y por qué no se aplican
| Resultado | Qué da | Por qué no basta aquí |
|---|---|---|
| BV para formas bilineales (Motohashi 1976; Bombieri–Friedlander–Iwaniec 1986, Thm. 0) [T] | nivel √N/L^B, con valores absolutos | K > √N/L^B por definición de R2 |
| Bombieri–Friedlander–Iwaniec I (Acta Math. 156, 1986), resultados de nivel 4/7 [T] | más allá de 1/2 en clase fija, con pesos **bien factorizables** y sin valores absolutos | μ(d) log d no es bien factorizable; la pieza con d primo (§7.3) es la antítesis |
| BFI II–III (1987, 1989) *(verificar enunciados)* | ganancias del orden de potencias de log por encima de √N, con valores absolutos | como mucho extienden R1 en factores logarítmicos; no alcanzan K = N^{1/2+δ} |
| Zhang (Annals 2014); Polymath 8a (2014) [T] | nivel 1/2 + 7/300 con valores absolutos para módulos **densamente divisibles** (lisos) | cubre solo la parte de d sin factor primo grande y K ≤ N^{1/2+7/300}; ni los d primos ni los K mayores |
| Maynard, *Primes in APs to large moduli I–III* (2020–) *(verificar)* | nivel hasta 3/5 con pesos triplemente bien factorizables, y otros casos especiales | misma restricción estructural sobre los pesos del módulo |
| **[C: EH]** Elliott–Halberstam con clase fija, valores absolutos y nivel N^{1−ε} | controlaría B₃ en R2-II-a (C ≥ N^δ) mediante la suma en c y la ausencia de término principal de μ | es una hipótesis, no un teorema; y **no** cubre R2-II-b |

### 7.3 La pieza con módulo primo (presente en cada caja de R2)
Como μ(d) log d = −(μ*Λ)(d), la parte de β con d = p primo da

    P_N(C,K) = −Σ_{c~C, (c,N)=1} Σ_{p~K} log p·Λ(N−cp).

Cuenta las representaciones **N = q + c·p** con p > √N/L^B. Su tamaño esperado es ≍ CK·𝔖-promedio ≍ N: **no es pequeña** (§9: entre −0,28·N y −0,73·N en las cajas mayores). Evaluarla asintóticamente exige primos en la clase N mod p para **módulos primos p > √N**, sumados con signo. Es un enunciado de tipo Elliott–Halberstam para módulos primos, que no está demostrado. Ningún resultado bien factorizable se aplica a módulos primos.

### 7.4 El núcleo de paridad: C < N^δ
Para c ≤ N^δ hay como mucho N^δ términos por clase módulo d: ninguna afirmación de distribución, ni siquiera EH, aporta información. La **criba asintótica de Bombieri** (Rend. Accad. Naz. XL, 1976, *verificar*) muestra que la información de nivel de distribución no fija Σ μ(d)·(…) hasta o(N). Por el Lema 6.0, la caja c = 1 contiene un bloque diádico de Goldbach. **Esta es la ubicación exacta de la barrera de paridad en el marco de 04/05.**

### 7.5 Un intento por dispersión no gana nada
Si se aplica Cauchy–Schwarz en d al error de distribución Δ(d) de la sucesión de Tipo II en m ≡ N (mod d):
- el término diagonal (m₁ = m₂) vale ≍ N²L²/K;
- con Σ_d |β_d|² ≍ K L², se obtiene Σ_d |β_d||Δ(d)| ≪ N L², que es del orden de la cota trivial.

Una versión de Barban–Davenport–Halberstam con **clase fija** que mejore esto no está demostrada.

---

## 8. Mejores cotas demostrables en la zona no controlada

**8.1 Por caja [T + P].** Para c fijo, los n = N − cd con d ∈ (K, 2K] recorren una progresión mod c de longitud cK. Brun–Titchmarsh (Montgomery–Vaughan 1973: π(x+y;q,a) − π(x;q,a) ≤ 2y/(φ(q) log(y/q))) da a lo sumo 2cK/(φ(c) log K) primos. Las potencias propias aportan en total ≪ N^{1/2+ε}. Con log d ≤ log 2K, log N/log K ≤ 2 + o(1) en R2, y Σ_{C<c≤2C} c/φ(c) ≤ (ζ(2)ζ(3)/ζ(6) + o(1))·C ≈ 1,944·C (Landau):

    |B_N(C,K)| ≤ (3,9 + o(1))·CK·L + O_ε(N^{1/2+ε}).

En la diagonal CK ≍ N esto es ≍ N L. **(HB) pide N/L^A: la distancia es un factor L^{A+1}.**

**8.2 Para la cola completa [P-rev + T].** De (4.15), R(N) ≥ 0 y la cota superior de criba R(N) ≤ (4 + o(1))·N𝔖(N) (Selberg; mejorada a ≈ 3,92 por Chen 1978 y J. Wu 2004, *verificar*):

    −(3 + o(1))·N𝔖(N) ≤ T_N^prim(D) ≤ (1 + o(1))·N𝔖(N).

Goldbach para N grande equivale a T_N^prim ≤ (1−κ)·N𝔖(N) (04 §5.2; (HB+) de 05). **La cota superior demostrable está exactamente en el umbral 1, y hace falta llegar a 1−κ.** La barrera de paridad aparece aquí como una constante, no como una potencia de log.

---

## 9. Datos numéricos [N] (evidencia, no demostración)

Script: `numerico/cajas_bilineales.py`. Resultados: `numerico/resultados/20261123_215414_cajas_bilineales.json` (N = 10⁵ y 10⁶) y `20261123_215627_cajas_bilineales.json` (N = 10⁷). D = ⌊N^{1/3}⌋, como en 04/05. Cajas: c ∈ [2^i, 2^{i+1}), d ∈ [(D+1)2^j, (D+1)2^{j+1}). P es la parte con d primo y Q la parte con d compuesto (libre de cuadrados, ω ≥ 2).

| N | T^prim/N | T/(N𝔖) | ΣP/N | ΣQ/N | max_cajas \|B\|/N (P/N, Q/N en esa caja) |
|---:|---:|---:|---:|---:|---|
| 10⁵ | −0,01741 | −0,00989 | −7,756 | +7,739 | 0,0633 (−0,215; +0,152) |
| 10⁶ | +0,06425 | +0,03650 | −9,360 | +9,424 | 0,0337 (−0,255; +0,288) |
| 10⁷ | +0,08625 | +0,04899 | −10,913 | +10,999 | 0,0181 (−0,277; +0,295) |

**Fila c = 1 con N = 10⁷.** La suma de P/N en todas sus cajas es ≈ −1,761, frente a 𝔖(N) = 1,7604. **La parte prima de la fila c = 1 reproduce −R(N)/N**, como predice el Lema 6.0. Caja mayor: d ∈ [3 538 944; 7 077 887], con P/N = −0,6221, Q/N = +0,6178 y B/N = −0,0043.

**Lectura.**
1. P y Q crecen como ≍ L, igual que la cota trivial, y se cancelan casi por completo.
2. El máximo de |B|/N por caja **disminuye** con N (0,063 → 0,034 → 0,018). Encaja con (HB), pero no la prueba.
3. T/(N𝔖) es positivo y del orden de 0,05 para N ≤ 10⁷. Con estos N, los errores de la reducción de 04 (truncación en D = N^{1/3}) no son despreciables, así que estos datos **no** deciden si T/N → 0.

---

## 10. Qué está demostrado y qué falta

**Demostrado en este documento** (todo **[P]**, sin revisión independiente):
- Lema 6.0: estructura de la caja c = 1 ([T] elemental).
- Lema 6.1: (HB) en R1, K ≤ √N/L^B.
- B₀ y los términos no coprimos son despreciables.
- Lema 6.2: Tipo I de Vaughan en R2, con UV ≤ K^{1/2}/L^{B₂}.
- Equivalencia "(HB) en R2 ⟺ B₃ ≪ N/L^A".
- Cotas de §8.
- Informes R-1 a R-4 sobre 04 y 05.

**Falta (abierto):**
- **F1 (R2-II-a, N^δ ≤ C < √N L^{B+1}).** Una cota B₃ ≪ N/L^A. Basta [C: EH en la clase fija N mod d, con nivel K]. No hay resultado publicado que la dé; en particular falta la pieza P_N con módulo primo.
- **F2 (R2-II-b, C < N^δ).** Núcleo de paridad: contiene bloques diádicos de Goldbach (Lema 6.0). No lo alcanza ninguna técnica conocida, ni siquiera bajo EH.
- **F3 (efectividad).** Todas las cotas que dependen de Siegel–Walfisz son inefectivas. Para pasar de "N suficientemente grande" a "todo N" haría falta además un N₀ explícito.

**Recomendaciones para el equipo [P]:**
1. **Dejar de perseguir (HB) caja por caja en R2.** Por R-3, en las cajas con C pequeño equivale a Goldbach por bloques, así que no es un objetivo intermedio.
2. El blanco honesto es **la cota unilateral (HB+) para la suma conjunta sobre C < N^δ**, es decir, mejorar la constante 1 de §8.2 a 1−κ. Cualquier ataque tiene que traer información que distinga ω(d) par de ω(d) impar. La única fuente conocida de algo así en problemas binarios es el análogo de cuerpos de funciones (doc 03, Pellet).
3. **Formalizar en Lean:** la Prop. 4.1 de 04 (identidad exacta), el Lema 6.0 y la Prop. 4.2. Son finitos y elementales, y ya se pueden escribir (§11).

---

## 11. Lean 4 + Mathlib (instalado por Claude a petición de Dal, 2026-11-23)

- elan 4.4.0 en `C:\Users\dal_n\.elan`, añadido al PATH del usuario.
- Proyecto Mathlib en **`C:\Users\dal_n\lean\goldbach_lumina`**, fuera de OneDrive para no sincronizar unos 5 GB de caché. Toolchain `leanprover/lean4:v4.37.0-rc1`. `lake exe cache get` bajó los 9133 archivos y `lake build` terminó con éxito.
- Prueba de humo `GoldbachLumina/Prueba.lean`, compilada con exit 0:
  - el enunciado `GoldbachFuerte`, **sin demostrar**;
  - el caso 100 = 3 + 97, comprobado por el kernel;
  - Λ(n) = −Σ_{d|n} μ(d) log d, mediante `ArithmeticFunction.sum_moebius_mul_log_eq` de Mathlib.
- Uso: `cd C:\Users\dal_n\lean\goldbach_lumina` y luego `lake env lean <archivo>.lean`.
