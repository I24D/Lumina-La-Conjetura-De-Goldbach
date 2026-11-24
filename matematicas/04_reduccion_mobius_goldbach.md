# 04 — Reducción explícita a una cola primitiva de Möbius–von Mangoldt

## Procedencia y estado

- Agente: Codex. Modelo exacto: no expuesto por la interfaz.
- Fecha: 2026-11-23, America/New_York.
- Método: derivación paso a paso, comprobaciones simbólicas finitas y consulta de fuentes primarias.
- Estado: **[P]**. Pruebas escritas no revisadas independientemente y no formalizadas. Las identidades clásicas no se presentan como descubrimientos. §4.2 obtiene un error en potencias usando GRH; **§4.3 elimina esa dependencia para la reducción asintótica** mediante resultados clásicos incondicionales. La certificación externa de GRH se distingue de ejecución local en 03.
- Fuentes: [Sawin–Shusterman, Cor. 6.1 y §6.3](https://arxiv.org/pdf/1808.04001v2); [Philippe Michel, notas de Park City, §1.2 y §3.1, Bombieri–Vinogradov](https://www.epfl.ch/labs/tan/wp-content/uploads/2018/10/Parkcitylectures.pdf); [Lynnelle Ye, introducción: RH y la función de Mertens](https://arxiv.org/abs/1406.7326); [D. Mastrostefano, tesis, §7.1.1: cota incondicional de M](https://wrap.warwick.ac.uk/id/eprint/167258/1/WRAP_Theses_Mastrostefano_2021.pdf); documentos 01–03 del proyecto.
- Alcance: aislar una desigualdad pendiente que bastaría para Goldbach para pares suficientemente grandes. No se demuestra esa desigualdad ni se obtiene un umbral efectivo para cubrir todos los pares.

## 1. Definiciones

N≥4 es par; 1≤D≤N−1 es entero; x=N−1. Los logaritmos son naturales. μ es Möbius, Λ(1)=0 y Λ(p^k)=log p.

    R(N)=Σ_{n=1}^{N−1} Λ(n)Λ(N−n),
    R*(N)=Σ_{p+q=N; p,q primos} log p log q,
    ψ(x;d,a)=Σ_{1≤m≤x; m≡a (mod d)} Λ(m),
    A_N(D)=−Σ_{d≤D; (d,N)=1} μ(d)log d/φ(d),
    𝔖(N)=2 C₂ ∏_{p|N; p>2}(p−1)/(p−2),
    C₂=∏_{p>2}(1−1/(p−1)²)>0.

**[T, clásico]** La positividad de C₂ se sigue de la convergencia de Σ_{p>2}1/(p−1)² y de que cada factor es positivo. Por tanto 𝔖(N)≥2C₂>0 uniformemente en N par. Esta positividad por sí sola no cuenta representaciones.

## 2. Identidad exacta de hipérbola [P, demostración algebraica]

Definimos

    S_N(D)=−Σ_{d≤D} μ(d)log d · ψ(N−1;d,N),
    T_N(D)=Σ_{c≤(N−1)/(D+1)} Σ_{D<d≤(N−1)/c}
                         μ(d)log d · Λ(N−cd).

**Proposición 4.1.** Para todos esos N,D,

    R(N)=S_N(D)−T_N(D).                                  (4.1)

**Prueba.** La identidad clásica Λ(n)=−Σ_{d|n}μ(d)log d da

    R(N)=−Σ_{cd≤N−1} μ(d)log d · Λ(N−cd).

Separamos d≤D y d>D, sin aproximaciones. Para d fijo, los enteros N−cd positivos son exactamente los m<N congruentes con N módulo d; esto produce S_N(D). En el segundo dominio d≥D+1, de modo que c≤(N−1)/(D+1); es T_N(D), con signo menos. Los términos d=1 y N−cd=1 se anulan como deben. □

La notación T_N(D) de esta nota es una **suma escalar**; no es el operador T_N de Li–Liu.

## 3. Los términos no coprimos son pequeños [P, prueba incondicional]

Sea

    K_N = log N · Σ_{p|N} log p · Σ_{k≥1; p^k<N} τ(N−p^k),

donde τ(m) cuenta divisores positivos. Sea T_N^prim(D) la misma suma de T_N(D), restringida a (c,N)=1, y T_N^bad(D) la complementaria.

**Proposición 4.2.**

    |T_N^bad(D)| ≤ K_N,
    |Σ_{d≤D; (d,N)>1} μ(d)log d · ψ(N−1;d,N)| ≤ K_N,
    K_N ≪_ε N^ε (log N)³,  para todo ε>0.                 (4.2)

**Prueba.** Si (c,N)>1 y Λ(N−cd)≠0, escribimos N−cd=p^k. Todo primo común a c y N divide p^k, así que su base es p y p|N. Lo mismo ocurre con (d,N)>1. Para cada p^k<N, los pares positivos cd=N−p^k son exactamente τ(N−p^k). Cada término tiene valor absoluto a lo sumo log N·log p. Quitar restricciones sobre d y c solo aumenta esta cota positiva. Esto prueba las dos desigualdades.

Finalmente τ(m)≪_ε m^ε, y

    Σ_{p|N} log p · #{k≥1:p^k<N}
       ≤ ω(N)log N ≤ (log N)²/log 2.

La cota clásica de τ es elemental: para primos suficientemente grandes, k+1≤p^(εk); para los primos pequeños se absorben los supremos de (k+1)/p^(εk) en una constante finita. □

Así, la cancelación pendiente solo necesita tratar **(c,N)=1**. No necesitamos trasladar una estimación indiscriminada que incluya progresiones degeneradas.

## 4. Término de divisores pequeños [C: GRH-Dirichlet]

Usamos dos consecuencias clásicas y distinguibles de GRH:

    (AP)  ψ(x;d,a)=x/φ(d)+O(√x·log²(dx)), si (a,d)=1,
    (M)   M(y)=Σ_{n≤y}μ(n) ≪_ρ y^(1/2+ρ), para todo ρ>0.

Las constantes de AP son uniformes en d,a; la de M depende de ρ. AP figura en la aplicación a progresiones de las notas de Michel. M es la equivalencia clásica de Littlewood recogida en la introducción de Ye. No se sustituye una función Λ(N−cd) variable por un carácter fijo al invocar estas consecuencias.

De AP y (4.2), para D≤N,

    S_N(D)=(N−1)A_N(D)
           +O(D√N·log³N + K_N).                         (4.3)

En efecto, cada error coprimo se multiplica por |μ(d)|log d≤log N y hay a lo sumo D términos. Los demás son precisamente los controlados en (4.2).

### 4.1 Evaluación uniforme de A_N(D), con la dependencia en N expuesta

**Proposición 4.3 [C: M].** Para 0<ρ<1/2 y ε>0,

    A_N(D)=𝔖(N)+O_{ρ,ε}(N^ε D^(−1/2+ρ)log(2D)),       (4.4)

uniformemente en N par y D≥1.

**Prueba.** Definimos una función multiplicativa h_N por h_N(1)=1 y, para k≥1,

    h_N(p^k)=1/p^k                       si p|N,
    h_N(p^k)=−1/((p−1)p^k)               si p∤N.

La identidad de coeficientes de Euler, verificable prima por prima, es

    1_{(n,N)=1} μ(n)/φ(n)
       = Σ_{ab=n} μ(a)/a · h_N(b).                      (4.5)

Para p∤N, el factor generador de h_N es

    (1−p^(−s)/(p−1))/(1−p^(−1−s));

para p|N es (1−p^(−1−s))⁻¹. En particular,

    Σ_b h_N(b)=∏_{p|N}(1−1/p)⁻¹
                ·∏_{p∤N}(1−1/(p−1))/(1−1/p)
              =𝔖(N).                                  (4.6)

La serie es absolutamente convergente. Más precisamente,

    H_abs(N,ρ)=Σ_b |h_N(b)| b^(1/2−ρ) ≪_{ρ,ε}N^ε.      (4.7)

Para los primos fuera de N los factores son 1+O_ρ(p^(−3/2−ρ)), cuyo producto converge. Para p|N se reemplazan por (1−p^(−1/2−ρ))⁻¹. Para p suficientemente grande cada uno es ≤p^ε; el producto de los restantes se absorbe en una constante, y ∏_{p|N}p≤N. Así no se esconde una constante arbitraria dependiente de N.

Por sumación parcial de M, y el polo de ζ con residuo 1 en s=1,

    Σ_{a≤y} μ(a)/a = O_ρ(y^(−1/2+ρ)),
    Σ_{a≤y} μ(a)log a/a = −1+O_ρ(y^(−1/2+ρ)log(2y)).   (4.8)

Para justificar los valores límite: las series y su derivada convergen por M; para s>0 corresponden a 1/ζ(1+s) y su derivada, y 1/ζ(1+s)=s+O(s²) cuando s baja a cero.

Insertamos (4.5) en la suma que define −A_N(D). Para cada b≤D, las sumas interiores de (4.8), con y=D/b, dan

    Σ_{a≤D/b} μ(a)/a · (log a+log b)
       = −1+O_ρ((D/b)^(−1/2+ρ)log(2D)).

El error total queda acotado por D^(−1/2+ρ)log(2D)H_abs(N,ρ). Reemplazar Σ_{b≤D}h_N(b) por (4.6) cuesta a lo sumo D^(−1/2+ρ)H_abs(N,ρ). Usando (4.7) obtenemos (4.4). □

### 4.2 Forma final de la parte controlada

Combinando (4.1)–(4.4),

    R(N)=N𝔖(N)−T_N^prim(D)
      +O_{ρ,ε}(N^(1+ε)D^(−1/2+ρ)log(2D)
                 +D√N log³N+N^ε log³N).                (4.9)

Cambiar N−1 por N cuesta O_ε(N^ε): (4.6)–(4.7) también acotan 𝔖(N).

Con D=⌊N^(1/3)⌋, los exponentes de las dos cotas se equilibran en 5/6. Para cualquier ν>0 suficientemente pequeño, elegimos ε=ν/4 y ρ=3ν/4 y absorbemos los logaritmos restantes en N^ν:

    R(N)=N𝔖(N)−T_N^prim(⌊N^(1/3)⌋)
                      +O_ν(N^(5/6+ν)).                (4.10)

**[C: GRH-Dirichlet]** Esta fórmula sí tiene una derivación en esta nota. Tomar ν<1/6 hace o(N) el error escrito. **No** dice que la cola sea pequeña ni que 5/6 sea una cota del conjunto excepcional.

### 4.3 Retirar GRH de la reducción asintótica [P, derivación con entradas clásicas incondicionales]

Esta mejora evita depender del estado de la noticia de Alpöge para aislar la cola. No conserva el error en potencias de (4.10).

**Entradas clásicas [T].** La región libre de ceros de ζ da, para cada J>0, M(y)≪_J y/(log(2y))^J. La cota de Vinogradov–Korobov, recogida en §7.1.1 de la tesis citada, es más fuerte. Bombieri–Vinogradov da, para cada J>1 y un B_J, un error agregado ≪_J x/(log x)^J en módulos hasta √x/(log x)^B_J. Se invocan estos resultados publicados; no se vuelve a demostrar aquí el teorema analítico completo.

**a) Evaluar A_N sin RH.** La primera entrada y sumación parcial dan, para cada J,

    Σ_{a≤y} μ(a)/a = O_J((log(2y))^(−J)),
    Σ_{a≤y} μ(a)log a/a = −1+O_J((log(2y))^(−J)).

Los valores límite se justifican como en (4.8). Además, por los factores de h_N,

    Σ_b |h_N(b)| ≤ C ∏_{p|N}(1−1/p)⁻¹ ≤ C(ω(N)+1)≪log(2N).

La desigualdad del medio es elemental: si p_1<…<p_r son los divisores primos de N, p_i≥i+1 y el producto de p_i/(p_i−1) es ≤∏(i+1)/i=r+1.

Dividimos (4.5) en b≤√D y b>√D. En el primer rango D/b≥√D y las dos cotas anteriores hacen el error interior O_J((log(2D))^(−J)); el factor log b se absorbe eligiendo previamente una potencia J mayor. Su suma cuesta O_J(log(2N)/(log(2D))^J).

En el segundo rango se usa la cota trivial O(log²(2D)) para la suma interior. Por (4.7) con ρ=1/4,

    Σ_{b>√D}|h_N(b)|≪_ε N^ε D^(−1/8).

Esto también cubre la cola de (4.6). Así,

    A_N(D)=𝔖(N)+O_J(log(2N)/(log(2D))^J)
                    +O_ε(N^ε D^(−1/8)log²(2D)).

Para D=⌊N^(1/3)⌋, tomamos ε=1/48. El último error es O(N^(−1/48)log²N); por tanto, para cada A>0,

    A_N(D)=𝔖(N)+O_A((log N)^(−A)).                     (4.14)

**b) Evaluar S_N mediante promedio en módulos, manteniendo cada N fijo.** Para cada J, D≤√(N−1)/(log(N−1))^B_J cuando N es suficientemente grande. Bombieri–Vinogradov incluye un máximo en la clase reducida, de modo que puede usarse la clase a=N mod d aunque N varíe. Multiplicar por |μ(d)|log d cuesta a lo sumo log N. Resulta E_N(D)≪_A N/(log N)^A para cualquier A, después de aumentar J.

La formulación de Michel centra en ψ(x)/φ(d), en vez de x/φ(d). Pasar a x cuesta lo mismo: PNT con ahorro de cualquier potencia logarítmica y Σ_{d≤D}1/φ(d)≪log(2D). Esta última cota se obtiene de

    1/φ(n)=n⁻¹Σ_{a|n} μ(a)²/φ(a),

sumando la serie armónica y usando la convergencia de ∏_p(1+1/(p(p−1))). Así el cambio de centro no introduce una pérdida polinómica oculta.

**Proposición 4.6 [P, sin premisa GRH].** Para D=⌊N^(1/3)⌋ y cada A>0,

    R(N)=N𝔖(N)−T_N^prim(D)+O_A(N/(log N)^A),           (4.15)

uniformemente para N par suficientemente grande. Se deduce de (4.1), (4.2), (4.14) y b). Las constantes de Bombieri–Vinogradov en esta versión no se hacen efectivas; no se proporciona un umbral numérico.

El promedio utilizado es **en d**, con un máximo sobre a. No es un promedio sobre N, y por ello no se ha cambiado «cada N» por «casi todo N».

## 5. Una estimación adicional suficiente y falsable [C: HμΛ]

Para L_c=⌊(N−1)/c⌋ y D≤t≤L_c ponemos

    F_{N,c,D}(t)=Σ_{D<d≤t} μ(d)Λ(N−cd),
    U_N(D)=Σ_{c≤(N−1)/(D+1); (c,N)=1}
                        sup_{D≤t≤L_c}|F_{N,c,D}(t)|.

Los supremos son sobre una función escalonada finita. Sumación parcial da, para cada c,

    Σ_{D<d≤L_c} μ(d)Λ(N−cd)log d
       = F_{N,c,D}(L_c)log L_c − ∫_D^{L_c}F_{N,c,D}(t)dt/t.

Por tanto

    |T_N^prim(D)| ≤ 2 log N · U_N(D).                   (4.11)

**Hipótesis adicional HμΛ [H, NO DEMOSTRADA].** Con D=⌊N^(1/3)⌋,

    U_N(D)=o(N/log N), uniformemente para N par→∞.       (HμΛ)

**Proposición 4.4 [C: HμΛ].** Bajo HμΛ y los resultados clásicos incondicionales usados en §4.3,

    R(N)=N𝔖(N)+o(N),

y todos los pares suficientemente grandes son suma de dos primos.

**Prueba.** (4.11) hace o(N) la cola de (4.15), cuyo error también es o(N). Alternativamente, GRH permite usar (4.10) con ν<1/6. Para separar potencias de primos de primos, sin recurrir a la constante aún no verificada de 01, usamos la cota elemental

    0≤R(N)−R*(N)≤2⌊√N⌋⌊log₂N⌋(log N)²=o(N).          (4.12)

En efecto, cada potencia propia de primo ≤N es algún a^k con 2≤a≤√N y 2≤k≤log₂N; esta enumeración permite sobreconteo. Cada sumando tiene peso ≤(log N)² y hay dos posiciones posibles de la potencia propia. Como 𝔖(N)≥2C₂>0, resulta R*(N)>0 para N par suficientemente grande. □

### 5.1 Una hipótesis puntual más fuerte, con todos los cuantificadores

Fijemos η∈(0,1), B≥0 y K>0. Bastaría que existiera N₀ tal que, **para cada N par≥N₀**, cada c≤(N−1)/(D+1) con (c,N)=1 y cada D≤t≤L_c,

    |F_{N,c,D}(t)|≤K (N/c)^(1−η)(log N)^B.              (Hpuntual)

Aquí η,B,K,N₀ deben ser independientes de N,c,t. Sumando c^(η−1),

    U_N(D)≪_{η}K N D^(−η)(log N)^B,
    |T_N^prim(D)|≪_{η}K N^(1−η/3)(log N)^(B+1)=o(N).

Un resultado para c fijo con constantes sin controlar, para casi todo N o para un desplazamiento fijo no acredita Hpuntual. HμΛ es más débil: permite excepciones en c si su contribución agregada es pequeña.

### 5.2 Para positividad basta aún menos, pero no se da por probado

Si existiese κ>0 fijo tal que

    T_N^prim(D)≤(1−κ)N𝔖(N)

para cada N par suficientemente grande, (4.15) y (4.12) también darían Goldbach en ese rango. Los valores negativos grandes de la cola no perjudican la cota inferior. Esta condición es una reformulación suficiente del obstáculo: **no se ofrece como una solución ni como hipótesis más fácil por definición**.

### 5.3 Balance de exponentes para una transferencia más fiel [P / C: entradas indicadas]

El Corolario 6.1 de Sawin–Shusterman no cubre coeficientes de tamaño arbitrario respecto de la variable: exige m≤αd, α<1, para sus grados. En nuestra elección D=N^(1/3), c puede llegar a N^(2/3) y los d de la cola pueden ser solo de tamaño N^(1/3). La semejanza de las sumas **no** autoriza a aplicar literalmente ese corolario.

**Proposición 4.5, balance de tamaños.** Sea 0<σ<1/2 y D=⌊N^(1/2+σ)⌋. Si d>D y cd≤N−1, entonces

    c<N^(1/2−σ),
    c<d^α,   α=(1−2σ)/(1+2σ)<1.                       (4.13)

**Prueba.** d>D implica d>N^(1/2+σ), incluso con el redondeo. De cd<N se obtiene c<N/d<N^(1/2−σ). Elevando la primera desigualdad a α, d^α>N^(1/2−σ)>c. □

Por ejemplo, σ=1/10 corresponde a D=N^(3/5), c<N^(2/5) y α=2/3. Esta es solo la compatibilidad de tamaños: permanecen las otras hipótesis de las formas y de los conductores.

El coste está en el otro lado de la hipérbola. Definamos el error agregado

    E_N(D)=Σ_{d≤D; (d,N)=1} |μ(d)|log d
                ·|ψ(N−1;d,N)−(N−1)/φ(d)|.

La identidad exacta y las pruebas de §§3–4 dan

    R(N)=N𝔖(N)−T_N^prim(D)
      +O_{ρ,ε}(N^(1+ε)D^(−1/2+ρ)log(2D)+E_N(D)+K_N+𝔖(N)).

Con esa elección de D, bastarían **ambas** entradas:

1. E_N(D)=o(N), para cada N par: distribución agregada más allá de 1/2.
2. U_N(D)=o(N/log N), o una cota unilateral suficiente para T_N^prim(D): cancelación sensible a la paridad.

AP individual solo produce aquí E_N(D)≪N^(1+σ)log³N, que no proporciona la primera entrada. Una hipótesis del tipo Elliott–Halberstam con suficiente ahorro logarítmico sí suministraría ese error agregado, pero no se presupone demostrada y **no sustituye** la segunda entrada.

Las variantes publicadas de distribución más allá de 1/2 para una clase **fija** o pesos **bien factorizables** tampoco se aplican automáticamente: aquí la clase es N mod d y los pesos vienen dados por μ(d)log d. Esas compatibilidades se deben demostrar, no asumir.

**Conclusión operativa.** Hay dos rutas explícitas: controlar la cola amplia con D<N^(1/2), o aportar distribución adicional y controlar una cola de coeficiente relativamente pequeño con D>N^(1/2). Esta última reproduce mejor el reparto de trabajo del método en polinomios. No se ha demostrado ninguna de las nuevas estimaciones requeridas sobre los enteros.

## 6. Intento de cerrar la cola: dónde se detiene la prueba [P]

No podemos aplicar M directamente a μ(d)Λ(N−cd). Para un peso b_d, sumación parcial utiliza la variación de b_d, no solo su tamaño. Λ(N−cd) no es un peso suave de variación pequeña.

Un contraejemplo general al atajo «M pequeño implica cancelación contra cualquier peso acotado no negativo» es b_d=(1+μ(d))/2:

    Σ_{d≤y} μ(d)b_d = (M(y)+Σ_{d≤y}μ(d)²)/2.

Como Σ μ(d)²∼6y/π², el lado derecho es ∼3y/π² aunque M(y)=o(y). Esto no afirma que Λ(N−cd) sea ese peso; identifica exactamente el uso inválido de sumación parcial sin control de variación.

Tampoco basta Cauchy–Schwarz: usando |μ|≤1 y Λ≤log N se obtiene solo un tamaño O((N/c)log N) para F. Tras sumar c, no produce HμΛ. Hace falta información aritmética sobre la correlación, no volver a aplicar una cota de una sola función.

**Resultado del intento:** se han demostrado las reducciones anteriores, pero no HμΛ, Hpuntual ni la desigualdad unilateral. No hay aquí una demostración incondicional de Goldbach.

## 7. Comprobaciones, aceptación y continuidad

**[N] Verificador:** `python numerico/puente_paridad_check.py --hasta 600`. Sin numpy, instalaciones, servicios, modelos ni cálculos grandes. Las 1193 comprobaciones de (4.1) comparan polinomios formales en log p con coeficientes enteros: el resultado no depende de tolerancias flotantes. Las 3000 comprobaciones de (4.5) usan Fraction. Se comprueba por separado la restricción de bases de potencias propias en la cola no coprima. Los diagnósticos de U y T usan logaritmos flotantes y son solamente muestras, sin validar un o(·).

**Revisión adversarial de esta nota:**

| Comprobación | Resultado |
|---|---|
| Signo de la cola, límites d>D y exclusión de Λ(0) | Resueltos algebraicamente y contrastados en (4.1). |
| Primos que dividen N y clases no reducidas | Separados en (4.2); no se invoca AP coprima sobre ellos. |
| Constantes dependientes de N en producto de Euler | Se explicita y prueba el factor N^ε en (4.7). |
| «Casi todos» frente a «todos» | HμΛ exige cada N par; ningún promedio en N la sustituye. |
| Eliminación de potencias propias | Cota elemental (4.12), sin la constante pendiente 2,9. |
| Identificación de condición con teorema | HμΛ permanece [H]; no se afirma una mejora de literatura. |
| Formalización / revisión independiente | Pendientes. Lean y lake no están instalados; no se simula un PASS. |
| Tramo finito para Goldbach completa | Falta un N₀ efectivo y una verificación certificada que cubra hasta ese N₀. |

**Siguiente problema real:** demostrar una cota agregada para las correlaciones primitivas F, o una cota unilateral para su suma logarítmica, usando un ingrediente que no sea AP/M por sí solos. Corolario 6.1 es el punto de comparación preciso en el cuerpo de funciones. Si se busca en la prueba automorfa de Alpöge, habrá que obtener un lema con esos cuantificadores, demostrarlo y revisar sus dependencias; citar RH/GRH no lo reemplaza.
