# Documento 107 — Prueba escrita de L105(a) y descomposición por niveles de ℳ_Ω

**Veredicto en una frase.** **L105(a) queda DEMOSTRADO**: la cadena
cota puntual (106.6) + conteo trivial del marco cierra
B_Δ ≤ C_Δ·T(N) con C_Δ = 10⁻³ explícita para todo N ≥ 2,36·10¹⁴ y
B_Δ/T(N) → 0, sin ningún ingrediente de distribución; de paso se
registra una **errata al doc 97**: su cota impresa |Δ_r| ≤ ω(r)·log N
(sin φ(r)) es falsa en general (contraejemplo medido y recomputado
por dos vías: r = 458815, factor 1846). Y ℳ_Ω por niveles muestra que
el hueco de L105(b) **no se concentra en un nivel**: la masa de peso
vive en los niveles conmutados S2–S5 (≥ 97%) y en módulos grandes
(d > D_A^{2/3}: 80–97%), pero la cota clásica por nivel sobrepasa
T(N) en **todos** los niveles (exceso mínimo por nivel ≥ 3,2·10³;
total 1,1·10⁷ en z = 11 creciendo a 1,8·10⁹ en z = 29): el enunciado
mínimo que cerraría L105(b) es un BDH de clase única sobre el
conjunto delgado 𝒟 con ganancia polinómica N^{2c}, que no es
consecuencia de BV, de BDH ni de los teoremas de dispersión de BFI
tal como están publicados — es un lema genuinamente nuevo, con la
descomposición por niveles (106.2) como su insumo estructural
plausible.

Cadena: doc 100 (Δ_r, T(N)) → doc 101 (conteo A100.1) → doc 105
(Lema L105: B_M ≤ B̂₅ [T]; faltan (a) B_Δ ≤ C_Δ·T(N) y
(b) D̂^C ≤ C_Ω·T(N)) → doc 106 (L105(b) reducido a ℳ_Ω (106.5);
B_Δ medido ≤ 2,16·10⁻⁴·T con prueba elemental esbozada) → **este
documento escribe la prueba de (a) completa y descompone ℳ_Ω por
niveles para dictaminar (b)**.

Script: `107_prueba_l105a_niveles_momega.py` (solo stdlib; importa
sin cambios la maquinaria verificada de los docs 102–106). Salida:
`salida_doc107_local.txt`. Tiempo total: 9,9 s. Chequeos K0–K6:
**PASS** (detalle en el Apéndice A).

**Disciplina.** [T] = teorema demostrado aquí o en docs anteriores
citados; [M] = medición reproducible del script (no es prueba);
[P-rev] = afirmación de literatura no reprobada aquí. La prueba de
§1 no usa ninguna medición: las mediciones son calibración
independiente de la cadena. No se usa GRH. Goldbach sigue abierta.
Este documento no modifica los docs 85–106.

**Notación (docs 105/106).** Marco: A = 3 fijo; z = N^{1/10},
y = N^{1/3}, D_A = N^{1/2}; H el soporte (cierre por divisores de
m = s·d); ℋ = {r ∈ H : (r, N) = 1, r > R0}, R0 = (ln N)^κ, κ = 1;
𝒟 = {d ≥ 2 : d | r para algún r ∈ ℋ}; S_r = M_r + Δ_r con Δ_r la
corrección imprimitiva (docs 97–100); B_Δ = B_Δ^P + B_Δ^C =
Σ_{r∈ℋ} Δ_r²/φ(r); Ω_d^C y U_d^C los pesos de conteo (106.2);
ℳ_Ω(N) = Σ_{d∈𝒟} Ω_d^C φ(d)²E(N; d, N)² (106.5);
T(N) = N^{2−2c}/(0,0225·ln⁴N) con c = 0,039884, 2−2c = 1,920232
(doc 100 §5.3). Nivel de d ∈ 𝒟: el par (tipo, k) con tipo **L**
(liso: todo primo de d es < z) o **S** (conmutado: exactamente un
primo ≥ z, el primo de conmutación s₀ — d no puede tener dos, pues
r ∈ ℋ tiene como mucho uno) y k = ω(d) = |supp(d)|, de modo que
τ(d) = 2^k; es la estructura de factorización de (106.2):
d = d₀ (parte lisa) × (s₀ opcional).

---

## 1. Tarea 1 — Prueba completa de L105(a): B_Δ ≤ C_Δ·T(N)

### 1.1. Paso 1: la forma cerrada de Δ_r [T, docs 97/100]

Para r > 1 con (r, N) = 1, el Teorema 97.1 (doc 97 §3.3; verificado
por dos rutas independientes — ortogonalidad y cálculo directo por
caracteres — en el doc 100 §2, acuerdo a 1,43·10⁻⁸) da

Δ_r = Σ_{p | r, p ∤ r/p} Σ*_{χ* mod r/p} χ̄*(N)·g_p(χ*),
g_p(χ*) := Σ_{j ≥ 1, p^j ≤ N} χ*(p)^j · log p,

donde Σ* suma sobre los caracteres primitivos mod r/p (con la
convención de que si r/p = 1 hay un único carácter, el principal
mod 1). La condición es p ∥ r (p divide a r exactamente una vez).
Equivalentemente Δ_r = Σ_{p∥r} P_{r/p} con P_e la forma de Möbius
del doc 100 §2.1.

### 1.2. Paso 2: la cota puntual (106.6) [T]

**Lema (cota puntual).** |Δ_r| ≤ φ(r)·ω(r)·ln N para todo r ∈ ℋ.

*Demostración.* Cuatro hechos elementales:

(i) |g_p(χ*)| ≤ Σ_{j: p^j≤N} log p = ⌊log_p N⌋·log p ≤ ln N,
    pues |χ*(p)| ≤ 1.
(ii) El número de caracteres primitivos mod e es ≤ φ(e): son un
    subconjunto de los φ(e) caracteres mod e.
(iii) Si p ∥ r, entonces φ(r/p) ≤ φ(r): escribiendo r = p·m con
    p ∤ m, φ(r) = (p−1)·φ(m) ≥ φ(m) = φ(r/p).
(iv) El número de primos p con p ∥ r es ≤ ω(r).

Aplicando la desigualdad triangular a la forma cerrada de §1.1 y
luego (i), (ii), (iii), (iv) en ese orden:

|Δ_r| ≤ Σ_{p∥r} φ(r/p)·ln N ≤ ω(r)·φ(r)·ln N. ∎

**Errata registrada (doc 97).** El doc 97 §3.3 imprime la cota
«|Δ_r| ≤ ω(r)·log N», **sin** el factor φ(r). Esa forma no se sigue
de la demostración: exigiría que la suma ortogonal
Σ*_{χ* mod e} χ̄*(N)χ*(p)^j fuese ≤ 1 en valor absoluto para todo
e, y no lo es — P_e es en general del orden de una fracción de
e·ln N, no de ln N. La forma impresa se cumple en los conductores
pequeños del doc 97 (donde la suma colapsa), pero como cota general
es **falsa**:

- *Contraejemplo [M].* r = 458815 ∈ ℋ(S23) (ω(r) = 3,
  φ(r) = 314592), N = 999982: |Δ_r| = 76522,07 frente a
  ω(r)·ln N = 41,4 — la cota impresa falla por un factor
  **1846,3**. El valor de Δ_r fue recomputado por dos
  implementaciones independientes (la `delta_r_orth` de la serie y
  una implementación desde cero de la fórmula P_e del doc 100 con
  factorización por división de prueba): acuerdo a 6 decimales
  (76522,068094 en ambas). Su contribución a B_Δ es
  Δ_r²/φ(r) ≈ 1,86·10⁴, el 21% del B_Δ total de S23: no es un
  artefacto numérico.

La prueba de L105(a) usa (106.6), demostrada arriba; no usa la
forma impresa del doc 97 en ningún paso. Este es el único punto de
la cadena donde un enunciado anterior era incorrecto, y queda
corregido con el contraejemplo.

*Calibración [M] (K1).* En R1, S11, S13, S17, S23 (ambos N):
máx_{r∈ℋ} |Δ_r|/(φ(r)ω(r) ln N) = **0,0513**, alcanzado en el primo
r = 17 (Δ_17 = (log 17)⌊log_17 N⌋ = 11,3329): la cota (106.6) se
cumple en ventana con una holgura de ×19,5. El ratio sin φ alcanza
1846,3 (el contraejemplo de arriba).

### 1.3. Paso 3: de la cota puntual a B_Δ [T]

B_Δ = Σ_{r∈ℋ} Δ_r²/φ(r) ≤ (ln N)² · Σ_{r∈ℋ} φ(r)·ω(r)²
=: (ln N)² · M_ω,

sustituyendo (106.6) término a término: Δ_r²/φ(r) ≤
φ(r)²ω(r)²ln²N/φ(r) = φ(r)ω(r)²ln²N. Álgebra pura. [M: la
desigualdad se verifica en los 12 puntos con holgura ≥ 2·10⁴ (K1).]

### 1.4. Paso 4: el conteo del marco [T] y los momentos auxiliares

Todo r ∈ ℋ divide a algún m = s·d < D_A por construcción del
soporte (d < D_A/s en cada rama). Luego:

(i) r < D_A = N^{1/2}, y φ(r) ≤ r < N^{1/2};
(ii) ω(r) ≤ ln r/ln 2 < (ln N)/(2 ln 2) (cada primo ≥ 2; en el
    soporte de hecho ≥ 5, pues (r, 6) = 1 — no se usa);
(iii) |ℋ| ≤ D_A = N^{1/2} (ℋ ⊂ [1, D_A), enteros distintos).

Por tanto M_ω ≤ N^{1/2}·(ln N)²/(4(ln 2)²)·N^{1/2} y

**B_Δ ≤ N·(ln N)⁴/(4(ln 2)²) = 0,520342·N·(ln N)⁴.**   (107.1)

*Sobre la crudeza de (iii).* Es deliberada: la prueba no necesita
más. La fórmula exacta de A100.1 (doc 101) da |H| ≈ 0,033·N^{0,4208}
en la familia escalada [M], lo que bajaría la escala de B_Δ a
polilog·N^{0,9208} (la cifra citada en el doc 106 §4.2); pero
convertir esa fórmula en una cota analítica |H| ≪ N^{0,43} exige
acotar |V⁺(L)| mediante las desigualdades de Nathanson, que es la
dificultad registrada en el doc 101 (α_W = 0,4289 es medición con
intervalo, no cota probada). La cadena [T] de este doc usa solo el
conteo trivial (iii) y aun así cierra con margen N^{0,92} (§1.5).

*Momentos auxiliares estilo A100.1 [T reindexación + M].* Los
momentos que la cadena consume, M_ω = Σ_{r∈ℋ} φ(r)ω(r)² y
M_φ = Σ_{r∈ℋ} φ(r), se calculan sin enumerar pares por la
descomposición en ramas del doc 101 (ℋ = 𝒱₀ ⊔ ⨆_s s·𝒱_s,
disjunta): en la rama s, φ(s·h) = (s−1)φ(h) y ω(s·h) = ω(h)+1
(s ∤ h siempre: los primos de h son < z ≤ s). La fórmula por ramas
reproduce la enumeración directa con **error relativo 0,00** en los
11 puntos comprobados (K2). Valores (familia gen): M_ω = 7,14·10⁶
(S11), 3,95·10⁷ (S13), 5,64·10⁸ (S17), 1,73·10⁹ (S19), 1,07·10¹⁰
(S23); M_φ = 1,64·10⁶, 7,43·10⁶, 8,95·10⁷, 2,52·10⁸, 1,37·10⁹.

### 1.5. Paso 5: comparación con T(N) y constante explícita [T]

De (107.1) y la definición de T(N):

B_Δ/T(N) ≤ [N(ln N)⁴/(4(ln 2)²)] · [0,0225·ln⁴N/N^{2−2c}]
= (0,0225/(4(ln 2)²))·(ln N)⁸·N^{−(1−2c)}
= **0,0117077·(ln N)⁸·N^{−0,920232} → 0**.   (107.2)

La identidad es algebraica (verificada numéricamente, K6) y el
cociente es decreciente para N > 5970. Constantes explícitas por
bisección sobre (107.2) [cómputo, K6]:

| C_Δ | N₀: la cota vale para todo N ≥ N₀ |
|---|---|
| 1 | 3,82·10⁹ |
| 10⁻² | 7,40·10¹² |
| 10⁻³ | 2,36·10¹⁴ |

Y en la ventana medida (12 puntos, docs 105/106 reproducidos, K3):
B_Δ/T(N) ≤ **2,16·10⁻⁴** (máximo en gen/S23), así que C_Δ = 10⁻³ ya
cubre la ventana con factor ≥ 4,6 — coherente con N₀(10⁻³) del
marco: los soportes S29/S31/S37 del marco tienen
N_f = z¹⁰ ≥ 4,2·10¹⁴ > N₀.

*Nota de ventana (la del doc 106 §4.2, con mis números).* Evaluada
en el N del modelo (≈ 10⁶) con soportes escalados al N del marco,
la cota cruda (ln N)²·M_ω sobrepasa T(N): crude/T va de 3,36
(gen/S11) a 5,03·10³ (gen/S23). Es el artefacto ya declarado: T se
evalúa en el N equivocado. En el marco (N = N_f = z¹⁰) la cadena
(107.2) da, p. ej. para z = 23, ratio ≈ 3,2·10⁻³. Ninguna de las
dos lecturas es una obstrucción: la del marco es la del enunciado.

### 1.6. Teorema y cadena de dependencias

**Teorema (L105(a)).** En el marco de la serie (A = 3,
z = N^{1/10}, y = N^{1/3}, D_A = N^{1/2}, κ = 1):

B_Δ ≤ C_Δ·T(N) con C_Δ = 10⁻³ para todo N ≥ 2,36·10¹⁴;

y B_Δ/T(N) → 0, luego todo C_Δ > 0 vale desde algún N₀(C_Δ).
**L105(a) queda DEMOSTRADO.**

Cadena de dependencias, paso a paso:

1. Forma cerrada de Δ_r (Thm 97.1) — [T], docs 97/100, verificada
   por dos rutas (ortogonalidad y caracteres directos).
2. Conteo de caracteres primitivos ≤ φ(e) y φ(r/p) ≤ φ(r) para
   p ∥ r — elemental [T], §1.2.
3. Cota puntual (106.6) — [T], §1.2 (con la errata del doc 97
   corregida: la forma sin φ(r) es falsa, contraejemplo en §1.2).
4. B_Δ ≤ (ln N)²·M_ω — [T], álgebra, §1.3.
5. Conteo del marco: r | m < D_A ⇒ φ(r) < N^{1/2},
   ω(r) < (ln N)/(2 ln 2), |ℋ| ≤ N^{1/2} — [T] por construcción
   del soporte, §1.4.
6. Definición de T(N) (doc 100 §5.3) y comparación (107.2) — [T].

Ningún paso usa BV, BDH, criba grande ni GRH; ningún paso es una
medición. Las mediciones del script (K1–K3, K6) son calibración
independiente de cada eslabón: la cota puntual se cumple con
holgura ×19,5, los momentos por ramas coinciden exactamente con
la enumeración, y B_Δ/T ≤ 2,16·10⁻⁴ en ventana.

**Estado del Lema L105 tras §1:** la pata (a) está cerrada. Queda
exactamente la pata (b), reducida por el doc 106 a ℳ_Ω (106.5).
La §2 dictamina qué falta allí.

---

## 2. Tarea 2 — ℳ_Ω por niveles

### 2.1. (a) La masa de peso por niveles [M, conteo exacto]

Masa de peso del nivel ℓ: Σ_{d ∈ 𝒟_ℓ} Ω_d^C·φ(d)², normalizada
por el total Σ_{d∈𝒟} Ω_d^C φ(d)² (computable sin ψ por (106.2);
verificado contra enumeración, K4: las fracciones suman 1 con
error < 10⁻¹²; |𝒟_C| = 577 en S11 y 47102 en S29).

| z | L (total) | S1 | S2 | S3 | S4 | S5 | d > D_A^{2/3} |
|---|---|---|---|---|---|---|---|
| 11 | 0,0000 | 0,019 | 0,381 | 0,600 | — | — | 0,798 |
| 13 | 0,0001 | 0,016 | 0,329 | 0,636 | 0,019 | — | 0,857 |
| 17 | 0,0003 | 0,014 | 0,282 | 0,652 | 0,042 | 0,011 | 0,916 |
| 19 | 0,0003 | 0,011 | 0,239 | 0,644 | 0,070 | 0,037 | 0,931 |
| 23 | 0,0001 | 0,007 | 0,165 | 0,660 | 0,100 | 0,068 | 0,954 |
| 29 | 0,0001 | 0,004 | 0,118 | 0,643 | 0,124 | 0,111 | 0,973 |

Lectura:

1. **Los niveles lisos no existen a efectos de masa**: ≤ 0,03% en
   todo z. La masa es íntegramente conmutada (d contiene el primo
   de conmutación s₀): 𝒟 hereda la estructura del soporte, cuyos
   elementos genéricos son s·h.
2. **El nivel S3 manda en conteo** (0,60–0,66, estable), seguido
   de S2, que cede masa a S4 y S5 al crecer z (S2: 0,381 → 0,118;
   S4+S5: 0 → 0,235). S1 es marginal y decreciente (0,019 → 0,004).
3. **La masa vive en módulos grandes**: la fracción con
   d > D_A^{2/3} crece de 0,798 a 0,973; la fracción con
   d ≤ D_A^{1/3} es 0,000 en todo z. ℳ_Ω no es un momento de
   módulos pequeños disfrazado.
4. *Calibración con ψ [M]*: la masa real medida de ℳ_Ω (= D̂^C,
   N = 999982) por nivel sigue el mismo patrón: S2+S3 = 0,906,
   0,908, 0,877, 0,846, 0,809 para z = 11…23 (S3 el mayor desde
   z = 17: 0,448–0,481; S4+S5 creciendo hasta 0,145; niveles L
   ≤ 0,002). La descomposición de conteo predice bien dónde está
   el momento real: lo que haya que acotar está en S2–S5.

### 2.2. (b) Cotas clásicas por nivel contra T(N) [M sobre fórmulas clásicas]

A cada nivel ℓ se le aplica, con sus pesos medidos
W*_ℓ = máx_{d∈ℓ} Ω_d^C φ(d) y su techo de módulo Q_ℓ = máx 𝒟_ℓ,
la mejor de las dos cotas clásicas disponibles, en el marco
N_f = z¹⁰ (D_A = z⁵):

- **Gran criba (vía Parseval del doc 106):**
  ℳ_Ω,ℓ ≤ W*_ℓ·(N_f + Q_ℓ²)·Λ₂(N_f), con Λ₂(N) ≈ N ln N.
- **BDH pagando el peso:** como BDH controla la varianza sobre
  todas las clases y ℳ_Ω carga φ(d)² en una sola clase, la ruta
  BDH paga φ(d) ≤ Q_ℓ: ℳ_Ω,ℓ ≤ W*_ℓ·Q_ℓ·(Q_ℓ N_f ln N_f +
  N_f²/(ln N_f)^A), con A = 2 representativo (la conclusión no
  depende de A: la brecha es polinómica, §2.3).

Resultados (exceso = cota/T(N_f)):

| z | exceso total (Σ mín) | exceso solo-LS | nivel de mayor exceso LS | exceso mínimo por nivel |
|---|---|---|---|---|
| 11 | 1,14·10⁷ | 3,39·10⁷ | S3: 8,5·10⁶ | L1 (BDH): 3,2·10³ |
| 13 | 5,79·10⁷ | 1,25·10⁸ | S4: 2,9·10⁷ | L1 (BDH): 8,1·10³ |
| 17 | 3,28·10⁸ | 5,47·10⁸ | S5: 1,2·10⁸ | L1 (BDH): 1,8·10⁴ |
| 19 | 4,65·10⁸ | 7,76·10⁸ | S5: 1,6·10⁸ | L1 (BDH): 3,2·10⁴ |
| 23 | 8,04·10⁸ | 1,36·10⁹ | S5: 2,6·10⁸ | L1 (BDH): 5,5·10⁴ |
| 29 | 1,84·10⁹ | 2,72·10⁹ | S5: 4,5·10⁸ | L1 (BDH): 1,1·10⁵ |

Y la identidad de consistencia con el doc 106: el exceso solo-LS
total es, a precisión de redondeo,

exceso_LS ≈ 0,0225·(Σ_ℓ W*_ℓ)·N_f^{2c}·(ln N_f)⁵,

con Σ_ℓ W*_ℓ = 28,0 (z = 11) creciendo a 189,6 (z = 29) — porque
Q_ℓ² ≪ N_f en casi todo nivel y cada nivel **repaga** la escala
completa (N_f + Q_ℓ²)·Λ₂ ≈ N_f² ln N_f. Es la cuenta global del
doc 106 (0,045·W*·N^{2c}·ln⁵N) descompuesta: mismo orden, ahora
con Σ_ℓ W*_ℓ en el papel de 2W*.

**Respuesta a (b).** El exceso N^{2c} **no se concentra en un
nivel: se reparte por todos, porque es estructural**. Cada nivel
— incluido el nivel liso L1, de masa ~10⁻⁵, tratado con BDH —
sobrepasa T(N_f) por ≥ 3,2·10³ individualmente: la gran criba no
ve ni la delgadez de 𝒟 (|𝒟| ~ N^{0,42}) ni la estructura de
niveles, pues su término (N + Q²) ya vale ≈ N con Q = 0. Si se
pregunta qué nivel manda *en la suma de cotas*, es el nivel
conmutado superior S_{k_max} (S5 desde z = 17; W*_ℓ = W* global
= 31 y Q_ℓ el más cercano a D_A), con el 16% del total solo-LS en
z = 29, seguido de S4 y S3 (~43% entre los tres); pero eliminar
ese nivel dejaría el exceso en ~10⁹ igualmente. En la masa (§2.1)
manda S3; en la cota manda S_{k_max}; en **ningún** nivel la
tecnología clásica cierra. La conclusión del doc 106 se afina:
el hueco no es «de distribución en general», es de distribución
*uniformemente en todos los niveles del soporte*.

### 2.3. (c) Veredicto sobre L105(b): el enunciado mínimo y su lugar en la literatura

**Enunciado mínimo que cerraría ℳ_Ω (y con él L105(b)).**

**(DIST-𝒟).** Existe una constante absoluta C tal que para todo N
del marco,
Σ_{d∈𝒟} Ω_d^C·φ(d)²·E(N; d, N)² ≤ C·T(N),
con 𝒟 y Ω_d^C los objetos de conteo explícitos de (106.2).

Equivalentemente: una cota de varianza de **clase única**
(a ≡ N mod d) sobre el **conjunto delgado** 𝒟 (|𝒟| ~ N^{0,42},
divisores de un soporte de criba con factorización por niveles
explícita), con pesos Ω_d^C φ(d) ≤ W* ≤ 32,4 medido [M]
(polilogarítmicos), que gana el factor **polinómico** N^{2c}
(2c = 0,079768) sobre la escala N²·polilog que dan gran criba y
BDH, a nivel D_A = N^{1/2}. Por §2.1 el contenido efectivo está en
los niveles S2–S5 y en módulos d > D_A^{2/3}; por §2.2 la ganancia
se necesita en todos los niveles a la vez, así que (DIST-𝒟) no se
debilita restringiéndolo a un nivel: esa es su forma mínima útil.

**Comparación con los tres resultados clásicos [P-rev en la forma
estándar; las diferencias son de enunciado, verificables].**

- **Bombieri–Vinogradov.** Promedia máx_a |E(x; q, a)| sobre
  q ≤ x^{1/2}/(log x)^B: es un resultado de primer momento con
  máximo sobre clases. No controla el segundo momento ponderado
  Σ Ω_d φ(d)²E², y su escala natural en segundo momento es la de
  la gran criba (x²·polilog): no contiene ninguna ganancia
  polinómica x^{2c} en esa escala. Además promedia sobre un
  intervalo de módulos, no sobre un conjunto delgado con pesos.
- **Barban–Davenport–Halberstam.** Σ_{q≤Q} Σ_{(a,q)=1} |E(x;q,a)|²
  ≪ Qx log x + x²/(log x)^A para Q ≤ x/(log x)^B. Difiere de
  (DIST-𝒟) en tres puntos exactos: (1) suma sobre **todas** las
  clases a y **todos** los q ≤ Q — restringir a una clase y a un
  conjunto delgado no mejora su demostración, porque el término
  x²/(log x)^A nace del término diagonal de la gran criba, que es
  óptimo en general; (2) su ahorro en el término de escala x² es
  solo (log x)^A, y (DIST-𝒟) exige x^{2c} polinómico — para
  cualquier A fijo, x²/(log x)^A / T(x) → ∞; (3) no admite pesos
  φ(q)²Ω_q sin pagarlos como en §2.2, lo que empeora la cota.
  La cuenta por niveles de §2.2 es la verificación numérica de
  (1)–(3): ni siquiera nivel a nivel, con W*_ℓ y Q_ℓ propios,
  BDH baja de la escala x²·polilog.
- **Bombieri–Friedlander–Iwaniec (dispersión).** Es el pariente
  cercano en espíritu: el método de dispersión sí obtiene
  ganancias polinómicas más allá de la barrera de la gran criba
  explotando estructura bilineal (p. ej. teoremas tipo BV con
  nivel 4/7 para pesos bien factorizables). Pero (DIST-𝒟) no es
  un caso de sus teoremas: (1) es un **segundo momento** con pesos
  φ(d)²Ω_d^C, no un promedio de primer momento tipo BV;
  (2) el conjunto de módulos es el conjunto de divisores de un
  soporte de criba — delgado y con factorización por niveles
  dada por (106.2) — no un peso bien factorizable sobre un
  intervalo en el sentido de BFI, aunque la factorización por
  niveles es exactamente el *tipo* de estructura que la
  dispersión explota; (3) la clase es la diagonal a ≡ N (mod d),
  con N variando con el problema, no una clase fija ni un promedio
  sobre clases con cancelación del término principal.

**Veredicto.** (DIST-𝒟) **no se deduce de BV, de BDH ni de los
teoremas publicados de BFI**: es un enunciado de distribución
genuinamente nuevo — un «BDH para el conjunto 𝒟 con ganancia
polinómica N^{2c}» — cuya única ruta plausible visible es una
dispersión construida sobre la descomposición por niveles
(106.2) de Ω_d^C y 𝒟, trabajo que no existe en la literatura
citada. L105(b) queda, tras este doc, reducido a ese único
enunciado, con su contenido localizado (niveles S2–S5, módulos
grandes) y su distancia a lo clásico medida nivel a nivel.

---

## 3. Qué NO se ha demostrado

- **L105(b) sigue abierto**: (DIST-𝒟) no está probado; §2 mide su
  distancia a lo clásico, no lo cierra. El Lema L105 está por
  tanto demostrado en su pata (a) y reducido a un enunciado
  explícito en su pata (b).
- La prueba de L105(a) es en el marco asintótico (familia de
  soportes escalada, N → ∞); la constante C_Δ = 10⁻³ está probada
  para N ≥ 2,36·10¹⁴ y medida en ventana (≤ 2,16·10⁻⁴), no
  probada para N pequeño — donde tampoco la necesita el balance
  de Chen, que es asintótico.
- Las leyes y fracciones de §2 son mediciones [M] en
  z ∈ [11, 29] (conteo) y z ∈ [11, 23] (con ψ), no asintóticas
  probadas; Λ₂(N) ≈ N ln N y la forma de BDH con A = 2 son las
  formas clásicas estándar usadas como escala [P-rev], y las
  conclusiones de §2 son de escala (órdenes con N^{2c}), no de
  constantes.
- Aun con L105 completo faltaría el puente de los docs 94/96
  (problema independiente) y el resto de la cadena de Chen en el
  marco Q-twisted. **Goldbach sigue abierta.**

## Apéndice A. Chequeos del script (todos PASS)

| Chequeo | Contenido | Resultado |
|---|---|---|
| K0 | Maquinaria: R2 \|H\| = 880; S11 == Model97; build_H106 == build_H en z = 11, 13, 17 | PASS |
| K1 | Cota puntual (106.6) en R1/S11/S13: máx \|Δ_r\|/(φ(r)ω(r) ln N) = 0,0513 ≤ 1; B_Δ ≤ (ln N)²·M_ω en los 12 puntos | PASS |
| K2 | M_ω por ramas (A100.1) == enumeración directa, error relativo 0,00 en 11 puntos | PASS |
| K3 | máx B_Δ/T(N) en 12 puntos = 2,159·10⁻⁴ (doc 106: 2,16·10⁻⁴); C_Δ = 10⁻³ cubre la ventana con factor ≥ 4,6 | PASS |
| K4 | Masas por nivel suman el total (rel < 10⁻¹²) en z = 11…29; \|𝒟_C\| = 577 (S11) … 47102 (S29) | PASS |
| K5 | W* por niveles == doc 106 (9,52/15,00/31,00 en z = 11/13/17); Q_ℓ ≤ D_A; la cota clásica total NO alcanza T (exceso > 1 en todo z) | PASS |
| K6 | Identidad (107.2) verificada; N₀ por bisección: 3,82·10⁹ (C_Δ=1), 7,40·10¹² (10⁻²), 2,36·10¹⁴ (10⁻³) | PASS |

Reproducibilidad: `python3 107_prueba_l105a_niveles_momega.py`
(solo stdlib, ~9,9 s), salida en `salida_doc107_local.txt`. El
script reutiliza sin cambios la maquinaria de los docs 102–106
(`Model97`, `build_H106`, `SegPsi`, `delta_r_orth`, `Counting`,
`measure105`, `T_of_N`) y añade: la medida por conductor de Δ_r
con los ratios de la cota puntual, los momentos M_ω/M_φ, el conteo
M_ω por ramas, la cadena (107.1)–(107.2) con bisección de N₀, y
la descomposición por niveles de la masa de ℳ_Ω con las cotas
clásicas por nivel. El contraejemplo de §1.2 (r = 458815) fue
además recomputado con una implementación independiente de la
fórmula P_e del doc 100 (factorización por división de prueba):
Δ_r = 76522,068094 por ambas vías.

## Apéndice B. Qué no es este doc

No es una demostración de la cuestión (b) ni de Goldbach: es el
cierre probado de la pata elemental del Lema L105 y la localización
exacta — por niveles, con la distancia a BV/BDH/BFI enunciada sin
inventar referencias — del único enunciado de distribución que
falta. La desigualdad B_Δ ≤ C_Δ·T(N) sí es un teorema [T] dentro
del marco, con constante y N₀ explícitos.
