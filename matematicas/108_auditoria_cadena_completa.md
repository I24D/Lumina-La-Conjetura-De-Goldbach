# Documento 108 — Auditoría de la cadena completa: de (DIST-𝒟) al balance final

**Veredicto en una frase.** Asumiendo (DIST-𝒟) **tal como está enunciado**
en el doc 107 (pesos Ω^C), la cadena **NO cierra** hasta el balance
final sin más huecos bloqueantes: quedan exactamente **dos** huecos
analíticos además de (DIST-𝒟) — **G1, la pieza prima de B_M**, que el
Lema L105 no acota (D̂^P es una identidad, no una cota; acotarla a la
escala T(N) exige la misma ganancia N^{2c} sobre un conjunto delgado,
es decir (DIST-𝒟) extendido a los pesos Ω^P) — y **G2, el primer factor
del puente de Cauchy**, la cota de Gram sup Q_s ≪ 1/log z_A (en su
forma normalizada), que sigue **medida y no probada**: sin ella el
puente entrega R_grande = O(N/(A ln²N)) con constante
√(G_grande·C_B/0,0225), no el o(·) que (89.8) exige, y la constante
0,0225 dentro de T(N) es ella misma una medición. Si (DIST-𝒟) se
extiende a (DIST-𝒟⁺) (ambos sistemas de pesos) **y** se concede G2,
entonces **sí**: la cadena cierra eslabón a eslabón hasta el balance,
y el teorema condicional resultante es **Chen en el marco Q-twisted
con uniformidad en el selector** — positividad de la **suma**
C_k + C_{k+1} vía I_k > 0 — que **no** separa las capas (el problema
central del doc 81 queda intacto, fuera de la cadena) y **no**
demuestra Goldbach. Todo lo demás abierto es deuda de escritura o
verificación [P-rev] (docs 81/83/94/96 no están en local), listada en
§5 con su tipo.

Cadena auditada: doc 85–88 (SW-Q) → doc 89 (uniformidad, UNIF-U) →
doc 93 (agregado, problema central) → docs 94/96 de Lumina (Gram,
puente) vía docs 95/97 → docs 98–100 (Δ_r, segundo factor, T(N)) →
docs 101/102 (BV muerta, holgura) → doc 103 (criba restringida muerta)
→ doc 104 (mayorantes agotados) → doc 105 (L105) → doc 106 (conteo Ω,
ℳ_Ω) → doc 107 (L105(a) [T]; (DIST-𝒟)) → **este documento audita las
implicaciones entre esos bloques, no rehace ninguno**.

Script: `108_auditoria_cadena_completa.py` (solo stdlib). Salida:
`salida_doc108_local.txt`. Chequeos K1–K7: **PASS** (Apéndice A): el
script verifica las identidades computables de la cadena (álgebra de
T(N) y del factor del puente, la cadena (107.2), la factorización del
precio de la pieza C, la ausencia de cancelación en la pieza prima,
el orden de umbrales, el exceso clásico, la extrapolación de |H|).

**Disciplina y etiquetas.** DEMOSTRADO = teorema [T] en la serie;
REDUCIDO = implicación exacta que traslada el problema a otro
enunciado nombrado; MEDIDO = medición reproducible [M], no prueba;
ABIERTO = sin cota ni reducción suficiente. [P-rev] = depende del
texto íntegro de un documento no disponible en local (docs 81, 83 —
viven en OneDrive — y docs 94, 96, 98, 99 de Lumina, conocidos por las
transcripciones y veredictos de los docs 95, 97, 100 y el chat). No se
usa GRH. **Goldbach sigue abierta; este documento no la demuestra ni
afirma que la cadena la demuestre.**

**Notación (docs 100–107).** Marco: A fijo (en la línea, A = 3 en
toda verificación); z = N^{1/10}, y = N^{1/3}, D_A = N^{1/2};
ℋ = {r ∈ H : (r, N) = 1, r > R0}, R0 = (ln N)^κ, κ = 1;
𝒟 = {d ≥ 2 : d | r para algún r ∈ ℋ};
S_r = M_r + Δ_r, M_r = Σ_{d|r} μ(r/d)φ(d)E(N; d, N);
B_M = Σ_{r∈ℋ} |S_r|²/φ(r), con piezas P (r primo) y C (r compuesto);
D̂^X = Σ_{d∈𝒟} Ω_d^X φ(d)²E(N; d, N)² el mayorante diagonal de la
pieza X (doc 105), con Ω_d^C dado por (106.2) y
Ω_d^P = [d primo ∈ ℋ]/φ(d);
ℳ_Ω(N) = Σ_{d∈𝒟} Ω_d^C φ(d)²E(N; d, N)² (106.5) — es decir,
ℳ_Ω ≡ D̂^C por definición;
T(N) = N^{2−2c}/(0,0225·ln⁴N), c = 0,039884 (doc 100 §5.3);
R_A = resto primitivo/bilineal del marco SW-Q con módulo exterior A
(docs 89, 94); G_A = norma de Gram (doc 94);
G_grande = f·G_A ≈ 0,05·0,45 = 0,0225, la energía de Gram del bloque
grande r > R0 con los valores medidos del doc 100 §5.3.

---

## 1. Eslabón 1: (DIST-𝒟) ⇒ ℳ_Ω ≤ C_Ω·T(N) ⇒ L105(b) ⇒ cota de B_M

### 1.1. El primer paso es definicional [REDUCIDO, exacto]

(DIST-𝒟) en la forma mínima del doc 107 §2.3 es, palabra por palabra,

Σ_{d∈𝒟} Ω_d^C·φ(d)²·E(N; d, N)² ≤ C·T(N),

y el miembro izquierdo es ℳ_Ω(N) de (106.5), que es D̂^C por
definición (doc 106 §3.2: «sustituyendo Ω_d^C por su fórmula de
conteo (106.2), D̂^C ≤ C_Ω·T(N) es la desigualdad concreta ℳ_Ω ≤
C_Ω·T(N)»). Luego **(DIST-𝒟) ⇒ D̂^C ≤ C_Ω·T(N) no es una implicación
con pérdida: es el mismo enunciado**, apoyado en que Ω_d^C y 𝒟 son
objetos de conteo explícitos con fórmulas exactas verificadas
(106.2)–(106.3) [T la reindexación; M la verificación ≤ 1,6·10⁻¹⁴].
Las reducciones intermedias del doc 106 (D̂^C ≤ ρ_C·Σ_d U_d^C F_d²;
ℳ_Ω ≤ W*·Σ_{d∈𝒟}Σ*_χ |S(χ)|² vía Parseval) son rutas alternativas
hacia (DIST-𝒟), no hipótesis adicionales de este eslabón; las leyes
«ρ_C y W* polilogarítmicos» que las acompañan son [M] (W* ≤ 32,4 y
ρ_C ≤ 31 medidos en z ≤ 29) y solo harían falta si se ataca (DIST-𝒟)
por esa ruta (hueco G8 de §5, no bloqueante para la implicación).

### 1.2. El ensamblaje de L105 [DEMOSTRADO el andamiaje; ABIERTO un sumando]

El doc 105 §1.2 demuestra [T] el andamiaje completo:

B_M ≤ Σ_{X∈{P,C}} (√D̂^X + √B_Δ^X)²,   D̂^P = B_M^{P,(M)} exactamente.

Con L105(a) [T, doc 107] — B_Δ = B_Δ^P + B_Δ^C ≤ C_Δ·T(N),
C_Δ = 10⁻³ explícita para N ≥ 2,36·10¹⁴, cadena elemental sin
BV/BDH/criba/GRH, y además B_Δ/T(N) → 0 — y con §1.1, la **pieza
compuesta cierra**:

B_M^{C,(S)} ≤ (√C_Ω + √C_Δ)²·T(N), con el √C_Δ absorbible en el o(1).

Pero la conclusión impresa del Lema L105 (doc 105 §5.4),
«B_M ≤ (√C_Ω + √C_Δ + o(1))²·T(N)», **no se sigue de (a) + (b) tal
como (b) está enunciado**: falta acotar el sumando de la pieza prima,

**B_M^{P,(M)} = D̂^P = Σ_{r primo ∈ ℋ} φ(r)·E(N; r, N)²**,

para el que la serie no contiene ninguna cota a la escala T(N). La
pieza prima está «cerrada por identidad» (doc 105 §5.5, ingrediente
4) en el sentido de que el mayorante no pierde nada en ella
(D̂^P = B_M^{P,(M)} sin factor), pero eso es exactitud del mayorante,
no una cota del objeto. Es el hueco **G1**:

- **Por qué no hay cancelación que lo abarate.** Para r primo,
  M_r = Σ_{d|r} μ(r/d)F_d = F_r exactamente (divisores 1 y r;
  F_1 = 0): la pieza prima es un **segundo momento crudo** de errores
  de progresión en módulos primos, sin estructura de Möbius
  (verificado como identidad computable, K4 del script: M_r = F_r en
  primos; M_{pq} = F_{pq} − F_p − F_q en semiprimos).
- **Qué módulos son.** Los primos de ℋ son los primos < z que dividen
  a algún d ∈ V⁺ (parte lisa) y los primos de conmutación
  s ∈ [z, y); el conteo lo dominan los de conmutación,
  ~π(y) ~ N^{1/3}/ln N. Nivel efectivo ≤ y = N^{1/3} < D_A.
- **Por qué lo clásico no basta.** Aun a nivel N^{1/3}, BDH sumando
  sobre todas las clases da la escala QN ln N + N²/(ln N)^A, y la
  gran criba de clase única la escala N²·ln N: ambas sobrepasan
  T(N) = N^{1,920232}/(0,0225 ln⁴N) por el factor
  ~N^{2c}·polilog — **la misma ganancia polinómica N^{2c} que define
  (DIST-𝒟)** (la cuenta global es la del doc 106 §3.2 / doc 107
  §2.2; K6 del script la reproduce en orden). La delgadez del
  conjunto de primos y la clase única a ≡ N son, también aquí, lo
  único que puede cerrar la brecha, y ningún teorema citado lo hace.
- **Por qué es un hueco «de la misma familia», no un muro nuevo.**
  Los pesos son más simples que los de (DIST-𝒟): Ω_d^P·φ(d) =
  [d primo ∈ ℋ] ≤ 1, frente a Ω_d^C·φ(d) ≤ W* ≈ 32 medido. El
  enunciado que falta es (DIST-𝒟) con los pesos Ω^P:

  **(DIST-𝒟_P).** Σ_{d∈𝒟} Ω_d^P·φ(d)²·E(N; d, N)² ≤ C_P·T(N).

  Nada en los docs 105–107 lo enuncia ni lo mide por separado (en
  ventana, la pieza prima es el 22–58% de B_M [M, docs 105/106:
  share_C = 0,42–0,78], así que tampoco es despreciable frente a la
  compuesta). Se registra como **bloqueante bajo la lectura literal
  de (DIST-𝒟)** y como **extensión natural (DIST-𝒟⁺) = (DIST-𝒟) +
  (DIST-𝒟_P)** en el veredicto de §6.

### 1.3. Las constantes y el factor 1/(A log²N): el álgebra encaja; la calibración es medida

El puente (§2) necesita R_grande = √(G_grande·B_M) ≪ N/(A·ln²N).
El doc 100 §5.3 construyó T(N) exactamente para eso: con
A = N^c y G_grande = 0,0225,

√(G_grande · T(N)) = √(0,0225 · N^{2−2c}/(0,0225·ln⁴N))
= N^{1−c}/ln²N = N/(A·ln²N)   **exacto**

(verificado a precisión de máquina, K1 del script, que también
reproduce la tabla de objetivos del doc 100 §5.3 a N = 10¹², 10²⁴,
10¹⁰⁰). Es decir: **la cota de B_M que sale de L105 es la que el
puente necesita, con el factor 1/(A log²N) correcto**, y con
B_M ≤ C_B·T(N) se obtiene R_grande ≤ √(G_grande·C_B/0,0225) ·
N/(A ln²N). Tres precisiones de auditoría:

1. **c dentro de T(N) es c\* = 0,039884** (el cero del balance del
   doc 89). Para c < c\* el objetivo correcto es mayor
   (N^{2−2c} crece al bajar c), así que usar T(N) con c\* es la
   elección conservadora: una cota a la escala T_{c\*}(N) implica
   la cota a la escala T_c(N) para todo c < c\*. ✓.
2. **La constante 0,0225 es [M], no [T].** Es f·G_A con f ≈ 0,05
   (fracción de energía del bloque grande según Lumina; en los
   modelos del doc 97 §4.4, f = 0,0036–0,0041 a κ = 2,7 — la
   discrepancia es la regla de atribución h frente a ah, sin fijar:
   hueco G3) y G_A ≈ 0,45 (norma medida del doc 94 a parámetros de
   modelo). Si el G_grande asintótico verdadero fuese mayor que
   0,0225, T(N) estaría calibrado demasiado pequeño y L105(b)
   pediría más de lo que el puente necesita; si es menor (lo
   esperado si G_A ~ 1/ln z, §2), T(N) es conservador. La cadena
   no tiene hoy ningún [T] que fije G_grande: su valor entra como
   medición en la definición misma del objetivo. Esto es parte del
   hueco G2.
3. **El o(·) no sale de este eslabón.** Con C_B y G_grande
   constantes absolutas, R_grande = O(N/(A ln²N)) con constante
   √(G_grande·C_B/0,0225) ≠ o(1). La forma (89.8) de UNIF-U exige
   o(·): ese factor tiene que venir del primer factor del puente
   (G_grande = o(1)), no de B_M. Véase §2.3.

**Estado del eslabón 1:** (DIST-𝒟) ⇒ pieza compuesta de B_M a la
escala T(N): **DEMOSTRADO como implicación** (definicional + L105(a)
[T] + andamiaje [T]). Cota de B_M completo: **ABIERTO por G1**
(pieza prima) y con la calibración de T(N) apoyada en [M] (parte
de G2).

---

## 2. Eslabón 2: B_M acotado ⇒ R_A — el puente de Cauchy, término a término

### 2.1. La forma del puente

R_A = Σ_{m∈M_A} b_A(m)·E(N; Am, N) (doc 95 §1). La ruta de los
docs 94/96: ortogonalidad agrupa los caracteres por conductor
primitivo r | Am; el coeficiente de los caracteres con la misma
parte nueva h = r/(r, A) es K_A(h) (doc 94, con K_A(r) = K_A(h)
[T]); Cauchy sobre el bloque grande da, en la forma que el doc 97
§6 confirma como la estructura correcta,

|R_grande| ≤ √(G_A^{grande}) · √(B_M),

con G_A^{grande} = Σ_{h: bloque grande} φ(h)|K_A(h)|² el primer
factor y B_M el segundo. **La forma exacta de la desigualdad de
Cauchy del doc 94 queda [P-rev]** (doc 97 §8.5: no verificable sin
el doc 94 íntegro; la transcripción que revisó el doc 95 no la
contenía, §4 de ese doc, y Lumina precisó después que sí existe en
el doc 94 — ambas afirmaciones son compatibles: la desigualdad
existe, su forma exacta no está auditada). La vía ingenua sin esta
estructura pierde 10⁷–10¹⁰ (doc 95 §4 [M]).

### 2.2. Inventario de términos del resto

| # | Término | Estado | Fuente |
|---|---|---|---|
| T1 | Bloque pequeño (conductores ≤ R0): Siegel–Walfisz da R_pequeño = o(N/(A log²N)) para M > 2κ+3 | **DEMOSTRADO condicional a SW** (teorema clásico, constante inefectiva), según el avance de Lumina; fijar κ y la regla de atribución (h vs ah) queda [P-rev] | doc 94 vía resumen; docs 95 §5(d2), 97 §4.1/§7.3 |
| T2 | Bloque grande: Cauchy √(G_grande·B_M) | **REDUCIDO**: factor 2 = eslabón 1; factor 1 = G2 (ABIERTO); forma exacta [P-rev] | docs 94/96 vía 95 §4, 97 §6 |
| T3 | Correcciones inducidas (carácter inducido mod Am frente a su primitivo), agregado O(D_A log²N log log N) | **Acotado según el avance de Lumina** [P-rev texto íntegro]; tamaño: D_A·polilog = N^{1/2}·polilog ≪ N/(A ln²N) para A ≤ N^{c\*} ✓ no bloqueante en tamaño | doc 94 vía resumen |
| T4 | Pieza del carácter principal: término O(\|ψ(N) − N\|·\|K_A(1)\|) | **Identificado y acotado según Lumina** [P-rev]; ψ(N) − N = o(N) por el TNP basta si el coeficiente es el declarado; \|K_A(1)\|² es el 16–34% de G_A [M], no un residuo despreciable en el factor 1 | doc 94 vía resumen; doc 95 §2.2 (H6) |
| T5 | Error imprimitivo del paso de agrupación R_A → conductores: por carácter ≤ Σ_{p\|Am, p∤r} log N; total ≪ (log N)(log log N)·Σ_m \|b_A(m)\| | **ABIERTO como escritura**: no aparece en la transcripción del doc 94 (doc 95 §2.1); su tamaño esperado es inofensivo (Σ_m\|b_A(m)\| = O(D_A)), pero **no está verificado que el Δ_r de los docs 97–100 lo absorba**: Δ_r corrige la identidad de Möbius S_r ↔ E(N; d, N) (paso posterior), no este paso previo de agrupación | doc 95 §2.1; heredado [P-rev] en doc 97 §8.6 |
| T6 | Hipótesis de soporte del paso de ortogonalidad: supp b_A ⊆ {m : (m, N) = 1} (y (A, N) = 1) | **Declarada y cierta en el marco** (es el régimen de la errata E93.1); debe escribirse en el enunciado del puente [P-rev como escritura] | doc 95 §2.1; doc 93 E93.1 |
| T7 | Δ_r dentro de S_r (identidad de Möbius corregida) | **DEMOSTRADO** [T]: identidad exacta por dos rutas (docs 97/100); su coste en el segundo factor es B_Δ, cerrado por L105(a) [T, doc 107]. La versión sin Δ_r del doc 96 fue retirada por Lumina | docs 97–100, 107 |
| T8 | Segunda reducción de K_A (límite fijo frente a variable) | Límite fijo [T]; variable falsa tal como estaba enunciada y **reparada** con la cola s′>s [T, doc 95 §2.3]; qué versión usa el doc 94: [P-rev] | doc 95 §§2.3, 7.2 |

**Respuesta a la pregunta del encargo:** el puente **no** está
completo solo con B_M y G_A como cotas sueltas: está completo como
*estructura* (T2 con T1, T3, T4 acotados o demostrados según los
avances citados), pero quedan (i) el factor 1 sin cota probada
(G2), (ii) el término T5 sin escribir ni absorber verificado, y
(iii) la forma exacta de T2 y las fijaciones de T1 en estado
[P-rev]. Ningún término adicional desconocido aparece en el
inventario: los ocho términos agotan la descomposición tal como la
serie la ha exhibido.

### 2.3. El punto fino: O(·) frente a o(·)

UNIF-U en la forma que el marco necesita (doc 89 (89.8)) es
R_B(A) ≪ (N/A)/log²N · **o(1)**. Del eslabón 1, el bloque grande
aporta √(G_grande·C_B/0,0225)·N/(A ln²N): para que el prefactor
sea o(1) hace falta **G_grande = o(1)** (C_B es una constante
absoluta de (DIST-𝒟⁺), no una cantidad que decaiga). La única vía
en la serie para G_grande = o(1) es la cota de Gram del doc 94:

sup_{z_A ≤ s < y_A} Q_s ≪ 1/log z_A,

en su **forma corregida** (doc 97 §7.2): normalizada por el
producto local, sup Q_s/Q_completo(z_A; N) ≪ 1 con
Q_completo ≍ 1/ln z_A (constante local en N: 1,7225 para
A = 3, N = 2^e [T el producto; M la cota de V⁺]). Su estado:
**MEDIDO, no probado** — Q_s·ln z ∈ [1,63, 2,83] en los modelos
(docs 95/97), con la diagonal dando ya el orden y la conservación
del orden por la restricción V⁺ sin demostrar. Este es el hueco
**G2: bloqueante** para el o(·) de (89.8) e independiente de
(DIST-𝒟) (vive en el lado de los coeficientes K_A, no en el de los
errores E). Sin G2, la cadena entrega el resto a la escala
correcta pero con constante absoluta no controlada frente a
0,0225 — y como la constante 0,0225 de T(N) es ella misma el valor
medido de G_grande (§1.3.2), el argumento sería circular si se
usara para cerrar el o(·): G2 es también lo que convierte la
calibración de T(N) en teorema.

**Estado del eslabón 2:** estructura **REDUCIDA** y completa en
inventario; cierre **ABIERTO por G2** (y por G1 vía el eslabón 1),
con deuda [P-rev] en T2/T5 y fijaciones de T1.

---

## 3. Eslabón 3: R_A controlado ⇒ UNIF-U-COND ⇒ su lugar en el marco SW-Q

### 3.1. Qué es UNIF-U-COND, exactamente

UNIF-U-COND es la conclusión (89.8) del doc 89 — el resto bilineal
del Lema SW-Q acotado por o((N/A)/log²N), es decir menor que el
término principal y no solo ≪ N/log^{A'} N — obtenida **por la ruta
de conductores primitivos** (docs 94–107) en lugar de por BV/BDH
directos. Si los eslabones 1–2 cierran, UNIF-U-COND sustituye a
UNIF-U **para el módulo exterior A del marco**, con el resto
controlado término a término (§2.2) en vez de por un teorema de
valor medio global. No es un enunciado más débil en su conclusión;
es el mismo (89.8) por otra vía.

### 3.2. Rango de u y de c: lo que cubre y lo que no

- **En c: cubre todo c fijo con c < c\* = 0,039884; no cubre
  c = c\*.** El balance del doc 89, Bal(c) = f(s_lower(c)) −
  (1/2)J(c) − (1/2)(0,4 e^γ c_{1/10})(φ(Q)/Q)/(1−c), se anula en
  c\* [T, doc 89 §3(b)]: en c = c\* el término principal del
  balance es cero y ningún resto, por pequeño que sea, produce
  positividad. Para cada c < c\* fijo, T(N) se recalibra con
  2−2c (§1.3.1) y la cadena es la misma. Los otros dos umbrales de
  la serie no son la restricción activa de esta ruta: geometría
  c ≤ 2/17 ≈ 0,1176 [T, doc 89] y umbral de Cauchy corregido
  c < 1/3 [T, doc 93 Prop. 93.5, para la ruta 2M-R] quedan ambos
  por encima de c\* (verificado, K5 del script).
- **En u (A = Qu): la línea está construida y verificada a A
  fijo.** Toda la maquinaria de los docs 94–107 (K_A, Gram, L105,
  (DIST-𝒟)) está enunciada para A fijo y verificada numéricamente
  solo a A = 3. El álgebra sigue la pista de φ(A) (la identidad de
  Gram lleva 1/φ(A)² [T]), y el doc 95 §3.3 ya advirtió que «en la
  cadena hacia R_A, el factor φ(A)² multiplica a G_A, y si A crece
  con N el objetivo para G_A debe leerse con ese factor». La
  calibración de T(N) a A = N^{c\*} (doc 100 §5.3) presupone esa
  lectura, pero **la uniformidad en A de las constantes (C_Ω, C_Δ,
  la cota de Gram, el bloque pequeño de Siegel–Walfisz con
  conductores ah donde a | A) no está demostrada ni auditada en la
  línea**: es el hueco G6 — no bloqueante para A fijo (que es lo
  que el balance de los docs 88/89 consume por selector), bloqueante
  para afirmar (89.8) uniformemente en u ≤ N^c con constantes
  independientes de u.
- **La familia de selectores Q ≤ N^θ (doc 81) es un paso
  posterior, no parte de UNIF-U-COND.** El doc 93 §6(a) lo dice:
  la uniformidad en el selector «sigue necesitando UNIF-U /
  UNIF-U-2M-R — esa es la línea viva». UNIF-U-COND da el resto por
  selector; sumar sobre la familia de selectores del doc 81 (con
  su conteo y acumulación de errores) es un eslabón distinto cuyo
  enunciado exacto depende del doc 81: **[P-rev]** (hueco G7).

**Estado del eslabón 3:** la implicación R_A controlado ⇒ (89.8)
para A fijo es **REDUCIDA/directa** (es la definición de la ruta);
su extensión uniforme en u ≤ N^c y a la familia de selectores está
**ABIERTO** en la forma de G6/G7, sin que ello afecte al balance
por selector fijo de §4.

---

## 4. Eslabón 4: del marco SW-Q + UNIF-U-COND al balance final — y el problema central

### 4.1. Qué demostraría la cadena completa

Con los eslabones 1–3 cerrados (bajo (DIST-𝒟⁺) + G2, §6), para
cada selector A = Qu del rango (§3.2) la máquina SW-Q de los
docs 86–88 corre con el resto en la forma (89.8), y el balance
(89.6) es positivo para c < c\*: se obtiene la incidencia del
doc 81 en su forma uniforme en el selector,

**I_k > 0, y en particular C_k + C_{k+1} > 0**,

con la identidad I_k = k·C_k + binom(k+1, 2)·C_{k+1} [P-rev: doc 81
no está en local; enunciado tomado de las citas de los docs 85–93].

Conviene medir el avance contra lo ya demostrado **sin**
UNIF-U-COND, para no vender dos veces el mismo teorema:

- El doc 88 ya demuestra [T] C_k + C_{k+1} > 0 para la
  sub-incidencia u = 1 (balance (88.5) positivo para todo Q;
  J_k(N; Q_N) > 0 para N grande), con la salvedad [P-rev] de la
  reescritura del Lema 32 de Bordignon con Q (Prop. 88.8; hueco
  G9, de escritura: «es escritura, no un hueco matemático», doc 88
  §6).
- Los docs 89 (Prop. 89.2) y 93 (Prop. 93.3) ya demuestran [T] la
  cota inferior sumada de I_k para Q fijo (cabeza polilogarítmica +
  no-negatividad de la cola + Rankin), sin control por-u en el
  rango potencia.

Lo que la cadena completa **añade** es, por tanto, exactamente:
**la uniformidad** — el resto (89.8) por selector en el rango
potencia u ≤ N^c, c < c\*, por la ruta de conductores — y con ella
la forma del doc 81 con selectores Q ≤ N^θ [P-rev G7]. Es el teorema
de Chen en el marco Q-twisted (tipo Lu–Cai, Chen en progresiones)
**con uniformidad en el selector**: más uniforme que el doc 88, no
más fuerte en contenido de capas que Chen.

### 4.2. Qué NO demostraría: la separación de capas queda intacta

El problema central del doc 81 (doc 93 §6(b2), [N] allí y suscrito
aquí) no es un eslabón de esta cadena: **ningún paso de los
docs 85–107 produce una segunda relación lineal independiente
entre C_k y C_{k+1}, ni una cota inferior directa de la capa C_k,
ni un argumento de signo que la fuerce**. La identidad (81) es una
ecuación en dos incógnitas no negativas; (C_k, C_{k+1}) =
(0, I_k/binom(k+1, 2)) es siempre consistente con I_k > 0, y los
pesos favorecen a la segunda capa (binom(k+1,2) > k; medido: la
capa P2 lleva el 56–61% de J_k, doc 93 H3). Es además la barrera de
paridad del switching (doc 88 Prop. 88.3 [T]: el peso W da P1+P2
conjuntamente y no puede separar). Por tanto:

- La positividad final es **de la suma** C_k + C_{k+1}, no de C_k.
- Ni siquiera con (DIST-𝒟⁺) + G2 + toda la deuda [P-rev] saldada
  la cadena toca Goldbach: Goldbach exigiría la capa P1 sola
  (C_k > 0 para el selector correspondiente), que es el problema
  central abierto, **sin ruta asignada** en la serie (doc 93 §8).

**Estado del eslabón 4:** la implicación (SW-Q + (89.8)) ⇒ balance
es **DEMOSTRADA** en su forma por-selector (docs 88/89, con la
salvedad G9); la separación de capas está **ABIERTO** y fuera de
la cadena por construcción.

---

## 5. Inventario final de huecos

Ordenados: primero los que bloquean el balance final **asumiendo
(DIST-𝒟)**; después los que no bloquean el tamaño pero impiden
declarar la cadena formalmente cerrada; al final el que no bloquea
el balance pero bloquea todo lo demás.

### 5.1. Bloqueantes asumiendo (DIST-𝒟)

| # | Hueco | Tipo | Dónde |
|---|---|---|---|
| G1 | Pieza prima B_M^{P,(M)} = Σ_{r primo ∈ ℋ} φ(r)E(N; r, N)² ≤ C_P·T(N): no enunciada ni acotada en la línea; L105(b) solo cubre D̂^C. Exige (DIST-𝒟_P): la misma ganancia N^{2c} de clase única sobre el conjunto delgado de primos (≤ y, pesos ≤ 1) | Lema nuevo (misma familia que (DIST-𝒟)) | §1.2 |
| G2 | Cota de Gram: sup Q_s/Q_completo(z_A; N) ≪ 1 (⇒ G_A ≍ 1/ln z_A, G_grande = o(1)): medida (Q_s·ln z ∈ [1,63, 2,83]), no probada; sin ella el puente da O(·) y no el o(·) de (89.8), y la constante 0,0225 de T(N) queda como medición | Demostración pendiente (criba sobre V⁺; la diagonal ya da el orden) | §§1.3, 2.3 |

### 5.2. No bloqueantes en tamaño; deuda de escritura/verificación [P-rev]

| # | Hueco | Tipo | Dónde |
|---|---|---|---|
| G3 | κ y regla de atribución de masa (h frente a ah ≤ R0) sin fijar: mueve f entre 0,004 y 0,05 y con ello la calibración de T(N) y el reparto pequeño/grande | [P-rev] + fijación de parámetros | §§1.3, 2.2 (T1) |
| G4 | Forma exacta de la desigualdad de Cauchy del doc 94 (Gram ↔ resto primitivo) no auditada | [P-rev] (doc 94 íntegro) | §2.1 |
| G5 | Error imprimitivo del paso de agrupación (T5) sin escribir en la transcripción; no verificado que Δ_r lo absorba (Δ_r pertenece al paso posterior de Möbius) | [P-rev] / escritura pendiente; tamaño esperado inofensivo | §2.2 (T5) |
| G6 | Uniformidad en A = Qu: constantes independientes de A ≤ N^{c\*} no demostradas (factor φ(A)² en el objetivo de G_A; bloque pequeño con a \| A); toda verificación es a A = 3 | Demostración pendiente (uniformidad) — bloqueante solo para la forma uniforme en u de (89.8), no para el balance por selector fijo | §3.2 |
| G7 | Suma sobre la familia de selectores Q ≤ N^θ del doc 81: enunciado y acumulación de errores | [P-rev] (doc 81) | §3.2 |
| G8 | Leyes de conteo ρ_C, W\* polilogarítmicos: [M] en z ≤ 29/37, no asintóticas probadas; solo necesarias si (DIST-𝒟) se ataca por la ruta Parseval/W\* del doc 106 | Medición sin prueba (ruta-dependiente) | §1.1 |
| G9 | Reescritura completa del Lema 32 de Bordignon con Q insertado (resto bilineal de la cola SW-Q, doc 88 Prop. 88.8): transferido por el mismo argumento [T la transferencia]; la reescritura es escritura pendiente | [P-rev] / escritura | §4.1 |
| G10 | Identificación exacta de las capas C_k/C_{k+1} del doc 81 con las capas P1/P2 de J_k (doc 88); enunciados de los docs 81/83 tomados de citas | [P-rev] (docs 81/83) | §4.1 |
| G11 | Versión (fija/variable) de la segunda reducción en el doc 94; la variable está reparada [T] en el doc 95, pero qué versión usa el doc 94 es [P-rev]; análogo: si el doc 96 ya contenía Δ_r (doc 97 §8.2, resuelto de facto por la retirada de Lumina en el doc 98) | [P-rev] | §2.2 (T8) |

### 5.3. El hueco que no bloquea el balance y bloquea todo lo demás

| # | Hueco | Tipo | Dónde |
|---|---|---|---|
| G0 | **Separación de capas**: pasar de C_k + C_{k+1} > 0 a C_k > 0. Una ecuación, dos incógnitas no negativas (doc 93 §6); barrera de paridad del switching (doc 88 Prop. 88.3). Ningún eslabón de la cadena lo toca | Problema abierto, sin ruta asignada | §4.2 |

Y el hueco de fondo, ya nombrado por el doc 107 y presupuesto por
esta auditoría: **(DIST-𝒟) mismo** — varianza de clase única sobre
𝒟 con ganancia N^{2c} a nivel N^{1/2}, no deducible de BV, BDH ni
BFI (doc 107 §2.2–2.3: el exceso clásico se reparte por **todos**
los niveles, ≥ 3,2·10³ por nivel; K6 del script). La pregunta del
encargo era si es el *único* hueco: la respuesta es que no — G1 es
de su misma familia y G2 es de otra (coeficientes, no errores) —,
pero sí es el único hueco de *distribución en la pieza compuesta*,
y tras esta auditoría el mapa completo de lo que falta es finito y
está nombrado: (DIST-𝒟⁺) + G2 para el balance; G0 para Goldbach.

---

## 6. Veredicto y teorema condicional

**¿Asumiendo (DIST-𝒟), la cadena cierra hasta el balance final sin
más huecos bloqueantes? NO**, en la lectura literal: quedan G1 y
G2 (§5.1). G1 no es un descuido de redacción del Lema L105 sino un
sumando real (la pieza prima es el 22–58% medido de B_M) cuya cota
a la escala T(N) es un enunciado de distribución del mismo tipo
que (DIST-𝒟) — por eso la forma honesta de cerrar la pregunta es
extender la hipótesis, no esconder el sumando. G2 es
independiente: vive en los coeficientes K_A y su estado es
medición consistente, no prueba.

**Sí cierra**, eslabón a eslabón y sin más huecos bloqueantes que
la deuda [P-rev] de §5.2, bajo la hipótesis extendida. El
enunciado exacto del teorema condicional que la línea ha
construido es:

> **Teorema condicional (línea UNIF-U-COND).** Supóngase:
> **(H1) (DIST-𝒟⁺)**: para el conjunto delgado 𝒟 del marco y ambos
> sistemas de pesos de conteo Ω^C (106.2) y
> Ω^P (Ω_d^P = [d primo ∈ ℋ]/φ(d)), vale
> Σ_{d∈𝒟} Ω_d·φ(d)²·E(N; d, N)² ≤ C·T(N) con C absoluta —
> la varianza de clase única a ≡ N sobre 𝒟 con ganancia
> polinómica N^{2c} sobre la escala de la gran criba/BDH, a nivel
> D_A = N^{1/2} —; y
> **(H2) la cota de Gram normalizada**: sup_{z_A≤s<y_A} Q_s ≪
> Q_completo(z_A; N) con Q_completo ≍ 1/log z_A, de donde
> G_A^{grande} = o(1).
> Entonces, para todo c fijo con 0 < c < c\* = 0,039884 y cada
> selector A = Qu con u ≤ N^c (con las constantes uniformes en A
> en el sentido de G6, y la suma sobre selectores del doc 81 en el
> sentido de G7), el resto bilineal del Lema SW-Q satisface
> **R_B(A) = o((N/A)/log²N)** — esto es UNIF-U-COND, la forma
> (89.8) de UNIF-U por la ruta de conductores primitivos —; el
> balance de Chen en el marco Q-twisted es positivo, **I_k > 0 y
> en particular C_k + C_{k+1} > 0** para N grande.
> El teorema **no** afirma C_k > 0 por separado (G0, problema
> central del doc 81) y **no** demuestra la conjetura de Goldbach.

Sin (H2), la misma hipótesis (H1) entrega el resto a la escala
exacta N/(A log²N) pero con constante absoluta
√(G_grande·C_B/0,0225) en lugar de o(1): la cadena llegaría «a un
factor constante del balance», que no es el balance. Sin la parte
Ω^P de (H1), entrega la cota solo de la pieza compuesta de B_M.
Ambas precisiones son parte del veredicto, no notas al pie.

---

## 7. Qué NO es este documento

No es una demostración de (DIST-𝒟), de G2, de UNIF-U-COND ni de
Goldbach: es la auditoría de las implicaciones entre los bloques
85–107, con cada eslabón clasificado y cada hueco nombrado con su
tipo. Las únicas afirmaciones nuevas son de auditoría: (i) la
conclusión impresa del Lema L105 no se sigue de sus patas (a)+(b)
sin acotar la pieza prima (G1, con la identidad M_r = F_r como
razón estructural); (ii) el o(·) de (89.8) no puede venir del lado
B_M y exige G_grande = o(1), es decir G2, y la constante 0,0225 de
T(N) es una medición hasta que G2 se pruebe; (iii) el inventario
de términos del puente (T1–T8) es exhaustivo sobre lo exhibido por
la serie, con T5 como único término ni escrito ni absorbido
verificado. Todo enunciado que depende de los docs 81/83/94/96/
98/99 está marcado [P-rev] y deberá re-verificarse contra los
textos íntegros cuando estén en local.

## Apéndice A. Chequeos del script (todos PASS)

| Chequeo | Contenido | Resultado |
|---|---|---|
| K1 | Álgebra del objetivo (doc 100 §5.3): T(N) reproduce la tabla (8,414·10¹⁸ / 5,803·10⁴⁰ / 1,668·10¹⁸⁴ a N = 10¹²/10²⁴/10¹⁰⁰, rel < 10⁻⁴) y √(0,0225·T(N)) = N/(A ln²N) a precisión de máquina (rel ≤ 1,4·10⁻¹⁶) | PASS |
| K2 | Cadena (107.2) de L105(a): coeficiente 0,0225/(4 ln²2) = 0,0117077; cociente decreciente para N > 5970; N₀(C_Δ = 10⁻³) = 2,361·10¹⁴ por bisección (doc 107: 2,36·10¹⁴) | PASS |
| K3 | Factorización del precio de la pieza C (doc 105 §4): Diag^C × precio τ = D̂^C/B_M^{C,(M)} en los seis puntos (rel < 2·10⁻³) | PASS |
| K4 | Pieza prima sin cancelación: M_r = F_r exacto en primos (F arbitrarios, F_1 = 0); M_{pq} = F_{pq} − F_p − F_q; Σ_{d\|r} μ(r/d) = 0 para 2 ≤ r ≤ 600 (cancelación del término principal en S_r, doc 97 paso 5) | PASS |
| K5 | Orden de umbrales: c\* = 0,039884 < 2/17 ≈ 0,1176 < 1/3; Bal(0,02) = 0,032934 > 0 > Bal(0,04) = −0,000208 (tabla del doc 89) | PASS |
| K6 | Exceso clásico (docs 106/107): 0,045·W\*·N^{2c}·ln⁵N en el marco z = 11 (N = 11¹⁰, W\* = 9,52) = 2,30·10⁷, mismo orden que el exceso medido del doc 107 (1,14·10⁷–3,39·10⁷) | PASS |
| K7 | Extrapolación del doc 102: \|H\| ≈ 0,033·N^{0,4208} a N = 10²⁴ da 4,15·10⁸ (doc 102: ~4·10⁸) | PASS |

Reproducibilidad: `python3 108_auditoria_cadena_completa.py`
(solo stdlib, < 1 s), salida en `salida_doc108_local.txt`. El
script no reutiliza la maquinaria de los docs anteriores: verifica
solo el álgebra y las identidades de esta auditoría (K4 usa una
criba de μ propia hasta 600 y valores F arbitrarios deterministas,
porque la identidad auditada es combinatoria, no numérica sobre ψ).

## Apéndice B. La cadena de un vistazo

| Eslabón | Estado | Lo que falta |
|---|---|---|
| (DIST-𝒟) ⇒ ℳ_Ω = D̂^C ≤ C_Ω·T(N) | REDUCIDO (definicional; conteo (106.2) [T]) | (DIST-𝒟) mismo: ABIERTO, lema nuevo |
| L105(a): B_Δ ≤ C_Δ·T(N) | **DEMOSTRADO** [T] (doc 107; C_Δ = 10⁻³, N ≥ 2,36·10¹⁴) | — |
| Andamiaje B_M ≤ Σ_X (√D̂^X + √B_Δ^X)² | **DEMOSTRADO** [T] (doc 105) | — |
| Pieza prima D̂^P ≤ C_P·T(N) | **ABIERTO (G1)** | (DIST-𝒟_P), misma familia que (DIST-𝒟) |
| B_M ≤ C_B·T(N) (completo) | REDUCIDO a (DIST-𝒟⁺) | G1 |
| Puente: factor 1, G_A^{grande} = o(1) | **ABIERTO (G2)** — medido, no probado | cota de Gram normalizada |
| Puente: términos T1, T3, T4, T7 | DEMOSTRADO (T7, T1 cond. SW) / acotados según Lumina (T3, T4) [P-rev] | saldar [P-rev] |
| Puente: término T5 (imprimitivo de agrupación) | ABIERTO como escritura; tamaño inofensivo | escribirlo o verificar absorción por Δ_r |
| R_A ⇒ UNIF-U-COND (89.8), A fijo | REDUCIDO (es la ruta) | G1 + G2 |
| Uniformidad en u ≤ N^c / selectores Q ≤ N^θ | ABIERTO (G6/G7) | uniformidad en A; doc 81 [P-rev] |
| SW-Q + (89.8) ⇒ I_k > 0, C_k + C_{k+1} > 0 | **DEMOSTRADO** por selector (docs 88/89; salvedad G9 de escritura) | — |
| Separación C_k > 0 (G0) | **ABIERTO**, fuera de la cadena | segunda relación independiente; sin ruta |
| Goldbach | **ABIERTO** | G0 y todo lo anterior |
