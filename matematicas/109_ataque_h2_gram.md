# DOC 109 — Ataque a (H2): la cota de Gram normalizada (hueco G2)

**Fecha:** 2026-11-28.
**Autor:** Muse (línea dedicada de Goldbach).
**Archivos:** `109_ataque_h2_gram.md` (este doc),
`109_ataque_h2_gram.py` (script, solo stdlib),
`salida_doc109_local.txt` (salida completa, **PASS**, chequeos K0–K7, 2,8 s).

**Convención de marcas:** [T] teorema / verificado (prueba escrita o cálculo
exacto por dos rutas); [N] no tenido; [M] medición; [P-rev] pendiente de
re-verificación contra un doc no disponible en local. MEDICIÓN ≠ PRUEBA en
todo el documento.

## §0. Alcance y método

La auditoría del doc 108 dejó el teorema condicional de la línea en la forma
(DIST-𝒟⁺) + (H2) ⇒ UNIF-U-COND ⇒ C_k + C_{k+1} > 0 para c < c* = 0,039884,
con dos huecos bloqueantes de familias distintas: G1/(DIST-𝒟⁺) (distribución)
y G2/(H2) (coeficientes). Este doc ataca (H2): lo enuncia con precisión (§1),
lo mide finamente (§2), intenta su prueba y la lleva hasta el paso exacto que
falta (§3), y confirma su independencia de (DIST-𝒟⁺) (§4).

Trabajo en el **modo de límite fijo** (la segunda reducción del doc 94,
certificada [T] en el doc 95 §2.3): cada d tiene un único primo de
conmutación s. En este modo el soporte V_s = V⁺(D_A/s) (Rosser–Iwaniec,
β = 2, Teorema 9.3 de Nathanson, como en los docs 95/97) depende de s solo
a través del nivel L = D_A//s, luego **Q_s = Q^{(L)}** [T por construcción;
verificado numéricamente en K0 contra el Model97 del doc 97, error 0,00e+00]:
el sup sobre s es un sup sobre niveles L ∈ [D_A/y_A, D_A/z_A], y el peor caso
es el nivel más pequeño. La versión de límite variable exigiría el mismo
enunciado para Q_s^full con la cola s′>s del doc 95; qué versión usa el
doc 94 es [P-rev] (doc 108, G11) y no cambia el objeto: la cota pedida es la
misma forma de Gram.

## §1. (H2) enunciada con precisión

### §1.1. Los objetos

Para un primo de conmutación s (z_A ≤ s < y_A, (s, AN) = 1):

- T_s(h0) = Σ_{d∈V⁺(D_A/s), h0|d} μ(d)/φ(d), y
  **Q_s = Σ_{h0} φ(h0)|T_s(h0)|²** (doc 94, vía doc 95 §2).
- Identidad de Gram [T, doc 95 §2]: Q_s = Σ_{d,e} μ(d)μ(e)(d,e)/(φ(d)φ(e)).
- Diagonalización [T, doc 97 §2.1]: **Q_s = Σ_r φ(r) C_s(r)²**, con
  **C_s(r) = Σ_{d∈V⁺(D_A/s), r|d} μ(d)/φ(d)**, la suma sobre r libres de
  cuadrados con primos < z_A y (r, AN) = 1.
- **Q_completo(z; N)** = ∏_{p<z, p∤AN} (1 − 1/(p−1) + 1/(p−1)²): el valor
  de Q_s si se quita la restricción V⁺ (identidad exacta [T, doc 97 §2.3]).
  Por Mertens, Q_completo(z;N)·ln z → C_loc(N) :=
  e^{−γ}·∏_{p|AN}(1−1/p)^{−1}·C₂(AN), con C_loc = 1,7225 para (A=3, N=2^e)
  [T, doc 97 §2.3]: Q_completo ≍ C_loc(N)/ln z, con constante **dependiente
  de N** a través de sus factores locales.

Aquí A es impar y fijo (la línea está construida a A fijo; la uniformidad
en A es el hueco G6 del doc 108, separado), N es par, (A, N) = 1, y N entra
en (H2) **únicamente** por el conjunto de primos excluidos {p : p | AN}.

### §1.2. Las dos formulaciones y cuál adopto

- **Forma del doc 108 (G2):** sup_{z_A≤s<y_A} Q_s ≪ 1/log z_A, con constante
  absoluta.
- **Forma del doc 97 (§7.2):** sup_{z_A≤s<y_A} Q_s ≪ Q_completo(z_A; N),
  i.e. sup_s Q_s/Q_completo(z_A; N) ≪ 1.

**Adopto la forma del doc 97**, por dos razones demostradas en este doc:

1. **La constante absoluta es falsa uniformemente en N** [T, K7]. Si todo
   primo p < z divide a AN (N primorial), el soporte se vacía: V_s = {1},
   Q_s = 1, y también Q_completo(z;N) = 1 (producto vacío), luego la forma
   normalizada da R = 1 exactamente, mientras que 1 ≫ C/ln z para toda
   constante absoluta C cuando z crece (a z=29 ya 3/ln 29 = 0,891 < 1).
   El propio doc 108 (§2.3) registra la forma corregida; su «1/log z_A» es
   una abreviatura válida solo cuando C_loc(N) está acotada (p. ej. N = 2^e).
2. **Sin condición de nivel, el enunciado es falso incluso a N fijo**
   [T, K4]. Si D_A/s es pequeño, V⁺(D_A/s) = {1} y Q_s = 1 ≫ Q_completo.
   La forma correcta lleva la hipótesis de nivel del marco,
   u_min := ln(D_A/y_A)/ln z_A ≥ u₀. En el marco auditado
   (D_A ≍ N^{1/2}, y_A = N^{1/3}, doc 108 §1.2) D_A/y_A ≍ N^{1/6}, así que
   la hipótesis es parte del marco, no una petición nueva; pero debe
   escribirse: es exactamente el análogo en el lado de coeficientes de la
   condición de nivel de cualquier lema fundamental.

**Enunciado adoptado (H2).** Existen constantes absolutas C y u₀ tales que
para todo marco (A impar fijo, N par, (A,N)=1, z_A < y_A ≤ D_A) con
ln(D_A/y_A)/ln z_A ≥ u₀:

  sup_{s primo, z_A ≤ s < y_A, (s,AN)=1} Q_s ≤ C · Q_completo(z_A; N).

### §1.3. Qué tasa necesita realmente el puente

Menos de lo que (H2) ofrece. El doc 108 (§2.3) es explícito: el o(·) de
(89.8) exige **G_A^{grande} = o(1)**, y nada más — cualquier tasa de
decaimiento de G_A sirve. La reducción del doc 94/95 convierte
sup_s Q_s ≪ Q_completo en G_A ≪ Q_completo(z_A;N) ≍ C_loc/ln z_A (la cota
implicada para G_A queda a un factor 1,05–1,41 del valor real [M, doc 95
§2.5]; la forma exacta de la reducción es [P-rev], doc 108 G4), y eso es
o(1) en el marco (z_A → ∞ con C_loc = o(ln z_A)). Recíprocamente, la tasa
de (H2) es **aguda**: no hay o(1/ln z_A) disponible (Q_s·ln z no decae
[M, doc 95 §2.5 y §2 de este doc]), así que (H2) es exactamente el
enunciado que el puente puede consumir: ni más débil (o(1) abstracto sin
mecanismo) ni más fuerte (tasas mejores son falsas).

## §2. Medición fina [M]

Script `109_ataque_h2_gram.py`, chequeos K0–K7, todo PASS (2,8 s).
K0 valida el modelo ligero contra el Model97 del doc 97 (S1 y S4: error
0,00e+00) y reproduce la identidad sin criba Q = Q_completo (relerr
≤ 3,8e-16). El estadístico central es **R(z) := sup_s Q_s/Q_completo(z;N)**
sobre los niveles muestreados (≤ 60 primos s log-uniformes por marco; en
modo fijo Q_s solo depende del nivel).

### §2.1. Barrido en z: R(z) es plano donde el nivel acompaña (K1/K2)

Marco A (D = z⁵, y = z^{10/3}, u_min = 5/3, el marco del doc 95 §2.5) y
marco B (D = z⁶, y = z³, u_min = 3), A=3, E0 (N = 2^e):

| z | R(z) marco A | R(z) marco B | Q_completo(z) |
|---|---|---|---|
| 11 | 1,429 | 1,000 | 0,6997 |
| 13 | 1,571 | 1,010 | 0,6367 |
| 17 | 1,701 | 1,032 | 0,5880 |
| 19 | 1,468 | 1,062 | 0,5536 |
| 23 | 1,549 | 1,064 | 0,5245 |
| 29 | 1,619 | 1,055 | 0,5018 |
| 31 | 1,677 | 1,069 | 0,4845 |
| 37 | 1,492 | 1,060 | 0,4689 |

- **Marco B (u ≥ 3): R(z) ∈ [1,000, 1,069]** en todo z = 11…37, sin
  tendencia (desde z=17: [1,032, 1,069]). Es la medición más limpia de la
  serie sobre (H2): con nivel suficiente, sup_s Q_s = (1,00–1,07)·Q_completo.
- **Marco A (u ≥ 5/3): R(z) ∈ [1,43, 1,70]**, y el sup se alcanza en la
  **meseta de niveles pequeños**: a z=11, 13, 17 el sup es Q_s = 1,000000
  exacto (el nivel mínimo L = D/y es tan pequeño que V⁺(L) = {1}: a z=11,
  L ≈ 54 < 5³, ni el primo 5 entra en el soporte); a z=19…31 el sup es
  0,8125 y a z=37 es 0,6997. Los rangos de Q_s·ln z reproducen los del
  doc 95 §2.5 (mín ≈ 1,63–1,70; máx 2,39–2,83) — validación cruzada ✓.
- Ajustes de tasa (mínimos cuadrados sobre ln ln z): γ̂ = 0,82 (A) y 0,80 (B)
  para sup Q_s ~ (ln z)^{−γ}; δ̂ = 0,14 (A) y 0,16 (B) para R(z) ~ (ln z)^δ.
  **Lectura honesta:** γ̂ < 1 no es evidencia contra γ = 1: el ajuste está
  sesgado a la baja por la meseta Q_s = 1 a z pequeño (marco A) y por la
  convergencia lenta de Mertens de Q_completo·ln z (1,678 a z=11 →
  1,7225 asintótico [T, doc 97]). El estadístico decisivo es R(z) mismo:
  plano. La medición es **compatible con la tasa requerida** y no sugiere
  que (H2) — en la forma adoptada de §1.2 — sea falsa; sugiere exactamente
  lo contrario en el marco B.

### §2.2. Dependencia en A y en el factor local de N (K3, z=29, marco A)

| caso | Q_completo | sup Q_s | R | sup Q_s·ln z |
|---|---|---|---|---|
| A=3, E0 | 0,5018 | 0,8125 | 1,619 | 2,736 |
| A=15, E0 | 0,6176 | 1,0000 | 1,619 | 3,367 |
| A=105, E0 | 0,7172 | 1,0000 | 1,394 | 3,367 |
| A=3, E1 (5,7,11,13 \| N) | 0,8533 | 1,0000 | 1,172 | 3,367 |

sup Q_s·ln z varía de 2,74 a 3,37 con el factor local (la «constante» de
la forma sin normalizar no existe), mientras R permanece en [1,17, 1,62]:
**la normalización por Q_completo(z;N) es la correcta**, como el doc 97
anticipó y el K7 demuestra.

### §2.3. Desglose por niveles: R(L) decrece de 1/Q_completo a 1 (K4)

z=17, A=3, D=2·10⁵, exhaustivo en L (u = ln L/ln z):

| banda u | 0–1,5 | 2,0 | 2,5 | 3,0 | 3,5 |
|---|---|---|---|---|---|
| máx R(L) | 1,701 | 1,382 | 1,190 | 1,032 | 1,018 |

El valor 1,701 de las bandas bajas es exactamente 1/Q_completo(17) = 1/0,588:
la meseta Q^{(L)} = 1 de los niveles con V = {1}. El umbral empírico para
R(L) ≤ 1,5 en todo L mayor es u > 1,70, y el peor infractor es L = 125 = 5³:
el nivel en el que el primo 5 entra en el soporte V⁺. Moraleja cuantitativa:
la condición de nivel de (H2) no es un tecnicismo — R(L) es una función
decreciente y explícita de u, con R → 1 cuando u crece. Un u₀ = 3 basta en
ventana para R ≤ 1,07 (§2.1, marco B).

## §3. Intento de prueba: reducción exacta a (H2′), y las rutas muertas

### §3.1. La factorización y la descomposición del defecto [T]

**Paso 1 [T].** Sea r | P(z_A) libre de cuadrados. Todo d con r | d se
escribe d = r·d′ con (d′, r) = 1 (d libre de cuadrados), y
μ(d)/φ(d) = μ(r)μ(d′)/(φ(r)φ(d′)). Luego

  C_s(r) = (μ(r)/φ(r)) · S_s(r),
  S_s(r) := Σ_{d′∈V_s^{(r)}} μ(d′)/φ(d′),

donde V_s^{(r)} = {d′ : (d′,r)=1, r·d′ ∈ V⁺(D_A/s), primos(d′) < z_A,
(d′, AN) = 1} es el **soporte desplazado por r**: la condición de Nathanson
se impone sobre el producto r·d′, no sobre d′ solo.

**Paso 2 [T].** Sin restricción V⁺, S_s(r) sería
Π(r) := ∏_{p<z_A, p∤rAN} (1 − 1/(p−1)), y Σ_r Π(r)²/φ(r) = Q_completo
(doc 97 §2.3). Definiendo el **defecto** δ_s(r) := S_s(r) − Π(r) se tiene
φ(r)C_s(r)² = (Π(r) + δ_s(r))²/φ(r) (verificado a 6,9e-18 en K5), y por
(a+b)² ≤ 2a² + 2b²:

  **Q_s ≤ 2·Q_completo(z_A;N) + 2·Σ_r δ_s(r)²/φ(r).** [T]

Por tanto **(H2) ⇐ (H2′)**, donde el sub-enunciado restante es:

> **(H2′).** Existen C, u₀ absolutas tales que, bajo la condición de nivel
> ln(D_A/y_A)/ln z_A ≥ u₀,
> sup_{z_A≤s<y_A} Σ_{r|P(z_A), (r,AN)=1} δ_s(r)²/φ(r) ≤ C·Q_completo(z_A;N).

Medición de (H2′) [M, K5]: el máximo sobre niveles de Σ_r δ²/φ es
0,947·Q_completo (z=17), 0,738·Q_completo (z=29), 0,553·Q_completo (z=37)
en el marco A: es O(Q_completo) y decrece con z. (H2′) es compatible con
la medición; probarla es el paso que falta, y §3.2–§3.4 delimitan por qué
no es un paso rutinario.

### §3.2. La ruta del triángulo está muerta [T, con números en K6]

δ_s(r) = −Σ_{d′∉V_s^{(r)}} μ(d′)/φ(d′) (la cola fuera del soporte). Por
triángulo, |C_s(r)| ≤ (1/φ(r))·Σ_{d′} 1/φ(d′) ≤ (1/φ(r))·P⁺ con
P⁺ := ∏_{p<z, p∤AN}(1 + 1/(p−1)) ≍ e^γ·ln z·(factor local), y entonces
Q_s ≤ (P⁺)³ ≍ (ln z)³. La pérdida contra el objetivo:
(P⁺)³/Q_completo = 4,43 (z=11) → 22,12 (z=37), creciendo como (ln z)⁴/C
[K6]. En el peor nivel medido (z=37): cota de triángulo por nivel 2,77 y
global 10,37, frente a Q_s = 0,70. **Ningún mayorante absoluto sobre
C_s(r) puede dar (H2)**: la pérdida crece polinómicamente en ln z. Es el
mismo fenómeno que agotó los mayorantes absolutos en la base E (docs
103–104), ahora en la base de coeficientes.

### §3.3. La ruta diagonal d=e está muerta [T] — y una errata al doc 95

La diagonal de la forma de Gram es Σ_{d∈V_s} d/φ(d)². Sobre el soporte
**sin restricción** vale exactamente, por multiplicatividad,

  Σ_{d|P(z)} d/φ(d)² = ∏_{p<z, p∤AN} (1 + p/(p−1)²) ≍ C₃·ln z,

con C₃ = 0,667 medido a z=37 (producto = 2,4068) [T la identidad; M la
constante]. Y sobre cualquier V_s, diag ≥ 1 siempre (el término d=1).
Como Q_completo ≍ C_loc/ln z → 0, **ninguna cota de la forma Q_s ≤ diag
puede implicar (H2)**, ni siquiera concediendo sin prueba que los términos
cruzados d≠e son netos negativos: la diagonal es Θ(ln z) en el peor caso
y ≥ 1 en todo caso, a un factor ~(ln z)² del objetivo.

**Errata (al doc 95 §2.5, con efecto en el doc 97 §7.2).** El doc 95
afirmó: «La diagonal sola ya es O(1/ln z): diag·ln z ≈ 3,8–7,5», y el
doc 97 §7.2 propuso como primer paso «probar la cota para la diagonal
d=e (que ya da el orden)». La lectura es inconsistente con la propia
tabla del doc 95: diag·ln z **crece** de 3,76 (z=11) a 7,54 (z=29),
siguiendo 0,667·(ln z)² (predicción: 3,84 y 7,56) — es decir
diag ≍ 0,667·ln z, no O(1/ln z). Lo que el doc 95 midió como
«diag/Q_s ≈ 2,2–4,3» es un cociente en ventana: asintóticamente
diag/Q_completo ~ (0,667/C_loc)·(ln z)² → ∞. La cancelación fuera de la
diagonal no es una mejora de constante de un factor 2–4: **es el orden
entero de (H2)**. La sugerencia del doc 97 §7.2 queda retirada por este
doc; la conclusión del doc 95 §2.5 («el objetivo sobrevive como orden de
magnitud, constante C ≈ 2,5») se mantiene, pero su justificación vía la
diagonal no.

### §3.4. La forma puntual de (H2′) es falsa en ventana [M, K5]

El análogo puntual del lema fundamental,
|δ_s(r)| ≤ Π(r)·ρ(u) con ρ(u) < 1 (p. ej. ρ = e^{−u/2}),
**no se sostiene**: máx_r |δ_s(r)|/Π(r) = 1,000 (z=17), 1,069 (z=29),
1,295 (z=37) en el marco A. El valor 1,000 es estructural: para todo r
que no divide a ningún d ∈ V_s se tiene C_s(r) = 0, i.e. δ_s(r) = −Π(r)
exactamente. Esos r son la cola de la propia suma de Q_completo, y su
masa Π(r)²/φ(r) es la que el promedio L² de (H2′) sí controla. Conclusión:
**(H2′) en forma L² es el sub-enunciado mínimo viable**; cualquier intento
de prueba debe ser en promedio sobre r, no uniforme en r.

### §3.5. El paso exacto que falta, y por qué el lema fundamental clásico no basta

Lo que falta es un **lema fundamental uniforme sobre los soportes
desplazados** {V_s^{(r)}}: cancelación en Σ μ(d′)/φ(d′) sobre el
complemento de V_s^{(r)}, uniforme en r | P(z_A) y en el nivel, en la
norma L²(Σ_r ·/φ(r)). El lema fundamental clásico (Halberstam–Richert;
Iwaniec–Kowalski) no lo da literalmente, por tres razones concretas:

1. Controla sumas **ponderadas por λ^±** (Σ_d λ_d/d y análogas) o medias
   de μ sobre d | P(z) con d < D; aquí la suma es μ(d′)/φ(d′) **sin
   pesos** sobre el soporte de λ⁺ — el soporte y los pesos no son
   intercambiables sin un argumento nuevo.
2. El soporte V_s^{(r)} está definido por la condición de Nathanson
   sobre el producto r·d′ (los primos de r entran en la ordenación
   decreciente y en los productos parciales): no es un soporte de
   Rosser–Iwaniec de parámetros (D_A/(sr), z_A) en ninguna forma
   verificada, y la dependencia en r es exactamente lo que hay que
   controlar uniformemente.
3. La conclusión necesaria es un promedio L² sobre ~2^{π(z_A)} valores
   de r, no una cota para un único conjunto de cribado.

**Estado de la prueba:** (H2) está reducida, por pasos [T], a (H2′);
(H2′) está apoyada por la medición (K5) y delimitada por dos rutas
probadamente muertas (§§3.2–3.3) y una forma puntual probadamente falsa
en ventana (§3.4). No declaro (H2) demostrada.

## §4. Independencia de (DIST-𝒟⁺) [T]

Confirmada, por inspección de las definiciones: Q_s, C_s(r), S_s(r),
Π(r) y δ_s(r) se construyen solo con μ, φ, el soporte combinatorio V⁺ y
los parámetros (A, z_A, y_A, D_A) más el conjunto de primos que dividen
a AN. **E(N; d, N) = ψ(N; d, N) − ψ(N)/φ(d) no aparece en ningún paso**;
ψ no aparece en absoluto (el script de este doc no contiene ninguna
suma sobre primos en progresiones: N entra solo como conjunto de primos
excluidos). (H2)/(H2′) es una propiedad de la forma de Gram de los
coeficientes del selector b_A; (DIST-𝒟⁺) es una propiedad de la
distribución de los primos en la clase a ≡ N sobre el conjunto delgado
𝒟. Ninguno de los dos enunciados implica formalmente al otro: (H2) es,
para cada tupla de parámetros, una afirmación finita sobre sumas de
Möbius en soportes de criba, decidible sin ningún input de distribución
de primos. Esto confirma la clasificación del doc 108: G2 vive en el
lado de los coeficientes; G1/(DIST-𝒟⁺), en el de los errores.

## §5. Veredicto

1. **(H2) queda fijada con precisión** en la forma del doc 97
   (sup_s Q_s ≪ Q_completo(z_A;N)), con dos reparaciones necesarias y
   demostradas: la normalización por el producto local (la forma
   1/log z_A es falsa uniformemente en N: contraejemplo primorial, K7) y
   la condición de nivel u_min ≥ u₀ (sin ella Q_s = 1 en niveles
   pequeños, K4). Al puente le basta G_A^{grande} = o(1); la tasa de
   (H2) es la aguda disponible.
2. **La medición apoya (H2) en la forma adoptada**: R(z) plano, ∈
   [1,00, 1,07] con u ≥ 3 y ≤ 1,70 con u ≥ 5/3 (z = 11…37, K1–K4,
   validación cruzada exacta contra el doc 97 en K0). Ninguna señal de
   falsedad; el sup vive en la meseta de niveles pequeños y R(L) → 1 al
   crecer u.
3. **La prueba avanza hasta una reducción exacta [T]**: (H2) ⇐ (H2′),
   la cota L² de los defectos δ_s(r) = S_s(r) − Π(r), con
   Q_s ≤ 2·Q_completo + 2·Σ_r δ_s(r)²/φ(r).
4. **Delimitación negativa completa**: triángulo y diagonal d=e están
   probadamente muertos (pérdidas (ln z)⁴/C y (ln z)²; errata al
   doc 95 §2.5 y retirada de la sugerencia del doc 97 §7.2), y la forma
   puntual de (H2′) es falsa en ventana (máx |δ|/Π = 1,00–1,30, K5).
   El paso que falta es un lema fundamental uniforme sobre los soportes
   desplazados V_s^{(r)}, en promedio L²: la FL clásica no aplica
   literalmente (soportes acoplados por r; sumas μ/φ sin pesos λ^±).
5. **Independencia confirmada**: (H2)/(H2′) no involucran E(N;·,·) ni
   ψ; es un hueco de coeficientes, independiente de (DIST-𝒟⁺).
   (H2) sigue abierta, reducida a (H2′). Goldbach no se declara
   demostrada.

## Apéndice. Cómo correr la verificación

```
cd ~/workspace/goldbach
python3 109_ataque_h2_gram.py | tee salida_doc109_local.txt
```

Solo stdlib (math, sys, time, importlib para la validación cruzada con
el script del doc 97). Tiempo: ~2,8 s. Salida: PASS, chequeos:

- K0: validación cruzada del modelo ligero contra Model97 del doc 97
  (S1: error 0,00e+00 en los 13 s; S4: 0,00e+00 en muestra) e identidad
  sin criba Q = Q_completo (relerr ≤ 3,8e-16, 3 casos).
- K1: barrido en z del marco A (D=z⁵, y=z^{10/3}): sup Q_s, R(z),
  Q_s·ln z (cf. doc 95 §2.5: mín 1,63–1,70, máx 2,39–2,83 ✓).
- K2: marco B (D=z⁶, y=z³) y ajustes γ̂, δ̂ sobre ln ln z.
- K3: A ∈ {3, 15, 105} y E0/E1 a z=29: R estable, sup Q_s·ln z no.
- K4: desglose exhaustivo por niveles a z=17: bandas de R(L) por u,
  umbral u > 1,70 para R ≤ 1,5, Q^{(2)} = 1 exacto.
- K5: defectos: Σ_r δ²/φ = 0,947/0,738/0,553 × Q_completo (z=17/29/37);
  máx |δ|/Π = 1,000/1,069/1,295; identidad de descomposición a 6,9e-18.
- K6: rutas absolutas a z=37: (P⁺)³/Q_completo = 4,43 (z=11) → 22,12
  (z=37); diag sin restricción = ∏(1+p/(p−1)²) = 0,667·ln z.
- K7: contraejemplo primorial a z=29: Q_s = Q_completo = 1, R = 1,
  3/ln z = 0,891 < 1.

---
*Goldbach sigue abierta. (H2) queda enunciada en su forma correcta,
medida en su favor, y reducida a un único sub-enunciado (H2′): un lema
fundamental en promedio L² sobre los soportes desplazados por r. Es un
hueco de coeficientes, independiente del muro de distribución
(DIST-𝒟⁺).*
