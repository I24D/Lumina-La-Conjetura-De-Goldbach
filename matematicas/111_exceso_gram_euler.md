# DOC 111 — El exceso de Gram, conductor a conductor: la masa positiva es el contenido mínimo de (H2)

**Fecha:** 2026-11-28.
**Autor:** Muse (línea dedicada de Goldbach).
**Archivos:** `111_exceso_gram_euler.md` (este doc),
`111_exceso_gram_euler.py` (script, solo stdlib),
`salida_doc111_local.txt` (salida completa, **PASS**, chequeos K0–K6, 2,3 s).

**Convención de marcas:** [T] teorema / verificado (prueba escrita o cálculo
exacto por dos rutas); [N] no tenido; [M] medición; [P-rev] pendiente de
re-verificación contra un doc no disponible en local. MEDICIÓN ≠ PRUEBA en
todo el documento.

## §0. Alcance y método

El doc 110 disolvió (H2′): por la identidad (110.2),
Def2 = (Q^{(L)} − Q_completo) + 2·Π_full·M_cola con el término de conteo
cerrado [T], (H2′) es (H2) disfrazado, y el contenido mínimo restante es
**(H2″)**: sup_L[Q^{(L)} − Q_completo] ≤ C·Q_completo — un lema fundamental
para la forma de Gram misma. Este doc ataca (H2″) por los dos lados que el
doc 110 dejó indicados: la **diagonalización** (docs 96/97), escribiendo el
exceso conductor a conductor y factorizando la diferencia de cuadrados, y
la **factorización de Euler** (doc 110: Q_completo = Π_full·M_full),
preguntando de qué bloque de primos viene el exceso.

Trabajo en el **modo de límite fijo** del doc 109 (Q_s = Q^{(L)},
L = D_A//s), con las definiciones exactas de los docs 96/97/109/110,
reutilizadas por import (el script importa el del doc 110, que importa el
del doc 109: son literalmente las mismas funciones):

- 𝒫 = primos permitidos: impares, < z_A, ∤ AN. a(d) := μ(d)/φ(d).
- V = V⁺(L) (Nathanson, Teorema 9.3, β = 2).
- C^{(L)}(r) = Σ_{d∈V, r|d} a(d); Q^{(L)} = Σ_r φ(r)C^{(L)}(r)²
  [T, doc 97 §2.1].
- C^full(r) = (μ(r)/φ(r))·Π(r), con Π(r) = ∏_{p∈𝒫∖primos(r)}(1 − 1/(p−1)):
  el coeficiente sin criba [T, doc 97 §2.3], y
  Q_completo = Σ_r Π(r)²/φ(r) = ∏_{p∈𝒫}(1 − 1/(p−1) + 1/(p−1)²).
- S(r) := φ(r)C^{(L)}(r)/μ(r); δ(r) = S(r) − Π(r); Def2 = Σ_r δ(r)²/φ(r).

Defino el **exceso** y su término por conductor:

  E^{(L)} := Q^{(L)} − Q_completo,
  **t(r) := φ(r)·[ C^{(L)}(r)² − C^full(r)² ].**

## §1. La expresión del exceso y su estructura de signos

### §1.1. La expresión y la factorización [T]

Por la diagonalización aplicada a las dos formas (la cribada y la
completa, que comparten el índice r y la base φ(r)):

  **E^{(L)} = Σ_r t(r) = Σ_r [ S(r)² − Π(r)² ] / φ(r).** [T]

La diferencia de cuadrados se factoriza con δ(r) = S(r) − Π(r):

  **t(r) = δ(r)·( S(r) + Π(r) ) / φ(r) = δ(r)·( 2Π(r) + δ(r) ) / φ(r).** [T]

Verificación [T, K1]: Σ_r t(r) reproduce E^{(L)} con error 0,00e+00 en los
niveles que deciden (y ≤ 2e-16 en todos los examinados); la factorización
vale por r con error ≤ 3,5e-18; y la expresión es consistente con (110.2):
E^{(L)} = Def2 − 2·Π_full·M_cola (error ≤ 8,3e-16), es decir el exceso es
también la **diferencia de dos términos ≥ 0**: el momento de los defectos
menos la cola de conteo.

### §1.2. No hay signo definido por conductor [T + M]

Como Π(r) > 0 siempre (cada factor 1 − 1/(p−1) ∈ [1/2, 1) para p ≥ 3), el
signo de t(r) es sign(δ(r))·sign(S(r) + Π(r)), y ninguna de las dos cosas
tiene signo fijo en general. La cancelación entre conductores es real y
del orden del objetivo. Descomponiendo

  P^{(L)} := Σ_{t(r)>0} t(r),   N^{(L)} := −Σ_{t(r)<0} t(r),
  E^{(L)} = P^{(L)} − N^{(L)},

se tienen dos hechos exactos:

- **El bloque invisible es exacto y no positivo [T].** Si r no divide a
  ningún d ∈ V (r *invisible*), C^{(L)}(r) = 0, S(r) = 0 y
  **t(r) = −Π(r)²/φ(r)**, sin error (K2: error 0,0e+00). Este bloque solo
  ayuda: resta.
- **La masa negativa total está acotada por el propio normalizador [T].**
  Si t(r) < 0 entonces |t(r)| = (Π(r)² − S(r)²)/φ(r) ≤ Π(r)²/φ(r), luego
  **0 ≤ N^{(L)} ≤ Σ_r Π(r)²/φ(r) = Q_completo**, sin hipótesis alguna.

Y un hecho medido que simplifica la lectura en ventana:

- **S(r) ≥ 0 en toda la ventana [M, K2b].** Sobre los 240 niveles
  muestreados de los marcos A (z = 17, 29, 37) y B (z = 29)
  —46 000 instancias de r en total— S(r) < 0 **no ocurre nunca**
  (mínimo medido: −0,0000, redondeo). En consecuencia, en ventana
  sign(t(r)) = sign(δ(r)): t(r) > 0 ⟺ δ(r) > 0, y el canal δ < −2Π está
  vacío (K6: la parte de P con δ < −2Π es 0,000 en los cuatro marcos).
  Los defectos pequeños negativos (0 < S(r) < Π(r)) existen pero son
  raros en el marco A (16/52/90 instancias) y frecuentes en el B (858).
  Esto es MEDICIÓN, no teorema: no afirmo S(r) ≥ 0 en general.

### §1.3. Dónde está el exceso en los niveles que deciden [M]

| marco | z | L | Q^{(L)} | E/Q_c | P/Q_c (nº de r) | N/Q_c | invisible/Q_c |
|---|---|---|---|---|---|---|---|
| A | 17 | 112 | 1,0000 | +0,701 | 1,248 (1) | 0,548 | −0,548 |
| A | 29 | 273 | 0,8125 | +0,619 | 1,073 (2) | 0,453 | −0,453 |
| A | 37 | 410 | 0,6997 | +0,492 | 0,862 (4) | 0,370 | −0,370 |
| B | 29 | 24399 | 0,5292 | +0,055 | 0,159 (28) | 0,104 | −0,019 |
| A/B | — | L máx. | — | +0,000…0,003 | ≤ 0,017 | ≤ 0,013 | ≈ 0 |

En los niveles que deciden el marco A el soporte es minúsculo y todo se
puede escribir a mano (verificado por la ruta de C y por la de los
desplazados): a z=17, L=112, V = {1} (la meseta del doc 109: ni el primo 5
supera p³ < 112), Q^{(L)} = 1, y P es un único término,
**t(1) = S(1)² − Π_full² = 1 − Π_full² = 1,248·Q_c**; toda la masa
negativa es el bloque invisible. A z=29, L=273, V = {1, 5} y
P = t(1) + t(5) = 0,743 + 0,330 (en Q_c). A z=37, L=410, V = {1, 5, 7, 35}
y P = t(1) + t(5) + t(7) + t(35) = 0,481 + 0,214 + 0,116 + 0,051.
En los tres casos **todo r visible tiene δ(r) > 0 y todo r con t(r) < 0
es invisible**: la masa negativa es, bloque a bloque de P⁺(r) y de ω(r),
exactamente el Def2 de los invisibles (K3/K4: neto t = −neto Def2 en cada
bloque negativo). En el marco B (nivel real, u ≥ 3) el cuadro cambia:
P se reparte entre 28 conductores, N = 0,104·Q_c del que solo 0,019 es
invisible — hay cancelación genuina entre conductores visibles — y el
exceso total es 0,055·Q_c.

**Advertencia de ventana.** La concentración de P en r=1 a u ≈ 5/3 es un
artefacto de soporte diminuto: a nivel L = z^{5/3} un primo entra en V
solo si p < z^{5/9}, y a z ≤ 37 eso deja fuera a casi todo 𝒫. Al crecer z
el soporte se puebla y P/Q_c *decrece* en el marco A (1,248 → 1,073 →
0,862, K6, máximo sobre los 60 niveles de cada z). El fenómeno asintótico
no es «r=1 manda» sino el de §3.

## §2. La cola de Euler: de qué primos viene el exceso (y de cuáles no)

La pregunta tenía tres lecturas posibles; las respondo por separado.

### §2.1. La cola de Euler del normalizador no es el exceso [T]

Q_completo se parte en cualquier corte w como
Q_completo = Q_completo^{≤w} · Q_completo^{>w}. Con w = z^{1/2} [M, K5]:

| z | Q_c^{≤√z} | Q_c^{>√z} | cuota de log de la cola |
|---|---|---|---|
| 17 | 1,0000 | 0,5880 | 100% (𝒫 empieza en 5 > √17) |
| 29 | 0,8125 | 0,6176 | 69,9% |
| 37 | 0,8125 | 0,5771 | 72,6% |

El factor de cola es ≈ 0,58–0,62, **acotado lejos de 1**: es el cociente
de Mertens (ln √z)/(ln z) = 1/2 disfrazado. El producto de Q_completo es
condicionalmente convergente (tipo ∏(1 − 1/p)): su «cola de primos
grandes» no es 1 + o(1) ni un resto despreciable — es parte del
normalizador, no del exceso. En cambio el otro factor de Euler del
doc 110 sí converge absolutamente, con cota explícita [T]:

  M_full^{>w} = ∏_{p>w}(1 + 1/((p−1)(p−2))) ≤ exp( Σ_{p>w} 1/((p−1)(p−2)) )
              ≤ **exp(1/(⌊w⌋−1))**,

porque Σ_{p>w} 1/((p−1)(p−2)) ≤ Σ_{m>w} 1/((m−1)(m−2)) = 1/(w−1)
(telescopio sobre enteros). Medido [M, K5]: M_full^{>√z} = 1,1405 /
1,0629 / 1,0655 a z = 17/29/37, bajo la cota 1,396 / 1,284 / 1,221.
Y su resto a nivel, M_cola, ya está cerrado por conteo [T, doc 110]:
por bloques de P⁺(d) (z=29, L=273) M_V vive en P⁺(d) ∈ {1, 5} —los d de V—
y M_cola vive entero en P⁺(d) ≥ 7 (0,0361 + 0,0124 + 0,0086 + ⋯,
decreciente): es exactamente la cola de los d que la condición de
Nathanson excluye a ese nivel, y está acotada entre 0 y M_full. **Ninguna
parte abierta del exceso es una cola de Euler.**

### §2.2. Los primos del soporte, por bloques [M]

Por primo mayor P⁺(r) (K3) y por ω(r) (K4), en los niveles de §1.3:

- Marco A: el bloque P⁺(r) = 1 (es decir r = 1) aporta +1,248 / +0,743 /
  +0,481 de Q_c — iguala o excede el exceso neto —; los demás bloques
  positivos son los de los primos visibles pequeños (5, y 7 a z=37);
  todos los bloques con P⁺(r) ≥ 7 (z=29) o ≥ 11 (z=37) son negativos y
  coinciden con el Def2 invisible. Por ω: el exceso positivo es ω = 0 y
  parte de ω = 1; ω ≥ 2 es neto negativo en el marco A.
- Marco B: ningún bloque manda; los netos por P⁺(r) son ±0,01–0,03 de Q_c
  con cancelación intra-bloque (pos y neg del mismo orden desde p=11).

Es decir: el exceso **no viene de los primos grandes del soporte** — sus
bloques restan (son los invisibles) — ni de una cola de Euler. En ventana
viene de los conductores visibles más pequeños, con r=1 a la cabeza en el
marco A: t(1) = S(1)² − Π_full², donde S(1) = Σ_{d∈V⁺(L)} a(d) es la suma
de tipo lema fundamental del propio soporte (medido: S(1) = 1 / 0,75 /
0,625 frente a Q_c^{1/2} = 0,767 / 0,708 / 0,685).

### §2.3. Los primos del nivel (z ≤ s < y_A): el sup lo fijan ellos [M]

Los primos entre z_A e y_A no dividen a r: entran solo por el nivel
L = D_A//s. El perfil del exceso en el nivel es el del doc 109 leído al
revés: E/Q_c es máximo en el nivel mínimo (s máximo, u = u_min) y
≈ 0,000–0,003 en el nivel máximo de cada marco (K1/K2); con u ≥ 3
(marco B) E/Q_c ≤ 0,055 en todo el barrido. **El sup de (H2″) lo deciden
los primos de conmutación grandes vía niveles pequeños**; es la condición
de nivel de (H2), no un bloque de primos del soporte, lo que gobierna el
exceso.

## §3. Intento de cota por bloques: lo cerrado y lo que se resiste

### §3.1. El libro mayor [T]

De E^{(L)} = P^{(L)} − N^{(L)} y de (110.2), las piezas cerradas por
conteo o por estructura exacta son:

1. **Bloque invisible:** Σ_{r invisible} t(r) = −Σ_{r invisible} Π(r)²/φ(r)
   ≤ 0, exacto [T, §1.2]. Cerrado.
2. **Masa negativa:** 0 ≤ N^{(L)} ≤ Q_completo [T, §1.2]. Cerrada como
   cota (no hace falta conocerla: solo acota lo que P debe compensar).
3. **Cola de conteo:** 0 ≤ 2·Π_full·M_cola ≤ 2·Q_completo [T, doc 110],
   con M_full ≤ 1,1618 absoluto a A=3. Cerrada.
4. **Cola de Euler de M:** M_full^{>w} ≤ exp(1/(⌊w⌋−1)) [T, §2.1].
   Cerrada.

### §3.2. La pieza que se resiste: la masa positiva P^{(L)}

Todo lo demás es **P^{(L)} = Σ_{t(r)>0} t(r)**, en ventana
P^{(L)} = Σ_{δ(r)>0} (2Π(r)δ(r) + δ(r)²)/φ(r) [M, §1.2]. Es la única
pieza sin cota: los intentos disponibles fallan por las razones ya
delimitadas en los docs 109–110 —

- P^{(L)} ≤ Σ_{r visible} S(r)²/φ(r) = Q^{(L)}: circular, es (H2) misma.
- Mayorantes absolutos término a término: pierden (ln z)²–(ln z)⁴
  (triángulo, doc 109 §3.2; diagonal d=e ≍ 0,667·ln z, doc 109 §3.3;
  diagonal de la cola, doc 110 §3.1).
- La parte lineal positiva Σ_{δ>0} 2Π(r)δ(r)/φ(r) por separado: la suma
  lineal completa Σ_r 2Π(r)δ(r)/φ(r) = −2·Π_full·M_cola ∈ [−2Q_c, 0] está
  cerrada [T, doc 110 §1.2], pero cerrarla *por partes* exigiría controlar
  la cancelación entre los conductores con δ > 0 y con δ < 0, que es
  exactamente el contenido no contable.

El sub-enunciado mínimo restante, equivalente a (H2″):

> **(H2‴) — masa positiva del exceso.** Existen constantes absolutas C
> y u₀ tales que para todo marco (A impar fijo, N par, (A,N)=1,
> z_A < y_A ≤ D_A) con ln(D_A/y_A)/ln z_A ≥ u₀:
>
>   sup_L P^{(L)} = sup_L Σ_{r : t(r)>0} φ(r)·[ C^{(L)}(r)² − C^full(r)² ]
>                ≤ C · Q_completo(z_A; N),
>
> con L ∈ [D_A/y_A, D_A/z_A] en el modo fijo (y su análogo de límite
> variable [P-rev], como en los docs 109/110).

**Equivalencia [T].** (H2″) con constante C ⟹ (H2‴) con C+1, porque
P = E + N ≤ C·Q_c + Q_c. (H2‴) con constante C ⟹ (H2″) con C, porque
E ≤ P. Las constantes difieren en ≤ 1 aditivo, como en la equivalencia
(H2) ⟺ (H2′) del doc 110: (H2‴) no es más débil que (H2″); es su
contenido tras retirar las tres piezas cerradas de §3.1. Lo que un lema
fundamental para la forma de Gram debe probar, despojado de todo lo que
el conteo ya da, es exactamente esto: **la suma de las desviaciones
positivas de φ(r)C^{(L)}(r)² sobre su valor completo es O(Q_completo)**.

**Medición de (H2‴) [M, K6].** Máximo de P/Q_c sobre los 60 niveles de
cada marco: **1,248 (A, z=17) → 1,073 (A, z=29) → 0,862 (A, z=37)**,
decreciente en z, y **0,159 (B, z=29)**. Compatible con C absoluto
pequeño; no es una prueba.

**El núcleo r=1, y por qué no basta [M].** En el marco A el término
mayor de P es t(1) = S(1)² − Π_full² con S(1) = Σ_{d∈V⁺(L)} μ(d)/φ(d);
como Π_full² ≍ (const/ln z)² ≪ Q_c, la parte r=1 de (H2‴) es en esencia
la cota de tipo lema fundamental |Σ_{d∈V⁺(L)} μ(d)/φ(d)| ≪ Q_c^{1/2}
para la suma sin pesos sobre el soporte — medido S(1)/Q_c^{1/2} =
1,30 / 1,06 / 0,91 a z = 17/29/37 en los niveles peores. Pero en el
marco B P se reparte entre 28 conductores y t(1) aporta solo 0,028 de
0,159: **ningún enunciado sobre r=1 solo cierra (H2‴) uniformemente**;
la masa positiva es un fenómeno de todos los conductores visibles.

**Estado de la prueba:** (H2″) no queda demostrada. Queda reducida a su
contenido mínimo (H2‴) con las piezas de conteo y de signo cerradas por
separado, y con la localización medida del exceso: niveles pequeños
(primos de conmutación grandes), conductores visibles pequeños, nada de
cola de Euler.

## §4. Coherencia con el puente: qué constante necesita G2

El doc 108 (G2, §§1.3/2.3 y tabla de huecos) pide exactamente:
sup Q_s/Q_completo(z_A; N) ≪ 1 con constante absoluta, para obtener
G_A ≍ C_loc/ln z_A y G_grande = o(1), que es lo que convierte el
R_grande = O(N/(A ln²N)) del puente en el o(·) de (89.8).

- **Sí es la constante que G2 necesita [T, por composición].** (H2″) con
  constante C da Q^{(L)} ≤ Q_completo + C·Q_completo = (C+1)·Q_completo,
  es decir (H2) con constante C+1 — la forma que el doc 109 §1.2 adoptó
  y cuya suficiencia para el puente el doc 109 §1.3 ya estableció: al
  puente le basta G_A^{grande} = o(1) con *cualquier* tasa, y aquí la
  tasa es Q_completo ≍ C_loc(N)/ln z_A. El valor de C es irrelevante
  para el o(·): ninguna constante absoluta concreta (ni 1, ni 3, ni la
  que salga de (H2‴)) cambia el cierre del eslabón 1 del doc 108.
- **Lo que (H2″) no calibra.** La constante numérica 0,0225 dentro de
  T(N) = N^{2−2c}/(0,0225·ln⁴N) es f·G_A con f ≈ 0,05 la fracción del
  bloque grande según la regla de atribución (doc 108 §2.3): seguirá
  siendo una medición mientras no se fijen la forma exacta de la
  reducción sup Q_s → G_A del doc 94 (hueco G4, [P-rev]; factor medido
  1,05–1,41, doc 95 §2.5) y la regla f. (H2″) cierra la *tasa* de G2,
  no su constante impresa en T(N): si (H2″) se probara mañana, el
  eslabón 1 cerraría como o(·), pero el 0,0225 seguiría marcado [M]
  hasta saldar G4.
- **Condiciones de marco que acompañan (no son parte de (H2″)).**
  (i) C_loc(N) = o(ln z_A), para que Q_completo = o(1) (doc 109 §1.3):
  hipótesis sobre el marco (falla si N absorbe demasiados primos
  pequeños, cf. el contraejemplo primorial del doc 109 K7, donde
  Q_completo = 1 y no hay o(1) que pedir). (ii) La versión de límite
  variable de (H2″)/(H2‴) (cola s′>s del doc 95) es la misma afirmación
  módulo [P-rev] (doc 108, G11).

## §5. Veredicto

1. **Expresión exacta [T]:** E^{(L)} = Σ_r φ(r)[C^{(L)}(r)² − C^full(r)²]
   = Σ_r δ(r)(2Π(r) + δ(r))/φ(r), verificada a precisión de máquina
   (K1); el exceso **no tiene signo definido por conductor** — hay
   cancelación entre conductores del orden de Q_c —, pero el bloque
   invisible es exacto, t(r) = −Π(r)²/φ(r) ≤ 0, y la masa negativa total
   cumple 0 ≤ N^{(L)} ≤ Q_completo [T]. En ventana S(r) ≥ 0 siempre
   [M, K2b: 46 000 instancias, ningún S(r) < 0], luego t(r) > 0 ⟺ δ(r) > 0.
2. **Origen del exceso:** no es una cola de Euler (la única cola
   absolutamente convergente, la de M, está cerrada [T] con
   M_full^{>w} ≤ exp(1/(⌊w⌋−1)); la de Q_completo es un factor de
   Mertens ≈ 1/2, parte del normalizador) ni viene de los primos grandes
   del soporte (sus bloques son los invisibles: restan). Lo gobiernan
   los **primos del nivel** (s grandes ⇒ L pequeño); dentro de la suma,
   en el marco A la masa positiva se concentra en r=1 y los conductores
   visibles del soporte minúsculo, y con nivel real (u ≥ 3) se reparte
   y el exceso cae a 0,055·Q_c [M, K2–K4].
3. **Cota por bloques:** cerradas [T] las piezas invisible, negativa y
   de conteo (§3.1); **se resiste la masa positiva P^{(L)}**, y el
   sub-enunciado mínimo restante es **(H2‴)**: sup_L Σ_{t(r)>0} t(r)
   ≤ C·Q_completo, **equivalente a (H2″) con constantes que difieren en
   ≤ 1** [T]. Medición: máx P/Q_c = 1,248 → 0,862 (marco A, decreciente
   en z) y 0,159 (marco B) [M, K6]. (H2″) no queda demostrada.
4. **Coherencia con el puente:** sí — cualquier constante absoluta de
   (H2″)/(H2‴) da (H2) con C+1 y por tanto el G_grande = o(1) que G2
   del doc 108 necesita; el valor de C no interviene. Lo que queda
   fuera de (H2″) es la calibración numérica 0,0225 de T(N), que sigue
   siendo [M] hasta fijar la reducción del doc 94 (G4, [P-rev]) y la
   regla de atribución f.
5. **Estado:** (H2) sigue abierta, reducida ahora a su contenido mínimo
   (H2‴) — la masa positiva del exceso de Gram — con las piezas de
   signo y de conteo cerradas por separado y la localización del exceso
   medida. Goldbach no se declara demostrada.

## Apéndice. Cómo correr la verificación

```
cd ~/workspace/goldbach
python3 111_exceso_gram_euler.py | tee salida_doc111_local.txt
```

Solo stdlib (math, sys, time, importlib para reutilizar las definiciones
de los docs 109/110 y garantizar que son las mismas). Tiempo: ~2,3 s.
Salida: PASS, chequeos:

- K0: validación cruzada contra Level110 del doc 110 (Q y la forma
  cerrada (110.2), acuerdo ≤ 3,9e-16 en los niveles peores).
- K1: Σ_r t(r) = Q − Q_completo (0,00e+00 en los niveles que deciden);
  factorización t(r) = δ(2Π+δ)/φ por conductor (≤ 3,5e-18); consistencia
  con (110.2) (≤ 8,3e-16).
- K2: bloque invisible exacto (error 0,0e+00); 0 ≤ N ≤ Q_c; E = P − N;
  conteos de signos en los niveles representativos.
- K2b: clasificación de S(r) en los 240 niveles: S(r) < 0 nunca;
  0 < S < Π raro en el marco A (16/52/90) y presente en el B (858);
  S > Π: 175/384/557/2205 instancias por marco.
- K3: bloques de t(r) y de Def2 por primo mayor P⁺(r): en el marco A
  cada bloque negativo coincide con el Def2 invisible; el bloque
  P⁺ = 1 (r=1) iguala o excede el exceso neto.
- K4: bloques por ω(r): marco A, ω=0 positivo y dominante; ω ≥ 2 neto
  negativo; marco B repartido.
- K5: partición de Euler de Q_completo en √z (factor de cola
  0,577–0,618, lejos de 1); cota telescópica M_full^{>w} ≤ exp(1/(⌊w⌋−1))
  verificada; M_cola/M_V por P⁺(d): M_cola vive en P⁺(d) ≥ 7 a z=29.
- K6: libro mayor de la cota: máximo de P/Q_c sobre los 60 niveles de
  cada marco (1,248/1,073/0,862 en A; 0,159 en B); P proviene entero de
  δ > 0 (el canal δ < −2Π aporta 0,000).

---
*Goldbach sigue abierta. El exceso de Gram, abierto conductor a
conductor, no esconde ninguna pieza de conteo más: cerrados el bloque
invisible, la masa negativa y la cola de Euler, todo lo que queda de
(H2) es su masa positiva — (H2‴) — concentrada en los niveles pequeños
que fijan los primos de conmutación grandes y, en ventana, en los
conductores visibles más pequeños. Un lema fundamental para la forma de
Gram es, desde este doc, un lema sobre desviaciones positivas.*
