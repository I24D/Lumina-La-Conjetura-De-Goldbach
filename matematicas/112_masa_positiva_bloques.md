# DOC 112 — La masa positiva por bloques de conductores: r=1 cae por el lema fundamental; el resto es (H2⁗)

**Fecha:** 2026-11-28.
**Autor:** Muse (línea dedicada de Goldbach).
**Archivos:** `112_masa_positiva_bloques.md` (este doc),
`112_masa_positiva_bloques.py` (script, solo stdlib),
`salida_doc112_local.txt` (salida completa, **PASS**, chequeos K0–K5, 0,8 s).

**Convención de marcas:** [T] teorema / verificado (prueba escrita o cálculo
exacto por dos rutas); [N] no tenido; [M] medición; [P-rev] pendiente de
re-verificación contra un doc no disponible en local. MEDICIÓN ≠ PRUEBA en
todo el documento.

## §0. Alcance y método

El doc 111 redujo (H2″) a su contenido mínimo **(H2‴)**: acotar la masa
positiva del exceso de Gram,

  P^{(L)} := Σ_{r : t(r)>0} t(r),   t(r) := δ(r)·(2Π(r) + δ(r))/φ(r),

con las piezas de signo y de conteo cerradas por separado (bloque
invisible ≤ 0 exacto; 0 ≤ N^{(L)} ≤ Q_completo; cola de conteo cerrada).
Este doc ataca (H2‴) **por bloques de conductores**, con las definiciones
exactas de los docs 109/110/111 reutilizadas por import (el script
importa el del doc 111, que importa el del 110, que importa el del 109:
son literalmente las mismas funciones). Trabajo en el modo de límite
fijo del doc 109 (Q_s = Q^{(L)}, L = D_A//s):

- 𝒫 = primos permitidos: impares, < z_A, ∤ AN. a(d) := μ(d)/φ(d).
- V = V⁺(L) (Nathanson, Teorema 9.3, β = 2); r es **visible** si divide
  a algún d ∈ V, e **invisible** en caso contrario.
- S(r), Π(r), δ(r) = S(r) − Π(r), Def2, Q_completo: docs 109/110.
- g(d) := ∏_{p|d} 1/(p−2) (doc 110); se usa la identidad
  Π(r)/φ(r) = Π_full·g(r) [T, doc 110 §1.2].

**Umbral.** Fijo **R₀ := z²**, y justifico la elección: los primos de 𝒫
son < z, luego si ω(r) ≤ 2 entonces r < z² (un primo: r < z; dos primos:
r = pq < z·z), es decir

  **r ≥ R₀ ⟹ ω(r) ≥ 3** [T; verificado por enumeración en K3:
  el máximo r con ω ≤ 2 es 143/437/899 < z² = 289/841/1369].

Es el primer umbral en el que la cola del producto tipo M sobre
{r ≥ R₀} se controla por la parte ω ≥ 3 de las sumas elementales de
{1/(p−2)} (§3). También tabulo R₀ = z en K1: la frontera (ii)/(iii) es
sensible al umbral (en el marco B casi todo el residuo cae en (iii) con
R₀ = z), pero la unión (ii)+(iii) no depende de él — y es la unión lo
que sobrevive como enunciado (§4).

**Los cuatro bloques.**

  (i)   r = 1;
  (ii)  r visible, 1 < r < R₀;
  (iii) r visible, r ≥ R₀;
  (iv)  r invisible.

## §1. La partición y las masas [T + M]

**Partición [T].** Los invisibles no aportan a P: si r es invisible,
t(r) = −Π(r)²/φ(r) ≤ 0 exactamente (doc 111 §1.2). Luego, nivel a nivel,

  **P^{(L)} = t(1)⁺ + bloque(ii) + bloque(iii),**

verificado con error ≤ 5,6e-17 en los 240 niveles muestreados y con
máx t(r) = 0,00e+00 sobre los invisibles (K0). Los máximos de P/Q_c
reproducen el doc 111 (K6): 1,248 / 1,073 / 0,862 (marco A, z = 17/29/37)
y 0,159 (marco B), diferencias ≤ 0,0005.

**Masas [M, K1/K5].** En el nivel de máximo P de cada marco (R₀ = z²):

| marco | z | L | P/Q_c | (i) t(1)⁺ | (ii) | (iii) |
|---|---|---|---|---|---|---|
| A | 17 | 112 | 1,248 | 1,248 (100%) | 0,000 (0 r) | 0,000 (0 r) |
| A | 29 | 273 | 1,073 | 0,743 (69%) | 0,330 (1 r) | 0,000 (0 r) |
| A | 37 | 410 | 0,862 | 0,481 (56%) | 0,381 (3 r) | 0,000 (0 r) |
| B | 29 | 24399 | 0,159 | 0,028 (18%) | 0,129 (26 r) | 0,002 (1 r) |

Y los máximos de cada bloque sobre los 60 niveles de cada marco:
(i) = 1,248/0,743/0,481/0,028; (ii) = 0,224/0,330/0,381/0,129;
(iii) = 0,001/0,004/0,007/0,011 (en unidades de Q_c). El bloque (iii)
es minúsculo en toda la ventana, y en los niveles que deciden el
marco A es literalmente cero: no hay conductores visibles ≥ z² en el
soporte minúsculo. El bloque (i) es el término dominante del marco A;
el (ii) manda en el marco B. Sobre este reparto actúan §§2–4.

## §2. Bloque (i): CERRADO [T] — el lema fundamental clásico sí aplica a r=1

### §2.1. La identificación: S(1) es una suma de pesos λ⁺

Por el Teorema 9.3 de Nathanson en la forma en que esta línea lo usa
(docs 95/97/109), V⁺(L) es **el soporte de los pesos superiores λ⁺ del
β-sieve (β = 2) de nivel L**, y sobre ese soporte **λ⁺(d) = μ(d)**.
Por tanto, con g(d) := 1/φ(d) (multiplicativa, g(p) = 1/(p−1)):

  **S(1) = Σ_{d∈V⁺(L)} μ(d)/φ(d) = Σ_{d | P_𝒫(z)} λ⁺(d)·g(d),** [T]

exactamente: la obstrucción del doc 109 §3.5 (los soportes desplazados
V^{(r)} no son soportes de Rosser–Iwaniec de parámetros verificados)
**no existe para r = 1**, porque V^{(1)} = V es el soporte sin desplazar.
Y el término principal de la suma es el correcto:
∏_{p∈𝒫}(1 − g(p)) = ∏_{p∈𝒫}(1 − 1/(p−1)) = Π_full = Π(1), luego
δ(1) = S(1) − Π_full es exactamente el error del lema fundamental
para esta suma de pesos.

### §2.2. El resultado clásico y la verificación de sus hipótesis

**Lema fundamental (forma usada).** Sea 𝒫 un conjunto de primos,
g multiplicativa con 0 ≤ g(p) < 1 verificando la condición de
dimensión (Ω_κ) con κ = 1: existe K absoluto tal que para
2 ≤ η ≤ ξ, ∏_{η≤p<ξ, p∈𝒫}(1 − g(p))^{−1} ≤ (ln ξ/ln η)·(1 + K/ln η).
Sean λ^± los pesos del β-sieve (β = 2) de nivel D = z^u. Entonces

  Σ_{d | P_𝒫(z)} λ^±(d) g(d)
      = ∏_{p∈𝒫, p<z} (1 − g(p)) · {1 + ρ(u)},
  con ρ(u) → 0 cuando u → ∞

(en la formulación estándar de Tenenbaum del lema fundamental de la
criba combinatoria, ρ(u) = O(u^{−u/2}); Halberstam–Richert, *Sieve
Methods*, Cap. 3, págs. 82–93, donde se formula el lema fundamental
del β-sieve). En la aplicación clásica a una sucesión cribada hay
además un término de restos O(Σ_{d<D}|R_d|); **aquí no hay sucesión**:
la suma de pesos es el objeto mismo y las densidades son exactas,
luego ese término no aparece.

**Verificación de hipótesis para g(p) = 1/(p−1)** [T]:

1. *Rango:* 0 ≤ g(p) = 1/(p−1) ≤ 1/2 < 1 para p ≥ 3 (en 𝒫, p ≥ 5 si
   3 | A). ✓
2. *Dimensión κ = 1:* para η ≤ p < ξ,
   −ln(1 − 1/(p−1)) = ln(1 + 1/(p−2)) ≤ 1/(p−2)
   = 1/p + 1/(p(p−2)).
   Sumando y usando Mertens clásico,
   Σ_{η≤p<ξ} 1/p ≤ ln ln ξ − ln ln η + C/η y
   Σ_{η≤p<ξ} 1/(p(p−2)) ≤ C′/η, luego
   ∏_{η≤p<ξ, p∈𝒫}(1 − g(p))^{−1}
     ≤ exp(ln ln ξ − ln ln η + O(1/η))
     = (ln ξ/ln η)·(1 + K/ln η)
   con K absoluto. Retirar primos de 𝒫 (el 2, los que dividen a AN)
   solo reduce el producto. ✓
3. *Nivel:* D = L = z^u con u = ln L/ln z. El lema da ρ(u) → 0; en
   particular existe un **u₀\* absoluto** (depende solo de κ = 1,
   β = 2 y de la constante implícita clásica) tal que u ≥ u₀\*
   implica |ρ(u)| ≤ 1/2. ✓ — con la precisión de §2.4.

### §2.3. La cota del bloque

Con ρ̃ := δ(1)/Π_full (con signo), la factorización del doc 111 da la
identidad exacta [T, verificada en K2 a 1,1e-16]:

  t(1) = ρ̃·(2 + ρ̃)·Π_full².

Y **Π_full² ≤ Q_completo** [T]: Q_completo = Π_full·M_full (doc 110)
con Π_full ≤ 1 (cada factor 1 − 1/(p−1) < 1) y M_full ≥ 1 (cada factor
1 + 1/((p−1)(p−2)) > 1); verificado en K2: Π_full²/Q_c =
0,452/0,379/0,352 a z = 17/29/37. Por tanto, para todo nivel con
u ≥ u₀\*:

  **t(1)⁺ ≤ ρ(u)·(2 + ρ(u))·Q_completo ≤ (5/4)·Q_completo.** [T]

**Bloque (i) CERRADO [T]**, con constante 5/4 (con ρ ≤ 1/2; con la
ρ(u) = O(u^{−u/2}) citada, la constante tiende a 0 con u).

### §2.4. Alcance exacto del cierre, y la medición

El cierre vale bajo la condición de nivel u ≥ u₀\* — la misma forma
de hipótesis que (H2)/(H2″)/(H2‴) ya llevan («existe u₀ absoluto»):
basta tomar u₀ := máx(u₀ del marco, u₀\*), que sigue siendo absoluto.
No afirmo que el u_min = 5/3 del marco A supere u₀\*: el marco A es el
marco de estrés de la línea (doc 109: el sup vive en la meseta de
niveles pequeños), no el marco del puente (doc 108 §1.2: D_A ≍ N^{1/2},
y_A = N^{1/3}). La medición [M, K2] es exactamente la que la FL
predice: |ρ₁(u)| = |δ(1)|/Π_full por bandas de u sobre los 240
niveles — u < 2: ≤ 0,939; 2 ≤ u < 3: ≤ 0,455; 3 ≤ u < 4: ≤ 0,037;
u ≥ 4: ≤ 0,002 — decae con u y es ≤ 0,037 en todo el marco B
(u ≥ 3). Compatible con ρ(u) → 0; no es la prueba.

## §3. Bloque (iii): la ruta absoluta NO cierra [T la pérdida]

### §3.1. La reducción exacta y los dos factores

De t(r) ≤ (2Π(r) + |δ(r)|)·|δ(r)|/φ(r) [T] y de
Π(r)/φ(r) = Π_full·g(r) [T, doc 110 §1.2; verificado en K3 a 2,8e-17]:

  bloque(iii) ≤ Λ₃ + Θ₃,
  Λ₃ := 2·Π_full·Σ_{r≥R₀} g(r)·|δ(r)|,
  Θ₃ := Σ_{r≥R₀} δ(r)²/φ(r).

La ruta prescrita acota |δ(r)| por la cola absoluta. En efecto
δ(r) = −T(r) (doc 110 §1.1), luego [T, verificado con 𝒯^{(r)}
explícito en K3, error 0,00e+00]:

  |δ(r)| ≤ D⁻(r) := Σ_{d′∈𝒯^{(r)}} 1/φ(d′) ≤ Σ_{d′} 1/φ(d′)
         = ∏_{p∈𝒫}(1 + 1/(p−1)) =: P⁺.

(La variante por el soporte, |δ(r)| ≤ D⁺(r) = Σ_{d′∈V^{(r)}} 1/φ(d′)
≤ Ψ(L/r), no decae en r y es peor; no se usa.) Y la cola de g se
controla por conteo [T]:

  Σ_{r≥z²} g(r) ≤ Σ_{ω(r)≥3} g(r) = Ĝ − 1 − G₁ − e₂,

con Ĝ := ∏_{p∈𝒫}(1 + 1/(p−2)), G₁ := Σ_{p∈𝒫} 1/(p−2),
e₂ := Σ_{p<q} 1/((p−2)(q−2)) (identidad exacta, verificada en K3 por
enumeración y por sumas de potencias: 0,0195/0,0645/0,0938 a
z = 17/29/37). El mayorante resultante para la parte lineal es

  Λ₃/Q_completo ≤ 2·P⁺·Σ_{r≥z²} g(r)/M_full
     = **0,060 / 0,128 / 0,177** (z = 17/29/37). [cálculo exacto, K3]

### §3.2. Por qué no es una constante absoluta [T]

El mayorante es engañosamente bueno en ventana y **no** es O(1):

- P⁺ = ∏_{p∈𝒫}(1 + 1/(p−1)) ≍ e^γ·(ln z)·(factor local de 𝒫) → ∞
  (Mertens clásico) [T].
- Σ_{ω≥3} g(r) ≥ e₃({1/(p−2)}) y e₃ → ∞: en efecto,
  e₃ = (G₁³ − 3G₁G₂ + 2G₃)/6 ≥ (G₁³ − 3G₁G₂)/6 con
  G₂ = Σ_{p∈𝒫} 1/(p−2)² ≤ Σ_{m≥3} 1/(m−2)² = π²/6 − 1 < 0,65
  absoluto, mientras G₁ = Σ_{p∈𝒫} 1/(p−2) ~ ln ln z + O(1) → ∞
  (Mertens clásico) [T].

Luego el factor 2·P⁺·Σ_{ω≥3}g(r)/M_full crece al menos como
(ln z)·(ln ln z)³/(6·M_full) — en la práctica, calculado con
Σ_{ω≥3} g(r) exacto: **0,615 (z=53), 0,868 (z=71), 1,234 (z=101)**,
creciente (K3). La ruta absoluta para Λ₃ **pierde
(ln z)(log log z)³**: no acota el bloque (iii) por C·Q_completo con
C absoluta. Es la misma muerte de todos los mayorantes absolutos de
la línea (triángulo, doc 109 §3.2; diagonal de la cola, doc 110
§3.1), ahora con la pérdida medida factor a factor: un ln z viene de
P⁺ (la cola absoluta de δ) y el (log log z)³ de la cola ω ≥ 3 del
producto g.

La parte cuadrática no rescata ni condena: Θ₃ ≤
(máx_{r≥R₀, visible} |δ(r)|/Π(r))·Λ₃/2 [T], y ese factor puntual no
está acotado por ninguna constante probada (la forma puntual del
lema fundamental es falsa en general, doc 109 §3.4; sobre los
visibles r ≥ R₀ es [M]: 0 en los niveles peores del marco A — no hay
tales conductores — y 0,452 en el marco B, K3). La ruta muere ya en
Λ₃.

### §3.3. Estado del bloque (iii)

**NO cerrado** por la ruta absoluta con productos tipo M. Su masa
real, en cambio, es despreciable en ventana [M, §1]: máximos
0,001/0,004/0,007 (marco A) y 0,011 (marco B) de Q_c sobre todos los
niveles, y 0 en los niveles que deciden el marco A. El bloque (iii)
es un problema de la prueba, no de la masa: lo que lo cerraría es una
cota |δ(r)| ≤ ρ·Π(r) **en promedio contra el peso g(r)** sobre los
visibles r ≥ R₀, que no está disponible (y cuya forma puntual es
falsa ya por los invisibles, que no pertenecen al bloque pero fijan
el fenómeno).

## §4. Bloque (ii): ni FL ni análisis finito — el residuo (H2⁗)

### §4.1. Por qué no cierra

**No es una familia finita uniformemente** [M, K4]. El censo de
conductores visibles con 1 < r < z² crece con z: a nivel máximo hay
10 (z=17), 34 (z=29), 60 (z=37) conductores — π_𝒫 primos,
C(π_𝒫, 2) semiprimos y unos pocos ω ≥ 3 pequeños (385, 5005, … <
z²) — y la familia contiene **todos** los primos de 𝒫. Un análisis
finito uniforme necesitaría una cota por conductor cuya suma sobre
~π(z)²/2 conductores fuera O(Q_completo); la cota individual
disponible es t(r) ≤ S(r)²/φ(r) [T], cuyo mayorante absoluto
(P⁺)²/φ(r) pierde (ln z)² (docs 109 §3.2, 110 §3.1). Y el máximo
t(r) individual medido en el bloque no es pequeño: 0,330·Q_c
(r = 5, z=29) y 0,214·Q_c (r = 5, z=37) en los niveles que deciden.

**La FL de §2 no aplica.** Para r = p primo, S(p) es la suma
μ(d′)/φ(d′) sobre el soporte **desplazado** V^{(p)} = {d′ : (d′,p)=1,
p·d′ ∈ V}: la condición de Nathanson se impone sobre el producto
p·d′, no sobre d′ — V^{(p)} no es un soporte de Rosser–Iwaniec de
parámetros verificados, y la suma no lleva pesos λ^± (doc 109 §3.5,
razones 1–2). Lo mismo para semiprimos y demás r > 1. El bloque (ii)
es exactamente el territorio de la obstrucción ya localizada.

### §4.2. El sub-enunciado residual: (H2⁗)

Como el bloque (iii) tampoco cierra por la ruta absoluta (§3), el
residuo honesto es la unión (ii)+(iii) — que no depende del umbral
R₀. Por la partición de §1 y el cierre de §2, es equivalente a
(H2‴):

> **(H2⁗) — masa positiva sobre conductores visibles r > 1.**
> Existen constantes absolutas C y u₀ tales que para todo marco
> (A impar fijo, N par, (A,N)=1, z_A < y_A ≤ D_A) con
> ln(D_A/y_A)/ln z_A ≥ u₀:
>
>   sup_L R^{(L)} := sup_L Σ_{r visible, r > 1, t(r)>0} t(r)
>                  ≤ C · Q_completo(z_A; N),
>
> donde t(r) = φ(r)·[ C^{(L)}(r)² − C^full(r)² ], r es visible si
> divide a algún d ∈ V⁺(L), y L ∈ [D_A/y_A, D_A/z_A] en el modo fijo
> (y su análogo de límite variable [P-rev], como en los docs
> 109–111).

**Equivalencia [T].** P^{(L)} = t(1)⁺ + R^{(L)} nivel a nivel (§1:
los invisibles no aportan y r=1 es el único conductor restante).
(H2⁗) con constante C ⟹ (H2‴) con C + 5/4 (por §2.3, para
u ≥ u₀\*). (H2‴) con constante C ⟹ (H2⁗) con C (pues
R^{(L)} ≤ P^{(L)}). Es decir, (H2⁗) ⟺ (H2‴) módulo la pieza ya
demostrada: no es un enunciado más débil disfrazado, es (H2‴) sin
su término ya cerrado.

**Tamaño medido [M, K4].** máx_L R^{(L)}/Q_c = **0,224 / 0,330 /
0,381** (marco A, z = 17/29/37) y **0,130** (marco B), contra
máx P/Q_c = 1,248/1,073/0,862/0,159. Cada término del residuo tiene
además φ(r) ≥ 4 y peso de comparación acotado:
Σ_{r>1} Π(r)²/φ(r) = Q_completo − Π_full² [T, verificado en K4].

## §5. Veredicto global

1. **(H2‴) NO queda demostrada por bloques.** El inventario exacto
   de lo que falta es un único sub-enunciado: **(H2⁗)**, la masa
   positiva sobre conductores visibles r > 1 (§4.2), equivalente a
   (H2‴) módulo las piezas cerradas.
2. **Cerrados [T]:** el bloque (iv) — invisibles, t(r) =
   −Π(r)²/φ(r) ≤ 0 exacto (doc 111) — y el bloque **(i) r=1**
   (este doc, §2): S(1) es la suma de pesos λ⁺ sin desplazar, el
   lema fundamental clásico (Halberstam–Richert; formulación de
   Tenenbaum, error O(u^{−u/2})) se le aplica con las hipótesis
   verificadas una a una (dimensión κ=1 para g(p) = 1/(p−1) vía
   Mertens; nivel u), y da t(1)⁺ ≤ ρ(2+ρ)·Q_completo con
   ρ(u) → 0 — en particular ≤ (5/4)·Q_completo para u ≥ u₀\*
   absoluto.
3. **No cerrados:** el bloque (iii) por la ruta absoluta prescrita
   — el mayorante con productos tipo M pierde
   (ln z)(log log z)³ [T, §3.2: P⁺ ≍ e^γ ln z y
   Σ_{ω≥3} g(r) ≥ e₃ → ∞; valores 0,060 → 1,234 entre z=17 y
   z=101] — pese a que su masa real es ≤ 0,011·Q_c en ventana [M];
   y el bloque (ii), cuya familia crece con z (contiene todos los
   primos de 𝒫) y cuyos conductores viven sobre soportes
   desplazados, fuera del alcance de la FL clásica.
4. **¿Es lo que falta más pequeño que (H2‴)?** En **rango de r**,
   sí, y sustancialmente: fuera r=1 — el 100%/69%/56% de P en los
   niveles que deciden el marco A, y el único conductor cuyo
   soporte no está desplazado — y fuera los invisibles. En **masa
   medida**, sí: el residuo es ≤ 0,381·Q_c frente a P ≤ 1,248·Q_c.
   En **niveles y uniformidad de marco**, no: el sup_L y las
   hipótesis son idénticos. En **tipo**, el residuo es homogéneo:
   toda su masa vive sobre soportes genuinamente desplazados
   V^{(r)}, r ≥ 5 — es decir, (H2⁗) es exactamente *el lema
   fundamental sobre soportes desplazados, en su parte positiva y
   sin el conductor trivial*, que es donde el doc 109 §3.5 ya había
   localizado la obstrucción. La descomposición por bloques no ha
   hecho el problema más fácil en su núcleo; lo ha dejado desnudo:
   cerradas la pieza trivial (r=1, por la FL clásica) y las piezas
   de signo y conteo (docs 110–111), todo (H2) es ya (H2⁗).
5. **Estado:** (H2) sigue abierta, reducida a (H2⁗). Goldbach no
   se declara demostrada.

## Apéndice. Cómo correr la verificación

```
cd ~/workspace/goldbach
python3 112_masa_positiva_bloques.py | tee salida_doc112_local.txt
```

Solo stdlib (math, sys, time, importlib para reutilizar las
definiciones de los docs 109/110/111 y garantizar que son las
mismas). Tiempo: ~0,8 s. Salida: PASS, chequeos:

- K0: validación cruzada contra el doc 111 (máx P/Q_c =
  1,248/1,073/0,862/0,159, dif. ≤ 0,0005); partición
  P = (i)+(ii)+(iii) con error ≤ 5,6e-17 en los 240 niveles;
  máx t(r) sobre invisibles = 0,00e+00.
- K1: masas de bloque por marco en el nivel de máximo P y máximos
  por bloque sobre los 60 niveles (R₀ = z²; comparación con R₀ = z).
- K2: Π_full ≤ 1 ≤ M_full y Π_full² ≤ Q_c (cocientes
  0,452/0,379/0,352); identidad t(1) = ρ(2+ρ)Π_full² (error
  1,1e-16); |ρ₁(u)| por bandas de u: 0,939 / 0,455 / 0,037 / 0,002.
- K3: r ≥ z² ⟹ ω ≥ 3 (máx r con ω ≤ 2: 143/437/899 < z²);
  Π(r)/φ(r) = Π_full·g(r) (error ≤ 2,8e-17);
  Σ_{r≥z²} g(r) ≤ Σ_{ω≥3} g(r) = Ĝ − 1 − G₁ − e₂ (exacto);
  mayorante de Λ₃: 0,060/0,128/0,177 (z ≤ 37) y
  0,615/0,868/1,234 (z = 53/71/101, solo primos);
  |δ(r)| ≤ D⁻(r) ≤ P⁺ con 𝒯^{(r)} explícito (error 0,00e+00);
  máx |δ|/Π sobre visibles r ≥ z²: 0 (marco A) y 0,452 (marco B) [M].
- K4: censo de visibles 1<r<z² por ω a nivel máximo
  (10/34/60 conductores); máx t(r) en (ii): 0,330/0,214·Q_c en los
  niveles que deciden; máx_L R^{(L)}/Q_c =
  0,224/0,330/0,381/0,130; Σ_{r>1} Π(r)²/φ(r) = Q_c − Π_full².
- K5: libro mayor de la cuota de r=1 en el máximo de P:
  100%/69%/56%/18%.

---
*Goldbach sigue abierta. La masa positiva del exceso de Gram,
partida por bloques de conductores, deja exactamente una pieza
demostrada nueva — r=1, por el lema fundamental clásico aplicado
al único soporte no desplazado — y un residuo homogéneo, (H2⁗):
las desviaciones positivas sobre los soportes desplazados
V^{(r)}, r > 1. Ni la ruta absoluta de los productos tipo M
(pierde (ln z)(log log z)³) ni un análisis finito (la familia
crece con z) lo tocan. (H2) es, desde este doc, (H2⁗).*
