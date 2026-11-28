# DOC 110 — El segundo momento de los defectos, expandido: (H2′) es (H2) disfrazado

**Fecha:** 2026-11-28.
**Autor:** Muse (línea dedicada de Goldbach).
**Archivos:** `110_segundo_momento_defectos.md` (este doc),
`110_segundo_momento_defectos.py` (script, solo stdlib),
`salida_doc110_local.txt` (salida completa, **PASS**, chequeos K0–K7, 0,9 s).

**Convención de marcas:** [T] teorema / verificado (prueba escrita o cálculo
exacto por dos rutas); [N] no tenido; [M] medición; [P-rev] pendiente de
re-verificación contra un doc no disponible en local. MEDICIÓN ≠ PRUEBA en
todo el documento.

## §0. Alcance y método

El doc 109 redujo (H2) — la cota de Gram normalizada
sup_{z_A≤s<y_A} Q_s ≪ Q_completo(z_A;N) — al sub-enunciado (H2′): la cota
L² de los defectos δ_s(r) = S_s(r) − Π(r),

  Σ_r δ_s(r)²/φ(r) ≤ C·Q_completo(z_A;N),

y delimitó el paso que faltaba como «un lema fundamental uniforme sobre
los soportes desplazados V_s^{(r)} en promedio L²». Este doc ataca (H2′)
por la vía que el doc 109 no intentó: **expandir el segundo momento
exactamente** y reindexarlo, al estilo de las reindexaciones A100.1
(doc 101) y (106.2) (doc 106). El resultado es una identidad cerrada que
lo cambia todo sobre (H2′): el momento de los defectos es, término a
término, el exceso de la forma de Gram sobre su producto completo más
un término de conteo puro y explícito. En consecuencia (H2′) no es un
sub-enunciado más débil que (H2): **es (H2) disfrazado**, y la vía L² de
los defectos queda agotada — no por un fracaso numérico, sino por una
identidad.

Trabajo en el **modo de límite fijo** del doc 109 (cada d tiene un único
primo de conmutación s; Q_s = Q^{(L)} con L = D_A//s), con las
definiciones exactas de los docs 96/97/109:

- 𝒫 = primos permitidos: impares, < z_A, ∤ AN. Todo d, r es libre de
  cuadrados sobre 𝒫. a(d) := μ(d)/φ(d).
- V = V⁺(L) (Nathanson, Teorema 9.3, β = 2): d ∈ V si d < L y para todo
  m impar de la factorización decreciente, p_m²·(p_1···p_m) < L.
- C(r) = Σ_{d∈V, r|d} a(d); Q^{(L)} = Σ_r φ(r)C(r)² [T, doc 97 §2.1].
- S(r) := φ(r)C(r)/μ(r) = Σ_{d′∈V^{(r)}} a(d′), con el **soporte
  desplazado** V^{(r)} = {d′ : (d′,r)=1, r·d′ ∈ V} [T, doc 109 §3.1].
- Π(r) := ∏_{p∈𝒫∖primos(r)}(1 − 1/(p−1)) = Σ_{(d′,r)=1} a(d′) (suma
  completa, sin criba) [T, doc 97 §2.3].
- δ(r) = S(r) − Π(r); Def2 := Σ_r δ(r)²/φ(r) — el momento de (H2′).
- Q_completo = ∏_{p∈𝒫}(1 − 1/(p−1) + 1/(p−1)²) [T, doc 97 §2.3].

## §1. La expansión exacta y su reindexación [T]

### §1.1. La suma doble sobre las colas de los desplazados

Definiendo la **cola del desplazado**
𝒯^{(r)} = {d′ : (d′,r)=1, r·d′ ∉ V}, la partición del conjunto de los
d′ coprimos con r en V^{(r)} ⊔ 𝒯^{(r)} da S(r) + T(r) = Π(r) con
T(r) := Σ_{d′∈𝒯^{(r)}} a(d′), luego δ(r) = −T(r) y la expansión literal
del momento es la suma doble sobre los soportes desplazados:

  **Def2 = Σ_r (1/φ(r)) · Σ_{d′,e′ ∈ 𝒯^{(r)}} a(d′)a(e′).**

Verificación [T, K1/K4]: con V^{(r)} y 𝒯^{(r)} construidos desde cero
(sin pasar por C), la partición S(r) + T(r) = Π(r) vale para todo r con
error ≤ 1,11e-15, y la suma doble reproduce Def2 con error absoluto
≤ 5,0e-16 en los niveles examinados. K0 valida Q y Def2 contra el
LevelData del doc 109: acuerdo 0,00e+00.

### §1.2. La reindexación del término cruzado

Expandiendo el cuadrado por δ = S − Π hay tres sumas:

  Def2 = Σ_r S(r)²/φ(r) − 2·Σ_r S(r)Π(r)/φ(r) + Σ_r Π(r)²/φ(r).

La primera es Q^{(L)} (pues S(r)²/φ(r) = φ(r)C(r)²) y la tercera es
Q_completo [T, doc 97 §2.3; verificado en K2b a ≤ 3,8e-16]. La suma
cruzada es la que se reindexa. Como S(r)/φ(r) = C(r)/μ(r) = μ(r)C(r):

  Σ_r S(r)Π(r)/φ(r) = Σ_r μ(r)Π(r)C(r)
                   = Σ_{d∈V} a(d) · Σ_{r|d} μ(r)Π(r)   (intercambio).

La suma interna colapsa. Escribiendo Π_full := ∏_{p∈𝒫}(1 − 1/(p−1)),
se tiene Π(r) = Π_full·∏_{p|r}(p−1)/(p−2), y con y_p := (p−1)/(p−2):

  Σ_{r|d} μ(r)Π(r) = Π_full·∏_{p|d}(1 − y_p)
                  = Π_full·∏_{p|d}(−1/(p−2))
                  = Π_full·μ(d)·g(d),

donde **g(d) := ∏_{p|d} 1/(p−2)** (g(1) = 1). Como a(d)·μ(d) = g(d)/φ(d)
·… precisamente a(d)·μ(d)·g(d) = g(d)/φ(d) (μ² = 1), la suma cruzada es
Π_full·M_V con

  **M_V := Σ_{d∈V} g(d)/φ(d)** — una suma de conteo sobre el soporte.

Resultado, la **primera forma cerrada**:

  **(110.1)  Def2 = Q^{(L)} + Q_completo − 2·Π_full·M_V.** [T]

### §1.3. La segunda forma cerrada: el término de la cola de M

La función g(d)/φ(d) = ∏_{p|d} 1/((p−1)(p−2)) es multiplicativa con
producto de Euler convergente:

  M_full := Σ_{d : primos ⊆ 𝒫} g(d)/φ(d)
          = ∏_{p∈𝒫}(1 + 1/((p−1)(p−2))),

y un cálculo local exacto —
(1 − 1/(p−1))·(1 + 1/((p−1)(p−2))) = (p²−3p+3)/(p−1)²
= 1 − 1/(p−1) + 1/(p−1)² — da la identidad

  **Q_completo = Π_full · M_full.** [T]

Definiendo la **cola de conteo** M_cola := Σ_{d∉V} g(d)/φ(d)
= M_full − M_V ≥ 0, la sustitución en (110.1) da la **segunda forma
cerrada**, la identidad central de este doc:

  **(110.2)  Def2 = ( Q^{(L)} − Q_completo ) + 2·Π_full·M_cola.** [T]

Verificación [T, K2/K3]: (110.1) y (110.2) reproducen Def2 en **240
niveles** (marcos A con z = 17, 29, 37 y marco B con z = 29, 60 niveles
muestreados por marco) con error ≤ 1,3e-15 relativo a la escala Q + Q_c;
el intercambio Σ_r μ(r)Π(r)C(r) = Π_full·M_V y la partición de M
(M_cola directo = M_full − M_V, chequeo absoluto por la cancelación
catastrófica de la diferencia) verificados por separado. Los máximos de
Def2/Q_c reproducen los del doc 109 (K5): 0,947 / 0,738 / 0,553
(z = 17/29/37, marco A) — validación cruzada ✓.

**Lectura.** El segundo momento de los defectos no es un objeto nuevo:
es el exceso de la forma de Gram cribada sobre su valor sin criba,
Q^{(L)} − Q_completo, más el término 2·Π_full·M_cola, que es **conteo
puro**: una cola de un producto de Euler convergente, sin ψ, sin
caracteres, sin el nivel más allá de la pertenencia d ∈ V.

## §2. Descomposición: diagonal contra cruzada, y las dos piezas de (110.2)

### §2.1. En la expansión de la cola, (H2′) vive de la cancelación cruzada [M]

Separando la suma doble de §1.1 en su diagonal d′=e′ y su cruzada
d′≠e′ (como el doc 103 separó Diag^C), en los niveles peores de cada
marco (los de menor L, que son los que deciden el sup):

| marco | L | Def2 | D_cola | X_cola | D_cola/Q_c | X_cola/Q_c |
|---|---|---|---|---|---|---|
| A, z=17 | 112 | 0,5568 | 0,8973 | −0,3405 | 1,53 | −0,58 |
| A, z=29 | 273 | 0,3701 | 0,9312 | −0,5611 | 1,86 | −1,12 |
| A, z=37 | 410 | 0,2591 | 0,8391 | −0,5800 | 1,79 | −1,24 |
| B, z=29 | 24399 | 0,0274 | 0,0317 | −0,0043 | 0,06 | −0,01 |

- En el marco A, la diagonal de la cola **sola ya excede el objetivo**:
  D_cola = 1,53–1,86 × Q_completo, más que el propio Def2
  (0,55–0,95 × Q_c). La cruzada es negativa y cancela el 38–69% de la
  diagonal. Sin cancelación cruzada no hay (H2′): la cota vive de ella.
- En el marco B (nivel alto) todo es pequeño y la diagonal casi basta
  (0,06 × Q_c): el fenómeno es de niveles bajos, coherente con el
  R(L) → 1 del doc 109 §2.3. En los niveles máximos muestreados las
  colas se vacían (a z=17, L=83521: V = enumeración completa, δ ≡ 0,
  Def2 = 0 exacto).

### §2.2. Tras reindexar: las dos piezas de (110.2) [M]

En los mismos niveles, la descomposición de (110.2):

| marco | L | (Q−Q_c)/Q_c | 2Π_full·M_cola/Q_c | Def2/Q_c |
|---|---|---|---|---|
| A, z=17 | 112 | +0,701 | 0,246 | 0,947 |
| A, z=29 | 273 | +0,619 | 0,118 | 0,738 |
| A, z=37 | 410 | +0,492 | 0,060 | 0,553 |
| B, z=29 | 24399 | +0,055 | 0,000 | 0,055 |

El grueso del momento en los niveles peores es el **exceso de Gram**
Q^{(L)} − Q_completo (positivo, 0,49–0,70 × Q_c); el término de conteo
es menor y decrece con z (0,246 → 0,060). Es decir: ni la diagonal de
la cola basta (§2.1), ni tras la reindexación aparece ninguna pieza
diagonal suficiente — el contenido no contable del momento es, entero,
el exceso de la forma de Gram sobre su producto de Euler.

## §3. Hasta dónde llega el conteo: equivalencia (H2) ⟺ (H2′), y (H2″)

### §3.1. Lo que el conteo cierra por completo [T]

M_cola es una subsuma de términos positivos de M_full, luego
M_cola ≤ M_full, y M_full está acotado por un producto convergente
explícito: a A=3, M_full(z) = 1,1405 / 1,1514 / 1,1543
(z = 17/29/37) y M_full(∞) ≤ ∏_{3<p≤10⁵}(1 + 1/((p−1)(p−2)))
= 1,1618. Por (110.2) y Q_completo = Π_full·M_full:

  **0 ≤ 2·Π_full·M_cola ≤ 2·Q_completo,** [T]

sin condición de nivel, sin lema fundamental, sin nada: conteo puro.
(Medido en los 240 niveles: el máximo de 2Π_full·M_cola/Q_c es
0,246/0,118/0,060 en el marco A y < 5e-4 en el marco B [M, K5/K6].)
En cambio el conteo NO acota la diagonal de la cola al orden correcto:
la cota pura D_cola ≤ P⁺·∏(1 + 1/(p−1)²), con
P⁺ = ∏_{p∈𝒫}(1 + 1/(p−1)) ≍ e^γ·ln z·(factor local), da
1,930 / 2,284 / 2,450 = 3,3 / 4,6 / 5,2 × Q_c a z = 17/29/37 [M, K5b]:
como múltiplo de Q_c **crece como (ln z)²** — la misma pérdida que el
doc 109 encontró en la ruta del triángulo, ahora en la diagonal de la
cola. El conteo no separa D_cola de X_cola; solo la identidad (110.2)
las disuelve a las dos en una sola pieza contable y el exceso.

### §3.2. La equivalencia [T]

De (110.2), nivel a nivel: Def2 = (Q^{(L)} − Q_completo) + τ_L con
0 ≤ τ_L ≤ 2·Q_completo. Por tanto, si C₀ es la mejor constante de (H2)
y C′ la de (H2′) (bajo la misma condición de nivel):

- Def2 ≤ C′·Q_c ⇒ Q^{(L)} = Def2 + Q_c − τ_L ≤ (C′+1)·Q_c;
- Q^{(L)} ≤ C₀·Q_c ⇒ Def2 ≤ (C₀−1)·Q_c + 2·Q_c = (C₀+1)·Q_c;

es decir **C′ ≤ C₀ + 1 y C₀ ≤ C′ + 1: (H2) ⟺ (H2′) como enunciados
O(·) con constantes absolutas.** Esto afina la reducción del doc 109
(que perdía un factor 2 por (a+b)² ≤ 2a² + 2b² y solo daba
(H2) ⇐ (H2′)): la implicación contraria también vale, y con constante
aditiva 1, no multiplicativa. La moraleja es incómoda pero exacta:
**ningún progreso sobre (H2′) puede ser más fácil que (H2)**, porque
son el mismo enunciado salvo un término de conteo acotado. En
particular, el «lema fundamental en promedio L² sobre los soportes
desplazados» que el doc 109 §3.5 describió como el paso que falta no es
un paso intermedio hacia (H2): es (H2) con otro nombre.

### §3.3. El término que se resiste y el sub-enunciado mínimo (H2″)

El término exacto que el conteo no toca es **Q^{(L)} − Q_completo**:
el exceso de la forma de Gram cribada sobre su producto de Euler.
Medido [M]: positivo y del orden del objetivo en los niveles peores
del marco A (+0,70/+0,62/+0,49 × Q_c a z = 17/29/37), pequeño en el
marco B (+0,055 × Q_c); no es un resto que decaiga en la ventana — es
Θ(Q_completo) donde el sup se decide. Ninguna ganancia asintótica hace
falta (no se pide o(Q_c), y el doc 109 ya mostró que no la hay); lo que
falta es solo la cota uniforme. El sub-enunciado mínimo restante es por
tanto:

> **(H2″) — exceso de Gram.** Existen constantes absolutas C y u₀
> tales que para todo marco (A impar fijo, N par, (A,N)=1,
> z_A < y_A ≤ D_A) con ln(D_A/y_A)/ln z_A ≥ u₀:
>
>   sup_L [ Q^{(L)} − Q_completo(z_A; N) ] ≤ C · Q_completo(z_A; N),
>
> donde Q^{(L)} = Σ_r φ(r)C^{(L)}(r)² es la forma de Gram al nivel
> L ∈ [D_A/y_A, D_A/z_A] del modo fijo (y su análogo de límite
> variable con la cola s′>s del doc 95 [P-rev], como en el doc 109).

Qué suma: la diferencia término a término
Σ_r [φ(r)C(r)² − Π(r)²/φ(r)] entre la forma cribada y el producto
completo. Qué rango: los niveles del marco con la condición de nivel
de (H2). Qué ganancia: ninguna asintótica — una cota O(Q_completo)
uniforme sobre un exceso que ya es Θ(Q_completo) positivo en medición.
Por §3.2, (H2″) ⟺ (H2) ⟺ (H2′): (H2″) no es un enunciado más débil,
es el **contenido mínimo** de los tres, despojado de las dos
reformulaciones (defectos L², colas de desplazados) que este doc
demuestra equivalentes. Lo que cerraría (H2″) es un lema fundamental
**para la forma de Gram misma** — control del exceso de
Σ_r φ(r)C(r)² sobre su producto de Euler bajo la restricción V⁺ — no
un lema sobre los defectos. La vía de los defectos está agotada.

**Estado de la prueba:** (H2′) no queda demostrada — queda disuelta.
(H2) sigue abierta, ahora con la certeza de que su formulación de
exceso (H2″) es la mínima y de que no existe atajo por el segundo
momento de los defectos.

## §4. Relación con Δ_r: objetos independientes

La pregunta era si δ_s(r) y la corrección imprimitiva Δ_r de los
docs 97–100 (S_r = Σ_{d|r} μ(r/d)φ(d)E(N;d,N) + Δ_r) están conectados
por una identidad o una cota, dado que ambos «viven sobre soportes
desplazados». Respuesta: **no hay identidad ni cota en ninguna
dirección; son objetos independientes**, y la expansión de §1 es la
evidencia principal.

1. **La forma cerrada no contiene a Δ_r.** (110.2) expresa Σ_r δ²/φ
   usando solo Q^{(L)}, Q_completo, Π_full y M_cola — objetos del lado
   de los coeficientes, construidos con μ, φ, V⁺ y (A, z_A, L, primos
   de AN). En ningún paso de la derivación (partición, intercambio,
   Euler) aparece E(N;d,N), ni un carácter, ni una potencia de primo
   p^j ≤ N. Recíprocamente, la identidad de Möbius corregida
   (docs 97 §3.3, 100 §2) no contiene V, ni s, ni Π(r): Δ_r
   = Σ_{p∥r} Σ*_{χ* mod r/p} χ̄*(N)·Σ_{p^j≤N} χ*(p)^j log p
   (doc 97, Teorema 97.1) vive en los datos (r, N, caracteres).
2. **Los conjuntos de índices difieren.** δ_s(r) solo está definido
   para r libre de cuadrados (la diagonalización es sobre r|P(z_A));
   Δ_r está definido para todo r, y su rasgo estructural es la
   condición p∥r. En particular Δ_r = 0 si ningún p∥r aporta — p. ej.
   r = p^k, k ≥ 2 (medido: r = 9, 25, 27 → Δ_r = 0, doc 97 §3.4) —
   un régimen que ni siquiera pertenece al índice de δ.
3. **Patrones de anulación opuestos.** δ_s(r) es **máximo**,
   δ_s(r) = −Π(r) exactamente, cuando r es invisible para el soporte
   (r no divide a ningún d ∈ V; doc 109 §3.4); Δ_r es **cero** en sus
   casos «invisibles» (sin p∥r). Uno mide la ausencia del soporte; el
   otro, la ausencia de caracteres inducidos mal contados.
4. **Dependencias disjuntas, demostradas [M/T, K7].** A (z=17, A=3)
   fijos, δ^{(L)}(35) salta de −0,825000 = −Π(35) exacto
   (L = 112, 138, 179, 250: 35 no divide a ningún d ∈ V) a +0,175000
   (L = 404, 954: S(35) = 1, solo d′=1 pasa el desplazamiento), con
   Π(35) = 0,825 fijo: δ es función del nivel L, un dato que la
   definición de Δ_r ni contiene. Y Δ_5 = (log 5)·⌊log_5 N⌋ crece con
   N (12,88 a N=10⁶; 27,36 a N=10¹², aritmética del Teorema 97.1)
   mientras δ_5 en el marco es O(Π(5)) = 0,688 y no recibe a N en
   absoluto (solo su soporte primo, vía excl). Fijado el marco, δ no
   ve variar Δ_r; fijado N, Δ_r no ve variar δ. Ninguna cota
   |δ| ≤ F(Δ) ni |Δ| ≤ F(δ) puede existir con F independiente de los
   datos ajenos.
5. **Lo único compartido es taxonomía, no matemática.** Ambos son
   correcciones de segundo orden a una identidad limpia (Π(r) → S_s(r)
   por restricción de soporte; inversión de Möbius sobre conductores
   por caracteres inducidos), ambos son efectos de borde (δ máximo
   fuera del soporte; Δ_r soportado en el borde p∥r de la
   factorización), y ambos deben llevarse explícitamente en el puente:
   Δ_r en el segundo factor (base E; su energía B_Δ es el objeto de
   L105(a), docs 105–107) y δ en el primero (base de coeficientes;
   docs 109–110). La premisa de que «viven sobre los mismos soportes
   desplazados» es solo una semejanza de notación: los objetos
   desplazados de Δ_r son los módulos r/p entre divisores; los de δ,
   los soportes V_s^{(r)} desplazados por r. Ahí termina el parecido.

## §5. Veredicto

1. **La expansión es exacta [T]:** Def2 = Σ_r (1/φ(r))Σ_{d′,e′∈𝒯^{(r)}}
   a(d′)a(e′), y su reindexación cierra en dos formas:
   (110.1) Def2 = Q + Q_completo − 2Π_full·M_V y
   (110.2) Def2 = (Q − Q_completo) + 2Π_full·M_cola, verificadas en
   240 niveles con error ≤ 1,3e-15 (K0–K3).
2. **Descomposición [M]:** en la expansión de la cola, (H2′) vive de
   la cancelación cruzada (D_cola = 1,53–1,86 × Q_c en los niveles
   peores del marco A; X_cola cancela el 38–69%); tras reindexar, el
   grueso del momento es el exceso de Gram (0,49–0,70 × Q_c) y el
   término de conteo es menor y decreciente (0,246 → 0,060).
3. **El conteo cierra M_cola por completo [T]:** 0 ≤ 2Π_full·M_cola
   ≤ 2·Q_completo, con M_full ≤ 1,1618 absoluto (A=3); y no cierra la
   diagonal de la cola (su cota pura pierde (ln z)²: 3,3–5,2 × Q_c y
   creciendo).
4. **Equivalencia [T]:** (H2) ⟺ (H2′) con constantes que difieren en
   ≤ 1 aditivo — la reducción del doc 109 era unidireccional y con
   pérdida ×2; la identidad muestra que la vía L² de los defectos no
   es un atajo sino el mismo enunciado. El término que se resiste es
   Q^{(L)} − Q_completo, y el sub-enunciado mínimo restante es
   **(H2″)**, el exceso de Gram, equivalente a (H2): lo que falta es
   un lema fundamental para la forma de Gram misma.
5. **δ_s(r) y Δ_r son independientes [T/M]:** datos disjuntos
   (nivel/soporte contra N/caracteres), índices distintos, anulación
   opuesta; la forma cerrada (110.2) no contiene a Δ_r en ningún paso.
   (H2) sigue abierta, reducida a su contenido mínimo (H2″).
   Goldbach no se declara demostrada.

## Apéndice. Cómo correr la verificación

```
cd ~/workspace/goldbach
python3 110_segundo_momento_defectos.py | tee salida_doc110_local.txt
```

Solo stdlib (math, sys, time, importlib para reutilizar las
definiciones del doc 109 y garantizar que son las mismas). Tiempo:
~0,9 s. Salida: PASS, chequeos:

- K0: validación cruzada contra LevelData del doc 109 (Q y Def2,
  acuerdo 0,00e+00 en 3 niveles).
- K1: partición S(r) + T(r) = Π(r) para todo r con soportes
  desplazados explícitos (error ≤ 1,11e-15); δ = −T.
- K2: intercambio Σ_r μ(r)Π(r)C(r) = Π_full·M_V e identidad (110.1);
  K2b: Σ_r Π(r)²/φ(r) = Q_completo (≤ 3,8e-16).
- K3: Q_completo = Π_full·M_full (Euler, suma y producto);
  M_cola directo = M_full − M_V (chequeo absoluto); identidad (110.2)
  en 240 niveles (error ≤ 1,3e-15 relativo a Q + Q_c).
- K4: descomposición D_cola/X_cola con 𝒯^{(r)} explícito: D + X = Def2
  (error absoluto ≤ 5,0e-16); tamaños contra Q_c por nivel.
- K5: 0 ≤ 2Π_full·M_cola ≤ 2Q_c en los 240 niveles; K5b: cota
  absoluta M_full ≤ 1,1618 y cota de conteo de D_cola (pierde (ln z)²).
- K6: equivalencia cuantitativa nivel a nivel:
  Def2/Q_c − (R−1) = 2·M_cola/M_full ∈ [0, 2] en los 240 niveles.
- K7: δ^{(L)}(35) contra L a Π(35) fijo (salto −Π → +0,175) y
  aritmética de Δ_5(N) del Teorema 97.1 (12,88 → 27,36): datos
  disjuntos.

---
*Goldbach sigue abierta. (H2′) queda disuelta por una identidad:
el segundo momento de los defectos es el exceso de la forma de Gram
sobre su producto de Euler, más una cola de conteo acotada entre 0 y
2·Q_completo. No hay atajo L² hacia (H2); el lema fundamental que
falta es para la forma de Gram misma, y su contenido mínimo es (H2″).*
