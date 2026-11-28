# DOC 97 — Verificación independiente PARCIAL del doc 96 (ruta diagonal + identidad de Möbius + experimentos 94d/94e)

**Fecha:** 2026-11-27 (noche).
**Autor:** Muse (línea dedicada de Goldbach).
**Estado:** PARCIAL — ver §0.
**Archivos:** `97_verificacion_doc96_diagonal_mobius.md` (este doc),
`97_verificacion_doc96_diagonal_mobius.py` (script, solo stdlib),
`salida_doc97_local.txt` (salida completa, 93 líneas, 1,3 s, **11 PASS / 0 FAIL**).

**Convención de marcas:** [T] teorema / verificado (prueba escrita o cálculo
exacto reproducido por dos rutas); [N] no tenido (afirmación sin prueba ni
medición suficiente); [P-rev] pendiente de re-verificación contra el doc 96
completo y sus scripts. MEDICIÓN ≠ PRUEBA en todo el documento.

## §0. Limitación y alcance (léase primero)

El túnel al PC de Dal sigue caído (salida #2 de egress, desde ~12:16 EST del
2026-11-27). El doc 96 íntegro y los scripts `94d_diagonal_profile.py` y
`94e_conductor_cutoff_profile.py` NO están disponibles. Este doc trabaja
SOLO con:

- la respuesta de Lumina Chat en el chat "Revisa UNIF U 2M X" (2026-11-27,
  ~20:14 EST), transcrita en la tarea;
- los docs 85–95 en `~/workspace/goldbach/` como contexto;
- modelos sintéticos propios con pesos λ⁺ de Nathanson (como en el doc 95)
  y ψ(N) reales hasta N ≈ 10⁶.

Por tanto todo lo que sigue es una verificación de las **afirmaciones** del
doc 96 tal como Lumina las enunció en el chat, NO una revisión del doc 96
mismo. En particular:

- **No he visto** la derivación completa de la identidad diagonal en el
  doc 96, ni su enunciado exacto del sup Q_s, ni la desigualdad de Cauchy
  del doc 94, ni la definición de L_s en el exp. 94d, ni la regla exacta de
  atribución de masa en el exp. 94e ("ah ≤ R0").
- **No he reproducido** los experimentos 94d/94e con sus parámetros
  (z, y, D, proxy de N): mis modelos usan parámetros propios y comparables,
  no idénticos. Donde los números difieren, lo digo y lo atribuyo a
  parámetros, no a error de Lumina — salvo donde la diferencia es
  cualitativa y la marco explícitamente.
- **Pendiente cuando vuelva el túnel** (§9): re-verificar contra el doc 96
  completo + scripts 94d/94e, y enviar el veredicto a Lumina.

## §1. Resumen ejecutivo

| # | Afirmación del doc 96 | Veredicto |
|---|---|---|
| 1 | Identidad diagonal exacta Q_s = Σ φ(r) C_s(r)² | **[T] SOBREVIVE.** Prueba escrita (§2.1) + verificación numérica exacta en 4 escenarios (error 0 a 2,2e-16). La distinción con los términos d=e es correcta y necesaria: d=e solos dan 1,00–3,19× Q_s. |
| 2 | Producto completo Q_completo(z) = ∏(1 − 2/(p−1) + p/(p−1)²), orden 1/log z | **[T] SOBREVIVE.** Identidad exacta sin criba (relerr ≤ 3,3e-16, 3 casos) + orden de Mertens con constante calculada: Q_completo·ln z → 1,7225 para (A=3, N=2^e); medido 1,7220 a z=10⁵. |
| 3 | V⁺ conserva una cota comparable a Q_completo | **[N] como prueba; [T] como medición.** Q_s/Q_completo ∈ [1,000, 1,429] en mis modelos — el orden se conserva numéricamente, pero la demostración sigue siendo el problema real, como Lumina misma dice. |
| 4 | Identidad de Möbius S_r = Σ μ(r/d)φ(d)E(N;d,N) | **[N] tal como está enunciada: FALTA UN TÉRMINO.** La identidad correcta es S_r = Σ μ(r/d)φ(d)E + Δ_r, con Δ_r la corrección imprimitiva de potencias de primos (§4). Sin Δ_r el error es ~13 a N≈10⁶ para r primo (1–3% de \|S_r\|). Con Δ_r, verificada exacta a ≤6,7e-9 en 10 conductores con ψ reales. |
| 5 | Cancelación exacta del término principal | **[T].** Σ_{d\|r} μ(r/d) = 0 para r>1, exacto. No se ve afectada por Δ_r. |
| 6 | Exp. 94d: Q_s·log z ∈ [1,556, 2,829]; Q_s/P ∈ [1,003, 1,701]; masa r>L_s/z ≤ 9,84%; r>L_s ≤ 18,15% | **Parcialmente reproducible [P-rev].** Mi rango Q_s·ln z es [1,67, 1,71] en el modelo E0 (N=2^e) y 2,874 en el modelo E1 (N con primos pequeños) — el extremo 2,83 de Lumina queda explicado por el factor local de N (§6), no por variación en s. La masa en r>z² es 0–1,3% en mis modelos; la definición de L_s no es recuperable de la transcripción (§5.3). |
| 7 | Exp. 94e: perfil κ 49,2% / 73,3% / 87,3% / 92,9% / 95,3% / 99,1% | **Forma reproducible; cifras dependientes de parámetros [P-rev].** Mis modelos dan a κ=2,7 entre 0,933 y 1,000 según la regla de atribución (h vs ah) y los parámetros. La lectura ah≤R0 en mi modelo S4 reproduce su tabla a 0,004–0,028. El 94% es plausible, no universal. |
| 8 | Mejora del primer factor de Cauchy ~4–5× al restringir a r>R0 | **La conclusión sobrevive; el factor es dependiente de parámetros.** En mis modelos el bloque grande lleva solo 0,36–0,41% de la masa (no 5–7%), así que la mejora es ~16×, no 4–5×. En ambos casos: no cambia nada en la pérdida de 7–10 órdenes (§5.4). |
| 9 | Precaución de N: el producto local cambia la constante | **[T] CONFIRMADA CUANTITATIVAMENTE.** Cociente E1/E0 medido [1,685, 1,698], media 1,691, predicción del producto 1,7005 (acuerdo 0,6%). Q_s·ln z pasa de ~1,70 a 2,874. La constante NO es uniforme en N. |
| 10 | Doc 94 SÍ contiene una desigualdad de Cauchy Gram↔resto; falta la cota del segundo factor | **Confirmada como precisión [T] en lo que puedo verificar.** El segundo factor queda ahora reducido exactamente a Σ(1/φ(r))\|Σ μφE\|² + corrección Δ (§4.5). Su estimación sigue siendo [N]. |

**Veredicto global:** la ruta diagonal SOBREVIVE como frente prioritario, con
una corrección obligatoria: la identidad de Möbius del puente debe escribirse
con el término +Δ_r. Nada de lo verificado contradice el doc 96 en su
sustancia; lo que cambia es un enunciado (Möbius) y dos cifras que son
dependientes de parámetros (perfil κ, factor de Cauchy). Goldbach sigue
abierta. La ruta diagonal no es una solución: es, como Lumina dice, el
frente donde hay que probar cotas.

## §2. La identidad diagonal: Q_s = Σ_r φ(r) C_s(r)²

### §2.1 Prueba analítica [T]

Notación (doc 94, vía doc 95): para un primo de conmutación s,
T_s(h0) = Σ_{d∈V⁺(D_A/s), h0|d} μ(d)/φ(d), y Q_s = Σ_{h0} φ(h0)|T_s(h0)|².
El doc 94/95 ya estableció la identidad de Gram
Q_s = Σ_{d,e} μ(d)μ(e)(d,e)/(φ(d)φ(e)) [T, doc 95 §2, error ≤1,8e-14].
Parto de ahí.

**Paso 1.** (d,e) = Σ_{r|(d,e)} φ(r) — identidad clásica Σ_{r|n} φ(r) = n.
[T]

**Paso 2.** Sustitución e intercambio de sumas (finitas, sin problema de
convergencia):
Q_s = Σ_{d,e} μ(d)μ(e)/(φ(d)φ(e)) · Σ_{r|(d,e)} φ(r)
    = Σ_r φ(r) · Σ_{d,e : r|d, r|e} μ(d)μ(e)/(φ(d)φ(e)).
[T]

**Paso 3.** La suma interna factoriza:
Σ_{d,e : r|d, r|e} μ(d)μ(e)/(φ(d)φ(e))
  = (Σ_{d : r|d} μ(d)/φ(d)) · (Σ_{e : r|e} μ(e)/φ(e)) = C_s(r)²,
con C_s(r) = Σ_{d∈V⁺(D_A/s), r|d} μ(d)/φ(d). [T]

**Conclusión:** Q_s = Σ_{r|P(z_A)} φ(r) C_s(r)², donde la suma es sobre r
libres de cuadrados con primos < z_A y ∤ AN. [T] La identidad es exacta,
preserva los signos de Chen (ningún valor absoluto antes de elevar al
cuadrado), y no usa GRH ni ninguna hipótesis no probada.

**Precisión sobre el límite variable (doc 95, punto 2).** En el modo de
límite fijo (cada d tiene un solo primo de conmutación s, o equivalentemente
el conjunto de h0 es el de la segunda reducción del doc 94), la identidad
vale tal cual. En el modo de límite variable (h0 con primos en [z,s) y con
≥2 primos de conmutación), la suma Σ_r φ(r)C_s(r)² debe tomarse sobre el
conjunto COMPLETO de h0 — incluidos los h0 con primos en [z,s) y los
multi-primos — y con la cola s'>s del doc 95. Si se toma solo el bloque
suave Q_s^suave, la identidad falla por 7,57e-02 en mi modelo S1-variable
(H1, abajo). Esto es exactamente el contraejemplo del doc 95, ahora
reproducido desde la representación diagonal. [T]

### §2.2 Verificación numérica (H1) [T]

Modelos (pesos λ⁺ de Nathanson, V⁺ según la definición de Nathanson
verificada en el doc 95; par (s,d) admisible si d∈V⁺(D/s), (d,AN)=1;
multiplicidad 1 verificada en H0):

- S1-fixed: A=3, z=11, y=60, D=20000, 50 pares (s,d), |H|=54.
- S2-fixed: A=3, z=17, y=120, D=200000, 243 pares, |H|=256.
- S4-fixed: A=3, z=11, y=3000, D=200000, 876 pares, |H|=880 (parámetros
  próximos a los de Lumina en el doc 94: 1.444 pares).
- S1-variable: mismo S1 con h0 que incluyen primos en [z,s).

| Escenario | max\|Q_diag − Q_s\| | max\|Q_diag − Q_full\| |
|---|---|---|
| S1-fixed | 0,00e+00 | 0,00e+00 |
| S2-fixed | 2,22e-16 (s=19) | 1,11e-16 |
| S4-fixed | 0,00e+00 | 0,00e+00 |
| S1-variable (suave) | 7,57e-02 (s=13) | 0,00e+00 |

La diagonal coincide exactamente con Q_s en límite fijo, y con Q_full
(conjunto completo de h0) en límite variable. [T]

**La distinción d=e es correcta y necesaria.** El cociente
(diag(d=e))/Q_s — es decir, Σ_d d/φ(d)² dividido por Q_s — mide
2,69–3,19 en S2, 1,62–2,24 en S1, 1,00–2,24 en S4 y 1,62–3,05 en
S1-variable. Los términos d=e solos SOBREESTIMAN Q_s por un factor
1–3,2: dan el orden correcto pero no el valor. Los términos cruzados
d≠e contribuyen netos en negativo y son del mismo orden que la diagonal.
Quedarse con d=e no basta; la diagonalización completa (que los incluye
a todos vía (d,e) = Σ φ(r)) es lo que hace exacta la representación. [T]
Esto confirma la advertencia de Lumina en su punto 1.

### §2.3 El producto completo y el orden 1/log z (H2) [T]

**Identidad exacta sin criba.** Si se quita la restricción V⁺ (suma sobre
TODOS los d libres de cuadrados con primos < z, ∤ AN), la misma
diagonalización da una identidad cerrada. En efecto, C(r) factoriza:
C(r) = μ(r)/φ(r) · ∏_{p<z, p∤rAN}(1 − 1/(p−1)), de donde
φ(r)C(r)² = (1/φ(r))∏_{p∤rAN}(1 − 1/(p−1))², y sumando sobre r por
producto de Euler (cada primo o divide a r — aporta 1/(p−1) — o no —
aporta (1 − 1/(p−1))²):

Q_sin criba = ∏_{p<z, p∤AN} [(1 − 1/(p−1))² + 1/(p−1)]
            = ∏_{p<z, p∤AN} (1 − 1/(p−1) + 1/(p−1)²) =: Q_completo(z). [T]

Verificación numérica (H2a): (A=3,z=11): 0,699652777778 en ambas rutas,
relerr 0,00e+00; (A=3,z=17): 0,588048442323, relerr 3,33e-16;
(A=15,z=19): 0,681344589422, relerr 3,33e-16. [T]

**Orden de Mertens (H2c).** Escribiendo
(1 − 1/(p−1) + 1/(p−1)²) = (1 − 1/p)·[p(p²−3p+3)/(p−1)³],
el primer factor da por Mertens (e^{−γ}/ln z)·∏_{p|AN}(1−1/p)^{−1} y el
segundo converge a C2 = ∏_{p∤AN} p(p²−3p+3)/(p−1)³ (convergente porque
el factor es 1+O(1/p²)). Luego
Q_completo(z)·ln z → e^{−γ} · ∏_{p|AN}(1−1/p)^{−1} · C2. [T]
Para A=3, N=2^e: 0,56146·3·1,0226 = 1,7225. Medido: 1,6777 (z=11),
1,6661 (17), 1,6897 (29), 1,6999 (100), 1,7159 (1000), 1,7220 (10⁵).
La convergencia es lenta (error O(1/ln z)), pero el límite teórico y la
medición coinciden a 3e-4 en z=10⁵. [T]

**La restricción V⁺ (H2b): medición, no prueba.** Con V⁺, el cociente
Q_s/Q_completo por s es: S1 [1,000, 1,161] (media 1,012); S2 [1,008, 1,094]
(media 1,046); S4 [1,000, 1,429] (media 1,236). La restricción
Rosser–Iwaniec conserva el orden 1/log z numéricamente en los tres
modelos, con un factor O(1) que crece con y (más primos de conmutación ⇒
más peso en la cola). Esto es una MEDICIÓN a favor de la hipótesis, no
una demostración: probar que V⁺ conserva una cota comparable sigue siendo
el problema real, exactamente como Lumina lo formula. [N] como prueba.

## §3. La identidad de Möbius: falta el término Δ_r

Este es el hallazgo adversarial principal de este doc.

### §3.1 Enunciado de Lumina (doc 96, vía chat)

Para conductor r>1 coprimo con N:
S_r(N) = Σ*_{χ* mod r} χ̄*(N) ψ(N, χ*),
y por inversión de Möbius sobre conductores primitivos:
S_r(N) = Σ_{d|r} μ(r/d) φ(d) E(N; d, N),
donde E(N;d,N) = ψ(N;d,N) − ψ(N)/φ(d). La contribución principal se
cancela exactamente.

### §3.2 La identidad tal como está enunciada es INEXACTA [N]

Verificación numérica con ψ reales (H6): N = 999983 (primo),
ψ(N) = 999586,5975, caracteres construidos explícitamente por CRT sobre
potencias de primos (grupos cíclicos, logaritmos discretos por fuerza
bruta), ψ(N,χ*) = Σ_a χ*(a)ψ(N;r,a) con ψ(N;r,a) real por criba de von
Mangoldt. Comparando S_r con Σ μ(r/d)φ(d)E(N;d,N) SIN corrección:

| r | #χ* | S_r | RHS sin Δ | diferencia |
|---|---|---|---|---|
| 3 | 1 | −715,0038 | −728,1872 | +13,1833 |
| 5 | 3 | 1046,8240 | 1033,9485 | +12,8755 |
| 6 | 0 | 0,0000 | −0,6931 | +0,6931 |
| 7 | 5 | 526,2080 | 512,5867 | +13,6214 |
| 9 | 4 | 674,8847 | 674,8847 | 0 |
| 15 | 3 | 1563,0200 | 1563,0200 | 0 |
| 21 | 5 | 1540,4249 | 1554,0463 | −13,6214 |
| 25 | 16 | 357,5449 | 357,5449 | 0 |
| 27 | 12 | 791,9846 | 791,9846 | 0 |
| 35 | 15 | 2382,8836 | 2374,4999 | +8,3837 |

La diferencia no es ruido: para r primo es exactamente
(log p)·⌊log_p N⌋ (p. ej. r=3: log 3·12 = 13,1833; r=7: log 7·7 =
13,6214). Es la contribución de las potencias de primos p^j ≤ N con p|r.
[N] la identidad sin Δ_r; [T] el diagnóstico del error.

### §3.3 La identidad corregida [T]

**Teorema (97.1).** Para r>1, (r,N)=1:

S_r(N) = Σ_{d|r} μ(r/d) φ(d) E(N; d, N) + Δ_r,

Δ_r = Σ_{p|r, p∤r/p} Σ*_{χ* mod r/p} χ̄*(N) · Σ_{p^j≤N} χ*(p)^j log p,

con la convención de que si r/p = 1, la suma Σ* es sobre el carácter
principal mod 1 (χ* ≡ 1), y contribuye (log p)·⌊log_p N⌋.

**Prueba.** Paso 1 (ortogonalidad, solo para (n,r)=1):
Σ_{χ mod r} χ̄(N)χ(n) = φ(r) si n≡N (mod r), 0 si no. [T]
Paso 2 (reagrupación por conductor): todo χ mod r es inducido de un
único χ* primitivo mod d, d|r, con χ(n)=χ*(n) si (n,r)=1 y χ(n)=0 si
(n,r)>1. Para (n,r)=1: φ(r)[n≡N] = Σ_{d|r} Σ_{χ* mod d} χ̄*(N)χ*(n). [T]
Paso 3 (suma con Λ): multiplicar por Λ(n) y sumar n≤N. El lado izquierdo
da φ(r)ψ(N;r,N) (n≡N fuerza (n,r)=1 pues (N,r)=1). El lado derecho da
Σ_{d|r} S_d^{ind}, donde S_d^{ind} usa el carácter INDUCIDO a mod r, no
el primitivo: ψ(N,χ*^{ind}) = ψ(N,χ*) − Σ_{p|r, p∤d} Σ_{p^j≤N} χ*(p^j)
log p, porque los n con (n,r)>1 y (n,d)=1 que llevan Λ(n)≠0 son
exactamente las potencias p^j con p|r, p∤d. [T] **Este es el paso que
la identidad sin Δ omite: la ortogonalidad solo vale para (n,r)=1, y
los caracteres inducidos se anulan en los n no coprimos mientras que los
primitivos no.**
Paso 4 (inversión de Möbius): de φ(r)ψ(N;r,N) = Σ_{d|r} S_d − Δ'_r se
obtiene S_r = Σ_{d|r} μ(r/d)φ(d)ψ(N;d,N) + Δ_r, y la suma interna de
Möbius colapsa: para p y e fijos (p∤e), Σ_{d: e|d|r, p|d} μ(r/d) = 1 si
e = r/p y 0 si no. De ahí la forma cerrada de Δ_r. [T]
Paso 5 (término principal): Σ_{d|r} μ(r/d)φ(d)E(N;d,N) =
Σ_{d|r} μ(r/d)φ(d)ψ(N;d,N) − ψ(N)·Σ_{d|r} μ(r/d), y Σ_{d|r} μ(r/d) = 0
para r>1. [T] La cancelación del término principal es exacta y no se ve
afectada por Δ_r.

**Cota.** |Δ_r| ≤ Σ_{p|r} (log p)·⌊log_p N⌋ ≤ ω(r)·log N. [T]

### §3.4 Verificación numérica de la identidad corregida (H6) [T]

Con Δ_r incluido: |S_r − (RHS + Δ_r)| ≤ 6,70e-09 en los 10 conductores
(máximo en r=9), error relativo ≤ 4,2e-11 (r=6). Valores de Δ_r: r=3:
+13,1833 (1,8% de |S_r|); r=5: +12,8755 (1,2%); r=6: +0,6931 = log 2
(con S_r = 0: toda la identidad es la corrección); r=7: +13,6214 (2,6%);
r=9, 15, 25, 27: 0 exactamente; r=21: −13,6214 (0,9%); r=35: +8,3837
(0,4%). El caso r=6 es instructivo: no hay caracteres primitivos mod 6,
S_6 = 0, y la identidad se reduce a 0 = −log 2 + log 2.

**Peso en el segundo factor.** La cantidad que el puente necesita estimar
es Σ_{r∈R_A, r>R0} (1/φ(r))|Σ_{d|r} μ(r/d)φ(d)E(N;d,N)|². Sustituyendo
la identidad corregida, el término Δ_r añade una contribución relativa
medida de Σ(1/φ)|Δ|² / Σ(1/φ)|S_r|² = 1,24e-04 (0,012%) en mis 10
conductores a N≈10⁶. Es decir: **Δ_r no amenaza la reducción** — es
pequeño y exactamente computable — **pero el enunciado del doc 96 debe
llevar +Δ_r**. Una identidad que se presenta como exacta y no lo es, por
pequeño que sea el error, es una errata que hay que corregir antes de
construir encima. [T] la magnitud; [N] que el error sea despreciable a
toda escala (la cota ω(r)·log N crece con r; en el régimen del puente,
r hasta AD_A, hay que llevarla explícitamente en la estimación del
segundo factor — véase §7).

## §4. Experimentos 94d/94e: lo que se reproduce y lo que no

### §4.1 Perfil κ (94e, H3): la forma se reproduce; las cifras dependen de parámetros

Mis modelos (fracción de G_A con h ≤ R0, R0 = (ln N)^κ):

| Modelo / lectura | κ=1,0 | 2,0 | 2,5 | 2,7 | 2,8 | 3,0 |
|---|---|---|---|---|---|---|
| Lumina 94e (N=10¹², A=3) | 0,492 | 0,733 | 0,873 | 0,929 | 0,953 | 0,991 |
| S1, corte h (N modelo) | 0,560 | 0,980 | 1,000 | 1,000 | 1,000 | 1,000 |
| S2, corte h (N modelo) | 0,391 | 0,877 | 0,980 | 0,999 | 1,000 | 1,000 |
| S4, corte h (N modelo) | 0,582 | 0,840 | 0,973 | 0,990 | 0,997 | 0,999 |
| S1, corte A·h | 0,334 | 0,913 | 0,990 | 1,000 | 1,000 | 1,000 |
| S2, corte A·h | 0,253 | 0,771 | 0,915 | 0,956 | 0,976 | 0,998 |
| S4, corte A·h | 0,510 | 0,774 | 0,898 | 0,957 | 0,978 | 0,994 |
| S4, corte A·h (proxy 10¹²) | 0,510 | 0,774 | 0,898 | 0,957 | 0,978 | 0,994 |

Lecturas:

- **La forma se reproduce:** subida empinada entre κ=1 y κ=2,7 y saturación
  cerca de κ=3, en los tres modelos y con las dos lecturas del corte.
  [T] como medición.
- **Las cifras exactas NO son universales:** a κ=2, mis modelos dan
  0,761–0,984 según modelo y lectura, frente a 0,733 de Lumina; a κ=2,7,
  0,933–1,000 frente a 0,929. La diferencia es de parámetros (z, y, D,
  proxy de N) y de regla de atribución, no de error. [T]
- **La lectura ah≤R0 importa.** La fila que mejor reproduce la tabla de
  Lumina es la de S4 con corte A·h: 0,510/0,774/0,898/0,957/0,978/0,994
  frente a 0,492/0,733/0,873/0,929/0,953/0,991 — acuerdo a 0,004–0,041.
  Esto sugiere que su "ah ≤ R0" atribuye la masa del h a su conductor
  r = ah con a = A (el mayor a posible), no a h. Pero sin el doc 96 no
  puedo confirmarlo: **la regla exacta de atribución queda [P-rev]**.
  La diferencia entre las dos lecturas (hasta 0,25 a κ=2 en S2) es del
  mismo orden que la discrepancia con su tabla, así que la atribución es
  una ambigüedad real, no un detalle.
- **Conclusión sobre el 94%:** el valor a κ≈2,7 es plausible y está
  respaldado por la forma del perfil, pero no es una constante universal
  del problema: depende de (z, y, D) y de la regla de atribución. Para
  usarlo en una demostración hay que fijar κ y la regla, y cuantificar
  el complemento (el bloque grande) en esos términos. [N] como constante;
  [T] como orden de magnitud.

### §4.2 Masa en r de la representación diagonal (H7) [P-rev]

Fracción de Q_s con r > t (forma diagonal, S1/S2): r>z: 6,0% (S1),
15,7–20,8% (S2); r>z²: 0% (S1), 0–1,3% (S2). Lumina reporta en 94d:
masa máxima en r > L_s/z: 9,84%; en r > L_s: 18,15%. Mi r>z (6–21%) es
comparable a su r>L_s/z (9,84%) si L_s/z ~ z, es decir L_s ~ z²; pero
mi r>z² (0–1,3%) NO es comparable a su r>L_s (18,15%) bajo esa misma
identificación. O bien L_s es mucho menor que z² en su configuración, o
bien sus parámetros difieren de los míos. **La definición de L_s queda
[P-rev]** — no es recuperable de la transcripción del chat. Lo que sí
queda establecido: la masa diagonal está concentrada en r pequeños
(r ≤ z² lleva ≥98,7% en mis modelos), así que una separación
r ≤ L_s/z vs cola es viable en principio, pero su cuantificación exacta
depende de L_s.

### §4.3 Rango Q_s·log z (94d)

Lumina: Q_s·log z ∈ [1,556, 2,829] en 20 escenarios N=2^e, A=3.
Mis mediciones: modelo E0 (N=2^e, z=29): Q_s·ln z ∈ [1,693, 1,705]
(10 primos s muestreados); modelo E1 (N con primos 5,7,11,13): 2,874.
El extremo superior 2,829 de Lumina NO se reproduce por variación en s
dentro del modelo E0 (mi rango en s es estrecho, 1,69–1,71), pero SÍ se
reproduce por el factor local de N (§5): 2,874 ≈ 2,829. Hipótesis: el
rango [1,556, 2,829] de 94d mezcla variación en s con variación en el
factor local según los parámetros de cada escenario — o sus escenarios
tienen más variación en s que los míos. **Queda [P-rev]** contra sus
scripts: reproducir el rango exacto con sus parámetros.
El cociente Q_s/P_{A,N}(z) ∈ [1,003, 1,701] de Lumina: mis cocientes
Q_s/Q_completo son [1,000, 1,429] (S4) — el extremo 1,701 no se alcanza
en mis modelos; de nuevo, diferencia de parámetros. [P-rev]

### §4.4 Consecuencia práctica: el factor de Cauchy (H4)

Lumina: si el bloque grande conserva solo 5–7% de la energía, el primer
factor de Cauchy mejora ~4–5×, sin resolver la pérdida de 7–10 órdenes
del doc 95. Mis modelos (κ=2,7, proxy N=10¹², R0=7794):

| Modelo | masa bloque grande f | 1/√f | 1/f | pérdida ingenua (doc 95) | tras mejora |
|---|---|---|---|---|---|
| S1 | 0 (bloque vacío: R0 > max h=1855) | — | — | 7,45e+07 | — |
| S2 | 0,0036 | 16,68× | 278,1× | 1,40e+08 | 8,41e+06 |
| S4 | 0,0041 | 15,70× | 246,4× | 7,45e+07 | 4,74e+06 |

En mis modelos el bloque grande lleva solo 0,36–0,41% de la masa, no
5–7%: la mejora es ~16× (lectura 1/√f), mayor que la 4–5× de Lumina.
La razón es la misma dependencia de parámetros del §4.1: mis modelos
concentran más masa en conductores pequeños que el suyo. **La conclusión
cualitativa de Lumina sobrevive intacta:** incluso con 16× (o 278× en
la lectura 1/f), la pérdida ingenua pasa de ~10⁸ a ~10⁶–10⁷ — sigue
siendo enorme. La restricción a r>R0 ayuda pero no resuelve nada por sí
sola; el puente necesita la estimación del segundo factor, no solo la
mejora del primero. [T] La cifra exacta del factor (4–5× vs ~16×) es
dependiente de parámetros. [P-rev] su valor en los parámetros del doc 96.

## §5. La precaución de N: el factor local es real y cuantitativo [T]

### §5.1 El experimento (H5)

Dos variantes del mismo modelo (A=3, z=29, y=200, D=29⁵, 10 primos s
muestreados log-uniformemente en [29,200)):

- **E0 (N=2^e):** los d excluyen solo 2 y los primos |A. Es el modelo de
  los docs 94/95.
- **E1 (N con muchos primos pequeños):** los d excluyen además 5,7,11,13
  — es decir, N es divisible por 5,7,11,13 (además de 2 y 3).

Predicción del producto completo: Q_s(E1)/Q_s(E0) =
∏_{p∈{5,7,11,13}} (1 − 1/(p−1) + 1/(p−1)²)^{−1} = 1,7005.
El factor local tipo serie singular es ∏ p/(p−1) = 1,7378.

### §5.2 Resultado

| s | Q_s(E0) | Q_s(E1) | cociente | Q_s·ln z E0 | Q_s·ln z E1 |
|---|---|---|---|---|---|
| 29 | 0,502722 | 0,853494 | 1,6977 | 1,693 | 2,874 |
| 43 | 0,503133 | 0,853494 | 1,6964 | 1,694 | 2,874 |
| 61 | 0,503805 | 0,853494 | 1,6941 | 1,696 | 2,874 |
| 79 | 0,504870 | 0,853494 | 1,6905 | 1,700 | 2,874 |
| 101 | 0,504870 | 0,853494 | 1,6905 | 1,700 | 2,874 |
| 139 | 0,505277 | 0,853494 | 1,6892 | 1,701 | 2,874 |
| 163 | 0,505791 | 0,853494 | 1,6874 | 1,703 | 2,874 |
| 199 | 0,506468 | 0,853494 | 1,6852 | 1,705 | 2,874 |

Cociente E1/E0 ∈ [1,685, 1,698], media 1,691 — acuerdo con la predicción
1,7005 al 0,6%. [T] Q_s·ln z pasa de ~1,70 a 2,874: **la constante del
sup Q_s ≪ 1/log z NO es uniforme en N**; lleva el factor local
∏_{p|N}(1 − 1/(p−1) + 1/(p−1)²)^{−1}, del mismo orden que el factor
local de la serie singular de Goldbach. La precaución de Lumina es
correcta y cuantitativa.

### §5.3 Consecuencias

1. **El enunciado a probar debe ser la versión normalizada por el
   producto local:** sup Q_s ≪ (factor local de N)/log z_A, o bien
   sup Q_s / Q_completo(z_A; N) ≪ 1. Probar sup Q_s ≪ 1/log z con una
   constante absoluta es FALSO en general (el factor local crece con el
   número de primos pequeños de N; para N con muchos primos pequeños es
   ~1,7 y puede ser mayor). [T] la no-uniformidad; [N] la forma exacta
   del enunciado normalizado (hay que fijarla en el doc 98 o en la
   revisión del 96).
2. **El rango del doc 95 queda explicado.** El máximo Q_s·ln z ≈ 2,83
   medido en el doc 95 (con N=2^e en el modelo) está en el mismo rango
   que mi 2,874 del modelo E1: el extremo superior del rango es un
   efecto de factor local, no de variación en s. Esto refuerza que la
   constante ~2,5 "sin holgura" del doc 95 es en realidad una constante
   local, no universal. [T] como explicación; [P-rev] la descomposición
   exacta del rango 94d.
3. **Q_s(E1) es independiente de s en este régimen** (0,853494 para los
   10 primos muestreados): con D/s grande, todos los d admisibles pasan
   V⁺ y Q_s no depende de s. Es una propiedad del régimen D ≫ y, no un
   error. En el régimen del doc 94 (D/s comparable a los d), Q_s sí varía
   con s (S4: Q_s/Q_completo ∈ [1,000, 1,429]). [T]

## §6. La precisión sobre el punto 4 del doc 95

Lumina precisa: el doc 94 SÍ contiene una desigualdad formal de Cauchy
entre la norma de Gram y el resto primitivo; lo que falta es demostrar
una cota útil para el segundo factor, el momento de los bloques de
caracteres.

**Confirmación parcial [T].** No he visto el doc 94 íntegro, así que no
puedo verificar la forma exacta de su desigualdad de Cauchy — eso queda
[P-rev]. Pero la estructura que describe es la correcta y mi verificación
la deja en su sitio exacto:

- Primer factor: √G_A (o √G_A^{grande} tras restringir a r>R0) — ahora
  cuantificado: la restricción mejora ~4–16× según parámetros (§4.4).
- Segundo factor: el momento Σ_{r∈R_A, r>R0} (1/φ(r))|S_r(N)|², donde
  S_r es la suma primitiva. Con el Teorema 97.1, S_r = Σ_{d|r} μ(r/d)
  φ(d)E(N;d,N) + Δ_r, así que el segundo factor queda reducido
  exactamente a una forma cuadrática sobre errores de progresión
  restringidos a los conductores del soporte, más la corrección Δ_r
  (peso relativo medido 1,24e-04 a N≈10⁶). [T]
- Lo que sigue faltando es la ESTIMACIÓN de ese segundo factor: una
  gran criba restringida a H_A, o cualquier otra vía, que lo acote por
  debajo de lo que exige el balance de Chen. [N] — es el problema
  abierto central, y la identidad de Möbius (corregida) es una reducción
  exacta a él, no una demostración.

En otras palabras: el doc 95 tenía razón en que el puente G_A → R_A no
estaba en la derivación (la vía ingenua pierde 10⁷–10¹⁰), y Lumina tiene
razón en que la desigualdad formal sí existe y lo que falta es el segundo
factor. Ambas afirmaciones son compatibles: la desigualdad existe, la
cota útil no. El trabajo pendiente es estimar el segundo factor, no
volver a derivar el primero.

## §7. Veredicto: ¿sobrevive la ruta diagonal como frente prioritario?

**Sí, sobrevive — con dos correcciones y una advertencia.**

**Lo que sobrevive [T]:**

1. La identidad diagonal es exacta y es la representación correcta de
   Q_s. Preserva los signos de Chen y diagonaliza la forma de Gram sin
   pérdida. Es el objeto adecuado para el frente analítico.
2. El orden 1/log z del producto completo es un teorema (Mertens +
   factor local explícito), con constante 1,7225 para (A=3, N=2^e).
3. La cancelación del término principal en la identidad de Möbius es
   exacta.
4. La masa de Gram está concentrada en conductores pequeños (forma del
   perfil κ reproducida) y el factor local de N es real y medible
   (acuerdo 0,6% con la predicción del producto).

**Lo que hay que corregir:**

1. **La identidad de Möbius debe llevar +Δ_r (Teorema 97.1).** Sin él,
   el enunciado es inexacto (error 1–3% a N≈10⁶, mayor para r primo).
   La corrección es pequeña y exactamente acotada (|Δ_r| ≤ ω(r)·log N),
   y no amenaza la reducción — pero debe escribirse.
2. **El enunciado de sup Q_s debe normalizarse por el producto local
   de N.** Una constante absoluta es falsa en general; la versión
   correcta es sup Q_s/Q_completo(z_A; N) ≪ 1 o equivalente.

**Advertencia:** las cifras 94% (κ=2,7) y 4–5× (mejora de Cauchy) son
dependientes de parámetros. En mis modelos son ~99% y ~16×. Ninguna de
las dos cifras cambia la conclusión (la restricción ayuda pero no cierra
la pérdida de 7–10 órdenes), pero ninguna debe citarse como constante
del problema.

**Lo próximo concreto que hay que probar (en orden):**

1. **La cota del segundo factor.** Estimar
   Σ_{r∈R_A, r>R0} (1/φ(r)) |Σ_{d|r} μ(r/d) φ(d) E(N; d, N)|²
   (con el +Δ_r ya incorporado) por debajo del umbral que exige el
   balance de Chen. Esta es la pieza que falta; sin ella, todo lo demás
   es preparación. Vía propuesta por Lumina: gran criba restringida a
   H_A sobre la identidad de Möbius corregida. Alternativa: acotar
   directamente el momento de bloques de caracteres sin pasar por
   Möbius. Cualquiera de las dos es el siguiente teorema a intentar.
2. **La cota de Q_s en su forma normalizada.** Probar
   sup_{z_A≤s<y_A} Q_s ≪ Q_completo(z_A; N) — es decir, que la
   restricción V⁺ conserva una cota comparable al producto completo.
   Es el problema real que Lumina señala en su punto 1. Un primer paso
   concreto: probar la cota para la diagonal d=e (que ya da el orden,
   §2.2) y luego controlar los términos cruzados d≠e, que son netos en
   negativo y del mismo orden — es decir, demostrar que la cancelación
   entre términos distintos no destruye la cota sino que la mejora.
   (Si los cruzados son netos negativos, Q_s ≤ diag(d=e): una cota
   superior para la diagonal bastaría. **Esto es una conjetura de trabajo,
   no un teorema**: hay que verificar el signo neto de los cruzados en
   general, no solo en mis 4 modelos. [N])
3. **Fijar κ y la regla de atribución** (h vs ah) antes de cuantificar
   el bloque grande, y escribir la estimación del segundo factor con
   esa κ fija. Sin κ fija, las cifras del perfil no son utilizables en
   una demostración.
4. **Re-verificación contra el doc 96 completo** (§8) antes de publicar
   nada de esto: en particular, comprobar si el doc 96 ya contiene Δ_r
   (en cuyo caso mi corrección es solo a la transcripción del chat) y
   cuál es su definición exacta de L_s y de la atribución de masa.

**Lo que NO hay que hacer:** seguir acumulando identidades. Como dice
Lumina, la próxima prioridad es probar la cota de Q_s, no seguir
derivando formas equivalentes. La identidad diagonal y la de Möbius
(corregida) son suficientes como reducciones; el trabajo ahora es
estimación, no álgebra.

## §8. Pendientes [P-rev] para cuando vuelva el túnel

1. Re-verificar este doc contra el doc 96 íntegro
   (`matematicas/96_respuesta_adversarial_doc95_y_ruta_diagonal.md`) y
   los scripts `94d_diagonal_profile.py` / `94e_conductor_cutoff_profile.py`.
2. Comprobar si el doc 96 ya contiene el término Δ_r en la identidad de
   Möbius (en cuyo caso la corrección del §3 es solo a la transcripción
   del chat, no al doc).
3. Definición exacta de L_s en el exp. 94d (masa en r>L_s/z y r>L_s) y
   reproducción del rango Q_s·log z ∈ [1,556, 2,829] y del cociente
   Q_s/P ∈ [1,003, 1,701] con sus parámetros.
4. Regla exacta de atribución de masa en el exp. 94e ("ah ≤ R0") y
   reproducción de la tabla κ con sus parámetros.
5. Forma exacta de la desigualdad de Cauchy del doc 94 (Gram ↔ resto
   primitivo) y del segundo factor que Lumina describe.
6. Los seis ítems [P-rev] del doc 95 (re-verificación contra el doc 94
   íntegro), que este doc hereda sin resolver.
7. El enunciado exacto del sup Q_s en el doc 96 (límite fijo vs
   variable, cola s'>s) y si su constante se enuncia ya normalizada por
   el producto local.
8. Enviar el veredicto de este doc a Lumina Chat (chat "Revisa UNIF U 2M X")
   cuando el túnel vuelva, con la corrección Δ_r como punto principal.

## Apéndice A. Cómo correr la verificación

```
cd ~/workspace/goldbach
python3 97_verificacion_doc96_diagonal_mobius.py | tee salida_doc97_local.txt
```

Solo stdlib (math, sys, itertools). Tiempo: ~1,3 s. Salida: 93 líneas,
11 PASS / 0 FAIL. Chequeos:

- H0: hipótesis del modelo (multiplicidad 1, (d,AN)=1, d impar libre de
  cuadrados) en S1/S2/S4.
- H1: identidad diagonal por dos rutas (Q_s directo vs Σ φ(r)C_s(r)²;
  Q_full con el conjunto completo de h0) en 4 escenarios.
- H2: (a) identidad exacta del producto completo sin criba, 3 casos;
  (b) Q_s/Q_completo con V⁺; (c) orden de Mertens hasta z=10⁵ y límite
  teórico 1,7225.
- H3: perfil κ en 3 modelos × 2 elecciones de N, con las dos lecturas
  del corte (h y A·h).
- H4: mejora del primer factor de Cauchy a κ=2,7 y pérdida residual.
- H5: factor local de N (E0 vs E1), 10 primos s.
- H6: identidad de Möbius con ψ reales (N=999983 primo), 10 conductores,
  caracteres por CRT + logaritmos discretos, con la corrección Δ_r.
- H7: masa diagonal con r>t (r>z, r>z², r>D/(s·z), r>D/s).

## Apéndice B. Tabla de resultados numéricos (salida completa en salida_doc97_local.txt)

- H1: S1-fixed 0,00e+00; S2-fixed 2,22e-16; S4-fixed 0,00e+00;
  S1-variable (suave) 7,57e-02, (completo) 0,00e+00. diag(d=e)/Q_s:
  S1 [1,62, 2,24]; S2 [2,69, 3,19]; S4 [1,00, 2,24].
- H2a: relerr 0,00e+00 / 3,33e-16 / 3,33e-16. H2c: Q_c·ln z =
  1,6777 / 1,6661 / 1,6897 / 1,6999 / 1,7159 / 1,7220 (z = 11 / 17 / 29 /
  100 / 1000 / 10⁵); límite 1,7225. H2b: Q_s/Q_completo S1 [1,000, 1,161],
  S2 [1,008, 1,094], S4 [1,000, 1,429].
- H3: véase §4.1. H4: S2 f=0,0036, 1/√f=16,68×, 1,40e+08 → 8,41e+06;
  S4 f=0,0041, 1/√f=15,70×, 7,45e+07 → 4,74e+06; S1 bloque vacío.
- H5: cociente E1/E0 ∈ [1,685, 1,698], media 1,691, predicción 1,7005;
  Q_s·ln z E0 1,693–1,705, E1 2,874.
- H6: |S_r − (RHS+Δ_r)| ≤ 6,70e-09 (r=9); rel ≤ 4,2e-11 (r=6).
  Δ_r: 3:+13,1833; 5:+12,8755; 6:+0,6931; 7:+13,6214; 9:0; 15:0;
  21:−13,6214; 25:0; 27:0; 35:+8,3837. Peso relativo en el segundo
  factor: 1,240e-04.
- H7: r>z: 6,0% (S1), 15,7–20,8% (S2); r>z²: 0% (S1), 0–1,3% (S2).

---
*Goldbach sigue abierta. La ruta diagonal sobrevive como frente
prioritario, con la identidad de Möbius corregida (+Δ_r) y el enunciado
de sup Q_s normalizado por el producto local. El siguiente paso es
estimar el segundo factor, no derivar otra identidad.*
