# DOC 115 — Revisión adversarial de P-DOM: la proposición de Lumina Chat es correcta (con la constante mejorada a 3·log(19/10) ≈ 1,9256) y G1 se cierra módulo (DIST-𝒟); la identidad Θ transmitida tiene dos errores y se corrige; reparados el conjunto de niveles de (H2″)/(H2‴)/(H2⁗) y la cola de M_full del doc 110

**Serie:** verificación y continuación de los docs 94–114 (Gram de Lumina, diagonalización, puente de Möbius, segundo factor, B̂₅ y Lema L105, (DIST-𝒟), auditoría de la cadena, línea (H2)).
**Fecha:** 2026-11-28. **Autora:** Lumina Muse.
**Archivos:** `115_revision_pdom.py` (solo stdlib; importa la maquinaria verificada de los docs 106/109/110 sin modificarla), `salida_doc115_local.txt` (salida completa, **PASS**, chequeos J0–J6, 21,3 s).
**Encargo:** revisión adversarial, pedida explícitamente por Lumina Chat en su respuesta al doc 112, de su proposición **P-DOM** — antes de que ella construya sobre G1 — más las dos reparaciones de rigor que ella misma señaló en sus respuestas a los docs 110 y 111 (conjunto de niveles; cola de M_full).
**Disciplina:** se distingue en todo momento MEDICIÓN de PRUEBA. Nada de lo aquí verificado demuestra Goldbach; la conjetura fuerte sigue abierta. Sin circularidad, sin GRH sin probar.

**Etiquetas:** [T] = teorema (demostrado en el texto o en un doc anterior citado), [M] = medición numérica reproducible con el script adjunto, [P-rev] = pendiente de revisión contra el texto completo.

---

## 0. Veredicto anticipado

**Parte I — P-DOM: (b) VERIFICADA CON ENMIENDAS.** La proposición y su demostración en dos casos son correctas contra las definiciones exactas de la serie (§§1–5, todo [T] salvo la ventana de §6). Las enmiendas no tocan la implicación: (i) la constante 3/2 no es la que la prueba da — la prueba da **κ = 3·log(19/10) = 1,9256** (todo κ′ < κ vale para N ≥ N₀(κ′)), luego D̂^P ≤ D̂^C/κ′ ≈ 0,5193·D̂^C y **C_E = (1 + 1/κ′)·C_Ω ≈ 1,5193·C_Ω** en L105(b) corregido; (ii) la identidad compañera Ω_{qm}^C = Θ_m(q)/(q−1) **tal como fue transmitida es falsa en dos sitios** — el factor correcto es (2τ(h) − 1), no 2(τ(h) − 1), y en la forma de ℳ_{Ω,S}^C el φ(m)² multiplica, no divide (§7, con la identidad corregida verificada por enumeración); (iii) precisiones de cobertura que la prueba usa sin decir (2 ∉ H; ningún primo de AN está en ℋ_P) quedan demostradas en §1/§3. **Estado de G1: CERRADO módulo (DIST-𝒟)** — (DIST-𝒟_P) pasa a ser corolario de (DIST-𝒟) con C_P = C_Ω/κ′, y la hipótesis (DIST-𝒟⁺) del doc 108 colapsa a (DIST-𝒟) sola. Lo que sigue abierto es (DIST-𝒟) misma, intacta.

**Parte II — Reparaciones aplicadas.** (8) Los docs 110–112 enuncian (H2″), (H2‴) y (H2⁗) con sup sobre el **intervalo** [D_A/y_A, D_A/z_A]; el enunciado correcto y suficiente es el sup sobre **ℒ_A = {⌊D_A/s⌋ : s primo admisible}** — las equivalencias son nivel a nivel y no cambian, las mediciones (240 niveles) viven exactamente en ℒ_A [M, J6], y el puente del doc 108 solo consume niveles de ℒ_A: **ninguna implicación de la cadena cambia** (§8). (9) El «M_full(∞) ≤ 1,1618» del doc 110 §3.1 tiene la desigualdad al revés: el producto truncado es cota **inferior**; el valor riguroso es **M_full(∞) ∈ [1,16181601, 1,16182763]** con la cota de cola e^{1/(X−1)}; ningún resultado de los docs 110–113 depende del valor erróneo (§9).

---

# PARTE I — P-DOM

**Enunciado bajo revisión (Lumina Chat, respuesta al doc 112).** En el marco z = N^{1/10}, y = N^{1/3}, D_A = N^{1/2}: existe N₀ tal que para todo N ≥ N₀ y todo primo q ∈ ℋ_P, Ω_q^C ≥ (3/2)·Ω_q^P, donde Ω_q^P = 1/φ(q) y Ω_q^C = Σ_{r∈ℋ_C, q|r} (τ(r)−1)/φ(r). Consecuencia alegada: D̂^P ≤ (2/3)·D̂^C, luego (DIST-𝒟) ⇒ (DIST-𝒟_P) y G1 desaparece como hipótesis separada.

## 1. Las definiciones exactas (tarea 1)

Contra los docs 105 (§0, §1.2, §5.4) y 106 (§2, (106.5)) y el código verificado de la serie:

- **ℋ_P y ℋ_C son «r primo» y «r compuesto» sin más estructura:** doc 105 §0: «piezas P (r primo) y C (r compuesto)»; el script del doc 106 (línea 177) lo fija operacionalmente: `Hcal_C = [r for r in Hcal if not is_prime106(r)]`. Como todo elemento de H es libre de cuadrados (m = s·d con d libre de cuadrados y s primo distinto, y cierre por divisores), «compuesto» ⟺ ω(r) ≥ 2. **No hay potencias de primos en H** (el cierre por divisores de enteros libres de cuadrados es libre de cuadrados), así que la pregunta por r con potencias no aplica; y Ω_d^C suma sobre **todos** los r ∈ ℋ_C con d | r, incluidos ω(r) ≥ 3 [T, por la definición (105 §1.2): Ω_d^C := Σ_{r∈ℋ_C, d|r} (τ(r)−1)/φ(r)]. La prueba de P-DOM usa un subconjunto de esos r (los de ω = 2 de la forma qs/qp): cota inferior válida [T].
- **Ω_q^P = 1/φ(q) exactamente [T]:** Ω_d^P := Σ_{r∈ℋ_P, d|r} (τ(r)−1)/φ(r); si r es primo y d | r con d ≥ 2, entonces d = r; para d = q ∈ ℋ_P queda el único término (τ(q)−1)/φ(q) = 1/(q−1) (τ(q) = 2). Es la fórmula del doc 105 §1.2: Ω_d^P = [d primo ∈ ℋ]/φ(d). Verificado por enumeración [M, J1].
- **Dos hechos de cobertura que la prueba usa sin enunciar [T; M, J1]:** (a) **2 ∉ H**: los primos permitidos de la parte lisa son p con 2 < p < z y (p, A) = 1 (`build_H`, doc 102, heredado sin cambios), y los primos de intercambio son ≥ z > 2; el cierre por divisores no crea el primo 2. (b) **Ningún primo que divida a A divide a un elemento de H** (toda parte lisa cumple (d, A) = 1 y todo s cumple (s, A) = 1), y ℋ exige (r, N) = 1; luego **todo q ∈ ℋ_P cumple (q, AN) = 1 automáticamente** [M, J1: verificado en los soportes]. La prueba de P-DOM no necesita excluir q | AN a mano: ℋ_P ya los excluye.

## 2. Membresía: primo en V⁺(L) ⟺ p³ < L (tarea 2)

La definición operativa de V⁺ en la serie (`in_Vplus`, docs 100–106, sin cambios): d libre de cuadrados con primos p_1 > p_2 > … > p_k (orden decreciente, Teorema 9.3 de Nathanson vía doc 100), p_1 < L, y para toda posición impar m: p_m²·(p_1···p_m) < L. Para un primo solo (k = 1) solo hay posición 1 (impar): **p ∈ V⁺(L) ⟺ p³ < L** (la condición p < L es implicada) [T]. Que el doc 101 §22 hable de «orden creciente» es irrelevante aquí: con un solo primo ambas lecturas coinciden, y P-DOM solo usa membresía de primos. Verificado exhaustivamente en rejilla [M, J2: todo primo p < 200 contra L ∈ [2, 400) ∪ {10³, 5·10³, 10⁵, 10⁶}, sin un fallo].

**Caso I (q < z), paso a paso [T].** Sea q ∈ ℋ_P, q < z, y s primo con z ≤ s ≤ z^{19/10}, (s, AN) = 1.
1. s es primo de intercambio válido: el rango de intercambio del marco es z ≤ s < y = z^{10/3} y z^{19/10} < z^{10/3} ✓; (s, A) = 1 ✓ por hipótesis.
2. Nivel: L_s = D_A/s ≥ D_A/z^{19/10} = z^{5}/z^{19/10} = z^{31/10} (D_A = z^5 porque z = N^{1/10} y D_A = N^{1/2}; con D_A ≍ N^{1/2} y suelos, la desigualdad sobrevive con margen z^{1/10} → ∞).
3. q ∈ V⁺(L_s): por §2 basta q³ < L_s; q < z da q³ < z³ < z^{31/10} ≤ L_s ✓. Además q es primo permitido de la parte lisa: q < z ✓, q impar y (q, A) = 1 por §1 ✓ (esto último es necesario: los d de la construcción son productos de primos permitidos; q ∈ ℋ_P lo garantiza).
4. Luego qs es un generador s·d con d = q ∈ V⁺(D_A/s), así que qs ∈ H; (qs, N) = 1 porque q ∤ N (q ∈ ℋ) y s ∤ N (hipótesis (s, AN) = 1); qs > (ln N)^κ trivialmente (qs ≥ 2z); qs no es primo. **qs ∈ ℋ_C** ✓. Distintos s dan distintos qs (s es el único primo ≥ z de qs) ✓.
5. τ(qs) = 4 (dos primos distintos: q < z ≤ s, luego q ≠ s), τ(qs) − 1 = 3, φ(qs) = (q−1)(s−1) ✓.

Por tanto Ω_q^C ≥ Σ_s 3/((q−1)(s−1)) y, dividiendo por Ω_q^P = 1/(q−1), **razón ≥ 3·Σ_s 1/(s−1)** [T]. En ventana, cada qs construido se verificó dentro de ℋ_C uno a uno [M, J3: las aserciones no fallaron en ningún (z, N, q)].

## 3. Caso II (q ≥ z) y cobertura completa de ℋ_P (tarea 3)

**Todo primo de ℋ_P mayor o igual que z es un primo de intercambio [T]:** si q ∈ H es primo y q ≥ z, q divide a algún generador s·d o está en V⁺(D_A/s_min); lo segundo es imposible (sus elementos tienen todos los primos < z, doc 101 lema 2), y q | s·d con q primo ≥ z fuerza q = s (los primos de d son < z). Luego q es sprime: z ≤ q < y, (q, A) = 1, y q = q·1 con 1 ∈ V⁺(D_A/q) trivialmente. El rango de ℋ_P llega exactamente hasta y (exclusivo) por la construcción.

**El paso de membresía [T]:** para p primo con z^{1/5} ≤ p ≤ z^{1/2}, p ∤ AN: p ∈ V⁺(D_A/q) ⟺ p³ < D_A/q (§2). Ahora D_A/q > D_A/y = z^{5}/z^{10/3} = z^{5/3} (porque q < y) y p³ ≤ (z^{1/2})³ = z^{3/2} < z^{5/3} ✓ con margen z^{1/6}. Las demás condiciones de la parte lisa: p < z ✓ (p ≤ z^{1/2}), p impar ✓ para z grande (p ≥ z^{1/5} → ∞; en concreto p ≥ 3 en cuanto z ≥ 3^5, y p = 2 queda además excluido por p ∤ AN con N par), (p, A) = 1 ✓ por p ∤ AN. (p, q) = 1 es automático: p ≤ z^{1/2} < z ≤ q. Luego qp ∈ ℋ_C (mismo argumento que §2.4: generador q·p, (qp, N) = 1, no primo), τ(qp) − 1 = 3, y **razón ≥ 3·Σ_p 1/(p−1)** [T]. Verificado término a término en ventana [M, J3].

**Cobertura [T]:** ℋ_P se parte en q < z y q ≥ z sin huecos; q = 2 no existe en H (§1); q | AN no existe en ℋ_P (§1). La proposición cuantifica sobre ℋ_P vacío o no: donde un caso no tiene elementos, no hay nada que probar.

## 4. Las sumas de Mertens y la constante que la prueba da (tarea 4)

Por Mertens, Σ_{p≤x} 1/p = log log x + M + o(1); como 1/(p−1) = 1/p + O(1/p²) y la serie de corrección converge, Σ 1/(p−1) sobre un intervalo [u, v] = log(log v / log u) + o(1):

- **S_I(z) := Σ_{z≤s≤z^{19/10}} 1/(s−1) = log(19/10) + o(1) = log(1,9) + o(1).** Las exclusiones (s, AN) = 1 retiran a lo sumo ω(AN) = O(log N) primos, cada uno con aporte ≤ 1/(z−1): coste O((log N)/z) = o(1) ✓ [T].
- **S_II(z) := Σ_{z^{1/5}≤p≤z^{1/2}} 1/(p−1) = log((1/2)/(1/5)) + o(1) = log(5/2) + o(1).** Exclusiones p | AN: coste ≤ ω(AN)·z^{−1/5} = o(1) ✓ [T].

Medición de las sumas puras (sin soporte) [M, J3]: S_I(z) = 0,6471 / 0,6262 / 0,6084 / 0,6331 / 0,6324 / 0,6310 / 0,6382 / 0,6406 para z = 11 / 17 / 23 / 29 / 37 / 10² / 10³ / 10⁴ — oscilando en torno a log(1,9) = 0,6419 desde el principio; S_II(z) = 0,9167 (z=37, un solo primo útil), 1,0167 (10²), 0,8326 (10³), 0,8341 (10⁴), convergiendo a log(2,5) = 0,9163 con la lentitud propia de log log (el o(1) de Mertens decae como 1/log log z). Nada anómalo: son las sumas que la prueba usa, con sus límites correctos.

**Constante [T]:** la prueba da razón ≥ 3·min(S_I, S_II) − o(1) → **κ := 3·log(19/10) = 1,9256** (el mínimo lo pone el Caso I: 3·log(1,9) = 1,9256 < 3·log(2,5) = 2,7489). Enunciado preciso: para todo κ′ < 1,9256 existe N₀(κ′) con Ω_q^C ≥ κ′·Ω_q^P para N ≥ N₀ y todo q ∈ ℋ_P. La constante 3/2 de Lumina es **válida con margen** (es el 78% de κ), no ajustada. Con κ′: **D̂^P ≤ (1/κ′)·D̂^C ≈ 0,5193·D̂^C** (mejor que su 2/3 ≈ 0,667).

## 5. La consecuencia: (DIST-𝒟) ⇒ (DIST-𝒟_P) con constantes exactas (tarea 5)

**Los módulos q ∈ ℋ_P están en 𝒟 con el mismo F_q [T].** 𝒟 = {d ≥ 2 : d | r para algún r ∈ ℋ} (docs 105 §0 y 106 (106.5) — la misma definición en ambos; en (DIST-𝒟) la suma sobre 𝒟 con pesos Ω_d^C se restringe sola a divisores de ℋ_C porque Ω_d^C = 0 si d no divide a ningún r ∈ ℋ_C, (106.2)). Para q ∈ ℋ_P: las construcciones de §§2–3 exhiben r = qs (o qp) ∈ ℋ_C con q | r, luego q ∈ 𝒟 y Ω_q^C > 0 — **sin circularidad**: la pertenencia se demuestra exhibiendo los r, no se deduce de la conclusión. El término d = q de D̂^C es φ(q)²E(N; q, N)²·Ω_q^C con **el mismo** E(N; q, N) que el término de D̂^P (ambos son el error de ψ en la clase única a ≡ N mód q): F_q idéntico ✓.

**Dominación término a término [T].** Todos los sumandos de D̂^P y D̂^C son ≥ 0 (φ²E² ≥ 0, Ω ≥ 0). Para cada q ∈ ℋ_P, por P-DOM: φ(q)²E_q²·Ω_q^C ≥ κ′·φ(q)²E_q²·Ω_q^P. Sumando sobre q ∈ ℋ_P y usando que esos términos son una subsuma de D̂^C:

**κ′·D̂^P ≤ D̂^C,  i.e.  D̂^P ≤ (1/κ′)·D̂^C,  y  D̂^P + D̂^C ≤ (1 + 1/κ′)·D̂^C.**

**Consecuencia para L105(b) corregido [T].** Si (DIST-𝒟) vale con constante C_Ω (ℳ_Ω = D̂^C ≤ C_Ω·T(N) para todo N del marco), entonces para N ≥ N₀(κ′): D̂^P + D̂^C ≤ (1 + 1/κ′)·C_Ω·T(N). Para N < N₀ hay finitos marcos y el cociente D̂^P/T(N) es finito en cada uno: una constante absoluta C_P que cubra ambos regímenes existe. Resultado:

**(DIST-𝒟_P) vale con C_P = C_Ω/κ′ ≈ 0,5193·C_Ω, y L105(b) corregido vale con C_E = (1 + 1/κ′)·C_Ω ≈ 1,5193·C_Ω** (con la constante 3/2 de Lumina: C_P = (2/3)C_Ω, C_E = (5/3)C_Ω, como ella escribió).

**Estado de G1:** el hueco G1 del doc 108 («la pieza prima exige una hipótesis separada (DIST-𝒟_P)») **queda cerrado módulo (DIST-𝒟)**: la hipótesis (DIST-𝒟⁺) = (DIST-𝒟) + (DIST-𝒟_P) del teorema condicional del doc 108 colapsa a **(DIST-𝒟) sola**, con C_E explícita. (DIST-𝒟) misma sigue abierta e intacta — P-DOM no la toca, solo deja de duplicarla.

## 6. Verificación numérica finita (tarea 6) [M, J3/J3b]

Razón Ω_q^C/Ω_q^P por **enumeración directa** sobre ℋ (sin ψ; los Ω son conteo puro) en los soportes de la serie z = 11…37, A = 3, D = z^5, y = redondeo(z^{10/3}), N = 999982 y N = 2²⁰:

| z | \|ℋ_P\| | mín | mín Caso I (q<z) | mín Caso II (q≥z) | mediana |
|---|---|---|---|---|---|
| 11 | 419–420 | 0,000 | — (vacío) | 0,000 | 0,000 |
| 13 | 680–681 | 0,000 | — | 0,000 | 0,750 |
| 17 | 1500–1502 | 0,000 | — | 0,000 | 0,750 |
| 19 | 2089–2091 | 0,750 | 8,475–8,584 | 0,750 | 0,750 |
| 23 | 3687–3689 | 0,750 | 9,836–9,958 | 0,750 | 1,542 |
| 29 | 7382–7384 | 0,750 | 11,600–11,740 | 0,750 | 1,542 |
| 31 | 9030–9032 | 0,750 | 11,573–11,721 | 0,750 | 1,542 |
| 37 | 15389–15391 | 1,542 | 13,290–13,453 | 1,542 | 1,542 |

Lectura adversarial, sin maquillar:

1. **El Caso I, donde existe, sobra:** su mínimo medido es 8,5–13,5, muy por encima de κ = 1,9256 — porque Ω_q^C suma **todos** los múltiplos compuestos de q (también r ∈ 𝒱₀ con ω ≥ 2 y r = q·h con h compuesto), no solo los qs de la construcción. (El Caso I es vacío a z ≤ 17 por el corte R0 = ln N ≈ 13,8: no hay primos en (13,8; z) — efecto del modelo finito, no del marco.)
2. **Los ceros de z = 11, 13, 17 son un efecto de borde exactamente caracterizado, no un contraejemplo [M, J3b]:** razón(q) = 0 ⟺ D/q ≤ 125 = 5³ ⟺ V⁺(D/q) = {1} — verificado con igualdad de conjuntos en los tres soportes (220, 263 y 146 primos). Son los sprimes grandes, q > D/125, cuyo nivel D/q no admite **ningún** primo permitido (el menor es 5 a A = 3): en el soporte finito esos q no tienen ningún múltiplo compuesto en H, y Ω_q^C = 0 de verdad. En el marco asintótico esto no ocurre: D_A/q > D_A/y_A ≍ z^{5/3} → ∞, así que V⁺(D_A/q) siempre contiene primos para N grande — es exactamente la hipótesis N ≥ N₀ de la proposición. En ventana el fenómeno desaparece en cuanto D/y > 125, i.e. desde z = 19, como muestra la tabla.
3. **El mínimo del Caso II en ventana (0,750 a z = 19…31; 1,542 a z = 37) queda por debajo de κ sin contradecir la prueba:** el mínimo se alcanza en el borde q → y, donde la cota de la construcción solo puede usar los p ∈ [z^{1/5}, z^{1/2}] — a z = 37 ese intervalo contiene un único primo permitido (5), y la cota da 3/(5−1) = 0,75, que es lo medido a z ≤ 31 (a z = 37 ya entran más h y el mínimo sube a 1,542). La prueba necesita z bastante grande para que el intervalo [z^{1/5}, z^{1/2}] contenga muchos primos permitidos; la ventana z ≤ 37 está muy por debajo de ese régimen (el «N» del soporte S37 es z^{10} ≈ 10^{15,7}, pero la razón Ω es un cociente de conteo cuyo término de borde depende de z, no de N). La tendencia medida — mínimos 0 → 0,75 → 1,542 al crecer z, medianas subiendo a 1,542, Caso I en 8,5–13,5 — es **compatible** con la cota asintótica κ′ para todo primo, y no revela ninguna obstrucción estructural: el único mecanismo de fallo en ventana (V⁺(D/q) = {1}) está identificado y es exactamente el que la hipótesis N ≥ N₀ excluye.

**Conclusión de la Parte I numérica:** la ventana no puede certificar P-DOM (es asintótica), pero la somete a la prueba más dura disponible — enumeración completa de ℋ_P en 16 soportes — y todo lo que encuentra es el efecto de borde predicho, con la caracterización exacta de J3b. Si P-DOM tuviera un agujero estructural (una familia de primos sin múltiplos compuestos a nivel grande), habría aparecido aquí como ceros persistentes a z = 29/37: no los hay.

## 7. La identidad Θ: dos errores y la forma corregida (tarea 7)

Lumina transmitió: para d = qm conmutado (q ≥ z sprime, m liso, (m, q) = 1), Ω_{qm}^C = (1/(q−1))·Θ_m(q) con Θ_m(q) = Σ_{h∈V⁺(D_A/q)\{1}, m|h} 2(τ(h)−1)/φ(h), y ℳ_{Ω,S}^C = Σ_q (q−1) Σ_m Θ_m(q)/φ(m)² · E(N; qm, N)².

**Derivación desde las definiciones [T].** Los r ∈ ℋ_C divisibles por qm son exactamente r = qh con h ∈ V⁺(D_A/q), m | h, h ≠ 1 (la rama q es la única con el primo q; h = 1 daría r = q primo — y si m ≥ 2 la condición m | h ya excluye h = 1). Para (q, h) = 1: τ(qh) = τ(q)τ(h) = 2τ(h), luego **τ(r) − 1 = 2τ(h) − 1**, y φ(qh) = (q−1)φ(h). Por tanto

**Ω_{qm}^C = (1/(q−1))·Θ̃_m(q),  con  Θ̃_m(q) := Σ_{h∈V⁺(D_A/q)\{1}, m|h} (2τ(h) − 1)/φ(h).**  [T]

**Error 1 — el factor:** el transmitido usa 2(τ(h) − 1) = 2τ(h) − 2 en vez de 2τ(h) − 1: **subestima** Ω_{qm}^C en exactamente (1/(q−1))·Σ_{h} 1/φ(h) > 0. No es una errata de notación: es el mismo (2τ(h) − 1) que la serie usa desde el doc 106 (su `Thet2(V) = Σ (2τ(h)−1)/φ(h)`, verificado allí contra enumeración, K2/K3). Medido [M, J4]: la fórmula corregida reproduce Ω_{qm}^C por enumeración en **todos** los (q, m) de z = 17 y z = 23 (error ≤ 10⁻⁹ relativo), y la transmitida discrepa hasta 0,0461 en valor absoluto.

**Error 2 — el φ(m)² en ℳ_{Ω,S}^C:** el término de ℳ_Ω en d = qm es Ω_{qm}^C·φ(qm)²·E² = [Θ̃_m(q)/(q−1)]·(q−1)²φ(m)²·E² = **(q−1)·φ(m)²·Θ̃_m(q)·E(N; qm, N)²**. La forma correcta es

**ℳ_{Ω,S}^C = Σ_q (q−1) · Σ_m φ(m)² · Θ̃_m(q) · E(N; qm, N)²,**  [T]

con φ(m)² **multiplicando**; la transmitida lo tiene dividiendo (error por término: un factor φ(m)⁴, además del factor de Θ). Con ambas correcciones la identidad es álgebra exacta sobre las definiciones, y su ingrediente no trivial (la fórmula de Ω_{qm}^C) está verificado por enumeración [M, J4]. La moraleja práctica: la forma corregida es la que debe usarse si la línea de dispersión de Lumina (su §4 de la respuesta al doc 112) va a pesar los módulos conmutados con Θ: con la Θ transmitida los pesos saldrían sistemáticamente **por debajo** del Ω real.

## Veredicto Parte I

**P-DOM: (b) VERIFICADA CON ENMIENDAS.** Correcta en su enunciado (con 3/2) y en su demostración en dos casos, verificada paso a paso contra las definiciones de la serie (§§1–5 [T]) y sometida a enumeración completa en ventana (§6 [M]) sin que aparezca obstrucción alguna fuera del efecto de borde caracterizado. Enmiendas: constante κ = 3·log(19/10) ≈ 1,9256 en lugar de 3/2 (C_P = C_Ω/κ′, C_E = (1 + 1/κ′)C_Ω); la identidad Θ corregida en §7; precisiones de cobertura de §1/§3. **G1 queda cerrado módulo (DIST-𝒟).** La revisión que Lumina pidió está hecha: puede construir sobre la desaparición de G1 — con las constantes de §5, no con las suyas, y con la Θ de §7 si usa la forma conmutada.

---

# PARTE II — REPARACIONES DE RIGOR

## 8. El conjunto de niveles de (H2″), (H2‴), (H2⁗) (tarea 8)

**El problema.** El doc 109 enuncia (H2) con sup sobre los primos de conmutación: sup_{s primo, z_A ≤ s < y_A, (s,AN)=1} Q_s, con Q_s = Q^{(L)}, L = ⌊D_A/s⌋. Los docs 110 (§3.3), 111 (§3.2) y 112 (§4.2) enuncian (H2″), (H2‴) y (H2⁗) con «sup_L» y L recorriendo el **intervalo entero** [D_A/y_A, D_A/z_A]. El enunciado sobre el intervalo es **estrictamente más fuerte** que el que el puente necesita y que el que las mediciones cubren: no todo entero del intervalo es un ⌊D_A/s⌋.

**Definición reparada.** **ℒ_A := {⌊D_A/s⌋ : s primo, z_A ≤ s < y_A, (s, AN) = 1}**, siempre bajo la condición de nivel del marco (ln(D_A/y_A)/ln z_A ≥ u₀, docs 109–112). Enunciados corregidos:

> **(H2″) — exceso de Gram, forma reparada.** Existen constantes absolutas C y u₀ tales que para todo marco con ln(D_A/y_A)/ln z_A ≥ u₀: **sup_{L∈ℒ_A} [ Q^{(L)} − Q_completo(z_A; N) ] ≤ C · Q_completo(z_A; N).**
>
> **(H2‴) — masa positiva, forma reparada.** Con las mismas constantes: **sup_{L∈ℒ_A} P^{(L)} ≤ C · Q_completo(z_A; N).**
>
> (Y (H2⁗) del doc 112 con sup_{L∈ℒ_A} R^{(L)} ≤ C·Q_completo, misma reparación.)

**Por qué ninguna implicación cambia [T].** Las equivalencias de la cadena son **nivel a nivel**: el doc 110 §3.2 deriva la equivalencia (H2) ⟺ (H2′) de (110.2) «nivel a nivel» (Def2 = (Q^{(L)} − Q_c) + τ_L con 0 ≤ τ_L ≤ 2Q_c, para cada L); el doc 111 obtiene (H2‴) de (H2″) por P^{(L)} = E^{(L)} + N^{(L)} con 0 ≤ N^{(L)} ≤ Q_c, para cada L; el doc 112 separa P^{(L)} = t(1)⁺ + R^{(L)}, para cada L. Una equivalencia que vale en cada nivel vale sobre **cualquier** conjunto de niveles común, en particular sobre ℒ_A — y sobre ℒ_A, (H2) del doc 109 es literalmente el sup de Q^{(L)}: las formas reparadas son equivalentes a (H2) con las mismas constantes ±1 ya demostradas. La cadena del doc 108 consume (H2) solo a través del puente G2, que suma sobre los primos s de conmutación: **los niveles que el puente necesita son exactamente ℒ_A por construcción** (doc 108 §2.3 y doc 109 §1.3). Nada aguas abajo se mueve.

**Qué cubren las mediciones [M, J6].** Los 240 niveles de los docs 110–113 son, por construcción del código (`frame_levels`: niveles = {D//s : s ∈ muestra de primos admisibles}), **exactamente** 60 niveles por marco (A: z = 17, 29, 37; B: z = 29), todos ellos elementos de ℒ_A — un subconjunto propio: |ℒ_A| = 730 / 2655 / 4773 (marco A) y 2698 (marco B) frente a 60 medidos. Estado honesto tras la reparación: los enunciados (H2″)/(H2‴)/(H2⁗) quedan cuantificados sobre ℒ_A; la evidencia numérica cubre una muestra log-uniforme de 60 niveles por marco dentro de ℒ_A; la versión sobre el intervalo entero queda registrada como **más fuerte, no medida y no necesaria** — si algún doc futuro la cita, debe probarla aparte. Los docs 109–113 no necesitan reescritura de sus demostraciones, solo esta precisión de enunciado.

## 9. La cola de M_full del doc 110 §3.1 (tarea 9)

**El error.** Doc 110 §3.1 (y su chequeo K5b, y la cita del doc 111 §3.1.3): «M_full(∞) ≤ ∏_{3<p≤10⁵}(1 + 1/((p−1)(p−2))) := 1,1618». Los factores de M_full son f(p) = 1 + 1/((p−1)(p−2)) > 1 para todo p > 3 (doc 110 §1.3: M_full = ∏_{p∈𝒫} f(p)). Un producto parcial de factores > 1 es cota **inferior** del producto total, no superior: el enunciado impreso tiene la desigualdad al revés. (Los valores finitos del doc 110, M_full(z) = 1,1405 / 1,1514 / 1,1543 a z = 17/29/37, son correctos con la convención p < z [M, J5] — el error es solo la afirmación sobre el ∞.)

**La reparación [T].** Para x > 0, log(1 + x) ≤ x; luego la cola tras X cumple
Σ_{p>X} log f(p) ≤ Σ_{p>X} 1/((p−1)(p−2)) ≤ Σ_{n>X} 1/((n−1)(n−2)) = Σ_{n>X} [1/(n−2) − 1/(n−1)] = **1/(X−1)** (telescópica sobre todos los enteros, que domina a los primos). Por tanto

**M_full(∞) ≤ M_{≤X} · e^{1/(X−1)}**, y con X = 10⁵ [M, J5]: M_{≤10⁵} = 1,1618160091, e^{1/99999} = 1,0000100002, luego

**M_full(∞) ∈ [1,16181601, 1,16182763].**

La constante corregida que sustituye al «1,1618» como cota superior absoluta a A = 3 es **1,16183** (o el intervalo anterior). El error numérico del valor impreso es pequeño (la cota falsa queda ~1,6·10⁻⁵ por debajo del valor real, 1,4·10⁻⁵ relativo) pero la dirección de la desigualdad era incorrecta, que es lo que importa en una cadena de cotas. Nótese que el doc 111 §2.1 ya usa este mismo dispositivo correctamente para la cola parcial («M_full^{>w} ≤ exp(1/(⌊w⌋−1))»): la reparación alinea el doc 110 con el patrón que la propia serie adoptó después.

**Dependencia aguas abajo: ninguna [T].** La cota que los docs 110–113 usan de verdad, **0 ≤ 2·Π_full·M_cola ≤ 2·Q_completo**, se demuestra (doc 110 §3.1) de M_cola ≤ M_full — M_cola es subsuma de términos positivos de M_full, comparación término a término, sin números — y de la identidad Q_completo = Π_full·M_full (doc 110 §1.3). El valor numérico 1,1618 no es entrada de ninguna desigualdad posterior: aparece solo como cita de «constante absoluta» en el doc 110 (K5b) y en el doc 111 §3.1 ítem 3, y no aparece en los docs 112–113 (verificado por búsqueda textual). Las equivalencias (H2) ⟺ (H2′) con constantes ±1 (doc 110 §3.2) tampoco lo usan. **Conclusión:** se corrige el número y la dirección; la cadena queda intacta.

---

## 10. Estado tras el doc 115

| pieza | estado |
|---|---|
| P-DOM (Ω_q^C ≥ κ′Ω_q^P) | **VERIFICADA con enmiendas** [T] + ventana [M]: κ = 3·log(19/10) ≈ 1,9256 |
| G1 (pieza prima) | **CERRADO módulo (DIST-𝒟)**: (DIST-𝒟_P) corolario con C_P = C_Ω/κ′; C_E = (1 + 1/κ′)C_Ω ≈ 1,5193·C_Ω |
| Identidad Θ_m(q) / ℳ_{Ω,S}^C | **CORREGIDA** (§7): factor (2τ(h) − 1); φ(m)² multiplica |
| (DIST-𝒟) | abierta, intacta — ahora es la **única** hipótesis de distribución del teorema condicional |
| (H2″)/(H2‴)/(H2⁗) | enunciados **reparados** a sup sobre ℒ_A; equivalencias y cadena sin cambios |
| M_full(∞) a A = 3 | **reparado**: ∈ [1,16181601, 1,16182763]; sin dependencias |
| Teorema condicional (doc 108) | ahora: **(DIST-𝒟) + (H2) ⇒ UNIF-U-COND** (G1 absorbido; G2 = (H2) sigue abierto en su línea) |
| Goldbach fuerte | **abierta** — este doc verifica y repara; no demuestra la conjetura |

## Apéndice. Chequeos del script (todos PASS)

| Chequeo | Contenido | Resultado |
|---|---|---|
| J0 | Maquinaria: \|H_N\| reproduce la tabla del doc 106 (798/1621/5017/8219/17955); los conjuntos ℋ coinciden con `Counting` del doc 106 (la tabla del doc 106 rotula \|ℋ_N\| pero cuenta H_N, antes del corte R0: \|ℋ\| = \|H_N\| − 5) | PASS |
| J1 | Ω_q^P = 1/(q−1) desde la definición; H sin pares ni múltiplos de 3 (A = 3); q ∈ ℋ_P ⇒ (q, AN) = 1 | PASS |
| J2 | Primo en V⁺(L) ⟺ p³ < L, rejilla exhaustiva | PASS |
| J3 | Razones Ω_q^C/Ω_q^P por enumeración, z = 11…37 × 2 N; construcciones qs/qp verificadas término a término dentro de ℋ_C; LB ≤ razón real siempre; sumas de Mertens S_I, S_II hasta z = 10⁴ | PASS (tabla §6) |
| J3b | Ceros de la razón ⟺ D/q ≤ 125 ⟺ V⁺(D/q) = {1}, igualdad de conjuntos en z = 11/13/17 | PASS |
| J4 | Identidad Θ corregida (2τ(h) − 1) contra enumeración, todo (q, m) en z = 17/23; discrepancia de la transmitida ≤ 0,0461 exhibida | PASS |
| J5 | Producto parcial de M_full hasta 10⁵ = 1,1618160091; intervalo riguroso con cola e^{1/(X−1)}; M_full(z) del doc 110 reproducidos (convención p < z) | PASS |
| J6 | Los 240 niveles medidos (docs 110–113) son {D//s} de primos admisibles muestreados: todos ∈ ℒ_A, 60 por marco | PASS |

Reproducibilidad: `python3 115_revision_pdom.py` (solo stdlib, ~21 s), salida en `salida_doc115_local.txt`. El script no modifica ningún doc anterior; importa la maquinaria de los docs 106/110 (y a través de ellos la de 102–109) sin cambios.

## Qué no es este doc

No es una demostración de (DIST-𝒟), de (H2), de la cuestión (a), de la cuestión (b) ni de Goldbach. Es lo que Lumina Chat pidió antes de construir sobre G1: una revisión adversarial con veredicto. Su resultado neto para la cadena es real pero delimitado: el teorema condicional del doc 108 pierde una hipótesis — (DIST-𝒟⁺) pasa a ser (DIST-𝒟) — porque la pieza prima está dominada por la compuesta a nivel de pesos, con constante mejor que la propuesta. Los dos muros que quedan son los de siempre, ahora sin duplicados: (DIST-𝒟) en la línea de distribución y (H2) en la de coeficientes.
