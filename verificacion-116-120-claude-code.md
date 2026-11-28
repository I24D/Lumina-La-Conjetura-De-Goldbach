# Verificación de las afirmaciones [T] y [M] de los docs 116–120a (línea (H2))

**Agente:** Claude Code (Anthropic, modelo `claude-opus-5-5`), brazo verificador de la línea Goldbach.
**Fecha:** 2026-11-28.
**Encargo:** de Muse, por orden de Dal: el lote propuesto en el §4 de `verificacion-116-121-claude-code.md`, en orden de dependencia 116 → 117 → 118 → 119 → 120a.

**Método.**
- Leí íntegros los docs 116, 116a, 116b, 116c, 117, 117a, 117b, 118, 118a, 118b1, 118b2, 119, 119a, 119b, 120 y 120a. También sus dependencias de definición: docs 109, 110, 113 y 114.
- Reimplementé desde las definiciones escritas, **sin importar código de nadie**. Del script del doc 109 solo tomé, leyéndola, la regla de muestreo de niveles (`sample_sprimes`): no está escrita en los docs y hace falta para reproducir sus conteos exactos.
- Los contraejemplos usan predicados propios y aritmética exacta (`Fraction`), y los niveles se comprueban contra la definición de ℒ_A.
- Script: `numerico/verificacion_116_120_claude_code.py` (solo stdlib, 4 s). Salida: `numerico/resultados/20261128_verificacion_116_120_claude_code.txt`.
- Resultado: **76 chequeos; 69 confirman el texto y 7 no lo confirman.** Los 7 son hallazgos y se explican abajo.

**Estado:** PROPUESTA hasta el contraste de Muse.
**Etiquetas:** [T] teorema, [M] medición, [P] pendiente.

**Aviso (regla 4 del README):** el código y el método son independientes, pero **no es una revisión independiente en sentido estricto**: corro en el mismo usuario de Windows que los demás brazos.

**Goldbach sigue abierta.** Nada de este informe la demuestra ni cierra (H2).

---

## §0. Veredicto en breve

**El grueso del lote se sostiene.** Se confirman con código propio:
- todos los teoremas centrales: (REM-PERS-small), la falsedad puntual de (R-FLOW-71+), la falsedad de (BAND) cruda, (REG-Π), (REG-Π-LOW), el Teorema A y la reindexación de (SLAB-V);
- todas las cifras decisivas que pediste: el margen ×28, 0,17–0,46·Q_c, la mediana 62,7 % y las tablas de 117b, 118a, 118b2, 119a, 119b y 120a.

**Hallazgos que cambian algo** (detalle en §1–§6; los que piden reclasificación están en la tabla del §6):

1. **Doc 116b §4 (contraejemplo R_R):** el nivel del ejemplo **no está en ℒ_A**. Usa s = 5 < z = 6052, y u_L = 4,815 > 4. El enunciado "R_R ≤ 1 es falso en marco" **sobrevive con un ejemplo propio legítimo** (§1.2).
2. **Docs 116 §2 y 116b §0 ("R_A > 1 solo en la configuración libre: y primo…"):** es **falso**. Construí un **y compuesto** con R_A(y) = 1,00220 > 1 dentro de ℒ_A (§1.2).
3. **Docs 117 y 117a: (FC-2) no es "estrictamente más débil" que (FLOW-COMP), es estrictamente más fuerte.** Es una condición suficiente. Además, el bloqueador "entrada de fuerza PNT" **no existe** para la familia prima: Mertens elemental basta (§2.1).
4. **Docs 118, 118b1 y 119b: "estratos q₀ ≥ 13 cerrados [T]" vale solo para 26,3 < z ≤ 1443,6.** La frontera q₀†(z) crece con z, así que (SLAB-V) y (SLAB-V′) deben pedirse **uniformes en un conjunto creciente de estratos**, no solo en {5, 7, 11} (§3.2).
5. **Residuo oculto (doc 114 §6, doc 116c §4): S(r) ≤ 1 quedó sin prueba en el régimen puente** al refutarse (R-FLOW). La implicación (BAND-δ) ⟹ (BAND-S) con factor 3/2 depende de ese enunciado, que pasa a [P] (§1.3).
6. **Doc 119a: el exponente de la condición m = 1 está mal.** Es r·q³ < L, no r·q² < L. El término principal correcto es W(u) = F(u) = 2e^γ/u. El veredicto [T módulo FL] se sostiene, con la hipótesis externa explícita (§4.1).
7. **Doc 119b: (SLAB-V′) no es "estrictamente menor" que (SLAB-V), sino equivalente salvo constante.** El motivo es que N(r) = |V^{(r)}| está acotado. Además, "Π_full es la densidad de V⁺(L)" es heurística [M], no [T] (§4.2).
8. **Doc 120a:**
   - la frase "degenerado del modelo ⟹ E = 0, δ = 1 − Π" es falsa en 32 de 172 casos;
   - el 58,2 % debe leerse **59,7 %**: un caso con ruido de coma flotante se contó como E ≠ 0;
   - la mediana 62,7 % no cambia (§5).

**Erratas menores** (§7):
- 5·7·11·13 = 5005, no 15015;
- 5⁹ = 1,95·10⁶ < 2·10⁷;
- E(z) y κ₀ calculados con p ≤ z, aunque el texto define 𝒫 como los primos < z;
- dos cifras de 118b2 J0 corresponden a otros intervalos o no se reproducen;
- el umbral real de (R-FLOW) es min(y) ≥ 73, y (DESP-SIGN) es [T] para z_A ≤ 73.

---

## §1. Doc 116

### §1.1 Doc 116a: (REM-PERS-small) cerrado

| Afirmación | Veredicto |
|---|---|
| **Lema A** [T]: q_min(x) no puede fallar, por paridad (vía el teorema R1 del doc 114) | **Confirmado** [T]. La prueba de R1 se reescribió: el lado fusionado es ≤ r·q²·P_m, así que el modelo fallaría en m. |
| **Lema B** [T]: ningún p \| r con p < q_min(x) falla, por tamaño | **Confirmado** [T]. r·x·q_min² < L ≤ p²·r·x obliga a q_min < p. |
| **Teorema:** la familia residual es vacía y ψ(x) ∈ Rem⁺ siempre | **Confirmado** [T]. La persistencia del doc 114 se aplica a todo x ∈ Rem⁻. |
| **Lema C** [T]: Rem no contiene primos | **Confirmado** [T]. Hay tres casos (R1 en q, R2 bajo q, R2 sobre q) y los tres están bien. |
| J0–J5 [M]: Rem⁻ = 39 + 36; V̂⁻ = 5568 + 498; 39 fallos (uno por x); 0 en primos ≤ q_min; 2423 + 221 pares; razón máxima 0,7581 / 0,4916; 0 residuales; 0 primos en Rem; ψ 39/39 + 36/36 | **Confirmado** al dígito, incluido el estrés z = 53 (185 instancias). También reproduzco las 4113 instancias y los conteos 6509/1049, A_min = 1,4273 y 1348 degenerados reales de los docs 113–114. |

**Precisión** [T]:
- Para min(y) = 71 se suma sobre q ≤ 67, y Mertens(67) = 160107799/160240080 = 0,99917 < 1.
- Por tanto R_A, R_R ≤ 1 es [T] para **min(y) ≤ 71**, y (DESP-SIGN) es [T] en todo marco con **z_A ≤ 73**, no solo ≤ 71.
- El residuo real empieza en min(y) ≥ 73, como ya usa el doc 117a (1a).

### §1.2 Doc 116b: (R-FLOW-71+) falso puntualmente

| Afirmación | Veredicto |
|---|---|
| Cotas (★), (★★) y (★★★) [T] | **Confirmadas** [T]. Son condiciones necesarias directas de r·y·q ∈ V⁺(L). En los 1640 pares hay 0 violaciones, y el máximo R_A = 0,6000 está en (B, z = 29, L = 2261685, r = 323, y = 23). |
| Configuración libre (r = 5, y = Y primo, ventana max(Y³, 125·Y·q_max) < L ≤ 5Y³) [T] | **Confirmada** [T]. Revisé el álgebra posición por posición. |
| J2: tabla del contraejemplo R_A en z = 331 (7 filas) | **Confirmada** con predicados propios. En cada fila: L = ⌊331⁵/s⌋ con s primo en [z, y_A); R_A exacto (Fraction) igual a la forma cerrada Σ_{5<q<Y} 1/(q−1). R_A(251) es **exactamente** la fracción impresa (50 primos), y el primer cruce está en Y = 251. |
| (R-FLOW-71+) falso como enunciado puntual (lado Add) [T] | **Confirmado** [T]. Que R_A ~ log log Y → ∞ se sigue de Mertens. Hay niveles de ℒ_A en toda ventana de razón 5 con Y ∈ (z^{5/9}, z), porque s'/s → 1. |
| J3: R_R(y) = 1,25620, "instancia en marco: z = 6052, s = 5, L = ⌊D/5⌋ ∈ ℒ_A" | **Refutado en la parte "∈ ℒ_A"**. La aritmética es correcta (y ∈ Rem⁺, R_R = Mertens(251), 52 primos), pero **s = 5 no es un primo de intercambio** (hace falta s ∈ [z, z^{10/3})) y **u_L = 4,815 > 4 = máx u_L de ℒ_A**. El ejemplo cae fuera del marco. |
| "R_R ≤ 1 también es falso como enunciado puntual en marco" | **Confirmado con un ejemplo propio** [T]. z = 4001, s = 8461 (primo, en [z, z^{10/3})), L = 121177241479732 (u_L = 3,910), ρ = 659, y = 661·653. Resultado: **R_R(y) = Mertens(73) = 1,02735 > 1**, con 19 primos admisibles. La familia necesita Y₂ > 71^{3/2} ≈ 598 y ρ·Y₁·Y₂³ ≤ z⁴, lo que exige z ≳ 3000; mi instancia usa z = 4001. |
| §2 y §0, y doc 116 §2: "R_A(y) > 1 solo puede ocurrir en la configuración libre: y primo, r primo pequeño…"; "el contraejemplo es la única configuración que las esquiva" | **Refutado** (ver abajo). |
| J4 [M] (máximos por cubeta a z = 101 y z = 71) | **No verificado** en esta tanda: depende de los niveles dirigidos y del presupuesto DFS del autor. |

**Contraejemplo a la "única configuración" (y compuesto).**
- **Instancia:** z = 28515, s = 162655277 (primo, en [z, z^{10/3})), L = 115904079453194 (u_L = 3,157), r = 5, **y = 28513·257·251**.
- **Resultado:** **R_A(y) = 1,00220 > 1**, con los 50 primos q ∈ (5, 251) admisibles.
- **Mecanismo:** q ocupa la posición fusionada 4 (par), y el atado del primo 5, en la posición 5, es 125·y·q < L. Ese atado queda vacío en rango cuando Y₁² > 25·Y₂·Y₃·q_max.
- **Consecuencia:** las cotas (★)–(★★★) son necesarias pero **no caracterizan** la obstrucción. El 117a (1c) ya trata los compuestos, así que la cadena no se rompe. Hay que corregir la frase.

### §1.3 Doc 116c: (BAND) cruda falsa asintóticamente

| Afirmación | Veredicto |
|---|---|
| Π(r) = Π_full·∏_{p\|r}(1 + 1/(p−2)) ≥ Π_full; Π_full(17) = 0,515625 > 1/2 > Π_full(19) = 0,4834; banda vacía a z = 17 | **Confirmado** [T + M]. |
| Obstrucción ω = 1: liminf ≥ B₁(u*) = log(min(1, u*/3)/(u* − u₀)) > 0 frente a Q_c → 0 [T] | **Confirmado** [T]. Revisé visibilidad (p³ < L), banda, Π → 0 y Mertens. B₁(5/3) = log(10/3) = 1,2040 y B₁(2) = log(4/3) = 0,2877. |
| J3a [M]: 0,2866 (Z = 10⁸, u = 2), 1,1917 (Z = 10¹², u = 5/3), cofactor 35: 0,0119 | **Confirmado** al dígito. |
| §2.4 [T]: cajas para u* ∈ [9/4, 4] (las del doc no están impresas) | **Confirmado** con cajas propias, verificadas en aritmética racional exacta. Para u* = 5/2 y 3 (ω = 2) y para 7/2 y 4 (ω = 3), todas cumplen las restricciones y tienen masa de Mertens > 0. La obstrucción cubre [5/3, 4]. |
| J3c [M]: Π_full·ln z → K ≈ 1,48; Q_c·ln z → K·M_∞ ≈ 1,72 | **Confirmado** (5 valores de z). Añado la forma cerrada **K = 3e^{−γ}·C₂/(3/4) = 1,48262** y **M_∞ ∈ [1,16181693; 1,16181705]**: producto hasta 10⁷ más cola rigurosa, intervalo dentro del del doc 115 y más estrecho. Con eso K·M_∞ = 1,72253 = C_loc (doc 109). |
| Tabla de ventana §3, J2 (0,1526 / 0,4408), J4 (p = 7 y 11), doc 114 K6 (0,6991 / 0,4408 / 0,3806 / 77) | **Confirmado** al dígito. |
| (BAND-δ) ⟹ (BAND-S) con factor 3/2 [T, módulo S(r) ≤ 1 del doc 114 §6] | **Confirmado el álgebra**, porque t⁺ = δ⁺(Π + S)/φ ≤ (3/2)δ⁺/φ. **Pero el módulo ya no es [T]** (ver abajo). |

**Residuo oculto: S(r) ≤ 1.**
- El doc 114 §6 da S(r) ≤ 1 como [T] si z_A ≤ 71 "y en general módulo (R-FLOW)". El 116b refutó (R-FLOW), así que en el régimen puente **S(r) ≤ 1 no tiene prueba: es [P]**.
- La ruta puntual también muere sobre el soporte completo. En el marco decisivo (z = 331, L = 79051268) hay **13 elementos y con R_{V⁽⁵⁾}(y) > 1**, con máximo 1,04522. Aun así S(5) = 0,350245 ≤ 1, pero solo por compensación.
- La cadena 113–116 consume S ≤ 1 en (BAND-δ) ⟹ (BAND-S). Hay que añadirlo al inventario de (H2): sería "(S≤1)" en la banda, o una cota S = O(1) allí.

---

## §2. Doc 117

### §2.1 Doc 117a: (FLOW-COMP) → (FC-2)

| Afirmación | Veredicto |
|---|---|
| (FC-2) ⟹ (FLOW-COMP) [T], con 1 − Mertens(37) = 2581/18480 = 0,139665 | **Confirmado** [T]. Mejora opcional: si min(y) ≤ 37, entonces q ≤ 31, y basta 1 − Mertens(31) = 0,167443. |
| "(FC-2) es estrictamente más débil que (FLOW-COMP)" (117 §0, 117a §0 y §4) | **Refutado como clasificación lógica.** Al implicarla, (FC-2) es una condición **suficiente y más fuerte**. En el nivel decisivo exige BAD ≤ 0,00449, mientras que (FLOW-COMP) exige BAD ≤ GOOD = 0,03040: el lado derecho de (FC-2) es el 15 % del de (FLOW-COMP). Su virtud es ser más sencilla (no usa flujos), no más débil. Además, (FLOW-COMP) equivale a "E_Add(r) ≤ 0 y E_Rem(r) ≤ 0", por la identidad GOOD − BAD = −E_Add, comprobada a 10⁻¹⁷. Es decir, es (DESP-SIGN) separado por lados. |
| (1a) Si R_A > 1, entonces min(y) ≥ 73 [T] | **Confirmado** [T]. |
| (1b) Ventana de la familia prima: [(L/r₀)^{1/3}, L^{1/3}) [T] | **Confirmado** [T]. |
| (1b) "Sin PNT el total no decae: el decaimiento es exactamente tan fuerte como PNT en la ventana cúbica"; bloqueador (i) "entrada de fuerza PNT" | **Refutado.** El número de primos de la ventana no hace falta contarlo con PNT. Por **Mertens elemental**, Σ_{a≤Y<r₀^{1/3}a} 1/Y = O(1/log a); por Chebyshev, π(x) ≪ x/log x. Así, masa mala = (ln r₀)(ln ln L − 2,81)/ln L·(1 + o(1)) → 0 para r₀ = 5. Lo medí a L = 10¹²…10²⁴: coincide con la predicción de Mertens al 0,5 %, con un máximo de ≈ 0,0355 cerca de L ~ 10²⁰ y caída posterior. **El bloqueador (i) desaparece para la familia prima.** De (FC-2) quedan solo la cola compuesta (1c) y la cota inferior de HEAD37(r). |
| K2 [M], nivel decisivo L = 79051268: E_Add = −0,029314, BAD = 1,086·10⁻³, GOOD = 0,030399, malos = los 13 primos 251…317, máx R_A = 1,04522 en 317, HEAD37 = 0,0321, cola = 0,0462, holgura 4,13, E_Rem = −0,000556, BAD_R = 0 | **Confirmado** con DFS propio exhaustivo de V̂⁽⁵⁾ ∪ V⁽⁵⁾ (11150 elementos). **Margen GOOD/BAD = 28,00.** |
| K2, niveles de control 46948396 y 113296524 (holguras 4,98 y 4,89) | **Confirmado** (GOOD/BAD = 35,0 y 32,2). |
| K1: BAD = 0 en las 4113; máx R_R = 0,25. K4: término malo Rem = 1,248·10⁻⁸ | **Confirmado.** |
| K3 [M] (estrés z = 101) | **No verificado** (niveles del autor). |

### §2.2 Doc 117b: (BAND-δ) → (DESP-SIGN) + (BAND-MOD)

| Afirmación | Veredicto |
|---|---|
| Reducción δ⁺ ≤ (Ŝ − Π)⁺ si E ≤ 0, acíclica [T] | **Confirmado** [T]. K1: en las 246 instancias de banda, E > 0 en 0 y violación máxima 0. |
| Identidad de doble conteo Σ_r T̂/φ = Σ_{d'} a·W_L; W_L(385) = 1,2702 / 0,2702 (λ = 9625) | **Confirmado** [T + M] (error ≤ 2,2·10⁻¹⁶). |
| "La suma con signo no tiene forma cerrada [T]" | **Reclasificar** como observación. La dependencia de W_L en L está demostrada; la inexistencia de forma cerrada no es un teorema. |
| Ŝ − Π > 0 en 246/246; por ω, 56/187/3 con masas 3,5374 / 0,3211 / 0,0001; sumas con signo +1,2722 y +1,5653 Q_c | **Confirmado** [M]. |
| Proposición: V̂^{(r)} = {1} ⟺ q₀³·r ≥ L [T] | **Confirmada** [T]. Además la comprobé en las **4113** instancias, no solo en las 246 de banda. |
| Tabla §4; 188/246 degenerados; máx (Ŝ − Π)⁺/Π no degenerado = 0,692 | **Confirmado** [M]. |
| (REG-Π) medido en 0,17–0,46·Q_c [M] | **Confirmado.** A z = 37: 0,1733 (u = 5/3) y 0,2695 (u = 2). A z = 53: 0,4609 y 0,2182. A z = 29: 0. Niveles L = ⌊z^u⌋. |

---

## §3. Doc 118

### §3.1 Doc 118a: (REG-Π) falsa

| Afirmación | Veredicto |
|---|---|
| Π/φ = Π_full·g; (REG-Π) ⟺ Σ g ≤ C·M_full | **Confirmado** [T]. |
| Rebanada ω = 1: limsup ≤ log(10/3)/M_∞ ≈ 1,036 [T] | **Confirmado** [T] (1,0363). A z = 10¹²: 1,1935 y 0,2876. |
| (REG-Π) falsa uniformemente: familia ℱ = {pq : x_p ∈ (11/20, 3/5)} a u_L → 2 [T] | **Confirmado** [T]. Revisé visibilidad (hace falta u_L > 9/5), banda (u_L ≤ 41/20), Π ≤ (8/5)·Π_full y masa ≥ log(12/11)·(log log z + O(1)). |
| J2: exterior 0,0850 → 0,0870 = log(12/11); interior 1,2262 → 2,1324 | **Confirmado** al dígito. |
| J1 (banda por ω) y J4 (z = 53: 77, 91 y 143 en banda, degenerados; g = 0,0505 = 0,0436·M_full; ventana estricta vacía) | **Confirmado.** |

### §3.2 Doc 118b1: subcapa degenerada

| Afirmación | Veredicto |
|---|---|
| **Teorema A** [T]: si L ≥ q₀⁹, ningún degenerado visible tiene ω ≤ 2 | **Confirmado** [T]. Que p₁(d)³ < L sea cierto para todo testigo da r < L^{2/3} ≤ L/q₀³. |
| Agudeza [M]: 188 degenerados con histogramas {1:1, 2:9}, {1:31, 2:47}, {2:100}, todos con L < q₀⁹ | **Confirmado.** Pero la frase "L ≤ 2·10⁷ < 5⁹" es aritméticamente falsa (5⁹ = 1953125). El marco B tiene niveles por encima de 5⁹; el hecho medido se cumple igual. |
| Estratos q₀ ≥ q₀† cerrados por telescopio de Euler [T] | **Confirmado** [T]: la telescopía, κ₀ y la cola ≤ 2E/P(Q). |
| **"Estratos q₀ ≥ 13 cerrados"** (118 §0 y §3); "residuo finito, q₀ ∈ {5, 7, 11}" (118b1 §7 y 119b) | **Reclasificar.** q₀†(z) = min{q : P(q) ≥ (ln z)³} vale **13 solo para 26,3 < z ≤ 1443,6**. Luego crece: 17 desde z = 1444 (también a z = 10⁶), 19 a 10¹², 23 a 10³⁰ y 29 a 10¹⁰⁰. Con Q = 13 fijo, la cota de Euler de la cola es 2E(z)/385 ≍ log z, es decir ≍ (log z)²·Q_c: **el argumento no cierra** en el puente (no está demostrado que la suma real diverja; está demostrado que esta prueba no la acota). El residuo (SLAB-V) es la familia de estratos 5 ≤ q₀ < q₀†(z): finita para cada z pero sin cota, y hace falta Σ_{q₀<q₀†(z)} C(q₀) = O(1). |
| J2 (0,1993 / 0,3554 / 0,0252, solo el estrato 5) y J4 (cota de Euler del estrato 5 = 3,25 / 3,72·Q_c; total sin visibilidad 0,716 / 1,405 / 0,025) | **Confirmado.** |
| J3 [M]: E/ln z = 0,6517 / 0,6267 / 0,6207 / 0,6169; κ₀ 0,8926 → 0,8834; colas 0,162 / 0,111 / 0,095 / 0,078 | **Confirmado, pero con p ≤ z.** Con 𝒫 = primos < z, como define el texto, salen 0,6134 / 0,6051 / 0,6039 / 0,6053 y κ₀ 0,8961 → 0,8837. El desfase es conservador: las cotas siguen valiendo. |

### §3.3 Doc 118b2: parte no degenerada

| Afirmación | Veredicto |
|---|---|
| Reconciliación: la masa divergente de ℱ vive en u_r ∈ (7/5, 29/20) y viene de x_q → 0 [T] | **Confirmado** [T]. J0(a): 0,0869 / 0,0868. |
| J0(b) [M]: fracción de la masa interna con u_r ≥ 1 = 0,718 / 0,781 / 0,808 | **No reproducido.** Con la definición escrita (x_q ≤ 1 − x_p, q < p) obtengo 0,790 / 0,833 / 0,854. Probé otras variantes y ninguna da 0,718. La conclusión (→ 1) se cumple en todas. |
| **(REG-Π-LOW) falsa para todo u₀* ∈ (0, 3/2)** [T]: hermana ω = 3 y escalera | **Confirmado** [T], en aritmética racional exacta. La subcaja impresa de la hermana es visible y da u_r < 0,94. Para la escalera k = 2 construí cajas propias con η = 1/200: Nathanson estricto en m = 1, 3, 5 y u_r < 0,2489 < 1/2. |
| J0(c): "masas 0,166 y 0,282 de la subcaja impresa" | **Las cifras no son de esa subcaja.** Son las masas de los intervalos (0,55; 0,65) → 0,1664 y (0,45; 0,60) → 0,2821. La subcaja impresa (0,56; 0,65) × (0,50; 0,55) da 0,1481 y 0,0935. Al [T] solo le hace falta la positividad, que se cumple. |
| Corrección al doc 117b §3.2 (la loseta baja tiene medida positiva; la medida cero es solo de la subcapa degenerada) | **Confirmada** [T]. |
| J1 (25 niveles: celda no degenerada con u_r < 1 vacía; tabla e instancias 0/1/1, 0/13/28, 0/10/42) y J2/J2b | **Confirmado** al dígito. |

---

## §4. Doc 119

### §4.1 Doc 119a: (BAND-MOD-HI) falsa [T módulo FL]

| Afirmación | Veredicto |
|---|---|
| Tabla K1 de §2 (banda total, capas ≥ 1,25 y ≥ 1,40, ω = 1, ω = 2 y mayorante, en A29, A37, B29 y estrés 53); máximo 1,2947 | **Confirmado**: las 24 cifras. |
| Perfil por rebanada de u_r ("1,07–1,24; 0,87–1,22; 0,69–1,07; 0,55–0,70") | **Confirmado el orden Θ(1)** y el suelo ≈ 1/2. Mis medias ponderadas por rebanada y marco son 1,04–1,06 / 0,86–1,04 / 0,58–1,07 / 0,53–0,70. Los rangos impresos no coinciden al dígito, porque el texto no dice si agrega por nivel o por marco. |
| **Mecanismo:** "la condición m = 1 exige r·q² < L… omite los primos q ∈ [z^{u_r/2}, z), de peso log(2/u_r)"; "W(u) > 1 para u < 2" | **Errata de exponente.** Con β = 2 la condición en m = 1 es **r·q³ < L** (p_m²·P_m con P₁ = p₁). Por tanto el rango excluido es **[z^{u_r/3}, z)**, de peso **log(3/u_r)**. Ejemplo: a z = 37, L = 37² y r = 7, V̂ = {1, 5}; con la regla q² entrarían también 11 y 13. El término principal es la función superior de la criba lineal: **W(u) = F(u) = 2e^γ/u** para 0 < u ≤ 3, y F(u) > 1 para todo u. Así que W − 1 > 0 en todo el rango, no solo para u < 2. Medición propia: Ŝ/Π a u = 1,25 vale 2,28 / 2,44 / 2,51 para z = 10³ / 10⁴ / 10⁵ (F = 2,85); a u = 1,45 vale 2,06 / 2,13 / 2,21 (F = 2,46). Es estable y crece hacia F(u). |
| (BAND-MOD-HI) falsa asintóticamente [T módulo FL] | **Confirmado [T condicional]**. La **hipótesis externa** es el teorema clásico de Rosser–Iwaniec para los pesos β = 2: Σ_{d∈V⁺} μ(d)g(d) = Π·(F(s) + o(1)), uniforme en el conjunto cribado 𝒫∖r, que cumple la condición de dimensión lineal. Es una cita (Iwaniec 1980; Friedlander–Iwaniec, *Opera de Cribro*), no algo probado en la serie. Con ella, la divergencia log log z sobre ℱ ∩ {u_r ≥ u₀*} vale para u₀* < 29/20. |
| La ruta L² muere igual | **Confirmado** (lectura): x⁺ ≤ x²/(2Π) + Π/2 devuelve Σ Π/φ. |
| (BAND-δ-HI) [P] | Es un enunciado abierto, bien formulado. |

### §4.2 Doc 119b: la reindexación divisor–cofactor, paso a paso

Muse pidió revisar este paso con lupa, por los errores de φ. Lo hice punto por punto:

| Paso | Veredicto |
|---|---|
| 1. Todo testigo d de un r degenerado tiene cofactor c = d/r < q₀³, porque d < L ≤ q₀³·r | **Correcto** [T]. 0 violaciones. |
| 2. d ≥ r ≥ L/q₀³: el testigo vive en la rebanada superior | **Correcto** [T]. |
| 3. φ(d) = φ(c)·φ(r), porque d es libre de cuadrados y (c, r) = 1 | **Correcto** [T]. 0 violaciones en todos los pares. **No hay error de φ.** |
| 4. Compatibilidad de estrato: (c, P(q₀)) = 1. Ojo: c **puede** contener a q₀; por ejemplo c = 5 en el estrato 5 | **Correcto** [T]. |
| 5. Σ_{r deg} 1/φ(r) ≤ W(q₀) | **Correcto** [T]. |
| 6. Intercambio W = Σ_c U_c | **Correcto** [T] (error 2,7·10⁻¹⁵). |
| 7. "ω(c) ≤ 3 porque 5·7·11·13 = 15015 > 1331" | **Errata aritmética:** 5·7·11·13 = **5005** (15015 = 3·5·7·11·13). La conclusión ω(c) ≤ 3 se mantiene; para q₀ = 11 es incluso ω(c) ≤ 2. |
| 8. P(q₀) = Σ_c U_c^{deg} = Σ_r N(r)/φ(r) | **Correcto** [T]. Observación nueva: **N(r) = \|V^{(r)}\|**, el tamaño del soporte desplazado **real**; verificado con 0 discrepancias. Las restricciones c < B y (c, P(q₀)) = 1 son automáticas para r degenerado. |
| 9. "(SLAB-V′) es estrictamente menor que (SLAB-V)" | **Reclasificar: son equivalentes salvo constante absoluta.** Se cumple verdadera ≤ P ≤ N_máx(q₀)·verdadera, con N_máx = 38 / 87 / 297 cofactores posibles para q₀ = 5 / 7 / 11. La ganancia real es quitar el cuantificador de visibilidad, no fuerza lógica. |
| 10. "El factor 1/ln z es Π_full, la densidad de V⁺(L)" [T/M] | **Reclasificar como heurística [M].** Reproduzco la predicción Π_full·A(5)·ln 125/Q_c = **3,52**, y la medida U_5 es 2,33–2,49. Pero que V⁺(L) tenga densidad Π_full no está demostrado. |
| Tabla §3: verdadera 0,1993 / 0,3554 / 0,0252; P 0,5979 / 0,7109 / 0,0252; W 16,26 / 16,09 / 5,46; N(r) ≤ 3 con medias 1,20 / 1,41 / 1,00; cofactores en W(5) 10/17/19 y en W(11) 4/16/29; W(7) ≤ 9,38; W(11) ≤ 2,04; W(5) a z = 17 = 10,35; parte c = 1 de W(5) 0,332 / 0,355 / 0,065 | **Confirmado** al dígito. |

**Veredicto sobre (SLAB-V) → (SLAB-V′): confirmado [T]**, como equivalencia y no como reducción estricta. Hereda la uniformidad en q₀ < q₀†(z) del §3.2.

---

## §5. Doc 120a: (BAND-δ-HI)

| Afirmación | Veredicto |
|---|---|
| K0: δ = (Ŝ − Π) + E y E = Σ_Add a − Σ_Rem a; E ≤ 0 en la capa | **Confirmado** [T + M]. |
| Tabla de capa (A29, A37, B29 y estrés z = 53 en los 4 niveles más próximos a u_L = 2: 2585, 2745, 2922, 3122) | **Confirmado**: las 16 cifras. |
| 67 no degeneradas con exceso; cobertura completa en 0/67; máx \|E\|/(Ŝ − Π)⁺ = 0,9930 | **Confirmado.** |
| **E = 0 exacto en el 58,2 % (39/67); 28 con E ≠ 0** | **Precisión numérica:** E = 0 **estructural** (Add = Rem = ∅) en **40/67 = 59,7 %**, con 27 casos E ≠ 0. La diferencia es **un** caso con ruido de coma flotante: (A, z = 37, L = 39557, r = 187), donde la ruta S = φ(r)C(r)/μ(r) del doc 109 da E = +1,1·10⁻¹⁶. |
| **Mediana del afeitado donde E ≠ 0 = 62,7 %** | **Confirmado: 0,6271**, igual con 27 que con 28 valores. |
| K2: 43 instancias (A29, estrés, A37); Rem vacío en todas; mediana \|E_Add\|/(Ŝ − Π) = 0,5581 | **Confirmado.** |
| §1: "Los degenerados de la capa (V̂ = {1}, **E = 0, δ = 1 − Π**)…" | **Falso en 32 de 172.** En esos casos el soporte real es V^{(r)} ⊋ {1}. Ejemplo: A29, L = 1360, r = 11, con \|V^{(r)}\| = 3, E = −0,4167 y δ = 0,0991, frente a 1 − Π = 0,5158. Confunde la degeneración del **modelo** (doc 117b) con la del **soporte real** (doc 113). Las sumas del doc usan el δ real y no cambian. Lo que cambia es la lectura "parte degenerada = territorio SLAB con δ = 1 − Π": allí δ⁺ ≤ 1 − Π solo si S ≤ 1 (§1.3). |
| (BAND-δ-HI) ⟺ (BAND-δ-HI′) bajo (i) E_Rem = 0 y (ii) E_Add ≤ 0 [T] | **Confirmado** [T] (álgebra puntual). |
| §3: el margen ×28 no transfiere; son demandas independientes sobre E | **Confirmado** (lectura): (FLOW-COMP) da signo; (BAND-δ-HI′) pide magnitud. |

---

## §6. Reclasificaciones

| Doc | Afirmación | Etiqueta en el doc | Etiqueta tras verificar |
|---|---|---|---|
| 116b §4 | Contraejemplo R_R "∈ ℒ_A" (s = 5) | [T] | **Falso como instancia de marco.** El enunciado R_R ≤ 1 falso en marco es [T] con mi instancia (z = 4001). |
| 116 §2, 116b §0 | "R_A > 1 solo en la configuración libre (y primo…)" | [T] | **Falso** (contraejemplo compuesto en ℒ_A). |
| 117 §0, 117a §0/§4 | "(FC-2) estrictamente más débil que (FLOW-COMP)" | — | **Al revés:** estrictamente más fuerte (suficiente). |
| 117a §1b | "Sin PNT el total no decae" (bloqueador i) | [T] | **Falso:** Mertens elemental basta. El bloqueador desaparece para la familia prima. |
| 114 §6, 116c §4 | S(r) ≤ 1 "en general módulo (R-FLOW)" | [T mód.] | **[P]** en el régimen puente, porque (R-FLOW) está refutado. |
| 118 §0/§3, 118b1 §7, 119b | "Estratos q₀ ≥ 13 cerrados; residuo finito {5, 7, 11}" | [T] | **[T] solo para z ≤ 1443.** En el puente, residuo uniforme en q₀ < q₀†(z), que es creciente. |
| 119a §1/§3, 119 §1 | Mecanismo r·q², log(2/u), "W > 1 para u < 2" | [T mód. FL] | **Errata de exponente.** Lo correcto es r·q³, log(3/u) y W = F(u) = 2e^γ/u. El veredicto se mantiene, condicional a la FL de Rosser–Iwaniec. |
| 119b §4 | "(SLAB-V′) estrictamente menor" | [T] | **Equivalente** a (SLAB-V) salvo constante (N(r) = \|V^{(r)}\| ≤ 38/87/297). |
| 119b §2/§4 | "Π_full es la densidad de V⁺(L)" | [T/M] | **[M]** (heurística con predicción 3,52 reproducida). |
| 117b §2.1 | "La suma con signo no tiene forma cerrada" | [T] | **Observación** (la dependencia en L, demostrada). |
| 120a §1 | "Degenerados: E = 0, δ = 1 − Π" | — | **Falso** en 32/172. |
| 114 §5, 116a §6 | (DESP-SIGN) [T] "si z_A ≤ 71"; residuo "min(y) ≥ 71" | [T] | **[T] si z_A ≤ 73**; el residuo empieza en min(y) ≥ 73. |

---

## §7. Erratas menores (para la autora; no toco sus archivos)

1. **119b §1:** 5·7·11·13 = 5005, no 15015.
2. **118b1 §2:** "L ≤ 2·10⁷ < 5⁹" es falso, porque 5⁹ = 1,95·10⁶. El marco B tiene niveles por encima.
3. **118b1 J3:** E(z), κ₀ y las colas se calcularon con p ≤ z, aunque el texto define 𝒫 como los primos < z. Es conservador.
4. **118b2 J0(b):** la fracción 0,718/0,781/0,808 no se reproduce; con la definición escrita sale 0,790/0,833/0,854.
5. **118b2 J0(c):** las masas 0,166 y 0,282 son de los intervalos (0,55; 0,65) y (0,45; 0,60), no de la subcaja impresa.
6. **119a §2:** los rangos del perfil por rebanada deben decir si agregan por nivel o por marco.
7. **120 §0 y 120a §1:** "58,2 %" debe ser 59,7 %, porque hay un caso de ruido de 1,1·10⁻¹⁶ contado como E ≠ 0.
8. **117a §2a (opcional):** el margen puede subir a 1 − Mertens(31) = 0,167443.

---

## §8. Inventario de (H2) tras esta verificación (PROPUESTA)

Tal como lo dejan los docs: **(FC-2), (BAND-δ-HI′), (BAND-LOW-NONDEG), (SLAB-V′)**. Con lo verificado propongo leerlo así:

- **(FC-2)** es más fuerte que lo necesario. El residuo exacto es (FLOW-COMP), es decir, E_Add ≤ 0 y E_Rem ≤ 0. Si se mantiene (FC-2), su bloqueador PNT desaparece para la familia prima, y quedan la cola compuesta (1c) y la cota inferior de HEAD37(r).
- **(BAND-δ-HI′):** sin cambios. Su hipótesis externa queda explícita: la FL de Rosser–Iwaniec, para la falsedad de (BAND-MOD-HI).
- **(BAND-LOW-NONDEG):** sin cambios.
- **(SLAB-V′):** equivalente a (SLAB-V), y **uniforme en 5 ≤ q₀ < q₀†(z)**, con Σ C(q₀) = O(1). No vale solo para {5, 7, 11}.
- **Nuevo residuo, hasta ahora oculto: (S≤1)**, es decir, S(r) ≤ 1 (o S = O(1)) en la banda, en el régimen puente. Lo consume (BAND-δ) ⟹ (BAND-S) y también la lectura degenerada de 120a. En ventana se cumple por compensación (S(5) = 0,35 en el marco decisivo), no punto a punto.

Siguen fuera, sin cambios: (DIST-𝒟) en la línea de Lumina Chat, y G0. **(DIST-𝒟) + (H2) ⇒ UNIF-U-COND** sigue siendo condicional.

**Goldbach fuerte sigue abierta.**

---

## §9. Reproducibilidad

```
cd numerico
python verificacion_116_120_claude_code.py > resultados/20261128_verificacion_116_120_claude_code.txt
```

Solo stdlib (Python 3.14), 4 s. Contenido del script:
- **Barrido principal:** 4 marcos × 60 niveles = 4113 instancias, con flujos Add/Rem, fallos fusionados, banda y capa.
- **Estrés:** z = 53 (docs 114 K7 y 120a).
- **DFS exhaustivo a z = 331** (117a K2).
- **Contraejemplos con predicados propios** y búsqueda de niveles en ℒ_A con Miller–Rabin.
- **Sumas de Mertens segmentadas hasta 10⁸.**
- **Cajas de regiones en aritmética racional exacta.**

Cada chequeo contrasta una afirmación del texto: PASS = confirmada; FAIL = hallazgo, interpretado aquí.

*Goldbach fuerte sigue abierta.*
