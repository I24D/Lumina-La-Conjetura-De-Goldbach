# ESTADO del proyecto (registro vivo)

Cada entrada lleva: etiqueta, enunciado, fuente y quién la añadió. Etiquetas definidas en `matematicas/01_marco_barreras_y_rutas.md`, sección 0.

## Objetivo principal
- **[H] Goldbach fuerte:** todo par n ≥ 4 es p + q. **ABIERTO.**

## Resultados conocidos relevantes
| Etiqueta | Enunciado | Fuente | Añadido por |
|---|---|---|---|
| [T] | Goldbach débil: todo impar n ≥ 7 es suma de tres primos | Helfgott 2013 | Claude, 2026-11-23 |
| [T] | Chen: todo par suficientemente grande es p + P₂ | Chen 1973 | Claude, 2026-11-23 |
| [T] | E(N) < N^{0,72} | Pintz, arXiv:1804.09084 (abstract comprobado) | Claude, 2026-11-23 |
| [T]? | E(N) = O(N^{7/10}), con constante inefectiva | Zhao, arXiv:2511.05631 (abstract comprobado; sin arbitrar) | Claude, 2026-11-23 |
| [C: GRH-Dirichlet] | E(N) ≪ N^{1/2+ε} (Goldston: N^{1/2}(log N)^c) | Hardy–Littlewood 1924; Goldston 1992 (vía el survey arXiv:2607.27282) | Claude, 2026-11-23 |
| [C: semiplano común libre de ceros] | E(N) ≪ N^{7/12+ε} | Pintz–Althöfer, arXiv:2610.15559 | Claude, 2026-11-23 |
| [T]? / [C: GRH] | Todo par grande es p + q + 6 potencias de 2 (incondicional); 5 bajo GRH | Huixi Li y Zhixin Liu, arXiv:2610.15784 (sin arbitrar) | Claude, 2026-11-23 |
| **[T]** | **Goldbach binaria en F_q[T]** para p impar y q > 685090 p² | Sawin–Shusterman, Annals 2022, Remark 1.2 (comprobado literalmente) | Claude, 2026-11-23 |
| [T-Lean]* | RH y GRH para caracteres de Dirichlet | Alpöge, 20-11-2026. Certificado por Lean FRO/Mathlib el 21-11-2026 | Claude, 2026-11-23 |
| [N] | Goldbach verificada hasta 4·10¹⁸ (sigue siendo el récord reconocido) | Oliveira e Silva, Herzog y Pardi 2014 | Claude, 2026-11-23 |
| [T] | Reducción 1.1: basta R(n) > 2,9 √n log n | elemental, `matematicas/01` §1 *(verificar la constante de Rosser–Schoenfeld)* | Claude, 2026-11-23 |

\* La verificación formal de Alpöge está certificada. **No** hay prueba legible ni aceptación matemática general. Si se acepta, se vuelven incondicionales la fila [C: GRH-Dirichlet] y el "5 bajo GRH" de Linnik.

## Trabajos analizados
| Trabajo | Veredicto (propio, [P]) | Archivo |
|---|---|---|
| Li & Liu, *A Conditional Proof of the Goldbach Conjecture* (13-11-2026, sin arbitrar) | **Reformulación, no reducción.** Por sus ecuaciones (2.28) y (2.30), las "hipótesis espectrales" S y R equivalen a cotas de ‖T_N‖ al menos tan fuertes como la conclusión. GRH sola da ‖T_N‖ ≪ (log N)²; hace falta o(1). Lo útil: la identidad exacta (2.39) y que basta acotar λ⁺_max(T_N), es decir, una cota de un solo lado | `matematicas/02_analisis_Li_Liu_2026.md` |

## Resultados numéricos propios [N]
| Fecha | Qué | Resultado | Archivo |
|---|---|---|---|
| 2026-11-23 | `selftest` | OK (criba, Miller–Rabin y p_min contra fuerza bruta hasta 20 000; r(n) hasta 1000) | — |
| 2026-11-23 | `verify` hasta 10⁸ | 49 999 999 pares, todos OK. Máximo p_min = 1093 (n = 60 119 912). Los 29 récords coinciden con OEIS A025018/A025019 | `numerico/resultados/20261123_203349_verify_100000000.json` |
| 2026-11-23 | `verify` hasta 10¹⁰ | 4 999 999 999 pares, todos OK, ningún contraejemplo. Máximo p_min = 1877 (n = 7 473 202 036). Los 36 récords coinciden con OEIS (`contrastar`). Tardó 1218 s. Usó la versión del script anterior a `contrastar`; la lógica de `verify` es la misma | `numerico/resultados/20261123_205419_verify_10000000000.json` |
| 2026-11-23 | `reps` hasta 10⁷ | Ningún r(n) = 0. r(n)/Hardy–Littlewood entre 0,947 y 1,025 por tramos (muestreado) | `numerico/resultados/20261123_203438_reps_10000000.json` |
| 2026-11-23 | T-LL2: `li_liu_check.py` hasta 10⁷ | El error con pesos de Li–Liu cumple max \|err\|/δ_N entre 0,14 y 0,23, sin crecimiento. Coherente con la conclusión de su Teorema 1.2; no contrasta su hipótesis | `numerico/resultados/20261123_204322_li_liu_check_10000000.json` |

## En curso
*(nada)*

## Pendientes y bloqueos
- ~~**Lean 4 + Mathlib sin instalar.** Bloquea la línea L2 (formalización).~~ **Resuelto el 2026-11-23** (Claude, a petición de Dal): ver la entrada de Claude al final y `matematicas/06` §11.
- **MCP `lumina_shared_memory` caído** (2026-11-23). Quedan pendientes de registrar como propuesta: la nota RH/Alpöge, `matematicas/01` y `matematicas/02`.
- Tareas abiertas T-LL1 a T-LL4 (ver `matematicas/02` §4).
- Lecturas pendientes:
  - survey de Bhowmik–Grimmelt (arXiv:2607.27282);
  - arXiv:2610.26756 (cuartos momentos de sumas exponenciales de primos);
  - Sawin–Shusterman, §§2–3 (mecanismo de paridad);
  - Sutherland, https://math.mit.edu/~drew/grh/
- Comprobar las referencias que siguen marcadas *(verificar)* en `matematicas/01`.

## Registro de afirmaciones refutadas [FALSO]
*(vacío)*


## Aporte de Codex — 2026-11-23 (entrada posterior; preserva las anteriores)

**Procedencia:** Codex, modelo exacto no expuesto; lectura primaria, derivaciones escritas y comprobaciones finitas. Estado global **[P]**, sin revisión independiente ni Lean local. Esta entrada amplía el registro previo: no lo sobrescribe ni eleva su estado por acuerdo entre agentes.

| Etiqueta | Resultado / alcance | Evidencia |
|---|---|---|
| [P] | Mecanismo de paridad, hipótesis de tamaños y comparación con la interfaz formal de Alpöge; no se afirma transferencia a los enteros | `matematicas/03_sawin_shusterman_alpoge_puente.md` |
| [P] | Identidad exacta R=S−T; términos no coprimos acotados por potencias propias; evaluación uniforme del coeficiente de Euler | `matematicas/04_reduccion_mobius_goldbach.md`, §§2–4 |
| [P; derivación con resultados clásicos incondicionales] | Para D=⌊N^(1/3)⌋, R(N)=N𝔖(N)−T_N^prim(D)+O_A(N/(log N)^A). No acota T_N^prim | `04`, §4.3, Prop. 4.6 |
| [C: GRH-Dirichlet; P] | Variante con error O_ν(N^(5/6+ν)); 5/6 NO es una cota del conjunto excepcional | `04`, (4.10) |
| [C: HμΛ; P] | La hipótesis adicional U_N(D)=o(N/log N) bastaría para Goldbach para pares suficientemente grandes; **HμΛ no demostrada** | `04`, §5 |
| [FALSO] | Atajo universal basado solo en soporte impar, densidad positiva y raíz cuadrada en los modos no triviales: contraejemplo denso por caracteres cuadráticos. No es un contraejemplo sobre primos | `03`, §4; amplía el registro de refutaciones que antes estaba vacío |
| [N] | 1193 identidades simbólicas, 3000 coeficientes racionales, 1092 casos de Pellet, 81 casos de fibras, 178657 verificaciones de tamaños y 16 ejemplos densos; PASS en 0,5244 s | `numerico/resultados/20261124_021251_764252_puente_paridad.json` (nombre fechado en UTC) |

**Lecturas realizadas / límites:** revisados el mecanismo algebraico de §§2–3 de Sawin–Shusterman, enunciados de distribución, Cor. 6.1 y la descomposición de §6.3; no se certificó de nuevo toda la prueba geométrica. Leída la página primaria de Sutherland y las declaraciones de Main/Statements/Challenge; no se ejecutaron verificadores ni se descargó la exportación grande. El survey y el trabajo de cuartos momentos previamente pendientes no se dan por auditados.

**Pendientes concretos:**
- Demostrar la estimación primitiva con cuantificadores uniformes, o una cota unilateral suficiente para la cola. La noticia de GRH no la acredita.
- Revisar de manera independiente las pruebas de 03–04 y formalizar sus pasos centrales. `lean` y `lake` no están disponibles en PATH.
- Obtener un umbral efectivo y una comprobación certificada que lo alcance antes de concluir Goldbach completa. La versión usada de Bombieri–Vinogradov lleva constantes no efectivas.
- Registrar este aporte como propuesta en memoria compartida: **PENDIENTE**, herramientas `lumina_shared_memory`/`lumina_propose` no expuestas. Respaldo con procedencia y hashes: `notas/2026-11-23_CODEX_aporte_memoria_pendiente.json`. No se eludió el contrato por REST/service_role.

**Ejecuciones:** terminadas; ningún cálculo del verificador queda en segundo plano. No se iniciaron modelos, builds, instalaciones ni cambios de servicios. **Goldbach no está demostrada por este aporte.**

## Aporte de LUMINA/OpenClaw - 2026-11-23

| Etiqueta | Resultado / alcance | Evidencia |
|---|---|---|
| [P] | Reindexacion exacta de la cola por `n=N-cd`, que permite enumerar divisores de `N-n` | `matematicas/05_diagnostico_cola_y_objetivo_bilineal.md`, (5.1) |
| [N] | Implementacion literal y reindexada coinciden para todos los pares `N<=300` | `numerico/cola_mobius_scan.py` |
| [N] | `U log(N)/N` vale 3.2367, 3.4482, 3.2752 y 3.1691 para `N=10^3,10^4,10^5,10^6`; HmuLambda no recibe apoyo numerico y queda despriorizada, no refutada | `numerico/resultados/20261124_023525_cola_mobius_1000000.json` |
| [P] | Una cota uniforme por cajas bilineales firmadas (HB) implica `T_N^prim=o(N)` conservando cancelacion entre fibras | `matematicas/05_diagnostico_cola_y_objetivo_bilineal.md`, seccion 4 |

**Estado honesto:** Goldbach sigue abierta. El avance elimina un blanco numericamente sospechoso y lo reemplaza por una desigualdad bilineal precisa. El siguiente cuello es una descomposicion Tipo I/II valida en todos los rangos y uniforme para cada N par.

## Aporte de Claude Code — 2026-11-23 (tarea HB; entrada posterior, preserva las anteriores)

**Procedencia:** Claude Code, modelo `claude-opus-5-5`. Derivación escrita y diagnóstico numérico propio. Estado **[P]**, sin revisión independiente.

| Etiqueta | Resultado / alcance | Evidencia |
|---|---|---|
| [P] | Lema 6.1: (HB) se cumple en las cajas con K ≤ √N/(log N)^B, por BV directo en c. Equivale a subir D de 04 a √N/(log N)^B | `matematicas/06_descomposicion_bilineal.md` §5 y R-2 |
| [P] | Lema 6.2: con Vaughan y UV ≤ K^{1/2}/L^B, los términos de Tipo I y los degenerados son ≪ N/L^A (BV para μ en d). En R2, B_N = B₃ (Tipo II) + O(N/L^A) | `06` §6 |
| **ABIERTO** | B₃ (Tipo II) en todas las cajas con K > √N/L^B. Con C ≥ N^δ basta [C: EH en clase fija]; con C < N^δ es el núcleo de paridad. Ningún resultado publicado lo controla (BFI, Zhang–Polymath y Maynard requieren módulos factorizables o lisos) | `06` §7 |
| [T] elemental | Lema 6.0: en la caja c = 1, la parte con d primo es −Σ_{p~K} log p·Λ(N−p), un bloque diádico de Goldbach. (HB) en las cajas con C pequeño equivale en dificultad a Goldbach por bloques | `06` §4 y R-3 |
| [P-rev + T] | Mejor cota: −(3+o(1))N𝔖(N) ≤ T^prim ≤ (1+o(1))N𝔖(N). Goldbach necesita ≤ (1−κ)N𝔖(N) | `06` §8 |
| [N] | Comprobación cruzada independiente: T^prim/N = 0,0642495 en N = 10⁶, idéntico al de OpenClaw. ΣP/N y ΣQ/N valen ≈ ∓10,9 en N = 10⁷ con T/N = 0,086. La fila c = 1 reproduce −R(N)/N ≈ −1,761 frente a 𝔖 = 1,7604. max\|B\|/N por caja: 0,063 → 0,034 → 0,018 | `numerico/cajas_bilineales.py`; `numerico/resultados/20261123_215414_cajas_bilineales.json` y `20261123_215627_cajas_bilineales.json` |
| infraestructura | **Lean 4 + Mathlib instalado.** elan 4.4.0; proyecto en `C:\Users\dal_n\lean\goldbach_lumina` (fuera de OneDrive), toolchain v4.37.0-rc1. La prueba de humo compila: enunciado de Goldbach, 100 = 3+97 e identidad Λ = −Σ μ log (Mathlib) | `06` §11 |

**Informes sobre 04 y 05 (ninguno es un error):**
- R-1: confirmación de (HB) ⇒ T = o(N), con la precisión de truncar en d > D la caja que contiene a D.
- R-2: la Prop. 4.6 vale con D = √N/L^B.
- R-3: la parte compuesta de la cola es ≈ N𝔖(N) incondicionalmente, así que "T = o(N)" es la asintótica de Hardy–Littlewood.

**Recomendación:** no perseguir (HB) caja por caja en R2. El blanco honesto es la cota unilateral sobre C < N^δ (mejorar la constante 1 a 1−κ), y requiere información sensible a la paridad. Formalizar ya en Lean la Prop. 4.1, la Prop. 4.2 y el Lema 6.0.

**Memoria compartida:** registro como propuesta **PENDIENTE**. El MCP `lumina_shared_memory` no está disponible en esta sesión de Claude. **Goldbach no está demostrada por este aporte.**

## Aporte de Claude Code — 2026-11-23, segunda ronda (entrada posterior, preserva las anteriores)

**Procedencia:** Claude Code, `claude-opus-5-5`. Experimentos FFT, derivación escrita y Lean. Detalle en `matematicas/07_sesgo_mobius_goldbach_y_casi_todo_N.md`.

| Etiqueta | Resultado / alcance | Evidencia |
|---|---|---|
| [P] (entradas clásicas) | **Teorema 7.1:** (HB) se cumple en todas las cajas para todo N par de (X/2, X] salvo O(X/(log X)^B), por Parseval sobre N y Davenport 1937. Recupera E(X) ≪ X/log^B en el marco de cajas. Lección: (HB) es fácil en promedio sobre N; la dificultad está en cada N | `07` §1 |
| [N] | **Fenómeno:** la correlación Möbius–Goldbach M0(N) = Σ μ(n)Λ(N−n) tiene varianza sub-aleatoria (≈ 0,90) y un **sesgo negativo cuando N tiene factores primos pequeños**: media de Z = −0,219 si 3 \| N y corr(Z, 𝔖(N)) = −0,49 en N ≤ 10⁶; −0,094 y −0,26 en N ≤ 10⁷ | `numerico/mobius_goldbach_aleatoriedad.py`; `resultados/20261123_220731_*.json` |
| [H + N] | **Fórmula explicativa (arcos mayores):** P_Q(N) = Σ_{q≤Q} μ(q)/φ(q)·Σ_{n<N} μ(n)c_q(N−n); el término q = 1 es la función de Mertens. Con Q = 6 desaparece el sesgo de 3; con Q = 30, todos los sesgos de clase y la correlación con 𝔖 (−0,49 → −0,045 en 10⁶; −0,26 → −0,017 en 10⁷). R² = 0,13 / 0,29 / 0,45 con Q = 30 / 210 / 1000 | `numerico/arcos_mayores_mobius.py`; `resultados/20261123_2209*_*.json`, `…221209…`, `…221311…` |
| [T-Lean] | **Verificado por el kernel** (solo propext, Classical.choice y Quot.sound): Prop. 4.1 de 04 (R = S − T para todo N y D), coprimalidad (N−cd, N) = 1 ⇒ (c, N) = 1, núcleo de la Prop. 4.2, Lema 6.0 y R-3 | `lean/README.md`; canónico en `C:\Users\dal_n\lean\goldbach_lumina\GoldbachLumina\Reduccion.lean` |

**Estado honesto:** el sesgo explicado es de arcos mayores (tamaño √N; bajo GRH es inofensivo). **La parte de arcos menores del núcleo de paridad sigue abierta** (06, F1 y F2). No se reivindica novedad bibliográfica sin buscarla. Memoria compartida: **PENDIENTE** (MCP no disponible). **Goldbach no está demostrada.**


## Codex - revision 06/07 en curso (2026-11-23)

Procedencia: Codex; modelo exacto no expuesto; estado [P]/[N]. Lectura del encargo completada. Se reproduciran secuencialmente los dos scripts a 10^6, y despues se ampliara el diagnostico hasta 10^7, sin modelos ni calculos pesados simultaneos. Se prepara matematicas/08_revision_codex_06_07.md. La entrada final indicara el cierre y los resultados.

## Cierre de Codex — revisión 06/07, 2026-11-23 ET

**Procedencia:** Codex; identificador exacto de modelo no expuesto. Método: lectura del encargo y de 04–07, derivaciones escritas, búsqueda de fuentes primarias, experimentos secuenciales hasta 10^7 y pruebas nuevas en Lean. **[P-rev] NO independiente**, al compartir usuario de Windows. Se preservaron los documentos y fuentes de otros agentes; esta entrada posterior matiza expresamente afirmaciones anteriores.

Informe principal: [08_revision_codex_06_07.md](matematicas/08_revision_codex_06_07.md).

| Tarea / estado | Resultado comprobado y alcance |
|---|---|
| T1 — revisión entregada [P-rev] | Lema 6.1 válido en el rango previsto, explicando K≥N^δ y 2K dentro de BV; se corrige el par de exponentes de Rankin. Lema 6.2: (q,c)=1 correcto, pero la extensión de BV para μ a clases no reducidas queda como **HUECO**, no teorema acreditado. R-2 se mantiene. |
| T1 — error explícito R-3 | Falta **−Tbad** en la identidad exacta de la parte compuesta. N=32,D=3 da Tbad=log 2·log 36≠0. Su asintótica N𝔖+o(N) sobrevive. HB+ implica una cota cuantitativa de representaciones; no se ha probado equivalencia con mera existencia. |
| T1 — Teorema 7.1 [P-rev] | Parseval requiere ≤ al restringir N≤X; Davenport uniforme es aplicable. Para cajas móviles faltaba un argumento maximal. 08 §4 aporta la reparación escrita con descomposición diádica de prefijos y control de su energía. Sigue siendo para **casi todo** N y no se formaliza aquí en Lean. |
| T2 — búsqueda entregada | Λ_Q=Σμ(q)c_q/φ(q) es un modelo aditivo conocido, explícito en Bhowmik–Grimmelt, arXiv:2607.27282v2 (2.6). Su convolución P_Q con μ no se presenta como fórmula nueva. No se localizó publicación del sesgo estadístico exacto −0,22/−0,49; búsqueda limitada, sin reivindicar novedad. |
| T3 — PASS reproducción [N] | Ambos comandos pedidos a 10^6 salieron 0 y reproducen 07. Error FFT/directo 3,65·10⁻¹⁵. Se detecta muestreo por pasos restringido a múltiplos de 5 en la correlación original; −0,49 no representa por sí solo a todos los pares. |
| T3 — PASS ampliación [N] | X=10^5,10^6,10^7; residuo Ramanujan Q=30 con medias −0,0005285,+0,0000536,−0,0000632 y asimetrías −0,17217,−0,08260,−0,02637. No emerge un signo estable para HB+. P_Q no es automáticamente una proyección continua de arcos mayores. |
| T3 — cola completa [N] | Todas las c, **25 N por escala**, no todos los N. Medias de (T−T_30)/(N𝔖): −0,00534,+0,00211,+0,00693. La identidad Tprim=S−R−Tbad se contrasta con bucle literal; usar R para evaluarla no prueba HB+ independientemente. |
| T4 — PASS [T-Lean] | `Reduccion2.lean`: reindexación (c,d), desigualdad explícita (4.12) y primera cota de Prop. 4.2 \|Tbad\|≤K_N. Fuente final compila con exit 0, sin errores ni advertencias. Seis `#print axioms`: únicamente propext, Classical.choice, Quot.sound. No se formalizan la segunda cota de Prop. 4.2, la estimación asintótica de K_N o las cotas analíticas de 06/07. |
| Memoria — PENDIENTE | Sin herramientas `lumina_shared_memory`/`lumina_propose` expuestas. Tres propuestas locales: aportes de Claude 06/07 con sus observaciones de revisión, y el aporte de Codex. Sin ID remoto ni persistencia MCP afirmada. |

**Archivos nuevos:** `matematicas/08_revision_codex_06_07.md`; `numerico/codex_residuos_06_07.py`; `lean/Reduccion2.lean`; `lean/2026-11-23_codex_Reduccion2_verificacion.json`; `notas/2026-11-23_CODEX_revision_06_07_memoria_pendiente.json`; `notas/2026-11-23_CODEX_revision_06_07_cierre.json`.

**JSON de experimentos:** `numerico/resultados/20261123_223206_mobius_goldbach_aleatoriedad_1000000.json`, `20261123_223326_arcos_mayores_mobius_1000000.json` y `20261124_033802_732515_codex_residuos_06_07.json`. Los dos primeros nombres usan hora local; el tercero, UTC. El último contiene hashes, semilla 20261123, parámetros, estadísticas por clases y bloques y muestra de cola completa.

**Lean canónico:** `C:\Users\dal_n\lean\goldbach_lumina\GoldbachLumina\Reduccion2.lean`, SHA-256 `2568abd7e5ccd6a125cc0fbef40c9dc80b66cd5593ea0129764f9cb5e77096e6`; copia idéntica UTF-8 sin BOM en OneDrive. `Reduccion.lean` permanece preservado. Se añaden resultados a `lean/README.md` y un enlace al README general.

**Pendiente matemático:** justificar la entrada no reducida de Lema 6.2; revisión independiente de la reparación maximal; control unilateral uniforme de Tipo II sensible a la paridad, con constantes y cuantificadores; umbral efectivo y tramo finito. Los resultados numéricos, el teorema para casi todos y las identidades formales no resuelven esa obligación. **Goldbach sigue abierta.**

## Verificacion LUMINA/OpenClaw - 2026-11-24

| Etiqueta | Resultado | Evidencia |
|---|---|---|
| [T-Lean, reproducido] | `Reduccion2.lean` compila con exit 0; seis declaraciones sin `sorryAx`; hash de copia y fuente canonica identico | `lean/2026-11-23_codex_Reduccion2_verificacion.json` |
| [P-rev] | Confirmados los tres errores documentados en 08 y el hueco de clases no reducidas; ninguno cierra el Tipo II | `matematicas/08_revision_codex_06_07.md` |
| [P] | Identidad exacta (9.1) reduce una clase no reducida de `mu` a clases reducidas mediante `g=(r,q)` y un parametro `g`-liso `k` | `matematicas/09_reduccion_clases_no_reducidas_mu.md` |
| [ABIERTO] | Falta demostrar la estimacion agregada WNR; aun con ella, el termino Tipo II `B_3` permanece como barrera principal | `matematicas/09_reduccion_clases_no_reducidas_mu.md`, secciones 3-5 |

**Veredicto:** la revision de Codex es tecnicamente util y su artefacto Lean fue reproducido. La investigacion avanzo al convertir el hueco vago de clases arbitrarias en una desigualdad explicita, pero no existe todavia demostracion de Goldbach.

## Continuacion autonoma LUMINA/OpenClaw - 2026-11-24

| Etiqueta | Resultado | Evidencia |
|---|---|---|
| [N] | 3000 verificaciones exactas de (9.1); `E_all/E_reduced`=1.1026,1.0668,1.0186 para x=10^4,10^5,10^6 | `matematicas/10_diagnostico_WNR_y_prioridad_tipo_II.md` |
| [P]+[N] | Descomposicion `T=sum(-1)^j A_j`; en N=10^6 las masas par/impar son 17.934M/17.870M y su diferencia es solo 64,249 | `matematicas/11_capas_paridad_y_funcion_generadora.md` |
| [P] | Formula exacta `H_m(z)=z(1+z)^(omega(m)-1)log(rad(m))`; version primitiva (12.4) y soporte ralo en z=-1 | `matematicas/12_polynomial_divisores_y_soporte_ralo.md` |
| [N] | Prop. 12.1 contrastada en 10000 casos aleatorios con z entero; PASS | registro de ejecucion LUMINA 2026-11-24 |

**Nuevo frente:** expandir (12.6) por divisores cuadrados libres `a|rad(N)` y buscar una forma positiva/unilateral que conserve el termino `a=1`. El objetivo Goldbach permanece activo; no hay demostracion completa.

## RH/GRH y transferencia a Goldbach - 2026-11-24

| Etiqueta | Resultado | Evidencia |
|---|---|---|
| [T, formal] | RH y GRH para caracteres de Dirichlet: enunciados canonicos aceptados por Lean y nanoda; recompilaciones independientes reportadas | `matematicas/13_riemann_alpoge_transferencia_goldbach.md`, seccion 1 |
| [P] | Identidad Riesz puntual `R_1(N)=R(N)/2`; la suavizacion conserva positividad y el coeficiente N | `13`, seccion 3 |
| [P, descarte] | GRH lineal no controla automaticamente el Tipo II aditivo; separar sectores par/impar pierde la cancelacion | `13`, secciones 4-5 |
| [PROGRAMA] | Buscar un complejo Z/2 y una desigualdad de pendientes que controle el supertrazo de Mobius mediante ganancias locales positivas | `13`, seccion 6 |

**Direccion activa:** estudiar el motor de positividad de Alpoge, no solo usar RH como caja negra. La transferencia necesita un puente multiplicativo-aditivo nuevo; Goldbach sigue abierta.

## Colapso de la expansion por divisores de N - 2026-11-24

| Etiqueta | Resultado | Evidencia |
|---|---|---|
| [P] | En (12.8), `G_a(N)=0` exactamente si `omega(a)>=2`; para `a=ell|N`, todo termino no nulo satisface `N/ell=p^k+ell^(j-1)` y es `O((log N)^3)` | `matematicas/14_colapso_divisores_de_N.md`, seccion 1 |
| [P] | La suma completa colapsa a `-G_1` mas correcciones sobre los primos divisores de N; ademas `G_1=R(N)+O(omega(N)log^3 N)` | `14`, ecuaciones (14.1)-(14.2) y corolario |
| [N] | Ambos lados de (14.1) coinciden en 549 pares; error maximo `1.14e-13`; soporte prohibido exactamente cero | `numerico/resultados/20261124_colapso_divisores_N.json` |
| [FALSO como frente] | La inclusion-exclusion grande sobre divisores impares propuesta en 12/13 no existe: los grados >=2 se anulan antes. Un puente nuevo debe actuar directamente sobre `G_1` | `14`, seccion 3 |

**Direccion activa corregida:** estudiar kernels aditivos Fejer/Riesz o una formula espectral unilateral que controle el coeficiente individual `G_1(N)`. Goldbach sigue abierta.

## Deficit unilateral puntual - 2026-11-24 (LUMINA/OpenClaw)

| Etiqueta | Resultado | Evidencia |
|---|---|---|
| [P] | Identificacion: el supertrazo Z/2 de 11 (ec. 11.1) ES la cola T_N^prim(D) de 04; la desigualdad de pendientes de 13 s6 y la H_eps de 04 s5 son la misma condicion, escrita puntual y normalizada como ho(N)=T_prim(D)/(N*G(N)) <= 1-eps | matematicas/15_deficit_unilateral_puntual.md, s1-s2 |
| [N] | Escaneo literal en 2299 N pares (todos <=2000 exhaustivos + muestras hasta 10^6, semilla 20261124): sup rho = 0.4201 (N=30916), sin deriva ascendente; banda alta (300 muestras) sup = 0.2535, mediana -0.013. Veredicto pre-registrado (umbral sup<=0.5): H_eps plausible con margen amplio (eps >= 0.58 medido) | 
umerico/deficit_unilateral.py; esultados/20261124_224341_...json y 20261124_224653_..._dense_...json |
| [N] | Estructura de los peores N: series singulares casi minimas (2^a*p o dos primos medianos): 30916=2^2*59*131, 31718=2*15859, 30658=2*15329, 9224=2^3*1153... Los N "dificiles" de Goldbach concentran el peor ratio unilateral | 15, s4 |
| [N; honesto] | La cadena (4.15) punto a punto: gap = 1-r-rho pequeno en N grandes tipicos (|gap| <~ 0.1), pero grande co-localizado con rho muy negativo (N=443190: rho=-0.561, r~1.000, gap~0.561); T_bad solo esta acotada por K_N (Prop. 4.2), no es pequena a priori | 15, s4 |

**Estado honesto:** no hay demostracion. La calibracion dice cuanto margen real tiene la cota unilateral y donde vive el peor caso (familias delgadas de serie singular minima). Barrido de familias EJECUTADO (2026-11-24, `deficit_familias.py`, 15 458 filas exactas): VEREDICTO PRE-REGISTRADO CONFIRMADO — las familias delgadas capturan el supremo (2p: 0.4816 @ 2218; 4p: 0.4798 @ 34324; ventana exhaustiva (25000,40000]: 0.4841 @ 32728=2^3*4091; todas > 0.4201 global). Nuevo sup rho medido 0.4841 → ε ≥ 0.52. Hallazgo estructural: medianas de rho en 2p/4p = +0.23 (cola primitiva cronicamente sobredimensionada en familias de G minima, frente a mediana global −0.013); minimo −0.878 @ 30030 (primorial). El sup familiar no crece con N dentro del rango. Siguiente pre-registro: 2p y 4p densos hasta 10^6 para acotar si sup rho puede crecer mas alla de ~0.5.

## Bucle LUMINA/OpenClaw - 2026-11-24 18:07 ET - revision adversarial de 17

- **[P-rev]** Confirmada contra Lichtman Def. 1.4 la obstruccion C7: `mu(q) log q` no es triplemente bien factorizable (ni bien factorizable), por su soporte en primos grandes; esto basta para cerrar la aplicacion directa de Prop. 3.1 a B3.
- **[P-rev, correccion de alcance]** C10 no demuestra que el termino principal de B3 sea `O_A(N/log^A N)`: el termino real contiene un peso movil `W(d)`, y la suma no ponderada de `mu(d) log d` no lo controla sin un lema de variacion/suavizado.
- Evidencia: `matematicas/18_revision_adversarial_doc17_Lichtman.md`. Goldbach sigue abierta.
## Bucle LUMINA/OpenClaw - 2026-11-24 19:17 ET - extension pre-registrada 2p/4p

- **[N]** 600 evaluaciones estratificadas de `rho(N)` en las familias `N=2p,4p`, hasta `10^6`, con hipotesis y refutacion registradas antes de ejecutar.
- Maximo `rho=0.364285` en `N=111292`; derivas del supremo entre primera y ultima banda: `-0.0902` (2p), `-0.1013` (4p). El criterio (`rho>=0.75` o deriva `>0.15`) no fue refutado.
- Evidencia: `matematicas/20_extension_prereg_familias_2p_4p.md`. Es evidencia finita, no cota uniforme; B3 sigue abierto y Goldbach no esta demostrada.
## Bucle LUMINA/OpenClaw - 2026-11-24 19:37 ET - Lean Hito 3

- **[T-Lean]** `lambdaFlat U V n` se expandio exactamente como suma finita sobre `n=m*d` y `m=a*b`, conservando los cortes `U,V` explicitos.
- Verificacion: `lake env lean GoldbachLumina/VaughanFlatExpansion.lean` exit 0; `lake build` verde (18299 jobs); axiomas unicamente `propext`, `Classical.choice`, `Quot.sound`; 0 `sorry`/axiomas propios.
- Evidencia: `matematicas/21_lean_expansion_tipoII_lambdaFlat.md`, `lean/VaughanFlatExpansion.lean`. No aporta cotas analiticas ni resuelve `B3`; Goldbach sigue abierta.
- **Memoria compartida:** propuesta pendiente; MCP `lumina_shared_memory` no expuesto en esta ejecucion. Respaldo local en el doc 21 y su recibo.
