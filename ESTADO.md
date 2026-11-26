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


## Bucle LUMINA/OpenClaw - 2026-11-24 20:12 ET - Lean Hito 4

- **[T-Lean]** Cuatro lemas puntuales verifican exactamente los soportes de las diferencias truncadas de `lambdaFlat`: se anulan para `a<=U`, `d<=V` y coinciden con `Lambda(a)`, `mu(d)` por encima de los cortes.
- Verificacion: `lake env lean GoldbachLumina/VaughanFlatSupport.lean` exit 0; `lake build` verde (18301 jobs); axiomas unicamente `propext`, `Classical.choice`, `Quot.sound`; 0 `sorry`/axiomas propios.
- Evidencia: `matematicas/22_lean_soporte_tipoII.md`, `lean/VaughanFlatSupport.lean`. Falta filtrar la suma completa y definir cajas diadicas; no hay cota analitica de `B3` y Goldbach sigue abierta.
- **Memoria compartida:** propuesta pendiente; MCP no expuesto en esta ejecucion.


## Bucle LUMINA/OpenClaw - 2026-11-24 20:37 ET - puente exacto `lambdaFlat` / `B3`

- **[P-rev]** Auditada la notacion de los docs 19--22 contra el doc 06: el `B3` clasico corresponde a `lambdaFlat V U`, no a `lambdaFlat U V`; los nombres de los dos cortes estan intercambiados entre ambas interfaces.
- **[T] elemental** Para `n>=1`, `lambdaFlat V U(n) = sum_{ar=n, a>V} Lambda(a) nu_U(r)` por expansion y reagrupacion finita; la variable clasica `r` absorbe el factor `zeta` y la variable interna de `mu`.
- Evidencia: `matematicas/23_puente_lambdaFlat_a_B3.md`. No hay cota analitica nueva; `B3` y Goldbach siguen abiertos.

## Bucle LUMINA/OpenClaw - 2026-11-25 04:39 ET - Lean Hito 5

- **[T-Lean]** Se verifico e integro el puente exacto `lambdaFlat V U` / `B3`: cuatro lemas compilan y el corte `nuAbove U r = 0` para `r <= U` queda probado por kernel.
- Verificacion: `lake env lean GoldbachLumina/VaughanB3Bridge.lean` y `lake build`, ambos exit 0; axiomas unicamente `propext`, `Classical.choice`, `Quot.sound`; 0 `sorry` propios.
- Evidencia: `matematicas/24_lean_puente_lambdaFlat_B3.md`, `lean/VaughanB3Bridge.lean`. Es una identidad finita, no una cota de Tipo II; Goldbach sigue abierta.

- **Memoria compartida:** propuesta pendiente; MCP `lumina_shared_memory` / `lumina_propose` no expuesto en esta ejecucion. Sin ID remoto ni persistencia afirmada.

## Bucle LUMINA/OpenClaw - 2026-11-25 05:15 ET - Lean Hito 6

- **[T-Lean]** Se definio una caja diadica finita con cortes `C<c<=2C`, `K<d<=2K`, `cd<N` y se probo su igualdad exacta con la expansion clasica de `B3`, para pesos exteriores arbitrarios.
- Verificacion: `lake env lean GoldbachLumina/VaughanB3DyadicBox.lean` y `lake build`, ambos exit 0; axiomas unicamente `propext`, `Classical.choice`, `Quot.sound`; 0 `sorry` propios.
- Evidencia: `matematicas/25_lean_caja_diadica_B3.md`, `lean/VaughanB3DyadicBox.lean`. Es empaquetado algebraico finito, no una cota de Tipo II; Goldbach sigue abierta.
- **Memoria compartida:** propuesta pendiente; MCP no expuesto, sin ID remoto ni persistencia afirmada.

## Bucle LUMINA/OpenClaw - 2026-11-25 05:39 ET - Lean Hito 7

- **[T-Lean]** Se especializo la caja diadica formal de `B3` a los pesos literales del doc 06: `1_{(c,N)=1}` y `mu(d) log d`.
- Verificacion: `lake env lean GoldbachLumina/VaughanB3Doc06Weights.lean` y `lake build`, ambos exit 0; axiomas unicamente `propext`, `Classical.choice`, `Quot.sound`; 0 `sorry` propios.
- Evidencia: `matematicas/26_lean_pesos_doc06_B3.md`, `lean/VaughanB3Doc06Weights.lean`. Es una identidad finita, no una cota de Tipo II; Goldbach sigue abierta.
- **Memoria compartida:** propuesta pendiente; MCP no expuesto, sin ID remoto ni persistencia afirmada.

## Bucle LUMINA/OpenClaw - 2026-11-25 06:10 ET - Lean Hito 8

- **[T-Lean]** Se formalizo el soporte del peso exterior mu(d) log d: se anula en d=0,1 y fuera de los libres de cuadrados; peso no nulo implica d libre de cuadrados y d>=2.
- Verificacion: lake env lean GoldbachLumina/VaughanB3WeightSupport.lean y lake build, ambos exit 0; axiomas unicamente propext, Classical.choice, Quot.sound; 0 sorry propios.
- Evidencia: matematicas/27_lean_soporte_peso_moebius_B3.md, lean/VaughanB3WeightSupport.lean. Es soporte algebraico finito, no una cota de Tipo II; Goldbach sigue abierta.
- **Memoria compartida:** propuesta pendiente; MCP no expuesto, sin ID remoto ni persistencia afirmada.
## Direccion LUMINA/Muse - 2026-11-25 06:32 ET - ataque analitico, no empaquetado

- Orden de Dal: lo publicado no es descubrir. Buscar avances no publicados o la respuesta.
- Balance: docs 18-27 son formalizacion Lean de identidades conocidas + un experimento finito (doc 20). No hay cota analitica nueva. B3 sigue abierto.
- Prioridad: (1) C10 / peso movil W(d): formular M_B3 exacto y acotar su variacion. (2) Uniones C_k+C_{k+1} con a=-N. (3) Variante Lema 4 Pintz-Althofer con sup en sigma > 5/7. (4) Idea nueva para B3 con pre-registro.
- Lean pasa a apoyo: solo formalizar lo que sirva al frente analitico.
- Detalle: notas/2026-11-25_DIRECCION_ataque_analitico.md

## Bucle LUMINA/OpenClaw - 2026-11-25 06:36 ET - termino principal B3 y coprimalidad

- **[T] elemental** El termino principal movil se expandio exactamente por inclusion-exclusion: `M=sum_{r,e} mu(e)log(re)W_r(re)/(phi(r)phi(e))`, con `r` libre de cuadrados y `(r,e)=1`.
- **[T] elemental** Para `r` fijo, la variacion total debida al movimiento de los extremos se acota por dos veces la masa absoluta de `F=Lambda*nu_U` sobre multiplos de `r`; la movilidad geometrica no es el obstaculo restante.
- Evidencia: `matematicas/28_descomposicion_coprimalidad_main_B3.md`; 500 verificaciones racionales exactas PASS. Quedan la uniformidad Mertens en `e` y la cola `r` grande; no hay cota de B3 ni demostracion de Goldbach.
- **Memoria compartida:** propuesta pendiente; MCP canonico no expuesto, sin ID remoto ni persistencia afirmada.

## Bucle LUMINA/OpenClaw - 2026-11-25 07:06 ET - cola de divisores grandes del main B3

- **[T] elemental** Para cualquier `0<delta<1/2`, la parte `r>N^delta` de (28.3) es `N^{1-delta+o(1)}`: el intercambio de sumas da `|M_{>R}| << log K ||F||_1 max_n sum_{r|n,r>R}1/phi(r)`, y para `F=Lambda_{>V}*nu_U` se usa `||F||_1<<N log^3 N`.
- Evidencia: `matematicas/29_cola_divisores_grandes_main_B3.md`. El cuello del termino principal queda reducido a `r<=N^delta`, con sumas de `e` de longitud al menos `N^(1/2-delta)/log^B N`; falta auditar la uniformidad Mertens en los primos de `rN`.
- No se controla el error de dispersion de `B3`; Goldbach sigue abierta. Memoria compartida pendiente: MCP no expuesto, sin persistencia remota afirmada.

## Bucle LUMINA/OpenClaw - 2026-11-25 07:36 ET - serie de Dirichlet del cuello Mertens

- **[T] elemental** Para `q` par, `A_q(s)=sum_{(e,q)=1} mu(e)/(phi(e)e^s)=zeta(1+s)^(-1)H_q(s)` tiene cero simple en `s=0`, con `H_q(0)=2 C_2 prod_{p|q,p>2}(p-1)/(p-2)`.
- **[T] elemental** La dependencia local es `O(log log q)`: excluir los primos de `rN` no produce por si solo una perdida exponencial en `omega(rN)`; la constante se cancela en diferencias diadicas.
- Evidencia: `matematicas/30_serie_Dirichlet_Mertens_main_B3.md`. Falta demostrar la cota uniforme de colas (30.8), incluyendo derivadas de factores locales; el error de dispersion de `B3` sigue abierto y Goldbach no esta demostrada.
- **Memoria compartida:** propuesta pendiente; MCP no expuesto, sin persistencia remota afirmada.

## Bucle LUMINA/OpenClaw - 2026-11-25 08:06 ET - derivadas locales del cuello Mertens

- **[T]** Para `q` par, el multiplicador local `G_q(s)` y cada derivada de orden fijo son `(log log(3q))^{O_j(1)}` uniformemente en la franja estrecha usada por un contorno de Vinogradov--Korobov.
- Evidencia: `matematicas/31_derivadas_factores_locales_main_B3.md`. Esto descarta una perdida `exp(omega(rN))` fuera de `s=0`; falta aun Perron truncado, segmentos horizontales y retirada del suavizado para demostrar (30.8).
- El error de dispersion de `B3` sigue abierto y Goldbach no esta demostrada. Memoria compartida pendiente: MCP no expuesto, sin persistencia remota afirmada.

## Bucle LUMINA/OpenClaw - 2026-11-25 08:36 ET - Perron suavizado del cuello Mertens

- **[T] derivado de teoremas clasicos** Para `q<=x^L`, las sumas exponencialmente suavizadas de `mu(e)/phi(e)` y `mu(e)log(e)/phi(e)`, con `(e,q)=1`, tienen ahorro arbitrario en logaritmos uniformemente en `q`; la segunda se centra exactamente en `-H_q(0)`.
- Evidencia: `matematicas/32_perron_suavizado_uniforme_main_B3.md`. El argumento usa Mellin, la region libre de ceros de Vinogradov--Korobov y el control local del doc 31.
- Falta insertar una particion suave directamente en (28.3) o justificar el desuavizado; el error de dispersion de `B3` sigue abierto y Goldbach no esta demostrada. Memoria compartida pendiente: MCP no expuesto, sin persistencia remota afirmada.

## Bucle LUMINA/OpenClaw - 2026-11-25 09:06 ET - desuavizado y termino principal B3

- **[T] derivado de teoremas clasicos** Una transicion suave de ancho relativo `log(N)^(-B)` se retira por masa absoluta; esto prueba las colas bruscas uniformes (30.8) para `mu(e)/phi(e)` y `mu(e)log(e)/phi(e)`.
- **[T]** Combinado con la cola `r>N^delta` del doc 29, el termino principal movil `M` de (28.1) es `O_A(N/log(N)^A)` para todo `A`.
- Evidencia: `matematicas/33_desuavizado_uniforme_main_B3.md`. El error de dispersion de `B3` no esta controlado; la barrera de paridad y Goldbach siguen abiertas.
- **Memoria compartida:** propuesta pendiente; MCP no expuesto, sin persistencia remota afirmada.

## Bucle LUMINA/OpenClaw - 2026-11-25 09:36 ET - error de dispersion B3

- **[T] elemental** Se escribio exactamente `B3=M+E`, con `E=sum mu(d)log(d) Delta_d`, y se obtuvo el split bilineal literal `E=sum rho(r)rho(e)log(re)Delta_{re}`, `rho=mu/2^omega`, sin multiplicidad espuria.
- **[T]/[P]** Se aislaron una norma bilineal suficiente y su version unilateral con los signos fijos de Mobius; no se probo ninguna ganancia sobre Cauchy ni una cota nueva de dispersion.
- Evidencia: `matematicas/34_error_dispersion_split_mobius.md`; `python numerico/verifica_dispersion_split.py` PASS en 500 configuraciones racionales. El error Tipo II y Goldbach siguen abiertos.
- **Memoria compartida:** propuesta pendiente; MCP no expuesto, sin persistencia remota afirmada.

## Bucle LUMINA/OpenClaw - 2026-11-25 10:14 ET - conservacion de la diagonal del split Mobius

- **[T] elemental** Para cualquier particion en rectangulos de las factorizaciones coprimas `d=re`, los coeficientes normalizados tienen el mismo signo y suma absoluta 1; la suma de las normas diagonales por bloque es al menos la diagonal original.
- Consecuencia: factorizar con `rho=mu/2^omega` y aplicar Cauchy bloque a bloque no produce por si solo ganancia; hace falta controlar productos distintos o una forma unilateral global.
- Evidencia: `matematicas/35_conservacion_diagonal_split_mobius.md`; `python numerico/verifica_energia_split_mobius.py` PASS en 12160 enteros libres de cuadrados. No hay cota nueva de `B3`; Goldbach sigue abierta.
- **Memoria compartida:** propuesta pendiente; MCP no expuesto, sin persistencia remota afirmada.

## Bucle LUMINA/OpenClaw - 2026-11-25 10:44 ET - productos distintos del split Mobius

- **[T] elemental** La forma global factorizada colapsa fibra a fibra exactamente al modulo original: `sum_{re=d} rho(r)rho(e)=mu(d)`.
- **[P-rev]** Los productos distintos no tienen signo favorable automatico y pueden superar la diagonal por un factor `m-1`; descartarlos o dominarlos por Cauchy no cierra `B3` sin una correlacion aritmetica nueva entre los residuos `Delta_d(F)`.
- Evidencia: `matematicas/36_productos_distintos_split_mobius.md`; `python numerico/verifica_productos_distintos_global.py` PASS. Se abandona esta version del split y el siguiente frente sera `C_k+C_{k+1}`. Goldbach sigue abierta.
- **Memoria compartida:** propuesta pendiente; MCP no expuesto, sin persistencia remota afirmada.

## Bucle LUMINA/OpenClaw - 2026-11-25 11:14 ET - auditoria de capas consecutivas

- **[T]/[P-rev]** La sustitucion `a=-N` en los teoremas de shifted primes mira `p>N` y `p-N>0`, mientras las capas de Goldbach requieren `p<N` y `N-p>0`; la cota inferior de Timofeev para `p+2` no se transfiere directamente a `C_2+C_3`.
- Evidencia: `matematicas/37_obstruccion_dominio_uniones_capas_shift.md`, contrastado con Ford Cor. 1. El siguiente frente es formular el tamiz directamente sobre `A_N={N-p:p<N}` y auditar su nivel uniforme.
- No hay cota inferior nueva de una capa, y Goldbach sigue abierta. Memoria compartida pendiente: MCP no expuesto, sin persistencia remota afirmada.

## Bucle LUMINA/OpenClaw - 2026-11-25 11:44 ET - nivel de criba reflejada

- **[T]** Para `A_N={N-p:p<N}`, Bombieri--Vinogradov da nivel uniforme `D=N^(1/2)/log(N)^B` en los modulos coprimos con `N`, con densidad local `1/phi(d)`; los primos que dividen `N` solo aportan excepciones puntuales.
- **[T]/[P-rev]** El tamiz lineal desnudo rompe exactamente en `s=2`: forzar `P_3` pide `z>N^(1/4)` y entonces `s=log D/log z<2`, donde la funcion inferior es cero; en el rango positivo solo se obtiene `P_4`.
- Evidencia: `matematicas/38_nivel_distribucion_criba_reflejada_s2.md`. Hace falta auditar pesos de Chen/Richert y switching; no hay cota inferior nueva de capas y Goldbach sigue abierta.
- **Memoria compartida:** propuesta pendiente; MCP no expuesto, sin persistencia remota afirmada.

## Bucle LUMINA/OpenClaw - 2026-11-25 12:14 ET - perdida de nivel al fijar q

- **[T]** Al fijar un primo q|N-p, la sucesion residual tiene nivel D/q; si q=N^delta, su nivel relativo es (1/2-delta)/(1-delta)<1/2.
- **[T]/[P-rev]** La criba lineal residual da s<2 para P2 y s<=2 para P3; el paso "aplicar Chen/BV" del doc 16 ocultaba una suma de switching bilineal sum_{r,ell} a_r b_ell Lambda(N-qrell), uniforme en q.
- Evidencia: matematicas/39_auditoria_fijar_q_switching.md. No hay cota inferior nueva de capas; B3 y Goldbach siguen abiertos. Memoria compartida pendiente: MCP no expuesto.

## Bucle LUMINA/OpenClaw - 2026-11-25 12:44 ET - geometria de cajas de switching

- **[T] elemental** Si `q=N^delta` y `r<=ell`, las cajas `q r ell~N` fuera de Bombieri--Vinogradov ocupan exactamente `1/2-delta<rho<=(1-delta)/2`, una cunha de anchura `delta/2`.
- **[T]/[P-rev]** Intercambiar `r` y `ell` no recupera BV: tanto `qr` como `qell` superan `N^(1/2)` en toda la region critica. Una adaptacion de Chen/Timofeev debe cubrir la cunha completa, no solo la caja balanceada.
- Evidencia: `matematicas/40_geometria_cajas_switching.md`; verificacion racional exacta PASS en 3540 cajas. No hay estimacion nueva de switching, cota inferior de capas ni demostracion de Goldbach. Memoria compartida pendiente: MCP no expuesto.


## Bucle LUMINA/OpenClaw - 2026-11-25 13:14 ET - peso acoplado del switching

- **[T] elemental** El peso real de una extraccion de Buchstab contiene `K(r,ell)=1_{P^-(ell)>r}`; en primos intercalados su matriz es triangular de rango arbitrario, por lo que no puede absorberse en `a_r b_ell` ni en un numero fijo de productos independientes.
- **[P-rev]** La forma (39.9) era una abreviatura geometrica, no aun un lema aplicable: una adaptacion de Chen debe conservar el nucleo de orden o introducir otra variable de switching y controlar la forma trilineal resultante.
- Evidencia: `matematicas/41_auditoria_peso_switching.md`; `python numerico/verifica_peso_switching.py` PASS para rangos `h=1..25`. No hay cota nueva de switching, capas o `B3`; Goldbach sigue abierta. Memoria compartida pendiente: MCP no expuesto.

## Bucle LUMINA/OpenClaw - 2026-11-25 13:44 ET - extraccion Buchstab y cofactor primo

- **[T]** La primera extraccion de Buchstab es exacta: `ell=s t`, `s=P^-(ell)>r`, `P^-(t)>=s`.
- **[T]/[P-rev]** En la cunha critica, si `q<=N^(1/4)` el cofactor `ell` debe ser primo; para `q>N^(1/4)` la rama compuesta solo ocupa `1/2-delta<rho<=(1-delta)/3`, de anchura `(4delta-1)/6`. Esto restringe la obstruccion abstracta de rango del doc 41.
- Evidencia: `matematicas/42_cuna_switching_fuerza_cofactor_primo.md`; verificacion exacta de identidad y geometria PASS. Falta estimar la rama bilineal sobre dos primos; no hay cota inferior de capas y Goldbach sigue abierta. Memoria compartida pendiente: MCP no expuesto.

## Bucle LUMINA/OpenClaw - 2026-11-25 14:14 ET - norma de la rama prima de switching

- **[T]** La norma bilineal absoluta de la rama prima contiene necesariamente toda su masa principal positiva; pedir cancelacion de esa norma antes de centrar es una meta imposible.
- **[T]/[P-rev]** Cauchy expande el problema en correlaciones de dos formas primas; Bombieri--Vinogradov no controla esos productos, y fijar cualquier variable deja modulos supra-raiz.
- Evidencia: `matematicas/43_norma_rama_prima_switching.md`; `python numerico/verifica_norma_rama_prima.py` PASS en 500 configuraciones. Falta una mayorante local con serie singular y el balance real de pesos Chen/Richert; no hay cota inferior de capas y Goldbach sigue abierta. Memoria compartida pendiente: MCP no expuesto.
## Bucle LUMINA/OpenClaw - 2026-11-25 14:44 ET - serie singular local de la rama prima

- **[T] elemental** Para `L_1(x)=x` y `L_2(x)=N-qrx`, la densidad local es `nu_p=1` si `p|Nqr` y `nu_p=2` en otro caso; por tanto `S(N;q,r)=2 C_2 prod_{p|Nqr,p>2}(p-1)/(p-2)`.
- **[T]/[P-rev]** La dependencia en `q,r` es positiva, separable y `1+O(1/q+1/r)`: una mayorante local aislada no aporta cancelacion; falta insertar los pesos reales de Chen/Richert y probar que su balance de constantes sea positivo.
- Evidencia: `matematicas/44_serie_singular_local_rama_prima.md`; `python numerico/verifica_serie_singular_switching.py` PASS en 400 ternas y 10000 factores locales. No hay cota inferior de capas y Goldbach sigue abierta. Memoria compartida pendiente: MCP no expuesto.

## Bucle LUMINA/OpenClaw - 2026-11-25 15:16 ET - serie singular bajo pesos firmados

- **[T] elemental** La extraccion de `S_N` de una caja con coeficientes `c(q,r)` tiene error acotado por `S_N O(1/Q+1/R) sum |c(q,r)|`, no por el tamano de la suma firmada.
- **[T]/[P-rev]** Dos coeficientes `+1,-1` dan suma base cero pero correccion local no nula: la casi constancia puntual de la serie singular no justifica comparar constantes Chen/Richert tras cancelacion sin control L1/margen.
- Evidencia: `matematicas/45_perturbacion_serie_singular_pesos_firmados.md`; `python numerico/verifica_perturbacion_serie_singular.py` PASS en 1000 matrices racionales. No hay cota inferior de capas y Goldbach sigue abierta. Memoria compartida pendiente: MCP no expuesto.

## Bucle LUMINA/OpenClaw - 2026-11-25 15:44 ET - cancelacion marginal de la serie singular

- **[T] elemental** La perturbacion de primer orden de `S(N;q,r)=S_N(1+a_q)(1+b_r)` depende solo de los margenes del peso firmado `c(q,r)`.
- **[T]/[P]** Si todos los margenes en `q` y `r` se anulan, los terminos de orden `1/Q+1/R` desaparecen exactamente y queda una cota de segundo orden `S_N sum|c|/((Q-2)(R-2))`; esto da una condicion auditable de diseno, no construye los pesos.
- Evidencia: `matematicas/46_cancelacion_marginal_serie_singular.md`; `python numerico/verifica_cancelacion_marginal.py` PASS en 1000 matrices racionales. Siguiente: extraer del primer peso explicito de Chen/Richert sus margenes y medir el defecto. No hay cota inferior de capas y Goldbach sigue abierta.
- **Memoria compartida:** propuesta pendiente; MCP no expuesto, sin persistencia remota afirmada.


## Bucle LUMINA/OpenClaw - 2026-11-25 16:14 ET - margenes del peso de orden de switching

- **[T] elemental** Toda matriz de pesos se descompone de forma canonica en una interaccion con margenes nulos y una proyeccion aditiva que contiene exactamente toda la perturbacion de primer orden de la serie singular.
- **[T]/[P-rev]** Para el nucleo explicito de Buchstab `K_ij=1_{i<j}`, los margenes son `h-i` y `j-1`, y su masa marginal total coincide con `||K||_1`: el criterio del doc 46 falla a escala principal; centrar el peso no permite descartar la proyeccion aditiva.
- Evidencia: `matematicas/47_descomposicion_marginal_peso_switching.md`; `python numerico/verifica_descomposicion_marginal.py` PASS en 1000 matrices racionales y `h=2..60`. Siguiente: calcular margenes de dos ramas consecutivas de Buchstab con sus signos. No hay cota inferior de capas y Goldbach sigue abierta.
- **Memoria compartida:** propuesta pendiente; MCP no expuesto, sin persistencia remota afirmada.

## Bucle LUMINA/OpenClaw - 2026-11-25 16:44 ET - dos ramas consecutivas de Buchstab

- **[T] elemental** En el modelo ordenado unitario, las dos primeras incidencias firmadas dan C_ij=1_{i<j}(2+i-j), con margenes cerrados R_i=m(3-m)/2, S_j=j(3-j)/2.
- **[T]/[P-rev]** La norma L1 de los margenes sigue asintoticamente a escala completa (sum|R_i|/||C||_1 -> 1): anadir una segunda rama no cancela automaticamente la perturbacion local de primer orden. Esto audita solo el modelo unitario, no un peso Chen/Richert real.
- Evidencia: matematicas/48_dos_ramas_buchstab_margenes.md; python numerico/verifica_dos_ramas_buchstab.py PASS para h=2..300. Siguiente: extraer coeficientes y cortes de una identidad publicada de Chen/Richert. No hay cota inferior de capas y Goldbach sigue abierta.
- **Memoria compartida:** propuesta pendiente; MCP no expuesto, sin persistencia remota afirmada.

## Bucle LUMINA/OpenClaw - 2026-11-25 17:18 ET - cancelacion marginal en Lean

- **[T-Lean]** Se formalizo para matrices finitas sobre un anillo conmutativo que margenes nulos cancelan exactamente los terminos constante y lineales de `(1+a_i)(1+b_j)`, dejando solo la interaccion mixta.
- Evidencia: `matematicas/49_cancelacion_marginal_Lean.md`; `lake build` verde (18311 jobs) y `#print axioms` solo `propext`, `Classical.choice`, `Quot.sound`. No se construyen pesos Chen/Richert ni una cota de switching; Goldbach sigue abierta.
- **Memoria compartida:** propuesta pendiente; MCP no expuesto, sin persistencia remota afirmada.


## Bucle LUMINA/OpenClaw - 2026-11-25 17:44 ET - proyeccion marginal en Lean

- **[T-Lean]** Se formalizo sobre cualquier anillo conmutativo la version sin denominadores de la descomposicion del doc 47: la interaccion escalada tiene ambos margenes exactamente nulos y reconstruye la matriz original junto con su proyeccion aditiva.
- Evidencia: `matematicas/50_proyeccion_marginal_Lean.md`; `lake build` verde (18313 jobs) y `#print axioms` solo `propext`, `Classical.choice`, `Quot.sound`. No se construyen pesos Chen/Richert ni una cota de switching; Goldbach sigue abierta.
- **Memoria compartida:** propuesta pendiente; MCP no expuesto, sin persistencia remota afirmada.

## Bucle LUMINA/OpenClaw - 2026-11-25 18:14 ET - auditoria del peso publicado de switching

- **[T-lit]/[T]** En el marco publicado de Matomaki--Zuniga Alterman, la rama penalizada usa `|mu(m)| W^+(p_1m)` (y una correccion positiva por bloque), de modo que sus coeficientes internos son no negativos.
- **[P-rev]** La matriz firmada y la cancelacion marginal de los docs 46--48 no representan ese peso Kuhn/Richert real: la resta ocurre externamente como `Sigma_1-Sigma_2`, tras acotar cada rama positiva. Se desprioriza seguir centrando matrices.
- Evidencia: `matematicas/51_revision_peso_publicado_switching.md`; `python numerico/verifica_peso_publicado_switching.py` PASS en 4000 pruebas racionales. Siguiente: instanciar la sucesion switched `b_{P,ell}` para `A_N={N-p}` y auditar A4--A5; no hay cota inferior nueva de capas y Goldbach sigue abierta.
- **Memoria compartida:** propuesta pendiente; MCP no expuesto, sin persistencia remota afirmada.


## LUMINA/Muse - 2026-11-25 18:34 ET - margenes del peso real de Chen contra el criterio del doc 46

- **[T-lit]/[T]** Identidad publicada de Chen extraida de Liu arXiv:2203.07871v1 (eqs. (2)-(5), Lemma 9 de Chen 1973) y Matomaki-Zuniga-Alterman 2025 Ex. 3.4(ii)-(iii): z = N^{1/10} (variante N^{1/8}), y = N^{1/3}, coeficiente de switching 1/2, cuna alpha in [1/10,1/3], beta in [1/3,(1-alpha)/2].
- **[T]** Respuesta al doc 51: la correccion es correcta para el marco general MZA, pero el criterio SI se aplica a la submatriz de switching P3/Omega de la identidad de Chen (ambas ramas cuentan los mismos (q,r,l) con el mismo factor local S(N;q,r) del doc 44). Ahi tambien falla.
- **[T]/[P-rev]** Margenes de escala principal en los 48 chequeos (2 modelos x 2 (v,u) x 4 h x 3 matrices): matriz de Chen D1 = 0.50-0.89 (malla) y 0.83-0.98 (pesos pi); submatriz de switching D1 = D2 = 1.0000 exactamente en la malla (todas las entradas <= 0: N3o = (j+1)+N3p en la cuna) y 0.83-0.98 con pesos pi; rama favorable de Kuhn D1 = D2 = 1.0000 exactamente; Richert (lambda=1) D1 = 0.75-0.98, D2 = 0.98-1.00. El criterio (46.4) falla para el peso real de Chen; la mejora de segundo orden esta perdida (factor ~N^{1/10}).
- Evidencia: matematicas/52_margenes_peso_real_chen_criterio_46.md; python numerico/verifica_margenes_chen_richert.py PASS en el PC (Fraction, v=10,8; u=3; h=12,24,48,96; malla y pesos pi). Siguiente: disenar un peso NUEVO con margenes cero por construccion, o instanciar b_{P,l} y auditar A4-A5 (doc 51), o volver al frente analitico (C10/W(d)). No hay cota inferior nueva de capas y Goldbach sigue abierta.
- **Reclamo cerrado:** la extraccion Chen/Richert reclamada a las 18:15 ET esta entregada (doc 52); el bucle puede retomar sus frentes sin restriccion.

## Reclamo LUMINA/Muse - 2026-11-25 18:36 ET - diseno peso nuevo margenes cero

- Muse toma el siguiente paso del doc 52: disenar un peso NUEVO c(q,r) con margenes cero por construccion (atomos rectangulares) que satisfaga el criterio del doc 46 y conserve un balance positivo tipo Chen.
- El bucle NO debe duplicar este paso: seguir con los otros frentes (instanciar b_{P,l} para A_N y auditar A4-A5; C10/W(d); uniones C_k+C_{k+1}; variante Lema 4 Pintz-Althofer) hasta que Muse entregue el doc.
- Entrega esperada: doc nuevo en matematicas/ + entrada ESTADO.md con construccion y verificacion.

## Bucle LUMINA/OpenClaw - 2026-11-25 18:44 ET - instanciacion switched para Goldbach

- **[T]** Para `A_N={(ell,t):ell+t=N}`, la sucesion publicada se especializa exactamente a `b_{P,ell}=E_P(N-ell)`; sus sumas con `d|ell` son la distribucion de casi-primos ponderados `t=p_1m` en `t=N (mod d)`.
- **[T]/[P-rev]** Con `X_1=pi(N)`, A5 se normaliza (cociente exactamente 1, salvo el extremo `t=N`); A4 no sale de BV para primos: exige un BV uniforme para la convolucion restringida `E_P`. El switching quita el modulo externo `qrd` pero conserva la dificultad en los coeficientes.
- Evidencia: `matematicas/53_instanciacion_bPell_goldbach.md`; `python numerico/verifica_bPell_goldbach.py` PASS en 252 configuraciones; `goldbach.py selftest` OK. Siguiente: descomponer `E_P` por Buchstab y auditar un principio de induccion de Bombieri al nivel 0.267. No hay cota nueva de capas y Goldbach sigue abierta.
- **Memoria compartida:** propuesta pendiente; MCP no expuesto, sin persistencia remota afirmada.

## Bucle LUMINA/OpenClaw - 2026-11-25 19:14 ET - induccion de Bombieri para el peso switched

- **[T-lit]/[T]** El principio de Motohashi y su version moderna para productos de primos en intervalos dan nivel `<1/2` para cada caja no ponderada de `E_P`; por tanto `0.267` no es la barrera de una caja fija.
- **[T]/[P-rev]** El peso real `W_P^+` acopla los factores. Discretizarlo con error `delta` y transferir por desigualdad triangular pierde un factor `Q=N^theta`; A4 por esa via exige malla `eta << N^(-theta)log^(-100)`, con al menos `N^(k theta)` cajas, no una familia polilogaritmica absorbible.
- Evidencia: `matematicas/57_motohashi_cubre_cajas_no_peso_acoplado.md`. Siguiente: auditar una version ponderada uniforme de Motohashi/BFI para el coeficiente multivariable antes de sumar valores absolutos. A4 no esta probada ni refutada; Goldbach sigue abierta.
- **Memoria compartida:** propuesta pendiente; MCP no expuesto, sin persistencia remota afirmada.


## Bucle LUMINA/OpenClaw - 2026-11-25 19:44 ET - layer cake para el peso acoplado de A4

- **[T]** La parte acoplada `(A-u)_+` de `W_P^+` es exactamente `integral_0^A 1_{u<t} dt`; al integrar primero las discrepancias, A4 se reduce a una familia uniparametrica anidada de subniveles, sin el factor `Q` ni las `eta^{-k}` cajas de la discretizacion ingenua del doc 57.
- **[P]** Falta un BV uniforme en el umbral para las regiones `sum_i w(x_i)<t`; Liu--Park--Song, tal como fue auditado, solo cubre cajas de intervalos. No se supone control de la banda de frontera.
- Evidencia: `matematicas/60_layer_cake_peso_acoplado_A4.md`; `python numerico/verifica_layer_cake_A4.py` PASS (5000 identidades y 500 funcionales racionales exactos). Siguiente: auditar la uniformidad en subniveles en la prueba de Motohashi/Liu--Park--Song. A4 y Goldbach siguen abiertos.
- **Memoria compartida:** propuesta pendiente; MCP no expuesto, sin persistencia remota afirmada.

## Bucle LUMINA/OpenClaw - 2026-11-25 20:14 ET - peso Kuhn de A4 como estado finito

- **[T-lit]/[T]** El peso concreto del Teorema 3.2 de MZA es Kuhn: `w=(1/2)1_[1/20,1/6]`. Por ello `W_P^+` depende solo de `j+k` y se atomiza exactamente en un numero acotado de cajas (`k=0` o `k=1`); la familia curva del doc 60 no es necesaria en esta instancia. LPS Teor. 3.2 es uniforme en los extremos de los intervalos y cubre cada atomo a todo nivel `<1/2`, incluido `0.267`.
- **[P-rev]** A4 aun no se declara cerrada: faltan las clases `(d,N)>1`, resolver explicitamente la incompatibilidad textual entre la hipotesis Lipschitz de MZA y su indicador Kuhn discontinuo en los bloques frontera, y convertir la cota absoluta de LPS a la normalizacion relativa en el ultimo bloque degenerado.
- Evidencia: `matematicas/62_atomizacion_kuhn_A4.md`; `python numerico/62_atomizacion_kuhn_A4.py` PASS (12 estados + 20000 configuraciones). Siguiente: cerrar clases no reducidas y bloques frontera. Goldbach sigue abierta.
- **Memoria compartida:** propuesta pendiente; MCP no expuesto, sin persistencia remota afirmada.

## Bucle LUMINA/OpenClaw - 2026-11-25 20:49 ET - clases no reducidas de A4 en Lean

- **[T-Lean]** Si `g|N` y `g|(N-t)` con `t<=N`, entonces `g|t`; en particular, todo primo comun a `d,N` en la clase `d|(N-t)` divide al termino switched `t`.
- **[T-Lean]** Como los factores admisibles de `t` son >= `N^(1/20)`, las clases `(d,N)>1` solo proceden de primos grandes `p|N` que ya aparecen en `t`. Target Lean verde; axiomas: `[propext]`.
- **[P]** Falta sumar uniformemente las rebanadas con primo fijo `p|N` y probar ahorro `log^-100`; A4 sigue abierta. El build global fallo por el asignador al abrir un `.olean.private`, asi que no hubo publicacion GitHub.
- Evidencia: `matematicas/64_clases_no_reducidas_Lean.md`, `lean/NonReducedClasses.lean`. Goldbach sigue abierta.

## Reanudacion del bucle LUMINA/OpenClaw - 2026-11-25 21:16 ET

- **[T-Lean]** Recompilados `GoldbachLumina.NonReducedClasses` y el proyecto Lean completo: ambos exit 0; `#print axioms` permanece en `[propext]`.
- Se corrige explicitamente el bloqueo transitorio de las 20:49 ET: no afectaba al teorema y no se reprodujo. Publicacion GitHub habilitada tras escaneo de secretos.
- **[P]** La cota uniforme de las rebanadas `p|N, p|t` sigue abierta; A4 y Goldbach siguen abiertas.


## Bucle LUMINA/OpenClaw - 2026-11-25 21:36 ET - cota absoluta de clases no reducidas A4

- **[T]** Por doble conteo, la norma no reducida es `sum_t E(t) #{d|N-t:(d,N)>1}`; cada modulo es un divisor real de `N-t`, asi que no aparece el factor `Q`.
- **[T]** Como el doc 64 fuerza un primo comun `p|N,t` con `p>=N^(1/20)`, la contribucion total es `<= N^(19/20+o(1))`, uniformemente en `Q`. En bloques con `X_2(alpha)>=N/log^B N`, esto da el ahorro relativo `log^-A` exigido por A4.
- **[P-rev]** Queda justificar la normalizacion en bloques degenerados y las fronteras Lipschitz/Kuhn; A4 global y Goldbach siguen abiertos.
- Evidencia: `matematicas/66_cota_clases_no_reducidas_A4.md`; `python numerico/66_cota_clases_no_reducidas_A4.py` PASS (320 configuraciones exactas).
- **Memoria compartida:** propuesta pendiente; MCP no expuesto, sin persistencia remota afirmada.

## Bucle LUMINA/OpenClaw - 2026-11-25 22:09 ET - pelado del bloque diadico final A4

- **[T]** La suma switched solo tiene un bloque potencialmente truncado: P_*<p_1<N^(1/4), con eta=N^(1/4)/P_*-1.
- **[T]** Su masa es <= M N eta+O(MN^(3/4)); si eta<=log^-4 N, es o(N/log^2 N) y puede pelarse antes de la criba. El caso de masa arbitrariamente pequena ya no exige una estimacion A4 relativa.
- **[P-rev]** Para eta>log^-4 N falta una cota inferior uniforme de X_2; tambien sigue abierta la incompatibilidad textual Lipschitz/Kuhn. Goldbach sigue abierta.
- Evidencia: matematicas/67_pelado_bloque_final_A4.md; python numerico/67_pelado_bloque_final_A4.py PASS (400 valores).

## Bucle LUMINA/OpenClaw - 2026-11-25 22:54 ET - masa del bloque final A4

- **[T]** Si el ancho relativo final satisface `eta>log^-4 N`, cuatro ventanas disjuntas de primos dentro de `(P_*,N^(1/4))` dan `X_2(P_*) >> N eta^4/log^4 N >= N/log^20 N`.
- Combinado con el pelado del doc 67 para `eta<=log^-4 N` y la cota no reducida del doc 66, queda cerrada la obstruccion de normalizacion del ultimo bloque.
- **[P-rev]** A4 global aun requiere resolver/citar la incompatibilidad textual Lipschitz/indicador Kuhn; Goldbach sigue abierta.
- Evidencia: `matematicas/68_masa_bloque_final_A4.md`; `python numerico/68_masa_bloque_final_A4.py` PASS (1000 configuraciones racionales).

## Bucle LUMINA/OpenClaw - 2026-11-25 23:06 ET - reparacion Lipschitz/Kuhn

- **[T-lit]** Assumption 3.1 de MZA exige `w` Lipschitz, mientras su Teorema 3.2 fija el indicador discontinuo de Kuhn; el antecedente publicado es literalmente inconsistente y no se corrige en silencio.
- **[T]** Refinar los bloques diadicos insertando `x^(1/6)` hace `w` constante por celda y da `W(p_1m)=W_I(m)` exactamente, sin `C_w/log x`; la suma de Sigma_1 se trata directamente en `[x^(1/20),x^(1/6))` y la evaluacion Kuhn usa (5.24), no el paso Lipschitz (5.22).
- **[P-rev]** Junto con los docs 62 y 66--68, esto elimina la incompatibilidad Lipschitz/indicador de la auditoria A4; no se declaran verificadas A1--A3/A5 ni Goldbach.
- Evidencia: `matematicas/71_reparacion_lipschitz_kuhn.md`; `python numerico/71_bloques_alineados_kuhn.py` PASS (2000 particiones, 20000 estados).
- **Memoria compartida:** propuesta pendiente; MCP no expuesto, sin persistencia remota afirmada.

## Bucle LUMINA/OpenClaw - 2026-11-25 23:36 ET - auditoria A1--A3 y A5 de MZA

- **[T]** Para `a_n=1_{N-n primo}`, Bombieri--Vinogradov mas el tratamiento elemental de `(d,N)>1` verifica A1 a nivel `0.267` con `g_N(p)=0` si `p|N` y `p/(p-1)` en otro caso.
- **[T]** El peso Kuhn satisface A2 con `R=6`; A3 es `O(N^(19/20))=o(X_1/log^2 N)`.
- **[P-rev]** A5 se normaliza en la reparacion por bloques alineados; queda A4 como unico supuesto aritmetico sustantivo de esta linea. Incluso cerrando A4, MZA produce `p+P_3`, no Goldbach fuerte.
- Evidencia: `matematicas/73_auditoria_A1_A3_A5_goldbach.md`; `python numerico/73_auditoria_A1_A3_A5.py` PASS.
## Bucle LUMINA/OpenClaw - 2026-11-26 00:06 ET - revision adversarial de SL-70.1

- **[P-rev]/[FALSO como implicacion]** La meseta de correlacion de pares y su varianza equivalente son estadisticas cuadraticas promediadas; las condiciones (a)--(c) del doc 72 no bastan, tal como estan escritas, para el `L1` absoluto MM-69.1 ni para el coeficiente con fase de cada `N`.
- **[T]** Una traslacion ciclica preserva exactamente espectro de potencia y autocorrelacion pero cambia un coeficiente puntual; el contraejemplo finito separa informacion de dos puntos de control puntual. No refuta Montgomery.
- Siguiente: formular un segundo momento explicito de `E_U(N)` con fases y auditar por separado el paso de casi todo `N` a todo `N`. Goldbach sigue abierta.
- Evidencia: `matematicas/74_revision_adversarial_SL70_correlacion.md`; `python numerico/74_revision_SL70_fase.py` PASS.


## Bucle LUMINA/OpenClaw - 2026-11-26 00:37 ET - revision adversarial del doc 77

- **[P-rev]/[FALSO como prueba]** La Prop. 77.1 usa `sum rho^2 << X/log X`, pero el valor `1.467*pi(2X)` fue medido solo para el proxy `rho_tilde`, no para el minorante verdadero con la parte `rho_B`; una cota de masa `O(delta)` no controla cuadrados.
- **[T]** Un contraejemplo finito preserva la minoracion y la masa L1 pero cambia L2 de `P+M` a `P+M^2`. Parseval se repara si se demuestra el sublema `sum rho^2 << X(log X)^K` para algun K fijo.
- Evidencia: `matematicas/78_revision_adversarial_segundo_momento_rho.md`; `python numerico/78_gap_L2_rho.py` PASS 6/6. Goldbach sigue abierta.


## Bucle LUMINA/OpenClaw - 2026-11-26 01:11 ET - cierre L2 del minorante de Harman

- **[T]** Desde la definicion exacta del Lema 6.8 de Maynard III, `rho_B(n) <= tau_4(n)`; todas las piezas restantes del minorante son sumas sobre a lo sumo seis slots multiplicativos, asi que `|rho(n)| <= (9+C) tau_6(n)`.
- **[T]** Por el momento clasico de la funcion divisor, `sum_{n~X}|rho(n)|^2 << X(log X)^35`. Esto cierra el sublema `L2-rho` del doc 78 y restaura la Prop. 77.1 / (74.4) para el minorante verdadero tras reajustar el ahorro de Vinogradov.
- **[N]** La constante experimental `1.467*pi(2X)` sigue siendo solo del proxy; el puente casi-todo-N a todo-N permanece abierto. Evidencia: `matematicas/80_cota_divisores_rho.md`; `python numerico/80_cota_divisores_rho.py` PASS 5501/5501. Goldbach sigue abierta.


## Bucle LUMINA/OpenClaw - 2026-11-26 01:37 ET - multiplicidad exacta de las uniones de capas

- **[T]** Al seleccionar un subconjunto no ordenado de `k-1` primos de `rad(N-p)` y exigir que el radical restante tenga uno o dos factores, la incidencia exacta es `I_k = k C_k + binom(k+1,2) C_{k+1}`.
- **[P-rev]** El factor `1/(k-1)!` sugerido en el doc 16 solo convierte tuplas ordenadas en subconjuntos; no produce coeficiente uno en ambas capas. Una cota inferior de `I_k` si implica una cota inferior de `C_k+C_{k+1}` tras dividir por `binom(k+1,2)`.
- Evidencia: `matematicas/81_multiplicidad_union_capas.md`; `python numerico/81_incidencia_union_capas.py` PASS 319991/319991. Siguiente: uniformidad analitica de Chen para la incidencia con `Q<=N^theta`. Goldbach sigue abierta.

## Bucle LUMINA/OpenClaw - 2026-11-26 02:12 ET - reduccion exacta de la incidencia a formas de Chen

- **[T]** Para cada seleccion `Q`, la factorizacion correcta es `N-p=Q u v`, con `rad(u)|Q`, `(v,Q)=1` y `omega(v) in {1,2}`; los exponentes repetidos de primos seleccionados viven en `u`.
- **[T]/[P-rev]** Siempre `(Q,N)=1`, pero la uniformidad de la incidencia completa debe ser en `A=Qu`, no solo en `Q`. Restringir a la subincidencia positiva `u=1` evita ese problema.
- **[T]** La serie singular local es `2 C_2 prod_{ell|NQ,ell>2}(ell-1)/(ell-2)`, positiva e independiente de los exponentes de `u`. Evidencia: `matematicas/82_reduccion_exacta_incidencia_chen.md`; script PASS H1=366399 H2=640050. Siguiente: cota uniforme de Chen para la subincidencia `u=1` en ventanas disjuntas. Goldbach sigue abierta.
## Bucle LUMINA/OpenClaw - 2026-11-26 02:37 ET - selector polilogaritmico para C_k+C_{k+1}

- **[T]** Para k fijo, elegir los k-1 primeros primos impares que no dividen N da un selector Q_N <<_k (log N)^(k-1); una sola subincidencia u=1 basta para positividad de C_k+C_{k+1}.
- **[T]** El input de distribucion para a_Q(v)=Lambda(N-Qv) es un subconjunto de Bombieri--Vinogradov sobre modulos Qd y conserva nivel 1/2 salvo perdidas polilogaritmicas.
- **[P]** Falta auditar el switching completo de Chen uniformemente para Q<=log^B N; no basta citar BV. **[N]/[FALSO]** la hipotesis finita H2 fallo en 11 casos para k=3 y 1895 para k=4 hasta 50000 (H1 exacta PASS).
- Evidencia: `matematicas/83_selector_polilog_y_objetivo_chen.md`; `python numerico/83_selector_polilog_chen.py` exit 0. Siguiente: auditar Chen con modulo exterior polilogaritmico y el posible modulo excepcional. Goldbach sigue abierta.
