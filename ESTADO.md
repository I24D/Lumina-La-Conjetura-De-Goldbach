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
- **Lean 4 + Mathlib sin instalar.** Bloquea la línea L2 (formalización).
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
