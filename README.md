# Proyecto Goldbach — LUMINA

**Objetivo:** trabajar hacia una demostración de la conjetura fuerte de Goldbach: todo par n ≥ 4 es suma de dos primos.
**Agentes:** Codex, Muse, OpenClaw, Cline y Claude Code, coordinados por el usuario.
**Estado honesto (2026-11-23):** la conjetura sigue **abierta**. Nadie sabe cómo demostrarla. Las técnicas conocidas, incluso suponiendo GRH, no bastan (ver `matematicas/01_marco_barreras_y_rutas.md`). Este proyecto busca avances reales y verificables, y no da por demostrado nada que no lo esté.

## Por dónde empezar (lectura obligatoria para cada agente)
1. `ESTADO.md`: qué está demostrado, qué es condicional y qué es solo evidencia.
2. `matematicas/01_marco_barreras_y_rutas.md`: barreras conocidas, trampas típicas y líneas de trabajo L1–L4.
3. `notas/2026-11-23_RH_Alpoge_relevancia_Goldbach.md`: la noticia de RH/GRH (Alpöge, noviembre de 2026) y lo que cambia para Goldbach.
4. `C:\I24D_WhatsApp\Claude_Lumina\CONTRATO_MEMORIA_COMPARTIDA_LUMINA.md`: el contrato de memoria compartida.

## Estructura
| Carpeta | Contenido |
|---|---|
| `matematicas/` | Documentos matemáticos numerados (`NN_tema.md`) |
| `notas/` | Notas fechadas: noticias, lecturas, ideas sueltas (`AAAA-MM-DD_tema.md`) |
| `numerico/` | `goldbach.py` (verificador y estadísticas) y `resultados/`, con un JSON con procedencia por ejecución |
| `lean/` | *(pendiente: falta instalar Lean 4)* Formalización |

## Reglas para todos los agentes
1. **Cabecera de procedencia en todo archivo nuevo:** agente, modelo, fecha, método, estado (PROPUESTA / revisado) y fuentes.
2. **Etiquetas de estado** en cada afirmación matemática: [T], [T-Lean], [C: H], [N], [H], [P], [P-rev], [FALSO] (definidas en `matematicas/01…`, sección 0).
3. **No sobrescribir el trabajo de otro agente.** Para discrepar, añadir una sección "Objeción (agente, fecha)" o crear un archivo nuevo. Los lemas refutados se marcan [FALSO], no se borran.
4. **La evidencia numérica no es demostración.** Ninguna revisión entre agentes que comparten usuario de Windows cuenta como independiente (contrato, punto 5).
5. **Cada lema propuesto** debe tener un enunciado preciso, pasar una prueba numérica cuando sea posible, someterse a la lista de "trampas típicas" (`matematicas/01…`, sección 5) y, si es central, formalizarse en Lean.
6. **Memoria compartida:** registrar los hallazgos como propuesta en `lumina-mcp`. Si el MCP falla, dejar aquí una nota local y anotarla como pendiente en `ESTADO.md`.
7. **Recursos:** no lanzar a la vez cálculos pesados de varios agentes. Anotar en `ESTADO.md` lo que se está ejecutando.

## Herramienta numérica
```
python numerico/goldbach.py selftest          # pruebas internas (contra fuerza bruta)
python numerico/goldbach.py verify --hasta 1e9
python numerico/goldbach.py reps --hasta 1e7  # r(n) exacto y comparación con Hardy–Littlewood
```
Requiere Python 3 y numpy. Los récords de p_min se contrastan con OEIS A025018/A025019 (`numerico/referencia_oeis_pmin.json`).


## Aporte de Codex (2026-11-23)

Procedencia: Codex, modelo exacto no expuesto; lectura primaria y derivación, estado [P].
- [03: mecanismo de paridad y comparación con Alpöge](matematicas/03_sawin_shusterman_alpoge_puente.md).
- [04: reducción a una cola primitiva, pruebas y desigualdad pendiente](matematicas/04_reduccion_mobius_goldbach.md).
- Verificador ligero: `python numerico/puente_paridad_check.py --hasta 600`.
- Ver la entrada fechada de Codex en ESTADO para resultados, límites y memoria pendiente. Se conserva el trabajo anterior.

## Aporte de LUMINA/OpenClaw (2026-11-23)

- Reindexacion exacta y verificador de `T_N^prim`.
- Calculo literal de `U_N(D)` hasta `N=10^6`: la normalizacion exigida por HmuLambda permanece cerca de 3.2, por lo que esa hipotesis queda despriorizada como blanco posiblemente demasiado fuerte.
- Nuevo blanco que conserva cancelacion: cajas bilineales firmadas `B_N(C,K)`.
- Ver `matematicas/05_diagnostico_cola_y_objetivo_bilineal.md`.

## Aporte de Claude Code (2026-11-23, tarea HB)

- [06: descomposición bilineal de B_N(C,K)](matematicas/06_descomposicion_bilineal.md). Qué rangos se controlan (BV, Vaughan Tipo I), cuáles no (Tipo II con K > √N) y por qué (clase fija más allá de 1/2; núcleo de paridad).
- Diagnóstico por cajas: `python numerico/cajas_bilineales.py --N 1e6`.

## Lean 4 + Mathlib (instalado el 2026-11-23)

- Proyecto: `C:\Users\dal_n\lean\goldbach_lumina`. Está **fuera de OneDrive a propósito**, porque la caché de Mathlib ocupa unos 5 GB.
- Compilar un archivo: `cd C:\Users\dal_n\lean\goldbach_lumina` y luego `lake env lean GoldbachLumina/<Archivo>.lean`.
- Ejemplo que compila: `GoldbachLumina/Prueba.lean` (enunciado de Goldbach, sin demostrar; un caso concreto; identidad Λ = −Σ μ log de Mathlib).
- Notación: usar `ArithmeticFunction.moebius` y `ArithmeticFunction.vonMangoldt` con nombre completo, porque las notaciones μ/Λ dan ambigüedad.

## Aporte de Claude Code (2026-11-23, segunda ronda)

- [07: sesgo Möbius–Goldbach, fórmula de arcos mayores, (HB) para casi todo N y Lean](matematicas/07_sesgo_mobius_goldbach_y_casi_todo_N.md).
- Experimentos: `python numerico/mobius_goldbach_aleatoriedad.py --hasta 1e6` y `python numerico/arcos_mayores_mobius.py --hasta 1e6 --Q 1 6 30`.
- Lean: ver `lean/README.md`, con 6 teoremas verificados por el kernel.

- **Codex, 2026-11-23:** [revisión adversarial 08 de 06/07](matematicas/08_revision_codex_06_07.md), bibliografía y experimentos de residuos hasta 10^7; [Reduccion2.lean](lean/Reduccion2.lean) verifica reindexación, (4.12) y \|Tbad\|≤K_N. Ver cierre en ESTADO y registro en `lean/README.md`. [P-rev]/[N]/[T-Lean], alcances separados; **Goldbach sigue abierta**.

## Verificacion LUMINA/OpenClaw posterior a Codex (2026-11-24)

- `lean/Reduccion2.lean` recompilado localmente con exit 0; su SHA-256 coincide con la fuente canonica y con el recibo de Codex.
- La objecion sobre clases no reducidas es valida. `matematicas/09_reduccion_clases_no_reducidas_mu.md` demuestra la reduccion algebraica exacta a clases reducidas y aisla la estimacion agregada WNR que aun falta.
- El termino Tipo II `B_3` sigue abierto incluso si WNR se cerrara. Goldbach no esta demostrada.

## Aporte autonomo continuo de LUMINA/OpenClaw (2026-11-24)

- [10: diagnostico finito de WNR](matematicas/10_diagnostico_WNR_y_prioridad_tipo_II.md): 3000 comprobaciones exactas de (9.1); la penalizacion de clases no reducidas cae de 10.3% a 1.9% entre 10^4 y 10^6. WNR parece plausible, pero no esta demostrada.
- [11: capas de paridad](matematicas/11_capas_paridad_y_funcion_generadora.md): la cola es `F(-1)`; las masas pares e impares son individualmente enormes y se cancelan casi por completo.
- [12: polinomio exacto de divisores](matematicas/12_polynomial_divisores_y_soporte_ralo.md): formula cerrada (12.4) y soporte ralo de la suma completa primitiva; Prop. 12.1 contrastada en 10000 casos aleatorios.
- Scripts: `mobius_clases_no_reducidas.py` y `capas_paridad_cola.py`. Goldbach sigue abierta; el objetivo activo continua.

## Transferencia desde RH/GRH formal (2026-11-24)

- [13: que aporta la prueba Alpoge--Claude a Goldbach](matematicas/13_riemann_alpoge_transferencia_goldbach.md).
- RH y GRH-Dirichlet se tratan como [T] formalmente verificadas por dos kernels.
- La localizacion de ceros sola no controla el Tipo II aditivo. La pista transferible es la positividad unilateral en ventanas, realizada por ganancias locales de reticulos.
- Nueva identidad puntual de Riesz: `R_1(N)=R(N)/2`; suaviza una variable sin perder el coeficiente individual N.
- [14: colapso exacto de divisores de N](matematicas/14_colapso_divisores_de_N.md): corrige el frente de 12; `G_a=0` para `omega(a)>=2`, las piezas `a=ell|N` son poli-logaritmicas y la unica pieza grande es `G_1`, esencialmente la propia correlacion de Goldbach. Verificacion: 549 casos, PASS.
