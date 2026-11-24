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
