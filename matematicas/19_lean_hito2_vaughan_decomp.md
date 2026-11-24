# 19 — Lean Hito 2: identidad algebraica de Vaughan

## Procedencia y estado

- **Agente:** LUMINA/OpenClaw (Codex; modelo exacto no expuesto al runtime).
- **Fecha:** 2026-11-24, America/New_York.
- **Método:** preflight de la formalización pública de Mellendijk, port mínimo a Mathlib y verificación local con Lean 4.
- **Contexto:** bucle 30 min, iteración de las 18:36 ET.
- **Estado:** **[T-Lean]** para la identidad algebraica; no contiene ninguna estimación analítica.
- **Fuente portada:** `amellendijk/lean-bombieri-vinogradov`, commit `b1b0904cdaf0ddf03e410c86ae4bebdb0eace131`, `BV/Defs.lean` (`Lambda_decomp`).

## Resultado [T-Lean]

Se creó `lean/VaughanDecomp.lean`, copia del archivo compilado en
`C:\Users\dal_n\lean\goldbach_lumina\GoldbachLumina\VaughanDecomp.lean`.
El port elimina la infraestructura global `ProofData` del repositorio fuente y
hace explícitos los cortes naturales `U,V`. Define

- `lambdaLE U = 1_{[1,U]} Λ`;
- `moebiusLE V = 1_{[1,V]} μ`;
- `lambdaSharp U V = μ_{≤V} * log - (Λ_{≤U} * μ_{≤V}) * ζ`;
- `lambdaFlat U V = (Λ-Λ_{≤U}) * ζ * (μ-μ_{≤V})`.

El kernel verifica, para todo `U,V,n : ℕ`,

```lean
theorem vaughan_decomp (U V n : ℕ) :
    ArithmeticFunction.vonMangoldt n =
      lambdaSharp U V n + lambdaFlat U V n + lambdaLE U n
```

La prueba usa únicamente el álgebra de convolución de funciones aritméticas y
las identidades ya presentes en Mathlib
`vonMangoldt_mul_zeta` y `coe_zeta_mul_coe_moebius`.

## Verificación

- `lake env lean GoldbachLumina/VaughanDecomp.lean`: **exit 0**.
- `lake build`: **Build completed successfully (7 jobs)**.
- `#print axioms GoldbachLumina.vaughan_decomp`:
  `[propext, Classical.choice, Quot.sound]`.
- Sin `sorry`, `axiom` propio ni hipótesis matemática añadida.

## Alcance honesto

Este hito empaqueta exactamente la descomposición que origina las piezas Tipo I
y Tipo II, pero no formaliza Bombieri–Vinogradov, la gran criba ni la cota
unilateral de `B3`. Por tanto no rompe la barrera de paridad y no demuestra
Goldbach. Su valor es eliminar una deuda formal concreta del doc 17 y dejar una
interfaz con `U,V` explícitos para formalizaciones posteriores.

## Siguiente paso

Relacionar `lambdaFlat U V` con la suma diádica `B3` de `06`, empezando por un
lema puramente finito que expanda la convolución en divisores y conserve todos
los cortes; las estimaciones analíticas deben permanecer separadas.
