# 21 - Lean Hito 3: expansion finita exacta del termino Tipo II `lambdaFlat`

## Procedencia y estado

- **Agente:** LUMINA/OpenClaw (modelo exacto no expuesto por el runtime).
- **Fecha:** 2026-11-24, America/New_York.
- **Metodo:** bucle 30 min; formalizacion Lean 4 + Mathlib y compilacion local.
- **Estado:** **[T-Lean]** para la expansion algebraica finita. No contiene cotas analiticas y no demuestra Goldbach.

## Resultado [T-Lean]

Se formalizo `GoldbachLumina.lambdaFlat_apply_expanded`. Para cortes `U,V` y `n`, expande exactamente

`lambdaFlat U V = ((Lambda-Lambda_{<=U}) * zeta) * (mu-mu_{<=V})`

como una suma finita exterior sobre factorizaciones `n=m*d` y una suma finita interior sobre `m=a*b`:

```lean
theorem lambdaFlat_apply_expanded (U V n : Nat) :
    lambdaFlat U V n =
      sum md in n.divisorsAntidiagonal,
        (sum ab in md.1.divisorsAntidiagonal,
          (ArithmeticFunction.vonMangoldt ab.1 - lambdaLE U ab.1) *
            ArithmeticFunction.zeta ab.2) *
        ((ArithmeticFunction.moebius : ArithmeticFunction Real) md.2 -
          moebiusLE V md.2)
```

La igualdad conserva las tres variables finitas que despues originan las cajas diadicas `B3` (`a,b,d`) y conserva explicitamente ambos cortes. Solo desarrolla dos convoluciones de Dirichlet; no intercambia series infinitas ni usa cancelacion esperada.

## Verificacion

- Fuente compilada: `C:\Users\dal_n\lean\goldbach_lumina\GoldbachLumina\VaughanFlatExpansion.lean`.
- `lake env lean GoldbachLumina/VaughanFlatExpansion.lean`: **exit 0**.
- `lake build`: **Build completed successfully (18299 jobs)**.
- `#print axioms`: `[propext, Classical.choice, Quot.sound]`.
- `sorry`/axiomas propios: **0**.
- SHA-256: `eaaef8687a43f86b85c50de79317f8c360790098efeed97838ae6beb366f7819`.

## Alcance honesto

Este hito cierra unicamente la expansion algebraica solicitada al final del doc 19. Todavia no identifica formalmente cada soporte con los cortes estrictos `a>U`, `d>V`, no introduce particiones diadicas y no aporta la cota unilateral de `B3`. La barrera de paridad permanece intacta y Goldbach fuerte sigue abierta.

## Siguiente paso

Formalizar los lemas de soporte que reemplazan las diferencias truncadas por indicadores `a>U` y `d>V` bajo `n>0`; despues empaquetar una caja diadica finita separada de cualquier estimacion analitica.
