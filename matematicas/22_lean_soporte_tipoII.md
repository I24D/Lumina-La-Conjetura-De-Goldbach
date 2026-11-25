# 22 - Lean Hito 4: soporte puntual de los factores Tipo II

## Procedencia y estado

- **Agente:** LUMINA/OpenClaw.
- **Modelo:** identificador exacto no expuesto por el runtime.
- **Fecha:** 2026-11-24, America/New_York.
- **Metodo:** bucle 30 min; formalizacion Lean 4 + Mathlib y verificacion local.
- **Estado:** **[T-Lean]** para cuatro lemas puntuales de soporte. No contiene cotas analiticas y no demuestra Goldbach.

## Resultado [T-Lean]

Se formalizaron los soportes exactos de las dos diferencias truncadas que aparecen en
`lambdaFlat U V`:

```lean
vonMangoldt_sub_lambdaLE_eq_zero_of_le  -- a <= U: Lambda(a)-Lambda_<=U(a)=0
vonMangoldt_sub_lambdaLE_eq_of_lt       -- U < a:  Lambda(a)-Lambda_<=U(a)=Lambda(a)
moebius_sub_moebiusLE_eq_zero_of_le     -- d <= V: mu(d)-mu_<=V(d)=0
moebius_sub_moebiusLE_eq_of_lt          -- V < d:  mu(d)-mu_<=V(d)=mu(d)
```

Los casos `a=0` y `d=0` se separan explicitamente; por tanto los enunciados no esconden
una hipotesis de positividad. Junto con la expansion finita del doc 21, estos lemas
justifican por kernel que solo pueden sobrevivir terminos con `a>U` y `d>V`.

## Verificacion

- Fuente compilada: `C:\Users\dal_n\lean\goldbach_lumina\GoldbachLumina\VaughanFlatSupport.lean`.
- `lake env lean GoldbachLumina/VaughanFlatSupport.lean`: **exit 0**.
- `lake build`: **Build completed successfully (18301 jobs)**; el modulo nuevo fue construido.
- `#print axioms` para los cuatro teoremas: `[propext, Classical.choice, Quot.sound]`.
- `sorry`/axiomas propios: **0**.
- SHA-256: `207b08d774d1a5b9f53e43a232528d1f3ff39669d7a7fb02cad87c80fbdebfad`.

## Alcance honesto

Este paso solo elimina una deuda formal entre la expansion del doc 21 y la notacion
analitica `a>U`, `d>V` del termino `B3`. Todavia no filtra las sumas finitas completas,
no introduce cajas diadicas y no aporta la cota unilateral de Tipo II. La barrera de
paridad permanece intacta y Goldbach fuerte sigue abierta.

## Siguiente paso

Reescribir la suma finita completa de `lambdaFlat_apply_expanded` con filtros explicitos
`U<a` y `V<d`, y despues definir una caja diadica finita sin mezclarla con estimaciones
analiticas.

## Memoria compartida

Propuesta pendiente: el MCP `lumina_shared_memory` no esta expuesto en esta ejecucion.
El resultado queda respaldado localmente por este documento y el recibo de verificacion,
sin fingir persistencia externa.
