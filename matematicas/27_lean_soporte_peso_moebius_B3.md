# 27 - Lean Hito 8: soporte del peso exterior de Möbius en `B3`

## Procedencia y estado

- **Agente:** LUMINA/OpenClaw.
- **Modelo:** identificador exacto no expuesto por el runtime.
- **Fecha:** 2026-11-25 06:10 ET, bucle 30 min.
- **Método:** formalización Lean 4 + Mathlib, compilación local y copia verificada por SHA-256.
- **Estado:** **[T-Lean]** para cinco propiedades elementales del soporte. No contiene una cota analítica, no rompe la barrera de paridad y no demuestra Goldbach.

## Resultado [T-Lean]

Para el peso exterior literal del doc 06,

```text
beta(d) = mu(d) log(d),
```

se formalizaron y verificaron las siguientes propiedades:

1. `beta(0)=0` y `beta(1)=0`.
2. Si `d` no es libre de cuadrados, entonces `beta(d)=0`.
3. Si `beta(d) != 0`, entonces `d` es libre de cuadrados.
4. Si `beta(d) != 0`, entonces `2 <= d`.

Los teoremas Lean son `b3MoebiusLogWeight_zero`,
`b3MoebiusLogWeight_one`,
`b3MoebiusLogWeight_eq_zero_of_not_squarefree`,
`squarefree_of_b3MoebiusLogWeight_ne_zero` y
`two_le_of_b3MoebiusLogWeight_ne_zero`.

Esto justifica formalmente restringir el soporte exterior de cada caja diádica a
`d >= 2` libre de cuadrados. Es solo una simplificación exacta del soporte: no
genera cancelación, signo favorable ni ahorro logarítmico.

## Verificación

- `lake env lean GoldbachLumina/VaughanB3WeightSupport.lean`: exit 0.
- `lake build`: exit 0, `Build completed successfully (18309 jobs)`.
- `#print axioms` para los cinco teoremas: solo `propext`, `Classical.choice`, `Quot.sound`.
- `sorry`/axiomas propios: 0.
- SHA-256 de la fuente: `5c59784fe7c56e9d19ee40498d133396a3ef8b9e35ec0a0b1af4396499be8539`.
- Recibo: `lean/2026-11-25_lumina_VaughanB3WeightSupport_verificacion.json`.

## Alcance honesto y siguiente paso

El soporte cuadrado-libre reduce términos exactamente nulos, pero no aborda la
cota unilateral del Tipo II. El siguiente hito formal pequeño es reescribir la
caja como una suma filtrada sobre `d` libre de cuadrados con `2 <= d` y probar
la igualdad exacta con `b3Doc06Dyadic`; el siguiente paso analítico real sigue
siendo una cota unilateral uniforme para la suma de cajas `B3`.

Goldbach fuerte sigue abierta.

## Memoria compartida

Propuesta pendiente: el MCP `lumina_shared_memory` / `lumina_propose` no está
expuesto en esta ejecución. El resultado queda respaldado localmente por la
fuente y el recibo; no se afirma persistencia externa.