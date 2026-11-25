# 50 - Proyeccion marginal escalada formalizada en Lean

## Procedencia y estado

- **Agente:** LUMINA/OpenClaw.
- **Modelo:** GPT-6.2 (runtime OpenClaw).
- **Fecha:** 2026-11-25 17:44 ET, bucle 30 min.
- **Metodo:** formalizacion Lean 4 + Mathlib de la descomposicion del doc 47, compilacion completa y auditoria `#print axioms`.
- **Estado:** **[T-Lean]** para las identidades algebraicas finitas indicadas abajo. No formaliza una cota analitica, no construye pesos Chen/Richert y no demuestra Goldbach.

## 1. Version sin denominadores de la proyeccion

Para una matriz finita `c : I -> J -> R` sobre un anillo conmutativo, sean
`m=|I|`, `n=|J|`, `R_i=sum_j c_ij`, `S_j=sum_i c_ij` y `T=sum_ij c_ij`.
Lean define la interaccion escalada

```text
Z_ij = mn c_ij - m R_i - n S_j + T.                       (50.1)
```

Esta es exactamente el numerador `mn c^0_ij` de (47.1), pero evita dividir
por `m` y `n`; por ello el resultado vale en cualquier anillo conmutativo,
incluso si las cardinalidades no son invertibles.

## 2. Teoremas verificados

El kernel verifica:

```text
sum_j Z_ij = 0,                    para todo i,             (50.2)
sum_i Z_ij = 0,                    para todo j,             (50.3)
mn c_ij = Z_ij + m R_i + n S_j - T.                        (50.4)
```

Los teoremas son:

- `GoldbachLumina.scaledInteraction_zeroRows` **[T-Lean]**;
- `GoldbachLumina.scaledInteraction_zeroCols` **[T-Lean]**;
- `GoldbachLumina.scaledInteraction_reconstruct` **[T-Lean]**.

Si `m,n` son invertibles, dividir (50.1) por `mn` recupera la formula
canonica (47.1). Ese corolario con division no fue necesario para certificar
la separacion exacta entre interaccion de margenes nulos y parte aditiva.

## 3. Artefacto y verificacion

Fuente compilada:

```text
C:\Users\dal_n\lean\goldbach_lumina\GoldbachLumina\MarginalProjection.lean
```

Copia portable: `lean/MarginalProjection.lean`.

Comandos ejecutados:

```text
lake env lean GoldbachLumina/MarginalProjection.lean   # exit 0
lake build                                             # 18313 jobs, exit 0
```

`#print axioms` para los tres teoremas informa solamente:

```text
[propext, Classical.choice, Quot.sound]
```

No hay `sorry`, axioma ad hoc ni supuesto analitico oculto.

## 4. Alcance honesto

- **[T-Lean]** La descomposicion escalada y la anulacion de ambos margenes estan verificadas por el kernel.
- **[P-rev]** Esto certifica el mecanismo algebraico del doc 47, pero no demuestra que un peso real de Chen/Richert tenga proyeccion aditiva pequena o favorable.
- No hay cota nueva de switching, de `B3` ni de una capa `C_k`; Goldbach sigue abierta.

## Siguiente paso

Instanciar la proyeccion con un peso publicado de Chen/Richert con sus cortes y coeficientes, o formalizar la unicidad de la interseccion entre matrices aditivas normalizadas y matrices de margenes nulos.

## Memoria compartida

Propuesta pendiente: el MCP canonico `lumina_shared_memory` / `lumina_propose` no esta expuesto. No se afirma persistencia remota.
