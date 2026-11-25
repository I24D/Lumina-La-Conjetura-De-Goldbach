# 49 - Cancelacion marginal de primer orden formalizada en Lean

## Procedencia y estado

- **Agente:** LUMINA/OpenClaw.
- **Modelo:** GPT-6.2 (runtime OpenClaw).
- **Fecha:** 2026-11-25, bucle 30 min.
- **Metodo:** formalizacion finita en Lean 4 + Mathlib, compilacion completa y auditoria `#print axioms`.
- **Estado:** **[T-Lean]** para la identidad algebraica finita del doc 46. No formaliza ninguna cota analitica, no construye pesos Chen/Richert y no demuestra Goldbach.

## 1. Enunciado formalizado

Sean `I,J` tipos finitos, `R` un anillo conmutativo, `c : I -> J -> R` una matriz de pesos y `a,b` perturbaciones locales. Bajo las hipotesis de margenes nulos

```text
forall i, sum_j c(i,j)=0,
forall j, sum_i c(i,j)=0,
```

Lean verifica exactamente

```text
sum_i sum_j c(i,j)(1+a(i))(1+b(j))
  = sum_i sum_j c(i,j)a(i)b(j).                    (49.1)
```

Esto es la version kernel-verificada de la cancelacion de orden cero y primer orden usada en el doc 46. La igualdad conserva solo la interaccion mixta; no usa estimaciones asintoticas ni supuestos probabilisticos.

## 2. Artefacto y prueba

Fuente compilada:

```text
C:\Users\dal_n\lean\goldbach_lumina\GoldbachLumina\MarginalCancellation.lean
```

Copia portable del proyecto:

```text
lean/MarginalCancellation.lean
```

Teoremas:

- `GoldbachLumina.marginalCancellation` **[T-Lean]**: identidad (49.1).
- `GoldbachLumina.zeroTotal_of_zeroRows` **[T-Lean]**: margenes de fila nulos implican masa total nula.

Verificacion ejecutada:

```text
lake env lean GoldbachLumina/MarginalCancellation.lean   # exit 0
lake build                                               # 18311 jobs, exit 0
```

`#print axioms` informa solo:

```text
[propext, Classical.choice, Quot.sound]
```

Son los axiomas estandar ya presentes en Mathlib; no hay `sorry`, axioma ad hoc ni hipotesis analitica oculta.

## 3. Alcance para el frente de switching

**[T-Lean]** Si un peso firmado real tiene ambos margenes exactamente nulos, la variacion separable de la serie singular `S_N(1+a_q)(1+b_r)` queda reducida de forma exacta al termino mixto. Este paso ya no depende de algebra informal.

**[P-rev]** La formalizacion no demuestra que los pesos de Chen/Richert tengan margenes nulos. Los docs 47--48 muestran que los primeros modelos de incidencia de Buchstab no satisfacen esa condicion. Por tanto, (49.1) es una interfaz verificada que un peso futuro debe cumplir, no una cota inferior de criba.

## 4. Siguiente paso

Portar a Lean la descomposicion canonica de una matriz en interaccion de margenes nulos mas proyeccion aditiva (doc 47), o bien instanciar (49.1) con coeficientes y cortes extraidos de una identidad publicada de Chen/Richert.

## Alcance honesto

- **[T-Lean]** Identidad algebraica finita verificada por el kernel.
- No hay nueva cota de `B3`, de switching ni de una capa `C_k`.
- Goldbach sigue abierta.

## Memoria compartida

Propuesta pendiente: el MCP canonico `lumina_shared_memory` / `lumina_propose` no esta expuesto. No se afirma persistencia remota.
