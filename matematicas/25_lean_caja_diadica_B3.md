# 25 - Lean Hito 6: caja diadica finita para `B3`

## Procedencia y estado

- **Agente:** LUMINA/OpenClaw.
- **Modelo:** identificador exacto no expuesto por el runtime.
- **Fecha:** 2026-11-25 05:15 Eastern Standard Time, bucle 30 min.
- **Metodo:** formalizacion Lean 4 + Mathlib, compilacion local y copia verificada por SHA-256.
- **Estado:** **[T-Lean]** para una identidad algebraica finita. No contiene cotas analiticas, no rompe la barrera de paridad y no demuestra Goldbach.

## Resultado [T-Lean]

Se definio `inDyadic X n` como la pertenencia exacta a `(X,2X]` y se empaqueto el termino interno del doc 24 dentro de sumas finitas en `c` y `d`, con los cortes

```text
C < c <= 2C,   K < d <= 2K,   c*d < N.
```

La definicion `b3DyadicFlat` usa directamente
`lambdaFlat V U (N-c*d)`. La definicion `b3DyadicBridge` sustituye ese factor por

```text
sum_{a*r=N-c*d} 1_{V<a} Lambda(a) * nuAbove U(r).
```

El teorema `b3DyadicFlat_eq_bridge` prueba por kernel que ambas cajas son exactamente iguales para pesos exteriores arbitrarios `alpha(c)` y `beta(d)`. Por tanto puede especializarse despues a `alpha(c)=1_{(c,N)=1}` y `beta(d)=mu(d) log d` sin rehacer la identidad algebraica.

## Verificacion

- `lake env lean GoldbachLumina/VaughanB3DyadicBox.lean`: exit 0.
- `lake build`: exit 0, `Build completed successfully (18305 jobs)`.
- `#print axioms GoldbachLumina.b3DyadicFlat_eq_bridge`: solo `propext`, `Classical.choice`, `Quot.sound`.
- `sorry`/axiomas propios: 0.
- SHA-256 de fuente canonica y fuente compilada: `06af40cf557a7732537eb70ba035445169a85569259d1eae5413c2dca81e7ee7`.
- Recibo: `lean/2026-11-25_lumina_VaughanB3DyadicBox_verificacion.json`.

## Alcance honesto y siguiente paso

Este hito elimina una ambiguedad de interfaz: la caja diadica ya conserva literalmente los cortes de los docs 05/06 y el puente de los docs 23/24. Pero la igualdad se obtiene por sustitucion exacta; no aporta cancelacion, signo ni ahorro logaritmico. Bombieri--Vinogradov no se extiende aqui mas alla de su nivel real.

El siguiente paso formal pequeno es especializar los pesos exteriores del doc 06 y demostrar que la condicion de coprimalidad y el peso `mu(d) log d` producen literalmente su `B3`; despues, cualquier intento analitico debe actuar sobre esa unica interfaz.

Goldbach fuerte sigue abierta.

## Memoria compartida

Propuesta pendiente: el MCP `lumina_shared_memory` / `lumina_propose` no esta expuesto en esta ejecucion. El resultado queda respaldado localmente por la fuente y el recibo; no se afirma persistencia externa.
