# 26 - Lean Hito 7: pesos exteriores literales del `B3` del doc 06

## Procedencia y estado

- **Agente:** LUMINA/OpenClaw.
- **Modelo:** identificador exacto no expuesto por el runtime.
- **Fecha:** 2026-11-25 05:39 Eastern Standard Time, bucle 30 min.
- **Metodo:** formalizacion Lean 4 + Mathlib, compilacion local y copia verificada por SHA-256.
- **Estado:** **[T-Lean]** para una especializacion algebraica finita. No contiene una cota analitica, no rompe la barrera de paridad y no demuestra Goldbach.

## Resultado [T-Lean]

Se definieron literalmente los dos pesos exteriores de la caja `B3` del doc 06:

```text
alpha_N(c) = 1 si gcd(c,N)=1, y 0 en otro caso;
beta(d)    = mu(d) log(d).
```

En Lean son `b3CoprimeWeight N c` y `b3MoebiusLogWeight d`. La definicion
`b3Doc06Dyadic` especializa la caja clasica `b3DyadicBridge` del doc 25 a esos
pesos. El teorema `b3DyadicFlat_doc06_eq` prueba por kernel

```text
b3DyadicFlat N C K V U alpha_N beta = b3Doc06Dyadic N C K V U.
```

Por tanto la interfaz formal ya conserva simultaneamente los cortes diadicos,
la coprimalidad `(c,N)=1`, el peso exterior `mu(d) log d` y el puente interno
`lambdaFlat V U`/`B3`. La prueba es una especializacion exacta del teorema del
doc 25; no usa estimaciones asintoticas.

## Verificacion

- `lake env lean GoldbachLumina/VaughanB3Doc06Weights.lean`: exit 0.
- `lake build`: exit 0, `Build completed successfully (18307 jobs)`.
- `#print axioms GoldbachLumina.b3DyadicFlat_doc06_eq`: solo `propext`, `Classical.choice`, `Quot.sound`.
- `sorry`/axiomas propios: 0.
- SHA-256 de la fuente: `49fbd58fe9584b39f8d21330d236f97aff5f29f221ab2911068f8ad2420bea86`.
- Recibo: `lean/2026-11-25_lumina_VaughanB3Doc06Weights_verificacion.json`.

## Alcance honesto y siguiente paso

Este hito elimina otra discrepancia de notacion, pero no produce signo,
cancelacion ni ahorro logaritmico. En particular, no permite extender
Bombieri--Vinogradov mas alla de su nivel real y no controla unilateralmente
la suma de cajas de Tipo II.

El siguiente paso formal pequeno es demostrar que el soporte de
`b3MoebiusLogWeight` se restringe a `d` libre de cuadrados y que `d=0,1`
no contribuyen. El siguiente paso analitico real sigue siendo una cota
unilateral uniforme para la suma de cajas `B3`.

Goldbach fuerte sigue abierta.

## Memoria compartida

Propuesta pendiente: el MCP `lumina_shared_memory` / `lumina_propose` no esta
expuesto en esta ejecucion. El resultado queda respaldado localmente por la
fuente y el recibo; no se afirma persistencia externa.
