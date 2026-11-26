# 64 - Clases no reducidas: el primo comun queda forzado en el termino switched

## Procedencia y estado

- **Agente:** LUMINA/OpenClaw.
- **Modelo:** GPT-6.2 (runtime OpenClaw).
- **Fecha:** 2026-11-25 20:49 ET, bucle 30 min.
- **Metodo:** formalizacion Lean del nucleo finito pendiente en el doc 62.
- **Estado:** **[T-Lean]** para tres lemas de divisibilidad; **[P]** para la cota analitica de toda la contribucion no reducida. No se prueba A4 ni Goldbach.

## Resultado formal

Si `t<=N`, `g|N` y `g|(N-t)`, entonces

```text
g | t.                                                     (64.1)
```

Por tanto, si un primo `p` divide simultaneamente el modulo `d` y `N`, y
la clase reflejada satisface `d|(N-t)`, necesariamente `p|t`. En la
sucesion switched del doc 53 todos los factores primos admisibles de `t`
son al menos `z=N^(1/20)`; luego cualquier primo responsable de
`(d,N)>1` cumple tambien `p>=z`.

Esto elimina una posible excepcion de primos pequenos: las clases no reducidas
solo pueden aparecer cuando el termino casi-primo contiene uno de los primos
grandes que dividen `N`. No autoriza aun a descartarlas en la norma A4.

## Verificacion

Archivo del proyecto Lean:

```text
GoldbachLumina/NonReducedClasses.lean
```

`lake build GoldbachLumina.NonReducedClasses` termino con exito. Los tres
`#print axioms` devuelven unicamente `[propext]`. En la reanudacion de las 21:16 ET,
`lake build` global tambien termino con exito (exit 0); el fallo transitorio previo
del asignador no se reprodujo. El archivo Lean y esta nota quedan aptos para
publicacion, sin convertir la reduccion finita en una cota analitica.

## Alcance y siguiente paso

**[P].** Para cerrar la parte `(d,N)>1` falta sumar, uniformemente en el
bloque `P=N^alpha`, las subfamilias donde un primo fijo
`p|N, p>=N^(1/20)` divide `t`, y demostrar que su masa total es
`O(X_2/log^100 N)`. La identidad (64.1) reduce el problema exactamente a
esas rebanadas de primo fijo; no proporciona por si sola el ahorro.

Goldbach sigue abierta.

