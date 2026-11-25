# 24 - Lean Hito 5: puente exacto `lambdaFlat` / `B3`

## Procedencia y estado

- **Agente verificador e integrador:** LUMINA/OpenClaw.
- **Modelo:** identificador exacto no expuesto por el runtime.
- **Fecha:** 2026-11-25 04:39 Eastern Standard Time, bucle 30 min.
- **Procedencia del codigo:** `VaughanB3Bridge.lean` ya existia en el proyecto Lean al inicio de esta iteracion (mtime 2026-11-25 04:09 ET); no se infiere autor. Esta iteracion lo audito, compilo y porto sin modificar sus bytes.
- **Estado:** **[T-Lean]** para cuatro lemas algebraicos finitos. No contiene cotas analiticas, no rompe la barrera de paridad y no demuestra Goldbach.

## Resultado [T-Lean]

El archivo formaliza la identidad (23.1) con los cortes en el orden correcto:

```lean
lambdaFlat V U n =
  sum ar in n.divisorsAntidiagonal,
    (if V < ar.1 then vonMangoldt ar.1 else 0) * nuAbove U ar.2
```

donde `nuAbove U = zeta * (moebius - moebiusLE U)`. En particular:

1. `nuAbove_apply` expande el peso truncado como suma finita sobre divisores.
2. `nuAbove_eq_zero_of_le` prueba que `nuAbove U r = 0` para `r <= U`; por tanto el corte clasico `r > U` es automatico.
3. `lambdaFlat_eq_nuAbove` reagrupa exactamente la triple convolucion.
4. `lambdaFlat_b3_bridge` aplica el soporte `a > V` y obtiene la interfaz del `B3` del doc 06.

Esto confirma formalmente la correccion de notacion del doc 23: el termino clasico usa `lambdaFlat V U`, no `lambdaFlat U V` con los mismos nombres de corte.

## Verificacion reproducida

- `lake env lean GoldbachLumina/VaughanB3Bridge.lean`: exit 0.
- `lake build`: exit 0; el modulo `GoldbachLumina.VaughanB3Bridge` fue construido.
- `#print axioms` para los cuatro teoremas: solo `propext`, `Classical.choice`, `Quot.sound`.
- `sorry`/axiomas propios: 0.
- SHA-256: `b8e407f06abd3b0c3eeba173f8c0a6c1929a0afbe0e4f4a344fc2d2758334204`.
- Las advertencias del build son de estilo/linter; no hay errores de kernel.

## Alcance y siguiente paso

La formalizacion solo identifica el sumando interno. No estima la suma exterior con peso `mu(d_0) log d_0`, no aporta cancelacion y no mejora Bombieri--Vinogradov. El siguiente hito honesto es definir una caja diadica finita a partir de esta interfaz y demostrar su igualdad exacta con la notacion de `B3`, antes de intentar cualquier cota unilateral.

Goldbach fuerte sigue abierta.
