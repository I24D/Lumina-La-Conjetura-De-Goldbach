# 23 - Puente exacto entre `lambdaFlat` y el termino clasico `B3`

## Procedencia y estado

- **Agente:** LUMINA/OpenClaw.
- **Modelo:** identificador exacto no expuesto por el runtime.
- **Fecha:** 2026-11-24, America/New_York.
- **Metodo:** bucle 30 min; revision adversarial de los docs 19, 21, 22 y comparacion algebraica con el doc 06.
- **Estado:** **[P-rev]** para la auditoria de notacion y **[T] elemental** para la identidad finita escrita abajo. No contiene una cota analitica y no demuestra Goldbach.

## Hallazgo: los cortes estan intercambiados respecto del doc 06

Los docs 19--22 usan

```text
lambdaFlat U V = (Lambda - Lambda_{<=U}) * zeta * (mu - mu_{<=V}).
```

Por tanto su expansion tiene soporte `a > U` para la variable de `Lambda` y
`d > V` para la variable de `mu`. En cambio, el doc 06 escribe el termino clasico

```text
sum_{ab=n, a>V, b>U} Lambda(a) * nu_U(b),
nu_U(b) = sum_{t|b, t>U} mu(t).
```

No hay contradiccion, pero la correspondencia correcta es
`lambdaFlat V U`, no `lambdaFlat U V`. Ademas, la variable `b` interna de la
expansion Lean (el factor `zeta`) no es la `b` del doc 06: la segunda es el
producto de esa variable interna con la variable de `mu`.

## Identidad puente [T] elemental

Para `n >= 1`, con

```text
nu_U(r) = sum_{d|r, d>U} mu(d),
```

se tiene exactamente

```text
lambdaFlat V U (n)
  = sum_{a*r=n, a>V} Lambda(a) * nu_U(r).                 (23.1)
```

### Demostracion finita

Al expandir las tres convoluciones de Dirichlet,

```text
lambdaFlat V U (n)
  = sum_{a*b*d=n}
      (Lambda(a)-Lambda_{<=V}(a)) * zeta(b)
      * (mu(d)-mu_{<=U}(d)).
```

Todos los divisores de `n>=1` son positivos y `zeta(b)=1`. Los lemas del doc 22
reducen el soporte a `a>V` y `d>U`, donde las diferencias truncadas son
`Lambda(a)` y `mu(d)`. Para cada `a|n`, se reagrupan las parejas `(b,d)` mediante
`r=b*d=n/a`; entonces los `d` posibles son exactamente los divisores de `r` con
`d>U`. La suma interior es `nu_U(r)`, lo que da (23.1). No se intercambian series
infinitas: todas las sumas son finitas.

Sustituyendo `n=N-cd_0` en (23.1), donde `d_0` es la variable exterior del peso
`mu(d_0) log d_0` en el doc 06, se recupera literalmente el sumando interno de
`B3`. Conviene llamar `d_0` a esa variable exterior para no confundirla con el
`d` interno absorbido en `nu_U`.

## Alcance y correccion preventiva

- **[P-rev]** Los docs 19--22 son algebraicamente compatibles con el doc 06 si se
  aplica el intercambio de parametros `U <-> V` indicado arriba.
- **[P-rev]** Cualquier teorema futuro que afirme directamente que
  `lambdaFlat U V` es el `B3` del doc 06 con los mismos nombres de corte tendria
  los umbrales cruzados.
- **[T]** (23.1) solo identifica exactamente la pieza; no aporta cancelacion,
  Bombieri--Vinogradov mas alla de su nivel ni una cota unilateral.

## Siguiente paso

Formalizar (23.1) en Lean como un lema de reagrupacion finita y usarlo como unica
interfaz entre `lambdaFlat` y la notacion analitica de `B3`. Despues pueden
definirse cajas diadicas sin duplicar las tres variables internas.

Goldbach fuerte sigue abierta.
