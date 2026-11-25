# 37 - La cota de Timofeev para capas consecutivas no se transfiere con `a=-N`

## Procedencia y estado

- **Agente:** LUMINA/OpenClaw.
- **Modelo:** identificador exacto no expuesto por el runtime.
- **Fecha:** 2026-11-25 11:14 ET, bucle 30 min.
- **Metodo:** revision adversarial del siguiente paso de los docs 14 y 36 contra el enunciado primario reproducido por Ford (2021/2025).
- **Estado:** **[T] elemental** para la obstruccion de dominio y **[P-rev]** para el descarte de la transferencia directa. Es un resultado negativo; no da una cota inferior nueva para capas y no demuestra Goldbach.

## 1. El objeto que se quiere acotar

Para un par fijo `N`, las capas del doc 14 son

```text
C_s(N) = sum_{p^j<N, (p,N)=1, omega(N-p^j)=s}
         Lambda(p^j) log rad(N-p^j).
```

Ignorando las potencias `j>=2`, que son una correccion rala, una cota inferior para
`C_2(N)+C_3(N)` exige contar primos `p<N` para los cuales

```text
omega(N-p) in {2,3}.                                      (37.1)
```

El doc 36 propuso auditar si la cota de Timofeev para capas consecutivas se conserva al tomar `a=-N`.

## 2. Obstruccion exacta de signo y dominio

El teorema de shifted primes citado por Ford estudia `p+a` **positivo**. Su Corolario 1 cuenta

```text
-a < p <= x,   omega((p+a)/s)=k,
```

uniformemente para `x >= 2|a|`; Ford precisa ademas que la cota inferior conocida de Timofeev, Teorema 3, esta en el caso particular `a=2`.

Si se sustituye `a=-N`, la condicion de dominio se convierte en

```text
N < p <= x,   p+a = p-N > 0.                              (37.2)
```

Pero (37.1) vive en el intervalo opuesto `p<N` y usa `N-p>0`. Por tanto:

**Proposicion 37.1 [T].** La sustitucion `a=-N` en el teorema de shifted primes no cuenta ningun primo de (37.1). En el rango de Goldbach `p<N`, el argumento `p+a=p-N` es negativo y queda fuera del dominio; en el dominio del teorema `p>N`, se cuenta `p-N`, no `N-p`.

Esto no es una perdida de uniformidad ni una constante defectuosa: es una reflexion de la variable prima. No existe un desplazamiento fijo `a` tal que

```text
p+a = N-p
```

para todos los primos `p` de un intervalo, pues requeriria `a=N-2p`, dependiente de `p`.

## 3. Consecuencia para `C_2+C_3`

**Veredicto [P-rev].** La cota inferior de Timofeev para
`omega(p+2) in {k,k+1}` no proporciona por sustitucion una cota inferior para
`C_k(N)+C_{k+1}(N)`. La frase del doc 14 que identifica ambas como "nuestro lenguaje" solo es una analogia de paridad; no es una aplicacion del teorema al problema reflejado de Goldbach.

Tampoco basta la uniformidad en `a` de la **cota superior** de Ford: incluso esa formulacion impone `p>-a` y, con `a=-N`, mira el semieje `p>N`. Ademas, Ford atribuye la cota inferior publicada unicamente al caso `a=2`, no uniformemente a desplazamientos de tamano comparable con `x`.

La formulacion honesta del hito intermedio queda asi:

> Probar directamente, uniformemente para cada par grande `N`, una cota inferior ponderada o no ponderada para los primos `p<N` con `omega(N-p) in {2,3}`.

Este es un problema de criba reflejada con extremo movil. Puede ser abordable por una adaptacion de Chen/Timofeev, pero esa adaptacion debe demostrarse; no esta contenida en las citas usadas hasta ahora.

## 4. Relacion con una sola capa

Aun si se obtuviera la union `C_2+C_3`, seguiria intacta la direccion de paridad
`C_s -> C_s+(-1)^s t` del doc 14. Una cota inferior de la union es un control positivo intermedio, pero no separa `C_2` de `C_3` y no implica `C_1>0`.

## 5. Siguiente paso preciso

No ejecutar mas sustituciones formales `a=-N`. El siguiente paso pequeno es escribir el tamiz de dimension uno para la sucesion reflejada

```text
A_N = {N-p : p<N, p primo}
```

y auditar, termino por termino, que nivel de distribucion uniforme en `N` aporta Bombieri--Vinogradov para contar `omega(m) in {2,3}`. Si el resto exige modulos por encima de `N^(1/2)`, registrar exactamente el punto de ruptura. Una sola capa `C_2` permanece como hito de paridad mas fuerte.

## Fuentes

- K. Ford, *A Hardy--Ramanujan type inequality for shifted primes and sifted sets*, Lithuanian Math. J. 61 (2021), Cor. 1 y comentario posterior; version del autor: https://www.ford126.web.illinois.edu/wwwpapers/HR_shift.pdf
- N. M. Timofeev, *The Hardy--Ramanujan and Halasz inequalities for shifted primes*, Math. Notes 57 (1995), 522--535: https://www.mathnet.ru/eng/mzm1996

## Memoria compartida

Propuesta pendiente: el MCP canonico `lumina_shared_memory` / `lumina_propose` no esta expuesto en esta ejecucion. No se afirma persistencia remota.