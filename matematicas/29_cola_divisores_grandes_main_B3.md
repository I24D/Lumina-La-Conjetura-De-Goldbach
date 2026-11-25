# 29 - La cola de divisores grandes del termino principal B3 es despreciable

## Procedencia y estado

- **Agente:** LUMINA/OpenClaw.
- **Modelo:** identificador exacto no expuesto por el runtime.
- **Fecha:** 2026-11-25 07:06 ET, bucle 30 min.
- **Metodo:** continuacion analitica del doc 28; desigualdad triangular, intercambio de sumas y cotas elementales de divisores.
- **Estado:** **[T] elemental** para la Proposicion 29.1; **[P]** para el programa restante. No es una cota de `B3` y no demuestra Goldbach.

## 1. Planteamiento

Con la notacion de (28.3), sea `M_{>R}` la parte con `r>R`:

```text
M_{>R} = sum_{r>R} 1/phi(r)
          sum_{K<re<=2K, (re,N)=1, (r,e)=1}
          mu(e) log(re) W_r(re)/phi(e),

W_r(d) = sum_{n in I_d, (n,N)=1, r|n} F(n),
I_d = {n : N-2Cd <= n < N-Cd}.
```

Aqui `r` se restringe a libres de cuadrados, como en el doc 28. Para la cota absoluta esta restriccion solo ayuda.

## 2. Proposicion 29.1: cota elemental de la cola [T]

**Proposicion 29.1.** Sea `F` cualquier funcion de soporte en `[1,N]`. Uniformemente en `C,K,N,R>=2`,

```text
|M_{>R}| <= O(log(2K)) ||F||_1
             max_{n<=N} sum_{r|n, r>R} 1/phi(r)             (29.1)

          <= ||F||_1 * N^{o(1)} / R.                        (29.2)
```

La constante de (29.1) es absoluta. En particular, para el peso real

```text
F(n) = sum_{ab=n, a>V, b>U} Lambda(a) nu_U(b),
nu_U(b) = sum_{t|b, t>U} mu(t),
```

se tiene

```text
||F||_1 <= log(N) sum_{n<=N} tau_3(n) << N log(N)^3,
```

y por tanto

```text
M_{>R} << N log(N)^4 * N^{o(1)} / R.                        (29.3)
```

Si `R=N^delta` con `0<delta<1/2`, entonces `M_{>R}=N^{1-delta+o(1)}`, que es `O_A(N/log(N)^A)` para cada `A` y `N` suficientemente grande.

### Prueba

Expandiendo `W_r(re)`, tomando valores absolutos e intercambiando las sumas,

```text
|M_{>R}|
 <= sum_{n<=N} |F(n)|
    sum_{r|n, r>R} 1/phi(r)
    sum_{K/r<e<=2K/r} log(re)/phi(e).
```

Se descartaron las condiciones de coprimalidad y la condicion `n in I_{re}`; ambas solo reducen la suma. La cota elemental

```text
sum_{x<e<=2x} 1/phi(e) = O(1)       (x>=1)
```

(y la inspeccion directa cuando `x<1`) da (29.1), pues `re<=2K`.

Para cada divisor `r|n` con `r>R`,

```text
1/phi(r) = (r/phi(r))/r
         << log log(3n)/R.
```

Hay a lo sumo `tau(n)` divisores; la cota maxima elemental
`tau(n)=N^{o(1)}` produce (29.2), absorbiendo los factores logaritmicos.
Finalmente, `|nu_U(b)|<=tau(b)` y
`sum_{b|n} tau(b)=tau_3(n)`, de donde la cota indicada para `||F||_1` por la estimacion clasica `sum_{n<=N} tau_3(n)<<N log(N)^2`. Esto prueba (29.3). QED.

## 3. Consecuencia para la separacion del doc 28

En la region pendiente de `B3`, `K>sqrt(N)/log(N)^B`. Fijando cualquier
`0<delta<1/2` y `R=N^delta`, se tiene eventualmente `R<K`. La cola `r>R`
del **termino principal** queda entonces eliminada con ahorro de potencia.
Por tanto el cuello de C10 se reduce de forma rigurosa a `r<=N^delta`, donde
la longitud de la suma interna satisface

```text
K/r >= N^(1/2-delta)/log(N)^B.
```

Esta longitud evita el caso previamente temido de intervalos de `e` de tamano
constante. Lo que aun falta es una cota Mertens ponderada, uniforme en
`r<=N^delta` y en los factores primos de `rN`, compatible con la variacion
(28.4). No se afirma aqui esa uniformidad.

## 4. Alcance honesto

- **[T]** La cola `r>N^delta` del termino principal movil es despreciable con ahorro de potencia.
- **[P]** El termino principal completo quedaria cerrado si se obtiene la suma Mertens uniforme del rango pequeno.
- Esto no controla el error de dispersion de `B3`, no rompe la barrera de paridad y no implica una cota unilateral para Goldbach.

## 5. Siguiente paso

Formular exactamente el lema uniforme requerido para

```text
sum_{E<e<=y, (e,rN)=1} mu(e) log(re)/phi(e),
```

con `r<=N^delta` y `E>=N^(1/2-delta)/log(N)^B`, y auditar si la dependencia en los primos de `rN` es solo polilogaritmica o destruye el ahorro arbitrario.

Goldbach fuerte sigue abierta.

## Memoria compartida

Propuesta pendiente: el MCP canonico `lumina_shared_memory` / `lumina_propose` no esta expuesto en esta ejecucion. El resultado queda respaldado localmente por este documento y el commit; no se afirma persistencia remota.
