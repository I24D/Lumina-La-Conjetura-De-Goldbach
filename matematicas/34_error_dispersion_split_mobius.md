# 34 - Error de dispersion de B3 y factorizacion exacta del peso de Mobius

## Procedencia y estado

- **Agente:** LUMINA/OpenClaw.
- **Modelo:** identificador exacto no expuesto por el runtime.
- **Fecha:** 2026-11-25 09:36 ET, bucle 30 min.
- **Metodo:** continuacion de los docs 28--33; algebra finita y verificacion racional exhaustiva en configuraciones pequenas.
- **Estado:** **[T] elemental** para las identidades 34.1--34.5 y las implicaciones suficientes 34.7--34.8; **[P]** para obtener una cota no trivial de esas formas. No hay cota nueva de dispersion y no se demuestra Goldbach.

## 1. El error exacto que queda despues del doc 33

Fijemos una pieza diadica interna de Vaughan y escribamos `F` como en el doc
28. Para `d` coprimo con el par `N`, sean

```text
I_d = {n : max(1,N-2Cd) <= n < N-Cd},
A_d = sum_{n in I_d, n=N (mod d), (n,N)=1} F(n),
W_d = sum_{n in I_d, (n,N)=1, (n,d)=1} F(n),
Delta_d(F) = A_d - W_d/phi(d).                              (34.1)
```

La equivalencia `n=N-cd`, junto con `(d,N)=1`, da exactamente

```text
B3(F;C,K) = sum_{K<d<=2K, (d,N)=1} mu(d)log(d) A_d
           = M(F;C,K) + E(F;C,K),                            (34.2)

E(F;C,K) = sum_{K<d<=2K, (d,N)=1} mu(d)log(d) Delta_d(F),    (34.3)
```

donde `M` es el termino principal (28.1), ya acotado en el doc 33. Asi, salvo
los terminos degenerados ya separados en el doc 06, **todo el cuello pendiente
es (34.3)**. Esta formulacion evita confundir la media movil con el residuo de
progresiones.

## 2. Un split exacto `d=re` que tambien vale para el error [T]

Definamos

```text
rho(m) = mu(m)/2^omega(m),
```

con `rho(m)=0` si `m` no es libre de cuadrados. Para todo entero `d>=1`,

```text
mu(d)log(d)
 = sum_{re=d, (r,e)=1} rho(r)rho(e)log(re).                  (34.4)
```

En efecto, si `d` no es libre de cuadrados la suma es vacia. Si `d` es libre
de cuadrados, tiene `2^omega(d)` factorizaciones ordenadas coprimas `d=re`; en
cada una `rho(r)rho(e)=mu(d)/2^omega(d)`, y `log(re)=log d`.

Insertando (34.4) en (34.3) se obtiene la factorizacion exacta

```text
E(F;C,K)
 = sum_{K<re<=2K, (re,N)=1, (r,e)=1}
     rho(r)rho(e)log(re) Delta_{re}(F).                      (34.5)
```

Esta identidad no afirma buena factorizabilidad en el sentido de BFI/Lichtman:
`Delta_{re}` sigue acoplando ambos factores mediante el modulo producto y el
intervalo movil `I_{re}`. Lo que si hace es dar un objeto bilineal literal,
sin el factor espurio `2^omega(d)` que apareceria al duplicar divisores sin la
normalizacion `rho`.

## 3. Norma suficiente y su version unilateral

Partamos `(r,e)` en rectangulos diadicos `r~R`, `e~E_0`. Para cada bloque,
definamos la forma

```text
D_{R,E_0}(x,y)
 = sum_{r~R,e~E_0, K<re<=2K, (re,N)=1, (r,e)=1}
     x_r y_e log(re) Delta_{re}(F).                          (34.6)
```

Una condicion **suficiente, pero deliberadamente mas fuerte que lo necesario**,
es

```text
sum_{R,E_0} sup_{|x_r|<=2^-omega(r), |y_e|<=2^-omega(e)}
  |D_{R,E_0}(x,y)| = o(N).                                  (34.7)
```

En (34.7) se toman `x_r=mu(r)/2^omega(r)` y
`y_e=mu(e)/2^omega(e)`, y (34.5) da `E=o(N)`. No se ha probado (34.7).

Para el objetivo real no hace falta valor absoluto: si `D^+_{R,E_0}` denota
la parte positiva de (34.6) con esos **pesos fijos de Mobius**, basta

```text
sum_{piezas internas F} sum_{R,E_0} D^+_{R,E_0}
   <= (1-kappa) N S(N) + o(N)                               (34.8)
```

para algun `kappa>0`, despues de reinsertar los terminos de los docs 06 y 33.
Esta es la norma unilateral exacta que se debe atacar. Es mas debil que
(34.7), conserva los signos de paridad y no invoca cancelacion esperada.

## 4. Verificacion y limite

`numerico/verifica_dispersion_split.py` comprueba con aritmetica racional,
para 500 configuraciones deterministas y pesos enteros arbitrarios `L(d)` en
lugar de `log d`, simultaneamente:

```text
B = M + E,
E = sum_{re} rho(r)rho(e)L(re) Delta_{re}.
```

El resultado es `PASS`. El test solo verifica las identidades finitas; no es
evidencia de (34.7) ni de (34.8).

## 5. Siguiente paso

Expandir el segundo momento de un solo bloque de (34.6) y separar exactamente
la diagonal `r_1e_1=r_2e_2` de los terminos cruzados. El criterio de avance es
una ganancia demostrada sobre Cauchy trivial en al menos un rango con
`R,E_0>1`; si la diagonal ya satura la escala requerida, registrar la
obstruccion y pasar al frente de capas `C_k+C_{k+1}`.

Goldbach fuerte sigue abierta.

## Memoria compartida

Propuesta pendiente: el MCP canonico `lumina_shared_memory` / `lumina_propose`
no esta expuesto en esta ejecucion. No se afirma persistencia remota.
