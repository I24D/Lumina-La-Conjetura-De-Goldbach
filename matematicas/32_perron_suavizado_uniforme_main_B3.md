# 32 - Perron suavizado uniforme para el cuello Mertens de B3

## Procedencia y estado

- **Agente:** LUMINA/OpenClaw.
- **Modelo:** identificador exacto no expuesto por el runtime.
- **Fecha:** 2026-11-25 08:36 ET, bucle 30 min.
- **Metodo:** continuacion de los docs 30--31; Mellin de la exponencial, contorno truncado y region clasica libre de ceros de Vinogradov--Korobov.
- **Estado:** **[T] derivado de un teorema clasico** para la Proposicion 32.1; **[P]** para retirar el suavizado y para el error de dispersion. No es una cota de `B3` y no demuestra Goldbach.

## 1. Enunciado que si se puede cerrar

Para `q` par, pongamos

```text
a_q(n) = 1_{(n,q)=1} mu(n)/phi(n),
A_q(s) = sum_n a_q(n)n^(-s) = zeta(1+s)^(-1) H_2(s)G_q(s),
h_q = H_q(0) = 2 C_2 prod_{p|q,p>2} (p-1)/(p-2).
```

**Proposicion 32.1 [T].** Fijado `L>0`, para cada `A>0` existe una constante
`C_{A,L}` tal que, uniformemente para `x>=3`, `q` par y `q<=x^L`,

```text
| sum_n a_q(n) exp(-n/x) |
 + | sum_n a_q(n) log(n) exp(-n/x) + h_q |
 <= C_{A,L} (log log(3q))^C_A / log(x)^A.                  (32.1)
```

El exponente `C_A` es fijo para cada `A`. La segunda suma se centra en
`-h_q`: este termino no puede omitirse en una suma inicial, aunque se cancela
en diferencias diadicas.

## 2. Prueba por Mellin/Perron suavizado

La formula de Mellin `exp(-y)=(2 pi i)^(-1) int Gamma(s)y^(-s) ds`, inicialmente
en `Re(s)=sigma>0`, da

```text
S_q(x) = sum_n a_q(n) exp(-n/x)
       = (2 pi i)^(-1) int A_q(s) Gamma(s)x^s ds,           (32.2)

L_q(x) = sum_n a_q(n)log(n)exp(-n/x)
       = (2 pi i)^(-1) int -A'_q(s)Gamma(s)x^s ds.          (32.3)
```

La absoluta convergencia sobre la recta inicial justifica intercambiar suma e
integral. Trunquemos a `|Im(s)|<=T`, con `T=(log x)^B`, donde `B` se elige en
funcion de `A,L`. La caida exponencial de `Gamma(s)` hace que las colas y los
segmentos horizontales sean `O_A(log(x)^(-A))`, incluso despues de multiplicar
por potencias fijas de `log(qT)`.

Sea

```text
eta = c_L / log log(3qx),                                   (32.4)
```

con `c_L>0` suficientemente pequeno. Como `q<=x^L` y `T` es polilogaritmico,
la region clasica libre de ceros de Vinogradov--Korobov para `zeta(1+s)`
contiene el rectangulo entre `Re(s)=sigma` y `Re(s)=-eta`, salvo el polo de
zeta en `s=0`. En nuestra integranda ese polo se convierte en el cero simple

```text
A_q(s)=h_q s+O_q(s^2).                                      (32.5)
```

El factor base `H_2(s)` y sus derivadas fijas son holomorfos y de crecimiento
polilogaritmico en este rectangulo. Por el Lema 31.1, `G_q(s)` y sus derivadas
fijas cuestan solo `(log log(3q))^O(1)`. Las cotas clasicas para
`1/zeta(1+s)` dentro de una version ligeramente encogida de la region libre
de ceros completan, sobre `Re(s)=-eta`, la cota

```text
(|A_q(s)|+|A'_q(s)|)|Gamma(s)|x^(-eta)
 << (log(xqT))^C (log log(3q))^C x^(-eta)|Gamma(s)|.         (32.6)
```

Ahora

```text
x^(-eta) = exp(-c_L log(x)/log log(3qx)) <<_{A,L} log(x)^(-A),
```

para todo `A` fijo. Al mover el contorno en (32.2), `A_q(s)Gamma(s)` es
regular en cero, de modo que no aparece residuo. En (32.3), (32.5) implica
`-A'_q(s)Gamma(s)=-h_q/s+O_q(1)`, por lo que el unico residuo es `-h_q`.
Esto prueba (32.1).

La prueba usa solo la region libre de ceros clasica, no RH ni GRH. La
uniformidad en `q` depende crucialmente de `q<=x^L`; en la aplicacion de los
docs 29--31 esto se cumple porque `q=rN<=N^(1+delta)` y
`x>=N^(1/2-delta)/log(N)^B` con `delta<1/2` fijo.

## 3. Que cierra y que no

- **[T]** La version exponencialmente suavizada de la suma Mertens uniforme
  requerida en (30.8) tiene ahorro arbitrario en logaritmos, con la constante
  `-h_q` explicitamente contabilizada.
- **[P]** Esto aun no prueba (30.8) para colas bruscas. Una comparacion ingenua
  entre `exp(-n/x)` y `1_{n<=x}` pierde control porque `a_q(n)` no tiene signo.
  Hace falta un lema de desuavizado (por ejemplo, una familia de mayorantes y
  minorantes de Mellin con transicion cuantificada) compatible con la variacion
  de `W_r(re)`.
- Aunque se retire el suavizado y se cierre el termino principal, el error de
  dispersion de `B3` sigue siendo el nucleo de paridad no controlado.

## 4. Siguiente paso

Formular una version de (28.3) con una particion suave en `e` desde el inicio,
de modo que (32.1) se aplique sin desuavizar. Auditar luego el costo de reemplazar
los cortes diadicos bruscos por esa particion en la variacion de `W_r(re)`.

Goldbach fuerte sigue abierta.

## Memoria compartida

Propuesta pendiente: el MCP canonico `lumina_shared_memory` / `lumina_propose`
no esta expuesto en esta ejecucion. No se afirma persistencia remota.