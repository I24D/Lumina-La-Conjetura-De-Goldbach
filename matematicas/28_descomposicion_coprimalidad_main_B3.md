# 28 - Descomposicion exacta de la coprimalidad en el termino principal de B3

## Procedencia y estado

- **Agente:** LUMINA/OpenClaw.
- **Modelo:** identificador exacto no expuesto por el runtime.
- **Fecha:** 2026-11-25 06:36 ET, bucle 30 min.
- **Metodo:** derivacion algebraica, revision adversarial de los docs 17--18 y verificacion finita con pesos arbitrarios.
- **Estado:** **[T] elemental** para las identidades y la cota de variacion; **[P]** para la lectura analitica. No es una cota de Tipo II y no demuestra Goldbach.

## 1. El termino principal exacto de una pieza diadica

Fijemos una pieza diadica de las variables internas de Vaughan y escribamos

```text
F(n) = sum_{ab=n, a en A, b en B, a>V, b>U} Lambda(a) nu_U(b).
```

Para `d` coprimo con el par `N`, sea

```text
I_d = {n : N-2Cd <= n < N-Cd},
W_d = sum_{n in I_d, (n,N)=1, (n,d)=1} F(n).
```

La condicion `(n,N)=1` es exactamente la condicion exterior `(c,N)=1`, pues
`n=N-cd` y `(d,N)=1` implican `(n,N)=(c,N)`. El termino principal natural de
esta pieza (despues de reemplazar la clase reducida por su densidad) es

```text
M = sum_{K<d<=2K, (d,N)=1} mu(d) log(d) W_d / phi(d).       (28.1)
```

Esto precisa el peso movil esquematico del doc 18: ademas del intervalo movil,
`W_d` contiene la coprimalidad `(n,d)=1`. Por tanto no es correcto tratarlo
directamente como una funcion de variacion puramente geometrica.

## 2. Expansion exacta del factor de coprimalidad [T]

Definamos, para `r|d`,

```text
W_r(d) = sum_{n in I_d, (n,N)=1, r|n} F(n).
```

Por inclusion-exclusion,

```text
W_d = sum_{r|d} mu(r) W_r(d).                               (28.2)
```

El factor `mu(d)` en (28.1) permite restringir a `d` libre de cuadrados. Al
escribir `d=re`, se tiene `(r,e)=1`, `mu(d)mu(r)=mu(e)` y
`phi(d)=phi(r)phi(e)`. En consecuencia,

```text
M = sum_{r,e: K<re<=2K, (re,N)=1, (r,e)=1, r squarefree}
      mu(e) log(re) W_r(re) / (phi(r) phi(e)).               (28.3)
```

La condicion `r` libre de cuadrados es necesaria porque el cambio de variables
parte del soporte `mu(d)!=0`; `mu(e)` ya fuerza que `e` sea libre de cuadrados.
La igualdad es finita y exacta. Se verifico tambien en 500 configuraciones
deterministas con una funcion `F` arbitraria y un peso exterior arbitrario
`L(d)` en lugar de `log d`; ver
`numerico/verifica_main_B3_coprimalidad.py`.

## 3. La variacion geometrica si esta controlada [T]

Para `r` fijo, al aumentar `e` en una unidad los dos extremos de `I_{re}` se
desplazan monotonamente en `Cr` y `2Cr`. Cada entero cruza cada extremo a lo
sumo una vez. Por desigualdad triangular,

```text
sum_e |W_r(r(e+1))-W_r(re)|
 <= 2 sum_{n en la union de los I_re, (n,N)=1, r|n} |F(n)|. (28.4)
```

Asi, el intervalo movil por si solo no es el obstaculo de C10: despues de
expandir la coprimalidad, cada rebanada `r` tiene variacion total geometrica,
no un costo igual al numero de modulos.

## 4. Que falta: una suma uniforme en e y la cola de r grandes [P]

La sumacion parcial aplicada a (28.3) pide controlar, uniformemente en `r,N`,

```text
sum_{E<e<=y, (e,rN)=1} mu(e) log(re)/phi(e),                (28.5)
```

y luego sumar el resultado con el costo `1/phi(r)` y la masa de `F` sobre
multiplos de `r`. Para `r=1`, (28.4) convierte el problema en una estimacion
Mertens ponderada estandar. Para `r>1`, (28.3) muestra exactamente la deuda
que el argumento del doc 17 omitio. En particular, cuando `r` es comparable
con `K`, el intervalo de `e` es corto y no se obtiene cancelacion de Mobius
por una PNT ordinaria.

**Conclusion honesta [P]:** C10 no queda cerrado, pero ahora se separa en dos
subproblemas precisos: (i) sumas uniformes (28.5) para `r` pequeno, y (ii) una
cota de la cola de divisores `r` grandes usando la escasez de multiplos de `r`
en `F`. La movilidad de los extremos ya no es una barrera independiente; la
coprimalidad/divisor grande es el cuello de botella restante del termino
principal. Esto no controla el error de dispersion ni la caja B3 completa.

## 5. Siguiente paso

Fijar `R`, separar `r<=R` y `r>R` en (28.3), y buscar una eleccion de `R` para
la cual una cota Mertens uniforme cierre la primera parte y una cota de
divisores para `F=Lambda*nu_U` ahorre en la segunda. Si la cola grande es del
mismo orden que B3, registrar explicitamente la ruta como otra forma de la
barrera de paridad.

Goldbach fuerte sigue abierta.
