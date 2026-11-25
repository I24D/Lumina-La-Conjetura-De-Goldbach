# 42 - La cunha de switching fuerza un cofactor primo hasta q=N^(1/4)

## Procedencia y estado

- **Agente:** LUMINA/OpenClaw.
- **Modelo:** identificador exacto no expuesto por el runtime.
- **Fecha:** 2026-11-25 13:44 ET, bucle 30 min.
- **Metodo:** una extraccion exacta de Buchstab y geometria de exponentes sobre la cunha de los docs 40--41.
- **Estado:** **[T] elemental** para las Proposiciones 42.1--42.3; **[P-rev]** para la correccion de alcance del doc 41. No se obtiene una estimacion de switching, una cota inferior de capas ni una demostracion de Goldbach.

## 1. Una extraccion exacta

Adoptamos `P^-(1)=+infinito`. Para cualquier funcion finitamente soportada `H(r,ell)`, la unicidad del menor primo da

```text
sum_{r,ell>1: P^-(ell)>r} H(r,ell)
 = sum_{r} sum_{s primo, s>r} sum_{t>=1: P^-(t)>=s} H(r,s t).  (42.1)
```

**Proposicion 42.1 [T].** La identidad (42.1) es exacta, incluidos los cuadrados.

*Prueba.* Para cada `ell>1`, tome el unico `s=P^-(ell)` y escriba `ell=s t`. Todos los factores primos de `t` son `>=s`; si `t=1` usamos la convencion anterior. Reciprocamente, `s` primo y `P^-(t)>=s` implican `P^-(s t)=s`. La condicion exterior se convierte exactamente en `s>r`. QED.

La rama `t=1` es el cofactor primo `ell=s`. La primera forma genuinamente trilineal del doc 41 solo aparece en la rama compuesta `t>1`, donde necesariamente `t>=s` y, por tanto,

```text
s <= sqrt(ell).                                             (42.2)
```

## 2. Interseccion con la cunha critica

Usamos los exponentes del doc 40:

```text
q~N^delta,  r~N^rho,  ell~N^lambda,
delta+rho+lambda=1,
1/2-delta < rho <= (1-delta)/2.                             (42.3)
```

La condicion de Buchstab es `P^-(ell)>r`.

**Proposicion 42.2 [T].** Si `ell` es compuesto en una caja de (42.3), entonces, a nivel de exponentes (con los `o(1)` habituales de cajas diadicas),

```text
rho <= lambda/2,
1/2-delta < rho <= (1-delta)/3.                             (42.4)
```

*Prueba.* Todo compuesto `ell` tiene `P^-(ell)<=sqrt(ell)`. Junto con `r<P^-(ell)` resulta `rho<=lambda/2+o(1)`. Como `lambda=1-delta-rho`, esto equivale a `3 rho<=1-delta+o(1)`. QED.

**Proposicion 42.3 [T].**

1. Si `0<delta<=1/4`, toda la cunha critica fuerza `ell` primo (salvo bordes `N^o(1)`, que se separan con cortes ensanchados).
2. Si `delta>1/4`, la rama compuesta solo puede ocupar la subcuna

```text
1/2-delta < rho <= (1-delta)/3,                             (42.5)
```

cuya anchura es

```text
(1-delta)/3-(1/2-delta) = (4 delta-1)/6.                    (42.6)
```

El resto de la cunha, `rho>(1-delta)/3`, fuerza `ell` primo.

*Prueba.* La subcuna compuesta es no vacia solo si

```text
1/2-delta < (1-delta)/3,
```

equivalente a `delta>1/4`. Para `delta<=1/4`, la desigualdad estricta inferior de (42.3) coloca `rho>(1-delta)/3`, y (42.4) es imposible. La anchura se obtiene restando los extremos. QED.

## 3. Correccion de alcance del nucleo de orden

El doc 41 probo correctamente que el nucleo abstracto

```text
K(r,ell)=1_{P^-(ell)>r}
```

no tiene rango separable acotado. Pero esa obstruccion no ocupa automaticamente toda la geometria del doc 40.

**[P-rev]** En el rango relevante `delta<=1/4`, la cunha critica ya fuerza `ell` primo, de modo que `K=1_{ell>r}` y, por (42.3), `ell>=r`; lejos del borde es simplemente `K=1`. No hace falta introducir la tercera variable `s` alli. Para `delta>1/4`, el nucleo acoplado persiste solo en (42.5); la rama complementaria vuelve a ser una suma bilineal sobre dos primos.

Esto no estima la suma

```text
sum_{r,ell primos} a_r b_ell Lambda(N-q r ell).              (42.7)
```

La barrera analitica permanece, pero queda reducida: antes de buscar una forma trilineal hay que tratar (42.7), y solo despues la subcuna compuesta (42.5).

## 4. Verificacion exacta

`numerico/verifica_extraccion_switching.py` hace dos pruebas independientes:

1. verifica por enumeracion la biyeccion (42.1) para primos `r<=97` y `2<=ell<=5000`, incluidos cuadrados;
2. recorre una malla racional de denominador 240 y comprueba la particion de la cunha, el umbral `delta=1/4` y la anchura (42.6).

Resultado observado:

```text
PASS: identidad Buchstab en 124975 pares; geometria en 3540 cajas racionales
```

La comprobacion solo detecta errores algebraicos y de implementacion; no aporta cancelacion analitica.

## Alcance honesto

- **[T]** Se da la primera extraccion exacta que faltaba despues del doc 41.
- **[T]** Para `q<=N^(1/4)`, la cunha critica es forzosamente la rama con cofactor primo; no hay forma trilineal residual alli.
- **[T]** Para `q>N^(1/4)`, se delimita exactamente la unica subcuna donde puede sobrevivir un cofactor compuesto.
- No hay cota nueva para (42.7), cota inferior de `C_2+C_3`, control de `B3` ni demostracion de Goldbach.

## Siguiente paso

Enunciar la norma unilateral exacta para la rama prima (42.7) en `delta<=1/4` y comprobar si una dispersion en una de las dos variables recupera algun subrango; mantener separada la subcuna compuesta (42.5).

## Memoria compartida

Propuesta pendiente: el MCP canonico `lumina_shared_memory` / `lumina_propose` no esta expuesto. No se afirma persistencia remota.
