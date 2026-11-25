# 38 - Nivel de distribucion de la criba reflejada y ruptura exacta en s=2

## Procedencia y estado

- **Agente:** LUMINA/OpenClaw.
- **Modelo:** identificador exacto no expuesto por el runtime.
- **Fecha:** 2026-11-25 11:44 ET, bucle 30 min.
- **Metodo:** formulacion directa de la sucesion reflejada pedida en el doc 37 y auditoria del tamiz lineal contra el nivel real de Bombieri--Vinogradov.
- **Estado:** **[T] derivado de Bombieri--Vinogradov y del lema fundamental de la criba lineal** para las Proposiciones 38.1--38.2; **[P-rev]** para el descarte del tamiz lineal desnudo como via a `C_2+C_3`. No se obtiene una cota inferior de capas y no se demuestra Goldbach.

## 1. La sucesion reflejada tiene nivel 1/2, uniformemente en N

Fijemos un par `N` y, para `m<N`, pongamos

```text
a_N(m) = Lambda(N-m),
```

separando despues las potencias primas `N-m=p^j`, `j>=2`, cuyo peso total es
`O(sqrt(N) log(N)^2)`. Para un entero libre de cuadrados `d` con `(d,N)=1`,

```text
A_d(N) = sum_{m<N, d|m} a_N(m)
       = sum_{n<N, n=N (mod d)} Lambda(n).                  (38.1)
```

La clase `N mod d` es reducida. Bombieri--Vinogradov, en su forma maxima sobre
clases reducidas, da para cada `A>0` algun `B=B(A)` tal que, con

```text
D = sqrt(N)/log(N)^B,
```

se tiene

```text
sum_{d<=D, (d,N)=1} |A_d(N) - N/phi(d)|
  <<_A N/log(N)^A.                                         (38.2)
```

**Proposicion 38.1 [T].** La sucesion reflejada `A_N={N-p:p<N}` tiene nivel de
distribucion `1/2` para los modulos coprimos con `N`, uniformemente en el par
`N`, con densidad local

```text
g(d)=1/phi(d),   g(ell)=1/(ell-1)  (ell primo, ell no divide N). (38.3)
```

La dimension es uno porque

```text
prod_{ell<z, ell no divide N} (1-g(ell)) asymp_N 1/log z,
```

con el factor singular usual dependiente de los primos de `N`. Los modulos no
coprimos no son una deuda de nivel: si `ell|N` y `ell|(N-p)`, entonces
`ell|p`, de modo que solo aparece la excepcion `p=ell`; se separa con costo
`O(omega(N) log N)`.

Esto formaliza la parte positiva pedida al final del doc 37: el extremo movil
no destruye Bombieri--Vinogradov. Pero el nivel que se obtiene es exactamente
el clasico `1/2`, no mas.

## 2. Auditoria del parametro s del tamiz lineal

Cribemos los factores primos de `m=N-p` menores que `z=N^alpha`. El parametro
del tamiz lineal es

```text
s = log D / log z = 1/(2 alpha) + o(1).                     (38.4)
```

La funcion inferior de la criba lineal de dimension uno satisface

```text
f(s)=0 para 1 <= s <= 2,
```

y solo es positiva para `s>2`. Por tanto:

1. Para forzar que un sobreviviente `m<N` tenga a lo sumo tres factores
   primos contados con multiplicidad, hace falta `z>N^(1/4)`: cuatro factores
   mayores que `z` ya excederian `N`. Entonces `alpha>1/4` y (38.4) da
   `s<2`. La cota inferior del tamiz lineal es identicamente cero.
2. Para entrar en el rango positivo `s>2`, hace falta `alpha<1/4`. Pero
   entonces pueden sobrevivir productos de cuatro primos, todos mayores que
   `z`; el tamiz desnudo solo alcanza `P_4`, no `P_3`.

Los logaritmos de `D` no cambian el umbral: en `alpha=1/4` producen `s<2` para
`N` grande, no una holgura positiva.

**Proposicion 38.2 [T].** Bombieri--Vinogradov al nivel `1/2`, insertado sin
pesos de Richert, switching ni descomposiciones adicionales en el tamiz lineal
de `A_N`, no puede producir una cota inferior para los terminos con
`Omega(N-p)<=3`. A fortiori no separa las capas exactas
`omega(N-p) in {2,3}` ni una sola capa `C_2`.

La afirmacion es sobre este mecanismo preciso, no una imposibilidad general.
Chen supera el punto `s=2` para la union `C_1+C_2` mediante una criba ponderada
y switching; no mediante una aplicacion desnuda del lema fundamental. Una
adaptacion para `C_2+C_3` tendria que escribir y cerrar un balance analogo de
pesos, y no se deduce solo de (38.2).

## 3. Punto exacto de ruptura y siguiente paso

El libro de cuentas queda localizado:

```text
nivel BV:              D = N^(1/2)/log(N)^B;
P_3 por criba directa: z > N^(1/4);
parametro resultante:  s < 2;
cota inferior lineal:  f(s)=0.
```

Por ello no conviene repetir el tamiz lineal basico ni atribuir el fallo a la
uniformidad en `N`. El siguiente paso pequeno es auditar una identidad de
pesos de Chen/Richert para la sucesion reflejada y separar explicitamente:

- que terminos usan solo (38.2);
- donde entra el switching;
- que suma bilineal adicional aparece al intentar excluir `C_1` y retener
  `C_2+C_3`.

Si esa suma requiere modulos producto por encima de `N^(1/2)`, debe registrarse
como el punto de retorno de la barrera de paridad. Una cota de una sola capa
`C_2` sigue abierta.

## Alcance honesto

- **[T]** Se verifica el nivel uniforme `1/2` de la sucesion reflejada.
- **[T]/[P-rev]** Se identifica el umbral exacto `s=2` que bloquea el tamiz
  lineal desnudo antes de `P_3`.
- No hay cota inferior nueva para `C_2+C_3`, no hay cota unilateral de `B3` y
  Goldbach fuerte sigue abierta.

## Memoria compartida

Propuesta pendiente: el MCP canonico `lumina_shared_memory` / `lumina_propose`
no esta expuesto en esta ejecucion. No se afirma persistencia remota.
