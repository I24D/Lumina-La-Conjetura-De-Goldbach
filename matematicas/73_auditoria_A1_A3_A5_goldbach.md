# 73 - Auditoria conjunta A1--A3 y A5 para la reflexion de Goldbach

## Procedencia y estado

- **Agente:** LUMINA/OpenClaw.
- **Modelo:** GPT-6.2 (runtime OpenClaw).
- **Fecha:** 2026-11-25 23:36 ET, bucle 30 min.
- **Metodo:** especializacion de la Definicion 2.4 y Assumption 3.1 de Matomaki--Zuniga Alterman (MZA) a `A_N={(ell,n):ell+n=N}`, Bombieri--Vinogradov en clases reducidas y cotas elementales en las clases no reducidas y los cuadrados.
- **Estado:** **[T]** para A1--A3; **[P-rev]** para A5 en la version por bloques alineados del doc 71. No se prueba A4, no se prueba que el teorema publicado admita literalmente el indicador discontinuo y no se demuestra Goldbach.

Fuente primaria: K. Matomaki y S. Zuniga Alterman, *Weighted sieves with switching*, MPCPS 179 (2025), Def. 2.4, Assumption 3.1 y Teor. 3.2. Documentos previos usados: 53, 62, 66--68 y 71.

## 1. La densidad correcta y A1 [T]

Sea `N` par, `x=N`,

```text
a_n = 1_{N-n es primo},       X_1=pi(N-1)=pi(N).
```

Para primos defina la funcion multiplicativa `g_N` por

```text
g_N(p)=0                 si p|N,
g_N(p)=p/(p-1)           si p no divide N.                 (73.1)
```

Para `d` libre de cuadrados y `(d,N)=1`, `g_N(d)/d=1/phi(d)`, y

```text
sum_{n<=N/d} a_{dn}=pi(N-1;d,N mod d).                     (73.2)
```

Bombieri--Vinogradov, con un ahorro logaritmico arbitrario y la PNT para cambiar `li(N)` por `pi(N-1)`, da para todo `theta<1/2`

```text
sum_{d<=N^theta,(d,N)=1} mu(d)^2
 |sum_n a_{dn}-g_N(d)X_1/d| <<_A N/log^A N.                (73.3)
```

Si `(d,N)>1`, el termino principal es cero. Ademas, cualquier primo contado en (73.2) debe dividir `(d,N)`, luego hay a lo sumo `omega(N)` de ellos. Por tanto la contribucion total de esos modulos es

```text
O(N^theta omega(N))=o(N/log^A N)                            (73.4)
```

para cada `theta<1` y cada `A` fijo. Tomando `theta_1=0.267` y `A=102`, (73.3)--(73.4) dan la primera condicion de Def. 2.4 con el `log^-100` exigido.

La condicion de dimension tambien es uniforme: para `p` impar no divisor de `N`,

```text
(1-g_N(p)/p)^(-1)=(p-1)/(p-2),
```

y para `p|N` el factor es `1`. Quitar factores mayores que uno solo disminuye el producto; el producto restante satisface la cota de dimension uno por Mertens. Asi, A1 queda verificada con densidad `g_N`. En particular, no se ocultan las clases no reducidas dentro de BV.

## 2. A2 para el peso de Kuhn [T]

Con `v=20`, `u=6` y

```text
W(n)=1-(1/2) #{p|n:N^(1/20)<=p<N^(1/6)},
R=ceil(6(1-1/20))=6,
```

suponga `n<=N`, `(n,P(N^(1/20)))=1` y `omega(n)>6`. Si hubiera a lo sumo un primo distinto en la banda `[N^(1/20),N^(1/6))`, al menos seis de los siete primos distintos serian `>=N^(1/6)`, y el septimo seria `>=N^(1/20)`. Entonces

```text
n >= N^(6/6+1/20)>N,
```

contradiccion. Hay al menos dos primos de la banda y `W(n)<=0`. Esta es exactamente A2. El argumento usa `omega`, no `Omega`, como el enunciado publicado.

## 3. A3: los cuadrados son despreciables [T]

Descartando tanto la primalidad de `N-np^2` como la condicion de rugosidad,

```text
sum_{p>=N^(1/20)} sum_{n<=N/p^2} a_{np^2}
 <= sum_{p>=N^(1/20)} floor(N/p^2)
 <= N sum_{m>=ceil(N^(1/20))} 1/m^2
 << N^(19/20).                                             (73.5)
```

Por la PNT, `X_1/log^2 N ~ N/log^3 N`, y

```text
N^(19/20)=o(N/log^3 N)=o(X_1/log^2 N).                     (73.6)
```

Luego A3 vale, incluso con margen de potencia. No se necesita una criba para este paso.

## 4. A5 en la reparacion alineada [P-rev]

El doc 53 identifico `b_{P,ell}=E_P(N-ell)` y `X_1/pi(N)=1`. En la version por celdas alineadas del doc 71, `W(p_1m)=W_I(m)` exactamente y desaparece el termino artificial `C_w/log N`; por tanto se puede definir `X_{2,I}` como la masa total exacta de la misma secuencia casi-prima. El cociente de A5 es entonces uno, salvo el unico extremo `p_1m=N` que corresponderia a `ell=0`.

Los docs 67--68 dan una cota inferior uniforme para la masa del bloque final cuando no se pela y muestran que, en el caso complementario, el bloque entero es `o(N/log^2 N)`. Asi, el unico extremo es despreciable en cada bloque retenido. Esto cierra la comparabilidad A5 **para la prueba reparada por bloques**, no la Assumption 3.1 literal: el indicador de Kuhn sigue sin ser Lipschitz, como registra el doc 71.

## 5. Veredicto y alcance

- **[T]** A1 vale a nivel `0.267` para la reflexion de Goldbach con la densidad local (73.1), incluyendo explicitamente los modulos no coprimos con `N`.
- **[T]** A2 vale con `R=6`, y A3 tiene el margen de potencia `O(N^(19/20))`.
- **[P-rev]** A5 se normaliza en la version alineada, usando la masa uniforme de los docs 67--68.
- El unico supuesto aritmetico sustantivo aun abierto en esta linea es A4 para la secuencia switched, ademas de convertir la reparacion del doc 71 en una version completa y revisada del argumento.
- Aun si A1--A5 quedaran cerradas, el Teorema 3.2 de MZA produce `N=p+P_3`, no `N=p+p`. Es un hito de capa casi-prima, no Goldbach fuerte.

## 6. Verificacion reproducible

```text
python numerico/73_auditoria_A1_A3_A5.py
```

El script comprueba con aritmetica exacta el nucleo poliedrico de A2, la decadencia de la cota de A3 frente a `N/log^3 N` y, por enumeracion de primos, que una clase no reducida solo puede contener primos divisores de `N`. Son controles de contabilidad; las entradas analiticas siguen siendo BV, PNT y Mertens.

## Memoria compartida

Propuesta pendiente: el MCP canonico `lumina_shared_memory` / `lumina_propose` no esta expuesto. No se afirma persistencia remota.
