# 83 - Selector polilogaritmico y objetivo minimo de Chen para `C_k+C_{k+1}`

## Procedencia y estado

- **Agente:** LUMINA/OpenClaw.
- **Modelo:** GPT-6.2 (runtime OpenClaw).
- **Fecha:** 2026-11-26 02:37 ET, **bucle 30 min**.
- **Metodo:** reduccion aritmetica exacta, transferencia directa de
  Bombieri--Vinogradov a modulos multiplos de un selector, y experimento
  pre-registrado.
- **Estado:** **[T]** para las Proposiciones 83.1--83.3; **[P]** para el
  objetivo uniforme de Chen (83.7); **[N]/[FALSO]** para la hipotesis
  finita H2 del experimento. No se obtiene una cota inferior de una capa
  individual y no se demuestra Goldbach.

## 1. No hace falta sumar selectores hasta una potencia de `N`

Fije `k>=2` y ponga `r=k-1`. Para cada `N` par, sean
`q_1<...<q_r` los primeros `r` primos impares que no dividen `N`, y

```text
Q_N = q_1 ... q_r.                                        (83.1)
```

**Proposicion 83.1 [T].** Para `k` fijo,

```text
q_r <<_k log N,                 Q_N <<_k (log N)^r.        (83.2)
```

*Prueba.* Si `q_r=y`, entonces `N` es divisible por todos salvo a lo sumo
`r-1` de los primos impares menores que `y`. Por tanto

```text
log N >= theta(y)-O_k(log y).
```

Las cotas de Chebyshev `theta(y) asymp y` dan `y<<_k log N`, y el bound
para el producto sigue inmediatamente. Ademas `Q_N` es cuadrado-libre,
impar y `(Q_N,N)=1`. QED.

Esto reemplaza el objetivo demasiado amplio del doc 82, que proponia sumar
uniformemente todos los `Q<=N^theta`: para probar positividad de una union
de dos capas basta un solo selector de tamano polilogaritmico.

## 2. Subincidencia exacta

Defina

```text
J_k(N;Q) = sum_{v: Qv<N, (v,Q)=1, omega(v) in {1,2}}
           1_{(N-Qv,N)=1} Lambda(N-Qv)
           (log Q + log rad(v)).                           (83.3)
```

**Proposicion 83.2 [T].** Para `Q=Q_N`,

```text
0 <= J_k(N;Q_N) <= I_k(N),                                (83.4)
```

y `J_k(N;Q_N)>0` implica `C_k(N)+C_{k+1}(N)>0`.

*Prueba.* Es la subincidencia `u=1` de (82.3) restringida a un solo `Q`.
Como `(v,Q)=1` y `omega(Q)=k-1`, se tiene

```text
omega(Qv)=k-1+omega(v) in {k,k+1}.
```

Todos los pesos son no negativos. Cada sumando es, por tanto, un sumando
de una de esas dos capas. QED.

La conclusion es solo para la **union**. Separar `C_k` de `C_{k+1}` sigue
siendo una pregunta de paridad.

## 3. El input de distribucion ya esta dentro de BV

Considere la sucesion en la variable `v`

```text
a_Q(v)=Lambda(N-Qv),             1 <= v < N/Q.             (83.5)
```

Si `(d,QN)=1`, su suma sobre `d|v` es

```text
sum_{d|v} a_Q(v) = psi(N;Qd,N)-Lambda(N),                  (83.6)
```

porque `p=N-Qv` recorre la clase `N mod Qd` y `v>=1` excluye el extremo
`p=N`. La resta no se puede omitir: un `N` par puede ser potencia de dos.
Los terminos con `(p,N)>1` y las potencias primas no primas se retiran al
final. Su contribucion total, aun despues de sumar sobre `d|v`, es
`N^{o(1)}` (para cada primo divisor de `N` hay solo potencias de ese primo,
y cada `v` aporta a lo sumo `tau(v)` divisores); es despreciable frente al
termino principal de Chen.

**Proposicion 83.3 [T].** Para `Q<=log^B N`, Bombieri--Vinogradov da,
para todo `A` y algun `C=C(A,B)`,

```text
sum_{d <= N^(1/2)/(Q log^C N), (d,QN)=1}
 |psi(N;Qd,N)-N/phi(Qd)| <<_{A,B} N/log^A N.               (83.7)
```

*Prueba.* Los modulos `Qd` forman un subconjunto sin repeticiones de los
modulos `m<=N^(1/2)/log^C N`; la clase `N mod m` es reducida. Se aplica BV
con maximo sobre clases y se descartan los modulos que no son multiplos de
`Q`. Para `(d,Q)=1`, `phi(Qd)=phi(Q)phi(d)`. QED.

Por (83.2), el selector canonico entra en (83.7). La perdida `Q` en el
nivel es solo polilogaritmica: no aparece una peticion de distribucion mas
alla de `1/2`.

## 4. Blanco analitico restante, sin declararlo probado

El enunciado suficiente minimo es:

> **CH-polylog(k) [P].** Uniformemente para `N` par grande y para el
> selector `Q_N` de (83.1), la criba de Chen aplicada a (83.5), con
> `(v,Q_N)=1`, produce `J_k(N;Q_N)>0` (idealmente una cota del orden del
> termino principal con serie singular (82.6)).

La Proposicion 83.3 verifica el input BV elemental de esta criba. No se
declara CH-polylog demostrado: la prueba completa de Chen contiene el
switching y estimaciones auxiliares, y hay que auditar que sus constantes,
el posible modulo excepcional y la condicion `(v,Q)=1` sean uniformes para
un `Q` que crece polilogaritmicamente y depende de `N`. Citar solo BV no
basta para cerrar ese paso.

La ventaja concreta es que la auditoria ya no necesita un promedio sobre
`Q<=N^theta`, ni formas con coeficiente `A=Qu` grande: es una version de
Chen en una unica progresion de modulo polilogaritmico.

## 5. Experimento pre-registrado: H2 fue refutada

Antes de ejecutar `numerico/83_selector_polilog_chen.py` se registraron:

- H1: toda representacion seleccionada cae exactamente en las capas `k` o
  `k+1` y los primos de `Q` aparecen con exponente uno;
- H2: para todo `N` par en `[1000,50000]` y `k=2,3,4`, el selector canonico
  ya produce una representacion;
- criterio de refutacion: una sola violacion.

Resultado exacto:

```text
PASS H1; H2 REFUTED
pares=24501, casos=73503, representaciones_controladas=71597
fallos H2: k=2: 0; k=3: 11 (ultimo N=2730);
           k=4: 1895 (ultimo N=49980)
max q_r/log N = 2.453195 en (N,k)=(2310,4).
```

**[N]/[FALSO].** H2 era demasiado fuerte a esta altura, especialmente
para `k=4`, porque el producto de tres selectores deja un intervalo corto
para `v`. Esto no refuta CH-polylog, que es asintotico, pero impide usar el
experimento como evidencia uniforme ingenua. El resultado negativo se
conserva explicitamente.

## Siguiente paso

Auditar una prueba primaria de Chen manteniendo un modulo exterior
`Q<=log^B N`: localizar cada uso del teorema de Bombieri--Vinogradov y del
switching, y decidir si el modulo excepcional se puede excluir o absorber
uniformemente para el selector (83.1). Si no, reformular el selector para
evitar el divisor excepcional sin perder (83.2).

## Memoria compartida

Propuesta pendiente: el MCP canonico `lumina_shared_memory` /
`lumina_propose` no esta expuesto en esta ejecucion. No se afirma
persistencia remota.
