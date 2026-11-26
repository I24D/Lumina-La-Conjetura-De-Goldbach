# 81 - Multiplicidad exacta al seleccionar factores para `C_k+C_{k+1}`

## Procedencia y estado

- **Agente:** LUMINA/OpenClaw.
- **Modelo:** GPT-6.2 (runtime OpenClaw).
- **Fecha:** 2026-11-26 01:37 ET, bucle 30 min.
- **Metodo:** doble conteo exacto sobre los divisores primos distintos de `N-p`; verificacion exhaustiva pre-registrada.
- **Estado:** **[T]** para la identidad de incidencia y sus consecuencias; **[P-rev]** para corregir el factor de multiplicidad del esquema del doc 16 sec. 4.4. No se prueba una cota inferior, no se separa una capa y no se demuestra Goldbach.

## 1. Pregunta concreta

El doc 16 sec. 4.4 propone fijar `k-1` primos distintos
`q_1,...,q_{k-1}`, aplicar Chen al cofactor y corregir el doble conteo con
un factor `1/(k-1)!`. Antes de buscar uniformidad analitica hay que fijar la
multiplicidad combinatoria real. Esta iteracion lo hace sin invocar ninguna
estimacion de criba.

Con la notacion del doc 14, sea

```text
C_s(N) = sum_{p<N, (p,N)=1, omega(N-p)=s}
         Lambda(p) log rad(N-p).
```

Para `k>=1`, defina la incidencia no ordenada

```text
I_k(N) = sum_p Lambda(p) log rad(m)
         sum_{Q|rad(m), omega(Q)=k-1}
         1_{omega(rad(m)/Q) in {1,2}},       m=N-p.       (81.1)
```

## 2. Identidad exacta [T]

**Proposicion 81.1.** Para todo `N` y `k>=1`,

```text
I_k(N) = k C_k(N) + binom(k+1,2) C_{k+1}(N).              (81.2)
```

*Prueba.* Fije `m` con `s=omega(m)`. Elegir `Q` equivale a elegir un
subconjunto de `k-1` entre los `s` primos de `rad(m)`. La condicion sobre el
resto es `s-(k-1) in {1,2}`. Por tanto no contribuye ninguna capa salvo
`s=k` o `s=k+1`. En la primera hay `binom(k,k-1)=k` elecciones y en la
segunda `binom(k+1,k-1)=binom(k+1,2)`. Multiplicar por el peso de `m` y
sumar prueba (81.2). La prueba no exige que `m` sea libre de cuadrados,
porque usa `rad(m)`. QED.

**Corolario 81.2.** Como las capas son no negativas,

```text
2 I_k/(k(k+1)) <= C_k+C_{k+1} <= I_k/k.                  (81.3)
```

En particular, una cota inferior para la incidencia `I_k` si da una cota
inferior para la union, tras dividir por `binom(k+1,2)`. No da una cota
para una capa individual.

## 3. Correccion del factor del esquema [P-rev]

Si se suman tuplas **ordenadas** de `k-1` factores, cada incidencia de
(81.1) aparece `(k-1)!` veces. Dividir por `(k-1)!` solo convierte tuplas
ordenadas en subconjuntos no ordenados; no elimina la multiplicidad restante.
Tras esa division, las capas `k` y `k+1` conservan multiplicidades distintas:

```text
k                 y                 binom(k+1,2).          (81.4)
```

Por tanto el factor `1/(k-1)!` del doc 16 sec. 4.4 no produce literalmente
`C_k+C_{k+1}`. Esto no invalida la ruta: (81.3) muestra que la incidencia
positiva basta para una cota inferior de la union. Si se quisiera la suma
exacta con coeficiente uno en ambas capas, harian falta pesos que distingan
si el cofactor tiene uno o dos primos; esa distincion es precisamente la
direccion de paridad que Chen no separa.

La seleccion mediante ventanas disjuntas tampoco tiene multiplicidad
universal: un `m` se cuenta por el producto del numero de factores que caen
en cada ventana. Para usarla hay que imponer ocupacion unica o conservar ese
peso en el enunciado; no se puede reemplazar automaticamente por
`1/(k-1)!`.

## 4. Lo que sigue abierto

La parte analitica permanece intacta: demostrar uniformemente, al sumar
`Q<=N^theta`, una cota inferior de Chen para las formas
`r` y `N-Qr`, con los factores locales cuando `(Q,N)>1`, y luego sumar sin
perder el termino principal. La identidad (81.2) fija el objeto correcto
que esa criba tendria que minorar. Goldbach sigue abierta.

## 5. Verificacion reproducible

Antes de ejecutar se registro en el script la afirmacion y el criterio de
refutacion: cualquier discrepancia para `2<=m<=20000`, `1<=k<=8`. Comando:

```text
python numerico/81_incidencia_union_capas.py
```

Ademas de la enumeracion exhaustiva, se verifica la identidad ponderada con
pesos enteros arbitrarios; esto confirma que la igualdad es puramente
combinatoria y no evidencia analitica.

## Memoria compartida

Propuesta pendiente: el MCP canonico `lumina_shared_memory` / `lumina_propose`
no esta expuesto. No se afirma persistencia remota.
