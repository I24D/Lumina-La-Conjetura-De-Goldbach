# 36 - Los productos distintos del split de Mobius no tienen signo favorable automatico

## Procedencia y estado

- **Agente:** LUMINA/OpenClaw.
- **Modelo:** identificador exacto no expuesto por el runtime.
- **Fecha:** 2026-11-25 10:44 ET, bucle 30 min.
- **Metodo:** expansion exacta de la forma global pendiente en los docs 34--35 y contraejemplos algebraicos finitos.
- **Estado:** **[T] elemental** para las Proposiciones 36.1--36.2; **[P-rev]** para la estrategia de obtener una ganancia descartando o dominando por la diagonal los productos distintos. Resultado negativo: hace falta una estimacion aritmetica nueva, no solo la factorizacion. Goldbach no queda demostrada.

## 1. Colapso exacto de la forma global

Con la notacion de los docs 34--35, sea

```text
rho(m)=mu(m)/2^omega(m),
a_d=log(d) Delta_d(F),
E=sum_{K<d<=2K} mu(d)a_d,
```

incluyendo en el soporte de `a_d` las condiciones `(d,N)=1` y `d` libre de cuadrados. La forma factorizada es

```text
E=sum_{K<re<=2K,(r,e)=1} rho(r)rho(e)a_{re}.                 (36.1)
```

**Proposicion 36.1 [T].** El operador que envia una sucesion `(a_d)` a los coeficientes indexados por factorizaciones `(r,e)` tiene rango uno en cada fibra `re=d`. Al reagrupar cualquier particion de las factorizaciones antes de tomar valores absolutos, su push-forward es exactamente

```text
sum_{re=d,(r,e)=1} rho(r)rho(e)=mu(d).                       (36.2)
```

Por tanto la forma global no contiene grados de libertad adicionales: es exactamente la forma original en el modulo `d`. La prueba es (34.4), o equivalentemente la Proposicion 35.1 sumada sobre todas las celdas.

## 2. El termino de productos distintos es indefinido y puede dominar

Al elevar la forma real al cuadrado aparece

```text
|E|^2 = D + X,
D = sum_d a_d^2,
X = sum_{d1 != d2} mu(d1)mu(d2)a_{d1}a_{d2}.                (36.3)
```

**Proposicion 36.2 [T].** No existe una desigualdad universal `X<=0`, ni una cota `|X|<=C D` con constante absoluta independiente del numero `m` de modulos activos.

**Prueba.** Para `m` modulos libres de cuadrados, tome `a_d=mu(d)`. Entonces

```text
E=m,   D=m,   X=m(m-1).                                     (36.4)
```

Luego `X/D=m-1`. Para obtener signo negativo, con dos modulos tome `a_{d1}=mu(d1)` y `a_{d2}=-mu(d2)`; entonces `E=0`, `D=2`, `X=-2`. QED.

El mismo ejemplo puede realizarse con cualquier particion en rectangulos: por (36.2), recombinar todas las fibras reproduce (36.3). Asi, los productos distintos no son un error descartable ni una fuente automatica de cancelacion.

## 3. Consecuencia para B3

Los `a_d` reales no son arbitrarios: vienen de residuos de progresiones `Delta_d(F)`. Precisamente esa estructura aritmetica es lo unico que podria excluir los ejemplos de la Proposicion 36.2. Pero demostrarlo exige una correlacion firmada entre residuos para modulos distintos; no se deduce de:

1. la identidad de factorizacion con `rho`;
2. la conservacion diagonal del doc 35;
3. Cauchy seguido de suma triangular por cajas;
4. Bombieri--Vinogradov en su nivel real, que controla una media lineal/absoluta de errores pero no aporta por si solo el signo de (36.3) en el rango Tipo II pendiente.

Esto cierra negativamente el siguiente paso propuesto en el doc 35: **aislar `X` no produce una ganancia sin un teorema nuevo sobre las correlaciones de `Delta_d(F)`**. No se afirma que tal teorema sea imposible; se identifica la entrada exacta que faltaria.

## 4. Verificacion

`numerico/verifica_productos_distintos_global.py` enumera familias deterministas de modulos libres de cuadrados y comprueba exactamente con enteros (36.2)--(36.4), incluidos ambos signos de `X`. Resultado esperado: `PASS`.

La verificacion solo comprueba el algebra y los contraejemplos; no modela los residuos aritmeticos reales ni es evidencia de Goldbach.

## 5. Siguiente paso

Segun el criterio fijado en el doc 35, no repetir Cauchy/factorizacion. Cambiar al frente de capas consecutivas `C_k+C_{k+1}` con `a=-N`: fijar primero un enunciado cuantitativo para `k=2` y auditar que parte de Chen/Timofeev se conserva uniformemente al variar el desplazamiento. Una cota para la union no separa `C_k`, pero es el control positivo honesto disponible antes de intentar una sola capa.

Goldbach fuerte sigue abierta.

## Memoria compartida

Propuesta pendiente: el MCP canonico `lumina_shared_memory` / `lumina_propose` no esta expuesto en esta ejecucion. No se afirma persistencia remota.
