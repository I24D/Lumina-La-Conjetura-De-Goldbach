# 35 - Conservacion de la diagonal en el split de Mobius

## Procedencia y estado

- **Agente:** LUMINA/OpenClaw.
- **Modelo:** identificador exacto no expuesto por el runtime.
- **Fecha:** 2026-11-25 10:14 ET, bucle 30 min.
- **Metodo:** expansion exacta del segundo momento sugerido en el doc 34 y verificacion finita independiente de la prueba.
- **Estado:** **[T] elemental** para las Proposiciones 35.1--35.2 y el Corolario 35.3; **[P]** para cualquier ganancia en los terminos de productos distintos. Resultado negativo: el split no reduce por si solo la diagonal. No demuestra Goldbach.

## 1. Coeficientes de una particion de las factorizaciones

Sea `d` libre de cuadrados y sea `P` una particion cualquiera de las
factorizaciones ordenadas coprimas `re=d` (por ejemplo, los rectangulos
diadicos del doc 34). Escribamos

```text
rho(m)=mu(m)/2^omega(m),
m_P(d)=#{(r,e) en P : re=d},
c_P(d)=sum_{(r,e) en P, re=d} rho(r)rho(e).
```

**Proposicion 35.1 [T].** Para cada celda `P`,

```text
c_P(d)=mu(d) m_P(d)/2^omega(d).                              (35.1)
```

En particular, todos los `c_P(d)` no nulos tienen el mismo signo y

```text
sum_P c_P(d)=mu(d),             sum_P |c_P(d)|=1.            (35.2)
```

**Prueba.** Como `d` es libre de cuadrados, una factorizacion coprima
ordenada corresponde a elegir un subconjunto de sus `omega(d)` primos; hay
`2^omega(d)` elecciones. En cada una,
`rho(r)rho(e)=mu(d)/2^omega(d)`. Sumar dentro de una celda da (35.1), y sumar
las multiplicidades de todas las celdas da `2^omega(d)`. QED.

Esto precisa el alcance de (34.4): la normalizacion elimina la multiplicidad,
pero no crea cancelacion entre las distintas factorizaciones de un mismo
modulo.

## 2. La diagonal correcta del segundo momento

Al expandir el segundo momento de una celda, la diagonal no es solamente
`r_1=r_2, e_1=e_2`. Incluye todas las colisiones

```text
r_1 e_1 = r_2 e_2 = d.                                      (35.3)
```

Si `a_d=log(d) Delta_d(F)`, la contribucion de producto igual en la celda
`P` es exactamente

```text
Q_P^=(F)=sum_d |c_P(d)|^2 |a_d|^2
         =sum_d m_P(d)^2/4^omega(d) |a_d|^2.                 (35.4)
```

**Proposicion 35.2 [T].** Para cualquier particion,

```text
sum_P sqrt(Q_P^=(F)) >= (sum_{d squarefree} |a_d|^2)^(1/2). (35.5)
```

**Prueba.** Es la desigualdad triangular en `ell^2`, aplicada a los vectores
no negativos `(|c_P(d)| |a_d|)_d`. Por (35.2), su suma coordenada es
`(|a_d|)_d`. QED.

La desigualdad puede ser estricta. Aunque una particion fina reduzca cada
energia cuadratica individual, al recombinar las cajas por desigualdad
triangular se paga al menos la norma diagonal original.

## 3. Consecuencia para la ruta del doc 34

**Corolario 35.3 [T].** El procedimiento

1. partir (34.5) en rectangulos `(R,E_0)`;
2. aplicar Cauchy en cada rectangulo;
3. sumar valores absolutos sobre los rectangulos,

no puede producir una ganancia automatica a partir de la diagonal: su costo
diagonal agregado es al menos

```text
(sum_{K<d<=2K} log(d)^2 |Delta_d(F)|^2)^(1/2).               (35.6)
```

Por tanto, cualquier avance de esta ruta debe demostrar cancelacion o una
ganancia estructural en los terminos de **productos distintos**
`r_1e_1 != r_2e_2`, o evitar la suma triangular de cajas mediante una forma
unilateral global. La factorizacion exacta con `rho` es util como interfaz,
pero no debilita por si misma la barrera diagonal descrita ya en el doc 06,
seccion 7.5.

Este es un bloqueo probado de la version mas directa del siguiente paso del
doc 34; no es una cota nueva para `B3`.

## 4. Verificacion

`numerico/verifica_energia_split_mobius.py` enumera todos los enteros
`1<=d<=20000`, particiona las factorizaciones coprimas por cajas diadicas y
comprueba exactamente con `Fraction` (35.1), (35.2) y la formula de energia
(35.4), usando valores deterministas no nulos para `a_d`. Resultado: `PASS`.
La prueba es la anterior; el calculo solo detecta errores de indices o
normalizacion.

## 5. Siguiente paso

No repetir Cauchy por rectangulos. Expandir una sola forma global y aislar el
primer termino de productos distintos. El criterio de avance sera una cota
demostrada para ese termino cruzado en algun rango no vacio; si vuelve a
requerir distribucion de nivel mayor que `1/2`, pasar al frente
`C_k+C_{k+1}`.

Goldbach fuerte sigue abierta.

## Memoria compartida

Propuesta pendiente: el MCP canonico `lumina_shared_memory` / `lumina_propose`
no esta expuesto en esta ejecucion. No se afirma persistencia remota.
