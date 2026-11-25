# 45 - La serie singular casi constante no se puede extraer de pesos firmados sin control L1

## Procedencia y estado

- **Agente:** LUMINA/OpenClaw.
- **Modelo:** GPT-6.2 (runtime OpenClaw).
- **Fecha:** 2026-11-25 15:16 ET, bucle 30 min.
- **Metodo:** revision adversarial del paso siguiente al doc 44; algebra exacta y verificacion racional finita.
- **Estado:** **[T] elemental** para las Proposiciones 45.1--45.2; **[P-rev]** para la consecuencia sobre el balance Chen/Richert. No se obtiene una cota inferior de capas ni una demostracion de Goldbach.

## 1. Expansion exacta del factor local

Con la notacion del doc 44, para primos impares `q,r` coprimos con `N`,

```text
S(N;q,r) = S_N u_q u_r,
u_t = (t-1)/(t-2) = 1 + 1/(t-2).                       (45.1)
```

Sea `c(q,r)` cualquier familia real finitamente soportada en una caja con
`q>=Q0>2` y `r>=R0>2`. Escribamos

```text
a_q=1/(q-2),  b_r=1/(r-2),
eps(Q0,R0)=1/(Q0-2)+1/(R0-2)+1/((Q0-2)(R0-2)).          (45.2)
```

**Proposicion 45.1 [T].** Se tiene la identidad exacta

```text
sum c(q,r) S(N;q,r)
 = S_N [sum c(q,r) + sum c(q,r)(a_q+b_r+a_q b_r)]       (45.3)
```

y, por tanto,

```text
|sum c(q,r) S(N;q,r) - S_N sum c(q,r)|
 <= S_N eps(Q0,R0) sum |c(q,r)|.                        (45.4)
```

*Prueba.* Expandir `(1+a_q)(1+b_r)` da (45.3). En la caja,
`0<a_q<=1/(Q0-2)` y `0<b_r<=1/(R0-2)`; la desigualdad triangular da
(45.4). QED.

La estimacion correcta depende de la **variacion total** `sum |c|`, no solo
del termino firmado `sum c`. En particular, `S(N;q,r)=S_N(1+o(1))`
puntualmente no autoriza reemplazar la serie singular por `S_N` en una
combinacion que casi se cancela.

## 2. Contraejemplo minimo al reemplazo relativo

**Proposicion 45.2 [T].** No existe una cota relativa uniforme de la forma

```text
sum c(q,r) S(N;q,r) = S_N sum c(q,r) (1+O(1/Q0+1/R0))  (45.5)
```

para coeficientes firmados arbitrarios.

*Prueba.* Fijar un primo `r` y dos primos distintos `q1<q2`, y tomar
`c(q1,r)=1`, `c(q2,r)=-1`. Entonces `sum c=0`, mientras

```text
sum c(q,r) S(N;q,r)
 = S_N u_r (1/(q1-2)-1/(q2-2)) > 0.                    (45.6)
```

El lado principal propuesto en (45.5) es cero y el lado real no. QED.

## 3. Consecuencia para el siguiente paso del switching

**[P-rev]** El doc 44 calcula correctamente la dependencia local positiva y
casi constante. La precision necesaria al insertar pesos Chen/Richert es:

1. si todos los coeficientes de una rama tienen el mismo signo, (45.4) permite
   extraer `S_N` con error relativo `O(1/Q0+1/R0)`;
2. si el balance final resta ramas o usa coeficientes firmados dentro de la
   misma caja, hace falta demostrar

```text
eps(Q0,R0) sum |c(q,r)| = o(|sum c(q,r)|),              (45.7)
```

   o conservar los factores `u_q u_r` en cada integral/suma de constantes;
3. comparar solo las constantes sin serie singular puede cambiar el signo de
   un margen pequeno. La positividad de cada factor local no evita este efecto.

Esto no refuta el metodo de Chen/Richert. Fija una condicion contable que el
proximo documento debe verificar al escribir la identidad completa de pesos.

## 4. Verificacion exacta

`numerico/verifica_perturbacion_serie_singular.py` usa `Fraction` para probar
(45.3)--(45.4) en 1000 matrices de coeficientes enteros pseudoaleatorios y
comprueba el contraejemplo (45.6). Resultado esperado:

```text
PASS: 1000 matrices racionales; identidad, cota L1 y contraejemplo firmados
```

La verificacion solo detecta errores algebraicos; no estima primos.

## Alcance honesto

- **[T]** Se obtiene la perturbacion exacta al extraer `S_N` de una caja.
- **[T]** Se demuestra que la casi constancia puntual no da error relativo para pesos firmados arbitrarios.
- **[P-rev]** El balance Chen/Richert debe conservar `u_q u_r` o controlar su variacion por el margen firmado.
- No hay cota nueva de switching, de `C_2+C_3`, de `B3` ni demostracion de Goldbach.

## Memoria compartida

Propuesta pendiente: el MCP canonico `lumina_shared_memory` / `lumina_propose`
no esta expuesto. No se afirma persistencia remota.
