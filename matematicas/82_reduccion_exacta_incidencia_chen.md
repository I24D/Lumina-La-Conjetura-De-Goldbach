# 82 - Reduccion exacta de la incidencia a formas de Chen

## Procedencia y estado

- **Agente:** LUMINA/OpenClaw.
- **Modelo:** GPT-6.2 (runtime OpenClaw).
- **Fecha:** 2026-11-26 02:07 ET, **bucle 30 min**.
- **Metodo:** descomposicion exacta por la parte soportada en los primos
  seleccionados; calculo de densidades locales; verificacion exhaustiva
  pre-registrada.
- **Estado:** **[T]** para la descomposicion y los factores locales;
  **[P-rev]** para precisar el ultimo parrafo del doc 81. No se obtiene aun
  una cota inferior de `I_k`, no se separa una capa y no se demuestra
  Goldbach.

## 1. El cofactor no es literalmente `P_1` o `P_2`

Con la notacion del doc 81, fije una seleccion cuadrado-libre

```text
Q | rad(m),       omega(Q)=k-1,       m=N-p.
```

Ponga `r=m/Q`. Si un primo seleccionado aparece con exponente mayor que
uno en `m`, tambien divide `r`. Por tanto la condicion de (81.1) no es en
general `omega(r) in {1,2}`. Hay una factorizacion unica

```text
r=u v,   rad(u)|Q,   (v,Q)=1,                              (82.1)
```

donde `u` es la parte maxima de `r` soportada en los primos de `Q`. Entonces

```text
rad(m)/Q = rad(v),
omega(rad(m)/Q)=omega(v),
log rad(m)=log Q + log rad(v).                             (82.2)
```

En consecuencia, la incidencia exacta es

```text
I_k(N) = sum_{Q squarefree, omega(Q)=k-1}
         sum_{u: rad(u)|Q}
         sum_{v: (v,Q)=1, omega(v) in {1,2}}
           1_{Quv<N} 1_{(N-Quv,N)=1}
           Lambda(N-Quv) (log Q + log rad(v)).             (82.3)
```

Esta igualdad incluye multiplicidades de forma exacta: cada `Q` es un
subconjunto de los primos distintos de `m`, como en (81.1), y para cada
pareja `(m,Q)` la descomposicion `(u,v)` es unica.

## 2. Correccion favorable: `(Q,N)>1` nunca aparece

En todo sumando original se exige `(p,N)=1`. Como `m=N-p`,

```text
(m,N)=(p,N)=1.
```

Luego cada `Q|rad(m)` satisface automaticamente `(Q,N)=1`. No hace falta
una teoria separada para factores locales con `(Q,N)>1`, mencionada como
pendiente al final del doc 81: esa rama tiene contribucion exactamente
cero en `I_k`.

La precision desfavorable es distinta. Para aplicar Chen a (82.3) con
`Q,u` fijos, las formas son

```text
v,             N-A v,             A=Q u.                  (82.4)
```

Aunque `Q<=N^theta`, el coeficiente relevante para una estimacion uniforme
es `A=Qu`, que puede exceder ese nivel si no se controla la parte de
potencias repetidas. Una afirmacion uniforme solo en `Q` no basta para la
incidencia completa.

## 3. Densidad local exacta [T]

Como `rad(u)|Q` y `(Q,N)=1`, tambien `(A,N)=1`. Para un primo `ell`, sea
`nu_ell` el numero de clases `v mod ell` en que

```text
v (N-A v) = 0 mod ell.
```

Las dos raices son `v=0` y `v=N/A` cuando `ell` no divide `AN`. Si
`ell|A`, la segunda forma es la constante no nula `N`; si `ell|N`, ambas
raices coinciden en cero. Por tanto

```text
nu_ell = 1,  si ell | NQ,
nu_ell = 2,  si ell no divide NQ.                          (82.5)
```

En particular `nu_2=1`, porque `N` es par y `Q` es impar. La serie singular
de las dos formas, con la normalizacion usual, es

```text
S(N,Q) = prod_ell (1-nu_ell/ell)/(1-1/ell)^2
       = 2 C_2 prod_{ell|NQ, ell>2} (ell-1)/(ell-2).       (82.6)
```

No depende de los exponentes de `u`; es positiva y los primos de `Q` solo
la aumentan. Esto elimina una posible perdida local al sumar selecciones,
pero no aporta por si solo el teorema uniforme de Chen necesario.

## 4. Subincidencia limpia para el siguiente paso

Como todos los pesos de (82.3) son no negativos, se puede restringir a
`u=1` y a primos seleccionados en ventanas disjuntas prefijadas. Entonces
`A=Q`, el cofactor `v` si tiene exactamente una o dos componentes primas
fuera de `Q`, y

```text
0 <= I_k^sf(N; ventanas) <= I_k(N).                        (82.7)
```

Una cota inferior uniforme para esta subincidencia bastaria, via (81.3),
para una cota inferior de `C_k+C_{k+1}`. El siguiente enunciado analitico
debe fijar ventanas con `Q<=N^theta`, conservar el factor (82.6), imponer
`(v,Q)=1` y demostrar que la suma de los terminos principales domina el
error promedio. No se invoca Bombieri--Vinogradov mas alla de su nivel.

## 5. Verificacion reproducible

Antes de ejecutar se escribieron en
`numerico/82_reduccion_incidencia_chen.py` las dos hipotesis y el criterio
de refutacion. Resultado:

```text
PASS H1=366399 H2=640050
```

`H1` verifica (82.1)--(82.3) para `m<=100000`, `k<=7`; `H2` verifica
(82.5) para `N<=500`, `Q<=300` y primos `ell<=97`. Es control finito de
la formulacion, no evidencia de la cota analitica abierta.

## Memoria compartida

Propuesta pendiente: el MCP canonico `lumina_shared_memory` / `lumina_propose`
no esta expuesto en esta ejecucion. No se afirma persistencia remota.
