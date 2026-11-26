# 55 - Balance del peso nuevo del doc 54: positivo pero de segundo orden (Etapa 2)

## Procedencia y estado

- **Agente:** LUMINA/Muse.
- **Fecha:** 2026-11-25.
- **Metodo:** calculo exacto con aritmetica racional (`Fraction`) del balance
  `B = sum c(q,r) a_q b_r` (46.5); verificacion por la identidad de atomos
  (Prop 46.2). Marco: docs 46 (46.1)-(46.6), 52 (constante de Chen),
  54 (peso por atomos).
- **Estado:** **[T]** para el calculo de B y el veredicto de signo y escala.
  No hay cota inferior de capas, no hay nivel de distribucion, no hay
  demostracion de Goldbach.

## 1. Definicion y modelo

Por (46.5), con margenes cero (46.4),

```text
sum_{q,r} c(q,r) S(N;q,r) = S_N * B,   B = sum_{q,r} c(q,r) a_q b_r,
```

con `a_q = 1/(q-2)` (46.1). B es la suma que multiplica a la serie
singular (y al termino `S(N)*N/log^2 N`): es el balance del peso.

- Modelo: `N = 10^12` (doc 52, modelo A); primos representativos
  `q_i = menor primo >= N^{alpha_i}` (convencion del doc 52; el float
  solo elige el primo, toda la aritmetica despues es exacta).
  Malla del doc 54: `v = 10`, control `v = 8`; `h = 12, 24, 48, 96`.
- Diagonal `q = r` (caso `q^2 | n`, doc 54 sec. 1): pares ordenados,
  la diagonal entra UNA vez. Su parte `D = sum_i c(i,i) a_{q_i}^2`
  se reporta aparte; la identidad de atomos la incluye automaticamente.
- Verificacion interna: B directo == B por atomos exactamente en todos
  los casos; margenes cero re-verificados; `sum c = 0` exacto.
  `||c||_1` reproduce la tabla del doc 54 (0.1763, 0.2905, 0.1830, 0.2918).

## 2. Proposicion 55.1 [T]: B > 0 en los 16 casos

*Prueba.* Por (46.7), cada atomo aporta
`w_{i,j} (a_{q_i} - a_{q_{i+1}})(b_{r_j} - b_{r_{j+1}})`, con
`w_{i,j} = psi(alpha_i, alpha_j) > 0` y `a_q` estrictamente decreciente
en `q`. Los pasos degenerados (dos `alpha` seguidos con el mismo primo
representativo: 3 en `v=10,h=96`; 1 en `v=8,h=96`; hasta 26 a `N=10^7`)
aportan exactamente 0. Luego cada termino es `>= 0` y al menos uno es
`> 0`. QED.

Tabla a `N = 10^12` (B en decimal; valores exactos en la salida del script):

| v | h | B | B/\|\|c\|\|_1 | B/cota (46.6) | diag/B | B/0.67 | Q0 |
|---|---|---|---|---|---|---|---|
| 10 | 12 | 4.291e-07 | 2.434e-06 | 3.703e-03 | 0.491 | 6.40e-07 | 41 |
| 10 | 24 | 7.461e-07 | 3.190e-06 | 2.325e-03 | 0.368 | 1.11e-06 | 29 |
| 10 | 48 | 1.044e-06 | 3.865e-06 | 1.705e-03 | 0.318 | 1.56e-06 | 23 |
| 10 | 96 | 1.404e-06 | 4.835e-06 | 1.397e-03 | 0.302 | 2.09e-06 | 19 |
| 10 | 192* | 1.651e-06 | 5.481e-06 | 1.233e-03 | 0.300 | 2.46e-06 | 17 |
| 8 | 12 | 2.313e-07 | 1.264e-06 | 6.016e-03 | 0.546 | 3.45e-07 | 71 |
| 8 | 24 | 3.729e-07 | 1.565e-06 | 4.070e-03 | 0.426 | 5.56e-07 | 53 |
| 8 | 48 | 5.470e-07 | 2.007e-06 | 3.052e-03 | 0.391 | 8.16e-07 | 41 |
| 8 | 96 | 6.474e-07 | 2.219e-06 | 2.718e-03 | 0.361 | 9.65e-07 | 37 |
| 8 | 192* | 6.629e-07 | 2.195e-06 | 2.689e-03 | 0.329 | 9.88e-07 | 37 |

(* control de convergencia.) Referencia de Chen [L] (doc 52 sec. 1):
`0.6706 = 2.6408 - 3.9404/2`.

## 3. Limite continuo y decaimiento con N [P/T]

- B crece con `h` hacia un limite continuo `B_inf`: medido `h=192`,
  `v=10: 1.651e-06`, `v=8: 6.629e-07`. Estimacion integral del limite
  (la suma de atomos tiende a `(log N)^2` por la integral doble de
  `psi(alpha,beta) N^{-(alpha+beta)}`): `~1.6e-06` (v=10) y `~7.5e-07`
  (v=8) a `N=10^12` — concuerda en escala con lo medido.
- B decae con N: a `N=10^7`, `B(v=10,h=96) = 3.689e-05`; a `N=10^12`,
  `1.404e-06` (factor x26). En los 8 casos el factor de `N=10^7` a
  `N=10^12` es x26-64. La esquina `(alpha,beta ~ 1/v)` domina y
  `Q0 ~ N^{1/v}` crece, asi que B se achica con N: no es una constante
  final, depende de N.

## 4. Veredicto

- **Signo: B > 0 en todos los casos.** El balance es positivo. La
  diagonal aporta el 30-55% de B a `N=10^12` (contada una vez; si se
  duplicara por error, B crecería en D, sin cambiar el orden).
- **Escala: B es de segundo orden.** `B/0.67 ~ 10^{-6}` a `N=10^12`
  (maximo 2.09e-06 en la tabla; limite continuo ~2.5e-06). B usa solo
  el 0.14-0.60% de la cota (46.6). B decae con N. Esta seis ordenes de
  magnitud por debajo de la constante `0.6706` de Chen.
- **Razon estructural [T]:** (46.4) fuerza `sum c = 0` (doc 46 (46.3)),
  asi que ningun peso de margenes exactamente cero puede llevar un
  termino principal de orden 1 por (46.5). El diseno del doc 54 elimina
  el error de primer orden y, con el, el termino principal: el balance
  positivo que queda es de segundo orden por construccion.
- **Respuesta:** el peso nuevo conserva un balance positivo en signo,
  pero NO a escala de Chen. No hay constante final comparable con 0.67
  en esta etapa.

## 5. Que ajuste lo podria salvar

- (a) **Mas peso en la esquina:** `psi = (alpha*beta) + (alpha*beta)^2`
  (dos monomios no proporcionales, mismo rango 2): B medido x34-47 mayor
  (`v=10,h=96: 6.570e-05`; `B/0.67 = 9.80e-05`; `B/cota = 1.58e-02`).
  Mejora la constante, no el orden: sigue segundo orden y << 0.67.
  (De esta variante solo se midio B con el mismo script; sus
  propiedades 2D no se re-verificaron aqui.)
- (b) **El ajuste estructural:** margenes PEQUENOS pero no cero.
  Maximizar `sum c` sujeto a `|sum a_q R_q + sum b_r C_r| << sum c`,
  con el trade cuantificado por (46.3). Es el siguiente diseno: un peso
  con `sum c` de orden 1 y defecto marginal controlado, en vez de
  margenes exactamente cero. (Doc 56 candidato.)
- Ademas sigue pendiente el nivel de distribucion de este peso
  (doc 54 sec. 6) y la integracion con los pesos analiticos del Lemma 4:
  sin eso, ni siquiera un B de orden 1 seria una constante de Goldbach.

## 6. Verificacion reproducible

```text
python numerico/balance_peso_nuevo_doc54.py
```

PASS en el PC: B > 0 exacto en los 16 casos + 2 de convergencia;
B directo == B por atomos; margenes cero; `sum c = 0`; diagonal una vez.

## Alcance honesto

- Lo que este doc SI dice: B calculado exactamente, positivo, de segundo
  orden, con constante explicita (`B_inf ~ 1.65e-06` a `N=10^12, v=10`)
  y comparacion con el 0.67 de Chen. El peso del doc 54 tiene senal
  positiva a escala principal en el sentido de (46.5), pero esa escala
  es de segundo orden.
- Lo que NO dice: no hay cota inferior de capas, no hay nivel de
  distribucion, no hay demostracion de Goldbach.

## Memoria compartida

Propuesta pendiente: el MCP canonico `lumina_shared_memory` /
`lumina_propose` no esta expuesto. No se afirma persistencia remota.
