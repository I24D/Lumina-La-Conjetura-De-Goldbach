# 54 - Atomos rectangulares de un peso nuevo con margenes cero por construccion (Etapa 1)

## Procedencia y estado

- **Agente:** LUMINA/Muse.
- **Fecha:** 2026-11-25 ~18:52 ET. Reclamo de ESTADO.md 18:36 ET (este paso; el bucle no lo duplica).
- **Metodo:** construccion directa por atomos rectangulares; verificacion exacta con aritmetica racional (`Fraction`).
- **Estado:** **[T]** para la construccion y la prueba de margenes cero; **[T]** para la verificacion numerica (criterio (46.4) SATISFECHO por primera vez en la serie). No hay balance ni cota inferior: eso es la etapa 2. No hay demostracion de Goldbach.

## 1. El dominio: el cuadrado completo, no la cuna

- Pares ordenados `(q,r)` de primos con tamanos logaritmicos `alpha, beta` en `[1/v, 1/2]`, `v = 10` (control `v = 8`). Malla del doc 52: `alpha_t = lo + (t+1)*step`, `lo = 1/v`, `hi = 1/2`, `step = (hi-lo)/(h+1)`, `t = 0..h-1`; `h = 12, 24, 48, 96`.
- Por que el cuadrado: la cuna triangular del doc 52 (sec. 6.2) es lo que impedia los atomos rectangulares en el peso de Chen. El peso nuevo usa pares ordenados sobre el cuadrado completo; el factor local (46.1) es simetrico en `(q,r)`, asi que la suma sobre pares ordenados esta bien definida. La diagonal `q = r` es el caso `q^2 | n`, que (46.1) admite con `a_q = b_r`; la etapa 2 lo contara con su multiplicidad correcta.
- Sin corte `u`: el peso nuevo no tiene truncacion P2/P3; su dominio es todo el cuadrado.

## 2. La construccion

- Atomo rectangular `A_{i,j}` en la celda `(i,j)`, `0 <= i,j <= h-2`: `+1` en `(i,j)`, `-1` en `(i,j+1)`, `-1` en `(i+1,j)`, `+1` en `(i+1,j+1)`. Es el atomo del doc 46 sec. 2, en celdas adyacentes.
- Funcion potencial `psi(alpha,beta) = (alpha*beta)^2 + (alpha*beta)^3 = alpha^2 beta^2 + alpha^3 beta^3`. Simetrica, suave, positiva en el cuadrado. Su parcial mixta es `d^2 psi/d alpha d beta = 4 alpha beta + 9 alpha^2 beta^2 = alpha beta (4 + 9 alpha beta)`, genuinamente bidimensional: no es aditiva (la parcial mixta no es cero) ni separable (si `alpha beta (4+9 alpha beta) = f(alpha) g(beta)`, en `alpha = beta = 0` daria `16/13 = 4`, contradiccion).
- Coeficiente de atomo `w_{i,j} = psi(alpha_i, alpha_j)`. Peso: `c = suma_{i,j} w_{i,j} A_{i,j}`, es decir `c(i,j) = w_{i,j} - w_{i-1,j} - w_{i,j-1} + w_{i-1,j-1}` (con `w = 0` fuera de rango).
- Formula interior exacta (`1 <= i,j <= h-2`, con `s = step`):

```text
c(i,j) = s^2 [(2 alpha_i - s)(2 alpha_j - s)
            + (3 alpha_i^2 - 3 alpha_i s + s^2)(3 alpha_j^2 - 3 alpha_j s + s^2)].
```

Es la diferencia mixta discreta de cada monomio: `Delta(alpha^2 beta^2) = s^2 (2 alpha_i - s)(2 alpha_j - s)` y `Delta(alpha^3 beta^3) = s^2 (3 alpha_i^2 - 3 alpha_i s + s^2)(3 alpha_j^2 - 3 alpha_j s + s^2)`.
- Signo: `c(i,j) > 0` en el interior (ambos factores positivos); las celdas del borde compensan para que los margenes sean cero (ver conteos en sec. 5).

## 3. Prueba de margenes cero

**Proposicion 54.1 [T].** Todos los margenes de `c` son cero: `R_i = 0` para todo `i` y `C_j = 0` para todo `j`. Luego el criterio (46.4) se satisface.

*Prueba.* Cada atomo `A_{i,j}` tiene margenes cero: en la fila `i` aporta `+w` en la columna `j` y `-w` en la columna `j+1`; en la fila `i+1` aporta `-w` en `j` y `+w` en `j+1`; en las demas filas, `0`. Lo mismo por columnas. Una suma de atomos hereda margenes cero. QED.

Consecuencia (doc 46 (46.5)-(46.6)): `sum c(q,r) S(N;q,r) = S_N sum c(q,r) a_q b_r`, con `|sum c S| <= S_N / ((Q0-2)(R0-2)) * sum |c|`. El error de primer orden desaparece por construccion.

## 4. Genuinamente bidimensional (no aditivo, no rango 1)

**Lema (trampa separable) [T].** Si `psi(alpha,beta) = f(alpha) g(beta)`, entonces `c` es de rango 1 en el interior: la diferencia mixta de un producto separable es `Delta f * Delta g`. Por eso un solo termino separable no basta: `psi` lleva dos monomios no proporcionales.

**Proposicion 54.2 [T].** (a) `c` no es aditiva: la diferencia mixta sobre cualquier rectangulo interior `(i,i') x (j,j')` es `s^2 (alpha_i - alpha_{i'})(alpha_j - alpha_{j'}) [4 + 9 (alpha_i + alpha_{i'} - s)(alpha_j + alpha_{j'} - s)] != 0`. En particular, para celdas adyacentes vale `s^4 (4 + 36 alpha_i alpha_j) > 0`. (b) `c` no es de rango 1: todo menor 2x2 interior es no nulo.

*Prueba de (a).* La formula interior de sec. 2 es suma de dos piezas de rango 1; la diferencia mixta de cada una es `(f_i - f_{i'})(g_j - g_{j'})`; sumando da la expresion indicada, cuyo corchete es `>= 4` (pues `alpha, beta >= 1/10 > s`). QED. *Prueba de (b).* `c` es suma de dos piezas de rango 1 con factores `(2 alpha - s)` y `(3 alpha^2 - 3 alpha s + s^2)`; el cociente `phi(alpha) = (2 alpha - s)/(3 alpha^2 - 3 alpha s + s^2)` es estrictamente decreciente en `[1/10, 1/2]` (su derivada es `-(6 alpha^2 - 6 alpha s + s^2)/(...)^2 < 0`, pues `6 alpha (alpha - s) + s^2 > 0`), asi que los dos factores no son proporcionales y el determinante `(f_i f'_{i'} - f'_i f_{i'})(g_j g'_{j'} - g_{j'} g'_j)` es no nulo. QED.

(a) y (b) se verifican ademas exactamente por el script en todas las celdas interiores adyacentes (sec. 5).

## 5. Verificacion reproducible

```text
python numerico/verifica_atomos_rectangulares.py
```

PASS en el PC (Fraction, `v = 10, 8`; `h = 12, 24, 48, 96`). El script comprueba: (a) todos los margenes cero; (b) `||c||_1 > 0`; (c) la formula interior exacta de sec. 2; (d) todas las diferencias mixtas interiores adyacentes no nulas; (e) todos los menores 2x2 interiores adyacentes no nulos.

| v | h | margenes | \|\|c\|\|_1 | entradas + / - | dif. mixtas | menores 2x2 |
|---|---|---|---|---|---|---|
| 10 | 12 | cero (46.4) | 0.1763 | 122 / 22 | 81/81 | 81/81 |
| 10 | 96 | cero (46.4) | 0.2905 | 9026 / 190 | 8649/8649 | 8649/8649 |
| 8 | 12 | cero (46.4) | 0.1830 | 122 / 22 | 81/81 | 81/81 |
| 8 | 96 | cero (46.4) | 0.2918 | 9026 / 190 | 8649/8649 | 8649/8649 |

(Los 8 casos pasan; `||c||_1` en decimal, valores exactos en la salida del script. `D1 = D2 = 0` exactamente en todos los casos.)

## 6. Alcance honesto y que sigue

- Lo que este doc SI dice: existe un peso `c(q,r)` explicito, genuinamente bidimensional, con margenes exactamente cero por construccion: el primero de la serie que satisface (46.4). La via de los pesos publicados esta cerrada (docs 47, 48, 52); la via de disenar el peso esta abierta.
- Lo que NO dice: no hay balance (etapa 2), no hay nivel de distribucion (este peso ya no es el de Chen: su nivel hay que probarlo desde cero, como advierte el doc 52 sec. 8.1), no hay cota inferior de capas ni demostracion de Goldbach.
- Etapa 2: balance. Calcular `sum c(q,r)` por el termino principal y compararlo con el objetivo; decidir si el peso tiene senal positiva a escala principal. Tambien: la diagonal `q = r` y el conteo de pares ordenados.

## Memoria compartida

Propuesta pendiente: el MCP canonico `lumina_shared_memory` / `lumina_propose` no esta expuesto. No se afirma persistencia remota.
