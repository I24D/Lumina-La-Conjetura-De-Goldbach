# 52 - Margenes del peso real de Chen contra el criterio del doc 46: el switching tampoco cancela

## Procedencia y estado

- **Agente:** LUMINA/Muse.
- **Fecha:** 2026-11-25, ~18:30 ET.
- **Metodo:** extraccion de la identidad publicada de Chen (fuentes leidas, no de memoria) y medicion exacta de margenes con aritmetica racional (`Fraction`).
- **Estado:** **[T-lit]** para la identidad extraida; **[T]** para las mediciones (exactas en dos modelos); **[P-rev]** para el criterio del doc 46 aplicado al peso real de Chen. No hay cota inferior nueva de capas ni demostracion de Goldbach.
- **Relacion con el doc 51:** este doc acepta la correccion adversarial del doc 51 y la precisa: el criterio no se aplica al marco general MZA, pero SI se aplica a la submatriz de switching P3/Omega de la identidad de Chen, donde ambas ramas cuentan los mismos objetos con el mismo factor local. Ahi tambien falla, a escala principal.

## 1. La identidad publicada de Chen (extraccion)

Fuente principal (leida completa): Z. Liu, "A Corrected Simplified Proof of Chen's Theorem", arXiv:2203.07871v1 (2022), que corrige a Pan et al. y sigue a Chen [2] = J.-R. Chen, *Sci. Sinica* 16 (1973) 157-176, Lemma 9. Corroborada por K. Matomaki y S. Zuniga-Alterman, "Weighted sieves with switching", *Math. Proc. Camb. Phil. Soc.* 179 (2025) 351-372, Example 3.4(ii)-(iii), que da la forma general del peso y atribuye los pesos de Kuhn con switching a Chen.

La identidad (ecuaciones como en [L]):

```text
W(N) = #{p : N-p sin factores primos <= N^{1/10}
              y a lo sumo un factor primo en (N^{1/10}, N^{1/3}]},
W(N) >= 2.6408 * S(N) * N / log^2 N.                          (2)
N-p contado por W(N) tiene exactamente tres factores primos ssi
N-p = p1 p2 p3, N^{1/10} < p1 <= N^{1/3} < p2 < p3.            (3)
Omega = sum_a sum_{a p3 <= N, N - a p3 primo} f(a),             (4)
f(a) = 1 ssi a = p1 p2, N^{1/10} < p1 <= N^{1/3}
       < p2 <= (N/p1)^{1/2}.
Omega <= 3.9404 * S(N) * N / log^2 N.                          (5)
#{p que cumplen (3)} <= Omega/2  (simetria), luego
#{1,2} >= W(N) - Omega/2 >= (2.6408 - 3.9404/2) * S(N) * N / log^2 N
       > 0.67 * S(N) * N / log^2 N.
```

Coeficientes y cortes extraidos:

- Cortes: rugosidad `z = N^{1/10}` (variante publicada: `z = N^{1/8}`, R. Baxter, charla Columbia 2024, ecs. (1)-(3)); ventana del peso `(N^{1/10}, N^{1/3}]`, es decir `(v,u) = (10,3)` y `(8,3)`.
- Coeficiente de switching: `1/2` (el conteo de P3 es a lo sumo `Omega/2`, por simetria: cada terna aparece dos veces en `Omega`).
- Peso de Kuhn/Chen [M, Ex. 3.4(ii)]: `w(alpha) = (1/2) 1_{[1/v,1/u]}`, `W_Kuhn(n) = 1 - (1/2) #{p|n : x^{1/v} <= p < x^{1/u}}`. "Kuhn type weights with switching have been used for detecting P2 numbers, starting from the work of Chen [2]."
- Peso de Richert [M, Ex. 3.4(iii)]: `w(alpha) = lambda (1 - u alpha) 1_{[1/v,1/u]}`, `W_Richert(n) = 1 - lambda sum_{p|n} (1 - u log p/log x)`; (A2) vale para `R >= 1/lambda + u - 1`. "Switching with Richert's weights was first utilised by Vaughan [14]."
- La cuna (q,r) de `Omega` es exactamente la region de integracion del Lemma 4 de [L]: `alpha = log p1/log N in [1/10,1/3]`, `beta = log p2/log N in [1/3,(1-alpha)/2]`, con `sum_a f(a)/(a log(N/a)) <= 0.49254/log N`. Su area es `49/3600` sobre el cuadrado `[1/10,1/2]^2` de area `0.16`.
- Teorema 3.2 de [M]: misma forma de peso de Kuhn con `(v,u) = (20,6)`, `S = 3`, nivel `0.267`. No se calcula aqui porque `S = 3` cambia la estructura de switching (P4); se calculan las instancias de Chen `(10,3)` y `(8,3)`.

## 2. Respuesta al doc 51: donde SI se aplica el criterio

El doc 51 tiene razon en el marco general de [M]: `Sigma_1` y `Sigma_2` son sucesiones distintas con factores locales distintos, la rama penalizada usa `W^+ >= 0`, la resta ocurre fuera de la rama, y no hay matriz firmada a la que aplicar el doc 46. Esa correccion se acepta.

Pero no cierra la pregunta para la identidad especifica de Chen, `W(N) - Omega/2`. Ahi hay una parte donde ambas ramas cuentan los mismos objetos con el mismo factor local del doc 44:

- La parte P3 de `W(N)`: fijado `(q,r) = (p1,p2)`, cuenta `l = p3` primo con `N - q r l` primo. Es exactamente el marco del doc 44: formas `L1(x) = x`, `L2(x) = N - q r x`, factor local `S(N;q,r)`.
- La parte `Omega`: fijado `(q,r)`, cuenta `l = p3` primo con `a l <= N` y `N - a l` primo. Mismo marco, mismo factor local `S(N;q,r)`.

La parte P2 de `W(N)` tiene otro factor local ("N - q r primo", una sola forma lineal) y queda fuera de este dominio. Por eso la medicion correcta es la **submatriz de switching** `cSW` (P3/Omega), no solo la matriz completa. El doc 46 se aplica exactamente a `cSW`.

## 3. El modelo (q,r) y sus identificaciones

Modelo ordenado del bucle (docs 47-48), truncacion P2/P3 (cofactor primo):

- Indices `i,j,k` en `{0..h-1}` para primos ordenados con tamanos logaritmicos `alpha_t` en `[1/v, 1/2]`, ventana `[1/v, 1/u)`, `b` = numero de indices en la ventana (`b/h = 7/12` para `(10,3)`).
- Par `(q,r) = (p_i, p_j)`, `i<j`; `l = p_k` es el cofactor primo.
- `N3p(i,j)` = #{k>j : alpha_i+alpha_j+alpha_k <= 1} (cofactores `l > r`).
- `N3o(i,j)` = #{k : alpha_i+alpha_j+alpha_k <= 1} (todos los cofactores).
- Cuna (ec. (4)): `i<b`, `j>=b`, `alpha_j <= (1-alpha_i)/2`.
- Matriz favorable de Kuhn: `c1[i][j] = (1 + N3p) * 1_{j>=b}` (P2: +1 si `j>=b`; P3: +1 por cada `k>j`).
- Matriz conmutada: `cOm[i][j] = N3o * 1_{cuna}`.
- Matriz diferencia de Chen: `c = c1 - (1/2) cOm`.
- Submatriz de switching: `cSW[i][j] = N3p * 1_{j>=b} - (1/2) cOm[i][j]`.
- Matriz favorable de Richert (`lambda = 1`, `u = 3`): `cR[i][j] = W_R(p_i p_j) + sum_{k>j} W_R(p_i p_j p_k)`, peso firmado sin exclusion dura.
- Dos ponderaciones: (A) modelo unitario (malla, `N = 10^12`, cada par/l cuenta 1); (B) pesos pi (conteo real de primos, `N2 = 10^7`, `N3p_pi = pi(N2/(p_i p_j)) - pi(p_j)`, `N3o_pi = pi(N2/(p_i p_j))`, exactos por criba a `10^6`).
- Medidas (docs 46-48): `D1 = sum_i |R_i| / ||c||_1`, `D2 = sum_j |C_j| / ||c||_1`; perturbacion de primer orden `P1 = sum_i a_i R_i + sum_j b_j C_j` con `a_i = 1/(p_i-2)` (primos representativos); cota de segundo orden `B2 = ||c||_1 / ((Q0-2)(R0-2))`; razon `|P1|/B2`. Criterio (46.4): `R_i = 0` para todo `i` Y `C_j = 0` para todo `j`.

Chequeos estructurales (exactos, verificados por el script):

- `N3p = 0` para `i>=b`: tres factores `> N^{1/3}` exceden `N`; es el "ssi" de la ec. (3). Luego `cSW` vive solo en `i<b<=j`.
- Fuera de la cuna, `N3p = 0` tambien: si `alpha_j > (1-alpha_i)/2` entonces `alpha_i + 2 alpha_j > 1`, luego ningun `k>j` cabe. Asi `cSW` esta soportada exactamente en la cuna triangular.
- Conteos de pares en la cuna: `(10,3)`: 12, 49, 196, 784 para `h = 12, 24, 48, 96`; `(8,3)`: 9, 49, 169, 729. Coherentes con el area `49/3600` del Lemma 4.

## 4. Resultados medidos

### 4.1 Modelo unitario (malla): matriz de Chen `c`

| (v,u) | h | D1 | D2 | \|P1\|/B2 |
|---|---|---|---|---|
| (10,3) | 12 | 63/125 = 0.504 | 77/125 = 0.616 | 4.60 |
| (10,3) | 24 | 505/807 = 0.626 | 565/807 = 0.700 | 5.34 |
| (10,3) | 48 | 2122/2787 = 0.761 | 751/929 = 0.808 | 4.58 |
| (10,3) | 96 | 17697/20482 = 0.864 | 480/539 = 0.891 | 3.49 |
| (8,3) | 12 | 32/59 = 0.542 | 37/59 = 0.627 | 7.45 |
| (8,3) | 24 | 218/389 = 0.560 | 281/389 = 0.722 | 8.88 |
| (8,3) | 48 | 1768/2451 = 0.721 | 1969/2451 = 0.803 | 8.30 |
| (8,3) | 96 | 15206/18241 = 0.834 | 16213/18241 = 0.889 | 8.21 |

El defecto crece con la resolucion `h` y se acerca a 1: escala principal, no pequeno.

### 4.2 Modelo unitario: submatriz de switching `cSW` (el dominio exacto del doc 46)

**D1 = D2 = 1.0000 exactamente en los 8 casos** (v = 10 y 8; h = 12, 24, 48, 96). `T/||cSW||_1 = -1.0000` exactamente: todas las entradas son `<= 0`. |P1|/B2: (10,3): 15.30, 9.34, 6.12, 4.05; (8,3): 33.31, 17.96, 12.22, 9.96. Filas no nulas (v=10): 6/12, 13/24, 27/48, 55/96; columnas no nulas: 3/12, 7/24, 14/48, 28/96.

Mecanismo (lema, exacto en la malla): para `(i,j)` en la cuna, todos los `k <= j` satisfacen la restriccion de tamano, porque `alpha_k <= alpha_j <= 1 - alpha_i - alpha_j` es exactamente la condicion de la cuna. Luego `N3o = (j+1) + N3p`, y

```text
cSW[i][j] = N3p - (1/2) N3o = (N3p - (j+1)) / 2 <= 0,
```

pues en la malla uniforme `N3p <= j+1`: la ventana en alpha por encima de `r` mide `1 - alpha_i - 2 alpha_j`, y la de por debajo mide `alpha_j - 1/v`; la primera es `<=` la segunda porque `alpha_j >= 1/3` y `alpha_i >= 1/v`. Es decir: **para cada par, la suma conmutada cuenta mas cofactores pequenos (j+1) que la rama favorable cuenta grandes (N3p)**. El termino Omega domina en cada par, la matriz es toda no positiva, y sus margenes son su masa total.

### 4.3 Pesos pi (conteo real de primos, N2 = 10^7)

Matriz de Chen `c`:

| (v,u) | h | D1 | D2 | \|P1\|/B2 |
|---|---|---|---|---|
| (10,3) | 12 | 9571/11513 = 0.831 | 11019/11513 = 0.957 | 1.04 |
| (10,3) | 24 | 4547/5398 = 0.842 | 2557/2699 = 0.947 | 0.28 |
| (10,3) | 48 | 107953/117474 = 0.919 | 112967/117474 = 0.962 | 0.31 |
| (10,3) | 96 | 493174/509683 = 0.968 | 492642/509683 = 0.967 | 0.23 |
| (8,3) | 12 | 5713/6150 = 0.929 | 1193/1230 = 0.970 | 1.80 |
| (8,3) | 24 | 109/125 = 0.872 | 2797/2900 = 0.965 | 1.03 |
| (8,3) | 48 | 117153/123751 = 0.947 | 120609/123751 = 0.975 | 1.08 |
| (8,3) | 96 | 131458/134613 = 0.977 | 131458/134613 = 0.977 | 0.75 |

Submatriz de switching `cSW` (pesos pi): D1 = 0.834, 0.844, 0.920, 0.968 (v=10); 0.931, 0.875, 0.948, 0.977 (v=8). D2 = 0.957, 0.946, 0.961, 0.967 (v=10); 0.969, 0.964, 0.974, 0.977 (v=8). Con pesos pi el signo ya no es uniforme (D < 1), porque `pi(x) ~ x/log x` hace que el conteo por encima de `r` pueda superar a `pi(r)` para pares pequenos; pero el defecto sigue siendo de escala principal (0.83-0.98) y crece con `h`. El veredicto no depende de la ponderacion.

### 4.4 Controles: rama favorable sola y Richert

- Rama favorable de Kuhn `c1`: **D1 = D2 = 1.0000 exactamente** en los 8 casos (todas las entradas `>= 0`, luego los margenes son la masa total). |P1|/B2 (malla): (10,3): 12.17, 8.31, 6.33, 4.61; (8,3): 20.70, 16.51, 12.37, 11.57. Es el fenomeno del doc 47 en el peso real: la rama favorable no tiene ninguna cancelacion interna.
- Rama favorable de Richert `cR` (lambda = 1, u = 3, malla): D1 = 0.915, 0.818, 0.776, 0.750 (v=10); 0.982, 0.940, 0.902, 0.876 (v=8). D2 = 1.000, 0.989, 0.983, 0.980 (v=10); 1.000, 0.995, 0.994, 0.992 (v=8). El peso continuo (que se anula en el borde de la ventana) no arregla los margenes: siguen a escala principal.

## 5. Veredicto

**El criterio (46.4) FALLA para el peso real de Chen, incluso en el dominio restringido donde se aplica exactamente.** Los margenes de la submatriz de switching son su masa total en el modelo unitario (D1 = D2 = 1 exactamente) y 0.83-0.98 con pesos pi reales; los de la matriz completa son 0.50-0.89 (malla) y 0.83-0.98 (pi), creciendo con la resolucion. La perturbacion de primer orden supera la cota de segundo orden por un factor 3.5-8.9 (matriz completa, malla, N = 10^12) y 4.1-33.3 (switching, malla). El switching (factor 1/2) reduce el defecto respecto a la rama favorable sola (de 1.00 a ~0.86 en la matriz completa), pero no lo elimina ni lo hace pequeno: la cancelacion es parcial, no estructural.

Sobre |P1|/B2 en el modelo pi (0.23-1.80, del orden de 1): es cancelacion firmada en la suma `P1` mas el decaimiento de `a_i = 1/(p_i-2)`; no es cancelacion estructural de los margenes (D ~ 0.97 ahi). Ademas, por la cota `|P1|/B2 <= N^{1/10} (D1+D2)`, el factor crece como `N^{1/10}` con `N`: a escala real el termino de primer orden domina ampliamente la cota de segundo orden. La mejora de segundo orden del doc 46 esta perdida para este peso.

## 6. Por que falla (mecanismo)

Tres razones, todas visibles en la identidad extraida:

1. **Asimetria del rango de `l`.** `Omega` cuenta todos los cofactores primos `l <= N/(qr)`; la rama favorable solo cuenta `l > r`. El conteo pequeno (`j+1` indices, o `pi(r)` primos) es comparable o mayor que el grande, asi que el neto por par es negativo o mixto sin cancelacion sistematica por fila.
2. **La cuna es triangular, no rectangular.** La condicion `beta <= (1-alpha)/2` hace que el soporte sea un triangulo; no hay simetria que produzca atomos rectangulares `(+1,-1,-1,+1)` del doc 46.
3. **Cada rama es positiva.** La rama favorable tiene D = 1 exactamente; el switching resta solo la mitad de una masa comparable. Ninguna de las dos ramas esta disenada para anular margenes: Chen las estima por separado y compara constantes escalares (2.6408 contra 3.9404/2).

## 7. Alcance honesto

- El total firmado `T` es negativo y `|T|/||c||_1 -> 0.78-0.98` en ambos modelos: el termino `Omega` domina en el conteo unitario/pi. **Esto NO contradice el teorema de Chen (0.67 > 0):** es un artefacto de contar cada (par, `l`) una vez, sin los pesos analiticos `1/(a log(N/a))` del Lemma 4, que son los que dan `2.6408 - 3.9404/2 = 0.6706 > 0`. Las medidas `D1, D2` son cocientes invariantes de escala y son la lectura estructural; el signo de `T` en estos modelos no es significativo.
- La truncacion P2/P3 (cofactor primo) es conservadora: los terminos P4+ omitidos son todos favorables (W(N) los cuenta), asi que solo pueden aumentar los margenes.
- El resultado `D = 1` exacto para `cSW` es un hecho de la malla uniforme; con pesos pi el defecto es 0.83-0.98, todavia principal. El veredicto no depende de la ponderacion.
- Lo que este doc NO dice: no refuta el teorema de Chen (que es correcto y usa las ramas por separado), ni dice que el switching no sirva para `{1,2}`. Dice que la ruta de cancelacion marginal del doc 46 no es la ruta por la que el peso de Chen funciona.

## 8. Que sigue

1. **Disenar un peso NUEVO con margenes cero por construccion** (creacion, no verificacion): coeficientes `c(q,r)` genuinamente bidimensionales, con atomos rectangulares, que anulen `R_q` y `C_r` por diseno en vez de heredarlos de un peso que solo depende de `n` por sus factores. Es la via que pide el mayor anhelo de Dal: crear, no verificar. El precio es que un peso asi ya no seria el de Chen y habria que probar su nivel de distribucion desde cero.
2. **Instanciar `b_{P,l}` para `A_N = {N-p}` y auditar A4-A5** (siguiente paso que propone el doc 51): determinar que parte de la comparabilidad de terminos principales se sigue de Bombieri-Vinogradov y donde aparece el modulo supra-raiz.
3. **Volver al frente analitico** (C10/peso movil W(d), uniones C_k+C_{k+1}, variante Lema 4 Pintz-Althofer), que es donde el bucle ya esta y donde la cancelacion marginal no es el cuello de botella.

## Verificacion reproducible

```text
python numerico/verifica_margenes_chen_richert.py
```

PASS exacto (`Fraction`) para `(v,u) = (10,3)` y `(8,3)`; `h = 12, 24, 48, 96`; modelo unitario (malla, N = 10^12) y pesos pi (N2 = 10^7). 48 chequeos del criterio (46.4), los 48 FALLAN. El script verifica ademas `N3p = 0` para `i >= b` (ec. (3)) y los conteos de la cuna contra el area del Lemma 4.
