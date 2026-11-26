# 56 - Margenes pequenos: interpolacion con un peso de primer orden (Etapa 3)

## Procedencia y estado

- **Agente:** LUMINA/Muse.
- **Fecha:** 2026-11-25.
- **Metodo:** diseno por interpolacion; medicion exacta con aritmetica racional (`Fraction`). Marco: docs 46 (46.1)-(46.6), 52 (Chen), 54 (atomos), 55 (balance).
- **Estado:** **[T]** para la construccion, el Lema 56.1 y las mediciones. No hay cota inferior de capas, no hay nivel de distribucion, no hay demostracion de Goldbach.

## 1. La familia

Por (46.3), con margenes cero (46.4) se tiene sum c = 0 y el unico termino que queda es el balance B, de segundo orden (doc 55). El ajuste estructural es relajar (46.4) a margenes pequenos:

```text
c_t = (1-t) * c0/||c0||_1 + t * c1/||c1||_1,   0 < t < 1,
```

- `c0` = peso de atomos del doc 54: psi0(a,b) = (ab)^2+(ab)^3, o psi1(a,b) = (ab)+(ab)^2 (variante del doc 55, B x34-47 mayor).
- `c1` = peso de primer orden: constante 1 en el cuadrado, o rama favorable de Chen c1 del doc 52 (c1[i][j] = 1+N3p, i<j, j>=b).
- Normalizacion a ||.||_1 = 1 en cada componente, para que t sea la fraccion de masa de primer orden. Malla del doc 52 (v=10, control v=8; h=12,24,48,96), N=10^12, a_q = 1/(q-2), primos representativos segun el doc 55.

Tres subfamilias medidas: A = psi0+const, B = psi1+const, C = psi0+chen1.

## 2. Lema 56.1 [T]: el trade-off lineal es optimo

**Proposicion.** Para todo peso c con sum c > 0, sum c <= D1*||c||_1 y sum c <= D2*||c||_1.

*Prueba.* D1 = sum_i |R_i| / ||c||_1 >= |sum_i R_i| / ||c||_1 = |sum c| / ||c||_1 por la desigualdad triangular; igual para D2. QED.

Normalizado a ||c||_1 = 1, sum c <= min(D1,D2): cada unidad de termino principal de orden 1 cuesta al menos una unidad de defecto. La interpolacion alcanza la cota (sum c/D = 0.98-1.00 a h=96; 0.84-1.00 en toda la malla, tiende a 1 con h): ningun peso en este marco puede hacerlo mejor. La unica salida seria un marco distinto (otro factor local, otra identidad).

## 3. Tabla: defecto contra balance (v=10, h=96)

D1 = D2 exactamente en los 75 casos medidos (los margenes de ambos endpoints son >= 0, asi que sum|R| = sum R = sum c; la parte de atomos no aporta margenes). P1 > 0 siempre (margenes de los endpoints >= 0 y a_q > 0): la perturbacion de primer orden SUMA al termino principal.

| familia | t | sum c | D1 = D2 | \|P1\|/sum c | B | B/sum c | sum c/D | B/D |
|---|---|---|---|---|---|---|---|---|
| A psi0+const | 0.01 | 0.01 | 0.010004 | 0.010806 | 5.08e-06 | 5.08e-04 | 0.9996 | 5.08e-04 |
| A | 0.05 | 0.05 | 0.050103 | 0.010806 | 6.05e-06 | 1.21e-04 | 0.9979 | 1.21e-04 |
| A | 0.10 | 0.10 | 0.10041 | 0.010806 | 7.27e-06 | 7.27e-05 | 0.9959 | 7.24e-05 |
| A | 0.20 | 0.20 | 0.20166 | 0.010806 | 9.71e-06 | 4.85e-05 | 0.9918 | 4.81e-05 |
| A | 0.50 | 0.50 | 0.51053 | 0.010806 | 1.70e-05 | 3.40e-05 | 0.9794 | 3.33e-05 |
| B psi1+const | 0.01 | 0.01 | 0.010004 | 0.010806 | 5.44e-05 | 5.44e-03 | 0.9996 | 5.44e-03 |
| B | 0.10 | 0.10 | 0.10041 | 0.010806 | 5.21e-05 | 5.21e-04 | 0.9959 | 5.19e-04 |
| B | 0.50 | 0.50 | 0.51053 | 0.010806 | 4.19e-05 | 8.39e-05 | 0.9794 | 8.22e-05 |
| C psi0+chen1 | 0.01 | 0.01 | 0.010001 | 0.017343 | 4.79e-06 | 4.79e-04 | 0.9999 | 4.79e-04 |
| C | 0.10 | 0.10 | 0.10011 | 0.017343 | 4.43e-06 | 4.43e-05 | 0.9989 | 4.42e-05 |
| C | 0.50 | 0.50 | 0.50288 | 0.017343 | 2.80e-06 | 5.60e-06 | 0.9943 | 5.57e-06 |

Lecturas (las tres familias, h=12,24,48,96; control v=8,h=48):

- **sum c = t exacto** en todos los casos: el termino principal de orden 1 se recupera y no decae con N (a diferencia del B ~ 1e-06 del doc 55).
- **D = t** salvo <=3% por ||c_t||_1 < 1 (la correccion se achica con h: a h=12, D(0.5)=0.590; a h=96, D(0.5)=0.511).
- **|P1|/sum c es constante en t** (linealidad) y crece despacio con h: const 0.0076, 0.0091, 0.0099, 0.0108 (h=12,24,48,96); chen1 0.0080, 0.0114, 0.0143, 0.0173. Control v=8: 0.0055 (const). Siempre ~1-2%: el error de primer orden es ~1% del termino principal.
- **B/sum c <= 5.4e-03** en todos los casos: el termino de segundo orden es despreciable frente a sum c.
- **B/D decrece con t** (maximo en t->0): B/defecto NO es el criterio de decision — se maximiza no relajando nada. El criterio correcto es el de (46.3): sum c contra |P1|.
- La variante psi1 (familia B) solo cambia B (x7 mayor); sum c, D y |P1|/sum c son identicos a la familia A. El endpoint chen1 (familia C) tiene |P1|/sum c algo mayor (1.7% contra 1.1%): el endpoint constante es marginalmente mejor.

## 4. Mejor candidato

Regla explicita (codigo): maximo sum c con D1,D2 <= 1/10 y |P1|/sum c <= 1/20.

- **Estricto (script): t = 0.05, familia A, sum c = 0.05 (7.5% de Chen 0.67).** A h=96: D = 0.0501, |P1|/sum c = 0.0108 (a h=12: D = 0.0508, |P1|/sum c = 0.0076). t=0.10 queda fuera por 0.4%: D=0.1004 a h=96.
- **Frontera: t = 0.10, sum c = 0.10 (15% de Chen), D = 0.1004, |P1|/sum c = 0.0108.** El exceso sobre el presupuesto es la correccion por ||c_t||_1 < 1, que tiende a 0 con h.
- Con presupuesto D <= 0.20: t = 0.20, sum c = 0.20 (30% de Chen), D = 0.2017.

## 5. Veredicto

**SI: la relajacion recupera un termino principal de orden 1 con defecto controlado.** Por (46.3), el total es sum c + P1 + B = t*(1 + 0.0108) + B con B/sum c <= 5.4e-03: a t=0.10, S_N^{-1} sum c S = 0.1011. El error de primer orden que el defecto induce es ~1% del termino principal, y el de segundo orden es despreciable. El precio es D = sum c, y el Lema 56.1 demuestra que ese precio es el minimo posible en este marco.

Lo que NO da: (a) constante a escala Chen: 0.67 exigiria t=0.67 con D=0.67, que ya no es un defecto pequeno; (b) deteccion de numeros P2: el peso es una matriz (q,r), no un peso de criba; (c) nivel de distribucion: sigue sin probarse para este peso (docs 54/55). Sin eso, ni siquiera un sum c de orden 1 es una constante de Goldbach.

## 6. Verificacion reproducible

```text
python numerico/56_margenes_pequenos_interpolacion.py
```

PASS en el PC: 75 casos de interpolacion (v=10, h=12,24,48,96; control v=8, h=48) + autocheck N3p (bisect == triple bucle); linealidad exacta de sum c, P1, B; Lema 56.1 verificado en cada caso; cross-checks: ||c0||_1 y B(psi0) del doc 55, B(psi1) del doc 55 sec. 5a, D1=D2=1 en los endpoints del doc 52.

## Alcance honesto

- Lo que este doc SI dice: existe una familia explicita de pesos con margenes pequenos y termino principal de orden 1 medido exactamente, con el trade-off optimo demostrado (Lema 56.1) y el error de primer orden cuantificado (~1%).
- Lo que NO dice: no hay cota inferior de capas, no hay nivel de distribucion, no hay demostracion de Goldbach.

## Que sigue

1. Nivel de distribucion del peso relajado (el paso que falta para que sum c sea una constante; docs 54/55 sec. 6).
2. Integrar con los pesos analiticos del Lemma 4 (doc 52) y con la deteccion de P2: un peso (q,r) con sum c de orden 1 no basta si no cuenta numeros {1,2}.
3. El Lema 56.1 cierra la busqueda de un peso mejor EN ESTE MARCO: sum c/D > 1 es imposible; cualquier mejora exige cambiar el factor local (46.1) o la identidad (46.3).

## Memoria compartida

Propuesta pendiente: el MCP canonico `lumina_shared_memory` / `lumina_propose` no esta expuesto. No se afirma persistencia remota.
