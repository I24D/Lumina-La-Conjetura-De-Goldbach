# 07 — Sesgo de la correlación Möbius–Goldbach, fórmula de arcos mayores, (HB) para casi todo N y formalización en Lean

## Procedencia y estado

- **Agente:** Claude Code. **Modelo:** `claude-opus-5-5`.
- **Fecha:** 2026-11-23 (America/New_York).
- **Método:** experimentos numéricos nuevos (FFT para todo N ≤ X), derivación escrita y formalización en Lean 4 + Mathlib.
- **Estado:**
  - §1: **[P]** con entradas clásicas.
  - §§2–3: **[N] + [H]** (fenómeno numérico y explicación heurística).
  - §4: **[T-Lean]**.
- **Dependencias:** 04, 05, 06.
- **Goldbach no queda demostrada.** Este documento aporta un teorema "para casi todo N" en el marco de cajas, un fenómeno numérico con su fórmula explicativa y las identidades base verificadas por máquina.

---

## 1. (HB) se cumple para casi todo N, en todas las cajas [P]

**Teorema 7.1.** Para cada A, B > 0 y X suficientemente grande, todos los N pares de (X/2, X] **salvo O_{A,B}(X/(log X)^B)** cumplen

    |B_N(C,K)| ≤ X/(log X)^A  en todas las cajas,   y por tanto   T_N^prim(D) ≪ X/(log X)^{A−2}.

Con la Prop. 4.6 de 04: R(N) = N𝔖(N) + O(N/(log N)^{A−2}) para esos N.

**Prueba.**
1. **Quitar la dependencia en N de los pesos.** Por `coprime_of_coprime_sub` (formalizado en §4), la condición (c, N) = 1 solo cambia términos con (N−cd, N) > 1. Por `primePow_of_not_coprime`, esos son potencias de primos que dividen a N, y su aporte total es ≪ N^ε (04, Prop. 4.2). Sea B̃_N(C,K) la caja **sin** la condición de coprimalidad.
2. **Cajas pequeñas en d.** Las cajas con K ≤ √X/(log X)^{B₁} se controlan para **todo** N por el Lema 6.1 de 06. Las demás cumplen K ≥ X^{1/3} y log K ≍ log X.
3. **Parseval sobre N.** Sea a(m) = Σ_{cd=m, c∈(C,2C], d∈(K,2K]} μ(d) log d, restringida a m < X. Sean A(α) = Σ a(m) e(mα) y S(α) = Σ_{n≤X} Λ(n) e(nα). Para N ≤ X, B̃_N es el coeficiente N-ésimo de A·S. Por tanto

        Σ_{N≤X} |B̃_N|² = ∫₀¹ |A|²|S|² ≤ sup_α|A(α)|² · Σ_{n≤X} Λ(n)² ≪ sup|A|² · X log X.

4. **Davenport (1937) [T].** Uniformemente en θ, Σ_{d≤y} μ(d) e(dθ) ≪_J y (log y)^{−J}. Por sumación parcial (peso log d) y aplicándolo con θ = cα para cada c, se obtiene sup_α|A(α)| ≪ C·K/(log X)^{J−1} ≤ X/(log X)^{J−1}.
5. **Chebyshev.** #{N : |B̃_N| > X/(log X)^A} ≪ X/(log X)^{2J−3−2A}. Se elige J grande y se toma la unión sobre las O(log² X) cajas. □

**Consecuencias y lectura:**
- Se recupera el teorema clásico E(X) ≪ X/(log X)^B (Chudakov, van der Corput y Estermann, 1937–38) **dentro del marco de cajas de 05**. No es un resultado nuevo de la literatura; lo nuevo es la ubicación.
- **Lección precisa:** (HB) es *fácil en promedio sobre N*, porque solo usa la cota uniforme de Davenport. Toda la dificultad está en **cada N por separado**. Cualquier demostración de (HB) tiene que usar información de un N concreto, que Parseval sobre N no ve. Es lo mismo que en 01 §2.1, ahora en el marco de 04–06.
- **[C: GRH] Variante.** Baker–Harman (1991, J. London Math. Soc., *verificar*) dan, bajo GRH, Σ_{d≤y} μ(d)e(dθ) ≪ y^{3/4+ε} uniformemente. Con K ≥ X^{1/2−ε}, eso da sup|A| ≪ X^{7/8+ε} y un conjunto excepcional ≪ X^{3/4+ε}. Es **peor** que el X^{1/2+ε} de Hardy–Littlewood, porque el marco de cajas ignora la estructura de arcos mayores de S. Es una comprobación de coherencia, no una mejora.

---

## 2. Fenómeno numérico: sesgo y varianza sub-aleatoria de la correlación Möbius–Goldbach [N]

Se define, para todo N ≤ X (por FFT, con error relativo 3,7·10⁻¹⁵ frente a sumas directas),

    M0(N) = Σ_{n<N} μ(n) Λ(N−n),         V0(N) = Σ_{n<N} μ(n)² Λ(N−n)²,
    M1(N) = Σ_{n<N} μ(n) log n Λ(N−n),   V1(N) = Σ μ(n)² log² n Λ(N−n)².

M1, cortada en n > D, es la fila c = 1 de la cola T de 04–06, es decir, el núcleo de paridad. En el **modelo aleatorio** (μ como signos independientes sobre los libres de cuadrados), Z = M/√V tendría media 0 y varianza 1, sin dependencia de N.

**Datos.** Ventana N par en (X/2, X]. Script `numerico/mobius_goldbach_aleatoriedad.py`, resultado `numerico/resultados/20261123_220731_mobius_goldbach_aleatoriedad_1000000.json`.

| X = 10⁶ | Var(Z), N par | Var(Z), N impar | media Z si 3 \| N | media Z si 3 ∤ N | corr(Z, 𝔖(N)) |
|---|---|---|---|---|---|
| M1 (μ log) | 0,875 | 0,716 | −0,167 | +0,040 / +0,051 | −0,427 |
| M0 (μ) | 0,896 | 0,746 | −0,219 | +0,055 / +0,066 | −0,489 |

Dos desviaciones claras del modelo aleatorio:
1. **Varianza sub-aleatoria**, alrededor de 0,87–0,90 para N par.
2. **Sesgo negativo cuando N tiene factores primos pequeños.** También cuando 5 | N (−0,12) y 7 | N (−0,10). En consecuencia, Z se correlaciona negativamente con la serie singular 𝔖(N).

Con X = 10⁷ el sesgo es menor (media de Z si 3 | N: −0,094; corr(Z, 𝔖) = −0,264): **decrece con N**.

---

## 3. Fórmula de arcos mayores que explica el sesgo [H + N]

**Heurística.** El desarrollo de Ramanujan–Hardy de la densidad de primos, Λ(m) ≈ Σ_{q≤Q} μ(q) c_q(m)/φ(q), donde c_q es la suma de Ramanujan, equivale a la contribución de los arcos mayores de nivel Q en el método del círculo. Sustituyéndolo en M0:

    M0(N) ≈ P_Q(N) := Σ_{q≤Q} (μ(q)/φ(q)) · Σ_{1≤n<N} μ(n) c_q(N−n).

- **El término q = 1 es la función de Mertens M(N−1).**
- El término q = 2 da 2Σ_{n<N, n impar} μ(n) para N par.
- En general, P_Q es una combinación lineal explícita de **sumas de Möbius en progresiones aritméticas de módulo q ≤ Q**, gobernada por los ceros de las funciones L de Dirichlet.

**Prueba numérica.** Script `numerico/arcos_mayores_mobius.py`. Resultados: `numerico/resultados/20261123_220959_…_1000000.json`, `…221209_…_1000000.json` y `…221311_…_10000000.json`.

X = 10⁶ (250 001 valores de N par):

| Q | corr(M0, P_Q) | R² | Var(M0−P_Q)/⟨V0⟩ | sesgo 3\|N | sesgo 5\|N | sesgo 7\|N | corr(residuo, 𝔖) |
|---:|---:|---:|---:|---:|---:|---:|---:|
| (sin) | — | 0 | 0,910 | −0,219 | −0,122 | −0,102 | −0,489 |
| 1 | 0,067 | 0,004 | 0,906 | −0,208 | −0,110 | −0,091 | −0,491 |
| 2 | 0,076 | 0,006 | 0,905 | −0,182 | −0,083 | −0,064 | −0,492 |
| 6 | 0,187 | 0,035 | 0,878 | **+0,004** | −0,084 | −0,069 | −0,393 |
| 30 | 0,355 | 0,126 | 0,796 | +0,001 | **+0,001** | **+0,011** | **−0,045** |
| 60 | 0,422 | 0,178 | 0,748 | 0,000 | 0,001 | 0,000 | −0,038 |
| 210 | 0,541 | 0,292 | 0,644 | 0,000 | 0,001 | −0,001 | +0,002 |
| 1000 | 0,668 | 0,446 | 0,505 | 0,000 | 0,001 | −0,001 | +0,011 |

X = 10⁷: sin corregir, los sesgos son −0,094, −0,046 y −0,041, con corr(Z, 𝔖) = −0,264. **Con Q = 30 pasan a +0,0006, +0,0003 y +0,005, y la correlación a −0,017.**

**Conclusiones [N].**
1. **El sesgo de clase queda explicado** por la fórmula explícita P_Q:
   - el módulo q = 3 aparece en P_Q a partir de Q = 6, y su sesgo desaparece exactamente ahí;
   - con Q = 30 desaparecen los sesgos de 5 y 7.

   Esto encaja con una causa de arcos mayores y no con un efecto aleatorio.
2. El sesgo es un **término secundario de tamaño ≍ √N**, con origen en las sumas de μ en progresiones. Bajo GRH (Alpöge), P_Q(N) ≪_Q N^{1/2+ε} *(derivación estándar; verificar)*, así que es inofensivo para Goldbach.
3. **El residuo M0 − P_Q no tiene sesgo de clase**, y su varianza sigue siendo sub-aleatoria (0,80 con Q = 30).

   **Explicación de la varianza sub-aleatoria [H].** Por Parseval, Σ_N M0(N)² = ∫|A_μ|²|S|², mientras que Σ_N V0(N) ≈ ∫|A_μ|²·∫|S|². El cociente es menor que 1 si |A_μ|² es pequeña justo donde |S|² es grande, es decir, en los arcos mayores. Es lo que ocurre: μ no se correlaciona con caracteres de módulo pequeño (Davenport) y es "más pequeña que aleatoria" en 0. La media de M(x)²/x es ≈ 0,05, frente al 6/π² ≈ 0,61 del modelo aleatorio (Ng 2004, bajo RH y hipótesis adicionales, *verificar*).
4. **Para Goldbach [P]:** el núcleo de paridad (fila c = 1 de la cola) se descompone en
   - una **parte de arcos mayores, explícita y calculable**, que bajo GRH es pequeña;
   - una **parte de arcos menores**, sin sesgo, que es la que contiene la barrera de paridad.

   La fórmula P_Q **no** controla esta segunda parte. **Esto no es una demostración.**

**Novedad.** No conozco una publicación que documente este sesgo de clase para la correlación binaria μ×Λ y su eliminación exacta por P_Q. Es el tipo de fenómeno (término secundario de tipo Chebyshev) que la heurística de arcos mayores predice, así que **no se reivindica como nuevo** sin una búsqueda bibliográfica. Queda registrado como observación verificable y reproducible.

---

## 4. Formalización en Lean 4 + Mathlib [T-Lean]

Archivo canónico: `C:\Users\dal_n\lean\goldbach_lumina\GoldbachLumina\Reduccion.lean`, con copia en `lean/Reduccion.lean`. Compila con `lake env lean` (exit 0). `#print axioms` da para todos **solo [propext, Classical.choice, Quot.sound]**: sin `sorry` y sin axiomas propios.

| Teorema | Enunciado | Documento |
|---|---|---|
| `R_eq_neg_sum` | R(N) = −Σ_{0<n<N} Σ_{d∣n} μ(d) log d · Λ(N−n) | 04 §2 |
| `R_eq_S_sub_T` | **R(N) = S_N(D) − T_N(D) para todo N y todo D** (forma de divisores) | 04, Prop. 4.1 |
| `coprime_of_coprime_sub` | gcd(N−cd, N) = 1 ⇒ gcd(c, N) = 1 | 06, Lema 6.1 paso 1; §1 de este documento |
| `primePow_of_not_coprime` | Λ(m) ≠ 0 y gcd(m,N) ≠ 1 ⇒ m = p^k con p ∣ N | 04, Prop. 4.2 (núcleo) |
| `mu_log_prime` | μ(p) log p = −log p | 06, Lema 6.0 |
| `divisors_prime_gt` | p primo y 1 ≤ D < p ⇒ {d ∣ p : d > D} = {p} | 06, R-3 |

Con estos teoremas, las identidades exactas de las que parten 04, 05 y 06 ya no dependen de la revisión entre agentes: las comprobó el kernel.

---

## 5. Qué queda y siguiente paso

- **Sigue abierto:** la parte de arcos menores del núcleo de paridad (06, F1 y F2). Ninguno de los resultados de este documento la toca.
- **Siguiente paso numérico [P]:** comprobar si el residuo de arcos menores de la fila c = 1 tiene una distribución **simétrica**, o si conserva algún sesgo de signo que pueda servir para la cota unilateral (HB+). Hasta Q = 1000, la media por clases del residuo es ≈ 0, así que no se ve sesgo aprovechable.
- **Siguiente paso formal [P]:** formalizar en Lean la reindexación (c, d) de 04 y la desigualdad (4.12) (potencias propias), que son finitas y elementales.
- **Siguiente paso teórico [P]:** escribir con detalle la cota P_Q(N) ≪ N^{1/2+ε} bajo GRH, con la dependencia en Q explícita.
