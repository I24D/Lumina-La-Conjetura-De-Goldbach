# 66. Minorante ρ de Harman (SL-65.2, etapa 1)

Fecha: 2026-11-26. Etapa 1 del proyecto SL-65.2 (doc 65 §8): sustituir Λ por el
minorante ρ de Harman. Solo etapa 1: formulación de ρ, por qué evita el bloqueo
de paridad del Teor. 1.3, y mediciones del término principal y el defecto.
Script: `numerico/66_minorante_harman_rho.py` — PASS 20/20, verificado local y
en el PC con números idénticos. No se afirma ninguna prueba de Goldbach.

## 1. Formulación de ρ

Fuente: J. Maynard, "Primes in arithmetic progressions to large moduli III:
Uniform residue classes", arXiv:2006.08250v1 (2020); Mem. Amer. Math. Soc.
306(1544), 2025. Sección 9 (prueba del Teor. 1.3), siguiendo la criba de
Harman ([21] en el paper: G. Harman, "Prime-Detecting Sieves", London Math.
Soc. Monographs 33, Princeton Univ. Press, 2007). Todo lo citado aquí fue
leído en la copia ORA del paper.

Notación: z₁ = x^{1/7}, z₂ = x^{3/7}, z₃ = x^{4/7}; z₀ = x^{1/14−4δ};
δ < 1/1000. ρ(n,z) := 1 si P⁻(n) ≥ z, 0 en otro caso. Para n ∼ x,
ρ(n,√(2x)) es el indicador de los primos.

Conjuntos bueno/malo (§9):
G = [x^{3/7+8δ}, x^{1/2−4δ}] ∪ [x^{1/2+4δ}, x^{4/7−8δ}],
B = [x^{3/7}, x^{3/7+8δ}] ∪ [x^{1/2−4δ}, x^{1/2+4δ}] ∪ [x^{4/7−8δ}, x^{4/7}].

Descomposición de Buchstab (§9): ρ(n,√(2x)) = S₁ − S₂ − S₃; S₂ se descompone
en S_{2,1} − S_{2,2} − S_{2,3} − S_{2,4} − S_{2,5}; S_{2,2} se vuelve a
descomponer (S_{2,2,1} − S_{2,2,2}); S₃, S_{2,4} y S_{2,2,2} se parten en
partes G y B. Resultado:

ρ(n,√(2x)) = S^triv(n) + S^G(n) + S^B(n) + S_{2,3}(n) + S_{2,5}(n).

Definición del minorante (§9): ρ(n) := S^triv(n) + S^G(n) − C·ρ_B(n), con C
constante absoluta grande y ρ_B la función del Lema 6.8.

Teorema 1.3 (leído): (1) ρ(n) ≤ 1 si n es primo, 0 en otro caso;
(2) Σ_{n≤x} ρ(n) ≥ π(x)/8; (3) ρ equidistribuido en progresiones
aritméticas para Q₁ ∈ [x^{2/5+5δ}, x^{3/7}], Q₂ = x^{1/2+δ}/Q₁, con error
≪ x/(φ(q₁q₂)(log x)^A), uniforme en a.

Los dos términos descartados (regiones según el cálculo de densidad de §9,
ver §4):
S_{2,3}(n) = Σ ρ(n/(p₁p₂), p₂) sobre z₁ ≤ p₂ < p₁ < z₂, p₁p₂ ≤ z₂,
  p₁p₂² ≥ z₃;
S_{2,5}(n) = Σ ρ(n/(p₁p₂), p₂) sobre z₁ ≤ p₂ < p₁ < z₂, p₁p₂ ≥ z₃.

## 2. Por qué ρ evita el bloqueo de paridad

S_{2,3}(n), S_{2,5}(n) ≥ 0: son sumas de ρ(n/(p₁p₂), p₂) ≥ 0, es decir,
términos de CONTEO (números con dos factores primos pequeños en rangos
dados). Como son ≥ 0, se pueden DESCARTAR en la definición del minorante:

ρ(n) = ρ(n,√(2x)) − S^B(n) − Cρ_B(n) − S_{2,3}(n) − S_{2,5}(n) ≤ ρ(n,√(2x)),

usando S^B(n) ≥ −Cρ_B(n) (§9). El precio es la densidad (se pierde la masa
de S_{2,3}+S_{2,5}, medida en §4), no la equidistribución: S^triv, S^G y ρ_B
equidistribuyen (Prop. 5.1 y Lema 6.8), luego ρ también (Teor. 1.3(3)).

Para Λ, la identidad de Heath-Brown (Lema 6.1, leído) es EXACTA y FIRMADA:
Λ(n) = Σ_j (−1)^j Σ μ(m₁)···μ(m_k) log n₁. Los términos análogos a
S_{2,3}, S_{2,5} llevan coeficientes μ firmados y NO se pueden descartar:
descartarlos cambiaría la suma, porque no hay desigualdad de un solo lado.
Esa es la obstrucción de paridad del Teor. 1.3 para Λ (doc 65): la prueba de
Maynard descarta S_{2,3}, S_{2,5} ≥ 0, lo que solo es legítimo para un
minorante construido con términos de conteo.

En resumen: ρ evita el bloqueo porque renuncia a la exactitud (es un
minorante, no una identidad) y porque los términos que descarta son conteos
no negativos. El bloqueo de paridad es un bloqueo para identidades exactas
con coeficientes firmados, no para minorantes.

## 3. Estructura exacta de ρ (Props. 66.1–66.2, derivadas de §9)

Prop. 66.1: ρ(n) = 1 para n primo, n ∼ x. Prueba: para n primo,
S_{2,3}(n) = S_{2,5}(n) = 0 (sin dos factores primos), S^B(n) = 0 y
ρ_B(n) = 0 (sin factor en B: los únicos divisores son 1 y n ∼ x, fuera de
B), S^triv(n) = 1 y S^G(n) = 0; luego ρ(n) = 1. Nótese: el Teor. 1.3(1)
solo afirma ρ(n) ≤ 1 en primos; la igualdad es consecuencia de §9, no del
enunciado.

Prop. 66.2: ρ(n) ≤ 0 para n compuesto, n ∼ x. Prueba: ρ(n,√(2x)) = 0 para
n compuesto, y ρ = ρ(n,√(2x)) − S^B − Cρ_B − S_{2,3} − S_{2,5} con
S^B ≥ −Cρ_B y S_{2,3}, S_{2,5} ≥ 0. Es decir, en el rango diádico
ρ = 1_p − D con D(n) = S^B(n) + Cρ_B(n) + S_{2,3}(n) + S_{2,5}(n) ≥ 0
soportado en compuestos. En particular ρ(n) < 0 cuando
S_{2,3}(n) + S_{2,5}(n) ≥ 1.

Corolario (corrección del doc 65 §8): la desigualdad Σ_n ρ(n)ρ(N−n) ≤ R(N)
es FALSA como está escrita, porque ρ es negativa en (muchos) compuestos y
el producto de dos valores negativos contribuye positivamente. El reemplazo
válido de un solo lado es Σ_{p<N} ρ(N−p) ≤ R(N) (término a término, pues
ρ(m) ≤ 1_p(m) para todo m), con constante del término principal c_ρ
(no c_ρ²). Ver §4 para el valor de c_ρ.

## 4. Medición: término principal y defecto

Pérdida L_S = (ΣS_{2,3} + ΣS_{2,5})/π(x). El paper (§9) la calcula con las
integrales de Buchstab (u = log p₁/log x, v = log p₂/log x, integrando
ω((1−u−v)/v)/(u·v²), ω la función de Buchstab del Lema 6.5):

I_a (S_{2,3}): u ∈ [4/21, 2/7], v ∈ [2/7−u/2, min(u, 3/7−u)] → 0.127;
I_b (S_{2,5}): u ∈ [2/7, 3/7], v ∈ [4/7−u, min(u, (1−u)/2)] → 0.409;
L_S = I_a + I_b = 0.535 (asintótico).

Verificaciones (script, PASS 20/20):
- ω verificada numéricamente (ω(2) = 1/2, ω(3) = (1+ln 2)/3,
  ω(10) → e^{−γ}) y contra conteos exactos de la criba (Lema 6.5:
  error rel. 0.076 a z=7, 0.009 a z=101, 0.029 a z=502, X=2×10⁶).
- I_a, I_b convergen <1% entre mallas 1600² y 3200²; doble verificación con
  scipy.quad: I_a = 0.127, I_b = 0.408, suma 0.535.
- Cota cruda del paper (ω ≤ 1): I_a + I_b = 0.839, 1 − 0.839 = 0.161 ≥ 1/8,
  consistente con Teor. 1.3(2). (El paper enuncia la conclusión como ≥ 1/8;
  nuestro cálculo da 0.161.)
- Criba exacta a X = 2×10⁶: L_S = 0.325 (S_a/π = 0.066, S_b/π = 0.259;
  π(X) = 148933). La diferencia con 0.535 es efecto de escala finita
  (O(1/log X) con constante grande: la densidad de primos en escala
  logarítmica es pobre a X = 2×10⁶).
- Convergencia Buchstab por pares: L_S = 0.418 (2×10⁶) → 0.462 (10⁸) →
  0.499 (10¹⁰) → 0.513 (10¹²) → 0.535 (integral). A 10¹² está a <5% del
  integral.
- Fracción con ρ < 0: ≥ 0.0189 a X = 2×10⁶ (n con S_{2,3}+S_{2,5} ≥ 1).
- Masa B: fracción de n ≤ X con factor primo en B = 0.0431 (O(δ),
  Lema 6.8(3)).

Densidad y defecto. De la identidad exacta
ρ(n,√(2x)) − ρ(n) = S^B + Cρ_B + S_{2,3} + S_{2,5}, sumando y usando
|S^B| ≤ Cρ_B (Lema 6.8(3)):

c_ρ := lim Σρ/π = 1 − L_S − O(δ) = 0.465 − O(δ).

La cota cruda implica O(δ) ≤ 0.036 (pues 1 − 0.839 − O(δ) ≥ 1/8), luego
c_ρ ∈ [0.43, 0.47] asintóticamente; cota rigurosa del teorema: c_ρ ≥ 1/8.
Defecto: 1 − c_ρ ≈ 0.535 asintótico; ≤ 7/8 riguroso.

Constante unilateral: el conteo Σ_{p<N} ρ(N−p) ≤ R(N) tiene constante del
término principal c_ρ ∈ [0.43, 0.47] (rigurosa ≥ 1/8). Esto corrige el
doc 65 §8: la constante NO es ~1/64 = 0.0156, sino ≥ 1/8 (8× mayor) y
plausiblemente ~0.465 (30× mayor).

Nota de verificación (fracaso como experimento): la primera lectura de
S_{2,3} como {z₂ ≤ p₁p₂ < z₃} daba L_S = 1.17 por integrales, incompatible
con Teor. 1.3(2) (daría c_ρ < 0). Las integrales de densidad del propio §9
fijan las regiones correctas (I_a, I_b arriba; la restricción u ≥ 4/21 que
el paper explica en §9 confirma la región de S_{2,3}: el corte es en
p₁p₂² ≷ z₃, no en p₁p₂ ≷ z₂). La criba exacta y la convergencia por pares
confirman los valores corregidos.

## 5. Veredicto y qué falta para la etapa 2

Veredicto: ρ SÍ evita el bloqueo de paridad del Teor. 1.3 (los términos
descartados son conteos ≥ 0; para Λ son firmados y no se pueden descartar).
El precio es el defecto: ρ captura solo c_ρ ≈ 0.46 de los primos
asintóticamente (≥ 1/8 riguroso). El conteo unilateral Σ_p ρ(N−p) ≤ R(N)
tiene constante c_ρ ∈ [0.43, 0.47], no c_ρ².

Etapa 2 (no hecha aquí): rehacer el análisis de criba / método del círculo
con ρ en lugar de Λ para el conteo unilateral Σ_{p<N} ρ(N−p). Término
principal esperado: c_ρ·𝔖(N)·N/log²N (heurístico: la serie singular es la
misma, porque ρ se anula en múltiplos de primos pequeños y es uniforme en
clases reducidas, igual que 1_p; a probar en la etapa 2). Control del error
con el nivel 1/2+δ del Teor. 1.3(3) en la ventana Q₁ ∈ [x^{2/5+5δ}, x^{3/7}].
Si el error queda bajo el término principal, se obtiene una cota inferior
de Goldbach con constante c_ρ. La esquina DUA (doc 65) no es el cuello de
botella para ρ: ρ equidistribuye en toda la ventana, no solo en la esquina.

Referencias (leídas, no de memoria):
- J. Maynard, "Primes in arithmetic progressions to large moduli III:
  Uniform residue classes", arXiv:2006.08250v1 (2020), §6 (Lemas 6.1, 6.2,
  6.5, 6.8) y §9 (prueba del Teor. 1.3). Copia ORA.
- G. Harman, "Prime-Detecting Sieves", London Math. Soc. Monographs 33,
  Princeton Univ. Press, 2007 — ref. [21] de Maynard III (la construcción
  de §9 sigue su criba).
