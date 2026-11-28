# Contraste numérico del doc 107 (L105(a): B_Δ ≤ C_Δ·T(N))

**Agente:** Claude Code (Anthropic, modelo `claude-opus-5-5`), brazo verificador de la línea Goldbach.
**Fecha:** 2026-11-28.
**Encargo:** de Muse, por orden de Dal (encargo 3 de 4). Blindar el teorema insignia recomputando con código propio sus calibraciones:
- los chequeos K0–K6;
- B_Δ ≤ 2,16·10⁻⁴·T(N);
- la descomposición por niveles (S2–S5 ≥ 97 %; d > D_A^{2/3} entre 80 y 97 %);
- el contraejemplo de la errata al doc 97 (r = 458815, factor 1846).

**Método:** reimplementación propia desde las definiciones de los docs 97, 98, 100, 101, 106 y 107, sin importar código de nadie.
- El contraejemplo se recalcula por **dos rutas independientes**: la forma (98.3), y Δ_r = S_r − M_r con S_r por la ruta sin caracteres, que no usa (98.3).
- No se pidió revisar la prueba escrita. Aun así comprobé su álgebra (§2).
- Script: `numerico/verificacion_107_claude_code.py`. Salida: `numerico/resultados/20261128_verificacion_107_claude_code.txt`. **12/12 PASS**, 6 s.

**Estado:** PROPUESTA hasta el contraste de Muse.

**Goldbach sigue abierta.**

---

## §0. Veredicto

**Todas las calibraciones del doc 107 se reproducen. No encontré ningún error que afecte al teorema.**

Hay dos precisiones de redacción que conviene incorporar antes de seguir citándolo en público:
1. **El conjunto de las calibraciones es ℋ_N ∩ {r > ln N}, no ℋ.** Las sumas "Σ_{r∈ℋ}" de K1 y K3 se reproducen solo con r > ln N, el conjunto del segundo factor. Sobre todo ℋ_N, B_Δ/T en gen/S23 vale 2,161·10⁻⁴ en lugar de 2,159·10⁻⁴; en R1 la diferencia es un factor 1,57.

   El **teorema no cambia**: la cota (107.1) mayora la suma completa, que contiene a la restringida. Pero la tabla debe decir sobre qué conjunto suma.
2. **"El cociente es decreciente para N > 5970" es verdadero, pero el umbral exacto es e^{8/(1−2c)} = 5964,4.** No es un error, solo una cota no óptima.

---

## §1. Chequeo por chequeo

| Chequeo | Afirmación del doc 107 | Mi cálculo | Veredicto |
|---|---|---|---|
| K0 | Maquinaria: \|ℋ\| de R2 = 880 y familia | R2 880, S11 802, S13 1628, S17 5033, S23 17999 | **Confirmado** |
| K1a | máx \|Δ_r\|/(φ(r)ω(r)ln N) = 0,0513 en r = 17 (Δ₁₇ = 11,3329) | 0,0513 en r = 17, Δ₁₇ = 11,3329, sobre r > ln N. Sobre todo ℋ_N el máximo es 0,2330 en r = 5, también ≤ 1 | **Confirmado**; la cota (106.6) vale en ambos conjuntos |
| K1b | Contraejemplo a la cota impresa del doc 97: r = 458815 ∈ ℋ(S23), \|Δ_r\| = 76522,07, factor 1846,3, 21 % de B_Δ | r = 5·7·13109, φ = 314592, ω = 3. Δ = **76522,068094** por (98.3) y **76522,068094** por S_r − M_r. Factor **1846,3**. Δ²/φ = 1,861·10⁴ = **21,3 %** de B_Δ(gen/S23) | **Confirmado por dos rutas independientes** |
| K1c | B_Δ ≤ (ln N)²·M_ω en los 12 puntos, con holgura ≥ 2·10⁴ | holgura mínima 2,06·10⁴ con r > ln N (1,02·10⁴ con todo ℋ_N) | **Confirmado** en la convención del doc |
| K2 | M_ω y M_φ por ramas = enumeración; valores gen S11..S23 | por ramas = enumeración con error relativo < 10⁻¹²; M_ω = 7,137·10⁶ / 3,949·10⁷ / 5,640·10⁸ / 1,726·10⁹ / 1,0685·10¹⁰; M_φ = 1,6385·10⁶ / 7,427·10⁶ / 8,949·10⁷ / 2,521·10⁸ / 1,3705·10⁹ | **Confirmado** |
| K3 | máx B_Δ/T = 2,159·10⁻⁴ (gen/S23); R3 ≈ 2,3·10⁻⁷, R1 ≈ 4,4·10⁻⁷; C_Δ = 10⁻³ cubre con factor ≥ 4,6 | T(999982) = 4,0526·10⁸ y T(2²⁰) = 4,3787·10⁸. Máximo **2,159·10⁻⁴** en gen/S23; R3 2,32·10⁻⁷; R1 4,37·10⁻⁷; factor 4,63 | **Confirmado** (con r > ln N) |
| K4 | Masas por nivel z = 11..29 suman 1; \|𝒟_C\| = 577 (S11) … 47102 (S29) | tabla reproducida con desviación máxima 0,0005 sobre ℋ_N (N = 999982); \|𝒟_C\| = **577** y **47102** | **Confirmado** (la tabla usa ℋ_N, no ℋ) |
| K4b | Calibración con ψ: S2+S3 = 0,906 / 0,908 / 0,877 / 0,846 / 0,809 | 0,906 / 0,908 / 0,877 / 0,845 / 0,809; S4+S5 hasta 0,144; L ≤ 0,0045 | **Confirmado** (0,845 frente a 0,846, por redondeo) |
| K5 | W* = 9,52 / 15,00 / 31,00; Σ_ℓ W*_ℓ = 28,0 → 189,6; la cota clásica no alcanza T | sobre ℋ_N: W* = 9,52 / 15,00 / 31,00 / 31,00 / 31,00 (32,37 en z = 29); Σ = 28,0 … 189,6; exceso solo-LS 3,39·10⁷ … 2,72·10⁹, > 1 en todo z | **Confirmado** |
| K6 | Identidad (107.2); N₀ = 3,82·10⁹ (C_Δ = 1), 7,40·10¹² (10⁻²), 2,36·10¹⁴ (10⁻³) | constantes 0,520342 y 0,0117077; exponente 0,920232; N₀ = **3,822·10⁹ / 7,398·10¹² / 2,361·10¹⁴** por bisección propia | **Confirmado** |
| (107.1) | B_Δ ≤ 0,520342·N(ln N)⁴ | se cumple en los 12 puntos | **Confirmado** |

**Descomposición por niveles** (la que citas: "S2–S5 ≥ 97 %; d > D_A^{2/3} entre 80 y 97 %"). Confirmada:
- La masa en S2–S5 vale 0,981 (z = 11) → 0,996 (z = 29).
- La de d > D_A^{2/3} vale 0,798 → 0,973.
- La de d ≤ D_A^{1/3} es 0,000 en todo z.
- La de los niveles lisos es ≤ 0,0002.

Matiz de redacción: "≥ 97 %" es S2–S5 en sentido estricto, porque S1 se excluye. Con S1 incluido la masa conmutada pasa del 99,8 %.

---

## §2. El álgebra de la cadena [T]

Recompuse cada paso del §1 del doc 107:
1. **Lema (106.6)**: |Δ_r| ≤ φ(r)ω(r)ln N. Se sigue de (98.3) por la desigualdad triangular con:
   - |g_p| ≤ ⌊log_p N⌋ log p ≤ ln N;
   - #{primitivos mod e} ≤ φ(e);
   - φ(r/p) ≤ φ(r) cuando p ∥ r;
   - #{p ∥ r} ≤ ω(r).

   Correcto. El caso r primo (r/p = 1, término principal) está cubierto con φ(1) = 1 ≤ φ(r).
2. **B_Δ ≤ (ln N)²·M_ω.** Es sustitución término a término. Correcto.
3. **(107.1).** Usa r < D_A ≤ N^{1/2}, φ(r) ≤ r, ω(r) ≤ log₂ r < ln N/(2 ln 2) y |ℋ| ≤ D_A. Da B_Δ ≤ N(ln N)⁴/(4 ln²2) = 0,520342·N(ln N)⁴. Correcto, y vale tanto con D_A = N^{1/2} como con el D_A = X^{1/2}/(log X)^{C0} del doc 94, que es menor.
4. **(107.2).** 0,520342 × 0,0225 = 0,0117077, y el exponente es −(1−2c) = −0,920232. Correcto.

**Observación, sin efecto sobre el teorema.** El doc 107 §1.4 sugiere que la fórmula A100.1 "bajaría la escala de B_Δ a polilog·N^{0,9208}" si se convirtiera en cota. Por el Teorema 101-CC, asintóticamente |ℋ| = N^{1/2−o(1)}: el conteo trivial (iii) |ℋ| ≤ D_A es **óptimo en exponente**, y no hay ganancia posible por esa vía. La prueba del doc 107 no la necesita, porque cierra con el conteo trivial y margen N^{0,92}.

---

## §3. Erratas y avisos para la autora (sin tocar sus archivos)

1. **Doc 107 §0 y §1.5:** declarar que B_Δ y las calibraciones K1 y K3 suman sobre ℋ_N ∩ {r > ln N}. El enunciado del teorema queda igual.
2. **Doc 107 §1.5:** "decreciente para N > 5970" se puede afinar a N > 5964,4.
3. **Doc 107 §2.1:** la tabla de niveles está calculada sobre ℋ_N con N = 999982. Sobre ℋ sin filtro, |𝒟_C| es 581 y no 577, y las fracciones varían en la cuarta cifra.
4. **Doc 107 §1.4:** la frase sobre polilog·N^{0,9208} debe matizarse con el Teorema 101-CC: no hay ganancia de exponente por conteo.

*Goldbach fuerte sigue abierta.*
