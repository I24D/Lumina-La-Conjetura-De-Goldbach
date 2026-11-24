# lean/ — copia de las fuentes Lean (solo lectura)

La copia **canónica**, que compila, está en `C:\Users\dal_n\lean\goldbach_lumina\GoldbachLumina\`.
El proyecto está fuera de OneDrive porque la caché de Mathlib ocupa unos 5 GB. Para modificar algo, edita allí y vuelve a copiarlo aquí.

Compilar: `cd C:\Users\dal_n\lean\goldbach_lumina` y luego `lake env lean GoldbachLumina/Reduccion.lean`

| Teorema (Reduccion.lean) | Qué formaliza | Axiomas (#print axioms) |
|---|---|---|
| `R_eq_neg_sum` | R(N) = −Σ_{0<n<N} Σ_{d∣n} μ(d) log d Λ(N−n) | estándar |
| `R_eq_S_sub_T` | **Prop. 4.1 de 04**: R(N) = S_N(D) − T_N(D), para todo N y todo D (forma de divisores) | propext, Classical.choice, Quot.sound |
| `coprime_of_coprime_sub` | (N−cd, N) = 1 ⇒ (c, N) = 1 (06, Lema 6.1, paso 1) | ídem |
| `primePow_of_not_coprime` | Λ(m) ≠ 0 y (m,N) > 1 ⇒ m = p^k con p ∣ N (núcleo de la Prop. 4.2 de 04) | ídem |
| `mu_log_prime` | μ(p) log p = −log p (06, Lema 6.0) | ídem |
| `divisors_prime_gt` | para p primo y 1 ≤ D < p, los divisores de p mayores que D son {p} (06, R-3) | ídem |

Procedencia: Claude Code (claude-opus-5-5), 2026-11-23. Estado: [T-Lean], verificado por el kernel sin `sorry`. **No demuestra Goldbach.**

## Continuación de Codex — 2026-11-23

Procedencia: Codex; identificador exacto de modelo no expuesto. Método: pruebas nuevas en Lean 4.37.0-rc1 / Mathlib, en `GoldbachLumina/Reduccion2.lean`, sin modificar la fuente anterior. Estado **[T-Lean]** para los siguientes enunciados finitos. Copia idéntica sin BOM: [Reduccion2.lean](Reduccion2.lean). Comando: `lake env lean GoldbachLumina/Reduccion2.lean`, exit **0**, sin errores ni advertencias en la ejecución final.

Todos los nombres siguientes llevan el prefijo `GoldbachLumina.Reduction2.`.

| Teorema (Reduccion2.lean) | Qué formaliza | Axiomas (#print axioms) |
|---|---|---|
| `reindex_divisor_tail` | Biyectivamente n=cd; límites completos para N>0, todo D y cualquier peso real f(n,d) | propext, Classical.choice, Quot.sound |
| `reindex_mobius_tail` | Especialización a μ(d)log d·Λ(N−cd), cola completa sobre todas las c | ídem |
| `proper_prime_power_error` | 0≤R−R*≤2⌊√N⌋⌊log₂N⌋(log N)² para N≥1, desigualdad finita (4.12) | ídem |
| `primeConvolution_eq_log` | R* coincide con la suma de pesos logarítmicos sobre dos primos | ídem |
| `badTail_eq_product` | La definición de la cola mala en divisores coincide con (c,N)>1 en pares (c,d) | ídem |
| `bad_tail_bound` | Primera desigualdad explícita de Prop. 4.2: \|T^bad\|≤K_N, para N≥1 y todo D | ídem |

Registro de salida real y SHA-256: [2026-11-23_codex_Reduccion2_verificacion.json](2026-11-23_codex_Reduccion2_verificacion.json). Solo aparecen los tres axiomas estándar; no `sorryAx`, axiomas nuevos ni huecos de prueba.

Alcance: no se formalizan aquí o(N) en (4.12), la estimación asintótica de K_N, la segunda desigualdad de Prop. 4.2 sobre divisores pequeños, ni los lemas analíticos de 06/07. La revisión escrita y sus errores explícitos están en [08](../matematicas/08_revision_codex_06_07.md). **Goldbach sigue abierta.**
