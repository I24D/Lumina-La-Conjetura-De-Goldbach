/-
Lean Hito 5 (bucle Goldbach): puente exacto entre `lambdaFlat` y la
notación clásica del término de Tipo II `B3` del doc 06.

`lambdaFlat V U = (Λ - Λ_{≤V}) * ζ * (μ - μ_{≤U})` se reagrupa
exactamente como

  `lambdaFlat V U n = ∑_{a * r = n, a > V} Λ(a) * ν_U(r)`,

donde `ν_U(r) = ∑_{t ∣ r, t > U} μ(t)` (doc 23, (23.1); doc 06, §6).
La condición `b > U` del doc 06 es automática porque `ν_U` se anula
en `[0, U]`.

Este archivo no aporta ninguna cota analítica y no demuestra la
conjetura de Goldbach.
-/

import GoldbachLumina.VaughanFlatSupport

open ArithmeticFunction
open scoped ArithmeticFunction Moebius zeta

namespace GoldbachLumina

/-- The classical truncated Möbius weight `ν_U = ζ * (μ - μ_{≤U})` of
Vaughan's Type-II term (doc 06, §6; doc 23, (23.1)). -/
noncomputable def nuAbove (U : ℕ) : ArithmeticFunction ℝ :=
  ArithmeticFunction.zeta *
    ((ArithmeticFunction.moebius : ArithmeticFunction ℝ) - moebiusLE U)

/-- Pointwise form: `ν_U(r) = ∑_{t ∣ r} (μ - μ_{≤U})(t)`, i.e. only the
divisors `t > U` of `r` contribute `μ(t)`. -/
theorem nuAbove_apply (U r : ℕ) :
    nuAbove U r =
      ∑ t ∈ r.divisors,
        ((ArithmeticFunction.moebius : ArithmeticFunction ℝ) t - moebiusLE U t) :=
  ArithmeticFunction.coe_zeta_mul_apply

/-- `ν_U` vanishes on `[0, U]`: the extra `b > U` condition written in
doc 06 for the second variable of `B3` is automatic. -/
theorem nuAbove_eq_zero_of_le (U r : ℕ) (hr : r ≤ U) :
    nuAbove U r = 0 := by
  rw [nuAbove_apply]
  refine Finset.sum_eq_zero fun t ht => ?_
  rcases Nat.mem_divisors.mp ht with ⟨hdvd, hr0⟩
  have htr : t ≤ r := Nat.le_of_dvd (Nat.pos_of_ne_zero hr0) hdvd
  have htu : t ≤ U := Nat.le_trans htr hr
  exact moebius_sub_moebiusLE_eq_zero_of_le U t htu

/-- Step 1 of the bridge: exact re-grouping of the triple Dirichlet
convolution defining `lambdaFlat V U` around `ν_U`. -/
theorem lambdaFlat_eq_nuAbove (V U n : ℕ) :
    lambdaFlat V U n =
      ∑ ar ∈ n.divisorsAntidiagonal,
        (ArithmeticFunction.vonMangoldt ar.1 - lambdaLE V ar.1) * nuAbove U ar.2 := by
  have h : lambdaFlat V U =
      (ArithmeticFunction.vonMangoldt - lambdaLE V) * nuAbove U := by
    simp only [lambdaFlat, nuAbove, mul_assoc]
  simp [h, ArithmeticFunction.mul_apply]

/-- The exact bridge (doc 23, (23.1)):
`lambdaFlat V U n = ∑_{a * r = n, a > V} Λ(a) * ν_U(r)`.
This is the inner summand of the classical Type-II term `B3`. -/
theorem lambdaFlat_b3_bridge (V U n : ℕ) :
    lambdaFlat V U n =
      ∑ ar ∈ n.divisorsAntidiagonal,
        (if V < ar.1 then ArithmeticFunction.vonMangoldt ar.1 else 0) * nuAbove U ar.2 := by
  rw [lambdaFlat_eq_nuAbove]
  refine Finset.sum_congr rfl fun ar _ => ?_
  by_cases hV : V < ar.1
  · rw [vonMangoldt_sub_lambdaLE_eq_of_lt V ar.1 hV]
    simp [hV]
  · rw [vonMangoldt_sub_lambdaLE_eq_zero_of_le V ar.1 (Nat.not_lt.mp hV)]
    simp [hV]

#print axioms GoldbachLumina.nuAbove_apply
#print axioms GoldbachLumina.nuAbove_eq_zero_of_le
#print axioms GoldbachLumina.lambdaFlat_eq_nuAbove
#print axioms GoldbachLumina.lambdaFlat_b3_bridge

end GoldbachLumina
