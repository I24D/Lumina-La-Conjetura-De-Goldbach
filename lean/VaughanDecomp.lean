/-
Lean Hito 2 (GB-L19): algebraic core of Vaughan's identity.

This is a Mathlib-only packaging of `Lambda_decomp` from
amellendijk/lean-bombieri-vinogradov (commit b1b0904, BV/Defs.lean).
The original theorem uses project-global cutoffs U,V and infrastructure
unrelated to the identity.  Here the cutoffs are explicit natural numbers.

No analytic estimate is formalized: this is only the exact arithmetic-
function decomposition.  In particular it does not prove Goldbach.
-/

import Mathlib

open ArithmeticFunction
open scoped ArithmeticFunction Moebius zeta

namespace GoldbachLumina

/-- Restrict an arithmetic function to a set. -/
noncomputable def ArithmeticFunction.restrictTo {R : Type*} [Zero R]
    (s : Set ℕ) (f : ArithmeticFunction R) : ArithmeticFunction R :=
  ⟨s.indicator f, by simp⟩

/-- The small von Mangoldt piece `1_{1 ≤ n ≤ U} · Λ(n)`. -/
noncomputable def lambdaLE (U : ℕ) : ArithmeticFunction ℝ :=
  ArithmeticFunction.restrictTo (Set.Icc 1 U) ArithmeticFunction.vonMangoldt

/-- The small Moebius piece `1_{1 ≤ n ≤ V} · μ(n)`. -/
noncomputable def moebiusLE (V : ℕ) : ArithmeticFunction ℝ :=
  ArithmeticFunction.restrictTo (Set.Icc 1 V) (ArithmeticFunction.moebius : ArithmeticFunction ℝ)

/-- Vaughan's Type-I (sharp) arithmetic-function piece. -/
noncomputable def lambdaSharp (U V : ℕ) : ArithmeticFunction ℝ :=
  moebiusLE V * ArithmeticFunction.log
    - lambdaLE U * moebiusLE V * ArithmeticFunction.zeta

/-- Vaughan's Type-II (flat) arithmetic-function piece. -/
noncomputable def lambdaFlat (U V : ℕ) : ArithmeticFunction ℝ :=
  (ArithmeticFunction.vonMangoldt - lambdaLE U) * ArithmeticFunction.zeta
    * ((ArithmeticFunction.moebius : ArithmeticFunction ℝ) - moebiusLE V)

/-- Exact Vaughan decomposition, with no size assumptions on the cutoffs.

This is the formal algebraic identity
`Λ = Λ♯ + Λ♭ + Λ_{≤U}`.  The Type-I/Type-II estimates needed
for Bombieri--Vinogradov or Goldbach are separate analytic statements. -/
theorem vaughan_decomp (U V n : ℕ) :
    ArithmeticFunction.vonMangoldt n =
      lambdaSharp U V n + lambdaFlat U V n + lambdaLE U n := by
  simp_rw [← add_apply]
  congr 1
  simp [lambdaFlat, lambdaSharp, lambdaLE,
    ← ArithmeticFunction.vonMangoldt_mul_zeta]
  have hzmu :
      (ArithmeticFunction.zeta *
        (ArithmeticFunction.moebius : ArithmeticFunction ℝ)) = 1 :=
    ArithmeticFunction.coe_zeta_mul_coe_moebius
  grind

#print axioms GoldbachLumina.vaughan_decomp

end GoldbachLumina
