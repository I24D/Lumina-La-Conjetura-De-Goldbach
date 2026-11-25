/-
Lean Hito 3 (bucle Goldbach): exact finite expansion of Vaughan's Type-II
piece.  This file adds no analytic estimate and does not prove Goldbach.
-/

import GoldbachLumina.VaughanDecomp

open ArithmeticFunction
open scoped ArithmeticFunction Moebius zeta

namespace GoldbachLumina

/-- Exact double-convolution expansion of Vaughan's Type-II term.

The outer factorisation is `n = m*d`; the inner one is `m = a*b`.
Thus the formula retains all three finite variables which later become the
dyadic `B3` variables.  The cutoffs remain explicit in `lambdaLE U` and
`moebiusLE V`; no estimate or rearrangement of an infinite sum is used. -/
theorem lambdaFlat_apply_expanded (U V n : ℕ) :
    lambdaFlat U V n =
      ∑ md ∈ n.divisorsAntidiagonal,
        (∑ ab ∈ md.1.divisorsAntidiagonal,
          (ArithmeticFunction.vonMangoldt ab.1 - lambdaLE U ab.1) *
            ArithmeticFunction.zeta ab.2) *
        ((ArithmeticFunction.moebius : ArithmeticFunction ℝ) md.2 -
          moebiusLE V md.2) := by
  simp [lambdaFlat, ArithmeticFunction.mul_apply]

#print axioms GoldbachLumina.lambdaFlat_apply_expanded

end GoldbachLumina
