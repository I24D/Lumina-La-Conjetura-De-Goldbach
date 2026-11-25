/-
Lean Hito 6 (bucle Goldbach): finite dyadic boxes for the Type-II term B3.

This file only packages the exact bridge from `VaughanB3Bridge` inside
finite outer sums.  It proves no analytic estimate and does not prove
Goldbach.
-/

import GoldbachLumina.VaughanB3Bridge

open ArithmeticFunction

namespace GoldbachLumina

/-- Membership in the half-open dyadic interval `(X, 2X]`. -/
def inDyadic (X n : Nat) : Prop := X < n ∧ n ≤ 2 * X

instance (X n : Nat) : Decidable (inDyadic X n) := by
  unfold inDyadic
  infer_instance

/-- A finite dyadic Type-II box written with `lambdaFlat V U`.

The cutoff `c*d < N` makes `N-c*d` the positive additive complement.
The weights are kept arbitrary, so the identity below can later be
specialized to the coprimality and Mobius-log weights of doc 06. -/
noncomputable def b3DyadicFlat
    (N C K V U : Nat) (alpha beta : Nat → Real) : Real :=
  ∑ c ∈ Finset.range N,
    if inDyadic C c then
      alpha c * ∑ d ∈ Finset.range N,
        if inDyadic K d ∧ c * d < N then
          beta d * lambdaFlat V U (N - c * d)
        else 0
    else 0

/-- The same finite box after expanding its Type-II inner factor through
the exact `lambdaFlat`/`B3` bridge. -/
noncomputable def b3DyadicBridge
    (N C K V U : Nat) (alpha beta : Nat → Real) : Real :=
  ∑ c ∈ Finset.range N,
    if inDyadic C c then
      alpha c * ∑ d ∈ Finset.range N,
        if inDyadic K d ∧ c * d < N then
          beta d *
            (∑ ar ∈ (N - c * d).divisorsAntidiagonal,
              (if V < ar.1 then ArithmeticFunction.vonMangoldt ar.1 else 0) *
                nuAbove U ar.2)
        else 0
    else 0

/-- Exact equality of the finite dyadic box with its classical B3
expansion.  No cancellation or asymptotic input is used. -/
theorem b3DyadicFlat_eq_bridge
    (N C K V U : Nat) (alpha beta : Nat → Real) :
    b3DyadicFlat N C K V U alpha beta =
      b3DyadicBridge N C K V U alpha beta := by
  simp only [b3DyadicFlat, b3DyadicBridge, lambdaFlat_b3_bridge]

#print axioms GoldbachLumina.b3DyadicFlat_eq_bridge

end GoldbachLumina
