/-
Lean Hito 8 (Goldbach loop): elementary support properties of the
outer Mobius-logarithm weight in the finite dyadic B3 box.

This file proves no analytic estimate and does not prove Goldbach.
-/

import GoldbachLumina.VaughanB3Doc06Weights

open ArithmeticFunction

namespace GoldbachLumina

/-- The Mobius-logarithm weight vanishes at zero. -/
@[simp] theorem b3MoebiusLogWeight_zero : b3MoebiusLogWeight 0 = 0 := by
  simp [b3MoebiusLogWeight]

/-- The Mobius-logarithm weight vanishes at one. -/
@[simp] theorem b3MoebiusLogWeight_one : b3MoebiusLogWeight 1 = 0 := by
  simp [b3MoebiusLogWeight]

/-- A nonsquarefree outer variable has zero Mobius-logarithm weight. -/
theorem b3MoebiusLogWeight_eq_zero_of_not_squarefree {d : Nat}
    (hd : ¬Squarefree d) : b3MoebiusLogWeight d = 0 := by
  simp [b3MoebiusLogWeight,
    ArithmeticFunction.moebius_eq_zero_of_not_squarefree hd]

/-- Nonzero outer weight forces squarefree support. -/
theorem squarefree_of_b3MoebiusLogWeight_ne_zero {d : Nat}
    (hd : b3MoebiusLogWeight d ≠ 0) : Squarefree d := by
  by_contra hsq
  exact hd (b3MoebiusLogWeight_eq_zero_of_not_squarefree hsq)

/-- Nonzero outer weight also excludes the exceptional indices 0 and 1. -/
theorem two_le_of_b3MoebiusLogWeight_ne_zero {d : Nat}
    (hd : b3MoebiusLogWeight d ≠ 0) : 2 ≤ d := by
  have h0 : d ≠ 0 := by
    intro h
    subst d
    exact hd b3MoebiusLogWeight_zero
  have h1 : d ≠ 1 := by
    intro h
    subst d
    exact hd b3MoebiusLogWeight_one
  omega

#print axioms GoldbachLumina.b3MoebiusLogWeight_zero
#print axioms GoldbachLumina.b3MoebiusLogWeight_one
#print axioms GoldbachLumina.b3MoebiusLogWeight_eq_zero_of_not_squarefree
#print axioms GoldbachLumina.squarefree_of_b3MoebiusLogWeight_ne_zero
#print axioms GoldbachLumina.two_le_of_b3MoebiusLogWeight_ne_zero

end GoldbachLumina
