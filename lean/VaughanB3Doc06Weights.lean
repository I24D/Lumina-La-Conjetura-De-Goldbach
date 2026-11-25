/-
Lean Hito 7 (Goldbach loop): specialize the finite dyadic B3 box to
the two outer weights used in project document 06.

This file proves only an exact finite identity. It contains no analytic
estimate, does not overcome the parity barrier, and does not prove
Goldbach.
-/

import GoldbachLumina.VaughanB3DyadicBox

open ArithmeticFunction

namespace GoldbachLumina

/-- The outer c-weight from doc 06: the indicator of `(c,N)=1`. -/
noncomputable def b3CoprimeWeight (N c : Nat) : Real :=
  if Nat.Coprime c N then 1 else 0

/-- The outer d-weight from doc 06: `mu(d) log(d)`. -/
noncomputable def b3MoebiusLogWeight (d : Nat) : Real :=
  (ArithmeticFunction.moebius d : Real) * Real.log d

/-- The doc-06 dyadic box, using its literal outer weights and the
classical inner `B3` bridge. -/
noncomputable def b3Doc06Dyadic (N C K V U : Nat) : Real :=
  b3DyadicBridge N C K V U (b3CoprimeWeight N) b3MoebiusLogWeight

/-- Specialization of the exact finite dyadic identity to the two
outer weights occurring in doc 06. -/
theorem b3DyadicFlat_doc06_eq
    (N C K V U : Nat) :
    b3DyadicFlat N C K V U (b3CoprimeWeight N) b3MoebiusLogWeight =
      b3Doc06Dyadic N C K V U := by
  exact b3DyadicFlat_eq_bridge N C K V U
    (b3CoprimeWeight N) b3MoebiusLogWeight

#print axioms GoldbachLumina.b3DyadicFlat_doc06_eq

end GoldbachLumina
