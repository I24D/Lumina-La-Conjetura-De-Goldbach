import Mathlib.Algebra.BigOperators.Group.Finset.Basic
import Mathlib.Tactic

open scoped BigOperators

namespace GoldbachLumina

/-!
# Exact cancellation of first-order local perturbations

This is the finite algebraic identity used in document 46 of the Goldbach
workspace.  It is deliberately independent of any analytic estimate.
-/

theorem marginalCancellation
    {ι κ R : Type*} [Fintype ι] [Fintype κ] [CommRing R]
    (c : ι → κ → R) (a : ι → R) (b : κ → R)
    (hrow : ∀ i, ∑ j, c i j = 0)
    (hcol : ∀ j, ∑ i, c i j = 0) :
    (∑ i, ∑ j, c i j * (1 + a i) * (1 + b j)) =
      ∑ i, ∑ j, c i j * a i * b j := by
  have htotal : (∑ i, ∑ j, c i j) = 0 := by
    simp [hrow]
  have ha : (∑ i, ∑ j, c i j * a i) = 0 := by
    simp_rw [← Finset.sum_mul]
    simp [hrow]
  have hb : (∑ i, ∑ j, c i j * b j) = 0 := by
    rw [Finset.sum_comm]
    simp_rw [← Finset.sum_mul]
    simp [hcol]
  have hpoint : ∀ i j,
      c i j * (1 + a i) * (1 + b j) =
        c i j + c i j * a i + c i j * b j + c i j * a i * b j := by
    intro i j
    ring
  simp_rw [hpoint, Finset.sum_add_distrib]
  rw [htotal, ha, hb]
  simp

theorem zeroTotal_of_zeroRows
    {ι κ R : Type*} [Fintype ι] [Fintype κ] [AddCommMonoid R]
    (c : ι → κ → R) (hrow : ∀ i, ∑ j, c i j = 0) :
    (∑ i, ∑ j, c i j) = 0 := by
  simp [hrow]

#print axioms GoldbachLumina.marginalCancellation
#print axioms GoldbachLumina.zeroTotal_of_zeroRows

end GoldbachLumina
