import Mathlib.Algebra.BigOperators.Group.Finset.Basic
import Mathlib.Tactic

open scoped BigOperators

namespace GoldbachLumina

/-!
# Scaled marginal projection of a finite matrix

This is a denominator-free form of the canonical decomposition in document 47.
It formalizes only a finite algebraic identity, not an analytic estimate.
-/

def matrixRowSum
    {ι κ R : Type*} [Fintype κ] [CommRing R]
    (c : ι → κ → R) (i : ι) : R :=
  ∑ j, c i j

def matrixColSum
    {ι κ R : Type*} [Fintype ι] [CommRing R]
    (c : ι → κ → R) (j : κ) : R :=
  ∑ i, c i j

def matrixTotal
    {ι κ R : Type*} [Fintype ι] [Fintype κ] [CommRing R]
    (c : ι → κ → R) : R :=
  ∑ i, ∑ j, c i j

/-- The numerator `mn * c⁰` of the doubly centered matrix from (47.1). -/
def scaledInteraction
    {ι κ R : Type*} [Fintype ι] [Fintype κ] [CommRing R]
    (c : ι → κ → R) (i : ι) (j : κ) : R :=
  ((Fintype.card ι * Fintype.card κ : ℕ) : R) * c i j
    - (Fintype.card ι : R) * matrixRowSum c i
    - (Fintype.card κ : R) * matrixColSum c j
    + matrixTotal c

theorem sum_matrixColSum_eq_matrixTotal
    {ι κ R : Type*} [Fintype ι] [Fintype κ] [CommRing R]
    (c : ι → κ → R) :
    (∑ j, matrixColSum c j) = matrixTotal c := by
  simp only [matrixColSum, matrixTotal]
  rw [Finset.sum_comm]

theorem scaledInteraction_reconstruct
    {ι κ R : Type*} [Fintype ι] [Fintype κ] [CommRing R]
    (c : ι → κ → R) (i : ι) (j : κ) :
    ((Fintype.card ι * Fintype.card κ : ℕ) : R) * c i j =
      scaledInteraction c i j
        + (Fintype.card ι : R) * matrixRowSum c i
        + (Fintype.card κ : R) * matrixColSum c j
        - matrixTotal c := by
  simp [scaledInteraction, Nat.cast_mul]
  ring

theorem scaledInteraction_zeroRows
    {ι κ R : Type*} [Fintype ι] [Fintype κ] [CommRing R]
    (c : ι → κ → R) (i : ι) :
    (∑ j, scaledInteraction c i j) = 0 := by
  simp_rw [scaledInteraction, Finset.sum_add_distrib, Finset.sum_sub_distrib]
  simp_rw [← Finset.mul_sum]
  simp [matrixRowSum, sum_matrixColSum_eq_matrixTotal, Nat.cast_mul]
  ring

theorem scaledInteraction_zeroCols
    {ι κ R : Type*} [Fintype ι] [Fintype κ] [CommRing R]
    (c : ι → κ → R) (j : κ) :
    (∑ i, scaledInteraction c i j) = 0 := by
  simp_rw [scaledInteraction, Finset.sum_add_distrib, Finset.sum_sub_distrib]
  simp_rw [← Finset.mul_sum]
  simp [matrixRowSum, matrixColSum, matrixTotal, Nat.cast_mul]
  ring

#print axioms GoldbachLumina.scaledInteraction_reconstruct
#print axioms GoldbachLumina.scaledInteraction_zeroRows
#print axioms GoldbachLumina.scaledInteraction_zeroCols

end GoldbachLumina
