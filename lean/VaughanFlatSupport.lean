import GoldbachLumina.VaughanFlatExpansion

/-! # Pointwise support for Vaughan's Type-II factors

Lean Hito 4 (Goldbach loop): pointwise support of the two truncated
differences occurring in Vaughan's Type-II term.

This file adds no analytic estimate and does not prove Goldbach.
-/

open ArithmeticFunction
open scoped ArithmeticFunction Moebius zeta

namespace GoldbachLumina

/-- The large von Mangoldt factor vanishes at and below its cutoff. -/
theorem vonMangoldt_sub_lambdaLE_eq_zero_of_le
    (U a : Nat) (ha : a <= U) :
    ArithmeticFunction.vonMangoldt a - lambdaLE U a = 0 := by
  by_cases hzero : a = 0
  · subst a
    simp [lambdaLE, ArithmeticFunction.restrictTo]
  · have hone : 1 <= a := Nat.one_le_iff_ne_zero.mpr hzero
    simp [lambdaLE, ArithmeticFunction.restrictTo, Set.mem_Icc, hone, ha]

/-- Strictly above the cutoff, the large von Mangoldt factor is unchanged. -/
theorem vonMangoldt_sub_lambdaLE_eq_of_lt
    (U a : Nat) (ha : U < a) :
    ArithmeticFunction.vonMangoldt a - lambdaLE U a =
      ArithmeticFunction.vonMangoldt a := by
  simp [lambdaLE, ArithmeticFunction.restrictTo, Set.mem_Icc,
    Nat.not_le.mpr ha]

/-- The large Moebius factor vanishes at and below its cutoff. -/
theorem moebius_sub_moebiusLE_eq_zero_of_le
    (V d : Nat) (hd : d <= V) :
    (ArithmeticFunction.moebius : ArithmeticFunction Real) d -
        moebiusLE V d = 0 := by
  by_cases hzero : d = 0
  · subst d
    simp [moebiusLE, ArithmeticFunction.restrictTo]
  · have hone : 1 <= d := Nat.one_le_iff_ne_zero.mpr hzero
    simp [moebiusLE, ArithmeticFunction.restrictTo, Set.mem_Icc, hone, hd]

/-- Strictly above the cutoff, the large Moebius factor is unchanged. -/
theorem moebius_sub_moebiusLE_eq_of_lt
    (V d : Nat) (hd : V < d) :
    (ArithmeticFunction.moebius : ArithmeticFunction Real) d -
        moebiusLE V d =
      (ArithmeticFunction.moebius : ArithmeticFunction Real) d := by
  simp [moebiusLE, ArithmeticFunction.restrictTo, Set.mem_Icc,
    Nat.not_le.mpr hd]

#print axioms GoldbachLumina.vonMangoldt_sub_lambdaLE_eq_zero_of_le
#print axioms GoldbachLumina.vonMangoldt_sub_lambdaLE_eq_of_lt
#print axioms GoldbachLumina.moebius_sub_moebiusLE_eq_zero_of_le
#print axioms GoldbachLumina.moebius_sub_moebiusLE_eq_of_lt

end GoldbachLumina
