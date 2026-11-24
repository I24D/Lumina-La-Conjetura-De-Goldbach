/-
Goldbach-LUMINA: finite reindexing and elementary error bounds.
Provenance: Codex (exact model identifier unavailable), 2026-11-23 ET.
Method: Lean 4.37.0-rc1 / Mathlib; kernel checked, no custom axioms.
Status: [T-Lean] applies only to declarations that compile and pass #print axioms.
This file does not prove Goldbach or any uniform analytic cancellation estimate.
Canonical source outside OneDrive; does not modify Reduccion.lean.
-/
import Mathlib

open Finset
noncomputable section

namespace GoldbachLumina.Reduction2

def divisorPairs (N D : ℕ) : Finset (Σ _ : ℕ, ℕ) :=
  (Ioo 0 N).sigma fun n => n.divisors.filter (fun d => D < d)

def productPairs (N D : ℕ) : Finset (Σ _ : ℕ, ℕ) :=
  (Icc 1 ((N - 1) / (D + 1))).sigma fun c => Icc (D + 1) ((N - 1) / c)

theorem product_mem_iff {N D c d : ℕ} (hN : 0 < N) :
    (⟨c, d⟩ : Σ _ : ℕ, ℕ) ∈ productPairs N D ↔
      0 < c ∧ D < d ∧ c * d < N := by
  simp only [productPairs, mem_sigma, mem_Icc]
  constructor
  · rintro ⟨⟨hc, hcmax⟩, hd, hdmax⟩
    have hmul := (Nat.le_div_iff_mul_le hc).mp hdmax
    have : c * d ≤ N - 1 := by simpa [Nat.mul_comm] using hmul
    omega
  · rintro ⟨hc, hd, hcd⟩
    have hcd' : c * d ≤ N - 1 := by omega
    have hcmax : c ≤ (N - 1) / (D + 1) := by
      apply (Nat.le_div_iff_mul_le (by omega : 0 < D + 1)).mpr
      exact le_trans (Nat.mul_le_mul_left c (by omega : D + 1 ≤ d)) hcd'
    have hdmax : d ≤ (N - 1) / c := by
      apply (Nat.le_div_iff_mul_le hc).mpr
      simpa [Nat.mul_comm] using hcd'
    exact ⟨⟨by omega, hcmax⟩, by omega, hdmax⟩

theorem mem_divisorPairs {N D n d : ℕ} :
    (⟨n, d⟩ : Σ _ : ℕ, ℕ) ∈ divisorPairs N D ↔
      0 < n ∧ n < N ∧ d ∣ n ∧ D < d := by
  simp only [divisorPairs, mem_sigma, mem_Ioo, mem_filter, Nat.mem_divisors]
  omega

/-- Exact change of variables n = c*d, including every positive c. [T-Lean] -/
theorem reindex_divisor_tail (N D : ℕ) (hN : 0 < N) (f : ℕ → ℕ → ℝ) :
    (∑ n ∈ Ioo 0 N, ∑ d ∈ n.divisors.filter (fun d => D < d), f n d) =
      ∑ c ∈ Icc 1 ((N - 1) / (D + 1)),
        ∑ d ∈ Icc (D + 1) ((N - 1) / c), f (c * d) d := by
  rw [sum_sigma', sum_sigma']
  change (∑ x ∈ divisorPairs N D, f x.1 x.2) =
    ∑ x ∈ productPairs N D, f (x.1 * x.2) x.2
  refine sum_bij (fun x _ => ⟨x.1 / x.2, x.2⟩) ?_ ?_ ?_ ?_
  · intro x hx
    obtain ⟨hn, hnN, hdvd, hD⟩ := mem_divisorPairs.mp hx
    apply (product_mem_iff hN).mpr
    have hd : 0 < x.2 := by omega
    have hc : 0 < x.1 / x.2 := Nat.div_pos (Nat.le_of_dvd hn hdvd) hd
    exact ⟨hc, hD, by simpa [Nat.div_mul_cancel hdvd] using hnN⟩
  · intro a ha b hb hab
    have da := (mem_divisorPairs.mp ha).2.2.1
    have db := (mem_divisorPairs.mp hb).2.2.1
    have hq : a.1 / a.2 = b.1 / b.2 := congrArg Sigma.fst hab
    have hd : a.2 = b.2 := congrArg Sigma.snd hab
    have hn : a.1 = b.1 := by
      calc
        a.1 = (a.1 / a.2) * a.2 := (Nat.div_mul_cancel da).symm
        _ = (b.1 / b.2) * b.2 := by rw [hq, hd]
        _ = b.1 := Nat.div_mul_cancel db
    exact Sigma.ext hn (by simpa using hd)
  · intro b hb
    obtain ⟨hc, hD, hmul⟩ := (product_mem_iff hN).mp hb
    have hd : 0 < b.2 := by omega
    have hmem : (⟨b.1 * b.2, b.2⟩ : Σ _ : ℕ, ℕ) ∈ divisorPairs N D := by
      apply mem_divisorPairs.mpr
      exact ⟨Nat.mul_pos hc hd, hmul, dvd_mul_left _ _, hD⟩
    refine ⟨⟨b.1 * b.2, b.2⟩, hmem, ?_⟩
    simp [Nat.mul_div_cancel _ hd]
  · intro x hx
    have hdvd := (mem_divisorPairs.mp hx).2.2.1
    simp [Nat.div_mul_cancel hdvd]

#print axioms reindex_divisor_tail

/-- An overcounting cover of proper prime powers at most N. -/
def properPowerCover (N : ℕ) : Finset ℕ :=
  ((Icc 2 (Nat.sqrt N)).product (Icc 2 (Nat.log 2 N))).image
    (fun pk => pk.1 ^ pk.2)

theorem properPowerCover_card_le (N : ℕ) :
    (properPowerCover N).card ≤ Nat.sqrt N * Nat.log 2 N := by
  calc
    (properPowerCover N).card ≤
        ((Icc 2 (Nat.sqrt N)).product (Icc 2 (Nat.log 2 N))).card := card_image_le
    _ = (Nat.sqrt N + 1 - 2) * (Nat.log 2 N + 1 - 2) := by simp
    _ ≤ Nat.sqrt N * Nat.log 2 N := by
      apply Nat.mul_le_mul <;> omega

theorem proper_prime_power_mem {m N : ℕ}
    (hm : ArithmeticFunction.vonMangoldt m ≠ 0) (hnp : ¬m.Prime) (hmN : m ≤ N) :
    m ∈ properPowerCover N := by
  obtain ⟨p, k, hp, hk, hpk⟩ :=
    (isPrimePow_nat_iff _).mp (ArithmeticFunction.vonMangoldt_ne_zero_iff.mp hm)
  have hk2 : 2 ≤ k := by
    by_contra h
    have hk1 : k = 1 := by omega
    rw [hk1, pow_one] at hpk
    exact hnp (hpk ▸ hp)
  have hpmax : p ≤ Nat.sqrt N := by
    apply Nat.le_sqrt'.mpr
    exact le_trans (Nat.pow_le_pow_right hp.pos hk2) (hpk ▸ hmN)
  have hkmax : k ≤ Nat.log 2 N := by
    apply Nat.le_log_of_pow_le (by decide : 1 < 2)
    exact le_trans (Nat.pow_le_pow_left hp.two_le k) (hpk ▸ hmN)
  apply mem_image.mpr
  exact ⟨(p, k), by simp [hp.two_le, hpmax, hk2, hkmax], hpk⟩

abbrev weight (N n : ℕ) : ℝ :=
  ArithmeticFunction.vonMangoldt n * ArithmeticFunction.vonMangoldt (N - n)

def convolution (N : ℕ) : ℝ := ∑ n ∈ Ioo 0 N, weight N n

/-- Prime-only contribution, using equivalent Mangoldt weights. -/
def primeConvolution (N : ℕ) : ℝ :=
  ∑ n ∈ Ioo 0 N, if n.Prime ∧ (N - n).Prime then weight N n else 0

def badPowerIndices (N : ℕ) : Finset ℕ :=
  (Ioo 0 N).filter (fun n => ¬(n.Prime ∧ (N - n).Prime) ∧ weight N n ≠ 0)

theorem weight_nonneg (N n : ℕ) : 0 ≤ weight N n :=
  mul_nonneg ArithmeticFunction.vonMangoldt_nonneg ArithmeticFunction.vonMangoldt_nonneg

theorem weight_le_log_sq {N n : ℕ} (hN : 1 ≤ N) (hn : n ∈ Ioo 0 N) :
    weight N n ≤ (Real.log N) ^ 2 := by
  have hn' := mem_Ioo.mp hn
  have hother : 0 < N - n := Nat.sub_pos_of_lt hn'.2
  have log_bound : ∀ m : ℕ, 0 < m → m ≤ N →
      ArithmeticFunction.vonMangoldt m ≤ Real.log N := by
    intro m hm hmN
    exact le_trans ArithmeticFunction.vonMangoldt_le_log
      (Real.log_le_log (by exact_mod_cast hm) (by exact_mod_cast hmN))
  rw [pow_two]
  exact mul_le_mul (log_bound n hn'.1 (le_of_lt hn'.2))
    (log_bound (N - n) hother (Nat.sub_le _ _))
    ArithmeticFunction.vonMangoldt_nonneg
    (Real.log_nonneg (by exact_mod_cast hN))

theorem badPowerIndices_card_le (N : ℕ) :
    (badPowerIndices N).card ≤ 2 * (properPowerCover N).card := by
  have hsub : badPowerIndices N ⊆
      properPowerCover N ∪ (properPowerCover N).image (fun m => N - m) := by
    intro n hn
    obtain ⟨hnI, hnp, hw⟩ := mem_filter.mp hn
    have hm1 : ArithmeticFunction.vonMangoldt n ≠ 0 := (mul_ne_zero_iff.mp hw).1
    have hm2 : ArithmeticFunction.vonMangoldt (N - n) ≠ 0 := (mul_ne_zero_iff.mp hw).2
    by_cases hnprime : n.Prime
    · apply mem_union_right
      have hother : ¬(N - n).Prime := fun hp => hnp ⟨hnprime, hp⟩
      apply mem_image.mpr
      exact ⟨N - n, proper_prime_power_mem hm2 hother (Nat.sub_le _ _),
        Nat.sub_sub_self (le_of_lt (mem_Ioo.mp hnI).2)⟩
    · exact mem_union_left _ (proper_prime_power_mem hm1 hnprime
        (le_of_lt (mem_Ioo.mp hnI).2))
  calc
    (badPowerIndices N).card ≤
        (properPowerCover N ∪ (properPowerCover N).image (fun m => N - m)).card :=
      card_le_card hsub
    _ ≤ (properPowerCover N).card +
        ((properPowerCover N).image (fun m => N - m)).card := card_union_le _ _
    _ ≤ 2 * (properPowerCover N).card := by
      have h := card_image_le (s := properPowerCover N) (f := fun m => N - m)
      omega

theorem convolution_sub_prime_eq (N : ℕ) :
    convolution N - primeConvolution N = ∑ n ∈ badPowerIndices N, weight N n := by
  unfold convolution primeConvolution badPowerIndices
  rw [← sum_sub_distrib, sum_filter]
  apply sum_congr rfl
  intro n _
  by_cases hp : n.Prime ∧ (N - n).Prime <;>
    by_cases hw : weight N n = 0 <;> simp [hp, hw]

/-- The elementary explicit estimate (04, Eq. 4.12). [T-Lean]
Nat.log 2 N is the floor of log base 2 N. No asymptotic or analytic hypothesis.
-/
theorem proper_prime_power_error (N : ℕ) (hN : 1 ≤ N) :
    0 ≤ convolution N - primeConvolution N ∧
    convolution N - primeConvolution N ≤
      2 * (Nat.sqrt N : ℝ) * (Nat.log 2 N : ℝ) * (Real.log N) ^ 2 := by
  rw [convolution_sub_prime_eq]
  constructor
  · exact sum_nonneg (fun n _ => weight_nonneg N n)
  · calc
      (∑ n ∈ badPowerIndices N, weight N n) ≤
          ∑ _n ∈ badPowerIndices N, (Real.log N) ^ 2 := by
        apply sum_le_sum
        intro n hn
        exact weight_le_log_sq hN (mem_filter.mp hn).1
      _ = ((badPowerIndices N).card : ℝ) * (Real.log N) ^ 2 := by simp
      _ ≤ ((2 * (Nat.sqrt N * Nat.log 2 N) : ℕ) : ℝ) * (Real.log N) ^ 2 := by
        apply mul_le_mul_of_nonneg_right _ (sq_nonneg _)
        exact_mod_cast le_trans (badPowerIndices_card_le N)
          (Nat.mul_le_mul_left 2 (properPowerCover_card_le N))
      _ = _ := by push_cast; ring

#print axioms proper_prime_power_error

theorem primeConvolution_eq_log (N : ℕ) :
    primeConvolution N = ∑ n ∈ Ioo 0 N,
      if n.Prime ∧ (N - n).Prime then Real.log n * Real.log ((N - n : ℕ) : ℝ) else 0 := by
  apply sum_congr rfl
  intro n _
  by_cases hp : n.Prime ∧ (N - n).Prime
  · simp [hp, weight, ArithmeticFunction.vonMangoldt_apply_prime hp.1,
      ArithmeticFunction.vonMangoldt_apply_prime hp.2]
  · simp [hp]

abbrev mobiusTerm (N n d : ℕ) : ℝ :=
  (ArithmeticFunction.moebius d : ℝ) * Real.log d *
    ArithmeticFunction.vonMangoldt (N - n)

theorem reindex_mobius_tail (N D : ℕ) (hN : 0 < N) :
    (∑ n ∈ Ioo 0 N, ∑ d ∈ n.divisors.filter (fun d => D < d), mobiusTerm N n d) =
      ∑ c ∈ Icc 1 ((N - 1) / (D + 1)), ∑ d ∈ Icc (D + 1) ((N - 1) / c),
        (ArithmeticFunction.moebius d : ℝ) * Real.log d *
          ArithmeticFunction.vonMangoldt (N - c * d) :=
  reindex_divisor_tail N D hN (mobiusTerm N)

/-- A local copy of the coprimality argument from Reduccion.lean, to avoid a build
dependency on its .olean. Its source is preserved. -/
theorem coprime_coefficient {N c d : ℕ} (hc : Nat.Coprime (N - c * d) N) :
    Nat.Coprime c N := by
  have h1 : Nat.gcd c N ∣ N - c * d :=
    Nat.dvd_sub (Nat.gcd_dvd_right c N) (Dvd.dvd.mul_right (Nat.gcd_dvd_left c N) d)
  have h2 : Nat.gcd c N ∣ Nat.gcd (N - c * d) N :=
    Nat.dvd_gcd h1 (Nat.gcd_dvd_right c N)
  rw [Nat.Coprime.gcd_eq_one hc] at h2
  exact Nat.dvd_one.mp h2

def badIntegers (N : ℕ) : Finset ℕ :=
  (Ioo 0 N).filter (fun m => ¬Nat.Coprime m N ∧ ArithmeticFunction.vonMangoldt m ≠ 0)

/-- Divisor form of the bad tail, with (c,N)>1 and c=n/d. -/
def badTail (N D : ℕ) : ℝ :=
  ∑ n ∈ Ioo 0 N, ∑ d ∈ n.divisors.filter (fun d => D < d),
    if ¬Nat.Coprime (n / d) N then mobiusTerm N n d else 0

theorem badTail_eq_product (N D : ℕ) (hN : 0 < N) :
    badTail N D =
      ∑ c ∈ Icc 1 ((N - 1) / (D + 1)), ∑ d ∈ Icc (D + 1) ((N - 1) / c),
        if ¬Nat.Coprime c N then mobiusTerm N (c * d) d else 0 := by
  rw [badTail, reindex_divisor_tail N D hN]
  apply sum_congr rfl
  intro c _
  apply sum_congr rfl
  intro d hd
  have hdpos : 0 < d := by have := (mem_Icc.mp hd).1; omega
  simp [Nat.mul_div_cancel _ hdpos]

theorem sum_reflection (N : ℕ) (f : ℕ → ℝ) :
    (∑ n ∈ Ioo 0 N, f (N - n)) = ∑ m ∈ Ioo 0 N, f m := by
  refine sum_bij (fun n _ => N - n) ?_ ?_ ?_ ?_
  · intro n hn
    simp only [mem_Ioo] at hn ⊢
    omega
  · intro a ha b hb hab
    have := mem_Ioo.mp ha
    have := mem_Ioo.mp hb
    omega
  · intro m hm
    refine ⟨N - m, ?_, ?_⟩
    · simp only [mem_Ioo] at hm ⊢
      omega
    · exact Nat.sub_sub_self (le_of_lt (mem_Ioo.mp hm).2)
  · intro _ _
    rfl

theorem abs_mobiusTerm_le {N n d : ℕ} (_hN : 1 ≤ N) (hn : n ∈ Ioo 0 N)
    (hd : d ∈ n.divisors) :
    |mobiusTerm N n d| ≤ Real.log N * ArithmeticFunction.vonMangoldt (N - n) := by
  have hdpos : 0 < d := Nat.pos_of_mem_divisors hd
  have hdN : d ≤ N := le_trans (Nat.divisor_le hd) (le_of_lt (mem_Ioo.mp hn).2)
  have hlog0 : 0 ≤ Real.log d := Real.log_nonneg (by exact_mod_cast hdpos)
  have hlog : Real.log d ≤ Real.log N :=
    Real.log_le_log (by exact_mod_cast hdpos) (by exact_mod_cast hdN)
  have hmu : |(ArithmeticFunction.moebius d : ℝ)| ≤ 1 := by
    exact_mod_cast (ArithmeticFunction.abs_moebius_le_one (n := d))
  simp only [mobiusTerm, abs_mul, abs_of_nonneg hlog0,
    abs_of_nonneg ArithmeticFunction.vonMangoldt_nonneg]
  calc
    _ ≤ 1 * Real.log d * ArithmeticFunction.vonMangoldt (N - n) :=
      mul_le_mul_of_nonneg_right (mul_le_mul_of_nonneg_right hmu hlog0)
        ArithmeticFunction.vonMangoldt_nonneg
    _ ≤ Real.log N * ArithmeticFunction.vonMangoldt (N - n) := by
      simpa using mul_le_mul_of_nonneg_right hlog ArithmeticFunction.vonMangoldt_nonneg

theorem bad_term_bound {N n d : ℕ} (hN : 1 ≤ N) (hn : n ∈ Ioo 0 N)
    (hd : d ∈ n.divisors) :
    |if ¬Nat.Coprime (n / d) N then mobiusTerm N n d else 0| ≤
      if N - n ∈ badIntegers N then
        Real.log N * ArithmeticFunction.vonMangoldt (N - n) else 0 := by
  by_cases hc : Nat.Coprime (n / d) N
  · have hz : (if ¬Nat.Coprime (n / d) N then mobiusTerm N n d else 0) = 0 :=
      ite_eq_right (not_not.mpr hc)
    rw [hz, abs_zero]
    split_ifs
    · exact mul_nonneg (Real.log_nonneg (by exact_mod_cast hN))
        ArithmeticFunction.vonMangoldt_nonneg
    · exact le_rfl
  · simp only [hc, not_false_eq_true, ite_true]
    by_cases hm : ArithmeticFunction.vonMangoldt (N - n) = 0
    · simp [mobiusTerm, hm]
    · have hnot : ¬Nat.Coprime (N - n) N := by
        intro h
        have h' : Nat.Coprime (N - (n / d) * d) N := by
          simpa [Nat.div_mul_cancel (Nat.mem_divisors.mp hd).1] using h
        exact hc (coprime_coefficient h')
      have hmem : N - n ∈ badIntegers N := by
        apply mem_filter.mpr
        exact ⟨mem_Ioo.mpr ⟨Nat.sub_pos_of_lt (mem_Ioo.mp hn).2, by
          have := (mem_Ioo.mp hn).1; omega⟩, hnot, hm⟩
      simp only [hmem, ite_true]
      exact abs_mobiusTerm_le hN hn hd

theorem badTail_le_mass (N D : ℕ) (hN : 1 ≤ N) :
    |badTail N D| ≤ Real.log N *
      ∑ m ∈ badIntegers N, ArithmeticFunction.vonMangoldt m * ((N - m).divisors.card : ℝ) := by
  let f : ℕ → ℝ := fun m => if m ∈ badIntegers N then
    Real.log N * (ArithmeticFunction.vonMangoldt m * ((N - m).divisors.card : ℝ)) else 0
  have hinner : ∀ n ∈ Ioo 0 N,
      |∑ d ∈ n.divisors.filter (fun d => D < d),
          if ¬Nat.Coprime (n / d) N then mobiusTerm N n d else 0| ≤ f (N - n) := by
    intro n hn
    let b : ℝ := if N - n ∈ badIntegers N then
      Real.log N * ArithmeticFunction.vonMangoldt (N - n) else 0
    have hb : 0 ≤ b := by
      dsimp [b]
      split_ifs
      · exact mul_nonneg (Real.log_nonneg (by exact_mod_cast hN))
          ArithmeticFunction.vonMangoldt_nonneg
      · exact le_rfl
    calc
      _ ≤ ∑ d ∈ n.divisors.filter (fun d => D < d),
          |if ¬Nat.Coprime (n / d) N then mobiusTerm N n d else 0| := abs_sum_le_sum_abs _ _
      _ ≤ ∑ _d ∈ n.divisors.filter (fun d => D < d), b := by
        apply sum_le_sum
        intro d hd
        exact bad_term_bound hN hn (mem_filter.mp hd).1
      _ ≤ ∑ _d ∈ n.divisors, b :=
        sum_le_sum_of_subset_of_nonneg (filter_subset _ _) (fun _ _ _ => hb)
      _ = f (N - n) := by
        have hnn : N - (N - n) = n := Nat.sub_sub_self (le_of_lt (mem_Ioo.mp hn).2)
        simp only [sum_const, nsmul_eq_mul, f, b, hnn]
        split_ifs <;> ring
  calc
    |badTail N D| ≤ ∑ n ∈ Ioo 0 N,
        |∑ d ∈ n.divisors.filter (fun d => D < d),
          if ¬Nat.Coprime (n / d) N then mobiusTerm N n d else 0| := abs_sum_le_sum_abs _ _
    _ ≤ ∑ n ∈ Ioo 0 N, f (N - n) := sum_le_sum hinner
    _ = ∑ m ∈ Ioo 0 N, f m := sum_reflection N f
    _ = Real.log N * ∑ m ∈ badIntegers N,
        ArithmeticFunction.vonMangoldt m * ((N - m).divisors.card : ℝ) := by
      rw [mul_sum]
      change (∑ m ∈ Ioo 0 N, if m ∈ badIntegers N then _ else 0) = _
      rw [← sum_filter]
      congr 1
      ext m
      simp only [mem_filter, badIntegers]
      tauto

/-- Finite indexing of exactly the p^k<N with p prime dividing N. -/
def badPrimePowers (N : ℕ) : Finset (ℕ × ℕ) :=
  (((Icc 2 N).filter (fun p => p.Prime ∧ p ∣ N)).product (Icc 1 (Nat.log 2 N))).filter
    (fun pk => pk.1 ^ pk.2 < N)

/-- The explicit K_N of Proposition 4.2, with finite exponent cutoff. -/
def badPrimePowerBound (N : ℕ) : ℝ :=
  Real.log N * ∑ pk ∈ badPrimePowers N,
    Real.log pk.1 * ((N - pk.1 ^ pk.2).divisors.card : ℝ)

theorem badIntegers_subset_powers (N : ℕ) (hN : 0 < N) :
    badIntegers N ⊆ (badPrimePowers N).image (fun pk => pk.1 ^ pk.2) := by
  intro m hm
  obtain ⟨hmI, hnot, hm0⟩ := mem_filter.mp hm
  obtain ⟨p, k, hp, hk, hpk⟩ :=
    (isPrimePow_nat_iff _).mp (ArithmeticFunction.vonMangoldt_ne_zero_iff.mp hm0)
  have hpN : p ∣ N := by
    by_contra h
    exact hnot (hpk ▸ Nat.Coprime.pow_left k ((Nat.Prime.coprime_iff_not_dvd hp).mpr h))
  have hkmax : k ≤ Nat.log 2 N := by
    apply Nat.le_log_of_pow_le (by decide : 1 < 2)
    exact le_trans (Nat.pow_le_pow_left hp.two_le k) (hpk ▸ le_of_lt (mem_Ioo.mp hmI).2)
  apply mem_image.mpr
  refine ⟨(p, k), ?_, hpk⟩
  change (p, k) ∈ (((Icc 2 N).filter (fun p => p.Prime ∧ p ∣ N)).product
    (Icc 1 (Nat.log 2 N))).filter (fun pk => pk.1 ^ pk.2 < N)
  apply Finset.mem_filter.mpr
  refine ⟨Finset.mem_product.mpr ⟨?_, ?_⟩, hpk ▸ (mem_Ioo.mp hmI).2⟩
  · exact Finset.mem_filter.mpr ⟨mem_Icc.mpr ⟨hp.two_le, Nat.le_of_dvd hN hpN⟩, hp, hpN⟩
  · exact mem_Icc.mpr ⟨by omega, hkmax⟩

/-- Proposition 4.2: the full bad tail, every c and d, is bounded by K_N. [T-Lean]
This is the explicit finite inequality, not the asymptotic divisor bound for K_N.
-/
theorem bad_tail_bound (N D : ℕ) (hN : 1 ≤ N) :
    |badTail N D| ≤ badPrimePowerBound N := by
  apply le_trans (badTail_le_mass N D hN)
  unfold badPrimePowerBound
  apply mul_le_mul_of_nonneg_left _ (Real.log_nonneg (by exact_mod_cast hN))
  let g : ℕ → ℝ := fun m =>
    ArithmeticFunction.vonMangoldt m * ((N - m).divisors.card : ℝ)
  have hg : ∀ m, 0 ≤ g m := fun _ =>
    mul_nonneg ArithmeticFunction.vonMangoldt_nonneg (Nat.cast_nonneg _)
  calc
    (∑ m ∈ badIntegers N, g m) ≤
        ∑ m ∈ (badPrimePowers N).image (fun pk => pk.1 ^ pk.2), g m :=
      sum_le_sum_of_subset_of_nonneg (badIntegers_subset_powers N (by omega))
        (fun m _ _ => hg m)
    _ ≤ ∑ pk ∈ badPrimePowers N, g (pk.1 ^ pk.2) :=
      sum_image_le_of_nonneg (fun m _ => hg m)
    _ = ∑ pk ∈ badPrimePowers N,
        Real.log pk.1 * ((N - pk.1 ^ pk.2).divisors.card : ℝ) := by
      apply sum_congr rfl
      intro pk hpk
      have hpk' : pk ∈ (((Icc 2 N).filter (fun p => p.Prime ∧ p ∣ N)).product
          (Icc 1 (Nat.log 2 N))).filter (fun pk => pk.1 ^ pk.2 < N) := hpk
      obtain ⟨hp, hk⟩ := mem_product.mp (mem_filter.mp hpk').1
      have hpprime := (mem_filter.mp hp).2.1
      have hkpos : 0 < pk.2 := (mem_Icc.mp hk).1
      simp only [g, ArithmeticFunction.vonMangoldt_apply_pow (by omega : pk.2 ≠ 0),
        ArithmeticFunction.vonMangoldt_apply_prime hpprime]

#print axioms reindex_mobius_tail
#print axioms primeConvolution_eq_log
#print axioms badTail_eq_product
#print axioms bad_tail_bound

end GoldbachLumina.Reduction2
