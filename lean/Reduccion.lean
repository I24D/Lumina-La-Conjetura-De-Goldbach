/-
Formalización de las identidades base del proyecto Goldbach-LUMINA
(matematicas/04, Prop. 4.1 y núcleo de la Prop. 4.2; matematicas/06, Lemas 6.0 y 6.1, paso 1).

Procedencia: Claude Code (claude-opus-5-5), 2026-11-23.
Estado: lo que compile sin `sorry` queda verificado por el kernel de Lean [T-Lean].
Nada de esto demuestra Goldbach: son las identidades exactas sobre las que se apoya la reducción.
-/
import Mathlib

open Finset

noncomputable section

namespace GoldbachLumina

/-- Abreviaturas sin notación ambigua. -/
abbrev mu (n : ℕ) : ℝ := (ArithmeticFunction.moebius n : ℝ)
abbrev vM (n : ℕ) : ℝ := ArithmeticFunction.vonMangoldt n

/-- R(N) = Σ_{0<n<N} Λ(n) Λ(N-n). -/
def R (N : ℕ) : ℝ := ∑ n ∈ Ioo 0 N, vM n * vM (N - n)

/-- Sumando común: μ(d) log d Λ(N-n). -/
def term (N n d : ℕ) : ℝ := mu d * Real.log d * vM (N - n)

/-- Λ(n) = −Σ_{d|n} μ(d) log d (Mathlib), con `Real.log`. -/
theorem sum_mu_log (n : ℕ) : ∑ d ∈ n.divisors, mu d * Real.log d = -vM n :=
  ArithmeticFunction.sum_moebius_mul_log_eq

/-- Identidad de partida (04, §2): R(N) = −Σ_{0<n<N} Σ_{d|n} μ(d) log d Λ(N−n). -/
theorem R_eq_neg_sum (N : ℕ) :
    R N = -∑ n ∈ Ioo 0 N, ∑ d ∈ n.divisors, term N n d := by
  unfold R term
  rw [← Finset.sum_neg_distrib]
  refine Finset.sum_congr rfl (fun n _ => ?_)
  rw [← Finset.sum_mul, sum_mu_log]
  ring

/-- Lado S (divisores pequeños d ≤ D), en forma de divisores. -/
def S (N D : ℕ) : ℝ := -∑ n ∈ Ioo 0 N, ∑ d ∈ n.divisors.filter (· ≤ D), term N n d

/-- Cola T (divisores grandes d > D), en forma de divisores. -/
def T (N D : ℕ) : ℝ := ∑ n ∈ Ioo 0 N, ∑ d ∈ n.divisors.filter (fun d => ¬ d ≤ D), term N n d

/-- **Prop. 4.1 de 04** (forma de divisores): R(N) = S_N(D) − T_N(D), para todo N y todo D. -/
theorem R_eq_S_sub_T (N D : ℕ) : R N = S N D - T N D := by
  rw [R_eq_neg_sum]
  unfold S T
  have h : ∀ n : ℕ, ∑ d ∈ n.divisors, term N n d =
      ∑ d ∈ n.divisors.filter (· ≤ D), term N n d +
      ∑ d ∈ n.divisors.filter (fun d => ¬ d ≤ D), term N n d :=
    fun n => (Finset.sum_filter_add_sum_filter_not _ _ _).symm
  simp_rw [h, Finset.sum_add_distrib]
  ring

/-- **06, Lema 6.1, paso 1:** si N − cd es coprimo con N, entonces c es coprimo con N.
Así la restricción (c,N)=1 de T^prim es automática salvo cuando (N−cd, N) > 1. -/
theorem coprime_of_coprime_sub {N c d : ℕ} (hc : Nat.Coprime (N - c * d) N) :
    Nat.Coprime c N := by
  have h1 : Nat.gcd c N ∣ N - c * d :=
    Nat.dvd_sub (Nat.gcd_dvd_right c N) (Dvd.dvd.mul_right (Nat.gcd_dvd_left c N) d)
  have h2 : Nat.gcd c N ∣ Nat.gcd (N - c * d) N := Nat.dvd_gcd h1 (Nat.gcd_dvd_right c N)
  rw [Nat.Coprime.gcd_eq_one hc] at h2
  exact Nat.dvd_one.mp h2

/-- **Núcleo de la Prop. 4.2 de 04:** si Λ(m) ≠ 0 y m no es coprimo con N, entonces m es
potencia de un primo que divide a N. (Esos son los únicos términos "malos".) -/
theorem primePow_of_not_coprime {m N : ℕ} (hm : vM m ≠ 0) (h : ¬ Nat.Coprime m N) :
    ∃ p k : ℕ, p.Prime ∧ 0 < k ∧ p ^ k = m ∧ p ∣ N := by
  rw [ArithmeticFunction.vonMangoldt_ne_zero_iff] at hm
  obtain ⟨p, k, hp, hk, rfl⟩ := (isPrimePow_nat_iff _).mp hm
  refine ⟨p, k, hp, hk, rfl, ?_⟩
  by_contra hpN
  exact h (Nat.Coprime.pow_left k ((Nat.Prime.coprime_iff_not_dvd hp).mpr hpN))

/-- **06, Lema 6.0 (paso clave):** para p primo, μ(p) log p = −log p. Por eso, en la caja c = 1,
la parte con d primo es exactamente −Σ log p · Λ(N−p), un bloque de la suma de Goldbach. -/
theorem mu_log_prime {p : ℕ} (hp : p.Prime) : mu p * Real.log p = -Real.log p := by
  unfold mu
  rw [ArithmeticFunction.moebius_apply_prime hp]
  push_cast
  ring

/-- **06, R-3:** un primo p con D < p contribuye a la cola T únicamente a través de d = p
(el divisor d = 1 está en el lado S y además log 1 = 0). Formalmente: los divisores de p
mayores que D son exactamente {p}. -/
theorem divisors_prime_gt {p D : ℕ} (hp : p.Prime) (hD : D < p) (hD1 : 1 ≤ D) :
    p.divisors.filter (fun d => ¬ d ≤ D) = {p} := by
  rw [Nat.Prime.divisors hp]
  ext x
  simp only [Finset.mem_filter, Finset.mem_insert, Finset.mem_singleton, not_le]
  constructor
  · rintro ⟨h | h, hx⟩
    · omega
    · exact h
  · rintro rfl
    exact ⟨Or.inr rfl, hD⟩

end GoldbachLumina
