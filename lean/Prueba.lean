/-
Prueba de humo de la instalación Lean 4 + Mathlib para el proyecto Goldbach-LUMINA.
Procedencia: Claude Code (claude-opus-5-5), 2026-11-23.
-/
import Mathlib

open ArithmeticFunction
open scoped ArithmeticFunction.Moebius ArithmeticFunction.vonMangoldt

/-- Conjetura fuerte de Goldbach (solo el enunciado; NO está demostrada). -/
def GoldbachFuerte : Prop :=
  ∀ n : ℕ, 4 ≤ n → Even n → ∃ p q : ℕ, p.Prime ∧ q.Prime ∧ p + q = n

/-- Un caso concreto, comprobado por el kernel. -/
example : ∃ p q : ℕ, p.Prime ∧ q.Prime ∧ p + q = 100 :=
  ⟨3, 97, by norm_num, by norm_num, by norm_num⟩

/-- Identidad base de la doc 04 (Prop. 4.1): Λ(n) = −Σ_{d|n} μ(d) log d. -/
example (n : ℕ) :
    ∑ d ∈ n.divisors, (ArithmeticFunction.moebius d : ℝ) * Real.log d
      = -ArithmeticFunction.vonMangoldt n :=
  ArithmeticFunction.sum_moebius_mul_log_eq
