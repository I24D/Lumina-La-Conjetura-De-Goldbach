import Mathlib.Data.Nat.Prime.Basic

/-!
# Non-reduced residue classes in the reflected Goldbach sequence

Finite divisibility facts used to isolate the `(d,N)>1` part of the switched
sequence.  These statements contain no analytic estimate.
-/

namespace GoldbachLumina

/-- If `g` divides both `N` and `N-t`, then it divides `t`. -/
theorem dvd_reflected_term {g N t : ℕ} (ht : t ≤ N)
    (hN : g ∣ N) (hNt : g ∣ N - t) : g ∣ t := by
  have h : g ∣ N - (N - t) := Nat.dvd_sub hN hNt
  simpa [Nat.sub_sub_self ht] using h

/-- A prime common to the modulus and `N` is forced into every reflected
term `t` lying in the residue class `N (mod d)`. -/
theorem common_prime_forced_into_term {p d N t : ℕ} (ht : t ≤ N)
    (hpd : p ∣ d) (hpN : p ∣ N) (hd : d ∣ N - t) : p ∣ t := by
  exact dvd_reflected_term ht hpN (dvd_trans hpd hd)

/-- Consequently, if every prime factor of an admissible switched term is at
least `z`, then a common prime coming from a non-reduced class is at least
`z` as well. -/
theorem common_prime_is_large {p d N t z : ℕ} (ht : t ≤ N)
    (hp : Nat.Prime p) (hpd : p ∣ d) (hpN : p ∣ N)
    (hd : d ∣ N - t) (hrough : ∀ q, Nat.Prime q → q ∣ t → z ≤ q) : z ≤ p := by
  exact hrough p hp (common_prime_forced_into_term ht hpd hpN hd)

#print axioms GoldbachLumina.dvd_reflected_term
#print axioms GoldbachLumina.common_prime_forced_into_term
#print axioms GoldbachLumina.common_prime_is_large

end GoldbachLumina
