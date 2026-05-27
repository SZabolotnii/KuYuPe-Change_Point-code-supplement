import GSA.Part2.BasisApprox

namespace GSA.Part2

open scoped BigOperators
open Filter

/-!
# 2.6. Теорема 2 (інформаційний функціонал J(s)) — каркас (HilbertBasis)

У `Part2_review.md`:
- J(s) ≤ J
- J(s) монотонно зростає
- за повного базису J(s) → J

Концептуально J(s) — "енергія" ортогональної проєкції LLR на вкладений підпростір.

Тут використовуємо канонічну реалізацію через `HilbertBasis`:
`J(s)` — часткова сума квадратів коефіцієнтів розкладу `z` по базису.
-/

variable {H : Type*} [NormedAddCommGroup H] [InnerProductSpace ℝ H]

/-- Information functional J(s): partial sum of squared Hilbert-basis coefficients of `z`.

  **Interpretation (paper § 2.7, C2 reformulation).**
  When `z` is the score vector (derivative of log-density) evaluated at a small perturbation δ,
  J(s) equals the s-term projection of the χ²-divergence onto the basis subspace, and — to
  leading order in δ — approximates twice the KL divergence via Fisher information in basis
  coordinates (Le Cam QMD / local asymptotic normality).  In the orthonormal case the equality
  J(s) = ‖Pₛ z‖² makes this a Parseval-type projection of χ², explaining why J(s) ≤ J = ‖z‖²
  (Theorem 2a) and J(s) → J as s → ∞ (Theorem 2c). -/
noncomputable def J (b : HilbertBasis ℕ ℝ H) (z : H) (s : ℕ) : ℝ :=
  (Finset.range s).sum fun i => (b.repr z i)^2

/-- Допоміжна рівність: повна сума квадратів коефіцієнтів дорівнює ‖z‖². -/
theorem tsum_repr_sq_eq_norm_sq (b : HilbertBasis ℕ ℝ H) (z : H) :
    (∑' i : ℕ, (b.repr z i)^2) = ‖z‖^2 := by
  have h := b.tsum_inner_mul_inner z z
  -- Приводимо терми до (b.repr z i)^2 і ⟪z,z⟫ = ‖z‖².
  have hterm :
      (fun i : ℕ => inner ℝ z (b i) * inner ℝ (b i) z) =
        (fun i : ℕ => (b.repr z i)^2) := by
    funext i
    -- b.repr z i = ⟪b i, z⟫, а для ℝ inner симетричний.
    simp [HilbertBasis.repr_apply_apply, real_inner_comm, pow_two]
  simpa [hterm, real_inner_self_eq_norm_sq] using h

/- Theorem 2 (a/b/c): the three structural properties of J(s) — upper bound by ‖z‖²,
   monotone growth with s, and convergence to ‖z‖² — are the Hilbert-space formalization of
   the local-Fisher / Parseval-projection picture described in paper § 2.7. -/
/-- (a) Обмеженість зверху: J(s) ≤ ‖z‖². -/
theorem theorem2_a_upper_bound
    (b : HilbertBasis ℕ ℝ H) (z : H) (s : ℕ) :
    J b z s ≤ ‖z‖^2 := by
  have hsum : Summable (fun i : ℕ => (b.repr z i)^2) := by
    -- Витікає з сумовності скалярних добутків для ортонормованої системи.
    have hsum' := (b.orthonormal.inner_products_summable (x := z))
    refine hsum'.congr ?_
    intro i
    -- ‖⟪b i, z⟫‖^2 = (b.repr z i)^2 для ℝ.
    simp [HilbertBasis.repr_apply_apply, Real.norm_eq_abs, pow_two]
  have hle :
      (Finset.range s).sum (fun i => (b.repr z i)^2) ≤
        ∑' i : ℕ, (b.repr z i)^2 := by
    refine Summable.sum_le_tsum (s := Finset.range s) ?_ hsum
    intro i hi
    exact sq_nonneg (b.repr z i)
  simpa [J, tsum_repr_sq_eq_norm_sq] using hle

/-- (b) Монотонність J(s). -/
theorem theorem2_b_monotone
    (b : HilbertBasis ℕ ℝ H) (z : H) :
    Monotone (J b z) := by
  classical
  intro s t hst
  have hsubset : Finset.range s ⊆ Finset.range t := by
    intro i hi
    exact Finset.mem_range.mpr (lt_of_lt_of_le (Finset.mem_range.mp hi) hst)
  refine Finset.sum_le_sum_of_subset_of_nonneg hsubset ?_
  intro i hi hnot
  exact sq_nonneg (b.repr z i)

/-- (c) Збіжність J(s) → ‖z‖². -/
theorem theorem2_c_tendsto
    (b : HilbertBasis ℕ ℝ H) (z : H) :
    Filter.Tendsto (J b z) Filter.atTop (nhds (‖z‖^2)) := by
  have hsum : Summable (fun i : ℕ => (b.repr z i)^2) := by
    have hsum' := (b.orthonormal.inner_products_summable (x := z))
    refine hsum'.congr ?_
    intro i
    simp [HilbertBasis.repr_apply_apply, Real.norm_eq_abs, pow_two]
  have ht :
      Tendsto (fun s : ℕ => (Finset.range s).sum fun i => (b.repr z i)^2)
        atTop (nhds (∑' i : ℕ, (b.repr z i)^2)) :=
    Summable.tendsto_sum_tsum_nat hsum
  simpa [J, tsum_repr_sq_eq_norm_sq] using ht

end GSA.Part2
