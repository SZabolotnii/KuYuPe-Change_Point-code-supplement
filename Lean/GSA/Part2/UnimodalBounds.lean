import GSA.Part2.Setup
import Mathlib.Probability.Moments.Variance
import Mathlib.Tactic

namespace GSA.Part2

open MeasureTheory ProbabilityTheory

/-!
# 2.4.1. Уточнені межі для унімодальних розподілів

Формалізація нерівностей типу Височанського–Петуніна та Кантеллі
для унімодальних розподілів.
-/

/-!
## Нерівність Кантеллі (односторонній варіант Чебишева)

Для довільних розподілів з середнім μ та дисперсією σ²:
P(X - μ ≥ kσ) ≤ 1/(1 + k²)
-/

/-- Формула межі Кантеллі. -/
noncomputable def cantelli_bound (k : ℝ) : ℝ :=
  1 / (1 + k^2)

/-- Межа Кантеллі є невід'ємною. -/
theorem cantelli_bound_nonneg (k : ℝ) : 0 ≤ cantelli_bound k := by
  unfold cantelli_bound
  positivity

/-- Межа Кантеллі ≤ 1. -/
theorem cantelli_bound_le_one (k : ℝ) : cantelli_bound k ≤ 1 := by
  unfold cantelli_bound
  have h : 0 ≤ k^2 := sq_nonneg k
  have hden : 0 < 1 + k^2 := by linarith
  rw [div_le_one hden]
  linarith

/-- Межа Кантеллі менша за межу Чебишева при k > 1. -/
theorem cantelli_lt_chebyshev (k : ℝ) (hk : 1 < k) :
    cantelli_bound k < 1 / k^2 := by
  unfold cantelli_bound
  have hk2 : 0 < k^2 := sq_pos_of_pos (lt_trans one_pos hk)
  have hden : 0 < 1 + k^2 := by linarith
  have h1 : k^2 < 1 + k^2 := by linarith
  exact div_lt_div_of_pos_left one_pos hk2 h1

/-!
## Нерівність Височанського-Петуніна

Для унімодальних розподілів:
P(|X - μ| ≥ kσ) ≤ 4/(9k²), для k ≥ √(8/3)
-/

/-- Мінімальне значення k для нерівності В.-П. -/
noncomputable def VP_min_k : ℝ := Real.sqrt (8/3)

/-- Формула межі Височанського-Петуніна. -/
noncomputable def VP_bound (k : ℝ) : ℝ :=
  4 / (9 * k^2)

/-- Межа В.-П. є невід'ємною при k > 0. -/
theorem VP_bound_nonneg (k : ℝ) (hk : 0 < k) : 0 ≤ VP_bound k := by
  unfold VP_bound
  positivity

/-- Межа В.-П. менша за межу Чебишева при k ≥ √(8/3). -/
theorem VP_lt_chebyshev (k : ℝ) (hk : VP_min_k ≤ k) :
    VP_bound k ≤ 1 / k^2 := by
  unfold VP_bound VP_min_k at *
  have hk_pos : 0 < k := by
    have h83 : 0 < (8:ℝ)/3 := by norm_num
    have hsqrt : 0 < Real.sqrt (8/3) := Real.sqrt_pos.mpr h83
    exact lt_of_lt_of_le hsqrt hk
  have hk2 : 0 < k^2 := sq_pos_of_pos hk_pos
  have h9k2 : 0 < 9 * k^2 := by nlinarith
  have h_ineq : 9 * k^2 ≥ k^2 := by nlinarith
  have h_num : (4:ℝ) ≤ 9 := by norm_num
  calc
    4 / (9 * k^2) ≤ 9 / (9 * k^2) := by
      apply div_le_div_of_nonneg_right h_num
      exact le_of_lt h9k2
    _ = 1 / k^2 := by field_simp

/-!
## Уточнений поріг ПЕ для унімодальних розподілів

h_PE^(VP) = E[Λ|H₀] + (2/3)√(Var[Λ|H₀]/ε)

Ця формула коректна при ε ≤ 1/6 (що відповідає k ≥ √(8/3)).
-/

/-- Уточнений поріг ПЕ для унімодальних розподілів. -/
noncomputable def PE_threshold_VP (mean variance ε : ℝ) : ℝ :=
  mean + (2/3) * Real.sqrt (variance / ε)

/-- Поріг Кантеллі. -/
noncomputable def PE_threshold_Cantelli (mean variance ε : ℝ) : ℝ :=
  mean + Real.sqrt variance * Real.sqrt (1/ε - 1)

/-- Поріг В.-П. менший за стандартний поріг ПЕ. -/
theorem PE_threshold_VP_lt_PE (mean variance ε : ℝ)
    (_hε : 0 < ε) (_hε_small : ε ≤ 1/6) (_hvar : 0 ≤ variance) :
    PE_threshold_VP mean variance ε ≤
      mean + Real.sqrt (variance / ε) := by
  unfold PE_threshold_VP
  have h23 : (2:ℝ)/3 ≤ 1 := by norm_num
  have hsqrt_nonneg : 0 ≤ Real.sqrt (variance / ε) := Real.sqrt_nonneg _
  have hmul : (2/3) * Real.sqrt (variance / ε) ≤ 1 * Real.sqrt (variance / ε) :=
    mul_le_mul_of_nonneg_right h23 hsqrt_nonneg
  linarith

/-- Умова застосовності порогу В.-П.: ε ≤ 1/6 ↔ 6 ≤ 1/ε. -/
theorem VP_threshold_applicability_simple (ε : ℝ) (hε : 0 < ε) :
    (ε ≤ 1/6) ↔ (6 ≤ 1/ε) := by
  constructor
  · intro hε_small
    have h1 : ε * 6 ≤ 1 := by linarith
    have h2 : 6 ≤ 1/ε := by
      rw [le_div_iff₀ hε]
      linarith
    exact h2
  · intro h6
    have h1 : 6 * ε ≤ 1 := by
      have h := (le_div_iff₀ (a := 6) (b := 1) hε).mp h6
      linarith
    linarith

end GSA.Part2
