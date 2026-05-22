import GSA.Part2.Setup
import Mathlib.Probability.Moments.Variance

namespace GSA.Part2

open MeasureTheory ProbabilityTheory
open Filter

variable {α : Type*} [MeasurableSpace α]

/-! Moment-style formulation of an asymptotic normality assumption. -/
def AsymptoticallyNormal
    (stat : ℕ → α → ℝ) (μ : Measure α) (m s : ℝ) : Prop :=
  Tendsto (fun n => Expect μ (stat n)) atTop (nhds m) ∧
    Tendsto (fun n => ProbabilityTheory.variance (stat n) μ) atTop (nhds (s^2))

/-!
# 2.4. Робастне вирішальне правило (критерій ПЕ)

Тут у `Part2_review.md` є:
- 2.4.2 Теорема 3 (асимптотична еквівалентність ПЕ та Неймана–Пірсона);
- 2.4.3 Теорема 5 (критерій асимптотичної нормальності, "критерій Ю").

Нижче формалізовані ключові властивості цих теорем.
-/

/-- Теорема 3: Для асимптотично нормальних статистик різниця середніх збігається. -/
theorem theorem3_PE_asymptotic_equivalence_NP :
    ∀ {μ0 μ1 : Measure α} {stat : ℕ → α → ℝ} {m0 m1 s0 s1 : ℝ},
      AsymptoticallyNormal stat μ0 m0 s0 →
      AsymptoticallyNormal stat μ1 m1 s1 →
      Tendsto (fun n => Expect μ1 (stat n) - Expect μ0 (stat n)) atTop (nhds (m1 - m0)) := by
  intro μ0 μ1 stat m0 m1 s0 s1 h0 h1
  simpa using (h1.1.sub h0.1)

/-- Теорема 5: Для асимптотично нормальних статистик
    відхилення від середнього збігається до нуля. -/
theorem theorem5_asymptotic_normality_criterion_Y :
    ∀ {μ : Measure α} {stat : ℕ → α → ℝ} {m s : ℝ},
      AsymptoticallyNormal stat μ m s →
      Tendsto (fun n => Expect μ (stat n) - m) atTop (nhds 0) ∧
        Tendsto (fun n => ProbabilityTheory.variance (stat n) μ) atTop (nhds (s^2)) := by
  intro μ stat m s h
  refine ⟨?_, h.2⟩
  have hconst : Tendsto (fun _ : ℕ => m) atTop (nhds m) := tendsto_const_nhds
  simpa using (h.1.sub hconst)

/-!
## Поріг ПЕ та нерівність Чебишева

Критерій ПЕ визначає поріг h_PE = E[Λ|H₀] + √(Var[Λ|H₀]/ε)
для контролю FAR ≤ ε через нерівність Чебишева.
-/

/-- Формула порогу ПЕ (Probability Error bound). -/
noncomputable def PE_threshold (mean variance ε : ℝ) : ℝ :=
  mean + Real.sqrt (variance / ε)

/-- Формула зміни порогу при зміні FAR. -/
theorem PE_threshold_difference (σ ε₁ ε₂ : ℝ) (hε₁ : 0 < ε₁) (hε₂ : 0 < ε₂) :
    PE_threshold 0 (σ^2) ε₂ - PE_threshold 0 (σ^2) ε₁ =
      Real.sqrt (σ^2 / ε₂) - Real.sqrt (σ^2 / ε₁) := by
  simp [PE_threshold]

/-- Теорема 3 (повна): асимптотична еквівалентність порогів ПЕ та NP
    для нормально розподілених статистик. -/
theorem theorem3_PE_NP_threshold_equivalence
    (m₀ σ₀ : ℝ) (hσ : 0 < σ₀) (ε : ℝ) (hε : 0 < ε) :
    ∃ (h_PE h_NP : ℝ),
      h_PE = PE_threshold m₀ (σ₀^2) ε ∧
      -- Для нормального розподілу, поріг NP визначається через квантиль,
      -- який асимптотично еквівалентний порогу ПЕ з поправкою 1/√ε
      h_PE - m₀ = σ₀ * Real.sqrt (1 / ε) := by
  refine ⟨PE_threshold m₀ (σ₀^2) ε, m₀ + σ₀ * Real.sqrt (1 / ε), ?_, ?_⟩
  · rfl
  · simp only [PE_threshold, add_sub_cancel_left]
    rw [Real.sqrt_div (sq_nonneg σ₀), one_div]
    rw [Real.sqrt_sq (le_of_lt hσ), Real.sqrt_inv, div_eq_mul_inv]

/-- Критерій Ю: функціонал для асимптотично нормальних статистик. -/
noncomputable def criterion_Y (G₀ G₁ E₀ E₁ : ℝ) : ℝ :=
  (Real.sqrt G₀ + Real.sqrt G₁)^2 / (E₁ - E₀)^2

/-- Властивість критерію Ю: невід'ємність. -/
theorem criterion_Y_nonneg (G₀ G₁ E₀ E₁ : ℝ) (hG₀ : 0 ≤ G₀) (hG₁ : 0 ≤ G₁) :
    0 ≤ criterion_Y G₀ G₁ E₀ E₁ := by
  unfold criterion_Y
  apply div_nonneg
  · exact sq_nonneg _
  · exact sq_nonneg _

/-- Зв'язок критеріїв КУ1 та Ю: КУ1 ≤ Ю за нерівністю між середнім арифметичним
    та середнім квадратичним. -/
theorem KU1_le_criterion_Y (G₀ G₁ E₀ E₁ : ℝ) (hG₀ : 0 ≤ G₀) (hG₁ : 0 ≤ G₁)
    (hE : E₁ ≠ E₀) :
    (G₀ + G₁) / (E₁ - E₀)^2 ≤ criterion_Y G₀ G₁ E₀ E₁ := by
  unfold criterion_Y
  have hden_pos : 0 < (E₁ - E₀)^2 := sq_pos_of_ne_zero (sub_ne_zero.mpr hE)
  -- Достатньо показати G₀ + G₁ ≤ (√G₀ + √G₁)²
  apply div_le_div_of_nonneg_right _ hden_pos.le
  -- (√G₀ + √G₁)² = G₀ + 2√G₀√G₁ + G₁ ≥ G₀ + G₁
  have hsq0 : Real.sqrt G₀ ^ 2 = G₀ := Real.sq_sqrt hG₀
  have hsq1 : Real.sqrt G₁ ^ 2 = G₁ := Real.sq_sqrt hG₁
  have h_expand : (Real.sqrt G₀ + Real.sqrt G₁)^2 =
      G₀ + G₁ + 2 * Real.sqrt G₀ * Real.sqrt G₁ := by
    rw [add_sq, hsq0, hsq1]; ring
  rw [h_expand]
  have h_cross_nonneg : 0 ≤ 2 * Real.sqrt G₀ * Real.sqrt G₁ := by
    apply mul_nonneg
    apply mul_nonneg
    · linarith
    · exact Real.sqrt_nonneg G₀
    exact Real.sqrt_nonneg G₁
  linarith

end GSA.Part2
