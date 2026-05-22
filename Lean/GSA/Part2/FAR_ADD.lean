import GSA.Part2.Setup
import Mathlib.Probability.Moments.Variance
import Mathlib.Tactic

open MeasureTheory ProbabilityTheory
open scoped BigOperators

namespace GSA.Part2

/-!
# 2.4.4. FAR–ADD компроміс та Теорема 6

У файлі Теорема 6 дає аналітичний поріг (через середнє та дисперсію)
для контролю FAR на заданому рівні ε.

Нижче — лема-скелет, яка зводить це до Chebyshev (`meas_ge_le_variance_div_sq`).
-/

variable {Ω : Type*} [MeasurableSpace Ω]
variable (μ : Measure Ω) [MeasureTheory.IsFiniteMeasure μ]
variable (X : Ω → ℝ)

/-- Подія перевищення порогу `m + t`, де `m = E[X]`. -/
def exceed (t : ℝ) : Set Ω :=
  {ω | (∫ x, X x ∂ μ) + t ≤ X ω}

theorem exceed_le_by_chebyshev (hX : MeasureTheory.MemLp X 2 μ) (t : ℝ) (ht : 0 < t) :
    μ (exceed (μ:=μ) X t) ≤ ENNReal.ofReal (ProbabilityTheory.variance X μ / t^2) := by
  -- 1) exceed ⊆ {ω | t ≤ |X ω - E[X]|}
  have hsubset : exceed (μ:=μ) X t ⊆ {ω | t ≤ |X ω - μ[X]|} := by
    intro ω hω
    have hω' : (∫ x, X x ∂ μ) + t ≤ X ω := by
      simpa [exceed] using hω
    have h1 : t ≤ X ω - (∫ x, X x ∂ μ) := by linarith
    have h2 : 0 ≤ X ω - (∫ x, X x ∂ μ) := by linarith
    have habs :
        |X ω - (∫ x, X x ∂ μ)| = X ω - (∫ x, X x ∂ μ) := by
      exact abs_of_nonneg h2
    have : t ≤ |X ω - (∫ x, X x ∂ μ)| := by
      simpa [habs] using h1
    -- Переписуємо `μ[X]` як інтеграл.
    simpa using this
  -- 2) застосувати Chebyshev
  have hcheb :=
    ProbabilityTheory.meas_ge_le_variance_div_sq (μ:=μ) (X:=X) hX (c:=t) ht
  exact (measure_mono hsubset).trans hcheb

/-- Зручний вибір порогу t = sqrt(Var/ε) (дає bound ≤ ε). -/
theorem exceed_le_eps (hX : MeasureTheory.MemLp X 2 μ) (ε : ℝ) (hε : 0 < ε)
    (hvar : 0 < ProbabilityTheory.variance X μ) :
    μ (exceed (μ:=μ) X (Real.sqrt (ProbabilityTheory.variance X μ / ε))) ≤ ENNReal.ofReal ε := by
  have ht : 0 < Real.sqrt (ProbabilityTheory.variance X μ / ε) := by
    refine Real.sqrt_pos.2 ?_
    exact div_pos hvar hε
  have hle := exceed_le_by_chebyshev (μ:=μ) (X:=X) (hX:=hX) (t:=_) ht
  -- спрощуємо праву частину до ε
  have hsq : (Real.sqrt (ProbabilityTheory.variance X μ / ε)) ^ 2 =
      ProbabilityTheory.variance X μ / ε := by
    simp [Real.sq_sqrt, div_nonneg (ProbabilityTheory.variance_nonneg _ _) (le_of_lt hε)]
  have hvar' : (ProbabilityTheory.variance X μ) ≠ 0 := by linarith
  have hcalc : ProbabilityTheory.variance X μ / (ProbabilityTheory.variance X μ / ε) = ε := by
    field_simp [hvar', hε.ne']
  -- підстановка та алгебра
  simpa [hsq, hcalc] using hle

/-!
Corollary 1/2 з тексту можна отримати, підставляючи відповідні статистики `G_k` та беручи `sup`/`max`.
-/

end GSA.Part2
