import GSA.Part2.Setup

namespace GSA.Part2

open MeasureTheory ProbabilityTheory InformationTheory

variable {α : Type*} [MeasurableSpace α]

/-!
# 2.1.1. Малий відносний розлад

У `Part2_review.md`:

ρ = D_KL(f₁ ∥ f₀) / sqrt(Var_{f₀}( log(f₁/f₀) )).

У термінах мір: `D_KL(μ1 ∥ μ0) = klDiv μ1 μ0`, а `log(f₁/f₀)` — `llr μ1 μ0`.
-/

/-- Параметр відносного розладу ρ (в `ℝ`).  
Зауваження: `klDiv` має тип `ℝ≥0∞`, тому тут використовується `ENNReal.toReal`.
Для строгих тверджень варто додавати умови `klDiv μ1 μ0 < ⊤`. -/
noncomputable def relDisorder (μ0 μ1 : Measure α) : ℝ :=
  (InformationTheory.klDiv μ1 μ0).toReal /
    Real.sqrt (ProbabilityTheory.variance (LLR μ1 μ0) μ0)

inductive DisorderRegime
  | small   -- ρ < 1
  | medium  -- 1 ≤ ρ < 3
  | large   -- 3 ≤ ρ
  deriving DecidableEq, Repr

/-- Класифікація режимів за ρ. -/
noncomputable def regime (μ0 μ1 : Measure α) : DisorderRegime :=
  let ρ := relDisorder (α:=α) μ0 μ1
  if _h1 : ρ < 1 then .small
  else if _h3 : ρ < 3 then .medium
  else .large

/-- Невід'ємність ρ (загальний наслідок невід'ємності `klDiv` і `variance`). -/
theorem relDisorder_nonneg (μ0 μ1 : Measure α) : 0 ≤ relDisorder (α:=α) μ0 μ1 := by
  unfold relDisorder
  have hnum : 0 ≤ (InformationTheory.klDiv μ1 μ0).toReal := ENNReal.toReal_nonneg
  have hden : 0 ≤ Real.sqrt (ProbabilityTheory.variance (LLR μ1 μ0) μ0) := by
    exact Real.sqrt_nonneg _
  exact div_nonneg hnum hden

/-!
Подальші кроки: коректність інтерпретації режимів (умови на KL/дисперсію),
та зв'язок з припущеннями статті.
-/

end GSA.Part2
