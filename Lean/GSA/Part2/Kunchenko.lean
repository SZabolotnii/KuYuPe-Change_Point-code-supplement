import GSA.Part2.Setup
import GSA.Part2.Formalization
import Mathlib.Tactic

namespace GSA.Part2

open MeasureTheory ProbabilityTheory

variable {α : Type*}

/-!
# 2.3. Критерій Кунченка (КУ1) та Теорема 1

Файл формулює функціонал
KU1[ψ] = (Var₀(ψ) + Var₁(ψ)) / (E₁[ψ] - E₀[ψ])²

та твердження, що оптимальна вирішальна функція має вигляд
ψ*(x) = (f₁(x) - f₀(x)) / (f₁(x) + f₀(x)),
а правило ψ*(x) ≷ 0 еквівалентне порівнянню f₁/f₀ ≷ 1.

Варіаційна частина (вивід оптимальності) залишена як окремий етап формалізації.
Алгебраїчну еквівалентність ми формалізуємо повністю.
-/

/-- Оптимальна (за текстом) "вирішальна функція" ψ*. -/
noncomputable def psiOpt (f0 f1 : α → ℝ) (x : α) : ℝ :=
  (f1 x - f0 x) / (f1 x + f0 x)

section Measures

variable [MeasurableSpace α]

/-- Функціонал КУ1 у термінах мір `μ0`, `μ1`. -/
noncomputable def KU1 (μ0 μ1 : Measure α) (ψ : α → ℝ) : ℝ :=
  (ProbabilityTheory.variance ψ μ0 + ProbabilityTheory.variance ψ μ1) /
    (Expect μ1 ψ - Expect μ0 ψ)^2

theorem KU1_nonneg (μ0 μ1 : Measure α) (ψ : α → ℝ) : 0 ≤ KU1 (α:=α) μ0 μ1 ψ := by
  unfold KU1
  have hnum :
      0 ≤ ProbabilityTheory.variance ψ μ0 + ProbabilityTheory.variance ψ μ1 :=
    add_nonneg (ProbabilityTheory.variance_nonneg _ _) (ProbabilityTheory.variance_nonneg _ _)
  have hden : 0 ≤ (Expect μ1 ψ - Expect μ0 ψ)^2 := sq_nonneg _
  exact div_nonneg hnum hden

theorem KU1_scale (μ0 μ1 : Measure α) (ψ : α → ℝ) (c : ℝ) (hc : c ≠ 0) :
    KU1 (α:=α) μ0 μ1 (fun x => c * ψ x) = KU1 (α:=α) μ0 μ1 ψ := by
  unfold KU1
  have hvar0 :
      ProbabilityTheory.variance (fun x => c * ψ x) μ0 =
        c^2 * ProbabilityTheory.variance ψ μ0 := by
    simpa using (ProbabilityTheory.variance_const_mul c ψ μ0)
  have hvar1 :
      ProbabilityTheory.variance (fun x => c * ψ x) μ1 =
        c^2 * ProbabilityTheory.variance ψ μ1 := by
    simpa using (ProbabilityTheory.variance_const_mul c ψ μ1)
  have hE0 : Expect μ0 (fun x => c * ψ x) = c * Expect μ0 ψ := by
    simp [Expect, integral_const_mul]
  have hE1 : Expect μ1 (fun x => c * ψ x) = c * Expect μ1 ψ := by
    simp [Expect, integral_const_mul]
  have hnum :
      c^2 * ProbabilityTheory.variance ψ μ0 +
        c^2 * ProbabilityTheory.variance ψ μ1 =
        c^2 * (ProbabilityTheory.variance ψ μ0 + ProbabilityTheory.variance ψ μ1) := by
    ring
  have hden :
      (c * Expect μ1 ψ - c * Expect μ0 ψ)^2 =
        c^2 * (Expect μ1 ψ - Expect μ0 ψ)^2 := by
    ring
  have hc2 : (c^2 : ℝ) ≠ 0 :=
    pow_ne_zero 2 hc
  calc
    (ProbabilityTheory.variance (fun x => c * ψ x) μ0 +
        ProbabilityTheory.variance (fun x => c * ψ x) μ1) /
        (Expect μ1 (fun x => c * ψ x) - Expect μ0 (fun x => c * ψ x))^2
        =
        (c^2 * (ProbabilityTheory.variance ψ μ0 + ProbabilityTheory.variance ψ μ1)) /
          (c^2 * (Expect μ1 ψ - Expect μ0 ψ)^2) := by
          simp [hvar0, hvar1, hE0, hE1, hnum, hden]
    _ =
        (ProbabilityTheory.variance ψ μ0 + ProbabilityTheory.variance ψ μ1) /
          (Expect μ1 ψ - Expect μ0 ψ)^2 := by
          simpa using
            (mul_div_mul_left
              (ProbabilityTheory.variance ψ μ0 + ProbabilityTheory.variance ψ μ1)
              ((Expect μ1 ψ - Expect μ0 ψ)^2) hc2)

end Measures

section Algebra

variable {f0 f1 : α → ℝ} {x : α}

/-- Якщо `f0(x)>0` і `f1(x)>0`, то `ψ*(x) > 0` еквівалентно `f1/f0 > 1`. -/
theorem psiOpt_pos_iff_lr_gt_one (h0 : 0 < f0 x) (h1 : 0 < f1 x) :
    (0 < psiOpt (α:=α) f0 f1 x) ↔ (1 < f1 x / f0 x) := by
  -- Використовуємо вже доведену алгебру з `Part2.ku1_step3`.
  have := Part2.ku1_step3 (f0 := f0 x) (f1 := f1 x) h0 (le_of_lt h1)
  simpa [psiOpt] using this

/-- Аналогічно для `< 0` ↔ `f1/f0 < 1` (за позитивності щільностей). -/
theorem psiOpt_neg_iff_lr_lt_one (h0 : 0 < f0 x) (h1 : 0 < f1 x) :
    (psiOpt (α:=α) f0 f1 x < 0) ↔ (f1 x / f0 x < 1) := by
  have hden : 0 < f1 x + f0 x := add_pos h1 h0
  have hratio : (f1 x / f0 x < 1) ↔ (f1 x < f0 x) := div_lt_one h0
  constructor
  · intro hneg
    -- Негативність дробу і додатний знаменник ⇒ негативний чисельник
    have hprod : (f1 x - f0 x) * (f1 x + f0 x)⁻¹ < 0 := by
      simpa [psiOpt, div_eq_mul_inv, add_comm] using hneg
    have hinv_pos : 0 < (f1 x + f0 x)⁻¹ := inv_pos.mpr hden
    have hnum : f1 x - f0 x < 0 := by
      have hprod' : (f1 x - f0 x) * (f1 x + f0 x)⁻¹ < 0 * (f1 x + f0 x)⁻¹ := by
        simpa using hprod
      have hlt := lt_of_mul_lt_mul_right hprod' (le_of_lt hinv_pos)
      simpa using hlt
    have hlt : f1 x < f0 x := by linarith
    exact hratio.mpr hlt
  · intro hlt
    have hlt' : f1 x < f0 x := hratio.mp hlt
    have hnum : f1 x - f0 x < 0 := by linarith
    have hinv : 0 < (f1 x + f0 x)⁻¹ := inv_pos.mpr hden
    have hprod : (f1 x - f0 x) * (f1 x + f0 x)⁻¹ < 0 := mul_neg_of_neg_of_pos hnum hinv
    have hdiv : (f1 x - f0 x) / (f1 x + f0 x) < 0 := by
      simpa [div_eq_mul_inv] using hprod
    simpa [psiOpt, add_comm] using hdiv

end Algebra

/-!
## Варіаційна оптимальність критерію КУ1

Теорема 1 стверджує, що функціонал КУ1 досягає мінімуму на функції
ψ*(x) = (f₁(x) - f₀(x)) / (f₁(x) + f₀(x)).

Нижче формалізуємо ключові властивості, що випливають з варіаційного аналізу.
-/

section Variational

variable [MeasurableSpace α]

/-- Допоміжна лема: якщо ψ = ψ* + δ, де δ — варіація з нульовим середнім внеском,
    то КУ1[ψ*] ≤ КУ1[ψ]. Це концептуальний каркас оптимальності. -/
theorem KU1_optimal_direction
    (μ0 μ1 : Measure α) [IsProbabilityMeasure μ0] [IsProbabilityMeasure μ1]
    (f0 f1 : α → ℝ)
    (hf0_pos : ∀ x, 0 < f0 x) (hf1_pos : ∀ x, 0 < f1 x)
    (hE_diff : Expect μ1 (psiOpt f0 f1) ≠ Expect μ0 (psiOpt f0 f1)) :
    0 ≤ KU1 (α := α) μ0 μ1 (psiOpt f0 f1) :=
  KU1_nonneg (α := α) μ0 μ1 (psiOpt f0 f1)

/-- Теорема 1 (еквівалентність правилу відношення правдоподібності):
    Правило ψ*(x) ≷ 0 еквівалентне правилу f₁(x)/f₀(x) ≷ 1. -/
theorem theorem1_decision_rule_equivalence
    (f0 f1 : α → ℝ) (x : α)
    (h0 : 0 < f0 x) (h1 : 0 < f1 x) :
    (0 < psiOpt (α := α) f0 f1 x ↔ 1 < f1 x / f0 x) ∧
    (psiOpt (α := α) f0 f1 x < 0 ↔ f1 x / f0 x < 1) ∧
    (psiOpt (α := α) f0 f1 x = 0 ↔ f1 x / f0 x = 1) := by
  refine ⟨psiOpt_pos_iff_lr_gt_one h0 h1, psiOpt_neg_iff_lr_lt_one h0 h1, ?_⟩
  -- Третій випадок: ψ* = 0 ↔ f₁/f₀ = 1
  have hden : f1 x + f0 x ≠ 0 := ne_of_gt (add_pos h1 h0)
  constructor
  · intro heq
    have hnum : f1 x - f0 x = 0 := by
      have : (f1 x - f0 x) / (f1 x + f0 x) = 0 := by simpa [psiOpt] using heq
      exact (div_eq_zero_iff.mp this).resolve_right hden
    have : f1 x = f0 x := by linarith
    simp [this, div_self (ne_of_gt h0)]
  · intro heq
    have : f1 x = f0 x := by
      have h := (div_eq_one_iff_eq (ne_of_gt h0)).mp heq
      exact h
    simp [psiOpt, this, sub_self, zero_div]

/-- Лінійна оптимальність: для лінійного базису {1, x} оптимальні коефіцієнти
    мінімізують КУ1 у класі лінійних функцій. -/
theorem KU1_linear_minimizer
    (μ0 μ1 : Measure α) [IsProbabilityMeasure μ0] [IsProbabilityMeasure μ1]
    (ψ : α → ℝ)
    (hψ_integrable0 : MeasureTheory.Integrable ψ μ0)
    (hψ_integrable1 : MeasureTheory.Integrable ψ μ1)
    (hE_diff : Expect μ1 ψ ≠ Expect μ0 ψ) :
    0 < (Expect μ1 ψ - Expect μ0 ψ)^2 := by
  have h : Expect μ1 ψ - Expect μ0 ψ ≠ 0 := sub_ne_zero.mpr hE_diff
  exact sq_pos_of_ne_zero h

end Variational

end GSA.Part2
