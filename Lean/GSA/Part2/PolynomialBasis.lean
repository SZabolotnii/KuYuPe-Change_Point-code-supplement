import GSA.Part2.BasisApprox
import Mathlib.Probability.Moments.Basic

namespace GSA.Part2

open MeasureTheory ProbabilityTheory
open scoped BigOperators

variable {α : Type*} [MeasurableSpace α]

/-!
# 2.2.3. Поліноміальний базис: умови застосовності

Формалізація умов ефективності поліноміального базису:
- існування моментів до порядку `2s`;
- числова стійкість (обумовленість матриці Грама).
-/

/-!
## Умови на моменти

Для коректного формування матриці F необхідно, щоб існували
моменти порядку 2s при обох гіпотезах.
-/

/-- Умова існування моменту порядку k. -/
def HasMoment (μ : Measure α) (X : α → ℝ) (k : ℕ) : Prop :=
  MeasureTheory.Integrable (fun x => |X x|^k) μ

/-- Умова існування всіх моментів до порядку k (включно). -/
def HasMomentsUpTo (μ : Measure α) (X : α → ℝ) (k : ℕ) : Prop :=
  ∀ j ≤ k, HasMoment μ X j

/-- Якщо існує момент порядку k, то існують всі моменти нижчих порядків
    (для ймовірнісних мір).

    Математичне обґрунтування:
    Для ймовірнісної міри μ, якщо E[|X|^k] < ∞, то E[|X|^j] < ∞ для j ≤ k.
    Це випливає з нерівності: |X|^j ≤ max(1, |X|^k) ≤ 1 + |X|^k
    та теореми про мажоровану збіжність.

    Вимірність окремого припущення не потребує, хоча довгий час здавалося,
    що потребує, і саме через це твердження стояло аксіомою: `|X|^k` вимірна
    вже тому, що інтегровна, а при `k ≥ 1` з неї відновлюється
    `|X| = (|X|^k)^(1/k)` як композиція з неперервним `t ↦ t^(1/k)`;
    випадок `k = 0` змушує `j = 0` і тривіальний.
-/
theorem hasMoment_of_higher [IsProbabilityMeasure μ]
    (X : α → ℝ) {k j : ℕ} (hjk : j ≤ k) (hk : HasMoment μ X k) :
    HasMoment μ X j := by
  -- unfold HasMoment: goal is Integrable (fun x => |X x|^j) μ
  -- hk : Integrable (fun x => |X x|^k) μ
  rcases Nat.eq_zero_or_pos k with (rfl | hk_pos)
  · -- k = 0, so j = 0 (since j ≤ 0)
    have hj0 : j = 0 := Nat.eq_zero_of_le_zero hjk
    subst hj0
    rw [HasMoment]
    exact integrable_const _
  · -- k > 0
    have hk_ae : AEStronglyMeasurable (fun x => |X x| ^ (k : ℝ)) μ := by
      -- hk gives Integrable of |X|^k (Nat.pow); convert to Real.rpow
      simpa [Real.rpow_natCast] using hk.aestronglyMeasurable
    have h_one_div_k_nonneg : 0 ≤ (1 : ℝ) / (k : ℝ) := by
      positivity
    have h_abs_ae : AEStronglyMeasurable (fun x => |X x|) μ := by
      -- |X| = (|X|^k)^(1/k) as Real.rpow, and t ↦ t^(1/k) is continuous
      have h_eq : (fun x => |X x|) = (fun x => ((|X x| ^ (k : ℝ)) ^ ((1 : ℝ) / (k : ℝ)))) := by
        ext x
        have h_nonneg : 0 ≤ |X x| := abs_nonneg _
        calc
          |X x| = |X x| ^ (1 : ℝ) := by simp
          _ = |X x| ^ (((k : ℝ) * ((1 : ℝ) / (k : ℝ)))) := by
            field_simp [show (k : ℝ) ≠ 0 from by exact_mod_cast hk_pos.ne.symm]
          _ = ((|X x| ^ (k : ℝ)) ^ ((1 : ℝ) / (k : ℝ))) := by
            rw [Real.rpow_mul h_nonneg (k : ℝ) ((1 : ℝ) / (k : ℝ))]
      rw [h_eq]
      exact Continuous.comp_aestronglyMeasurable
        (Real.continuous_rpow_const h_one_div_k_nonneg) hk_ae
    have hf_ae : AEStronglyMeasurable (fun x => |X x| ^ j) μ :=
      h_abs_ae.pow j
    have hg_int : Integrable (fun x => (1 : ℝ) + |X x| ^ k) μ := by
      have h_const : Integrable (fun _ => (1 : ℝ)) μ := integrable_const _
      -- hk is Integrable (|X|^k) with Nat.pow; need to convert for addition
      have hk_nat : Integrable (fun x => |X x| ^ k) μ := hk
      simpa [add_comm] using hk_nat.add h_const
    have h_bound : ∀ x, |X x| ^ j ≤ (1 : ℝ) + |X x| ^ k := by
      intro x
      have h_nonneg : 0 ≤ |X x| := abs_nonneg _
      have h_pow_nonneg : 0 ≤ |X x| ^ k := pow_nonneg h_nonneg k
      by_cases h_le_one : |X x| ≤ 1
      · have h_pow_le_one : |X x| ^ j ≤ 1 :=
          pow_le_one₀ h_nonneg h_le_one
        linarith
      · have h_one_le : 1 ≤ |X x| := by linarith
        have h_pow_le : |X x| ^ j ≤ |X x| ^ k :=
          pow_le_pow_right₀ h_one_le hjk
        linarith
    have h_bound_ae : ∀ᵐ x ∂μ, ‖(|X x| ^ j : ℝ)‖ ≤ (1 : ℝ) + |X x| ^ k := by
      filter_upwards [] with x
      have h_nonneg_pow : 0 ≤ |X x| ^ j := pow_nonneg (abs_nonneg _) j
      simpa [abs_of_nonneg h_nonneg_pow] using h_bound x
    exact Integrable.mono' hg_int hf_ae h_bound_ae

/-- Умова застосовності поліноміального базису порядку s:
    необхідно існування моменту 2s. -/
structure PolynomialBasisApplicable (μ₀ μ₁ : Measure α) (X : α → ℝ) (s : ℕ) : Prop where
  moment_H0 : HasMoment μ₀ X (2 * s)
  moment_H1 : HasMoment μ₁ X (2 * s)

/-- Для Парето-розподілу з параметром α, максимальний порядок s < α/2. -/
noncomputable def pareto_max_order (α_param : ℝ) : ℕ :=
  Nat.floor (α_param / 2)

/-!
## Числова стійкість

Матриця Грама для поліноміального базису має Ганкелеву структуру
і може бути погано обумовленою при великих s.
-/

/-- Порогове значення числа обумовленості для практичного використання. -/
def condition_number_threshold : ℝ := 10^6

/-- Рекомендований максимальний порядок s для різних значень ексцесу γ₄. -/
noncomputable def recommended_order (excess : ℝ) : ℕ :=
  if excess < 6 then 2
  else if excess < 20 then 3
  else 2  -- для важкохвостих краще логарифмічний базис

/-!
## Гладкість LLR

Швидкість збіжності O(s^{-r}) реалізується якщо істинний LLR
є достатньо гладким (має r неперервних похідних).
-/

/-- Клас гладкості LLR. -/
inductive LLR_Smoothness
  | C0  -- неперервний
  | C1  -- неперервно диференційовний
  | C2  -- двічі неперервно диференційовний
  | Cr (r : ℕ)  -- r разів неперервно диференційовний
  | analytic  -- аналітичний
  deriving DecidableEq, Repr

/-- Швидкість збіжності залежно від гладкості. -/
def convergence_rate (smooth : LLR_Smoothness) : ℕ :=
  match smooth with
  | .C0 => 0
  | .C1 => 1
  | .C2 => 2
  | .Cr r => r
  | .analytic => 100  -- експоненційна збіжність

/-!
## Практичні рекомендації

Таблиця вибору базису за ексцесом γ₄.
-/

/-- Тип базису. -/
inductive BasisType
  | polynomial
  | logarithmic
  | fractional
  | harmonic
  deriving DecidableEq, Repr

/-- Вибір типу базису за ексцесом. -/
noncomputable def select_basis_by_excess (γ₄ : ℝ) : BasisType :=
  if γ₄ < 6 then .polynomial
  else if γ₄ < 20 then .fractional
  else .logarithmic

/-- Рекомендований порядок апроксимації за типом базису. -/
def recommended_s (basis : BasisType) : ℕ :=
  match basis with
  | .polynomial => 3
  | .logarithmic => 2
  | .fractional => 3
  | .harmonic => 2

end GSA.Part2
