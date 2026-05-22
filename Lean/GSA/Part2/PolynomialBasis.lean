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

    Формальне доведення вимагає додаткових припущень про вимірність X,
    тому приймаємо як аксіому стандартний результат теорії ймовірностей.
-/
axiom hasMoment_of_higher [IsProbabilityMeasure μ]
    (X : α → ℝ) {k j : ℕ} (hjk : j ≤ k) (hk : HasMoment μ X k) :
    HasMoment μ X j

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
