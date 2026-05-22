import Mathlib.LinearAlgebra.Matrix.Defs
import Mathlib.LinearAlgebra.Matrix.ToLinearEquiv
import Mathlib.Probability.Moments.Covariance
import GSA.Part2.Setup

open scoped BigOperators
open MeasureTheory ProbabilityTheory

namespace GSA.Part2

/-!
# 2.3.2–2.3.3. Система F K = Y

У тексті:
- маємо базисні функції φᵢ, i=1..s;
- статистика Λ^{(s)}(x) = ∑ Kᵢ φᵢ(x);
- варіаційний вивід дає лінійну систему F K = Y,
  де Fᵢⱼ = cov₀(φᵢ,φⱼ)+cov₁(φᵢ,φⱼ),
      Yⱼ  = E₁[φⱼ] - E₀[φⱼ].

Нижче — каркас в термінах матриць `Matrix (Fin s) (Fin s) ℝ`.
-/

variable {α : Type*} [MeasurableSpace α]
variable (μ0 μ1 : Measure α) [MeasureTheory.IsFiniteMeasure μ0] [MeasureTheory.IsFiniteMeasure μ1]
variable {s : ℕ} (φ : Fin s → α → ℝ)

/-- Компонента Yⱼ = E₁[φⱼ] - E₀[φⱼ]. -/
noncomputable def Yvec : Fin s → ℝ := fun j =>
  (∫ x, φ j x ∂ μ1) - (∫ x, φ j x ∂ μ0)

/-- Коваріація (в mathlib є `ProbabilityTheory.covariance`, але для загальної міри потрібні умови).
Тут залишаємо окрему дефініцію як проміжний каркас. -/
noncomputable def cov (μ : Measure α) (X Y : α → ℝ) : ℝ :=
  (∫ x, (X x - (∫ t, X t ∂ μ)) * (Y x - (∫ t, Y t ∂ μ)) ∂ μ)

/-- Матриця Fᵢⱼ = cov₀(φᵢ,φⱼ)+cov₁(φᵢ,φⱼ). -/
noncomputable def Fmat : Matrix (Fin s) (Fin s) ℝ :=
  fun i j => cov μ0 (φ i) (φ j) + cov μ1 (φ i) (φ j)

omit [IsFiniteMeasure μ0] [IsFiniteMeasure μ1] in
/-- Лінійна система (існування/єдиність розв'язку залежить від невиродженості F). -/
theorem has_solution_FK_eq_Y
    (hInv : IsUnit (Matrix.det (Fmat (μ0:=μ0) (μ1:=μ1) φ))) :
    ∃ K : Fin s → ℝ, (Fmat (μ0:=μ0) (μ1:=μ1) φ).mulVec K = Yvec (μ0:=μ0) (μ1:=μ1) φ := by
  classical
  let A : Matrix (Fin s) (Fin s) ℝ := Fmat (μ0:=μ0) (μ1:=μ1) φ
  let Y : Fin s → ℝ := Yvec (μ0:=μ0) (μ1:=μ1) φ
  -- Інстанс оберненості з IsUnit det.
  let _ := Matrix.invertibleOfIsUnitDet (A := A) hInv
  refine ⟨A⁻¹.mulVec Y, ?_⟩
  calc
    A.mulVec (A⁻¹.mulVec Y)
        = (A * A⁻¹).mulVec Y := by
            simp [Matrix.mulVec_mulVec]
    _ = (1 : Matrix (Fin s) (Fin s) ℝ).mulVec Y := by
            simp [Matrix.mul_inv_of_invertible]
    _ = Y := by
            simp

end GSA.Part2
