import Mathlib.Data.Matrix.Basic
import Mathlib.Data.Matrix.Mul
import Mathlib.Data.Fin.VecNotation
import Mathlib.Tactic.NormNum
import Mathlib.Tactic.FinCases
import GSA.Part2.LinearSystem

open scoped BigOperators
open Matrix

namespace GSA.Part2.MAry

/-!
# Багатоальтернативний випадок: попарні правила КУ1 та їх узгодженість

Розширення §2.3 на `M` гіпотез. У бінарній задачі система одна:
`F K = Y`, `F = C⁰ + C¹`, `Y = μ¹ - μ⁰`.

За `M` гіпотез кожна пара `(m,n)` дає **свою** матрицю `F_{mn} = C^m + C^n`.
Вектори `Y` адитивні завжди, а коефіцієнти `K` — лише коли `F` спільна.
Звідси два результати:

* **Т2 (звідність).** За гомоскедастичності `C^m = C₀` попарні статистики є
  різницями значень однієї скалярної функції (`lam_eq_half_g_sub`), тому
  попарні рішення переходові (`transitive_of_potential`).
* **Т1 (незвідність).** У загальному випадку переходовості немає:
  `exists_intransitive` дає явний контрприклад на `M = 3`, `s = 2`.

Контрприклад знайдено чисельно (`Zhila-diss/verification/p2_intransitivity.py`)
і тут подано в точних раціональних числах.
-/

variable {s : ℕ} {ι : Type*}

/-- `Y_{mn} = μ^n - μ^m` — права частина нормальної системи для пари. -/
def Ypair (μ : ι → (Fin s → ℝ)) (m n : ι) : Fin s → ℝ := μ n - μ m

/-- `F_{mn} = C^m + C^n` — матриця нормальної системи для пари. -/
def Fpair (C : ι → Matrix (Fin s) (Fin s) ℝ) (m n : ι) : Matrix (Fin s) (Fin s) ℝ :=
  C m + C n

/-- **L1.** Адитивність правих частин: `Y_{mn} + Y_{np} = Y_{mp}`. -/
theorem Ypair_add (μ : ι → (Fin s → ℝ)) (m n p : ι) :
    Ypair μ m n + Ypair μ n p = Ypair μ m p := by
  simp only [Ypair]
  abel

/-- **L3.** Якщо матриця нормальної системи **спільна** для двох пар, то їхні
розв'язки додаються: сума є розв'язком для складеної пари.

Це і є точна умова звідності: адитивність `K` успадковується від адитивності `Y`
лише за спільної `A`. За `A = F_{mn}`, що залежить від пари, висновок не діє. -/
theorem Ksol_add (μ : ι → (Fin s → ℝ)) (A : Matrix (Fin s) (Fin s) ℝ)
    {m n p : ι} {K₁ K₂ : Fin s → ℝ}
    (h₁ : A *ᵥ K₁ = Ypair μ m n) (h₂ : A *ᵥ K₂ = Ypair μ n p) :
    A *ᵥ (K₁ + K₂) = Ypair μ m p := by
  rw [Matrix.mulVec_add, h₁, h₂, Ypair_add]

/-- За гомоскедастичності матриця нормальної системи одна для всіх пар. -/
theorem Fpair_const {C : ι → Matrix (Fin s) (Fin s) ℝ} {C₀ : Matrix (Fin s) (Fin s) ℝ}
    (hC : ∀ m, C m = C₀) (m n : ι) : Fpair C m n = C₀ + C₀ := by
  simp [Fpair, hC]

/-- Симетрична матриця дає симетричну білінійну форму. -/
theorem dot_mulVec_comm {B : Matrix (Fin s) (Fin s) ℝ} (hB : Bᵀ = B)
    (u v : Fin s → ℝ) : (B *ᵥ u) ⬝ᵥ v = (B *ᵥ v) ⬝ᵥ u := by
  have hsym : ∀ i j, B j i = B i j := fun i j => by
    simpa [Matrix.transpose_apply] using congrFun (congrFun hB i) j
  simp only [Matrix.mulVec, dotProduct, Finset.sum_mul]
  rw [Finset.sum_comm]
  refine Finset.sum_congr rfl fun i _ => Finset.sum_congr rfl fun j _ => ?_
  rw [hsym]; ring

/-- Потенціал узгодженого правила: `g_m(x) = ⟪B μ^m, x⟫ - ½⟪B μ^m, μ^m⟫`,
де `B` грає роль `C₀⁻¹`. -/
noncomputable def gpot (B : Matrix (Fin s) (Fin s) ℝ) (μ : ι → (Fin s → ℝ)) (m : ι)
    (x : Fin s → ℝ) : ℝ :=
  (B *ᵥ μ m) ⬝ᵥ x - (1 / 2) * ((B *ᵥ μ m) ⬝ᵥ μ m)

/-- **L4 — ядро Т2.** За гомоскедастичності попарна статистика КУ1 дорівнює
половині різниці потенціалів: `Λ_{mn}(x) = ½ (g_n(x) - g_m(x))`.

Отже всі попарні статистики породжені **однією** скалярною функцією. -/
theorem lam_eq_half_g_sub {B : Matrix (Fin s) (Fin s) ℝ} (hB : Bᵀ = B)
    (μ : ι → (Fin s → ℝ)) (m n : ι) (x : Fin s → ℝ) :
    ((1 / 2 : ℝ) • (B *ᵥ Ypair μ m n)) ⬝ᵥ (x - (1 / 2 : ℝ) • (μ m + μ n))
      = (1 / 2) * (gpot B μ n x - gpot B μ m x) := by
  have hcross : (B *ᵥ μ n) ⬝ᵥ μ m = (B *ᵥ μ m) ⬝ᵥ μ n := dot_mulVec_comm hB _ _
  simp only [Ypair, gpot, Matrix.mulVec_sub, Matrix.mulVec_add,
    Matrix.smul_mulVec_assoc, sub_dotProduct, add_dotProduct, smul_dotProduct,
    dotProduct_sub, dotProduct_add, dotProduct_smul, smul_eq_mul]
  rw [hcross]; ring

/-- **Переходовість.** Якщо попарні статистики є різницями значень однієї
функції, то попарні рішення транзитивні: цикл неможливий. -/
theorem transitive_of_potential {lam : ι → ι → ℝ} {f : ι → ℝ}
    (h : ∀ m n, lam m n = f n - f m) {m n p : ι}
    (h₁ : 0 < lam m n) (h₂ : 0 < lam n p) : 0 < lam m p := by
  rw [h] at h₁ h₂ ⊢
  linarith

/-! ### Т1: контрприклад до переходовості

`μ⁰ = (0,0)`, `μ¹ = (2,0)`, `μ² = (1,2)`;
`C⁰ = I`, `C¹ = [[1,-1],[-1,2]]`, `C² = [[1,1],[1,3]]` — усі додатно визначені.
У точці `x = (1, ½)` попарні рішення утворюють цикл `1 ≻ 0`, `2 ≻ 1`, `0 ≻ 2`.

Коефіцієнти подано як розв'язки систем `F K = Y`, тому обертати матриці не треба. -/

theorem exists_intransitive :
    ∃ (μ₀ μ₁ μ₂ : Fin 2 → ℝ) (C₀ C₁ C₂ : Matrix (Fin 2) (Fin 2) ℝ)
      (K₀₁ K₁₂ K₀₂ x : Fin 2 → ℝ),
      (C₀ + C₁) *ᵥ K₀₁ = μ₁ - μ₀ ∧
      (C₁ + C₂) *ᵥ K₁₂ = μ₂ - μ₁ ∧
      (C₀ + C₂) *ᵥ K₀₂ = μ₂ - μ₀ ∧
      0 < K₀₁ ⬝ᵥ (x - (1 / 2 : ℝ) • (μ₀ + μ₁)) ∧
      0 < K₁₂ ⬝ᵥ (x - (1 / 2 : ℝ) • (μ₁ + μ₂)) ∧
      K₀₂ ⬝ᵥ (x - (1 / 2 : ℝ) • (μ₀ + μ₂)) < 0 := by
  refine ⟨![0, 0], ![2, 0], ![1, 2],
    !![1, 0; 0, 1], !![1, -1; -1, 2], !![1, 1; 1, 3],
    ![6 / 5, 2 / 5], ![-1 / 2, 2 / 5], ![2 / 7, 3 / 7], ![1, 1 / 2],
    ?_, ?_, ?_, ?_, ?_, ?_⟩
  · funext i; fin_cases i <;>
      simp [Matrix.mulVec, dotProduct, Fin.sum_univ_two] <;> norm_num
  · funext i; fin_cases i <;>
      simp [Matrix.mulVec, dotProduct, Fin.sum_univ_two] <;> norm_num
  · funext i; fin_cases i <;>
      simp [Matrix.mulVec, dotProduct, Fin.sum_univ_two] <;> norm_num
  · norm_num [dotProduct, Fin.sum_univ_two]
  · norm_num [dotProduct, Fin.sum_univ_two]
  · norm_num [dotProduct, Fin.sum_univ_two]

end GSA.Part2.MAry
