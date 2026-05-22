import Mathlib.Analysis.InnerProductSpace.GramSchmidtOrtho
import Mathlib.Analysis.InnerProductSpace.PiL2
import GSA.Part2.Setup

namespace GSA.Part2

open scoped BigOperators

variable {α : Type*} {n : ℕ}

abbrev EmpVec (n : ℕ) := PiLp 2 (fun _ : Fin n => ℝ)

/-!
# 2.2.4. Ортогоналізація (Грам–Шмідт)

У `Part2_review.md` описана "емпірична" ортогоналізація відносно скалярного добутку
⟪u,v⟫ₙ = (1/n) ∑ u(xᵢ) v(xᵢ).

У mathlib є загальний інструмент `gramSchmidt` для внутрішніх добутків.
Нижче фіксуємо емпіричний скалярний добуток як середнє значення по вибірці
та явний перехід до вектора в `Fin n → ℝ`.
-/

/-- Емпіричний вектор значень функції на вибірці. -/
noncomputable def empVec (x : Fin n → α) (u : α → ℝ) : EmpVec n :=
  WithLp.toLp 2 (fun i => u (x i))

/-- Емпіричний скалярний добуток на вибірці. -/
noncomputable def empInner (x : Fin n → α) (u v : α → ℝ) : ℝ :=
  (1 / (n:ℝ)) * ∑ i, u (x i) * v (x i)

lemma empInner_eq_empVec (x : Fin n → α) (u v : α → ℝ) :
    empInner x u v = (1 / (n:ℝ)) * ∑ i, empVec x u i * empVec x v i := by
  rfl

/-!
Працюємо з емпіричними векторами у просторі `Fin n → ℝ` і використовуємо
`gramSchmidtNormed` над стандартним inner product. Масштабування `1/√n`
дозволяє перенести ортонормованість на емпіричний скалярний добуток.
-/

/-- Емпіричний inner product на векторах. -/
noncomputable def empInnerVec (v w : EmpVec n) : ℝ :=
  (1 / (n:ℝ)) * inner ℝ v w

lemma empInner_eq_empInnerVec (x : Fin n → α) (u v : α → ℝ) :
    empInner x u v = empInnerVec (empVec x u) (empVec x v) := by
  unfold empInner empInnerVec empVec
  simp [PiLp.inner_apply, mul_comm]

/-- Масштабована емпірична репрезентація (щоб `empInner` збігався зі стандартним inner). -/
noncomputable def empVecScaled (x : Fin n → α) (u : α → ℝ) : EmpVec n :=
  (1 / Real.sqrt (n:ℝ)) • empVec x u

/-- Визначення ортонормованості відносно `empInnerVec`. -/
def empOrthonormalVec {ι : Type*} [DecidableEq ι] (v : ι → EmpVec n) : Prop :=
  ∀ i j, empInnerVec (n:=n) (v i) (v j) = (if i = j then 1 else 0)

lemma empInnerVec_smul_smul (hn : 0 < n) (v w : EmpVec n) :
    empInnerVec (n:=n) ((Real.sqrt (n:ℝ)) • v) ((Real.sqrt (n:ℝ)) • w) = inner ℝ v w := by
  have hn' : (n:ℝ) ≠ 0 := by
    exact_mod_cast (ne_of_gt hn)
  have hpos : 0 ≤ (n:ℝ) := by
    exact_mod_cast (Nat.cast_nonneg n)
  have hsqrt : Real.sqrt (n:ℝ) * Real.sqrt (n:ℝ) = (n:ℝ) := by
    simp [Real.mul_self_sqrt hpos]
  unfold empInnerVec
  have hinner :
      inner ℝ ((Real.sqrt (n:ℝ)) • v) ((Real.sqrt (n:ℝ)) • w) =
        (n:ℝ) * inner ℝ v w := by
    calc
      inner ℝ ((Real.sqrt (n:ℝ)) • v) ((Real.sqrt (n:ℝ)) • w)
          = (Real.sqrt (n:ℝ)) * (Real.sqrt (n:ℝ) * inner ℝ v w) := by
              simp [inner_smul_left, inner_smul_right]
      _ = (Real.sqrt (n:ℝ) * Real.sqrt (n:ℝ)) * inner ℝ v w := by
            simpa using
              (mul_assoc (Real.sqrt (n:ℝ)) (Real.sqrt (n:ℝ)) (inner ℝ v w)).symm
      _ = (n:ℝ) * inner ℝ v w := by
            simp [hsqrt]
  calc
    (1 / (n:ℝ)) * inner ℝ ((Real.sqrt (n:ℝ)) • v) ((Real.sqrt (n:ℝ)) • w)
        = (1 / (n:ℝ)) * (n:ℝ) * inner ℝ v w := by
          simp [hinner, mul_assoc]
    _ = (1:ℝ) * inner ℝ v w := by
          simp [div_eq_mul_inv, hn']
    _ = inner ℝ v w := by
          simp

/-- Емпірична ортонормалізація (векторна форма). -/
noncomputable def empGramSchmidt {ι : Type*} [LinearOrder ι] [LocallyFiniteOrderBot ι]
    [WellFoundedLT ι]
    (x : Fin n → α) (φ : ι → α → ℝ) : ι → EmpVec n :=
  fun i =>
    (Real.sqrt (n:ℝ)) •
      InnerProductSpace.gramSchmidtNormed ℝ (fun i => empVecScaled (n:=n) x (φ i)) i

theorem empGramSchmidt_orthonormal {ι : Type*} [LinearOrder ι] [LocallyFiniteOrderBot ι]
    [WellFoundedLT ι]
    (x : Fin n → α) (φ : ι → α → ℝ)
    (hlin : LinearIndependent ℝ (fun i => empVecScaled (n:=n) x (φ i)))
    (hn : 0 < n) :
    empOrthonormalVec (n:=n) (empGramSchmidt (n:=n) x φ) := by
  classical
  intro i j
  let g : ι → EmpVec n :=
    InnerProductSpace.gramSchmidtNormed ℝ (fun i => empVecScaled (n:=n) x (φ i))
  have horth : Orthonormal ℝ g :=
    InnerProductSpace.gramSchmidtNormed_orthonormal hlin
  by_cases h : i = j
  · subst h
    have hinner : inner ℝ (g i) (g i) = 1 := by
      have hnorm : ‖g i‖ = 1 := horth.1 i
      simp [hnorm]
    calc
      empInnerVec (n:=n) (empGramSchmidt (n:=n) x φ i) (empGramSchmidt (n:=n) x φ i)
          = inner ℝ (g i) (g i) := by
              simpa [empGramSchmidt, g] using
                (empInnerVec_smul_smul (n:=n) hn (v := g i) (w := g i))
      _ = 1 := hinner
      _ = (if i = i then 1 else 0) := by simp
  · have hinner : inner ℝ (g i) (g j) = 0 := horth.2 h
    calc
      empInnerVec (n:=n) (empGramSchmidt (n:=n) x φ i) (empGramSchmidt (n:=n) x φ j)
          = inner ℝ (g i) (g j) := by
              simpa [empGramSchmidt, g] using
                (empInnerVec_smul_smul (n:=n) hn (v := g i) (w := g j))
      _ = 0 := hinner
      _ = (if i = j then 1 else 0) := by simp [h]

/-! Спеціалізація на `Fin s` (умови на тип індексів виконуються автоматично). -/

noncomputable def empGramSchmidtFin {s : ℕ}
    (x : Fin n → α) (φ : Fin s → α → ℝ) : Fin s → EmpVec n :=
  empGramSchmidt (n:=n) x φ

theorem empGramSchmidtFin_orthonormal {s : ℕ}
    (x : Fin n → α) (φ : Fin s → α → ℝ)
    (hlin : LinearIndependent ℝ (fun i => empVecScaled (n:=n) x (φ i)))
    (hn : 0 < n) :
    empOrthonormalVec (n:=n) (empGramSchmidtFin (n:=n) x φ) := by
  simpa [empGramSchmidtFin] using
    (empGramSchmidt_orthonormal (n:=n) x φ hlin hn)

end GSA.Part2
