/-
GSA.Part2.Formalization.lean

Фрагменти формальної верифікації тверджень із Part2_review.md:
- Теорема 1 (КУ1), крок 3: (f1-f0)/(f1+f0) > 0 ↔ f1/f0 > 1
- Теорема 4(a): L²-збіжність часткових сум розкладу в гільбертовому просторі
- Узгодження коефіцієнтів із inner product
-/

import Mathlib

open scoped BigOperators Topology
open Filter

namespace Part2

/- =========================
   Теорема 1 (крок 3): алгебра
   ========================= -/

-- У файлі: ψ_opt(x) = (f1-f0)/(f1+f0), і ψ_opt(x) > 0 ↔ f1/f0 > 1.
-- Доводимо як чисту лему над ℝ за умов f0>0, f1≥0 (щільності невід'ємні; f0>0 для ділення).
lemma ku1_step3 {f0 f1 : ℝ} (hf0 : 0 < f0) (hf1 : 0 ≤ f1) :
    ((f1 - f0) / (f1 + f0) > 0) ↔ (f1 / f0 > 1) := by
  have hden : 0 < f1 + f0 := add_pos_of_nonneg_of_pos hf1 hf0
  have hne : (f1 + f0) ≠ 0 := ne_of_gt hden

  have hfrac : ((f1 - f0) / (f1 + f0) > 0) ↔ (0 < f1 - f0) := by
    constructor
    · intro h
      -- множимо нерівність на додатний знаменник
      have : 0 < ((f1 - f0) / (f1 + f0)) * (f1 + f0) := by
        simpa [zero_mul] using (mul_lt_mul_of_pos_right h hden)
      -- скорочуємо (…/(f1+f0))*(f1+f0) = (…)
      simpa [div_eq_mul_inv, hne, mul_assoc] using this
    · intro h
      exact div_pos h hden

  -- 0 < f1 - f0 ↔ f0 < f1 (лема sub_pos)
  -- f0 < f1 ↔ 1 < f1/f0 (домноження на додатний f0)
  have : ((f1 - f0) / (f1 + f0) > 0) ↔ (1 < f1 / f0) := by
    calc
      ((f1 - f0) / (f1 + f0) > 0) ↔ (0 < f1 - f0) := hfrac
      _ ↔ (f0 < f1) := by simp [sub_pos]
      _ ↔ (1 < f1 / f0) := by
            -- one_lt_div hf0 : 1 < f1/f0 ↔ f0 < f1
            simpa using (one_lt_div (a := f1) (b := f0) hf0).symm

  -- Переписуємо “1 < …” як “… > 1” (це те саме визначення)
  simpa using this

/- =======================================
   Теорема 4(a): збіжність у гільбертовому просторі
   ======================================= -/

section HilbertApprox

variable {E : Type*} [NormedAddCommGroup E] [InnerProductSpace ℝ E]

-- Нехай b : HilbertBasis ℕ ℝ E — повний ортонормований базис (mathlib-структура).
variable (b : HilbertBasis ℕ ℝ E) (z : E)

-- Часткова сума (0..s-1). Це Lean-версія Λ^(s)(x) = Σ_{i=1..s} k_i φ_i(x),
-- тільки з нульовою індексацією.
noncomputable def approx (s : ℕ) : E :=
  (Finset.range s).sum fun i => (b.repr z i) • (b i)

-- “Рісс–Фішер” / збіжність розкладу: HasSum … = z (готова теорема mathlib).
-- Див. HilbertBasis.hasSum_repr.
theorem approx_tendsto :
    Tendsto (fun s : ℕ => approx b z s) atTop (𝓝 z) := by
  have hsum : HasSum (fun i : ℕ => (b.repr z i) • (b i)) z := by
    simpa using (HilbertBasis.hasSum_repr (b := b) z)

  have hsumm : Summable (fun i : ℕ => (b.repr z i) • (b i)) :=
    hsum.summable

  have htsum : (∑' i : ℕ, (b.repr z i) • (b i)) = z :=
    hsum.tsum_eq

  -- Summable.tendsto_sum_tsum_nat: часткові суми → tsum
  have ht :
      Tendsto (fun s : ℕ => (Finset.range s).sum fun i => (b.repr z i) • (b i))
        atTop (𝓝 (∑' i : ℕ, (b.repr z i) • (b i))) :=
    Summable.tendsto_sum_tsum_nat hsumm

  simpa [approx, htsum] using ht

-- Узгодження коефіцієнта з inner product:
-- (b.repr z i) = ⟪b i, z⟫  (готова лема repr_apply_apply)
theorem coeff_eq_inner (i : ℕ) :
    (b.repr z i) = inner ℝ (b i) z := by
  simpa using (HilbertBasis.repr_apply_apply (b := b) z i)

end HilbertApprox

end Part2
