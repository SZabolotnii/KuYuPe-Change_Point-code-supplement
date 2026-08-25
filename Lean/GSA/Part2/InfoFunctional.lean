import GSA.Part2.BasisApprox

namespace GSA.Part2

open scoped BigOperators
open Filter

/-!
# Часткові суми Парсеваля по **ортонормованому** базису

## Що цей модуль доводить

Фіксуємо `HilbertBasis ℕ ℝ H` `b` і вектор `z : H`. Величина

  `J b z s = ∑_{i < s} ⟪b i, z⟫²`

— часткова сума Парсеваля. Три теореми нижче кажуть рівно ось що:

* `theorem2_a_upper_bound` — `J b z s ≤ ‖z‖²`;
* `theorem2_b_monotone`    — `s ↦ J b z s` монотонна;
* `theorem2_c_tendsto`     — `J b z s → ‖z‖²` при `s → ∞`.

Це чисті гільбертові факти, і нічого понад них тут не доведено.

## Що цей модуль НЕ доводить

1. **Базис ортонормований.** Відвантажені словники (`Φ_poly`, `Φ_log`, `Φ_frac`)
   ортонормованими не є. Неортогональний об'єкт рукопису `J = Yᵀ F⁻¹ Y` живе в
   `Kernel.lean`; його варіаційна характеризація — `BridgeGap.lean`.
2. **`z` — довільний вектор `H`.** Ніде в цьому файлі `z` не ототожнюється з
   логарифмом відношення правдоподібності, а `‖z‖²` — з жодною дивергенцією.
   Ані Фішер, ані χ², ані KL, ані Джеффріс тут не фігурують і не можуть
   фігурувати: у формулюваннях немає ні мір, ні щільностей.

## ⚠ Імена `theorem2_*` НЕ засвідчують Теорему 2 рукопису

Теорема 2 рукопису (arXiv:2605.23419, `paper/preprint_en/02_theory.tex`)
стверджує `J(s) ≤ D_J` і `J(s) → D_J`, де `D_J` — дивергенція Джеффріса. Це
статистичне ототожнення **хибне** і в цій формалізації не доводиться — див.
`erratum/ERRATUM_theorem2c_2026-08-23.md` і перелік `-- NOT FORMALISED` у
`GSA/Part2/BridgeGap.lean`. Правильна границя `J` рукопису за повного базису —
`(1/c)·2Δ/(2−Δ)`, де `Δ` — трикутна (Vincze–Le Cam) дискримінація.

Імена декларацій збережено заради неперервності аудиторського сліду. Читати їх
слід як «гільбертове твердження, з якого списано 2(a)/(b)/(c)», і ніколи як
«Теорема 2 рукопису Lean-verified».
-/

variable {H : Type*} [NormedAddCommGroup H] [InnerProductSpace ℝ H]

/-- Часткова сума Парсеваля `J b z s = ∑_{i<s} ⟪b i, z⟫²` для ортонормованого
базису `b` і довільного вектора `z : H`.

Це гільбертова величина — і тільки вона. Зокрема:

* це **не** інформаційний функціонал рукопису `J = Yᵀ F⁻¹ Y`; той визначено для
  довільного (можливо неортогонального) словника як `GSA.Part2.Jof` у
  `Kernel.lean`, а його варіаційну форму доведено в `BridgeGap.lean`;
* жодна теорема нижче не приписує `z`, `‖z‖²` чи границі `J b z s` статистичного
  змісту — ні інформації Фішера, ні χ², ні KL, ні дивергенції Джеффріса. -/
noncomputable def J (b : HilbertBasis ℕ ℝ H) (z : H) (s : ℕ) : ℝ :=
  (Finset.range s).sum fun i => (b.repr z i)^2

/-- Допоміжна рівність: повна сума квадратів коефіцієнтів дорівнює ‖z‖². -/
theorem tsum_repr_sq_eq_norm_sq (b : HilbertBasis ℕ ℝ H) (z : H) :
    (∑' i : ℕ, (b.repr z i)^2) = ‖z‖^2 := by
  have h := b.tsum_inner_mul_inner z z
  -- Приводимо терми до (b.repr z i)^2 і ⟪z,z⟫ = ‖z‖².
  have hterm :
      (fun i : ℕ => inner ℝ z (b i) * inner ℝ (b i) z) =
        (fun i : ℕ => (b.repr z i)^2) := by
    funext i
    -- b.repr z i = ⟪b i, z⟫, а для ℝ inner симетричний.
    simp [HilbertBasis.repr_apply_apply, real_inner_comm, pow_two]
  simpa [hterm, real_inner_self_eq_norm_sq] using h

/- Три властивості нижче — факти Парсеваля про `J b z s` та `‖z‖²`. І межа
   зверху, і границя — це `‖z‖²`, квадрат норми довільного вектора `H`, а не
   дивергенція. Див. заголовок модуля. -/
/-- (a) Обмеженість зверху квадратом норми: `J b z s ≤ ‖z‖²`.
Межа — саме `‖z‖²`, а не дивергенція Джеффріса. -/
theorem theorem2_a_upper_bound
    (b : HilbertBasis ℕ ℝ H) (z : H) (s : ℕ) :
    J b z s ≤ ‖z‖^2 := by
  have hsum : Summable (fun i : ℕ => (b.repr z i)^2) := by
    -- Витікає з сумовності скалярних добутків для ортонормованої системи.
    have hsum' := (b.orthonormal.inner_products_summable (x := z))
    refine hsum'.congr ?_
    intro i
    -- ‖⟪b i, z⟫‖^2 = (b.repr z i)^2 для ℝ.
    simp [HilbertBasis.repr_apply_apply, Real.norm_eq_abs, pow_two]
  have hle :
      (Finset.range s).sum (fun i => (b.repr z i)^2) ≤
        ∑' i : ℕ, (b.repr z i)^2 := by
    refine Summable.sum_le_tsum (s := Finset.range s) ?_ hsum
    intro i hi
    exact sq_nonneg (b.repr z i)
  simpa [J, tsum_repr_sq_eq_norm_sq] using hle

/-- (b) Монотонність `s ↦ J b z s`. -/
theorem theorem2_b_monotone
    (b : HilbertBasis ℕ ℝ H) (z : H) :
    Monotone (J b z) := by
  classical
  intro s t hst
  have hsubset : Finset.range s ⊆ Finset.range t := by
    intro i hi
    exact Finset.mem_range.mpr (lt_of_lt_of_le (Finset.mem_range.mp hi) hst)
  refine Finset.sum_le_sum_of_subset_of_nonneg hsubset ?_
  intro i hi hnot
  exact sq_nonneg (b.repr z i)

/-- (c) Збіжність `J b z s → ‖z‖²`.
Границя — квадрат норми `z`. Ототожнення цієї границі з дивергенцією
Джеффріса (Теорема 2(c) рукопису) хибне і тут не доводиться:
див. `GSA/Part2/BridgeGap.lean`. -/
theorem theorem2_c_tendsto
    (b : HilbertBasis ℕ ℝ H) (z : H) :
    Filter.Tendsto (J b z) Filter.atTop (nhds (‖z‖^2)) := by
  have hsum : Summable (fun i : ℕ => (b.repr z i)^2) := by
    have hsum' := (b.orthonormal.inner_products_summable (x := z))
    refine hsum'.congr ?_
    intro i
    simp [HilbertBasis.repr_apply_apply, Real.norm_eq_abs, pow_two]
  have ht :
      Tendsto (fun s : ℕ => (Finset.range s).sum fun i => (b.repr z i)^2)
        atTop (nhds (∑' i : ℕ, (b.repr z i)^2)) :=
    Summable.tendsto_sum_tsum_nat hsum
  simpa [J, tsum_repr_sq_eq_norm_sq] using ht

end GSA.Part2
