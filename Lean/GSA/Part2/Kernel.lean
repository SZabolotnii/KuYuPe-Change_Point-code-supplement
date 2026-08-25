import Mathlib.Analysis.InnerProductSpace.Basic
import Mathlib.Analysis.Matrix.PosDef
import Mathlib.LinearAlgebra.Matrix.NonsingularInverse

open scoped BigOperators
open Finset

namespace GSA.Part2

/-!
# Theorem kernel системи `F K = B`

Модуль формалізує алгебраїчне ядро узагальненої стохастичної апроксимації —
ту частину §2.4, що **не залежить** ні від DGP, ні від вибору словника, ні від
скінченної вибірки.

## Чому окремо від `InfoFunctional.lean`

`InfoFunctional.lean` доводить властивості `J(s)` в **ортонормованому** базисі
(`HilbertBasis`), де `J(s)` — часткова сума квадратів коефіцієнтів. Це не той
об'єкт, з яким працює рукопис: відвантажені словники (`Φ_poly`, `Φ_log`,
`Φ_frac`) **не ортогональні**, і `J = Bᵀ F⁻¹ B` з невиродженою, але сильно
неодиничною матрицею Грама `F`. Саме неортогональний випадок породжує всі
питання про обумовленість (§2.5), і саме він тут формалізується.

## Ключове рішення: варіаційне означення

`J` означується не формулою `Bᵀ F⁻¹ B`, а через **ризик апроксимації**
`R(K) = ‖q − Σ Kᵢ vᵢ‖²`:

  `J(K) = ‖q‖² − R(K)`.

Формула `Bᵀ F⁻¹ B` виникає як теорема (`Jof_eq_dot_of_normalEq`) для розв'язку
нормальних рівнянь. Виграш від такого порядку означень у тому, що інформаційна
межа, монотонність за вкладеними словниками та інваріантність відносно заміни
базису стають наслідками **однієї** геометричної тотожності
(`normalEq_iff_orthogonal`: залишок ортогональний до span), а не окремих
матричних викладок.

## Відповідність пунктам програми формалізації

| # | Твердження | Теорема |
|---|---|---|
| 1 | існування та єдиність розв'язку `F K = B` | `existsUnique_normalEq_of_posDef` |
| 2 | `K* = F⁻¹B` мінімізує квадратичний ризик | `normalEq_isMinimizer` |
| 3 | residual identity | `risk_eq_risk_add_normSq` |
| 4 | інформаційна межа `J ≤ E[q²]` | `Jof_le_normSq` |
| 5 | монотонність за вкладеними словниками | `Jof_init_le_Jof_of_isMinimizer` |
| 6 | Schur-приріст від однієї фічі | `Jof_snoc_general` (частковий випадок: `Jof_snoc_orthogonal`) |
| 7 | інваріантність span при заміні базису | `risk_eq_of_range_eq`, `range_combo_eq_of_isUnit` |
| 8 | межа ортогоналізації (те саме `J`, інша `F`) | `isMinimizer_transfer` |
| 9 | вироджена `F`: `J` визначена коректно | `risk_eq_of_normalEq_of_normalEq` |

Пункт 10 програми (виродження PATP при `α = 1/2`) стосується іншого рукопису
й тут не формалізується.

## ⚠ Перетин із `Ku-Projection-Framework` (виявлено 2026-08-13, після написання)

Трек `Ku-Projection-Framework` (подано в Archives of Control Sciences, OJS 2649)
має власну формалізацію `lean/GEP/Core.lean`, 162 рядки, у **тій самій**
загальності (словник не ортогональний і не обов'язково лінійно незалежний):

| Там | Тут | Хто сильніший |
|---|---|---|
| `normalSystem_iff_starProjection` | `normalEq_iff_orthogonal` | там — через `Submodule.starProjection`, тут — елементарно; еквівалентні |
| `energy_eq_dotProduct` | `Jof_eq_dot_of_normalEq` | еквівалентні |
| `energy_le` | `Jof_le_normSq` | еквівалентні |
| `energy_mono` | `Jof_init_le_Jof_of_isMinimizer` | **там сильніше**: довільне вкладення span, а не лише префікс |
| — (наслідок `energy_eq_dotProduct`) | `risk_eq_of_normalEq_of_normalEq` | **там сильніше**: енергія від вибору розв'язку не залежить за побудовою |

Важливо: `GEP/Core.lean` прямо називає відповідність `FK=Y ↔ L²-проєкція`
питанням, **поставленим як відкрите в цьому самому GSA-рукописі**. Тобто на
відкрите питання parent'а вже відповів сусідній трек — і revision має це
процитувати, а не переоткривати.

**Що додає цей модуль понад `GEP/Core.lean`:**

1. **варіаційний шар** — `risk`, `Jof`, `IsMinimizer` і residual identity
   `risk_eq_risk_add_normSq`. GEP працює з нормами проєкцій і поняття
   квадратичного ризику не має взагалі; пункти 2–3 §12.2 там не виражені;
2. **`Jof_snoc_general`** — Schur-приріст від однієї фічі. У GEP нічого
   подібного немає, і саме він відповідає на питання «скільки дає ще одна
   фіча»;
3. **`existsUnique_normalEq_of_posDef`** — єдиність за `PosDef`. GEP цей факт
   лише згадує в docstring, не доводить.

Решта — незалежний перевивід. Це прийнятно (крос-репозиторна Lean-залежність
між приватними треками непрактична) і навіть корисне як незалежна перевірка
двома різними стилями, але **видавати це за нове не можна**.

**Чого це ядро не доводить:** скінченно-вибіркової ефективності, коректності
плагін-оцінок кумулянтів, адекватності розподілу, новизни. Формалізація та
експериментальна перевірка лишаються двома окремими шарами доказів.
-/

variable {H : Type*} [NormedAddCommGroup H] [InnerProductSpace ℝ H]
variable {s : ℕ}

/-! ## 1. Об'єкти -/

/-- Матриця Грама словника: `Fᵢⱼ = ⟪vᵢ, vⱼ⟫`.

Статистично `vᵢ` — центрована фіча `φ̃ᵢ` як елемент `L²(μ)`, тож
`Fᵢⱼ = E[φ̃ᵢ φ̃ⱼ]` — саме матриця з §2.4. -/
noncomputable def gram (v : Fin s → H) : Matrix (Fin s) (Fin s) ℝ :=
  fun i j => inner ℝ (v i) (v j)

/-- Права частина: `Bᵢ = ⟪vᵢ, q⟫`; статистично `Bᵢ = E[φ̃ᵢ q]`, де `q` — score
напряму зміни. -/
noncomputable def bvec (v : Fin s → H) (q : H) : Fin s → ℝ :=
  fun i => inner ℝ (v i) q

/-- Наближення `Σ Kᵢ vᵢ` — GSA-increment із коефіцієнтами `K`. -/
noncomputable def combo (v : Fin s → H) (K : Fin s → ℝ) : H :=
  ∑ i, K i • v i

/-- Ризик апроксимації `R(K) = ‖q − Σ Kᵢ vᵢ‖²`. -/
noncomputable def risk (v : Fin s → H) (q : H) (K : Fin s → ℝ) : ℝ :=
  ‖q - combo v K‖ ^ 2

/-- Захоплена інформація `J(K) = ‖q‖² − R(K)`.

Це варіаційне означення; зв'язок із `Bᵀ F⁻¹ B` — `Jof_eq_dot_of_normalEq`. -/
noncomputable def Jof (v : Fin s → H) (q : H) (K : Fin s → ℝ) : ℝ :=
  ‖q‖ ^ 2 - risk v q K

/-- Нормальні рівняння `F K = B`. -/
def NormalEq (v : Fin s → H) (q : H) (K : Fin s → ℝ) : Prop :=
  ∀ i, ∑ j, gram v i j * K j = bvec v q i

/-- `K` мінімізує ризик серед усіх коефіцієнтів. -/
def IsMinimizer (v : Fin s → H) (q : H) (K : Fin s → ℝ) : Prop :=
  ∀ L, risk v q K ≤ risk v q L

/-! ## 2. Лінійність наближення і геометрія нормальних рівнянь -/

lemma gram_symm (v : Fin s → H) (i j : Fin s) : gram v i j = gram v j i := by
  simp [gram, real_inner_comm]

lemma combo_sub (v : Fin s → H) (K L : Fin s → ℝ) :
    combo v (K - L) = combo v K - combo v L := by
  simp [combo, sub_smul, Finset.sum_sub_distrib]

lemma combo_add (v : Fin s → H) (K L : Fin s → ℝ) :
    combo v (K + L) = combo v K + combo v L := by
  simp [combo, add_smul, Finset.sum_add_distrib]

lemma combo_smul (v : Fin s → H) (c : ℝ) (K : Fin s → ℝ) :
    combo v (c • K) = c • combo v K := by
  simp [combo, Finset.smul_sum, smul_smul]

lemma inner_combo_left (v : Fin s → H) (K : Fin s → ℝ) (x : H) :
    inner ℝ (combo v K) x = ∑ i, K i * inner ℝ (v i) x := by
  simp [combo, sum_inner, real_inner_smul_left]

lemma inner_v_combo (v : Fin s → H) (K : Fin s → ℝ) (i : Fin s) :
    inner ℝ (v i) (combo v K) = ∑ j, gram v i j * K j := by
  simp only [combo, inner_sum, real_inner_smul_right, gram]
  exact Finset.sum_congr rfl fun j _ => by ring

/-- **Геометричний зміст нормальних рівнянь:** `F K = B` рівносильне тому, що
залишок `q − Σ Kᵢ vᵢ` ортогональний до кожної фічі, тобто до всього span.

Це та єдина тотожність, з якої далі виводиться все інше. -/
theorem normalEq_iff_orthogonal (v : Fin s → H) (q : H) (K : Fin s → ℝ) :
    NormalEq v q K ↔ ∀ i, inner ℝ (v i) (q - combo v K) = (0 : ℝ) := by
  simp only [NormalEq, inner_sub_right, inner_v_combo, bvec, sub_eq_zero]
  exact ⟨fun h i => (h i).symm, fun h i => (h i).symm⟩

/-- Залишок оптимального наближення ортогональний до будь-якого елемента span. -/
lemma inner_combo_resid_eq_zero (v : Fin s → H) (q : H) {K : Fin s → ℝ}
    (hK : NormalEq v q K) (L : Fin s → ℝ) :
    inner ℝ (combo v L) (q - combo v K) = (0 : ℝ) := by
  rw [inner_combo_left]
  refine Finset.sum_eq_zero fun i _ => ?_
  rw [(normalEq_iff_orthogonal v q K).mp hK i, mul_zero]

/-! ## 3. Residual identity, оптимальність, коректність при виродженні -/

/-- **(3) Residual identity.** Для розв'язку нормальних рівнянь `K*`

  `R(K) = R(K*) + ‖Σ (K − K*)ᵢ vᵢ‖²`,

де другий доданок — це рівно квадратична форма `(K − K*)ᵀ F (K − K*)`.
Невиродженість `F` не потрібна. -/
theorem risk_eq_risk_add_normSq (v : Fin s → H) (q : H) {Ks : Fin s → ℝ}
    (hKs : NormalEq v q Ks) (K : Fin s → ℝ) :
    risk v q K = risk v q Ks + ‖combo v (K - Ks)‖ ^ 2 := by
  have hdecomp : q - combo v K = (q - combo v Ks) - combo v (K - Ks) := by
    rw [combo_sub]; abel
  have horth : inner ℝ (q - combo v Ks) (combo v (K - Ks)) = (0 : ℝ) := by
    rw [real_inner_comm]
    exact inner_combo_resid_eq_zero v q hKs _
  rw [risk, hdecomp, norm_sub_sq_real, horth, risk]
  ring

/-- Квадратична форма Грама у явному вигляді — для звірки з матричним записом
рукопису. -/
theorem normSq_combo_eq_quadForm (v : Fin s → H) (K : Fin s → ℝ) :
    ‖combo v K‖ ^ 2 = ∑ i, ∑ j, K i * K j * gram v i j := by
  rw [← real_inner_self_eq_norm_sq, inner_combo_left]
  refine Finset.sum_congr rfl fun i _ => ?_
  rw [inner_v_combo, Finset.mul_sum]
  exact Finset.sum_congr rfl fun j _ => by ring

/-- **(2) Оптимальність.** Розв'язок нормальних рівнянь мінімізує ризик.
Додатної визначеності не потрібно: `F` — матриця Грама, тож форма невід'ємна
автоматично. -/
theorem normalEq_isMinimizer (v : Fin s → H) (q : H) {Ks : Fin s → ℝ}
    (hKs : NormalEq v q Ks) : IsMinimizer v q Ks := by
  intro L
  rw [risk_eq_risk_add_normSq v q hKs L]
  have : (0 : ℝ) ≤ ‖combo v (L - Ks)‖ ^ 2 := sq_nonneg _
  linarith

/-- **(9) Коректність `J` при виродженій `F`.** Два різні розв'язки нормальних
рівнянь дають той самий ризик, тож `J` визначена однозначно навіть для лінійно
залежного словника, коли `K*` не єдиний. -/
theorem risk_eq_of_normalEq_of_normalEq (v : Fin s → H) (q : H) {K₁ K₂ : Fin s → ℝ}
    (h₁ : NormalEq v q K₁) (h₂ : NormalEq v q K₂) :
    risk v q K₁ = risk v q K₂ := by
  have ha := risk_eq_risk_add_normSq v q h₂ K₁
  have hb := risk_eq_risk_add_normSq v q h₁ K₂
  have hx : (0 : ℝ) ≤ ‖combo v (K₁ - K₂)‖ ^ 2 := sq_nonneg _
  have hy : (0 : ℝ) ≤ ‖combo v (K₂ - K₁)‖ ^ 2 := sq_nonneg _
  linarith

/-! ## 4. Інформаційна межа і формула `Bᵀ F⁻¹ B` -/

/-- **(4) Інформаційна межа.** `J ≤ ‖q‖²` — для будь-яких коефіцієнтів, не лише
оптимальних. Статистично: жоден словник не захоплює більше інформації, ніж є в
напрямі зміни. -/
theorem Jof_le_normSq (v : Fin s → H) (q : H) (K : Fin s → ℝ) :
    Jof v q K ≤ ‖q‖ ^ 2 := by
  have : (0 : ℝ) ≤ risk v q K := by rw [risk]; positivity
  simp only [Jof]
  linarith

/-- Захоплена інформація дорівнює квадрату норми проєкції. Це та сама
Parseval-картина, що в `InfoFunctional.lean`, але **без** припущення
ортонормованості словника. -/
theorem Jof_eq_normSq_combo (v : Fin s → H) (q : H) {Ks : Fin s → ℝ}
    (hKs : NormalEq v q Ks) : Jof v q Ks = ‖combo v Ks‖ ^ 2 := by
  have horth : inner ℝ (combo v Ks) (q - combo v Ks) = (0 : ℝ) :=
    inner_combo_resid_eq_zero v q hKs Ks
  rw [inner_sub_right, real_inner_self_eq_norm_sq, sub_eq_zero] at horth
  have hqc : inner ℝ q (combo v Ks) = ‖combo v Ks‖ ^ 2 := by
    rw [real_inner_comm]; exact horth
  rw [Jof, risk, norm_sub_sq_real, hqc]
  ring

/-- Для розв'язку нормальних рівнянь `J = Bᵀ K`, тобто `Bᵀ F⁻¹ B` за
невиродженої `F`. -/
theorem Jof_eq_dot_of_normalEq (v : Fin s → H) (q : H) {Ks : Fin s → ℝ}
    (hKs : NormalEq v q Ks) :
    Jof v q Ks = ∑ i, Ks i * bvec v q i := by
  have horth : inner ℝ (combo v Ks) (q - combo v Ks) = (0 : ℝ) :=
    inner_combo_resid_eq_zero v q hKs Ks
  rw [inner_sub_right, real_inner_self_eq_norm_sq, sub_eq_zero] at horth
  rw [Jof_eq_normSq_combo v q hKs, ← horth, inner_combo_left]
  rfl

/-! ## 5. Існування та єдиність -/

/-- **(1) Існування та єдиність.** Якщо матриця Грама додатно визначена
(словник лінійно незалежний), нормальні рівняння мають рівно один розв'язок. -/
theorem existsUnique_normalEq_of_posDef (v : Fin s → H) (q : H)
    (hF : (gram v).PosDef) : ∃! K : Fin s → ℝ, NormalEq v q K := by
  classical
  have hunit : IsUnit (gram v).det :=
    isUnit_iff_ne_zero.mpr (ne_of_gt (Matrix.PosDef.det_pos hF))
  have _ := Matrix.invertibleOfIsUnitDet (A := gram v) hunit
  have hmv : ∀ K : Fin s → ℝ, NormalEq v q K ↔ (gram v).mulVec K = bvec v q := by
    intro K
    constructor
    · intro h; funext i; simpa [Matrix.mulVec, dotProduct] using h i
    · intro h i; simpa [Matrix.mulVec, dotProduct] using congrFun h i
  refine ⟨(gram v)⁻¹.mulVec (bvec v q), ?_, ?_⟩
  · show NormalEq v q ((gram v)⁻¹.mulVec (bvec v q))
    rw [hmv, Matrix.mulVec_mulVec, Matrix.mul_inv_of_invertible, Matrix.one_mulVec]
  · intro K hK
    show K = (gram v)⁻¹.mulVec (bvec v q)
    rw [hmv] at hK
    rw [← hK, Matrix.mulVec_mulVec, Matrix.inv_mul_of_invertible, Matrix.one_mulVec]

/-- Оптимальний `K*` існує за додатної визначеності `F`. -/
theorem exists_isMinimizer_of_posDef (v : Fin s → H) (q : H)
    (hF : (gram v).PosDef) :
    ∃ K : Fin s → ℝ, IsMinimizer v q K ∧ NormalEq v q K := by
  obtain ⟨K, hK, _⟩ := existsUnique_normalEq_of_posDef v q hF
  exact ⟨K, normalEq_isMinimizer v q hK, hK⟩

/-! ## 6. Монотонність за вкладеними словниками -/

/-- Доповнення коефіцієнта нулем не змінює наближення. -/
lemma combo_snoc_zero (v : Fin (s + 1) → H) (K : Fin s → ℝ) :
    combo v (Fin.snoc K (0 : ℝ)) = combo (Fin.init v) K := by
  simp [combo, Fin.sum_univ_castSucc, Fin.init]

lemma risk_snoc_zero (v : Fin (s + 1) → H) (q : H) (K : Fin s → ℝ) :
    risk v q (Fin.snoc K (0 : ℝ)) = risk (Fin.init v) q K := by
  simp [risk, combo_snoc_zero]

/-- **(5) Монотонність за вкладеними словниками.** Додавання фічі не може
зменшити захоплену інформацію.

Це те твердження, яким числовий артефакт відрізняється від результату:
спостережене падіння `κ` зі зростанням `s` у подвійній точності **не може** бути
властивістю словника (див. `reports/basis_information_ceiling_20260813.md`,
§6.1). -/
theorem Jof_init_le_Jof_of_isMinimizer (v : Fin (s + 1) → H) (q : H)
    {Kbig : Fin (s + 1) → ℝ} (hopt : IsMinimizer v q Kbig) (K : Fin s → ℝ) :
    Jof (Fin.init v) q K ≤ Jof v q Kbig := by
  have h1 : risk v q Kbig ≤ risk v q (Fin.snoc K (0 : ℝ)) := hopt _
  rw [risk_snoc_zero] at h1
  simp only [Jof]
  linarith

/-! ## 7. Інваріантність відносно заміни базису того самого span -/

/-- **`J` залежить лише від span, а не від словника.** Якщо два словники —
можливо, різної довжини — породжують ту саму множину досяжних наближень, їхні
оптимальні ризики збігаються.

Це найзагальніша форма пунктів 7–8: вона не потребує ні матриці переходу, ні
однакової кількості фіч. Ортогоналізація, перенормування і будь-яка інша заміна
базису того самого span потрапляють сюди як окремі випадки. -/
theorem risk_eq_of_range_eq {m n : ℕ} {v : Fin m → H} {w : Fin n → H} (q : H)
    (hrange : Set.range (combo v) = Set.range (combo w))
    {Kv : Fin m → ℝ} {Kw : Fin n → ℝ}
    (hv : IsMinimizer v q Kv) (hw : IsMinimizer w q Kw) :
    risk v q Kv = risk w q Kw := by
  have h₁ : risk w q Kw ≤ risk v q Kv := by
    obtain ⟨L, hL⟩ : combo v Kv ∈ Set.range (combo w) := by
      rw [← hrange]; exact ⟨Kv, rfl⟩
    calc risk w q Kw ≤ risk w q L := hw L
      _ = risk v q Kv := by rw [risk, risk, hL]
  have h₂ : risk v q Kv ≤ risk w q Kw := by
    obtain ⟨L, hL⟩ : combo w Kw ∈ Set.range (combo v) := by
      rw [hrange]; exact ⟨Kw, rfl⟩
    calc risk v q Kv ≤ risk v q L := hv L
      _ = risk w q Kw := by rw [risk, risk, hL]
  linarith

/-- Заміна словника `wᵢ = Σⱼ Aᵢⱼ vⱼ` переносить коефіцієнти через `Aᵀ`. -/
lemma combo_transform (A : Matrix (Fin s) (Fin s) ℝ) (v : Fin s → H) (K : Fin s → ℝ) :
    combo (fun i => ∑ j, A i j • v j) K = combo v (Matrix.vecMul K A) := by
  simp only [combo, Matrix.vecMul, dotProduct, Finset.smul_sum, Finset.sum_smul]
  rw [Finset.sum_comm]
  exact Finset.sum_congr rfl fun j _ => Finset.sum_congr rfl fun i _ => by
    rw [smul_smul]

/-- **(7) Інваріантність span.** Невироджена заміна базису не змінює множину
досяжних наближень. Звідси однакові оптимальні ризики, тобто однакове `J` —
попри те, що матриці Грама (і їхня обумовленість) різні. -/
theorem range_combo_eq_of_isUnit (A : Matrix (Fin s) (Fin s) ℝ) (v : Fin s → H)
    (hA : IsUnit A.det) :
    Set.range (combo (fun i => ∑ j, A i j • v j)) = Set.range (combo v) := by
  classical
  have _ := Matrix.invertibleOfIsUnitDet (A := A) hA
  apply Set.Subset.antisymm
  · rintro x ⟨K, rfl⟩
    exact ⟨Matrix.vecMul K A, (combo_transform A v K).symm⟩
  · rintro x ⟨L, rfl⟩
    refine ⟨Matrix.vecMul L A⁻¹, ?_⟩
    rw [combo_transform, Matrix.vecMul_vecMul, Matrix.inv_mul_of_invertible,
      Matrix.vecMul_one]

/-- **(8) Межа ортогоналізації.** Мінімізатор переноситься між словниками того
самого span зі збереженням ризику. Це і є твердження «ортогоналізація змінює
обумовленість, але не популяційний оптимум». -/
theorem isMinimizer_transfer (A : Matrix (Fin s) (Fin s) ℝ) (v : Fin s → H) (q : H)
    (hA : IsUnit A.det) {K : Fin s → ℝ} (hK : IsMinimizer v q K) :
    IsMinimizer (fun i => ∑ j, A i j • v j) q (Matrix.vecMul K A⁻¹) ∧
      risk (fun i => ∑ j, A i j • v j) q (Matrix.vecMul K A⁻¹) = risk v q K := by
  classical
  have _ := Matrix.invertibleOfIsUnitDet (A := A) hA
  have hback : combo (fun i => ∑ j, A i j • v j) (Matrix.vecMul K A⁻¹) = combo v K := by
    rw [combo_transform, Matrix.vecMul_vecMul, Matrix.inv_mul_of_invertible,
      Matrix.vecMul_one]
  have hrisk : risk (fun i => ∑ j, A i j • v j) q (Matrix.vecMul K A⁻¹) = risk v q K := by
    simp [risk, hback]
  refine ⟨?_, hrisk⟩
  intro L
  rw [hrisk]
  have hL : risk (fun i => ∑ j, A i j • v j) q L = risk v q (Matrix.vecMul L A) := by
    simp [risk, combo_transform]
  rw [hL]
  exact hK _

/-! ## 8. Приріст від ортогональної фічі -/

lemma combo_snoc (v : Fin s → H) (u : H) (K : Fin s → ℝ) (c : ℝ) :
    combo (Fin.snoc v u) (Fin.snoc K c) = combo v K + c • u := by
  simp [combo, Fin.sum_univ_castSucc]

lemma combo_snoc_apply (v : Fin s → H) (u : H) (L : Fin (s + 1) → ℝ) :
    combo (Fin.snoc v u) L = combo v (Fin.init L) + (L (Fin.last s)) • u := by
  simp [combo, Fin.sum_univ_castSucc, Fin.init]

/-- **(6) Ортогональний приріст.** Якщо нова фіча ортогональна до всіх наявних,
приріст інформації дорівнює `⟪u, q⟫² / ‖u‖²` — рівно одновимірна проєкція на неї.
У парі з уже формалізованою ортогоналізацією (`GramSchmidt.lean`) це дає покрокову
форму приросту. -/
theorem Jof_snoc_orthogonal (v : Fin s → H) (q u : H) {K : Fin s → ℝ}
    (hK : NormalEq v q K) (horth : ∀ i, inner ℝ (v i) u = (0 : ℝ)) (hu : u ≠ 0) :
    NormalEq (Fin.snoc v u) q (Fin.snoc K (inner ℝ u q / ‖u‖ ^ 2)) ∧
      Jof (Fin.snoc v u) q (Fin.snoc K (inner ℝ u q / ‖u‖ ^ 2))
        = Jof v q K + (inner ℝ u q : ℝ) ^ 2 / ‖u‖ ^ 2 := by
  classical
  set c : ℝ := inner ℝ u q / ‖u‖ ^ 2 with hc
  have hun : (0 : ℝ) < ‖u‖ ^ 2 := by positivity
  have hcombo := combo_snoc v u K c
  -- залишок нового наближення
  have hres : q - combo (Fin.snoc v u) (Fin.snoc K c) = (q - combo v K) - c • u := by
    rw [hcombo]; abel
  -- фіча `u` ортогональна до старого наближення
  have hu_combo : inner ℝ u (combo v K) = (0 : ℝ) := by
    rw [combo, inner_sum]
    refine Finset.sum_eq_zero fun i _ => ?_
    rw [real_inner_smul_right, real_inner_comm, horth i, mul_zero]
  have hcombo_u : inner ℝ (combo v K) u = (0 : ℝ) := by
    rw [real_inner_comm]; exact hu_combo
  have hnormal : NormalEq (Fin.snoc v u) q (Fin.snoc K c) := by
    rw [normalEq_iff_orthogonal]
    intro i
    rw [hres]
    rcases Fin.eq_castSucc_or_eq_last i with ⟨i₀, rfl⟩ | rfl
    · rw [Fin.snoc_castSucc, inner_sub_right, real_inner_smul_right, horth i₀,
        (normalEq_iff_orthogonal v q K).mp hK i₀]
      ring
    · rw [Fin.snoc_last, inner_sub_right, inner_sub_right, hu_combo,
        real_inner_smul_right, real_inner_self_eq_norm_sq, hc]
      field_simp
      ring
  refine ⟨hnormal, ?_⟩
  rw [Jof_eq_normSq_combo _ q hnormal, Jof_eq_normSq_combo v q hK, hcombo]
  rw [norm_add_sq_real, real_inner_smul_right, hcombo_u, norm_smul,
    Real.norm_eq_abs, mul_pow, sq_abs, hc]
  field_simp
  ring

/-- Приєднання фічі `u` і приєднання її залишку `r = u − proj u` дають той самий
span, тож і той самий оптимальний ризик. -/
lemma range_combo_snoc_resid (v : Fin s → H) (u : H) (Ku : Fin s → ℝ) :
    Set.range (combo (Fin.snoc v u)) =
      Set.range (combo (Fin.snoc v (u - combo v Ku))) := by
  apply Set.Subset.antisymm
  · rintro x ⟨L, rfl⟩
    refine ⟨Fin.snoc (Fin.init L + (L (Fin.last s)) • Ku) (L (Fin.last s)), ?_⟩
    rw [combo_snoc, combo_add, combo_smul, combo_snoc_apply, smul_sub]
    abel
  · rintro x ⟨L, rfl⟩
    refine ⟨Fin.snoc (Fin.init L - (L (Fin.last s)) • Ku) (L (Fin.last s)), ?_⟩
    rw [combo_snoc, combo_sub, combo_smul, combo_snoc_apply, smul_sub]
    abel

/-- **(6) Загальний приріст (Schur-complement у геометричній формі).**

Додавання довільної фічі `u` збільшує захоплену інформацію рівно на
`⟪r, q⟫² / ‖r‖²`, де `r = u − Σ (Ku)ⱼ vⱼ` — частина `u`, ортогональна до
наявного словника.

Це кількісна відповідь на питання «скільки дає ще одна фіча»: внесок визначає не
сама `u`, а **її новизна відносно наявного словника**. Зокрема, майже колінеарна
фіча (мале `‖r‖`) не додає інформації, зате псує обумовленість `F` — рівно той
компроміс, який рукопис обговорює в §2.5. -/
theorem Jof_snoc_general (v : Fin s → H) (q u : H) {K Ku : Fin s → ℝ}
    (hK : NormalEq v q K) (hKu : NormalEq v u Ku)
    (hr : u - combo v Ku ≠ 0)
    {Kbig : Fin (s + 1) → ℝ} (hbig : IsMinimizer (Fin.snoc v u) q Kbig) :
    Jof (Fin.snoc v u) q Kbig
      = Jof v q K + (inner ℝ (u - combo v Ku) q : ℝ) ^ 2 / ‖u - combo v Ku‖ ^ 2 := by
  have horth : ∀ i, inner ℝ (v i) (u - combo v Ku) = (0 : ℝ) :=
    (normalEq_iff_orthogonal v u Ku).mp hKu
  obtain ⟨hnormal, hJ⟩ := Jof_snoc_orthogonal v q (u - combo v Ku) hK horth hr
  have hmin : IsMinimizer (Fin.snoc v (u - combo v Ku)) q _ :=
    normalEq_isMinimizer _ q hnormal
  have hrisk := risk_eq_of_range_eq (v := Fin.snoc v u)
    (w := Fin.snoc v (u - combo v Ku)) q (range_combo_snoc_resid v u Ku) hbig hmin
  simp only [Jof, hrisk]
  exact hJ

end GSA.Part2
