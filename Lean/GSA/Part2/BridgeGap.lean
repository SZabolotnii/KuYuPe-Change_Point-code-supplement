import GSA.Part2.Kernel

open scoped BigOperators
open Finset

namespace GSA.Part2

/-!
# Місток `J = Yᵀ F⁻¹ Y = sup Rayleigh` — і те, чого в ньому немає

Модуль робить дві різні речі, і їх не можна плутати.

## Частина 1 (доведено): варіаційна характеризація `J`

Для скінченновимірного простору ознак із матрицею Грама `F` і вектором
різниці середніх `Y` величина рукопису `J = Yᵀ F⁻¹ Y` дорівнює максимуму
відношення Релея

  `J = sup_{K} (Kᵀ Y)² / (Kᵀ F K)`,

і максимум досягається на `K* = F⁻¹ Y`. Це `isGreatest_rayleigh` та
`isGreatest_rayleigh_eq_quadForm_inv` нижче. Твердження скінченновимірне,
доведене без `sorry` і без власних аксіом.

Це — чесний зміст `J`. Воно не потребує ортонормованості словника (на відміну
від `InfoFunctional.lean`) і не потребує невиродженості `F`
(`isGreatest_rayleigh` бере довільний розв'язок нормальних рівнянь).

## Частина 2 (НЕ доведено): статистична ідентифікація границі

Нижче, у секції «Інвентар `NOT FORMALISED`», перелічено рівно ті кроки, які
рукопис робить понад Частину 1 і які **тут не формалізовані**. Головний із них:
ототожнення границі `J` за повного базису з якоюсь конкретною дивергенцією —
і, зокрема, **з дивергенцією Джеффріса. Це ототожнення хибне.**

Джерела: `erratum/ERRATUM_theorem2c_2026-08-23.md`,
`erratum/THEOREM2C_AUDIT_2026-08-23.md` §5, `erratum/verify_theorem2c.py`.
-/

variable {H : Type*} [NormedAddCommGroup H] [InnerProductSpace ℝ H]
variable {s : ℕ}

/-! ## 1. Відношення Релея -/

/-- Відношення Релея словника: `(Kᵀ Y)² / (Kᵀ F K)`, де `Y = bvec v q` —
вектор кореляцій фіч із напрямом зміни, а `F = gram v` — матриця Грама.

Ділення в Lean тотальне: при `Kᵀ F K = 0` значення дорівнює `0`. Це не
патч, а зручність — у виродженому напрямі `Kᵀ Y` теж нуль
(див. `rayleigh_le_Jof`), тож нуль є коректним значенням, і формулювання
`IsGreatest` не потребує окремої умови невиродженості. -/
noncomputable def rayleigh (v : Fin s → H) (q : H) (K : Fin s → ℝ) : ℝ :=
  (∑ i, K i * bvec v q i) ^ 2 / (∑ i, ∑ j, K i * K j * gram v i j)

/-- Чисельник відношення Релея — це скалярний добуток наближення з `q`, а для
розв'язку нормальних рівнянь — з оптимальним наближенням. -/
lemma dot_bvec_eq_inner_combo (v : Fin s → H) (q : H) {Ks : Fin s → ℝ}
    (hKs : NormalEq v q Ks) (K : Fin s → ℝ) :
    (∑ i, K i * bvec v q i) = inner ℝ (combo v K) (combo v Ks) := by
  have h0 : inner ℝ (combo v K) (q - combo v Ks) = (0 : ℝ) :=
    inner_combo_resid_eq_zero v q hKs K
  rw [inner_sub_right, sub_eq_zero] at h0
  rw [← h0, inner_combo_left]
  rfl

/-- **Межа Релея.** Для будь-якого `K` відношення Релея не перевищує захопленої
інформації `J`. Це нерівність Коші–Буняковського в span словника. -/
theorem rayleigh_le_Jof (v : Fin s → H) (q : H) {Ks : Fin s → ℝ}
    (hKs : NormalEq v q Ks) (K : Fin s → ℝ) :
    rayleigh v q K ≤ Jof v q Ks := by
  have hnum := dot_bvec_eq_inner_combo v q hKs K
  have hden : (∑ i, ∑ j, K i * K j * gram v i j) = ‖combo v K‖ ^ 2 :=
    (normSq_combo_eq_quadForm v K).symm
  have hJ : Jof v q Ks = ‖combo v Ks‖ ^ 2 := Jof_eq_normSq_combo v q hKs
  simp only [rayleigh]
  rw [hnum, hden, hJ]
  rcases eq_or_lt_of_le (sq_nonneg ‖combo v K‖) with h | h
  · rw [← h, div_zero]
    positivity
  · rw [div_le_iff₀ h]
    have h1 : |inner ℝ (combo v K) (combo v Ks)| ≤ ‖combo v K‖ * ‖combo v Ks‖ :=
      abs_real_inner_le_norm _ _
    have h2 : (0 : ℝ) ≤ |inner ℝ (combo v K) (combo v Ks)| := abs_nonneg _
    have h3 : |inner ℝ (combo v K) (combo v Ks)| ^ 2
        = (inner ℝ (combo v K) (combo v Ks) : ℝ) ^ 2 := sq_abs _
    nlinarith [h1, h2, h3, norm_nonneg (combo v K), norm_nonneg (combo v Ks)]

/-- Нульовий вектор коефіцієнтів дає нульове відношення Релея. -/
lemma rayleigh_zero (v : Fin s → H) (q : H) : rayleigh v q (0 : Fin s → ℝ) = 0 := by
  simp [rayleigh]

/-- **Досяжність.** На розв'язку нормальних рівнянь `K* = F⁻¹ Y` відношення
Релея дорівнює `J` (за `J ≠ 0`; при `J = 0` максимум досягається на `K = 0`,
див. `isGreatest_rayleigh`). -/
theorem rayleigh_eq_Jof_of_normalEq (v : Fin s → H) (q : H) {Ks : Fin s → ℝ}
    (hKs : NormalEq v q Ks) (hJ : Jof v q Ks ≠ 0) :
    rayleigh v q Ks = Jof v q Ks := by
  have hnum : (∑ i, Ks i * bvec v q i) = Jof v q Ks :=
    (Jof_eq_dot_of_normalEq v q hKs).symm
  have hden : (∑ i, ∑ j, Ks i * Ks j * gram v i j) = Jof v q Ks := by
    rw [← normSq_combo_eq_quadForm]
    exact (Jof_eq_normSq_combo v q hKs).symm
  simp only [rayleigh]
  rw [hnum, hden, sq, mul_div_assoc, div_self hJ, mul_one]

/-- **Місток (основна форма).** Захоплена інформація `J` — це рівно максимум
відношення Релея по всіх векторах коефіцієнтів:

  `J = max_K (Kᵀ Y)² / (Kᵀ F K)`.

Невиродженість `F` не потрібна: достатньо будь-якого розв'язку нормальних
рівнянь `F K = Y`, а `J` від його вибору не залежить
(`risk_eq_of_normalEq_of_normalEq`). -/
theorem isGreatest_rayleigh (v : Fin s → H) (q : H) {Ks : Fin s → ℝ}
    (hKs : NormalEq v q Ks) :
    IsGreatest (Set.range (rayleigh v q)) (Jof v q Ks) := by
  constructor
  · rcases eq_or_ne (Jof v q Ks) 0 with h | h
    · exact ⟨0, by rw [rayleigh_zero, h]⟩
    · exact ⟨Ks, rayleigh_eq_Jof_of_normalEq v q hKs h⟩
  · rintro r ⟨K, rfl⟩
    exact rayleigh_le_Jof v q hKs K

/-- Явний розв'язок `K* = F⁻¹ Y` задовольняє нормальні рівняння за додатно
визначеної `F`. -/
lemma normalEq_inv_mulVec (v : Fin s → H) (q : H) (hF : (gram v).PosDef) :
    NormalEq v q ((gram v)⁻¹.mulVec (bvec v q)) := by
  classical
  have hunit : IsUnit (gram v).det :=
    isUnit_iff_ne_zero.mpr (ne_of_gt (Matrix.PosDef.det_pos hF))
  have _ := Matrix.invertibleOfIsUnitDet (A := gram v) hunit
  have hmv : (gram v).mulVec ((gram v)⁻¹.mulVec (bvec v q)) = bvec v q := by
    rw [Matrix.mulVec_mulVec, Matrix.mul_inv_of_invertible, Matrix.one_mulVec]
  intro i
  simpa [Matrix.mulVec, dotProduct] using congrFun hmv i

/-- **Місток (матрична форма).** Для `F ≻ 0`

  `Yᵀ F⁻¹ Y = max_{K} (Kᵀ Y)² / (Kᵀ F K)`,

і максимум досягається на `K = F⁻¹ Y`. Це рівно те твердження, яке аудит
(`THEOREM2C_AUDIT_2026-08-23.md`, §9 P2.1) називає відсутнім містком «від
скінченновимірної коваріаційної системи до відношення Релея». -/
theorem isGreatest_rayleigh_eq_quadForm_inv (v : Fin s → H) (q : H)
    (hF : (gram v).PosDef) :
    IsGreatest (Set.range (rayleigh v q))
      (∑ i, ((gram v)⁻¹.mulVec (bvec v q)) i * bvec v q i) := by
  have hK := normalEq_inv_mulVec v q hF
  have h := isGreatest_rayleigh v q hK
  rwa [Jof_eq_dot_of_normalEq v q hK] at h

/-! ## 2. Арифметика правильної границі

Наступна лема — **арифметична половина** виправленої границі. Вона показує:
якщо для оптимального напряму `h*` виконуються дві моментні тотожності

  `E₁h* − E₀h* = Δ`   і   `Var₀(h*) + Var₁(h*) = Δ − Δ²/2`,

то значення відношення Релея при нормуванні `F = c(C₀ + C₁)` дорівнює
`(1/c)·2Δ/(2−Δ)`.

Самі моментні тотожності тут **не доведені** — див. інвентар нижче. -/

/-- Значення відношення Релея на оптимальному напрямі: `Δ²/(c(Δ − Δ²/2))`
згортається до `(1/c)·2Δ/(2−Δ)`.

За `c = 1/2` (нормування arXiv v2) це `4Δ/(2−Δ)`; за `c = 1` (arXiv v1 і код) —
`2Δ/(2−Δ)`. Для `P₀ = N(0,1)`, `P₁ = N(1,1)` маємо `Δ ≈ 0.4081085`, звідки
`4Δ/(2−Δ) ≈ 1.0254682`, тоді як дивергенція Джеффріса дорівнює рівно `1`.
Розбіжність у третьому знаку — це і є спростування Теореми 2(c). -/
theorem rayleighValue_of_triangular {c Δ : ℝ}
    (hc : c ≠ 0) (hΔ : 0 < Δ) (hΔ2 : Δ < 2) :
    Δ ^ 2 / (c * (Δ - Δ ^ 2 / 2)) = (1 / c) * (2 * Δ / (2 - Δ)) := by
  have h1 : Δ ≠ 0 := ne_of_gt hΔ
  have h2 : (2 : ℝ) - Δ ≠ 0 := by linarith
  have h3 : Δ - Δ ^ 2 / 2 ≠ 0 := by
    have : Δ - Δ ^ 2 / 2 = Δ * (2 - Δ) / 2 := by ring
    rw [this]
    positivity
  field_simp

/-! ## 3. Інвентар `NOT FORMALISED`

Нижче — вичерпний перелік кроків, які рукопис робить понад доведене вище і
які **в цій формалізації не доводяться**. Перелік навмисно записаний прозою:
його призначення — щоб жоден майбутній читач не міг прийняти Частину 1 за
підтвердження Теореми 2 рукопису.

-- NOT FORMALISED (1) — **ототожнення границі з дивергенцією Джеффріса.**
   Ніде не доведено, що границя `J` за повного базису дорівнює
   `D_KL(f₁‖f₀) + D_KL(f₀‖f₁)`. Понад те: це твердження **хибне**.
   Контрприклад (точний, скінченний носій): для Bernoulli `p₀ = 0.2`,
   `p₁ = 0.8` при `c = 1/2` границя дорівнює `2.25`, а `D_J ≈ 1.6635532`.
   Той самий контрприклад одночасно спростовує заявлену верхню межу
   Теореми 2(a). Джерело: `erratum/ERRATUM_theorem2c_2026-08-23.md` §4.2.

-- NOT FORMALISED (2) — **ототожнення границі з будь-якою іншою дивергенцією.**
   Правильна границя — `(1/c)·2Δ/(2−Δ)`, де
   `Δ = ∫ (f₁−f₀)²/(f₁+f₀)` — трикутна (Vincze–Le Cam) дискримінація.
   Формалізовано лише арифметику цього виразу (`rayleighValue_of_triangular`),
   але не його виведення.

-- NOT FORMALISED (3) — **оптимальний напрям.**
   Не доведено, що супремум відношення Релея по **всіх** вимірних `h`
   (а не по скінченновимірному span) досягається на
   `h*(x) = (f₁−f₀)/(f₁+f₀) = tanh(ℓ(x)/2)`, де `ℓ = log(f₁/f₀)`.
   Причина: потрібен вимірно-теоретичний шар (RN-похідні, `L²(f₀+f₁)`,
   факторизація за константами), непропорційний до виграшу.

-- NOT FORMALISED (4) — **дві моментні тотожності.**
   `E₁h* − E₀h* = Δ` та `Var₀(h*) + Var₁(h*) = Δ − Δ²/2` взято як задані
   в `rayleighValue_of_triangular`, а не доведено.

-- NOT FORMALISED (5) — **перехід від скінченного `s` до повного базису.**
   `isGreatest_rayleigh` — твердження при фіксованому скінченному `s`.
   Щільність об'єднання вкладених `V_s` у відповідному коваріаційному
   гільбертовому просторі, а отже й збіжність `J_s` до глобального супремуму,
   тут не формалізована. `theorem2_c_tendsto` в `InfoFunctional.lean` дає
   збіжність лише для **ортонормованого** базису і лише до `‖z‖²`.

-- NOT FORMALISED (6) — **ототожнення `q` з LLR.**
   У всій формалізації `q : H` — довільний вектор. Ніде не сказано і не
   використано, що це логарифм відношення правдоподібності. Саме через це
   Частина 1 не має і не може мати статистичного змісту сама по собі.

-- NOT FORMALISED (7) — **узгодження нормування.**
   Стаття v2 бере `F = ½(C₀+C₁)`, а відвантажений код — `F = C₀+C₁`.
   Це різниця вдвічі у всіх значеннях `J`. Формалізація не фіксує жодного з
   двох нормувань: `gram` — абстрактна матриця Грама.

### Що з цього випливає для формулювань

Дозволено писати: «варіаційна характеризація `J = Yᵀ F⁻¹ Y` як максимуму
відношення Релея формалізована в Lean і перевірена на стандартних аксіомах
Mathlib».

Заборонено писати: «Теорема 2 Lean-verified», «Lean підтверджує збіжність до
дивергенції Джеффріса», «захоплена частка дивергенції Джеффріса». Останнє
формулювання слід замінити на частку від правильної стелі
`(1/c)·2Δ/(2−Δ)` або зняти нормовану інтерпретацію взагалі.
-/

end GSA.Part2
