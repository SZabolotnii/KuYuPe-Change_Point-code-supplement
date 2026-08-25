import Mathlib.Probability.Moments.Variance
import Mathlib.MeasureTheory.Measure.Typeclasses.Probability
import Mathlib.Tactic.Linarith
import GSA.Part2.RobustPE
import GSA.Part2.FAR_ADD

open MeasureTheory ProbabilityTheory
open scoped ENNReal BigOperators

namespace GSA.Part2.MAry

/-!
# Межа ймовірності помилки для багатоальтернативного правила (критерій ПЕ)

Продовження `GSA.Part2.MAry` (`MultiAlternative.lean`) з боку ймовірностей.
Там показано, коли попарні правила КУ1 узгоджені (Т2) і коли ні (Т1).
Тут — скільки коштує помилка, якщо рішення все одно ухвалюється попарно.

Ланцюг такий:

1. `copeland_of_wins_all` — комбінаторна частина. Якщо гіпотеза `m` виграла
   всі `M-1` дуелей, вона переможець за Коуплендом **однозначно**. Отже цикли
   з Т1 не можуть зіпсувати рішення там, де правильна гіпотеза бездоганна.
2. `error_subset_duels` — звідси помилка вкладена в об'єднання подій «дуель `n`
   програна», а їх лише `M-1`.
3. `PE_duel_bound` — за порогом ПЕ `h = E[Λ] + √(Var[Λ]/ε)` дуель програється
   з імовірністю не більшою за `ε`. Це бінарний результат `exceed_le_eps`
   з `FAR_ADD.lean`; тут лише показано, що подія `{Λ ≥ h_ПЕ}` — це та сама
   `exceed`, тож переносити доведення не треба.
4. `PE_error_bound` — об'єднавча межа: `P(помилка | H_m) ≤ Σ_{n≠m} ε_n`,
   а за рівномірного розподілу бюджету `ε_n = ε/(M-1)` — просто `ε`
   (`PE_error_bound_uniform`).

Ціна багатоальтернативності — саме множник `M-1` у кроці 4: щоб утримати
сумарну похибку на рівні `ε`, кожна дуель мусить бути жорсткішою в `M-1` разів,
що піднімає поріг у `√(M-1)` разів.

Жодних припущень про незалежність дуелей не робиться — тому межа груба, але
чинна за будь-якої залежності між `Λ_{mn}`.

Крок 3 успадковує від критерію ПЕ саме нерівність Чебишева. Однобічна межа
Кантеллі дала б тугішу оцінку (`cantelli_lt_chebyshev` в `UnimodalBounds.lean`),
але там вона доведена лише як арифметична нерівність між межами, без
міро-теоретичної версії, тож підставити її сюди поки нема чого.
-/

variable {Ω : Type*} [MeasurableSpace Ω] {μ : Measure Ω}
variable {ι : Type*} [Fintype ι] [DecidableEq ι]

/-! ### 1. Комбінаторна частина: бездоганна гіпотеза перемагає за Коуплендом -/

/-- Кількість попарних перемог гіпотези `a` у турнірі `beats`. -/
def score (beats : ι → ι → Bool) (a : ι) : ℕ :=
  (Finset.univ.filter fun b => beats a b = true).card

/-- **L5.1.** Якщо `m` виграла всі дуелі, то її рахунок Коупленда строго вищий
за рахунок будь-якої іншої гіпотези.

Доведення не спирається на переходовість: суперник `n` програв щонайменше одну
дуель — саме дуель проти `m`, — тож його рахунок не досягає максимуму `M-1`.
Тому контрприклад `exists_intransitive` з `MultiAlternative.lean` не заважає:
цикли виникають лише серед гіпотез, які вже комусь програли. -/
theorem copeland_of_wins_all (beats : ι → ι → Bool)
    (hirr : ∀ a, beats a a = false)
    (hasym : ∀ a b, beats a b = true → beats b a = false)
    {m : ι} (hm : ∀ n, n ≠ m → beats m n = true) {n : ι} (hn : n ≠ m) :
    score beats n < score beats m := by
  classical
  have hcard : 2 ≤ Fintype.card ι := Finset.one_lt_card.mpr
    ⟨n, Finset.mem_univ _, m, Finset.mem_univ _, hn⟩
  -- `m` перемагає рівно всіх, крім себе
  have hM : score beats m = Fintype.card ι - 1 := by
    have : (Finset.univ.filter fun b => beats m b = true) = Finset.univ.erase m := by
      ext b
      simp only [Finset.mem_filter, Finset.mem_univ, true_and, Finset.mem_erase, and_true]
      constructor
      · intro hb hbm; rw [hbm, hirr] at hb; exact Bool.noConfusion hb
      · intro hb; exact hm b hb
    rw [score, this, Finset.card_erase_of_mem (Finset.mem_univ _), Finset.card_univ]
  -- `n` не перемагає ні себе, ні `m`
  have hN : score beats n ≤ Fintype.card ι - 2 := by
    have hsub : (Finset.univ.filter fun b => beats n b = true)
        ⊆ (Finset.univ.erase m).erase n := by
      intro b hb
      simp only [Finset.mem_filter, Finset.mem_univ, true_and] at hb
      refine Finset.mem_erase.mpr ⟨?_, Finset.mem_erase.mpr ⟨?_, Finset.mem_univ _⟩⟩
      · rintro rfl; rw [hirr] at hb; exact Bool.noConfusion hb
      · rintro rfl; rw [hasym _ _ (hm n hn)] at hb; exact Bool.noConfusion hb
    calc score beats n ≤ ((Finset.univ.erase m).erase n).card := Finset.card_le_card hsub
      _ = Fintype.card ι - 2 := by
          rw [Finset.card_erase_of_mem (Finset.mem_erase.mpr ⟨hn, Finset.mem_univ _⟩),
            Finset.card_erase_of_mem (Finset.mem_univ _), Finset.card_univ]
          omega
  omega

/-! ### 2. Подія помилки вкладена в об'єднання програних дуелей -/

/-- **L5.2.** Якщо правило `choose` повертає `m` щоразу, коли `m` виграє всі
дуелі, то помилка можлива лише через програш якоїсь однієї дуелі.

Гіпотеза `hchoose` — це саме те, що дає `copeland_of_wins_all` для агрегатора
Коупленда; але твердження сформульовано для будь-якого агрегатора з цією
властивістю (Коупленд, зважений Коупленд, MFAS — усі її мають). -/
theorem error_subset_duels (win : ι → Ω → Prop) (choose : Ω → ι) (m : ι)
    (hchoose : ∀ ω, (∀ n ∈ Finset.univ.erase m, win n ω) → choose ω = m) :
    {ω | choose ω ≠ m} ⊆ ⋃ n ∈ Finset.univ.erase m, {ω | ¬ win n ω} := by
  intro ω hω
  by_contra hcon
  simp only [Set.mem_iUnion, Set.mem_setOf_eq, not_exists, not_not] at hcon
  exact hω (hchoose ω fun n hn => hcon n hn)

/-! ### 3. Поріг ПЕ дає `ε` на одну дуель -/

/-- **L5.3.** За порогом ПЕ `h = E[Λ] + √(Var[Λ]/ε)` дуель програється з
імовірністю не більшою за `ε`.

Це бінарний результат `GSA.Part2.exceed_le_eps` (`FAR_ADD.lean`): подія
`{Λ ≥ h_ПЕ}` збігається з `exceed μ Λ √(Var[Λ]/ε)` за означенням, бо
`PE_threshold mean v ε = mean + √(v/ε)`. Переформульовано лише заради того,
щоб далі говорити мовою порогів, а не відступів. -/
theorem PE_duel_bound [IsFiniteMeasure μ] {Λ : Ω → ℝ} (hΛ : MemLp Λ 2 μ)
    {ε : ℝ} (hε : 0 < ε) (hv : 0 < variance Λ μ) :
    μ {ω | GSA.Part2.PE_threshold (μ[Λ]) (variance Λ μ) ε ≤ Λ ω} ≤ ENNReal.ofReal ε :=
  GSA.Part2.exceed_le_eps μ Λ hΛ ε hε hv

/-! ### 4. Об'єднавча межа для `M` гіпотез -/

/-- **L5 (головна).** Імовірність помилки багатоальтернативного правила під
`H_m` не перевищує суми бюджетів окремих дуелей.

Незалежність дуелей не потрібна — об'єднавча межа чинна за будь-якої
залежності між статистиками `Λ_{mn}`. -/
theorem PE_error_bound [IsProbabilityMeasure μ] {m : ι}
    (Lam : ι → Ω → ℝ) (hL : ∀ n, MemLp (Lam n) 2 μ)
    (hv : ∀ n, 0 < variance (Lam n) μ)
    (εs : ι → ℝ) (hε : ∀ n, 0 < εs n)
    (choose : Ω → ι)
    (hchoose : ∀ ω, (∀ n ∈ Finset.univ.erase m,
        Lam n ω < GSA.Part2.PE_threshold (μ[Lam n]) (variance (Lam n) μ) (εs n)) →
      choose ω = m) :
    μ {ω | choose ω ≠ m} ≤ ENNReal.ofReal (∑ n ∈ Finset.univ.erase m, εs n) := by
  classical
  set thr : ι → ℝ := fun n =>
    GSA.Part2.PE_threshold (μ[Lam n]) (variance (Lam n) μ) (εs n) with hthr
  have hsub := error_subset_duels (fun n ω => Lam n ω < thr n) choose m hchoose
  have hstep : μ {ω | choose ω ≠ m}
      ≤ ∑ n ∈ Finset.univ.erase m, μ {ω | ¬ Lam n ω < thr n} := by
    refine le_trans (measure_mono hsub) ?_
    exact measure_biUnion_finset_le _ _
  refine le_trans hstep ?_
  have hterm : ∀ n ∈ Finset.univ.erase m,
      μ {ω | ¬ Lam n ω < thr n} ≤ ENNReal.ofReal (εs n) := by
    intro n _
    have : {ω | ¬ Lam n ω < thr n} = {ω | thr n ≤ Lam n ω} := by
      ext ω; simp [not_lt]
    rw [this]
    exact PE_duel_bound (hL n) (hε n) (hv n)
  refine le_trans (Finset.sum_le_sum hterm) ?_
  rw [← ENNReal.ofReal_sum_of_nonneg fun n _ => (hε n).le]

/-- **Наслідок.** Рівномірний розподіл бюджету `ε_n = ε/(M-1)` утримує сумарну
похибку на рівні `ε`. Це і є ціна багатоальтернативності: кожна дуель мусить
бути в `M-1` разів жорсткішою, тобто поріг зростає в `√(M-1)` разів. -/
theorem PE_error_bound_uniform [IsProbabilityMeasure μ] {m : ι}
    (Lam : ι → Ω → ℝ) (hL : ∀ n, MemLp (Lam n) 2 μ)
    (hv : ∀ n, 0 < variance (Lam n) μ)
    {ε : ℝ} (hε : 0 < ε) (hcard : 2 ≤ Fintype.card ι)
    (choose : Ω → ι)
    (hchoose : ∀ ω, (∀ n ∈ Finset.univ.erase m,
        Lam n ω < GSA.Part2.PE_threshold (μ[Lam n]) (variance (Lam n) μ)
          (ε / ((Fintype.card ι : ℝ) - 1))) → choose ω = m) :
    μ {ω | choose ω ≠ m} ≤ ENNReal.ofReal ε := by
  classical
  have hN : (0 : ℝ) < (Fintype.card ι : ℝ) - 1 := by
    have : (2 : ℝ) ≤ (Fintype.card ι : ℝ) := by exact_mod_cast hcard
    linarith
  have hcards : ((Finset.univ.erase m).card : ℝ) = (Fintype.card ι : ℝ) - 1 := by
    rw [Finset.card_erase_of_mem (Finset.mem_univ _), Finset.card_univ]
    have : 1 ≤ Fintype.card ι := by omega
    push_cast [Nat.cast_sub this]
    ring
  refine le_trans (PE_error_bound Lam hL hv (fun _ => ε / ((Fintype.card ι : ℝ) - 1))
    (fun _ => div_pos hε hN) choose hchoose) ?_
  rw [Finset.sum_const, nsmul_eq_mul, hcards, mul_comm, div_mul_cancel₀ _ hN.ne']

end GSA.Part2.MAry
