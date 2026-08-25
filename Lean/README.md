# Lean формалізація (KuYuPe-Change_Point)

Файли у `Lean/` і підкаталогах формалізують частину 2 `paper/Part2_review.md`.

## Статус (13/08/2026)

| Метрика | Значення |
|---------|----------|
| Компіляція | ✅ Успішна (`lake build GSA`, 7764 jobs) |
| `sorry` | **0** |
| `error` | **0** |
| `axiom` | **0** (лише стандартні `propext`, `Classical.choice`, `Quot.sound`) |
| Покриття теорем | див. застереження нижче |

Усі теореми перевірено через `#print axioms` (`GSA/Audit.lean`): залежать лише
від `propext`, `Classical.choice`, `Quot.sound`.

## ⚠ Чиста збірка ≠ підтверджена теорема рукопису

`lake build` і аудит аксіом засвідчують, що **Lean-твердження** доведені. Вони
не засвідчують, що кожне Lean-твердження збігається з однойменним твердженням
препринта. Один такий розрив відомий і задокументований:

> **Теорема 2 препринта (arXiv:2605.23419) хибна в частині (a) і (c).** Вона
> ототожнює `J(s)` та його границю з дивергенцією Джеффріса
> `D_KL(f₁‖f₀) + D_KL(f₀‖f₁)`. Правильна границя за повного базису —
> `(1/c)·2Δ/(2−Δ)`, де `Δ` — трикутна (Vincze–Le Cam) дискримінація; при
> нормуванні препринта `c = 1/2` це `4Δ/(2−Δ)`. `InfoFunctional.lean` доводить
> **не це**, а факт Парсеваля про ортонормований базис. Жодна теорема цієї
> формалізації не пов'язує `J` з дивергенцією Джеффріса.

Джерела: `erratum/ERRATUM_theorem2c_2026-08-23.md`,
`erratum/THEOREM2C_AUDIT_2026-08-23.md` §5,
`GSA/Part2/BridgeGap.lean` (секція `Інвентар NOT FORMALISED`).

Шляхи `erratum/…`, які трапляються тут і в докстрінгах, вказують на приватне
дослідницьке репо треку. **Публічне джерело того самого твердження — сам
препринт: arXiv:2605.23419**, починаючи з **v3** (24.08.2026), де виправлені
формулювання, доведення й контрприклад Bernoulli(0.2, 0.8) стоять у тексті, а
`correction_note` — одразу після анотації.

**Не писати** «Теорема 2 Lean-verified». Писати можна: «варіаційна
характеризація `J = Yᵀ F⁻¹ Y` як максимуму відношення Релея формалізована».

## Ключові точки входу
- `GSA.lean` — головний модуль, що імпортує `GSA.Part2` разом із базовими лемами/перевірками.
- `GSA/Part2.lean` — точка входу розділу 2: послідовно імпортує всі тематичні підмодулі (`Setup`, `Kunchenko`, `Algorithm` тощо) після `Formalization` і `Checks`.

## Модулі під `GSA/Part2/`

### Базові модулі
- `Formalization.lean` — базові леми: `ku1_step3` для КУ1 та `approx_tendsto`/`coeff_eq_inner` про гільбертовий розклад.
- `Checks.lean` — прості `#check`/`#print`, щоб бачити типи ключових лем.
- `Setup.lean` — загальна термінологія: LLR, дивергенція Джеффріса, інтерфейс мір.
  `JeffreysENN` — лише означення: **жодна теорема з ним не працює**.

### Теореми (повністю формалізовані)
- `Kunchenko.lean` — Теорема 1 (КУ1): `psiOpt_pos_iff_lr_gt_one`, `theorem1_decision_rule_equivalence`, варіаційна секція
- `Convergence.lean` — Теорема 4 (L²-збіжність, оцінка хвоста, швидкість ~260 рядків)
- `InfoFunctional.lean` — часткові суми Парсеваля по **ортонормованому** базису:
  `J b z s ≤ ‖z‖²`, монотонність, `J b z s → ‖z‖²`. Імена `theorem2_*` — данина
  аудиторському сліду; Теорему 2 препринта це **не** підтверджує (див. вище).
- `BridgeGap.lean` — місток «скінченновимірна коваріаційна система → відношення
  Релея»: `J = Yᵀ F⁻¹ Y = max_K (KᵀY)²/(KᵀFK)`, максимум на `K = F⁻¹Y`
  (`isGreatest_rayleigh`, `isGreatest_rayleigh_eq_quadForm_inv`). Той самий файл
  містить явний інвентар того, що **не** формалізовано, — насамперед ототожнення
  границі з дивергенцією Джеффріса.
- `RobustPE.lean` — Теорема 3, 5: `PE_threshold`, `criterion_Y`, `KU1_le_criterion_Y`, `theorem3_PE_NP_threshold_equivalence`
- `FAR_ADD.lean` — Теорема 6: `exceed_le_by_chebyshev`, `exceed_le_eps`
- `UnimodalBounds.lean` — нерівності Кантеллі, В.-П.: `cantelli_bound`, `VP_bound`, `PE_threshold_*`
- `LinearSystem.lean` — система FK=Y: `has_solution_FK_eq_Y`
- `Kernel.lean` — **theorem kernel системи `F K = B`** (доданий 13/08/2026): алгебраїчне
  ядро §2.4 для **неортогонального** словника. `InfoFunctional.lean` доводить властивості
  `J(s)` лише в ортонормованому `HilbertBasis`; відвантажені словники (`Φ_poly`, `Φ_log`,
  `Φ_frac`) не ортогональні, і саме звідти беруться всі питання §2.5 про обумовленість.
  `J` означено варіаційно (`Jof = ‖q‖² − risk`), формула `BᵀF⁻¹B` — теорема.
  Усе виводиться з однієї геометричної тотожності `normalEq_iff_orthogonal`
  (нормальні рівняння ⟺ залишок ортогональний до span):
  | § | Твердження | Теорема |
  |---|---|---|
  | 1 | існування та єдиність `F K = B` | `existsUnique_normalEq_of_posDef` |
  | 2 | `K*` мінімізує квадратичний ризик | `normalEq_isMinimizer` |
  | 3 | residual identity | `risk_eq_risk_add_normSq` |
  | 4 | інформаційна межа `J ≤ ‖q‖²` | `Jof_le_normSq` |
  | 5 | монотонність за вкладеними словниками | `Jof_init_le_Jof_of_isMinimizer` |
  | 6 | Schur-приріст від однієї фічі | `Jof_snoc_general` |
  | 7 | `J` залежить лише від span | `risk_eq_of_range_eq` |
  | 8 | межа ортогоналізації | `isMinimizer_transfer` |
  | 9 | вироджена `F`: `J` коректно визначена | `risk_eq_of_normalEq_of_normalEq` |
- `GramSchmidt.lean` — емпірична ортонормалізація: `empGramSchmidt_orthonormal`
- `PolynomialBasis.lean` — умови на моменти: `HasMoment`, `PolynomialBasisApplicable`, `BasisType`, `select_basis_by_excess`

### Допоміжні модулі
- `Problem.lean`, `RelativeDisorder.lean`, `Algorithm.lean` — блоки, близькі до тексту статті
- `BasisApprox.lean` — дефініції `coeff`/`approx` для ортонормованих систем
- `GaussianLimit.lean` — гаусівський ліміт
- `Universality.lean`, `Autocorrelation.lean` — допоміжні оцінки/функціонали

## Збірка
```bash
# З кореня репозиторія
lake build GSA
```

- Lean 4.26.0 / Mathlib 4.26.0 (налаштовано в `lean-toolchain` і `lakefile.lean`)
- Формалізація розділу 2 збирається без `sorry`
- Власних аксіом немає: `hasMoment_of_higher` доведено в `PolynomialBasis.lean`
- Аудит аксіом: `lake env lean Lean/GSA/Audit.lean`

## Див. також
- `reports/Part2_analysis_report.md` — детальний звіт аналізу теорії та формалізації
