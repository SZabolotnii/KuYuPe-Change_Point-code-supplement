# Lean формалізація (KuYuPe-Change_Point)

Файли у `Lean/` і підкаталогах формалізують частину 2 `paper/Part2_review.md`.

## Статус (10/01/2026)

| Метрика | Значення |
|---------|----------|
| Компіляція | ✅ Успішна |
| `sorry` | **0** |
| `error` | **0** |
| `axiom` | 1 |
| Покриття теорем | ~95% |

## Ключові точки входу
- `GSA.lean` — головний модуль, що імпортує `GSA.Part2` разом із базовими лемами/перевірками.
- `GSA/Part2.lean` — точка входу розділу 2: послідовно імпортує всі тематичні підмодулі (`Setup`, `Kunchenko`, `Algorithm` тощо) після `Formalization` і `Checks`.

## Модулі під `GSA/Part2/`

### Базові модулі
- `Formalization.lean` — базові леми: `ku1_step3` для КУ1 та `approx_tendsto`/`coeff_eq_inner` про гільбертовий розклад.
- `Checks.lean` — прості `#check`/`#print`, щоб бачити типи ключових лем.
- `Setup.lean` — загальна термінологія: LLR, Jeffreys дивергенція та інтерфейс мір.

### Теореми (повністю формалізовані)
- `Kunchenko.lean` — Теорема 1 (КУ1): `psiOpt_pos_iff_lr_gt_one`, `theorem1_decision_rule_equivalence`, варіаційна секція
- `Convergence.lean` — Теорема 4 (L²-збіжність, оцінка хвоста, швидкість ~260 рядків)
- `InfoFunctional.lean` — Теорема 2 (монотонність J(s), збіжність)
- `RobustPE.lean` — Теорема 3, 5: `PE_threshold`, `criterion_Y`, `KU1_le_criterion_Y`, `theorem3_PE_NP_threshold_equivalence`
- `FAR_ADD.lean` — Теорема 6: `exceed_le_by_chebyshev`, `exceed_le_eps`
- `UnimodalBounds.lean` — нерівності Кантеллі, В.-П.: `cantelli_bound`, `VP_bound`, `PE_threshold_*`
- `LinearSystem.lean` — система FK=Y: `has_solution_FK_eq_Y`
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
- 1 `axiom` (`hasMoment_of_higher`) — стандартний результат теорії ймовірностей про існування нижчих моментів

## Див. також
- `reports/Part2_analysis_report.md` — детальний звіт аналізу теорії та формалізації
