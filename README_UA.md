# GSA Детектор Точок Зміни

Програмний код до статті:

> **"Узагальнена стохастична апроксимація логарифму відношення правдоподібності для робастного послідовного виявлення точок зміни"**
>
> Сергій Заболотній

## Основні переваги

- **Апроксимація LLR на основі моментів** -- не потребує знання щільності розподілу, лише моменти до порядку 2s
- **Три типи базисів:** поліноміальний, дробовий, логарифмічний (та Hermite)
- **PE-критерій Кунченка** для гарантованого контролю FAR через нерівності Чебишева / Височанського-Петуніна / Кантеллі
- **Працює на важкохвостих даних**, де класичні методи відмовляють (kurtosis > 20)
- **O(s) на відлік** -- придатний для edge/embedded пристроїв
- **Формальні докази на Lean 4** (Mathlib) для основних теорем
- **Модульні правила зупинки:** CUSUM, GRSh (баєсівське), Shiryaev-Roberts

## Встановлення

```bash
# З вихідного коду (рекомендовано для експериментів)
git clone https://github.com/SZabolotnii/KuYuPe-Change_Point-code-supplement.git
cd KuYuPe-Change_Point-code-supplement
pip install -e "."

# З залежностями для експериментів (matplotlib, pandas, ruptures)
pip install -e ".[experiments]"

# Все разом (експерименти + тестування)
pip install -e ".[all]"
```

**Вимоги:** Python >= 3.9, NumPy >= 1.21, SciPy >= 1.7.

## Швидкий старт

```python
import numpy as np
from gsa_cpd import GSADetector, BasisType, ThresholdType

# 1. Калібрування на даних нормального режиму
calibration = np.random.normal(0, 1, size=2000)

detector = GSADetector(
    basis=BasisType.FRAC,       # дробовий базис -- найкращий для важких хвостів
    degree=2,                    # порядок апроксимації s=2
    epsilon=0.01,                # цільовий рівень FAR
    threshold_type=ThresholdType.CHEBYSHEV,
)
detector.fit(calibration, delta=0.3)

# 2. Моніторинг потоку даних
stream = np.concatenate([
    np.random.normal(0, 1, 500),   # нормальний режим
    np.random.normal(1, 1, 200),   # після зміни (зсув середнього)
])

for t, x in enumerate(stream):
    if detector.predict(x):
        print(f"Зміну виявлено на t={detector.alarm_time}")
        break

# 3. Діагностика
diag = detector.diagnostics
print(f"Поріг: {diag.threshold:.3f}")
print(f"J(s): {diag.J_s:.4f}")
print(f"cond(F): {diag.condition_number:.1f}")
print(f"Метод розв'язку: {diag.solver_method}")
```

## Відтворення результатів статті

> Експерименти запускайте як модулі з кореня репозиторію (після
> `pip install -e .`), щоб коректно розв'язувалися і встановлений пакет
> `gsa_cpd`, і локальний пакет `experiments`.

### Monte Carlo моделювання (Розділ 4)

```bash
python -m experiments.monte_carlo.run_all --quick    # швидко (~5 хв)
python -m experiments.monte_carlo.run_all            # повний (~2 год)

# Або окремий експеримент, напр. перевірка гаусівської границі (Теорема 1):
python -m experiments.monte_carlo.exp01_gaussian_limit --quick
```

Результати записуються в каталоги з мітками часу під `results/`.

### Бенчмарки на реальних даних (Розділ 5)

Результати Розділу 5 (NASA IMS Bearing, NSL-KDD, SKAB, TCPD, макроряди FRED,
PhysioNet 2019 тощо) спираються на великі зовнішні датасети, частина яких
потребує реєстрації або окремого завантаження. **Ці скрипти не входять до
цього додатку.** Репозиторій зосереджений на повністю самодостатніх,
відтворюваних частинах статті: Monte Carlo (Розділ 4), формальні докази
(Розділ 2) та пакет `gsa_cpd` з набором тестів. Джерела датасетів і
наведені результати на реальних даних див. у Розділі 5 статті та
[docs/DATASETS.md](docs/DATASETS.md).

### Формальні докази (Розділ 2)

```bash
cd Lean && lake build GSA
```

Основні файли: `InfoFunctional.lean` (Теорема 2), `Convergence.lean` (Теорема 4), `FAR_ADD.lean` (Теорема 6).

### Запуск тестів

```bash
pytest tests/ -v
```

## Основні результати

| Сценарій | Перевага GSA | Відтворюється тут |
|---|---|---|
| Гаусівська границя | S=1 poly = класичний CUSUM (точний збіг, перевірено) | Так (MC + тест) |
| Негаусові дані (gamma_3 >= 8) | Зменшення ADD на 30--36% порівняно з класичним CUSUM | Так (MC) |
| Важкі хвости (kurtosis > 20) | Єдиний працюючий метод (датасет NASA IMS Bearing) | Розділ 5 статті |
| Кібербезпека (NSL-KDD) | FAR = 0%, DetRate = 100% | Розділ 5 статті |

Рядки з позначкою «Розділ 5 статті» — бенчмарки на реальних даних, наведені
у статті; синтетичні (Monte Carlo) результати та гаусівська границя
відтворюються безпосередньо з цього репозиторію (див. вище).

## Структура пакету

```
src/gsa_cpd/
    __init__.py             # Публічний API: GSADetector, BasisType, ThresholdType
    core/
        detector.py         # GSADetector -- основний клас детектора
        basis.py            # BasisType enum + обчислення базисних функцій
        threshold.py        # ThresholdType enum + пороги за PE-критерієм
        solver.py           # Розв'язувач системи FK=Y (direct / ridge / SVD)
        moments.py          # Утиліти обчислення моментів
        diagnostics.py      # GSADiagnostics dataclass
    stopping_rules/
        cusum.py            # CUSUMRule -- мінімаксне (Lorden)
        grsh.py             # GRShRule -- баєсівське (Гіршик-Рубін-Ширяєв)
        srp.py              # SRPRule -- квазімінімаксне (Ширяєв-Робертс)
    baselines/
        oracle_cusum.py     # OracleCUSUM -- верхня межа (відомий LLR)
        sign_cusum.py       # SignCUSUM -- непараметричний знаковий
        mad_cusum.py        # MADCUSUM -- MAD-нормалізований робастний
        ewma.py             # EWMA -- експоненційно зважене ковзне середнє
    data/
        schema.py           # TimeSeriesData, DatasetInfo контейнери
        preprocessing.py    # Завантаження та попередня обробка даних
    utils/
        distributions.py    # Фабрика розподілів та генерація вибірок
        metrics.py          # Метрики FAR, ADD, J(s)
experiments/
    monte_carlo/            # Розділ 4: Monte Carlo моделювання (запуск через -m)
    real_data/              # Розділ 5: лише нотатки (датасети зовнішні; див. статтю)
results/                    # Згенеровані результати (створюється при першому запуску)
tests/                      # Модульні тести
Lean/                       # Формальні докази на Lean 4 (Mathlib)
```

## Цитування

```bibtex
@article{zabolotnii2025gsa,
  title   = {Generalized Stochastic Approximation of Log-Likelihood Ratio
             for Robust Sequential Change-Point Detection},
  author  = {Zabolotnii, Serhii},
  year    = {2025},
  note    = {Software available at https://github.com/SZabolotnii/KuYuPe-Change_Point-code-supplement}
}
```

## Ліцензія

MIT License. Деталі див. у файлі [LICENSE](LICENSE).
