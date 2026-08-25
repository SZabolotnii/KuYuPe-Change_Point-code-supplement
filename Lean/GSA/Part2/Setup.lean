import Mathlib.Probability.Notation
import Mathlib.Probability.Moments.Variance
import Mathlib.InformationTheory.KullbackLeibler.Basic
import Mathlib.Analysis.InnerProductSpace.Projection.Basic

open scoped BigOperators
open MeasureTheory ProbabilityTheory InformationTheory

namespace GSA.Part2

/-!
# Setup (спільні означення)

Мета: мати мінімальні спільні визначення для частини 2 (`Part2_review.md`):
- дві гіпотези як два ймовірнісні мірі `μ0` і `μ1`;
- LLR як `MeasureTheory.llr μ1 μ0`;
- `klDiv` як KL-дивергенція (значення в `ℝ≥0∞`);
- `variance` як дисперсія реальної випадкової величини.
-/

variable {α : Type*} [MeasurableSpace α]

/-- Позначення очікування (Lebesgue integral) відносно міри `μ`. -/
noncomputable def Expect (μ : Measure α) (f : α → ℝ) : ℝ :=
  ∫ x, f x ∂ μ

/-- Логарифм відношення правдоподібності: `z = log (dμ1/dμ0)` (в термінах RN-похідної). -/
noncomputable def LLR (μ1 μ0 : Measure α) : α → ℝ :=
  MeasureTheory.llr μ1 μ0

/-- Симетрична дивергенція Джеффріса (в `ℝ≥0∞`).

⚠ Це означення, а не результат. **Жодна теорема цієї формалізації не пов'язує
`JeffreysENN` з інформаційним функціоналом `J`** — ні з `GSA.Part2.J`
(суми Парсеваля, `InfoFunctional.lean`), ні з `GSA.Part2.Jof`
(`Yᵀ F⁻¹ Y`, `Kernel.lean`). Ототожнення границі `J` з `JeffreysENN` —
Теорема 2(c) рукопису — хибне; правильна границя виражається через трикутну
дискримінацію. Див. `GSA/Part2/BridgeGap.lean` і
`erratum/ERRATUM_theorem2c_2026-08-23.md`. `JeffreysENN` лишається тут лише як
термінологічний орієнтир і ніде не використовується. -/
noncomputable def JeffreysENN (μ0 μ1 : Measure α) : ENNReal :=
  InformationTheory.klDiv μ1 μ0 + InformationTheory.klDiv μ0 μ1

end GSA.Part2
