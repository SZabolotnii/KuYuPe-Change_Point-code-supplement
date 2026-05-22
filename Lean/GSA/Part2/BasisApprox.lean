import GSA.Part2.Setup
import GSA.Part2.Formalization
import Mathlib.Analysis.InnerProductSpace.Orthonormal
import Mathlib.Analysis.InnerProductSpace.Projection.Basic

open scoped BigOperators
open Real

namespace GSA.Part2

/-!
# 2.2. Узагальнена стохастична апроксимація LLR в базисі

Ми формалізуємо "апроксимацію в базисі" максимально абстрактно:
нехай `H` — гільбертовий простір, `φ : ℕ → H` — ортонормована система,
і `z : H` — цільовий елемент (LLR як елемент `L2` — це окремий інстанс).

Тоді часткова сума
`z_s = ∑ i < s, ⟪z, φ i⟫ • φ i`
є ортогональною проєкцією на `span {φ 0, …, φ (s-1)}`.
-/

variable {H : Type*} [NormedAddCommGroup H] [InnerProductSpace ℝ H]

/-- Коефіцієнти розкладу (узгоджено з `Part2_review.md`: `k_i`). -/
noncomputable def coeff (φ : ℕ → H) (z : H) (i : ℕ) : ℝ :=
  inner ℝ z (φ i)

/-- Часткова сума (апроксимація порядку `s`). -/
noncomputable def approx (φ : ℕ → H) (z : H) (s : ℕ) : H :=
  (Finset.range s).sum fun i => (coeff φ z i) • (φ i)

/-!
Зв'язок із формалізацією через `HilbertBasis` з `GSA.Part2.Formalization`.
Якщо `φ i := b i`, то визначення коефіцієнтів і часткових сум збігаються
з `b.repr` і `Part2.approx`.
-/

lemma coeff_hilbertBasis_eq_repr (b : HilbertBasis ℕ ℝ H) (z : H) (i : ℕ) :
    coeff (φ := fun i => b i) z i = b.repr z i := by
  have h := Part2.coeff_eq_inner (b := b) (z := z) (i := i)
  have h' : b.repr z i = inner ℝ z (b i) := by
    simpa [real_inner_comm] using h
  unfold coeff
  simp [h']

@[simp] lemma approx_hilbertBasis_eq (b : HilbertBasis ℕ ℝ H) (z : H) (s : ℕ) :
    _root_.Part2.approx b z s = approx (φ := fun i => b i) z s := by
  simp [approx, _root_.Part2.approx, coeff_hilbertBasis_eq_repr (b := b) (z := z)]

end GSA.Part2
