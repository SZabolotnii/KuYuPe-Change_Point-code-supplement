import GSA.Part2.BasisApprox
import GSA.Part2.InfoFunctional
import Mathlib.Analysis.InnerProductSpace.l2Space
import Mathlib.Analysis.PSeries
import Mathlib.Tactic

namespace GSA.Part2

open scoped BigOperators Topology
open Filter

variable {H : Type*} [NormedAddCommGroup H] [InnerProductSpace ℝ H] [CompleteSpace H]
variable (φ : ℕ → H) (z : H)

/-!
# 2.2.2. Теорема 4 (каркас)

`Part2_review.md` формулює:

(a) ‖z - z_s‖² → 0 (збіжність в середньоквадратичному сенсі / L²);
(b) ‖z - z_s‖² = ∑_{i>s} k_i²;
(c) якщо k_i = O(i^{-r-1/2}), то ‖z - z_s‖² = O(s^{-2r}).

У гільбертовому просторі це випливає з Parseval / повноти ортонормованого базису.
У mathlib це зазвичай робиться через `HilbertBasis`/`OrthonormalBasis` + `tsum` API.

Нижче — точка входу для доведень.
-/

omit [CompleteSpace H] in
/-- Спеціалізація частини (a) на випадок гільбертового базису.
    Використовує готову збіжність `Part2.approx_tendsto` та рівність
    визначень `approx` з `Part2.approx` для `φ i := b i`. -/
theorem theorem4_a_L2_convergence_basis
    (b : HilbertBasis ℕ ℝ H) (z : H) :
    Filter.Tendsto (fun s : ℕ => ‖z - _root_.Part2.approx b z s‖^2)
      Filter.atTop (nhds 0) := by
  -- Перенесення збіжності з `Part2.approx_tendsto`.
  have happrox : Tendsto (fun s : ℕ => _root_.Part2.approx b z s) atTop (𝓝 z) :=
    _root_.Part2.approx_tendsto (b := b) (z := z)
  have hconst : Tendsto (fun _ : ℕ => z) atTop (𝓝 z) := tendsto_const_nhds
  have hsub :
      Tendsto (fun s : ℕ => z - _root_.Part2.approx b z s) atTop (𝓝 0) := by
    simpa [sub_self] using (hconst.sub happrox)
  have hnorm :
      Tendsto (fun s : ℕ => ‖z - _root_.Part2.approx b z s‖) atTop (𝓝 0) := by
    simpa using hsub.norm
  have hmul :
      Tendsto (fun s : ℕ => ‖z - _root_.Part2.approx b z s‖ *
        ‖z - _root_.Part2.approx b z s‖) atTop (𝓝 (0 * 0)) :=
    hnorm.mul hnorm
  simpa [pow_two] using hmul

theorem theorem4_a_L2_convergence
    (hφ : Orthonormal ℝ φ)
    (hcomplete : Dense ((Submodule.span ℝ (Set.range φ) : Submodule ℝ H) : Set H)) :
    Filter.Tendsto (fun s : ℕ => ‖z - approx φ z s‖^2) Filter.atTop (nhds 0) := by
  classical
  have htop :
      (Submodule.span ℝ (Set.range φ) : Submodule ℝ H).topologicalClosure = ⊤ := by
    simpa using (Submodule.dense_iff_topologicalClosure_eq_top).1 hcomplete
  have hsp :
      (⊤ : Submodule ℝ H) ≤
        (Submodule.span ℝ (Set.range φ) : Submodule ℝ H).topologicalClosure := by
    simp [htop]
  let b : HilbertBasis ℕ ℝ H := HilbertBasis.mk (v:=φ) hφ hsp
  have hb : (fun i => b i) = φ := by
    simp [b]
  simpa [hb, approx_hilbertBasis_eq (b := b) (z := z)] using
    (theorem4_a_L2_convergence_basis (b:=b) (z:=z))

theorem theorem4_b_error_as_tail_sum
    (hφ : Orthonormal ℝ φ)
    (hcomplete : Dense ((Submodule.span ℝ (Set.range φ) : Submodule ℝ H) : Set H))
    (s : ℕ) :
    ‖z - approx φ z s‖^2 = ∑' i : ℕ, (if s ≤ i then (coeff φ z i)^2 else 0) := by
  classical
  -- Побудова HilbertBasis з (φ, dense).
  have htop :
      (Submodule.span ℝ (Set.range φ) : Submodule ℝ H).topologicalClosure = ⊤ := by
    simpa using (Submodule.dense_iff_topologicalClosure_eq_top).1 hcomplete
  have hsp :
      (⊤ : Submodule ℝ H) ≤
        (Submodule.span ℝ (Set.range φ) : Submodule ℝ H).topologicalClosure := by
    simp [htop]
  let b : HilbertBasis ℕ ℝ H := HilbertBasis.mk (v:=φ) hφ hsp
  have hb : (fun i => b i) = φ := by
    simp [b]
  have hinner :
      ∀ i,
        inner ℝ (approx (φ := fun i => b i) z s) (b i) =
          if i < s then coeff (φ := fun i => b i) z i else 0 := by
    intro i
    classical
    have hsum :
        inner ℝ (approx (φ := fun i => b i) z s) (b i) =
          ∑ x ∈ Finset.range s,
            inner ℝ ((coeff (φ := fun i => b i) z x) • b x) (b i) := by
      simpa [approx] using
        (sum_inner (s := Finset.range s)
          (f := fun x => (coeff (φ := fun i => b i) z x) • b x) (x := b i))
    have hsum' :
        inner ℝ (approx (φ := fun i => b i) z s) (b i) =
          ∑ x ∈ Finset.range s,
            (coeff (φ := fun i => b i) z x) * inner ℝ (b x) (b i) := by
      simpa [inner_smul_left] using hsum
    have hsum'' :
        inner ℝ (approx (φ := fun i => b i) z s) (b i) =
          ∑ x ∈ Finset.range s,
            (coeff (φ := fun i => b i) z x) * (if x = i then 1 else 0) := by
      simpa [orthonormal_iff_ite.mp b.orthonormal] using hsum'
    have hsum''' :
        inner ℝ (approx (φ := fun i => b i) z s) (b i) =
          ∑ x ∈ Finset.range s, if x = i then coeff (φ := fun i => b i) z x else 0 := by
      simpa [mul_boole] using hsum''
    have hsum'''' :
        inner ℝ (approx (φ := fun i => b i) z s) (b i) =
          if i ∈ Finset.range s then coeff (φ := fun i => b i) z i else 0 := by
      simpa [Finset.sum_ite_eq'] using hsum'''
    simpa [Finset.mem_range] using hsum''''
  -- Працюємо в координатах базису b.
  have hcoeff :
      ∀ i, coeff (φ := fun i => b i) (z - approx (φ := fun i => b i) z s) i =
        if s ≤ i then coeff (φ := fun i => b i) z i else 0 := by
    intro i
    by_cases hi : i < s
    · have hi' : ¬ s ≤ i := not_le.mpr hi
      have hinner' :
          inner ℝ (approx (φ := fun i => b i) z s) (b i) =
            coeff (φ := fun i => b i) z i := by
        simpa [hi] using hinner i
      simp [coeff, inner_sub_left, hinner', hi']
    · have hi' : s ≤ i := le_of_not_gt hi
      have hinner' :
          inner ℝ (approx (φ := fun i => b i) z s) (b i) = 0 := by
        simpa [hi] using hinner i
      simp [coeff, inner_sub_left, hinner', hi']
  have hcoeff_sq :
      ∀ i,
        (coeff (φ := fun i => b i) (z - approx (φ := fun i => b i) z s) i)^2 =
          if s ≤ i then (coeff (φ := fun i => b i) z i)^2 else 0 := by
    intro i
    by_cases hi : s ≤ i <;> simp [hcoeff, hi]
  have hparseval :
      (∑' i : ℕ,
          (coeff (φ := fun i => b i) (z - approx (φ := fun i => b i) z s) i)^2) =
        ‖z - approx (φ := fun i => b i) z s‖^2 := by
    simpa [coeff_hilbertBasis_eq_repr] using
      (tsum_repr_sq_eq_norm_sq (b := b) (z := z - approx (φ := fun i => b i) z s))
  have hsum :
      (∑' i : ℕ,
          (coeff (φ := fun i => b i) (z - approx (φ := fun i => b i) z s) i)^2) =
        ∑' i : ℕ, (if s ≤ i then (coeff (φ := fun i => b i) z i)^2 else 0) := by
    refine tsum_congr ?_
    intro i
    simpa using hcoeff_sq i
  -- Повертаємося до початкового φ.
  simpa [hb] using hparseval.symm.trans hsum

theorem theorem4_c_rate
    (hφ : Orthonormal ℝ φ)
    (hcomplete : Dense ((Submodule.span ℝ (Set.range φ) : Submodule ℝ H) : Set H))
    (r : ℝ) (hr : 0 < r)
    (hdecay : ∃ C > 0, ∀ i, |coeff φ z i| ≤ C * ( (i+1:ℝ) ^ (-r - (1/2:ℝ)) )) :
    ∃ C' > 0, ∀ s,
      ‖z - approx φ z s‖^2 ≤
        C' * (∑' i : ℕ, if s ≤ i then ((i+1:ℝ) ^ (-r - (1/2:ℝ)))^2 else 0) := by
  rcases hdecay with ⟨C, hCpos, hdecay⟩
  refine ⟨C^2, by nlinarith, ?_⟩
  intro s
  -- Використовуємо попередню лему та оцінку коефіцієнтів.
  have htail :
      ‖z - approx φ z s‖^2 = ∑' i : ℕ, (if s ≤ i then (coeff φ z i)^2 else 0) :=
    theorem4_b_error_as_tail_sum (φ := φ) (z := z) hφ hcomplete s
  -- Порівняння поелементно.
  let f : ℕ → ℝ := fun i => if s ≤ i then (coeff φ z i)^2 else 0
  let g : ℕ → ℝ :=
    fun i => if s ≤ i then (C^2) * ((i+1:ℝ) ^ (-r - (1/2:ℝ)))^2 else 0
  have hle_point : ∀ i, f i ≤ g i := by
    intro i
    by_cases hi : s ≤ i
    · have hbound := hdecay i
      have hnonneg :
          0 ≤ C * ((i+1:ℝ) ^ (-r - (1/2:ℝ))) := by
        exact mul_nonneg (le_of_lt hCpos) (by positivity)
      have h1 : |coeff φ z i|^2 ≤ (C * ((i+1:ℝ) ^ (-r - (1/2:ℝ))))^2 := by
        have hnonneg_left : 0 ≤ |coeff φ z i| := abs_nonneg _
        have hmul := mul_le_mul hbound hbound hnonneg_left hnonneg
        simpa [pow_two] using hmul
      -- |a|^2 = a^2 для ℝ
      have hsq : (coeff φ z i)^2 ≤ (C^2) * ((i+1:ℝ) ^ (-r - (1/2:ℝ)))^2 := by
        simpa [pow_two, sq_abs, mul_comm, mul_left_comm, mul_assoc] using h1
      simpa [f, g, hi] using hsq
    · simp [f, g, hi]
  have hsum_base :
      Summable (fun i : ℕ => ((i+1:ℝ) ^ (-r - (1/2:ℝ)))^2) := by
    -- p-series for squared terms
    have hp : (-2*r - 1) < (-1:ℝ) := by linarith
    have hsum0 : Summable (fun i : ℕ => (i:ℝ) ^ (-2*r - 1)) :=
      (Real.summable_nat_rpow.mpr hp)
    have hsum' : Summable (fun i : ℕ => (i+1:ℝ) ^ (-2*r - 1)) := by
      simpa [Nat.cast_add, Nat.cast_one] using
        (summable_nat_add_iff (f := fun i : ℕ => (i:ℝ) ^ (-2*r - 1)) 1).2 hsum0
    refine hsum'.congr ?_
    intro i
    have hx : 0 ≤ (i+1:ℝ) := by positivity
    calc
      (i+1:ℝ) ^ (-2*r - 1)
          = (i+1:ℝ) ^ ((-r - (1/2:ℝ)) * 2) := by ring_nf
      _ = ((i+1:ℝ) ^ (-r - (1/2:ℝ)))^2 := by
            simpa using (Real.rpow_mul_natCast hx (-r - (1/2:ℝ)) 2)
  have hsum_g0 :
      Summable (fun i : ℕ => (C^2) * ((i+1:ℝ) ^ (-r - (1/2:ℝ)))^2) :=
    hsum_base.mul_left (C^2)
  have hnonneg_g : ∀ i, 0 ≤ g i := by
    intro i
    by_cases hi : s ≤ i
    · have h0 :
        0 ≤ (C^2) * ((i+1:ℝ) ^ (-r - (1/2:ℝ)))^2 :=
        mul_nonneg (sq_nonneg C) (sq_nonneg _)
      simpa [g, hi] using h0
    · have hgi : g i = 0 := by simp [g, hi]
      simp [hgi]
  have hsum_g : Summable g := by
    refine hsum_g0.of_nonneg_of_le hnonneg_g ?_
    intro i
    by_cases hi : s ≤ i
    · simp [g, hi]
    · have h0 :
        0 ≤ (C^2) * ((i+1:ℝ) ^ (-r - (1/2:ℝ)))^2 :=
        mul_nonneg (sq_nonneg C) (sq_nonneg _)
      simpa [g, hi] using h0
  have hnonneg_f : ∀ i, 0 ≤ f i := by
    intro i
    by_cases hi : s ≤ i
    · simpa [f, hi] using (sq_nonneg (coeff φ z i))
    · simp [f, hi]
  have hsum_f : Summable f := by
    exact hsum_g.of_nonneg_of_le hnonneg_f hle_point
  have hle_inside : (∑' i : ℕ, f i) ≤ ∑' i : ℕ, g i :=
    hsum_f.tsum_le_tsum hle_point hsum_g
  have hmul :
      (∑' i : ℕ, g i) =
        C^2 * (∑' i : ℕ,
          if s ≤ i then ((i+1:ℝ) ^ (-r - (1/2:ℝ)))^2 else 0) := by
    have hrewrite :
        (fun i : ℕ => g i) =
          fun i : ℕ => C^2 * (if s ≤ i then ((i+1:ℝ) ^ (-r - (1/2:ℝ)))^2 else 0) := by
      funext i
      by_cases hi : s ≤ i <;> simp [g, hi]
    simpa [hrewrite] using
      (tsum_mul_left (a := C^2)
        (f := fun i : ℕ => if s ≤ i then ((i+1:ℝ) ^ (-r - (1/2:ℝ)))^2 else 0))
  have hle :
      (∑' i : ℕ, f i) ≤
        C^2 * (∑' i : ℕ,
          if s ≤ i then ((i+1:ℝ) ^ (-r - (1/2:ℝ)))^2 else 0) := by
    simpa [hmul] using hle_inside
  -- Завершуємо оцінку
  simpa [f, htail] using hle

end GSA.Part2
