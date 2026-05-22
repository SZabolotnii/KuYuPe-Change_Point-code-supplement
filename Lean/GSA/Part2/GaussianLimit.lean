import GSA.Part2.BasisApprox

namespace GSA.Part2

/-!
# 2.2.5. Гаусівський ліміт при S=1

Твердження з `Part2_review.md` (ідея):
для `H0 : N(μ0,σ²)`, `H1 : N(μ1,σ²)` LLR має вигляд `ℓ(x) = a*x + b`.
Тому для базису порядку 1 (підпростір `span{1,x}`) ортогональна проєкція дає точний LLR.

Щоби не "прибиватися" до конкретного API ортогональної проєкції (яке залежить від деталей mathlib),
ми формулюємо загальну лему про будь-який проєктор `projL`, що є тотожністю на `L`.
-/

theorem projection_eq_self_of_mem
    {H : Type*} [NormedAddCommGroup H] [InnerProductSpace ℝ H]
    (L : Submodule ℝ H) (projL : H →ₗ[ℝ] H)
    (hproj : ∀ z : H, z ∈ L → projL z = z) {z : H} (hz : z ∈ L) :
    projL z = z :=
  hproj z hz

end GSA.Part2
