import GSA.Part2.Formalization
import GSA.Part2.Checks
import GSA.Part2

/-!
# Аудит аксіом

`#print axioms` для кожної теореми шару GSA. Очікуваний вивід — стандартні
аксіоми Lean/Mathlib (`propext`, `Classical.choice`, `Quot.sound`) і нічого
понад них.

Навіщо це окремим файлом. `lake build GSA` завершується з кодом 0 навіть тоді,
коли в дереві лишився `sorry`, і взагалі мовчить про те, на що теорема
спирається. Тобто «збирається» — не твердження про доведення. Цей файл робить
залежності частиною виводу збірки.

## Власних аксіом немає

Раніше шар містив одну — `GSA.Part2.hasMoment_of_higher` (монотонність
існування моментів). Її доведено в `GSA/Part2/PolynomialBasis.lean`, і тепер
кожна теорема нижче звітує рівно `[propext, Classical.choice, Quot.sound]`.

## Чиста збірка ≠ підтверджена теорема рукопису

Цей файл засвідчує, що перелічені **Lean-твердження** доведені без `sorry` і
без нестандартних аксіом. Він нічого не каже про те, чи збігається кожне з них
із однойменним твердженням рукопису. Один такий розрив відомий і задокументований:
Теорема 2 препринта ототожнює границю `J(s)` з дивергенцією Джеффріса, і це
ототожнення **хибне** — див. `GSA/Part2/BridgeGap.lean`, секція
`Інвентар NOT FORMALISED`, та `erratum/ERRATUM_theorem2c_2026-08-23.md`.

Це твердження перевіряється щоразу при збірці, а не тримається на слові:
варто комусь ввести нову аксіому чи `sorry`, і вона з'явиться у списку
теореми, яка на неї спирається.
-/

/-! ### `GSA/Part2/BasisApprox.lean` — 1 -/

#print axioms GSA.Part2.coeff_hilbertBasis_eq_repr

/-! ### `GSA/Part2/BridgeGap.lean` — 7

Місток «скінченновимірна коваріаційна система → відношення Релея».
Що саме доведено і що навмисно НЕ доведено — див. заголовок і секцію
`Інвентар NOT FORMALISED` у самому файлі. Зокрема тут **немає** теореми, яка
пов'язує границю `J` з дивергенцією Джеффріса: таке ототожнення хибне.
-/

#print axioms GSA.Part2.dot_bvec_eq_inner_combo
#print axioms GSA.Part2.rayleigh_le_Jof
#print axioms GSA.Part2.rayleigh_zero
#print axioms GSA.Part2.rayleigh_eq_Jof_of_normalEq
#print axioms GSA.Part2.isGreatest_rayleigh
#print axioms GSA.Part2.normalEq_inv_mulVec
#print axioms GSA.Part2.isGreatest_rayleigh_eq_quadForm_inv
#print axioms GSA.Part2.rayleighValue_of_triangular

/-! ### `GSA/Part2/Convergence.lean` — 4 -/

#print axioms GSA.Part2.theorem4_a_L2_convergence_basis
#print axioms GSA.Part2.theorem4_a_L2_convergence
#print axioms GSA.Part2.theorem4_b_error_as_tail_sum
#print axioms GSA.Part2.theorem4_c_rate

/-! ### `GSA/Part2/FAR_ADD.lean` — 2 -/

#print axioms GSA.Part2.exceed_le_by_chebyshev
#print axioms GSA.Part2.exceed_le_eps

/-! ### `GSA/Part2/Formalization.lean` — 3 -/

#print axioms Part2.ku1_step3
#print axioms Part2.approx_tendsto
#print axioms Part2.coeff_eq_inner

/-! ### `GSA/Part2/GaussianLimit.lean` — 1 -/

#print axioms GSA.Part2.projection_eq_self_of_mem

/-! ### `GSA/Part2/GramSchmidt.lean` — 5 -/

#print axioms GSA.Part2.empInner_eq_empVec
#print axioms GSA.Part2.empInner_eq_empInnerVec
#print axioms GSA.Part2.empInnerVec_smul_smul
#print axioms GSA.Part2.empGramSchmidt_orthonormal
#print axioms GSA.Part2.empGramSchmidtFin_orthonormal

/-! ### `GSA/Part2/InfoFunctional.lean` — 4

⚠ Імена `theorem2_*` названо за Теоремою 2 рукопису, але доведено рівно
гільбертові факти про часткові суми Парсеваля по **ортонормованому** базису:
межа зверху і границя — це `‖z‖²`, а не дивергенція Джеффріса. Див. заголовок
`InfoFunctional.lean` і `BridgeGap.lean`.
-/

#print axioms GSA.Part2.tsum_repr_sq_eq_norm_sq
#print axioms GSA.Part2.theorem2_a_upper_bound
#print axioms GSA.Part2.theorem2_b_monotone
#print axioms GSA.Part2.theorem2_c_tendsto

/-! ### `GSA/Part2/Kernel.lean` — 29 -/

#print axioms GSA.Part2.gram_symm
#print axioms GSA.Part2.combo_sub
#print axioms GSA.Part2.combo_add
#print axioms GSA.Part2.combo_smul
#print axioms GSA.Part2.inner_combo_left
#print axioms GSA.Part2.inner_v_combo
#print axioms GSA.Part2.normalEq_iff_orthogonal
#print axioms GSA.Part2.inner_combo_resid_eq_zero
#print axioms GSA.Part2.risk_eq_risk_add_normSq
#print axioms GSA.Part2.normSq_combo_eq_quadForm
#print axioms GSA.Part2.normalEq_isMinimizer
#print axioms GSA.Part2.risk_eq_of_normalEq_of_normalEq
#print axioms GSA.Part2.Jof_le_normSq
#print axioms GSA.Part2.Jof_eq_normSq_combo
#print axioms GSA.Part2.Jof_eq_dot_of_normalEq
#print axioms GSA.Part2.existsUnique_normalEq_of_posDef
#print axioms GSA.Part2.exists_isMinimizer_of_posDef
#print axioms GSA.Part2.combo_snoc_zero
#print axioms GSA.Part2.risk_snoc_zero
#print axioms GSA.Part2.Jof_init_le_Jof_of_isMinimizer
#print axioms GSA.Part2.risk_eq_of_range_eq
#print axioms GSA.Part2.combo_transform
#print axioms GSA.Part2.range_combo_eq_of_isUnit
#print axioms GSA.Part2.isMinimizer_transfer
#print axioms GSA.Part2.combo_snoc
#print axioms GSA.Part2.combo_snoc_apply
#print axioms GSA.Part2.Jof_snoc_orthogonal
#print axioms GSA.Part2.range_combo_snoc_resid
#print axioms GSA.Part2.Jof_snoc_general

/-! ### `GSA/Part2/Kunchenko.lean` — 7 -/

#print axioms GSA.Part2.KU1_nonneg
#print axioms GSA.Part2.KU1_scale
#print axioms GSA.Part2.psiOpt_pos_iff_lr_gt_one
#print axioms GSA.Part2.psiOpt_neg_iff_lr_lt_one
#print axioms GSA.Part2.KU1_optimal_direction
#print axioms GSA.Part2.theorem1_decision_rule_equivalence
#print axioms GSA.Part2.KU1_linear_minimizer

/-! ### `GSA/Part2/LinearSystem.lean` — 1 -/

#print axioms GSA.Part2.has_solution_FK_eq_Y

/-! ### `GSA/Part2/PolynomialBasis.lean` — 1 -/

#print axioms GSA.Part2.hasMoment_of_higher

/-! ### `GSA/Part2/RelativeDisorder.lean` — 1 -/

#print axioms GSA.Part2.relDisorder_nonneg

/-! ### `GSA/Part2/RobustPE.lean` — 6 -/

#print axioms GSA.Part2.theorem3_PE_asymptotic_equivalence_NP
#print axioms GSA.Part2.theorem5_asymptotic_normality_criterion_Y
#print axioms GSA.Part2.PE_threshold_difference
#print axioms GSA.Part2.theorem3_PE_NP_threshold_equivalence
#print axioms GSA.Part2.criterion_Y_nonneg
#print axioms GSA.Part2.KU1_le_criterion_Y

/-! ### `GSA/Part2/UnimodalBounds.lean` — 7 -/

#print axioms GSA.Part2.cantelli_bound_nonneg
#print axioms GSA.Part2.cantelli_bound_le_one
#print axioms GSA.Part2.cantelli_lt_chebyshev
#print axioms GSA.Part2.VP_bound_nonneg
#print axioms GSA.Part2.VP_lt_chebyshev
#print axioms GSA.Part2.PE_threshold_VP_lt_PE
#print axioms GSA.Part2.VP_threshold_applicability_simple

/-! ### `GSA/Part2/MultiAlternative.lean` — 7 -/

#print axioms GSA.Part2.MAry.Ypair_add
#print axioms GSA.Part2.MAry.Ksol_add
#print axioms GSA.Part2.MAry.Fpair_const
#print axioms GSA.Part2.MAry.dot_mulVec_comm
#print axioms GSA.Part2.MAry.lam_eq_half_g_sub
#print axioms GSA.Part2.MAry.transitive_of_potential
#print axioms GSA.Part2.MAry.exists_intransitive

/-! ### `GSA/Part2/MultiAlternativePE.lean` — 5 -/

#print axioms GSA.Part2.MAry.copeland_of_wins_all
#print axioms GSA.Part2.MAry.error_subset_duels
#print axioms GSA.Part2.MAry.PE_duel_bound
#print axioms GSA.Part2.MAry.PE_error_bound
#print axioms GSA.Part2.MAry.PE_error_bound_uniform
