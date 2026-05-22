import Lake
open Lake DSL

package KuYuPeChangePoint where
  srcDir := "Lean"

require mathlib from git
  "https://github.com/leanprover-community/mathlib4" @ "v4.26.0"

lean_lib KuYuPeChangePoint where
  roots := #[`KuYuPeChangePoint, `GSA]
