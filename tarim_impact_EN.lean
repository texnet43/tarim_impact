/-
  Tarim Impact Hypothesis — formalization skeleton (Lean 4 + Mathlib).
  Machine-checkable statements of the axis mechanics, sections 10.1-10.3.
  `sorry` marks open proof obligations; the numeric constants are from
  the preprint (section 2). Build:  lake build  (requires elan).
  License: CC BY 4.0, texnet43
-/

import Mathlib

namespace TarimImpact

/-- Angular momentum of the impactor: L = m v R sin(theta). -/
noncomputable def LImpactor (m v R theta : ℝ) : ℝ := m * v * R * Real.sin theta

/-- Earth's angular momentum: L = C omega. -/
noncomputable def LEarth (C omega : ℝ) : ℝ := C * omega

/-- TPW inertia requirement of section 10.2: a figure rotation by eps
    costs dI = (C - A) tan(eps) — no angular momentum. -/
noncomputable def dIRequired (C A eps : ℝ) : ℝ := (C - A) * Real.tan eps

/-- The "188" bound of section 10.2: rotating the momentum vector by the
    full obliquity requires L_E sin(eps) — about 188 impactor momenta. -/
noncomputable def boundFullRotation (LE eps : ℝ) : ℝ := LE * Real.sin eps

/-- Conservation ceiling (section 10.2): the maximal tilt of the total
    spin axis attainable by ANY internal redistribution is bounded by
    arctan(L_imp / L_E); with the preprint constants it is below
    2.2e-3 rad (~0.125 deg). Open proof obligation. -/
theorem ceiling_internal
    (Limp LE : ℝ) (hpos : 0 < Limp ∧ 0 < LE)
    (himp : Limp = 1.26e31) (hle : LE = 5.86e33) :
    Real.arctan (Limp / LE) < 2.2e-3 := by
  sorry

/-- Inner-core tier (section 10.2): the impactor momentum exceeds the
    inner core's spin angular momentum (~4.3e30 kg m²/s) by a factor
    of 2-3, so a full re-orientation of its spin is energetically
    available. Open proof obligation. -/
theorem inner_core_tier
    (Limp LIC : ℝ) (himp : Limp = 1.26e31) (hlic : 4.28e30 ≤ LIC)
    (hle : LIC ≤ 5.4e30) :
    2 ≤ Limp / LIC ∧ Limp / LIC ≤ 3.5 := by
  sorry

end TarimImpact
