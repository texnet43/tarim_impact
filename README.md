# Tarim Impact Hypothesis — Impact Paradigm

**Author:** texnet43
**Version:** 1.0 (November 2026)
**License:** Creative Commons Attribution 4.0 International (CC BY 4.0)
**DOI:** [10.5281/zenodo.24405368](https://doi.org/10.5281/zenodo.24405368)
**Contact:** texnet43@gmail.com · [@texnet43](https://t.me/texnet43)

---

## What is this?

The **Tarim Impact Hypothesis** proposes that the impact of a 300-km
iron-nickel body on the Tarim Basin (39°N, 83°E) 12,800 cal BP (the
Allerød/Younger Dryas boundary; the GISP2 platinum spike dated against
annually layered ice) at 37° against Earth's rotation, with an energy of
5.0×10²⁸ J, is the **primary event that shaped modern Earth**: the fracture
of a single lithosphere into radially diverging plates, mountain systems as
the raised edge of the wound, the mid-ocean ridge system as a re-accretion
seam, Chicxulub as an antipodal structure, the LLSVP pair, the South
Atlantic Anomaly, and the geochemical boundary layer of the Younger Dryas.

The work is presented as an **alternative scientific paradigm**: where the
proposed interpretation conflicts with the standard model, the conflict is
stated **explicitly as an open question** (§17), with 13 falsifiable
predictions (§18).

## Key canonical numbers

| Parameter | Value |
|---|---|
| Impactor | Fe-Ni (84% Fe, 16% Ni), ⌀ 300 km, 30 km/s, 37° against rotation |
| Mass / energy | 1.113×10²⁰ kg (0.002% of Earth) / 5.0×10²⁸ J |
| Event date (markers) | 12,800 cal BP (GISP2 Pt peak, annual layers) |
| Plate relaxation τ | 462 years; today ~8 cm/yr |
| Instantaneous tilt impulse | 0.125°; the rest of the reorientation is TPW |
| Obliquity | 24.04° at 12.8 ka → 23.44° (Milankovitch, pre-existing) |
| Wave deployment | 90° tangential circle, ~12 km crust spall, peel line 95–105° (§4.2–4.4) |
| Antipodal focusing | Chicxulub in the 84–92°W band; gain ~4×10⁵× |
| Geochemistry | ¹⁸⁷Os/¹⁸⁸Os = 0.53 (~57% extraterrestrial Os); Pt ×100, superchondritic Pt/Ir |

## Deposit contents

| File | Content |
|---|---|
| `tarim_impact_RU.txt` / `tarim_impact_EN.txt` | Full parallel texts (primary source; 19 sections, 94 sources) |
| `tarim_impact_RU.md` / `.html` / `.docx` (+EN) | Derived formats |
| `tarim_impact_RU.pdf` / `tarim_impact_EN.pdf` | Print versions (from .docx) |
| `tarim_impact_RU.png` / `tarim_impact_EN.png` | Title previews |
| `tarim_impact_RU.py` / `tarim_impact_EN.py` | Self-contained verification bundles: 9 checks embedded, no external dependencies |
| `tarim_impact_RU.lean` / `tarim_impact_EN.lean` | Lean 4 (Mathlib) formalization skeletons |
| `README.md` | This overview |

## Verification

All checks are reproducible without network access:

- `tarim_impact_RU.py` / `tarim_impact_EN.py` — self-contained bundles,
  9/9 checks passed, exit 0 (run: `python tarim_impact_RU.py`)
- `verify_tarim` — 24 numerical checks (mass/energy/geometry/osmium/axis)
- `axis_mechanism_check` — 17 checks (TPW budget, conservation ceiling,
  inner-core fingerprint, Berger obliquity 24.04° at 12.8 ka)
- `all_tests` — 11 scripts passed; 14 of 14 key numerical predictions confirmed
- Cross-reference consistency — 94 literature entries in both languages,
  all cited in the text, continuous numbering

## Honest scope

The preprint explicitly separates the three axis-related processes
(precession 25,771 yr — always existed; TPW — damped, impact-initiated;
Milankovitch obliquity — not reachable internally), reports the conflicts
with the standard model as open questions (§17: Atlantic floor ages,
LLSVP age, impactor identity), and lists the data state of every
still-open test (§18).

## Keywords

impact paradigm; Tarim Basin; lithospheric plates; antipodal focusing;
osmium anomaly; platinum anomaly; PGE; impact glass; Younger Dryas;
wave deployment; true polar wander; inner core; Fe-Ni impactor
