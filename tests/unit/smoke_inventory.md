# Smoke check inventory (honest revamp)

Units stay on `combat_sim_unit_tests`. Chrono orchestration is `combat_sim_physics_smoke`.
This suite is **test-only**: checks match the claim. They are allowed to fail until a later engine plan.

Default smoke mesh is still coarse (`meshResolution=0.25`, `dynamicResolution=off`) except:
- `editor` and motion jobs M2/M3 call `keepEditorDefaults()` → **1.0 / dynres on**, matching [`scenario.hpp`](../../src/physics_core/scenario.hpp) / [`main.cpp`](../../src/main.cpp)
- `remesh` keeps dynamic resolution on

| Group | Jobs | Kind | What would catch the 21-46-13 clips |
|---|---|---|---|
| mesh | tet fill, size field, cube/sphere/bar/tube | faithful | no |
| motion | M1 airborne 1 g; M2 hollow editor drop; M3 solid editor drop; M4 path log; M5 roll; M6 slide-to-stop; M7 rest no skate | faithful kinematics | **M2** hollow still on a corner (span ~0.69). **M3** solid ΔCOM.xz ≈ 0 while rotating |
| ductile | tension E/σy/UTS; compression; Poisson; 3pt bend; shear G; torsion; weldment tension+bend; Charpy-J; 20-cycle fatigue | faithful continuum + **feature-gap** Charpy/S–N | no |
| brittle | 3pt MOR; crush; punch (not Voronoi); thermal shock; Vickers HV | faithful Rankine + **feature-gap** HV | no |
| compliant | uniaxial λ; biaxial vs uniaxial; compression set; trouser; DMA tanδ; Shore; rebound CoR | Hencky faithful + **feature-gap** set/DMA/Shore | no |
| joints | weld throat 0.3×/1.2×; screw proof | faithful | no |
| heat | handbook retention tables; torch gradient; thermal bow | faithful | no |
| remesh | mass and COM while refining | faithful | no |
| editor | hollow + solid, keepEditorDefaults, 1 s | **editor identity** | same metrics as M2 and M3, last in the sweep |

## Handbook windows (assertion targets, not loosened)

Ductile: Al 6061-T6 E=68.9 GPa σy=276 MPa UTS=310 MPa A≈12%; A36 E=200 GPa σy=250 MPa UTS=400 MPa A≈20%; AR500 σy=1.30 GPa UTS=1.60 GPa A≈8%; G=E/2(1+ν). Charpy-V Al ~10–40 J (feature-gap).

Brittle: soda-lime MOR 40–70 MPa (preset 45), crush ~800 MPa, K_IC 0.75 MPa√m; tempered MOR ~120 MPa. HV ~500–550 (feature-gap). Shatter is Rankine bisection, not Voronoi.

Compliant: preset E=15 MPa ν=0.49 tearStretch=4 (stiff TPU analog). Compression set, DMA, Shore are feature-gaps and must fail if residual/E″/hardness cannot be produced.

## Cut from the sweep (merged, not “too hard”)

`tets-rest`, `materials`, `dent`, `malleable`, `warp`, `drop`, `alive`, `shapes`, `failure`, `midair`, `glass-punch`, `glass-impact`, `contact-settle`, `contact-solid`, `depenetrate`, `welds`, `solid`, `rigid`, `screws`, `weld-mixed`, `mixed-drop`, `chips`, `arena`, `aero`, `respan`, `loads`, `tempered`, `conserve`, `elastic`, `plastic`, `stretch`, `hits`, `thermal`, `stress`, `fracture`, `extreme`, `refine`, `friction` — claims live in the groups above. `--only=catalog` still runs the old coupon file for debugging.

## Aliases

`tets` → mesh. `catalog` → old catalog. `drop`/`alive`/`midair`/`rigid`/`friction` → motion. `glass` → brittle.
