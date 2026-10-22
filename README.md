# Robot Combat Sim

[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](LICENSE)
[![C++17](https://img.shields.io/badge/C%2B%2B-17-00599C.svg)](https://isocpp.org/)
[![CMake](https://img.shields.io/badge/CMake-3.18+-064F8C.svg?logo=cmake)](https://cmake.org/)
[![Platform](https://img.shields.io/badge/platform-Windows%20x64-0078D6.svg?logo=windows)](https://www.microsoft.com/windows)
[![Physics](https://img.shields.io/badge/physics-Project%20Chrono%20FEA-orange.svg)](https://projectchrono.org/)
[![Solver](https://img.shields.io/badge/solver-Intel%20oneMKL%20Pardiso-0071C5.svg?logo=intel)](https://www.intel.com/content/www/us/en/developer/tools/oneapi/onemkl.html)
[![Graphics](https://img.shields.io/badge/graphics-Magnum%20%2B%20OpenGL-8A2BE2.svg)](https://magnum.graphics/)
[![Status](https://img.shields.io/badge/status-early%20development-critical.svg)](#disclaimer)

Desktop structural simulator for combat-robot style assemblies: independent parts, explicit welds and screws, continuum FEA under gravity, floor contact, impact loads, heat, and fracture. You can author a scenario, then press **Simulate**, then scrub the finished run with free camera and stress / material / temperature colouring. Pretty cool, right?

This is a batch FEA tool, not yet a full SIL (ROS / firmware) stack.

## Disclaimer

This is **very, very, *AHEM*, very early development**. The solver, contact, joints, heat, and fracture paths are still being hardened. Expect crashes, stalled runs, wrong rebound, and material behaviour that does not match a production FE code. You probably shouldn't treat a pretty playback as a validated result. Unit and smoke targets below catch some failures; they do not mean the product is reliable.

A few of the larger known holes:

- **Heated welds / torch.** A 2000 K torch on a seam is supposed to release the bead once both sides reach the solidus. Aluminium is the path that usually works. Other continuum modes (glass, wood, rubber, steel) can keep the weld, inject energy, or leave the Newton loop entirely.
- **Heat and mechanics are only loosely coupled.** Thermal strain can jump in one increment; the energy gate that rejects non-physical contact steps does not apply the same test to heat-only steps. A part that starts hot will cool, but torch-on-seam is not signed off.
- **Elastomers are incomplete.** A slender rubber patch is a layered hyperelastic solid, not a Reissner shell, so large stretch no longer wrecks the element. Rebound is in a plausible range. Compression set, DMA, Shore hardness, and several tear/stretch handbook checks are still feature gaps or failing labs.
- **Implicit seating can be very slow.** Mixed shell–solid rest (for example rubber on glass) can take many minutes for a fraction of a second of simulated time. The drop / rest policy tests pass; they are not cheap.
- **Importing CAD files** Just don't try it, I need time for that one.

## Features

- **Assemblies of parts** — each panel or imported body is its own FEA mesh, not one fused solid.
- **Manual joints** — welds (edge beads with shear / peel / cohesive opening) and screws; links can fail progressively under load. Volume nodes that occupy the same point get a zero-length tie. A shell seam keeps the rest gap the mesh already has (midsurface offset of half the thickness). Mixed shell–solid uses a frame link on the side that carries rotation.
- **Continuum materials** — aluminium 6061-T6, mild steel, AR500, soda-lime and tempered glass, rubber, hardwood. Thermal softening follows handbook-style retention curves keyed off solidus and soften temperature, not the material’s name.
- **Mesh strategy** — shell or solid from thickness over span, for every material. Reissner or BST shells for slender metal and brittle walls. TetGen fills stocky regions and finite-strain walls, including a slender rubber patch extruded as layered hyperelastic tets. `throughThickness` sets solid layer count, not the element family.
- **Contact and loads** — gravity, optional floor and arena planes, part–part node-to-face contact, blunt applied forces, optional simple aerodynamics. Deformable contact is a softened penalty with Hunt–Crossley (Flores) damping from `floorRestitution` (0.15 for concrete). Rigid chips still bounce with impulses.
- **Batch HHT integrate** — HHT (α = −0.05) with Pardiso (Intel oneMKL). First touch is limited to about a contact half-period; a seated body near the floor keeps a contact-sized step instead of throwing `h_min`. Plastic return mapping commits only on an accepted increment. Record frames, then play / loop / scrub.
- **Editor viewport** — Magnum + ImGui: select / move / rotate / splice / weld tools, stress and temperature maps, weld utilisation on playback.

## How it works

```text
CAD / test cube  →  mesh (tet / shell)  →  Chrono FEA system
                                            │
                         welds / screws as node constraints
                                            │
                         HHT + Pardiso step  →  Recording
                                            │
                         Magnum playback (stress / material / temp)
```

1. **Definition** — parts carry material, thickness, and either an authored panel, a surface mesh (OBJ/STL), or a volume mesh produced at simulate time.
2. **Meshing** — `ensurePartMesh` decides shell vs solid from thickness over span. `throughThickness` only sets how many layers a solid gets. Slender metal and brittle walls are midsurface shells. A finite-strain material has no linear shell, so a slender rubber patch is extruded as a layered hyperelastic solid rather than refused.
3. **Assembly** — coincident volume nodes share a zero-length tie. A shell seam keeps the separation the mesh already has, including a midsurface offset of half the thickness. Mixed shell–solid uses a frame link where one side carries rotation. Idealization and joints key off geometry and material fields (thickness over span, continuum mode, strengths, release rate), not the material’s name.
4. **Contact** — every live node against the floor / arena is a Chrono load with a consistent tangent. Closing uses Hunt–Crossley `1 − D v_n` with `D` sized at the latched approach speed; a Newton overshoot cannot evaluate that law at an unbounded gap or speed. Part–part pairs use the same law on a node-to-face load. Intact welds skip proximity between the two parents.
5. **Solve** — each user timestep may cut back on contact onset. Newton / HHT residuals drive rejects. A millimetre hover after bounce is still treated as contact, not free flight. A contact increment that injects more than 2% of the energy scale is rejected; heat-only increments are not held to that gate. Failed weld links disable without aborting the whole run when possible. A weld that is actually hot releases at the solidus; a link whose node has left the live mesh also releases.
6. **Record** — node poses, stress peaks, temperatures, and joint state are stored for offline scrubbing so the viewport stays interactive during the worker solve.

The built-in test cube stands in for a CAD robot: six independent shell panels, twelve welds, drop or rest on a floor. That is the reference path for “does a welded box stay together.” The welded box is a fixture of the geometry rules above, not a privileged code path. A brittle crack is a face separation inside the existing mesh: both sides stay in the FEA solve. It is not an XFEM enrichment.

## Requirements

| Dependency | Role |
| --- | --- |
| **Windows x64** + MSVC (VS 2022 / 2026 Build Tools or full IDE) | Primary supported platform |
| **CMake** ≥ 3.18 | Configure / build |
| **Intel oneMKL** (+ oneTBB on Windows) | Pardiso FEA solver; required at configure time |
| **Git** with submodule support | Chrono, Magnum, Assimp, ImGui, … |
| **TetGen sources** under `external/tetgen/` | Volume fill; AGPL, not stored in this repo |
| **OpenGL** driver | Viewport |
| **SDL2** | Fetched automatically if missing |

Linux / Ninja is partially wired in CMake but Windows + Visual Studio is the path that is exercised day to day.

Place Hang Si’s TetGen (`tetgen.h`, `tetgen.cxx`, `predicates.cxx`) in `external/tetgen/` before configure. Those files are gitignored because TetGen is AGPL; this repository does not redistribute them.

## Install and build

### 1. Clone with submodules

Chrono and Magnum are large; the first clone takes a while.

```bash
git clone --recurse-submodules https://github.com/KnobNA/robot_combat_sim.git
cd robot_combat_sim
```

If you already cloned without submodules:

```bash
git submodule update --init --recursive
```

### 2. Install Intel oneMKL

```powershell
winget install Intel.oneMKL
```

Reconfigure after installing MKL if CMake ran before that.

### 3. Configure

```powershell
cmake -S . -B build -G "Visual Studio 18 2026" -A x64
```

Use `"Visual Studio 17 2022"` if that is what you have installed.

If CMake cannot find MKL:

```powershell
cmake -S . -B build -G "Visual Studio 18 2026" -A x64 `
  -DMKL_DIR="C:/Program Files (x86)/Intel/oneAPI/mkl/latest/lib/cmake/mkl"
```

The first configure also sets up the Chrono tree under `build/chrono` with `CH_ENABLE_MODULE_PARDISO_MKL=ON`.

### 4. Build and run

```powershell
cmake --build build --config Release --target combat_sim
.\build\bin\Release\combat_sim.exe
```

MKL / TBB / Chrono / SDL2 runtime DLLs are copied next to the executable on build.

**Ninja (optional):**

```powershell
cmake -S . -B build -G Ninja -DCMAKE_BUILD_TYPE=Release
cmake --build build --target combat_sim
```

## Using the app

1. **Scenario** — duration, integrator timestep, capture rate, drop height / orientation / launch velocity, floor, arena, ambient temperature.
2. **Assembly** — start from the default welded box or import OBJ / STL; assign materials; place welds / screws with the tools.
3. **Loads** — optional blunt hits (part, point, direction, force, time window). The UI warns when a hit is sharper than the timestep can resolve.
4. **Heat** — optional point sources (torch). See the disclaimer: heated seams are not reliable yet.
5. **Simulate** — runs on a worker thread with progress; the viewport stays live.
6. **Playback** — play, loop, scrub. Orbit freely at any time. Colour by material, von Mises stress, or temperature; weld utilisation updates with the frame.

**Camera / tools**

- Right-drag: orbit · Middle-drag: pan · Scroll: zoom · Space: play / pause
- Select / Move / Rotate are exclusive tools before a run exists; editing the assembly invalidates the previous recording

Welds hold until filtered reaction / peel / opening exceeds the bead capacity for a short hold time; a seam can unzip link by link. A bead that reaches the solidus is meant to drop out the same way; that path is still buggy outside aluminium.

## Tests

```powershell
# Headless physics smoke (free fall, rest, drop, weld unzip, bad loads)
cmake --build build --config Release --target combat_sim_physics_smoke
.\build\bin\Release\combat_sim_physics_smoke.exe

# GoogleTest unit suite (mesh, welds, fracture, solver policy, …)
cmake --build build --config Release --target combat_sim_unit_tests
.\build\bin\Release\combat_sim_unit_tests.exe
```

`--only=motion` and `--only=compliant` on the smoke binary are the honest kinematics and elastomer labs. Several compliant checks are allowed to fail until those constitutive gaps close. Policy unit tests (`SolverPolicy`, `MixedJoints`, `IdealizationPolicy`, …) pin geometry and seating behaviour; do not loosen them to make a run look green.

## Repository layout

```text
robot_combat_sim/
├── CMakeLists.txt
├── LICENSE
├── README.md
├── src/
│   ├── main.cpp           # Magnum application entry
│   ├── physics_core/      # meshing, parts, joints, contact, HHT sim, recording
│   ├── render_core/       # camera, stress-shaded meshes, gizmos
│   └── editor/            # ImGui scenario, assembly, loads, playback
├── tests/unit/            # GoogleTest targets
├── tools/                 # helper scripts
└── external/              # Git submodules (Chrono, Magnum, Assimp, ImGui, …)
```

`external/tetgen/` is expected at build time and is not part of the git tree (AGPL). Chrono is a fork with FEA hooks this project needs (Reissner plastic slots, weld rest offset).

## License

This project is released under the [MIT License](LICENSE). Copyright (c) 2026 Nirav Ramdhanie.

Third-party code under `external/` (Project Chrono, Magnum, Assimp, ImGui, and their dependencies) retains its own licenses; see those trees for details. TetGen, if you add it locally, remains AGPL and is not redistributed from this repository.

## Roadmap (not in this tree yet)

- Finish heat–weld release so a torch on a seam works for every continuum mode
- Close elastomer handbook gaps (compression set, tear, DMA)
- Deeper CAD / STEP pipelines and automatic joint inference
- ROS 2 sensor / actuator traffic and Renode firmware-in-the-loop
- Broader platform support and packaging
