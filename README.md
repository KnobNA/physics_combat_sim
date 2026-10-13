# Robot Combat Sim

[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](LICENSE)
[![C++17](https://img.shields.io/badge/C%2B%2B-17-00599C.svg)](https://isocpp.org/)
[![CMake](https://img.shields.io/badge/CMake-3.18+-064F8C.svg?logo=cmake)](https://cmake.org/)
[![Platform](https://img.shields.io/badge/platform-Windows%20x64-0078D6.svg?logo=windows)](https://www.microsoft.com/windows)
[![Physics](https://img.shields.io/badge/physics-Project%20Chrono%20FEA-orange.svg)](https://projectchrono.org/)
[![Solver](https://img.shields.io/badge/solver-Intel%20oneMKL%20Pardiso-0071C5.svg?logo=intel)](https://www.intel.com/content/www/us/en/developer/tools/oneapi/onemkl.html)
[![Graphics](https://img.shields.io/badge/graphics-Magnum%20%2B%20OpenGL-8A2BE2.svg)](https://magnum.graphics/)
[![Status](https://img.shields.io/badge/status-early%20development-critical.svg)](#status)

Desktop structural simulator for combat-robot style assemblies: independent parts, explicit welds and screws, continuum FEA under gravity, floor contact, impact loads, heat, and fracture. Author a scenario, press **Simulate**, then scrub the finished run with free camera and stress / material / temperature colouring.

This is a batch FEA tool, not a game physics sandbox and not yet a full SIL (ROS / firmware) stack.

## Status

Early development. The default scene is a **400 mm hollow aluminium box**: six 8 mm panels meshed as tetrahedra and joined by twelve edge welds. OBJ / STL import, weld tools, plasticity, brittle / ductile / elastomer paths, and shell / solid classification are in place and still being hardened. Expect rough edges; unit and smoke targets below catch many of them.

## Features

- **Assemblies of parts** — each panel or imported body is its own FEA mesh, not one fused solid.
- **Manual joints** — welds (edge beads with shear / peel / cohesive opening) and screws; links can fail progressively under load.
- **Continuum materials** — aluminium 6061-T6, mild steel, AR500, soda-lime and tempered glass, rubber, hardwood; handbook-style thermal softening where applicable.
- **Mesh strategy** — TetGen volume fills for stocky / authored solids; Reissner or BST shells for slender walls; mixed classification from local thickness vs span.
- **Contact and loads** — gravity, optional floor plane, blunt applied forces, optional simple aerodynamics.
- **Batch HHT integrate** — Pardiso (Intel oneMKL) sparse direct solve; record frames, then play / loop / scrub.
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
2. **Meshing** — `ensurePartMesh` decides shell vs solid from geometry (and authored through-thickness for ductile solids like the test panels). TetGen fills volumes; slender closed sheets stay midsurface shells.
3. **Assembly** — kissing skins get coincident node ties (`ChLinkNodeNode`); shell seams keep rest gap (`OffsetNodeWeld`); mixed shell–solid uses frame links where needed.
4. **Solve** — each user timestep may cut back on contact onset; Newton / HHT residuals drive rejects. Failed weld links disable without aborting the whole run when possible.
5. **Record** — node poses, stress peaks, temperatures, and joint state are stored for offline scrubbing so the viewport stays interactive during the worker solve.

The built-in test cube stands in for a CAD robot: six independent solids, twelve welds, drop or rest on a floor. That is the reference path for “does a welded box stay together.”

## Requirements

| Dependency | Role |
| --- | --- |
| **Windows x64** + MSVC (VS 2022 / 2026 Build Tools or full IDE) | Primary supported platform |
| **CMake** ≥ 3.18 | Configure / build |
| **Intel oneMKL** (+ oneTBB on Windows) | Pardiso FEA solver; required at configure time |
| **Git** with submodule support | Chrono, Magnum, Assimp, ImGui, … |
| **OpenGL** driver | Viewport |
| **SDL2** | Fetched automatically if missing |

Linux / Ninja is partially wired in CMake but Windows + Visual Studio is the path that is exercised day to day.

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
4. **Simulate** — runs on a worker thread with progress; the viewport stays live.
5. **Playback** — play, loop, scrub. Orbit freely at any time. Colour by material, von Mises stress, or temperature; weld utilisation updates with the frame.

**Camera / tools**

- Right-drag: orbit · Middle-drag: pan · Scroll: zoom · Space: play / pause
- Select / Move / Rotate are exclusive tools before a run exists; editing the assembly invalidates the previous recording

Welds hold until filtered reaction / peel / opening exceeds the bead capacity for a short hold time; a seam can unzip link by link.

## Tests

```powershell
# Headless physics smoke (free fall, rest, drop, weld unzip, bad loads)
cmake --build build --config Release --target combat_sim_physics_smoke
.\build\bin\Release\combat_sim_physics_smoke.exe

# GoogleTest unit suite (mesh, welds, fracture, …)
cmake --build build --config Release --target combat_sim_unit_tests
.\build\bin\Release\combat_sim_unit_tests.exe
```

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

TetGen is built as part of the project for volume meshing.

## License

This project is released under the [MIT License](LICENSE). Copyright (c) 2026 Nirav Ramdhanie.

Third-party code under `external/` (Project Chrono, Magnum, Assimp, ImGui, TetGen, and their dependencies) retains its own licenses; see those trees for details.

## Roadmap (not in this tree yet)

- Deeper CAD / STEP pipelines and automatic joint inference
- ROS 2 sensor / actuator traffic and Renode firmware-in-the-loop
- Broader platform support and packaging

Contributions and issue reports are welcome once the repo is public; please include the scenario settings and, when possible, a short clip of the failure.
