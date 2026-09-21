# Robot Combat Sim

Software-in-the-loop arena for designing, stressing, and fighting virtual combat robots. CAD meshes go in through Assimp/Magnum, Project Chrono runs rigid-body and FEA physics, ROS 2 carries sensor/actuator traffic, and Renode can run unmodified robot firmware.

This slice is a batch structural simulator. You author an assembly of independent parts joined by manual welds and screws, set a duration, press Simulate, and then scrub the finished run with the camera free and von Mises stress on the surfaces. The built-in test assembly is a 400 mm aluminium box built from six separate 8 mm panels welded along their twelve shared edges; it stands in until CAD import lands.

## Clone

Submodules are required. Chrono in particular is large, so the first clone takes a while.

```bash
git clone --recurse-submodules <your-github-repo-url>
```

If you already cloned without submodules:

```bash
git submodule update --init --recursive
```

## Build

Intel oneMKL is required. Chrono FEA uses its Pardiso sparse direct solver (the same path later tearing/plasticity will keep). On Windows, Pardiso is threaded with oneTBB so it does not fight Chrono's MSVC OpenMP. MUMPS + AOCL is intentionally off for now.

```powershell
winget install Intel.oneMKL
```

Configure out of source. The first configure also configures the Chrono submodule into `build/chrono` with `CH_ENABLE_MODULE_PARDISO_MKL=ON`. The first build compiles Chrono (including `Chrono_pardisomkl`), Magnum, SDL2, and the app and can take a while. Reconfigure after installing MKL if CMake ran before that.

If CMake cannot find MKL, pass the config directory:

```powershell
cmake -S . -B build -G "Visual Studio 18 2026" -A x64 `
  -DMKL_DIR="C:/Program Files (x86)/Intel/oneAPI/mkl/latest/lib/cmake/mkl"
```

**Windows (Visual Studio):**

```powershell
cmake -S . -B build -G "Visual Studio 18 2026" -A x64
cmake --build build --config Release --target combat_sim
.\build\bin\Release\combat_sim.exe
```

Use `"Visual Studio 17 2022"` if that is the generator installed on your machine.

**Ninja / single-config:**

```powershell
cmake -S . -B build -G Ninja -DCMAKE_BUILD_TYPE=Release
cmake --build build --target combat_sim
```

SDL2 is fetched automatically if it is not already installed. A working OpenGL driver is required. MKL runtime DLLs are copied next to `combat_sim.exe`.

## Using it

1. **Scenario** — set the run duration, the integrator timestep, and how many frames per simulated second to record. Set the drop height, orientation, and launch velocity.
2. **Applied loads** — add a blunt hit: pick the part, the point, the direction, the force, and the window it acts over. The panel reports the per-step velocity jump the hit implies and warns when it is sharper than the timestep can resolve.
3. **Simulate** — runs on a worker thread with a progress bar; the viewport stays live.
4. **Playback** — the bottom bar plays, loops, and scrubs the finished run. Orbit the camera freely at any point in the timeline. The stress colour map and the per-joint weld utilisation update with the frame.

Welds hold two parts together rigidly until the reaction in a bead exceeds its shear capacity, at which point that link fails; a bead can unzip link by link before the parts fully separate.

- **Select / Move / Rotate** — exclusive tools, available before a run exists. Moving a part edits the assembly definition and invalidates the previous recording.
- Right-drag: orbit. Middle-drag: pan. Scroll: zoom. Space: play/pause.

## Headless check

`combat_sim_physics_smoke` runs the same engine without a window and asserts free fall tracks 1 g, the assembly is unstressed at rest, sound welds survive a corner drop, weakened welds fail progressively without destabilising the run, and an unresolvable load is reported rather than silently producing garbage.

```powershell
cmake --build build --config Release --target combat_sim_physics_smoke
.\build\bin\Release\combat_sim_physics_smoke.exe
```

## Layout

```
robot_combat_sim/
├── CMakeLists.txt
├── src/
│   ├── main.cpp
│   ├── physics_core/     # parts, joints, contact, batch simulator, recording
│   ├── render_core/      # Magnum camera, stress-shaded meshes, gizmos
│   └── editor/           # scenario, load, and playback panels
├── firmware/             # Robot MCU code (not added yet)
├── renode/               # Virtual MCU configs (not added yet)
└── external/             # Git submodules
```

## Publish later

This project is local-only until you create a GitHub repo and push:

```bash
gh repo create robot_combat_sim --source=. --remote=origin --push
```
