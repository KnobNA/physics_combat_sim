# Robot Combat Sim

Software-in-the-loop arena for designing, stressing, and fighting virtual combat robots. CAD meshes go in through Assimp/Magnum, Project Chrono runs rigid-body and FEA physics, ROS 2 carries sensor/actuator traffic, and Renode can run unmodified robot firmware.

This repo is early scaffolding: CMake and vendor submodules are in place. Application source (`src/`), firmware, and Renode configs will land next.

## Clone

Submodules are required. Chrono in particular is large, so the first clone takes a while.

```bash
git clone --recurse-submodules <your-github-repo-url>
```

If you already cloned without submodules:

```bash
git submodule update --init --recursive
```

## Layout

```
robot_combat_sim/
├── CMakeLists.txt      # Physics/render engine (Chrono, Magnum, ROS 2)
├── src/                # Chrono, Magnum, ROS bridge (not added yet)
├── firmware/           # Robot MCU code (not added yet)
├── renode/             # Virtual MCU configs (not added yet)
└── external/           # Git submodules: Chrono, Magnum, Corrade, Assimp
```

## Publish later

This project is local-only until you create a GitHub repo and push:

```bash
gh repo create robot_combat_sim --private --source=. --remote=origin --push
```

Omit `--private` if you want it public. After that, `git push` is enough for new commits.
