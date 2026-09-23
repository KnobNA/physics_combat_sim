#pragma once

#include "chrono/core/ChVector3.h"

#include <cstddef>
#include <string>
#include <vector>

namespace rcs {

/// Everything about the world the assembly sits in.
struct Environment {
    chrono::ChVector3d gravity{0.0, -9.81, 0.0};

    /// Sea-level standard atmosphere.
    double airDensity = 1.225;  ///< kg/m^3
    /// Pressure drag coefficient for flow normal to a flat panel.
    double dragCoefficient = 1.28;
    bool aerodynamicsEnabled = true;

    bool floorEnabled = true;
    double floorHeight = 0.0;  ///< world Y of the floor surface
    double floorHalfSize = 2.5;

    /// Penalty contact, per node. Stiffness and damping are clamped per node
    /// to whatever the timestep can integrate stably.
    ///
    /// The floor is fixed concrete: stiff, and it does not give the impact
    /// back. `floorRestitution` is the fraction of the normal impact speed
    /// returned on the way out (0 is a dead stop, 1 is elastic). Aluminum on
    /// dry concrete is about 0.15; a trampoline is near 1. Friction is the
    /// dry-concrete range.
    double contactStiffness = 1.0e6;  ///< N/m
    double contactDamping = 2.0e4;    ///< N*s/m, clamped to the stable limit
    double floorRestitution = 0.15;
    double contactFriction = 0.65;
    /// Relative tangential speed below which friction ramps in, m/s.
    double frictionRegularization = 0.02;
};

/// A manually authored load: a blunt hit, a thruster, a clamp, etc. The force
/// is spread over every node within `radius` of `point` (resolved against the
/// undeformed assembly) so it behaves like a finite contact patch rather than
/// a point singularity.
struct ImpactLoad {
    std::string name = "Impact";
    std::size_t partIndex = 0;
    chrono::ChVector3d point{0.0, 0.0, 0.0};
    chrono::ChVector3d direction{0.0, 0.0, -1.0};
    double force = 5000.0;  ///< N, total over the patch
    double radius = 0.05;   ///< m
    double startTime = 0.0;
    double duration = 0.01;
    bool enabled = true;
};

/// Initial placement / launch conditions for the whole assembly.
struct InitialState {
    chrono::ChVector3d position{0.0, 0.55, 0.0};
    /// Intrinsic XYZ Euler angles, radians.
    chrono::ChVector3d orientation{0.0, 0.0, 0.5};
    chrono::ChVector3d linearVelocity{0.0, 0.0, 0.0};
    chrono::ChVector3d angularVelocity{0.0, 0.0, 0.0};
};

struct SimulationSettings {
    double duration = 1.0;      ///< total simulated time, s
    /// Shell stiffness is integrated implicitly, but contact and user loads
    /// are applied as explicit nodal forces, so a sharp hit needs a short step.
    /// The loads panel reports the per-step velocity jump it implies.
    double timestep = 5.0e-4;   ///< integrator step, s
    double captureRate = 120.0; ///< playback frames recorded per simulated second
    /// Newton convergence tolerance. A step that cannot get inside this is
    /// aborted, which is what stops a very large applied load. Looser values
    /// accept a less exact step so that load can finish.
    double solverTolerance = 1.0e-4;
};

struct Scenario {
    SimulationSettings settings;
    Environment environment;
    InitialState initial;
    std::vector<ImpactLoad> loads;
};

}  // namespace rcs
