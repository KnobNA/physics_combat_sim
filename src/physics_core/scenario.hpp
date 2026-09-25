#pragma once

#include "chrono/core/ChVector3.h"

#include <algorithm>
#include <cmath>
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
    /// Closed arena, metres, centred on X and Z and standing on the floor.
    /// The floor grid is 1 m inside it, and the walls use the floor's contact.
    chrono::ChVector3d boundsSize{10.0, 10.0, 10.0};

    /// Grid spacing that stays readable as the arena grows. Default is 1 m.
    static double gridStep(double span) {
        if (span <= 30.0)
            return 1.0;
        if (span <= 80.0)
            return 2.0;
        return 5.0;
    }

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

    /// Still air and parts that have no other initial temperature start here.
    double ambientTemperature = 293.15;  ///< K
};

/// How force moves from its initial value to its final value.
enum class LoadRamp {
    Linear,
    Logarithmic,
    Exponential,
    Polynomial
};

/// A manually authored load: a blunt hit, a thruster, a clamp, etc. The force
/// is spread over every node within `radius` of `point` (resolved against the
/// undeformed assembly) so it behaves like a finite contact patch rather than
/// a point singularity. Over `rampTime` it moves from `initialForce` to
/// `finalForce`, then holds until `duration` ends. A ramp time of zero is a
/// step at the final force.
struct ImpactLoad {
    std::string name = "Impact";
    std::size_t partIndex = 0;
    chrono::ChVector3d point{0.0, 0.0, 0.0};
    chrono::ChVector3d direction{0.0, 0.0, -1.0};
    double initialForce = 0.0;  ///< N, total over the patch at the start
    double finalForce = 5000.0;  ///< N, total over the patch at the end of the ramp
    double rampTime = 0.0;       ///< s, zero applies finalForce for the whole duration
    LoadRamp ramp = LoadRamp::Linear;
    /// Used only for LoadRamp::Polynomial. 1 matches linear; higher stays near
    /// the initial force and then rises.
    int polynomialDegree = 2;
    double radius = 0.05;  ///< m
    double startTime = 0.0;
    double duration = 0.01;
    bool enabled = true;
};

/// Blend in [0, 1] for a ramp parameter `u` that is also in [0, 1].
inline double loadRampBlend(LoadRamp ramp, int degree, double u) {
    u = std::clamp(u, 0.0, 1.0);
    switch (ramp) {
    case LoadRamp::Logarithmic: {
        // Steep at the start, then flattening toward the final force.
        const double e = std::exp(1.0);
        return std::log(1.0 + u * (e - 1.0));
    }
    case LoadRamp::Exponential: {
        // Flat at the start, then climbing toward the final force.
        const double e = std::exp(1.0);
        return (std::exp(u) - 1.0) / (e - 1.0);
    }
    case LoadRamp::Polynomial:
        return std::pow(u, std::clamp(degree, 1, 8));
    case LoadRamp::Linear:
    default:
        return u;
    }
}

/// Largest slope of the blend, d(blend)/du, on [0, 1].
inline double loadRampMaxSlope(LoadRamp ramp, int degree) {
    switch (ramp) {
    case LoadRamp::Logarithmic:
        return std::exp(1.0) - 1.0;
    case LoadRamp::Exponential:
        return std::exp(1.0) / (std::exp(1.0) - 1.0);
    case LoadRamp::Polynomial:
        return static_cast<double>(std::clamp(degree, 1, 8));
    case LoadRamp::Linear:
    default:
        return 1.0;
    }
}

/// Patch force at `time`, in newtons. Zero outside the load's window.
inline double loadMagnitude(const ImpactLoad& load, double time) {
    if (!load.enabled || load.duration <= 0.0)
        return 0.0;
    if (time < load.startTime || time >= load.startTime + load.duration)
        return 0.0;
    const double elapsed = time - load.startTime;
    if (load.rampTime <= 1.0e-9 || elapsed >= load.rampTime)
        return load.finalForce;
    const double u = elapsed / load.rampTime;
    return load.initialForce +
           (load.finalForce - load.initialForce) * loadRampBlend(load.ramp, load.polynomialDegree, u);
}

/// Largest change in patch force from one step to the next, including the
/// switch-on. The drop at the end of the duration is not included.
inline double loadMaxStepChange(const ImpactLoad& load, double dt) {
    if (!load.enabled || load.duration <= 0.0)
        return 0.0;
    const double turnOn = std::abs(loadMagnitude(load, load.startTime));
    if (load.rampTime <= 1.0e-9)
        return std::max(turnOn, std::abs(load.finalForce));
    const double slope = std::abs(load.finalForce - load.initialForce) *
                         loadRampMaxSlope(load.ramp, load.polynomialDegree) / load.rampTime;
    return std::max(turnOn, slope * std::max(dt, 0.0));
}

/// A spherical reservoir. It has no collision. While it is on, nodes and air
/// cells inside it relax toward `temperature`; heat outside the sphere then
/// moves only by conduction and still-air diffusion.
struct HeatSource {
    std::string name = "Heat";
    chrono::ChVector3d point{0.0, 0.5, 0.0};
    double radius = 0.15;          ///< m
    double temperature = 900.0;    ///< K
    double startTime = 0.0;
    double duration = 1.0;
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
    /// 1 keeps each part's element size. Higher cuts that size down for every
    /// part, so the same knob adds tetrahedra everywhere.
    double meshResolution = 1.0;
    /// When set, a hot or loaded tetrahedron is split finer during the run.
    /// When clear, the mesh stays at the resolution chosen above.
    bool dynamicResolution = true;
};

struct Scenario {
    SimulationSettings settings;
    Environment environment;
    InitialState initial;
    std::vector<ImpactLoad> loads;
    std::vector<HeatSource> heatSources;
};

}  // namespace rcs
