// Issue #8: contact / HHT / Newton robustness.
//
// The tests step the Simulator the same way Simulator::run does (reset, then
// step(dt, time, ...) with time += dt) so they can watch energy and contact
// after every user step. Nothing here is scenario-specific: the same drop is
// repeated per continuum mode.
//
// Energy is E = KE + M g h_com + contact elastic energy, which is what the
// reject gate in Simulator::step compares (e1 <= e0 + 0.02 * M g delta).
// Without loads, damping and plasticity only remove energy, so any growth
// is energy the integrator or contact injected.
//
// Accept paths today: a converged HHT step, a terminated Newton accepted at
// h_min or when "seated", and one explicit slice at h_min after HHT fails.
// All three bump hhtCommitCount() and none is counted separately, so a run
// cannot tell them apart. The target-API block at the end asks for that.

#include "physics_core/cad.hpp"
#include "physics_core/material.hpp"
#include "physics_core/scenario.hpp"
#include "physics_core/simulator.hpp"

#include <gtest/gtest.h>

#include <algorithm>
#include <cmath>
#include <functional>
#include <string>
#include <vector>

namespace {

struct Trace {
    double e0 = 0.0;
    double worstGain = 0.0;  ///< max over steps of E(t) - E(0), J
    double worstTime = 0.0;
    double startHeight = 0.0;
    double maxHeightAfterContact = -1.0e30;
    double firstContact = -1.0;
    double finalKinetic = 0.0;
    double mass = 0.0;
    bool finite = true;
    int steps = 0;
    std::vector<double> finalPositions;
};

double energyOf(const rcs::Assembly& a, double g) {
    return a.kineticEnergy() + a.mass() * g * a.massHeight() + a.contactElasticEnergy();
}

rcs::Scenario dropScene(double height, double duration) {
    rcs::Scenario scene;
    scene.environment.gravity = {0.0, -9.81, 0.0};
    scene.environment.floorEnabled = true;
    scene.environment.aerodynamicsEnabled = false;
    scene.settings.duration = duration;
    scene.settings.timestep = 1.0e-3;
    scene.settings.captureRate = 20.0;
    scene.settings.meshResolution = 0.25;
    scene.initial.position = {0.0, height, 0.0};
    scene.initial.orientation = {0.0, 0.0, 0.0};
    return scene;
}

Trace drop(const rcs::AssemblyDefinition& definition, const rcs::Scenario& scene) {
    rcs::Simulator sim;
    sim.reset(definition, scene);
    const double g = std::abs(scene.environment.gravity.y());
    Trace t;
    t.mass = sim.assembly().mass();
    t.e0 = energyOf(sim.assembly(), g);
    t.startHeight = sim.assembly().massHeight();
    const double dt = scene.settings.timestep;
    const auto steps = static_cast<long long>(std::ceil(scene.settings.duration / dt));
    double time = 0.0;
    for (long long i = 0; i < steps; ++i) {
        if (!sim.step(dt, time, scene.environment, scene.loads)) {
            t.finite = false;
            break;
        }
        time += dt;
        ++t.steps;
        const rcs::Assembly& a = sim.assembly();
        const double gain = energyOf(a, g) - t.e0;
        if (gain > t.worstGain) {
            t.worstGain = gain;
            t.worstTime = time;
        }
        if (t.firstContact < 0.0 && a.activeContactCount() > 0)
            t.firstContact = time;
        if (t.firstContact >= 0.0)
            t.maxHeightAfterContact = std::max(t.maxHeightAfterContact, a.massHeight());
    }
    t.finalKinetic = sim.assembly().kineticEnergy();
    for (const auto& part : sim.assembly().parts()) {
        for (std::size_t n = 0; n < part.nodeCount(); ++n) {
            const auto p = part.nodePosition(n);
            t.finalPositions.insert(t.finalPositions.end(), {p.x(), p.y(), p.z()});
        }
    }
    return t;
}

rcs::AssemblyDefinition boxOf(const rcs::Material& material) {
    rcs::TestCubeOptions options;
    options.material = material;
    return rcs::makeTestCube(options);
}

struct Mode {
    std::string label;
    rcs::Material material;
};

std::vector<Mode> modes() {
    return {{"ductile", rcs::aluminum6061()}, {"compliant", rcs::rubber()}, {"brittle", rcs::glass()}};
}

}  // namespace

// Goal 3 / Phase A: no artificial energy on impact. A 15 cm drop, one rule
// for every continuum mode. The allowance is 1% of the drop energy; the
// existing per-step gate allows 2% of M g delta per step, which can add up.
// MAY FAIL TODAY: this is the "energy spike on first floor hit" symptom.
// Brittle panels that shatter also exercise Assembly::kineticEnergy and
// massHeight, which add shard.mass on top of the shard nodes' own mass.
TEST(SolverPolicy, DropDoesNotInjectEnergy) {
    for (const Mode& mode : modes()) {
        const Trace t = drop(boxOf(mode.material), dropScene(0.35, 0.3));
        ASSERT_TRUE(t.finite) << mode.label << " went non-finite after " << t.steps << " steps";
        ASSERT_GT(t.firstContact, 0.0) << mode.label << " never reached the floor";
        const double dropEnergy = t.mass * 9.81 * (t.startHeight - 0.2);
        EXPECT_LE(t.worstGain, 0.01 * dropEnergy)
            << mode.label << ": energy rose by " << t.worstGain << " J at t=" << t.worstTime
            << " s (drop energy " << dropEnergy << " J, first contact at " << t.firstContact << " s)";
    }
}

// Goal 3: "stabilization that does not glue a seated body" and that does not
// launch it either. A box set down on the floor should stay there: its centre
// of mass may settle, never rise, and motion should die out.
// PASSES OR FAILS TODAY depending on the seated / ring-down heuristics
// (KE <= 0.05 M g delta, approach <= 2 mm/s, bounce <= 5 cm).
TEST(SolverPolicy, SeatedBoxStaysSeated) {
    for (const Mode& mode : modes()) {
        const Trace t = drop(boxOf(mode.material), dropScene(0.201, 0.2));
        ASSERT_TRUE(t.finite) << mode.label;
        ASSERT_GE(t.firstContact, 0.0) << mode.label << " never touched the floor";
        EXPECT_LE(t.maxHeightAfterContact - t.startHeight, 1.0e-3)
            << mode.label << ": resting box rose " << (t.maxHeightAfterContact - t.startHeight) * 1000.0 << " mm";
        EXPECT_LE(t.finalKinetic, t.mass * 9.81 * 1.0e-4)
            << mode.label << ": still moving after 0.2 s, KE=" << t.finalKinetic << " J";
    }
}

// Phase A: "reproducible metrics ... no random spike outliers". The same run
// twice on the same machine must give the same answer. (Thread count comes
// from hardware_concurrency, so this is per-machine, not across machines.)
TEST(SolverPolicy, RepeatRunsAreIdentical) {
    const rcs::Scenario scene = dropScene(0.35, 0.2);
    const Trace a = drop(rcs::makeTestCube(), scene);
    const Trace b = drop(rcs::makeTestCube(), scene);
    ASSERT_EQ(a.finalPositions.size(), b.finalPositions.size());
    double worst = 0.0;
    for (std::size_t i = 0; i < a.finalPositions.size(); ++i)
        worst = std::max(worst, std::abs(a.finalPositions[i] - b.finalPositions[i]));
    EXPECT_LE(worst, 1.0e-9) << "identical runs differ by " << worst << " m";
}

// Phase A: "inconsistent bounce across similar drops". A 0.1 mm change in
// drop height is physically nothing; the rebound must not change by more
// than a few millimetres. A large difference means a predictor or a
// threshold (onset-only contactH, seated / ring-down switches) flipped.
// MAY FAIL TODAY.
TEST(SolverPolicy, NearbyDropsBounceAlike) {
    for (const Mode& mode : modes()) {
        const Trace a = drop(boxOf(mode.material), dropScene(0.35, 0.3));
        const Trace b = drop(boxOf(mode.material), dropScene(0.3501, 0.3));
        ASSERT_TRUE(a.finite && b.finite) << mode.label;
        const double reboundA = a.maxHeightAfterContact - 0.2;
        const double reboundB = b.maxHeightAfterContact - 0.2;
        EXPECT_NEAR(reboundA, reboundB, 3.0e-3)
            << mode.label << ": rebound " << reboundA * 1000.0 << " mm vs " << reboundB * 1000.0 << " mm";
    }
}

// ---------------------------------------------------------------------------
// Target API (Goal 5 / Phase A instrumentation). Compiles once the simulator
// exposes per-run solver statistics, so it does not break the build today.
// Suggested shape:
//
//   // physics_core/solver_stats.hpp
//   struct SolverStats {
//       int converged = 0;          // HHT steps that met tolerance
//       int terminatedAccepts = 0;  // AcceptTerminatedStep commits
//       int explicitSlices = 0;     // explicit fallback commits
//       int cutbacks = 0;
//       int rejects = 0;
//       double minStep = 0.0;
//       double maxEnergyGain = 0.0;     // J, max of E(t) - E(t-1) over commits
//       double maxPenetration = 0.0;    // m
//   };
//   // Simulator gains `const SolverStats& stats() const;`
// ---------------------------------------------------------------------------
#if __has_include("physics_core/solver_stats.hpp")
#include "physics_core/solver_stats.hpp"

TEST(SolverStats, EveryCommitHasAKind) {
    rcs::Simulator sim;
    const rcs::Scenario scene = dropScene(0.35, 0.2);
    sim.run(rcs::makeTestCube(), scene, nullptr);
    const rcs::SolverStats& s = sim.stats();
    EXPECT_EQ(s.converged + s.terminatedAccepts + s.explicitSlices, sim.hhtCommitCount());
    EXPECT_GT(s.minStep, 0.0);
}

TEST(SolverStats, RestingBoxNeedsNoEmergencyAccepts) {
    rcs::Simulator sim;
    sim.run(rcs::makeTestCube(), dropScene(0.201, 0.1), nullptr);
    EXPECT_EQ(sim.stats().explicitSlices, 0);
    EXPECT_EQ(sim.stats().terminatedAccepts, 0);
}
#endif
