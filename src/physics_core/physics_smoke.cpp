// Headless check of the batch simulator: build the welded test box, drop it
// on a corner, punch it, and assert the run stays physical.
#include "physics_core/cad.hpp"
#include "physics_core/simulator.hpp"

#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <limits>

namespace {

struct Bounds {
    double lowest = std::numeric_limits<double>::max();
    double highest = -std::numeric_limits<double>::max();
};

/// Mean node height. The mesh is uniform, so this tracks the centre of mass
/// closely enough to test the free-fall trajectory without elastic ringing.
double meanHeight(const rcs::Frame& frame) {
    double sum = 0.0;
    std::size_t count = 0;
    for (const auto& part : frame.parts) {
        for (const auto& p : part.positions) {
            sum += double(p.y);
            ++count;
        }
    }
    return count ? sum / double(count) : 0.0;
}

Bounds verticalBounds(const rcs::Frame& frame) {
    Bounds bounds;
    for (const auto& part : frame.parts) {
        for (const auto& p : part.positions) {
            bounds.lowest = std::min(bounds.lowest, double(p.y));
            bounds.highest = std::max(bounds.highest, double(p.y));
        }
    }
    return bounds;
}

int failures = 0;

void check(bool condition, const char* what) {
    std::printf("%-46s %s\n", what, condition ? "ok" : "FAILED");
    if (!condition)
        ++failures;
}

/// The reported failure: a box dropped from 2.81 m, tumbling, got shoved at
/// about 1.09 s and 3.04 s while nothing was touching it. While the lowest
/// node is clear of the floor, the centre-of-mass acceleration has to stay
/// on gravity (drag on this box is a fraction of a m/s^2).
void checkNoMidairKick() {
    const rcs::AssemblyDefinition definition = rcs::makeTestCube();
    rcs::Scenario scenario;
    scenario.settings.timestep = 5.0e-4;
    scenario.initial.position = {0.0, 2.810, 0.0};
    scenario.initial.orientation = {0.350, 0.0, 0.450};
    scenario.environment.aerodynamicsEnabled = true;

    rcs::Simulator simulator;
    simulator.reset(definition, scenario);
    const double mass = simulator.assembly().mass();
    const double dt = scenario.settings.timestep;
    const chrono::ChVector3d gravity = scenario.environment.gravity;

    chrono::ChVector3d prevV{0, 0, 0};
    double t = 0.0;
    double worst = 0.0;
    double worstTime = 0.0;
    double landedAt = -1.0;
    double maxClearance = 0.0;
    bool finite = true;
    for (int step = 0; step < 7200; ++step) {
        if (!simulator.step(dt, t, scenario.environment, scenario.loads)) {
            finite = false;
            break;
        }
        t += dt;

        chrono::ChVector3d momentum{0, 0, 0};
        double minY = 1.0e9;
        for (const auto& part : simulator.assembly().parts()) {
            for (std::size_t i = 0; i < part.nodes().size(); ++i) {
                const auto& node = part.nodes()[i];
                momentum += part.nodeMass()[i] * node->GetPosDt();
                minY = std::min(minY, node->GetPos().y());
            }
        }
        if (landedAt < 0.0 && minY < 0.01)
            landedAt = t;
        // How high the lowest point gets after the first hit. A trampoline
        // sends it back up by a large fraction of the drop; concrete should
        // leave a short hop.
        if (landedAt > 0.0 && t > landedAt + 0.01)
            maxClearance = std::max(maxClearance, minY);

        const chrono::ChVector3d v = momentum / mass;
        const double stray = ((v - prevV) / dt - gravity).Length();
        // Clear of the floor: a bouncing corner sits a few millimetres up for
        // one step after contact lets go, and that step is still the bounce.
        if (t > 0.02 && minY > 0.05 && stray > worst) {
            worst = stray;
            worstTime = t;
        }
        prevV = v;
    }

    std::printf("airborne stray acceleration %.2f m/s^2 at t=%.3f s (ran to %.2f s)\n",
                worst, worstTime, t);
    std::printf("first landing at %.3f s, highest clearance afterwards %.3f m\n",
                landedAt, maxClearance);
    check(finite, "3.6 s tumble stays finite");
    check(simulator.assembly().brokenLinkCount() == 0, "sound welds survive the tumble");
    check(worst < 3.0, "no mid-air impulses while clear of the floor");
    check(landedAt > 0.0 && maxClearance < 0.40,
          "concrete floor keeps the rebound under 0.40 m");
}

}  // namespace

int main(int argc, char** argv) {
    const rcs::AssemblyDefinition definition = rcs::makeTestCube();

    const double hitForce = argc > 1 ? std::atof(argv[1]) : 6000.0;
    const double stepSize = argc > 2 ? std::atof(argv[2]) : 5.0e-4;
    const double hitRadius = argc > 3 ? std::atof(argv[3]) : 0.12;

    rcs::Scenario scenario;
    scenario.settings.duration = 1.0;
    scenario.settings.timestep = stepSize;
    scenario.settings.captureRate = 120.0;
    scenario.initial.position = {0.0, 0.55, 0.0};
    scenario.initial.orientation = {0.35, 0.0, 0.45};

    rcs::ImpactLoad hit;
    hit.name = "Hammer";
    hit.partIndex = 2;  // front panel
    hit.point = {0.0, 0.55, 0.196};
    hit.direction = {0.0, 0.0, -1.0};
    hit.force = hitForce;
    hit.radius = hitRadius;
    hit.startTime = 0.6;
    hit.duration = 0.02;
    scenario.loads.push_back(hit);

    std::printf("hit %.0f N over %.0f ms, patch r=%.3f m, dt=%.5f s\n", hitForce,
                hit.duration * 1000.0, hitRadius, stepSize);

    rcs::Simulator simulator;
    const rcs::Recording recording = simulator.run(definition, scenario, nullptr);

    std::printf("%s\n", recording.message.c_str());
    std::printf("mass %.2f kg, frames %zu, peak stress %.1f MPa, broken links %d\n",
                simulator.assembly().mass(), recording.frames.size(),
                recording.maxStress * 1.0e-6f, recording.brokenLinks);

    check(!recording.diverged, "run completed without diverging");
    check(recording.frames.size() > 100, "recorded a full playback timeline");
    check(std::abs(simulator.assembly().mass() - 19.9) < 2.0, "assembly mass is plausible");

    const rcs::Frame& first = recording.frames.front();
    const rcs::Frame& last = recording.frames.back();
    const Bounds start = verticalBounds(first);
    const Bounds end = verticalBounds(last);

    std::printf("start y [%.3f, %.3f] -> end y [%.3f, %.3f]\n", start.lowest, start.highest,
                end.lowest, end.highest);

    check(start.lowest > 0.1, "box starts clear of the floor");
    check(end.lowest > -0.02, "box never sinks through the floor");
    check(end.highest < start.highest, "box fell and settled");
    check(first.maxStress < 5.0e6f, "no pre-stress at t = 0");
    check(recording.maxStress > 1.0e6f, "impact and landing produced stress");
    check(recording.brokenLinks == 0, "sound welds survive a 0.5 m drop");

    // Free fall for the first 0.1 s should track 0.5*g*t^2; anything else
    // means gravity is being double-counted, damped, or lost in the mass matrix.
    const rcs::Frame* early = recording.sample(0.1);
    if (early) {
        const double drop = meanHeight(first) - meanHeight(*early);
        const double expected = 0.5 * 9.81 * early->time * early->time;
        std::printf("free fall at %.3f s: %.4f m (expected %.4f m)\n", early->time, drop,
                    expected);
        check(std::abs(drop - expected) < 0.05 * expected, "free fall matches 1 g");
    }

    // A 0.4 m box spans 0.69 m corner to corner, 0.57 m on an edge and 0.40 m
    // flat. Landing on a corner must start it tipping rather than balancing.
    const rcs::Frame* settled = recording.sample(recording.duration);
    if (settled) {
        const Bounds bounds = verticalBounds(*settled);
        const double span = bounds.highest - bounds.lowest;
        std::printf("vertical span %.3f m at t = %.2f s (corner 0.69, edge 0.57, flat 0.40)\n",
                    span, settled->time);
        check(span < 0.66, "box tips off its corner instead of balancing");
    }

    // Same drop with a deliberately weak bead: the failure path, the solver
    // refresh after a topology change, and progressive weld failure.
    rcs::AssemblyDefinition weak = definition;
    for (auto& weld : weak.welds)
        weld.strengthScale = 2.0e-4;

    rcs::Simulator weakRun;
    const rcs::Recording weakRecording = weakRun.run(weak, scenario, nullptr);
    std::printf("\nweak welds: %s, broken links %d, peak stress %.1f MPa\n",
                weakRecording.diverged ? "diverged" : "stable", weakRecording.brokenLinks,
                weakRecording.maxStress * 1.0e-6f);
    check(!weakRecording.diverged, "weak-weld run stays stable after breaks");
    check(weakRecording.brokenLinks > 0, "overloaded welds actually break");
    check(weakRecording.brokenLinks < 52 || weakRecording.frames.size() > 100,
          "panels keep simulating after the box comes apart");

    // A load far too sharp for the timestep must be reported, not quietly
    // turned into garbage geometry the user might trust.
    scenario.loads.front().force = 4.0e6;
    rcs::Simulator overRun;
    const rcs::Recording overRecording = overRun.run(definition, scenario, nullptr);
    std::printf("\nover-driven load: %s\n", overRecording.message.c_str());
    check(overRecording.diverged, "an unresolvable load is flagged as diverged");

    std::printf("\n");
    checkNoMidairKick();

    std::printf(failures == 0 ? "\nall checks passed\n" : "\n%d check(s) failed\n", failures);
    return failures == 0 ? 0 : 1;
}
