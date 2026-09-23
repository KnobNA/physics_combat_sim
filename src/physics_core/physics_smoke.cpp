// Headless check of the batch simulator: build the welded test box, drop it
// on a corner, punch it, and assert the run stays physical.
#include "physics_core/cad.hpp"
#include "physics_core/simulator.hpp"

#include <algorithm>
#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <fstream>
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
        for (std::size_t i = 0; i < part.nodeCount(); ++i) {
            sum += double(part.position(i).y);
            ++count;
        }
    }
    return count ? sum / double(count) : 0.0;
}

Bounds verticalBounds(const rcs::Frame& frame) {
    Bounds bounds;
    for (const auto& part : frame.parts) {
        for (std::size_t i = 0; i < part.nodeCount(); ++i) {
            const double y = double(part.position(i).y);
            bounds.lowest = std::min(bounds.lowest, y);
            bounds.highest = std::max(bounds.highest, y);
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
            for (std::size_t i = 0; i < part.nodeCount(); ++i) {
                momentum += part.nodeMass()[i] * part.nodeVelocity(i);
                minY = std::min(minY, part.nodePosition(i).y());
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

void checkMaterialsAndImport() {
    const rcs::Material glass = rcs::glass();
    const rcs::Material rubber = rcs::rubber();
    check(glass.fracture == rcs::FractureMode::Brittle && glass.fractureStrength > 3.0e7 &&
              !glass.allowPlasticity && glass.tearStretch < 1.01,
          "glass cracks in tension before it can yield");
    check(rubber.fracture == rcs::FractureMode::Compliant && !rubber.allowPlasticity &&
              rubber.tearStretch > 3.0,
          "elastomer stretches without plasticity or shatter");

    const char* path = "import_smoke.obj";
    {
        std::ofstream file(path);
        file << "o Plate\n"
             << "v 0 0 0\nv 0.2 0 0\nv 0.2 0.1 0\nv 0 0.1 0\n"
             << "f 1 2 3\nf 1 3 4\n";
    }
    rcs::AssemblyDefinition loaded;
    std::string error;
    const bool imported = rcs::importAssembly(path, loaded, error);
    check(imported && loaded.parts.size() == 1 && !loaded.parts.front().surface.empty(),
          "OBJ import builds one triangle shell");
    check(loaded.welds.empty(), "import does not invent welds");

    rcs::Simulator preview;
    rcs::Scenario scenario;
    scenario.settings.duration = 0.05;
    scenario.settings.timestep = 5.0e-4;
    if (imported) {
        preview.reset(loaded, scenario);
        check(preview.assembly().parts().size() == 1 &&
                  !preview.assembly().parts().front().rotational(),
              "imported surface is a triangle shell");
    }
}

/// A soft yield should leave a dent (plastic strain) on a short drop without
/// shedding elements. Aluminium's real yield is far above this.
void checkPermanentDent() {
    rcs::TestCubeOptions options;
    options.material.yieldStress = 5.0e6;
    options.material.failPlasticStrain = 2.0;
    options.material.tearStretch = 1.5;
    const rcs::AssemblyDefinition definition = rcs::makeTestCube(options);

    rcs::Scenario scenario;
    scenario.settings.timestep = 5.0e-4;
    scenario.initial.position = {0.0, 0.35, 0.0};
    scenario.initial.orientation = {0.2, 0.0, 0.15};
    scenario.environment.aerodynamicsEnabled = false;

    rcs::Simulator simulator;
    simulator.reset(definition, scenario);
    bool finite = true;
    double t = 0.0;
    const double dt = scenario.settings.timestep;
    for (int step = 0; step < 800 && finite; ++step) {
        finite = simulator.step(dt, t, scenario.environment, scenario.loads);
        t += dt;
    }
    rcs::Frame frame;
    simulator.assembly().writeFrame(frame);
    float plastic = 0.0f;
    for (const auto& part : frame.parts)
        for (float value : part.plasticStrain)
            plastic = std::max(plastic, value);
    std::printf("soft drop plastic strain %.5f at t=%.3f s, torn=%d\n", plastic, t,
                simulator.assembly().parts().empty()
                    ? 0
                    : simulator.assembly().parts().front().torn() ? 1 : 0);
    check(finite, "soft-yield drop stays finite");
    check(plastic > 0.0f, "over-yield impact leaves plastic strain");
    check(!simulator.assembly().parts().front().torn(), "a dent does not tear the panel");
}

/// One clamped panel. The patch load is shared by the interior nodes so the
/// step stays inside what HHT can converge, and the boundary is what makes
/// the face bend instead of translating.
rcs::AssemblyDefinition clampedPlate(const rcs::Material& material, int divisions) {
    rcs::AssemblyDefinition definition;
    definition.name = "Plate";
    rcs::PartDefinition plate;
    plate.name = "Plate";
    plate.material = material;
    plate.thickness = 0.008;
    plate.clampBoundary = true;
    plate.patch.origin = {-0.2, 0.0, -0.2};
    plate.patch.uAxis = {1.0, 0.0, 0.0};
    plate.patch.vAxis = {0.0, 0.0, 1.0};
    plate.patch.uSize = 0.4;
    plate.patch.vSize = 0.4;
    plate.patch.divU = divisions;
    plate.patch.divV = divisions;
    definition.parts.push_back(plate);
    return definition;
}

/// A clamped glass face under a patch load must crack, and the crack has to
/// leave the panel less stressed than it was at the peak. Rubber, pushed far
/// more gently relative to its modulus, must stay in one piece.
void checkGlassShatters() {
    rcs::Scenario scenario;
    scenario.settings.timestep = 1.0e-4;
    scenario.initial.position = {0.0, 0.3, 0.0};
    scenario.environment.aerodynamicsEnabled = false;
    scenario.environment.floorEnabled = false;

    rcs::ImpactLoad hit;
    hit.name = "Punch";
    hit.partIndex = 0;
    hit.point = {0.0, 0.3, 0.0};
    hit.direction = {0.0, -1.0, 0.0};
    hit.force = 12000.0;
    hit.radius = 0.15;
    hit.startTime = 0.002;
    hit.duration = 0.003;
    scenario.loads.push_back(hit);

    rcs::Simulator simulator;
    simulator.reset(clampedPlate(rcs::glass(), 4), scenario);
    const auto advice = simulator.assembly().adviseLoad(0, scenario.loads, scenario.settings.timestep);
    std::printf("glass load on %zu nodes, %.2f m/s per step\n", advice.nodeCount, advice.deltaVPerStep);
    const std::size_t before = simulator.assembly().parts()[0].quads().size();
    double peakDuring = 0.0;
    double endTension = 0.0;
    bool finite = true;
    double t = 0.0;
    const double dt = scenario.settings.timestep;
    for (int step = 0; step < 500 && finite; ++step) {
        finite = simulator.step(dt, t, scenario.environment, scenario.loads);
        t += dt;
        const double tension = simulator.assembly().parts()[0].peakTension();
        if (simulator.assembly().parts()[0].quads().size() == before)
            peakDuring = std::max(peakDuring, tension);
        endTension = tension;
    }
    const std::size_t after = simulator.assembly().parts()[0].quads().size();
    const auto& chips = simulator.assembly().parts()[0].shards();
    double worstAspect = 1.0;
    for (const auto& shard : chips) {
        const double edgeU = (shard.local[1] - shard.local[0]).Length();
        const double edgeV = (shard.local[3] - shard.local[0]).Length();
        const double shortest = std::min(edgeU, edgeV);
        if (shortest > 1.0e-6)
            worstAspect = std::max(worstAspect, std::max(edgeU, edgeV) / shortest);
    }
    std::printf(
        "glass quads %zu -> %zu, rigid shards %zu, worst aspect %.2f, tension during %.1f MPa, end %.1f MPa at t=%.3f\n",
        before, after, chips.size(), worstAspect, peakDuring * 1.0e-6, endTension * 1.0e-6, t);
    check(finite, "glass punch stays finite");
    check(after > before, "glass face shatters into more than one piece");
    check(after > before && peakDuring > endTension, "shatter relieves the peak tension");
    check(!chips.empty(), "detached chips leave the shell solve");
    check(worstAspect <= 3.0, "rigid shards are not long needles");

    scenario.loads.front().force = 5.0;
    scenario.loads.front().duration = 0.02;
    rcs::Simulator rubberRun;
    rubberRun.reset(clampedPlate(rcs::rubber(), 2), scenario);
    const std::size_t rubberBefore = rubberRun.assembly().parts()[0].quads().size();
    finite = true;
    t = 0.0;
    for (int step = 0; step < 80 && finite; ++step) {
        finite = rubberRun.step(dt, t, scenario.environment, scenario.loads);
        t += dt;
    }
    std::printf("elastomer finite=%d quads %zu -> %zu at t=%.3f\n", finite ? 1 : 0, rubberBefore,
                rubberRun.assembly().parts()[0].quads().size(), t);
    check(finite && rubberRun.assembly().parts()[0].quads().size() == rubberBefore,
          "elastomer does not shatter");
}

/// Drop a flat panel and, once its vertical motion has died, keep going. The
/// floor spring used to unload a resting sag and then catch the fall harder.
void checkContactSettles() {
    rcs::AssemblyDefinition definition;
    definition.name = "Plate";
    rcs::PartDefinition plate;
    plate.name = "Plate";
    plate.material = rcs::aluminum6061();
    plate.thickness = 0.008;
    plate.patch.origin = {-0.2, 0.0, -0.2};
    plate.patch.uAxis = {1.0, 0.0, 0.0};
    plate.patch.vAxis = {0.0, 0.0, 1.0};
    plate.patch.uSize = 0.4;
    plate.patch.vSize = 0.4;
    plate.patch.divU = 4;
    plate.patch.divV = 4;
    definition.parts.push_back(plate);

    rcs::Scenario scenario;
    scenario.settings.timestep = 5.0e-4;
    scenario.initial.position = {0.0, 0.15, 0.0};
    scenario.environment.aerodynamicsEnabled = false;
    scenario.environment.floorEnabled = true;

    rcs::Simulator simulator;
    simulator.reset(definition, scenario);
    const double dt = scenario.settings.timestep;
    double t = 0.0;
    double settledAt = -1.0;
    double settledKe = 0.0;
    double laterPeak = 0.0;
    bool finite = true;
    auto verticalKe = [&]() {
        double energy = 0.0;
        for (const auto& part : simulator.assembly().parts()) {
            for (std::size_t i = 0; i < part.nodeCount(); ++i) {
                const double vy = part.nodeVelocity(i).y();
                energy += 0.5 * part.nodeMass()[i] * vy * vy;
            }
        }
        return energy;
    };
    for (int step = 0; step < 4000 && finite; ++step) {
        finite = simulator.step(dt, t, scenario.environment, scenario.loads);
        t += dt;
        const double energy = verticalKe();
        if (settledAt < 0.0 && t > 0.35 && energy < 0.02) {
            settledAt = t;
            settledKe = energy;
        }
        if (settledAt > 0.0 && t > settledAt + 0.05)
            laterPeak = std::max(laterPeak, energy);
        if (settledAt > 0.0 && t > settledAt + 1.0)
            break;
    }
    std::printf("settled at %.3f s with %.4f J, later peak %.4f J\n", settledAt, settledKe, laterPeak);
    check(finite, "flat drop stays finite");
    check(settledAt > 0.0, "panel vertical speed dies down");
    check(settledAt > 0.0 && laterPeak < std::max(0.05, settledKe * 4.0),
          "settled panel does not gain bounce energy");
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
    checkMaterialsAndImport();
    checkPermanentDent();
    checkGlassShatters();
    checkContactSettles();
    std::printf("\n");
    checkNoMidairKick();

    std::printf(failures == 0 ? "\nall checks passed\n" : "\n%d check(s) failed\n", failures);
    return failures == 0 ? 0 : 1;
}
