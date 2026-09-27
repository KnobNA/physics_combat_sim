// Headless check of the batch simulator: build the welded test box, drop it
// on a corner, punch it, and assert the run stays physical.
#include "physics_core/cad.hpp"
#include "physics_core/mesh_split.hpp"
#include "physics_core/simulator.hpp"

#include <algorithm>
#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <fstream>
#include <string>
#include <limits>
#include <tuple>

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

/// A node held by the world, or one that has been frozen because its material
/// left the solve. Either way it carries no momentum.
bool nodeIsFixed(const rcs::Part& part, std::size_t index) {
    if (index < part.xyzNodes().size() && part.xyzNodes()[index])
        return part.xyzNodes()[index]->IsFixed();
    if (index < part.nodes().size() && part.nodes()[index])
        return part.nodes()[index]->IsFixed();
    return true;
}

int failures = 0;

void check(bool condition, const char* what) {
    std::printf("%-46s %s\n", what, condition ? "ok" : "FAILED");
    std::fflush(stdout);
    if (!condition)
        ++failures;
}

/// Finite is required but never enough: every physics case also needs a band.
void expectFinite(bool finite, const char* what) {
    check(finite, what);
}

/// Pass only when `value` sits in [lo, hi]. Prints the miss so a green run
/// cannot hide a clamped speed, a 48-tet peel, or a retuned yield.
void expectRange(double value, double lo, double hi, const char* what) {
    const bool ok = std::isfinite(value) && std::isfinite(lo) && std::isfinite(hi) &&
                    value >= lo && value <= hi;
    if (ok) {
        std::printf("%-46s ok\n", what);
    } else {
        std::printf("%-46s FAILED  (got %.6g, want [%.6g, %.6g])\n", what, value, lo, hi);
        ++failures;
    }
    std::fflush(stdout);
}

void expectRange(double value, double lo, double hi, const std::string& what) {
    expectRange(value, lo, hi, what.c_str());
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
    scenario.environment.aerodynamicsEnabled = std::getenv("RCS_NO_AERO") == nullptr;

    rcs::Simulator simulator;
    simulator.reset(definition, scenario);
    const double mass = simulator.assembly().mass();
    const double dt = scenario.settings.timestep;
    const chrono::ChVector3d gravity = scenario.environment.gravity;

    chrono::ChVector3d prevV{0, 0, 0};
    double prevMass = -1.0;
    int logged = 0;
    double peakStress = 0.0;
    double peakPlastic = 0.0;
    double t = 0.0;
    double worst = 0.0;
    double worstTime = 0.0;
    double worstClearance = 0.0;
    double worstSpeed = 0.0;
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
        double moving = 0.0;
        double minY = 1.0e9;
        for (const auto& part : simulator.assembly().parts()) {
            for (std::size_t i = 0; i < part.nodeCount(); ++i) {
                if (nodeIsFixed(part, i))
                    continue;
                momentum += part.nodeMass()[i] * part.nodeVelocity(i);
                moving += part.nodeMass()[i];
                minY = std::min(minY, part.nodePosition(i).y());
            }
            for (const auto& shard : part.shards()) {
                momentum += shard.mass * shard.velocity;
                moving += shard.mass;
            }
            std::vector<float> stress;
            part.evaluateStress(stress);
            for (float value : stress)
                peakStress = std::max(peakStress, static_cast<double>(value));
            std::vector<float> plastic;
            part.evaluatePlasticStrain(plastic);
            for (float value : plastic)
                peakPlastic = std::max(peakPlastic, static_cast<double>(value));
        }
        if (landedAt < 0.0 && minY < 0.01)
            landedAt = t;
        // How high the lowest point gets after the first hit. A trampoline
        // sends it back up by a large fraction of the drop; concrete should
        // leave a short hop.
        if (landedAt > 0.0 && t > landedAt + 0.01)
            maxClearance = std::max(maxClearance, minY);

        const chrono::ChVector3d v = moving > 0.0 ? momentum / moving : chrono::ChVector3d(0, 0, 0);
        const double stray = ((v - prevV) / dt - gravity).Length();
        // Clear of the floor: a bouncing corner sits a few millimetres up for
        // one step after contact lets go, and that step is still the bounce.
        // Steps where a chip left the mesh are skipped, because the sum then
        // covers different material and the change is not an impulse.
        if (t > 0.02 && minY > 0.05 && std::abs(moving - prevMass) < 1.0e-9 && stray > worst) {
            worst = stray;
            worstTime = t;
            worstClearance = minY;
            worstSpeed = v.Length();
        }
        if (std::getenv("RCS_TRACE") && t > 0.02 && stray > 100.0 && logged < 12) {
            ++logged;
            std::size_t nodes = 0;
            std::size_t shards = 0;
            for (const auto& part : simulator.assembly().parts()) {
                nodes += part.nodeCount();
                shards += part.shards().size();
            }
            std::printf("  stray %.0f m/s^2 t=%.4f minY=%.4f mass %.6f->%.6f nodes %zu shards %zu "
                        "v (%.3f %.3f %.3f) -> (%.3f %.3f %.3f)\n",
                        stray, t, minY, prevMass, moving, nodes, shards, prevV.x(), prevV.y(),
                        prevV.z(), v.x(), v.y(), v.z());
        }
        prevV = v;
        prevMass = moving;
    }

    std::printf("airborne stray acceleration %.2f m/s^2 at t=%.3f s, clearance %.3f m, "
                "speed %.2f m/s (ran to %.2f s)\n",
                worst, worstTime, worstClearance, worstSpeed, t);
    std::printf("first landing at %.3f s, highest clearance afterwards %.3f m\n",
                landedAt, maxClearance);
    std::printf("peak stress %.1f MPa, peak plastic strain %.4f, moving mass %.2f of %.2f kg\n",
                peakStress * 1.0e-6, peakPlastic, prevMass, mass);
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
    check(rcs::rankineCompressionLimit(glass) > 5.0 * glass.fractureStrength,
          "glass crush strength is well above its tensile strength");
    check(rubber.fracture == rcs::FractureMode::Compliant && !rubber.allowPlasticity &&
              rubber.tearStretch > 3.0,
          "elastomer stretches without plasticity or shatter");

    const rcs::Material aluminum = rcs::aluminum6061();
    const rcs::Material steel = rcs::mildSteel();
    check(rcs::thermalStrengthRetention(aluminum, 293.15) > 0.98,
          "room-temperature aluminum keeps its strength");
    check(rcs::thermalStrengthRetention(aluminum, 533.0) < 0.20 &&
              rcs::thermalStrengthRetention(aluminum, 533.0) > 0.08,
          "6061 yield has collapsed by 260 C");
    check(rcs::thermalStiffnessRetention(aluminum, 573.0) > 0.70 &&
              rcs::thermalStiffnessRetention(aluminum, 573.0) < 0.85,
          "aluminum modulus at 300 C follows the handbook");
    check(rcs::thermalStrengthRetention(aluminum, 3000.0) <= 2.0e-4 &&
              rcs::thermalStiffnessRetention(aluminum, 3000.0) <= 2.0e-4,
          "aluminum far above the solidus has essentially no solid strength");
    check(rcs::thermalDuctilityRetention(aluminum, 293.15) > 0.99, "room-temperature aluminum keeps its ductility");
    check(rcs::thermalDuctilityRetention(aluminum, 3000.0) < 0.08,
          "molten aluminum tears at a small mechanical strain");
    check(rcs::thermalStrengthRetention(steel, 673.0) > 0.95, "mild steel yield holds through 400 C");
    check(rcs::thermalStrengthRetention(steel, 873.0) < 0.55 && rcs::thermalStrengthRetention(steel, 873.0) > 0.35,
          "mild steel yield has fallen by 600 C");

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
    check(imported && loaded.parts.size() == 1 && !loaded.parts.front().source.empty(),
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

rcs::AssemblyDefinition clampedPlate(const rcs::Material& material, int divisions);

/// 6061-T6 (MMPDS / ASTM B209): σ_y = 276 MPa, elongation ~12%. The clamp is
/// the min-X edge, so this is a 200 mm cantilever, not a four-sided plate.
/// First-yield for a centre fibre is about 6 kN; 15 kN is ~2.5× that, so the
/// panel must dent and spring back only part way, without tearing.
void checkPermanentDent() {
    const rcs::Material aluminum = rcs::aluminum6061();
    rcs::Scenario scenario;
    scenario.settings.timestep = 2.0e-4;
    scenario.initial.position = {0.0, 0.3, 0.0};
    scenario.environment.gravity = {0.0, 0.0, 0.0};
    scenario.environment.floorEnabled = false;
    scenario.environment.aerodynamicsEnabled = false;
    rcs::ImpactLoad hit;
    hit.partIndex = 0;
    hit.point = {0.0, 0.3, 0.0};
    hit.direction = {0.0, -1.0, 0.0};
    hit.initialForce = 70000.0;
    hit.finalForce = 70000.0;
    hit.radius = 0.03;
    hit.duration = 0.02;
    scenario.loads.push_back(hit);

    rcs::Simulator simulator;
    simulator.reset(clampedPlate(aluminum, 6), scenario);
    rcs::Part& plate = simulator.assembly().parts().front();
    std::size_t centre = 0;
    double best = 1.0e9;
    for (std::size_t n = 0; n < plate.nodeCount(); ++n) {
        const double d = (plate.nodePosition(n) - hit.point).Length();
        if (d < best) {
            best = d;
            centre = n;
        }
    }
    const double startY = plate.nodePosition(centre).y();
    bool finite = true;
    double t = 0.0;
    const double dt = scenario.settings.timestep;
    double peakDent = 0.0;
    float plastic = 0.0f;
    for (int step = 0; step < 400 && finite; ++step) {
        finite = simulator.step(dt, t, scenario.environment, scenario.loads);
        t += dt;
        const rcs::Part& live = simulator.assembly().parts().front();
        peakDent = std::max(peakDent, startY - live.nodePosition(centre).y());
        std::vector<float> field;
        live.evaluatePlasticStrain(field);
        for (float value : field)
            plastic = std::max(plastic, value);
    }
    const rcs::Part& done = simulator.assembly().parts().front();
    const double residual = startY - done.nodePosition(centre).y();
    const double residualRatio = peakDent > 1.0e-9 ? residual / peakDent : 0.0;
    std::printf("6061 dent peak %.4f m residual %.4f m plastic %.5f torn=%d shards %zu\n",
                peakDent, residual, plastic, done.torn() ? 1 : 0, done.shards().size());

    expectFinite(finite, "6061 punch stays in the solver");
    expectRange(plastic, 0.002, 0.12, "6061 punch yields (MMPDS 276 MPa)");
    expectRange(peakDent, 2.0e-4, 0.08, "6061 punch leaves a millimetre-scale dent");
    expectRange(residualRatio, 0.30, 1.0, "6061 dent keeps at least 30% after unload");
    check(!done.torn() && done.shards().empty(), "6061 dent does not tear the panel");
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

/// A clamped glass face under a patch load must crack once Rankine tension
/// (~45 MPa, soda-lime handbook) is exceeded, and the crack has to leave the
/// panel less stressed than it was at the peak. Rubber under a light push
/// relative to its modulus must stay in one piece with no plastic strain.
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
    hit.initialForce = 60000.0;
    hit.finalForce = 60000.0;
    hit.radius = 0.15;
    hit.startTime = 0.002;
    hit.duration = 0.003;
    scenario.loads.push_back(hit);

    rcs::Simulator simulator;
    simulator.reset(clampedPlate(rcs::glass(), 4), scenario);
    const auto advice = simulator.assembly().adviseLoad(0, scenario.loads, scenario.settings.timestep);
    std::printf("glass load on %zu nodes, %.2f m/s per step\n", advice.nodeCount, advice.deltaVPerStep);
    const int before = simulator.assembly().parts()[0].liveTetCount();
    double peakDuring = 0.0;
    double endTension = 0.0;
    bool finite = true;
    double t = 0.0;
    const double dt = scenario.settings.timestep;
    for (int step = 0; step < 500 && finite; ++step) {
        finite = simulator.step(dt, t, scenario.environment, scenario.loads);
        t += dt;
        const double tension = simulator.assembly().parts()[0].peakTension();
        peakDuring = std::max(peakDuring, tension);
        endTension = tension;
    }
    const rcs::Part& glassPart = simulator.assembly().parts()[0];
    const int after = glassPart.liveTetCount();
    const auto& chips = glassPart.shards();
    double worstAspect = 1.0;
    for (const auto& shard : chips) {
        const double edgeU = (shard.local[1] - shard.local[0]).Length();
        const double edgeV = (shard.local[3] - shard.local[0]).Length();
        const double shortest = std::min(edgeU, edgeV);
        if (shortest > 1.0e-6)
            worstAspect = std::max(worstAspect, std::max(edgeU, edgeV) / shortest);
    }
    std::printf(
        "glass tets %d -> %d, rigid shards %zu, worst aspect %.2f, tension during %.1f MPa, end %.1f MPa at t=%.3f\n",
        before, after, chips.size(), worstAspect, peakDuring * 1.0e-6, endTension * 1.0e-6, t);

    expectFinite(finite, "glass punch stays in the solver");
    expectRange(peakDuring, 0.8 * rcs::glass().fractureStrength, 5.0 * rcs::glass().fractureStrength,
                "glass punch reaches Rankine tension");
    check(!chips.empty() || (before > 0 && after < before), "glass face shatters into more than one piece");
    expectRange(endTension / std::max(peakDuring, 1.0), 0.0, 0.70, "shatter relieves the peak tension");
    check(!chips.empty(), "detached chips leave the shell solve");
    expectRange(worstAspect, 1.0, 20.0, "glass shards are not long needles");
    if (before > 0)
        expectRange(double(after) / double(before), 0.05, 0.95, "glass punch does not delete the whole panel");

    scenario.loads.front().initialForce = 5.0;
    scenario.loads.front().finalForce = 5.0;
    scenario.loads.front().duration = 0.02;
    rcs::Simulator rubberRun;
    rubberRun.reset(clampedPlate(rcs::rubber(), 2), scenario);
    const int rubberBefore = rubberRun.assembly().parts()[0].liveTetCount();
    finite = true;
    t = 0.0;
    float rubberPlastic = 0.0f;
    for (int step = 0; step < 80 && finite; ++step) {
        finite = rubberRun.step(dt, t, scenario.environment, scenario.loads);
        t += dt;
        std::vector<float> field;
        rubberRun.assembly().parts()[0].evaluatePlasticStrain(field);
        for (float value : field)
            rubberPlastic = std::max(rubberPlastic, value);
    }
    const rcs::Part& rubberDone = rubberRun.assembly().parts()[0];
    std::printf("elastomer finite=%d tets %d -> %d plastic %.5f shards %zu at t=%.3f\n", finite ? 1 : 0, rubberBefore,
                rubberDone.liveTetCount(), rubberPlastic, rubberDone.shards().size(), t);
    expectFinite(finite, "elastomer punch stays in the solver");
    check(rubberDone.liveTetCount() == rubberBefore && rubberDone.shards().empty(),
          "elastomer does not shatter");
    expectRange(rubberPlastic, 0.0, 1.0e-4, "elastomer does not yield");
}

/// High-speed floor landings. Metals must yield (MMPDS 6061-T6 276 MPa).
/// Glass must fragment once Rankine tension is exceeded. Rubber must stretch.
/// Finite is required; a 48-tet peel or a 15 m/s speed clamp is not a pass.
void checkImpactFracture() {
    auto land = [](const rcs::Material& material, double speed, double duration) {
        const bool compliant = material.fracture == rcs::FractureMode::Compliant;
        rcs::PartDefinition cube = rcs::makeSolidCube(compliant ? 0.16 : 0.40);
        cube.material = material;
        if (compliant) {
            cube.sizing.maxElement = 0.05;
            cube.sizing.minElement = 0.02;
            cube.sizing.tetBudget = 300;
        }
        cube.sizing.cached = false;
        cube.volume = {};
        check(rcs::ensurePartMesh(cube), "impact cube meshes");

        rcs::AssemblyDefinition definition;
        definition.name = material.name;
        definition.parts.push_back(cube);

        rcs::Scenario scenario;
        const bool metal = material.fracture == rcs::FractureMode::Ductile;
        scenario.settings.timestep =
            (metal && speed > 0.0) ? (speed <= 20.0 ? 1.0e-4 : 2.0e-4) : 5.0e-4;
        scenario.settings.duration = duration;
        scenario.initial.position = {0.0, compliant ? 0.16 : 0.28, 0.0};
        scenario.initial.orientation = {0.0, 0.0, 0.0};
        scenario.initial.linearVelocity = {0.0, -speed, 0.0};
        scenario.environment.gravity = {0.0, -9.81, 0.0};
        scenario.environment.floorEnabled = true;
        scenario.environment.aerodynamicsEnabled = false;

        rcs::Simulator simulator;
        simulator.reset(definition, scenario);
        const int startLive = simulator.assembly().parts()[0].liveTetCount();
        std::printf("landing %s at %.0f m/s, %d tets, %.3f s\n", material.name.c_str(), speed, startLive,
                    duration);
        std::fflush(stdout);
        bool finite = true;
        double t = 0.0;
        const double dt = scenario.settings.timestep;
        const int steps = static_cast<int>(std::ceil(duration / dt));
        float plastic = 0.0f;
        double contactClosing = 0.0;
        bool seenHit = false;
        const double floor = scenario.environment.floorHeight;
        auto sampleContact = [&]() {
            const rcs::Part& live = simulator.assembly().parts()[0];
            const double offset = live.contactOffset();
            double minGap = 1.0e9;
            double closing = 0.0;
            for (std::size_t n = 0; n < live.nodeCount(); ++n) {
                if (!live.nodeOnLiveTet(n))
                    continue;
                const double gap = live.nodePosition(n).y() - (floor + offset);
                minGap = std::min(minGap, gap);
                if (gap < 0.10)
                    closing = std::max(closing, -live.nodeVelocity(n).y());
            }
            if (!seenHit && minGap < 0.10) {
                contactClosing = std::max(contactClosing, closing);
                if (minGap < 0.02)
                    seenHit = true;
            }
        };
        for (int step = 0; step < steps && finite; ++step) {
            sampleContact();
            finite = simulator.step(dt, t, scenario.environment, scenario.loads);
            t += dt;
            const rcs::Part& live = simulator.assembly().parts()[0];
            std::vector<float> field;
            live.evaluatePlasticStrain(field);
            for (float value : field)
                if (std::isfinite(value))
                    plastic = std::max(plastic, value);
            sampleContact();
        }
        const rcs::Part& done = simulator.assembly().parts()[0];
        struct Result {
            bool finite = false;
            int startLive = 0;
            int live = 0;
            std::size_t shards = 0;
            float plastic = 0.0f;
            bool torn = false;
            double lowest = 0.0;
            double contactClosing = 0.0;
            double worstAspect = 1.0;
        } result;
        result.finite = finite;
        result.startLive = startLive;
        result.live = done.liveTetCount();
        result.shards = done.shards().size();
        result.plastic = plastic;
        result.torn = done.torn();
        result.contactClosing = contactClosing;
        result.lowest = 1.0e9;
        for (std::size_t n = 0; n < done.nodeCount(); ++n)
            result.lowest = std::min(result.lowest, done.nodePosition(n).y());
        for (const auto& shard : done.shards()) {
            const double edgeU = (shard.local[1] - shard.local[0]).Length();
            const double edgeV = (shard.local[3] - shard.local[0]).Length();
            const double shortest = std::min(edgeU, edgeV);
            if (shortest > 1.0e-6)
                result.worstAspect = std::max(result.worstAspect, std::max(edgeU, edgeV) / shortest);
        }
        return result;
    };

    const auto aluminumFast = land(rcs::aluminum6061(), 50.0, 0.025);
    const auto aluminumSlow = land(rcs::aluminum6061(), 15.0, 0.04);
    const auto glassHit = land(rcs::glass(), 50.0, 0.025);
    const auto rubber = land(rcs::rubber(), 50.0, 0.025);
    const auto glassSoft = land(rcs::glass(), 0.0, 0.20);

    std::printf("\n15 m/s aluminum finite=%d live %d -> %d plastic %.4f contact %.1f m/s\n",
                aluminumSlow.finite ? 1 : 0, aluminumSlow.startLive, aluminumSlow.live, aluminumSlow.plastic,
                aluminumSlow.contactClosing);
    std::printf("50 m/s aluminum finite=%d live %d -> %d plastic %.4f torn=%d shards %zu lowest %.4f contact %.1f\n",
                aluminumFast.finite ? 1 : 0, aluminumFast.startLive, aluminumFast.live, aluminumFast.plastic,
                aluminumFast.torn ? 1 : 0, aluminumFast.shards, aluminumFast.lowest, aluminumFast.contactClosing);
    std::printf("50 m/s glass     finite=%d live %d -> %d torn=%d shards %zu aspect %.1f contact %.1f\n",
                glassHit.finite ? 1 : 0, glassHit.startLive, glassHit.live, glassHit.torn ? 1 : 0, glassHit.shards,
                glassHit.worstAspect, glassHit.contactClosing);
    std::printf("50 m/s rubber    finite=%d live %d -> %d torn=%d shards %zu contact %.1f\n",
                rubber.finite ? 1 : 0, rubber.startLive, rubber.live, rubber.torn ? 1 : 0, rubber.shards,
                rubber.contactClosing);
    std::printf("gentle glass     finite=%d live %d -> %d shards %zu lowest %.4f m\n",
                glassSoft.finite ? 1 : 0, glassSoft.startLive, glassSoft.live, glassSoft.shards, glassSoft.lowest);

    const double alSlowFrac =
        aluminumSlow.startLive > 0 ? double(aluminumSlow.live) / double(aluminumSlow.startLive) : 0.0;
    const double alFastFrac =
        aluminumFast.startLive > 0 ? double(aluminumFast.live) / double(aluminumFast.startLive) : 0.0;
    const double glassFrac = glassHit.startLive > 0 ? double(glassHit.live) / double(glassHit.startLive) : 0.0;
    const double rubberFrac = rubber.startLive > 0 ? double(rubber.live) / double(rubber.startLive) : 0.0;
    const double glassSoftFrac =
        glassSoft.startLive > 0 ? double(glassSoft.live) / double(glassSoft.startLive) : 0.0;

    expectFinite(aluminumSlow.finite, "aluminum 15 m/s landing stays in the solver");
    expectRange(aluminumSlow.plastic, 0.002, 0.12, "aluminum 15 m/s landing yields");
    expectRange(alSlowFrac, 0.80, 1.0, "aluminum 15 m/s landing keeps the mesh");
    expectRange(aluminumSlow.contactClosing, 12.0, 18.0, "aluminum 15 m/s still closing at contact");

    expectFinite(aluminumFast.finite, "aluminum 50 m/s landing stays in the solver");
    expectRange(aluminumFast.plastic, 0.002, 0.20, "aluminum 50 m/s landing yields");
    expectRange(alFastFrac, 0.40, 1.05, "aluminum 50 m/s landing keeps a surviving core");
    expectRange(aluminumFast.contactClosing, 40.0, 60.0, "aluminum 50 m/s still closing at contact");

    expectFinite(glassHit.finite, "glass 50 m/s landing stays in the solver");
    expectRange(glassHit.plastic, 0.0, 1.0e-4, "glass 50 m/s does not yield");
    expectRange(glassFrac, 0.20, 0.96, "glass 50 m/s shatters the contact zone, not the peel budget");
    expectRange(glassHit.worstAspect, 1.0, 20.0, "glass 50 m/s shards are not long needles");
    expectRange(glassHit.contactClosing, 40.0, 60.0, "glass 50 m/s still closing at contact");
    check(glassHit.shards > 0 || glassHit.live < glassHit.startLive, "glass 50 m/s detaches overloaded tets");

    expectFinite(rubber.finite, "rubber 50 m/s landing stays in the solver");
    expectRange(rubberFrac, 1.0, 1.0, "rubber 50 m/s landing does not lose tetrahedra");
    expectRange(double(rubber.shards), 0.0, 0.0, "rubber 50 m/s landing makes no shards");
    expectRange(rubber.plastic, 0.0, 1.0e-4, "rubber 50 m/s does not yield");
    expectRange(rubber.contactClosing, 40.0, 60.0, "rubber 50 m/s still closing at contact");

    expectFinite(glassSoft.finite, "a gentle glass drop stays in the solver");
    expectRange(glassSoftFrac, 1.0, 1.0, "a gentle glass drop does not shatter");
    expectRange(double(glassSoft.shards), 0.0, 0.0, "a gentle glass drop makes no shards");
    expectRange(glassSoft.lowest, -0.02, 0.20, "a gentle glass drop does not fall through the floor");
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
    // Level, so the whole face lands at once. The default scenario tips the
    // part half a radian, which is an edge landing and a different test.
    scenario.initial.orientation = {0.0, 0.0, 0.0};
    scenario.environment.aerodynamicsEnabled = false;
    scenario.environment.floorEnabled = true;

    rcs::Simulator simulator;
    simulator.reset(definition, scenario);
    const double dt = scenario.settings.timestep;
    const double touchY = simulator.assembly().parts()[0].contactOffset() + 1.0e-4;
    double t = 0.0;
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
    // A plate slapped onto rigid ground rings. The honest question is whether
    // that ringing decays, so compare the quarter second after first touch
    // with the quarter second at the end of the run.
    double impactAt = -1.0;
    double earlyPeak = 0.0;
    double latePeak = 0.0;
    double finalKe = 0.0;
    double lowest = 1.0e9;
    const double duration = 3.5;
    for (int step = 0; step * dt < duration && finite; ++step) {
        finite = simulator.step(dt, t, scenario.environment, scenario.loads);
        t += dt;
        const double energy = verticalKe();
        double minY = 1.0e9;
        for (const auto& part : simulator.assembly().parts())
            for (std::size_t i = 0; i < part.nodeCount(); ++i)
                minY = std::min(minY, part.nodePosition(i).y());
        lowest = std::min(lowest, minY);
        if (impactAt < 0.0 && minY < touchY)
            impactAt = t;
        if (impactAt > 0.0 && t <= impactAt + 0.25)
            earlyPeak = std::max(earlyPeak, energy);
        if (t > duration - 0.25)
            latePeak = std::max(latePeak, energy);
        finalKe = energy;
        if (std::getenv("RCS_TRACE") && step % 200 == 0)
            std::printf("  t=%.3f minY=%+.6f ke=%.5f\n", t, minY, energy);
    }
    std::printf("panel touched at %.3f s, ring %.4f J -> %.4f J, final %.5f J, lowest %.5f m\n", impactAt, earlyPeak,
                latePeak, finalKe, lowest);
    check(finite, "flat drop stays finite");
    check(lowest > -0.002, "a landed panel does not sink into the floor");
    check(impactAt > 0.0 && latePeak < earlyPeak * 0.2,
          "the landing ring decays instead of being fed by contact");
    check(finalKe < 0.02, "the panel comes to rest");
}

rcs::PartDefinition plate(const char* name, double x, double width) {
    rcs::PartDefinition part;
    part.name = name;
    part.material = rcs::aluminum6061();
    part.thickness = 0.004;
    part.clampBoundary = true;
    part.patch.origin = {x, 0.0, 0.0};
    part.patch.uAxis = {1.0, 0.0, 0.0};
    part.patch.vAxis = {0.0, 0.0, 1.0};
    part.patch.uSize = width;
    part.patch.vSize = 0.12;
    part.patch.divU = 6;
    part.patch.divV = 2;
    return part;
}

rcs::Scenario stillScenario() {
    rcs::Scenario scenario;
    scenario.settings.timestep = 1.0e-3;
    scenario.settings.duration = 0.5;
    scenario.initial.orientation = {0.0, 0.0, 0.0};
    scenario.environment.gravity = {0.0, 0.0, 0.0};
    scenario.environment.floorEnabled = false;
    scenario.environment.aerodynamicsEnabled = false;
    scenario.environment.ambientTemperature = 293.15;
    return scenario;
}

bool stepUntil(rcs::Simulator& simulator, const rcs::Scenario& scenario, double duration) {
    const double dt = scenario.settings.timestep;
    double t = 0.0;
    while (t < duration - 1.0e-9) {
        if (!simulator.step(dt, t, scenario.environment, scenario.loads, scenario.heatSources))
            return false;
        t += dt;
    }
    return simulator.assembly().finite();
}

void checkHeat() {
    rcs::AssemblyDefinition heated;
    heated.name = "Hot plate";
    heated.parts.push_back(plate("Plate", 0.0, 0.40));
    rcs::Scenario scenario = stillScenario();
    scenario.initial.position = {0.20, 0.40, 0.06};

    rcs::Simulator cold;
    cold.reset(heated, scenario);
    const bool coldFinite = stepUntil(cold, scenario, 0.2);
    const rcs::Recording coldFrame = cold.capture();
    bool coldRoom = true;
    for (const auto& part : cold.assembly().parts()) {
        for (std::size_t n = 0; n < part.nodeCount(); ++n) {
            if (std::abs(part.nodeTemperature(n) - 293.15) > 0.05)
                coldRoom = false;
        }
    }
    std::printf("isothermal stress %.2f MPa\n", coldFrame.maxStress * 1.0e-6f);
    check(coldFinite, "isothermal plate stays finite");
    check(coldRoom, "isothermal plate stays at 293.15 K");
    check(coldFrame.maxStress < 1.0e6f, "isothermal plate has no thermal stress");

    rcs::HeatSource source;
    source.temperature = 1100.0;
    source.radius = 0.07;
    source.startTime = 0.0;
    source.duration = 0.5;
    double hotX = 1.0e9;
    for (std::size_t n = 0; n < cold.assembly().parts()[0].nodeCount(); ++n) {
        const auto p = cold.assembly().parts()[0].nodePosition(n);
        if (p.x() < hotX) {
            hotX = p.x();
            source.point = p;
        }
    }
    scenario.heatSources.push_back(source);

    rcs::Simulator hot;
    hot.reset(heated, scenario);
    const bool hotFinite = stepUntil(hot, scenario, 0.5);
    const rcs::Recording hotFrame = hot.capture();
    double nearT = 0.0;
    double farT = 1.0e9;
    double farX = -1.0e9;
    hotX = 1.0e9;
    const auto& platePart = hot.assembly().parts()[0];
    for (std::size_t n = 0; n < platePart.nodeCount(); ++n) {
        const auto p = platePart.nodePosition(n);
        const double kelvin = platePart.nodeTemperature(n);
        if (p.x() < hotX + 1.0e-4) {
            if (p.x() < hotX - 1.0e-4)
                nearT = 0.0;
            hotX = std::min(hotX, p.x());
            nearT = std::max(nearT, kelvin);
        }
        if (p.x() > farX - 1.0e-4) {
            if (p.x() > farX + 1.0e-4)
                farT = 1.0e9;
            farX = std::max(farX, p.x());
            farT = std::min(farT, kelvin);
        }
    }
    float airMax = 293.15f;
    if (!hotFrame.frames.empty()) {
        for (float cell : hotFrame.frames.front().air.temperature)
            airMax = std::max(airMax, cell);
    }
    std::printf("hot edge %.1f K, far edge %.1f K, stress %.1f MPa, air %.1f K\n",
                nearT, farT, hotFrame.maxStress * 1.0e-6f, airMax);
    check(hotFinite, "heated plate stays finite");
    check(nearT > farT + 30.0, "the heated edge is hotter than the far edge");
    check(hotFrame.maxStress > coldFrame.maxStress + 1.0e6f, "the gradient stresses the plate");
    check(!hotFrame.frames.empty() && !hotFrame.frames.front().air.temperature.empty() && airMax > 300.0f,
          "still air inside the source warms");

    rcs::AssemblyDefinition welded;
    welded.name = "Welded plates";
    welded.parts.push_back(plate("Left", 0.0, 0.16));
    welded.parts.push_back(plate("Right", 0.16, 0.16));
    rcs::WeldDefinition weld;
    weld.partA = 0;
    weld.partB = 1;
    weld.tolerance = 3.0e-3;
    welded.welds.push_back(weld);

    rcs::Scenario weldScenario = stillScenario();
    weldScenario.initial.position = {0.16, 0.40, 0.06};
    rcs::Simulator weldCold;
    weldCold.reset(welded, weldScenario);
    const bool weldColdFinite = stepUntil(weldCold, weldScenario, 0.15);
    check(weldColdFinite && weldCold.capture().brokenLinks == 0, "a room-temperature weld holds");

    rcs::HeatSource torch;
    torch.temperature = 1600.0;
    torch.radius = 0.06;
    torch.duration = 0.45;
    const auto& left = weldCold.assembly().parts()[0];
    chrono::ChVector3d seam{0, 0, 0};
    int seamCount = 0;
    double seamX = -1.0e9;
    for (std::size_t n = 0; n < left.nodeCount(); ++n)
        seamX = std::max(seamX, left.nodePosition(n).x());
    for (std::size_t n = 0; n < left.nodeCount(); ++n) {
        const auto p = left.nodePosition(n);
        if (std::abs(p.x() - seamX) < 1.0e-4) {
            seam += p;
            ++seamCount;
        }
    }
    if (seamCount > 0)
        torch.point = seam / seamCount;
    weldScenario.heatSources.push_back(torch);

    rcs::Simulator melted;
    melted.reset(welded, weldScenario);
    const bool meltFinite = stepUntil(melted, weldScenario, 0.45);
    const int broken = melted.capture().brokenLinks;
    std::printf("melted weld links %d\n", broken);
    check(meltFinite, "a melted weld stays finite");
    check(broken > 0, "a weld through the aluminum solidus breaks");

    // A body that has already settled, then heated with no applied load, must
    // stay put. Expansion can change its size. It cannot start walking.
    rcs::PartDefinition cube = rcs::makeSolidCube(0.20);
    cube.sizing.maxElement = 0.08;
    cube.sizing.minElement = 0.04;
    cube.sizing.refineLimit = 0.1;
    check(rcs::ensurePartMesh(cube), "resting heat cube meshes");
    rcs::AssemblyDefinition resting;
    resting.parts.push_back(cube);
    rcs::Scenario restScene;
    restScene.settings.timestep = 1.0e-3;
    restScene.initial.position = {0.4, 0.12, -0.3};
    restScene.initial.orientation = {0.0, 0.0, 0.0};
    restScene.environment.gravity = {0.0, -9.81, 0.0};
    restScene.environment.floorEnabled = true;
    restScene.environment.aerodynamicsEnabled = false;
    rcs::Simulator restSim;
    restSim.reset(resting, restScene);
    const double dt = restScene.settings.timestep;
    double t = 0.0;
    bool restFinite = true;
    for (int step = 0; step < 400 && restFinite; ++step) {
        restFinite = restSim.step(dt, t, restScene.environment, {});
        t += dt;
    }
    auto comOf = [](const rcs::Part& part) {
        chrono::ChVector3d com(0, 0, 0);
        double mass = 0.0;
        for (std::size_t n = 0; n < part.nodeCount(); ++n) {
            if (part.nodeRigid(n))
                continue;
            const double m = part.nodeMass()[n];
            com += part.nodePosition(n) * m;
            mass += m;
        }
        return mass > 0.0 ? com / mass : com;
    };
    const chrono::ChVector3d settled = comOf(restSim.assembly().parts()[0]);
    double settledSpeed = 0.0;
    {
        chrono::ChVector3d momentum(0, 0, 0);
        double mass = 0.0;
        const auto& part = restSim.assembly().parts()[0];
        for (std::size_t n = 0; n < part.nodeCount(); ++n) {
            if (part.nodeRigid(n))
                continue;
            momentum += part.nodeVelocity(n) * part.nodeMass()[n];
            mass += part.nodeMass()[n];
        }
        settledSpeed = mass > 0.0 ? (momentum / mass).Length() : 0.0;
    }
    for (int step = 0; step < 300 && restFinite; ++step) {
        rcs::Part& part = restSim.assembly().parts()[0];
        for (std::size_t n = 0; n < part.nodeCount(); ++n) {
            if (part.nodeRigid(n))
                continue;
            part.setNodeTemperature(n, 500.0);
            part.setStrainTemperature(n, 500.0);
            part.setHeatActive(n, true);
        }
        part.syncThermalMechanics();
        restFinite = restSim.step(dt, t, restScene.environment, {});
        t += dt;
    }
    const chrono::ChVector3d heatedCom = comOf(restSim.assembly().parts()[0]);
    const double walked = std::hypot(heatedCom.x() - settled.x(), heatedCom.z() - settled.z());
    double heatedSpeed = 0.0;
    {
        chrono::ChVector3d momentum(0, 0, 0);
        double mass = 0.0;
        const auto& part = restSim.assembly().parts()[0];
        for (std::size_t n = 0; n < part.nodeCount(); ++n) {
            if (part.nodeRigid(n))
                continue;
            const chrono::ChVector3d v = part.nodeVelocity(n);
            momentum += chrono::ChVector3d(v.x(), 0.0, v.z()) * part.nodeMass()[n];
            mass += part.nodeMass()[n];
        }
        heatedSpeed = mass > 0.0 ? (momentum / mass).Length() : 0.0;
    }
    std::printf("resting heat finite=%d settled speed %.4f m/s, heated horizontal %.4f m/s, walked %.4f m\n",
                restFinite ? 1 : 0, settledSpeed, heatedSpeed, walked);
    check(restFinite, "a resting solid stays finite while it heats");
    check(settledSpeed < 0.05, "the solid is at rest before the heat");
    check(heatedSpeed < 0.02 && walked < 0.01, "heating a resting solid does not start it moving");
}

}  // namespace

/// A clamped aluminum plate at 3000 K must bow far more under the same ramp
/// than the same plate at room temperature. This is the malleability check.
void checkMalleable() {
    rcs::AssemblyDefinition def;
    auto sheet = plate("Plate", 0.0, 0.40);
    sheet.patch.vSize = 0.40;
    sheet.patch.divU = 4;
    sheet.patch.divV = 4;
    sheet.thickness = 0.008;
    def.name = "Hot sheet";
    def.parts.push_back(sheet);

    rcs::Scenario scenario = stillScenario();
    scenario.settings.timestep = 5.0e-4;
    scenario.initial.position = {0.20, 0.50, 0.20};

    rcs::Simulator probe;
    probe.reset(def, scenario);
    const auto& coldPart = probe.assembly().parts()[0];
    chrono::ChVector3d center = coldPart.nodePosition(0);
    for (std::size_t n = 1; n < coldPart.nodeCount(); ++n) {
        const auto p = coldPart.nodePosition(n);
        if ((p - chrono::ChVector3d(0.20, 0.50, 0.20)).Length2() <
            (center - chrono::ChVector3d(0.20, 0.50, 0.20)).Length2())
            center = p;
    }

    auto deflect = [&](double kelvin) {
        rcs::ImpactLoad hit;
        hit.name = "Press";
        hit.partIndex = 0;
        hit.point = center;
        hit.direction = {0.0, -1.0, 0.0};
        hit.initialForce = 0.0;
        hit.finalForce = 10000.0;
        hit.rampTime = 0.08;
        hit.duration = 0.12;
        hit.radius = 0.05;
        scenario.loads.clear();
        scenario.loads.push_back(hit);

        rcs::Simulator simulator;
        simulator.reset(def, scenario);
        auto& part = simulator.assembly().parts()[0];
        std::vector<chrono::ChVector3d> start;
        start.reserve(part.nodeCount());
        for (std::size_t n = 0; n < part.nodeCount(); ++n)
            start.push_back(part.nodePosition(n));
        if (kelvin > 290.0) {
            // Colour is already hot. Strain is still at the build temperature,
            // which is the state right after a sphere brings a face to 3000 K.
            for (std::size_t n = 0; n < part.nodeCount(); ++n)
                part.setNodeTemperature(n, kelvin);
            part.syncThermalMechanics();
        }

        const double dt = scenario.settings.timestep;
        double t = 0.0;
        bool finite = true;
        double peak = 0.0;
        double end = 0.0;
        for (int step = 0; step < 440 && finite; ++step) {
            finite = simulator.step(dt, t, scenario.environment, scenario.loads);
            t += dt;
            if (!finite)
                break;
            end = 0.0;
            for (std::size_t n = 0; n < part.nodeCount(); ++n) {
                const double moved = (part.nodePosition(n) - start[n]).Length();
                peak = std::max(peak, moved);
                end = std::max(end, moved);
            }
        }
        std::printf("  strain %.1f K, node %.1f K\n", part.strainTemperature(part.nodeCount() / 2),
                    part.nodeTemperature(part.nodeCount() / 2));
        std::printf("malleable T=%.0f K finite=%d torn=%d t=%.4f peak %.4f m, after release %.4f m\n", kelvin,
                    finite ? 1 : 0, part.torn() ? 1 : 0, t, peak, end);
        return std::tuple<double, double, bool, bool>(peak, end, finite, part.torn());
    };

    const auto room = deflect(293.15);
    const auto hot = deflect(3000.0);
    check(std::get<2>(hot), "3000 K aluminum stays in the solver under a 10 kN ramp");
    check(!std::get<3>(room), "room-temperature aluminum does not tear under a 10 kN ramp");
    check(std::get<3>(hot) || (std::get<0>(hot) > std::get<0>(room) * 5.0 && std::get<1>(hot) > 0.02),
          "3000 K aluminum tears or keeps a deep fold under a 10 kN ramp");
}

/// A weld has a rated shear strength, and a bead of a given throat over a
/// given seam length should therefore carry a predictable force and no more.
/// Two bars are butted end to end, the first is held, and the second is pulled
/// with a force swept from well under the rating to well over it. The bead has
/// to survive the low pull and let go somewhere near the rating.
void checkWeldStrength() {
    const rcs::Material metal = rcs::aluminum6061();
    const double width = 0.06;
    const double thick = 0.006;
    const double length = 0.12;
    // Throat defaults to the part thickness, so the seam carries
    // weldShearStrength over width*thickness.
    const double rated = metal.weldShearStrength * width * thick;

    auto pull = [&](double force) {
        rcs::PartDefinition left;
        left.name = "Held";
        left.material = metal;
        left.thickness = thick;
        left.fixed = true;
        left.patch.origin = {0.0, 0.0, -0.5 * width};
        left.patch.uAxis = {1.0, 0.0, 0.0};
        left.patch.vAxis = {0.0, 0.0, 1.0};
        left.patch.uSize = length;
        left.patch.vSize = width;
        left.patch.divU = 6;
        left.patch.divV = 3;

        rcs::PartDefinition right = left;
        right.name = "Pulled";
        right.fixed = false;
        right.patch.origin = {length, 0.0, -0.5 * width};

        rcs::WeldDefinition seam;
        seam.label = "Butt";
        seam.partA = 0;
        seam.partB = 1;
        seam.tolerance = 3.0e-3;

        rcs::AssemblyDefinition definition;
        definition.name = "Butt weld";
        definition.parts.push_back(left);
        definition.parts.push_back(right);
        definition.welds.push_back(seam);

        rcs::Scenario scenario;
        scenario.settings.timestep = 2.0e-4;
        scenario.initial.position = {0.0, 0.6, 0.0};
        scenario.initial.orientation = {0.0, 0.0, 0.0};
        scenario.environment.gravity = {0.0, 0.0, 0.0};
        scenario.environment.floorEnabled = false;
        scenario.environment.aerodynamicsEnabled = false;

        // Pull the free bar straight off the seam, along the bar axis.
        rcs::ImpactLoad tug;
        tug.name = "Tug";
        tug.partIndex = 1;
        tug.point = {2.0 * length, 0.6, 0.0};
        tug.direction = {1.0, 0.0, 0.0};
        tug.initialForce = 0.0;
        tug.finalForce = force;
        tug.rampTime = 0.08;
        tug.duration = 0.12;
        tug.radius = length;
        scenario.loads.push_back(tug);

        rcs::Simulator simulator;
        simulator.reset(definition, scenario);
        const int links = simulator.assembly().brokenLinkCount();
        bool finite = true;
        double t = 0.0;
        for (int step = 0; step < 800 && finite; ++step) {
            finite = simulator.step(2.0e-4, t, scenario.environment, scenario.loads);
            t += 2.0e-4;
        }
        const int broken = simulator.assembly().brokenLinkCount() - links;
        return std::pair<bool, int>(finite, broken);
    };

    const auto light = pull(0.3 * rated);
    const auto heavy = pull(2.5 * rated);
    std::printf("\nbutt weld rated %.0f N: light pull broke %d links, heavy pull broke %d\n", rated,
                light.second, heavy.second);
    check(light.first && heavy.first, "a weld under load stays in the solver");
    check(light.second == 0, "a weld holds well under its rated shear");
    check(heavy.second > 0, "a weld lets go well over its rated shear");
}

/// Every material, driven to its own failure mode by a load derived from its
/// own strength rather than a number picked to suit aluminium. A cantilever of
/// width w and thickness t loaded at its tip carries an outer-fibre stress of
/// 6*F*L/(w*t^2), so asking for a given fraction of a material's strength is a
/// closed-form force. Each bar is loaded once well inside its capacity, where
/// nothing may fail and it has to spring back, and once well past it, where
/// the mode of failure has to match what the material is.
void checkMaterialFailureModes() {
    struct Case {
        const char* name;
        rcs::Material material;
    };
    const std::vector<Case> cases = {
        {"aluminum 6061", rcs::aluminum6061()}, {"mild steel", rcs::mildSteel()},
        {"AR500", rcs::ar500()},                {"glass", rcs::glass()},
        {"tempered glass", rcs::temperedGlass()}, {"rubber", rcs::rubber()},
    };

    const double length = 0.20;
    const double width = 0.05;
    const double thick = 0.006;

    struct Result {
        bool finite = false;
        bool torn = false;
        double peakStress = 0.0;
        double peakPlastic = 0.0;
        double peakTip = 0.0;
        double residualTip = 0.0;
        std::size_t shards = 0;
    };

    for (const Case& item : cases) {
        rcs::PartDefinition bar;
        bar.name = item.name;
        bar.material = item.material;
        bar.thickness = thick;
        bar.clampBoundary = true;
        bar.patch.origin = {0.0, 0.0, -0.5 * width};
        bar.patch.uAxis = {1.0, 0.0, 0.0};
        bar.patch.vAxis = {0.0, 0.0, 1.0};
        bar.patch.uSize = length;
        bar.patch.vSize = width;
        bar.patch.divU = 10;
        bar.patch.divV = 3;

        rcs::AssemblyDefinition definition;
        definition.name = item.name;
        definition.parts.push_back(bar);

        // Brittle materials never reach a yield point, so their reference is
        // the stress that cracks them; everything else is measured against
        // first yield.
        const bool brittle = item.material.fracture == rcs::FractureMode::Brittle;
        const bool compliant = item.material.fracture == rcs::FractureMode::Compliant;
        const double reference = brittle ? item.material.fractureStrength
                                         : (compliant ? item.material.youngsModulus : item.material.yieldStress);
        auto forceFor = [&](double stress) {
            return stress * width * thick * thick / (6.0 * length);
        };
        const double rated = brittle ? item.material.fractureStrength : item.material.yieldStress;
        const double area = width * thick;
        const double serviceForce =
            compliant ? item.material.youngsModulus * 0.15 * area : forceFor(0.5 * rated);
        const double overloadForce =
            compliant ? item.material.youngsModulus * 0.90 * area : forceFor(3.0 * rated);
        const chrono::ChVector3d loadDir = compliant ? chrono::ChVector3d(1.0, 0.0, 0.0)
                                                     : chrono::ChVector3d(0.0, -1.0, 0.0);

        auto press = [&](double force) {
            rcs::Scenario scenario;
            scenario.settings.timestep = 2.0e-4;
            // Authored bar lives on x in [0, L]. The drop pose offsets the
            // origin, not the live COM, so keep X at zero.
            scenario.initial.position = {0.0, 0.6, 0.0};
            scenario.initial.orientation = {0.0, 0.0, 0.0};
            scenario.environment.gravity = {0.0, 0.0, 0.0};
            scenario.environment.floorEnabled = false;
            scenario.environment.aerodynamicsEnabled = false;

            rcs::ImpactLoad hit;
            hit.name = "Tip";
            hit.partIndex = 0;
            hit.point = {0.95 * length, 0.6, 0.0};
            hit.direction = loadDir;
            hit.initialForce = 0.0;
            hit.finalForce = force;
            hit.rampTime = 0.04;
            hit.duration = 0.06;
            hit.radius = 0.5 * width;
            scenario.loads.push_back(hit);

            rcs::Simulator simulator;
            simulator.reset(definition, scenario);
            rcs::Part& part = simulator.assembly().parts()[0];

            // Axial stretch is read at the free end; bending at mid-span.
            std::size_t tip = 0;
            double best = 1.0e9;
            const double tipX = compliant ? length : 0.6 * length;
            for (std::size_t n = 0; n < part.nodeCount(); ++n) {
                const auto rest = part.restPosition(n);
                const double score = std::abs(rest.x() - tipX) + std::abs(rest.z()) * 4.0;
                if (score < best) {
                    best = score;
                    tip = n;
                }
            }
            const double start = compliant ? part.nodePosition(tip).x() : part.nodePosition(tip).y();

            Result result;
            result.finite = true;
            double t = 0.0;
            const double dt = scenario.settings.timestep;
            // Load, then twice as long again unloaded so anything elastic has
            // time to come back and only a real permanent set is left.
            for (int step = 0; step < 900 && result.finite; ++step) {
                result.finite = simulator.step(dt, t, scenario.environment, scenario.loads);
                t += dt;
                if (!result.finite)
                    break;
                const rcs::Part& live = simulator.assembly().parts()[0];
                const double now = compliant ? live.nodePosition(tip).x() : live.nodePosition(tip).y();
                result.peakTip = std::max(result.peakTip, std::abs(now - start));
                std::vector<float> stresses;
                live.evaluateStress(stresses);
                for (float value : stresses)
                    result.peakStress = std::max(result.peakStress, double(value));
                std::vector<float> plastic;
                live.evaluatePlasticStrain(plastic);
                for (float value : plastic)
                    result.peakPlastic = std::max(result.peakPlastic, double(value));
            }
            const rcs::Part& done = simulator.assembly().parts()[0];
            const double end = compliant ? done.nodePosition(tip).x() : done.nodePosition(tip).y();
            result.residualTip = std::abs(end - start);
            result.torn = done.torn();
            result.shards = done.shards().size();
            return result;
        };

        const Result service = press(serviceForce);
        const Result overload = press(overloadForce);

        std::printf("\n%s (%s, ref %.0f MPa)\n", item.name,
                    brittle ? "brittle" : (compliant ? "compliant" : "ductile"), reference * 1.0e-6);
        std::printf("  service : finite=%d peak %.1f MPa plastic %.5f tip %.4f m residual %.5f m torn=%d shards %zu\n",
                    service.finite ? 1 : 0, service.peakStress * 1.0e-6, service.peakPlastic,
                    service.peakTip, service.residualTip, service.torn ? 1 : 0, service.shards);
        std::printf("  overload: finite=%d peak %.1f MPa plastic %.5f tip %.4f m residual %.5f m torn=%d shards %zu\n",
                    overload.finite ? 1 : 0, overload.peakStress * 1.0e-6, overload.peakPlastic,
                    overload.peakTip, overload.residualTip, overload.torn ? 1 : 0, overload.shards);

        const std::string tag = item.name;
        expectFinite(service.finite && overload.finite, (tag + " stays in the solver at both loads").c_str());
        check(!service.torn && service.shards == 0, (tag + " survives half its rated stress").c_str());
        if (!compliant)
            expectRange(service.peakStress / std::max(reference, 1.0), 0.25, 1.05,
                        tag + " service peak stress is near half its rating");
        expectRange(service.residualTip / std::max(service.peakTip, 1.0e-9), 0.0, 0.05,
                    tag + " springs back from a service load");

        if (brittle) {
            expectRange(overload.peakPlastic, 0.0, 1.0e-4, tag + " does not yield before it cracks");
            check(overload.torn || overload.shards > 0, (tag + " cracks under three times its strength").c_str());
        } else if (compliant) {
            expectRange(overload.peakPlastic, 0.0, 1.0e-4, tag + " stretches without yielding");
            expectRange(overload.peakTip / std::max(service.peakTip, 1.0e-9), 3.0, 80.0,
                        tag + " goes much further as the load grows");
        } else {
            expectRange(overload.peakPlastic, 0.002, 0.50, tag + " yields under three times its yield stress");
            expectRange(overload.residualTip, 1.0e-4, 0.25, tag + " keeps a permanent set");
            expectRange(double(overload.shards), 0.0, 0.0, tag + " bends instead of shattering");
        }
    }
}

/// A bar hotter on one face must curl, because the hot face is longer. This
/// is the bimetallic mechanism and it is the clearest test that thermal strain
/// is carried per element rather than as one number for the whole part.
///
/// For a linear gradient the curvature is alpha*dT/h, so a cantilever of
/// length L lifts its tip by alpha*dT*L^2/(2h). The band is 0.5-2x that
/// closed form, enough for shear and a coarse mesh but not a sign error.
void checkThermalWarp() {
    const double length = 0.30;
    const double thickness = 0.010;
    const double hot = 500.0;
    const double cold = 293.15;

    rcs::PartDefinition bar;
    bar.name = "Bar";
    bar.material = rcs::aluminum6061();
    bar.thickness = thickness;
    bar.clampBoundary = false;
    bar.patch.origin = {0.0, 0.0, -0.03};
    bar.patch.uAxis = {1.0, 0.0, 0.0};
    bar.patch.vAxis = {0.0, 0.0, 1.0};
    bar.patch.uSize = length;
    bar.patch.vSize = 0.06;
    bar.patch.divU = 12;
    bar.patch.divV = 2;

    rcs::AssemblyDefinition definition;
    definition.name = "Heated bar";
    definition.parts.push_back(bar);

    rcs::Scenario scenario;
    scenario.settings.timestep = 5.0e-4;
    scenario.initial.position = {0.15, 0.5, 0.0};
    scenario.initial.orientation = {0.0, 0.0, 0.0};
    scenario.environment.gravity = {0.0, 0.0, 0.0};
    scenario.environment.floorEnabled = false;
    scenario.environment.aerodynamicsEnabled = false;

    rcs::Simulator simulator;
    simulator.reset(definition, scenario);

    rcs::Part& part = simulator.assembly().parts()[0];
    auto pinRoot = [](rcs::Part& p) {
        double root = 1.0e9;
        for (std::size_t n = 0; n < p.nodeCount(); ++n)
            root = std::min(root, p.restPosition(n).x());
        auto pin = [&](auto& nodes) {
            for (std::size_t n = 0; n < nodes.size() && n < p.nodeCount(); ++n) {
                if (nodes[n] && p.restPosition(n).x() <= root + 1.0e-3)
                    nodes[n]->SetFixed(true);
            }
        };
        pin(p.nodes());
        pin(p.xyzNodes());
    };
    pinRoot(part);
    simulator.system().Setup();
    simulator.system().ForceUpdate();

    double minY = 1.0e9;
    double maxY = -1.0e9;
    double maxX = -1.0e9;
    for (std::size_t n = 0; n < part.nodeCount(); ++n) {
        minY = std::min(minY, part.restPosition(n).y());
        maxY = std::max(maxY, part.restPosition(n).y());
        maxX = std::max(maxX, part.restPosition(n).x());
    }
    const double span = std::max(maxY - minY, 1.0e-9);
    // The tip node on the midplane, so the measurement is bending and not
    // the through-thickness growth of the bar itself.
    std::size_t tip = 0;
    double best = 1.0e9;
    for (std::size_t n = 0; n < part.nodeCount(); ++n) {
        const auto rest = part.restPosition(n);
        const double score = (maxX - rest.x()) + std::abs(rest.y() - 0.5 * (minY + maxY)) * 10.0;
        if (score < best) {
            best = score;
            tip = n;
        }
    }
    const double tipStart = part.nodePosition(tip).y();
    const double alpha = rcs::aluminum6061().thermalExpansion;
    const double kappa = -alpha * (hot - cold) / thickness;

    // A linear through-thickness gradient is thermal curvature on a Reissner
    // midsurface. Driving that (and the node temperature) keeps this a test
    // of the thermo-mechanical coupling and not of the air model.
    auto soak = [&]() {
        rcs::Part& p = simulator.assembly().parts()[0];
        for (std::size_t n = 0; n < p.nodeCount(); ++n) {
            const double u = (p.restPosition(n).y() - minY) / span;
            const double kelvin = cold + (hot - cold) * u;
            p.setNodeTemperature(n, kelvin);
            p.setStrainTemperature(n, kelvin);
            p.setHeatActive(n, true);
        }
        p.syncThermalMechanics();
        p.setThermalCurvature(kappa);
    };

    bool finite = true;
    double t = 0.0;
    for (int step = 0; step < 1200 && finite; ++step) {
        soak();
        finite = simulator.step(5.0e-4, t, scenario.environment, scenario.loads);
        t += 5.0e-4;
    }

    const rcs::Part& warped = simulator.assembly().parts()[0];
    const double tipEnd = warped.nodePosition(tip).y();
    const double deflection = tipEnd - tipStart;
    const double expected = alpha * (hot - cold) * length * length / (2.0 * thickness);
    std::printf("thermal warp tip %+.4f m (beam theory %+.4f m), finite=%d\n", deflection, expected, finite ? 1 : 0);

    expectFinite(finite, "a bar with a thermal gradient stays in the solver");
    check(deflection < 0.0, "a bar hotter on top curls away from the heat");
    expectRange(std::abs(deflection) / std::max(expected, 1.0e-12), 0.50, 2.0,
                "the curl is within 0.5-2x of alpha*dT*L^2/2h");

    // The same bar heated uniformly must grow but not bend.
    rcs::Simulator even;
    even.reset(definition, scenario);
    pinRoot(even.assembly().parts()[0]);
    even.system().Setup();
    even.system().ForceUpdate();
    const double evenStart = even.assembly().parts()[0].nodePosition(tip).y();
    bool evenFinite = true;
    t = 0.0;
    for (int step = 0; step < 1200 && evenFinite; ++step) {
        rcs::Part& p = even.assembly().parts()[0];
        for (std::size_t n = 0; n < p.nodeCount(); ++n) {
            p.setNodeTemperature(n, hot);
            p.setStrainTemperature(n, hot);
            p.setHeatActive(n, true);
        }
        p.syncThermalMechanics();
        p.setThermalCurvature(0.0);
        evenFinite = even.step(5.0e-4, t, scenario.environment, scenario.loads);
        t += 5.0e-4;
    }
    const double evenDeflection = even.assembly().parts()[0].nodePosition(tip).y() - evenStart;
    std::printf("uniform soak tip %+.5f m\n", evenDeflection);
    expectFinite(evenFinite, "a uniformly heated bar stays in the solver");
    expectRange(std::abs(evenDeflection), 0.0, 0.10 * expected,
                "a uniformly heated bar grows without bending");
}

void checkSolidRigidMotion() {
    rcs::PartDefinition cube = rcs::makeSolidCube();
    cube.sizing.maxElement = 0.1;
    cube.sizing.minElement = 0.04;
    cube.sizing.refineLimit = 0.1;
    check(rcs::ensurePartMesh(cube), "coarse solid cube meshes");

    rcs::AssemblyDefinition solid;
    solid.parts.push_back(cube);

    rcs::Scenario spinScene;
    spinScene.settings.timestep = 1.0e-3;
    spinScene.initial.position = {0.0, 0.6, 0.0};
    spinScene.initial.orientation = {0.0, 0.0, 0.0};
    spinScene.initial.angularVelocity = {0.0, 2.0, 0.0};
    spinScene.environment.gravity = {0.0, 0.0, 0.0};
    spinScene.environment.floorEnabled = false;
    spinScene.environment.aerodynamicsEnabled = false;

    rcs::Simulator spin;
    spin.reset(solid, spinScene);
    const rcs::Part& spun = spin.assembly().parts()[0];
    chrono::ChVector3d center(0, 0, 0);
    for (std::size_t n = 0; n < spun.nodeCount(); ++n)
        center += spun.nodePosition(n);
    center /= static_cast<double>(std::max<std::size_t>(1, spun.nodeCount()));
    std::size_t marker = 0;
    double markerX = 0.0;
    for (std::size_t n = 0; n < spun.nodeCount(); ++n) {
        const double dx = spun.nodePosition(n).x() - center.x();
        if (std::abs(dx) > std::abs(markerX)) {
            markerX = dx;
            marker = n;
        }
    }
    const chrono::ChVector3d markerStart = spun.nodePosition(marker);
    const double radius = (markerStart - center).Length();
    bool finite = true;
    for (int step = 0; step < 200 && finite; ++step)
        finite = spin.step(1.0e-3, step * 1.0e-3, spinScene.environment, {});
    const rcs::Part& after = spin.assembly().parts()[0];
    chrono::ChVector3d centerAfter(0, 0, 0);
    for (std::size_t n = 0; n < after.nodeCount(); ++n)
        centerAfter += after.nodePosition(n);
    centerAfter /= static_cast<double>(std::max<std::size_t>(1, after.nodeCount()));
    const chrono::ChVector3d markerEnd = after.nodePosition(marker);
    const double turned = (markerEnd - centerAfter).z() - (markerStart - center).z();
    double worstStretch = 0.0;
    for (int tet = 0; tet < after.tetCount() && tet < 40; ++tet) {
        if (!after.tetLive(static_cast<std::size_t>(tet)))
            continue;
        const auto corners = after.tets()[static_cast<std::size_t>(tet)];
        for (int a = 0; a < 4; ++a) {
            const int b = (a + 1) % 4;
            const double rest = (after.restPosition(static_cast<std::size_t>(corners[a])) -
                                 after.restPosition(static_cast<std::size_t>(corners[b])))
                                    .Length();
            const double now = (after.nodePosition(static_cast<std::size_t>(corners[a])) -
                                after.nodePosition(static_cast<std::size_t>(corners[b])))
                                   .Length();
            if (rest > 1.0e-6)
                worstStretch = std::max(worstStretch, std::abs(now - rest) / rest);
        }
    }
    std::printf("solid spin finite=%d turn %.4f m (radius %.3f) stretch %.4f\n", finite ? 1 : 0, turned, radius,
                worstStretch);
    check(finite, "a spinning solid cube stays in the solver");
    check(std::abs(turned) > 0.03, "a solid cube given a spin actually rotates");
    check(worstStretch < 0.05, "a spinning solid cube keeps its shape");

    rcs::Scenario dropScene = spinScene;
    dropScene.initial.angularVelocity = {0.0, 0.0, 0.0};
    dropScene.initial.orientation = {0.6, 0.0, 0.0};
    dropScene.initial.position = {0.0, 0.55, 0.0};
    dropScene.environment.gravity = {0.0, -9.81, 0.0};
    dropScene.environment.floorEnabled = true;
    rcs::Simulator drop;
    drop.reset(solid, dropScene);
    const rcs::Part& resting = drop.assembly().parts()[0];
    chrono::ChVector3d dropCenter(0, 0, 0);
    for (std::size_t n = 0; n < resting.nodeCount(); ++n)
        dropCenter += resting.nodePosition(n);
    dropCenter /= static_cast<double>(resting.nodeCount());
    std::size_t top = 0;
    for (std::size_t n = 0; n < resting.nodeCount(); ++n)
        if (resting.nodePosition(n).y() > resting.nodePosition(top).y())
            top = n;
    const chrono::ChVector3d arm = resting.nodePosition(top) - dropCenter;
    bool dropFinite = true;
    double lowest = 1.0e9;
    for (int step = 0; step < 400 && dropFinite; ++step) {
        dropFinite = drop.step(1.0e-3, step * 1.0e-3, dropScene.environment, {});
        const rcs::Part& moving = drop.assembly().parts()[0];
        for (std::size_t n = 0; n < moving.nodeCount(); ++n)
            lowest = std::min(lowest, moving.nodePosition(n).y());
    }
    const rcs::Part& landed = drop.assembly().parts()[0];
    chrono::ChVector3d landedCenter(0, 0, 0);
    for (std::size_t n = 0; n < landed.nodeCount(); ++n)
        landedCenter += landed.nodePosition(n);
    landedCenter /= static_cast<double>(landed.nodeCount());
    const chrono::ChVector3d armAfter = landed.nodePosition(top) - landedCenter;
    const double alignment = arm.Length() > 1.0e-4 && armAfter.Length() > 1.0e-4
                                 ? arm.Dot(armAfter) / (arm.Length() * armAfter.Length())
                                 : 1.0;
    std::printf("solid tumble finite=%d alignment %.3f lowest %.4f m\n", dropFinite ? 1 : 0, alignment, lowest);
    check(dropFinite, "a tilted solid cube lands without diverging");
    check(lowest > -0.02, "a solid cube does not fall through the floor");
    check(alignment < 0.98, "a tilted solid cube rotates when it hits the floor");
}

void checkSolidAndMixed() {
    checkSolidRigidMotion();
    rcs::PartDefinition cube = rcs::makeSolidCube();
    check(cube.volume.empty(), "the solid cube waits to be meshed");
    check(rcs::ensurePartMesh(cube), "the solid cube meshes on demand");
    const int generation = cube.sizing.generation;
    const std::size_t tets = cube.volume.tets.size();
    rcs::ensurePartMesh(cube);
    check(cube.sizing.generation == generation, "a second open does not remesh the cube");
    check(cube.volume.tets.size() == tets && tets > 100, "solid cube is a tet fill");
    check(cube.surface.empty(), "solid cube has no shell skin");

    rcs::AssemblyDefinition solid;
    solid.name = "Solid cube";
    solid.parts.push_back(cube);
    solid.parts.back().clampBoundary = true;

    rcs::Scenario scenario;
    scenario.settings.timestep = 1.0e-3;
    scenario.settings.duration = 0.08;
    scenario.initial.position = {0.0, 0.30, 0.0};
    scenario.initial.orientation = {0.0, 0.0, 0.0};
    scenario.environment.gravity = {0.0, 0.0, 0.0};
    scenario.environment.aerodynamicsEnabled = false;
    scenario.environment.floorEnabled = false;
    rcs::ImpactLoad load;
    load.name = "Dent";
    load.partIndex = 0;
    load.point = {0.20, 0.30, 0.0};
    load.direction = {-1.0, 0.0, 0.0};
    load.initialForce = 0.0;
    load.finalForce = 2.5e6;
    load.rampTime = 0.04;
    load.duration = 0.06;
    load.radius = 0.02;
    scenario.loads = {load};

    rcs::Simulator dent;
    dent.reset(solid, scenario);
    const double mass = dent.assembly().mass();
    std::printf("solid cube mass %.1f kg, %d tets\n", mass, dent.assembly().parts()[0].tetCount());
    check(mass > 140.0 && mass < 200.0, "solid cube has the full aluminum mass");
    bool finite = true;
    double t = 0.0;
    const double dt = scenario.settings.timestep;
    double peakInward = 0.0;
    for (int step = 0; step < 80 && finite; ++step) {
        finite = dent.step(dt, t, scenario.environment, scenario.loads);
        t += dt;
        const rcs::Part& moving = dent.assembly().parts()[0];
        for (std::size_t n = 0; n < moving.nodeCount(); ++n) {
            const auto rest = moving.restPosition(n);
            if (rest.x() < 0.17 || (rest - load.point).Length() > load.radius)
                continue;
            peakInward = std::max(peakInward, rest.x() - moving.nodePosition(n).x());
        }
    }
    const rcs::Part& dented = dent.assembly().parts()[0];
    float plastic = 0.0f;
    std::vector<float> plasticField;
    dented.evaluatePlasticStrain(plasticField);
    for (float value : plasticField)
        plastic = std::max(plastic, value);
    double inward = 0.0;
    int faceCount = 0;
    for (std::size_t n = 0; n < dented.nodeCount(); ++n) {
        const auto rest = dented.restPosition(n);
        if (rest.x() < 0.17)
            continue;
        if ((rest - load.point).Length() > load.radius)
            continue;
        inward = std::max(inward, rest.x() - dented.nodePosition(n).x());
        ++faceCount;
    }
    std::printf("solid dent finite=%d plastic=%.4f peak %.5f m inward %.5f m (%d nodes)\n", finite ? 1 : 0, plastic,
                peakInward, inward, faceCount);
    check(finite, "solid cube dent stays finite");
    check(plastic > 1.0e-4, "loaded face of the solid cube yields");
    check(faceCount > 0 && inward > 1.0e-5 && inward < 0.05, "solid cube keeps a permanent set");

    rcs::Simulator soakSim;
    soakSim.reset(solid, scenario);
    rcs::HeatSource source;
    source.point = {0.22, 0.30, 0.0};
    source.radius = 0.05;
    source.temperature = 700.0;
    source.duration = 1.0;
    std::vector<rcs::HeatSource> heat = {source};
    bool soakFinite = true;
    t = 0.0;
    for (int step = 0; step < 100 && soakFinite; ++step) {
        soakFinite = soakSim.step(dt, t, scenario.environment, {}, heat);
        t += dt;
    }
    const rcs::Part& soaked = soakSim.assembly().parts()[0];
    double surface = 0.0;
    double core = 0.0;
    int surfaceCount = 0;
    int coreCount = 0;
    for (std::size_t n = 0; n < soaked.nodeCount(); ++n) {
        const auto p = soaked.nodePosition(n);
        if (p.x() > 0.15 && (p - source.point).Length() < source.radius + 0.03) {
            surface = std::max(surface, soaked.nodeTemperature(n));
            ++surfaceCount;
        } else if (std::abs(p.x()) < 0.04 && std::abs(p.y() - 0.30) < 0.04 && std::abs(p.z()) < 0.04) {
            core += soaked.nodeTemperature(n);
            ++coreCount;
        }
    }
    if (coreCount)
        core /= coreCount;
    std::printf("heat soak surface %.1f K (%d) core %.1f K (%d)\n", surface, surfaceCount, core, coreCount);
    check(soakFinite && surfaceCount > 0 && coreCount > 0, "heat soak has a surface and a core");
    check(surface > core + 30.0, "a heated face warms before the core");

    rcs::Simulator chipSim;
    chipSim.reset(solid, scenario);
    rcs::Part& hot = chipSim.assembly().parts()[0];
    int hotTet = -1;
    double hotDistance = 1.0e9;
    chrono::ChVector3d hotPoint = load.point;
    for (std::size_t tet = 0; tet < hot.tets().size(); ++tet) {
        int onFace = 0;
        chrono::ChVector3d centroid{0, 0, 0};
        int count = 0;
        for (int node : hot.tets()[tet]) {
            if (node < 0)
                continue;
            const auto p = hot.restPosition(static_cast<std::size_t>(node));
            centroid += p;
            ++count;
            if (p.x() > 0.17)
                ++onFace;
        }
        if (onFace < 3 || count != 4)
            continue;
        centroid *= 0.25;
        const double distance = (centroid - load.point).Length();
        if (distance < hotDistance) {
            hotDistance = distance;
            hotTet = static_cast<int>(tet);
            hotPoint = centroid;
        }
    }
    check(hotTet >= 0, "solid cube has a surface tet");
    if (hotTet >= 0) {
        for (int node : hot.tets()[static_cast<std::size_t>(hotTet)]) {
            hot.setNodeTemperature(static_cast<std::size_t>(node), 2000.0);
            hot.setStrainTemperature(static_cast<std::size_t>(node), 2000.0);
            hot.setHeatActive(static_cast<std::size_t>(node), true);
        }
        hot.syncThermalMechanics();
    }
    const int before = hot.liveTetCount();
    rcs::ImpactLoad poke = load;
    poke.point = hotPoint;
    poke.finalForce = 3.0e3;
    poke.rampTime = 0.02;
    poke.duration = 0.05;
    poke.radius = 0.03;
    bool chipFinite = true;
    t = 0.0;
    for (int step = 0; step < 40 && chipFinite; ++step) {
        chipFinite = chipSim.step(dt, t, scenario.environment, {poke});
        t += dt;
    }
    const rcs::Part& chipped = chipSim.assembly().parts()[0];
    int coreLive = 0;
    for (std::size_t tet = 0; tet < chipped.tets().size(); ++tet) {
        if (!chipped.tetLive(tet))
            continue;
        bool interior = true;
        for (int node : chipped.tets()[tet]) {
            if (node < 0)
                continue;
            const auto p = chipped.restPosition(static_cast<std::size_t>(node));
            if (std::abs(p.x()) > 0.05 || std::abs(p.y() - 0.30) > 0.05 || std::abs(p.z()) > 0.05)
                interior = false;
        }
        coreLive += interior ? 1 : 0;
    }
    std::printf("molten surface live %d -> %d, core live %d, finite=%d\n", before, chipped.liveTetCount(),
                coreLive, chipFinite ? 1 : 0);
    check(chipFinite, "a molten surface tet leaves without aborting");
    check(!chipped.shards().empty(), "a molten surface tet leaves the solve");
    check(coreLive > 0, "the core tet stays after the surface chip leaves");

    // A molten tet has to fall onto the floor and then leave the shard loop.
    rcs::PartDefinition looseCube = rcs::makeSolidCube(0.16);
    looseCube.sizing.maxElement = 0.08;
    looseCube.sizing.minElement = 0.04;
    looseCube.sizing.refineLimit = 0.1;
    check(rcs::ensurePartMesh(looseCube), "molten chip cube meshes");
    rcs::AssemblyDefinition loose;
    loose.parts.push_back(looseCube);
    rcs::Scenario fallScene;
    fallScene.settings.timestep = 1.0e-3;
    fallScene.initial.position = {0.0, 0.20, 0.0};
    fallScene.initial.orientation = {0.0, 0.0, 0.0};
    fallScene.environment.gravity = {0.0, -9.81, 0.0};
    fallScene.environment.floorEnabled = true;
    fallScene.environment.aerodynamicsEnabled = false;
    rcs::Simulator fallSim;
    fallSim.reset(loose, fallScene);
    rcs::Part& donor = fallSim.assembly().parts()[0];
    int donorTet = -1;
    double donorY = -1.0e9;
    for (std::size_t tet = 0; tet < donor.tets().size(); ++tet) {
        if (!donor.tetLive(tet))
            continue;
        double y = 0.0;
        int count = 0;
        bool exposed = false;
        for (int node : donor.tets()[tet]) {
            if (node < 0)
                continue;
            const auto p = donor.restPosition(static_cast<std::size_t>(node));
            y += p.y();
            ++count;
            if (p.y() > 0.25)
                exposed = true;
        }
        if (!exposed || count != 4)
            continue;
        y /= 4.0;
        if (y > donorY) {
            donorY = y;
            donorTet = static_cast<int>(tet);
        }
    }
    check(donorTet >= 0, "the cube has a tet that can fall to the floor");
    if (donorTet >= 0) {
        for (int node : donor.tets()[static_cast<std::size_t>(donorTet)]) {
            donor.setNodeTemperature(static_cast<std::size_t>(node), 2000.0);
            donor.setStrainTemperature(static_cast<std::size_t>(node), 2000.0);
            donor.setHeatActive(static_cast<std::size_t>(node), true);
        }
        donor.syncThermalMechanics();
    }
    bool fallFinite = true;
    double ft = 0.0;
    int asleep = 0;
    int peakShards = 0;
    for (int step = 0; step < 800 && fallFinite; ++step) {
        fallFinite = fallSim.step(1.0e-3, ft, fallScene.environment, {});
        ft += 1.0e-3;
        asleep = 0;
        const int live = static_cast<int>(fallSim.assembly().parts()[0].shards().size());
        peakShards = std::max(peakShards, live);
        for (const auto& shard : fallSim.assembly().parts()[0].shards())
            asleep += shard.asleep ? 1 : 0;
    }
    const auto& fallen = fallSim.assembly().parts()[0];
    double chipY = 1.0e9;
    for (const auto& shard : fallen.shards())
        chipY = std::min(chipY, shard.com.y());
    std::printf("molten chip finite=%d shards %zu peak %d asleep %d lowest %.4f m at t=%.3f\n", fallFinite ? 1 : 0,
                fallen.shards().size(), peakShards, asleep, fallen.shards().empty() ? -1.0 : chipY, ft);
    for (const auto& shard : fallen.shards()) {
        std::printf("  leftover y=%.4f speed=%.4f spin=%.3f\n", shard.com.y(), shard.velocity.Length(),
                    shard.angularVelocity.Length());
        // #region agent log
        {
            std::ofstream dbg("c:/Users/nirav/Documents/Projects/Robot Sim/robot_combat_sim/debug-b58308.log",
                              std::ios::app);
            if (dbg)
                dbg << "{\"sessionId\":\"b58308\",\"hypothesisId\":\"A\",\"location\":\"physics_smoke.cpp:molten\","
                       "\"message\":\"leftover shard\",\"data\":{\"y\":"
                    << shard.com.y() << ",\"speed\":" << shard.velocity.Length()
                    << ",\"spin\":" << shard.angularVelocity.Length() << ",\"asleep\":" << (shard.asleep ? 1 : 0)
                    << "},\"timestamp\":0}\n";
        }
        // #endregion
    }
    check(fallFinite, "a molten tet falls without aborting");
    check(peakShards > 0, "a molten tet becomes a chip");
    check(fallen.shards().empty(), "a molten chip on the floor leaves the active shard loop");

    rcs::PartDefinition fin;
    fin.name = "Fin block";
    fin.material = rcs::aluminum6061();
    fin.thickness = 0.003;
    fin.sizing.shellCutoff = 0.008;
    fin.sizing.minElement = 0.002;
    fin.sizing.maxElement = 0.02;
    fin.clampBoundary = true;
    const chrono::ChVector3d block[8] = {
        {0, 0, -0.02}, {0.08, 0, -0.02}, {0.08, 0.04, -0.02}, {0, 0.04, -0.02},
        {0, 0, 0.02},  {0.08, 0, 0.02},  {0.08, 0.04, 0.02},  {0, 0.04, 0.02}};
    fin.source.vertices.assign(block, block + 8);
    fin.source.vertices.push_back({0.02, 0.04, 0.0});
    fin.source.vertices.push_back({0.06, 0.04, 0.0});
    const int r0 = 8;
    const int r1 = 9;
    const int boxFaces[10][3] = {
        {0, 2, 1}, {0, 3, 2}, {4, 5, 6}, {4, 6, 7}, {0, 1, 5}, {0, 5, 4},
        {0, 4, 7}, {0, 7, 3}, {1, 2, 6}, {1, 6, 5}};
    for (const auto& face : boxFaces)
        fin.source.triangles.push_back({face[0], face[1], face[2]});
    const int top[6][3] = {{3, r1, 2}, {3, r0, r1}, {3, 7, r0}, {2, r1, 6}, {r1, 7, 6}, {r1, r0, 7}};
    for (const auto& face : top)
        fin.source.triangles.push_back({face[0], face[1], face[2]});
    fin.source.vertices.push_back({0.02, 0.07, 0.0015});
    fin.source.vertices.push_back({0.06, 0.07, 0.0015});
    fin.source.vertices.push_back({0.02, 0.07, -0.0015});
    fin.source.vertices.push_back({0.06, 0.07, -0.0015});
    const int t0 = 10, t1 = 11, u0 = 12, u1 = 13;
    fin.source.triangles.push_back({r0, r1, t1});
    fin.source.triangles.push_back({r0, t1, t0});
    fin.source.triangles.push_back({r0, u0, u1});
    fin.source.triangles.push_back({r0, u1, r1});
    check(rcs::ensurePartMesh(fin), "fin on a block meshes");
    check(!fin.volume.empty(), "the block is tetrahedra");
    check(fin.surface.empty(), "the fin is tetrahedra too");
    std::printf("mixed fin: %zu tets\n", fin.volume.tets.size());

    rcs::AssemblyDefinition mixed;
    mixed.parts.push_back(fin);
    chrono::ChVector3d centroid{0, 0, 0};
    for (const auto& vertex : fin.source.vertices)
        centroid += vertex;
    centroid /= static_cast<double>(fin.source.vertices.size());
    rcs::Scenario finRun;
    finRun.settings.timestep = 5.0e-4;
    finRun.initial.position = centroid + chrono::ChVector3d{0.0, 0.05, 0.0};
    finRun.initial.orientation = {0.0, 0.0, 0.0};
    finRun.environment.gravity = {0.0, 0.0, 0.0};
    finRun.environment.aerodynamicsEnabled = false;
    finRun.environment.floorEnabled = false;
    rcs::ImpactLoad tip;
    tip.partIndex = 0;
    tip.point = chrono::ChVector3d{0.04, 0.07, 0.0} + (finRun.initial.position - centroid);
    tip.direction = {0.0, 0.0, 1.0};
    tip.finalForce = 400.0;
    tip.duration = 0.02;
    tip.radius = 0.015;
    finRun.loads = {tip};
    rcs::Simulator finSim;
    finSim.reset(mixed, finRun);
    const rcs::Part& built = finSim.assembly().parts()[0];
    double minY = 1.0e9;
    for (std::size_t n = 0; n < built.nodeCount(); ++n)
        minY = std::min(minY, built.restPosition(n).y());
    bool bridge = false;
    double highest = -1.0e9;
    for (std::size_t tet = 0; tet < built.tets().size(); ++tet) {
        if (!built.tetLive(tet))
            continue;
        double ymin = 1.0e9;
        double ymax = -1.0e9;
        for (int node : built.tets()[tet]) {
            if (node < 0)
                continue;
            const double y = built.restPosition(static_cast<std::size_t>(node)).y();
            ymin = std::min(ymin, y);
            ymax = std::max(ymax, y);
        }
        highest = std::max(highest, ymax);
        if (ymax > minY + 0.06 && ymin <= minY + 0.045)
            bridge = true;
    }
    std::printf("fin highest %.3f m above base, bridge %d\n", highest - minY, bridge ? 1 : 0);
    check(bridge, "fin and block share the root nodes");
    bool finFinite = true;
    t = 0.0;
    for (int step = 0; step < 40 && finFinite; ++step) {
        finFinite = finSim.step(finRun.settings.timestep, t, finRun.environment, finRun.loads);
        t += finRun.settings.timestep;
    }
    std::vector<float> finStress;
    finSim.assembly().parts()[0].evaluateStress(finStress);
    float peak = 0.0f;
    for (float value : finStress)
        peak = std::max(peak, value);
    std::printf("fin load finite=%d peak %.2f MPa\n", finFinite ? 1 : 0, peak * 1.0e-6);
    check(finFinite, "a load through the fin into the block stays finite");
    check(peak > 1.0e4, "the fin load reaches the block");
}

void checkUnifiedTets() {
    rcs::PartDefinition cube = rcs::makeSolidCube();
    check(rcs::ensurePartMesh(cube), "solid cube meshes");
    check(!cube.volume.empty() && cube.surface.empty(), "solid cube is only tetrahedra");
    check(cube.volume.tets.size() < 8000, "resting solid stays under the tet budget");
    std::printf("solid cube tets %zu\n", cube.volume.tets.size());

    rcs::PartDefinition sheet;
    sheet.name = "Sheet";
    sheet.material = rcs::aluminum6061();
    sheet.thickness = 0.008;
    sheet.patch.origin = {-0.1, -0.1, 0.0};
    sheet.patch.uAxis = {1, 0, 0};
    sheet.patch.vAxis = {0, 1, 0};
    sheet.patch.uSize = 0.2;
    sheet.patch.vSize = 0.2;
    sheet.patch.divU = 2;
    sheet.patch.divV = 2;
    sheet.sizing.throughThickness = 2;
    check(rcs::ensurePartMesh(sheet), "thin panel meshes");
    check(!sheet.volume.empty() && sheet.surface.empty(), "thin panel is a tet volume");

    rcs::PartDefinition glass = cube;
    glass.material = rcs::glass();
    glass.sizing.cached = false;
    glass.volume = {};
    rcs::PartDefinition rubber = cube;
    rubber.material = rcs::rubber();
    rubber.sizing.cached = false;
    rubber.volume = {};
    check(rcs::ensurePartMesh(glass) && rcs::ensurePartMesh(rubber), "brittle and compliant solids mesh the same way");

    rcs::AssemblyDefinition box = rcs::makeTestCube();
    bool hollow = box.parts.size() == 6;
    for (auto& part : box.parts) {
        if (!rcs::ensurePartMesh(part) || part.volume.empty() || !part.surface.empty())
            hollow = false;
        std::printf("hollow wall %s tets %zu\n", part.name.c_str(), part.volume.tets.size());
        check(part.volume.tets.size() < 3000, "a hollow-cube wall stays well under 12k tets");
        const double inner = 0.5 * (0.40 - 0.008) - 0.004;
        for (const auto& tet : part.volume.tets) {
            chrono::ChVector3d center(0, 0, 0);
            for (int node : tet)
                center += part.volume.vertices[static_cast<std::size_t>(node)];
            center *= 0.25;
            if (std::abs(center.x()) < inner && std::abs(center.y()) < inner && std::abs(center.z()) < inner)
                hollow = false;
        }
    }
    check(hollow, "hollow six-panel box is thin tets with an empty cavity");

    rcs::Scenario rest;
    rest.settings.duration = 0.08;
    rest.settings.timestep = 5.0e-4;
    rest.settings.captureRate = 20.0;
    rest.initial.position = {0.0, 0.28, 0.0};
    rest.environment.floorEnabled = true;
    rest.environment.aerodynamicsEnabled = false;
    rcs::AssemblyDefinition solid;
    solid.parts.push_back(cube);
    rcs::Simulator sitting;
    const rcs::Recording sat = sitting.run(solid, rest, nullptr);
    check(!sat.diverged, "a resting solid cube does not diverge");

    rcs::Scenario boxRest = rest;
    boxRest.initial.position = {0.0, 0.30, 0.0};
    rcs::Simulator boxRun;
    const rcs::Recording boxed = boxRun.run(rcs::makeTestCube(), boxRest, nullptr);
    std::printf("hollow box rest %s broken %d\n", boxed.diverged ? "diverged" : "finite", boxed.brokenLinks);
    check(!boxed.diverged, "hollow welded box rests without diverging");
    check(boxed.brokenLinks == 0, "self-weight does not open a hollow-box seam");
}

/// Volume enclosed by a closed triangle soup, by the divergence theorem. Used
/// to check that the tet fill reproduces the volume of the surface it was
/// given, which is the one number an imported CAD part has to get right.
double surfaceVolume(const rcs::SurfaceMesh& mesh) {
    double sum = 0.0;
    for (const auto& triangle : mesh.triangles) {
        const auto& a = mesh.vertices[static_cast<std::size_t>(triangle[0])];
        const auto& b = mesh.vertices[static_cast<std::size_t>(triangle[1])];
        const auto& c = mesh.vertices[static_cast<std::size_t>(triangle[2])];
        sum += a.Dot(b.Cross(c));
    }
    return std::abs(sum) / 6.0;
}

/// Worst element shape in a tet fill, as volume divided by the cube of the
/// mean edge length. A regular tetrahedron scores about 0.117 and a sliver
/// tends to zero. Slivers are what collapse the integrator step, so this is
/// the number that decides whether a mesh is usable, not the tet count.
double worstTetQuality(const rcs::VolumeMesh& mesh) {
    double worst = 1.0;
    for (const auto& tet : mesh.tets) {
        const auto& a = mesh.vertices[static_cast<std::size_t>(tet[0])];
        const auto& b = mesh.vertices[static_cast<std::size_t>(tet[1])];
        const auto& c = mesh.vertices[static_cast<std::size_t>(tet[2])];
        const auto& d = mesh.vertices[static_cast<std::size_t>(tet[3])];
        const double volume = std::abs((b - a).Dot((c - a).Cross(d - a))) / 6.0;
        const double edges = ((b - a).Length() + (c - a).Length() + (d - a).Length() +
                              (c - b).Length() + (d - b).Length() + (d - c).Length()) /
                             6.0;
        if (edges <= 0.0)
            return 0.0;
        worst = std::min(worst, volume / (edges * edges * edges));
    }
    return worst;
}

/// A ring of points on a circle in the XZ plane at height y.
void addRing(rcs::SurfaceMesh& mesh, double radius, double y, int segments) {
    for (int i = 0; i < segments; ++i) {
        const double angle = 2.0 * 3.14159265358979323846 * double(i) / double(segments);
        mesh.vertices.push_back({radius * std::cos(angle), y, radius * std::sin(angle)});
    }
}

/// Stitch two equal-sized rings into a band of triangles wound outward.
void bridgeRings(rcs::SurfaceMesh& mesh, int lower, int upper, int segments, bool flip) {
    for (int i = 0; i < segments; ++i) {
        const int j = (i + 1) % segments;
        const int a = lower + i;
        const int b = lower + j;
        const int c = upper + j;
        const int d = upper + i;
        if (flip) {
            mesh.triangles.push_back({a, c, b});
            mesh.triangles.push_back({a, d, c});
        } else {
            mesh.triangles.push_back({a, b, c});
            mesh.triangles.push_back({a, c, d});
        }
    }
}

/// Close a ring onto a single point, as a cap or a cone apex.
void capRing(rcs::SurfaceMesh& mesh, int ring, int segments, int tip, bool flip) {
    for (int i = 0; i < segments; ++i) {
        const int j = (i + 1) % segments;
        if (flip)
            mesh.triangles.push_back({ring + j, ring + i, tip});
        else
            mesh.triangles.push_back({ring + i, ring + j, tip});
    }
}

rcs::SurfaceMesh sphereSurface(double radius, int segments, int rings) {
    rcs::SurfaceMesh mesh;
    for (int r = 1; r < rings; ++r) {
        const double phi = 3.14159265358979323846 * double(r) / double(rings);
        addRing(mesh, radius * std::sin(phi), radius * std::cos(phi), segments);
    }
    for (int r = 0; r + 2 < rings; ++r)
        bridgeRings(mesh, r * segments, (r + 1) * segments, segments, true);
    const int north = static_cast<int>(mesh.vertices.size());
    mesh.vertices.push_back({0.0, radius, 0.0});
    const int south = static_cast<int>(mesh.vertices.size());
    mesh.vertices.push_back({0.0, -radius, 0.0});
    capRing(mesh, 0, segments, north, false);
    capRing(mesh, (rings - 2) * segments, segments, south, true);
    return mesh;
}

rcs::SurfaceMesh cylinderSurface(double radius, double height, int segments) {
    rcs::SurfaceMesh mesh;
    addRing(mesh, radius, 0.0, segments);
    addRing(mesh, radius, height, segments);
    bridgeRings(mesh, 0, segments, segments, true);
    const int bottom = static_cast<int>(mesh.vertices.size());
    mesh.vertices.push_back({0.0, 0.0, 0.0});
    const int top = static_cast<int>(mesh.vertices.size());
    mesh.vertices.push_back({0.0, height, 0.0});
    capRing(mesh, 0, segments, bottom, false);
    capRing(mesh, segments, segments, top, true);
    return mesh;
}

rcs::SurfaceMesh coneSurface(double radius, double height, int segments) {
    rcs::SurfaceMesh mesh;
    addRing(mesh, radius, 0.0, segments);
    const int base = static_cast<int>(mesh.vertices.size());
    mesh.vertices.push_back({0.0, 0.0, 0.0});
    const int apex = static_cast<int>(mesh.vertices.size());
    mesh.vertices.push_back({0.0, height, 0.0});
    capRing(mesh, 0, segments, base, false);
    capRing(mesh, 0, segments, apex, true);
    return mesh;
}

/// An open-ended tube: a wall thin enough to exercise the thin-region path
/// while still enclosing a cavity the fill must leave empty.
rcs::SurfaceMesh tubeSurface(double outer, double inner, double height, int segments) {
    rcs::SurfaceMesh mesh;
    addRing(mesh, outer, 0.0, segments);     // 0: outer rim, bottom
    addRing(mesh, outer, height, segments);  // 1: outer rim, top
    addRing(mesh, inner, 0.0, segments);     // 2: bore, bottom
    addRing(mesh, inner, height, segments);  // 3: bore, top
    // The bore faces inward, so its winding is the opposite of the outer wall.
    bridgeRings(mesh, 0 * segments, 1 * segments, segments, true);
    bridgeRings(mesh, 2 * segments, 3 * segments, segments, false);
    bridgeRings(mesh, 2 * segments, 0 * segments, segments, true);
    bridgeRings(mesh, 1 * segments, 3 * segments, segments, true);
    return mesh;
}

/// Arbitrary CAD has to work with nothing but a material assigned to it. Each
/// shape is meshed, weighed against the volume of the surface it came from,
/// then dropped on the floor and required to settle without sinking. A curved
/// or pointed contact patch is a much harder case than a flat panel, because
/// the load arrives through one or two nodes.
void checkImportedShapes() {
    struct Shape {
        const char* name;
        rcs::SurfaceMesh surface;
        double analytic;
        double dropHeight;
        /// A body of revolution standing on a curved side should roll away
        /// under its own spin; a cone on its base or a tube on its rim should
        /// not travel far.
        bool rolls;
    };

    const double pi = 3.14159265358979323846;
    std::vector<Shape> shapes;
    shapes.push_back(
        {"sphere", sphereSurface(0.06, 24, 12), 4.0 / 3.0 * pi * 0.06 * 0.06 * 0.06, 0.10, true});
    shapes.push_back({"cylinder", cylinderSurface(0.05, 0.12, 24), pi * 0.05 * 0.05 * 0.12, 0.10, false});
    shapes.push_back({"cone", coneSurface(0.06, 0.14, 24), pi * 0.06 * 0.06 * 0.14 / 3.0, 0.10, false});
    shapes.push_back({"tube", tubeSurface(0.05, 0.042, 0.12, 24),
                      pi * (0.05 * 0.05 - 0.042 * 0.042) * 0.12, 0.10, false});

    const char* only = std::getenv("RCS_SHAPE");
    for (Shape& shape : shapes) {
        if (only && std::strcmp(only, shape.name) != 0)
            continue;
        rcs::PartDefinition part;
        part.name = shape.name;
        part.material = rcs::aluminum6061();
        part.source = shape.surface;
        part.sizing.shellCutoff = 0.008;
        part.sizing.minElement = 0.002;
        part.sizing.maxElement = 0.025;
        part.sizing.throughThickness = 3;

        const double faceted = surfaceVolume(shape.surface);
        const bool meshed = rcs::ensurePartMesh(part);
        const double quality = worstTetQuality(part.volume);
        std::printf("\n%s: %zu tets, worst element %.4f, faceted volume %.3e m^3 (analytic %.3e)\n",
                    shape.name, part.volume.tets.size(), quality, faceted, shape.analytic);
        check(meshed && !part.volume.empty(), std::string(std::string(shape.name) + " meshes into tets").c_str());
        check(part.surface.empty(), std::string(std::string(shape.name) + " needs no shell skin").c_str());
        // A 24-sided approximation loses about a percent of a true circle.
        check(faceted > 0.94 * shape.analytic && faceted < 1.01 * shape.analytic,
              std::string(std::string(shape.name) + " surface encloses the expected volume").c_str());
        // A regular tetrahedron scores 0.118. Delaunay fills of a curved or
        // pointed boundary bottom out three to four orders below that, and
        // pushing further costs elements without buying accuracy; what is not
        // acceptable is an inverted or numerically flat element, which is
        // what this rules out.
        check(quality > 1.0e-5,
              std::string(std::string(shape.name) + " fill has no degenerate elements").c_str());
        // The point of the sizing field is that an arbitrary import stays
        // affordable. A part that quietly meshes to a hundred thousand
        // elements will not run.
        check(part.volume.tets.size() < 45000,
              std::string(std::string(shape.name) + " stays inside the element budget").c_str());
        if (!meshed || part.volume.empty() || std::getenv("RCS_MESH_ONLY"))
            continue;

        rcs::AssemblyDefinition definition;
        definition.name = shape.name;
        definition.parts.push_back(part);

        rcs::Scenario scenario;
        scenario.settings.timestep = 2.0e-4;
        scenario.initial.position = {0.0, shape.dropHeight, 0.0};
        scenario.initial.orientation = {0.0, 0.0, 0.0};
        // A spin about Z rolls the shape along X once it is down, which is the
        // check that contact friction drives rotation instead of locking it.
        scenario.initial.angularVelocity = {0.0, 0.0, -6.0};
        scenario.environment.aerodynamicsEnabled = false;

        rcs::Simulator simulator;
        simulator.reset(definition, scenario);
        const double mass = simulator.assembly().mass();
        const double expected = rcs::aluminum6061().density * faceted;
        std::printf("  mass %.4f kg (surface volume gives %.4f kg)\n", mass, expected);
        check(std::abs(mass - expected) < 0.06 * expected,
              std::string(std::string(shape.name) + " fills its volume to within 6% by mass").c_str());

        double lowest = 1.0e9;
        double startX = 0.0;
        double endX = 0.0;
        bool finite = true;
        double t = 0.0;
        const double dt = scenario.settings.timestep;
        for (int step = 0; step < 3000 && finite; ++step) {
            finite = simulator.step(dt, t, scenario.environment, scenario.loads);
            t += dt;
            const rcs::Part& live = simulator.assembly().parts()[0];
            chrono::ChVector3d com(0, 0, 0);
            double total = 0.0;
            for (std::size_t i = 0; i < live.nodeCount(); ++i) {
                lowest = std::min(lowest, live.nodePosition(i).y());
                com += live.nodeMass()[i] * live.nodePosition(i);
                total += live.nodeMass()[i];
            }
            if (total > 0.0)
                com /= total;
            if (step == 0)
                startX = com.x();
            endX = com.x();
            if (std::getenv("RCS_TRACE") && step % 250 == 0)
                std::printf("    t=%.3f com y %.4f x %.4f lowest %.5f\n", t, com.y(), com.x(),
                            lowest);
        }

        const double travel = endX - startX;
        std::printf("  finite=%d lowest %+.5f m, com travelled %+.4f m in x\n", finite ? 1 : 0,
                    lowest, travel);
        check(finite, std::string(std::string(shape.name) + " lands without diverging").c_str());
        check(lowest > -0.003, std::string(std::string(shape.name) + " does not sink into the floor").c_str());
        if (shape.rolls) {
            // Spun backwards about Z, contact friction has to convert the spin
            // into travel along +x. A floor that locks rotation leaves it
            // where it fell.
            check(travel > 0.02, std::string(std::string(shape.name) + " rolls on its spin instead of locking").c_str());
        }
    }
}

/// The original end-to-end drop: a welded hollow box falls half a metre onto a
/// corner, takes a hammer blow, and has to settle. This is the one scenario
/// that exercises recording, welds, contact and the load path together.
void checkBoxDrop(double hitForce, double stepSize, double hitRadius) {
    const rcs::AssemblyDefinition definition = rcs::makeTestCube();

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
    hit.initialForce = hitForce;
    hit.finalForce = hitForce;
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

    // Judge the landing on the whole timeline rather than the last frame: the
    // box has to reach the floor, stay out of it, and give up its drop energy
    // instead of being thrown back into the air by the contact. The rebound is
    // only measured before the hammer lands, since a 6 kN blow is meant to
    // throw the box around.
    // Rebound is measured on the centre of mass, not the lowest node: a box
    // tumbling off a corner lifts its lowest node simply by rotating, which
    // costs no energy and is not a bounce. How far the centre of mass climbs
    // back compared with how far it fell is the restitution of the landing.
    double touchdown = -1.0;
    double comAtTouchdown = 0.0;
    double comRise = 0.0;
    double deepest = 1.0e9;
    const double comStart = meanHeight(first);
    for (const rcs::Frame& frame : recording.frames) {
        const Bounds bounds = verticalBounds(frame);
        deepest = std::min(deepest, bounds.lowest);
        if (touchdown < 0.0) {
            if (bounds.lowest < 0.01) {
                touchdown = frame.time;
                comAtTouchdown = meanHeight(frame);
            }
        } else if (frame.time < hit.startTime) {
            comRise = std::max(comRise, meanHeight(frame) - comAtTouchdown);
        }
    }
    const double comFall = std::max(comStart - comAtTouchdown, 1.0e-9);

    std::printf("start y [%.3f, %.3f] -> end y [%.3f, %.3f]\n", start.lowest, start.highest,
                end.lowest, end.highest);
    std::printf("touchdown %.3f s, deepest %.4f m, com fell %.3f m and recovered %.3f m\n",
                touchdown, deepest, comFall, comRise);

    check(start.lowest > 0.1, "box starts clear of the floor");
    check(touchdown > 0.0, "the box reaches the floor");
    check(deepest > -0.02, "box never sinks through the floor");
    // A quarter of the drop back is a restitution of one half, which is
    // already generous for a hollow box on concrete.
    check(comRise < 0.25 * comFall, "the box does not bounce back off the floor");
    check(end.lowest < 0.05, "the box is still on the floor when the run ends");
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
    scenario.loads.front().initialForce = 4.0e6;
    scenario.loads.front().finalForce = 4.0e6;
    rcs::Simulator overRun;
    const rcs::Recording overRecording = overRun.run(definition, scenario, nullptr);
    std::printf("\nover-driven load: %s\n", overRecording.message.c_str());
    check(overRecording.diverged, "an unresolvable load is flagged as diverged");
}

#include "physics_core/physics_catalog.inc"

int main(int argc, char** argv) {
    // `--only=<group>` runs one group. The suite is minutes long, so a single
    // failure is diagnosed against the group that produced it.
    if (argc > 1 && std::strncmp(argv[1], "--only=", 7) == 0) {
        const std::string only = argv[1] + 7;
        if (only == "tets")
            checkUnifiedTets();
        else if (only == "materials")
            checkMaterialsAndImport();
        else if (only == "dent")
            checkPermanentDent();
        else if (only == "glass") {
            checkGlassShatters();
            checkImpactFracture();
        }
        else if (only == "contact")
            checkContactSettles();
        else if (only == "heat")
            checkHeat();
        else if (only == "malleable")
            checkMalleable();
        else if (only == "solid")
            checkSolidAndMixed();
        else if (only == "rigid")
            checkSolidRigidMotion();
        else if (only == "warp")
            checkThermalWarp();
        else if (only == "drop")
            checkBoxDrop(6000.0, 5.0e-4, 0.12);
        else if (only == "shapes")
            checkImportedShapes();
        else if (only == "failure")
            checkMaterialFailureModes();
        else if (only == "welds")
            checkWeldStrength();
        else if (only == "midair")
            checkNoMidairKick();
        else if (only == "conserve")
            checkConserve();
        else if (only == "elastic")
            checkElasticityCatalog();
        else if (only == "plastic")
            checkPlasticityCatalog();
        else if (only == "brittle")
            checkBrittleCatalog();
        else if (only == "stretch")
            checkCompliantCatalog();
        else if (only == "hits")
            checkHitsCatalog();
        else if (only == "thermal")
            checkThermalCatalog();
        else if (only == "stress")
            checkStressCatalog();
        else if (only == "fracture")
            checkFractureCatalog();
        else if (only == "extreme")
            checkExtremeCatalog();
        else if (only == "refine")
            checkRefineCatalog();
        else if (only == "friction")
            checkFrictionCatalog();
        else if (only == "catalog")
            checkCatalog();
        else {
            std::printf("unknown group '%s'\n", only.c_str());
            std::printf("groups: tets materials dent glass contact heat malleable solid rigid warp drop "
                        "shapes failure welds midair conserve elastic plastic brittle stretch hits thermal "
                        "stress fracture extreme refine friction catalog\n");
            return 2;
        }
        std::printf(failures == 0 ? "\nall checks passed\n" : "\n%d check(s) failed\n", failures);
        return failures == 0 ? 0 : 1;
    }

    checkUnifiedTets();

    const double hitForce = argc > 1 ? std::atof(argv[1]) : 6000.0;
    const double stepSize = argc > 2 ? std::atof(argv[2]) : 5.0e-4;
    const double hitRadius = argc > 3 ? std::atof(argv[3]) : 0.12;
    checkBoxDrop(hitForce, stepSize, hitRadius);

    std::printf("\n");
    checkMaterialsAndImport();
    checkPermanentDent();
    checkGlassShatters();
    checkImpactFracture();
    checkContactSettles();
    checkHeat();
    checkThermalWarp();
    checkMalleable();
    checkSolidRigidMotion();
    checkSolidAndMixed();
    checkImportedShapes();
    checkMaterialFailureModes();
    checkWeldStrength();
    std::printf("\n");
    checkNoMidairKick();
    checkCatalog();

    std::printf(failures == 0 ? "\nall checks passed\n" : "\n%d check(s) failed\n", failures);
    return failures == 0 ? 0 : 1;
}

