// Issue #6: coupled thermal-stress validation.
//
// The chain under test: node temperature -> (a) thermal strain alpha * dT and
// (b) stiffness / strength retention from the handbook tables in
// material.hpp -> element stress -> joints release at the solidus.
//
// What the code does today:
//   Part::syncThermalMechanics hands each tet setThermalState(alpha * dT,
//   stiffnessRetention, strengthRetention); expansion uses the lagged strain
//   temperature, retention uses the node temperature.
//   Part::evaluateStress then multiplies tet von Mises by
//   thermalStiffnessRetention(part-mean temperature) AGAIN, although the tet
//   computed that stress with the reduced modulus already (part.cpp ~975 and
//   ~1008). Hot solid parts therefore display r^2 instead of r.
//   material.hpp picks a retention table by thresholds: Brittle -> glass,
//   Compliant -> rubber, meltTemperature >= 1400 K -> Eurocode steel, else
//   6061; and strength lookups run 1.45x faster when yieldStress >= 1 GPa
//   (AR500). These are preset choices disguised as numbers, and each one is
//   a cliff for any material near the threshold.
//   HeatField::advance returns early when no source is on and nothing is
//   marked dirty, so a part that starts hot may never cool.
//   Welds release when the cooler end is at or above the lower solidus and
//   either end is heat-active (assembly.cpp ~1786).

#include "physics_core/material.hpp"
#include "physics_core/part.hpp"
#include "physics_core/physics_units.hpp"
#include "physics_core/scenario.hpp"
#include "physics_core/simulator.hpp"
#include "physics_core/thermal_tet.hpp"

#include "chrono/core/ChVector3.h"
#include "chrono/fea/ChNodeFEAxyz.h"

#include <gtest/gtest.h>

#include <algorithm>
#include <array>
#include <cmath>
#include <memory>
#include <string>
#include <vector>

namespace {

using chrono::ChVector3d;

struct Mode {
    std::string label;
    rcs::Material material;
};

/// One material per element law: J2 (ductile), Mooney-Rivlin (compliant),
/// linear + Rankine (brittle).
std::vector<Mode> modes() {
    return {{"ductile", rcs::aluminum6061()}, {"compliant", rcs::rubber()}, {"brittle", rcs::glass()}};
}

struct Tet {
    std::array<std::shared_ptr<chrono::fea::ChNodeFEAxyz>, 4> nodes;
    rcs::ThermalTet4 element;
};

/// One 0.1 m tet configured the way Part::buildVolume does it.
std::unique_ptr<Tet> makeTet(const rcs::Material& m) {
    auto t = std::make_unique<Tet>();
    const double L = 0.1;
    const ChVector3d at[4] = {{0, 0, 0}, {L, 0, 0}, {0, L, 0}, {0, 0, L}};
    for (int i = 0; i < 4; ++i)
        t->nodes[i] = chrono_types::make_shared<chrono::fea::ChNodeFEAxyz>(at[i]);
    t->element.SetNodes(t->nodes[0], t->nodes[1], t->nodes[2], t->nodes[3]);
    t->element.configure(m.youngsModulus, m.poissonRatio, m.allowPlasticity ? m.yieldStress : 0.0, m.density, 0.0);
    if (m.fracture == rcs::FractureMode::Brittle)
        t->element.setRankine(rcs::rankineTensionLimit(m), rcs::rankineCompressionLimit(m));
    t->element.SetupInitial(nullptr);
    return t;
}

void evaluate(Tet& t) {
    chrono::ChVectorDynamic<> fi(12);
    fi.setZero();
    t.element.ComputeInternalForces(fi);
}

/// Mean stress of a fully restrained linear solid heated by dT:
/// sigma_m = -E alpha dT / (1 - 2 nu).
double restrainedPressure(const rcs::Material& m, double strain) {
    return -m.youngsModulus * strain / std::max(1.0 - 2.0 * m.poissonRatio, 0.02);
}

/// The test_fracture.cpp single-tet part (2 cm edges).
rcs::PartDefinition tetPart(const rcs::Material& material) {
    rcs::PartDefinition part;
    part.name = material.name;
    part.material = material;
    part.thickness = 0.02;
    const double L = 0.02;
    part.volume.vertices = {{0, 0, 0}, {L, 0, 0}, {0, L, 0}, {0, 0, L}};
    part.volume.tets = {{0, 1, 2, 3}};
    part.sizing.cached = true;
    return part;
}

float peakDisplayed(const rcs::Part& part) {
    std::vector<float> out;
    part.evaluateStress(out);
    float peak = 0.0f;
    for (float v : out)
        peak = std::max(peak, v);
    return peak;
}

rcs::PartDefinition patch(const std::string& name, const rcs::Material& m, const ChVector3d& origin, double t) {
    rcs::PartDefinition p;
    p.name = name;
    p.material = m;
    p.thickness = t;
    p.patch.origin = origin;
    p.patch.uAxis = {1, 0, 0};
    p.patch.vAxis = {0, 0, 1};
    p.patch.uSize = 0.1;
    p.patch.vSize = 0.1;
    p.patch.divU = 5;
    p.patch.divV = 5;
    p.sizing.throughThickness = 0;
    p.sizing.cached = true;
    return p;
}

rcs::Scenario quiet(double duration, double dt) {
    rcs::Scenario scene;
    scene.environment.gravity = {0.0, 0.0, 0.0};
    scene.environment.floorEnabled = false;
    scene.environment.aerodynamicsEnabled = false;
    scene.initial.position = {0.0, 0.0, 0.0};
    scene.initial.orientation = {0.0, 0.0, 0.0};
    scene.settings.duration = duration;
    scene.settings.timestep = dt;
    scene.settings.captureRate = 20.0;
    return scene;
}

}  // namespace

// --- Element level: the thermal state reaches every element law -----------

// Free expansion carries no stress. Guard rail for every element law.
// PASSES TODAY (expected).
TEST(ThermalStress, FreeExpansionIsStressFree) {
    for (const Mode& mode : modes()) {
        auto t = makeTet(mode.material);
        const double e = 1.0e-4;
        t->element.setThermalState(e, 1.0, 1.0);
        for (auto& node : t->nodes)
            node->SetPos(node->GetX0() * (1.0 + e));
        evaluate(*t);
        const double scale = std::abs(restrainedPressure(mode.material, e));
        EXPECT_LE(std::abs(t->element.rawMean()), 1.0e-3 * scale) << mode.label;
        EXPECT_LE(t->element.vonMises(), 1.0e-3 * scale) << mode.label;
    }
}

// A fully restrained heated solid carries the textbook thermal pressure,
// -E alpha dT / (1 - 2 nu), and no shear. Exact for the linear law, first
// order for Mooney-Rivlin. PASSES TODAY (expected) -- this is the proof the
// issue asks for that expansion is not ignored by any element law.
TEST(ThermalStress, RestrainedExpansionGivesThermalPressure) {
    for (const Mode& mode : modes()) {
        auto t = makeTet(mode.material);
        const double e = 1.0e-4;
        t->element.setThermalState(e, 1.0, 1.0);
        evaluate(*t);
        const double expected = restrainedPressure(mode.material, e);
        const double tolerance = mode.material.fracture == rcs::FractureMode::Compliant ? 0.15 : 0.01;
        EXPECT_NEAR(t->element.rawMean() / expected, 1.0, tolerance)
            << mode.label << ": mean stress " << t->element.rawMean() << " Pa, expected " << expected << " Pa";
    }
}

// Stiffness retention must scale the element's response. Halving the scale
// halves the restrained thermal pressure. PASSES TODAY (expected).
TEST(ThermalStress, StiffnessRetentionReachesTheElement) {
    for (const Mode& mode : modes()) {
        auto full = makeTet(mode.material);
        auto soft = makeTet(mode.material);
        full->element.setThermalState(1.0e-4, 1.0, 1.0);
        soft->element.setThermalState(1.0e-4, 0.5, 1.0);
        evaluate(*full);
        evaluate(*soft);
        ASSERT_NE(full->element.rawMean(), 0.0) << mode.label;
        EXPECT_NEAR(soft->element.rawMean() / full->element.rawMean(), 0.5, 0.02) << mode.label;
    }
}

// --- Part level: what the user sees ----------------------------------------

// Displayed stress must be the stress the element carries. At 573 K
// aluminium keeps r ~ 0.76 of its modulus; the same small strain should read
// r times the cold stress. Today evaluateStress applies r a second time to
// tet stress, so the hot part reads ~ r^2 (~0.58).
// FAILS TODAY.
TEST(ThermalStress, DisplayedStressIsNotKnockedDownTwice) {
    const rcs::Material al = rcs::aluminum6061();
    const double hot = 573.0;
    auto loaded = [&](double kelvin) {
        rcs::Part part;
        part.build(tetPart(al));
        for (std::size_t n = 0; n < part.nodeCount(); ++n) {
            part.setNodeTemperature(n, kelvin);
            part.setStrainTemperature(n, part.initialTemperature());  // no expansion, only softening
        }
        part.syncThermalMechanics();
        auto& nodes = part.xyzNodes();
        ChVector3d p = nodes[1]->GetPos();
        p.x() *= 1.0001;  // elastic at both temperatures
        nodes[1]->SetPos(p);
        part.refreshTetTrials();
        return static_cast<double>(peakDisplayed(part));
    };
    const double cold = loaded(al.referenceTemperature);
    const double warm = loaded(hot);
    ASSERT_GT(cold, 0.0);
    const double r = rcs::thermalStiffnessRetention(al, hot);
    EXPECT_NEAR(warm / cold, r, 0.02) << "hot/cold displayed stress " << warm / cold << ", retention r = " << r
                                      << ", r^2 = " << r * r;
}

// --- Material tables: driven by fields, no hidden preset cliffs ------------

// Goal 2: degradation follows the material's own fields. Two materials that
// differ only by 2 K of solidus should soften alike. Today 1399 K selects the
// 6061 table and 1401 K the Eurocode steel table (material.hpp:165).
// FAILS TODAY.
TEST(ThermalMaterials, SolidusThresholdIsNotACliff) {
    rcs::Material below = rcs::aluminum6061();
    below.softenTemperature = 600.0;
    below.meltTemperature = 1399.0;
    rcs::Material above = below;
    above.meltTemperature = 1401.0;
    for (double kelvin : {700.0, 900.0, 1100.0}) {
        EXPECT_NEAR(rcs::thermalStiffnessRetention(below, kelvin), rcs::thermalStiffnessRetention(above, kelvin), 0.05)
            << "stiffness at " << kelvin << " K";
        EXPECT_NEAR(rcs::thermalStrengthRetention(below, kelvin), rcs::thermalStrengthRetention(above, kelvin), 0.05)
            << "strength at " << kelvin << " K";
    }
}

// Same idea for strength: yield 0.99 GPa vs 1.01 GPa should not change how
// fast a steel loses strength. Today >= 1 GPa runs the table 1.45x faster
// (material.hpp:217, the AR500 rule).
// FAILS TODAY.
TEST(ThermalMaterials, YieldThresholdIsNotACliff) {
    rcs::Material below = rcs::mildSteel();
    below.yieldStress = 0.99e9;
    rcs::Material above = below;
    above.yieldStress = 1.01e9;
    for (double kelvin : {700.0, 800.0, 900.0}) {
        EXPECT_NEAR(rcs::thermalStrengthRetention(below, kelvin), rcs::thermalStrengthRetention(above, kelvin), 0.05)
            << "strength at " << kelvin << " K";
    }
}

// Every preset softens monotonically from room temperature to well past its
// solidus and is nearly gone there. Guard rail. PASSES TODAY (expected).
TEST(ThermalMaterials, EveryPresetSoftensTowardItsSolidus) {
    for (const rcs::Material& m : rcs::materialPresets()) {
        double lastStiff = 1.0;
        double lastStrength = 1.0;
        EXPECT_NEAR(rcs::thermalStiffnessRetention(m, 293.15), 1.0, 1.0e-9) << m.name;
        EXPECT_NEAR(rcs::thermalStrengthRetention(m, 293.15), 1.0, 1.0e-9) << m.name;
        for (double kelvin = 300.0; kelvin <= m.meltTemperature + 100.0; kelvin += 10.0) {
            const double s = rcs::thermalStiffnessRetention(m, kelvin);
            const double y = rcs::thermalStrengthRetention(m, kelvin);
            EXPECT_LE(s, lastStiff + 1.0e-12) << m.name << " stiffness rises at " << kelvin << " K";
            EXPECT_LE(y, lastStrength + 1.0e-12) << m.name << " strength rises at " << kelvin << " K";
            lastStiff = s;
            lastStrength = y;
        }
        EXPECT_LE(rcs::thermalStrengthRetention(m, m.meltTemperature + 50.0), 0.1) << m.name;
    }
}

// --- Scenario level: the heat field drives the mechanics -------------------

// A part authored hot (initialTemperature) in a room-temperature arena must
// lose heat by convection even when no torch is on. HeatField::advance
// returns early when no source is on and nothing is marked dirty, so the part
// may never cool. A 2 mm aluminium plate loses ~1 K/s at this temperature
// difference; the test asks for 0.1 K over 1 s.
// MAY FAIL TODAY.
TEST(ThermalScenario, HotPartCoolsWithoutASource) {
    for (const Mode& mode : modes()) {
        rcs::PartDefinition plate = patch("Hot plate", mode.material, {-0.05, 0.5, -0.05}, 0.002);
        plate.fixed = true;
        plate.initialTemperature = std::min(600.0, mode.material.softenTemperature);
        rcs::AssemblyDefinition a;
        a.parts.push_back(plate);
        rcs::Simulator sim;
        sim.run(a, quiet(1.0, 5.0e-3), nullptr);
        const double end = sim.assembly().parts().front().temperature();
        EXPECT_LT(end, plate.initialTemperature - 0.1)
            << mode.label << ": started at " << plate.initialTemperature << " K, still " << end << " K after 1 s";
    }
}

// Acceptance: the scenario matrix, continuum modes x weld bridge under a
// torch. Two plates welded edge to edge; a 2000 K torch on the seam must
// release the weld for every material (each solidus is below 2000 K), and
// the same assembly without the torch must hold. PASSES TODAY for aluminium
// (smoke check "a weld through the aluminum solidus breaks"); the other
// modes are unproven, which is what the issue asks to fix.
TEST(ThermalScenario, TorchOnSeamReleasesWeldForEveryMode) {
    std::vector<Mode> matrix = modes();
    matrix.push_back({"ductile steel", rcs::mildSteel()});
    matrix.push_back({"brittle wood", rcs::wood()});
    for (const Mode& mode : matrix) {
        rcs::AssemblyDefinition a;
        a.parts.push_back(patch("Left", mode.material, {-0.1, 0.5, -0.05}, 0.004));
        a.parts.push_back(patch("Right", mode.material, {0.0, 0.5, -0.05}, 0.004));
        rcs::WeldDefinition seam;
        seam.partA = 0;
        seam.partB = 1;
        seam.tolerance = 0.004;
        seam.bead = {{0.0, 0.5, -0.05}, {0.0, 0.5, 0.05}};
        a.welds.push_back(seam);

        rcs::Scenario cold = quiet(0.8, 2.0e-3);
        rcs::Simulator control;
        control.run(a, cold, nullptr);
        EXPECT_EQ(control.assembly().brokenLinkCount(), 0) << mode.label << ": seam broke with no torch";

        rcs::Scenario torch = cold;
        rcs::HeatSource source;
        source.point = {0.0, 0.5, 0.0};
        source.radius = 0.03;
        source.temperature = 2000.0;
        source.startTime = 0.0;
        source.duration = 0.8;
        torch.heatSources.push_back(source);
        rcs::Simulator hot;
        const rcs::Recording rec = hot.run(a, torch, nullptr);
        EXPECT_FALSE(rec.diverged) << mode.label << ": " << rec.message;
        EXPECT_GT(hot.assembly().brokenLinkCount(), 0)
            << mode.label << ": seam under a 2000 K torch held (solidus " << mode.material.meltTemperature << " K)";
    }
}
