// Issue #3: dissimilar-material joints. Interface capacity should be a
// property of the interface (throat, weld/adhesive class, initiation stress,
// Gc), not min(parent bulk allowables) on coincident nodes.
//
// Where it lives today (assembly.cpp, Assembly::buildWeld):
//   allowable       = min(a.weldShearStrength, b.weldShearStrength) * strengthScale
//   breakPeel       = (peelStrength or allowable) * area
//   criticalOpening = cohesiveCriticalOpening(max(allowable, peelPa),
//                         std::max(Gc_a, Gc_b), throat)   <- named `weakerGc`
// Links that union-find defers are re-created by Assembly::respan through
// createLink with only breakForce / breakTorque, so their breakPeel and
// criticalOpening stay 0 and updateJoints skips their cohesive check.
//
// Tests marked FAILS TODAY pin a specific gap. Tests marked MAY FAIL TODAY
// encode the reported symptom (soft-hard seams popping under mild load) and
// should be run to confirm.

#include "physics_core/assembly.hpp"
#include "physics_core/cad.hpp"
#include "physics_core/joint.hpp"
#include "physics_core/material.hpp"
#include "physics_core/physics_units.hpp"
#include "physics_core/scenario.hpp"
#include "physics_core/simulator.hpp"

#include <gtest/gtest.h>

#include <algorithm>
#include <string>
#include <vector>

namespace {

rcs::Scenario restingScene(double duration) {
    rcs::Scenario scene;
    scene.environment.gravity = {0.0, -9.81, 0.0};
    scene.environment.floorEnabled = true;
    scene.environment.aerodynamicsEnabled = false;
    scene.settings.duration = duration;
    scene.settings.timestep = 2.0e-3;
    scene.settings.captureRate = 20.0;
    scene.settings.meshResolution = 0.25;
    scene.initial.position = {0.0, 0.201, 0.0};
    scene.initial.orientation = {0.0, 0.0, 0.0};
    return scene;
}

/// Default welded box with one panel swapped to `odd`. The swap is by field
/// values only; nothing downstream sees a preset name.
rcs::AssemblyDefinition boxWithOddPanel(const rcs::Material& odd, const std::string& panel = "Top") {
    rcs::AssemblyDefinition box = rcs::makeTestCube();
    for (auto& part : box.parts) {
        if (part.name != panel)
            continue;
        part.material = odd;
        part.sizing.cached = false;
        part.volume = {};
        part.surface = {};
    }
    return box;
}

struct WeakPartner {
    std::string label;
    rcs::Material material;
};

/// Weak partners parameterized by continuum mode and by strength ratio. The
/// last one is aluminium with a 100x weaker weld rating: same mode, same
/// stiffness, only the ratio changes.
std::vector<WeakPartner> weakPartners() {
    rcs::Material weakWeldAluminium = rcs::aluminum6061();
    weakWeldAluminium.name = "Aluminium, weld rating / 100";
    weakWeldAluminium.weldShearStrength /= 100.0;
    return {
        {"compliant", rcs::rubber()},
        {"brittle", rcs::glass()},
        {"ductile, 1/100 weld rating", weakWeldAluminium},
    };
}

}  // namespace

// Phase B building block. The cohesive opening of an interface cannot exceed
// what the weaker side can release. buildWeld takes std::max of the two G_c
// values (the tougher side) while naming the result `weakerGc`.
// FAILS TODAY whenever the two parents have different, non-zero toughness.
TEST(WeldCohesive, CriticalOpeningUsesTheWeakerToughness) {
    rcs::Material tough = rcs::aluminum6061();
    tough.name = "Tough metal";
    tough.fractureToughness = 30.0e6;  // G_c ~ 1.2e4 J/m^2
    const rcs::Material brittle = rcs::glass();  // G_c ~ 7.6 J/m^2
    ASSERT_GT(rcs::modeIReleaseRate(tough), 100.0 * rcs::modeIReleaseRate(brittle));

    rcs::TestCubeOptions options;
    options.material = tough;
    rcs::AssemblyDefinition box = rcs::makeTestCube(options);
    for (auto& part : box.parts) {
        if (part.name == "Top") {
            part.material = brittle;
            part.sizing.cached = false;
        }
    }
    rcs::Simulator sim;
    sim.reset(box, restingScene(0.01));

    const auto& parts = sim.assembly().parts();
    int checked = 0;
    for (const rcs::Joint& joint : sim.assembly().joints()) {
        if (joint.kind != rcs::JointKind::Weld)
            continue;
        const rcs::Material& a = parts[joint.partA].material();
        const rcs::Material& b = parts[joint.partB].material();
        if (a.fractureToughness == b.fractureToughness)
            continue;
        const double weaker = std::min(rcs::modeIReleaseRate(a), rcs::modeIReleaseRate(b));
        const double allowable = std::min(a.weldShearStrength, b.weldShearStrength);
        const double expected = rcs::cohesiveCriticalOpening(allowable, weaker, 0.0);
        for (const rcs::JointLink& link : joint.links) {
            ++checked;
            EXPECT_LE(link.criticalOpening, expected * 1.001)
                << joint.label << ": opening " << link.criticalOpening << " m is set by the tougher parent; "
                << "the weaker side allows only " << expected << " m";
        }
    }
    EXPECT_GT(checked, 0) << "no mixed-toughness weld links were built";
}

// Phase A/B. Every live weld link carries the same failure model. Links that
// respan() re-creates after a break come back through createLink without
// breakPeel / criticalOpening, so the seam's cohesive check silently skips
// them. FAILS TODAY once any link breaks and a deferred link is respanned.
TEST(WeldCohesive, RespannedLinksKeepPeelAndCohesiveData) {
    rcs::AssemblyDefinition box = rcs::makeTestCube();
    for (auto& weld : box.welds)
        weld.strengthScale = 1.0e-4;  // weak enough that seating breaks links
    rcs::Simulator sim;
    rcs::Scenario scene = restingScene(0.1);
    scene.initial.position = {0.0, 0.25, 0.0};
    sim.reset(box, scene);
    std::size_t initialLinks = 0;
    for (const auto& joint : sim.assembly().joints())
        initialLinks += joint.links.size();

    double time = 0.0;
    const double dt = scene.settings.timestep;
    for (int i = 0; i < 50; ++i, time += dt)
        ASSERT_TRUE(sim.step(dt, time, scene.environment, scene.loads)) << "step " << i;

    std::size_t finalLinks = 0;
    for (const auto& joint : sim.assembly().joints())
        finalLinks += joint.links.size();
    if (finalLinks <= initialLinks)
        GTEST_SKIP() << "no deferred link was respanned in this run; raise the drop or lower strengthScale";

    for (const rcs::Joint& joint : sim.assembly().joints()) {
        if (joint.kind != rcs::JointKind::Weld)
            continue;
        for (std::size_t i = 0; i < joint.links.size(); ++i) {
            const rcs::JointLink& link = joint.links[i];
            if (link.broken)
                continue;
            EXPECT_GT(link.breakPeel, 0.0) << joint.label << " link " << i;
            EXPECT_GT(link.criticalOpening, 0.0) << joint.label << " link " << i;
        }
    }
}

// Phase A acceptance: "mixed weak-strong pairs survive gravity seating".
// A 400 mm box resting on the floor puts a few tens of newtons through each
// seam; every partner below is rated in kN. A break here is a reaction spike
// from coincident ties between unlike continua, not a real overload.
// MAY FAIL TODAY: this is the symptom the issue reports.
TEST(MixedJoints, WeakStrongSeamsSurviveGravitySeating) {
    for (const WeakPartner& partner : weakPartners()) {
        rcs::Simulator sim;
        const rcs::Recording rec = sim.run(boxWithOddPanel(partner.material), restingScene(0.1), nullptr);
        EXPECT_FALSE(rec.diverged) << partner.label << ": " << rec.message;
        EXPECT_EQ(sim.assembly().brokenLinkCount(), 0)
            << partner.label << " top panel lost weld links while resting on the floor";
    }
}

// Control for the test above. "Homogeneous structural seams still meet
// existing drop survival regressions." PASSES TODAY.
TEST(MixedJoints, HomogeneousSeamsSurviveGravitySeating) {
    rcs::Simulator sim;
    const rcs::Recording rec = sim.run(rcs::makeTestCube(), restingScene(0.1), nullptr);
    EXPECT_FALSE(rec.diverged) << rec.message;
    EXPECT_EQ(sim.assembly().brokenLinkCount(), 0);
}

// Phase A: capacity is an interface property, optionally bounded by the
// parents. Today the shear allowable is min(parent weldShearStrength), so a
// weld's rating moves when an unrelated bulk property of one parent moves.
// This test changes only the stronger parent and expects no change.
// PASSES TODAY (min already ignores the stronger side); it is the guard rail
// that keeps the fix from coupling capacity to the stronger parent instead.
TEST(WeldCapacity, StrongerParentDoesNotSetTheRating) {
    auto topSeamCapacity = [](double strongShear) {
        rcs::Material strong = rcs::aluminum6061();
        strong.weldShearStrength = strongShear;
        rcs::TestCubeOptions options;
        options.material = strong;
        rcs::AssemblyDefinition box = rcs::makeTestCube(options);
        for (auto& part : box.parts) {
            if (part.name == "Top") {
                part.material = rcs::rubber();
                part.sizing.cached = false;
            }
        }
        rcs::Simulator sim;
        sim.reset(box, restingScene(0.01));
        double total = 0.0;
        const auto& parts = sim.assembly().parts();
        for (const rcs::Joint& joint : sim.assembly().joints()) {
            if (parts[joint.partA].name() != "Top" && parts[joint.partB].name() != "Top")
                continue;
            for (const rcs::JointLink& link : joint.links)
                total += link.breakForce;
        }
        return total;
    };
    const double base = topSeamCapacity(1.65e8);
    const double stronger = topSeamCapacity(3.30e8);
    ASSERT_GT(base, 0.0);
    EXPECT_NEAR(stronger / base, 1.0, 1.0e-9);
}

// ---------------------------------------------------------------------------
// Target API. Compiles only once an interface definition exists, so it does
// not break the build today. Suggested shape:
//
//   // physics_core/interface.hpp
//   struct InterfaceDefinition {
//       double shearStrength = 0.0;    // Pa, 0 = derive from parents
//       double peelStrength = 0.0;     // Pa
//       double releaseRate = 0.0;      // G_c, J/m^2, 0 = derive (weaker parent)
//       bool boundedByParents = true;  // cap at the weaker parent's allowable
//   };
//   // WeldDefinition gains `InterfaceDefinition interface;`
// ---------------------------------------------------------------------------
#if __has_include("physics_core/interface.hpp")
#include "physics_core/interface.hpp"

TEST(WeldCapacity, ExplicitInterfaceRatingIsUsedForAnyMaterialPair) {
    for (const WeakPartner& partner : weakPartners()) {
        rcs::AssemblyDefinition box = boxWithOddPanel(partner.material);
        for (auto& weld : box.welds) {
            weld.interface.shearStrength = 5.0e6;
            weld.interface.boundedByParents = false;
        }
        rcs::Simulator sim;
        sim.reset(box, restingScene(0.01));
        const auto& parts = sim.assembly().parts();
        for (const rcs::Joint& joint : sim.assembly().joints()) {
            const double throat = std::min(parts[joint.partA].thickness(), parts[joint.partB].thickness());
            for (const rcs::JointLink& link : joint.links) {
                // breakForce / (allowable * throat) is the tributary length, which
                // must be positive and no longer than the panel edge.
                const double length = link.breakForce / (5.0e6 * throat);
                EXPECT_GT(length, 0.0) << partner.label << " " << joint.label;
                EXPECT_LE(length, 0.41) << partner.label << " " << joint.label;
            }
        }
    }
}
#endif
