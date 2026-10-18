// Issue #4: brittle fracture toward cohesive / XFEM-class behaviour, not
// element deletion.
//
// What the code does today (part.cpp):
//   shells  shatterQuads -> detachQuad (RemoveElement) -> promoteLooseFragments
//           -> releaseShard: an over-stressed quad leaves the FEA solve and
//           becomes a RigidShard chip. This is element erosion.
//   tets    fractureCluster cuts a stress-driven plane and clones nodes, so
//           both sides can stay live (real face separation). But the positive
//           side is turned into rigid chips whenever
//             shatterToDust || piece volume < 0.125 * minimumShardSize^3
//           so a display-sized floor, not the mesh, decides what is simulated.
//
// Tests are keyed on FractureMode::Brittle plus fields. Presets are only used
// as data. Tests marked FAILS TODAY pin the deletion path; the others are
// guard rails the fix must keep green.

#include "physics_core/material.hpp"
#include "physics_core/part.hpp"
#include "physics_core/physics_units.hpp"

#include "chrono/core/ChVector3.h"

#include <gtest/gtest.h>

#include <string>
#include <utility>
#include <vector>

namespace {

/// Same fixture as test_fracture.cpp: one regular tet, 2 cm edges (V ~ 1.3 cm^3).
rcs::PartDefinition regularTet(const rcs::Material& material) {
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

void stretchTet(rcs::Part& part, double scale) {
    auto& nodes = part.xyzNodes();
    if (nodes.size() < 2 || !nodes[1])
        return;
    chrono::ChVector3d p = nodes[1]->GetPos();
    p.x() *= scale;
    nodes[1]->SetPos(p);
    part.refreshTetTrials();
}

/// Brittle by fields: past Rankine at 1.5x stretch, no Griffith gate, and no
/// display floor. Everything else is soda-lime.
rcs::Material crackable() {
    rcs::Material m = rcs::glass();
    m.name = "Brittle continuum";
    m.fractureStrength = 1.0e6;
    m.fractureToughness = 0.0;
    m.minimumShardSize = 0.0;
    m.shatterToDust = false;
    return m;
}

struct Outcome {
    int tets = 0;
    int liveTets = 0;
    std::size_t shards = 0;
    double massBefore = 0.0;
    double massAfter = 0.0;
};

Outcome crackOneTet(const rcs::Material& material) {
    rcs::Part part;
    part.build(regularTet(material));
    part.refreshTetTrials();
    Outcome out;
    out.massBefore = part.mass();
    stretchTet(part, 1.5);
    part.tearFailed();
    out.tets = part.tetCount();
    out.liveTets = part.liveTetCount();
    out.shards = part.shards().size();
    out.massAfter = part.mass();
    return out;
}

/// 0.2 m square brittle plate on a 4 x 4 Reissner patch.
rcs::PartDefinition brittlePlate(const rcs::Material& material) {
    rcs::PartDefinition plate;
    plate.name = material.name;
    plate.material = material;
    plate.thickness = 0.006;
    plate.patch.origin = {0.0, 0.0, 0.0};
    plate.patch.uAxis = {1.0, 0.0, 0.0};
    plate.patch.vAxis = {0.0, 0.0, 1.0};
    plate.patch.uSize = 0.2;
    plate.patch.vSize = 0.2;
    plate.patch.divU = 4;
    plate.patch.divV = 4;
    plate.sizing.throughThickness = 0;
    plate.sizing.cached = true;
    return plate;
}

}  // namespace

// Control. With no display floor the tet path already does what Phase B asks:
// the crack plane cuts the tet, nodes on the plane are cloned, and both sides
// stay in the FEA solve. PASSES TODAY.
TEST(BrittleFracture, CrackSeparatesFacesAndKeepsBothSidesLive) {
    const Outcome out = crackOneTet(crackable());
    EXPECT_GT(out.tets, 1) << "the crack plane did not cut the tet";
    EXPECT_EQ(out.liveTets, out.tets) << "a cut child left the solve";
    EXPECT_EQ(out.shards, 0u) << "a resolved fragment was turned into a rigid chip";
}

// Goal 2 / non-goal "shard floors dominate". Setting a display-scale floor
// (the glass preset's 2 cm) or shatterToDust must not change which resolved
// pieces are simulated. Today either one rigidizes the positive side
// (fractureCluster `erase`), so the same crack deletes half the body.
// FAILS TODAY for both variants.
TEST(BrittleFracture, ShardFloorDoesNotDecideWhatIsSimulated) {
    const Outcome reference = crackOneTet(crackable());

    rcs::Material floored = crackable();
    floored.minimumShardSize = rcs::glass().minimumShardSize;
    rcs::Material dust = crackable();
    dust.shatterToDust = true;

    for (const auto& variant : {std::make_pair(std::string("minimumShardSize=0.02"), floored),
                                std::make_pair(std::string("shatterToDust"), dust)}) {
        const Outcome out = crackOneTet(variant.second);
        EXPECT_EQ(out.liveTets, reference.liveTets)
            << variant.first << ": " << out.liveTets << " of " << out.tets << " tets still live, "
            << out.shards << " rigid chips (reference keeps " << reference.liveTets << " live)";
    }
}

// Phase B for shells. A brittle shell past Rankine should open a crack between
// elements (duplicate nodes on the crack line) and keep both sides as live
// shells. Today shatterQuads detaches every over-stressed quad and
// promoteLooseFragments turns it into a rigid chip.
// FAILS TODAY.
TEST(BrittleFracture, ShellCrackSeparatesInsteadOfDeleting) {
    rcs::Material m = rcs::glass();
    m.fractureToughness = 0.0;
    m.shatterToDust = false;
    rcs::Part plate;
    plate.build(brittlePlate(m));
    const int total = plate.totalElementCount();
    const std::size_t nodesBefore = plate.nodeCount();
    ASSERT_GT(total, 0);

    // 1% uniform stretch along u: ~0.7 GPa membrane stress, far past Rankine.
    for (std::size_t n = 0; n < plate.nodeCount(); ++n) {
        chrono::ChVector3d p = plate.nodePosition(n);
        p.x() *= 1.01;
        plate.setNodeKinematics(n, p, {0, 0, 0});
    }
    ASSERT_TRUE(plate.tearFailed()) << "plate did not fracture at 1% strain";
    EXPECT_GT(plate.nodeCount(), nodesBefore) << "no crack faces were created";
    EXPECT_EQ(plate.liveElementCount(), total)
        << total - plate.liveElementCount() << " of " << total << " shell elements were deleted";
    EXPECT_TRUE(plate.shards().empty()) << plate.shards().size() << " quads became rigid chips";
}

// Goal 3: fabric enters as fields, not as preset routines. A tempered pane
// renamed to anything must fracture identically. PASSES TODAY (there is no
// name check in part.cpp); keep it that way.
TEST(BrittleFracture, PresetsAreOnlyData) {
    rcs::Material tempered = rcs::temperedGlass();
    tempered.fractureStrength = 1.0e6;
    tempered.fractureToughness = 0.0;
    rcs::Material renamed = tempered;
    renamed.name = "Anything";
    const Outcome a = crackOneTet(tempered);
    const Outcome b = crackOneTet(renamed);
    EXPECT_EQ(a.tets, b.tets);
    EXPECT_EQ(a.liveTets, b.liveTets);
    EXPECT_EQ(a.shards, b.shards);
}

// "Multiple live fragments; stable solve" needs fracture not to create or
// destroy mass. Node cloning moves a quarter of the parent's mass to the clone.
// PASSES TODAY as far as Part::mass() goes; see the PR notes on Assembly
// kineticEnergy / massHeight, which add shard.mass on top of the shard's node
// masses.
TEST(BrittleFracture, FractureConservesMass) {
    rcs::Material floored = crackable();
    floored.minimumShardSize = 0.02;
    for (const rcs::Material& m : {crackable(), floored}) {
        const Outcome out = crackOneTet(m);
        ASSERT_GT(out.massBefore, 0.0);
        EXPECT_NEAR(out.massAfter / out.massBefore, 1.0, 1.0e-9)
            << "minimumShardSize=" << m.minimumShardSize;
    }
}
