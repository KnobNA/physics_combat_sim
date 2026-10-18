// Issue #5: shell-solid and seam coupling.
//
// How a seam is tied today (assembly.cpp):
//   createLink picks the constraint from the node types of the two parts:
//     xyz - xyz       (tets, BST triangle shells)  ChLinkNodeNode: forces the two
//                                                   nodes to the SAME point
//     xyz - xyzrot    (tet or BST to Reissner)     ChLinkNodeFrame: keeps the
//                                                   current offset (rigid arm)
//     xyzrot - xyzrot (Reissner to Reissner)       ChLinkMateGeneric at the
//                                                   midpoint: keeps the offset
//   So only xyz-xyz pairs can start with a gap the solver must close. To
//   avoid that, buildWeld "snaps" xyz-xyz pairs within snapReach together with
//   snapWeldNode, which moves the node AND its rest position (SetX0): the
//   part is permanently reshaped at that node. Pairs beyond snapReach, or
//   whose node was already snapped by another weld, keep their gap and are
//   pulled shut on the first step.
//   Inside one mixed part, the shell edge is tied to the nearest tet node
//   within one thickness with ChLinkNodeNode and no snap (assembly.cpp
//   ~260-275), so those ties can start up to a thickness apart.
//
// The commercial pattern (Abaqus shell-to-solid coupling, MPC/RBE links) is
// to keep both meshes where they are and tie them with a coupling that
// carries the offset. The tests below state that as three checks, run over
// shell-shell, shell-solid and solid-solid seams and over continuum modes:
//   1. no tie starts violated,
//   2. building a seam does not move any part's rest shape,
//   3. an unloaded assembly stays still (no constraint yank, no overload).

#include "physics_core/assembly.hpp"
#include "physics_core/cad.hpp"
#include "physics_core/joint.hpp"
#include "physics_core/material.hpp"
#include "physics_core/mesh_split.hpp"
#include "physics_core/part.hpp"
#include "physics_core/scenario.hpp"
#include "physics_core/simulator.hpp"

#include "chrono/core/ChVector3.h"

#include <gtest/gtest.h>

#include <algorithm>
#include <cmath>
#include <functional>
#include <string>
#include <vector>

namespace {

using chrono::ChVector3d;

rcs::Scenario quietScene() {
    rcs::Scenario scene;
    scene.environment.gravity = {0.0, 0.0, 0.0};
    scene.environment.floorEnabled = false;
    scene.environment.aerodynamicsEnabled = false;
    scene.initial.position = {0.0, 0.0, 0.0};
    scene.initial.orientation = {0.0, 0.0, 0.0};
    scene.settings.timestep = 1.0e-3;
    return scene;
}

rcs::PartDefinition patchPlate(const std::string& name, const rcs::Material& material, double t,
                               const ChVector3d& origin, const ChVector3d& u, const ChVector3d& v, double size,
                               int div) {
    rcs::PartDefinition plate;
    plate.name = name;
    plate.material = material;
    plate.thickness = t;
    plate.patch.origin = origin;
    plate.patch.uAxis = u;
    plate.patch.vAxis = v;
    plate.patch.uSize = size;
    plate.patch.vSize = size;
    plate.patch.divU = div;
    plate.patch.divV = div;
    return plate;
}

/// Closed 0.2 x t x 0.2 box surface; ensurePartMesh turns it into a BST
/// (xyz-node) triangle shell on its midsurface at y = t/2.
rcs::PartDefinition closedThinPlate(const rcs::Material& material, double t) {
    rcs::PartDefinition plate;
    plate.name = "Thin plate";
    plate.material = material;
    plate.thickness = t;
    const double L = 0.20;
    plate.source.vertices = {
        {0, 0, 0}, {L, 0, 0}, {L, t, 0}, {0, t, 0},
        {0, 0, L}, {L, 0, L}, {L, t, L}, {0, t, L},
    };
    const int tris[][3] = {
        {0, 2, 1}, {0, 3, 2}, {4, 5, 6}, {4, 6, 7}, {0, 1, 5}, {0, 5, 4},
        {3, 6, 2}, {3, 7, 6}, {0, 4, 7}, {0, 7, 3}, {1, 2, 6}, {1, 6, 5},
    };
    for (const auto& tri : tris)
        plate.source.triangles.push_back({tri[0], tri[1], tri[2]});
    plate.sizing.minElement = 0.002;
    plate.sizing.maxElement = 0.025;
    plate.sizing.tetBudget = 800;
    return plate;
}

rcs::WeldDefinition faceWeld(double y, double cx, double cz, double half) {
    rcs::WeldDefinition weld;
    weld.partA = 0;
    weld.partB = 1;
    weld.tolerance = 0.004;
    weld.throat = 0.006;
    weld.bead = {{cx - half, y, cz - half}, {cx + half, y, cz - half}, {cx + half, y, cz + half},
                 {cx - half, y, cz + half}};
    return weld;
}

struct Fixture {
    std::string label;
    rcs::AssemblyDefinition assembly;
};

/// Seam fixtures named by idealization pair, for one material. Every
/// fixture is built from fields only; nothing is keyed on a preset.
std::vector<Fixture> fixtures(const rcs::Material& m) {
    std::vector<Fixture> out;

    // solid - solid (or whatever the box panels classify as): the default box.
    {
        rcs::TestCubeOptions options;
        options.material = m;
        out.push_back({"welded box", rcs::makeTestCube(options)});
    }
    // shell (Reissner patch) - solid: pane under a cube, cube bottom on the
    // pane's top skin, 4 mm off the midsurface. Same layout as the
    // checkWeldMixed smoke test.
    {
        rcs::AssemblyDefinition a;
        a.parts.push_back(patchPlate("Pane", m, 0.008, {-0.16, 0.0, -0.16}, {1, 0, 0}, {0, 0, 1}, 0.32, 8));
        rcs::PartDefinition cube = rcs::makeSolidCube(0.08, {0.0, 0.044, 0.0});
        cube.material = rcs::aluminum6061();
        a.parts.push_back(cube);
        a.welds.push_back(faceWeld(0.004, 0.0, 0.0, 0.04));
        out.push_back({"patch shell - solid", a});
    }
    // shell (BST triangles, xyz nodes) - solid: the case that gets a
    // ChLinkNodeNode across the midsurface offset.
    {
        rcs::AssemblyDefinition a;
        a.parts.push_back(closedThinPlate(m, 0.008));
        rcs::PartDefinition cube = rcs::makeSolidCube(0.08, {0.1, 0.048, 0.1});
        cube.material = rcs::aluminum6061();
        a.parts.push_back(cube);
        a.welds.push_back(faceWeld(0.008, 0.1, 0.1, 0.04));
        out.push_back({"triangle shell - solid", a});
    }
    // shell - shell T-joint: a stem standing on the base midsurface, no bead
    // (same layout as catalog check S3).
    {
        rcs::AssemblyDefinition a;
        a.parts.push_back(patchPlate("Base", m, 0.004, {0.0, 0.0, 0.0}, {1, 0, 0}, {0, 0, 1}, 0.16, 8));
        a.parts.push_back(patchPlate("Stem", m, 0.004, {0.08, 0.0, 0.0}, {0, 1, 0}, {0, 0, 1}, 0.16, 8));
        rcs::WeldDefinition weld;
        weld.partA = 0;
        weld.partB = 1;
        weld.tolerance = 0.004;
        a.welds.push_back(weld);
        out.push_back({"shell T-joint", a});
    }
    // One part that is shell and solid at once (thin web on a boss): the
    // in-part shell-solid ties.
    {
        rcs::AssemblyDefinition a;
        rcs::PartDefinition stepped = rcs::makeSteppedSolid();
        stepped.material = m;
        stepped.sizing.cached = false;
        stepped.volume = {};
        stepped.surface = {};
        a.parts.push_back(stepped);
        out.push_back({"mixed shell-solid part", a});
    }
    return out;
}

std::vector<rcs::Material> modes() {
    return {rcs::aluminum6061(), rcs::rubber(), rcs::glass()};
}

/// Mesh exactly as Simulator::reset does, so a test can compare against
/// the parts before buildWeld touches them.
rcs::AssemblyDefinition meshLikeReset(rcs::AssemblyDefinition a) {
    rcs::seedWeldBeads(a);
    for (auto& part : a.parts)
        if (part.simulate)
            rcs::ensurePartMesh(part);
    return a;
}

double linkGap(const rcs::Assembly& assembly, const rcs::JointLink& link) {
    const auto& parts = assembly.parts();
    return (parts[link.partB].nodePosition(link.nodeB) - parts[link.partA].nodePosition(link.nodeA)).Length();
}

}  // namespace

// Check 1 (acceptance: "zero initial violation"). After reset, before any
// step, every node-to-node tie must already be satisfied. A gap here is a
// force the first Newton iterate has to produce out of nothing.
// FAILS TODAY wherever an xyz-xyz pair lies beyond snapReach or its node was
// already used by another weld (box corners are the usual suspects).
TEST(SeamCoupling, NoTieStartsViolated) {
    for (const rcs::Material& m : modes()) {
        for (Fixture& f : fixtures(m)) {
            rcs::Simulator sim;
            sim.reset(f.assembly, quietScene());
            int checked = 0;
            double worst = 0.0;
            std::string where;
            for (const rcs::Joint& joint : sim.assembly().joints()) {
                for (const rcs::JointLink& link : joint.links) {
                    if (!link.pointLink)
                        continue;  // frame and mate links keep their build-time offset
                    ++checked;
                    const double gap = linkGap(sim.assembly(), link);
                    if (gap > worst) {
                        worst = gap;
                        where = joint.label;
                    }
                }
            }
            EXPECT_LE(worst, 1.0e-6) << m.name << " / " << f.label << ": a coincident tie starts "
                                     << worst * 1000.0 << " mm open (" << where << ", " << checked
                                     << " node-node ties checked)";
        }
    }
}

// Check 2 (goal 2: "preserve rest gap"). Building a seam must not reshape the
// parts. snapWeldNode moves a node's rest position (SetX0), so the part's
// unstressed shape changes at every snapped node by up to the snap reach
// (a shell thickness for shell-solid). Compared against the same parts built
// standalone, relative to each part's centroid so the rigid placement
// cancels out.
// FAILS TODAY for the triangle shell - solid fixture, and wherever a snap moved
// a node by more than a micron.
TEST(SeamCoupling, BuildingSeamsDoesNotReshapeParts) {
    for (const rcs::Material& m : modes()) {
        for (Fixture& f : fixtures(m)) {
            const rcs::AssemblyDefinition meshed = meshLikeReset(f.assembly);
            rcs::Simulator sim;
            sim.reset(meshed, quietScene(), false);
            ASSERT_EQ(sim.assembly().parts().size(), meshed.parts.size()) << f.label;
            for (std::size_t i = 0; i < meshed.parts.size(); ++i) {
                rcs::Part reference;
                reference.build(meshed.parts[i]);
                const rcs::Part& welded = sim.assembly().parts()[i];
                const std::size_t n = std::min(reference.nodeCount(), welded.nodeCount());
                ASSERT_GT(n, 0u) << f.label;
                ChVector3d cRef{0, 0, 0};
                ChVector3d cWeld{0, 0, 0};
                for (std::size_t k = 0; k < n; ++k) {
                    cRef += reference.restPosition(k);
                    cWeld += welded.restPosition(k);
                }
                cRef /= static_cast<double>(n);
                cWeld /= static_cast<double>(n);
                double worst = 0.0;
                for (std::size_t k = 0; k < n; ++k) {
                    const ChVector3d a = reference.restPosition(k) - cRef;
                    const ChVector3d b = welded.restPosition(k) - cWeld;
                    worst = std::max(worst, (a - b).Length());
                }
                EXPECT_LE(worst, 1.0e-6) << m.name << " / " << f.label << " / " << meshed.parts[i].name
                                         << ": seam building moved a rest node by " << worst * 1000.0 << " mm";
            }
        }
    }
}

// Check 3 (acceptance: "no first-frame overload release"). With no gravity,
// no floor and no loads, nothing should move and no seam should carry load.
// A started-open tie (check 1) or an in-part shell-solid tie that starts a
// thickness apart shows up here as kinetic energy and weld utilisation.
// 20 ms is more than twice the reaction filter time constant (8 ms).
// MAY FAIL TODAY, in particular for the mixed shell-solid part.
TEST(SeamCoupling, UnloadedAssemblyStaysStill) {
    for (const rcs::Material& m : modes()) {
        for (Fixture& f : fixtures(m)) {
            const rcs::Scenario scene = quietScene();
            rcs::Simulator sim;
            sim.reset(f.assembly, scene);
            const double dt = scene.settings.timestep;
            double time = 0.0;
            bool finite = true;
            for (int i = 0; i < 20 && finite; ++i, time += dt)
                finite = sim.step(dt, time, scene.environment, scene.loads);
            ASSERT_TRUE(finite) << m.name << " / " << f.label;
            const rcs::Assembly& a = sim.assembly();
            EXPECT_LE(a.kineticEnergy(), 1.0e-6)
                << m.name << " / " << f.label << ": unloaded assembly picked up " << a.kineticEnergy() << " J";
            double peak = 0.0;
            for (const rcs::Joint& joint : a.joints())
                peak = std::max(peak, joint.peakUtilisation());
            EXPECT_LE(peak, 0.05) << m.name << " / " << f.label << ": an unloaded seam is at "
                                  << peak * 100.0 << "% of its capacity";
            EXPECT_EQ(a.brokenLinkCount(), 0) << m.name << " / " << f.label;
        }
    }
}
