// Issue #1: shell vs solid classification from geometry, for every material.
//
// These tests describe the policy the issue asks for, written against today's
// API. Several of them FAIL on the current tree on purpose: each failure
// message prints the classification table so it is clear which branch of
// ensurePartMesh made the call. They should go green as the policy lands,
// without any change to the expected values.
//
// Where the decisions live today (all compare thickness / span against 0.15,
// but each one measures "span" differently):
//   mesh_split.cpp faceIsThinWall    t / hypot(region u, v)  (+ t < 3*minElement)
//   mesh_split.cpp meshClosedBody    sorted AABB extents e0/e1, gated by throughThickness
//   mesh_split.cpp fillFromPatch     t / max(uSize, vSize)
//   mesh_split.cpp ensurePartMesh    open sheet: mean t / largest AABB extent
//   mesh_split.cpp seedWeldBeads     t / max(uSize, vSize)
// plus cad.cpp makeTestCube setting sizing.throughThickness = 2, which is the
// only reason the default box panels are tets.

#include "physics_core/cad.hpp"
#include "physics_core/material.hpp"
#include "physics_core/mesh_split.hpp"
#include "physics_core/part.hpp"
#include "physics_core/physics_units.hpp"

#include "chrono/core/ChVector3.h"

#include <gtest/gtest.h>

#include <algorithm>
#include <cmath>
#include <set>
#include <sstream>
#include <string>
#include <vector>

namespace {

enum class Family { Shell, Solid, Mixed, Refused };

const char* familyName(Family f) {
    switch (f) {
        case Family::Shell: return "shell";
        case Family::Solid: return "solid";
        case Family::Mixed: return "mixed";
        case Family::Refused: return "refused";
    }
    return "?";
}

/// What ensurePartMesh produced. A patch shell leaves both meshes empty and
/// Part::build makes a Reissner grid; a triangle shell fills `surface`.
Family familyOf(const rcs::PartDefinition& part, bool ok) {
    if (!ok)
        return Family::Refused;
    const bool solid = !part.volume.empty();
    const bool shell = !part.surface.empty();
    if (solid && shell)
        return Family::Mixed;
    return solid ? Family::Solid : Family::Shell;
}

Family mesh(rcs::PartDefinition& part) {
    rcs::applyResolutionSizing(part, 1.0);
    const bool ok = rcs::ensurePartMesh(part);
    return familyOf(part, ok);
}

/// L x L plate of thickness t as an authored Reissner patch.
rcs::PartDefinition patchPlate(const rcs::Material& material, double L, double t) {
    rcs::PartDefinition plate;
    plate.name = "patch plate";
    plate.material = material;
    plate.thickness = t;
    plate.patch.origin = {0.0, 0.0, 0.0};
    plate.patch.uAxis = {1.0, 0.0, 0.0};
    plate.patch.vAxis = {0.0, 0.0, 1.0};
    plate.patch.uSize = L;
    plate.patch.vSize = L;
    return plate;
}

/// The same plate as a closed 12-triangle box, which is what an OBJ/STL
/// import hands ensurePartMesh. Winding matches thinClosedPlate in
/// test_mesh_geom.cpp. The default 1 m patch is left alone on purpose so the
/// "authored patch" shortcuts do not apply.
rcs::PartDefinition importedPlate(const rcs::Material& material, double L, double t) {
    rcs::PartDefinition plate;
    plate.name = "imported plate";
    plate.material = material;
    plate.thickness = t;
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
    return plate;
}

/// Distinct vertex levels along the axis where the part is thinnest, minus
/// one. For an extruded or tetgen-filled slab this is the number of element
/// layers through the thickness (an upper bound if tetgen adds Steiner points).
int layersThroughThickness(const rcs::VolumeMesh& volume) {
    if (volume.vertices.empty())
        return 0;
    chrono::ChVector3d lo = volume.vertices.front();
    chrono::ChVector3d hi = lo;
    for (const auto& p : volume.vertices) {
        lo = {std::min(lo.x(), p.x()), std::min(lo.y(), p.y()), std::min(lo.z(), p.z())};
        hi = {std::max(hi.x(), p.x()), std::max(hi.y(), p.y()), std::max(hi.z(), p.z())};
    }
    const double extent[3] = {hi.x() - lo.x(), hi.y() - lo.y(), hi.z() - lo.z()};
    const int axis = static_cast<int>(std::min_element(extent, extent + 3) - extent);
    std::set<long long> levels;
    for (const auto& p : volume.vertices)
        levels.insert(std::llround(p[axis] * 1.0e6));
    return static_cast<int>(levels.size()) - 1;
}

double slenderness(const rcs::PartDefinition& part) {
    const auto& v = part.source.vertices;
    if (v.empty())
        return part.thickness / std::max(std::max(part.patch.uSize, part.patch.vSize), 1.0e-9);
    chrono::ChVector3d lo = v.front();
    chrono::ChVector3d hi = lo;
    for (const auto& p : v) {
        lo = {std::min(lo.x(), p.x()), std::min(lo.y(), p.y()), std::min(lo.z(), p.z())};
        hi = {std::max(hi.x(), p.x()), std::max(hi.y(), p.y()), std::max(hi.z(), p.z())};
    }
    double e[3] = {hi.x() - lo.x(), hi.y() - lo.y(), hi.z() - lo.z()};
    std::sort(e, e + 3);
    return e[0] / std::max(e[2], 1.0e-9);
}

/// Every preset whose continuum mode is compressible. Issue #2 owns the
/// near-incompressible ones; for these, geometry alone should pick the family.
std::vector<rcs::Material> compressiblePresets() {
    std::vector<rcs::Material> out;
    for (const rcs::Material& m : rcs::materialPresets()) {
        if (!rcs::volumeHyperelastic(m))
            out.push_back(m);
    }
    return out;
}

const std::vector<double> kThicknessSweep = {0.008, 0.02, 0.04, 0.055, 0.065, 0.075, 0.09, 0.12, 0.20};

}  // namespace

// Phase A. One plate, two representations. A preprocessor classifies the
// geometry, not the data structure it arrived in. Today the patch path uses
// t / max(u, v) while the closed-body path uses region diagonals and AABB
// extents, so the two disagree in a band of thicknesses around t/L = 0.15.
TEST(IdealizationPolicy, PatchAndImportedPlateClassifyTheSame) {
    const double L = 0.40;
    for (const rcs::Material& material : compressiblePresets()) {
        std::ostringstream table;
        bool agree = true;
        for (double t : kThicknessSweep) {
            rcs::PartDefinition patch = patchPlate(material, L, t);
            rcs::PartDefinition imported = importedPlate(material, L, t);
            const Family a = mesh(patch);
            const Family b = mesh(imported);
            table << "  t=" << t << " (t/L=" << t / L << "): patch=" << familyName(a)
                  << " imported=" << familyName(b) << "\n";
            agree = agree && a == b;
        }
        EXPECT_TRUE(agree) << material.name << ": same plate, different element family\n" << table.str();
    }
}

// Phase A. Thickening a wall can turn a shell into a solid, never the reverse.
// A non-monotone sweep means a budget or fallback branch, not geometry, made
// the call (e.g. fillFromPatch returning keepAsShell when the tet budget is hit).
TEST(IdealizationPolicy, ThickerNeverFlipsBackToShell) {
    const double L = 0.40;
    for (const rcs::Material& material : compressiblePresets()) {
        for (bool imported : {false, true}) {
            std::ostringstream table;
            bool seenSolid = false;
            bool monotone = true;
            for (double t : kThicknessSweep) {
                rcs::PartDefinition plate = imported ? importedPlate(material, L, t) : patchPlate(material, L, t);
                const Family f = mesh(plate);
                table << "  t=" << t << ": " << familyName(f) << "\n";
                if (f == Family::Solid)
                    seenSolid = true;
                else if (seenSolid && f == Family::Shell)
                    monotone = false;
            }
            EXPECT_TRUE(monotone) << material.name << (imported ? " (imported)" : " (patch)") << "\n"
                                  << table.str();
        }
    }
}

// Goal 3: "same policy for every material". Aluminium, steel, glass and wood
// are all compressible, so the 8 mm / 400 mm box panels (t/L = 0.02) should get
// the same family whatever they are made of. Today aluminium and steel come out
// as tets (authoredSolid via throughThickness = 2) and glass / wood come out as
// shells (the `plateFabric` brittle branch in meshClosedBody).
TEST(IdealizationPolicy, CompressibleMaterialsShareTheGeometryDecision) {
    std::ostringstream table;
    std::set<Family> families;
    for (const rcs::Material& material : compressiblePresets()) {
        rcs::TestCubeOptions options;
        options.material = material;
        rcs::AssemblyDefinition box = rcs::makeTestCube(options);
        for (auto& part : box.parts)
            rcs::applyResolutionSizing(part, 1.0);
        rcs::seedWeldBeads(box);
        rcs::PartDefinition& panel = box.parts.front();
        const Family f = familyOf(panel, rcs::ensurePartMesh(panel));
        families.insert(f);
        table << "  " << material.name << ": " << familyName(f) << "\n";
    }
    EXPECT_EQ(families.size(), 1u) << "8 mm panel, t/L = 0.02, classified by material:\n" << table.str();
}

// Goal 4: the default box passes because it obeys the policy, not because it
// is special. `sizing.throughThickness` is an element-size knob (MeshSizing);
// it must not also be the switch that picks shell vs solid. Today the same
// panel is a shell at 1 and a tet solid at 2.
TEST(IdealizationPolicy, ThroughThicknessKnobDoesNotPickTheFamily) {
    const rcs::AssemblyDefinition box = rcs::makeTestCube();
    rcs::PartDefinition one = box.parts.front();
    rcs::PartDefinition two = box.parts.front();
    one.sizing.throughThickness = 1;
    two.sizing.throughThickness = 2;
    const Family a = mesh(one);
    const Family b = mesh(two);
    EXPECT_EQ(a, b) << box.parts.front().name << ": throughThickness=1 -> " << familyName(a)
                    << ", throughThickness=2 -> " << familyName(b);
}

// Phase B quality gate. If policy does put a slender wall on linear tets,
// there must be enough layers to carry bending. One linear tet through the
// thickness has no bending mode (the code's own comment in fillFromPatch
// says so and requires 3); tetgen-filled panels have no such guard. Failing
// this is what "fail meshing loudly" should catch.
TEST(IdealizationPolicy, SlenderSolidsHaveLayersThroughThickness) {
    constexpr int kMinLayers = 3;
    rcs::AssemblyDefinition box = rcs::makeTestCube();
    for (auto& part : box.parts)
        rcs::applyResolutionSizing(part, 1.0);
    rcs::seedWeldBeads(box);
    for (auto& part : box.parts) {
        const double slender = slenderness(part);
        ASSERT_TRUE(rcs::ensurePartMesh(part)) << part.name;
        if (part.volume.empty() || slender >= 0.15)
            continue;
        EXPECT_GE(layersThroughThickness(part.volume), kMinLayers)
            << part.name << " is a tet solid with t/L = " << slender
            << " but only " << layersThroughThickness(part.volume) << " element layer(s) through the thickness";
    }
}
