#include "physics_core/physics_units.hpp"
#include "physics_core/material.hpp"
#include "physics_core/mesh_split.hpp"
#include "physics_core/cad.hpp"
#include "physics_core/assembly.hpp"
#include "physics_core/simulator.hpp"

#include "chrono/core/ChVector3.h"

#include <gtest/gtest.h>

#include <algorithm>
#include <array>
#include <cmath>
#include <vector>

namespace {

rcs::SurfaceMesh unitCube() {
    rcs::SurfaceMesh mesh;
    mesh.vertices = {
        {0, 0, 0}, {1, 0, 0}, {1, 1, 0}, {0, 1, 0},
        {0, 0, 1}, {1, 0, 1}, {1, 1, 1}, {0, 1, 1},
    };
    mesh.triangles = {
        {0, 2, 1}, {0, 3, 2}, {4, 5, 6}, {4, 6, 7},
        {0, 1, 5}, {0, 5, 4}, {1, 2, 6}, {1, 6, 5},
        {2, 3, 7}, {2, 7, 6}, {3, 0, 4}, {3, 4, 7},
    };
    return mesh;
}

}  // namespace

TEST(meshVolume, UnitCube) {
    EXPECT_NEAR(rcs::meshVolume(unitCube()), 1.0, 1.0e-12);
}

TEST(surfaceVolume, UnitCubeAbsolute) {
    EXPECT_NEAR(rcs::surfaceVolume(unitCube()), 1.0, 1.0e-12);
}

TEST(closedEnough, ClosedCubeAndOpenTriangle) {
    EXPECT_TRUE(rcs::closedEnough(unitCube()));
    rcs::SurfaceMesh open;
    open.vertices = {{0, 0, 0}, {1, 0, 0}, {0, 1, 0}};
    open.triangles = {{0, 1, 2}};
    EXPECT_FALSE(rcs::closedEnough(open));
}

TEST(faceArea, RightTriangleHalf) {
    rcs::SurfaceMesh mesh;
    mesh.vertices = {{0, 0, 0}, {1, 0, 0}, {0, 1, 0}};
    mesh.triangles = {{0, 1, 2}};
    EXPECT_NEAR(rcs::faceArea(mesh, mesh.triangles[0]), 0.5, 1.0e-12);
}

TEST(budgetEdge, CubeWithDefaultBudget) {
    rcs::MeshSizing sizing;
    sizing.tetBudget = 2500;
    sizing.minElement = 0.002;
    const double expected = std::max(std::cbrt(6.0 / 2500.0), 0.002);
    EXPECT_NEAR(rcs::budgetEdge(1.0, sizing), expected, 1.0e-12);
}

TEST(targetEdge, ScalesWithLengthAtFixedBudget) {
    rcs::MeshSizing sizing;
    sizing.tetBudget = 2500;
    sizing.minElement = 1.0e-9;
    const double h1 = rcs::targetEdge(1.0, {1.0, 1.0, 1.0}, sizing);
    const double h2 = rcs::targetEdge(8.0, {2.0, 2.0, 2.0}, sizing);
    EXPECT_NEAR(h2 / h1, 2.0, 1.0e-12);
    const double n1 = 6.0 / (h1 * h1 * h1);
    const double n2 = 48.0 / (h2 * h2 * h2);
    EXPECT_NEAR(n1 / n2, 1.0, 1.0e-9);
}

TEST(targetEdge, ElongatedBarTracksLongSpan) {
    rcs::MeshSizing sizing;
    sizing.tetBudget = 2500;
    sizing.minElement = 1.0e-9;
    const double h1 = rcs::targetEdge(8.0, {8.0, 1.0, 1.0}, sizing);
    const double h2 = rcs::targetEdge(16.0, {16.0, 2.0, 2.0}, sizing);
    EXPECT_NEAR(h2 / h1, 2.0, 1.0e-12);
}

TEST(targetEdge, CompliantDoesNotRefineIntoLock) {
    rcs::MeshSizing sizing;
    sizing.tetBudget = 2500;
    sizing.minElement = 1.0e-9;
    rcs::Material stiff = rcs::rubber();
    stiff.poissonRatio = 0.30;
    rcs::Material lock = rcs::rubber();
    lock.poissonRatio = 0.49;
    const double hStiff = rcs::targetEdge(1.0, {1.0, 1.0, 1.0}, sizing, stiff);
    const double hLock = rcs::targetEdge(1.0, {1.0, 1.0, 1.0}, sizing, lock);
    EXPECT_NEAR(hLock, hStiff, 1.0e-12);
}

TEST(tetBudgetCap, CompliantDoesNotInflateCount) {
    rcs::MeshSizing sizing;
    sizing.tetBudget = 400;
    rcs::Material metal = rcs::aluminum6061();
    rcs::Material lock = rcs::rubber();
    EXPECT_EQ(rcs::tetBudgetCap(sizing, metal), 400);
    EXPECT_EQ(rcs::tetBudgetCap(sizing, lock), rcs::tetBudgetCap(sizing, metal));
}

TEST(worstTetQuality, RegularTetNearPointOneOneSeven) {
    rcs::VolumeMesh mesh;
    mesh.vertices = {{1, 1, 1}, {1, -1, -1}, {-1, 1, -1}, {-1, -1, 1}};
    mesh.tets = {{0, 1, 2, 3}};
    EXPECT_NEAR(rcs::worstTetQuality(mesh), 0.117851, 2.0e-4);
}

TEST(worstTetQuality, SliverNearZero) {
    rcs::VolumeMesh mesh;
    mesh.vertices = {{0, 0, 0}, {1, 0, 0}, {0, 1, 0}, {0.5, 0.5, 1.0e-6}};
    mesh.tets = {{0, 1, 2, 3}};
    EXPECT_LT(rcs::worstTetQuality(mesh), 1.0e-5);
}

TEST(closestPointOnTriangle, ProjectsToInterior) {
    const chrono::ChVector3d q =
        rcs::closestPointOnTriangle({0.25, 0.25, 1.0}, {0, 0, 0}, {1, 0, 0}, {0, 1, 0});
    EXPECT_NEAR(q.x(), 0.25, 1.0e-12);
    EXPECT_NEAR(q.y(), 0.25, 1.0e-12);
    EXPECT_NEAR(q.z(), 0.0, 1.0e-12);
}

TEST(barycentric, CentroidIsThirds) {
    double wa = 0.0;
    double wb = 0.0;
    double wc = 0.0;
    rcs::barycentric({1.0 / 3.0, 1.0 / 3.0, 0.0}, {0, 0, 0}, {1, 0, 0}, {0, 1, 0}, wa, wb, wc);
    EXPECT_NEAR(wa, 1.0 / 3.0, 1.0e-12);
    EXPECT_NEAR(wb, 1.0 / 3.0, 1.0e-12);
    EXPECT_NEAR(wc, 1.0 / 3.0, 1.0e-12);
}

TEST(geometricEdgeRequest, ThicknessAndCurvatureClamp) {
    EXPECT_NEAR(rcs::geometricEdgeRequest(0.009, 0.0), 0.003, 1.0e-12);
    EXPECT_NEAR(rcs::geometricEdgeRequest(0.0, 10.0), 0.05, 1.0e-12);
    EXPECT_GT(rcs::geometricEdgeRequest(0.0, 0.0), 1.0e20);
    EXPECT_LT(rcs::geometricEdgeRequest(0.03, 100.0), rcs::geometricEdgeRequest(0.03, 0.0));
}

TEST(clampSizeField, KeepsLocalGeometryUnderCeiling) {
    std::vector<double> field = {0.01, 0.02, 0.03};
    rcs::clampSizeField(field, 0.05, 0.001, 0.05);
    EXPECT_NEAR(field[0], 0.01, 1.0e-12);
    EXPECT_NEAR(field[1], 0.02, 1.0e-12);
    EXPECT_NEAR(field[2], 0.03, 1.0e-12);
}

TEST(clampSizeField, BudgetCeilingClipsCoarseEntries) {
    std::vector<double> field = {0.01, 0.02, 0.03};
    rcs::clampSizeField(field, 0.015, 0.001, 0.05);
    EXPECT_NEAR(field[0], 0.01, 1.0e-12);
    EXPECT_NEAR(field[1], 0.015, 1.0e-12);
    EXPECT_NEAR(field[2], 0.015, 1.0e-12);
    std::vector<double> floor = {0.01, 0.02, 0.03};
    rcs::clampSizeField(floor, 0.001, 0.002, 0.05);
    EXPECT_GE(*std::min_element(floor.begin(), floor.end()), 0.002);
}

TEST(clampSizeField, UnknownBecomesTarget) {
    std::vector<double> field = {1.0e30, 0.0};
    rcs::clampSizeField(field, 0.01, 0.001, 0.05);
    EXPECT_NEAR(field[0], 0.01, 1.0e-12);
    EXPECT_NEAR(field[1], 0.01, 1.0e-12);
}

TEST(makeSteppedSolid, ClosedWatertight) {
    EXPECT_TRUE(rcs::closedEnough(rcs::makeSteppedSolid().source));
    EXPECT_GT(std::abs(rcs::meshVolume(rcs::makeSteppedSolid().source)), 1.0e-6);
}

TEST(makeSteppedSolid, WebFinerThanBossAndResolutionScalesBoth) {
    auto meanRegion = [](const rcs::VolumeMesh& volume, double xLo, double xHi) {
        double sum = 0.0;
        int count = 0;
        for (const auto& tet : volume.tets) {
            chrono::ChVector3d c{0, 0, 0};
            for (int i = 0; i < 4; ++i)
                c += volume.vertices[static_cast<std::size_t>(tet[static_cast<std::size_t>(i)])];
            c *= 0.25;
            if (c.x() < xLo || c.x() > xHi)
                continue;
            const int edges[6][2] = {{0, 1}, {0, 2}, {0, 3}, {1, 2}, {1, 3}, {2, 3}};
            for (const auto& e : edges) {
                sum += (volume.vertices[static_cast<std::size_t>(tet[static_cast<std::size_t>(e[0])])] -
                        volume.vertices[static_cast<std::size_t>(tet[static_cast<std::size_t>(e[1])])])
                           .Length();
                ++count;
            }
        }
        return count > 0 ? sum / static_cast<double>(count) : 0.0;
    };
    auto meshAt = [](double resolution) {
        rcs::PartDefinition part = rcs::makeSteppedSolid();
        rcs::applyResolutionSizing(part, resolution);
        EXPECT_TRUE(rcs::ensurePartMesh(part));
        return part;
    };
    const rcs::PartDefinition coarse = meshAt(1.0);
    const rcs::PartDefinition fine = meshAt(4.0);
    const double boss1 = meanRegion(coarse.volume, -1.0, 0.10);
    const double web1 = meanRegion(coarse.volume, 0.13, 1.0);
    const double boss4 = meanRegion(fine.volume, -1.0, 0.10);
    const double web4 = meanRegion(fine.volume, 0.13, 1.0);
    if (web1 > 1.0e-9 && boss1 > 1.0e-9) {
        EXPECT_LT(web1, boss1 * 0.85);
        EXPECT_LT(boss4, boss1 * 0.95);
        EXPECT_LT(web4, boss4 * 0.95);
        EXPECT_LE(web4, web1 * 1.05);
    } else {
        EXPECT_FALSE(coarse.surface.empty());
        EXPECT_FALSE(coarse.volume.empty());
        EXPECT_GT(fine.volume.tets.size() + fine.surface.triangles.size(),
                  coarse.volume.tets.size() + coarse.surface.triangles.size());
    }
}

TEST(makeTestCube, PanelsStayTetrahedralSolids) {
    rcs::AssemblyDefinition box = rcs::makeTestCube();
    for (auto& part : box.parts)
        rcs::applyResolutionSizing(part, 1.0);
    rcs::seedWeldBeads(box);
    const double inner = 0.5 * (0.40 - 0.008) - 0.004;
    for (auto& part : box.parts) {
        ASSERT_TRUE(rcs::ensurePartMesh(part)) << part.name;
        EXPECT_FALSE(part.volume.empty()) << part.name;
        EXPECT_LT(part.volume.tets.size(), 8000u) << part.name;
        for (const auto& tet : part.volume.tets) {
            chrono::ChVector3d center{0, 0, 0};
            for (int node : tet)
                center += part.volume.vertices[static_cast<std::size_t>(node)];
            center *= 0.25;
            EXPECT_FALSE(std::abs(center.x()) < inner && std::abs(center.y()) < inner &&
                         std::abs(center.z()) < inner)
                << part.name;
        }
    }
}

TEST(ensurePartMesh, SlenderHyperelasticStaysShell) {
    rcs::PartDefinition sheet;
    sheet.name = "Rubber sheet";
    sheet.material = rcs::rubber();
    sheet.thickness = 0.008;
    sheet.patch.origin = {-0.2, 0.0, -0.2};
    sheet.patch.uAxis = {1.0, 0.0, 0.0};
    sheet.patch.vAxis = {0.0, 0.0, 1.0};
    sheet.patch.uSize = 0.4;
    sheet.patch.vSize = 0.4;
    sheet.patch.divU = 4;
    sheet.patch.divV = 4;
    rcs::applyResolutionSizing(sheet, 1.0);
    ASSERT_TRUE(rcs::ensurePartMesh(sheet));
    EXPECT_TRUE(sheet.volume.empty());
    EXPECT_EQ(sheet.sizing.throughThickness, 0);
}

TEST(ensurePartMesh, BulkyHyperelasticStaysTets) {
    rcs::PartDefinition cube = rcs::makeSolidCube(0.08);
    cube.material = rcs::rubber();
    cube.sizing.cached = false;
    cube.volume = {};
    cube.surface = {};
    rcs::applyResolutionSizing(cube, 1.0);
    ASSERT_TRUE(rcs::ensurePartMesh(cube));
    EXPECT_FALSE(cube.volume.empty());
    EXPECT_TRUE(cube.surface.empty());
}

TEST(tetFill, PreservesImportedSurfaceVertices) {
    rcs::PartDefinition cube = rcs::makeSolidCube(0.08);
    const rcs::SurfaceMesh plc = cube.source;
    rcs::applyResolutionSizing(cube, 1.0);
    ASSERT_TRUE(rcs::ensurePartMesh(cube));
    ASSERT_FALSE(cube.volume.empty());
    chrono::ChVector3d plcLo = plc.vertices.front();
    chrono::ChVector3d plcHi = plcLo;
    for (const auto& p : plc.vertices) {
        plcLo = {std::min(plcLo.x(), p.x()), std::min(plcLo.y(), p.y()), std::min(plcLo.z(), p.z())};
        plcHi = {std::max(plcHi.x(), p.x()), std::max(plcHi.y(), p.y()), std::max(plcHi.z(), p.z())};
        bool found = false;
        for (const auto& q : cube.volume.vertices) {
            if ((p - q).Length() < 1.0e-9) {
                found = true;
                break;
            }
        }
        EXPECT_TRUE(found) << "tetgen moved or dropped a PLC vertex";
    }
    chrono::ChVector3d volLo = cube.volume.vertices.front();
    chrono::ChVector3d volHi = volLo;
    for (const auto& p : cube.volume.vertices) {
        volLo = {std::min(volLo.x(), p.x()), std::min(volLo.y(), p.y()), std::min(volLo.z(), p.z())};
        volHi = {std::max(volHi.x(), p.x()), std::max(volHi.y(), p.y()), std::max(volHi.z(), p.z())};
    }
    EXPECT_NEAR(volLo.x(), plcLo.x(), 1.0e-9);
    EXPECT_NEAR(volLo.y(), plcLo.y(), 1.0e-9);
    EXPECT_NEAR(volLo.z(), plcLo.z(), 1.0e-9);
    EXPECT_NEAR(volHi.x(), plcHi.x(), 1.0e-9);
    EXPECT_NEAR(volHi.y(), plcHi.y(), 1.0e-9);
    EXPECT_NEAR(volHi.z(), plcHi.z(), 1.0e-9);
}

TEST(makeSolidSphere, ClosedTetFillContainsCentre) {
    rcs::PartDefinition ball = rcs::makeSolidSphere(0.10);
    rcs::applyResolutionSizing(ball, 1.0);
    ASSERT_TRUE(rcs::ensurePartMesh(ball));
    EXPECT_FALSE(ball.volume.empty());
    EXPECT_TRUE(ball.surface.empty());
    EXPECT_TRUE(rcs::closedEnough(ball.source));
    bool centreHit = false;
    for (const auto& tet : ball.volume.tets) {
        const auto& a = ball.volume.vertices[static_cast<std::size_t>(tet[0])];
        const auto& b = ball.volume.vertices[static_cast<std::size_t>(tet[1])];
        const auto& c = ball.volume.vertices[static_cast<std::size_t>(tet[2])];
        const auto& d = ball.volume.vertices[static_cast<std::size_t>(tet[3])];
        if (rcs::tetContainsPoint(a, b, c, d, {0, 0, 0})) {
            centreHit = true;
            break;
        }
    }
    EXPECT_TRUE(centreHit);
}

TEST(tetContainsPoint, CentreInsideFarOutside) {
    const chrono::ChVector3d a{1, 1, 1};
    const chrono::ChVector3d b{1, -1, -1};
    const chrono::ChVector3d c{-1, 1, -1};
    const chrono::ChVector3d d{-1, -1, 1};
    EXPECT_TRUE(rcs::tetContainsPoint(a, b, c, d, {0, 0, 0}));
    EXPECT_FALSE(rcs::tetContainsPoint(a, b, c, d, {10, 0, 0}));
}

rcs::PartDefinition thinClosedPlate(const rcs::Material& material) {
    rcs::PartDefinition plate;
    plate.name = material.name;
    plate.material = material;
    plate.thickness = 0.008;
    const double L = 0.20;
    const double t = 0.008;
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

TEST(ensurePartMesh, ThinClosedSheetsStayShells) {
    const std::vector<rcs::Material> materials = {rcs::rubber(), rcs::aluminum6061(), rcs::glass(), rcs::wood()};
    for (const rcs::Material& material : materials) {
        rcs::PartDefinition plate = thinClosedPlate(material);
        rcs::applyResolutionSizing(plate, 1.0);
        ASSERT_TRUE(rcs::ensurePartMesh(plate)) << material.name;
        EXPECT_TRUE(plate.volume.empty()) << material.name;
        EXPECT_FALSE(plate.surface.empty()) << material.name;
        rcs::Part part;
        part.build(plate);
        EXPECT_EQ(part.tetCount(), 0) << material.name;
        EXPECT_GT(part.triangles().size(), 0u) << material.name;
    }
}

TEST(ensurePartMesh, BulkyCubeKeepsComparableTetEdges) {
    rcs::PartDefinition cube = rcs::makeSolidCube(0.08);
    rcs::applyResolutionSizing(cube, 1.0);
    ASSERT_TRUE(rcs::ensurePartMesh(cube));
    ASSERT_FALSE(cube.volume.empty());
    double longest = 0.0;
    double shortest = 1.0e9;
    const int edges[6][2] = {{0, 1}, {0, 2}, {0, 3}, {1, 2}, {1, 3}, {2, 3}};
    for (const auto& tet : cube.volume.tets) {
        for (const auto& edge : edges) {
            const double len = (cube.volume.vertices[static_cast<std::size_t>(tet[edge[0]])] -
                                cube.volume.vertices[static_cast<std::size_t>(tet[edge[1]])])
                                   .Length();
            longest = std::max(longest, len);
            shortest = std::min(shortest, len);
        }
    }
    ASSERT_GT(shortest, 1.0e-9);
    EXPECT_LT(longest / shortest, 12.0);
}

TEST(seedWeldBeads, MarksInterfaceFacets) {
    rcs::AssemblyDefinition assembly;
    assembly.parts.push_back(rcs::makeSolidCube(0.08, {-0.04, 0, 0}));
    assembly.parts.push_back(rcs::makeSolidCube(0.08, {0.04, 0, 0}));
    rcs::WeldDefinition weld;
    weld.partA = 0;
    weld.partB = 1;
    weld.peelStrength = 1.0e8;
    weld.bead = {{0, -0.04, -0.04}, {0, 0.04, -0.04}};
    assembly.welds.push_back(weld);
    rcs::seedWeldBeads(assembly);
    bool marked = false;
    for (int marker : assembly.parts[0].source.markers)
        marked = marked || marker == 2;
    EXPECT_TRUE(marked);
}

rcs::PartDefinition thinCurvedShell() {
    rcs::PartDefinition part;
    part.name = "Curved rubber";
    part.material = rcs::rubber();
    part.thickness = 0.008;
    const double L = 0.20;
    const double t = 0.008;
    auto bow = [](double x) {
        const double u = (x - 0.10) / 0.10;
        return 0.004 * (1.0 - u * u);
    };
    part.source.vertices = {
        {0, bow(0), 0}, {L, bow(L), 0}, {L, t + bow(L), 0}, {0, t + bow(0), 0},
        {0, bow(0), L}, {L, bow(L), L}, {L, t + bow(L), L}, {0, t + bow(0), L},
    };
    const int tris[][3] = {
        {0, 2, 1}, {0, 3, 2}, {4, 5, 6}, {4, 6, 7}, {0, 1, 5}, {0, 5, 4},
        {3, 6, 2}, {3, 7, 6}, {0, 4, 7}, {0, 7, 3}, {1, 2, 6}, {1, 6, 5},
    };
    for (const auto& tri : tris)
        part.source.triangles.push_back({tri[0], tri[1], tri[2]});
    part.sizing.minElement = 0.002;
    part.sizing.maxElement = 0.025;
    part.sizing.tetBudget = 400;
    return part;
}

TEST(ensurePartMesh, CurvedThinRubberStaysShellAndBends) {
    rcs::PartDefinition curved = thinCurvedShell();
    rcs::applyResolutionSizing(curved, 1.0);
    ASSERT_TRUE(rcs::ensurePartMesh(curved));
    EXPECT_TRUE(curved.volume.empty());
    EXPECT_FALSE(curved.surface.empty());

    rcs::PartDefinition sheet;
    sheet.name = "Rubber sheet";
    sheet.material = rcs::rubber();
    sheet.thickness = 0.008;
    sheet.clampBoundary = true;
    sheet.patch.origin = {-0.1, 0.05, -0.1};
    sheet.patch.uAxis = {1.0, 0.0, 0.0};
    sheet.patch.vAxis = {0.0, 0.0, 1.0};
    sheet.patch.uSize = 0.2;
    sheet.patch.vSize = 0.2;
    rcs::applyResolutionSizing(sheet, 0.25);
    ASSERT_TRUE(rcs::ensurePartMesh(sheet));
    EXPECT_TRUE(sheet.volume.empty());

    rcs::AssemblyDefinition assembly;
    assembly.parts.push_back(sheet);
    rcs::ImpactLoad load;
    load.partIndex = 0;
    load.point = {0.0, 0.05, 0.0};
    load.direction = {0.0, -1.0, 0.0};
    load.initialForce = 5.0;
    load.finalForce = 5.0;
    load.duration = 0.05;
    load.radius = 0.03;
    rcs::Scenario scene;
    scene.loads.push_back(load);
    scene.settings.duration = 0.05;
    scene.settings.timestep = 2.0e-3;
    scene.settings.captureRate = 20.0;
    scene.settings.meshResolution = 0.25;
    scene.environment.gravity = {0.0, 0.0, 0.0};
    scene.environment.floorEnabled = false;
    scene.environment.aerodynamicsEnabled = false;
    scene.initial.position = {0.0, 0.0, 0.0};
    scene.initial.orientation = {0.0, 0.0, 0.0};
    rcs::Simulator sim;
    const rcs::Recording rec = sim.run(assembly, scene, nullptr);
    EXPECT_FALSE(rec.diverged);
    ASSERT_FALSE(rec.frames.empty());
    double lowest = 1.0e9;
    const auto& part = sim.assembly().parts().front();
    for (std::size_t n = 0; n < part.nodeCount(); ++n)
        lowest = std::min(lowest, part.nodePosition(n).y());
    EXPECT_LT(lowest, 0.05 - 1.0e-3);
}
