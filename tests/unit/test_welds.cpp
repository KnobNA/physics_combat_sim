#include "physics_core/physics_units.hpp"
#include "physics_core/cad.hpp"
#include "physics_core/mesh_split.hpp"
#include "physics_core/simulator.hpp"

#include "chrono/core/ChVector3.h"

#include <gtest/gtest.h>

#include <algorithm>
#include <chrono>
#include <cmath>
#include <limits>
#include <string>
#include <vector>

TEST(closedBeadArea, UnitSquare) {
    const std::vector<chrono::ChVector3d> square{{0, 0, 0}, {1, 0, 0}, {1, 1, 0}, {0, 1, 0}};
    EXPECT_NEAR(rcs::closedBeadArea(square), 1.0, 1.0e-12);
}

TEST(beadPerimeter, UnitSquare) {
    const std::vector<chrono::ChVector3d> square{{0, 0, 0}, {1, 0, 0}, {1, 1, 0}, {0, 1, 0}};
    EXPECT_NEAR(rcs::beadPerimeter(square), 4.0, 1.0e-12);
}

TEST(distanceToBead, PointOnEdgeIsZero) {
    const std::vector<chrono::ChVector3d> square{{0, 0, 0}, {1, 0, 0}, {1, 1, 0}, {0, 1, 0}};
    EXPECT_NEAR(rcs::distanceToBead(square, {0.5, 0.0, 0.0}), 0.0, 1.0e-12);
}

TEST(distanceToBead, CentreOfUnitSquare) {
    const std::vector<chrono::ChVector3d> square{{0, 0, 0}, {1, 0, 0}, {1, 1, 0}, {0, 1, 0}};
    EXPECT_NEAR(rcs::distanceToBead(square, {0.5, 0.5, 0.0}), 0.5, 1.0e-12);
}

TEST(closestOnBead, ProjectsOntoEdge) {
    const std::vector<chrono::ChVector3d> square{{0, 0, 0}, {1, 0, 0}, {1, 1, 0}, {0, 1, 0}};
    const chrono::ChVector3d q = rcs::closestOnBead(square, {0.5, -0.25, 0.0});
    EXPECT_NEAR(q.x(), 0.5, 1.0e-12);
    EXPECT_NEAR(q.y(), 0.0, 1.0e-12);
    EXPECT_NEAR(q.z(), 0.0, 1.0e-12);
}

TEST(unite, AlreadyUnitedIsFalse) {
    std::vector<int> set{0, 1, 2};
    EXPECT_TRUE(rcs::unite(set, 0, 1));
    EXPECT_FALSE(rcs::unite(set, 0, 1));
    EXPECT_EQ(rcs::findRoot(set, 0), rcs::findRoot(set, 1));
}

TEST(findRoot, PathCompressedToRepresentative) {
    std::vector<int> set{0, 1, 2, 3};
    ASSERT_TRUE(rcs::unite(set, 0, 1));
    ASSERT_TRUE(rcs::unite(set, 1, 2));
    EXPECT_EQ(rcs::findRoot(set, 2), rcs::findRoot(set, 0));
    EXPECT_EQ(rcs::findRoot(set, 3), 3);
}

namespace {

struct Aabb {
    chrono::ChVector3d lo;
    chrono::ChVector3d hi;
};

Aabb partAabb(const rcs::PartDefinition& part) {
    Aabb box{part.source.vertices.front(), part.source.vertices.front()};
    for (const auto& p : part.source.vertices) {
        box.lo = {std::min(box.lo.x(), p.x()), std::min(box.lo.y(), p.y()), std::min(box.lo.z(), p.z())};
        box.hi = {std::max(box.hi.x(), p.x()), std::max(box.hi.y(), p.y()), std::max(box.hi.z(), p.z())};
    }
    return box;
}

double aabbOverlap(const Aabb& a, const Aabb& b, int axis) {
    const double lo = axis == 0 ? std::max(a.lo.x(), b.lo.x()) : axis == 1 ? std::max(a.lo.y(), b.lo.y())
                                                                          : std::max(a.lo.z(), b.lo.z());
    const double hi = axis == 0 ? std::min(a.hi.x(), b.hi.x()) : axis == 1 ? std::min(a.hi.y(), b.hi.y())
                                                                          : std::min(a.hi.z(), b.hi.z());
    return hi - lo;
}

bool pointOnAabb(const Aabb& box, const chrono::ChVector3d& p, double tol) {
    return p.x() >= box.lo.x() - tol && p.x() <= box.hi.x() + tol && p.y() >= box.lo.y() - tol &&
           p.y() <= box.hi.y() + tol && p.z() >= box.lo.z() - tol && p.z() <= box.hi.z() + tol;
}

}  // namespace

TEST(makeTestCube, SixPanelsKissWithoutOverlapping) {
    const rcs::AssemblyDefinition cube = rcs::makeTestCube();
    ASSERT_EQ(cube.parts.size(), 6u);
    ASSERT_EQ(cube.welds.size(), 12u);
    std::vector<Aabb> boxes;
    for (const auto& part : cube.parts) {
        ASSERT_EQ(part.source.vertices.size(), 8u);
        boxes.push_back(partAabb(part));
    }
    for (std::size_t i = 0; i < boxes.size(); ++i) {
        for (std::size_t j = i + 1; j < boxes.size(); ++j) {
            const double ox = aabbOverlap(boxes[i], boxes[j], 0);
            const double oy = aabbOverlap(boxes[i], boxes[j], 1);
            const double oz = aabbOverlap(boxes[i], boxes[j], 2);
            EXPECT_FALSE(ox > 1.0e-6 && oy > 1.0e-6 && oz > 1.0e-6)
                << cube.parts[i].name << " overlaps " << cube.parts[j].name;
        }
    }
    const chrono::ChVector3d cavity(0, 0, 0);
    for (std::size_t i = 0; i < boxes.size(); ++i)
        EXPECT_FALSE(pointOnAabb(boxes[i], cavity, -1.0e-6)) << cube.parts[i].name << " fills the cavity";
}

TEST(makeTestCube, TwelveWeldsLieOnKissingFaces) {
    const rcs::TestCubeOptions options;
    const rcs::AssemblyDefinition cube = rcs::makeTestCube(options);
    ASSERT_EQ(cube.welds.size(), 12u);
    std::vector<Aabb> boxes;
    for (const auto& part : cube.parts)
        boxes.push_back(partAabb(part));
    const double t = options.thickness;
    const double inner = options.outerSize - 2.0 * t;
    int backSeams = 0;
    for (const auto& weld : cube.welds) {
        ASSERT_LT(weld.partA, boxes.size());
        ASSERT_LT(weld.partB, boxes.size());
        ASSERT_EQ(weld.bead.size(), 4u) << weld.label;
        const double area = rcs::closedBeadArea(weld.bead);
        const double stripInner = t * inner;
        const double stripOuter = t * options.outerSize;
        EXPECT_TRUE(std::abs(area - stripInner) < 1.0e-9 || std::abs(area - stripOuter) < 1.0e-9)
            << weld.label << " area " << area;
        for (const auto& p : weld.bead) {
            EXPECT_TRUE(pointOnAabb(boxes[weld.partA], p, 1.0e-9)) << weld.label << " misses " << cube.parts[weld.partA].name;
            EXPECT_TRUE(pointOnAabb(boxes[weld.partB], p, 1.0e-9)) << weld.label << " misses " << cube.parts[weld.partB].name;
        }
        if (weld.label.find("Back - Right") != std::string::npos ||
            weld.label.find("Back - Left") != std::string::npos) {
            ++backSeams;
            const double z = 0.25 * (weld.bead[0].z() + weld.bead[1].z() + weld.bead[2].z() + weld.bead[3].z());
            EXPECT_NEAR(z, -(0.5 * options.outerSize - t), 1.0e-9) << weld.label;
            EXPECT_GT(std::abs(z + 0.5 * options.outerSize), 0.5 * t) << weld.label << " sat on the outer back face";
        }
    }
    EXPECT_EQ(backSeams, 2);
}

// Every seam needs nodes on the bead on both sides, whether the policy makes
// the panels shells or solids (issue #1). Checked on the built parts, so a
// Reissner patch counts as much as a tet fill.
TEST(makeTestCube, MeshedSeamsKeepKissingNodes) {
    rcs::AssemblyDefinition box = rcs::makeTestCube();
    for (auto& part : box.parts)
        rcs::applyResolutionSizing(part, 1.0);
    rcs::seedWeldBeads(box);
    const double inner = 0.5 * (0.40 - 0.008) - 0.004;
    std::vector<rcs::Part> built(box.parts.size());
    int elementsMin = std::numeric_limits<int>::max();
    int elementsMax = 0;
    for (std::size_t i = 0; i < box.parts.size(); ++i) {
        auto& part = box.parts[i];
        ASSERT_TRUE(rcs::ensurePartMesh(part)) << part.name;
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
        built[i].build(part);
        elementsMin = std::min(elementsMin, built[i].totalElementCount());
        elementsMax = std::max(elementsMax, built[i].totalElementCount());
    }
    EXPECT_GE(elementsMin, 16);
    EXPECT_LT(elementsMax, elementsMin * 25);

    const double capture = 1.5 * 0.008;
    for (const auto& weld : box.welds) {
        auto hits = [&](const rcs::Part& part) {
            int count = 0;
            for (std::size_t n = 0; n < part.nodeCount(); ++n)
                count += rcs::distanceToBead(weld.bead, part.nodePosition(n)) <= capture ? 1 : 0;
            return count;
        };
        EXPECT_GE(hits(built[weld.partA]), 4) << weld.label << " missed " << box.parts[weld.partA].name;
        EXPECT_GE(hits(built[weld.partB]), 4) << weld.label << " missed " << box.parts[weld.partB].name;
    }
}

TEST(makeTestCube, WeldSnapKeepsOuterSkin) {
    rcs::AssemblyDefinition cube = rcs::makeTestCube();
    for (auto& part : cube.parts)
        rcs::applyResolutionSizing(part, 1.0);
    rcs::seedWeldBeads(cube);
    struct Box {
        chrono::ChVector3d lo{0, 0, 0};
        chrono::ChVector3d hi{0, 0, 0};
    };
    std::vector<Box> before;
    for (auto& part : cube.parts) {
        ASSERT_TRUE(rcs::ensurePartMesh(part)) << part.name;
        Box box{part.volume.vertices.front(), part.volume.vertices.front()};
        for (const auto& p : part.volume.vertices) {
            box.lo = {std::min(box.lo.x(), p.x()), std::min(box.lo.y(), p.y()), std::min(box.lo.z(), p.z())};
            box.hi = {std::max(box.hi.x(), p.x()), std::max(box.hi.y(), p.y()), std::max(box.hi.z(), p.z())};
        }
        before.push_back(box);
    }

    rcs::Scenario scene;
    scene.environment.gravity = {0.0, 0.0, 0.0};
    scene.environment.floorEnabled = false;
    scene.environment.aerodynamicsEnabled = false;
    scene.initial.position = {0.0, 0.0, 0.0};
    scene.initial.orientation = {0.0, 0.0, 0.0};
    scene.settings.meshResolution = 1.0;
    rcs::Simulator sim;
    sim.reset(cube, scene, false);

    const double maxMove = 0.0015;
    ASSERT_EQ(sim.assembly().parts().size(), before.size());
    for (std::size_t i = 0; i < before.size(); ++i) {
        const auto& part = sim.assembly().parts()[i];
        ASSERT_GT(part.nodeCount(), 0u) << cube.parts[i].name;
        chrono::ChVector3d lo = part.nodePosition(0);
        chrono::ChVector3d hi = lo;
        for (std::size_t n = 1; n < part.nodeCount(); ++n) {
            const chrono::ChVector3d p = part.nodePosition(n);
            lo = {std::min(lo.x(), p.x()), std::min(lo.y(), p.y()), std::min(lo.z(), p.z())};
            hi = {std::max(hi.x(), p.x()), std::max(hi.y(), p.y()), std::max(hi.z(), p.z())};
        }
        EXPECT_NEAR(lo.x(), before[i].lo.x(), maxMove) << cube.parts[i].name;
        EXPECT_NEAR(lo.y(), before[i].lo.y(), maxMove) << cube.parts[i].name;
        EXPECT_NEAR(lo.z(), before[i].lo.z(), maxMove) << cube.parts[i].name;
        EXPECT_NEAR(hi.x(), before[i].hi.x(), maxMove) << cube.parts[i].name;
        EXPECT_NEAR(hi.y(), before[i].hi.y(), maxMove) << cube.parts[i].name;
        EXPECT_NEAR(hi.z(), before[i].hi.z(), maxMove) << cube.parts[i].name;
    }
    int welded = 0;
    for (const auto& joint : sim.assembly().joints()) {
        if (joint.kind != rcs::JointKind::Weld)
            continue;
        EXPECT_GE(joint.links.size(), 2u) << joint.label;
        ++welded;
    }
    EXPECT_EQ(welded, 12);
}

// Issue #2: the rubber box may be shells or solids under the compatibility
// policy; either way the welds must hold and the box must land on the floor.
TEST(makeTestCube, RubberBoxKeepsWeldsAndFloor) {
    rcs::AssemblyDefinition box = rcs::makeTestCube();
    for (auto& part : box.parts) {
        part.material = rcs::rubber();
        part.sizing.cached = false;
        part.volume = {};
        part.surface = {};
    }
    rcs::Scenario scene;
    scene.environment.gravity = {0.0, -9.81, 0.0};
    scene.environment.floorEnabled = true;
    scene.environment.aerodynamicsEnabled = false;
    scene.settings.duration = 0.08;
    scene.settings.timestep = 2.0e-3;
    scene.settings.captureRate = 20.0;
    scene.settings.meshResolution = 0.25;
    scene.initial.position = {0.0, 0.22, 0.0};
    scene.initial.orientation = {0.0, 0.0, 0.0};
    rcs::Simulator sim;
    const auto t0 = std::chrono::steady_clock::now();
    const rcs::Recording rec = sim.run(box, scene, nullptr);
    const double wall = std::chrono::duration<double>(std::chrono::steady_clock::now() - t0).count();
    EXPECT_FALSE(rec.diverged);
    EXPECT_LT(wall, 20.0) << "s/step crawl commits=" << sim.hhtCommitCount();
    int links = 0;
    int broken = 0;
    for (const auto& joint : sim.assembly().joints()) {
        if (joint.kind != rcs::JointKind::Weld)
            continue;
        links += static_cast<int>(joint.links.size());
        broken += joint.brokenCount();
    }
    EXPECT_GT(links, 12);
    EXPECT_LT(broken, links);
    double lowest = 1.0e9;
    for (const auto& part : sim.assembly().parts()) {
        for (std::size_t n = 0; n < part.nodeCount(); ++n)
            lowest = std::min(lowest, part.nodePosition(n).y());
    }
    EXPECT_GT(lowest, 0.0) << "fell through the floor, lowest=" << lowest;
    EXPECT_LT(lowest, 0.02) << "never reached the floor, lowest=" << lowest;
}

TEST(makeTestCube, RestingHhtStaysInteractive) {
    rcs::Scenario scene;
    scene.environment.gravity = {0.0, -9.81, 0.0};
    scene.environment.floorEnabled = true;
    scene.environment.aerodynamicsEnabled = false;
    scene.settings.duration = 0.002;
    scene.settings.timestep = 5.0e-4;
    scene.settings.captureRate = 20.0;
    scene.settings.meshResolution = 1.0;
    scene.initial.position = {0.0, 0.201, 0.0};
    scene.initial.orientation = {0.0, 0.0, 0.0};
    rcs::Simulator sim;
    const auto t0 = std::chrono::steady_clock::now();
    const rcs::Recording rec = sim.run(rcs::makeTestCube(), scene, nullptr);
    const double wall = std::chrono::duration<double>(std::chrono::steady_clock::now() - t0).count();
    const int commits = std::max(sim.hhtCommitCount(), 1);
    EXPECT_FALSE(rec.diverged);
    EXPECT_GT(sim.hhtCommitCount(), 0);
    EXPECT_LT(wall / static_cast<double>(commits), 2.0) << "s/step=" << wall / static_cast<double>(commits)
                                                        << " rejects=" << sim.newtonRejectCount();
}
