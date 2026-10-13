#include "physics_core/material.hpp"
#include "physics_core/part.hpp"
#include "physics_core/physics_units.hpp"

#include <gtest/gtest.h>

TEST(brittleShouldOpen, TensionBelowRankineIsFalse) {
    EXPECT_FALSE(rcs::brittleShouldOpen(rcs::glass(), 293.15, 4.4e7, 1.0, 0.0, 0.01));
}

TEST(brittleShouldOpen, TensionAtRankineOpens) {
    EXPECT_TRUE(rcs::brittleShouldOpen(rcs::glass(), 293.15, 5.0e7, 10.0, 0.0, 0.01));
}

TEST(brittleShouldOpen, MixedModeStillOpens) {
    EXPECT_TRUE(rcs::brittleShouldOpen(rcs::glass(), 293.15, 5.0e7, 10.0, 0.0, 0.01));
}

TEST(brittleShouldOpen, TemperedAddsResidualDrive) {
    EXPECT_EQ(rcs::temperedGlass().breakFabric, rcs::BreakFabric::Tempered);
    EXPECT_TRUE(rcs::brittleShouldOpen(rcs::temperedGlass(), 293.15, 1.3e8, 10.0, 0.0, 0.01));
    EXPECT_TRUE(rcs::brittleShouldOpen(rcs::temperedGlass(), 293.15, 8.5e7, 10.0, 0.0, 0.01));
    EXPECT_FALSE(rcs::brittleShouldOpen(rcs::temperedGlass(), 293.15, 7.0e7, 10.0, 0.0, 0.01));
}

TEST(brittleShouldOpen, GrainOpensAcrossFiberFirst) {
    const rcs::Material wood = rcs::wood();
    EXPECT_EQ(wood.breakFabric, rcs::BreakFabric::Grain);
    EXPECT_NEAR(rcs::grainOpeningScale(wood, {0, 1, 0}), 1.0, 1.0e-12);
    EXPECT_NEAR(rcs::grainOpeningScale(wood, {1, 0, 0}), 3.0, 1.0e-12);
    EXPECT_TRUE(rcs::brittleShouldOpen(wood, 293.15, 5.5e6, 10.0, 0.0, 0.01, 0.0, 0.0, {0, 1, 0}));
    EXPECT_FALSE(rcs::brittleShouldOpen(wood, 293.15, 5.5e6, 10.0, 0.0, 0.01, 0.0, 0.0, {1, 0, 0}));
    EXPECT_TRUE(rcs::brittleShouldOpen(wood, 293.15, 1.6e7, 10.0, 0.0, 0.01, 0.0, 0.0, {1, 0, 0}));
}

TEST(brittleShouldOpen, DuctileNeverOpens) {
    EXPECT_FALSE(rcs::brittleShouldOpen(rcs::aluminum6061(), 293.15, 1.0e9, 1.0e6, 1.0, 0.01));
}

TEST(brittleShouldOpen, CornerHitStillOpens) {
    EXPECT_TRUE(rcs::brittleShouldOpen(rcs::glass(), 293.15, 5.0e7, 10.0, 0.0, 0.01, 4.0e8));
}

TEST(brittleShouldOpen, TensionIndependentOfCrushCap) {
    EXPECT_TRUE(rcs::brittleShouldOpen(rcs::glass(), 293.15, 5.0e7, 10.0, 0.0, 0.01, 8.0e8));
}

TEST(brittleShouldOpen, TensionDominatedStillOpens) {
    EXPECT_TRUE(rcs::brittleShouldOpen(rcs::glass(), 293.15, 5.0e7, 10.0, 0.0, 0.01, 1.0e7, 1.7e7));
}

TEST(brittleShouldOpen, HoopAtRankineOpens) {
    EXPECT_TRUE(rcs::brittleShouldOpen(rcs::glass(), 293.15, 5.0e7, 10.0, 0.0, 0.01, 4.0e7, 5.0e6));
}

TEST(brittleShouldOpen, RankineIgnoresMeanStress) {
    EXPECT_TRUE(rcs::brittleShouldOpen(rcs::glass(), 293.15, 5.0e7, 10.0, 0.0, 0.01, 4.0e7, -8.0e7));
}

TEST(brittleShouldOpen, EnergyBelowGriffithStaysClosed) {
    EXPECT_FALSE(rcs::brittleShouldOpen(rcs::glass(), 293.15, 5.0e7, 1.0e-6, 0.0, 0.01));
}

TEST(fabricCrackNormal, GrainPlaneContainsTheFiber) {
    const chrono::ChVector3d normal = rcs::fabricCrackNormal(rcs::wood(), {0.0, 1.0, 0.0});
    EXPECT_NEAR(normal.Dot({1.0, 0.0, 0.0}), 0.0, 1.0e-9);
    EXPECT_NEAR(normal.Length(), 1.0, 1.0e-12);
}

TEST(interfaceShouldOpen, PlyBoundaryOpensBeforeTheParent) {
    rcs::Material laminate = rcs::aluminum6061();
    laminate.interfaceStrength = 2.0e7;
    laminate.interfaceToughness = 0.0;
    laminate.stackNormal[0] = 0.0;
    laminate.stackNormal[1] = 1.0;
    laminate.stackNormal[2] = 0.0;
    laminate.plies.push_back({});
    laminate.plies.push_back({});
    laminate.plies[1].direction[0] = 0.0;
    laminate.plies[1].direction[1] = 0.0;
    laminate.plies[1].direction[2] = 1.0;
    EXPECT_TRUE(rcs::interfaceShouldOpen(laminate, 293.15, 3.0e7, {0, 1, 0}, 1.0, 0.01));
    EXPECT_FALSE(rcs::interfaceShouldOpen(laminate, 293.15, 3.0e7, {1, 0, 0}, 1.0, 0.01));
    EXPECT_FALSE(rcs::brittleShouldOpen(laminate, 293.15, 3.0e7, 10.0, 0.0, 0.01));
}

TEST(brittleShouldCrush, HalfHandbookDoesNotCrush) {
    EXPECT_FALSE(rcs::brittleShouldCrush(rcs::glass(), 293.15, 0.5 * rcs::glass().compressiveStrength));
}

TEST(brittleShouldCrush, OverHandbookCrushes) {
    EXPECT_TRUE(rcs::brittleShouldCrush(rcs::glass(), 293.15, 1.2 * rcs::glass().compressiveStrength));
}

TEST(brittleShouldCrush, DuctileNeverCrushes) {
    EXPECT_FALSE(rcs::brittleShouldCrush(rcs::aluminum6061(), 293.15, 1.0e10));
}

TEST(brittleShouldOpen, RankineDoesNotNeedToughness) {
    rcs::Material plate = rcs::glass();
    plate.fractureToughness = 0.0;
    EXPECT_TRUE(rcs::brittleShouldOpen(plate, 293.15, 5.0e7, 0.0, 0.001, 0.01));
}

TEST(cohesiveDamage, OpeningPastTheThroatExhaustsTheBead) {
    const double critical = rcs::cohesiveCriticalOpening(1.65e8, 0.0, 0.008);
    const rcs::CohesiveDamage intact = rcs::cohesiveDamage(0.0, 1.0e-5, 0.0, 0.0, 100.0, 100.0, critical);
    EXPECT_FALSE(intact.exhausted);
    const rcs::CohesiveDamage peeled = rcs::cohesiveDamage(critical, 0.0, 0.0, 0.0, 100.0, 100.0, critical);
    EXPECT_TRUE(peeled.exhausted);
    const rcs::CohesiveDamage sheared = rcs::cohesiveDamage(0.0, 0.0, 0.0, 100.0, 100.0, 100.0, critical);
    EXPECT_TRUE(sheared.exhausted);
    EXPECT_LT(rcs::cohesiveDamage(0.0, 0.0, 10.0, 10.0, 100.0, 100.0, critical).damage, 1.0);
}

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
    if (nodes.empty() || !nodes[0])
        return;
    chrono::ChVector3d p = nodes[1]->GetPos();
    p.x() *= scale;
    nodes[1]->SetPos(p);
    part.refreshTetTrials();
}

TEST(fractureCluster, SharedPlaneCutsBothTetsWhenEnergyAllows) {
    rcs::Material glass = rcs::glass();
    glass.fractureStrength = 1.0e6;
    glass.fractureToughness = 0.0;
    glass.minimumShardSize = 10.0;
    glass.shatterToDust = false;
    rcs::Part open;
    open.build(regularTet(glass));
    open.refreshTetTrials();
    const int before = open.tetCount();
    const std::size_t nodesBefore = open.nodeCount();
    stretchTet(open, 1.5);
    EXPECT_TRUE(open.tearFailed());
    EXPECT_GT(open.tetCount(), before);
    EXPECT_GT(open.nodeCount(), nodesBefore);
    EXPECT_GT(open.liveTetCount(), 0);

    rcs::Material blocked = glass;
    blocked.fractureToughness = 1.0e12;
    EXPECT_FALSE(rcs::energyAllowsOpening(blocked, 1.0, 1.0e-4));
    rcs::Part shut;
    shut.build(regularTet(blocked));
    shut.refreshTetTrials();
    const int shutBefore = shut.tetCount();
    stretchTet(shut, 1.5);
    shut.tearFailed();
    EXPECT_EQ(shut.tetCount(), shutBefore);
}

TEST(fractureCluster, GrainAndPlyOpenWithoutDeletingTheParent) {
    rcs::Material wood = rcs::wood();
    wood.fractureStrength = 1.0e6;
    wood.fractureToughness = 0.0;
    wood.minimumShardSize = 10.0;
    rcs::Part grained;
    grained.build(regularTet(wood));
    grained.refreshTetTrials();
    const std::size_t grainNodes = grained.nodeCount();
    stretchTet(grained, 1.5);
    EXPECT_TRUE(grained.tearFailed());
    EXPECT_GT(grained.nodeCount(), grainNodes);
    EXPECT_GT(grained.liveTetCount(), 0);

    rcs::Material laminate = rcs::aluminum6061();
    laminate.interfaceStrength = 1.0e6;
    laminate.interfaceToughness = 0.0;
    laminate.stackNormal[0] = 1.0;
    laminate.stackNormal[1] = 0.0;
    laminate.stackNormal[2] = 0.0;
    laminate.plies = {rcs::Ply{}, rcs::Ply{}};
    laminate.plies[0].thickness = 0.02;
    laminate.plies[1].thickness = 0.02;
    laminate.plies[1].direction[0] = 0.0;
    laminate.plies[1].direction[2] = 1.0;
    rcs::Part stack;
    stack.build(regularTet(laminate));
    stack.refreshTetTrials();
    EXPECT_GT(stack.tetCount(), 1);
    const int liveBefore = stack.liveTetCount();
    stretchTet(stack, 1.01);
    EXPECT_TRUE(stack.tearFailed());
    EXPECT_GE(stack.liveTetCount(), 1);
    EXPECT_GE(stack.liveTetCount(), liveBefore > 4 ? liveBefore - 4 : 1);
}
