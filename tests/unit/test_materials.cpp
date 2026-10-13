#include "physics_core/material.hpp"

#include <gtest/gtest.h>

TEST(mapThermalAxis, RoomTemperatureIdentity) {
    EXPECT_DOUBLE_EQ(rcs::mapThermalAxis(293.15, 420.0, 855.0, 420.0, 855.0), 293.15);
}

TEST(mapThermalAxis, PastMeltScales) {
    const double mapped = rcs::mapThermalAxis(955.0, 420.0, 855.0, 420.0, 855.0);
    EXPECT_DOUBLE_EQ(mapped, 955.0);
}

TEST(sampleThermalKnots, InterpolatesMidSpan) {
    const double temperatures[] = {0.0, 10.0};
    const double values[] = {1.0, 3.0};
    EXPECT_DOUBLE_EQ(rcs::sampleThermalKnots(5.0, temperatures, values, 2), 2.0);
}

TEST(sampleThermalKnots, ClampsBelowFirstKnot) {
    const double temperatures[] = {293.15, 373.0};
    const double values[] = {1.0, 0.97};
    EXPECT_DOUBLE_EQ(rcs::sampleThermalKnots(200.0, temperatures, values, 2), 1.0);
}

TEST(thermalStiffnessRetention, AluminumRoomNearUnity) {
    EXPECT_NEAR(rcs::thermalStiffnessRetention(rcs::aluminum6061(), 293.15), 1.0, 1.0e-12);
}

TEST(thermalStiffnessRetention, AluminumAboveSolidusNearFloor) {
    EXPECT_NEAR(rcs::thermalStiffnessRetention(rcs::aluminum6061(), 5000.0), rcs::kThermalFloor, 1.0e-6);
}

TEST(thermalStrengthRetention, AluminumRoomNearUnity) {
    EXPECT_NEAR(rcs::thermalStrengthRetention(rcs::aluminum6061(), 293.15), 1.0, 1.0e-12);
}

TEST(thermalDuctilityRetention, AluminumRoomIsOne) {
    EXPECT_DOUBLE_EQ(rcs::thermalDuctilityRetention(rcs::aluminum6061(), 293.15), 1.0);
}

TEST(kThermalFloor, Residual) {
    EXPECT_DOUBLE_EQ(rcs::kThermalFloor, 1.0e-4);
}

TEST(rankineTensionLimit, DuctileIsZero) {
    EXPECT_DOUBLE_EQ(rcs::rankineTensionLimit(rcs::aluminum6061()), 0.0);
}

TEST(rankineTensionLimit, GlassUsesFractureStrength) {
    EXPECT_DOUBLE_EQ(rcs::rankineTensionLimit(rcs::glass()), 4.5e7);
}

TEST(rankineCompressionLimit, GlassUsesExplicitCrush) {
    EXPECT_DOUBLE_EQ(rcs::rankineCompressionLimit(rcs::glass()), 8.0e8);
}

TEST(rankineCompressionLimit, DerivedFifteenTimesTension) {
    rcs::Material brittle = rcs::glass();
    brittle.compressiveStrength = 0.0;
    EXPECT_DOUBLE_EQ(rcs::rankineCompressionLimit(brittle), 15.0 * 4.5e7);
}

TEST(aluminum6061, DefaultNameAndDensity) {
    EXPECT_EQ(rcs::aluminum6061().name, "Aluminum 6061-T6");
    EXPECT_DOUBLE_EQ(rcs::aluminum6061().density, 2700.0);
}

TEST(mildSteel, MeltTemperature) {
    EXPECT_DOUBLE_EQ(rcs::mildSteel().meltTemperature, 1780.0);
}

TEST(ar500, YieldStress) {
    EXPECT_DOUBLE_EQ(rcs::ar500().yieldStress, 1.30e9);
}

TEST(glass, BrittleFracture) {
    EXPECT_EQ(rcs::glass().fracture, rcs::FractureMode::Brittle);
}

TEST(temperedGlass, ShatterToDust) {
    EXPECT_TRUE(rcs::temperedGlass().shatterToDust);
    EXPECT_DOUBLE_EQ(rcs::temperedGlass().fractureStrength, 1.20e8);
}

TEST(rubber, CompliantFracture) {
    EXPECT_EQ(rcs::rubber().fracture, rcs::FractureMode::Compliant);
}

TEST(materialPresets, LibraryOrder) {
    const auto presets = rcs::materialPresets();
    ASSERT_EQ(presets.size(), 7u);
    EXPECT_EQ(presets[0].name, rcs::aluminum6061().name);
    EXPECT_EQ(presets[1].name, rcs::mildSteel().name);
    EXPECT_EQ(presets[2].name, rcs::ar500().name);
    EXPECT_EQ(presets[3].name, rcs::glass().name);
    EXPECT_EQ(presets[4].name, rcs::temperedGlass().name);
    EXPECT_EQ(presets[5].name, rcs::rubber().name);
    EXPECT_EQ(presets[6].name, rcs::wood().name);
}

TEST(presetIndex, MatchesLibraryOrder) {
    EXPECT_EQ(rcs::presetIndex(rcs::aluminum6061()), 0);
    EXPECT_EQ(rcs::presetIndex(rcs::glass()), 3);
    EXPECT_EQ(rcs::presetIndex(rcs::rubber()), 5);
    EXPECT_EQ(rcs::presetIndex(rcs::wood()), 6);
}

TEST(wood, GrainFabricAlongX) {
    EXPECT_EQ(rcs::wood().breakFabric, rcs::BreakFabric::Grain);
    EXPECT_DOUBLE_EQ(rcs::wood().grainDirection[0], 1.0);
    EXPECT_EQ(rcs::wood().fracture, rcs::FractureMode::Brittle);
}

TEST(glass, IsotropicFabric) {
    EXPECT_EQ(rcs::glass().breakFabric, rcs::BreakFabric::Isotropic);
}

TEST(aluminum6061, IsotropicFabric) {
    EXPECT_EQ(rcs::aluminum6061().breakFabric, rcs::BreakFabric::Isotropic);
}
