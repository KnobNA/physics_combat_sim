#include "physics_core/scenario.hpp"

#include <gtest/gtest.h>

#include <cmath>

TEST(loadRampBlend, LinearMidpoint) {
    EXPECT_DOUBLE_EQ(rcs::loadRampBlend(rcs::LoadRamp::Linear, 2, 0.5), 0.5);
}

TEST(loadRampBlend, ExponentialEndpoints) {
    EXPECT_DOUBLE_EQ(rcs::loadRampBlend(rcs::LoadRamp::Exponential, 2, 0.0), 0.0);
    EXPECT_DOUBLE_EQ(rcs::loadRampBlend(rcs::LoadRamp::Exponential, 2, 1.0), 1.0);
}

TEST(loadRampBlend, LogarithmicEndpoints) {
    EXPECT_DOUBLE_EQ(rcs::loadRampBlend(rcs::LoadRamp::Logarithmic, 2, 0.0), 0.0);
    EXPECT_DOUBLE_EQ(rcs::loadRampBlend(rcs::LoadRamp::Logarithmic, 2, 1.0), 1.0);
}

TEST(loadRampBlend, MidRampLogGreaterThanLinearGreaterThanExp) {
    const double linear = rcs::loadRampBlend(rcs::LoadRamp::Linear, 2, 0.5);
    const double logRamp = rcs::loadRampBlend(rcs::LoadRamp::Logarithmic, 2, 0.5);
    const double expRamp = rcs::loadRampBlend(rcs::LoadRamp::Exponential, 2, 0.5);
    EXPECT_GT(logRamp, linear);
    EXPECT_GT(linear, expRamp);
}

TEST(loadRampBlend, PolynomialDegreeTwo) {
    EXPECT_DOUBLE_EQ(rcs::loadRampBlend(rcs::LoadRamp::Polynomial, 2, 0.5), 0.25);
}

TEST(loadRampMaxSlope, KnownShapes) {
    EXPECT_DOUBLE_EQ(rcs::loadRampMaxSlope(rcs::LoadRamp::Linear, 2), 1.0);
    EXPECT_DOUBLE_EQ(rcs::loadRampMaxSlope(rcs::LoadRamp::Logarithmic, 2), std::exp(1.0) - 1.0);
    EXPECT_DOUBLE_EQ(rcs::loadRampMaxSlope(rcs::LoadRamp::Exponential, 2),
                     std::exp(1.0) / (std::exp(1.0) - 1.0));
    EXPECT_DOUBLE_EQ(rcs::loadRampMaxSlope(rcs::LoadRamp::Polynomial, 3), 3.0);
}

TEST(loadMagnitude, ZeroOutsideWindow) {
    rcs::ImpactLoad load;
    load.startTime = 1.0;
    load.duration = 2.0;
    load.rampTime = 0.5;
    load.initialForce = 0.0;
    load.finalForce = 100.0;
    load.ramp = rcs::LoadRamp::Linear;
    EXPECT_DOUBLE_EQ(rcs::loadMagnitude(load, 0.5), 0.0);
    EXPECT_DOUBLE_EQ(rcs::loadMagnitude(load, 3.0), 0.0);
}

TEST(loadMagnitude, LinearRampThenHold) {
    rcs::ImpactLoad load;
    load.startTime = 1.0;
    load.duration = 2.0;
    load.rampTime = 0.5;
    load.initialForce = 0.0;
    load.finalForce = 100.0;
    load.ramp = rcs::LoadRamp::Linear;
    EXPECT_DOUBLE_EQ(rcs::loadMagnitude(load, 1.0), 0.0);
    EXPECT_DOUBLE_EQ(rcs::loadMagnitude(load, 1.25), 50.0);
    EXPECT_DOUBLE_EQ(rcs::loadMagnitude(load, 1.5), 100.0);
    EXPECT_DOUBLE_EQ(rcs::loadMagnitude(load, 2.5), 100.0);
}

TEST(loadMagnitude, DisabledOrZeroDuration) {
    rcs::ImpactLoad load;
    load.enabled = false;
    load.finalForce = 100.0;
    EXPECT_DOUBLE_EQ(rcs::loadMagnitude(load, 0.0), 0.0);
    load.enabled = true;
    load.duration = 0.0;
    EXPECT_DOUBLE_EQ(rcs::loadMagnitude(load, 0.0), 0.0);
}

TEST(loadMaxStepChange, StepAndRamp) {
    rcs::ImpactLoad step;
    step.startTime = 0.0;
    step.duration = 0.01;
    step.rampTime = 0.0;
    step.finalForce = 5000.0;
    EXPECT_DOUBLE_EQ(rcs::loadMaxStepChange(step, 5.0e-4), 5000.0);

    rcs::ImpactLoad ramp;
    ramp.startTime = 0.0;
    ramp.duration = 1.0;
    ramp.rampTime = 0.5;
    ramp.initialForce = 0.0;
    ramp.finalForce = 100.0;
    ramp.ramp = rcs::LoadRamp::Linear;
    EXPECT_DOUBLE_EQ(rcs::loadMaxStepChange(ramp, 0.1), 20.0);
}

TEST(gridStep, ArenaSpans) {
    EXPECT_DOUBLE_EQ(rcs::Environment::gridStep(10.0), 1.0);
    EXPECT_DOUBLE_EQ(rcs::Environment::gridStep(50.0), 2.0);
    EXPECT_DOUBLE_EQ(rcs::Environment::gridStep(100.0), 5.0);
}
