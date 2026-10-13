#include "physics_core/physics_units.hpp"

#include <cmath>
#include <gtest/gtest.h>

TEST(continuumGains, ScalesWithModulusAndArea) {
    const rcs::ContactGains gains = rcs::continuumGains(69.0e9, 0.01, 1.0);
    EXPECT_NEAR(gains.stiffness, rcs::kContactPenaltyScale * 69.0e9 * 0.1, 1.0e3);
    EXPECT_NEAR(gains.damping, rcs::kContactDampingFraction * 2.0 * std::sqrt(gains.stiffness), 1.0e-6);
}

TEST(stableGains, KeepsMaterialStiffness) {
    const rcs::ContactGains gains = rcs::stableGains(1.0e12, 1.0e12, 1.0, 1.0e-3);
    EXPECT_NEAR(gains.stiffness, 1.0e12, 1.0);
    EXPECT_NEAR(gains.damping, rcs::kContactDampingFraction * 2.0 * std::sqrt(1.0e12), 1.0);
}

TEST(contactDashpot, StabilizationVanishesAtRest) {
    const double k = 1.0e6;
    const double m = 1.0;
    const double crit = 2.0 * std::sqrt(k * m);
    EXPECT_NEAR(rcs::contactDashpot(k, m, 0.0, 1.0e-3), rcs::kContactDampingFraction * crit, 1.0e-9);
    EXPECT_NEAR(rcs::contactDashpot(k, m, 0.5, 1.0e-3), rcs::kContactDampingFraction * crit, 1.0e-9);
    EXPECT_NEAR(rcs::contactDashpot(k, m, -1.0, 1.0e-3),
                (rcs::kContactDampingFraction + rcs::kContactStabFraction) * crit, 1.0e-9);
}

TEST(stableGains, ZeroMassKeepsK) {
    EXPECT_DOUBLE_EQ(rcs::stableGains(1.0e6, 1.0e3, 0.0, 1.0e-3).stiffness, 1.0e6);
    EXPECT_DOUBLE_EQ(rcs::stableGains(1.0e6, 1.0e3, 0.0, 1.0e-3).damping, 0.0);
    EXPECT_DOUBLE_EQ(rcs::stableGains(1.0e6, 1.0e3, 1.0, 0.0).stiffness, 1.0e6);
}

TEST(timeToClosingPlane, LinearWhenNoAcceleration) {
    EXPECT_NEAR(rcs::timeToClosingPlane(0.1, 2.0, 0.0), 0.05, 1.0e-12);
}

TEST(timeToClosingPlane, DropFromRestIsSqrtTwoHOverG) {
    const double g = 9.81;
    EXPECT_NEAR(rcs::timeToClosingPlane(0.08, 0.0, g), std::sqrt(2.0 * 0.08 / g), 1.0e-12);
}

TEST(timeToClosingPlane, GravityShortensFallWithSpeed) {
    const double g = 9.81;
    const double linear = 0.08 / 2.0;
    const double t = rcs::timeToClosingPlane(0.08, 2.0, g);
    EXPECT_LT(t, linear);
    EXPECT_GT(t, 0.0);
    EXPECT_NEAR(0.08 - 2.0 * t - 0.5 * g * t * t, 0.0, 1.0e-12);
}

TEST(timeToClosingPlane, ClosedOrReceding) {
    EXPECT_DOUBLE_EQ(rcs::timeToClosingPlane(0.0, 2.0, 9.81), 0.0);
    EXPECT_DOUBLE_EQ(rcs::timeToClosingPlane(-1.0e-4, 2.0, 9.81), 0.0);
    EXPECT_DOUBLE_EQ(rcs::timeToClosingPlane(0.1, -1.0, 0.0), rcs::kNoContactTime);
}
