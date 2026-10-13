#include "physics_core/recording.hpp"

#include <gtest/gtest.h>

#include <cmath>
#include <limits>

TEST(packHalf, MillimetreRoundTrip) {
    const float metres = 0.001f;
    EXPECT_NEAR(rcs::unpackHalf(rcs::packHalf(metres)), metres, 1.0e-3f);
}

TEST(unpackHalf, MillimetreScalePositions) {
    EXPECT_NEAR(rcs::unpackHalf(rcs::packHalf(0.123f)), 0.123f, 1.0e-3f);
    EXPECT_NEAR(rcs::unpackHalf(rcs::packHalf(-0.250f)), -0.250f, 1.0e-3f);
    EXPECT_NEAR(rcs::unpackHalf(rcs::packHalf(1.5f)), 1.5f, 1.0e-3f);
}

TEST(packHalf, NonFiniteStaysFinite) {
    EXPECT_TRUE(std::isfinite(rcs::unpackHalf(rcs::packHalf(std::numeric_limits<float>::quiet_NaN()))));
    EXPECT_TRUE(std::isfinite(rcs::unpackHalf(rcs::packHalf(std::numeric_limits<float>::infinity()))));
    EXPECT_TRUE(std::isfinite(rcs::unpackHalf(rcs::packHalf(-std::numeric_limits<float>::infinity()))));
    EXPECT_TRUE(std::isfinite(rcs::unpackHalf(rcs::packHalf(1.0e20f))));
}
