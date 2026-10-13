#include "physics_core/joint.hpp"

#include <gtest/gtest.h>

TEST(brokenCount, CountsFailedLinks) {
    rcs::Joint joint;
    joint.links.resize(3);
    joint.links[1].broken = true;
    EXPECT_EQ(joint.brokenCount(), 1);
}

TEST(severed, EmptyOrAllBroken) {
    rcs::Joint empty;
    EXPECT_TRUE(empty.severed());

    rcs::Joint joint;
    joint.links.resize(2);
    EXPECT_FALSE(joint.severed());
    joint.links[0].broken = true;
    joint.links[1].broken = true;
    EXPECT_TRUE(joint.severed());
}

TEST(peakUtilisation, IgnoresBrokenLinks) {
    rcs::Joint joint;
    joint.links.resize(3);
    joint.links[0].utilisation = 0.4;
    joint.links[1].utilisation = 1.2;
    joint.links[1].broken = true;
    joint.links[2].utilisation = 0.8;
    EXPECT_DOUBLE_EQ(joint.peakUtilisation(), 0.8);
}
