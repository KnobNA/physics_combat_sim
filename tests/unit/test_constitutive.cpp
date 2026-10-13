#include "physics_core/physics_units.hpp"
#include "physics_core/thermal_tet.hpp"

#include "chrono/core/ChMatrix.h"
#include "chrono/core/ChMatrix33.h"
#include "chrono/fea/ChNodeFEAxyz.h"

#include <cmath>
#include <memory>
#include <gtest/gtest.h>

TEST(vonMisesStress, UniaxialTensionEqualsAbsSigma) {
    EXPECT_NEAR(rcs::vonMisesStress(1.0e8, 0.0, 0.0, 0.0, 0.0, 0.0), 1.0e8, 1.0e-6);
    EXPECT_NEAR(rcs::vonMisesStress(-2.0e8, 0.0, 0.0, 0.0, 0.0, 0.0), 2.0e8, 1.0e-6);
}

TEST(vonMisesStress, ChVectorDynamicUniaxial) {
    chrono::ChVectorDynamic<> stress(6);
    stress.setZero();
    stress(0) = 5.0e7;
    EXPECT_NEAR(rcs::vonMisesStress(stress), 5.0e7, 1.0e-6);
}

TEST(vonMisesPlaneStress, UniaxialEqualsAbsSigma) {
    EXPECT_NEAR(rcs::vonMisesPlaneStress(3.0e7, 0.0, 0.0), 3.0e7, 1.0e-9);
    EXPECT_DOUBLE_EQ(rcs::vonMisesPlaneStress(0.0, 0.0, 0.0), 0.0);
}

TEST(equivalentPlastic, UniaxialDeviators) {
    const double strain[6] = {0.1, -0.05, -0.05, 0.0, 0.0, 0.0};
    EXPECT_NEAR(rcs::equivalentPlastic(strain), 0.1, 1.0e-12);
}

TEST(projectRankine, ClipsOnlyOverLimitPrincipal) {
    double stress[6] = {9.0e7, 1.0e7, -2.0e7, 0.0, 0.0, 0.0};
    ASSERT_TRUE(rcs::projectRankine(stress, 4.5e7, 8.0e8));
    EXPECT_NEAR(stress[0], 4.5e7, 1.0);
    EXPECT_NEAR(stress[1], 1.0e7, 1.0);
    EXPECT_NEAR(stress[2], -2.0e7, 1.0);
    EXPECT_NEAR(stress[3], 0.0, 1.0);
    EXPECT_NEAR(stress[4], 0.0, 1.0);
    EXPECT_NEAR(stress[5], 0.0, 1.0);
}

TEST(projectRankine, UntouchedInsideLimits) {
    double stress[6] = {1.0e7, 0.0, 0.0, 0.0, 0.0, 0.0};
    EXPECT_FALSE(rcs::projectRankine(stress, 4.5e7, 8.0e8));
    EXPECT_DOUBLE_EQ(stress[0], 1.0e7);
}

TEST(henckyVoigt, UniaxialLogStretch) {
    chrono::ChMatrix33<> U;
    U.setZero();
    U(0, 0) = std::exp(1.0);
    U(1, 1) = 1.0;
    U(2, 2) = 1.0;
    chrono::ChVectorDynamic<> strain;
    rcs::henckyVoigt(U, strain);
    EXPECT_NEAR(strain(0), 1.0, 1.0e-9);
    EXPECT_NEAR(strain(1), 0.0, 1.0e-9);
    EXPECT_NEAR(strain(2), 0.0, 1.0e-9);
    EXPECT_NEAR(strain(3), 0.0, 1.0e-9);
}

TEST(henckyVoigt, IdentityIsZero) {
    chrono::ChMatrix33<> U;
    U.setIdentity();
    chrono::ChVectorDynamic<> strain;
    rcs::henckyVoigt(U, strain);
    EXPECT_NEAR(strain(0), 0.0, 1.0e-12);
    EXPECT_NEAR(strain(1), 0.0, 1.0e-12);
    EXPECT_NEAR(strain(2), 0.0, 1.0e-12);
}

namespace {

void bindTet(rcs::ThermalTet4& tet,
             std::shared_ptr<chrono::fea::ChNodeFEAxyz> a,
             std::shared_ptr<chrono::fea::ChNodeFEAxyz> b,
             std::shared_ptr<chrono::fea::ChNodeFEAxyz> c,
             std::shared_ptr<chrono::fea::ChNodeFEAxyz> d) {
    tet.SetNodes(std::move(a), std::move(b), std::move(c), std::move(d));
    tet.configure(6.9e10, 0.33, 2.4e8, 2700.0, 0.0);
    tet.SetupInitial(nullptr);
}

}  // namespace

TEST(ThermalTet4, SkipsPolarOnSingularF) {
    using chrono::ChVector3d;
    auto n0 = chrono_types::make_shared<chrono::fea::ChNodeFEAxyz>(ChVector3d(0, 0, 0));
    auto n1 = chrono_types::make_shared<chrono::fea::ChNodeFEAxyz>(ChVector3d(0.1, 0, 0));
    auto n2 = chrono_types::make_shared<chrono::fea::ChNodeFEAxyz>(ChVector3d(0, 0.1, 0));
    auto n3 = chrono_types::make_shared<chrono::fea::ChNodeFEAxyz>(ChVector3d(0, 0, 0.1));
    rcs::ThermalTet4 tet;
    bindTet(tet, n0, n1, n2, n3);
    n3->SetPos(ChVector3d(0.05, 0.05, 0.0));
    chrono::ChVectorDynamic<> fi(12);
    fi.setZero();
    tet.ComputeInternalForces(fi);
    EXPECT_TRUE(tet.deformationSingular());
    EXPECT_FALSE(tet.inverted());
    ASSERT_EQ(fi.size(), 12);
    for (int i = 0; i < 12; ++i)
        EXPECT_DOUBLE_EQ(fi(i), 0.0);
}

TEST(ThermalTet4, SkipsForceOnInvertedF) {
    using chrono::ChVector3d;
    auto n0 = chrono_types::make_shared<chrono::fea::ChNodeFEAxyz>(ChVector3d(0, 0, 0));
    auto n1 = chrono_types::make_shared<chrono::fea::ChNodeFEAxyz>(ChVector3d(0.1, 0, 0));
    auto n2 = chrono_types::make_shared<chrono::fea::ChNodeFEAxyz>(ChVector3d(0, 0.1, 0));
    auto n3 = chrono_types::make_shared<chrono::fea::ChNodeFEAxyz>(ChVector3d(0, 0, 0.1));
    rcs::ThermalTet4 tet;
    bindTet(tet, n0, n1, n2, n3);
    n3->SetPos(ChVector3d(0, 0, -0.1));
    chrono::ChVectorDynamic<> fi(12);
    fi.setZero();
    tet.ComputeInternalForces(fi);
    EXPECT_TRUE(tet.inverted());
    EXPECT_FALSE(tet.deformationSingular());
    ASSERT_EQ(fi.size(), 12);
    for (int i = 0; i < 12; ++i)
        EXPECT_DOUBLE_EQ(fi(i), 0.0);
}

TEST(flowStress, OffsetYieldAndUts) {
    const rcs::Material al = rcs::aluminum6061();
    EXPECT_NEAR(rcs::flowStress(al.yieldStress, al.ultimateStress, al.failPlasticStrain, 0.0),
                0.90 * al.yieldStress, 1.0);
    EXPECT_NEAR(rcs::flowStress(al.yieldStress, al.ultimateStress, al.failPlasticStrain, 0.002),
                al.yieldStress, 1.0);
    EXPECT_NEAR(rcs::flowStress(al.yieldStress, al.ultimateStress, al.failPlasticStrain,
                                al.failPlasticStrain),
                al.ultimateStress, 1.0);
    EXPECT_GT(rcs::flowStress(al.yieldStress, al.ultimateStress, al.failPlasticStrain, 0.0),
              0.30 * al.yieldStress);
}

TEST(ductileTearStretch, TwiceHandbookElongation) {
    EXPECT_NEAR(rcs::ductileTearStretch(0.12, 1.0), 1.24, 1.0e-12);
    EXPECT_NEAR(rcs::ductileTearStretch(0.20, 0.5), 1.20, 1.0e-12);
}

TEST(ThermalTet4, OffsetYieldProducesPlastic) {
    using chrono::ChVector3d;
    const rcs::Material al = rcs::aluminum6061();
    const double L = 0.1;
    auto n0 = chrono_types::make_shared<chrono::fea::ChNodeFEAxyz>(ChVector3d(0, 0, 0));
    auto n1 = chrono_types::make_shared<chrono::fea::ChNodeFEAxyz>(ChVector3d(L, 0, 0));
    auto n2 = chrono_types::make_shared<chrono::fea::ChNodeFEAxyz>(ChVector3d(0, L, 0));
    auto n3 = chrono_types::make_shared<chrono::fea::ChNodeFEAxyz>(ChVector3d(0, 0, L));
    rcs::ThermalTet4 tet;
    tet.SetNodes(n0, n1, n2, n3);
    tet.configure(al.youngsModulus, al.poissonRatio, al.yieldStress, al.density, 0.0);
    tet.setFailureStrain(al.failPlasticStrain);
    tet.setHardening(al.ultimateStress);
    tet.SetupInitial(nullptr);
    const double nu = al.poissonRatio;
    const double eps = al.yieldStress / al.youngsModulus + 0.002;
    n1->SetPos(ChVector3d(L * (1.0 + eps), 0, 0));
    n2->SetPos(ChVector3d(0, L * (1.0 - nu * eps), 0));
    n3->SetPos(ChVector3d(0, 0, L * (1.0 - nu * eps)));
    chrono::ChVectorDynamic<> fi(12);
    fi.setZero();
    tet.ComputeInternalForces(fi);
    EXPECT_GE(tet.plasticStrain(), 0.001);
    EXPECT_LE(tet.plasticStrain(), 0.012);
}

TEST(ThermalTet4, BelowYieldStaysElastic) {
    using chrono::ChVector3d;
    const rcs::Material al = rcs::aluminum6061();
    const double L = 0.1;
    auto n0 = chrono_types::make_shared<chrono::fea::ChNodeFEAxyz>(ChVector3d(0, 0, 0));
    auto n1 = chrono_types::make_shared<chrono::fea::ChNodeFEAxyz>(ChVector3d(L, 0, 0));
    auto n2 = chrono_types::make_shared<chrono::fea::ChNodeFEAxyz>(ChVector3d(0, L, 0));
    auto n3 = chrono_types::make_shared<chrono::fea::ChNodeFEAxyz>(ChVector3d(0, 0, L));
    rcs::ThermalTet4 tet;
    tet.SetNodes(n0, n1, n2, n3);
    tet.configure(al.youngsModulus, al.poissonRatio, al.yieldStress, al.density, 0.0);
    tet.setFailureStrain(al.failPlasticStrain);
    tet.setHardening(al.ultimateStress);
    tet.SetupInitial(nullptr);
    const double nu = al.poissonRatio;
    const double eps = 0.30 * al.yieldStress / al.youngsModulus;
    n1->SetPos(ChVector3d(L * (1.0 + eps), 0, 0));
    n2->SetPos(ChVector3d(0, L * (1.0 - nu * eps), 0));
    n3->SetPos(ChVector3d(0, 0, L * (1.0 - nu * eps)));
    chrono::ChVectorDynamic<> fi(12);
    fi.setZero();
    tet.ComputeInternalForces(fi);
    EXPECT_LT(tet.plasticStrain(), 1.0e-6);
    EXPECT_NEAR(tet.vonMises() / (0.30 * al.yieldStress), 1.0, 0.20);
}

TEST(ThermalTet4, HealthyFRemainsUsable) {
    using chrono::ChVector3d;
    auto n0 = chrono_types::make_shared<chrono::fea::ChNodeFEAxyz>(ChVector3d(0, 0, 0));
    auto n1 = chrono_types::make_shared<chrono::fea::ChNodeFEAxyz>(ChVector3d(0.1, 0, 0));
    auto n2 = chrono_types::make_shared<chrono::fea::ChNodeFEAxyz>(ChVector3d(0, 0.1, 0));
    auto n3 = chrono_types::make_shared<chrono::fea::ChNodeFEAxyz>(ChVector3d(0, 0, 0.1));
    rcs::ThermalTet4 tet;
    bindTet(tet, n0, n1, n2, n3);
    chrono::ChVectorDynamic<> fi(12);
    fi.setZero();
    tet.ComputeInternalForces(fi);
    EXPECT_FALSE(tet.deformationSingular());
    EXPECT_FALSE(tet.inverted());
}

TEST(mooneyRivlinCauchy, SmallStrainMatchesYoungs) {
    const double E = 1.5e7;
    const double nu = 0.49;
    const double eps = 1.0e-4;
    chrono::ChMatrix33<> U;
    U.setZero();
    U(0, 0) = 1.0 + eps;
    U(1, 1) = 1.0 - nu * eps;
    U(2, 2) = 1.0 - nu * eps;
    chrono::ChVectorDynamic<> stress;
    rcs::mooneyRivlinCauchy(U, E, nu, 0.0, stress);
    EXPECT_NEAR(stress(0) / (E * eps), 1.0, 0.20);
}

TEST(mooneyRivlinCauchy, UniaxialEngineeringPRisesPastHenckyPeak) {
    const double E = 1.5e7;
    const double nu = 0.49;
    const double mu = E / (2.0 * (1.0 + nu));
    auto engineeringP = [&](double lam) {
        chrono::ChMatrix33<> U;
        U.setZero();
        U(0, 0) = lam;
        const double lat = 1.0 / std::sqrt(lam);
        U(1, 1) = lat;
        U(2, 2) = lat;
        chrono::ChVectorDynamic<> stress;
        rcs::mooneyRivlinCauchy(U, E, nu, 0.0, stress);
        return stress(0) / lam;
    };
    const double lamE = std::exp(1.0);
    const double pHencky = 2.0 * mu * std::log(lamE) / lamE;
    const double pMr = engineeringP(lamE);
    EXPECT_GT(pMr, pHencky);
    EXPECT_GT(engineeringP(lamE + 0.4), pMr);
}

TEST(mooneyRivlinCauchy, EquibiaxialDiffersFromUniaxial) {
    const double E = 1.5e7;
    const double nu = 0.49;
    const double lam = 1.4;
    chrono::ChMatrix33<> uni;
    uni.setZero();
    uni(0, 0) = lam;
    uni(1, 1) = 1.0 / std::sqrt(lam);
    uni(2, 2) = 1.0 / std::sqrt(lam);
    chrono::ChMatrix33<> bi;
    bi.setZero();
    bi(0, 0) = lam;
    bi(1, 1) = lam;
    bi(2, 2) = 1.0 / (lam * lam);
    chrono::ChVectorDynamic<> pUni;
    chrono::ChVectorDynamic<> pBi;
    rcs::mooneyRivlinCauchy(uni, E, nu, 0.0, pUni);
    rcs::mooneyRivlinCauchy(bi, E, nu, 0.0, pBi);
    EXPECT_GT(std::abs(pBi(0) - pUni(0)), 0.05 * std::max(std::abs(pUni(0)), 1.0));
}

TEST(mooneyRivlinCauchy, FreeThermalExpansionUnstressed) {
    const double th = 0.01;
    chrono::ChMatrix33<> U;
    U.setZero();
    U(0, 0) = U(1, 1) = U(2, 2) = 1.0 + th;
    chrono::ChVectorDynamic<> stress;
    rcs::mooneyRivlinCauchy(U, 1.5e7, 0.49, th, stress);
    EXPECT_NEAR(stress(0), 0.0, 1.0e3);
    EXPECT_NEAR(stress(1), 0.0, 1.0e3);
    EXPECT_NEAR(stress(2), 0.0, 1.0e3);
}

TEST(volumeHyperelastic, RubberIsVolumeNotPlate) {
    EXPECT_TRUE(rcs::volumeHyperelastic(rcs::rubber()));
    EXPECT_FALSE(rcs::volumeHyperelastic(rcs::aluminum6061()));
    EXPECT_FALSE(rcs::volumeHyperelastic(rcs::glass()));
}
