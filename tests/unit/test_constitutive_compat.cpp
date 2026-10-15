// Issue #2: constitutive model <-> element family compatibility.
//
// Today the constitutive law is chosen in three places that do not consult
// each other:
//   shells (part.cpp, Reissner quads / BST triangles): always linear
//       ThermalReissnerElasticity / ThermalKirchhoff, for every material
//   tets   (part.cpp configure + thermal_tet.cpp trialStress): Mooney-Rivlin
//       whenever the yield handed to configure() is 0 and Rankine is off
//   mesher (physics_units.cpp volumeHyperelastic): Compliant && yield <= 0
//       && nu >= 0.45, used to route parts in mesh_split.cpp
// So a rubber part is Mooney-Rivlin when it is bulky and linear elastic when it
// is slender, and two different predicates disagree about what "hyperelastic"
// means. Tests marked FAILS TODAY pin each of those down.

#include "physics_core/cad.hpp"
#include "physics_core/material.hpp"
#include "physics_core/mesh_split.hpp"
#include "physics_core/part.hpp"
#include "physics_core/physics_units.hpp"

#include <gtest/gtest.h>

#include <sstream>
#include <string>
#include <vector>

namespace {

/// The tet's rule, copied from Part::buildVolume (part.cpp, `configure(...,
/// allowPlasticity ? yieldStress : 0.0, ...)` then `setRankine`) and
/// ThermalTet4::trialStress (`hyper = _yield <= 0 && no Rankine`).
/// If this helper drifts from the source, that is the bug in miniature:
/// the decision should be one exported function both sides call.
bool tetUsesMooneyRivlin(const rcs::Material& m) {
    const double yield = m.allowPlasticity ? m.yieldStress : 0.0;
    return yield <= 0.0 && rcs::rankineTensionLimit(m) <= 0.0 && rcs::rankineCompressionLimit(m) <= 0.0;
}

struct Case {
    std::string label;
    rcs::PartDefinition part;
};

rcs::PartDefinition rubberPatch(double L, double t) {
    rcs::PartDefinition sheet;
    sheet.name = "rubber patch";
    sheet.material = rcs::rubber();
    sheet.thickness = t;
    sheet.patch.origin = {-0.5 * L, 0.0, -0.5 * L};
    sheet.patch.uAxis = {1.0, 0.0, 0.0};
    sheet.patch.vAxis = {0.0, 0.0, 1.0};
    sheet.patch.uSize = L;
    sheet.patch.vSize = L;
    return sheet;
}

}  // namespace

// Acceptance (b): a continuum mode with no valid shell formulation is never
// placed on a shell. Either the slender geometry is meshed as a (hybrid /
// adequately layered) solid, or ensurePartMesh refuses it. Until a
// finite-strain shell exists, "shell" is not an allowed answer for rubber.
// FAILS TODAY: every slender case below becomes a linear Reissner/BST shell.
TEST(ConstitutiveCompat, FiniteStrainModeNeverLandsOnSmallStrainShell) {
    std::vector<Case> cases;
    for (double t : {0.004, 0.008, 0.02, 0.05})
        cases.push_back({"patch 0.4 m, t=" + std::to_string(t), rubberPatch(0.4, t)});
    rcs::TestCubeOptions options;
    options.material = rcs::rubber();
    const rcs::AssemblyDefinition box = rcs::makeTestCube(options);
    cases.push_back({"test-cube panel " + box.parts.front().name, box.parts.front()});

    for (Case& c : cases) {
        rcs::applyResolutionSizing(c.part, 1.0);
        const bool ok = rcs::ensurePartMesh(c.part);
        if (!ok)
            continue;  // refusing is a valid outcome under the policy
        rcs::Part built;
        built.build(c.part);
        EXPECT_TRUE(built.hasVolume())
            << c.label << ": " << c.part.material.name
            << " became a shell. Shell elements use linear ThermalKirchhoff/ThermalReissnerElasticity, "
               "while the same material on tets uses mooneyRivlinCauchy.";
    }
}

// The control for the test above: bulky rubber already gets tets, so the
// difference between slender and bulky is purely the geometry branch.
// PASSES TODAY.
TEST(ConstitutiveCompat, BulkyRubberUsesHyperelasticTets) {
    rcs::PartDefinition cube = rcs::makeSolidCube(0.08);
    cube.material = rcs::rubber();
    rcs::applyResolutionSizing(cube, 1.0);
    ASSERT_TRUE(rcs::ensurePartMesh(cube));
    rcs::Part built;
    built.build(cube);
    EXPECT_TRUE(built.hasVolume());
    EXPECT_TRUE(tetUsesMooneyRivlin(cube.material));
}

// Goal 1: one dispatch, not two predicates. Whatever the mesher believes about
// a material must be what the element then does.
// FAILS TODAY for:
//   - aluminium with plasticity switched off in the editor: the mesher treats
//     it as a metal, the tet silently runs Mooney-Rivlin;
//   - an elastomer with nu = 0.40: volumeHyperelastic says no (nu < 0.45), the
//     tet says yes (yield 0, no Rankine).
TEST(ConstitutiveCompat, MesherAndTetAgreeOnHyperelastic) {
    std::vector<rcs::Material> materials = rcs::materialPresets();
    rcs::Material elasticAluminium = rcs::aluminum6061();
    elasticAluminium.name = "Aluminium, plasticity off";
    elasticAluminium.allowPlasticity = false;
    materials.push_back(elasticAluminium);
    rcs::Material softElastomer = rcs::rubber();
    softElastomer.name = "Elastomer nu=0.40";
    softElastomer.poissonRatio = 0.40;
    materials.push_back(softElastomer);

    for (const rcs::Material& m : materials) {
        EXPECT_EQ(rcs::volumeHyperelastic(m), tetUsesMooneyRivlin(m))
            << m.name << ": mesher volumeHyperelastic=" << rcs::volumeHyperelastic(m)
            << ", tet Mooney-Rivlin=" << tetUsesMooneyRivlin(m);
    }
}

// Goal 4 / non-goal "re-tuning one preset": a 0.001 change in Poisson's ratio
// is not a change of continuum mode, so it must not change the element family.
// FAILS TODAY: at nu = 0.449 the box panel is an authored tet solid; at 0.45
// volumeHyperelastic flips, the lockProne/plateFabric branches fire and the
// same panel becomes a shell.
TEST(ConstitutiveCompat, PoissonCliffDoesNotSwitchTheFamily) {
    auto panelHasVolume = [](double nu) {
        rcs::TestCubeOptions options;
        options.material = rcs::rubber();
        options.material.poissonRatio = nu;
        rcs::AssemblyDefinition box = rcs::makeTestCube(options);
        rcs::PartDefinition panel = box.parts.front();
        rcs::applyResolutionSizing(panel, 1.0);
        EXPECT_TRUE(rcs::ensurePartMesh(panel));
        return !panel.volume.empty();
    };
    const bool below = panelHasVolume(0.449);
    const bool above = panelHasVolume(0.450);
    EXPECT_EQ(below, above) << "nu=0.449 -> " << (below ? "solid" : "shell") << ", nu=0.450 -> "
                            << (above ? "solid" : "shell");
}

// Acceptance: "Metals' existing ductile shell/solid paths remain valid under
// the same table." Guard rail so the fix does not regress them.
// PASSES TODAY.
TEST(ConstitutiveCompat, DuctileMetalsStillMeshBothWays) {
    for (const rcs::Material& m : {rcs::aluminum6061(), rcs::mildSteel(), rcs::ar500()}) {
        rcs::PartDefinition sheet = rubberPatch(0.4, 0.008);
        sheet.material = m;
        rcs::applyResolutionSizing(sheet, 1.0);
        EXPECT_TRUE(rcs::ensurePartMesh(sheet)) << m.name << " thin sheet";

        rcs::PartDefinition cube = rcs::makeSolidCube(0.08);
        cube.material = m;
        rcs::applyResolutionSizing(cube, 1.0);
        ASSERT_TRUE(rcs::ensurePartMesh(cube)) << m.name << " cube";
        EXPECT_FALSE(cube.volume.empty()) << m.name << " cube";
        EXPECT_FALSE(tetUsesMooneyRivlin(m)) << m.name;
    }
}

// ---------------------------------------------------------------------------
// Target API. These compile only once the dispatch table exists, so they do
// not break the build today. Suggested shape:
//
//   // physics_core/element_policy.hpp
//   enum class ElementFamily { ReissnerShell, KirchhoffTriShell, LinearTet, HybridTet };
//   enum class ConstitutiveLaw { LinearElastic, J2Plastic, Rankine, MooneyRivlin };
//   struct Compatibility { bool allowed; ConstitutiveLaw law; std::string reason; };
//   /// Single source of truth. Keyed on Material fields, never Material::name.
//   Compatibility compatibility(const Material& material, ElementFamily family);
// ---------------------------------------------------------------------------
#if __has_include("physics_core/element_policy.hpp")
#include "physics_core/element_policy.hpp"

TEST(ElementPolicy, EveryPresetHasAnAnswerForEveryFamily) {
    const rcs::ElementFamily families[] = {rcs::ElementFamily::ReissnerShell, rcs::ElementFamily::KirchhoffTriShell,
                                           rcs::ElementFamily::LinearTet, rcs::ElementFamily::HybridTet};
    for (const rcs::Material& m : rcs::materialPresets()) {
        for (rcs::ElementFamily f : families) {
            const rcs::Compatibility c = rcs::compatibility(m, f);
            if (!c.allowed)
                EXPECT_FALSE(c.reason.empty()) << m.name << ": a refusal must say what to do instead";
        }
    }
}

TEST(ElementPolicy, RubberRefusedOnSmallStrainShells) {
    const rcs::Material m = rcs::rubber();
    EXPECT_FALSE(rcs::compatibility(m, rcs::ElementFamily::ReissnerShell).allowed);
    EXPECT_FALSE(rcs::compatibility(m, rcs::ElementFamily::KirchhoffTriShell).allowed);
    EXPECT_TRUE(rcs::compatibility(m, rcs::ElementFamily::HybridTet).allowed);
}

TEST(ElementPolicy, LawDependsOnFieldsNotName) {
    rcs::Material renamed = rcs::rubber();
    renamed.name = "Anything";
    for (rcs::ElementFamily f : {rcs::ElementFamily::LinearTet, rcs::ElementFamily::HybridTet}) {
        EXPECT_EQ(rcs::compatibility(renamed, f).allowed, rcs::compatibility(rcs::rubber(), f).allowed);
        EXPECT_EQ(rcs::compatibility(renamed, f).law, rcs::compatibility(rcs::rubber(), f).law);
    }
}

TEST(ElementPolicy, MetalsAllowedOnShellsAndTets) {
    for (const rcs::Material& m : {rcs::aluminum6061(), rcs::mildSteel(), rcs::ar500()}) {
        EXPECT_TRUE(rcs::compatibility(m, rcs::ElementFamily::ReissnerShell).allowed) << m.name;
        EXPECT_TRUE(rcs::compatibility(m, rcs::ElementFamily::LinearTet).allowed) << m.name;
        EXPECT_EQ(rcs::compatibility(m, rcs::ElementFamily::LinearTet).law, rcs::ConstitutiveLaw::J2Plastic);
    }
}
#endif
