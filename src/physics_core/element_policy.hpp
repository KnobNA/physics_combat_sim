#pragma once

#include "physics_core/material.hpp"

#include <string>

namespace rcs {

/// Element families the mesher can choose. Shells here are the small-strain
/// elements already in the solver. HybridTet is the near-incompressible solid
/// a hyperelastic law is allowed to use; that element is not built yet.
enum class ElementFamily { ReissnerShell, KirchhoffTriShell, LinearTet, HybridTet };

enum class ConstitutiveLaw { LinearElastic, J2Plastic, Rankine, MooneyRivlin };

struct Compatibility {
    bool allowed = false;
    ConstitutiveLaw law = ConstitutiveLaw::LinearElastic;
    /// Set whenever `allowed` is false, so a refusal says what to change.
    std::string reason;
};

/// The tet's Mooney–Rivlin switch: yield handed to the element is 0 when
/// plasticity is off, and the law runs when that yield is 0 and Rankine is off.
/// Mesher and element both use this. It does not look at Poisson's ratio or
/// at Material::name.
inline bool usesHyperelasticLaw(const Material& material) {
    const double yield = material.allowPlasticity ? material.yieldStress : 0.0;
    return yield <= 0.0 && rankineTensionLimit(material) <= 0.0 && rankineCompressionLimit(material) <= 0.0;
}

/// Material fields to the constitutive law for one element family.
/// A refusal carries a reason. The answer does not depend on Material::name.
inline Compatibility compatibility(const Material& material, ElementFamily family) {
    const bool shell = family == ElementFamily::ReissnerShell || family == ElementFamily::KirchhoffTriShell;
    if (usesHyperelasticLaw(material)) {
        if (shell) {
            return {false, ConstitutiveLaw::MooneyRivlin,
                    "No finite-strain shell for this material. Use a solid, or thicken the part."};
        }
        return {true, ConstitutiveLaw::MooneyRivlin, {}};
    }
    if (material.fracture == FractureMode::Brittle) {
        return {true, ConstitutiveLaw::Rankine, {}};
    }
    if (material.allowPlasticity && material.yieldStress > 0.0) {
        return {true, ConstitutiveLaw::J2Plastic, {}};
    }
    return {true, ConstitutiveLaw::LinearElastic, {}};
}

}  // namespace rcs
