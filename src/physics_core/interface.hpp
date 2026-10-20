#pragma once

namespace rcs {

/// Rating of a weld or adhesive seam. Capacity belongs to the interface.
/// A zero field means "derive it from the two parents."
struct InterfaceDefinition {
    double shearStrength = 0.0;    ///< Pa, 0 = derive from parents
    double peelStrength = 0.0;     ///< Pa
    double releaseRate = 0.0;      ///< G_c, J/m^2, 0 = derive (weaker parent)
    bool boundedByParents = true;  ///< cap at the weaker parent's allowable
};

}  // namespace rcs
