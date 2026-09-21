#pragma once

#include <string>

namespace rcs {

/// Isotropic structural material. Strength values are used for weld/fastener
/// failure checks and for scaling the stress colour map.
struct Material {
    std::string name = "Aluminum 6061-T6";
    double density = 2700.0;         ///< kg/m^3
    double youngsModulus = 6.89e10;  ///< Pa
    double poissonRatio = 0.33;
    double yieldStress = 2.76e8;     ///< Pa
    double ultimateStress = 3.10e8;  ///< Pa
    /// Allowable shear stress in a weld bead, fraction-of-parent convention.
    double weldShearStrength = 1.65e8;  ///< Pa
    /// Stiffness-proportional structural damping (Rayleigh beta), seconds.
    double dampingBeta = 6.0e-4;
};

inline Material aluminum6061() {
    return Material{};
}

inline Material mildSteel() {
    Material m;
    m.name = "Mild steel A36";
    m.density = 7850.0;
    m.youngsModulus = 2.0e11;
    m.poissonRatio = 0.29;
    m.yieldStress = 2.50e8;
    m.ultimateStress = 4.00e8;
    m.weldShearStrength = 2.10e8;
    m.dampingBeta = 4.0e-4;
    return m;
}

inline Material ar500() {
    Material m;
    m.name = "AR500 wear plate";
    m.density = 7850.0;
    m.youngsModulus = 2.05e11;
    m.poissonRatio = 0.29;
    m.yieldStress = 1.30e9;
    m.ultimateStress = 1.60e9;
    m.weldShearStrength = 4.00e8;
    m.dampingBeta = 3.0e-4;
    return m;
}

}  // namespace rcs
