#pragma once

#include <string>
#include <vector>

namespace rcs {

/// How a shell turns stress into a crack. Ductile metals yield and then tear.
/// Brittle solids crack on tensile principal stress. Compliant solids only
/// tear after a large stretch.
enum class FractureMode {
    Ductile,
    Brittle,
    Compliant
};

/// Isotropic structural material. Strength values drive weld failure, the
/// stress colour map, and fracture. Thermal fields are stored so a later
/// heat solver can write a temperature; with the part at the reference
/// temperature they change nothing.
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
    /// Edge stretch (current / rest) at which a ductile or compliant shell tears.
    double tearStretch = 1.12;
    /// Accumulated equivalent plastic strain at which a ductile shell tears.
    /// Zero disables this check.
    double failPlasticStrain = 0.12;
    /// False for materials with no plastic range (glass, rubber). Yield is
    /// ignored and the element stays elastic until it fractures or tears.
    bool allowPlasticity = true;

    FractureMode fracture = FractureMode::Ductile;
    /// Rankine tensile strength, Pa. Brittle shells crack above this.
    double fractureStrength = 0.0;
    /// Mode I toughness, Pa*sqrt(m). Zero skips the Griffith energy gate.
    double fractureToughness = 0.0;
    /// Tempered glass: the first crack runs down to the minimum shard size.
    bool shatterToDust = false;
    /// A detached piece smaller than this leaves the shell solve and becomes
    /// one rigid body. Metres. Annealed plate is about 2 cm; tempered is finer.
    double minimumShardSize = 0.02;

    double thermalExpansion = 23.6e-6;  ///< 1/K
    double specificHeat = 900.0;        ///< J/(kg*K)
    double conductivity = 167.0;        ///< W/(m*K)
    double referenceTemperature = 293.15;  ///< K
};

inline Material aluminum6061() {
    return Material{};
}

inline Material mildSteel() {
    Material m;
    m.name = "Mild steel A36";
    m.density = 7850.0;
    m.youngsModulus = 2.0e11;
    m.poissonRatio = 0.26;
    m.yieldStress = 2.50e8;
    m.ultimateStress = 4.00e8;
    m.weldShearStrength = 2.10e8;
    m.dampingBeta = 4.0e-4;
    m.tearStretch = 1.20;
    m.failPlasticStrain = 0.20;
    m.thermalExpansion = 12.0e-6;
    m.specificHeat = 490.0;
    m.conductivity = 50.0;
    return m;
}

inline Material glass() {
    Material m;
    m.name = "Soda-lime glass";
    m.density = 2500.0;
    m.youngsModulus = 7.0e10;
    m.poissonRatio = 0.23;
    m.yieldStress = 0.0;
    m.ultimateStress = 4.5e7;
    m.weldShearStrength = 1.0e7;
    m.dampingBeta = 2.0e-4;
    m.tearStretch = 1.002;
    m.failPlasticStrain = 0.0;
    m.allowPlasticity = false;
    m.fracture = FractureMode::Brittle;
    m.fractureStrength = 4.5e7;
    m.fractureToughness = 0.75e6;
    m.minimumShardSize = 0.02;
    m.thermalExpansion = 9.0e-6;
    m.specificHeat = 840.0;
    m.conductivity = 1.0;
    return m;
}

inline Material temperedGlass() {
    Material m = glass();
    m.name = "Tempered soda-lime";
    m.ultimateStress = 1.20e8;
    m.fractureStrength = 1.20e8;
    m.shatterToDust = true;
    m.minimumShardSize = 0.015;
    return m;
}

inline Material rubber() {
    Material m;
    m.name = "Hard elastomer";
    m.density = 1100.0;
    m.youngsModulus = 1.5e7;
    m.poissonRatio = 0.49;
    m.yieldStress = 0.0;
    m.ultimateStress = 2.5e7;
    m.weldShearStrength = 2.0e6;
    m.dampingBeta = 2.0e-3;
    m.tearStretch = 4.0;
    m.failPlasticStrain = 0.0;
    m.allowPlasticity = false;
    m.fracture = FractureMode::Compliant;
    m.thermalExpansion = 200.0e-6;
    m.specificHeat = 1500.0;
    m.conductivity = 0.15;
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
    m.tearStretch = 1.10;
    m.failPlasticStrain = 0.08;
    m.thermalExpansion = 12.0e-6;
    m.specificHeat = 490.0;
    m.conductivity = 50.0;
    return m;
}

/// Library order used by the assembly panel picker.
inline std::vector<Material> materialPresets() {
    return {aluminum6061(), mildSteel(), ar500(), glass(), temperedGlass(), rubber()};
}

inline int presetIndex(const Material& material) {
    const auto presets = materialPresets();
    for (int i = 0; i < static_cast<int>(presets.size()); ++i) {
        if (presets[i].name == material.name)
            return i;
    }
    return 0;
}

}  // namespace rcs
