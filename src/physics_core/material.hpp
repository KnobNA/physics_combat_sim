#pragma once

#include <algorithm>
#include <cmath>
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
/// stress colour map, and fracture. Expansion, heat capacity, and conductivity
/// feed the heat solver. Stiffness and yield follow handbook curves and keep
/// falling past the solidus, down to a small residual so the shell matrix
/// stays invertible. The solidus is where a weld lets go.
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
    /// Rankine compressive strength, Pa. Brittle solids crush above this.
    /// Zero means "derive from tensile strength" for brittle materials.
    double compressiveStrength = 0.0;
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
    /// Oxidized aluminum. A hot sphere heats a facing surface by radiation;
    /// polished metal would be lower, paint higher.
    double emissivity = 0.25;
    double referenceTemperature = 293.15;  ///< K, handbook room temperature
    /// Where rapid softening starts. Welds still use the solidus below.
    double softenTemperature = 420.0;  ///< K
    /// Solidus. A weld breaks when both parent nodes reach the lower solidus.
    double meltTemperature = 855.0;  ///< K
};

/// Numerical residual so a molten shell matrix is not exactly singular.
/// It is not a structural strength: an 8 mm aluminum sheet at this floor
/// yields under a few newtons. Handbook curves already sit here once the
/// metal is well above its solidus.
inline constexpr double kThermalFloor = 1.0e-4;

/// Slide a handbook curve onto this material's soften and solidus temperatures.
inline double mapThermalAxis(double kelvin, double soften, double melt, double refSoften, double refMelt) {
    constexpr double room = 293.15;
    soften = std::max(soften, room + 1.0);
    melt = std::max(melt, soften + 1.0);
    refSoften = std::max(refSoften, room + 1.0);
    refMelt = std::max(refMelt, refSoften + 1.0);
    if (kelvin <= soften) {
        const double u = (kelvin - room) / (soften - room);
        return room + u * (refSoften - room);
    }
    const double span = melt - soften;
    const double refSpan = refMelt - refSoften;
    if (kelvin >= melt)
        return refMelt + (kelvin - melt) * refSpan / span;
    const double u = (kelvin - soften) / span;
    return refSoften + u * refSpan;
}

inline double sampleThermalKnots(double kelvin, const double* temperatures, const double* values, int count) {
    if (count <= 0)
        return 1.0;
    if (kelvin <= temperatures[0])
        return values[0];
    for (int i = 1; i < count; ++i) {
        if (kelvin > temperatures[i])
            continue;
        const double span = temperatures[i] - temperatures[i - 1];
        const double u = span > 1.0 ? (kelvin - temperatures[i - 1]) / span : 0.0;
        return values[i - 1] + u * (values[i] - values[i - 1]);
    }
    const double span = std::max(temperatures[count - 1] - temperatures[0], 1.0);
    const double over = (kelvin - temperatures[count - 1]) / span;
    return values[count - 1] * std::exp(-3.5 * over);
}

/// Fraction of room-temperature Young's modulus still available at `kelvin`.
/// Aluminum follows Kaufman; steels follow Eurocode 3; glass and rubber hold
/// until their softening range and then collapse. Past the solidus the solid
/// modulus keeps falling toward `kThermalFloor`.
inline double thermalStiffnessRetention(const Material& material, double kelvin) {
    double axis = kelvin;
    const double* temperatures = nullptr;
    const double* values = nullptr;
    int count = 0;
    if (material.fracture == FractureMode::Brittle) {
        // Soda-lime stays stiff until the softening range, then slumps.
        static const double t[] = {293.15, 700.0, 750.0, 850.0, 1000.0};
        static const double v[] = {1.0, 0.96, 0.88, 0.35, 0.04};
        axis = mapThermalAxis(kelvin, material.softenTemperature, material.meltTemperature, 750.0, 1000.0);
        temperatures = t;
        values = v;
        count = 5;
    } else if (material.fracture == FractureMode::Compliant) {
        static const double t[] = {293.15, 350.0, 400.0, 450.0};
        static const double v[] = {1.0, 0.75, 0.30, 0.04};
        axis = mapThermalAxis(kelvin, material.softenTemperature, material.meltTemperature, 350.0, 450.0);
        temperatures = t;
        values = v;
        count = 4;
    } else if (material.meltTemperature >= 1400.0) {
        // Eurocode 3 Table 3.1, k_E. Structural steel is effectively gone
        // by 1200 °C, well before the liquidus.
        static const double t[] = {293.15, 373.0, 473.0, 573.0, 673.0, 773.0, 873.0,
                                   973.0, 1073.0, 1173.0, 1273.0, 1373.0, 1473.0};
        static const double v[] = {1.0, 1.0, 0.90, 0.80, 0.70, 0.60, 0.31, 0.13, 0.09, 0.0675, 0.045, 0.0225, 0.008};
        axis = mapThermalAxis(kelvin, material.softenTemperature, material.meltTemperature, 600.0, 1780.0);
        temperatures = t;
        values = v;
        count = 13;
    } else {
        // 6061 modulus, from Kaufman, Properties of Aluminum Alloys.
        static const double t[] = {293.15, 373.0, 473.0, 573.0, 673.0, 773.0, 855.0};
        static const double v[] = {1.0, 0.97, 0.91, 0.76, 0.52, 0.25, 0.08};
        axis = mapThermalAxis(kelvin, material.softenTemperature, material.meltTemperature, 420.0, 855.0);
        temperatures = t;
        values = v;
        count = 7;
    }
    return std::clamp(sampleThermalKnots(axis, temperatures, values, count), kThermalFloor, 1.0);
}

/// Fraction of room-temperature yield (or fracture strength) still available.
/// 6061-T6 loses most of its yield between 150 °C and 300 °C. Mild steel holds
/// until 400 °C (Eurocode 3 k_y) and then falls. Far above the solidus the
/// flow stress sits on `kThermalFloor`, so the shell warps instead of springing
/// back. This is still a solid: there is no liquid.
inline double thermalStrengthRetention(const Material& material, double kelvin) {
    double axis = kelvin;
    const double* temperatures = nullptr;
    const double* values = nullptr;
    int count = 0;
    if (material.fracture == FractureMode::Brittle) {
        static const double t[] = {293.15, 750.0, 900.0, 1000.0};
        static const double v[] = {1.0, 0.85, 0.40, 0.10};
        axis = mapThermalAxis(kelvin, material.softenTemperature, material.meltTemperature, 750.0, 1000.0);
        temperatures = t;
        values = v;
        count = 4;
    } else if (material.fracture == FractureMode::Compliant) {
        static const double t[] = {293.15, 350.0, 400.0, 450.0};
        static const double v[] = {1.0, 0.80, 0.35, 0.05};
        axis = mapThermalAxis(kelvin, material.softenTemperature, material.meltTemperature, 350.0, 450.0);
        temperatures = t;
        values = v;
        count = 4;
    } else if (material.meltTemperature >= 1400.0) {
        // Eurocode 3 k_y. Quenched wear plate tempers earlier than mild steel,
        // so a yield above 1 GPa walks the same table at a faster rate.
        static const double t[] = {293.15, 673.0, 773.0, 873.0, 973.0, 1073.0, 1173.0, 1273.0, 1373.0, 1473.0};
        static const double v[] = {1.0, 1.0, 0.78, 0.47, 0.23, 0.11, 0.06, 0.04, 0.02, 0.008};
        double lookup = kelvin;
        if (material.yieldStress >= 1.0e9)
            lookup = 293.15 + (kelvin - 293.15) * 1.45;
        axis = mapThermalAxis(lookup, material.softenTemperature, material.meltTemperature, 600.0, 1780.0);
        temperatures = t;
        values = v;
        count = 10;
    } else {
        // 6061-T6 yield, ASM. About 80% at 150 °C, 12% at 260 °C, a few
        // percent near the 582 °C solidus.
        static const double t[] = {293.15, 373.0, 422.0, 477.0, 533.0, 589.0, 700.0, 855.0};
        static const double v[] = {1.0, 0.95, 0.78, 0.38, 0.12, 0.06, 0.025, 0.012};
        axis = mapThermalAxis(kelvin, material.softenTemperature, material.meltTemperature, 420.0, 855.0);
        temperatures = t;
        values = v;
        count = 8;
    }
    return std::clamp(sampleThermalKnots(axis, temperatures, values, count), kThermalFloor, 1.0);
}

/// Fraction of the room-temperature tearing strain left at `kelvin`.
/// Ductility follows the strength curve: while the metal can still carry a
/// few percent of its yield it can still stretch, including through a thermal
/// buckle just past the solidus. Once the handbook strength has fallen to the
/// numerical floor, a small mechanical stretch tears the sheet. Uniform
/// thermal growth is not part of that stretch.
inline double thermalDuctilityRetention(const Material& material, double kelvin) {
    constexpr double floor = 0.04;
    const double strength = thermalStrengthRetention(material, kelvin);
    if (strength >= 0.05)
        return 1.0;
    if (strength <= kThermalFloor * 3.0)
        return floor;
    const double span = 0.05 - kThermalFloor * 3.0;
    const double u = (0.05 - strength) / span;
    return 1.0 + u * (floor - 1.0);
}

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
    m.emissivity = 0.80;
    m.softenTemperature = 600.0;
    m.meltTemperature = 1780.0;
    return m;
}

/// Tensile Rankine limit used by the constitutive cap and by shatter.
inline double rankineTensionLimit(const Material& material) {
    if (material.fracture != FractureMode::Brittle)
        return 0.0;
    if (material.fractureStrength > 0.0)
        return material.fractureStrength;
    return std::max(0.0, material.ultimateStress);
}

/// Compressive Rankine limit. Soda-lime handbook crush strength is about
/// 15–20× the tensile strength; an explicit `compressiveStrength` wins.
inline double rankineCompressionLimit(const Material& material) {
    if (material.fracture != FractureMode::Brittle)
        return 0.0;
    if (material.compressiveStrength > 0.0)
        return material.compressiveStrength;
    const double tension = rankineTensionLimit(material);
    return tension > 0.0 ? 15.0 * tension : 0.0;
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
    m.compressiveStrength = 8.0e8;
    m.fractureToughness = 0.75e6;
    m.minimumShardSize = 0.02;
    m.thermalExpansion = 9.0e-6;
    m.specificHeat = 840.0;
    m.conductivity = 1.0;
    m.emissivity = 0.92;
    m.softenTemperature = 750.0;
    m.meltTemperature = 1000.0;
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
    m.emissivity = 0.90;
    m.softenTemperature = 350.0;
    m.meltTemperature = 450.0;
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
    m.emissivity = 0.80;
    m.softenTemperature = 600.0;
    m.meltTemperature = 1780.0;
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
