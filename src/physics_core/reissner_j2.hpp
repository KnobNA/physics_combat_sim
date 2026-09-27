#pragma once

#include "physics_core/material.hpp"

#include "chrono/core/ChMatrix.h"
#include "chrono/fea/ChMaterialShellReissner.h"

#include <algorithm>

namespace rcs {

/// Plastic membrane strain, plastic curvature, and accumulated equivalent
/// plastic strain for one integration point of a Reissner shell.
struct ReissnerPlasticData : chrono::fea::ChShellReissnerInternalData {
    chrono::ChVector3d epsU;
    chrono::ChVector3d epsV;
    chrono::ChVector3d kurU;
    chrono::ChVector3d kurV;

    void Copy(const chrono::fea::ChShellReissnerInternalData& other) override;
};

/// Isotropic Reissner elasticity with a uniform thermal membrane strain.
/// `alpha * (T - T_ref)` is subtracted from the normal strains before the
/// elastic law. Zero at the reference temperature, so an isothermal run is
/// unchanged. A later heat solver only has to call `SetThermalStrain`.
class ThermalReissnerElasticity : public chrono::fea::ChElasticityReissnerIsothropic {
public:
    ThermalReissnerElasticity(double youngsModulus, double poissonRatio);

    void SetThermalStrain(double strain) { _thermal = strain; }
    double thermalStrain() const { return _thermal; }
    /// Uniform thermal curvature κ_T ≈ α (T_top − T_bot) / h. Subtracted from
    /// the bending strains so a through-thickness gradient curls the midplane.
    void SetThermalCurvature(double kappa) { _thermalKur = kappa; }
    double thermalCurvature() const { return _thermalKur; }
    /// Scales the tangent and the stress. The floor keeps the shell invertible
    /// once the handbook curve has fallen past the solidus.
    void SetStiffnessScale(double scale) { _stiffness = std::clamp(scale, kThermalFloor, 1.0); }
    double stiffnessScale() const { return _stiffness; }

    void ComputeStress(chrono::ChVector3d& nU,
                       chrono::ChVector3d& nV,
                       chrono::ChVector3d& mU,
                       chrono::ChVector3d& mV,
                       const chrono::ChVector3d& epsU,
                       const chrono::ChVector3d& epsV,
                       const chrono::ChVector3d& kurU,
                       const chrono::ChVector3d& kurV,
                       double zInf,
                       double zSup,
                       double angle) override;

    void ComputeStiffnessMatrix(chrono::ChMatrixRef mC,
                                const chrono::ChVector3d& epsU,
                                const chrono::ChVector3d& epsV,
                                const chrono::ChVector3d& kurU,
                                const chrono::ChVector3d& kurV,
                                double zInf,
                                double zSup,
                                double angle) override;

private:
    double _thermal = 0.0;
    double _thermalKur = 0.0;
    double _stiffness = 1.0;
};

/// J2 / Ilyushin return mapping on the shell resultants. The trial stress
/// comes from the layer's isotropic elasticity. Over yield, the elastic
/// strain is scaled back onto the von Mises surface at the outer fibres and
/// the difference is stored as plastic strain.
class ReissnerJ2 : public chrono::fea::ChPlasticityReissner {
public:
    ReissnerJ2(double yieldStress, double youngsModulus);

    void SetYieldScale(double scale) { _yieldScale = std::clamp(scale, kThermalFloor, 1.0); }

    bool ComputeStressWithReturnMapping(
        chrono::ChVector3d& nU,
        chrono::ChVector3d& nV,
        chrono::ChVector3d& mU,
        chrono::ChVector3d& mV,
        chrono::fea::ChShellReissnerInternalData& dataNew,
        const chrono::ChVector3d& epsUTrial,
        const chrono::ChVector3d& epsVTrial,
        const chrono::ChVector3d& kurUTrial,
        const chrono::ChVector3d& kurVTrial,
        const chrono::fea::ChShellReissnerInternalData& data,
        double zInf,
        double zSup,
        double angle) override;

    void CreatePlasticityData(
        int numpoints,
        std::vector<std::unique_ptr<chrono::fea::ChShellReissnerInternalData>>& plasticData) override;

    /// Temperature-scaled elastic tangent. Chrono's numerical perfect-plastic
    /// tangent is singular once the stress sits on a very low yield surface,
    /// and the step then aborts. The return mapping still holds the stress
    /// on the real yield surface, shear included.
    void ComputeStiffnessMatrixElastoplastic(
        chrono::ChMatrixRef K,
        const chrono::ChVector3d& epsU,
        const chrono::ChVector3d& epsV,
        const chrono::ChVector3d& kurU,
        const chrono::ChVector3d& kurV,
        const chrono::fea::ChShellReissnerInternalData& data,
        double zInf,
        double zSup,
        double angle) override;

private:
    double _yield = 0.0;
    double _yieldScale = 1.0;
    double _young = 1.0;
};

/// Rayleigh damping that follows the current thermal stiffness. Chrono's own
/// Rayleigh model freezes the stiffness matrix on first use, so a panel heated
/// after the first step stays damped as if it were still cold.
class ThermalRayleighDamping : public chrono::fea::ChDampingReissner {
public:
    ThermalRayleighDamping(std::shared_ptr<ThermalReissnerElasticity> elasticity, double beta);

    void ComputeStress(chrono::ChVector3d& nU,
                       chrono::ChVector3d& nV,
                       chrono::ChVector3d& mU,
                       chrono::ChVector3d& mV,
                       const chrono::ChVector3d& depsU,
                       const chrono::ChVector3d& depsV,
                       const chrono::ChVector3d& dkurU,
                       const chrono::ChVector3d& dkurV,
                       double zInf,
                       double zSup,
                       double angle) override;

    void ComputeDampingMatrix(chrono::ChMatrixRef R,
                              const chrono::ChVector3d& depsU,
                              const chrono::ChVector3d& depsV,
                              const chrono::ChVector3d& dkurU,
                              const chrono::ChVector3d& dkurV,
                              double zInf,
                              double zSup,
                              double angle) override;

private:
    void refresh(double zInf, double zSup, double angle);

    std::shared_ptr<ThermalReissnerElasticity> _elasticity;
    chrono::ChMatrixNM<double, 12, 12> _stiffness;
    double _beta = 0.0;
    double _cachedScale = -1.0;
};

}  // namespace rcs
