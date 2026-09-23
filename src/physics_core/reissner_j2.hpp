#pragma once

#include "chrono/fea/ChMaterialShellReissner.h"

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

private:
    double _thermal = 0.0;
};

/// J2 / Ilyushin return mapping on the shell resultants. The trial stress
/// comes from the layer's isotropic elasticity. Over yield, the elastic
/// strain is scaled back onto the von Mises surface at the outer fibres and
/// the difference is stored as plastic strain.
class ReissnerJ2 : public chrono::fea::ChPlasticityReissner {
public:
    ReissnerJ2(double yieldStress, double youngsModulus);

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

private:
    double _yield = 0.0;
    double _young = 1.0;
};

}  // namespace rcs
