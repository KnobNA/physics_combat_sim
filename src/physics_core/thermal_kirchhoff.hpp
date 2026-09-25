#pragma once

#include "physics_core/material.hpp"

#include "chrono/fea/ChMaterialShellKirchhoff.h"

#include <algorithm>
#include <cmath>

namespace rcs {

/// Kirchhoff membrane/bending law with the same handbook thermal strain and
/// J2 cutoff the Reissner quads use. Imported triangle shells go through this
/// instead of the temperature-blind isotropic law.
class ThermalKirchhoff : public chrono::fea::ChElasticityKirchhoffIsothropic {
public:
    ThermalKirchhoff(double youngs, double nu)
        : chrono::fea::ChElasticityKirchhoffIsothropic(youngs, nu), _youngs(youngs) {}

    void setThermalStrain(double strain) { _thermal = strain; }
    void setYield(double yield) { _yield = std::max(0.0, yield); }
    void setScales(double stiffness, double strength) {
        _stiffness = std::clamp(stiffness, kThermalFloor, 1.0);
        _strength = std::clamp(strength, kThermalFloor, 1.0);
    }

    double plastic() const { return _plastic; }
    double stiffnessScale() const { return _stiffness; }

    void ComputeStress(chrono::ChVector3d& n,
                       chrono::ChVector3d& m,
                       const chrono::ChVector3d& eps,
                       const chrono::ChVector3d& kur,
                       const double z_inf,
                       const double z_sup,
                       const double angle) override {
        chrono::ChVector3d mechanical = eps;
        mechanical.x() -= _thermal;
        mechanical.y() -= _thermal;
        chrono::fea::ChElasticityKirchhoffIsothropic::ComputeStress(
            n, m, mechanical, kur, z_inf, z_sup, angle);
        n *= _stiffness;
        m *= _stiffness;

        const double h = std::max(z_sup - z_inf, 1.0e-8);
        const double sxx = n.x() / h;
        const double syy = n.y() / h;
        const double sxy = n.z() / h;
        const double vm = std::sqrt(sxx * sxx - sxx * syy + syy * syy + 3.0 * sxy * sxy);
        const double yield = _yield * _strength;
        if (yield > 0.0 && vm > yield && std::isfinite(vm)) {
            const double scale = yield / vm;
            n *= scale;
            m *= scale;
            const double extra = (vm - yield) / std::max(_youngs * _stiffness, 1.0);
            if (extra > _plastic)
                _plastic = extra;
        }
    }

    void ComputeStiffnessMatrix(chrono::ChMatrixRef mC,
                                const chrono::ChVector3d& eps,
                                const chrono::ChVector3d& kur,
                                const double z_inf,
                                const double z_sup,
                                const double angle) override {
        chrono::ChVector3d mechanical = eps;
        mechanical.x() -= _thermal;
        mechanical.y() -= _thermal;
        chrono::fea::ChElasticityKirchhoffIsothropic::ComputeStiffnessMatrix(
            mC, mechanical, kur, z_inf, z_sup, angle);
        mC *= _stiffness;
    }

private:
    double _youngs = 0.0;
    double _thermal = 0.0;
    double _yield = 0.0;
    double _stiffness = 1.0;
    double _strength = 1.0;
    double _plastic = 0.0;
};

}  // namespace rcs
