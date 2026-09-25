#include "physics_core/reissner_j2.hpp"

#include <algorithm>
#include <cmath>

namespace rcs {
namespace {

double vonMisesPlane(double xx, double yy, double xy) {
    const double value = xx * xx - xx * yy + yy * yy + 3.0 * xy * xy;
    return std::sqrt(std::max(0.0, value));
}

/// Peak von Mises at the two outer fibres, from membrane forces and moments
/// per unit length. Section modulus of a plate is h^2/6.
double peakVonMises(const chrono::ChVector3d& nU,
                    const chrono::ChVector3d& nV,
                    const chrono::ChVector3d& mU,
                    const chrono::ChVector3d& mV,
                    double thickness)
{
    const double h = std::max(thickness, 1.0e-9);
    const double nxx = nU.x() / h;
    const double nyy = nV.y() / h;
    const double nxy = 0.5 * (nU.y() + nV.x()) / h;
    const double factor = 6.0 / (h * h);
    const double mxx = mU.y() * factor;
    const double myy = mV.x() * factor;
    const double mxy = 0.5 * (mU.x() + mV.y()) * factor;
    const double top = vonMisesPlane(nxx + mxx, nyy + myy, nxy + mxy);
    const double bottom = vonMisesPlane(nxx - mxx, nyy - myy, nxy - mxy);
    // Transverse shear and drilling are resultants too. Leaving them out of
    // the yield check lets a soft shell build a huge shear force while the
    // membrane still looks elastic, and the node then flies in one step.
    const double shear = std::sqrt(nU.z() * nU.z() + nV.z() * nV.z()) / h;
    const double drill = factor * std::sqrt(mU.z() * mU.z() + mV.z() * mV.z());
    const double plane = std::max(top, bottom);
    return std::sqrt(plane * plane + 3.0 * shear * shear + drill * drill);
}

const ReissnerPlasticData* asPlastic(const chrono::fea::ChShellReissnerInternalData& data) {
    return dynamic_cast<const ReissnerPlasticData*>(&data);
}

ReissnerPlasticData* asPlastic(chrono::fea::ChShellReissnerInternalData& data) {
    return dynamic_cast<ReissnerPlasticData*>(&data);
}

}  // namespace

void ReissnerPlasticData::Copy(const chrono::fea::ChShellReissnerInternalData& other) {
    chrono::fea::ChShellReissnerInternalData::Copy(other);
    if (const auto* src = asPlastic(other)) {
        epsU = src->epsU;
        epsV = src->epsV;
        kurU = src->kurU;
        kurV = src->kurV;
    }
}

ThermalReissnerElasticity::ThermalReissnerElasticity(double youngsModulus, double poissonRatio)
    : chrono::fea::ChElasticityReissnerIsothropic(youngsModulus, poissonRatio, 1.0, 0.1) {}

void ThermalReissnerElasticity::ComputeStress(chrono::ChVector3d& nU,
                                              chrono::ChVector3d& nV,
                                              chrono::ChVector3d& mU,
                                              chrono::ChVector3d& mV,
                                              const chrono::ChVector3d& epsU,
                                              const chrono::ChVector3d& epsV,
                                              const chrono::ChVector3d& kurU,
                                              const chrono::ChVector3d& kurV,
                                              double zInf,
                                              double zSup,
                                              double angle)
{
    chrono::ChVector3d eu = epsU;
    chrono::ChVector3d ev = epsV;
    eu.x() -= _thermal;
    ev.y() -= _thermal;
    chrono::fea::ChElasticityReissnerIsothropic::ComputeStress(nU, nV, mU, mV, eu, ev, kurU, kurV, zInf,
                                                               zSup, angle);
    if (_stiffness != 1.0) {
        nU *= _stiffness;
        nV *= _stiffness;
        mU *= _stiffness;
        mV *= _stiffness;
    }
}

void ThermalReissnerElasticity::ComputeStiffnessMatrix(chrono::ChMatrixRef mC,
                                                       const chrono::ChVector3d& epsU,
                                                       const chrono::ChVector3d& epsV,
                                                       const chrono::ChVector3d& kurU,
                                                       const chrono::ChVector3d& kurV,
                                                       double zInf,
                                                       double zSup,
                                                       double angle)
{
    chrono::fea::ChElasticityReissnerIsothropic::ComputeStiffnessMatrix(mC, epsU, epsV, kurU, kurV, zInf,
                                                                        zSup, angle);
    if (_stiffness != 1.0)
        mC *= _stiffness;
}

ReissnerJ2::ReissnerJ2(double yieldStress, double youngsModulus)
    : _yield(yieldStress), _young(std::max(youngsModulus, 1.0)) {}

bool ReissnerJ2::ComputeStressWithReturnMapping(
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
    double angle)
{
    const ReissnerPlasticData* oldState = asPlastic(data);
    ReissnerPlasticData* newState = asPlastic(dataNew);
    ReissnerPlasticData fallback;
    const ReissnerPlasticData& old = oldState ? *oldState : fallback;

    chrono::ChVector3d epsU = epsUTrial - old.epsU;
    chrono::ChVector3d epsV = epsVTrial - old.epsV;
    chrono::ChVector3d kurU = kurUTrial - old.kurU;
    chrono::ChVector3d kurV = kurVTrial - old.kurV;

    auto* elasticity = section ? section->GetElasticity().get() : nullptr;
    if (!elasticity) {
        nU = nV = mU = mV = chrono::ChVector3d(0, 0, 0);
        return false;
    }

    elasticity->ComputeStress(nU, nV, mU, mV, epsU, epsV, kurU, kurV, zInf, zSup, angle);
    const double thickness = std::abs(zSup - zInf);
    const double yield = _yield * _yieldScale;
    const double vm = peakVonMises(nU, nV, mU, mV, thickness);

    auto store = [&](const ReissnerPlasticData& state, bool yielded) {
        if (newState)
            *newState = state;
        else
            dataNew.p_strain_acc = state.p_strain_acc;
        return yielded;
    };

    if (!(vm > yield) || yield <= 0.0) {
        ReissnerPlasticData kept = old;
        kept.p_strain_acc = old.p_strain_acc;
        return store(kept, false);
    }

    double scale = std::clamp(yield / vm, 0.0, 1.0);
    double thermal = 0.0;
    if (const auto* law = dynamic_cast<const ThermalReissnerElasticity*>(elasticity))
        thermal = law->thermalStrain();
    const double halfH = 0.5 * thickness;
    const double trialShear = std::sqrt(epsU.z() * epsU.z() + epsV.z() * epsV.z());
    const double kurY = kurU.y();
    const double kurX = kurV.x();
    // Thermal strain is a stress offset, not plastic strain. Scale only the
    // mechanical part that remains after it is removed, or the next step
    // subtracts the full thermal strain again and the residual runs away.
    chrono::ChVector3d mechU = epsU;
    chrono::ChVector3d mechV = epsV;
    mechU.x() = thermal + scale * (epsU.x() - thermal);
    mechV.y() = thermal + scale * (epsV.y() - thermal);
    mechU.y() *= scale;
    mechU.z() *= scale;
    mechV.x() *= scale;
    mechV.z() *= scale;
    kurU *= scale;
    kurV *= scale;
    elasticity->ComputeStress(nU, nV, mU, mV, mechU, mechV, kurU, kurV, zInf, zSup, angle);
    const double vmBack = peakVonMises(nU, nV, mU, mV, thickness);
    if (vmBack > yield && vmBack > 0.0) {
        const double extra = yield / vmBack;
        nU *= extra;
        nV *= extra;
        mU *= extra;
        mV *= extra;
        scale *= extra;
    }

    ReissnerPlasticData updated;
    updated.epsU = epsUTrial - mechU;
    updated.epsV = epsVTrial - mechV;
    updated.kurU = kurUTrial - kurU;
    updated.kurV = kurVTrial - kurV;
    // Equivalent plastic strain from the geometric increment. Dividing the
    // trial stress by the cold modulus under-counts ductility once the
    // handbook has softened the metal, so a molten sheet never reached its
    // tear strain.
    const double dexx = (1.0 - scale) * (epsU.x() - thermal);
    const double deyy = (1.0 - scale) * (epsV.y() - thermal);
    const double dexy = (1.0 - scale) * 0.5 * (epsU.y() + epsV.x());
    const double membrane = std::sqrt(std::max(0.0, dexx * dexx - dexx * deyy + deyy * deyy + 3.0 * dexy * dexy));
    const double bkxx = (1.0 - scale) * kurY * halfH;
    const double bkyy = (1.0 - scale) * kurX * halfH;
    const double bend = std::sqrt(std::max(0.0, bkxx * bkxx - bkxx * bkyy + bkyy * bkyy));
    const double shearInc = (1.0 - scale) * trialShear;
    updated.p_strain_acc = old.p_strain_acc + std::max(membrane, std::max(bend, shearInc));
    return store(updated, true);
}

void ReissnerJ2::ComputeStiffnessMatrixElastoplastic(
    chrono::ChMatrixRef K,
    const chrono::ChVector3d& epsU,
    const chrono::ChVector3d& epsV,
    const chrono::ChVector3d& kurU,
    const chrono::ChVector3d& kurV,
    const chrono::fea::ChShellReissnerInternalData& data,
    double zInf,
    double zSup,
    double angle)
{
    (void)data;
    if (!section || !section->GetElasticity())
        return;
    // Chrono's numerical perfect-plastic tangent is singular on a molten yield
    // surface and the step aborts. The temperature-scaled elastic tangent
    // stays positive definite. The return mapping still caps the stress, so
    // the face cannot carry a real load.
    section->GetElasticity()->ComputeStiffnessMatrix(K, epsU, epsV, kurU, kurV, zInf, zSup, angle);
}

ThermalRayleighDamping::ThermalRayleighDamping(std::shared_ptr<ThermalReissnerElasticity> elasticity, double beta)
    : _elasticity(std::move(elasticity)), _beta(beta) {}

void ThermalRayleighDamping::refresh(double zInf, double zSup, double angle) {
    if (!_elasticity)
        return;
    const double scale = _elasticity->stiffnessScale();
    if (std::abs(scale - _cachedScale) <= 1.0e-6)
        return;
    const chrono::ChVector3d zero(0.0, 0.0, 0.0);
    _elasticity->ComputeStiffnessMatrix(_stiffness, zero, zero, zero, zero, zInf, zSup, angle);
    _cachedScale = scale;
}

void ThermalRayleighDamping::ComputeStress(chrono::ChVector3d& nU,
                                           chrono::ChVector3d& nV,
                                           chrono::ChVector3d& mU,
                                           chrono::ChVector3d& mV,
                                           const chrono::ChVector3d& depsU,
                                           const chrono::ChVector3d& depsV,
                                           const chrono::ChVector3d& dkurU,
                                           const chrono::ChVector3d& dkurV,
                                           double zInf,
                                           double zSup,
                                           double angle)
{
    refresh(zInf, zSup, angle);
    chrono::ChVectorN<double, 12> rate;
    rate.segment(0, 3) = depsU.eigen();
    rate.segment(3, 3) = depsV.eigen();
    rate.segment(6, 3) = dkurU.eigen();
    rate.segment(9, 3) = dkurV.eigen();
    const chrono::ChVectorN<double, 12> stress = _beta * _stiffness * rate;
    nU = stress.segment(0, 3);
    nV = stress.segment(3, 3);
    mU = stress.segment(6, 3);
    mV = stress.segment(9, 3);
}

void ThermalRayleighDamping::ComputeDampingMatrix(chrono::ChMatrixRef R,
                                                  const chrono::ChVector3d& depsU,
                                                  const chrono::ChVector3d& depsV,
                                                  const chrono::ChVector3d& dkurU,
                                                  const chrono::ChVector3d& dkurV,
                                                  double zInf,
                                                  double zSup,
                                                  double angle)
{
    (void)depsU;
    (void)depsV;
    (void)dkurU;
    (void)dkurV;
    refresh(zInf, zSup, angle);
    R = _beta * _stiffness;
}

void ReissnerJ2::CreatePlasticityData(
    int numpoints,
    std::vector<std::unique_ptr<chrono::fea::ChShellReissnerInternalData>>& plasticData)
{
    plasticData.resize(static_cast<std::size_t>(std::max(0, numpoints)));
    for (auto& slot : plasticData)
        slot = std::unique_ptr<chrono::fea::ChShellReissnerInternalData>(new ReissnerPlasticData());
}

}  // namespace rcs
