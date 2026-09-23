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
    return std::max(top, bottom);
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
    const double vm = peakVonMises(nU, nV, mU, mV, thickness);

    auto store = [&](const ReissnerPlasticData& state, bool yielded) {
        if (newState)
            *newState = state;
        else
            dataNew.p_strain_acc = state.p_strain_acc;
        return yielded;
    };

    if (!(vm > _yield) || _yield <= 0.0) {
        ReissnerPlasticData kept = old;
        kept.p_strain_acc = old.p_strain_acc;
        return store(kept, false);
    }

    const double scale = std::clamp(_yield / vm, 0.0, 1.0);
    epsU *= scale;
    epsV *= scale;
    kurU *= scale;
    kurV *= scale;
    elasticity->ComputeStress(nU, nV, mU, mV, epsU, epsV, kurU, kurV, zInf, zSup, angle);

    ReissnerPlasticData updated;
    updated.epsU = epsUTrial - epsU;
    updated.epsV = epsVTrial - epsV;
    updated.kurU = kurUTrial - kurU;
    updated.kurV = kurVTrial - kurV;
    updated.p_strain_acc = old.p_strain_acc + (1.0 - scale) * vm / _young;
    return store(updated, true);
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
