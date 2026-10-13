#include "physics_core/thermal_tet.hpp"
#include "physics_core/material.hpp"
#include "physics_core/physics_units.hpp"

#include "chrono/fea/ChMatrixCorotation.h"
#include "chrono/fea/ChPolarDecomposition.h"

#include <algorithm>
#include <cmath>

namespace rcs {
namespace {

using chrono::ChMatrix33;
using chrono::ChMatrixDynamic;
using chrono::ChVector3d;
using chrono::ChVectorDynamic;

void scaleHyperelasticTangent(ChMatrixDynamic<>& C, const ChMatrix33<>& stretch, double yield,
                              double rankineTension, double rankineCompression) {
    if (yield > 0.0 || rankineTension > 0.0 || rankineCompression > 0.0)
        return;
    const double j = std::max(stretch.determinant(), 1.0e-12);
    C *= 1.0 / std::max(std::cbrt(j), 1.0e-3);
}

}  // namespace

ThermalTet4::ThermalTet4() {
    _B.setZero(6, 12);
    _K.setZero(12, 12);
    _rotation.setIdentity();
    _stretch.setIdentity();
}

void ThermalTet4::SetNodes(std::shared_ptr<chrono::fea::ChNodeFEAxyz> a,
                           std::shared_ptr<chrono::fea::ChNodeFEAxyz> b,
                           std::shared_ptr<chrono::fea::ChNodeFEAxyz> c,
                           std::shared_ptr<chrono::fea::ChNodeFEAxyz> d) {
    _nodes[0] = std::move(a);
    _nodes[1] = std::move(b);
    _nodes[2] = std::move(c);
    _nodes[3] = std::move(d);
    std::vector<chrono::ChVariables*> vars;
    vars.reserve(4);
    for (int i = 0; i < 4; ++i)
        vars.push_back(&_nodes[i]->Variables());
    Kmatr.SetVariables(vars);
}

void ThermalTet4::configure(double youngs, double nu, double yield, double density, double beta) {
    _youngs = std::max(youngs, 1.0);
    _nu = std::clamp(nu, 0.0, 0.49);
    _yield = std::max(0.0, yield);
    _density = std::max(density, 1.0);
    _beta = std::max(0.0, beta);
}

void ThermalTet4::setThermalState(double thermalStrain, double stiffnessScale, double strengthScale) {
    _thermal = thermalStrain;
    _stiffness = std::clamp(stiffnessScale, kThermalFloor, 1.0);
    _strength = std::clamp(strengthScale, kThermalFloor, 1.0);
}

void ThermalTet4::setRankine(double tensionLimit, double compressionLimit) {
    _rankineTension = std::max(0.0, tensionLimit);
    _rankineCompression = std::max(0.0, compressionLimit);
}

void ThermalTet4::setFailureStrain(double strain) {
    _failPlastic = std::max(0.0, strain);
}

void ThermalTet4::setHardening(double ultimate) {
    _ultimate = std::max(0.0, ultimate);
}

double ThermalTet4::tangentScale() const {
    double scale = 1.0;
    if (_failPlastic > 0.0 && _plastic > 0.4 * _failPlastic) {
        const double span = std::max(0.6 * _failPlastic, 1.0e-12);
        const double u = std::clamp((_plastic - 0.4 * _failPlastic) / span, 0.0, 1.0);
        scale *= 1.0 - 0.95 * u;
    }
    // On the yield surface the elastic slope over-corrects a fast impact.
    // About 0.35 of the elastic matrix stays positive-definite; near the
    // solidus keep the full tangent so a molten yield surface does not go
    // singular. Rankine glass sheds with the same idea at half slope.
    if (_strength > 0.05 && _yield > 0.0 &&
        (_plastic > 1.0e-6 || _vonMises >= 0.99 * _yield * _strength))
        scale *= 0.35;
    if (_overstrength && _rankineTension > 0.0)
        scale *= 0.5;
    return scale;
}

void ThermalTet4::rebuildReference() {
    if (!_nodes[0])
        return;
    chrono::ChMatrixNM<double, 4, 4> rest;
    for (int i = 0; i < 4; ++i) {
        const ChVector3d x = _nodes[i]->GetX0();
        rest(0, i) = x.x();
        rest(1, i) = x.y();
        rest(2, i) = x.z();
        rest(3, i) = 1.0;
    }
    _shape = rest.inverse();

    // N = S * [x y z 1]. If that does not reproduce the nodes, S was built
    // transposed and the derivatives sit in the columns instead.
    bool rows = true;
    for (int i = 0; i < 4 && rows; ++i) {
        const ChVector3d x = _nodes[i]->GetX0();
        for (int j = 0; j < 4; ++j) {
            const double n = _shape(j, 0) * x.x() + _shape(j, 1) * x.y() + _shape(j, 2) * x.z() + _shape(j, 3);
            const double expect = i == j ? 1.0 : 0.0;
            if (std::abs(n - expect) > 1.0e-4)
                rows = false;
        }
    }
    ChVector3d grad[4];
    for (int i = 0; i < 4; ++i) {
        grad[i] = rows ? ChVector3d(_shape(i, 0), _shape(i, 1), _shape(i, 2))
                       : ChVector3d(_shape(0, i), _shape(1, i), _shape(2, i));
    }

    const ChVector3d e1 = _nodes[1]->GetX0() - _nodes[0]->GetX0();
    const ChVector3d e2 = _nodes[2]->GetX0() - _nodes[0]->GetX0();
    const ChVector3d e3 = _nodes[3]->GetX0() - _nodes[0]->GetX0();
    _volume = std::abs(e1.Dot(e2.Cross(e3))) / 6.0;
    if (_volume < 1.0e-16)
        _volume = 1.0e-16;

    _B.setZero(6, 12);
    for (int i = 0; i < 4; ++i) {
        const int c = 3 * i;
        _B(0, c + 0) = grad[i].x();
        _B(1, c + 1) = grad[i].y();
        _B(2, c + 2) = grad[i].z();
        _B(3, c + 0) = grad[i].y();
        _B(3, c + 1) = grad[i].x();
        _B(4, c + 1) = grad[i].z();
        _B(4, c + 2) = grad[i].y();
        _B(5, c + 0) = grad[i].z();
        _B(5, c + 2) = grad[i].x();
    }
    _ready = true;
}

void ThermalTet4::elasticLaw(ChMatrixDynamic<>& C) const {
    const double E = _youngs * _stiffness;
    const double nu = _nu;
    const double lambda = E * nu / ((1.0 + nu) * (1.0 - 2.0 * nu));
    const double mu = E / (2.0 * (1.0 + nu));
    C.setZero(6, 6);
    C(0, 0) = C(1, 1) = C(2, 2) = lambda + 2.0 * mu;
    C(0, 1) = C(0, 2) = C(1, 0) = C(1, 2) = C(2, 0) = C(2, 1) = lambda;
    C(3, 3) = C(4, 4) = C(5, 5) = mu;
}

void ThermalTet4::refreshRotation() {
    if (!_ready)
        return;
    chrono::ChMatrixNM<double, 4, 4> current;
    for (int i = 0; i < 4; ++i) {
        const ChVector3d x = _nodes[i]->GetPos();
        current(0, i) = x.x();
        current(1, i) = x.y();
        current(2, i) = x.z();
        current(3, i) = 1.0;
    }
    ChMatrix33<> F;
    for (int col = 0; col < 3; ++col) {
        for (int row = 0; row < 3; ++row) {
            double sum = 0.0;
            for (int k = 0; k < 4; ++k)
                sum += current(row, k) * _shape(k, col);
            F(row, col) = sum;
        }
    }
    const double detF = F.determinant();
    const double fNorm = F.norm();
    const double scale = std::max(1.0, fNorm * fNorm * fNorm);
    constexpr double kSingular = 1.0e-12;
    _collapsed = !std::isfinite(detF) || std::abs(detF) < kSingular * scale;
    _inverted = std::isfinite(detF) && detF < 0.0 && !_collapsed;
    // Chrono's polar decomposition printf's and writes Q = F when det == 0.
    // Keep the last good rotation instead of feeding it a planar F.
    if (_collapsed)
        return;
    const double det = chrono::fea::ChPolarDecomposition<>::Compute(F, _rotation, _stretch, 1.0e-6);
    if (det < 0.0)
        _rotation *= -1.0;
}

void ThermalTet4::GetStateBlock(ChVectorDynamic<>& state) {
    state.setZero(12);
    refreshRotation();
    for (int i = 0; i < 4; ++i) {
        const ChVector3d local = _rotation.transpose() * _nodes[i]->GetPos() - _nodes[i]->GetX0();
        state(3 * i + 0) = local.x();
        state(3 * i + 1) = local.y();
        state(3 * i + 2) = local.z();
    }
}

void ThermalTet4::deformationalSpeed(ChVectorDynamic<>& speed) const {
    speed.setZero(12);
    ChVector3d com(0, 0, 0);
    ChVector3d vcom(0, 0, 0);
    ChVector3d corner[4];
    ChVector3d velocity[4];
    for (int i = 0; i < 4; ++i) {
        if (!_nodes[i])
            return;
        corner[i] = _nodes[i]->GetPos();
        velocity[i] = _nodes[i]->GetPosDt();
        com += corner[i];
        vcom += velocity[i];
    }
    com *= 0.25;
    vcom *= 0.25;
    ChMatrix33<> inertia;
    inertia.setZero();
    ChVector3d angular(0, 0, 0);
    for (int i = 0; i < 4; ++i) {
        const ChVector3d r = corner[i] - com;
        const ChVector3d w = velocity[i] - vcom;
        angular += r.Cross(w);
        inertia(0, 0) += r.y() * r.y() + r.z() * r.z();
        inertia(1, 1) += r.x() * r.x() + r.z() * r.z();
        inertia(2, 2) += r.x() * r.x() + r.y() * r.y();
        inertia(0, 1) -= r.x() * r.y();
        inertia(0, 2) -= r.x() * r.z();
        inertia(1, 2) -= r.y() * r.z();
    }
    inertia(1, 0) = inertia(0, 1);
    inertia(2, 0) = inertia(0, 2);
    inertia(2, 1) = inertia(1, 2);
    ChVector3d omega(0, 0, 0);
    if (std::abs(inertia.determinant()) > 1.0e-16) {
        const ChMatrix33<> inverse(inertia.inverse());
        omega = inverse * angular;
    }
    for (int i = 0; i < 4; ++i) {
        const ChVector3d deform = velocity[i] - vcom - omega.Cross(corner[i] - com);
        const ChVector3d local = _rotation.transpose() * deform;
        speed(3 * i + 0) = local.x();
        speed(3 * i + 1) = local.y();
        speed(3 * i + 2) = local.z();
    }
}

void ThermalTet4::trialStress(ChVectorDynamic<>& stress, ChVectorDynamic<>& strain) {
    ChVectorDynamic<> local(12);
    GetStateBlock(local);
    const bool hyper = _yield <= 0.0 && _rankineTension <= 0.0 && _rankineCompression <= 0.0;
    if (hyper) {
        mooneyRivlinCauchy(_stretch, _youngs * _stiffness, _nu, _thermal, stress);
        henckyVoigt(_stretch, strain);
        strain(0) -= _thermal;
        strain(1) -= _thermal;
        strain(2) -= _thermal;
    } else {
        strain = _B * local;
        strain(0) -= _thermal;
        strain(1) -= _thermal;
        strain(2) -= _thermal;
    }
    for (int i = 0; i < 6; ++i) {
        strain(i) -= _plasticStrain[i];
        _plasticPending[i] = _plasticStrain[i];
    }
    if (!hyper) {
        ChMatrixDynamic<> C;
        elasticLaw(C);
        stress = C * strain;
    }

    const double vm = vonMisesStress(stress);
    _vonMises = std::isfinite(vm) ? vm : 0.0;
    const double ep0 = equivalentPlastic(_plasticStrain);
    const double flow0 = flowStress(_yield, _ultimate, _failPlastic, ep0, _strength);
    if (flow0 > 0.0 && std::isfinite(vm) && vm >= flow0) {
        const double E = _youngs * _stiffness;
        const double mu = E / (2.0 * (1.0 + _nu));
        const double twoMu = std::max(2.0 * mu, 1.0);
        const double threeMu = std::max(3.0 * mu, 1.0);
        const double mean = (stress(0) + stress(1) + stress(2)) / 3.0;
        double dep = 0.0;
        double flow = flow0;
        for (int iter = 0; iter < 12; ++iter) {
            const double ep = ep0 + dep;
            flow = flowStress(_yield, _ultimate, _failPlastic, ep, _strength);
            const double residual = vm - threeMu * dep - flow;
            if (std::abs(residual) <= 1.0e-6 * std::max(flow, vm))
                break;
            const double h = (flowStress(_yield, _ultimate, _failPlastic, ep + 1.0e-6, _strength) - flow) /
                             1.0e-6;
            dep = std::max(0.0, dep + residual / std::max(threeMu + h, 1.0));
        }
        flow = flowStress(_yield, _ultimate, _failPlastic, ep0 + dep, _strength);
        const double scale = vm > 0.0 ? std::clamp(flow / vm, 0.0, 1.0) : 1.0;
        for (int i = 0; i < 3; ++i) {
            const double dev = stress(i) - mean;
            const double kept = dev * scale;
            _plasticPending[i] = _plasticStrain[i] + (dev - kept) / twoMu;
            stress(i) = mean + kept;
        }
        for (int i = 3; i < 6; ++i) {
            const double kept = stress(i) * scale;
            _plasticPending[i] = _plasticStrain[i] + (stress(i) - kept) / std::max(mu, 1.0);
            stress(i) = kept;
        }
        _vonMises = flow;
        _plastic = equivalentPlastic(_plasticPending);
    } else {
        _plastic = ep0;
    }

    noteTension(stress);
    _rawTension = _tension;
    const double rawCompression = _rawCompression;
    _rawMean = (stress(0) + stress(1) + stress(2)) / 3.0;
    double strainPrincipal = 0.0;
    {
        chrono::ChVector3d v(1.0, 1.0, 1.0);
        v.Normalize();
        for (int iter = 0; iter < 8; ++iter) {
            const chrono::ChVector3d w(strain(0) * v.x() + 0.5 * strain(3) * v.y() + 0.5 * strain(5) * v.z(),
                                       0.5 * strain(3) * v.x() + strain(1) * v.y() + 0.5 * strain(4) * v.z(),
                                       0.5 * strain(5) * v.x() + 0.5 * strain(4) * v.y() + strain(2) * v.z());
            strainPrincipal = w.Dot(v);
            if (w.Length() < 1.0e-12)
                break;
            v = w / w.Length();
        }
    }
    _tensileStrain = std::isfinite(strainPrincipal) ? std::max(0.0, strainPrincipal) : 0.0;
    const double E = std::max(_youngs * _stiffness, 1.0);
    _trialEnergy = 0.5 * _rawTension * _rawTension / E * std::max(_volume, 0.0);

    _overstrength = projectRankine(stress, _rankineTension * _strength, _rankineCompression * _strength);
    if (_overstrength)
        _vonMises = vonMisesStress(stress);

    // Nearly-incompressible compliant solids lock under a shock. Cap the
    // hydrostatic part so lambda cannot send the residual to GPa.
    if (_yield <= 0.0 && _rankineTension <= 0.0 && _rankineCompression <= 0.0) {
        const double mean = (stress(0) + stress(1) + stress(2)) / 3.0;
        const double bulkLimit = std::max(_youngs * _stiffness, 1.0) * 3.0;
        if (mean < -bulkLimit) {
            const double shift = -bulkLimit - mean;
            stress(0) += shift;
            stress(1) += shift;
            stress(2) += shift;
            _overstrength = true;
        }
    }
    // Uniform elongation is near A. A load-controlled UTS hold equilibrates
    // around 0.05 if the coupon never necks; begin softening at 0.4 A so
    // constant UTS force can pull through 0.5 A without unloading at 0.03.
    if (_failPlastic > 0.0 && _plastic > 0.4 * _failPlastic) {
        const double span = std::max(0.6 * _failPlastic, 1.0e-12);
        const double u = std::clamp((_plastic - 0.4 * _failPlastic) / span, 0.0, 1.0);
        stress *= 1.0 - 0.95 * u;
    }
    noteTension(stress);
    // Fracture compares trial principals. Rankine (and the bulk cap) only
    // change the force; they must not rewrite the pre-projection compression
    // as ~0 and make a crush look tension-dominated.
    _rawCompression = rawCompression;
}

void ThermalTet4::noteTension(const ChVectorDynamic<>& stress) {
    auto principal = [&](double sign) {
        ChVector3d v(1.0, 1.0, 1.0);
        v.Normalize();
        double value = 0.0;
        for (int iter = 0; iter < 8; ++iter) {
            const ChVector3d w(sign * (stress(0) * v.x() + stress(3) * v.y() + stress(5) * v.z()),
                               sign * (stress(3) * v.x() + stress(1) * v.y() + stress(4) * v.z()),
                               sign * (stress(5) * v.x() + stress(4) * v.y() + stress(2) * v.z()));
            value = sign * w.Dot(v);
            if (w.Length() < 1.0e-12)
                break;
            v = w / w.Length();
        }
        return value;
    };
    const double tension = principal(1.0);
    const double compression = principal(-1.0);
    _tension = std::isfinite(tension) ? std::max(0.0, tension) : 0.0;
    _rawCompression = std::isfinite(compression) ? std::max(0.0, -compression) : 0.0;
    ChVector3d v(1.0, 1.0, 1.0);
    v.Normalize();
    for (int iter = 0; iter < 8; ++iter) {
        const ChVector3d w(stress(0) * v.x() + stress(3) * v.y() + stress(5) * v.z(),
                           stress(3) * v.x() + stress(1) * v.y() + stress(4) * v.z(),
                           stress(5) * v.x() + stress(4) * v.y() + stress(2) * v.z());
        if (w.Length() < 1.0e-12)
            break;
        v = w / w.Length();
    }
    _tensionAxis = _rotation * v;
    if (_tensionAxis.Length() > 1.0e-12)
        _tensionAxis.Normalize();
}

void ThermalTet4::commitPlastic() {
    for (int i = 0; i < 6; ++i)
        _plasticStrain[i] = _plasticPending[i];
    _plastic = equivalentPlastic(_plasticStrain);
}

void ThermalTet4::commitDamage(double value) {
    _damage = std::clamp(value, 0.0, 1.0);
}

void ThermalTet4::discardTrial() {
    for (int i = 0; i < 6; ++i)
        _plasticPending[i] = _plasticStrain[i];
    _plastic = equivalentPlastic(_plasticStrain);
    _overstrength = false;
    _rawTension = 0.0;
    _rawCompression = 0.0;
    _rawMean = 0.0;
    _trialEnergy = 0.0;
    _tensileStrain = 0.0;
}

void ThermalTet4::copyPlasticFrom(const ThermalTet4& other) {
    for (int i = 0; i < 6; ++i) {
        _plasticStrain[i] = other._plasticStrain[i];
        _plasticPending[i] = other._plasticStrain[i];
    }
    _plastic = other._plastic;
    _damage = other._damage;
}

void ThermalTet4::ComputeKRMmatricesGlobal(chrono::ChMatrixRef H, double Kfactor, double Rfactor, double Mfactor) {
    if (!_ready)
        rebuildReference();
    refreshRotation();
    ChMatrixDynamic<> C;
    elasticLaw(C);
    scaleHyperelasticTangent(C, _stretch, _yield, _rankineTension, _rankineCompression);
    C *= tangentScale();
    _K = _volume * (_B.transpose() * C * _B);
    for (int row = 0; row < 12; ++row)
        for (int col = row + 1; col < 12; ++col)
            _K(row, col) = _K(col, row);

    ChMatrixDynamic<> CK(12, 12);
    ChMatrixDynamic<> CKCt(12, 12);
    chrono::fea::ChMatrixCorotation::ComputeCK(_K, _rotation, 4, CK);
    chrono::fea::ChMatrixCorotation::ComputeKCt(CK, _rotation, 4, CKCt);
    H = (Kfactor + Rfactor * _beta) * CKCt;
    if (Mfactor != 0.0) {
        const double lump = _volume * _density / 4.0;
        for (int id = 0; id < 12; ++id)
            H(id, id) += Mfactor * lump;
    }
}

void ThermalTet4::ComputeInternalForces(ChVectorDynamic<>& Fi) {
    if (!_ready)
        rebuildReference();
    refreshRotation();
    if (_collapsed || _inverted) {
        Fi.setZero(12);
        return;
    }
    ChVectorDynamic<> stress(6);
    ChVectorDynamic<> strain(6);
    trialStress(stress, strain);
    ChMatrixDynamic<> C;
    elasticLaw(C);
    scaleHyperelasticTangent(C, _stretch, _yield, _rankineTension, _rankineCompression);
    C *= tangentScale();
    _K = _volume * (_B.transpose() * C * _B);

    // Linear B is conjugate to first PK on the reference tet. Cauchy in that
    // slot overstates axial force by ~λ, so a 0.9 E pull stalls short of tear.
    if (_yield <= 0.0 && _rankineTension <= 0.0 && _rankineCompression <= 0.0 &&
        std::abs(_stretch.determinant()) > 1.0e-12) {
        ChMatrix33<> sigma;
        sigma.setZero();
        sigma(0, 0) = stress(0);
        sigma(1, 1) = stress(1);
        sigma(2, 2) = stress(2);
        sigma(0, 1) = sigma(1, 0) = stress(3);
        sigma(1, 2) = sigma(2, 1) = stress(4);
        sigma(0, 2) = sigma(2, 0) = stress(5);
        const double J = std::max(_stretch.determinant(), 1.0e-12);
        const ChMatrix33<> pk = (J * sigma) * _stretch.inverse().transpose();
        stress(0) = pk(0, 0);
        stress(1) = pk(1, 1);
        stress(2) = pk(2, 2);
        stress(3) = 0.5 * (pk(0, 1) + pk(1, 0));
        stress(4) = 0.5 * (pk(1, 2) + pk(2, 1));
        stress(5) = 0.5 * (pk(0, 2) + pk(2, 0));
    }

    ChVectorDynamic<> localForce = _B.transpose() * stress * _volume;

    ChVectorDynamic<> speed(12);
    deformationalSpeed(speed);
    localForce += _beta * (_K * speed);
    if (_yield <= 0.0 && _rankineTension <= 0.0) {
        const ChVectorDynamic<> rate = _B * speed;
        const double tr = rate(0) + rate(1) + rate(2);
        if (tr < 0.0) {
            const double length = std::cbrt(std::max(6.0 * _volume, 1.0e-18));
            const double wave = std::sqrt(std::max(_youngs * _stiffness, 1.0) / _density);
            const double q = _density * wave * length * (-tr);
            ChVectorDynamic<> bulk(6);
            bulk.setZero(6);
            bulk(0) = bulk(1) = bulk(2) = q;
            localForce += _B.transpose() * bulk * _volume;
        }
    }
    localForce *= -1.0;
    chrono::fea::ChMatrixCorotation::ComputeCK(localForce, _rotation, 4, Fi);
}

void ThermalTet4::SetupInitial(chrono::ChSystem*) {
    rebuildReference();
    _rotation.setIdentity();
}

void ThermalTet4::Update() {
    chrono::fea::ChElementGeneric::Update();
    refreshRotation();
}

}  // namespace rcs
