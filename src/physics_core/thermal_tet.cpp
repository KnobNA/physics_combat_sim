#include "physics_core/thermal_tet.hpp"
#include "physics_core/material.hpp"

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

double vonMisesStress(const ChVectorDynamic<>& s) {
    const double mean = (s(0) + s(1) + s(2)) / 3.0;
    const double d0 = s(0) - mean;
    const double d1 = s(1) - mean;
    const double d2 = s(2) - mean;
    return std::sqrt(1.5 * (d0 * d0 + d1 * d1 + d2 * d2 + 2.0 * (s(3) * s(3) + s(4) * s(4) + s(5) * s(5))));
}

double equivalentPlastic(const double* strain) {
    double dev = strain[0] * strain[0] + strain[1] * strain[1] + strain[2] * strain[2];
    dev += 0.5 * (strain[3] * strain[3] + strain[4] * strain[4] + strain[5] * strain[5]);
    return std::sqrt(std::max(0.0, 2.0 / 3.0 * dev));
}

/// Clamp principal stresses onto the Rankine rectangle. Only the principals
/// that are actually over the limit move; the others keep carrying load. A
/// uniform scale of the whole tensor is what flattened tets on impact.
bool projectRankine(ChVectorDynamic<>& stress, double tension, double compression) {
    if (tension <= 0.0 && compression <= 0.0)
        return false;

    double A[3][3] = {{stress(0), stress(3), stress(5)},
                      {stress(3), stress(1), stress(4)},
                      {stress(5), stress(4), stress(2)}};
    double V[3][3] = {{1.0, 0.0, 0.0}, {0.0, 1.0, 0.0}, {0.0, 0.0, 1.0}};
    for (int iter = 0; iter < 16; ++iter) {
        int p = 0;
        int q = 1;
        double maxOff = std::abs(A[0][1]);
        if (std::abs(A[0][2]) > maxOff) {
            maxOff = std::abs(A[0][2]);
            p = 0;
            q = 2;
        }
        if (std::abs(A[1][2]) > maxOff) {
            maxOff = std::abs(A[1][2]);
            p = 1;
            q = 2;
        }
        const double scale = 1.0 + std::abs(A[0][0]) + std::abs(A[1][1]) + std::abs(A[2][2]);
        if (maxOff < 1.0e-12 * scale)
            break;

        const double app = A[p][p];
        const double aqq = A[q][q];
        const double apq = A[p][q];
        const double tau = (aqq - app) / (2.0 * apq);
        const double t = (tau >= 0.0 ? 1.0 : -1.0) / (std::abs(tau) + std::sqrt(1.0 + tau * tau));
        const double c = 1.0 / std::sqrt(1.0 + t * t);
        const double s = t * c;

        const double appNew = c * c * app - 2.0 * s * c * apq + s * s * aqq;
        const double aqqNew = s * s * app + 2.0 * s * c * apq + c * c * aqq;
        A[p][p] = appNew;
        A[q][q] = aqqNew;
        A[p][q] = A[q][p] = 0.0;
        for (int k = 0; k < 3; ++k) {
            if (k == p || k == q)
                continue;
            const double apk = A[p][k];
            const double aqk = A[q][k];
            A[p][k] = A[k][p] = c * apk - s * aqk;
            A[q][k] = A[k][q] = s * apk + c * aqk;
        }
        for (int k = 0; k < 3; ++k) {
            const double vip = V[k][p];
            const double viq = V[k][q];
            V[k][p] = c * vip - s * viq;
            V[k][q] = s * vip + c * viq;
        }
    }

    const double lo = compression > 0.0 ? -compression : -1.0e30;
    const double hi = tension > 0.0 ? tension : 1.0e30;
    double lambda[3] = {A[0][0], A[1][1], A[2][2]};
    bool clipped = false;
    for (double& value : lambda) {
        const double kept = std::clamp(value, lo, hi);
        clipped = clipped || kept != value;
        value = kept;
    }
    if (!clipped)
        return false;

    double R[3][3] = {};
    for (int i = 0; i < 3; ++i)
        for (int j = 0; j < 3; ++j)
            for (int k = 0; k < 3; ++k)
                R[i][j] += V[i][k] * lambda[k] * V[j][k];
    stress(0) = R[0][0];
    stress(1) = R[1][1];
    stress(2) = R[2][2];
    stress(3) = 0.5 * (R[0][1] + R[1][0]);
    stress(4) = 0.5 * (R[1][2] + R[2][1]);
    stress(5) = 0.5 * (R[0][2] + R[2][0]);
    return true;
}

}  // namespace

ThermalTet4::ThermalTet4() {
    _B.setZero(6, 12);
    _K.setZero(12, 12);
    _rotation.setIdentity();
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
    const double nu = (_yield <= 0.0 && _rankineTension <= 0.0) ? std::min(_nu, 0.45) : _nu;
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
    ChMatrix33<> stretch;
    const double det = chrono::fea::ChPolarDecomposition<>::Compute(F, _rotation, stretch, 1.0e-6);
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

void ThermalTet4::trialStress(ChVectorDynamic<>& stress, ChVectorDynamic<>& strain) {
    ChVectorDynamic<> local(12);
    GetStateBlock(local);
    strain = _B * local;
    strain(0) -= _thermal;
    strain(1) -= _thermal;
    strain(2) -= _thermal;
    for (int i = 0; i < 6; ++i) {
        strain(i) -= _plasticStrain[i];
        _plasticPending[i] = _plasticStrain[i];
    }
    ChMatrixDynamic<> C;
    elasticLaw(C);
    stress = C * strain;

    const double vm = vonMisesStress(stress);
    _vonMises = std::isfinite(vm) ? vm : 0.0;
    const double yield = _yield * _strength;
    if (yield > 0.0 && vm > yield && std::isfinite(vm)) {
        const double E = _youngs * _stiffness;
        const double mu = E / (2.0 * (1.0 + _nu));
        const double twoMu = std::max(2.0 * mu, 1.0);
        const double mean = (stress(0) + stress(1) + stress(2)) / 3.0;
        const double scale = yield / vm;
        for (int i = 0; i < 3; ++i) {
            const double dev = stress(i) - mean;
            const double kept = dev * scale;
            _plasticPending[i] += (dev - kept) / twoMu;
            stress(i) = mean + kept;
        }
        for (int i = 3; i < 6; ++i) {
            const double kept = stress(i) * scale;
            _plasticPending[i] += (stress(i) - kept) / std::max(mu, 1.0);
            stress(i) = kept;
        }
        _vonMises = yield;
        _plastic = equivalentPlastic(_plasticPending);
    } else {
        _plastic = equivalentPlastic(_plasticStrain);
    }

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
    noteTension(stress);
}

void ThermalTet4::noteTension(const ChVectorDynamic<>& stress) {
    ChVector3d v(1.0, 1.0, 1.0);
    v.Normalize();
    double tension = 0.0;
    for (int iter = 0; iter < 8; ++iter) {
        const ChVector3d w(stress(0) * v.x() + stress(3) * v.y() + stress(5) * v.z(),
                           stress(3) * v.x() + stress(1) * v.y() + stress(4) * v.z(),
                           stress(5) * v.x() + stress(4) * v.y() + stress(2) * v.z());
        tension = w.Dot(v);
        if (w.Length() < 1.0e-12)
            break;
        v = w / w.Length();
    }
    _tension = std::isfinite(tension) ? std::max(0.0, tension) : 0.0;
    _tensionAxis = _rotation * v;
    if (_tensionAxis.Length() > 1.0e-12)
        _tensionAxis.Normalize();
}

void ThermalTet4::commitPlastic() {
    for (int i = 0; i < 6; ++i)
        _plasticStrain[i] = _plasticPending[i];
    _plastic = equivalentPlastic(_plasticStrain);
}

void ThermalTet4::copyPlasticFrom(const ThermalTet4& other) {
    for (int i = 0; i < 6; ++i) {
        _plasticStrain[i] = other._plasticStrain[i];
        _plasticPending[i] = other._plasticStrain[i];
    }
    _plastic = other._plastic;
}

void ThermalTet4::ComputeKRMmatricesGlobal(chrono::ChMatrixRef H, double Kfactor, double Rfactor, double Mfactor) {
    if (!_ready)
        rebuildReference();
    refreshRotation();
    ChMatrixDynamic<> C;
    elasticLaw(C);
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
    ChVectorDynamic<> stress(6);
    ChVectorDynamic<> strain(6);
    trialStress(stress, strain);
    ChMatrixDynamic<> C;
    elasticLaw(C);
    _K = _volume * (_B.transpose() * C * _B);

    ChVectorDynamic<> localForce = _B.transpose() * stress * _volume;

    ChVector3d com(0, 0, 0);
    ChVector3d vcom(0, 0, 0);
    ChVector3d corner[4];
    ChVector3d velocity[4];
    for (int i = 0; i < 4; ++i) {
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

    ChVectorDynamic<> speed(12);
    for (int i = 0; i < 4; ++i) {
        const ChVector3d deform = velocity[i] - vcom - omega.Cross(corner[i] - com);
        const ChVector3d local = _rotation.transpose() * deform;
        speed(3 * i + 0) = local.x();
        speed(3 * i + 1) = local.y();
        speed(3 * i + 2) = local.z();
    }
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
