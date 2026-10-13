#include "physics_core/contact_load.hpp"
#include "physics_core/physics_units.hpp"

#include <algorithm>
#include <cmath>

namespace rcs {

using chrono::ChVector3d;

namespace {

void addOuter(chrono::ChMatrixDynamic<double>& m, const ChVector3d& a, const ChVector3d& b, double s) {
    const double ax = a.x(), ay = a.y(), az = a.z();
    const double bx = b.x(), by = b.y(), bz = b.z();
    m(0, 0) += s * ax * bx;
    m(0, 1) += s * ax * by;
    m(0, 2) += s * ax * bz;
    m(1, 0) += s * ay * bx;
    m(1, 1) += s * ay * by;
    m(1, 2) += s * ay * bz;
    m(2, 0) += s * az * bx;
    m(2, 1) += s * az * by;
    m(2, 2) += s * az * bz;
}

}  // namespace

NodeContact::NodeContact(std::shared_ptr<chrono::fea::ChNodeFEAxyz> node,
                         std::shared_ptr<const ContactWorld> world)
    : chrono::ChLoadCustom(node), _node(std::move(node)), _world(std::move(world)) {}

NodeContact::NodeContact(std::shared_ptr<chrono::fea::ChNodeFEAxyzrot> node,
                         std::shared_ptr<const ContactWorld> world)
    : chrono::ChLoadCustom(node), _rotNode(std::move(node)), _world(std::move(world)) {}

void NodeContact::configure(double offset, double stiffness, double dissipation, double mass) {
    (void)dissipation;
    _offset = offset;
    _mass = std::max(0.0, mass);
    _stiffness = std::max(0.0, stiffness);
    _damping = 0.0;
    if (_mass > 0.0 && _stiffness > 0.0)
        _damping = kContactDampingFraction * 2.0 * std::sqrt(_stiffness * _mass);
}

void NodeContact::addPlane(const ChVector3d& position,
                           const ChVector3d& velocity,
                           const ChVector3d& normal,
                           double planeOffset,
                           int slot,
                           ChVector3d& force,
                           bool jac)
{
    const double gap = planeOffset + _offset - position.Dot(normal);
    // A 10 µm open band keeps the penalty in the residual when Newton
    // overshoots the plane. Dropping k to zero mid-iteration is why a
    // second, slower landing would not converge and then froze.
    if (gap < -1.0e-5)
        return;
    const double vn = velocity.Dot(normal);
    const double vStick = _world ? std::max(_world->stickSpeed, 1.0e-6) : 1.0e-3;
    const double c = contactDashpot(_stiffness, _mass, vn, vStick);
    // Continuum penalty: F_n = k δ − c v_n. Stabilization is inside c and
    // is zero when the node is not closing. Unilateral: no tensile pull.
    double fn = _stiffness * gap - c * vn;
    if (fn <= 0.0)
        return;
    _normalForce = std::max(_normalForce, fn);
    force += normal * fn;
    const double mu = _world ? std::max(0.0, _world->friction) : 0.0;
    const ChVector3d vt = velocity - normal * vn;
    const double speed = vt.Length();
    const double ramp = std::min(1.0, speed / vStick);
    ChVector3d tHat(0, 0, 0);
    if (mu > 0.0 && speed > 1.0e-12) {
        tHat = vt / speed;
        force -= tHat * (mu * fn * ramp);
    }
    if (!jac && slot >= 0 && slot < 6) {
        _planeForce[slot] = fn;
        _planeNormal[slot] = normal;
    }
    if (!jac || !m_jacobians)
        return;
    // Chrono K = −dQ/dx, R = −dQ/dv. F = k δ − c v_n, δ decreasing along n.
    addOuter(m_jacobians->K, normal, normal, _stiffness);
    if (c > 0.0)
        addOuter(m_jacobians->R, normal, normal, c);
    if (mu <= 0.0 || fn <= 0.0)
        return;
    if (speed < vStick) {
        const double rStick = mu * fn / vStick;
        m_jacobians->R(0, 0) += rStick;
        m_jacobians->R(1, 1) += rStick;
        m_jacobians->R(2, 2) += rStick;
        addOuter(m_jacobians->R, normal, normal, -rStick);
    } else if (speed > 1.0e-12) {
        const double ft = mu * fn;
        const double s = ft / speed;
        m_jacobians->R(0, 0) += s;
        m_jacobians->R(1, 1) += s;
        m_jacobians->R(2, 2) += s;
        addOuter(m_jacobians->R, tHat, tHat, -s);
        addOuter(m_jacobians->R, normal, normal, -s);
        if (c > 0.0)
            addOuter(m_jacobians->R, tHat, normal, mu * c);
    }
}

void NodeContact::ComputeQ(chrono::ChState* state_x, chrono::ChStateDelta* state_w) {
    const int ndof = _rotNode ? 6 : 3;
    if (load_Q.size() != ndof)
        load_Q.setZero(ndof);
    load_Q.setZero();
    _normalForce = 0.0;
    if (!_fillJac) {
        for (double& fn : _planeForce)
            fn = 0.0;
    }
    if (!_active || !_world || _stiffness <= 0.0)
        return;

    const ChVector3d position = state_x ? ChVector3d((*state_x)(0), (*state_x)(1), (*state_x)(2))
                                        : (_node     ? _node->GetPos()
                                           : _rotNode ? _rotNode->GetPos()
                                                      : ChVector3d(0, 0, 0));
    const ChVector3d velocity = state_w ? ChVector3d((*state_w)(0), (*state_w)(1), (*state_w)(2))
                                        : (_node     ? _node->GetPosDt()
                                           : _rotNode ? _rotNode->GetPosDt()
                                                      : ChVector3d(0, 0, 0));

    ChVector3d force(0, 0, 0);
    const bool jac = _fillJac;
    if (_world->floorEnabled)
        addPlane(position, velocity, ChVector3d(0, 1, 0), _world->floorHeight, 0, force, jac);
    addPlane(position, velocity, ChVector3d(1, 0, 0), -_world->halfX, 1, force, jac);
    addPlane(position, velocity, ChVector3d(-1, 0, 0), -_world->halfX, 2, force, jac);
    addPlane(position, velocity, ChVector3d(0, 0, 1), -_world->halfZ, 3, force, jac);
    addPlane(position, velocity, ChVector3d(0, 0, -1), -_world->halfZ, 4, force, jac);
    addPlane(position, velocity, ChVector3d(0, -1, 0), -_world->ceiling, 5, force, jac);

    load_Q(0) = force.x();
    load_Q(1) = force.y();
    load_Q(2) = force.z();
}

void NodeContact::ComputeJacobian(chrono::ChState* state_x, chrono::ChStateDelta* state_w) {
    if (!m_jacobians)
        CreateJacobianMatrices();
    m_jacobians->K.setZero();
    m_jacobians->R.setZero();
    m_jacobians->M.setZero();
    _fillJac = true;
    ComputeQ(state_x, state_w);
    _fillJac = false;
}

}  // namespace rcs
