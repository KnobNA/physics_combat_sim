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
    // Force acts only on penetration. A dashpot reaching into an open gap
    // pushes a node that has not touched and makes the force jump at the
    // band edge, which Newton cannot follow.
    if (gap <= 0.0 || !std::isfinite(gap))
        return;
    const double gapLaw = std::min(gap, kContactGapCap);
    const double vn = velocity.Dot(normal);
    const double vStick = _world ? std::max(_world->stickSpeed, 1.0e-6) : 1.0e-3;
    // Softened contact: force and tangent are continuous at touch, so the
    // many nodes of a seated face do not chatter in and out of the active
    // set between Newton iterates.
    const double soft = std::max(kContactSoftening, 1.0e-9);
    const double depthScale = std::min(1.0, gapLaw / soft);
    const double spring = gapLaw < soft ? 0.5 * _stiffness * gapLaw * gapLaw / soft
                                        : _stiffness * (gapLaw - 0.5 * soft);
    // Hunt-Crossley hysteresis damping. It vanishes at touch, needs no node
    // mass, and returns a struck mass at the surface restitution. v⁻ is the
    // approach speed committed before this node touched.
    const double D = huntCrossleyFactor(_world ? _world->restitution : 0.15,
                                        std::max(0.0, -_approach.Dot(normal)));
    const HuntCrossleyScale hc =
        huntCrossleyScale(D, vn, std::max(0.0, -_approach.Dot(normal)));
    const double c = hc.damping * spring;
    const double fn = spring * hc.scale;
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
    // Chrono K = −dQ/dx, R = −dQ/dv. F = F_s(δ) s(v_n), δ decreasing along n.
    const double kn = _stiffness * depthScale * hc.scale;
    addOuter(m_jacobians->K, normal, normal, kn);
    if (c > 0.0)
        addOuter(m_jacobians->R, normal, normal, c);
    if (mu <= 0.0 || fn <= 0.0)
        return;
    // Friction scales with fn, so it inherits fn's gap and closing-speed
    // derivatives along the slip direction φ = vt / max(|vt|, vStick).
    const ChVector3d phi = vt / std::max(speed, vStick);
    addOuter(m_jacobians->K, phi, normal, -mu * kn);
    addOuter(m_jacobians->R, phi, normal, -mu * c);
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

void NodeContact::commitApproach() {
    if (!_world || (!_node && !_rotNode))
        return;
    const ChVector3d p = _node ? _node->GetPos() : _rotNode->GetPos();
    const ChVector3d v = _node ? _node->GetPosDt() : _rotNode->GetPosDt();
    auto open = [&](const ChVector3d& n, double plane) { return plane + _offset - p.Dot(n) <= 0.0; };
    const bool clear = (!_world->floorEnabled || open(ChVector3d(0, 1, 0), _world->floorHeight)) &&
                       open(ChVector3d(1, 0, 0), -_world->halfX) && open(ChVector3d(-1, 0, 0), -_world->halfX) &&
                       open(ChVector3d(0, 0, 1), -_world->halfZ) && open(ChVector3d(0, 0, -1), -_world->halfZ) &&
                       open(ChVector3d(0, -1, 0), -_world->ceiling);
    if (clear)
        _approach = v;
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

namespace {

void addBlock(chrono::ChMatrixDynamic<double>& m, int r, int c, const ChVector3d& a, const ChVector3d& b, double s) {
    for (int i = 0; i < 3; ++i)
        for (int j = 0; j < 3; ++j)
            m(r + i, c + j) += s * a[i] * b[j];
}

void addIdentity(chrono::ChMatrixDynamic<double>& m, int r, int c, double s) {
    for (int i = 0; i < 3; ++i)
        m(r + i, c + i) += s;
}

}  // namespace

NodeFaceContact::NodeFaceContact(std::vector<std::shared_ptr<chrono::ChLoadable>>& loadables,
                                 std::vector<Corner> corners)
    : chrono::ChLoadCustomMultiple(loadables), _corners(std::move(corners)) {
    int x = 0;
    int w = 0;
    for (const Corner& c : _corners) {
        const bool fixed = c.xyz ? c.xyz->IsFixed() : (!c.rot || c.rot->IsFixed());
        _xOffset.push_back(fixed ? -1 : x);
        _wOffset.push_back(fixed ? -1 : w);
        if (fixed)
            continue;
        x += c.xyz ? 3 : 7;
        w += c.xyz ? 3 : 6;
    }
}

void NodeFaceContact::configure(const ChVector3d& normal,
                                const std::vector<double>& shares,
                                double reach,
                                double stiffness,
                                double friction,
                                double stickSpeed,
                                double restitution) {
    _normal = normal;
    for (std::size_t c = 0; c < _corners.size() && c < shares.size(); ++c)
        _corners[c].share = shares[c];
    _reach = reach;
    _stiffness = std::max(0.0, stiffness);
    _friction = std::max(0.0, friction);
    _stickSpeed = std::max(stickSpeed, 1.0e-6);
    _restitution = restitution;
    if (!_latched) {
        ChVector3d v(0, 0, 0);
        for (std::size_t c = 0; c < _corners.size(); ++c)
            v += cornerVel(c, nullptr) * _corners[c].share;
        _approach = std::max(0.0, -v.Dot(_normal));
        _latched = true;
    }
}

ChVector3d NodeFaceContact::cornerPos(std::size_t c, const chrono::ChState* state_x) const {
    const int o = _xOffset[c];
    if (state_x && o >= 0)
        return ChVector3d((*state_x)(o), (*state_x)(o + 1), (*state_x)(o + 2));
    return _corners[c].xyz ? _corners[c].xyz->GetPos() : _corners[c].rot->GetPos();
}

ChVector3d NodeFaceContact::cornerVel(std::size_t c, const chrono::ChStateDelta* state_w) const {
    const int o = _wOffset[c];
    if (o < 0)
        return ChVector3d(0, 0, 0);
    if (state_w)
        return ChVector3d((*state_w)(o), (*state_w)(o + 1), (*state_w)(o + 2));
    return _corners[c].xyz ? _corners[c].xyz->GetPosDt() : _corners[c].rot->GetPosDt();
}

double NodeFaceContact::gap() const {
    ChVector3d d(0, 0, 0);
    for (std::size_t c = 0; c < _corners.size(); ++c)
        d += cornerPos(c, nullptr) * _corners[c].share;
    return _reach - d.Dot(_normal);
}

double NodeFaceContact::elasticEnergy() const {
    const double g = gap();
    if (g <= 0.0 || _stiffness <= 0.0)
        return 0.0;
    const double soft = kContactSoftening;
    return g < soft ? _stiffness * g * g * g / (6.0 * soft)
                    : _stiffness * (soft * soft / 6.0 + 0.5 * (g - soft) * (g - soft) + 0.5 * soft * (g - soft));
}

void NodeFaceContact::commitApproach() {
    if (gap() > 0.0)
        return;
    ChVector3d v(0, 0, 0);
    for (std::size_t c = 0; c < _corners.size(); ++c)
        v += cornerVel(c, nullptr) * _corners[c].share;
    _approach = std::max(0.0, -v.Dot(_normal));
}

void NodeFaceContact::evaluate(const chrono::ChState* state_x, const chrono::ChStateDelta* state_w, bool jac) {
    const int ndof = LoadGetNumCoordsVelLevel();
    if (load_Q.size() != ndof)
        load_Q.setZero(ndof);
    load_Q.setZero();
    if (_stiffness <= 0.0)
        return;
    ChVector3d d(0, 0, 0);
    ChVector3d vrel(0, 0, 0);
    for (std::size_t c = 0; c < _corners.size(); ++c) {
        d += cornerPos(c, state_x) * _corners[c].share;
        vrel += cornerVel(c, state_w) * _corners[c].share;
    }
    const ChVector3d& n = _normal;
    const double g = _reach - d.Dot(n);
    if (g <= 0.0 || !std::isfinite(g))
        return;
    const double gapLaw = std::min(g, kContactGapCap);
    const double vn = vrel.Dot(n);
    const double soft = std::max(kContactSoftening, 1.0e-9);
    const double depthScale = std::min(1.0, gapLaw / soft);
    const double spring = gapLaw < soft ? 0.5 * _stiffness * gapLaw * gapLaw / soft
                                        : _stiffness * (gapLaw - 0.5 * soft);
    const HuntCrossleyScale hc =
        huntCrossleyScale(huntCrossleyFactor(_restitution, _approach), vn, _approach);
    const double fn = spring * hc.scale;
    if (fn <= 0.0)
        return;
    const ChVector3d vt = vrel - n * vn;
    const double speed = vt.Length();
    ChVector3d force = n * fn;
    ChVector3d tHat(0, 0, 0);
    if (_friction > 0.0 && speed > 1.0e-12) {
        tHat = vt / speed;
        force -= tHat * (_friction * fn * std::min(1.0, speed / _stickSpeed));
    }
    for (std::size_t c = 0; c < _corners.size(); ++c) {
        const int o = _wOffset[c];
        if (o < 0)
            continue;
        const ChVector3d f = force * _corners[c].share;
        load_Q(o) = f.x();
        load_Q(o + 1) = f.y();
        load_Q(o + 2) = f.z();
    }
    if (!jac || !m_jacobians)
        return;
    // K = −dQ/dx and R = −dQ/dv. Every corner sees the same relative motion
    // scaled by its share, so each block is s_i s_j times the point-contact
    // block of the arena contact.
    const double kn = _stiffness * depthScale * hc.scale;
    const double cn = hc.damping * spring;
    double rStick = 0.0;
    double rSlide = 0.0;
    const ChVector3d phi = vt / std::max(speed, _stickSpeed);
    if (_friction > 0.0) {
        if (speed < _stickSpeed)
            rStick = _friction * fn / _stickSpeed;
        else
            rSlide = _friction * fn / speed;
    }
    for (std::size_t i = 0; i < _corners.size(); ++i) {
        const int ri = _wOffset[i];
        if (ri < 0)
            continue;
        for (std::size_t j = 0; j < _corners.size(); ++j) {
            const int cj = _wOffset[j];
            if (cj < 0)
                continue;
            const double s = _corners[i].share * _corners[j].share;
            addBlock(m_jacobians->K, ri, cj, n, n, s * kn);
            addBlock(m_jacobians->R, ri, cj, n, n, s * cn);
            if (_friction > 0.0) {
                addBlock(m_jacobians->K, ri, cj, phi, n, -s * _friction * kn);
                addBlock(m_jacobians->R, ri, cj, phi, n, -s * _friction * cn);
            }
            if (rStick > 0.0) {
                addIdentity(m_jacobians->R, ri, cj, s * rStick);
                addBlock(m_jacobians->R, ri, cj, n, n, -s * rStick);
            } else if (rSlide > 0.0) {
                addIdentity(m_jacobians->R, ri, cj, s * rSlide);
                addBlock(m_jacobians->R, ri, cj, tHat, tHat, -s * rSlide);
                addBlock(m_jacobians->R, ri, cj, n, n, -s * rSlide);
            }
        }
    }
}

void NodeFaceContact::ComputeQ(chrono::ChState* state_x, chrono::ChStateDelta* state_w) {
    evaluate(state_x, state_w, false);
}

void NodeFaceContact::ComputeJacobian(chrono::ChState* state_x, chrono::ChStateDelta* state_w) {
    if (!m_jacobians)
        CreateJacobianMatrices();
    m_jacobians->K.setZero();
    m_jacobians->R.setZero();
    m_jacobians->M.setZero();
    evaluate(state_x, state_w, true);
}

}  // namespace rcs
