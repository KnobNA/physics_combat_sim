#include "physics_core/contact_load.hpp"

#include <algorithm>
#include <cmath>
#include <cstdlib>

namespace rcs {

using chrono::ChVector3d;

NodeContact::NodeContact(std::shared_ptr<chrono::fea::ChNodeFEAxyz> node,
                         std::shared_ptr<const ContactWorld> world)
    : chrono::ChLoadCustom(node), _node(std::move(node)), _world(std::move(world)) {}

void NodeContact::configure(double offset, double stiffness, double dissipation) {
    _offset = offset;
    _stiffness = stiffness;
    _damping = dissipation;
}

void NodeContact::addPlane(const ChVector3d& position,
                           const ChVector3d& velocity,
                           const ChVector3d& normal,
                           double planeOffset,
                           ChVector3d& force)
{
    // Positive when the node has crossed the surface of the plane.
    const double gap = planeOffset + _offset - position.Dot(normal);
    if (gap <= 0.0)
        return;
    const double vn = velocity.Dot(normal);
    // A spring that switches on at full stiffness is not differentiable at
    // touchdown, and Newton then chatters across the surface instead of
    // converging. Ramp the stiffness in over `kRamp` so force and slope are
    // both continuous; beyond that it is the plain linear spring.
    constexpr double kRamp = 2.0e-5;
    const double stiff = gap < kRamp ? _stiffness * gap * gap / (2.0 * kRamp)
                                     : _stiffness * (gap - 0.5 * kRamp);
    // Hunt-Crossley dissipation: the stiffness itself is scaled by the closing
    // rate, so the surface pushes harder going in than coming out and the loop
    // encloses an area. A dashpot in parallel would instead have to be sized
    // against the mass being stopped, and the tributary mass of one surface
    // node is nowhere near the mass of the body arriving behind it, which is
    // what turns a floor into a trampoline. This form needs no mass at all:
    // the restitution falls off as roughly 1 - _damping * closing speed.
    // For small rates this is the plain `1 + a*closing` of Hunt-Crossley. The
    // exponential-of-tanh form keeps that slope near zero while staying smooth,
    // strictly positive and bounded: clamping the linear form at zero puts a
    // kink in the force that Newton has to step across, and letting it grow
    // without limit makes a fast impact arbitrarily stiff.
    // The cap matters as much as the coefficient. Unbounded, a fast impact
    // makes the contact several times stiffer than the material behind it, and
    // the Newton solve on a loose panel slamming the floor stops converging.
    // At 0.9 the surface never gets more than about two and a half times its
    // resting stiffness, which is enough hysteresis to hold restitution near a
    // half without wrecking the step.
    static const double kSaturation = [] {
        const char* override = std::getenv("RCS_CONTACT_SATURATION");
        return override ? std::atof(override) : 0.9;
    }();
    const double closing = -vn;
    const double hysteresis = std::exp(kSaturation * std::tanh(_damping * closing / kSaturation));
    const double fn = stiff * hysteresis;
    if (fn <= 0.0)
        return;
    _normalForce = std::max(_normalForce, fn);
    force += normal * fn;

    const ChVector3d vt = velocity - normal * vn;
    const double speed = vt.Length();
    if (speed <= 1.0e-12)
        return;
    // Regularised Coulomb: full mu once sliding, ramped below the stick speed
    // so the force stays differentiable through zero. Friction follows the
    // elastic pressure, not the total normal force: the hysteresis term is
    // energy going into the substrate, not extra load pressing the faces
    // together. Letting it drive friction makes the tangential force grow
    // with the normal approach rate, which couples the two directions and
    // feeds the sliding motion instead of opposing it.
    const double scale = std::min(1.0, speed / std::max(_world->stickSpeed, 1.0e-9));
    force -= vt * (_world->friction * stiff * scale / speed);
}

void NodeContact::ComputeQ(chrono::ChState* state_x, chrono::ChStateDelta* state_w) {
    if (load_Q.size() != 3)
        load_Q.setZero(3);
    load_Q.setZero();
    _normalForce = 0.0;
    if (!_active || !_world || _stiffness <= 0.0)
        return;

    const ChVector3d position = state_x ? ChVector3d((*state_x)(0), (*state_x)(1), (*state_x)(2))
                                        : _node->GetPos();
    const ChVector3d velocity = state_w ? ChVector3d((*state_w)(0), (*state_w)(1), (*state_w)(2))
                                        : _node->GetPosDt();

    ChVector3d force(0, 0, 0);
    if (_world->floorEnabled)
        addPlane(position, velocity, ChVector3d(0, 1, 0), _world->floorHeight, force);
    addPlane(position, velocity, ChVector3d(1, 0, 0), -_world->halfX, force);
    addPlane(position, velocity, ChVector3d(-1, 0, 0), -_world->halfX, force);
    addPlane(position, velocity, ChVector3d(0, 0, 1), -_world->halfZ, force);
    addPlane(position, velocity, ChVector3d(0, 0, -1), -_world->halfZ, force);
    addPlane(position, velocity, ChVector3d(0, -1, 0), -_world->ceiling, force);

    load_Q(0) = force.x();
    load_Q(1) = force.y();
    load_Q(2) = force.z();
}

}  // namespace rcs
