#pragma once

#include "chrono/fea/ChNodeFEAxyz.h"
#include "chrono/physics/ChLoad.h"

#include <memory>

namespace rcs {

/// Ground and arena state shared by every contact element. One instance per
/// assembly, rewritten each step, so the scenario can move the floor or change
/// the friction without rebuilding the loads.
struct ContactWorld {
    bool floorEnabled = true;
    double floorHeight = 0.0;
    double friction = 0.6;
    /// Sliding speed below which friction is ramped in linearly. A true
    /// Coulomb step is not differentiable and would stall the Newton loop.
    double stickSpeed = 1.0e-3;
    /// Half extents of the arena. The walls and ceiling use the same law.
    double halfX = 5.0;
    double halfZ = 5.0;
    double ceiling = 10.0;
};

/// One node pressed against the floor, the walls, or the ceiling, written as a
/// Chrono load so the implicit integrator sees its stiffness.
///
/// The explicit alternative cannot work here. A penalty force applied from the
/// state at the start of the step is stable only while `k < m/dt^2`, which for
/// a millimetre timestep is millions of times softer than metal. A body then
/// either sinks through the floor or is thrown off it. Inside the Jacobian the
/// same spring can carry the real `E*A/L` and the step stays stable.
class NodeContact : public chrono::ChLoadCustom {
public:
    NodeContact(std::shared_ptr<chrono::fea::ChNodeFEAxyz> node,
                std::shared_ptr<const ContactWorld> world);

    /// Contact geometry and gains for this node. Rewritten whenever the skin
    /// changes, so a refined or detached face contacts with its new area.
    /// `dissipation` is the Hunt-Crossley coefficient in s/m: restitution
    /// drops as roughly 1 - dissipation * closing speed.
    void configure(double offset, double stiffness, double dissipation);

    /// A node that has left the deformable body (a rigid chip, a frozen
    /// scrap) is driven directly and must not also be pushed by a contact.
    void setActive(bool active) { _active = active; }
    bool active() const { return _active; }

    /// Normal force carried on the last evaluation. Used for reporting and
    /// for the thermal contact conductance.
    double normalForce() const { return _normalForce; }

    NodeContact* Clone() const override { return new NodeContact(*this); }

    void ComputeQ(chrono::ChState* state_x, chrono::ChStateDelta* state_w) override;
    /// The load depends on position, so the integrator has to differentiate
    /// it. This is what makes the contact implicit.
    bool IsStiff() override { return true; }

private:
    /// Accumulate one half-space into `force`, given the node state.
    void addPlane(const chrono::ChVector3d& position,
                  const chrono::ChVector3d& velocity,
                  const chrono::ChVector3d& normal,
                  double planeOffset,
                  chrono::ChVector3d& force);

    std::shared_ptr<chrono::fea::ChNodeFEAxyz> _node;
    std::shared_ptr<const ContactWorld> _world;
    double _offset = 0.0;
    double _stiffness = 0.0;
    /// Hunt-Crossley dissipation coefficient, s/m.
    double _damping = 0.0;
    double _normalForce = 0.0;
    bool _active = true;
};

}  // namespace rcs
