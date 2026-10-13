#pragma once

#include "chrono/fea/ChNodeFEAxyz.h"
#include "chrono/fea/ChNodeFEAxyzrot.h"
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
    /// Integrator substep, seconds.
    double stepSize = 5.0e-4;
    /// Rigid-chip bounce only. Deformable residual uses continuum penalty
    /// (k δ); it does not read this.
    double restitution = 0.15;
};

/// One node pressed against the floor, the walls, or the ceiling, written as a
/// Chrono load so the implicit integrator sees its stiffness.
///
/// Penalty is continuum: F_n = k δ plus critical-fraction damping on that k,
/// and extra viscous stabilization only while the gap is closing. k is a
/// fraction of element stiffness E √A, not m/dt². Rebound comes from
/// continuum elasticity and plasticity, not a rigid-body χ.
class NodeContact : public chrono::ChLoadCustom {
public:
    NodeContact(std::shared_ptr<chrono::fea::ChNodeFEAxyz> node,
                std::shared_ptr<const ContactWorld> world);
    NodeContact(std::shared_ptr<chrono::fea::ChNodeFEAxyzrot> node,
                std::shared_ptr<const ContactWorld> world);

    /// Contact geometry and gains for this node. Rewritten whenever the skin
    /// changes, so a refined or detached face contacts with its new area.
    /// `stiffness` is the material penalty k = c E √A. `dissipation` is unused
    /// (kept so call sites stay stable). `mass` sets 0.2-of-critical damping
    /// on that k.
    void configure(double offset, double stiffness, double dissipation, double mass = 0.0);

    /// A node that has left the deformable body (a rigid chip, a frozen
    /// scrap) is driven directly and must not also be pushed by a contact.
    void setActive(bool active) { _active = active; }
    bool active() const { return _active; }

    /// Normal force carried on the last evaluation. Used for reporting and
    /// for the thermal contact conductance.
    double normalForce() const { return _normalForce; }
    double stiffness() const { return _stiffness; }
    double offset() const { return _offset; }

    NodeContact* Clone() const override { return new NodeContact(*this); }

    void ComputeQ(chrono::ChState* state_x, chrono::ChStateDelta* state_w) override;
    void ComputeJacobian(chrono::ChState* state_x, chrono::ChStateDelta* state_w) override;
    /// The load depends on position, so the integrator has to differentiate
    /// it. Analytical K = −dQ/dx matches Chrono's convention for ChLoadCustom.
    bool IsStiff() override { return true; }

private:
    /// Accumulate one half-space into `force`, given the node state.
    /// When `jac` is set, also add analytical K = −dQ/dx and R = −dQ/dv
    /// into the load Jacobians.
    void addPlane(const chrono::ChVector3d& position,
                  const chrono::ChVector3d& velocity,
                  const chrono::ChVector3d& normal,
                  double planeOffset,
                  int slot,
                  chrono::ChVector3d& force,
                  bool jac);

    bool _fillJac = false;

    std::shared_ptr<chrono::fea::ChNodeFEAxyz> _node;
    std::shared_ptr<chrono::fea::ChNodeFEAxyzrot> _rotNode;
    std::shared_ptr<const ContactWorld> _world;
    double _offset = 0.0;
    double _stiffness = 0.0;
    double _damping = 0.0;
    double _mass = 0.0;
    double _normalForce = 0.0;
    double _planeForce[6] = {};
    chrono::ChVector3d _planeNormal[6] = {};
    bool _active = true;
};

}  // namespace rcs
