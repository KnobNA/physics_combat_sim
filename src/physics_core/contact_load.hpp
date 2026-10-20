#pragma once

#include "chrono/fea/ChNodeFEAxyz.h"
#include "chrono/fea/ChNodeFEAxyzrot.h"
#include "chrono/physics/ChLoad.h"

#include <memory>
#include <vector>

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
    /// Restitution of the arena surfaces: the impulse for rigid chips and
    /// the Hunt-Crossley damping of the deformable penalty contact.
    double restitution = 0.15;
};

/// One node pressed against the floor, the walls, or the ceiling, written as a
/// Chrono load so the implicit integrator sees its stiffness.
///
/// Penalty is continuum: F_n = k δ (softened at touch) with Hunt-Crossley
/// hysteresis damping set by the surface restitution. k is a fraction of
/// element stiffness E √A, not m/dt².
class NodeContact : public chrono::ChLoadCustom {
public:
    NodeContact(std::shared_ptr<chrono::fea::ChNodeFEAxyz> node,
                std::shared_ptr<const ContactWorld> world);
    NodeContact(std::shared_ptr<chrono::fea::ChNodeFEAxyzrot> node,
                std::shared_ptr<const ContactWorld> world);

    /// Contact geometry and gains for this node. Rewritten whenever the skin
    /// changes, so a refined or detached face contacts with its new area.
    /// `stiffness` is the material penalty k = c E √A. `dissipation` and
    /// `mass` are unused (kept so call sites stay stable); damping comes from
    /// the surface restitution.
    void configure(double offset, double stiffness, double dissipation, double mass = 0.0);

    /// A node that has left the deformable body (a rigid chip, a frozen
    /// scrap) is driven directly and must not also be pushed by a contact.
    void setActive(bool active) { _active = active; }
    bool active() const { return _active; }

    /// Normal force carried on the last evaluation. Used for reporting and
    /// for the thermal contact conductance.
    double normalForce() const { return _normalForce; }
    double stiffness() const { return _stiffness; }
    /// Record the node velocity as its approach velocity while it is clear
    /// of every plane. Called on accepted increments only, so the
    /// Hunt-Crossley factor is fixed within a Newton solve.
    void commitApproach();
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
    chrono::ChVector3d _approach{0, 0, 0};
    bool _active = true;
};

/// One node against one triangle of another part (or a separated piece of the
/// same part), written as a Chrono load so the implicit integrator sees its
/// stiffness. The pair, its normal and the barycentric point are found once
/// per user step, as an implicit FE contact search does at the start of an
/// increment; the gap and the relative velocity are evaluated at every
/// Newton iterate. Normal force is the softened penalty with the same
/// Hunt-Crossley damping as the arena contact.
///
/// A fixed node has no solver variables and is not a loadable here: its
/// position is read directly and its reaction is not applied.
class NodeFaceContact : public chrono::ChLoadCustomMultiple {
public:
    struct Corner {
        std::shared_ptr<chrono::fea::ChNodeFEAxyz> xyz;
        std::shared_ptr<chrono::fea::ChNodeFEAxyzrot> rot;
        /// +1 for the contacting node, minus the barycentric weight for a
        /// triangle corner. The gap is reach − n · Σ s x.
        double share = 0.0;
    };

    /// `corners[0]` is the node; the rest are the triangle corners.
    /// `loadables` lists the free corners only, in `corners` order.
    NodeFaceContact(std::vector<std::shared_ptr<chrono::ChLoadable>>& loadables, std::vector<Corner> corners);

    /// Rewritten each user step while the pair persists.
    void configure(const chrono::ChVector3d& normal,
                   const std::vector<double>& shares,
                   double reach,
                   double stiffness,
                   double friction,
                   double stickSpeed,
                   double restitution);

    /// Latch the closing speed as the approach speed while the pair is open,
    /// on accepted increments only.
    void commitApproach();
    /// Gap and spring energy at the current node positions.
    double gap() const;
    double elasticEnergy() const;

    NodeFaceContact* Clone() const override { return new NodeFaceContact(*this); }
    void ComputeQ(chrono::ChState* state_x, chrono::ChStateDelta* state_w) override;
    void ComputeJacobian(chrono::ChState* state_x, chrono::ChStateDelta* state_w) override;
    bool IsStiff() override { return true; }

private:
    void evaluate(const chrono::ChState* state_x, const chrono::ChStateDelta* state_w, bool jac);
    chrono::ChVector3d cornerPos(std::size_t c, const chrono::ChState* state_x) const;
    chrono::ChVector3d cornerVel(std::size_t c, const chrono::ChStateDelta* state_w) const;

    std::vector<Corner> _corners;
    /// Offsets of each corner into the load state and Q, or −1 when fixed.
    std::vector<int> _xOffset;
    std::vector<int> _wOffset;
    chrono::ChVector3d _normal{0, 1, 0};
    double _reach = 0.0;
    double _stiffness = 0.0;
    double _friction = 0.0;
    double _stickSpeed = 1.0e-3;
    double _restitution = 0.15;
    double _approach = 0.0;
    bool _latched = false;
};

}  // namespace rcs
