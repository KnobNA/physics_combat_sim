#pragma once

#include "physics_core/joint.hpp"
#include "physics_core/part.hpp"
#include "physics_core/recording.hpp"
#include "physics_core/scenario.hpp"

#include "chrono/physics/ChSystem.h"

#include <cstddef>
#include <string>
#include <vector>

namespace rcs {

/// Authoring-time description of a whole machine: a bag of parts plus the
/// manual joints that hold them together. CAD import produces one of these.
struct AssemblyDefinition {
    std::string name = "Assembly";
    std::vector<PartDefinition> parts;
    std::vector<WeldDefinition> welds;
    std::vector<ScrewDefinition> screws;
};

/// Runtime assembly: owns the parts, the joint constraints, and the per-step
/// load assembly (gravity, aerodynamics, contact, user loads).
class Assembly {
public:
    void build(chrono::ChSystem& system, const AssemblyDefinition& definition);

    /// Place and launch the whole assembly as a rigid group.
    void place(const InitialState& state);

    /// Bind each user load to the nodes it acts on. Call after `place`.
    void bindLoads(const std::vector<ImpactLoad>& loads);

    /// What a load will actually do to the mesh it lands on. `deltaVPerStep`
    /// much above a metre per second means the step is too long to resolve
    /// the hit and the run will blow up.
    struct LoadAdvisory {
        std::size_t nodeCount = 0;
        double peakNodalForce = 0.0;
        double deltaVPerStep = 0.0;
    };
    LoadAdvisory adviseLoad(std::size_t loadIndex,
                            const std::vector<ImpactLoad>& loads,
                            double dt) const;

    /// Accumulate every external nodal force for the coming step.
    void applyLoads(const Environment& environment,
                    const std::vector<ImpactLoad>& loads,
                    double time,
                    double dt);

    /// Move rigid shards outside the shell solve, then sit them on the floor.
    void integrateShards(double dt, const Environment& environment);

    /// Update joint reactions and fail overloaded links.
    /// Returns true when the constraint topology changed.
    bool updateJoints(double dt);

    /// Detach shell elements past the material's tear limit. True if nodes were added.
    bool tearOverstretched();

    /// Promote a weld that was skipped as redundant, once a tree link has broken.
    /// True when a new constraint was added.
    bool respan();

    void writeLayout(RecordingLayout& layout) const;
    void writeFrame(Frame& frame) const;

    bool finite() const;
    chrono::ChVector3d center() const;
    chrono::ChVector3d centerOf(const std::vector<std::size_t>& partIndices) const;
    double mass() const;
    /// Parts reachable from `seed` through joints that are still intact.
    std::vector<std::size_t> connectedParts(std::size_t seed) const;

    const std::vector<Part>& parts() const { return _parts; }
    std::vector<Part>& parts() { return _parts; }
    const std::vector<Joint>& joints() const { return _joints; }
    int brokenLinkCount() const;

    /// `place` maps authoring coordinates to the world as
    /// world = placementRotation * (authored - authoredCenter) + placementPosition.
    /// Editor gizmos need the inverse to write edits back into the definition.
    const chrono::ChVector3d& authoredCenter() const { return _authoredCenter; }
    const chrono::ChQuaterniond& placementRotation() const { return _placementRotation; }
    const chrono::ChVector3d& placementPosition() const { return _placementPosition; }

private:
    void buildWeld(chrono::ChSystem& system,
                   const WeldDefinition& definition,
                   std::vector<int>& translationSet,
                   std::vector<int>& rotationSet);
    void buildScrew(chrono::ChSystem& system,
                    const ScrewDefinition& definition,
                    std::vector<int>& translationSet,
                    std::vector<int>& rotationSet);
    void reindexNodes();
    void createLink(chrono::ChSystem& system,
                    Joint& joint,
                    std::size_t partA,
                    std::size_t nodeA,
                    std::size_t partB,
                    std::size_t nodeB,
                    double breakForce,
                    double breakTorque,
                    bool wantRotationLock,
                    bool lockRotation);

    void addLink(chrono::ChSystem& system,
                 Joint& joint,
                 std::size_t partA,
                 std::size_t nodeA,
                 std::size_t partB,
                 std::size_t nodeB,
                 double breakForce,
                 double breakTorque,
                 bool wantRotationLock,
                 std::vector<int>& translationSet,
                 std::vector<int>& rotationSet);

    void applyGravity(const chrono::ChVector3d& gravity);
    void applyAerodynamics(const Environment& environment);
    void applyFloorContact(const Environment& environment, double dt);
    void applyShardContact(const Environment& environment, double dt);
    void separateShards();
    void applyPartContact(const Environment& environment, double dt);
    void applyUserLoads(const std::vector<ImpactLoad>& loads, double time);
    void refreshAdjacency();

    std::size_t globalNode(std::size_t part, std::size_t node) const {
        return _nodeOffset[part] + node;
    }

    struct DeferredLink {
        std::size_t partA = 0;
        std::size_t nodeA = 0;
        std::size_t partB = 0;
        std::size_t nodeB = 0;
        double breakForce = 0.0;
        double breakTorque = 0.0;
        bool wantRotation = false;
        std::string label;
        JointKind kind = JointKind::Weld;
    };

    std::vector<Part> _parts;
    std::vector<Joint> _joints;
    std::vector<DeferredLink> _deferred;
    chrono::ChSystem* _system = nullptr;
    chrono::ChVector3d _authoredCenter{0, 0, 0};
    chrono::ChQuaterniond _placementRotation{1, 0, 0, 0};
    chrono::ChVector3d _placementPosition{0, 0, 0};
    std::vector<std::size_t> _nodeOffset;
    /// Parts joined by an intact joint skip mutual contact: their meshes share
    /// a seam and would otherwise push each other apart forever.
    std::vector<char> _jointAdjacency;
    /// Per part pair: whether contact between them has been armed. See
    /// applyPartContact for why a pair starts disarmed.
    std::vector<char> _contactArmed;
    /// Deepest floor penetration of each node during a real inbound impact.
    /// Restitution reads it once on the way out. A resting sag does not.
    std::vector<std::vector<double>> _floorPeak;
    /// 0 resting, 1 a fast inbound impact is in progress.
    std::vector<std::vector<char>> _floorPhase;
    /// 1 when the node was clear of the floor on the previous step. An impact
    /// is latched only on entry, so a resting sag cannot turn the stiff spring
    /// back on.
    std::vector<std::vector<char>> _floorClear;
    /// Last part-pair contact normal, so a closest-point jump cannot kick sideways.
    std::vector<chrono::ChVector3d> _contactNormal;
    struct LoadTarget {
        std::size_t part = 0;
        std::size_t node = 0;
        double share = 0.0;  ///< fraction of the load's total force
    };
    /// Per load, the nodes it acts on and their share of the total force.
    std::vector<std::vector<LoadTarget>> _loadTargets;
};

}  // namespace rcs
