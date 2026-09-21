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

    /// Update joint reactions and fail overloaded links.
    /// Returns true when the constraint topology changed.
    bool updateJoints(double dt);

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
    void applyPartContact(const Environment& environment, double dt);
    void applyUserLoads(const std::vector<ImpactLoad>& loads, double time);
    void refreshAdjacency();

    std::size_t globalNode(std::size_t part, std::size_t node) const {
        return _nodeOffset[part] + node;
    }

    std::vector<Part> _parts;
    std::vector<Joint> _joints;
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
    /// Deepest floor penetration of each node during the current contact.
    /// The rebound spring is scaled from this so concrete does not return the
    /// energy it just absorbed. Cleared when the node leaves the floor.
    std::vector<std::vector<double>> _floorPeak;
    struct LoadTarget {
        std::size_t part = 0;
        std::size_t node = 0;
        double share = 0.0;  ///< fraction of the load's total force
    };
    /// Per load, the nodes it acts on and their share of the total force.
    std::vector<std::vector<LoadTarget>> _loadTargets;
};

}  // namespace rcs
