#pragma once

#include "physics_core/contact_load.hpp"
#include "physics_core/heat.hpp"
#include "physics_core/joint.hpp"
#include "physics_core/part.hpp"
#include "physics_core/recording.hpp"
#include "physics_core/scenario.hpp"

#include "chrono/fea/ChLinkNodeNode.h"
#include "chrono/physics/ChLoadContainer.h"
#include "chrono/physics/ChSystem.h"

#include <array>
#include <memory>
#include <utility>

#include <cstddef>
#include <cstdint>
#include <string>
#include <unordered_map>
#include <vector>

namespace rcs {

/// An empty folder in the assembly. It has no mesh. Translating it moves every
/// part that lists it as `group`, and the parts parented under those.
struct GroupDefinition {
    std::string name = "Group";
    /// Another folder. -1 is a root folder.
    int parent = -1;
    /// Copied onto every child when the folder checkbox is used.
    bool simulate = true;
    /// Sibling rank among folders that share `parent`.
    int order = 0;
};

/// Authoring-time description of a whole machine: a bag of parts plus the
/// manual joints that hold them together. CAD import produces one of these.
struct AssemblyDefinition {
    std::string name = "Assembly";
    std::vector<GroupDefinition> groups;
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

    /// Conduction, still air, and heat sources. Skipped entirely when the
    /// scene is still at its initial temperature and no source is on.
    void advanceHeat(double dt, double time, const Environment& environment,
                     const std::vector<HeatSource>& sources);
    /// Thermal strain or curvature is in the elements (a continuum load).
    bool hasThermalPrestress() const;
    /// Strain temperature or prescribed curvature changed this step.
    bool thermalLoadPending() const;
    void clearThermalLoadChanged();

    /// Accumulate every external nodal force for the coming step. True when
    /// the implicit contact set changed and the solver needs a fresh setup.
    bool applyLoads(const Environment& environment,
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
    /// Detach molten tetrahedra on the same slow clock as heat. The mesh
    /// itself is chosen at simulate start and does not grow during the run.
    bool maintainVolume(double dt);
    /// Split a fast floor approach so one substep does not jump through the
    /// contact patch (closing speed times dt vs sqrt of nodal area).
    int impactSubsteps(const Environment& environment, double dt) const;
    /// Contact onset uses this as the free-flight length of one substep.
    void setContactStep(double h);
    /// Throw away tet trial state from a rejected HHT step.
    void discardTrial();
    /// Keep tet plastic from the last accepted Newton iterate.
    void commitPlastic();
    void enforceRigidNodes();
    double peakPlastic() const;
    void captureKinematics();
    void restoreKinematics();
    /// Unused on the live HHT path. Kept so leftover call sites compile.
    void applyFloorBounce(const Environment& environment, bool newtonFailed = false);
    /// True when any live node is within `pad` of the floor contact plane.
    bool touchingFloor(const Environment& environment, double pad) const;
    /// Surface nodes currently penetrating any arena plane (floor, walls,
    /// ceiling) or another part's face.
    int activeContactCount() const;
    /// Fastest closing speed of a live contact node into an arena plane, m/s.
    double maxApproachSpeed() const;
    /// Soonest time for a live skin node to reach an arena plane under
    /// `gravity`. 0 if already on/in a plane; `kNoContactTime` if none close.
    double minTimeToContact(const chrono::ChVector3d& gravity) const;
    /// Shortest half period π√(m/k) over penetrating contact nodes, seconds.
    /// `kNoContactTime` when nothing touches.
    double contactHalfPeriod() const;
    /// Latch each clear node's velocity as its contact approach velocity.
    void commitContactApproach();
    /// Deepest live-node penetration into any arena plane, metres.
    double maxContactPenetration() const;
    /// Deepest floor overshoot of a free node, metres. 0 if none or no floor.
    double maxFloorOvershoot(const Environment& environment) const;
    /// True when any live brittle part is at or past Rankine on the trial pose.
    bool pastBrittleRankine() const;
    /// Mesh-scaled length a contact step may penetrate before it is rejected.
    double contactLengthScale() const;
    /// Unused on the live HHT path. Kept so leftover call sites compile.
    void ballisticAdvance(double dt, const chrono::ChVector3d& gravity);
    /// Unused on the live HHT path. Kept so leftover call sites compile.
    void rigidAdvance(double dt, const Environment& environment);
    /// Integrate the already-assembled nodal loads (gravity, user punches,
    /// part contact) when HHT did not converge, then the caller seats the floor.
    void explicitAdvance(double dt);
    void refreshTetTrials();
    /// Translational kinetic energy of live nodes and rigid chips, joules.
    double kineticEnergy() const;
    /// Penalty spring energy ½ k δ² on live arena contacts, joules.
    double contactElasticEnergy() const;
    /// Mass-weighted centre-of-mass height, metres.
    double massHeight() const;
    /// Project live (non-molten, non-chip) nodes out of other live volumes.
    void correctInterpenetration();
    /// Remove chips that have not settled. Prefer `sanitizeShards` on a
    /// rejected step; deleting flying chips is what made a pane vanish.
    bool discardAwakeShards();
    /// Reset non-finite rigid chips in place so they cannot abort the FEA step.
    void sanitizeShards();

    /// Promote a weld that was skipped as redundant, once a tree link has broken.
    /// True when a new constraint was added.
    bool respan();

    void writeLayout(RecordingLayout& layout) const;
    void writeFrame(Frame& frame) const;

    bool finite() const;
    /// Live tetrahedra whose current volume has flipped or flattened.
    int invertedTetCount() const;
    /// True when at least one part still has a free FEA node for HHT to solve.
    bool hasFreeSolve() const;
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
                   std::vector<int>& rotationSet,
                   std::vector<char>& snapped);
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

    bool addLink(chrono::ChSystem& system,
                 Joint& joint,
                 std::size_t partA,
                 std::size_t nodeA,
                 std::size_t partB,
                 std::size_t nodeB,
                 double breakForce,
                 double breakTorque,
                 bool wantRotationLock,
                 std::vector<int>& translationSet,
                 std::vector<int>& rotationSet,
                 double breakPeel,
                 double criticalOpening,
                 const chrono::ChVector3d& peelAxis);

    void applyGravity(const chrono::ChVector3d& gravity);
    void applyAerodynamics(const Environment& environment);
    /// Create or refresh one implicit contact element per surface node, and
    /// copy the current arena into the shared world. Called before the solve.
    /// True when elements were added or removed, so the caller has to re-run
    /// Setup and rebuild the sparsity pattern.
    bool syncContacts(const Environment& environment);
    void applyShardContact(const Environment& environment, double dt);
    void separateShards();
    void applyPartContact(const Environment& environment, double dt);
    void applyUserLoads(const std::vector<ImpactLoad>& loads, double time);
    void refreshAdjacency();
    void inheritWeldNodes();
    void classifySeams(std::vector<char>& state) const;

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
        double breakPeel = 0.0;
        double criticalOpening = 0.0;
        chrono::ChVector3d peelAxis{0.0, 1.0, 0.0};
        bool wantRotation = false;
        std::string label;
        JointKind kind = JointKind::Weld;
    };

    HeatField _heat;
    std::vector<WeldDefinition> _weldDefs;
    double _volumeWait = 0.0;
    std::vector<Part> _parts;
    std::vector<Joint> _joints;
    std::vector<std::shared_ptr<chrono::fea::ChLinkNodeNode>> _continuity;
    std::vector<DeferredLink> _deferred;
    chrono::ChSystem* _system = nullptr;
    chrono::ChVector3d _authoredCenter{0, 0, 0};
    chrono::ChQuaterniond _placementRotation{1, 0, 0, 0};
    chrono::ChVector3d _placementPosition{0, 0, 0};
    std::vector<std::size_t> _nodeOffset;
    /// Still used by connectedParts; part-part contact no longer blanks a
    /// whole pair just because a weld pin remains.
    std::vector<char> _jointAdjacency;
    /// Per part pair: 0 = not yet colliding, 1 = first contact step (overlap
    /// capped), 2 = ongoing. Overlapping pairs arm immediately so a broken
    /// weld cannot fold through, and so stacked volumes actually collide.
    std::vector<char> _contactArmed;
    /// Arena and friction as the contact elements see them.
    std::shared_ptr<ContactWorld> _contactWorld;
    std::shared_ptr<chrono::ChLoadContainer> _contactContainer;
    /// One element per node of each part, indexed the same way. Null where the
    /// node has no outside face.
    std::vector<std::vector<std::shared_ptr<NodeContact>>> _contacts;
    /// Node count, skin count, and free-surface-node count the contact
    /// elements were built against. A refine, a detach, or a node being
    /// frozen rebuilds them; nothing else does.
    std::vector<std::array<std::size_t, 3>> _contactShape;
    /// Last part-pair contact normal, so a closest-point jump cannot kick sideways.
    std::vector<chrono::ChVector3d> _contactNormal;
    /// Hunt-Crossley approach speed per shard penalty contact, latched when
    /// the contact opens. A contact not seen on the previous evaluation is
    /// new and re-latches.
    struct ApproachLatch {
        double speed = 0.0;
        unsigned stamp = 0;
    };
    std::unordered_map<std::uint64_t, ApproachLatch> _penaltyApproach;
    unsigned _penaltyStamp = 0;
    double penaltyDamping(std::uint64_t key, double vn);
    /// Implicit node-to-face contacts between parts, keyed by part pair and
    /// node. A pair not found on a user step is removed; a pair whose
    /// triangle changed is rebuilt. Either changes the solver structure.
    struct PairContact {
        std::shared_ptr<NodeFaceContact> load;
        std::array<int, 3> tri{-1, -1, -1};
        unsigned stamp = 0;
    };
    std::unordered_map<std::uint64_t, PairContact> _pairContacts;
    unsigned _pairStamp = 0;
    bool _pairsChanged = false;
    /// Pool the contact of `node` on `nodePart` against triangle `tri` of
    /// `facePart` for this user step.
    void requestPairContact(std::size_t nodeSide,
                            std::size_t node,
                            std::size_t faceSide,
                            const std::array<int, 3>& tri,
                            const chrono::ChVector3d& normal,
                            double wa,
                            double wb,
                            double wc,
                            double reach,
                            double stiffness,
                            double friction,
                            double stickSpeed);
    /// Drop the pair contacts not requested this user step, or all of them.
    void sweepPairContacts(bool all = false);
    struct LoadTarget {
        std::size_t part = 0;
        std::size_t node = 0;
        double share = 0.0;  ///< fraction of the load's total force
    };
    /// Per load, the nodes it acts on and their share of the total force.
    std::vector<std::vector<LoadTarget>> _loadTargets;
};

}  // namespace rcs
