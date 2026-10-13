#pragma once

#include "physics_core/material.hpp"
#include "physics_core/reissner_j2.hpp"
#include "physics_core/thermal_kirchhoff.hpp"
#include "physics_core/thermal_tet.hpp"

#include "chrono/core/ChQuaternion.h"
#include "chrono/core/ChVector3.h"
#include "chrono/fea/ChElementShellBST.h"
#include "chrono/fea/ChElementShellReissner4.h"
#include "chrono/fea/ChMesh.h"
#include "chrono/fea/ChNodeFEAxyz.h"
#include "chrono/fea/ChNodeFEAxyzrot.h"

#include <array>
#include <cstddef>
#include <map>
#include <memory>
#include <string>
#include <utility>
#include <vector>

namespace rcs {

/// A flat rectangular shell patch on the part's mid-surface. CAD import
/// produces one PlatePatch per planar face it recognises; the built-in test
/// cube uses six of them.
struct PlatePatch {
    chrono::ChVector3d origin{0.0, 0.0, 0.0};  ///< corner of the patch
    chrono::ChVector3d uAxis{1.0, 0.0, 0.0};
    chrono::ChVector3d vAxis{0.0, 0.0, 1.0};
    double uSize = 1.0;
    double vSize = 1.0;
    int divU = 4;
    int divV = 4;
};

/// A planar authored polygon. Triangle fans in `triangles` are derived from
/// these; the editor picks and splices faces, not tet faces.
struct AuthoredFace {
    std::vector<int> loop;
};

/// A triangulated surface. Non-empty means this part is a triangle shell
/// (ChElementShellBST) rather than a rectangular Reissner patch.
struct SurfaceMesh {
    std::vector<chrono::ChVector3d> vertices;
    std::vector<std::array<int, 3>> triangles;
    std::vector<AuthoredFace> faces;
    /// Per-triangle thickness, metres. Empty means the part thickness.
    std::vector<double> thickness;
    /// PLC facet marker. 2 is a weld or ply interface. Empty means every face is 1.
    std::vector<int> markers;

    bool empty() const { return vertices.empty() || triangles.empty(); }
};

/// Element-size knobs chosen at import. Changing them and clearing `cached`
/// rebuilds this part only.
struct MeshSizing {
    double shellCutoff = 0.008;  ///< m, thinner than this is extruded as one tet through the thickness
    double minElement = 0.002;   ///< m, finest edge runtime refinement may reach
    double maxElement = 0.025;   ///< m, finest edge the initial surface may use
    int throughThickness = 1;
    /// Resting tetrahedra. `meshResolution` scales this at simulate time.
    int tetBudget = 2500;
    /// Finest edge the heat clock may split a tetrahedron down to.
    double refineLimit = 0.025;
    bool cached = false;
    int generation = 0;
};

/// Tetrahedral fill of the thick region. Empty for a pure shell.
struct VolumeMesh {
    std::vector<chrono::ChVector3d> vertices;
    std::vector<std::array<int, 4>> tets;
    std::vector<std::array<int, 3>> boundary;
    /// Parallel to `boundary`. 2 is an interface face (weld or ply).
    std::vector<int> boundaryMarker;
    double elementSize = 0.025;
    bool empty() const { return tets.empty(); }
};

/// Authoring-time description of one rigidly-manufactured piece of the
/// assembly. Nothing downstream knows this came from a cube.
struct PartDefinition {
    std::string name = "Part";
    Material material;
    double thickness = 0.008;  ///< m
    PlatePatch patch;
    /// Authored triangles. The sized shell and solid meshes are derived from
    /// this once and then reused.
    SurfaceMesh source;
    /// Thin-region midsurface. Empty when the whole part is a solid or a plate.
    SurfaceMesh surface;
    VolumeMesh volume;
    MeshSizing sizing;
    bool fixed = false;
    /// Hold the patch boundary still so a face can be loaded in bending.
    bool clampBoundary = false;
    /// The authored mesh is unstressed at this temperature.
    double initialTemperature = 293.15;  ///< K
    /// Empty-folder parent in `AssemblyDefinition::groups`. -1 means this part
    /// is not inside a folder.
    /// When false the part is not built, not drawn, and not joined.
    bool simulate = true;
    int group = -1;
    /// Another part that owns this one. -1 means it is not parented to a part.
    /// A folder and a part parent are not used together; the folder wins.
    int parentPart = -1;
    /// Sibling rank inside the same parent or folder. Lower draws first.
    int order = 0;
    /// Intrinsic XYZ Euler, radians. Solids keep this as their editor frame;
    /// plates also store orientation in `patch.uAxis` / `patch.vAxis`.
    chrono::ChVector3d euler{0.0, 0.0, 0.0};
};

/// Runtime representation: a shell FEA mesh plus the bookkeeping the solver
/// loop needs (lumped nodal mass, quad topology, per-node load accumulator).
class Part {
public:
    void build(const PartDefinition& definition);

    /// Rigidly move the part. `alsoRest` moves the undeformed reference too,
    /// which is what you want while authoring: the part must not be born
    /// pre-stressed just because it was dragged into place.
    void transform(const chrono::ChVector3d& translation,
                   const chrono::ChQuaterniond& rotation,
                   const chrono::ChVector3d& pivot,
                   bool alsoRest);

    void setVelocity(const chrono::ChVector3d& linear,
                     const chrono::ChVector3d& angular,
                     const chrono::ChVector3d& pivot);

    void clearLoads();
    void addLoad(std::size_t nodeIndex, const chrono::ChVector3d& force);
    chrono::ChVector3d nodeForce(std::size_t nodeIndex) const;
    /// A load on a node whose shell quads have all torn follows the rigid
    /// chips that replaced it. A node that is still in a live element keeps
    /// the force itself.
    void addStructuralLoad(std::size_t nodeIndex, const chrono::ChVector3d& force);
    void commitLoads();
    /// Recompute tet trial stress from the current node pose. Used after an
    /// explicit floor seat so Rankine/J2 see the compressed configuration.
    void refreshTetTrials();
    /// Chip clones must stay out of HHT; Chrono Setup() can clear Fixed.
    void enforceRigidNodes();

    /// Per-quad peak von Mises stress across the shell thickness, in Pa.
    /// Plastic strain is subtracted, so a dent that has unloaded reads low.
    void evaluateStress(std::vector<float>& out) const;
    void evaluatePlasticStrain(std::vector<float>& out) const;
    double tetPlasticStrain(std::size_t index) const;
    /// Highest outer-fiber tensile principal stress, Pa. Brittle fracture uses this.
    double peakTension() const;
    /// Throw away tet state from a rejected Newton step.
    void discardTrial();
    /// Keep tet plastic from the last accepted Newton iterate.
    void commitPlastic();
    /// Forward the current HHT slice length into every live tet.
    void setImpactStep(double h);
    /// Copy free-node positions and velocities so a rejected HHT step can
    /// be rolled back if Chrono left a non-finite iterate.
    void captureKinematics();
    void restoreKinematics();
    /// Copy the live pose into the undeformed reference so a free thermal
    /// swell does not fight the original rest in the next implicit step.
    void bakeRest();

    /// Detach elements that have stretched, yielded, or cracked past the
    /// material limit. Returns true when the mesh gained nodes.
    bool tearFailed();
    /// On the heat clock: split tetrahedra that are heating, yielding, or
    /// loaded, then let molten ones leave as closed chips.
    bool maintainVolume();

    /// Mean node temperature, kelvin. A new part is uniform at its initial
    /// temperature, which is the stress-free temperature of the authored mesh.
    double temperature() const { return _temperature; }
    double initialTemperature() const { return _initialTemperature; }
    void setTemperature(double kelvin);
    double nodeTemperature(std::size_t index) const;
    void setNodeTemperature(std::size_t index, double kelvin);
    /// Temperature already written into the shell strain. It follows the node
    /// temperature, but slowly enough that one heat update does not shock the
    /// implicit step.
    double strainTemperature(std::size_t index) const;
    void setStrainTemperature(std::size_t index, double kelvin);
    /// True when this node is inside a source or out of equilibrium with a
    /// neighbour, so weld melting only has to look at these.
    bool heatActive(std::size_t index) const;
    void setHeatActive(std::size_t index, bool active);
    /// Write thermal strain, and stiffness or yield only when the knockdown
    /// has moved. Also grows rigid chips with alpha * deltaT.
    void syncThermalMechanics();
    /// Through-thickness thermal curvature for a Reissner midsurface, 1/m.
    /// Used when the mesh has no nodes through the thickness to fit dT/dn.
    void setThermalCurvature(double kappa);
    /// True when live elements carry thermal strain or curvature.
    bool thermalPrestress() const { return _thermalPrestress; }
    /// True when the last sync or curvature write changed the thermal load.
    bool thermalLoadChanged() const { return _thermalLoadChanged; }
    void clearThermalLoadChanged() { _thermalLoadChanged = false; }
    void evaluateTemperature(std::vector<float>& out) const;

    std::size_t nodeCount() const;
    chrono::ChVector3d nodePosition(std::size_t index) const;
    chrono::ChVector3d restPosition(std::size_t index) const { return nodeRest(index); }
    chrono::ChVector3d nodeVelocity(std::size_t index) const;
    /// Pose captured by `captureKinematics` (the start of the last accepted
    /// attempt). Falls back to the live pose when nothing has been saved.
    chrono::ChVector3d savedPosition(std::size_t index) const;
    chrono::ChVector3d savedVelocity(std::size_t index) const;
    /// Move a free node without touching rest. Used to eject a tunneled node.
    void setNodeKinematics(std::size_t index,
                           const chrono::ChVector3d& pos,
                           const chrono::ChVector3d& vel);
    /// Reissner nodes carry rotation. Triangle shells are translation only.
    bool rotational() const { return _xyzNodes.empty(); }

    struct TearTwin {
        int clone = 0;
        int original = 0;
    };
    bool torn() const { return _torn; }
    const std::vector<TearTwin>& tearTwins() const { return _tearTwins; }

    chrono::ChVector3d center() const;
    void bounds(chrono::ChVector3d& lo, chrono::ChVector3d& hi, double pad = 0.0) const;
    double mass() const;
    bool finite() const;

    const std::string& name() const { return _name; }
    const Material& material() const { return _material; }
    double thickness() const { return _thickness; }
    bool fixed() const { return _fixed; }

    const std::shared_ptr<chrono::fea::ChMesh>& mesh() const { return _mesh; }
    const std::vector<std::shared_ptr<chrono::fea::ChNodeFEAxyzrot>>& nodes() const { return _nodes; }
    const std::vector<std::shared_ptr<chrono::fea::ChNodeFEAxyz>>& xyzNodes() const { return _xyzNodes; }
    const std::vector<double>& nodeMass() const { return _nodeMass; }
    const std::vector<std::array<int, 4>>& quads() const { return _quads; }
    const std::vector<std::array<int, 3>>& triangles() const { return _triangles; }
    const std::vector<std::array<int, 4>>& tets() const { return _tets; }
    bool hasVolume() const { return !_tetElements.empty(); }
    /// False when every incident tet has left the solve, so floor contact
    /// must not keep pushing an under-massed orphan.
    bool nodeOnLiveTet(std::size_t index) const;
    bool nodeOnLiveShell(std::size_t index) const;
    int skinTriangleCount() const { return _skinTriangles; }
    int tetCount() const { return static_cast<int>(_tets.size()); }
    int liveTetCount() const;
    /// Live tets whose current signed volume has flipped or flattened.
    int invertedTetCount() const;
    /// Live tetrahedra, or live shell quads/triangles when there is no volume.
    int liveElementCount() const;
    /// True when a live brittle element is at or past handbook Rankine. Used
    /// so a terminated contact iterate is not rewound before shatter can run.
    bool pastBrittleRankine() const;
    /// Authored tet count, or shell quad/triangle slots (including released).
    int totalElementCount() const;
    bool tetLive(std::size_t index) const;
    /// False when every node is frozen (the whole solid has become chips).
    bool hasFreeNodes() const;
    /// Nodes on the outline of the patch, the only ones a weld can pair up.
    const std::vector<int>& boundaryNodes() const { return _boundaryNodes; }

    /// A chip that has left the shell solve. Four corners, one rigid motion.
    struct RigidShard {
        std::array<int, 4> nodes{};
        chrono::ChVector3d com;
        chrono::ChQuaterniond rotation{1, 0, 0, 0};
        chrono::ChVector3d velocity;
        chrono::ChVector3d angularVelocity;
        chrono::ChVector3d pendingForce;
        chrono::ChVector3d pendingTorque;
        std::array<chrono::ChVector3d, 4> local{};
        /// Corner offsets at the part's initial temperature. `local` is this
        /// grown by alpha * deltaT, so a later temperature does not compound.
        std::array<chrono::ChVector3d, 4> localRest{};
        double thermalFactor = 1.0;
        double mass = 0;
        double inertiaU = 1;
        double inertiaV = 1;
        double inertiaN = 1;
        bool asleep = false;
    };

    /// Tributary area of each node on the outside of the part: a third of
    /// every boundary triangle it corners. Zero for an interior node, which
    /// is what keeps contact on the real surface. Rebuilt when the skin
    /// changes, so a refine or a detach is picked up on the next step.
    const std::vector<double>& contactAreas() const;
    /// Half the offset a contact sits off a node. A tet volume carries its
    /// own surface, so the node is the surface. A shell is a midsurface and
    /// stands off by half its thickness. Used by floor/arena contact.
    double contactOffset() const { return hasVolume() ? 0.0 : 0.5 * _thickness; }
    /// Part-part penalty thickness. Volumes need a real skin (floor contact
    /// does not): a zero band skips every node-vs-face pair.
    double contactSkin() const {
        if (!hasVolume())
            return 0.5 * _thickness;
        return std::max(1.0e-3, std::min(0.25 * std::max(_thickness, 1.0e-3), 5.0e-3));
    }

    bool nodeRigid(std::size_t index) const;
    void clearLinkedNodes();
    void markLinkedNode(int index);
    std::vector<RigidShard>& shards() { return _shards; }
    const std::vector<RigidShard>& shards() const { return _shards; }
    /// Shell boundary node tied to a tet node so a mixed part stays one solid.
    const std::vector<std::pair<int, int>>& shellSolidTies() const { return _shellSolidTies; }
    /// Drop chips that are still moving, and the skin that belongs only to
    /// them. Settled chips stay. True when something was removed.
    bool discardAwakeShards();
    /// Reset a chip whose pose or velocity went non-finite. Does not delete it.
    void sanitizeShards();
    /// Free motion from gravity, contact reactions, and nodal loads.
    void advanceShards(double dt, const chrono::ChVector3d& gravity);
    /// Floor/arena contact for chips, then copy the rigid pose onto their nodes.
    /// A chip only leaves the solve when its motion has actually died; touching
    /// the floor is not enough.
    void finishShards(bool floorEnabled,
                      double floorHeight,
                      const chrono::ChVector3d& arenaSize,
                      const chrono::ChVector3d& gravity,
                      double restitution,
                      double friction);

private:
    struct QuadMechanics {
        double tension = 0.0;
        double theta = 0.0;
        double vonMises = 0.0;
        double energy = 0.0;
        double area = 0.0;
        double edgeU = 0.0;
        double edgeV = 0.0;
        bool ok = false;
    };

    void buildSurface(const PartDefinition& definition, bool append);
    void buildVolume(const PartDefinition& definition);
    bool releaseMoltenTets(bool includeBrittle = true);
    int xyzMidpoint(int a, int b);
    int xyzOnPlane(int a, int b, const chrono::ChVector3d& origin, const chrono::ChVector3d& normal);
    std::vector<std::size_t> splitTetOnTensionPlane(std::size_t index);
    std::vector<std::size_t> splitTetOnPlane(std::size_t index, const chrono::ChVector3d& origin,
                                             const chrono::ChVector3d& axis,
                                             std::map<std::pair<int, int>, int>* cuts);
    /// Open every tet the crack plane hits. Shared edges share nodes. A child
    /// that would be a sliver is released along its nearest face instead.
    bool fractureCluster(const std::vector<std::size_t>& seeds, std::vector<char>& shard);
    void releaseTetAlongFace(std::size_t index, const chrono::ChVector3d& normal);
    /// Cut tets on the planes between plies so a laminate has real interface faces.
    void stampPlyPlanes();
    void bisectTet(std::size_t index, int a, int b, int mid);
    int cloneNode(int original);
    int insertMidpoint(int a, int b);
    void rebuildTrianglesFromQuads();
    void detachQuad(std::size_t quad);
    bool splitQuad(std::size_t quad, bool alongU);
    /// Halve a quad into two visual pieces without adding a shell element.
    int splitQuadGeometry(std::size_t quad, bool alongU);
    void isolateQuad(int quadIndex);
    void releaseShard(int quadIndex);
    /// Freeze nodes that are off every live tet and every live shell on this
    /// part. A shattered weld partner must not pin the other part's live mesh.
    void freezeUnusedNodes();
    bool shatterQuads();
    /// A piece that is no longer held by a clamp or a weld leaves the implicit
    /// solve. `moltenOnly` keeps a room-temperature tear as a shell and only
    /// rigidizes a fragment that is already above its own solidus.
    void promoteLooseFragments(bool moltenOnly = false);
    chrono::ChVector3d nodeRest(std::size_t index) const;
    double maxEdgeStretch(const int* index, int count) const;
    bool averagePlastic(std::size_t element, ReissnerPlasticData& out) const;
    QuadMechanics quadMechanics(std::size_t quad) const;
    double thermalStrain() const;
    void seedTemperatures();
    void noteNewNode(int index, double kelvin);
    std::shared_ptr<chrono::fea::ChMaterialShellReissner> makeQuadMaterial(std::size_t quad);
    mutable std::vector<double> _contactArea;
    mutable std::size_t _contactAreaNodes = 0;
    mutable std::size_t _contactAreaFaces = 0;
    std::string _name;
    Material _material;
    double _thickness = 0.008;
    double _refineLimit = 0.025;
    double _minElement = 0.002;
    double _temperature = 293.15;
    double _initialTemperature = 293.15;
    double _sharedStiffness = 1.0;
    double _sharedStrength = 1.0;
    bool _fixed = false;
    bool _shattering = false;
    int _shatterLive0 = 0;
    std::shared_ptr<ThermalReissnerElasticity> _elasticity;
    std::shared_ptr<ReissnerJ2> _plasticity;
    std::shared_ptr<chrono::fea::ChMaterialShellReissner> _shellMaterial;
    struct QuadThermal {
        std::shared_ptr<ThermalReissnerElasticity> elasticity;
        std::shared_ptr<ReissnerJ2> plasticity;
        double stiffness = 1.0;
        double strength = 1.0;
    };
    std::vector<QuadThermal> _quadThermal;
    std::vector<double> _nodeTemperature;
    std::vector<double> _strainTemperature;
    std::vector<char> _heatActive;
    double _prescribedThermalKur = 0.0;
    double _thermalStamp = 0.0;
    bool _thermalPrestress = false;
    bool _thermalLoadChanged = false;

    std::shared_ptr<chrono::fea::ChMesh> _mesh;
    std::vector<std::shared_ptr<chrono::fea::ChNodeFEAxyzrot>> _nodes;
    std::vector<std::shared_ptr<chrono::fea::ChNodeFEAxyz>> _xyzNodes;
    std::vector<std::shared_ptr<chrono::fea::ChElementShellReissner4>> _elements;
    std::vector<std::shared_ptr<chrono::fea::ChElementShellBST>> _trianglesElements;
    std::vector<std::shared_ptr<ThermalKirchhoff>> _triangleThermal;
    std::vector<std::shared_ptr<ThermalTet4>> _tetElements;
    std::vector<std::array<int, 4>> _tets;
    std::vector<char> _tetLive;
    int _skinTriangles = 0;
    std::vector<TearTwin> _tearTwins;
    bool _torn = false;
    double _peakFractureTension = 0.0;
    std::vector<double> _nodeMass;
    std::vector<std::array<int, 4>> _quads;
    std::vector<std::array<int, 3>> _triangles;
    std::vector<int> _boundaryNodes;
    std::vector<chrono::ChVector3d> _loads;
    std::vector<char> _rigidNode;
    std::vector<char> _rigidQuad;
    std::vector<char> _linked;
    std::vector<RigidShard> _shards;
    std::vector<chrono::ChVector3d> _savedPos;
    std::vector<chrono::ChVector3d> _savedVel;
    std::vector<chrono::ChVector3d> _savedAcc;
    std::vector<chrono::ChQuaterniond> _savedRot;
    std::vector<chrono::ChQuaterniond> _savedRotDt;
    std::vector<chrono::ChVector3d> _savedAngAcc;
    std::vector<std::pair<int, int>> _shellSolidTies;
};

}  // namespace rcs
