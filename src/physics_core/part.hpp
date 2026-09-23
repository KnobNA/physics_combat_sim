#pragma once

#include "physics_core/material.hpp"
#include "physics_core/reissner_j2.hpp"

#include "chrono/core/ChQuaternion.h"
#include "chrono/core/ChVector3.h"
#include "chrono/fea/ChElementShellBST.h"
#include "chrono/fea/ChElementShellReissner4.h"
#include "chrono/fea/ChMesh.h"
#include "chrono/fea/ChNodeFEAxyz.h"
#include "chrono/fea/ChNodeFEAxyzrot.h"

#include <array>
#include <cstddef>
#include <memory>
#include <string>
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

/// A triangulated surface. Non-empty means this part is a triangle shell
/// (ChElementShellBST) rather than a rectangular Reissner patch.
struct SurfaceMesh {
    std::vector<chrono::ChVector3d> vertices;
    std::vector<std::array<int, 3>> triangles;

    bool empty() const { return vertices.empty() || triangles.empty(); }
};

/// Authoring-time description of one rigidly-manufactured piece of the
/// assembly. Nothing downstream knows this came from a cube.
struct PartDefinition {
    std::string name = "Part";
    Material material;
    double thickness = 0.008;  ///< m
    PlatePatch patch;
    SurfaceMesh surface;
    bool fixed = false;
    /// Hold the patch boundary still so a face can be loaded in bending.
    bool clampBoundary = false;
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
    void commitLoads();

    /// Per-quad peak von Mises stress across the shell thickness, in Pa.
    /// Plastic strain is subtracted, so a dent that has unloaded reads low.
    void evaluateStress(std::vector<float>& out) const;
    void evaluatePlasticStrain(std::vector<float>& out) const;
    /// Highest outer-fiber tensile principal stress, Pa. Brittle fracture uses this.
    double peakTension() const;

    /// Detach elements that have stretched, yielded, or cracked past the
    /// material limit. Returns true when the mesh gained nodes.
    bool tearFailed();

    /// Uniform part temperature, kelvin. Built at the material reference so
    /// the thermal strain is zero until a heat solver writes something else.
    double temperature() const { return _temperature; }
    void setTemperature(double kelvin);

    std::size_t nodeCount() const;
    chrono::ChVector3d nodePosition(std::size_t index) const;
    chrono::ChVector3d nodeVelocity(std::size_t index) const;
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
        double mass = 0;
        double inertiaU = 1;
        double inertiaV = 1;
        double inertiaN = 1;
        bool asleep = false;
    };

    bool nodeRigid(std::size_t index) const;
    void clearLinkedNodes();
    void markLinkedNode(int index);
    std::vector<RigidShard>& shards() { return _shards; }
    const std::vector<RigidShard>& shards() const { return _shards; }
    /// Free motion from gravity, contact reactions, and nodal loads.
    void advanceShards(double dt, const chrono::ChVector3d& gravity);
    /// Sit chips on the floor and copy the rigid pose back onto their nodes.
    void finishShards(bool floorEnabled, double floorHeight);

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

    void buildSurface(const PartDefinition& definition);
    int cloneNode(int original);
    int insertMidpoint(int a, int b);
    void rebuildTrianglesFromQuads();
    void detachQuad(std::size_t quad);
    bool splitQuad(std::size_t quad, bool alongU);
    /// Halve a quad into two visual pieces without adding a shell element.
    int splitQuadGeometry(std::size_t quad, bool alongU);
    void isolateQuad(int quadIndex);
    void releaseShard(int quadIndex);
    bool shatterQuads();
    void promoteLooseFragments();
    chrono::ChVector3d nodeRest(std::size_t index) const;
    double maxEdgeStretch(const int* index, int count) const;
    bool averagePlastic(std::size_t element, ReissnerPlasticData& out) const;
    QuadMechanics quadMechanics(std::size_t quad) const;
    double thermalStrain() const;
    std::string _name;
    Material _material;
    double _thickness = 0.008;
    double _temperature = 293.15;
    bool _fixed = false;
    bool _shattering = false;
    std::shared_ptr<ThermalReissnerElasticity> _elasticity;
    std::shared_ptr<chrono::fea::ChMaterialShellReissner> _shellMaterial;

    std::shared_ptr<chrono::fea::ChMesh> _mesh;
    std::vector<std::shared_ptr<chrono::fea::ChNodeFEAxyzrot>> _nodes;
    std::vector<std::shared_ptr<chrono::fea::ChNodeFEAxyz>> _xyzNodes;
    std::vector<std::shared_ptr<chrono::fea::ChElementShellReissner4>> _elements;
    std::vector<std::shared_ptr<chrono::fea::ChElementShellBST>> _trianglesElements;
    std::vector<TearTwin> _tearTwins;
    bool _torn = false;
    std::vector<double> _nodeMass;
    std::vector<std::array<int, 4>> _quads;
    std::vector<std::array<int, 3>> _triangles;
    std::vector<int> _boundaryNodes;
    std::vector<chrono::ChVector3d> _loads;
    std::vector<char> _rigidNode;
    std::vector<char> _rigidQuad;
    std::vector<char> _linked;
    std::vector<RigidShard> _shards;
};

}  // namespace rcs
