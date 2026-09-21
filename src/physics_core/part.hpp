#pragma once

#include "physics_core/material.hpp"

#include "chrono/core/ChQuaternion.h"
#include "chrono/core/ChVector3.h"
#include "chrono/fea/ChElementShellReissner4.h"
#include "chrono/fea/ChMesh.h"
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

/// Authoring-time description of one rigidly-manufactured piece of the
/// assembly. Nothing downstream knows this came from a cube.
struct PartDefinition {
    std::string name = "Part";
    Material material;
    double thickness = 0.008;  ///< m
    PlatePatch patch;
    bool fixed = false;
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
    void evaluateStress(std::vector<float>& out) const;

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
    const std::vector<double>& nodeMass() const { return _nodeMass; }
    const std::vector<std::array<int, 4>>& quads() const { return _quads; }
    const std::vector<std::array<int, 3>>& triangles() const { return _triangles; }
    /// Nodes on the outline of the patch, the only ones a weld can pair up.
    const std::vector<int>& boundaryNodes() const { return _boundaryNodes; }

private:
    std::string _name;
    Material _material;
    double _thickness = 0.008;
    bool _fixed = false;

    std::shared_ptr<chrono::fea::ChMesh> _mesh;
    std::vector<std::shared_ptr<chrono::fea::ChNodeFEAxyzrot>> _nodes;
    std::vector<std::shared_ptr<chrono::fea::ChElementShellReissner4>> _elements;
    std::vector<double> _nodeMass;
    std::vector<std::array<int, 4>> _quads;
    std::vector<std::array<int, 3>> _triangles;
    std::vector<int> _boundaryNodes;
    std::vector<chrono::ChVector3d> _loads;
};

}  // namespace rcs
