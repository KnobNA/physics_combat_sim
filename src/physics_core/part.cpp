#include "physics_core/part.hpp"
#include "physics_core/physics_units.hpp"
#include "physics_core/reissner_j2.hpp"

#include "chrono/core/ChMatrix33.h"
#include "chrono/fea/ChMaterialShellKirchhoff.h"
#include "chrono/fea/ChMaterialShellReissner.h"

#include <algorithm>
#include <cmath>
#include <limits>
#include <map>
#include <set>
#include <utility>

namespace rcs {
namespace {

using chrono::ChVector3d;

/// Orthonormal frame of a quad from its four corners: e1 runs along the
/// averaged 0->1 direction, n is the outward normal, e2 completes the triad.
bool quadFrame(const ChVector3d p[4], ChVector3d& e1, ChVector3d& e2, ChVector3d& n) {
    ChVector3d u = (p[1] + p[2]) - (p[0] + p[3]);
    ChVector3d v = (p[3] + p[2]) - (p[0] + p[1]);
    n = u.Cross(v);
    if (u.Length() < 1.0e-12 || n.Length() < 1.0e-12)
        return false;
    n.Normalize();
    u.Normalize();
    e1 = u;
    e2 = n.Cross(e1);
    return true;
}

/// Least-squares gradient of a 2-component nodal field over a quad:
/// solves G = (sum f_i a_i^T)(sum a_i a_i^T)^-1 with the centroid removed.
bool fitGradient(const double a[4][2], const double f[4][2], double g[2][2]) {
    double m00 = 0.0, m01 = 0.0, m11 = 0.0;
    double s[2][2] = {{0.0, 0.0}, {0.0, 0.0}};
    for (int i = 0; i < 4; ++i) {
        m00 += a[i][0] * a[i][0];
        m01 += a[i][0] * a[i][1];
        m11 += a[i][1] * a[i][1];
        for (int r = 0; r < 2; ++r) {
            s[r][0] += f[i][r] * a[i][0];
            s[r][1] += f[i][r] * a[i][1];
        }
    }
    const double det = m00 * m11 - m01 * m01;
    if (std::abs(det) < 1.0e-20)
        return false;
    const double i00 = m11 / det;
    const double i01 = -m01 / det;
    const double i11 = m00 / det;
    for (int r = 0; r < 2; ++r) {
        g[r][0] = s[r][0] * i00 + s[r][1] * i01;
        g[r][1] = s[r][0] * i01 + s[r][1] * i11;
    }
    return true;
}

double liveStretchLimit(const Material& material, double ductility) {
    if (material.fracture == FractureMode::Ductile && material.failPlasticStrain > 0.0)
        return ductileTearStretch(material.failPlasticStrain, ductility);
    return 1.0 + std::max(0.0, material.tearStretch - 1.0) * ductility;
}

}  // namespace

void Part::build(const PartDefinition& definition) {
    using namespace chrono;
    using namespace chrono::fea;

    _name = definition.name;
    _material = definition.material;
    _thickness = definition.thickness;
    if (!_material.plies.empty()) {
        double plySum = 0.0;
        for (const Ply& ply : _material.plies)
            plySum += std::max(ply.thickness, 0.0);
        if (plySum > 1.0e-9)
            _thickness = plySum;
    }
    _refineLimit = definition.sizing.refineLimit > 0.0 ? definition.sizing.refineLimit : definition.sizing.maxElement;
    _minElement = definition.sizing.minElement > 0.0 ? definition.sizing.minElement : 0.002;
    _initialTemperature = definition.initialTemperature;
    _temperature = _initialTemperature;
    _fixed = definition.fixed;
    _torn = false;
    _peakFractureTension = 0.0;
    _shattering = false;
    _shatterLive0 = 0;
    _elasticity.reset();
    _plasticity.reset();
    _shellMaterial.reset();
    _quadThermal.clear();
    _nodeTemperature.clear();
    _strainTemperature.clear();
    _heatActive.clear();
    _prescribedThermalKur = 0.0;
    _thermalStamp = 0.0;
    _thermalPrestress = false;
    _thermalLoadChanged = false;
    _tearTwins.clear();
    _xyzNodes.clear();
    _trianglesElements.clear();
    _triangleThermal.clear();
    _tetElements.clear();
    _tets.clear();
    _tetLive.clear();
    _skinTriangles = 0;

    _mesh = chrono_types::make_shared<ChMesh>();
    // Gravity is one of several nodal loads we assemble ourselves each step,
    // so the mesh must not add its own on top.
    _mesh->SetAutomaticGravity(false);

    if (!definition.volume.empty()) {
        buildVolume(definition);
        if (!definition.surface.empty())
            buildSurface(definition, true);
        _loads.assign(nodeCount(), ChVector3d(0, 0, 0));
        seedTemperatures();
        return;
    }

    if (!definition.surface.empty()) {
        buildSurface(definition, false);
        return;
    }

    // A plate shell is the patch midsurface (Reissner), even when the CAD
    // solid is still attached for display. Meshing that solid as BST triangles
    // leaves the weld as a point constraint Newton cannot satisfy.
    const bool patchShell = definition.sizing.throughThickness <= 0 && definition.patch.uSize > 1.0e-6 &&
                            definition.patch.vSize > 1.0e-6;
    if (!patchShell && !definition.source.empty()) {
        PartDefinition display = definition;
        display.surface = definition.source;
        if (display.surface.thickness.empty())
            display.surface.thickness.assign(display.surface.triangles.size(), _thickness);
        buildSurface(display, false);
        return;
    }

    _elasticity = chrono_types::make_shared<ThermalReissnerElasticity>(
        _material.youngsModulus, _material.poissonRatio);
    _elasticity->SetThermalStrain(thermalStrain());
    auto damping = chrono_types::make_shared<ThermalRayleighDamping>(
        _elasticity, _material.dampingBeta);
    _plasticity.reset();
    if (_material.allowPlasticity && _material.yieldStress > 0.0)
        _plasticity = chrono_types::make_shared<ReissnerJ2>(
            _material.yieldStress, _material.youngsModulus, _material.ultimateStress,
            _material.failPlasticStrain);
    _shellMaterial = chrono_types::make_shared<ChMaterialShellReissner>(
        _elasticity, _plasticity, damping);
    _shellMaterial->SetDensity(_material.density);

    const PlatePatch& patch = definition.patch;
    const int nu = std::max(1, patch.divU);
    const int nv = std::max(1, patch.divV);
    const int stride = nv + 1;

    ChVector3d uDir = patch.uAxis;
    ChVector3d vDir = patch.vAxis;
    uDir.Normalize();
    vDir.Normalize();
    ChVector3d nDir = uDir.Cross(vDir);
    nDir.Normalize();
    ChMatrix33d nodeRot;
    nodeRot.SetFromAxisZ(nDir, uDir);

    _nodes.clear();
    _elements.clear();
    _quads.clear();
    _triangles.clear();
    _boundaryNodes.clear();
    _nodeMass.assign(static_cast<std::size_t>(nu + 1) * stride, 0.0);

    for (int iu = 0; iu <= nu; ++iu) {
        for (int iv = 0; iv <= nv; ++iv) {
            const double su = static_cast<double>(iu) / nu;
            const double sv = static_cast<double>(iv) / nv;
            const ChVector3d pos =
                patch.origin + uDir * (su * patch.uSize) + vDir * (sv * patch.vSize);
            auto node = chrono_types::make_shared<ChNodeFEAxyzrot>(ChFrame<>(pos, nodeRot));
            // The Reissner element mass matrix already lumps translational and
            // rotary inertia; leaving nodal values at zero avoids counting twice.
            node->SetMass(0);
            node->GetInertia().fillDiagonal(0);
            node->SetFixed(_fixed);
            _mesh->AddNode(node);
            _nodes.push_back(node);
            if (iu == 0 || iu == nu || iv == 0 || iv == nv)
                _boundaryNodes.push_back(iu * stride + iv);
        }
    }

    const double du = patch.uSize / nu;
    const double dv = patch.vSize / nv;
    const double quadMass = _material.density * _thickness * du * dv;

    for (int iu = 1; iu <= nu; ++iu) {
        for (int iv = 1; iv <= nv; ++iv) {
            const int i00 = (iu - 1) * stride + (iv - 1);
            const int i10 = iu * stride + (iv - 1);
            const int i11 = iu * stride + iv;
            const int i01 = (iu - 1) * stride + iv;

            auto element = chrono_types::make_shared<ChElementShellReissner4>();
            element->SetNodes(_nodes[i00], _nodes[i10], _nodes[i11], _nodes[i01]);
            element->AddLayer(_thickness, 0, makeQuadMaterial(_elements.size()));
            _mesh->AddElement(element);
            _elements.push_back(element);

            _quads.push_back({i00, i10, i11, i01});
            _triangles.push_back({i00, i10, i11});
            _triangles.push_back({i00, i11, i01});

            for (int index : {i00, i10, i11, i01})
                _nodeMass[index] += quadMass * 0.25;
        }
    }

    if (definition.clampBoundary) {
        double minX = 1.0e9, maxX = -1.0e9;
        double minY = 1.0e9, maxY = -1.0e9;
        double minZ = 1.0e9, maxZ = -1.0e9;
        for (const auto& node : _nodes) {
            const ChVector3d p = node->GetPos();
            minX = std::min(minX, p.x());
            maxX = std::max(maxX, p.x());
            minY = std::min(minY, p.y());
            maxY = std::max(maxY, p.y());
            minZ = std::min(minZ, p.z());
            maxZ = std::max(maxZ, p.z());
        }
        const double extX = maxX - minX;
        const double extY = maxY - minY;
        const double extZ = maxZ - minZ;
        int thick = 0;
        double thickLen = extX;
        if (extY < thickLen) {
            thick = 1;
            thickLen = extY;
        }
        if (extZ < thickLen)
            thick = 2;
        const int axisA = (thick + 1) % 3;
        const int axisB = (thick + 2) % 3;
        const double inA = axisA == 0 ? extX : (axisA == 1 ? extY : extZ);
        const double inB = axisB == 0 ? extX : (axisB == 1 ? extY : extZ);
        const double aspect = std::max(inA, inB) / std::max(std::min(inA, inB), 1.0e-9);
        // Same rule as the tet clamp: a square face is a plate (all four
        // silhouette edges), a long bar is a cantilever (min-X root).
        const bool plate = aspect < 2.0;
        auto coord = [](const ChVector3d& p, int axis) {
            return axis == 0 ? p.x() : (axis == 1 ? p.y() : p.z());
        };
        auto lo = [&](int axis) { return axis == 0 ? minX : (axis == 1 ? minY : minZ); };
        auto hi = [&](int axis) { return axis == 0 ? maxX : (axis == 1 ? maxY : maxZ); };
        for (int index : _boundaryNodes) {
            if (index < 0 || static_cast<std::size_t>(index) >= _nodes.size())
                continue;
            const ChVector3d p = _nodes[static_cast<std::size_t>(index)]->GetPos();
            bool clamp = false;
            if (plate) {
                if (coord(p, axisA) <= lo(axisA) + 1.0e-4 || coord(p, axisA) >= hi(axisA) - 1.0e-4)
                    clamp = true;
                if (coord(p, axisB) <= lo(axisB) + 1.0e-4 || coord(p, axisB) >= hi(axisB) - 1.0e-4)
                    clamp = true;
            } else if (p.x() <= minX + 1.0e-4) {
                clamp = true;
            }
            if (clamp)
                _nodes[static_cast<std::size_t>(index)]->SetFixed(true);
        }
    }

    _loads.assign(_nodes.size(), ChVector3d(0, 0, 0));
    seedTemperatures();
}

void Part::transform(const chrono::ChVector3d& translation,
                     const chrono::ChQuaterniond& rotation,
                     const chrono::ChVector3d& pivot,
                     bool alsoRest)
{
    if (!_xyzNodes.empty()) {
        for (auto& node : _xyzNodes) {
            const chrono::ChVector3d rel = node->GetPos() - pivot;
            node->SetPos(pivot + rotation.Rotate(rel) + translation);
            if (alsoRest)
                node->SetX0(node->GetPos());
        }
        if (alsoRest) {
            for (const auto& tet : _tetElements)
                if (tet)
                    tet->refreshRest();
        }
        return;
    }
    for (auto& node : _nodes) {
        const chrono::ChVector3d rel = node->GetPos() - pivot;
        node->SetPos(pivot + rotation.Rotate(rel) + translation);
        node->SetRot(rotation * node->GetRot());
        if (alsoRest)
            node->SetX0(node->Frame());
        node->ForceToRest();
    }
}

void Part::setVelocity(const chrono::ChVector3d& linear,
                       const chrono::ChVector3d& angular,
                       const chrono::ChVector3d& pivot)
{
    if (!_xyzNodes.empty()) {
        for (auto& node : _xyzNodes)
            node->SetPosDt(linear + angular.Cross(node->GetPos() - pivot));
        return;
    }
    for (auto& node : _nodes) {
        node->SetPosDt(linear + angular.Cross(node->GetPos() - pivot));
        node->SetAngVelParent(angular);
    }
}

void Part::clearLoads() {
    std::fill(_loads.begin(), _loads.end(), chrono::ChVector3d(0, 0, 0));
}

void Part::addLoad(std::size_t nodeIndex, const chrono::ChVector3d& force) {
    if (nodeIndex < _loads.size())
        _loads[nodeIndex] += force;
}

chrono::ChVector3d Part::nodeForce(std::size_t nodeIndex) const {
    return nodeIndex < _loads.size() ? _loads[nodeIndex] : chrono::ChVector3d(0, 0, 0);
}

bool Part::nodeOnLiveShell(std::size_t index) const {
    const int node = static_cast<int>(index);
    for (std::size_t q = 0; q < _quads.size() && q < _elements.size(); ++q) {
        if (!_elements[q])
            continue;
        if (q < _rigidQuad.size() && _rigidQuad[q])
            continue;
        for (int corner : _quads[q]) {
            if (corner == node)
                return true;
        }
    }
    for (std::size_t t = 0; t < _trianglesElements.size(); ++t) {
        if (!_trianglesElements[t])
            continue;
        const std::size_t face = t + static_cast<std::size_t>(std::max(0, _skinTriangles));
        if (face >= _triangles.size())
            continue;
        for (int corner : _triangles[face]) {
            if (corner == node)
                return true;
        }
    }
    return false;
}

void Part::addStructuralLoad(std::size_t nodeIndex, const chrono::ChVector3d& force) {
    if (nodeOnLiveShell(nodeIndex)) {
        addLoad(nodeIndex, force);
        return;
    }
    std::vector<int> clones;
    for (const auto& twin : _tearTwins) {
        if (twin.original == static_cast<int>(nodeIndex))
            clones.push_back(twin.clone);
    }
    if (clones.empty()) {
        addLoad(nodeIndex, force);
        return;
    }
    const chrono::ChVector3d share = force / static_cast<double>(clones.size());
    for (int clone : clones)
        addLoad(static_cast<std::size_t>(clone), share);
}

const std::vector<double>& Part::contactAreas() const {
    const std::size_t triFaces = _trianglesElements.empty()
                                     ? static_cast<std::size_t>(std::max(0, _skinTriangles))
                                     : _triangles.size();
    const std::size_t faces = _quads.empty() ? triFaces : _quads.size();
    if (_contactAreaNodes == nodeCount() && _contactAreaFaces == faces && !_contactArea.empty())
        return _contactArea;
    _contactArea.assign(nodeCount(), 0.0);
    _contactAreaNodes = nodeCount();
    _contactAreaFaces = faces;
    auto share = [&](int index, double area) {
        if (index >= 0 && static_cast<std::size_t>(index) < _contactArea.size())
            _contactArea[static_cast<std::size_t>(index)] += area;
    };
    if (_quads.empty()) {
        for (std::size_t f = 0; f < faces && f < _triangles.size(); ++f) {
            const auto& tri = _triangles[f];
            if (tri[0] < 0 || tri[1] < 0 || tri[2] < 0)
                continue;
            const ChVector3d a = nodePosition(static_cast<std::size_t>(tri[0]));
            const ChVector3d b = nodePosition(static_cast<std::size_t>(tri[1]));
            const ChVector3d c = nodePosition(static_cast<std::size_t>(tri[2]));
            const double area = 0.5 * Vcross(b - a, c - a).Length() / 3.0;
            for (int corner : tri)
                share(corner, area);
        }
    } else {
        for (const auto& quad : _quads) {
            if (quad[0] < 0 || quad[1] < 0 || quad[2] < 0 || quad[3] < 0)
                continue;
            const ChVector3d a = nodePosition(static_cast<std::size_t>(quad[0]));
            const ChVector3d b = nodePosition(static_cast<std::size_t>(quad[1]));
            const ChVector3d c = nodePosition(static_cast<std::size_t>(quad[2]));
            const ChVector3d d = nodePosition(static_cast<std::size_t>(quad[3]));
            const double area =
                0.25 * (0.5 * Vcross(b - a, c - a).Length() + 0.5 * Vcross(c - a, d - a).Length());
            for (int corner : quad)
                share(corner, area);
        }
    }
    return _contactArea;
}

void Part::freezeUnusedNodes() {
    auto live = [&](std::size_t i) {
        if (hasVolume() && nodeOnLiveTet(i))
            return true;
        return nodeOnLiveShell(i);
    };
    const std::size_t count = nodeCount();
    if (_rigidNode.size() < count)
        _rigidNode.resize(count, 0);
    for (std::size_t i = 0; i < _xyzNodes.size(); ++i) {
        if (live(i) || !_xyzNodes[i] || _xyzNodes[i]->IsFixed())
            continue;
        _xyzNodes[i]->SetFixed(true);
        _xyzNodes[i]->SetPosDt(chrono::ChVector3d(0, 0, 0));
        _rigidNode[i] = 1;
    }
    for (std::size_t i = 0; i < _nodes.size(); ++i) {
        if (live(i) || !_nodes[i] || _nodes[i]->IsFixed())
            continue;
        _nodes[i]->SetFixed(true);
        _nodes[i]->SetPosDt(chrono::ChVector3d(0, 0, 0));
        _nodes[i]->SetAngVelParent(chrono::ChVector3d(0, 0, 0));
        _rigidNode[i] = 1;
    }
}

void Part::commitLoads() {
    if (!_xyzNodes.empty()) {
        for (std::size_t i = 0; i < _xyzNodes.size() && i < _loads.size(); ++i) {
            if (i < _rigidNode.size() && _rigidNode[i] && _xyzNodes[i]) {
                _xyzNodes[i]->SetFixed(true);
                _xyzNodes[i]->SetForce(chrono::ChVector3d(0, 0, 0));
                continue;
            }
            _xyzNodes[i]->SetForce(_loads[i]);
        }
        return;
    }
    for (std::size_t i = 0; i < _nodes.size(); ++i)
        _nodes[i]->SetForce(_loads[i]);
}

void Part::enforceRigidNodes() {
    for (std::size_t i = 0; i < _xyzNodes.size(); ++i) {
        if (i < _rigidNode.size() && _rigidNode[i] && _xyzNodes[i]) {
            _xyzNodes[i]->SetFixed(true);
            _xyzNodes[i]->SetPosDt(chrono::ChVector3d(0, 0, 0));
        }
    }
}

double Part::thermalStrain() const {
    return _material.thermalExpansion * (_temperature - _initialTemperature);
}

void Part::seedTemperatures() {
    _nodeTemperature.assign(nodeCount(), _initialTemperature);
    _strainTemperature.assign(nodeCount(), _initialTemperature);
    _heatActive.assign(nodeCount(), 0);
    _sharedStiffness = 1.0;
    _sharedStrength = 1.0;
    _temperature = _initialTemperature;
}

void Part::noteNewNode(int index, double kelvin) {
    if (index < 0)
        return;
    const std::size_t slot = static_cast<std::size_t>(index);
    if (_nodeTemperature.size() <= slot)
        _nodeTemperature.resize(slot + 1, _initialTemperature);
    if (_heatActive.size() <= slot)
        _heatActive.resize(slot + 1, 0);
    if (_strainTemperature.size() <= slot)
        _strainTemperature.resize(slot + 1, _initialTemperature);
    _nodeTemperature[slot] = kelvin;
    _strainTemperature[slot] = kelvin;
}

double Part::nodeTemperature(std::size_t index) const {
    if (index >= _nodeTemperature.size())
        return _initialTemperature;
    return _nodeTemperature[index];
}

void Part::setNodeTemperature(std::size_t index, double kelvin) {
    if (!std::isfinite(kelvin))
        return;
    if (index >= _nodeTemperature.size())
        _nodeTemperature.resize(index + 1, _initialTemperature);
    _nodeTemperature[index] = kelvin;
}

double Part::strainTemperature(std::size_t index) const {
    if (index >= _strainTemperature.size())
        return nodeTemperature(index);
    return _strainTemperature[index];
}

void Part::setStrainTemperature(std::size_t index, double kelvin) {
    if (!std::isfinite(kelvin))
        return;
    if (index >= _strainTemperature.size())
        _strainTemperature.resize(index + 1, _initialTemperature);
    _strainTemperature[index] = kelvin;
}

bool Part::heatActive(std::size_t index) const {
    return index < _heatActive.size() && _heatActive[index];
}

void Part::setHeatActive(std::size_t index, bool active) {
    if (index >= _heatActive.size())
        _heatActive.resize(index + 1, 0);
    _heatActive[index] = active ? 1 : 0;
}

void Part::setThermalCurvature(double kappa) {
    if (std::abs(kappa - _prescribedThermalKur) > 1.0e-12)
        _thermalLoadChanged = true;
    _prescribedThermalKur = kappa;
    if (std::abs(kappa) > 1.0e-12)
        _thermalPrestress = true;
    if (_elasticity)
        _elasticity->SetThermalCurvature(kappa);
    for (QuadThermal& slot : _quadThermal) {
        if (slot.elasticity)
            slot.elasticity->SetThermalCurvature(kappa);
    }
}

void Part::setTemperature(double kelvin) {
    _temperature = kelvin;
    for (double& value : _nodeTemperature)
        value = kelvin;
    for (double& value : _strainTemperature)
        value = kelvin;
    if (_elasticity)
        _elasticity->SetThermalStrain(thermalStrain());
    syncThermalMechanics();
}

std::shared_ptr<chrono::fea::ChMaterialShellReissner> Part::makeQuadMaterial(std::size_t quad) {
    using namespace chrono;
    using namespace chrono::fea;
    if (_quadThermal.size() <= quad)
        _quadThermal.resize(quad + 1);
    auto elasticity = chrono_types::make_shared<ThermalReissnerElasticity>(
        _material.youngsModulus, _material.poissonRatio);
    auto damping = chrono_types::make_shared<ThermalRayleighDamping>(
        elasticity, _material.dampingBeta);
    std::shared_ptr<ReissnerJ2> plasticity;
    if (_material.allowPlasticity && _material.yieldStress > 0.0)
        plasticity = chrono_types::make_shared<ReissnerJ2>(
            _material.yieldStress, _material.youngsModulus, _material.ultimateStress,
            _material.failPlasticStrain);
    auto material = chrono_types::make_shared<ChMaterialShellReissner>(
        elasticity, plasticity, damping);
    material->SetDensity(_material.density);
    _quadThermal[quad].elasticity = elasticity;
    _quadThermal[quad].plasticity = plasticity;
    _quadThermal[quad].stiffness = 1.0;
    _quadThermal[quad].strength = 1.0;
    return material;
}

void Part::syncThermalMechanics() {
    if (_nodeTemperature.empty())
        return;
    double sum = 0.0;
    for (double value : _nodeTemperature)
        sum += value;
    _temperature = sum / static_cast<double>(_nodeTemperature.size());

    const double alpha = _material.thermalExpansion;
    const double t0 = _initialTemperature;
    // Expansion follows the lagged strain temperature so a sudden heat-up
    // does not inject the full thermal strain in one step. Stiffness and
    // yield follow the node temperature the colour is already showing.
    const double firstStrain = strainTemperature(0);
    const double firstNode = nodeTemperature(0);
    bool uniform = true;
    for (std::size_t n = 0; n < _nodeTemperature.size(); ++n) {
        if (std::abs(strainTemperature(n) - firstStrain) > 0.05 ||
            std::abs(nodeTemperature(n) - firstNode) > 0.05) {
            uniform = false;
            break;
        }
    }

    auto applyScale = [&](ThermalReissnerElasticity& elasticity, ReissnerJ2* plasticity, double& storedStiffness,
                          double& storedStrength, double kelvin) {
        const double stiffness = thermalStiffnessRetention(_material, kelvin);
        const double strength = thermalStrengthRetention(_material, kelvin);
        auto drifted = [](double next, double previous) {
            return std::abs(next - previous) > 0.002 * std::max(previous, kThermalFloor);
        };
        if (!drifted(stiffness, storedStiffness) && !drifted(strength, storedStrength))
            return;
        elasticity.SetStiffnessScale(stiffness);
        if (plasticity)
            plasticity->SetYieldScale(strength);
        storedStiffness = stiffness;
        storedStrength = strength;
    };

    double thermalKur = 0.0;
    if (_elasticity && !_nodes.empty() && _thickness > 1.0e-9) {
        ChVector3d nDir(0.0, 1.0, 0.0);
        if (!_quads.empty() && _quads[0][0] >= 0 && _quads[0][1] >= 0 && _quads[0][3] >= 0) {
            const ChVector3d p0 = restPosition(static_cast<std::size_t>(_quads[0][0]));
            const ChVector3d p1 = restPosition(static_cast<std::size_t>(_quads[0][1]));
            const ChVector3d p3 = restPosition(static_cast<std::size_t>(_quads[0][3]));
            nDir = (p1 - p0).Cross(p3 - p0);
            if (nDir.Length() > 1.0e-12)
                nDir.Normalize();
            else
                nDir = {0.0, 1.0, 0.0};
        }
        double sumS = 0.0, sumT = 0.0, sumSS = 0.0, sumST = 0.0;
        double sMin = 1.0e30, sMax = -1.0e30;
        int samples = 0;
        for (std::size_t i = 0; i < _nodes.size(); ++i) {
            const double s = restPosition(i).Dot(nDir);
            const double kelvin = strainTemperature(i);
            sMin = std::min(sMin, s);
            sMax = std::max(sMax, s);
            sumS += s;
            sumT += kelvin;
            sumSS += s * s;
            sumST += s * kelvin;
            ++samples;
        }
        const double span = sMax - sMin;
        const double denom = static_cast<double>(samples) * sumSS - sumS * sumS;
        if (samples >= 3 && span > 0.2 * _thickness && std::abs(denom) > 1.0e-18)
            thermalKur = alpha * (static_cast<double>(samples) * sumST - sumS * sumT) / denom;
        else
            thermalKur = _prescribedThermalKur;
    }

    if (uniform && _elasticity) {
        _elasticity->SetThermalStrain(alpha * (firstStrain - t0));
        _elasticity->SetThermalCurvature(thermalKur);
        applyScale(*_elasticity, _plasticity.get(), _sharedStiffness, _sharedStrength, firstNode);
        for (QuadThermal& slot : _quadThermal) {
            if (!slot.elasticity)
                continue;
            slot.elasticity->SetThermalStrain(alpha * (firstStrain - t0));
            slot.elasticity->SetThermalCurvature(thermalKur);
            applyScale(*slot.elasticity, slot.plasticity.get(), slot.stiffness, slot.strength, firstNode);
        }
    } else if (_elasticity) {
        _elasticity->SetThermalStrain(0.0);
        _elasticity->SetThermalCurvature(thermalKur);
        for (std::size_t q = 0; q < _quads.size() && q < _elements.size(); ++q) {
            if (!_elements[q])
                continue;
            double strainKelvin = 0.0;
            double nodeKelvin = 0.0;
            int count = 0;
            for (int index : _quads[q]) {
                if (index < 0)
                    continue;
                const std::size_t node = static_cast<std::size_t>(index);
                strainKelvin += strainTemperature(node);
                nodeKelvin += nodeTemperature(node);
                ++count;
            }
            if (count == 0)
                continue;
            strainKelvin /= count;
            nodeKelvin /= count;
            if (q >= _quadThermal.size() || !_quadThermal[q].elasticity)
                continue;
            _quadThermal[q].elasticity->SetThermalStrain(alpha * (strainKelvin - t0));
            _quadThermal[q].elasticity->SetThermalCurvature(thermalKur);
            applyScale(*_quadThermal[q].elasticity, _quadThermal[q].plasticity.get(), _quadThermal[q].stiffness,
                       _quadThermal[q].strength, nodeKelvin);
        }
    }

    for (RigidShard& shard : _shards) {
        double kelvin = 0.0;
        int count = 0;
        for (int index : shard.nodes) {
            if (index < 0)
                continue;
            kelvin += strainTemperature(static_cast<std::size_t>(index));
            ++count;
        }
        if (count == 0)
            continue;
        kelvin /= count;
        const double factor = std::max(0.5, 1.0 + alpha * (kelvin - t0));
        if (std::abs(factor - shard.thermalFactor) < 1.0e-5)
            continue;
        double rest = 0.0;
        for (const auto& offset : shard.localRest)
            rest += offset.Length2();
        if (rest < 1.0e-16)
            continue;
        for (int k = 0; k < 4; ++k)
            shard.local[static_cast<std::size_t>(k)] = shard.localRest[static_cast<std::size_t>(k)] * factor;
        shard.thermalFactor = factor;
    }

    const double alphaScale = _material.thermalExpansion;
    for (std::size_t t = 0; t < _tetElements.size(); ++t) {
        if (!_tetElements[t] || t >= _tets.size())
            continue;
        double strainKelvin = 0.0;
        double nodeKelvin = 0.0;
        int count = 0;
        for (int index : _tets[t]) {
            if (index < 0)
                continue;
            strainKelvin += strainTemperature(static_cast<std::size_t>(index));
            nodeKelvin += nodeTemperature(static_cast<std::size_t>(index));
            ++count;
        }
        if (count == 0)
            continue;
        strainKelvin /= count;
        nodeKelvin /= count;
        _tetElements[t]->setThermalState(alphaScale * (strainKelvin - t0),
                                         thermalStiffnessRetention(_material, nodeKelvin),
                                         thermalStrengthRetention(_material, nodeKelvin));
    }
    for (std::size_t s = 0; s < _triangleThermal.size(); ++s) {
        if (!_triangleThermal[s])
            continue;
        const std::size_t tri = static_cast<std::size_t>(_skinTriangles) + s;
        if (tri >= _triangles.size())
            continue;
        double strainKelvin = 0.0;
        double nodeKelvin = 0.0;
        int count = 0;
        for (int index : _triangles[tri]) {
            if (index < 0)
                continue;
            strainKelvin += strainTemperature(static_cast<std::size_t>(index));
            nodeKelvin += nodeTemperature(static_cast<std::size_t>(index));
            ++count;
        }
        if (count == 0)
            continue;
        _triangleThermal[s]->setThermalStrain(alphaScale * ((strainKelvin / count) - t0));
        _triangleThermal[s]->setScales(thermalStiffnessRetention(_material, nodeKelvin / count),
                                       thermalStrengthRetention(_material, nodeKelvin / count));
    }

    double stamp = thermalKur;
    _thermalPrestress = std::abs(thermalKur) > 1.0e-12;
    for (std::size_t n = 0; n < _strainTemperature.size(); ++n) {
        stamp += strainTemperature(n);
        if (std::abs(strainTemperature(n) - t0) > 0.05)
            _thermalPrestress = true;
    }
    if (std::abs(stamp - _thermalStamp) > 1.0e-4) {
        _thermalLoadChanged = true;
        _thermalStamp = stamp;
    }
}

void Part::evaluateTemperature(std::vector<float>& out) const {
    const std::size_t count = _quads.empty() ? _triangles.size() : _quads.size();
    out.assign(count, static_cast<float>(_initialTemperature));
    auto mean = [&](const int* index, int n) {
        double sum = 0.0;
        int used = 0;
        for (int i = 0; i < n; ++i) {
            if (index[i] < 0)
                continue;
            sum += nodeTemperature(static_cast<std::size_t>(index[i]));
            ++used;
        }
        return used ? static_cast<float>(sum / used) : static_cast<float>(_initialTemperature);
    };
    if (!_quads.empty()) {
        for (std::size_t q = 0; q < _quads.size() && q < out.size(); ++q)
            out[q] = mean(_quads[q].data(), 4);
    } else {
        for (std::size_t t = 0; t < _triangles.size() && t < out.size(); ++t)
            out[t] = mean(_triangles[t].data(), 3);
    }
}

Part::QuadMechanics Part::quadMechanics(std::size_t quad) const {
    using namespace chrono;

    QuadMechanics out;
    if (quad >= _quads.size() || _nodes.empty())
        return out;
    const auto& idx = _quads[quad];

    ChVector3d rest[4], now[4];
    ChQuaterniond restRot[4], nowRot[4];
    for (int i = 0; i < 4; ++i) {
        const auto& node = _nodes[idx[i]];
        rest[i] = node->GetX0().GetPos();
        restRot[i] = node->GetX0().GetRot();
        now[i] = node->GetPos();
        nowRot[i] = node->GetRot();
    }
    out.edgeU = (rest[1] - rest[0]).Length();
    out.edgeV = (rest[3] - rest[0]).Length();
    out.area = 0.5 * (rest[2] - rest[0]).Cross(rest[3] - rest[1]).Length();

    ChVector3d e1r, e2r, nr, e1c, e2c, nc;
    if (!quadFrame(rest, e1r, e2r, nr) || !quadFrame(now, e1c, e2c, nc))
        return out;

    const double E = _material.youngsModulus;
    const double nu = _material.poissonRatio;
    if (E <= 0.0)
        return out;
    const double cm = E / (1.0 - nu * nu);
    const double shear = E / (1.0 + nu);
    const double zb = 0.5 * _thickness;

    ChVector3d restCentre = 0.25 * (rest[0] + rest[1] + rest[2] + rest[3]);
    ChVector3d nowCentre = 0.25 * (now[0] + now[1] + now[2] + now[3]);

    double a[4][2], x[4][2];
    for (int i = 0; i < 4; ++i) {
        const ChVector3d dr = rest[i] - restCentre;
        const ChVector3d dc = now[i] - nowCentre;
        a[i][0] = dr.Dot(e1r);
        a[i][1] = dr.Dot(e2r);
        x[i][0] = dc.Dot(e1c);
        x[i][1] = dc.Dot(e2c);
    }

    double F[2][2];
    if (!fitGradient(a, x, F))
        return out;

    // In-plane Green-Lagrange strain, so a large rigid rotation of the part
    // registers as zero stress. Subtract this quad's own thermal strain, the
    // one the solver has actually applied. Using the part average paints a
    // stress that the mesh has not expanded into yet.
    double cornerT = 0.0;
    int cornerCount = 0;
    for (int index : idx) {
        if (index < 0)
            continue;
        cornerT += strainTemperature(static_cast<std::size_t>(index));
        ++cornerCount;
    }
    const double eth = cornerCount
                           ? _material.thermalExpansion * (cornerT / cornerCount - _initialTemperature)
                           : 0.0;
    double exx = 0.5 * (F[0][0] * F[0][0] + F[1][0] * F[1][0] - 1.0) - eth;
    double eyy = 0.5 * (F[0][1] * F[0][1] + F[1][1] * F[1][1] - 1.0) - eth;
    double exy = 0.5 * (F[0][0] * F[0][1] + F[1][0] * F[1][1]);
    ReissnerPlasticData plastic;
    const bool yielded = averagePlastic(quad, plastic);
    if (yielded) {
        exx -= plastic.epsU.x();
        eyy -= plastic.epsV.y();
        exy -= 0.5 * (plastic.epsU.y() + plastic.epsV.x());
    }

    const double mxx = cm * (exx + nu * eyy);
    const double myy = cm * (eyy + nu * exx);
    const double mxy = shear * exy;

    const ChQuaterniond refNow = nowRot[0].GetConjugate();
    const ChQuaterniond refRest = restRot[0].GetConjugate();
    double w[4][2];
    for (int i = 0; i < 4; ++i) {
        const ChQuaterniond relative = (refNow * nowRot[i]) * (refRest * restRot[i]).GetConjugate();
        const ChVector3d omega = nowRot[0].Rotate(relative.GetRotVec());
        w[i][0] = omega.Dot(e1c);
        w[i][1] = omega.Dot(e2c);
    }

    double bxx = 0.0, byy = 0.0, bxy = 0.0;
    double G[2][2];
    if (fitGradient(a, w, G)) {
        double kxx = G[1][0];
        double kyy = -G[0][1];
        double kxy = 0.5 * (G[1][1] - G[0][0]);
        if (yielded) {
            kxx -= plastic.kurU.y();
            kyy -= plastic.kurV.x();
            kxy -= 0.5 * (plastic.kurU.x() + plastic.kurV.y());
        }
        bxx = cm * zb * (kxx + nu * kyy);
        byy = cm * zb * (kyy + nu * kxx);
        bxy = shear * zb * kxy;
    }

    auto principal = [](double sxx, double syy, double sxy) {
        const double avg = 0.5 * (sxx + syy);
        const double rad =
            std::sqrt(std::max(0.0, 0.25 * (sxx - syy) * (sxx - syy) + sxy * sxy));
        return std::pair<double, double>{avg + rad, 0.5 * std::atan2(2.0 * sxy, sxx - syy)};
    };
    const auto topP = principal(mxx + bxx, myy + byy, mxy + bxy);
    const auto botP = principal(mxx - bxx, myy - byy, mxy - bxy);
    const auto& face = topP.first >= botP.first ? topP : botP;
    const double top = vonMisesPlaneStress(mxx + bxx, myy + byy, mxy + bxy);
    const double bottom = vonMisesPlaneStress(mxx - bxx, myy - byy, mxy - bxy);
    out.vonMises = std::max(top, bottom);
    out.tension = std::max(0.0, face.first);
    out.theta = face.second;
    out.energy = 0.5 * out.tension * out.tension / E * out.area * _thickness;
    out.ok = std::isfinite(out.vonMises) && std::isfinite(out.tension);
    if (!out.ok)
        out.vonMises = out.tension = 0.0;
    return out;
}

void Part::evaluateStress(std::vector<float>& out) const {
    using namespace chrono;

    out.assign(_quads.size(), 0.0f);

    const double E = _material.youngsModulus;
    const double nu = _material.poissonRatio;
    const double cm = E / (1.0 - nu * nu);
    const double shear = E / (1.0 + nu);

    for (std::size_t q = 0; q < _quads.size(); ++q) {
        // A rigid chip has already spent its elastic energy. Colouring it from
        // the old shell strain leaves a resting pile looking loaded.
        if (q < _rigidQuad.size() && _rigidQuad[q]) {
            out[q] = 0.0f;
            continue;
        }
        const QuadMechanics fiber = quadMechanics(q);
        // The colour is computed with room-temperature modulus. The solver has
        // already scaled that modulus with temperature, so match it or a soft
        // hot face looks like it is carrying hundreds of MPa.
        const double scale = q < _quadThermal.size() ? _quadThermal[q].stiffness : 1.0;
        out[q] = static_cast<float>(fiber.vonMises * scale);
    }

    double tetPeak = 0.0;
    for (std::size_t t = 0; t < _tetElements.size(); ++t) {
        if (!_tetElements[t] || (t < _tetLive.size() && !_tetLive[t]))
            continue;
        const double scale = thermalStiffnessRetention(_material, _temperature);
        tetPeak = std::max(tetPeak, _tetElements[t]->vonMises() * scale);
    }
    if (!_quads.empty()) {
        if (!out.empty() && tetPeak > 0.0)
            out[0] = std::max(out[0], static_cast<float>(tetPeak));
        return;
    }
    if (_triangles.empty()) {
        if (tetPeak > 0.0)
            out.assign(1, static_cast<float>(tetPeak));
        return;
    }

    out.assign(_triangles.size(), 0.0f);
    std::map<std::array<int, 3>, std::size_t> faceHost;
    for (std::size_t tet = 0; tet < _tets.size(); ++tet) {
        if (tet >= _tetElements.size() || !_tetElements[tet] || (tet < _tetLive.size() && !_tetLive[tet]))
            continue;
        const int faces[4][3] = {{0, 1, 2}, {0, 1, 3}, {0, 2, 3}, {1, 2, 3}};
        for (const auto& face : faces) {
            std::array<int, 3> key = {_tets[tet][face[0]], _tets[tet][face[1]], _tets[tet][face[2]]};
            std::sort(key.begin(), key.end());
            faceHost.emplace(key, tet);
        }
    }
    auto paint = [&](std::size_t t) {
        std::array<int, 3> key = _triangles[t];
        std::sort(key.begin(), key.end());
        const auto host = faceHost.find(key);
        if (host == faceHost.end())
            return;
        const auto& element = _tetElements[host->second];
        const double scale = thermalStiffnessRetention(_material, _temperature);
        out[t] = static_cast<float>(element->vonMises() * scale);
    };
    for (std::size_t t = 0; t < _triangles.size(); ++t)
        paint(t);
    if (!_tetElements.empty())
        return;
    for (std::size_t t = static_cast<std::size_t>(_skinTriangles); t < _triangles.size(); ++t) {
        const auto& idx = _triangles[t];
        const ChVector3d rest[3] = {nodeRest(idx[0]), nodeRest(idx[1]), nodeRest(idx[2])};
        const ChVector3d now[3] = {nodePosition(idx[0]), nodePosition(idx[1]), nodePosition(idx[2])};
        const ChVector3d er1 = rest[1] - rest[0];
        const ChVector3d er2 = rest[2] - rest[0];
        const double restLen1 = er1.Length();
        const double restLen2 = er2.Length();
        if (restLen1 < 1.0e-9 || restLen2 < 1.0e-9)
            continue;
        const ChVector3d e1 = er1 / restLen1;
        ChVector3d n = er1.Cross(er2);
        if (n.Length() < 1.0e-12)
            continue;
        n.Normalize();
        const ChVector3d e2 = n.Cross(e1);
        const ChVector3d d1 = now[1] - now[0];
        const ChVector3d d2 = now[2] - now[0];
        const double f00 = d1.Dot(e1) / restLen1;
        const double f10 = d1.Dot(e2) / restLen1;
        const double a2 = er2.Dot(e1);
        const double b2 = er2.Dot(e2);
        if (std::abs(b2) < 1.0e-12)
            continue;
        const double f01 = (d2.Dot(e1) - f00 * a2) / b2;
        const double f11 = (d2.Dot(e2) - f10 * a2) / b2;
        double cornerT = 0.0;
        for (int index : idx)
            cornerT += strainTemperature(static_cast<std::size_t>(std::max(0, index)));
        const double eth = _material.thermalExpansion * (cornerT / 3.0 - _initialTemperature);
        const double exx = 0.5 * (f00 * f00 + f10 * f10 - 1.0) - eth;
        const double eyy = 0.5 * (f01 * f01 + f11 * f11 - 1.0) - eth;
        const double exy = 0.5 * (f00 * f01 + f10 * f11);
        const double mxx = cm * (exx + nu * eyy);
        const double myy = cm * (eyy + nu * exx);
        const double mxy = shear * exy;
        const double peak = vonMisesPlaneStress(mxx, myy, mxy);
        double scale = 1.0;
        if (t >= static_cast<std::size_t>(_skinTriangles)) {
            const std::size_t slot = t - static_cast<std::size_t>(_skinTriangles);
            if (slot < _triangleThermal.size() && _triangleThermal[slot])
                scale = _triangleThermal[slot]->stiffnessScale();
        }
        out[t] = std::isfinite(peak) ? static_cast<float>(peak * scale) : 0.0f;
    }
}

double Part::peakTension() const {
    double peak = 0.0;
    for (std::size_t q = 0; q < _quads.size(); ++q) {
        if (q < _rigidQuad.size() && _rigidQuad[q])
            continue;
        peak = std::max(peak, quadMechanics(q).tension);
    }
    for (const auto& element : _tetElements) {
        if (!element)
            continue;
        peak = std::max(peak, element->tension());
    }
    return peak;
}

void Part::discardTrial() {
    for (const auto& tet : _tetElements)
        if (tet)
            tet->discardTrial();
}

void Part::commitPlastic() {
    for (const auto& tet : _tetElements)
        if (tet)
            tet->commitPlastic();
}

void Part::refreshTetTrials() {
    if (_tetElements.empty() || _tets.empty())
        return;
    chrono::ChVectorDynamic<> fi(12);
    fi.setZero();
    for (std::size_t t = 0; t < _tetElements.size() && t < _tets.size(); ++t) {
        if (!_tetElements[t])
            continue;
        const auto& corners = _tets[t];
        chrono::ChVector3d p[4];
        bool ok = true;
        for (int k = 0; k < 4; ++k) {
            if (corners[static_cast<std::size_t>(k)] < 0) {
                ok = false;
                break;
            }
            p[k] = nodePosition(static_cast<std::size_t>(corners[static_cast<std::size_t>(k)]));
            if (!std::isfinite(p[k].x()) || !std::isfinite(p[k].y()) || !std::isfinite(p[k].z())) {
                ok = false;
                break;
            }
        }
        if (!ok)
            continue;
        const double vol = (p[1] - p[0]).Dot((p[2] - p[0]).Cross(p[3] - p[0])) / 6.0;
        if (vol <= 1.0e-12)
            continue;
        _tetElements[t]->ComputeInternalForces(fi);
    }
}

void Part::setImpactStep(double h) {
    for (const auto& tet : _tetElements)
        if (tet)
            tet->setImpactStep(h);
}

void Part::captureKinematics() {
    const std::size_t n = nodeCount();
    _savedPos.resize(n);
    _savedVel.resize(n);
    _savedAcc.resize(n);
    _savedRot.clear();
    _savedRotDt.clear();
    _savedAngAcc.clear();
    for (std::size_t i = 0; i < n; ++i) {
        _savedPos[i] = nodePosition(i);
        _savedVel[i] = nodeVelocity(i);
        if (!_xyzNodes.empty() && i < _xyzNodes.size() && _xyzNodes[i])
            _savedAcc[i] = _xyzNodes[i]->GetPosDt2();
        else if (i < _nodes.size() && _nodes[i])
            _savedAcc[i] = _nodes[i]->GetPosDt2();
        else
            _savedAcc[i] = {0, 0, 0};
    }
    if (_xyzNodes.empty()) {
        _savedRot.resize(n);
        _savedRotDt.resize(n);
        _savedAngAcc.resize(n);
        for (std::size_t i = 0; i < n && i < _nodes.size(); ++i) {
            if (!_nodes[i])
                continue;
            _savedRot[i] = _nodes[i]->GetRot();
            _savedRotDt[i] = _nodes[i]->GetRotDt();
            _savedAngAcc[i] = _nodes[i]->GetAngAccLocal();
        }
    }
}

void Part::restoreKinematics() {
    const std::size_t n = nodeCount();
    for (std::size_t i = 0; i < n && i < _savedPos.size(); ++i) {
        const chrono::ChVector3d vel = i < _savedVel.size() ? _savedVel[i] : nodeVelocity(i);
        setNodeKinematics(i, _savedPos[i], vel);
        const chrono::ChVector3d acc = i < _savedAcc.size() ? _savedAcc[i] : chrono::ChVector3d(0, 0, 0);
        if (!_xyzNodes.empty()) {
            if (i < _xyzNodes.size() && _xyzNodes[i] && !_xyzNodes[i]->IsFixed())
                _xyzNodes[i]->SetPosDt2(acc);
        } else if (i < _nodes.size() && _nodes[i] && !_nodes[i]->IsFixed()) {
            _nodes[i]->SetPosDt2(acc);
            if (i < _savedRot.size())
                _nodes[i]->SetRot(_savedRot[i]);
            if (i < _savedRotDt.size())
                _nodes[i]->SetRotDt(_savedRotDt[i]);
            if (i < _savedAngAcc.size())
                _nodes[i]->SetAngAccLocal(_savedAngAcc[i]);
        }
    }
}

void Part::bakeRest() {
    if (!_xyzNodes.empty()) {
        for (auto& node : _xyzNodes) {
            if (node)
                node->SetX0(node->GetPos());
        }
        for (const auto& tet : _tetElements) {
            if (tet)
                tet->refreshRest();
        }
        return;
    }
    for (auto& node : _nodes) {
        if (!node)
            continue;
        node->SetX0(node->Frame());
        node->ForceToRest();
    }
    if (_mesh) {
        for (const auto& element : _elements) {
            if (element)
                _mesh->InitializeElement(element);
        }
        for (const auto& element : _trianglesElements) {
            if (element)
                _mesh->InitializeElement(element);
        }
    }
}

chrono::ChVector3d Part::center() const {
    chrono::ChVector3d sum{0, 0, 0};
    const std::size_t count = nodeCount();
    if (count == 0)
        return sum;
    for (std::size_t i = 0; i < count; ++i)
        sum += nodePosition(i);
    return sum / static_cast<double>(count);
}

void Part::bounds(chrono::ChVector3d& lo, chrono::ChVector3d& hi, double pad) const {
    const double big = std::numeric_limits<double>::max();
    lo = {big, big, big};
    hi = {-big, -big, -big};
    for (std::size_t i = 0; i < nodeCount(); ++i) {
        const auto p = nodePosition(i);
        lo = {std::min(lo.x(), p.x()), std::min(lo.y(), p.y()), std::min(lo.z(), p.z())};
        hi = {std::max(hi.x(), p.x()), std::max(hi.y(), p.y()), std::max(hi.z(), p.z())};
    }
    const double margin = pad + 0.5 * _thickness;
    lo -= chrono::ChVector3d(margin, margin, margin);
    hi += chrono::ChVector3d(margin, margin, margin);
}

double Part::mass() const {
    double total = 0.0;
    for (double m : _nodeMass)
        total += m;
    return total;
}

bool Part::finite() const {
    for (std::size_t i = 0; i < nodeCount(); ++i) {
        const auto p = nodePosition(i);
        if (!std::isfinite(p.x()) || !std::isfinite(p.y()) || !std::isfinite(p.z()))
            return false;
        if (p.Length() > 1.0e3)
            return false;
        const auto v = nodeVelocity(i);
        if (!std::isfinite(v.x()) || !std::isfinite(v.y()) || !std::isfinite(v.z()))
            return false;
        if (v.Length() > 1.0e3)
            return false;
    }
    return true;
}

std::size_t Part::nodeCount() const {
    return _xyzNodes.empty() ? _nodes.size() : _xyzNodes.size();
}

chrono::ChVector3d Part::nodePosition(std::size_t index) const {
    if (!_xyzNodes.empty())
        return index < _xyzNodes.size() ? _xyzNodes[index]->GetPos() : chrono::ChVector3d(0, 0, 0);
    return index < _nodes.size() ? _nodes[index]->GetPos() : chrono::ChVector3d(0, 0, 0);
}

chrono::ChVector3d Part::nodeVelocity(std::size_t index) const {
    if (!_xyzNodes.empty())
        return index < _xyzNodes.size() ? _xyzNodes[index]->GetPosDt() : chrono::ChVector3d(0, 0, 0);
    return index < _nodes.size() ? _nodes[index]->GetPosDt() : chrono::ChVector3d(0, 0, 0);
}

chrono::ChVector3d Part::savedPosition(std::size_t index) const {
    if (index < _savedPos.size())
        return _savedPos[index];
    return nodePosition(index);
}

chrono::ChVector3d Part::savedVelocity(std::size_t index) const {
    if (index < _savedVel.size())
        return _savedVel[index];
    return nodeVelocity(index);
}

void Part::setNodeKinematics(std::size_t index,
                             const chrono::ChVector3d& pos,
                             const chrono::ChVector3d& vel) {
    if (!_xyzNodes.empty()) {
        if (index >= _xyzNodes.size() || !_xyzNodes[index] || _xyzNodes[index]->IsFixed())
            return;
        _xyzNodes[index]->SetPos(pos);
        _xyzNodes[index]->SetPosDt(vel);
        return;
    }
    if (index >= _nodes.size() || !_nodes[index] || _nodes[index]->IsFixed())
        return;
    _nodes[index]->SetPos(pos);
    _nodes[index]->SetPosDt(vel);
}

chrono::ChVector3d Part::nodeRest(std::size_t index) const {
    if (!_xyzNodes.empty())
        return index < _xyzNodes.size() ? _xyzNodes[index]->GetX0() : chrono::ChVector3d(0, 0, 0);
    return index < _nodes.size() ? _nodes[index]->GetX0().GetPos() : chrono::ChVector3d(0, 0, 0);
}

double Part::maxEdgeStretch(const int* index, int count) const {
    double stretch = 1.0;
    for (int i = 0; i < count; ++i) {
        const int a = index[i];
        const int b = index[(i + 1) % count];
        const double rest = (nodeRest(static_cast<std::size_t>(a)) - nodeRest(static_cast<std::size_t>(b))).Length();
        if (rest < 1.0e-9)
            continue;
        const double now =
            (nodePosition(static_cast<std::size_t>(a)) - nodePosition(static_cast<std::size_t>(b))).Length();
        // Thermal growth is free expansion, not a tear. Compare against the
        // length the edge would have at the strain temperature it has eased to.
        const double meanT = 0.5 * (strainTemperature(static_cast<std::size_t>(a)) +
                                    strainTemperature(static_cast<std::size_t>(b)));
        const double grown = rest * std::max(0.5, 1.0 + _material.thermalExpansion * (meanT - _initialTemperature));
        stretch = std::max(stretch, now / grown);
    }
    return stretch;
}

bool Part::nodeRigid(std::size_t index) const {
    return index < _rigidNode.size() && _rigidNode[index] != 0;
}

bool Part::nodeOnLiveTet(std::size_t index) const {
    if (_tetElements.empty())
        return true;
    const int node = static_cast<int>(index);
    for (std::size_t t = 0; t < _tets.size() && t < _tetElements.size(); ++t) {
        if (!_tetElements[t] || (t < _tetLive.size() && !_tetLive[t]))
            continue;
        for (int corner : _tets[t]) {
            if (corner == node)
                return true;
        }
    }
    return false;
}

void Part::clearLinkedNodes() {
    std::fill(_linked.begin(), _linked.end(), 0);
}

void Part::markLinkedNode(int index) {
    if (index < 0)
        return;
    const std::size_t count = std::max(_nodes.size(), _xyzNodes.size());
    if (_linked.size() < count)
        _linked.resize(count, 0);
    if (static_cast<std::size_t>(index) < _linked.size())
        _linked[static_cast<std::size_t>(index)] = 1;
}

bool Part::averagePlastic(std::size_t element, ReissnerPlasticData& out) const {
    if (element >= _elements.size() || !_elements[element])
        return false;
    ReissnerPlasticData sum;
    int count = 0;
    for (int ip = 0; ip < 4; ++ip) {
        const auto* data = dynamic_cast<const ReissnerPlasticData*>(_elements[element]->GetPlasticData(ip));
        if (!data)
            continue;
        sum.epsU += data->epsU;
        sum.epsV += data->epsV;
        sum.kurU += data->kurU;
        sum.kurV += data->kurV;
        sum.p_strain_acc += data->p_strain_acc;
        ++count;
    }
    if (count == 0)
        return false;
    const double scale = 1.0 / static_cast<double>(count);
    sum.epsU *= scale;
    sum.epsV *= scale;
    sum.kurU *= scale;
    sum.kurV *= scale;
    sum.p_strain_acc *= scale;
    out = sum;
    return sum.p_strain_acc > 0.0 || sum.epsU.Length2() + sum.kurU.Length2() > 0.0;
}

void Part::evaluatePlasticStrain(std::vector<float>& out) const {
    const std::size_t count = _quads.empty() ? _triangles.size() : _quads.size();
    out.assign(count, 0.0f);
    for (std::size_t i = 0; i < _elements.size() && i < out.size(); ++i) {
        ReissnerPlasticData plastic;
        if (averagePlastic(i, plastic))
            out[i] = static_cast<float>(plastic.p_strain_acc);
    }
    if (!_quads.empty())
        return;
    for (std::size_t t = 0; t < _triangles.size() && t < static_cast<std::size_t>(_skinTriangles); ++t) {
        for (std::size_t tet = 0; tet < _tets.size(); ++tet) {
            if (tet >= _tetElements.size() || !_tetElements[tet])
                continue;
            int hits = 0;
            for (int corner : _triangles[t])
                for (int node : _tets[tet])
                    hits += corner == node ? 1 : 0;
            if (hits == 3)
                out[t] = static_cast<float>(_tetElements[tet]->plasticStrain());
        }
    }
    if (!_tetElements.empty()) {
        const std::size_t base = out.size();
        out.resize(base + _tetElements.size(), 0.0f);
        for (std::size_t tet = 0; tet < _tetElements.size(); ++tet) {
            if (_tetElements[tet])
                out[base + tet] = static_cast<float>(_tetElements[tet]->plasticStrain());
        }
    }
    for (std::size_t s = 0; s < _triangleThermal.size(); ++s) {
        const std::size_t tri = static_cast<std::size_t>(_skinTriangles) + s;
        if (tri < out.size() && _triangleThermal[s])
            out[tri] = static_cast<float>(_triangleThermal[s]->plastic());
    }
}

double Part::tetPlasticStrain(std::size_t index) const {
    if (index >= _tetElements.size() || !_tetElements[index])
        return 0.0;
    return _tetElements[index]->plasticStrain();
}

int Part::cloneNode(int original) {
    const double share = std::max(1.0e-6, (original >= 0 && static_cast<std::size_t>(original) < _nodeMass.size())
                                              ? _nodeMass[static_cast<std::size_t>(original)] * 0.25
                                              : 1.0e-4);
    if (original >= 0 && static_cast<std::size_t>(original) < _nodeMass.size() &&
        _nodeMass[static_cast<std::size_t>(original)] > share)
        _nodeMass[static_cast<std::size_t>(original)] -= share;

    int clone = 0;
    if (!_xyzNodes.empty()) {
        const auto& src = _xyzNodes[static_cast<std::size_t>(original)];
        auto node = chrono_types::make_shared<chrono::fea::ChNodeFEAxyz>(src->GetPos());
        node->SetPosDt(src->GetPosDt());
        node->SetX0(src->GetX0());
        node->SetMass(0);
        node->SetFixed(_fixed);
        _mesh->AddNode(node);
        _xyzNodes.push_back(node);
        clone = static_cast<int>(_xyzNodes.size()) - 1;
    } else {
        const auto& src = _nodes[static_cast<std::size_t>(original)];
        auto node = chrono_types::make_shared<chrono::fea::ChNodeFEAxyzrot>(src->Frame());
        node->SetPosDt(src->GetPosDt());
        node->SetAngVelParent(src->GetAngVelParent());
        node->SetX0(src->GetX0());
        node->SetMass(0);
        node->GetInertia().fillDiagonal(0);
        node->SetFixed(_fixed);
        _mesh->AddNode(node);
        _nodes.push_back(node);
        clone = static_cast<int>(_nodes.size()) - 1;
    }
    _nodeMass.push_back(share);
    _loads.emplace_back(0, 0, 0);
    const double inherited = original >= 0 && static_cast<std::size_t>(original) < _nodeTemperature.size()
                                 ? _nodeTemperature[static_cast<std::size_t>(original)]
                                 : _initialTemperature;
    noteNewNode(clone, inherited);
    _tearTwins.push_back({clone, original});
    if (std::find(_boundaryNodes.begin(), _boundaryNodes.end(), original) != _boundaryNodes.end())
        _boundaryNodes.push_back(clone);
    return clone;
}

void Part::rebuildTrianglesFromQuads() {
    _triangles.clear();
    for (std::size_t q = 0; q < _quads.size(); ++q) {
        if (q < _rigidQuad.size() && _rigidQuad[q])
            continue;
        const auto& quad = _quads[q];
        _triangles.push_back({quad[0], quad[1], quad[2]});
        _triangles.push_back({quad[0], quad[2], quad[3]});
    }
}

int Part::insertMidpoint(int a, int b) {
    using namespace chrono;
    const std::size_t ia = static_cast<std::size_t>(a);
    const std::size_t ib = static_cast<std::size_t>(b);
    const double takeA = ia < _nodeMass.size() ? _nodeMass[ia] * 0.25 : 0.0;
    const double takeB = ib < _nodeMass.size() ? _nodeMass[ib] * 0.25 : 0.0;
    if (ia < _nodeMass.size())
        _nodeMass[ia] = std::max(0.0, _nodeMass[ia] - takeA);
    if (ib < _nodeMass.size())
        _nodeMass[ib] = std::max(0.0, _nodeMass[ib] - takeB);

    const auto& srcA = _nodes[ia];
    const auto& srcB = _nodes[ib];
    const ChVector3d mid = 0.5 * (srcA->GetPos() + srcB->GetPos());
    const ChVector3d midRest = 0.5 * (srcA->GetX0().GetPos() + srcB->GetX0().GetPos());
    auto node = chrono_types::make_shared<chrono::fea::ChNodeFEAxyzrot>(ChFrame<>(mid, srcA->GetRot()));
    node->SetX0(ChFrame<>(midRest, srcA->GetX0().GetRot()));
    node->SetPosDt(0.5 * (srcA->GetPosDt() + srcB->GetPosDt()));
    node->SetAngVelParent(0.5 * (srcA->GetAngVelParent() + srcB->GetAngVelParent()));
    node->SetMass(0);
    node->GetInertia().fillDiagonal(0);
    node->SetFixed(_fixed);
    _mesh->AddNode(node);
    _nodes.push_back(node);
    const int index = static_cast<int>(_nodes.size()) - 1;
    _nodeMass.push_back(std::max(1.0e-6, takeA + takeB));
    _loads.emplace_back(0, 0, 0);
    const double tempA = nodeTemperature(ia);
    const double tempB = nodeTemperature(ib);
    noteNewNode(index, 0.5 * (tempA + tempB));
    const bool boundaryA =
        std::find(_boundaryNodes.begin(), _boundaryNodes.end(), a) != _boundaryNodes.end();
    const bool boundaryB =
        std::find(_boundaryNodes.begin(), _boundaryNodes.end(), b) != _boundaryNodes.end();
    if (boundaryA && boundaryB)
        _boundaryNodes.push_back(index);
    return index;
}

void Part::detachQuad(std::size_t quad) {
    if (quad >= _quads.size())
        return;
    isolateQuad(static_cast<int>(quad));
    if (quad < _elements.size() && _elements[quad]) {
        _mesh->RemoveElement(_elements[quad]);
        _elements[quad].reset();
    }
}

bool Part::splitQuad(std::size_t quad, bool alongU) {
    using namespace chrono;
    using namespace chrono::fea;
    if (quad >= _quads.size() || quad >= _elements.size() || !_shellMaterial)
        return false;

    const int n0 = _quads[quad][0];
    const int n1 = _quads[quad][1];
    const int n2 = _quads[quad][2];
    const int n3 = _quads[quad][3];

    std::array<int, 4> kept;
    std::array<int, 4> opened;
    if (alongU) {
        const int left = insertMidpoint(n3, n0);
        const int right = insertMidpoint(n1, n2);
        kept = {n0, n1, right, left};
        opened = {cloneNode(left), cloneNode(right), n2, n3};
    } else {
        const int bottom = insertMidpoint(n0, n1);
        const int top = insertMidpoint(n2, n3);
        kept = {n0, bottom, top, n3};
        opened = {cloneNode(bottom), n1, n2, cloneNode(top)};
    }

    _quads[quad] = kept;

    auto element = chrono_types::make_shared<ChElementShellReissner4>();
    element->SetNodes(_nodes[opened[0]], _nodes[opened[1]], _nodes[opened[2]], _nodes[opened[3]]);
    element->AddLayer(_thickness, 0.0, makeQuadMaterial(_elements.size()));
    _mesh->AddElement(element);

    // InitializeElement snaps each node's rest pose to wherever it sits right
    // now. Park every node on its own rest pose first, so a shared corner keeps
    // the reference the neighbouring quads were built with, and the two new
    // pieces inherit that same unstressed shape.
    std::vector<int> involved = {kept[0], kept[1], kept[2], kept[3],
                                 opened[0], opened[1], opened[2], opened[3]};
    std::sort(involved.begin(), involved.end());
    involved.erase(std::unique(involved.begin(), involved.end()), involved.end());
    std::vector<ChFrame<>> parked;
    parked.reserve(involved.size());
    for (int index : involved) {
        auto& node = _nodes[static_cast<std::size_t>(index)];
        parked.push_back(node->Frame());
        node->SetPos(node->GetX0().GetPos());
        node->SetRot(node->GetX0().GetRot());
    }
    _elements[quad]->SetNodes(_nodes[kept[0]], _nodes[kept[1]], _nodes[kept[2]], _nodes[kept[3]]);
    _mesh->InitializeElement(_elements[quad]);
    _mesh->InitializeElement(element);
    for (std::size_t n = 0; n < involved.size(); ++n) {
        auto& node = _nodes[static_cast<std::size_t>(involved[n])];
        node->SetPos(parked[n].GetPos());
        node->SetRot(parked[n].GetRot());
    }

    _elements.push_back(element);
    _quads.push_back(opened);
    return true;
}

void Part::advanceShards(double dt, const chrono::ChVector3d& gravity) {
    using namespace chrono;
    if (dt <= 0.0)
        return;
    for (auto& shard : _shards) {
        if (shard.asleep || shard.mass <= 0.0) {
            shard.pendingForce = {0, 0, 0};
            shard.pendingTorque = {0, 0, 0};
            continue;
        }
        ChVector3d force = shard.mass * gravity + shard.pendingForce;
        ChVector3d torque = shard.pendingTorque;
        for (int k = 0; k < 4; ++k) {
            const std::size_t index = static_cast<std::size_t>(shard.nodes[static_cast<std::size_t>(k)]);
            if (index >= _loads.size())
                continue;
            const ChVector3d load = _loads[index];
            force += load;
            const ChVector3d arm = shard.rotation.Rotate(shard.local[static_cast<std::size_t>(k)]);
            torque += arm.Cross(load);
        }
        shard.pendingForce = {0, 0, 0};
        shard.pendingTorque = {0, 0, 0};
        shard.velocity += force * (dt / shard.mass);
        if (shard.velocity.Length() > 25.0)
            shard.velocity *= 25.0 / shard.velocity.Length();
        const ChVector3d torqueBody = shard.rotation.RotateBack(torque);
        ChVector3d spinBody = shard.rotation.RotateBack(shard.angularVelocity);
        spinBody.x() += torqueBody.x() / std::max(shard.inertiaU, 1.0e-9) * dt;
        spinBody.y() += torqueBody.y() / std::max(shard.inertiaV, 1.0e-9) * dt;
        spinBody.z() += torqueBody.z() / std::max(shard.inertiaN, 1.0e-9) * dt;
        shard.angularVelocity = shard.rotation.Rotate(spinBody);
        double radius = 0.0;
        for (int k = 0; k < 4; ++k)
            radius = std::max(radius, shard.local[static_cast<std::size_t>(k)].Length());
        const double spinCap = radius > 1.0e-4 ? 15.0 / radius : 40.0;
        const double spin = shard.angularVelocity.Length();
        if (spin > spinCap)
            shard.angularVelocity *= spinCap / spin;
        shard.com += shard.velocity * dt;
        const double turn = shard.angularVelocity.Length();
        if (turn > 1.0e-8) {
            ChQuaterniond step;
            step.SetFromAngleAxis(turn * dt, shard.angularVelocity / turn);
            shard.rotation = step * shard.rotation;
            shard.rotation.Normalize();
        }
    }
}

bool Part::discardAwakeShards() {
    std::vector<char> debris(nodeCount(), 0);
    bool any = false;
    std::vector<RigidShard> kept;
    kept.reserve(_shards.size());
    for (const RigidShard& shard : _shards) {
        if (shard.asleep || shard.mass <= 0.0) {
            kept.push_back(shard);
            continue;
        }
        any = true;
        for (int node : shard.nodes) {
            if (node >= 0 && static_cast<std::size_t>(node) < debris.size())
                debris[static_cast<std::size_t>(node)] = 1;
        }
    }
    if (!any)
        return false;
    _shards.swap(kept);

    std::vector<std::array<int, 3>> remain;
    remain.reserve(_triangles.size());
    int skin = 0;
    for (std::size_t i = 0; i < _triangles.size(); ++i) {
        const auto& tri = _triangles[i];
        bool onlyDebris = true;
        for (int corner : tri) {
            if (corner < 0 || static_cast<std::size_t>(corner) >= debris.size() ||
                !debris[static_cast<std::size_t>(corner)])
                onlyDebris = false;
        }
        if (onlyDebris)
            continue;
        if (static_cast<int>(i) < _skinTriangles)
            ++skin;
        remain.push_back(tri);
    }
    _triangles.swap(remain);
    _skinTriangles = skin;
    return true;
}

void Part::sanitizeShards() {
    auto write = [&](std::size_t index, const chrono::ChVector3d& pos, const chrono::ChVector3d& vel) {
        if (!_xyzNodes.empty()) {
            if (index >= _xyzNodes.size() || !_xyzNodes[index])
                return;
            _xyzNodes[index]->SetPos(pos);
            _xyzNodes[index]->SetPosDt(vel);
            return;
        }
        if (index >= _nodes.size() || !_nodes[index])
            return;
        _nodes[index]->SetPos(pos);
        _nodes[index]->SetPosDt(vel);
    };
    auto finiteVec = [](const chrono::ChVector3d& v) {
        return std::isfinite(v.x()) && std::isfinite(v.y()) && std::isfinite(v.z()) && v.Length() <= 1.0e3;
    };
    auto sanitizeNode = [&](std::size_t index) {
        chrono::ChVector3d p = nodePosition(index);
        chrono::ChVector3d v = nodeVelocity(index);
        if (finiteVec(p) && finiteVec(v))
            return;
        p = nodeRest(index);
        if (!finiteVec(p))
            p = {0, 0, 0};
        write(index, p, {0, 0, 0});
    };
    for (std::size_t i = 0; i < nodeCount(); ++i)
        sanitizeNode(i);
    for (auto& shard : _shards) {
        if (finiteVec(shard.com) && finiteVec(shard.velocity) && finiteVec(shard.angularVelocity))
            continue;
        chrono::ChVector3d com{0, 0, 0};
        int corners = 0;
        for (int k = 0; k < 4; ++k) {
            const int node = shard.nodes[static_cast<std::size_t>(k)];
            if (node < 0)
                continue;
            const chrono::ChVector3d p = savedPosition(static_cast<std::size_t>(node));
            if (finiteVec(p)) {
                com += p;
                ++corners;
            }
        }
        if (corners > 0)
            shard.com = com / static_cast<double>(corners);
        else if (!finiteVec(shard.com))
            shard.com = {0, 0, 0};
        shard.velocity = {0, 0, 0};
        shard.angularVelocity = {0, 0, 0};
        shard.asleep = true;
        for (int k = 0; k < 4; ++k) {
            const int node = shard.nodes[static_cast<std::size_t>(k)];
            if (node < 0)
                continue;
            const chrono::ChVector3d arm = shard.rotation.Rotate(shard.local[static_cast<std::size_t>(k)]);
            write(static_cast<std::size_t>(node), shard.com + arm, {0, 0, 0});
        }
    }
}

void Part::finishShards(bool floorEnabled, double floorHeight, const chrono::ChVector3d& arenaSize,
                       const chrono::ChVector3d& gravity, double restitution, double friction) {
    using namespace chrono;
    // A chip rests on the same surface as the body it came off.
    const double plane = floorHeight + contactOffset();
    const ChVector3d up{0, 1, 0};
    const double mu = std::max(0.0, friction);
    const double e = std::clamp(restitution, 0.0, 1.0);
    const bool noBodyForce = !floorEnabled && gravity.Length2() < 0.01;
    for (auto& shard : _shards) {
        if (!shard.asleep && noBodyForce && shard.mass > 0.0) {
            // No floor and no gravity: leftover Newton velocity would walk a
            // chip to y → −∞. The chip has already left the live solve.
            shard.velocity = {0, 0, 0};
            shard.angularVelocity = {0, 0, 0};
            shard.asleep = true;
        }
        if (!shard.asleep && floorEnabled && shard.mass > 0.0) {
            auto cornerAt = [&](int k) {
                return shard.com + shard.rotation.Rotate(shard.local[static_cast<std::size_t>(k)]);
            };
            auto applyImpulse = [&](const ChVector3d& impulse, const ChVector3d& arm) {
                shard.velocity += impulse / shard.mass;
                const ChVector3d moment = shard.rotation.RotateBack(arm.Cross(impulse));
                ChVector3d spin = shard.rotation.RotateBack(shard.angularVelocity);
                spin.x() += moment.x() / std::max(shard.inertiaU, 1.0e-9);
                spin.y() += moment.y() / std::max(shard.inertiaV, 1.0e-9);
                spin.z() += moment.z() / std::max(shard.inertiaN, 1.0e-9);
                shard.angularVelocity = shard.rotation.Rotate(spin);
            };
            double size = 0.0;
            for (int k = 0; k < 4; ++k)
                size = std::max(size, shard.local[static_cast<std::size_t>(k)].Length());
            const double skin = std::max({contactOffset(), contactSkin(), size * 0.05, 1.0e-4});
            // The floor pushes at the corner that is actually touching, so a
            // chip on its tip topples instead of being glued there.
            for (int pass = 0; pass < 4; ++pass) {
                for (int k = 0; k < 4; ++k) {
                    const ChVector3d point = cornerAt(k);
                    const double pen = plane - point.y();
                    if (pen <= 0.0)
                        continue;
                    const ChVector3d arm = point - shard.com;
                    const ChVector3d rxn = arm.Cross(up);
                    const ChVector3d rxnBody = shard.rotation.RotateBack(rxn);
                    const ChVector3d ang = shard.rotation.Rotate(
                        {rxnBody.x() / std::max(shard.inertiaU, 1.0e-9),
                         rxnBody.y() / std::max(shard.inertiaV, 1.0e-9),
                         rxnBody.z() / std::max(shard.inertiaN, 1.0e-9)});
                    const double eff = 1.0 / shard.mass + rxn.Dot(ang);
                    if (eff < 1.0e-12)
                        continue;
                    const double vn = (shard.velocity + shard.angularVelocity.Cross(arm)).Dot(up);
                    const double bias = std::min(pen / skin, 1.0);
                    const double j = vn < 0.0 ? -((1.0 + e) * vn + bias) / eff : -bias / eff;
                    if (j <= 0.0)
                        continue;
                    const double limited = std::min(j, shard.mass * std::max(1.0, std::abs(vn) + bias));
                    applyImpulse(up * limited, arm);
                    ChVector3d tangent = shard.velocity + shard.angularVelocity.Cross(arm);
                    tangent.y() = 0.0;
                    const double slide = tangent.Length();
                    if (slide > 1.0e-6 && mu > 0.0) {
                        const double jf = std::min(mu * limited, slide * shard.mass);
                        applyImpulse(tangent * (-jf / slide), arm);
                    }
                }
            }

            double lowest = std::numeric_limits<double>::max();
            for (int k = 0; k < 4; ++k)
                lowest = std::min(lowest, cornerAt(k).y());
            if (lowest < plane)
                shard.com.y() += plane - lowest;
            const double gAbs = std::abs(gravity.y());
            const double vSleep = std::sqrt(std::max(gAbs, 1.0e-6) * skin);
            const double wSleep = vSleep / skin;
            const bool slow = shard.velocity.Length() < 0.25 * vSleep &&
                              shard.angularVelocity.Length() < 0.25 * wSleep;
            if (slow && lowest <= plane + skin) {
                shard.velocity = {0, 0, 0};
                shard.angularVelocity = {0, 0, 0};
                shard.asleep = true;
            }

            const double hx = std::max(0.5, arenaSize.x() * 0.5);
            const double hz = std::max(0.5, arenaSize.z() * 0.5);
            const double ceiling = floorHeight + std::max(1.0, arenaSize.y());
            const double reach = 0.5 * _thickness;
            const ChVector3d wallOrigin[5] = {
                {-hx, 0.0, 0.0}, {hx, 0.0, 0.0}, {0.0, 0.0, -hz}, {0.0, 0.0, hz}, {0.0, ceiling, 0.0}};
            const ChVector3d wallNormal[5] = {
                {1.0, 0.0, 0.0}, {-1.0, 0.0, 0.0}, {0.0, 0.0, 1.0}, {0.0, 0.0, -1.0}, {0.0, -1.0, 0.0}};
            for (int pass = 0; pass < 2; ++pass) {
                for (int w = 0; w < 5; ++w) {
                    for (int k = 0; k < 4; ++k) {
                        const ChVector3d point = cornerAt(k);
                        const double pen = reach - (point - wallOrigin[w]).Dot(wallNormal[w]);
                        if (pen <= 1.0e-4)
                            continue;
                        const ChVector3d arm = point - shard.com;
                        const ChVector3d rxn = arm.Cross(wallNormal[w]);
                        const ChVector3d rxnBody = shard.rotation.RotateBack(rxn);
                        const ChVector3d ang = shard.rotation.Rotate(
                            {rxnBody.x() / std::max(shard.inertiaU, 1.0e-9),
                             rxnBody.y() / std::max(shard.inertiaV, 1.0e-9),
                             rxnBody.z() / std::max(shard.inertiaN, 1.0e-9)});
                        const double eff = 1.0 / shard.mass + rxn.Dot(ang);
                        if (eff < 1.0e-12)
                            continue;
                        const double vn = (shard.velocity + shard.angularVelocity.Cross(arm)).Dot(wallNormal[w]);
                        const double bias = std::min(std::max(0.0, pen - 1.0e-4) * 30.0, 1.0);
                        const double j = -(vn + bias) / eff;
                        if (j <= 0.0)
                            continue;
                        const double limited = std::min(j, shard.mass * 20.0);
                        applyImpulse(wallNormal[w] * limited, arm);
                        ChVector3d tangent = shard.velocity + shard.angularVelocity.Cross(arm);
                        tangent -= wallNormal[w] * tangent.Dot(wallNormal[w]);
                        const double slide = tangent.Length();
                        if (slide > 1.0e-6) {
                            const double jf = std::min(mu * limited, slide * shard.mass);
                            applyImpulse(tangent * (-jf / slide), arm);
                        }
                    }
                }
            }
        }
        for (int k = 0; k < 4; ++k) {
            const std::size_t index = static_cast<std::size_t>(shard.nodes[static_cast<std::size_t>(k)]);
            const ChVector3d arm = shard.rotation.Rotate(shard.local[static_cast<std::size_t>(k)]);
            const ChVector3d pos = shard.com + arm;
            const ChVector3d vel = shard.velocity + shard.angularVelocity.Cross(arm);
            if (!std::isfinite(pos.x()) || !std::isfinite(pos.y()) || !std::isfinite(pos.z()) ||
                !std::isfinite(vel.x()) || !std::isfinite(vel.y()) || !std::isfinite(vel.z()))
                continue;
            if (index < _xyzNodes.size() && _xyzNodes[index]) {
                _xyzNodes[index]->SetPos(pos);
                _xyzNodes[index]->SetPosDt(vel);
                _xyzNodes[index]->SetX0(pos);
                continue;
            }
            if (index >= _nodes.size())
                continue;
            auto& node = _nodes[index];
            node->SetPos(pos);
            node->SetRot(shard.rotation);
            node->SetPosDt(vel);
            node->SetAngVelParent(shard.angularVelocity);
            node->SetX0(node->Frame());
        }
    }
    sanitizeShards();
    // A chip that has actually come to rest is scenery. It is only marked
    // asleep after its speed dies, so a first floor touch still bounces.
    _shards.erase(std::remove_if(_shards.begin(), _shards.end(),
                                 [](const RigidShard& shard) { return shard.asleep || shard.mass <= 0.0; }),
                  _shards.end());
}

int Part::splitQuadGeometry(std::size_t quad, bool alongU) {
    if (quad >= _quads.size() || _nodes.empty())
        return -1;
    const int n0 = _quads[quad][0];
    const int n1 = _quads[quad][1];
    const int n2 = _quads[quad][2];
    const int n3 = _quads[quad][3];
    if (n0 < 0 || n1 < 0 || n2 < 0 || n3 < 0)
        return -1;
    const auto inRange = [&](int index) {
        return index >= 0 && static_cast<std::size_t>(index) < _nodes.size();
    };
    if (!inRange(n0) || !inRange(n1) || !inRange(n2) || !inRange(n3))
        return -1;

    std::array<int, 4> kept;
    std::array<int, 4> opened;
    if (alongU) {
        const int left = insertMidpoint(n3, n0);
        const int right = insertMidpoint(n1, n2);
        kept = {n0, n1, right, left};
        opened = {cloneNode(left), cloneNode(right), n2, n3};
    } else {
        const int bottom = insertMidpoint(n0, n1);
        const int top = insertMidpoint(n2, n3);
        kept = {n0, bottom, top, n3};
        opened = {cloneNode(bottom), n1, n2, cloneNode(top)};
    }

    _quads[quad] = kept;
    if (_elements.size() < _quads.size())
        _elements.resize(_quads.size());
    _quads.push_back(opened);
    _elements.push_back(nullptr);
    if (_rigidQuad.size() < _quads.size())
        _rigidQuad.resize(_quads.size(), 0);
    return static_cast<int>(_quads.size()) - 1;
}

void Part::isolateQuad(int quadIndex) {
    if (quadIndex < 0 || static_cast<std::size_t>(quadIndex) >= _quads.size())
        return;
    auto& quad = _quads[static_cast<std::size_t>(quadIndex)];
    for (int k = 0; k < 4; ++k) {
        const int node = quad[static_cast<std::size_t>(k)];
        bool shared = false;
        for (std::size_t q = 0; q < _quads.size() && !shared; ++q) {
            if (static_cast<int>(q) == quadIndex)
                continue;
            const auto& other = _quads[q];
            for (int j = 0; j < 4; ++j)
                shared = shared || other[static_cast<std::size_t>(j)] == node;
        }
        if (shared)
            quad[static_cast<std::size_t>(k)] = cloneNode(node);
    }
}

void Part::releaseShard(int quadIndex) {
    using namespace chrono;
    if (quadIndex < 0 || static_cast<std::size_t>(quadIndex) >= _quads.size())
        return;
    if (_rigidQuad.size() < _quads.size())
        _rigidQuad.resize(_quads.size(), 0);
    if (_rigidNode.size() < _nodes.size())
        _rigidNode.resize(_nodes.size(), 0);
    if (_rigidQuad[static_cast<std::size_t>(quadIndex)])
        return;

    const auto& quad = _quads[static_cast<std::size_t>(quadIndex)];
    RigidShard shard;
    shard.rotation = _nodes[quad[0]]->GetRot();
    for (int k = 0; k < 4; ++k) {
        const int index = quad[static_cast<std::size_t>(k)];
        if (index < 0 || static_cast<std::size_t>(index) >= _nodes.size())
            return;
        const auto& node = _nodes[static_cast<std::size_t>(index)];
        shard.nodes[static_cast<std::size_t>(k)] = index;
        shard.com += node->GetPos();
        shard.velocity += node->GetPosDt();
        shard.angularVelocity += node->GetAngVelParent();
        if (static_cast<std::size_t>(index) < _nodeMass.size())
            shard.mass += _nodeMass[static_cast<std::size_t>(index)];
    }
    shard.com *= 0.25;
    shard.velocity *= 0.25;
    shard.angularVelocity *= 0.25;
    if (shard.velocity.Length() > 25.0)
        shard.velocity *= 25.0 / shard.velocity.Length();
    if (shard.mass <= 1.0e-8)
        shard.mass = std::max(1.0e-4, _material.density * _thickness * 1.0e-4);
    for (int k = 0; k < 4; ++k) {
        const auto& node = _nodes[quad[static_cast<std::size_t>(k)]];
        shard.local[static_cast<std::size_t>(k)] = shard.rotation.RotateBack(node->GetPos() - shard.com);
    }
    double meanT = 0.0;
    for (int index : shard.nodes)
        meanT += nodeTemperature(static_cast<std::size_t>(std::max(0, index)));
    meanT *= 0.25;
    const double factor =
        std::max(0.5, 1.0 + _material.thermalExpansion * (meanT - _initialTemperature));
    for (int k = 0; k < 4; ++k)
        shard.localRest[static_cast<std::size_t>(k)] = shard.local[static_cast<std::size_t>(k)] / factor;
    shard.thermalFactor = factor;
    const double au = (shard.local[1] - shard.local[0]).Length();
    const double av = (shard.local[3] - shard.local[0]).Length();
    const double radius = 0.5 * std::sqrt(au * au + av * av);
    const double spinCap = radius > 1.0e-4 ? 15.0 / radius : 40.0;
    const double spin = shard.angularVelocity.Length();
    if (spin > spinCap)
        shard.angularVelocity *= spinCap / spin;
    const double h = _thickness;
    const double minArm = std::max(h, 1.0e-2);
    const double floorI = shard.mass * minArm * minArm / 12.0;
    shard.inertiaU = std::max(shard.mass * (av * av + h * h) / 12.0, floorI);
    shard.inertiaV = std::max(shard.mass * (au * au + h * h) / 12.0, floorI);
    shard.inertiaN = std::max(shard.mass * (au * au + av * av) / 12.0, floorI);

    if (_rigidNode.size() < _nodes.size())
        _rigidNode.resize(_nodes.size(), 0);
    for (int k = 0; k < 4; ++k) {
        const int index = quad[static_cast<std::size_t>(k)];
        bool shared = false;
        for (std::size_t q = 0; q < _quads.size() && !shared; ++q) {
            if (static_cast<int>(q) == quadIndex)
                continue;
            if (q < _rigidQuad.size() && _rigidQuad[q])
                continue;
            if (q >= _elements.size() || !_elements[q])
                continue;
            for (int corner : _quads[q])
                shared = shared || corner == index;
        }
        if (shared)
            continue;
        auto& node = _nodes[static_cast<std::size_t>(index)];
        node->SetX0(node->Frame());
        node->SetFixed(true);
        _rigidNode[static_cast<std::size_t>(index)] = 1;
    }
    if (static_cast<std::size_t>(quadIndex) < _elements.size() && _elements[static_cast<std::size_t>(quadIndex)]) {
        _mesh->RemoveElement(_elements[static_cast<std::size_t>(quadIndex)]);
        _elements[static_cast<std::size_t>(quadIndex)].reset();
    }
    _rigidQuad[static_cast<std::size_t>(quadIndex)] = 1;
    _shards.push_back(shard);
}

void Part::promoteLooseFragments(bool moltenOnly) {
    if (_quads.size() < 2 || _nodes.empty())
        return;

    const int count = static_cast<int>(_quads.size());
    std::vector<int> parent(static_cast<std::size_t>(count));
    for (int i = 0; i < count; ++i)
        parent[static_cast<std::size_t>(i)] = i;
    auto findRoot = [&](int i) {
        while (parent[static_cast<std::size_t>(i)] != i) {
            parent[static_cast<std::size_t>(i)] = parent[static_cast<std::size_t>(parent[static_cast<std::size_t>(i)])];
            i = parent[static_cast<std::size_t>(i)];
        }
        return i;
    };
    auto unite = [&](int a, int b) {
        a = findRoot(a);
        b = findRoot(b);
        if (a != b)
            parent[static_cast<std::size_t>(b)] = a;
    };

    std::map<int, int> owner;
    for (int q = 0; q < count; ++q) {
        for (int k = 0; k < 4; ++k) {
            const int node = _quads[static_cast<std::size_t>(q)][static_cast<std::size_t>(k)];
            const auto found = owner.find(node);
            if (found == owner.end())
                owner[node] = q;
            else
                unite(q, found->second);
        }
    }

    std::map<int, std::vector<int>> groups;
    for (int q = 0; q < count; ++q)
        groups[findRoot(q)].push_back(q);

    if (_rigidQuad.size() < _quads.size())
        _rigidQuad.resize(_quads.size(), 0);
    if (_rigidNode.size() < _nodes.size())
        _rigidNode.resize(_nodes.size(), 0);

    // A free fragment of any size leaves the shell solve. Cutting stops at
    // 64 shards so a face cannot become a thousand-body contact problem;
    // pieces already free are still released, just not subdivided further.
    constexpr int kMaxShards = 128;
    const double minEdge = _material.minimumShardSize > 1.0e-4 ? _material.minimumShardSize : 0.02;

    struct LiveGroup {
        std::vector<int> quads;
        bool anchored = false;
    };
    std::vector<LiveGroup> liveGroups;
    liveGroups.reserve(groups.size());
    for (const auto& group : groups) {
        LiveGroup live;
        for (int quadIndex : group.second) {
            if (quadIndex < 0 || _rigidQuad[static_cast<std::size_t>(quadIndex)])
                continue;
            live.quads.push_back(quadIndex);
            const auto& quad = _quads[static_cast<std::size_t>(quadIndex)];
            for (int k = 0; k < 4; ++k) {
                const int index = quad[static_cast<std::size_t>(k)];
                if (index < 0 || static_cast<std::size_t>(index) >= _nodes.size()) {
                    live.anchored = true;
                    continue;
                }
                const bool linked =
                    static_cast<std::size_t>(index) < _linked.size() && _linked[static_cast<std::size_t>(index)];
                live.anchored = live.anchored || _nodes[static_cast<std::size_t>(index)]->IsFixed() || linked;
            }
        }
        if (live.quads.empty())
            continue;
        liveGroups.push_back(std::move(live));
    }

    for (auto& live : liveGroups) {
        std::vector<int>& quads = live.quads;
        (void)live.anchored;

        auto quadMolten = [&](int quadIndex) {
            double kelvin = 0.0;
            int samples = 0;
            for (int corner : _quads[static_cast<std::size_t>(quadIndex)]) {
                if (corner < 0)
                    continue;
                kelvin += nodeTemperature(static_cast<std::size_t>(corner));
                ++samples;
            }
            return samples > 0 && kelvin / samples >= _material.meltTemperature;
        };

        std::vector<int> failed;
        failed.reserve(quads.size());
        for (int quadIndex : quads) {
            if (quadIndex < 0 || static_cast<std::size_t>(quadIndex) >= _quads.size())
                continue;
            const bool rigid = static_cast<std::size_t>(quadIndex) < _rigidQuad.size() &&
                               _rigidQuad[static_cast<std::size_t>(quadIndex)];
            const bool liveElem = static_cast<std::size_t>(quadIndex) < _elements.size() &&
                                  _elements[static_cast<std::size_t>(quadIndex)];
            const bool molten = quadMolten(quadIndex);
            if (moltenOnly) {
                if (molten)
                    failed.push_back(quadIndex);
                continue;
            }
            if (molten)
                failed.push_back(quadIndex);
            else if ((rigid || !liveElem) && _material.fracture == FractureMode::Brittle)
                failed.push_back(quadIndex);
        }
        if (failed.empty())
            continue;
        quads.swap(failed);

        bool splitting = !moltenOnly;
        while (splitting) {
            splitting = false;
            if (static_cast<int>(_shards.size()) + static_cast<int>(quads.size()) >= kMaxShards)
                break;
            for (std::size_t i = 0; i < quads.size(); ++i) {
                if (static_cast<int>(_shards.size()) + static_cast<int>(quads.size()) >= kMaxShards)
                    break;
                const auto corners = _quads[static_cast<std::size_t>(quads[i])];
                bool usable = true;
                for (int k = 0; k < 4; ++k)
                    usable = usable && corners[static_cast<std::size_t>(k)] >= 0 &&
                             static_cast<std::size_t>(corners[static_cast<std::size_t>(k)]) < _nodes.size();
                if (!usable)
                    continue;
                const double edgeU =
                    (_nodes[corners[1]]->GetPos() - _nodes[corners[0]]->GetPos()).Length();
                const double edgeV =
                    (_nodes[corners[3]]->GetPos() - _nodes[corners[0]]->GetPos()).Length();
                const double longest = std::max(edgeU, edgeV);
                const double shortest = std::max(1.0e-9, std::min(edgeU, edgeV));
                const double aspect = longest / shortest;
                if (aspect <= 2.0)
                    continue;
                if (longest * 0.5 < minEdge && aspect <= 3.0)
                    continue;
                if (longest * 0.5 < 1.0e-3)
                    continue;
                // alongU halves the V edge, so it shortens a long V side.
                const int opened = splitQuadGeometry(static_cast<std::size_t>(quads[i]), edgeV >= edgeU);
                if (opened < 0)
                    continue;
                quads.push_back(opened);
                splitting = true;
                break;
            }
        }
        for (int quadIndex : quads) {
            isolateQuad(quadIndex);
            releaseShard(quadIndex);
        }
    }
}

bool Part::shatterQuads() {
    using chrono::ChVector3d;
    if (_fixed || _elements.empty() || _material.fractureStrength <= 0.0)
        return false;

    const double minEdge = _material.minimumShardSize > 1.0e-4 ? _material.minimumShardSize : 0.02;

    struct Candidate {
        std::size_t quad = 0;
        double tension = 0.0;
        double theta = 0.0;
        double edgeU = 0.0;
        double edgeV = 0.0;
        bool over = false;
    };
    std::vector<Candidate> candidates;
    const std::size_t limit = std::min(_quads.size(), _elements.size());
    for (std::size_t q = 0; q < limit; ++q) {
        if (q < _rigidQuad.size() && _rigidQuad[q])
            continue;
        if (!_elements[q])
            continue;
        const QuadMechanics fiber = quadMechanics(q);
        if (!fiber.ok)
            continue;
        double kelvin = 0.0;
        int corners = 0;
        for (int index : _quads[q]) {
            if (index < 0)
                continue;
            kelvin += strainTemperature(static_cast<std::size_t>(index));
            ++corners;
        }
        kelvin = corners ? kelvin / corners : _initialTemperature;
        const double strength = _material.fractureStrength * thermalStrengthRetention(_material, kelvin);
        const double longest = std::max(fiber.edgeU, fiber.edgeV);
        // Handbook Rankine already includes typical flaw size. Requiring a
        // Griffith energy as well left 12×-overstrength panes sitting in the
        // floor: G_c·A is not a second gate on σ_t.
        const double crackArea = std::max(fiber.edgeU, fiber.edgeV) * std::max(_thickness, 1.0e-6);
        const bool over = fiber.tension >= strength && energyAllowsOpening(_material, fiber.energy, crackArea);
        const bool dust = _material.shatterToDust && _shattering && longest > minEdge;
        if (!over && !dust)
            continue;
        if (over)
            _peakFractureTension = std::max(_peakFractureTension, fiber.tension);
        candidates.push_back({q, fiber.tension, fiber.theta, fiber.edgeU, fiber.edgeV, over});
    }
    std::sort(candidates.begin(), candidates.end(),
              [](const Candidate& a, const Candidate& b) { return a.tension > b.tension; });

    // Rankine-failed elements leave the live shell this step. A 4-quad peel
    // plus an 8% contact remnant is what left a 12×-overstrength pane sitting
    // in the floor and pinning whatever was welded to it.
    const int budget = std::min(128, std::max(1, static_cast<int>(candidates.size())));
    bool changed = false;
    int made = 0;
    for (const Candidate& candidate : candidates) {
        if (made >= budget)
            break;
        if (candidate.quad >= _quads.size())
            continue;
        if (candidate.over) {
            detachQuad(candidate.quad);
        } else {
            const double longest = std::max(candidate.edgeU, candidate.edgeV);
            const double shortest = std::min(candidate.edgeU, candidate.edgeV);
            const double aspect = shortest > 1.0e-9 ? longest / shortest : 1.0;
            if (longest < minEdge || aspect > 8.0)
                detachQuad(candidate.quad);
            else {
                bool alongU = candidate.edgeV >= candidate.edgeU;
                if ((_material.breakFabric == BreakFabric::Grain || !_material.plies.empty()) &&
                    candidate.quad < _quads.size()) {
                    const auto& quad = _quads[candidate.quad];
                    bool finiteQuad = true;
                    for (int node : quad)
                        finiteQuad = finiteQuad && node >= 0 && static_cast<std::size_t>(node) < _nodes.size();
                    if (finiteQuad) {
                        const ChVector3d u = _nodes[static_cast<std::size_t>(quad[1])]->GetPos() -
                                             _nodes[static_cast<std::size_t>(quad[0])]->GetPos();
                        const ChVector3d v = _nodes[static_cast<std::size_t>(quad[3])]->GetPos() -
                                             _nodes[static_cast<std::size_t>(quad[0])]->GetPos();
                        ChVector3d fiber(_material.grainDirection[0], _material.grainDirection[1],
                                         _material.grainDirection[2]);
                        if (!_material.plies.empty())
                            fiber = {_material.plies.front().direction[0], _material.plies.front().direction[1],
                                     _material.plies.front().direction[2]};
                        if (fiber.Length2() > 0.25 && u.Length2() > 1.0e-12 && v.Length2() > 1.0e-12) {
                            fiber.Normalize();
                            alongU = std::abs(u.GetNormalized().Dot(fiber)) >= std::abs(v.GetNormalized().Dot(fiber));
                        }
                    }
                }
                if (!splitQuad(candidate.quad, alongU))
                    continue;
            }
        }
        changed = true;
        _shattering = true;
        ++made;
    }
    if (changed) {
        promoteLooseFragments();
        freezeUnusedNodes();
        rebuildTrianglesFromQuads();
    }
    return changed;
}

bool Part::tearFailed() {
    using chrono::ChVector3d;
    for (const auto& tet : _tetElements)
        if (tet)
            tet->commitPlastic();
    bool changed = false;
    bool tetChanged = false;
    if (!_tetElements.empty())
        tetChanged = releaseMoltenTets();
    if (!_elements.empty()) {
        if (_material.fracture == FractureMode::Brittle) {
            const bool shattered = shatterQuads();
            if (shattered)
                _torn = true;
            return shattered || tetChanged;
        }

        std::vector<std::size_t> failed;
        for (std::size_t q = 0; q < _quads.size() && q < _elements.size(); ++q) {
            const double stretch = maxEdgeStretch(_quads[q].data(), 4);
            ReissnerPlasticData plastic;
            const double accumulated = averagePlastic(q, plastic) ? plastic.p_strain_acc : 0.0;
            double kelvin = 0.0;
            int corners = 0;
            for (int index : _quads[q]) {
                if (index < 0)
                    continue;
                kelvin += nodeTemperature(static_cast<std::size_t>(index));
                ++corners;
            }
            kelvin = corners ? kelvin / corners : _initialTemperature;
            const double ductility = thermalDuctilityRetention(_material, kelvin);
            const double stretchLimit = liveStretchLimit(_material, ductility);
            const double plasticLimit = _material.failPlasticStrain * ductility;
            const bool stretchFail =
                _material.fracture != FractureMode::Brittle && stretch >= stretchLimit;
            const bool plasticFail = _material.fracture == FractureMode::Ductile &&
                                     plasticLimit > 0.0 && accumulated >= plasticLimit;
            if (stretchFail || plasticFail)
                failed.push_back(q);
        }
        for (std::size_t q : failed) {
            double kelvin = 0.0;
            int corners = 0;
            for (int index : _quads[q]) {
                if (index < 0)
                    continue;
                kelvin += nodeTemperature(static_cast<std::size_t>(index));
                ++corners;
            }
            kelvin = corners ? kelvin / corners : _initialTemperature;
            detachQuad(q);
            // A torn patch that can no longer carry stress leaves the implicit
            // solve as a rigid chip. Leaving it in as a molten shell is what
            // aborts the next step.
            if (kelvin >= _material.meltTemperature)
                releaseShard(static_cast<int>(q));
            changed = true;
        }
        if (changed) {
            // Room-temp fragments that the tear cut loose must leave too:
            // a disconnected Reissner shell under a 10x yield load NaNs (X1).
            // Molten-only promotion would keep those shells in HHT.
            promoteLooseFragments(false);
            freezeUnusedNodes();
            rebuildTrianglesFromQuads();
        }
    } else if (!_trianglesElements.empty()) {
        std::vector<float> stress;
        if (_material.fracture == FractureMode::Brittle)
            evaluateStress(stress);
        std::vector<std::size_t> failed;
        for (std::size_t t = static_cast<std::size_t>(_skinTriangles);
             t < _triangles.size() && (t - static_cast<std::size_t>(_skinTriangles)) < _trianglesElements.size(); ++t) {
            const double stretch = maxEdgeStretch(_triangles[t].data(), 3);
            const bool stretchFail = _material.fracture != FractureMode::Brittle && stretch >= _material.tearStretch;
            const bool brittleFail = _material.fracture == FractureMode::Brittle &&
                                     t < stress.size() &&
                                     stress[t] >= static_cast<float>(_material.fractureStrength);
            if (stretchFail || brittleFail) {
                const ChVector3d a = nodePosition(static_cast<std::size_t>(_triangles[t][0]));
                const ChVector3d b = nodePosition(static_cast<std::size_t>(_triangles[t][1]));
                const ChVector3d c = nodePosition(static_cast<std::size_t>(_triangles[t][2]));
                const double area = 0.5 * (b - a).Cross(c - a).Length();
                const double longest = std::max({(b - a).Length(), (c - b).Length(), (a - c).Length()});
                const double sigma = t < stress.size() ? static_cast<double>(stress[t]) : 0.0;
                const double energy = 0.5 * sigma * sigma / std::max(_material.youngsModulus, 1.0) * area *
                                      std::max(_thickness, 1.0e-6);
                if (energyAllowsOpening(_material, energy, longest * std::max(_thickness, 1.0e-6)))
                    failed.push_back(t);
            }
        }
        for (std::size_t t : failed) {
            const std::size_t slot = t - static_cast<std::size_t>(_skinTriangles);
            int cloned[3];
            for (int i = 0; i < 3; ++i)
                cloned[i] = cloneNode(_triangles[t][static_cast<std::size_t>(i)]);
            _trianglesElements[slot]->SetNodes(_xyzNodes[cloned[0]], _xyzNodes[cloned[1]], _xyzNodes[cloned[2]],
                                            nullptr, nullptr, nullptr);
            _mesh->InitializeElement(_trianglesElements[slot]);
            _triangles[t] = {cloned[0], cloned[1], cloned[2]};
            changed = true;
        }
    }
    if (changed)
        _torn = true;
    return tetChanged || changed;
}

int Part::liveTetCount() const {
    int live = 0;
    for (char flag : _tetLive)
        live += flag ? 1 : 0;
    return live;
}

int Part::invertedTetCount() const {
    int inverted = 0;
    for (std::size_t t = 0; t < _tets.size(); ++t) {
        if (t >= _tetLive.size() || !_tetLive[t])
            continue;
        const auto& corners = _tets[t];
        chrono::ChVector3d rest[4];
        chrono::ChVector3d live[4];
        bool ok = true;
        for (int k = 0; k < 4; ++k) {
            const int idx = corners[static_cast<std::size_t>(k)];
            if (idx < 0) {
                ok = false;
                break;
            }
            rest[k] = restPosition(static_cast<std::size_t>(idx));
            live[k] = nodePosition(static_cast<std::size_t>(idx));
            if (!std::isfinite(live[k].x()) || !std::isfinite(live[k].y()) || !std::isfinite(live[k].z())) {
                ok = false;
                break;
            }
        }
        if (!ok)
            continue;
        const double restVol = (rest[1] - rest[0]).Dot((rest[2] - rest[0]).Cross(rest[3] - rest[0]));
        const double liveVol = (live[1] - live[0]).Dot((live[2] - live[0]).Cross(live[3] - live[0]));
        if (!std::isfinite(liveVol) || liveVol * restVol <= 0.0)
            ++inverted;
        else if (std::abs(restVol) > 0.0 && std::abs(liveVol) < 1.0e-12 * std::abs(restVol))
            ++inverted;
    }
    return inverted;
}

bool Part::pastBrittleRankine() const {
    if (_material.fracture != FractureMode::Brittle || _material.fractureStrength <= 0.0)
        return false;
    const double strength = _material.fractureStrength * std::max(thermalStrengthRetention(_material, _temperature), 1.0e-4);
    for (std::size_t t = 0; t < _tetElements.size(); ++t) {
        if (!_tetElements[t] || (t < _tetLive.size() && !_tetLive[t]))
            continue;
        if (_tetElements[t]->overstrength())
            return true;
        // A collapsed or inverted tet never reports Rankine; restoring it
        // hid every corner hit. On glass that element has failed.
        if (_tetElements[t]->inverted() || _tetElements[t]->deformationSingular())
            return true;
        if (_tetElements[t]->rawTension() >= strength)
            return true;
        if (brittleShouldCrush(_material, _temperature, _tetElements[t]->rawCompression()))
            return true;
    }
    const std::size_t limit = std::min(_quads.size(), _elements.size());
    for (std::size_t q = 0; q < limit; ++q) {
        if (q < _rigidQuad.size() && _rigidQuad[q])
            continue;
        if (!_elements[q])
            continue;
        if (quadMechanics(q).tension >= strength)
            return true;
    }
    return false;
}

int Part::liveElementCount() const {
    if (!_tetElements.empty())
        return liveTetCount();
    int live = 0;
    for (const auto& element : _elements)
        if (element)
            ++live;
    for (const auto& element : _trianglesElements)
        if (element)
            ++live;
    return live;
}

int Part::totalElementCount() const {
    if (!_tetElements.empty())
        return tetCount();
    return static_cast<int>(_elements.size() + _trianglesElements.size());
}

bool Part::hasFreeNodes() const {
    for (const auto& node : _xyzNodes)
        if (node && !node->IsFixed())
            return true;
    for (const auto& node : _nodes)
        if (node && !node->IsFixed())
            return true;
    return false;
}

bool Part::tetLive(std::size_t index) const {
    return index < _tetLive.size() && _tetLive[index];
}

void Part::buildVolume(const PartDefinition& definition) {
    using namespace chrono;
    using namespace chrono::fea;

    const VolumeMesh& volume = definition.volume;
    _nodes.clear();
    _elements.clear();
    _quads.clear();
    _triangles.clear();
    _boundaryNodes.clear();
    _xyzNodes.clear();
    _nodeMass.assign(volume.vertices.size(), 0.0);
    _tets = volume.tets;
    _tetElements.clear();
    _tetLive.assign(_tets.size(), 1);
    _skinTriangles = 0;

    for (const ChVector3d& vertex : volume.vertices) {
        auto node = chrono_types::make_shared<ChNodeFEAxyz>(vertex);
        node->SetMass(0);
        node->SetFixed(_fixed);
        node->SetX0(vertex);
        _mesh->AddNode(node);
        _xyzNodes.push_back(node);
    }

    double minX = 1.0e30;
    double minY = 1.0e30;
    double minZ = 1.0e30;
    double maxX = -1.0e30;
    double maxY = -1.0e30;
    double maxZ = -1.0e30;
    for (const ChVector3d& vertex : volume.vertices) {
        minX = std::min(minX, vertex.x());
        minY = std::min(minY, vertex.y());
        minZ = std::min(minZ, vertex.z());
        maxX = std::max(maxX, vertex.x());
        maxY = std::max(maxY, vertex.y());
        maxZ = std::max(maxZ, vertex.z());
    }

    for (std::size_t t = 0; t < _tets.size(); ++t) {
        const auto& tet = _tets[t];
        bool valid = true;
        for (int index : tet)
            valid = valid && index >= 0 && static_cast<std::size_t>(index) < _xyzNodes.size();
        if (!valid)
            continue;
        auto element = chrono_types::make_shared<ThermalTet4>();
        element->configure(_material.youngsModulus, _material.poissonRatio,
                           _material.allowPlasticity ? _material.yieldStress : 0.0,
                           _material.density, _material.dampingBeta);
        element->setHardening(_material.ultimateStress);
        if (_material.fracture == FractureMode::Brittle)
            element->setRankine(rankineTensionLimit(_material), rankineCompressionLimit(_material));
        if (_material.fracture == FractureMode::Ductile)
            element->setFailureStrain(_material.failPlasticStrain);
        element->SetNodes(_xyzNodes[tet[0]], _xyzNodes[tet[1]], _xyzNodes[tet[2]], _xyzNodes[tet[3]]);
        _mesh->AddElement(element);
        _tetElements.push_back(element);
        const ChVector3d p0 = volume.vertices[static_cast<std::size_t>(tet[0])];
        const ChVector3d p1 = volume.vertices[static_cast<std::size_t>(tet[1])];
        const ChVector3d p2 = volume.vertices[static_cast<std::size_t>(tet[2])];
        const ChVector3d p3 = volume.vertices[static_cast<std::size_t>(tet[3])];
        const double tetVolume = std::abs((p1 - p0).Dot((p2 - p0).Cross(p3 - p0))) / 6.0;
        const double share = _material.density * tetVolume / 4.0;
        for (int index : tet)
            _nodeMass[static_cast<std::size_t>(index)] += share;
    }

    std::vector<char> boundary(volume.vertices.size(), 0);
    for (const auto& face : volume.boundary) {
        if (face[0] < 0 || face[1] < 0 || face[2] < 0)
            continue;
        if (static_cast<std::size_t>(face[0]) >= _xyzNodes.size() ||
            static_cast<std::size_t>(face[1]) >= _xyzNodes.size() ||
            static_cast<std::size_t>(face[2]) >= _xyzNodes.size())
            continue;
        _triangles.push_back(face);
        for (int index : face)
            boundary[static_cast<std::size_t>(index)] = 1;
    }
    _skinTriangles = static_cast<int>(_triangles.size());
    for (std::size_t i = 0; i < boundary.size(); ++i) {
        if (boundary[i])
            _boundaryNodes.push_back(static_cast<int>(i));
    }
    if (definition.clampBoundary) {
        const double extX = maxX - minX;
        const double extY = maxY - minY;
        const double extZ = maxZ - minZ;
        int thick = 0;
        double thickLen = extX;
        if (extY < thickLen) {
            thick = 1;
            thickLen = extY;
        }
        if (extZ < thickLen)
            thick = 2;
        const int axisA = (thick + 1) % 3;
        const int axisB = (thick + 2) % 3;
        const double inA = axisA == 0 ? extX : (axisA == 1 ? extY : extZ);
        const double inB = axisB == 0 ? extX : (axisB == 1 ? extY : extZ);
        const double aspect = std::max(inA, inB) / std::max(std::min(inA, inB), 1.0e-9);
        // Square plates clamp every in-plane silhouette node, matching the
        // Reissner outline. Long bars keep the min-X cantilever convention.
        const bool plate = aspect < 2.0;
        auto coord = [](const ChVector3d& p, int axis) {
            return axis == 0 ? p.x() : (axis == 1 ? p.y() : p.z());
        };
        auto lo = [&](int axis) { return axis == 0 ? minX : (axis == 1 ? minY : minZ); };
        auto hi = [&](int axis) { return axis == 0 ? maxX : (axis == 1 ? maxY : maxZ); };
        for (int index : _boundaryNodes) {
            if (index < 0 || static_cast<std::size_t>(index) >= _xyzNodes.size())
                continue;
            const ChVector3d p = _xyzNodes[static_cast<std::size_t>(index)]->GetPos();
            bool clamp = false;
            if (plate) {
                if (coord(p, axisA) <= lo(axisA) + 1.0e-4 || coord(p, axisA) >= hi(axisA) - 1.0e-4)
                    clamp = true;
                if (coord(p, axisB) <= lo(axisB) + 1.0e-4 || coord(p, axisB) >= hi(axisB) - 1.0e-4)
                    clamp = true;
            } else if (p.x() <= minX + 1.0e-4) {
                clamp = true;
            }
            if (clamp)
                _xyzNodes[static_cast<std::size_t>(index)]->SetFixed(true);
        }
    }
    stampPlyPlanes();
}

void Part::stampPlyPlanes() {
    if (_material.plies.size() < 2 || _tetElements.empty() || _xyzNodes.empty())
        return;
    using chrono::ChVector3d;
    ChVector3d normal(_material.stackNormal[0], _material.stackNormal[1], _material.stackNormal[2]);
    if (normal.Length2() < 0.25)
        normal = {0, 1, 0};
    else
        normal.Normalize();
    double lo = 1.0e30;
    double hi = -1.0e30;
    for (const auto& node : _xyzNodes) {
        if (!node)
            continue;
        const double s = node->GetPos().Dot(normal);
        lo = std::min(lo, s);
        hi = std::max(hi, s);
    }
    const double span = hi - lo;
    if (span < 1.0e-8)
        return;
    double thickSum = 0.0;
    for (const Ply& ply : _material.plies)
        thickSum += std::max(ply.thickness, 0.0);
    std::vector<double> planes;
    if (thickSum < 1.0e-8) {
        const double step = span / static_cast<double>(_material.plies.size());
        for (std::size_t i = 1; i < _material.plies.size(); ++i)
            planes.push_back(lo + step * static_cast<double>(i));
    } else {
        double cursor = lo;
        for (std::size_t i = 0; i + 1 < _material.plies.size(); ++i) {
            cursor += span * std::max(_material.plies[i].thickness, 0.0) / thickSum;
            planes.push_back(cursor);
        }
    }
    for (double plane : planes) {
        const ChVector3d origin = normal * plane;
        std::vector<std::size_t> hit;
        for (std::size_t t = 0; t < _tets.size() && t < _tetElements.size(); ++t) {
            if (!_tetElements[t] || (t < _tetLive.size() && !_tetLive[t]))
                continue;
            int pos = 0;
            int neg = 0;
            for (int node : _tets[t]) {
                if (node < 0 || static_cast<std::size_t>(node) >= _xyzNodes.size() || !_xyzNodes[static_cast<std::size_t>(node)])
                    continue;
                const double s = _xyzNodes[static_cast<std::size_t>(node)]->GetPos().Dot(normal) - plane;
                if (s > 1.0e-8)
                    ++pos;
                else if (s < -1.0e-8)
                    ++neg;
            }
            if (pos > 0 && neg > 0)
                hit.push_back(t);
        }
        std::map<std::pair<int, int>, int> shared;
        for (std::size_t tet : hit)
            splitTetOnPlane(tet, origin, normal, &shared);
    }
}

bool Part::maintainVolume() {
    if (_tetElements.empty())
        return false;
    return releaseMoltenTets(false);
}

int Part::xyzMidpoint(int a, int b) {
    using namespace chrono;
    using namespace chrono::fea;
    const std::size_t ia = static_cast<std::size_t>(a);
    const std::size_t ib = static_cast<std::size_t>(b);
    const auto& srcA = _xyzNodes[ia];
    const auto& srcB = _xyzNodes[ib];
    const ChVector3d mid = 0.5 * (srcA->GetPos() + srcB->GetPos());
    const ChVector3d midRest = 0.5 * (srcA->GetX0() + srcB->GetX0());
    auto node = chrono_types::make_shared<ChNodeFEAxyz>(mid);
    node->SetMass(0);
    node->SetFixed(_fixed);
    node->SetX0(midRest);
    node->SetPosDt(0.5 * (srcA->GetPosDt() + srcB->GetPosDt()));
    _mesh->AddNode(node);
    _xyzNodes.push_back(node);
    const int index = static_cast<int>(_xyzNodes.size()) - 1;
    _nodeMass.push_back(0.0);
    if (_loads.size() < _xyzNodes.size())
        _loads.resize(_xyzNodes.size(), ChVector3d(0, 0, 0));
    noteNewNode(index, 0.5 * (nodeTemperature(ia) + nodeTemperature(ib)));
    const bool boundaryA = std::find(_boundaryNodes.begin(), _boundaryNodes.end(), a) != _boundaryNodes.end();
    const bool boundaryB = std::find(_boundaryNodes.begin(), _boundaryNodes.end(), b) != _boundaryNodes.end();
    if (boundaryA && boundaryB)
        _boundaryNodes.push_back(index);
    return index;
}

int Part::xyzOnPlane(int a, int b, const ChVector3d& origin, const ChVector3d& normal) {
    using namespace chrono;
    using namespace chrono::fea;
    const std::size_t ia = static_cast<std::size_t>(a);
    const std::size_t ib = static_cast<std::size_t>(b);
    const auto& srcA = _xyzNodes[ia];
    const auto& srcB = _xyzNodes[ib];
    const ChVector3d pa = srcA->GetPos();
    const ChVector3d pb = srcB->GetPos();
    const double da = (pa - origin).Dot(normal);
    const double db = (pb - origin).Dot(normal);
    double t = 0.5;
    if (std::abs(da - db) > 1.0e-12)
        t = da / (da - db);
    t = std::clamp(t, 0.05, 0.95);
    const ChVector3d pos = pa + (pb - pa) * t;
    const ChVector3d rest = srcA->GetX0() + (srcB->GetX0() - srcA->GetX0()) * t;
    auto node = chrono_types::make_shared<ChNodeFEAxyz>(pos);
    node->SetMass(0);
    node->SetFixed(_fixed);
    node->SetX0(rest);
    node->SetPosDt(srcA->GetPosDt() + (srcB->GetPosDt() - srcA->GetPosDt()) * t);
    _mesh->AddNode(node);
    _xyzNodes.push_back(node);
    const int index = static_cast<int>(_xyzNodes.size()) - 1;
    _nodeMass.push_back(0.0);
    if (_loads.size() < _xyzNodes.size())
        _loads.resize(_xyzNodes.size(), ChVector3d(0, 0, 0));
    noteNewNode(index, nodeTemperature(ia) + (nodeTemperature(ib) - nodeTemperature(ia)) * t);
    const bool boundaryA = std::find(_boundaryNodes.begin(), _boundaryNodes.end(), a) != _boundaryNodes.end();
    const bool boundaryB = std::find(_boundaryNodes.begin(), _boundaryNodes.end(), b) != _boundaryNodes.end();
    if (boundaryA && boundaryB)
        _boundaryNodes.push_back(index);
    return index;
}

void Part::releaseTetAlongFace(std::size_t index, const chrono::ChVector3d& normal) {
    if (index >= _tets.size() || index >= _tetElements.size() || !_tetElements[index])
        return;
    if (index < _tetLive.size() && !_tetLive[index])
        return;
    using chrono::ChVector3d;
    const auto tet = _tets[index];
    const int faces[4][3] = {{0, 2, 1}, {0, 1, 3}, {0, 3, 2}, {1, 2, 3}};
    int best = -1;
    double bestDot = -1.0;
    for (int face = 0; face < 4; ++face) {
        bool valid = true;
        for (int corner : faces[face])
            valid = valid && tet[corner] >= 0 && static_cast<std::size_t>(tet[corner]) < _xyzNodes.size();
        if (!valid)
            continue;
        const ChVector3d a = _xyzNodes[static_cast<std::size_t>(tet[faces[face][0]])]->GetPos();
        const ChVector3d b = _xyzNodes[static_cast<std::size_t>(tet[faces[face][1]])]->GetPos();
        const ChVector3d c = _xyzNodes[static_cast<std::size_t>(tet[faces[face][2]])]->GetPos();
        ChVector3d n = (b - a).Cross(c - a);
        if (n.Length2() < 1.0e-24)
            continue;
        n.Normalize();
        const double align = std::abs(n.Dot(normal));
        if (align > bestDot) {
            bestDot = align;
            best = face;
        }
    }
    if (best < 0)
        return;
    int clones[3];
    for (int k = 0; k < 3; ++k)
        clones[k] = cloneNode(tet[static_cast<std::size_t>(faces[best][k])]);
    std::array<int, 4> next = tet;
    for (int k = 0; k < 3; ++k)
        next[static_cast<std::size_t>(faces[best][k])] = clones[k];
    _tets[index] = next;
    _tetElements[index]->SetNodes(_xyzNodes[static_cast<std::size_t>(next[0])],
                                  _xyzNodes[static_cast<std::size_t>(next[1])],
                                  _xyzNodes[static_cast<std::size_t>(next[2])],
                                  _xyzNodes[static_cast<std::size_t>(next[3])]);
    const int original[3] = {tet[static_cast<std::size_t>(faces[best][0])],
                             tet[static_cast<std::size_t>(faces[best][1])],
                             tet[static_cast<std::size_t>(faces[best][2])]};
    for (int face = 0; face < _skinTriangles && face < static_cast<int>(_triangles.size()); ++face) {
        int hits = 0;
        for (int corner : _triangles[static_cast<std::size_t>(face)])
            for (int node : original)
                hits += corner == node ? 1 : 0;
        if (hits == 3)
            _triangles[static_cast<std::size_t>(face)] = {clones[0], clones[1], clones[2]};
    }
}

bool Part::fractureCluster(const std::vector<std::size_t>& seeds, std::vector<char>& shard) {
    using chrono::ChVector3d;
    std::vector<char> in(_tets.size(), 0);
    std::vector<std::size_t> cluster;
    ChVector3d origin(0, 0, 0);
    ChVector3d axis(0, 0, 0);
    int samples = 0;
    for (std::size_t seed : seeds) {
        if (seed >= _tets.size() || seed >= _tetElements.size() || !_tetElements[seed])
            continue;
        if (seed < _tetLive.size() && !_tetLive[seed])
            continue;
        if (seed < in.size() && in[seed])
            continue;
        in[seed] = 1;
        cluster.push_back(seed);
        for (int node : _tets[seed]) {
            if (node < 0 || static_cast<std::size_t>(node) >= _xyzNodes.size())
                continue;
            origin += _xyzNodes[static_cast<std::size_t>(node)]->GetPos();
            ++samples;
        }
        axis += fabricCrackNormal(_material, _tetElements[seed]->tensionAxis());
    }
    if (cluster.empty() || samples == 0)
        return false;
    origin *= 1.0 / static_cast<double>(samples);
    if (axis.Length2() < 0.25)
        axis = {1, 0, 0};
    else
        axis.Normalize();

    auto shareFace = [](const std::array<int, 4>& a, const std::array<int, 4>& b) {
        int hits = 0;
        for (int i : a)
            for (int j : b)
                hits += i == j ? 1 : 0;
        return hits >= 3;
    };
    auto crosses = [&](std::size_t tet) {
        int pos = 0;
        int neg = 0;
        for (int node : _tets[tet]) {
            if (node < 0 || static_cast<std::size_t>(node) >= _xyzNodes.size())
                continue;
            const double s = (_xyzNodes[static_cast<std::size_t>(node)]->GetPos() - origin).Dot(axis);
            if (s > 1.0e-10)
                ++pos;
            else if (s < -1.0e-10)
                ++neg;
        }
        return pos > 0 && neg > 0;
    };
    for (std::size_t k = 0; k < cluster.size() && cluster.size() < 24; ++k) {
        for (std::size_t other = 0; other < _tets.size(); ++other) {
            if (other < in.size() && in[other])
                continue;
            if (other >= _tetElements.size() || !_tetElements[other])
                continue;
            if (other < _tetLive.size() && !_tetLive[other])
                continue;
            if (!shareFace(_tets[cluster[k]], _tets[other]) || !crosses(other))
                continue;
            if (other >= in.size())
                in.resize(other + 1, 0);
            in[other] = 1;
            cluster.push_back(other);
        }
    }

    auto usesEdge = [](const std::array<int, 4>& tet, int a, int b) {
        bool hasA = false;
        bool hasB = false;
        for (int node : tet) {
            hasA = hasA || node == a;
            hasB = hasB || node == b;
        }
        return hasA && hasB;
    };
    std::map<std::pair<int, int>, int> cuts;
    for (std::size_t tet : cluster) {
        if (!crosses(tet))
            continue;
        const auto corners = _tets[tet];
        double s[4];
        for (int k = 0; k < 4; ++k) {
            const int node = corners[static_cast<std::size_t>(k)];
            s[k] = node >= 0 && static_cast<std::size_t>(node) < _xyzNodes.size()
                       ? (_xyzNodes[static_cast<std::size_t>(node)]->GetPos() - origin).Dot(axis)
                       : 0.0;
        }
        const int edge[6][2] = {{0, 1}, {0, 2}, {0, 3}, {1, 2}, {1, 3}, {2, 3}};
        for (const auto& pair : edge) {
            if ((s[pair[0]] > 0.0) == (s[pair[1]] > 0.0))
                continue;
            int a = corners[static_cast<std::size_t>(pair[0])];
            int b = corners[static_cast<std::size_t>(pair[1])];
            if (a > b)
                std::swap(a, b);
            if (cuts.count({a, b}))
                continue;
            cuts.emplace(std::make_pair(a, b), xyzOnPlane(a, b, origin, axis));
        }
    }
    for (const auto& edge : cuts) {
        std::vector<std::size_t> users;
        for (std::size_t tet = 0; tet < _tets.size(); ++tet) {
            if (tet < in.size() && in[tet])
                continue;
            if (tet >= _tetElements.size() || !_tetElements[tet])
                continue;
            if (tet < _tetLive.size() && !_tetLive[tet])
                continue;
            if (usesEdge(_tets[tet], edge.first.first, edge.first.second))
                users.push_back(tet);
        }
        for (std::size_t tet : users) {
            if (tet < _tets.size() && usesEdge(_tets[tet], edge.first.first, edge.first.second))
                bisectTet(tet, edge.first.first, edge.first.second, edge.second);
        }
    }

    std::vector<std::size_t> opened;
    std::vector<std::size_t> fallback;
    for (std::size_t tet : cluster) {
        if (tet >= _tets.size() || tet >= _tetElements.size() || !_tetElements[tet])
            continue;
        if (tet < _tetLive.size() && !_tetLive[tet])
            continue;
        const auto kids = splitTetOnPlane(tet, origin, axis, &cuts);
        if (kids.size() < 2)
            fallback.push_back(tet);
        else
            opened.insert(opened.end(), kids.begin(), kids.end());
    }
    for (std::size_t tet : fallback) {
        releaseTetAlongFace(tet, axis);
        if (tet < shard.size())
            shard[tet] = 0;
    }
    if (opened.empty())
        return !fallback.empty();

    std::vector<char> child(std::max(_tets.size(), opened.back() + 1), 0);
    for (std::size_t tet : opened) {
        if (tet >= child.size())
            child.resize(tet + 1, 0);
        child[tet] = 1;
    }
    std::set<int> midpoints;
    for (const auto& edge : cuts)
        midpoints.insert(edge.second);
    std::set<int> front;
    for (std::size_t tet = 0; tet < _tets.size() && tet < _tetElements.size(); ++tet) {
        if (!_tetElements[tet] || (tet < _tetLive.size() && !_tetLive[tet]))
            continue;
        if (tet < child.size() && child[tet])
            continue;
        for (int node : _tets[tet])
            if (midpoints.count(node))
                front.insert(node);
    }
    std::map<int, int> cloned;
    std::vector<std::size_t> positive;
    for (std::size_t tet : opened) {
        if (tet >= _tets.size() || tet >= _tetElements.size() || !_tetElements[tet])
            continue;
        ChVector3d c(0, 0, 0);
        int corners = 0;
        for (int node : _tets[tet]) {
            if (node < 0 || static_cast<std::size_t>(node) >= _xyzNodes.size())
                continue;
            c += _xyzNodes[static_cast<std::size_t>(node)]->GetPos();
            ++corners;
        }
        if (corners == 0)
            continue;
        c *= 1.0 / static_cast<double>(corners);
        if ((c - origin).Dot(axis) < 0.0) {
            if (tet < shard.size())
                shard[tet] = 0;
            continue;
        }
        positive.push_back(tet);
        auto next = _tets[tet];
        bool moved = false;
        for (int& node : next) {
            if (!midpoints.count(node) || front.count(node))
                continue;
            if (!cloned.count(node))
                cloned.emplace(node, cloneNode(node));
            node = cloned[node];
            moved = true;
        }
        if (!moved)
            continue;
        _tets[tet] = next;
        _tetElements[tet]->SetNodes(_xyzNodes[static_cast<std::size_t>(next[0])],
                                    _xyzNodes[static_cast<std::size_t>(next[1])],
                                    _xyzNodes[static_cast<std::size_t>(next[2])],
                                    _xyzNodes[static_cast<std::size_t>(next[3])]);
    }

    double volume = 0.0;
    for (std::size_t tet : positive) {
        if (tet < _tetElements.size() && _tetElements[tet])
            volume += std::abs(_tetElements[tet]->volume());
    }
    const double speck = 0.125 * _material.minimumShardSize * _material.minimumShardSize * _material.minimumShardSize;
    const bool erase = _material.shatterToDust || (_material.minimumShardSize > 0.0 && volume < speck);
    for (std::size_t tet : positive) {
        if (shard.size() <= tet)
            shard.resize(tet + 1, 0);
        shard[tet] = erase ? 1 : 0;
        if (tet < _tetElements.size() && _tetElements[tet])
            _tetElements[tet]->commitDamage(1.0);
    }
    return true;
}

std::vector<std::size_t> Part::splitTetOnTensionPlane(std::size_t index) {
    using chrono::ChVector3d;
    if (index >= _tetElements.size() || !_tetElements[index])
        return {};
    ChVector3d centroid(0, 0, 0);
    if (index >= _tets.size())
        return {};
    for (int node : _tets[index]) {
        if (node < 0 || static_cast<std::size_t>(node) >= _xyzNodes.size())
            return {};
        centroid += _xyzNodes[static_cast<std::size_t>(node)]->GetPos();
    }
    centroid *= 0.25;
    return splitTetOnPlane(index, centroid, fabricCrackNormal(_material, _tetElements[index]->tensionAxis()), nullptr);
}

std::vector<std::size_t> Part::splitTetOnPlane(std::size_t index, const chrono::ChVector3d& origin,
                                              const chrono::ChVector3d& axisIn,
                                              std::map<std::pair<int, int>, int>* cuts) {
    using namespace chrono;
    using namespace chrono::fea;
    std::vector<std::size_t> children;
    if (index >= _tets.size() || index >= _tetElements.size() || !_tetElements[index])
        return children;
    if (index < _tetLive.size() && !_tetLive[index])
        return children;
    const auto parent = _tets[index];
    const auto parentElement = _tetElements[index];
    ChVector3d axis = axisIn;
    if (axis.Length2() < 0.25)
        axis = {1, 0, 0};
    else
        axis.Normalize();
    ChVector3d p[4];
    for (int k = 0; k < 4; ++k)
        p[k] = _xyzNodes[static_cast<std::size_t>(parent[static_cast<std::size_t>(k)])]->GetPos();
    double s[4];
    int pos = 0;
    int neg = 0;
    int isolated = -1;
    for (int k = 0; k < 4; ++k) {
        s[k] = (p[k] - origin).Dot(axis);
        if (s[k] > 1.0e-10) {
            ++pos;
            isolated = k;
        } else if (s[k] < -1.0e-10) {
            ++neg;
            isolated = k;
        }
    }
    if (pos == 0 || neg == 0)
        return children;
    if (!((pos == 1 && neg == 3) || (pos == 3 && neg == 1) || (pos == 2 && neg == 2)))
        return children;

    auto shapeOk = [](const ChVector3d& a, const ChVector3d& b, const ChVector3d& c, const ChVector3d& d) {
        const double volume = std::abs((b - a).Dot((c - a).Cross(d - a))) / 6.0;
        const ChVector3d q[4] = {a, b, c, d};
        double longest = 0.0;
        double shortest = 1.0e30;
        double sum = 0.0;
        for (int i = 0; i < 4; ++i) {
            for (int j = i + 1; j < 4; ++j) {
                const double len = (q[i] - q[j]).Length();
                longest = std::max(longest, len);
                shortest = std::min(shortest, len);
                sum += len;
            }
        }
        if (shortest < 1.0e-12 || longest > 20.0 * shortest)
            return false;
        const double mean = sum / 6.0;
        return mean > 0.0 && volume / (mean * mean * mean) >= 1.0e-4;
    };
    auto cutPoint = [&](int i, int j) {
        const double denom = s[i] - s[j];
        double t = std::abs(denom) > 1.0e-12 ? s[i] / denom : 0.5;
        t = std::clamp(t, 0.05, 0.95);
        return p[i] + (p[j] - p[i]) * t;
    };
    bool sound = true;
    if ((pos == 1 && neg == 3) || (pos == 3 && neg == 1)) {
        int others[3];
        int n = 0;
        for (int k = 0; k < 4; ++k)
            if (k != isolated)
                others[n++] = k;
        if (n == 3) {
            const ChVector3d m0 = cutPoint(isolated, others[0]);
            const ChVector3d m1 = cutPoint(isolated, others[1]);
            const ChVector3d m2 = cutPoint(isolated, others[2]);
            sound = shapeOk(p[isolated], m0, m1, m2) && shapeOk(m0, p[others[0]], p[others[1]], p[others[2]]) &&
                    shapeOk(m1, m0, p[others[1]], p[others[2]]) && shapeOk(m2, m0, m1, p[others[2]]);
        }
    } else {
        int a[2];
        int b[2];
        int na = 0;
        int nb = 0;
        for (int k = 0; k < 4; ++k) {
            if (s[k] > 0.0)
                a[na++] = k;
            else
                b[nb++] = k;
        }
        if (na == 2 && nb == 2) {
            const ChVector3d p02 = cutPoint(a[0], b[0]);
            const ChVector3d p03 = cutPoint(a[0], b[1]);
            const ChVector3d p12 = cutPoint(a[1], b[0]);
            const ChVector3d p13 = cutPoint(a[1], b[1]);
            sound = shapeOk(p[a[0]], p[a[1]], p02, p12) && shapeOk(p[a[0]], p02, p03, p12) &&
                    shapeOk(p[a[0]], p03, p13, p12) && shapeOk(p[b[0]], p[b[1]], p02, p12) &&
                    shapeOk(p[b[0]], p02, p03, p12) && shapeOk(p[b[0]], p03, p13, p12);
        }
    }
    if (!sound)
        return children;

    parentElement->commitPlastic();
    const ChVector3d p0 = p[0];
    const ChVector3d p1 = p[1];
    const ChVector3d p2 = p[2];
    const ChVector3d p3 = p[3];
    const double oldVolume = std::abs((p1 - p0).Dot((p2 - p0).Cross(p3 - p0))) / 6.0;
    const double oldShare = _material.density * oldVolume / 4.0;
    for (int node : parent)
        if (node >= 0 && static_cast<std::size_t>(node) < _nodeMass.size())
            _nodeMass[static_cast<std::size_t>(node)] =
                std::max(0.0, _nodeMass[static_cast<std::size_t>(node)] - oldShare);

    auto addChild = [&](int n0, int n1, int n2, int n3) -> std::size_t {
        const ChVector3d q0 = _xyzNodes[static_cast<std::size_t>(n0)]->GetPos();
        const ChVector3d q1 = _xyzNodes[static_cast<std::size_t>(n1)]->GetPos();
        const ChVector3d q2 = _xyzNodes[static_cast<std::size_t>(n2)]->GetPos();
        const ChVector3d q3 = _xyzNodes[static_cast<std::size_t>(n3)]->GetPos();
        double volume = (q1 - q0).Dot((q2 - q0).Cross(q3 - q0));
        if (std::abs(volume) < 1.0e-18)
            return static_cast<std::size_t>(-1);
        if (volume < 0.0) {
            std::swap(n2, n3);
            volume = -volume;
        }
        auto element = chrono_types::make_shared<ThermalTet4>();
        element->configure(_material.youngsModulus, _material.poissonRatio,
                           _material.allowPlasticity ? _material.yieldStress : 0.0, _material.density,
                           _material.dampingBeta);
        element->setHardening(_material.ultimateStress);
        if (_material.fracture == FractureMode::Brittle)
            element->setRankine(rankineTensionLimit(_material), rankineCompressionLimit(_material));
        if (_material.fracture == FractureMode::Ductile)
            element->setFailureStrain(_material.failPlasticStrain);
        element->SetNodes(_xyzNodes[static_cast<std::size_t>(n0)], _xyzNodes[static_cast<std::size_t>(n1)],
                          _xyzNodes[static_cast<std::size_t>(n2)], _xyzNodes[static_cast<std::size_t>(n3)]);
        element->copyPlasticFrom(*parentElement);
        _mesh->AddElement(element);
        const double share = _material.density * (volume / 6.0) / 4.0;
        for (int node : {n0, n1, n2, n3})
            if (static_cast<std::size_t>(node) < _nodeMass.size())
                _nodeMass[static_cast<std::size_t>(node)] += share;
        if (children.empty()) {
            _tets[index] = {n0, n1, n2, n3};
            _tetElements[index] = element;
            children.push_back(index);
            return index;
        }
        _tets.push_back({n0, n1, n2, n3});
        _tetElements.push_back(element);
        _tetLive.push_back(1);
        children.push_back(_tets.size() - 1);
        return _tets.size() - 1;
    };

    std::map<std::pair<int, int>, int> localCuts;
    std::map<std::pair<int, int>, int>& sharedCuts = cuts ? *cuts : localCuts;
    auto cutEdge = [&](int i, int j) {
        int a = parent[static_cast<std::size_t>(i)];
        int b = parent[static_cast<std::size_t>(j)];
        if (a > b)
            std::swap(a, b);
        const auto key = std::make_pair(a, b);
        const auto found = sharedCuts.find(key);
        if (found != sharedCuts.end())
            return found->second;
        const int mid = xyzOnPlane(a, b, origin, axis);
        sharedCuts.emplace(key, mid);
        return mid;
    };

    if ((pos == 1 && neg == 3) || (pos == 3 && neg == 1)) {
        int others[3];
        int n = 0;
        for (int k = 0; k < 4; ++k)
            if (k != isolated)
                others[n++] = k;
        if (n != 3)
            return children;
        const int m0 = cutEdge(isolated, others[0]);
        const int m1 = cutEdge(isolated, others[1]);
        const int m2 = cutEdge(isolated, others[2]);
        addChild(parent[static_cast<std::size_t>(isolated)], m0, m1, m2);
        addChild(m0, parent[static_cast<std::size_t>(others[0])], parent[static_cast<std::size_t>(others[1])],
                 parent[static_cast<std::size_t>(others[2])]);
        addChild(m1, m0, parent[static_cast<std::size_t>(others[1])], parent[static_cast<std::size_t>(others[2])]);
        addChild(m2, m0, m1, parent[static_cast<std::size_t>(others[2])]);
    } else if (pos == 2 && neg == 2) {
        int a[2];
        int b[2];
        int na = 0;
        int nb = 0;
        for (int k = 0; k < 4; ++k) {
            if (s[k] > 0.0)
                a[na++] = k;
            else
                b[nb++] = k;
        }
        if (na != 2 || nb != 2)
            return children;
        const int p02 = cutEdge(a[0], b[0]);
        const int p03 = cutEdge(a[0], b[1]);
        const int p12 = cutEdge(a[1], b[0]);
        const int p13 = cutEdge(a[1], b[1]);
        addChild(parent[static_cast<std::size_t>(a[0])], parent[static_cast<std::size_t>(a[1])], p02, p12);
        addChild(parent[static_cast<std::size_t>(a[0])], p02, p03, p12);
        addChild(parent[static_cast<std::size_t>(a[0])], p03, p13, p12);
        addChild(parent[static_cast<std::size_t>(b[0])], parent[static_cast<std::size_t>(b[1])], p02, p12);
        addChild(parent[static_cast<std::size_t>(b[0])], p02, p03, p12);
        addChild(parent[static_cast<std::size_t>(b[0])], p03, p13, p12);
    }
    if (children.size() < 2)
        return {};
    if (_mesh)
        _mesh->RemoveElement(parentElement);

    const int corner[6][2] = {{0, 1}, {0, 2}, {0, 3}, {1, 2}, {1, 3}, {2, 3}};
    for (const auto& e : corner) {
        const int ia = parent[static_cast<std::size_t>(e[0])];
        const int ib = parent[static_cast<std::size_t>(e[1])];
        if ((s[e[0]] > 0.0) == (s[e[1]] > 0.0))
            continue;
        for (int face = 0; face < _skinTriangles && face < static_cast<int>(_triangles.size()); ++face) {
            auto& tri = _triangles[static_cast<std::size_t>(face)];
            int hits = 0;
            for (int cornerNode : tri)
                hits += (cornerNode == ia || cornerNode == ib) ? 1 : 0;
            if (hits != 2)
                continue;
            int keep = -1;
            for (int cornerNode : tri)
                if (cornerNode != ia && cornerNode != ib)
                    keep = cornerNode;
            if (keep < 0)
                continue;
            int lo = ia;
            int hi = ib;
            if (lo > hi)
                std::swap(lo, hi);
            int mid = -1;
            const auto found = sharedCuts.find({lo, hi});
            if (found != sharedCuts.end())
                mid = found->second;
            else
                mid = xyzOnPlane(ia, ib, origin, axis);
            tri = {ia, mid, keep};
            _triangles.insert(_triangles.begin() + face + 1, {mid, ib, keep});
            ++_skinTriangles;
            ++face;
        }
    }
    return children;
}

void Part::bisectTet(std::size_t index, int a, int b, int mid) {
    if (index >= _tets.size() || index >= _tetElements.size() || !_tetElements[index])
        return;
    const auto parent = _tets[index];
    const auto parentElement = _tetElements[index];
    parentElement->commitPlastic();
    std::array<int, 4> first = parent;
    std::array<int, 4> second = parent;
    bool sawA = false;
    bool sawB = false;
    for (int i = 0; i < 4; ++i) {
        if (parent[static_cast<std::size_t>(i)] == a) {
            first[static_cast<std::size_t>(i)] = mid;
            sawA = true;
        }
        if (parent[static_cast<std::size_t>(i)] == b) {
            second[static_cast<std::size_t>(i)] = mid;
            sawB = true;
        }
    }
    if (!sawA || !sawB)
        return;

    auto addChild = [&](int n0, int n1, int n2, int n3) {
        auto element = chrono_types::make_shared<ThermalTet4>();
        element->configure(_material.youngsModulus, _material.poissonRatio,
                           _material.allowPlasticity ? _material.yieldStress : 0.0, _material.density,
                           _material.dampingBeta);
        element->setHardening(_material.ultimateStress);
        if (_material.fracture == FractureMode::Brittle)
            element->setRankine(rankineTensionLimit(_material), rankineCompressionLimit(_material));
        if (_material.fracture == FractureMode::Ductile)
            element->setFailureStrain(_material.failPlasticStrain);
        element->SetNodes(_xyzNodes[static_cast<std::size_t>(n0)], _xyzNodes[static_cast<std::size_t>(n1)],
                          _xyzNodes[static_cast<std::size_t>(n2)], _xyzNodes[static_cast<std::size_t>(n3)]);
        element->copyPlasticFrom(*parentElement);
        _mesh->AddElement(element);
        const ChVector3d p0 = _xyzNodes[static_cast<std::size_t>(n0)]->GetPos();
        const ChVector3d p1 = _xyzNodes[static_cast<std::size_t>(n1)]->GetPos();
        const ChVector3d p2 = _xyzNodes[static_cast<std::size_t>(n2)]->GetPos();
        const ChVector3d p3 = _xyzNodes[static_cast<std::size_t>(n3)]->GetPos();
        const double volume = std::abs((p1 - p0).Dot((p2 - p0).Cross(p3 - p0))) / 6.0;
        const double share = _material.density * volume / 4.0;
        for (int node : {n0, n1, n2, n3})
            if (static_cast<std::size_t>(node) < _nodeMass.size())
                _nodeMass[static_cast<std::size_t>(node)] += share;
        return std::pair<std::array<int, 4>, std::shared_ptr<ThermalTet4>>{{n0, n1, n2, n3}, element};
    };

    const ChVector3d p0 = _xyzNodes[static_cast<std::size_t>(parent[0])]->GetPos();
    const ChVector3d p1 = _xyzNodes[static_cast<std::size_t>(parent[1])]->GetPos();
    const ChVector3d p2 = _xyzNodes[static_cast<std::size_t>(parent[2])]->GetPos();
    const ChVector3d p3 = _xyzNodes[static_cast<std::size_t>(parent[3])]->GetPos();
    const double oldVolume = std::abs((p1 - p0).Dot((p2 - p0).Cross(p3 - p0))) / 6.0;
    const double oldShare = _material.density * oldVolume / 4.0;
    for (int node : parent)
        if (node >= 0 && static_cast<std::size_t>(node) < _nodeMass.size())
            _nodeMass[static_cast<std::size_t>(node)] = std::max(0.0, _nodeMass[static_cast<std::size_t>(node)] - oldShare);

    const auto childA = addChild(first[0], first[1], first[2], first[3]);
    const auto childB = addChild(second[0], second[1], second[2], second[3]);
    _tets[index] = childA.first;
    _tetElements[index] = childA.second;
    _tets.push_back(childB.first);
    _tetElements.push_back(childB.second);
    _tetLive.push_back(1);
    if (_mesh)
        _mesh->RemoveElement(parentElement);

    for (int face = 0; face < _skinTriangles && face < static_cast<int>(_triangles.size()); ++face) {
        auto& tri = _triangles[static_cast<std::size_t>(face)];
        int hits = 0;
        for (int corner : tri)
            hits += (corner == a || corner == b) ? 1 : 0;
        if (hits != 2)
            continue;
        int keep = -1;
        for (int corner : tri)
            if (corner != a && corner != b)
                keep = corner;
        if (keep < 0)
            continue;
        tri = {a, mid, keep};
        _triangles.insert(_triangles.begin() + face + 1, {mid, b, keep});
        ++_skinTriangles;
        ++face;
    }
}

bool Part::releaseMoltenTets(bool includeBrittle) {
    if (_tetElements.empty())
        return false;

    auto meanFaceArea = [&](const std::array<int, 4>& tet) {
        auto tri = [&](int i, int j, int k) {
            const chrono::ChVector3d a = nodePosition(static_cast<std::size_t>(tet[static_cast<std::size_t>(i)]));
            const chrono::ChVector3d b = nodePosition(static_cast<std::size_t>(tet[static_cast<std::size_t>(j)]));
            const chrono::ChVector3d c = nodePosition(static_cast<std::size_t>(tet[static_cast<std::size_t>(k)]));
            return 0.5 * (b - a).Cross(c - a).Length();
        };
        return 0.25 * (tri(0, 1, 2) + tri(0, 1, 3) + tri(0, 2, 3) + tri(1, 2, 3));
    };
    auto brittleOpen = [&](std::size_t t, double kelvin) {
        if (!includeBrittle)
            return false;
        return brittleShouldOpen(_material, kelvin, _tetElements[t]->rawTension(),
                                 _tetElements[t]->trialEnergy(), _tetElements[t]->tensileStrain(),
                                 meanFaceArea(_tets[t]), _tetElements[t]->rawCompression(),
                                 _tetElements[t]->rawMean(), _tetElements[t]->tensionAxis());
    };

    struct Candidate {
        std::size_t index = 0;
        double score = 0.0;
        int kind = 2;
    };
    std::vector<Candidate> candidates;
    for (std::size_t t = 0; t < _tets.size() && t < _tetElements.size(); ++t) {
        if (!_tetElements[t] || (t < _tetLive.size() && !_tetLive[t]))
            continue;
        const auto& tet = _tets[t];
        int corners = 0;
        bool molten = true;
        double kelvin = 0.0;
        for (int index : tet) {
            if (index < 0) {
                molten = false;
                continue;
            }
            ++corners;
            kelvin += nodeTemperature(static_cast<std::size_t>(index));
            if (nodeTemperature(static_cast<std::size_t>(index)) < _material.meltTemperature)
                molten = false;
        }
        if (corners != 4)
            continue;
        kelvin /= 4.0;
        const bool crush =
            includeBrittle && brittleShouldCrush(_material, kelvin, _tetElements[t]->rawCompression());
        const bool invertedFail = includeBrittle && _material.fracture == FractureMode::Brittle &&
                                 (_tetElements[t]->inverted() || _tetElements[t]->deformationSingular());
        const bool opened =
            includeBrittle && _material.fracture == FractureMode::Brittle &&
            (invertedFail || brittleOpen(t, kelvin));
        if (includeBrittle && _material.fracture == FractureMode::Brittle)
            _tetElements[t]->commitDamage((opened || crush) ? 1.0 : 0.0);
        const double ductility = thermalDuctilityRetention(_material, kelvin);
        const double grain = grainOpeningScale(_material, _tetElements[t]->tensionAxis());
        const double plasticLimit = _material.failPlasticStrain * ductility * grain;
        const double stretchLimit = 1.0 + (liveStretchLimit(_material, ductility) - 1.0) * grain;
        const double plastic = _tetElements[t]->plasticStrain();
        const double stretch = maxEdgeStretch(tet.data(), 4);
        const double tension = _tetElements[t]->tension();
        double softest = 1.0;
        for (int index : tet) {
            if (index < 0)
                continue;
            softest = std::min(softest, thermalStiffnessRetention(_material, nodeTemperature(static_cast<std::size_t>(index))));
        }
        double score = 0.0;
        int kind = 2;
        // Only a fully molten tet leaves as a chip. One hot corner used to
        // mark every neighbour as melt (softest <= thermal floor) and the
        // core vanished with the surface.
        if (molten) {
            score = 1.0e6 + tension;
            kind = 0;
        } else if (!includeBrittle) {
            continue;
        } else if (interfaceShouldOpen(_material, kelvin, tension, _tetElements[t]->tensionAxis(),
                                      _tetElements[t]->trialEnergy(), meanFaceArea(_tets[t]))) {
            score = tension;
            kind = 5;
        } else if (opened) {
            // Tensile Rankine chips leave even when the failed cluster is a
            // thin-sheet corner smaller than minimumShardSize³. That cube-volume
            // floor is what kept an over-strength pane in the floor.
            score = std::max({_tetElements[t]->trialEnergy(), _tetElements[t]->rawTension(),
                              _tetElements[t]->rawCompression()}) /
                    std::max(grain, 1.0e-6);
            kind = 1;
        } else if (crush) {
            // Cluster volume vs minimumShardSize³ is the crush-powder gate:
            // B3's tiny hoop zone stays, B4's crushed neighbourhood releases.
            score = std::max({_tetElements[t]->trialEnergy(), _tetElements[t]->rawTension(),
                              _tetElements[t]->rawCompression()});
            kind = 4;
        } else if (_material.fracture == FractureMode::Ductile) {
            if (((plasticLimit > 0.0 && plastic >= plasticLimit) || stretch >= stretchLimit) &&
                energyAllowsOpening(_material, _tetElements[t]->trialEnergy(), meanFaceArea(_tets[t])))
                score = plastic * 1.0e9 + tension;
        } else if (_material.fracture != FractureMode::Brittle && stretch >= stretchLimit &&
                   energyAllowsOpening(_material, _tetElements[t]->trialEnergy(), meanFaceArea(_tets[t]))) {
            score = stretch;
        }
        if (score > 0.0)
            candidates.push_back({t, score, kind});
    }
    if (!_material.shatterToDust && _material.minimumShardSize > 0.0) {
        const double minVol = 0.125 * _material.minimumShardSize * _material.minimumShardSize *
                              _material.minimumShardSize;
        std::vector<char> brittle(_tets.size(), 0);
        for (const auto& candidate : candidates)
            if (candidate.kind == 4 && candidate.index < brittle.size())
                brittle[candidate.index] = 1;
        std::vector<int> parent(_tets.size());
        for (int i = 0; i < static_cast<int>(_tets.size()); ++i)
            parent[static_cast<std::size_t>(i)] = i;
        auto shareFace = [](const std::array<int, 4>& a, const std::array<int, 4>& b) {
            int hits = 0;
            for (int i : a)
                for (int j : b)
                    hits += i == j ? 1 : 0;
            return hits >= 3;
        };
        for (std::size_t i = 0; i < _tets.size(); ++i) {
            if (!brittle[i])
                continue;
            for (std::size_t j = i + 1; j < _tets.size(); ++j) {
                if (!brittle[j])
                    continue;
                if (shareFace(_tets[i], _tets[j]))
                    unite(parent, static_cast<int>(i), static_cast<int>(j));
            }
        }
        std::vector<double> groupVol(_tets.size(), 0.0);
        for (std::size_t i = 0; i < _tets.size() && i < _tetElements.size(); ++i) {
            if (!brittle[i] || !_tetElements[i])
                continue;
            groupVol[static_cast<std::size_t>(findRoot(parent, static_cast<int>(i)))] += _tetElements[i]->volume();
        }
        candidates.erase(std::remove_if(candidates.begin(), candidates.end(),
                                        [&](const Candidate& c) {
                                            if (c.kind != 4 || c.index >= groupVol.size())
                                                return false;
                                            const int root = findRoot(parent, static_cast<int>(c.index));
                                            return groupVol[static_cast<std::size_t>(root)] < minVol;
                                        }),
                         candidates.end());
    }
    std::sort(candidates.begin(), candidates.end(),
              [](const Candidate& a, const Candidate& b) { return a.score > b.score; });

    std::vector<char> chosen(_tets.size(), 0);
    int meltTaken = 0;
    int otherTaken = 0;
    for (const auto& candidate : candidates) {
        if (candidate.kind == 0) {
            if (meltTaken >= 48)
                continue;
            ++meltTaken;
        } else if (candidate.kind == 1 || candidate.kind == 4) {
            // Brittle count follows the load, not a scene-size budget.
        } else {
            if (otherTaken >= 16)
                continue;
            ++otherTaken;
        }
        chosen[candidate.index] = 1;
    }

    bool changed = false;
    {
        std::vector<std::size_t> interfaces;
        std::vector<std::size_t> opening;
        for (const auto& candidate : candidates) {
            if (candidate.index >= chosen.size() || !chosen[candidate.index])
                continue;
            if (candidate.kind == 5)
                interfaces.push_back(candidate.index);
            else if (candidate.kind == 1 || candidate.kind == 2)
                opening.push_back(candidate.index);
        }
        chrono::ChVector3d stack(_material.stackNormal[0], _material.stackNormal[1], _material.stackNormal[2]);
        if (stack.Length2() < 0.25)
            stack = {0, 1, 0};
        else
            stack.Normalize();
        for (std::size_t tet : interfaces) {
            releaseTetAlongFace(tet, stack);
            if (tet < chosen.size())
                chosen[tet] = 0;
            changed = true;
        }
        if (!opening.empty()) {
            if (chosen.size() < _tets.size())
                chosen.resize(_tets.size(), 0);
            if (fractureCluster(opening, chosen))
                changed = true;
        }
    }

    const std::size_t shardsBefore = _shards.size();
    for (std::size_t t = 0; t < _tets.size() && t < _tetElements.size(); ++t) {
        if (t >= chosen.size() || !chosen[t])
            continue;
        if (!_tetElements[t] || (t < _tetLive.size() && !_tetLive[t]))
            continue;
        const auto& tet = _tets[t];
        int corners = 0;
        bool molten = true;
        double kelvin = 0.0;
        for (int index : tet) {
            if (index < 0) {
                molten = false;
                continue;
            }
            ++corners;
            const double nodeT = nodeTemperature(static_cast<std::size_t>(index));
            kelvin += nodeT;
            if (nodeT < _material.meltTemperature)
                molten = false;
        }
        if (corners != 4)
            continue;
        kelvin /= 4.0;
        const double ductility = thermalDuctilityRetention(_material, kelvin);
        const double grain = grainOpeningScale(_material, _tetElements[t]->tensionAxis());
        const double plasticLimit = _material.failPlasticStrain * ductility * grain;
        const double stretchLimit = 1.0 + (liveStretchLimit(_material, ductility) - 1.0) * grain;
        const double plastic = _tetElements[t]->plasticStrain();
        const double stretch = maxEdgeStretch(tet.data(), 4);
        bool failed = molten;
        if (!failed && _material.fracture == FractureMode::Brittle && _material.fractureStrength > 0.0)
            failed = _tetElements[t]->damage() >= 1.0;
        if (!failed && _material.fracture == FractureMode::Ductile)
            failed = (plasticLimit > 0.0 && plastic >= plasticLimit) || stretch >= stretchLimit;
        if (!failed && _material.fracture == FractureMode::Compliant)
            failed = stretch >= stretchLimit;
        if (!failed)
            continue;

        int clones[4];
        for (int k = 0; k < 4; ++k)
            clones[k] = cloneNode(tet[k]);
        const int faces[4][3] = {{0, 2, 1}, {0, 1, 3}, {0, 3, 2}, {1, 2, 3}};
        for (const auto& face : faces) {
            const int original[3] = {tet[face[0]], tet[face[1]], tet[face[2]]};
            int skin = -1;
            for (int f = 0; f < _skinTriangles && f < static_cast<int>(_triangles.size()); ++f) {
                int hits = 0;
                for (int corner : _triangles[static_cast<std::size_t>(f)])
                    for (int node : original)
                        hits += corner == node ? 1 : 0;
                if (hits == 3) {
                    skin = f;
                    break;
                }
            }
            if (skin >= 0) {
                _triangles[static_cast<std::size_t>(skin)] = {clones[face[0]], clones[face[1]], clones[face[2]]};
            } else {
                _triangles.insert(_triangles.begin() + _skinTriangles, {original[0], original[1], original[2]});
                ++_skinTriangles;
                _triangles.push_back({clones[face[0]], clones[face[1]], clones[face[2]]});
            }
        }
        const std::array<int, 4> originals = {tet[0], tet[1], tet[2], tet[3]};
        _tets[t] = {clones[0], clones[1], clones[2], clones[3]};
        if (_mesh)
            _mesh->RemoveElement(_tetElements[t]);
        _tetElements[t].reset();
        if (t < _tetLive.size())
            _tetLive[t] = 0;
        // The hole is not a 2000 K furnace. Neighbours that kept these
        // corners used to melt in the 0.4 s fall and spawn extra chips.
        const double belowMelt = _material.meltTemperature - 25.0;
        for (int orig : originals) {
            if (orig < 0)
                continue;
            const std::size_t node = static_cast<std::size_t>(orig);
            if (nodeTemperature(node) > belowMelt)
                setNodeTemperature(node, belowMelt);
            if (node < _strainTemperature.size() && strainTemperature(node) > belowMelt)
                setStrainTemperature(node, belowMelt);
        }

        using namespace chrono;
        RigidShard shard;
        for (int k = 0; k < 4; ++k) {
            const int index = clones[k];
            const auto& node = _xyzNodes[static_cast<std::size_t>(index)];
            shard.nodes[static_cast<std::size_t>(k)] = index;
            shard.com += node->GetPos();
            shard.velocity += node->GetPosDt();
            if (static_cast<std::size_t>(index) < _nodeMass.size())
                shard.mass += _nodeMass[static_cast<std::size_t>(index)];
        }
        shard.com *= 0.25;
        shard.velocity *= 0.25;
        if (shard.velocity.Length() > 25.0)
            shard.velocity *= 25.0 / shard.velocity.Length();
        if (shard.mass <= 1.0e-8)
            shard.mass = 1.0e-3;
        for (int k = 0; k < 4; ++k) {
            const auto& node = _xyzNodes[static_cast<std::size_t>(clones[k])];
            shard.local[static_cast<std::size_t>(k)] = node->GetPos() - shard.com;
            shard.localRest[static_cast<std::size_t>(k)] = shard.local[static_cast<std::size_t>(k)];
        }
        const double au = (shard.local[1] - shard.local[0]).Length();
        const double av = (shard.local[2] - shard.local[0]).Length();
        const double aw = (shard.local[3] - shard.local[0]).Length();
        const double floorI = shard.mass * 1.0e-6;
        shard.inertiaU = std::max(shard.mass * (av * av + aw * aw) / 12.0, floorI);
        shard.inertiaV = std::max(shard.mass * (au * au + aw * aw) / 12.0, floorI);
        shard.inertiaN = std::max(shard.mass * (au * au + av * av) / 12.0, floorI);
        if (_rigidNode.size() < _xyzNodes.size())
            _rigidNode.resize(_xyzNodes.size(), 0);
        for (int index : clones) {
            auto& node = _xyzNodes[static_cast<std::size_t>(index)];
            node->SetFixed(true);
            node->SetX0(node->GetPos());
            _rigidNode[static_cast<std::size_t>(index)] = 1;
        }
        _shards.push_back(shard);
        changed = true;
    }
    if (changed) {
        if (_shards.size() > shardsBefore)
            _torn = true;
        _contactArea.clear();
        _contactAreaNodes = 0;
        _contactAreaFaces = 0;
        std::vector<char> used(_xyzNodes.size(), 0);
        for (std::size_t i = 0; i < _tets.size() && i < _tetElements.size(); ++i) {
            if (!_tetElements[i] || (i < _tetLive.size() && !_tetLive[i]))
                continue;
            for (int node : _tets[i])
                if (node >= 0 && static_cast<std::size_t>(node) < used.size())
                    used[static_cast<std::size_t>(node)] = 1;
        }
        for (std::size_t q = 0; q < _quads.size() && q < _elements.size(); ++q) {
            if (!_elements[q])
                continue;
            if (q < _rigidQuad.size() && _rigidQuad[q])
                continue;
            for (int node : _quads[q]) {
                if (node >= 0 && static_cast<std::size_t>(node) < used.size())
                    used[static_cast<std::size_t>(node)] = 1;
            }
        }
        if (_rigidNode.size() < _xyzNodes.size())
            _rigidNode.resize(_xyzNodes.size(), 0);
        for (std::size_t i = 0; i < _xyzNodes.size(); ++i) {
            if (used[i] || !_xyzNodes[i] || _xyzNodes[i]->IsFixed())
                continue;
            _xyzNodes[i]->SetFixed(true);
            _xyzNodes[i]->SetPosDt(chrono::ChVector3d(0, 0, 0));
            _rigidNode[i] = 1;
        }
    }
    return changed;
}

void Part::buildSurface(const PartDefinition& definition, bool append) {
    using namespace chrono;
    using namespace chrono::fea;

    const SurfaceMesh& surface = definition.surface;
    std::vector<int> map(surface.vertices.size(), -1);
    const int volumeNodes = append ? static_cast<int>(_xyzNodes.size()) : 0;
    if (!append) {
        _nodes.clear();
        _elements.clear();
        _quads.clear();
        _triangles.clear();
        _boundaryNodes.clear();
        _xyzNodes.clear();
        _nodeMass.clear();
    }
    for (std::size_t v = 0; v < surface.vertices.size(); ++v) {
        const ChVector3d& vertex = surface.vertices[v];
        if (append) {
            double best = 1.0e-8;
            for (std::size_t n = 0; n < _xyzNodes.size(); ++n) {
                const double d = (_xyzNodes[n]->GetPos() - vertex).Length2();
                if (d < best) {
                    best = d;
                    map[v] = static_cast<int>(n);
                }
            }
            if (map[v] >= 0)
                continue;
        }
        auto node = chrono_types::make_shared<ChNodeFEAxyz>(vertex);
        node->SetMass(0);
        node->SetFixed(_fixed);
        node->SetX0(vertex);
        _mesh->AddNode(node);
        _xyzNodes.push_back(node);
        _nodeMass.push_back(0.0);
        map[v] = static_cast<int>(_xyzNodes.size()) - 1;
    }

    struct EdgeUse {
        int opposite = -1;
    };
    std::map<std::pair<int, int>, EdgeUse> edges;
    std::vector<std::array<int, 3>> faces;
    std::vector<double> faceThickness;
    faces.reserve(surface.triangles.size());
    for (std::size_t s = 0; s < surface.triangles.size(); ++s) {
        const auto& face = surface.triangles[s];
        if (face[0] == face[1] || face[1] == face[2] || face[2] == face[0])
            continue;
        if (face[0] < 0 || face[1] < 0 || face[2] < 0)
            continue;
        const int n = static_cast<int>(surface.vertices.size());
        if (face[0] >= n || face[1] >= n || face[2] >= n)
            continue;
        const int index = static_cast<int>(faces.size());
        std::array<int, 3> mapped = {map[static_cast<std::size_t>(face[0])],
                                     map[static_cast<std::size_t>(face[1])],
                                     map[static_cast<std::size_t>(face[2])]};
        if (mapped[0] < 0 || mapped[1] < 0 || mapped[2] < 0)
            continue;
        faces.push_back(mapped);
        faceThickness.push_back(s < surface.thickness.size() ? surface.thickness[s] : _thickness);
        (void)face;
        edges[{mapped[0], mapped[1]}] = {mapped[2]};
        edges[{mapped[1], mapped[2]}] = {mapped[0]};
        edges[{mapped[2], mapped[0]}] = {mapped[1]};
        (void)index;
    }

    const bool plastic = _material.allowPlasticity && _material.yieldStress > 0.0;

    auto neighbour = [&](int from, int to) -> std::shared_ptr<ChNodeFEAxyz> {
        const auto found = edges.find({from, to});
        if (found == edges.end() || found->second.opposite < 0)
            return nullptr;
        return _xyzNodes[static_cast<std::size_t>(found->second.opposite)];
    };

    for (const auto& face : faces) {
        const int i0 = face[0];
        const int i1 = face[1];
        const int i2 = face[2];
        if (i0 == i1 || i1 == i2 || i2 == i0)
            continue;
        const ChVector3d& a = _xyzNodes[static_cast<std::size_t>(i0)]->GetPos();
        const ChVector3d& b = _xyzNodes[static_cast<std::size_t>(i1)]->GetPos();
        const ChVector3d& c = _xyzNodes[static_cast<std::size_t>(i2)]->GetPos();
        const double area = 0.5 * (b - a).Cross(c - a).Length();
        if (area < 1.0e-12)
            continue;
        auto elasticity = chrono_types::make_shared<ThermalKirchhoff>(
            _material.youngsModulus, _material.poissonRatio);
        if (plastic)
            elasticity->setYield(_material.yieldStress);
        auto shell = chrono_types::make_shared<ChMaterialShellKirchhoff>(elasticity);
        shell->SetDensity(_material.density);
        const std::size_t sourceFace = static_cast<std::size_t>(&face - faces.data());
        const double thick = sourceFace < faceThickness.size() ? faceThickness[sourceFace] : _thickness;
        auto element = chrono_types::make_shared<ChElementShellBST>();
        element->SetNodes(_xyzNodes[i0], _xyzNodes[i1], _xyzNodes[i2], neighbour(i2, i1), neighbour(i0, i2),
                          neighbour(i1, i0));
        element->AddLayer(std::max(thick, 1.0e-5), 0.0, shell);
        _mesh->AddElement(element);
        _trianglesElements.push_back(element);
        _triangleThermal.push_back(elasticity);
        _triangles.push_back(face);

        const double share = _material.density * std::max(thick, 1.0e-5) * area / 3.0;
        _nodeMass[i0] += share;
        _nodeMass[i1] += share;
        _nodeMass[i2] += share;
    }

    for (const auto& edge : edges) {
        if (edges.find({edge.first.second, edge.first.first}) != edges.end())
            continue;
        _boundaryNodes.push_back(edge.first.first);
        _boundaryNodes.push_back(edge.first.second);
    }
    std::sort(_boundaryNodes.begin(), _boundaryNodes.end());
    _boundaryNodes.erase(std::unique(_boundaryNodes.begin(), _boundaryNodes.end()), _boundaryNodes.end());
    if (append && volumeNodes > 0) {
        _shellSolidTies.clear();
        const double reach = std::max(_thickness, 1.0e-4);
        const double reach2 = reach * reach;
        for (int index : _boundaryNodes) {
            if (index < volumeNodes || index < 0 || static_cast<std::size_t>(index) >= _xyzNodes.size())
                continue;
            const ChVector3d p = _xyzNodes[static_cast<std::size_t>(index)]->GetPos();
            int nearest = -1;
            double best = reach2;
            for (int solid = 0; solid < volumeNodes; ++solid) {
                if (solid < 0 || static_cast<std::size_t>(solid) >= _xyzNodes.size() || !_xyzNodes[static_cast<std::size_t>(solid)])
                    continue;
                const double d = (_xyzNodes[static_cast<std::size_t>(solid)]->GetPos() - p).Length2();
                if (d < best) {
                    best = d;
                    nearest = solid;
                }
            }
            if (nearest >= 0)
                _shellSolidTies.push_back({index, nearest});
        }
    }
    if (!append) {
        _loads.assign(nodeCount(), ChVector3d(0, 0, 0));
        seedTemperatures();
    }
}

}  // namespace rcs
