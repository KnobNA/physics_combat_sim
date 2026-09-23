#include "physics_core/part.hpp"
#include "physics_core/reissner_j2.hpp"

#include "chrono/core/ChMatrix33.h"
#include "chrono/fea/ChMaterialShellKirchhoff.h"
#include "chrono/fea/ChMaterialShellReissner.h"

#include <algorithm>
#include <cmath>
#include <limits>
#include <map>

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

double vonMisesPlaneStress(double sxx, double syy, double sxy) {
    const double v = sxx * sxx - sxx * syy + syy * syy + 3.0 * sxy * sxy;
    return v > 0.0 ? std::sqrt(v) : 0.0;
}

}  // namespace

void Part::build(const PartDefinition& definition) {
    using namespace chrono;
    using namespace chrono::fea;

    _name = definition.name;
    _material = definition.material;
    _thickness = definition.thickness;
    _temperature = _material.referenceTemperature;
    _fixed = definition.fixed;
    _torn = false;
    _shattering = false;
    _elasticity.reset();
    _shellMaterial.reset();
    _tearTwins.clear();
    _xyzNodes.clear();
    _trianglesElements.clear();

    _mesh = chrono_types::make_shared<ChMesh>();
    // Gravity is one of several nodal loads we assemble ourselves each step,
    // so the mesh must not add its own on top.
    _mesh->SetAutomaticGravity(false);

    if (!definition.surface.empty()) {
        buildSurface(definition);
        return;
    }

    _elasticity = chrono_types::make_shared<ThermalReissnerElasticity>(
        _material.youngsModulus, _material.poissonRatio);
    _elasticity->SetThermalStrain(thermalStrain());
    auto damping = chrono_types::make_shared<ChDampingReissnerRayleigh>(
        _elasticity, _material.dampingBeta);
    std::shared_ptr<chrono::fea::ChPlasticityReissner> plasticity;
    if (_material.allowPlasticity && _material.yieldStress > 0.0)
        plasticity = chrono_types::make_shared<ReissnerJ2>(_material.yieldStress, _material.youngsModulus);
    _shellMaterial = chrono_types::make_shared<ChMaterialShellReissner>(
        _elasticity, plasticity, damping);
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
            element->AddLayer(_thickness, 0, _shellMaterial);
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
        for (int index : _boundaryNodes) {
            if (index >= 0 && static_cast<std::size_t>(index) < _nodes.size())
                _nodes[static_cast<std::size_t>(index)]->SetFixed(true);
        }
    }

    _loads.assign(_nodes.size(), ChVector3d(0, 0, 0));
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

void Part::commitLoads() {
    if (!_xyzNodes.empty()) {
        for (std::size_t i = 0; i < _xyzNodes.size() && i < _loads.size(); ++i)
            _xyzNodes[i]->SetForce(_loads[i]);
        return;
    }
    for (std::size_t i = 0; i < _nodes.size(); ++i)
        _nodes[i]->SetForce(_loads[i]);
}

double Part::thermalStrain() const {
    return _material.thermalExpansion * (_temperature - _material.referenceTemperature);
}

void Part::setTemperature(double kelvin) {
    _temperature = kelvin;
    if (_elasticity)
        _elasticity->SetThermalStrain(thermalStrain());
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
    // registers as zero stress. Plastic strain and the uniform thermal strain
    // are removed before the elastic law.
    const double eth = thermalStrain();
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
        out[q] = static_cast<float>(fiber.vonMises);
    }

    if (!_quads.empty() || _triangles.empty())
        return;

    out.assign(_triangles.size(), 0.0f);
    for (std::size_t t = 0; t < _triangles.size(); ++t) {
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
        const double eth = thermalStrain();
        const double exx = 0.5 * (f00 * f00 + f10 * f10 - 1.0) - eth;
        const double eyy = 0.5 * (f01 * f01 + f11 * f11 - 1.0) - eth;
        const double exy = 0.5 * (f00 * f01 + f10 * f11);
        const double mxx = cm * (exx + nu * eyy);
        const double myy = cm * (eyy + nu * exx);
        const double mxy = shear * exy;
        const double peak = vonMisesPlaneStress(mxx, myy, mxy);
        out[t] = std::isfinite(peak) ? static_cast<float>(peak) : 0.0f;
    }
}

double Part::peakTension() const {
    double peak = 0.0;
    for (std::size_t q = 0; q < _quads.size(); ++q) {
        if (q < _rigidQuad.size() && _rigidQuad[q])
            continue;
        peak = std::max(peak, quadMechanics(q).tension);
    }
    return peak;
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
        stretch = std::max(stretch, now / rest);
    }
    return stretch;
}

bool Part::nodeRigid(std::size_t index) const {
    return index < _rigidNode.size() && _rigidNode[index] != 0;
}

void Part::clearLinkedNodes() {
    std::fill(_linked.begin(), _linked.end(), 0);
}

void Part::markLinkedNode(int index) {
    if (index < 0)
        return;
    if (_linked.size() < _nodes.size())
        _linked.resize(_nodes.size(), 0);
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
    const bool boundaryA =
        std::find(_boundaryNodes.begin(), _boundaryNodes.end(), a) != _boundaryNodes.end();
    const bool boundaryB =
        std::find(_boundaryNodes.begin(), _boundaryNodes.end(), b) != _boundaryNodes.end();
    if (boundaryA && boundaryB)
        _boundaryNodes.push_back(index);
    return index;
}

void Part::detachQuad(std::size_t quad) {
    if (quad >= _quads.size() || quad >= _elements.size())
        return;
    int cloned[4];
    for (int i = 0; i < 4; ++i)
        cloned[i] = cloneNode(_quads[quad][static_cast<std::size_t>(i)]);
    _elements[quad]->SetNodes(_nodes[cloned[0]], _nodes[cloned[1]], _nodes[cloned[2]], _nodes[cloned[3]]);
    _mesh->InitializeElement(_elements[quad]);
    _quads[quad] = {cloned[0], cloned[1], cloned[2], cloned[3]};
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
    element->AddLayer(_thickness, 0.0, _shellMaterial);
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

void Part::finishShards(bool floorEnabled, double floorHeight) {
    using namespace chrono;
    const double plane = floorHeight + 0.5 * _thickness;
    const ChVector3d up{0, 1, 0};
    constexpr double kFriction = 0.55;
    for (auto& shard : _shards) {
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
            // The floor pushes at the corner that is actually touching, so a
            // chip on its tip topples instead of being glued there.
            for (int pass = 0; pass < 4; ++pass) {
                for (int k = 0; k < 4; ++k) {
                    const ChVector3d point = cornerAt(k);
                    const double pen = plane - point.y();
                    if (pen <= 1.0e-4)
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
                    const double bias = std::min(std::max(0.0, pen - 1.0e-4) * 30.0, 1.0);
                    const double j = -(vn + bias) / eff;
                    if (j <= 0.0)
                        continue;
                    const double limited = std::min(j, shard.mass * 20.0);
                    applyImpulse(up * limited, arm);
                    ChVector3d tangent = shard.velocity + shard.angularVelocity.Cross(arm);
                    tangent.y() = 0.0;
                    const double slide = tangent.Length();
                    if (slide > 1.0e-6) {
                        const double jf = std::min(kFriction * limited, slide * shard.mass);
                        applyImpulse(tangent * (-jf / slide), arm);
                    }
                }
            }

            int support = 0;
            double lowest = std::numeric_limits<double>::max();
            ChVector3d normal = shard.rotation.Rotate(
                (shard.local[1] - shard.local[0]).Cross(shard.local[3] - shard.local[0]));
            if (normal.Length() > 1.0e-12)
                normal.Normalize();
            else
                normal = up;
            if (normal.y() < 0.0)
                normal = -normal;
            for (int k = 0; k < 4; ++k) {
                const double y = cornerAt(k).y();
                lowest = std::min(lowest, y);
                if (y <= plane + 1.5e-3)
                    ++support;
            }
            // Three corners on the plane means the face is down. A tip or an
            // edge is not a resting pose for a thin piece of glass.
            const bool faceDown = support >= 3 && normal.y() > 0.85;
            if (faceDown && shard.velocity.Length() < 0.08 && shard.angularVelocity.Length() < 1.0) {
                shard.velocity = {0, 0, 0};
                shard.angularVelocity = {0, 0, 0};
                shard.asleep = true;
                if (lowest < plane)
                    shard.com.y() += plane - lowest;
            } else if (support > 0 && normal.y() < 0.75 && shard.angularVelocity.Length() < 0.4) {
                const ChVector3d axis = up.Cross(normal);
                const double len = axis.Length();
                if (len > 1.0e-4)
                    shard.angularVelocity += axis * (0.6 / len);
            } else if (lowest < plane - 0.02) {
                shard.com.y() += (plane - 0.002) - lowest;
            }
        }
        for (int k = 0; k < 4; ++k) {
            const std::size_t index = static_cast<std::size_t>(shard.nodes[static_cast<std::size_t>(k)]);
            if (index >= _nodes.size())
                continue;
            const ChVector3d arm = shard.rotation.Rotate(shard.local[static_cast<std::size_t>(k)]);
            auto& node = _nodes[index];
            node->SetPos(shard.com + arm);
            node->SetRot(shard.rotation);
            node->SetPosDt(shard.velocity + shard.angularVelocity.Cross(arm));
            node->SetAngVelParent(shard.angularVelocity);
            node->SetX0(node->Frame());
        }
    }
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

void Part::promoteLooseFragments() {
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
    constexpr int kMaxShards = 64;
    const double minEdge = _material.minimumShardSize > 1.0e-4 ? _material.minimumShardSize : 0.02;

    for (const auto& group : groups) {
        std::vector<int> quads;
        bool anchored = false;
        for (int quadIndex : group.second) {
            if (quadIndex < 0 || _rigidQuad[static_cast<std::size_t>(quadIndex)])
                continue;
            quads.push_back(quadIndex);
            const auto& quad = _quads[static_cast<std::size_t>(quadIndex)];
            for (int k = 0; k < 4; ++k) {
                const int index = quad[static_cast<std::size_t>(k)];
                if (index < 0 || static_cast<std::size_t>(index) >= _nodes.size()) {
                    anchored = true;
                    continue;
                }
                const bool linked =
                    static_cast<std::size_t>(index) < _linked.size() && _linked[static_cast<std::size_t>(index)];
                anchored = anchored || _nodes[static_cast<std::size_t>(index)]->IsFixed() || linked;
            }
        }
        // A clamped corner or a live weld is still part of the plate.
        if (quads.empty() || anchored)
            continue;

        bool splitting = true;
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
    if (_fixed || _elements.empty() || _material.fractureStrength <= 0.0)
        return false;

    const double minEdge = _material.minimumShardSize > 1.0e-4 ? _material.minimumShardSize : 0.02;
    const int budget = _material.shatterToDust ? 24 : 4;
    const double strength = _material.fractureStrength;
    const double toughness = _material.fractureToughness;
    const double gc = (_material.youngsModulus > 0.0 && toughness > 0.0)
                          ? toughness * toughness / _material.youngsModulus
                          : 0.0;

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
        const double longest = std::max(fiber.edgeU, fiber.edgeV);
        const double crackLength = std::max(minEdge, std::min(fiber.edgeU, fiber.edgeV));
        const double crackArea = crackLength * _thickness;
        const bool energyOk = gc <= 0.0 || fiber.energy >= gc * crackArea;
        const bool over = fiber.tension >= strength && energyOk;
        const bool dust = _material.shatterToDust && _shattering && longest > minEdge;
        if (!over && !dust)
            continue;
        candidates.push_back({q, fiber.tension, fiber.theta, fiber.edgeU, fiber.edgeV, over});
    }
    std::sort(candidates.begin(), candidates.end(),
              [](const Candidate& a, const Candidate& b) { return a.tension > b.tension; });

    int live = 0;
    for (std::size_t q = 0; q < _elements.size() && q < _quads.size(); ++q) {
        if (_elements[q] && (q >= _rigidQuad.size() || !_rigidQuad[q]))
            ++live;
    }

    bool changed = false;
    int made = 0;
    for (const Candidate& candidate : candidates) {
        if (made >= budget)
            break;
        if (candidate.quad >= _quads.size())
            continue;
        const double longest = std::max(candidate.edgeU, candidate.edgeV);
        const double shortest = std::min(candidate.edgeU, candidate.edgeV);
        const double aspect = shortest > 1.0e-9 ? longest / shortest : 1.0;
        if (longest < minEdge) {
            if (!candidate.over)
                continue;
            detachQuad(candidate.quad);
        } else if (live >= 48) {
            // Further cracks leave the shell instead of refining it.
            detachQuad(candidate.quad);
        } else {
            // A needle is cut across its long side. Rankine still chooses
            // the direction while the quad is roughly square.
            const bool crossCut = aspect > 2.0 && (_shattering || candidate.over);
            const bool alongU = (crossCut || (!candidate.over && _material.shatterToDust))
                                    ? candidate.edgeV >= candidate.edgeU
                                    : std::abs(std::cos(candidate.theta)) >= std::abs(std::sin(candidate.theta));
            if (!splitQuad(candidate.quad, alongU))
                continue;
            ++live;
        }
        changed = true;
        _shattering = true;
        ++made;
    }
    if (changed) {
        promoteLooseFragments();
        rebuildTrianglesFromQuads();
    }
    return changed;
}

bool Part::tearFailed() {
    bool changed = false;
    if (!_elements.empty()) {
        if (_material.fracture == FractureMode::Brittle)
            return shatterQuads() ? (_torn = true) : false;

        std::vector<std::size_t> failed;
        for (std::size_t q = 0; q < _quads.size() && q < _elements.size(); ++q) {
            const double stretch = maxEdgeStretch(_quads[q].data(), 4);
            ReissnerPlasticData plastic;
            const double accumulated = averagePlastic(q, plastic) ? plastic.p_strain_acc : 0.0;
            const bool stretchFail =
                _material.fracture != FractureMode::Brittle && stretch >= _material.tearStretch;
            const bool plasticFail = _material.fracture == FractureMode::Ductile &&
                                     _material.failPlasticStrain > 0.0 &&
                                     accumulated >= _material.failPlasticStrain;
            if (stretchFail || plasticFail)
                failed.push_back(q);
        }
        for (std::size_t q : failed) {
            detachQuad(q);
            changed = true;
        }
        if (changed)
            rebuildTrianglesFromQuads();
    } else if (!_trianglesElements.empty()) {
        std::vector<float> stress;
        if (_material.fracture == FractureMode::Brittle)
            evaluateStress(stress);
        std::vector<std::size_t> failed;
        for (std::size_t t = 0; t < _triangles.size() && t < _trianglesElements.size(); ++t) {
            const bool stretchFail = _material.fracture != FractureMode::Brittle &&
                                     maxEdgeStretch(_triangles[t].data(), 3) >= _material.tearStretch;
            const bool brittleFail = _material.fracture == FractureMode::Brittle &&
                                     t < stress.size() &&
                                     stress[t] >= static_cast<float>(_material.fractureStrength);
            if (stretchFail || brittleFail)
                failed.push_back(t);
        }
        for (std::size_t t : failed) {
            int cloned[3];
            for (int i = 0; i < 3; ++i)
                cloned[i] = cloneNode(_triangles[t][static_cast<std::size_t>(i)]);
            _trianglesElements[t]->SetNodes(_xyzNodes[cloned[0]], _xyzNodes[cloned[1]], _xyzNodes[cloned[2]],
                                            nullptr, nullptr, nullptr);
            _mesh->InitializeElement(_trianglesElements[t]);
            _triangles[t] = {cloned[0], cloned[1], cloned[2]};
            changed = true;
        }
    }
    if (changed)
        _torn = true;
    return changed;
}

void Part::buildSurface(const PartDefinition& definition) {
    using namespace chrono;
    using namespace chrono::fea;

    const SurfaceMesh& surface = definition.surface;
    _nodes.clear();
    _elements.clear();
    _quads.clear();
    _triangles.clear();
    _boundaryNodes.clear();
    _nodeMass.assign(surface.vertices.size(), 0.0);

    for (const ChVector3d& vertex : surface.vertices) {
        auto node = chrono_types::make_shared<ChNodeFEAxyz>(vertex);
        node->SetMass(0);
        node->SetFixed(_fixed);
        node->SetX0(vertex);
        _mesh->AddNode(node);
        _xyzNodes.push_back(node);
    }

    struct EdgeUse {
        int opposite = -1;
    };
    std::map<std::pair<int, int>, EdgeUse> edges;
    std::vector<std::array<int, 3>> faces;
    faces.reserve(surface.triangles.size());
    for (const auto& face : surface.triangles) {
        if (face[0] == face[1] || face[1] == face[2] || face[2] == face[0])
            continue;
        if (face[0] < 0 || face[1] < 0 || face[2] < 0)
            continue;
        const int n = static_cast<int>(surface.vertices.size());
        if (face[0] >= n || face[1] >= n || face[2] >= n)
            continue;
        const int index = static_cast<int>(faces.size());
        faces.push_back(face);
        edges[{face[0], face[1]}] = {face[2]};
        edges[{face[1], face[2]}] = {face[0]};
        edges[{face[2], face[0]}] = {face[1]};
        (void)index;
    }

    auto elasticity = chrono_types::make_shared<ChElasticityKirchhoffIsothropic>(
        _material.youngsModulus, _material.poissonRatio);
    auto shell = chrono_types::make_shared<ChMaterialShellKirchhoff>(elasticity);
    shell->SetDensity(_material.density);

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
        auto element = chrono_types::make_shared<ChElementShellBST>();
        element->SetNodes(_xyzNodes[i0], _xyzNodes[i1], _xyzNodes[i2], neighbour(i2, i1), neighbour(i0, i2),
                          neighbour(i1, i0));
        element->AddLayer(_thickness, 0.0, shell);
        _mesh->AddElement(element);
        _trianglesElements.push_back(element);
        _triangles.push_back(face);

        const ChVector3d& a = surface.vertices[i0];
        const ChVector3d& b = surface.vertices[i1];
        const ChVector3d& c = surface.vertices[i2];
        const double area = 0.5 * (b - a).Cross(c - a).Length();
        const double share = _material.density * _thickness * area / 3.0;
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
    _loads.assign(nodeCount(), ChVector3d(0, 0, 0));
}

}  // namespace rcs
