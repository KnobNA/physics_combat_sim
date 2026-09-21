#include "physics_core/part.hpp"

#include "chrono/core/ChMatrix33.h"
#include "chrono/fea/ChMaterialShellReissner.h"

#include <algorithm>
#include <cmath>
#include <limits>

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
    _fixed = definition.fixed;

    _mesh = chrono_types::make_shared<ChMesh>();
    // Gravity is one of several nodal loads we assemble ourselves each step,
    // so the mesh must not add its own on top.
    _mesh->SetAutomaticGravity(false);

    auto elasticity = chrono_types::make_shared<ChElasticityReissnerIsothropic>(
        _material.youngsModulus, _material.poissonRatio, 1.0, 0.1);
    auto damping = chrono_types::make_shared<ChDampingReissnerRayleigh>(
        elasticity, _material.dampingBeta);
    auto shellMaterial = chrono_types::make_shared<ChMaterialShellReissner>(
        elasticity, nullptr, damping);
    shellMaterial->SetDensity(_material.density);

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
            element->AddLayer(_thickness, 0, shellMaterial);
            _mesh->AddElement(element);
            _elements.push_back(element);

            _quads.push_back({i00, i10, i11, i01});
            _triangles.push_back({i00, i10, i11});
            _triangles.push_back({i00, i11, i01});

            for (int index : {i00, i10, i11, i01})
                _nodeMass[index] += quadMass * 0.25;
        }
    }

    _loads.assign(_nodes.size(), ChVector3d(0, 0, 0));
}

void Part::transform(const chrono::ChVector3d& translation,
                     const chrono::ChQuaterniond& rotation,
                     const chrono::ChVector3d& pivot,
                     bool alsoRest)
{
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
    for (std::size_t i = 0; i < _nodes.size(); ++i)
        _nodes[i]->SetForce(_loads[i]);
}

void Part::evaluateStress(std::vector<float>& out) const {
    using namespace chrono;

    out.assign(_quads.size(), 0.0f);

    const double E = _material.youngsModulus;
    const double nu = _material.poissonRatio;
    const double cm = E / (1.0 - nu * nu);
    const double shear = E / (1.0 + nu);
    const double zb = 0.5 * _thickness;

    for (std::size_t q = 0; q < _quads.size(); ++q) {
        const auto& idx = _quads[q];
        ChVector3d rest[4], now[4];
        ChQuaterniond restRot[4], nowRot[4];
        for (int i = 0; i < 4; ++i) {
            const auto& node = _nodes[idx[i]];
            rest[i] = node->GetX0().GetPos();
            restRot[i] = node->GetX0().GetRot();
            now[i] = node->GetPos();
            nowRot[i] = node->GetRot();
        }

        ChVector3d e1r, e2r, nr, e1c, e2c, nc;
        if (!quadFrame(rest, e1r, e2r, nr) || !quadFrame(now, e1c, e2c, nc))
            continue;

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
            continue;

        // In-plane Green-Lagrange strain, so a large rigid rotation of the part
        // registers as zero stress.
        const double exx = 0.5 * (F[0][0] * F[0][0] + F[1][0] * F[1][0] - 1.0);
        const double eyy = 0.5 * (F[0][1] * F[0][1] + F[1][1] * F[1][1] - 1.0);
        const double exy = 0.5 * (F[0][0] * F[0][1] + F[1][0] * F[1][1]);

        const double mxx = cm * (exx + nu * eyy);
        const double myy = cm * (eyy + nu * exx);
        const double mxy = shear * exy;

        // Curvature from the shell director field. Each node's rotation is
        // taken relative to node 0 so that tumbling never wraps the rotation
        // vector and shows up as phantom bending.
        const ChQuaterniond refNow = nowRot[0].GetConjugate();
        const ChQuaterniond refRest = restRot[0].GetConjugate();
        double w[4][2];
        for (int i = 0; i < 4; ++i) {
            const ChQuaterniond relative =
                (refNow * nowRot[i]) * (refRest * restRot[i]).GetConjugate();
            const ChVector3d omega = nowRot[0].Rotate(relative.GetRotVec());
            w[i][0] = omega.Dot(e1c);
            w[i][1] = omega.Dot(e2c);
        }

        double G[2][2];
        double bxx = 0.0, byy = 0.0, bxy = 0.0;
        if (fitGradient(a, w, G)) {
            const double kxx = G[1][0];
            const double kyy = -G[0][1];
            const double kxy = 0.5 * (G[1][1] - G[0][0]);
            bxx = cm * zb * (kxx + nu * kyy);
            byy = cm * zb * (kyy + nu * kxx);
            bxy = shear * zb * kxy;
        }

        const double top = vonMisesPlaneStress(mxx + bxx, myy + byy, mxy + bxy);
        const double bottom = vonMisesPlaneStress(mxx - bxx, myy - byy, mxy - bxy);
        const double peak = std::max(top, bottom);
        out[q] = std::isfinite(peak) ? static_cast<float>(peak) : 0.0f;
    }
}

chrono::ChVector3d Part::center() const {
    chrono::ChVector3d sum{0, 0, 0};
    if (_nodes.empty())
        return sum;
    for (const auto& node : _nodes)
        sum += node->GetPos();
    return sum / static_cast<double>(_nodes.size());
}

void Part::bounds(chrono::ChVector3d& lo, chrono::ChVector3d& hi, double pad) const {
    const double big = std::numeric_limits<double>::max();
    lo = {big, big, big};
    hi = {-big, -big, -big};
    for (const auto& node : _nodes) {
        const auto& p = node->GetPos();
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
    for (const auto& node : _nodes) {
        const auto& p = node->GetPos();
        if (!std::isfinite(p.x()) || !std::isfinite(p.y()) || !std::isfinite(p.z()))
            return false;
        if (p.Length() > 1.0e3)
            return false;
    }
    return true;
}

}  // namespace rcs
