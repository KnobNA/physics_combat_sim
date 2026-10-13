#include "physics_core/physics_units.hpp"

#include <algorithm>
#include <cmath>
#include <map>
#include <utility>

namespace rcs {
namespace {

using chrono::ChVector3d;
using chrono::ChVectorDynamic;

bool projectRankineArray(double stress[6], double tension, double compression) {
    if (tension <= 0.0 && compression <= 0.0)
        return false;

    double A[3][3] = {{stress[0], stress[3], stress[5]},
                      {stress[3], stress[1], stress[4]},
                      {stress[5], stress[4], stress[2]}};
    double V[3][3] = {{1.0, 0.0, 0.0}, {0.0, 1.0, 0.0}, {0.0, 0.0, 1.0}};
    for (int iter = 0; iter < 16; ++iter) {
        int p = 0;
        int q = 1;
        double maxOff = std::abs(A[0][1]);
        if (std::abs(A[0][2]) > maxOff) {
            maxOff = std::abs(A[0][2]);
            p = 0;
            q = 2;
        }
        if (std::abs(A[1][2]) > maxOff) {
            maxOff = std::abs(A[1][2]);
            p = 1;
            q = 2;
        }
        const double scale = 1.0 + std::abs(A[0][0]) + std::abs(A[1][1]) + std::abs(A[2][2]);
        if (maxOff < 1.0e-12 * scale)
            break;

        const double app = A[p][p];
        const double aqq = A[q][q];
        const double apq = A[p][q];
        const double tau = (aqq - app) / (2.0 * apq);
        const double t = (tau >= 0.0 ? 1.0 : -1.0) / (std::abs(tau) + std::sqrt(1.0 + tau * tau));
        const double c = 1.0 / std::sqrt(1.0 + t * t);
        const double s = t * c;

        const double appNew = c * c * app - 2.0 * s * c * apq + s * s * aqq;
        const double aqqNew = s * s * app + 2.0 * s * c * apq + c * c * aqq;
        A[p][p] = appNew;
        A[q][q] = aqqNew;
        A[p][q] = A[q][p] = 0.0;
        for (int k = 0; k < 3; ++k) {
            if (k == p || k == q)
                continue;
            const double apk = A[p][k];
            const double aqk = A[q][k];
            A[p][k] = A[k][p] = c * apk - s * aqk;
            A[q][k] = A[k][q] = s * apk + c * aqk;
        }
        for (int k = 0; k < 3; ++k) {
            const double vip = V[k][p];
            const double viq = V[k][q];
            V[k][p] = c * vip - s * viq;
            V[k][q] = s * vip + c * viq;
        }
    }

    const double lo = compression > 0.0 ? -compression : -1.0e30;
    const double hi = tension > 0.0 ? tension : 1.0e30;
    double lambda[3] = {A[0][0], A[1][1], A[2][2]};
    bool clipped = false;
    for (double& value : lambda) {
        const double kept = std::clamp(value, lo, hi);
        clipped = clipped || kept != value;
        value = kept;
    }
    if (!clipped)
        return false;

    double R[3][3] = {};
    for (int i = 0; i < 3; ++i)
        for (int j = 0; j < 3; ++j)
            for (int k = 0; k < 3; ++k)
                R[i][j] += V[i][k] * lambda[k] * V[j][k];
    stress[0] = R[0][0];
    stress[1] = R[1][1];
    stress[2] = R[2][2];
    stress[3] = 0.5 * (R[0][1] + R[1][0]);
    stress[4] = 0.5 * (R[1][2] + R[2][1]);
    stress[5] = 0.5 * (R[0][2] + R[2][0]);
    return true;
}

}  // namespace

double vonMisesStress(double sxx, double syy, double szz, double sxy, double syz, double szx) {
    const double mean = (sxx + syy + szz) / 3.0;
    const double d0 = sxx - mean;
    const double d1 = syy - mean;
    const double d2 = szz - mean;
    return std::sqrt(1.5 * (d0 * d0 + d1 * d1 + d2 * d2 + 2.0 * (sxy * sxy + syz * syz + szx * szx)));
}

double vonMisesStress(const ChVectorDynamic<>& stress) {
    return vonMisesStress(stress(0), stress(1), stress(2), stress(3), stress(4), stress(5));
}

double vonMisesPlaneStress(double sxx, double syy, double sxy) {
    const double v = sxx * sxx - sxx * syy + syy * syy + 3.0 * sxy * sxy;
    return v > 0.0 ? std::sqrt(v) : 0.0;
}

double equivalentPlastic(const double* strain) {
    double dev = strain[0] * strain[0] + strain[1] * strain[1] + strain[2] * strain[2];
    dev += 0.5 * (strain[3] * strain[3] + strain[4] * strain[4] + strain[5] * strain[5]);
    return std::sqrt(std::max(0.0, 2.0 / 3.0 * dev));
}

double flowStress(double yield, double ultimate, double failPlastic, double equivalentPlastic,
                  double strengthScale) {
    const double scale = std::max(0.0, strengthScale);
    const double sy = std::max(0.0, yield) * scale;
    if (sy <= 0.0)
        return 0.0;
    const double su = std::max(ultimate, yield) * scale;
    constexpr double kOffset = 0.002;
    const double A = std::max(failPlastic, kOffset);
    const double ep = std::max(0.0, equivalentPlastic);
    const double s0 = 0.90 * sy;
    if (ep <= kOffset)
        return s0 + (sy - s0) * (ep / kOffset);
    // Yield plateau, then rise to UTS at A. A flat range after the 0.2%
    // offset is the handbook metal: holding UTS is past uniform elongation
    // and the bar must flow, not sit at a few millistrain on a steep H.
    const double plateau = 0.70 * A;
    if (ep <= plateau)
        return sy;
    const double span = std::max(A - plateau, 1.0e-12);
    const double u = std::min(1.0, (ep - plateau) / span);
    return sy + (su - sy) * u;
}

double ductileTearStretch(double failPlastic, double ductility) {
    return 1.0 + 2.0 * std::max(0.0, failPlastic) * std::max(0.0, ductility);
}

bool projectRankine(double stress[6], double tension, double compression) {
    return projectRankineArray(stress, tension, compression);
}

bool projectRankine(ChVectorDynamic<>& stress, double tension, double compression) {
    double values[6] = {stress(0), stress(1), stress(2), stress(3), stress(4), stress(5)};
    if (!projectRankineArray(values, tension, compression))
        return false;
    stress(0) = values[0];
    stress(1) = values[1];
    stress(2) = values[2];
    stress(3) = values[3];
    stress(4) = values[4];
    stress(5) = values[5];
    return true;
}

void mooneyRivlinCauchy(const chrono::ChMatrix33<>& stretch, double youngs, double nu,
                        double thermalStrain, ChVectorDynamic<>& stress) {
    const double E = std::max(youngs, 1.0);
    const double n = std::clamp(nu, 0.0, 0.49);
    const double mu = E / (2.0 * (1.0 + n));
    const double bulk = E / (3.0 * std::max(1.0 - 2.0 * n, 0.02));
    const double c10 = (8.0 / 9.0) * (0.5 * mu);
    const double c01 = (1.0 / 9.0) * (0.5 * mu);
    const double lamTh = std::max(1.0 + thermalStrain, 1.0e-6);
    chrono::ChMatrix33<> U = stretch;
    U *= 1.0 / lamTh;
    double J = U.determinant();
    J = std::max(J, 1.0e-12);
    const chrono::ChMatrix33<> b = U * U;
    const double j23 = std::pow(J, -2.0 / 3.0);
    const chrono::ChMatrix33<> bBar = b * j23;
    const double i1 = bBar(0, 0) + bBar(1, 1) + bBar(2, 2);
    const chrono::ChMatrix33<> bBar2 = bBar * bBar;
    chrono::ChMatrix33<> tau = 2.0 * c10 * bBar + 2.0 * c01 * (i1 * bBar - bBar2);
    const double tr = (tau(0, 0) + tau(1, 1) + tau(2, 2)) / 3.0;
    tau(0, 0) -= tr;
    tau(1, 1) -= tr;
    tau(2, 2) -= tr;
    tau *= 1.0 / J;
    const double p = bulk * (J - 1.0);
    tau(0, 0) += p;
    tau(1, 1) += p;
    tau(2, 2) += p;
    stress.setZero(6);
    stress(0) = tau(0, 0);
    stress(1) = tau(1, 1);
    stress(2) = tau(2, 2);
    stress(3) = 0.5 * (tau(0, 1) + tau(1, 0));
    stress(4) = 0.5 * (tau(1, 2) + tau(2, 1));
    stress(5) = 0.5 * (tau(0, 2) + tau(2, 0));
}

bool volumeHyperelastic(const Material& material) {
    return material.fracture == FractureMode::Compliant && material.yieldStress <= 0.0 &&
           material.poissonRatio >= 0.45;
}

void henckyVoigt(const chrono::ChMatrix33<>& stretch, ChVectorDynamic<>& strain) {
    chrono::ChMatrix33<> evec;
    chrono::ChVectorN<double, 3> evals;
    stretch.SelfAdjointEigenSolve(evec, evals);
    chrono::ChMatrix33<> logU;
    logU.setZero();
    for (int i = 0; i < 3; ++i) {
        const double lam = std::max(evals(i), 1.0e-12);
        const double le = std::log(lam);
        const double vx = evec(0, i);
        const double vy = evec(1, i);
        const double vz = evec(2, i);
        logU(0, 0) += le * vx * vx;
        logU(0, 1) += le * vx * vy;
        logU(0, 2) += le * vx * vz;
        logU(1, 0) += le * vy * vx;
        logU(1, 1) += le * vy * vy;
        logU(1, 2) += le * vy * vz;
        logU(2, 0) += le * vz * vx;
        logU(2, 1) += le * vz * vy;
        logU(2, 2) += le * vz * vz;
    }
    strain.setZero(6);
    strain(0) = logU(0, 0);
    strain(1) = logU(1, 1);
    strain(2) = logU(2, 2);
    strain(3) = logU(0, 1) + logU(1, 0);
    strain(4) = logU(1, 2) + logU(2, 1);
    strain(5) = logU(0, 2) + logU(2, 0);
}

bool interfaceShouldOpen(const Material& material, double kelvin, double tension,
                         const chrono::ChVector3d& tensionAxis, double energy, double area) {
    if (material.interfaceStrength <= 0.0 || material.plies.size() < 2)
        return false;
    chrono::ChVector3d normal(material.stackNormal[0], material.stackNormal[1], material.stackNormal[2]);
    chrono::ChVector3d axis = tensionAxis;
    if (normal.Length2() < 0.25 || axis.Length2() < 0.25)
        return false;
    normal.Normalize();
    axis.Normalize();
    const double retention = std::max(thermalStrengthRetention(material, kelvin), 1.0e-4);
    const double drive = tension * std::abs(axis.Dot(normal));
    if (drive < material.interfaceStrength * retention)
        return false;
    if (material.interfaceToughness <= 0.0)
        return true;
    Material gate = material;
    gate.fractureToughness = material.interfaceToughness;
    return energyAllowsOpening(gate, energy, area);
}

bool brittleShouldOpen(const Material& material, double kelvin, double tension, double trialEnergy,
                       double tensileStrain, double faceArea, double compression, double mean,
                       const chrono::ChVector3d& tensionAxis) {
    (void)tensileStrain;
    (void)compression;
    (void)mean;
    if (material.fracture != FractureMode::Brittle || material.fractureStrength <= 0.0)
        return false;
    const double retention = thermalStrengthRetention(material, kelvin);
    double strength = material.fractureStrength * std::max(retention, 1.0e-4);
    strength *= grainOpeningScale(material, tensionAxis);
    const double residual =
        material.breakFabric == BreakFabric::Tempered ? std::max(0.0, material.residualTension) : 0.0;
    // Rankine is max principal. Tempered residual is extra tensile drive, not
    // a new contact law. Grain raises the opening bar along the fiber.
    // Toughness spends G_c on the new area so a finer mesh does not open first.
    if (tension + residual < strength)
        return false;
    return energyAllowsOpening(material, trialEnergy, faceArea);
}

bool brittleShouldCrush(const Material& material, double kelvin, double compression) {
    if (material.fracture != FractureMode::Brittle)
        return false;
    const double fc =
        rankineCompressionLimit(material) * std::max(thermalStrengthRetention(material, kelvin), 1.0e-4);
    return fc > 0.0 && compression >= fc;
}

int findRoot(std::vector<int>& set, int i) {
    while (set[i] != i) {
        set[i] = set[set[i]];
        i = set[i];
    }
    return i;
}

bool unite(std::vector<int>& set, int a, int b) {
    const int ra = findRoot(set, a);
    const int rb = findRoot(set, b);
    if (ra == rb)
        return false;
    set[rb] = ra;
    return true;
}

double distanceToBead(const std::vector<ChVector3d>& bead, const ChVector3d& p) {
    if (bead.size() < 2)
        return 0.0;
    double best = 1.0e30;
    const std::size_t n = bead.size();
    for (std::size_t i = 0; i < n; ++i) {
        const ChVector3d a = bead[i];
        const ChVector3d b = bead[(i + 1) % n];
        const ChVector3d ab = b - a;
        const double len2 = ab.Length2();
        if (len2 < 1.0e-24)
            continue;
        const double t = std::clamp((p - a).Dot(ab) / len2, 0.0, 1.0);
        best = std::min(best, (p - (a + ab * t)).Length());
    }
    return best;
}

double beadPerimeter(const std::vector<ChVector3d>& bead) {
    if (bead.size() < 2)
        return 0.0;
    double len = 0.0;
    const std::size_t n = bead.size();
    const std::size_t last = n >= 3 ? n : n - 1;
    for (std::size_t i = 0; i < last; ++i) {
        const ChVector3d a = bead[i];
        const ChVector3d b = bead[(i + 1) % n];
        len += (b - a).Length();
    }
    return len;
}

double closedBeadArea(const std::vector<ChVector3d>& bead) {
    if (bead.size() < 3)
        return 0.0;
    ChVector3d acc{0, 0, 0};
    const std::size_t n = bead.size();
    for (std::size_t i = 0; i < n; ++i)
        acc += bead[i].Cross(bead[(i + 1) % n]);
    return 0.5 * acc.Length();
}

ChVector3d closestOnBead(const std::vector<ChVector3d>& bead, const ChVector3d& p) {
    if (bead.size() < 2)
        return p;
    double best = 1.0e30;
    ChVector3d closest = p;
    const std::size_t n = bead.size();
    for (std::size_t i = 0; i < n; ++i) {
        const ChVector3d a = bead[i];
        const ChVector3d b = bead[(i + 1) % n];
        const ChVector3d ab = b - a;
        const double len2 = ab.Length2();
        if (len2 < 1.0e-24)
            continue;
        const double t = std::clamp((p - a).Dot(ab) / len2, 0.0, 1.0);
        const ChVector3d q = a + ab * t;
        const double d = (p - q).Length();
        if (d < best) {
            best = d;
            closest = q;
        }
    }
    return closest;
}

ChVector3d closestPointOnTriangle(const ChVector3d& p, const ChVector3d& a, const ChVector3d& b,
                                  const ChVector3d& c) {
    const ChVector3d ab = b - a;
    const ChVector3d ac = c - a;
    const ChVector3d ap = p - a;
    const double d1 = ab.Dot(ap);
    const double d2 = ac.Dot(ap);
    if (d1 <= 0.0 && d2 <= 0.0)
        return a;

    const ChVector3d bp = p - b;
    const double d3 = ab.Dot(bp);
    const double d4 = ac.Dot(bp);
    if (d3 >= 0.0 && d4 <= d3)
        return b;

    const double vc = d1 * d4 - d3 * d2;
    if (vc <= 0.0 && d1 >= 0.0 && d3 <= 0.0)
        return a + ab * (d1 / (d1 - d3));

    const ChVector3d cp = p - c;
    const double d5 = ab.Dot(cp);
    const double d6 = ac.Dot(cp);
    if (d6 >= 0.0 && d5 <= d6)
        return c;

    const double vb = d5 * d2 - d1 * d6;
    if (vb <= 0.0 && d2 >= 0.0 && d6 <= 0.0)
        return a + ac * (d2 / (d2 - d6));

    const double va = d3 * d6 - d5 * d4;
    if (va <= 0.0 && (d4 - d3) >= 0.0 && (d5 - d6) >= 0.0)
        return b + (c - b) * ((d4 - d3) / ((d4 - d3) + (d5 - d6)));

    const double denom = 1.0 / (va + vb + vc);
    return a + ab * (vb * denom) + ac * (vc * denom);
}

void barycentric(const ChVector3d& p, const ChVector3d& a, const ChVector3d& b, const ChVector3d& c, double& wa,
                 double& wb, double& wc) {
    const ChVector3d v0 = b - a;
    const ChVector3d v1 = c - a;
    const ChVector3d v2 = p - a;
    const double d00 = v0.Dot(v0);
    const double d01 = v0.Dot(v1);
    const double d11 = v1.Dot(v1);
    const double d20 = v2.Dot(v0);
    const double d21 = v2.Dot(v1);
    const double denom = d00 * d11 - d01 * d01;
    if (std::abs(denom) < 1.0e-20) {
        wa = 1.0;
        wb = 0.0;
        wc = 0.0;
        return;
    }
    wb = (d11 * d20 - d01 * d21) / denom;
    wc = (d00 * d21 - d01 * d20) / denom;
    wa = 1.0 - wb - wc;
}

bool tetContainsPoint(const ChVector3d& a, const ChVector3d& b, const ChVector3d& c, const ChVector3d& d,
                      const ChVector3d& p) {
    const ChVector3d ab = b - a;
    const ChVector3d ac = c - a;
    const ChVector3d ad = d - a;
    const double den = ab.Dot(ac.Cross(ad));
    if (std::abs(den) < 1.0e-18)
        return false;
    const ChVector3d ap = p - a;
    const double u = ap.Dot(ac.Cross(ad)) / den;
    const double v = ab.Dot(ap.Cross(ad)) / den;
    const double w = ab.Dot(ac.Cross(ap)) / den;
    const double s = 1.0 - u - v - w;
    return u > 1.0e-4 && v > 1.0e-4 && w > 1.0e-4 && s > 1.0e-4;
}

ContactGains continuumGains(double youngs, double patchArea, double nodeMass) {
    const double patch = std::sqrt(std::max(patchArea, 0.0));
    const double k = kContactPenaltyScale * std::max(1.0e5, youngs) * std::max(patch, 1.0e-12);
    if (nodeMass <= 0.0 || k <= 0.0)
        return {k, 0.0};
    return {k, kContactDampingFraction * 2.0 * std::sqrt(k * nodeMass)};
}

ContactGains stableGains(double stiffness, double damping, double nodeMass, double dt) {
    (void)dt;
    const double k = std::max(0.0, stiffness);
    if (k <= 0.0)
        return {0.0, 0.0};
    if (nodeMass <= 0.0)
        return {k, 0.0};
    const double critical = 2.0 * std::sqrt(k * nodeMass);
    const double c = std::min(std::max(0.0, damping), kContactDampingFraction * critical);
    return {k, c};
}

double faceArea(const SurfaceMesh& mesh, const std::array<int, 3>& face) {
    const ChVector3d& a = mesh.vertices[static_cast<std::size_t>(face[0])];
    const ChVector3d& b = mesh.vertices[static_cast<std::size_t>(face[1])];
    const ChVector3d& c = mesh.vertices[static_cast<std::size_t>(face[2])];
    return 0.5 * (b - a).Cross(c - a).Length();
}

double meshVolume(const SurfaceMesh& mesh) {
    double volume = 0.0;
    for (const auto& face : mesh.triangles) {
        const ChVector3d& a = mesh.vertices[static_cast<std::size_t>(face[0])];
        const ChVector3d& b = mesh.vertices[static_cast<std::size_t>(face[1])];
        const ChVector3d& c = mesh.vertices[static_cast<std::size_t>(face[2])];
        volume += a.Dot(b.Cross(c)) / 6.0;
    }
    return volume;
}

double surfaceVolume(const SurfaceMesh& mesh) {
    return std::abs(meshVolume(mesh));
}

bool closedEnough(const SurfaceMesh& mesh) {
    std::map<std::pair<int, int>, int> edges;
    for (const auto& face : mesh.triangles) {
        for (int i = 0; i < 3; ++i)
            edges[{face[i], face[(i + 1) % 3]}]++;
    }
    for (const auto& edge : edges) {
        if (edges.find({edge.first.second, edge.first.first}) == edges.end())
            return false;
    }
    return !mesh.triangles.empty();
}

double budgetEdge(double volume, const MeshSizing& sizing) {
    const int budget = std::max(200, sizing.tetBudget);
    const double h = std::cbrt(std::max(volume, 1.0e-12) * 6.0 / static_cast<double>(budget));
    return std::max(h, sizing.minElement);
}

double targetEdge(double volume, const chrono::ChVector3d& span, const MeshSizing& sizing) {
    const double L = std::max({std::abs(span.x()), std::abs(span.y()), std::abs(span.z()),
                               std::cbrt(std::max(volume, 0.0))});
    const int budget = std::max(200, sizing.tetBudget);
    const double n = static_cast<double>(budget);
    const double hVol = std::cbrt(std::max(volume, 1.0e-12) * 6.0 / n);
    const double hSpan = L / std::cbrt(n);
    return std::max({hVol, hSpan, sizing.minElement});
}

double targetEdge(double volume, const chrono::ChVector3d& span, const MeshSizing& sizing,
                  const Material& material) {
    double h = targetEdge(volume, span, sizing);
    if (material.fracture == FractureMode::Brittle && material.minimumShardSize > 0.0)
        h = std::min(h, material.minimumShardSize);
    return std::max(h, sizing.minElement);
}

int tetBudgetCap(const MeshSizing& sizing, const Material& material) {
    (void)material;
    const int n = std::max(400, sizing.tetBudget);
    constexpr int kSolverCeiling = 2000;
    return std::min(n, kSolverCeiling);
}

double targetEdge(const PartDefinition& part) {
    chrono::ChVector3d span{part.patch.uSize, part.thickness, part.patch.vSize};
    double volume = std::abs(part.patch.uSize * part.patch.vSize * part.thickness);
    if (!part.source.vertices.empty()) {
        chrono::ChVector3d lo = part.source.vertices.front();
        chrono::ChVector3d hi = lo;
        for (const auto& p : part.source.vertices) {
            lo = {std::min(lo.x(), p.x()), std::min(lo.y(), p.y()), std::min(lo.z(), p.z())};
            hi = {std::max(hi.x(), p.x()), std::max(hi.y(), p.y()), std::max(hi.z(), p.z())};
        }
        span = hi - lo;
        volume = std::abs(meshVolume(part.source));
        if (volume < 1.0e-12)
            volume = std::abs(span.x() * span.y() * span.z());
    }
    return targetEdge(volume, span, part.sizing, part.material);
}

double worstTetQuality(const VolumeMesh& mesh) {
    double worst = 1.0;
    for (const auto& tet : mesh.tets) {
        const auto& a = mesh.vertices[static_cast<std::size_t>(tet[0])];
        const auto& b = mesh.vertices[static_cast<std::size_t>(tet[1])];
        const auto& c = mesh.vertices[static_cast<std::size_t>(tet[2])];
        const auto& d = mesh.vertices[static_cast<std::size_t>(tet[3])];
        const double volume = std::abs((b - a).Dot((c - a).Cross(d - a))) / 6.0;
        const double edges = ((b - a).Length() + (c - a).Length() + (d - a).Length() + (c - b).Length() +
                              (d - b).Length() + (d - c).Length()) /
                             6.0;
        if (edges <= 0.0)
            return 0.0;
        worst = std::min(worst, volume / (edges * edges * edges));
    }
    return worst;
}

}  // namespace rcs
