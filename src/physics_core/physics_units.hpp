#pragma once

#include "physics_core/material.hpp"
#include "physics_core/part.hpp"

#include "chrono/core/ChMatrix.h"
#include "chrono/core/ChMatrix33.h"
#include "chrono/core/ChVector3.h"

#include <algorithm>
#include <array>
#include <cmath>
#include <vector>

namespace rcs {

struct ContactGains {
    double stiffness = 0.0;
    double damping = 0.0;
};

double vonMisesStress(double sxx, double syy, double szz, double sxy, double syz, double szx);
double vonMisesStress(const chrono::ChVectorDynamic<>& stress);
double vonMisesPlaneStress(double sxx, double syy, double sxy);
double equivalentPlastic(const double* strain);

/// Isotropic flow stress. Handbook `yieldStress` is the 0.2% offset: σ(εp=0.002)
/// equals yield, and σ(A) tracks UTS. `strengthScale` is the thermal knockdown.
double flowStress(double yield, double ultimate, double failPlastic, double equivalentPlastic,
                  double strengthScale = 1.0);
/// Ductile stretch backup past handbook elongation A. Primary tear is plastic.
double ductileTearStretch(double failPlastic, double ductility);

bool projectRankine(double stress[6], double tension, double compression);
bool projectRankine(chrono::ChVectorDynamic<>& stress, double tension, double compression);

/// Hencky (log-stretch) strain in Voigt order with engineering shear.
/// Linear small-strain Poisson inverts a tet past λ ≈ 1+1/ν; this stays
/// defined for any positive stretch, which is what a ν→1/2 elastomer needs.
void henckyVoigt(const chrono::ChMatrix33<>& stretch, chrono::ChVectorDynamic<>& strain);
/// Mooney-Rivlin Cauchy stress in Voigt order from the corotated stretch U.
/// C10+C01 = μ/2 from E,ν with C01/C10 = 1/8 so equibiaxial differs from
/// uniaxial. Thermal strain scales U by 1/(1+αΔT). Used for yield=0 solids.
void mooneyRivlinCauchy(const chrono::ChMatrix33<>& stretch, double youngs, double nu,
                        double thermalStrain, chrono::ChVectorDynamic<>& stress);
/// True when tets use Mooney-Rivlin. Same predicate as ThermalTet4: the yield
/// passed in is 0 when plasticity is off, and Rankine is off. Shell placement
/// asks compatibility(), not this flag.
bool volumeHyperelastic(const Material& material);

/// Strength multiplier for opening: 1 across the grain, 3 along the fiber.
inline double grainOpeningScale(const Material& material, const chrono::ChVector3d& tensionAxis) {
    if (material.breakFabric != BreakFabric::Grain)
        return 1.0;
    chrono::ChVector3d grain(material.grainDirection[0], material.grainDirection[1],
                             material.grainDirection[2]);
    if (grain.Length2() < 0.25 || tensionAxis.Length2() < 0.25)
        return 1.0;
    grain.Normalize();
    chrono::ChVector3d axis = tensionAxis;
    axis.Normalize();
    return 1.0 + 2.0 * std::abs(axis.Dot(grain));
}

/// Griffith release rate G_c = K_Ic^2 (1-ν^2) / E, in J/m^2.
/// `fractureToughness` is K_Ic (Pa*sqrt(m)). Zero means the stress check stands alone.
inline double modeIReleaseRate(const Material& material) {
    if (material.fractureToughness <= 0.0 || material.youngsModulus <= 1.0)
        return 0.0;
    const double nu = std::clamp(material.poissonRatio, -0.9, 0.49);
    const double k = material.fractureToughness;
    return k * k * std::max(0.0, 1.0 - nu * nu) / material.youngsModulus;
}

/// Separation at which a linear cohesive law has spent its energy, metres.
/// `releaseRate` is G_c in J/m^2. With no toughness the throat itself is the scale.
inline double cohesiveCriticalOpening(double strengthPa, double releaseRate, double throat) {
    if (strengthPa > 1.0 && releaseRate > 0.0)
        return std::max(2.0 * releaseRate / strengthPa, 1.0e-6);
    return std::max({0.25 * std::max(throat, 0.0), 1.0e-3, 1.0e-6});
}

/// Bilinear cohesive damage in [0, 1]. Exhausted means the bead or ply face
/// may separate. The parent continuum is not part of this check.
struct CohesiveDamage {
    double damage = 0.0;
    bool exhausted = false;
};

inline CohesiveDamage cohesiveDamage(double peelOpening, double shearOpening, double peelForce, double shearForce,
                                     double peelCapacity, double shearCapacity, double criticalOpening) {
    CohesiveDamage out;
    const double peelUtil = peelCapacity > 1.0e-12 ? std::abs(peelForce) / peelCapacity : 0.0;
    const double shearUtil = shearCapacity > 1.0e-12 ? std::abs(shearForce) / shearCapacity : 0.0;
    const double openUtil = criticalOpening > 1.0e-12
                                ? std::max(std::abs(peelOpening), std::abs(shearOpening)) / criticalOpening
                                : 0.0;
    out.damage = std::clamp(std::max({peelUtil, shearUtil, openUtil}), 0.0, 1.0);
    out.exhausted = out.damage >= 1.0;
    return out;
}

/// True when the element stores at least G_c times the prospective crack area.
/// A material with no published toughness always passes.
inline bool energyAllowsOpening(const Material& material, double energy, double area) {
    const double gc = modeIReleaseRate(material);
    if (!(gc > 0.0))
        return true;
    return energy >= gc * std::max(area, 0.0);
}

/// Crack normal for a bulk opening. Grain and a composite ply keep the fiber
/// inside the crack plane, so the normal is the tension axis with the fiber removed.
inline chrono::ChVector3d fabricCrackNormal(const Material& material, const chrono::ChVector3d& tensionAxis) {
    chrono::ChVector3d axis = tensionAxis;
    if (axis.Length2() < 0.25)
        axis = {1, 0, 0};
    else
        axis.Normalize();
    chrono::ChVector3d fiber(material.grainDirection[0], material.grainDirection[1], material.grainDirection[2]);
    if (!material.plies.empty()) {
        fiber = {material.plies.front().direction[0], material.plies.front().direction[1],
                 material.plies.front().direction[2]};
    } else if (material.breakFabric != BreakFabric::Grain) {
        return axis;
    }
    if (fiber.Length2() < 0.25)
        return axis;
    fiber.Normalize();
    chrono::ChVector3d normal = axis - fiber * axis.Dot(fiber);
    if (normal.Length2() < 1.0e-8)
        return axis;
    normal.Normalize();
    return normal;
}

bool interfaceShouldOpen(const Material& material, double kelvin, double tension,
                         const chrono::ChVector3d& tensionAxis, double energy, double area);
bool brittleShouldOpen(const Material& material, double kelvin, double tension, double trialEnergy,
                       double tensileStrain, double faceArea, double compression = 0.0, double mean = 0.0,
                       const chrono::ChVector3d& tensionAxis = chrono::ChVector3d(0, 0, 0));
/// Rankine tensile open: max principal at/above handbook σ_t. Mixed-mode
/// contact (a corner hit) still opens; crush is the compressive cap below.
bool brittleShouldCrush(const Material& material, double kelvin, double compression);

int findRoot(std::vector<int>& set, int i);
bool unite(std::vector<int>& set, int a, int b);

double distanceToBead(const std::vector<chrono::ChVector3d>& bead, const chrono::ChVector3d& p);
double beadPerimeter(const std::vector<chrono::ChVector3d>& bead);
double closedBeadArea(const std::vector<chrono::ChVector3d>& bead);
chrono::ChVector3d closestOnBead(const std::vector<chrono::ChVector3d>& bead, const chrono::ChVector3d& p);

chrono::ChVector3d closestPointOnTriangle(const chrono::ChVector3d& p, const chrono::ChVector3d& a,
                                          const chrono::ChVector3d& b, const chrono::ChVector3d& c);
void barycentric(const chrono::ChVector3d& p, const chrono::ChVector3d& a, const chrono::ChVector3d& b,
                 const chrono::ChVector3d& c, double& wa, double& wb, double& wc);
bool tetContainsPoint(const chrono::ChVector3d& a, const chrono::ChVector3d& b, const chrono::ChVector3d& c,
                      const chrono::ChVector3d& d, const chrono::ChVector3d& p);

/// Smallest HHT increment the live path will take. Failed Newton steps cut
/// h down to this; impactSubsteps is bounded by ceil(dt / this).
constexpr double kHhtMinStep = 1.0e-5;
/// Sentinel: no live node is approaching an arena plane.
constexpr double kNoContactTime = 1.0e30;
/// Penalty scale on element stiffness E √A.
constexpr double kContactPenaltyScale = 0.1;
/// Viscous contact damping as a fraction of critical on the material k.
constexpr double kContactDampingFraction = 0.5;
/// Softened-contact depth, metres. Below it pressure rises quadratically
/// (k δ²/2δ₀) and damping ramps in, so force and tangent are continuous at
/// touch; beyond it the penalty is linear, k (δ − δ₀/2).
constexpr double kContactSoftening = 1.0e-5;
/// Largest gap the penalty may use. MixedJoints glass seating produced
/// Newton iterates at gap ~ 1e99 m (fn ~ 1e109) while committed energy
/// stayed millijoules; a 10 mm cap leaves real seating (microns) unchanged.
constexpr double kContactGapCap = 1.0e-2;
/// Extra critical damping while the gap is closing, Abaqus-style contact
/// stabilization. The ramp is 0 at rest so a seated body is not glued.
constexpr double kContactStabFraction = 1.0;

/// Slowest approach treated as an impact. Hunt-Crossley damping grows as
/// 1/v⁻, which is singular for resting contact; a node that arrives slower
/// than this is resting contact and gets the damping of an impact at it.
constexpr double kContactImpactSpeedFloor = 0.1;

/// Flores et al. (2011) Hunt-Crossley factor D in F = F_s (1 + D dδ/dt),
/// D = 8(1 - e) / (5 e v⁻). A struck mass rebounds at restitution e.
inline double huntCrossleyFactor(double restitution, double approachSpeed) {
    const double e = std::clamp(restitution, 0.02, 0.98);
    const double v = std::max(approachSpeed, kContactImpactSpeedFloor);
    return 8.0 * (1.0 - e) / (5.0 * e * v);
}

/// Hunt-Crossley multiplier on the spring force, F = F_s · s, and the
/// dashpot coefficient −F_s ds/dvn per unit F_s. Closing (vn ≤ 0) is the
/// linear 1 − D vn. Separating uses e^{−D vn}, which matches its value and
/// slope at vn = 0 but never reaches zero: a hard tensile cut at vn = 1/D
/// (about 1 cm/s for slow contact) is a kink Newton cycles across.
struct HuntCrossleyScale {
    double scale = 1.0;
    double damping = 0.0;
};
inline HuntCrossleyScale huntCrossleyScale(double D, double vn, double approachSpeed = kContactImpactSpeedFloor) {
    // D is sized at v⁻. A Newton overshoot with |vn| ≫ v⁻ must not
    // evaluate 1 + D|vn| at that iterate: MixedJoints seating hit fn=inf
    // at vn ~ 1e209 with vApp=0 (D from the 0.1 m/s floor).
    const double vRef = std::max(std::max(0.0, approachSpeed), kContactImpactSpeedFloor);
    if (vn <= 0.0) {
        const double vnLaw = std::max(vn, -vRef);
        return {1.0 - D * vnLaw, D};
    }
    const double decay = std::exp(-D * vn);
    return {decay, D * decay};
}

/// Penalty dashpot: 0.5-of-critical always, plus stabilization that tracks
/// closing speed and vanishes when the node is not approaching the plane.
inline double contactDashpot(double stiffness, double nodeMass, double vn, double stickSpeed) {
    if (stiffness <= 0.0 || nodeMass <= 0.0)
        return 0.0;
    const double critical = 2.0 * std::sqrt(stiffness * nodeMass);
    const double closing = std::max(0.0, -vn);
    const double vRef = std::max(stickSpeed, 1.0e-3);
    const double stab = std::min(1.0, closing / vRef) * kContactStabFraction;
    return (kContactDampingFraction + stab) * critical;
}

/// Time to close an open gap under constant closing speed and acceleration.
/// `openGap` is distance outside the plane. `closingSpeed` / `closingAcc` are
/// positive toward the plane (floor drop: closingAcc = g). Same residual as
/// integrating gravity; used only to size a free-flight HHT increment.
inline double timeToClosingPlane(double openGap, double closingSpeed, double closingAcc) {
    if (openGap <= 1.0e-12)
        return 0.0;
    if (closingAcc > 1.0e-12) {
        const double disc = closingSpeed * closingSpeed + 2.0 * closingAcc * openGap;
        if (disc < 0.0)
            return kNoContactTime;
        const double t = (-closingSpeed + std::sqrt(disc)) / closingAcc;
        return t >= 0.0 ? t : kNoContactTime;
    }
    if (closingSpeed > 1.0e-12)
        return openGap / closingSpeed;
    return kNoContactTime;
}

/// Continuum contact: k = 1.2 E √A, c = 0.2 of critical. No m/dt² cap —
/// if k is too stiff for the current step, the integrator cuts h.
ContactGains continuumGains(double youngs, double patchArea, double nodeMass);
/// Keep requested stiffness; cap damping at 0.2 of critical. `dt` unused.
ContactGains stableGains(double stiffness, double damping, double nodeMass, double dt);

double faceArea(const SurfaceMesh& mesh, const std::array<int, 3>& face);
double meshVolume(const SurfaceMesh& mesh);
double surfaceVolume(const SurfaceMesh& mesh);
bool closedEnough(const SurfaceMesh& mesh);
double budgetEdge(double volume, const MeshSizing& sizing);
/// Uniform target edge from bounding span, volume, and tetBudget.
double targetEdge(double volume, const chrono::ChVector3d& span, const MeshSizing& sizing);
/// Same geometric field, then FractureMode and handbook lengths scale it.
double targetEdge(double volume, const chrono::ChVector3d& span, const MeshSizing& sizing,
                  const Material& material);
double targetEdge(const PartDefinition& part);
/// Tet count the filler may emit: tetBudget, then a solver allocation ceiling.
int tetBudgetCap(const MeshSizing& sizing, const Material& material);
double worstTetQuality(const VolumeMesh& mesh);

/// Local geometric edge request: about three tets through the wall, and
/// enough samples on a curve of curvature κ. Huge when both are unknown.
inline double geometricEdgeRequest(double thickness, double curvature) {
    double h = 1.0e30;
    if (thickness > 1.0e-12)
        h = std::min(h, thickness / 3.0);
    if (curvature > 1.0e-12)
        h = std::min(h, 0.5 / curvature);
    return h;
}

/// Clamp a vertex size field. Known geometric requests stay, so a thin web is
/// not coarsened to the part-wide mean. Unknown entries become `target`.
/// `target` is a ceiling (budget / resolution): h = min(h_geom, target).
inline void clampSizeField(std::vector<double>& field, double target, double minH, double maxH) {
    const double lo = std::max(minH, 1.0e-9);
    const double hi = std::max(lo, maxH);
    const double goal = std::clamp(target, lo, hi);
    for (double& h : field) {
        if (h > 1.0e-12 && h < 1.0e29)
            h = std::min(h, goal);
        else
            h = goal;
        h = std::clamp(h, lo, hi);
    }
}

}  // namespace rcs
