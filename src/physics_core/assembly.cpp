#include "physics_core/assembly.hpp"

#include "chrono/core/ChMatrix33.h"
#include "chrono/core/ChRotation.h"

#include <algorithm>
#include <cmath>
#include <limits>
#include <unordered_map>

namespace rcs {
namespace {

using chrono::ChVector3d;

/// A link is failed only after it has been over its allowable for this many
/// consecutive steps, so a single solver spike cannot unzip a weld.
constexpr int kOverloadHoldSteps = 4;
/// Time constant of the reaction low-pass filter, seconds.
constexpr double kReactionFilterTau = 2.0e-3;

int findRoot(std::vector<int>& set, int i) {
    while (set[i] != i) {
        set[i] = set[set[i]];
        i = set[i];
    }
    return i;
}

/// Returns false when the two entries were already connected, i.e. the caller
/// is about to add a redundant constraint.
bool unite(std::vector<int>& set, int a, int b) {
    const int ra = findRoot(set, a);
    const int rb = findRoot(set, b);
    if (ra == rb)
        return false;
    set[rb] = ra;
    return true;
}

ChVector3d closestPointOnTriangle(const ChVector3d& p,
                                  const ChVector3d& a,
                                  const ChVector3d& b,
                                  const ChVector3d& c)
{
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

void barycentric(const ChVector3d& p,
                 const ChVector3d& a,
                 const ChVector3d& b,
                 const ChVector3d& c,
                 double& wa,
                 double& wb,
                 double& wc)
{
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

/// Contact is applied as an external nodal force, so the integrator sees it
/// explicitly and an over-stiff spring will shake itself apart. A node of mass
/// m stepped at dt is stable while k < 4m/dt^2 and c < 2m/dt; clamping to a
/// quarter of each keeps a comfortable margin and, more importantly, means the
/// user can dial in any stiffness without being able to detonate the run.
struct ContactGains {
    double stiffness;
    double damping;
};

/// Inbound speed above this is an impact. Slower motion is a resting contact
/// and must not be treated as a rebound.
constexpr double kRestSpeed = 0.05;

ContactGains stableGains(double stiffness, double damping, double nodeMass, double dt) {
    if (nodeMass <= 0.0 || dt <= 0.0)
        return {0.0, 0.0};
    const double k = std::min(stiffness, nodeMass / (dt * dt));
    const double criticalDamping = 2.0 * std::sqrt(k * nodeMass);
    const double c = std::min({damping, criticalDamping, 0.5 * nodeMass / dt});
    return {k, c};
}

/// Soft enough that critical damping fits inside the explicit stability cap
/// `c <= 0.5 m/dt`. A settled node then sits in an overdamped sag instead of
/// ringing at the timestep.
ContactGains restGains(double nodeMass, double dt) {
    if (nodeMass <= 0.0 || dt <= 0.0)
        return {0.0, 0.0};
    const double k = nodeMass / (16.0 * dt * dt);
    const double c = std::min(2.0 * std::sqrt(k * nodeMass), 0.5 * nodeMass / dt);
    return {k, c};
}

bool aabbOverlap(const ChVector3d& loA, const ChVector3d& hiA,
                 const ChVector3d& loB, const ChVector3d& hiB)
{
    return loA.x() <= hiB.x() && hiA.x() >= loB.x() &&
           loA.y() <= hiB.y() && hiA.y() >= loB.y() &&
           loA.z() <= hiB.z() && hiA.z() >= loB.z();
}

}  // namespace

void Assembly::build(chrono::ChSystem& system, const AssemblyDefinition& definition) {
    _parts.clear();
    _joints.clear();
    _deferred.clear();
    _system = &system;
    _nodeOffset.clear();
    _loadTargets.clear();
    _contactArmed.clear();
    _floorPeak.clear();

    _parts.resize(definition.parts.size());
    std::size_t offset = 0;
    for (std::size_t i = 0; i < definition.parts.size(); ++i) {
        _parts[i].build(definition.parts[i]);
        _nodeOffset.push_back(offset);
        offset += _parts[i].nodes().size();
        system.Add(_parts[i].mesh());
    }

    // Coincident nodes welded in a loop (three plates meeting at a box corner,
    // a bolt circle, ...) produce linearly dependent constraint rows, which a
    // direct solver cannot factor. Two union-finds keep the constraint graph a
    // spanning tree: one for the translational locks, one for rotational.
    std::vector<int> translationSet(offset);
    std::vector<int> rotationSet(offset);
    for (std::size_t i = 0; i < offset; ++i) {
        translationSet[i] = static_cast<int>(i);
        rotationSet[i] = static_cast<int>(i);
    }

    for (const auto& weld : definition.welds)
        buildWeld(system, weld, translationSet, rotationSet);
    for (const auto& screw : definition.screws)
        buildScrew(system, screw, translationSet, rotationSet);

    _authoredCenter = center();
    _placementRotation = chrono::QUNIT;
    _placementPosition = _authoredCenter;
    refreshAdjacency();
}

void Assembly::createLink(chrono::ChSystem& system,
                          Joint& joint,
                          std::size_t partA,
                          std::size_t nodeA,
                          std::size_t partB,
                          std::size_t nodeB,
                          double breakForce,
                          double breakTorque,
                          bool wantRotationLock,
                          bool lockRotation)
{
    using namespace chrono;

    JointLink entry;
    entry.partA = partA;
    entry.nodeA = nodeA;
    entry.partB = partB;
    entry.nodeB = nodeB;
    entry.breakForce = breakForce;
    entry.breakTorque = breakTorque;
    entry.locksRotation = lockRotation;

    if (!_parts[partA].rotational() && !_parts[partB].rotational()) {
        auto link = chrono_types::make_shared<fea::ChLinkNodeNode>();
        link->Initialize(_parts[partA].xyzNodes()[nodeA], _parts[partB].xyzNodes()[nodeB]);
        system.Add(link);
        entry.pointLink = link;
        entry.locksRotation = false;
    } else {
        const auto& a = _parts[partA].nodes()[nodeA];
        const auto& b = _parts[partB].nodes()[nodeB];
        auto link = chrono_types::make_shared<ChLinkMateGeneric>(
            true, true, true, lockRotation, lockRotation, lockRotation);
        link->Initialize(a, b, ChFrame<>(0.5 * (a->GetPos() + b->GetPos()), a->GetRot()));
        system.Add(link);
        entry.link = link;
    }
    (void)wantRotationLock;
    joint.links.push_back(entry);
}

void Assembly::addLink(chrono::ChSystem& system,
                       Joint& joint,
                       std::size_t partA,
                       std::size_t nodeA,
                       std::size_t partB,
                       std::size_t nodeB,
                       double breakForce,
                       double breakTorque,
                       bool wantRotationLock,
                       std::vector<int>& translationSet,
                       std::vector<int>& rotationSet)
{
    const int ga = static_cast<int>(globalNode(partA, nodeA));
    const int gb = static_cast<int>(globalNode(partB, nodeB));

    if (!unite(translationSet, ga, gb)) {
        DeferredLink skipped;
        skipped.partA = partA;
        skipped.nodeA = nodeA;
        skipped.partB = partB;
        skipped.nodeB = nodeB;
        skipped.breakForce = breakForce;
        skipped.breakTorque = breakTorque;
        skipped.wantRotation = wantRotationLock;
        skipped.label = joint.label;
        skipped.kind = joint.kind;
        _deferred.push_back(skipped);
        return;
    }

    const bool lockRotation = wantRotationLock && unite(rotationSet, ga, gb);
    createLink(system, joint, partA, nodeA, partB, nodeB, breakForce, breakTorque, wantRotationLock,
               lockRotation);
}

void Assembly::buildWeld(chrono::ChSystem& system,
                         const WeldDefinition& definition,
                         std::vector<int>& translationSet,
                         std::vector<int>& rotationSet)
{
    if (definition.partA >= _parts.size() || definition.partB >= _parts.size() ||
        definition.partA == definition.partB)
        return;

    const Part& a = _parts[definition.partA];
    const Part& b = _parts[definition.partB];

    auto insideRegion = [&](const ChVector3d& p) {
        if (!definition.useRegion)
            return true;
        return p.x() >= definition.regionMin.x() && p.x() <= definition.regionMax.x() &&
               p.y() >= definition.regionMin.y() && p.y() <= definition.regionMax.y() &&
               p.z() >= definition.regionMin.z() && p.z() <= definition.regionMax.z();
    };

    // Pair boundary nodes that already sit on top of each other: that is what
    // "these two edges are welded" means geometrically, and it works for any
    // CAD shape rather than just a cube.
    struct Pair {
        std::size_t nodeA;
        std::size_t nodeB;
        ChVector3d position;
    };
    std::vector<Pair> pairs;
    for (int ia : a.boundaryNodes()) {
        const ChVector3d pa = a.nodePosition(static_cast<std::size_t>(ia));
        if (!insideRegion(pa))
            continue;
        std::size_t best = b.nodeCount();
        double bestDistance = definition.tolerance;
        for (int ib : b.boundaryNodes()) {
            const ChVector3d pb = b.nodePosition(static_cast<std::size_t>(ib));
            if (!insideRegion(pb))
                continue;
            const double distance = (pa - pb).Length();
            if (distance <= bestDistance) {
                bestDistance = distance;
                best = static_cast<std::size_t>(ib);
            }
        }
        if (best < b.nodeCount())
            pairs.push_back({static_cast<std::size_t>(ia), best,
                             0.5 * (pa + b.nodePosition(best))});
    }
    if (pairs.empty())
        return;

    const double throat =
        definition.throat > 0.0 ? definition.throat : std::min(a.thickness(), b.thickness());
    const double allowable =
        std::min(a.material().weldShearStrength, b.material().weldShearStrength) *
        definition.strengthScale;

    Joint joint;
    joint.kind = JointKind::Weld;
    joint.label = definition.label;
    joint.partA = definition.partA;
    joint.partB = definition.partB;

    for (std::size_t i = 0; i < pairs.size(); ++i) {
        // Tributary bead length: distance to the nearest neighbouring link.
        // Along an edge weld this is the nodal spacing, which is what the
        // strength of one link should be based on.
        double spacing = 0.0;
        if (pairs.size() > 1) {
            double nearest = std::numeric_limits<double>::max();
            for (std::size_t j = 0; j < pairs.size(); ++j) {
                if (j == i)
                    continue;
                nearest = std::min(nearest, (pairs[i].position - pairs[j].position).Length());
            }
            spacing = nearest;
        } else {
            spacing = throat;
        }

        const double area = spacing * throat;
        const double breakForce = allowable * area;
        const double breakTorque = breakForce * throat;
        addLink(system, joint, definition.partA, pairs[i].nodeA, definition.partB, pairs[i].nodeB,
                breakForce, breakTorque, true, translationSet, rotationSet);
    }

    if (!joint.links.empty())
        _joints.push_back(std::move(joint));
}

void Assembly::buildScrew(chrono::ChSystem& system,
                          const ScrewDefinition& definition,
                          std::vector<int>& translationSet,
                          std::vector<int>& rotationSet)
{
    if (definition.partA >= _parts.size() || definition.partB >= _parts.size() ||
        definition.partA == definition.partB)
        return;

    auto nearestNode = [&](const Part& part) -> std::size_t {
        std::size_t best = part.nodeCount();
        double bestDistance = definition.clampRadius;
        for (std::size_t i = 0; i < part.nodeCount(); ++i) {
            const double distance = (part.nodePosition(i) - definition.position).Length();
            if (distance <= bestDistance) {
                bestDistance = distance;
                best = i;
            }
        }
        return best;
    };

    const std::size_t nodeA = nearestNode(_parts[definition.partA]);
    const std::size_t nodeB = nearestNode(_parts[definition.partB]);
    if (nodeA >= _parts[definition.partA].nodeCount() ||
        nodeB >= _parts[definition.partB].nodeCount())
        return;

    // Tensile stress area of a standard thread is about 75% of the shank.
    const double d = definition.shankDiameter;
    const double stressArea = 0.75 * 0.25 * chrono::CH_PI * d * d;
    const double breakForce = definition.proofStress * stressArea;

    Joint joint;
    joint.kind = JointKind::Screw;
    joint.label = definition.label;
    joint.partA = definition.partA;
    joint.partB = definition.partB;
    addLink(system, joint, definition.partA, nodeA, definition.partB, nodeB, breakForce,
            breakForce * d * 0.5, definition.locksRotation, translationSet, rotationSet);

    if (!joint.links.empty())
        _joints.push_back(std::move(joint));
}

void Assembly::place(const InitialState& state) {
    using namespace chrono;

    const ChQuaterniond rotation =
        QuatFromAngleZ(state.orientation.z()) *
        QuatFromAngleY(state.orientation.y()) *
        QuatFromAngleX(state.orientation.x());

    const ChVector3d pivot = _authoredCenter;
    _placementRotation = rotation;
    _placementPosition = state.position;
    for (auto& part : _parts)
        part.transform(state.position - pivot, rotation, pivot, true);

    const ChVector3d newPivot = center();
    for (auto& part : _parts)
        part.setVelocity(state.linearVelocity, state.angularVelocity, newPivot);

    // The links were initialised at the authoring pose; re-seat them so they
    // carry no constraint violation at t = 0.
    for (auto& joint : _joints) {
        for (auto& entry : joint.links) {
            if (entry.link && _parts[entry.partA].rotational() && _parts[entry.partB].rotational()) {
                const auto& a = _parts[entry.partA].nodes()[entry.nodeA];
                const auto& b = _parts[entry.partB].nodes()[entry.nodeB];
                entry.link->Initialize(
                    a, b, ChFrame<>(0.5 * (a->GetPos() + b->GetPos()), a->GetRot()));
            } else if (entry.pointLink) {
                entry.pointLink->Initialize(_parts[entry.partA].xyzNodes()[entry.nodeA],
                                            _parts[entry.partB].xyzNodes()[entry.nodeB]);
            }
        }
    }
}

void Assembly::bindLoads(const std::vector<ImpactLoad>& loads) {
    _loadTargets.assign(loads.size(), {});

    for (std::size_t i = 0; i < loads.size(); ++i) {
        const ImpactLoad& load = loads[i];
        if (load.partIndex >= _parts.size())
            continue;
        const Part& part = _parts[load.partIndex];

        double totalMass = 0.0;
        std::vector<LoadTarget> targets;
        for (std::size_t n = 0; n < part.nodeCount(); ++n) {
            if ((part.nodePosition(n) - load.point).Length() > load.radius)
                continue;
            targets.push_back({load.partIndex, n, part.nodeMass()[n]});
            totalMass += part.nodeMass()[n];
        }

        // Nothing inside the patch: fall back to the single closest node so a
        // mis-placed load still does something visible.
        if (targets.empty()) {
            std::size_t best = 0;
            double bestDistance = std::numeric_limits<double>::max();
            for (std::size_t n = 0; n < part.nodeCount(); ++n) {
                const double distance = (part.nodePosition(n) - load.point).Length();
                if (distance < bestDistance) {
                    bestDistance = distance;
                    best = n;
                }
            }
            if (part.nodeCount() > 0)
                targets.push_back({load.partIndex, best, 1.0});
            totalMass = 1.0;
        }

        if (totalMass <= 0.0)
            totalMass = static_cast<double>(targets.size());
        for (auto& target : targets)
            target.share /= totalMass;
        _loadTargets[i] = std::move(targets);
    }
}

Assembly::LoadAdvisory Assembly::adviseLoad(std::size_t loadIndex,
                                            const std::vector<ImpactLoad>& loads,
                                            double dt) const
{
    LoadAdvisory advisory;
    if (loadIndex >= loads.size() || loadIndex >= _loadTargets.size())
        return advisory;

    const double total = loads[loadIndex].force;
    for (const auto& target : _loadTargets[loadIndex]) {
        const double force = total * target.share;
        const double mass = _parts[target.part].nodeMass()[target.node];
        advisory.peakNodalForce = std::max(advisory.peakNodalForce, force);
        if (mass > 0.0)
            advisory.deltaVPerStep = std::max(advisory.deltaVPerStep, force * dt / mass);
    }
    advisory.nodeCount = _loadTargets[loadIndex].size();
    return advisory;
}

void Assembly::applyLoads(const Environment& environment,
                          const std::vector<ImpactLoad>& loads,
                          double time,
                          double dt)
{
    for (auto& part : _parts)
        part.clearLoads();

    applyGravity(environment.gravity);
    if (environment.aerodynamicsEnabled)
        applyAerodynamics(environment);
    if (environment.floorEnabled)
        applyFloorContact(environment, dt);
    applyPartContact(environment, dt);
    applyShardContact(environment, dt);
    applyUserLoads(loads, time);

    for (auto& part : _parts)
        part.commitLoads();
}

void Assembly::applyGravity(const chrono::ChVector3d& gravity) {
    for (auto& part : _parts) {
        for (std::size_t i = 0; i < part.nodeCount(); ++i) {
            if (part.nodeRigid(i))
                continue;
            part.addLoad(i, part.nodeMass()[i] * gravity);
        }
    }
}

void Assembly::applyAerodynamics(const Environment& environment) {
    const double factor = 0.5 * environment.airDensity * environment.dragCoefficient;
    if (factor <= 0.0)
        return;

    auto dragFace = [&](Part& part, const int* index, int count, const ChVector3d& areaVector) {
        const double area = areaVector.Length();
        if (area < 1.0e-12)
            return;
        const ChVector3d normal = areaVector / area;
        ChVector3d velocity{0, 0, 0};
        for (int i = 0; i < count; ++i)
            velocity += part.nodeVelocity(static_cast<std::size_t>(index[i]));
        velocity /= static_cast<double>(count);
        const double vn = velocity.Dot(normal);
        const ChVector3d force = normal * (-factor * area * std::abs(vn) * vn);
        for (int i = 0; i < count; ++i)
            part.addLoad(static_cast<std::size_t>(index[i]), force / static_cast<double>(count));
    };

    for (auto& part : _parts) {
        if (!part.quads().empty()) {
            for (const auto& quad : part.quads()) {
                const ChVector3d p0 = part.nodePosition(quad[0]);
                const ChVector3d p1 = part.nodePosition(quad[1]);
                const ChVector3d p2 = part.nodePosition(quad[2]);
                const ChVector3d p3 = part.nodePosition(quad[3]);
                dragFace(part, quad.data(), 4, 0.5 * (p2 - p0).Cross(p3 - p1));
            }
        } else {
            for (const auto& tri : part.triangles()) {
                const ChVector3d a = part.nodePosition(tri[0]);
                const ChVector3d b = part.nodePosition(tri[1]);
                const ChVector3d c = part.nodePosition(tri[2]);
                dragFace(part, tri.data(), 3, 0.5 * (b - a).Cross(c - a));
            }
        }
    }
}

void Assembly::applyShardContact(const Environment& environment, double dt) {
    using namespace chrono;
    (void)environment;
    for (auto& owner : _parts) {
        if (owner.shards().empty())
            continue;
        const double ownerReach = 0.5 * owner.thickness();
        for (auto& shard : owner.shards()) {
            if (shard.mass <= 0.0)
                continue;
            ChVector3d corner[4];
            ChVector3d lo{1.0e9, 1.0e9, 1.0e9};
            ChVector3d hi{-1.0e9, -1.0e9, -1.0e9};
            for (int k = 0; k < 4; ++k) {
                corner[k] = shard.com + shard.rotation.Rotate(shard.local[static_cast<std::size_t>(k)]);
                lo = {std::min(lo.x(), corner[k].x()), std::min(lo.y(), corner[k].y()), std::min(lo.z(), corner[k].z())};
                hi = {std::max(hi.x(), corner[k].x()), std::max(hi.y(), corner[k].y()), std::max(hi.z(), corner[k].z())};
            }
            const ChVector3d face = (corner[1] - corner[0]).Cross(corner[2] - corner[0]);
            if (face.Length() < 1.0e-12)
                continue;
            const ChVector3d normal = face.GetNormalized();
            for (auto& other : _parts) {
                const double reach = ownerReach + 0.5 * other.thickness();
                for (std::size_t i = 0; i < other.nodeCount(); ++i) {
                    if (other.nodeRigid(i))
                        continue;
                    const ChVector3d p = other.nodePosition(i);
                    if (p.x() < lo.x() - reach || p.x() > hi.x() + reach || p.y() < lo.y() - reach ||
                        p.y() > hi.y() + reach || p.z() < lo.z() - reach || p.z() > hi.z() + reach)
                        continue;
                    const double signedGap = (p - corner[0]).Dot(normal);
                    if (std::abs(signedGap) >= reach)
                        continue;
                    const ChVector3d projected = p - normal * signedGap;
                    double wa = 0, wb = 0, wc = 0, wa2 = 0, wb2 = 0, wc2 = 0;
                    barycentric(projected, corner[0], corner[1], corner[2], wa, wb, wc);
                    barycentric(projected, corner[0], corner[2], corner[3], wa2, wb2, wc2);
                    const bool inside = (wa > 1.0e-3 && wb > 1.0e-3 && wc > 1.0e-3) ||
                                        (wa2 > 1.0e-3 && wb2 > 1.0e-3 && wc2 > 1.0e-3);
                    if (!inside)
                        continue;
                    const ChVector3d outward = signedGap >= 0.0 ? normal : -normal;
                    const double overlap = reach - std::abs(signedGap);
                    const ChVector3d arm = projected - shard.com;
                    const ChVector3d shardVel = shard.velocity + shard.angularVelocity.Cross(arm);
                    const double vn = (other.nodeVelocity(i) - shardVel).Dot(outward);
                    const double mass = other.nodeMass()[i];
                    ContactGains gains = restGains(mass, dt);
                    if (vn < -kRestSpeed && mass > 0.0 && dt > 0.0) {
                        const double k = mass / (8.0 * dt * dt);
                        const double c = std::min(2.0 * std::sqrt(k * mass), 0.5 * mass / dt);
                        gains = {k, c};
                    }
                    const double support = std::min(overlap, std::max(0.0, -vn) * dt + 1.0e-4);
                    const double fn = gains.stiffness * support - gains.damping * vn;
                    if (fn <= 0.0)
                        continue;
                    other.addLoad(i, outward * fn);
                    // A sleeping chip is a fixed obstacle unless the shell is
                    // still arriving hard enough to knock it loose.
                    if (shard.asleep && vn > -0.2)
                        continue;
                    shard.asleep = false;
                    shard.pendingForce -= outward * fn;
                    shard.pendingTorque -= arm.Cross(outward * fn);
                }
            }
        }
    }
}

void Assembly::separateShards() {
    using namespace chrono;
    struct Ref {
        Part* part = nullptr;
        std::size_t index = 0;
    };
    std::vector<Ref> all;
    for (auto& part : _parts) {
        for (std::size_t s = 0; s < part.shards().size(); ++s)
            all.push_back({&part, s});
    }
    auto cornerOf = [](const Part::RigidShard& shard, int k) {
        return shard.com + shard.rotation.Rotate(shard.local[static_cast<std::size_t>(k)]);
    };
    for (std::size_t i = 0; i < all.size(); ++i) {
        for (std::size_t j = i + 1; j < all.size(); ++j) {
            auto& a = all[i].part->shards()[all[i].index];
            auto& b = all[j].part->shards()[all[j].index];
            if (a.asleep && b.asleep)
                continue;
            const double reach = 0.5 * (all[i].part->thickness() + all[j].part->thickness());
            auto push = [&](Part::RigidShard& mover, Part::RigidShard& wall, bool wakeWall) {
                ChVector3d face[4];
                for (int k = 0; k < 4; ++k)
                    face[k] = cornerOf(wall, k);
                const ChVector3d n = (face[1] - face[0]).Cross(face[2] - face[0]);
                if (n.Length() < 1.0e-12)
                    return;
                const ChVector3d normal = n.GetNormalized();
                for (int k = 0; k < 4; ++k) {
                    const ChVector3d p = cornerOf(mover, k);
                    const double signedGap = (p - face[0]).Dot(normal);
                    if (std::abs(signedGap) >= reach)
                        continue;
                    const ChVector3d projected = p - normal * signedGap;
                    double wa = 0, wb = 0, wc = 0, wa2 = 0, wb2 = 0, wc2 = 0;
                    barycentric(projected, face[0], face[1], face[2], wa, wb, wc);
                    barycentric(projected, face[0], face[2], face[3], wa2, wb2, wc2);
                    const bool inside = (wa > 1.0e-3 && wb > 1.0e-3 && wc > 1.0e-3) ||
                                        (wa2 > 1.0e-3 && wb2 > 1.0e-3 && wc2 > 1.0e-3);
                    if (!inside)
                        continue;
                    const ChVector3d outward = signedGap >= 0.0 ? normal : -normal;
                    const double overlap = reach - std::abs(signedGap);
                    mover.com += outward * overlap;
                    const double vn = mover.velocity.Dot(outward);
                    if (vn < 0.0)
                        mover.velocity -= outward * vn;
                    if (wakeWall && vn < -0.2)
                        wall.asleep = false;
                    return;
                }
            };
            if (!a.asleep)
                push(a, b, b.asleep);
            if (!b.asleep)
                push(b, a, a.asleep);
        }
    }
}

void Assembly::integrateShards(double dt, const Environment& environment) {
    bool any = false;
    for (const auto& part : _parts)
        any = any || !part.shards().empty();
    if (!any)
        return;
    for (auto& part : _parts)
        part.advanceShards(dt, environment.gravity);
    separateShards();
    for (auto& part : _parts)
        part.finishShards(environment.floorEnabled, environment.floorHeight);
}

void Assembly::applyFloorContact(const Environment& environment, double dt) {
    const double mu = environment.contactFriction;
    const double vStick = std::max(1.0e-4, environment.frictionRegularization);

    if (_floorPeak.size() != _parts.size())
        _floorPeak.assign(_parts.size(), {});
    if (_floorPhase.size() != _parts.size())
        _floorPhase.assign(_parts.size(), {});
    if (_floorClear.size() != _parts.size())
        _floorClear.assign(_parts.size(), {});

    for (std::size_t p = 0; p < _parts.size(); ++p) {
        auto& part = _parts[p];
        auto& peak = _floorPeak[p];
        auto& phase = _floorPhase[p];
        auto& clear = _floorClear[p];
        if (peak.size() != part.nodeCount())
            peak.assign(part.nodeCount(), 0.0);
        if (phase.size() != part.nodeCount())
            phase.assign(part.nodeCount(), 0);
        if (clear.size() != part.nodeCount())
            clear.assign(part.nodeCount(), 1);

        const double radius = 0.5 * part.thickness();
        // A resting node has to rise this far clear before the next arrival
        // counts as a new impact. One solver step of liftoff is not a bounce.
        constexpr double kClearGap = 1.0e-3;
        for (std::size_t i = 0; i < part.nodeCount(); ++i) {
            if (part.nodeRigid(i))
                continue;
            const double penetration = (environment.floorHeight + radius) - part.nodePosition(i).y();
            const ChVector3d v = part.nodeVelocity(i);
            const double vy = v.y();
            // phase: 0 resting, 1 inbound impact that started in clear air.
            // The support is never scaled down while the node is inside the
            // floor. A resting sag keeps the soft spring the whole time.
            if (penetration <= 0.0) {
                if (penetration < -kClearGap || phase[i] != 0) {
                    peak[i] = 0.0;
                    phase[i] = 0;
                    clear[i] = 1;
                }
                continue;
            }

            const bool entered = clear[i] != 0;
            if (entered && vy < -kRestSpeed)
                phase[i] = 1;
            clear[i] = 0;
            // The stiff spring only acts while the node is still arriving.
            // Keeping it on the way back up, or scaling it down and then
            // catching the fall, is what tips a panel that had nearly stopped.
            if (phase[i] == 1 && vy >= -kRestSpeed) {
                phase[i] = 0;
                peak[i] = 0.0;
            }

            // The usual explicit clamp k = m/dt^2 is applied a step late and
            // pumps a rocking mode. One eighth of that stops a fast arrival
            // and still leaves the contact overdamped.
            ContactGains gains = restGains(part.nodeMass()[i], dt);
            if (phase[i] == 1 && part.nodeMass()[i] > 0.0 && dt > 0.0) {
                const double nodeMass = part.nodeMass()[i];
                const double k = nodeMass / (8.0 * dt * dt);
                const double c = std::min(2.0 * std::sqrt(k * nodeMass), 0.5 * nodeMass / dt);
                gains = {k, c};
            }
            double support = penetration;
            if (entered)
                support = std::min(penetration, std::max(0.0, -vy) * dt + 1.0e-4);
            double fn = gains.stiffness * support - gains.damping * vy;
            if (fn <= 0.0)
                continue;

            ChVector3d force{0.0, fn, 0.0};
            const ChVector3d vt{v.x(), 0.0, v.z()};
            const double vtMagnitude = vt.Length();
            if (vtMagnitude > 1.0e-9) {
                double ft = mu * fn;
                if (vtMagnitude < vStick)
                    ft *= vtMagnitude / vStick;
                force -= vt * (ft / vtMagnitude);
            }
            part.addLoad(i, force);
        }
    }
}

void Assembly::applyPartContact(const Environment& environment, double dt) {
    const std::size_t count = _parts.size();
    if (count == 0)
        return;
    bool anyTorn = false;
    for (const auto& part : _parts)
        anyTorn = anyTorn || part.torn();
    if (count < 2 && !anyTorn)
        return;
    if (_contactArmed.size() != count * count)
        _contactArmed.assign(count * count, 0);
    if (_contactNormal.size() != count * count)
        _contactNormal.assign(count * count, ChVector3d(0, 0, 0));

    const double mu = environment.contactFriction;
    const double vStick = std::max(1.0e-4, environment.frictionRegularization);
    const double restitution = std::clamp(environment.floorRestitution, 0.0, 1.0);
    const double reboundScale = restitution * restitution;

    std::vector<ChVector3d> lo(count), hi(count);
    for (std::size_t i = 0; i < count; ++i)
        _parts[i].bounds(lo[i], hi[i]);

    // Node-against-current-triangle proximity: the contact surface is whatever
    // shape the mesh has right now, deformed or not. Returns true when the two
    // parts overlap at all, which is what arming uses.
    auto resolve = [&](std::size_t nodeSide, std::size_t faceSide, bool apply, bool self, bool fresh,
                       double forceScale) {
        bool touching = false;
        Part& nodePart = _parts[nodeSide];
        Part& facePart = _parts[faceSide];
        const double reach = 0.5 * (nodePart.thickness() + facePart.thickness());
        const auto& faces = facePart.triangles();
        if (faces.empty() || nodePart.nodeCount() == 0)
            return false;

        // Nodes that still share a shell are one piece. Self-contact is only
        // for a fragment that has fully separated; a crack-edge node sits on
        // the neighbouring quad and must not be kicked by a full-thickness penalty.
        std::vector<int> component;
        if (self) {
            const int n = static_cast<int>(nodePart.nodeCount());
            component.resize(static_cast<std::size_t>(n));
            for (int nIndex = 0; nIndex < n; ++nIndex)
                component[static_cast<std::size_t>(nIndex)] = nIndex;
            auto link = [&](int a, int b) {
                if (a >= 0 && b >= 0 && a < n && b < n)
                    unite(component, a, b);
            };
            if (!nodePart.quads().empty()) {
                for (const auto& quad : nodePart.quads())
                    for (int k = 1; k < 4; ++k)
                        link(quad[0], quad[k]);
            } else {
                for (const auto& tri : faces)
                    for (int k = 1; k < 3; ++k)
                        link(tri[0], tri[k]);
            }
        }

        const double cell = std::max(0.05, reach * 2.0);
        auto cellKey = [](int x, int y, int z) {
            return (x * 73856093) ^ (y * 19349663) ^ (z * 83492791);
        };
        std::unordered_map<int, std::vector<int>> grid;
        grid.reserve(faces.size() * 2);
        for (int t = 0; t < static_cast<int>(faces.size()); ++t) {
            ChVector3d a = facePart.nodePosition(faces[t][0]);
            ChVector3d b = facePart.nodePosition(faces[t][1]);
            ChVector3d c = facePart.nodePosition(faces[t][2]);
            const double minX = std::min({a.x(), b.x(), c.x()}) - reach;
            const double minY = std::min({a.y(), b.y(), c.y()}) - reach;
            const double minZ = std::min({a.z(), b.z(), c.z()}) - reach;
            const double maxX = std::max({a.x(), b.x(), c.x()}) + reach;
            const double maxY = std::max({a.y(), b.y(), c.y()}) + reach;
            const double maxZ = std::max({a.z(), b.z(), c.z()}) + reach;
            const int x0 = static_cast<int>(std::floor(minX / cell));
            const int y0 = static_cast<int>(std::floor(minY / cell));
            const int z0 = static_cast<int>(std::floor(minZ / cell));
            const int x1 = static_cast<int>(std::floor(maxX / cell));
            const int y1 = static_cast<int>(std::floor(maxY / cell));
            const int z1 = static_cast<int>(std::floor(maxZ / cell));
            for (int x = x0; x <= x1 && x < x0 + 6; ++x)
                for (int y = y0; y <= y1 && y < y0 + 6; ++y)
                    for (int z = z0; z <= z1 && z < z0 + 6; ++z)
                        grid[cellKey(x, y, z)].push_back(t);
        }

        std::vector<int> seen(faces.size(), 0);
        int stamp = 1;
        for (std::size_t i = 0; i < nodePart.nodeCount(); ++i) {
            if (nodePart.nodeRigid(i))
                continue;
            const ChVector3d p = nodePart.nodePosition(i);
            const int cx = static_cast<int>(std::floor(p.x() / cell));
            const int cy = static_cast<int>(std::floor(p.y() / cell));
            const int cz = static_cast<int>(std::floor(p.z() / cell));
            ++stamp;
            if (stamp == 0) {
                std::fill(seen.begin(), seen.end(), 0);
                stamp = 1;
            }
            for (int dx = -1; dx <= 1; ++dx) {
                for (int dy = -1; dy <= 1; ++dy) {
                    for (int dz = -1; dz <= 1; ++dz) {
                        const auto bucket = grid.find(cellKey(cx + dx, cy + dy, cz + dz));
                        if (bucket == grid.end())
                            continue;
                        for (int t : bucket->second) {
                            if (seen[static_cast<std::size_t>(t)] == stamp)
                                continue;
                            seen[static_cast<std::size_t>(t)] = stamp;
                            const auto& tri = faces[static_cast<std::size_t>(t)];
                            if (self) {
                                bool own = false;
                                for (int k = 0; k < 3; ++k)
                                    own = own || tri[k] == static_cast<int>(i);
                                if (own)
                                    continue;
                                for (const auto& twin : nodePart.tearTwins()) {
                                    const bool nodeIsClone = twin.clone == static_cast<int>(i);
                                    const bool nodeIsOriginal = twin.original == static_cast<int>(i);
                                    if (!nodeIsClone && !nodeIsOriginal)
                                        continue;
                                    const int other = nodeIsClone ? twin.original : twin.clone;
                                    for (int k = 0; k < 3; ++k) {
                                        if (tri[k] != other)
                                            continue;
                                        if ((p - nodePart.nodePosition(static_cast<std::size_t>(other)))
                                                .Length() < reach * 2.0)
                                            own = true;
                                    }
                                }
                                if (own)
                                    continue;
                                if (!component.empty() &&
                                    findRoot(component, static_cast<int>(i)) == findRoot(component, tri[0]))
                                    continue;
                            }
            const ChVector3d a = facePart.nodePosition(tri[0]);
            const ChVector3d b = facePart.nodePosition(tri[1]);
            const ChVector3d c = facePart.nodePosition(tri[2]);
                const ChVector3d face = (b - a).Cross(c - a);
                const double faceLen = face.Length();
                if (faceLen < 1.0e-12)
                    continue;
                const ChVector3d nHat = face / faceLen;
                // Face overlap only. A node sitting beside an edge, in the
                // same plane, is not inside the shell, and treating that
                // lateral gap as penetration kicks crack faces apart.
                const double planeGap = (p - a).Dot(nHat);
                if (std::abs(planeGap) >= reach)
                    continue;
                const ChVector3d projected = p - nHat * planeGap;
                double wa, wb, wc;
                barycentric(projected, a, b, c, wa, wb, wc);
                // On the edge means beside the shell, not through it.
                if (wa <= 1.0e-3 || wb <= 1.0e-3 || wc <= 1.0e-3)
                    continue;
                const double distance = std::abs(planeGap);
                touching = true;
                if (!apply)
                    continue;

                ChVector3d normal = planeGap >= 0.0 ? nHat : -nHat;
                if (!self) {
                    ChVector3d& previous = _contactNormal[nodeSide * count + faceSide];
                    if (previous.Length() > 0.5) {
                        if (normal.Dot(previous) < 0.0)
                            normal = -normal;
                        if (normal.Dot(previous) < 0.5)
                            normal = previous;
                    }
                    previous = normal;
                }

                const ChVector3d faceVelocity = facePart.nodeVelocity(tri[0]) * wa +
                                                facePart.nodeVelocity(tri[1]) * wb +
                                                facePart.nodeVelocity(tri[2]) * wc;
                const ChVector3d relative = nodePart.nodeVelocity(i) - faceVelocity;
                const double vn = relative.Dot(normal);

                // The lighter of the two sides sets the stability limit.
                const double effectiveMass =
                    std::min({nodePart.nodeMass()[i], facePart.nodeMass()[tri[0]],
                              facePart.nodeMass()[tri[1]], facePart.nodeMass()[tri[2]]});
                const bool resting = std::abs(vn) <= kRestSpeed;
                const ContactGains gains =
                    resting ? restGains(effectiveMass, dt)
                            : stableGains(environment.contactStiffness, environment.contactDamping,
                                          effectiveMass, dt);
                double overlap = reach - distance;
                if (fresh) {
                    const double closing = std::max(0.0, -vn) * dt + 1.0e-4;
                    overlap = std::min(overlap, closing);
                }
                // Separating contacts give the restitution fraction back.
                // Resting contacts keep their support.
                const double scale = (!resting && vn > 0.0) ? reboundScale : 1.0;
                double fn = forceScale * (gains.stiffness * scale * overlap - gains.damping * vn);
                if (fn <= 0.0)
                    continue;

                ChVector3d force = normal * fn;
                const ChVector3d vt = relative - normal * vn;
                const double vtMagnitude = vt.Length();
                if (vtMagnitude > 1.0e-9) {
                    double ft = mu * fn;
                    if (vtMagnitude < vStick)
                        ft *= vtMagnitude / vStick;
                    force -= vt * (ft / vtMagnitude);
                }

                nodePart.addLoad(i, force);
                facePart.addLoad(static_cast<std::size_t>(tri[0]), force * -wa);
                facePart.addLoad(static_cast<std::size_t>(tri[1]), force * -wb);
                facePart.addLoad(static_cast<std::size_t>(tri[2]), force * -wc);
                        }
                    }
                }
            }
        }
        return touching;
    };

    for (std::size_t i = 0; i < count; ++i) {
        for (std::size_t j = i + 1; j < count; ++j) {
            if (_jointAdjacency[i * count + j])
                continue;

            char& armed = _contactArmed[i * count + j];
            if (!aabbOverlap(lo[i], hi[i], lo[j], hi[j])) {
                armed = 1;
                continue;
            }

            // Two parts that were welded still share a seam the instant the
            // weld lets go. Pushing on that overlap would fire the whole panel
            // across the arena, so a pair only starts colliding once it has
            // been seen fully clear of the other at least one step.
            // armed == 1 means the pair has been clear and this is the first
            // step of contact. Cap that step's overlap so a gap that only
            // flickered open cannot dump a multi-metre-per-second kick.
            // Both directions run at half strength so the penalty is not applied twice.
            const bool apply = armed != 0;
            const bool fresh = armed == 1;
            const bool touchingA = resolve(i, j, apply, false, fresh, 0.5);
            const bool touchingB = resolve(j, i, apply, false, fresh, 0.5);
            if (apply && (touchingA || touchingB))
                armed = 2;
            if (!apply && !touchingA && !touchingB)
                armed = 1;
        }
    }

    for (std::size_t i = 0; i < count; ++i) {
        if (_parts[i].torn())
            resolve(i, i, true, true, false, 1.0);
    }
}

void Assembly::applyUserLoads(const std::vector<ImpactLoad>& loads, double time) {
    for (std::size_t i = 0; i < loads.size() && i < _loadTargets.size(); ++i) {
        const ImpactLoad& load = loads[i];
        if (!load.enabled || load.duration <= 0.0)
            continue;
        if (time < load.startTime || time >= load.startTime + load.duration)
            continue;

        ChVector3d direction = load.direction;
        if (direction.Length() < 1.0e-12)
            continue;
        direction.Normalize();

        for (const auto& target : _loadTargets[i])
            _parts[target.part].addLoad(target.node, direction * (load.force * target.share));
    }
}

bool Assembly::updateJoints(double dt) {
    const double alpha = dt / (kReactionFilterTau + dt);
    bool changed = false;

    for (auto& joint : _joints) {
        for (auto& entry : joint.links) {
            if (entry.broken || (!entry.link && !entry.pointLink))
                continue;

            double force = 0.0;
            double torque = 0.0;
            if (entry.link) {
                const auto reaction = entry.link->GetReaction1();
                force = reaction.force.Length();
                torque = reaction.torque.Length();
            } else if (entry.pointLink) {
                force = entry.pointLink->GetReactionOnNode1().Length();
            }
            if (!std::isfinite(force))
                force = entry.breakForce * 10.0;
            if (!std::isfinite(torque))
                torque = entry.breakTorque * 10.0;

            entry.filteredForce += alpha * (force - entry.filteredForce);
            entry.filteredTorque += alpha * (torque - entry.filteredTorque);

            double utilisation =
                entry.breakForce > 0.0 ? entry.filteredForce / entry.breakForce : 0.0;
            if (entry.locksRotation && entry.breakTorque > 0.0) {
                utilisation =
                    std::max(utilisation, entry.filteredTorque / entry.breakTorque);
            }
            entry.utilisation = utilisation;

            if (utilisation > 1.0)
                ++entry.overloadSteps;
            else
                entry.overloadSteps = 0;

            if (entry.overloadSteps >= kOverloadHoldSteps) {
                if (entry.link)
                    entry.link->SetBroken(true);
                if (entry.pointLink)
                    entry.pointLink->SetBroken(true);
                entry.broken = true;
                changed = true;
            }
        }
    }

    if (changed)
        refreshAdjacency();
    return changed;
}

void Assembly::reindexNodes() {
    _nodeOffset.resize(_parts.size());
    std::size_t offset = 0;
    for (std::size_t i = 0; i < _parts.size(); ++i) {
        _nodeOffset[i] = offset;
        offset += _parts[i].nodeCount();
    }
}

bool Assembly::tearOverstretched() {
    for (auto& part : _parts)
        part.clearLinkedNodes();
    for (const auto& joint : _joints) {
        for (const auto& link : joint.links) {
            if (link.broken)
                continue;
            if (link.partA < _parts.size())
                _parts[link.partA].markLinkedNode(static_cast<int>(link.nodeA));
            if (link.partB < _parts.size())
                _parts[link.partB].markLinkedNode(static_cast<int>(link.nodeB));
        }
    }
    for (const auto& link : _deferred) {
        if (link.partA < _parts.size())
            _parts[link.partA].markLinkedNode(static_cast<int>(link.nodeA));
        if (link.partB < _parts.size())
            _parts[link.partB].markLinkedNode(static_cast<int>(link.nodeB));
    }

    bool changed = false;
    for (auto& part : _parts)
        changed = part.tearFailed() || changed;
    if (changed)
        reindexNodes();
    return changed;
}

bool Assembly::respan() {
    if (!_system || _deferred.empty())
        return false;
    reindexNodes();
    std::size_t total = 0;
    for (const auto& part : _parts)
        total += part.nodeCount();
    if (total == 0)
        return false;

    std::vector<int> translationSet(total);
    std::vector<int> rotationSet(total);
    for (std::size_t i = 0; i < total; ++i) {
        translationSet[i] = static_cast<int>(i);
        rotationSet[i] = static_cast<int>(i);
    }
    for (const auto& joint : _joints) {
        for (const auto& entry : joint.links) {
            if (entry.broken)
                continue;
            const int ga = static_cast<int>(globalNode(entry.partA, entry.nodeA));
            const int gb = static_cast<int>(globalNode(entry.partB, entry.nodeB));
            if (ga < 0 || gb < 0 || static_cast<std::size_t>(ga) >= total ||
                static_cast<std::size_t>(gb) >= total)
                continue;
            unite(translationSet, ga, gb);
            if (entry.locksRotation)
                unite(rotationSet, ga, gb);
        }
    }

    bool added = false;
    std::vector<DeferredLink> still;
    for (const auto& skipped : _deferred) {
        if (skipped.partA >= _parts.size() || skipped.partB >= _parts.size())
            continue;
        if (skipped.nodeA >= _parts[skipped.partA].nodeCount() ||
            skipped.nodeB >= _parts[skipped.partB].nodeCount()) {
            still.push_back(skipped);
            continue;
        }
        const int ga = static_cast<int>(globalNode(skipped.partA, skipped.nodeA));
        const int gb = static_cast<int>(globalNode(skipped.partB, skipped.nodeB));
        if (!unite(translationSet, ga, gb)) {
            still.push_back(skipped);
            continue;
        }
        const bool lockRotation = skipped.wantRotation && unite(rotationSet, ga, gb);
        Joint* joint = nullptr;
        for (auto& existing : _joints) {
            if (existing.partA == skipped.partA && existing.partB == skipped.partB &&
                existing.kind == skipped.kind) {
                joint = &existing;
                break;
            }
        }
        if (!joint) {
            Joint created;
            created.kind = skipped.kind;
            created.label = skipped.label.empty() ? "Respanned" : skipped.label;
            created.partA = skipped.partA;
            created.partB = skipped.partB;
            _joints.push_back(std::move(created));
            joint = &_joints.back();
        }
        createLink(*_system, *joint, skipped.partA, skipped.nodeA, skipped.partB, skipped.nodeB,
                   skipped.breakForce, skipped.breakTorque, skipped.wantRotation, lockRotation);
        added = true;
    }
    _deferred.swap(still);
    if (added)
        refreshAdjacency();
    return added;
}

void Assembly::refreshAdjacency() {
    const std::size_t count = _parts.size();
    _jointAdjacency.assign(count * count, 0);
    for (const auto& joint : _joints) {
        if (joint.severed())
            continue;
        _jointAdjacency[joint.partA * count + joint.partB] = 1;
        _jointAdjacency[joint.partB * count + joint.partA] = 1;
    }
}

void Assembly::writeLayout(RecordingLayout& layout) const {
    layout.parts.clear();
    layout.joints.clear();
    for (const auto& part : _parts) {
        PartLayout entry;
        entry.name = part.name();
        entry.materialName = part.material().name;
        entry.thickness = static_cast<float>(part.thickness());
        entry.nodeCount = static_cast<int>(part.nodeCount());
        entry.quads = part.quads();
        entry.triangles = part.triangles();
        layout.parts.push_back(std::move(entry));
    }
    for (const auto& joint : _joints)
        layout.joints.push_back({joint.label, static_cast<int>(joint.links.size())});
}

void Assembly::writeFrame(Frame& frame) const {
    frame.parts.resize(_parts.size());
    frame.maxStress = 0.0f;
    frame.brokenLinks = 0;

    for (std::size_t i = 0; i < _parts.size(); ++i) {
        const Part& part = _parts[i];
        PartFrame& out = frame.parts[i];
        out.xyz.resize(part.nodeCount() * 3);
        for (std::size_t n = 0; n < part.nodeCount(); ++n) {
            const auto p = part.nodePosition(n);
            out.setPosition(n, p.x(), p.y(), p.z());
        }
        if (part.torn()) {
            out.quads = part.quads();
            out.triangles = part.triangles();
        }
        part.evaluateStress(out.stress);
        part.evaluatePlasticStrain(out.plasticStrain);
        for (float s : out.stress)
            frame.maxStress = std::max(frame.maxStress, s);
    }

    frame.joints.resize(_joints.size());
    for (std::size_t i = 0; i < _joints.size(); ++i) {
        const Joint& joint = _joints[i];
        JointFrame& out = frame.joints[i];
        out.utilisation.resize(joint.links.size());
        out.broken.resize(joint.links.size());
        for (std::size_t l = 0; l < joint.links.size(); ++l) {
            out.utilisation[l] = static_cast<float>(joint.links[l].utilisation);
            out.broken[l] = joint.links[l].broken ? 1u : 0u;
            if (joint.links[l].broken)
                ++frame.brokenLinks;
        }
    }
}

bool Assembly::finite() const {
    for (const auto& part : _parts) {
        if (!part.finite())
            return false;
    }
    return true;
}

chrono::ChVector3d Assembly::center() const {
    std::vector<std::size_t> all(_parts.size());
    for (std::size_t i = 0; i < all.size(); ++i)
        all[i] = i;
    return centerOf(all);
}

chrono::ChVector3d Assembly::centerOf(const std::vector<std::size_t>& partIndices) const {
    ChVector3d sum{0, 0, 0};
    std::size_t count = 0;
    for (std::size_t index : partIndices) {
        if (index >= _parts.size())
            continue;
        for (std::size_t n = 0; n < _parts[index].nodeCount(); ++n) {
            sum += _parts[index].nodePosition(n);
            ++count;
        }
    }
    if (count == 0)
        return {0, 0, 0};
    return sum / static_cast<double>(count);
}

double Assembly::mass() const {
    double total = 0.0;
    for (const auto& part : _parts)
        total += part.mass();
    return total;
}

std::vector<std::size_t> Assembly::connectedParts(std::size_t seed) const {
    std::vector<std::size_t> component;
    if (seed >= _parts.size())
        return component;

    std::vector<char> seen(_parts.size(), 0);
    std::vector<std::size_t> stack{seed};
    seen[seed] = 1;
    while (!stack.empty()) {
        const std::size_t i = stack.back();
        stack.pop_back();
        component.push_back(i);
        for (const auto& joint : _joints) {
            if (joint.severed())
                continue;
            std::size_t other = _parts.size();
            if (joint.partA == i)
                other = joint.partB;
            else if (joint.partB == i)
                other = joint.partA;
            if (other >= _parts.size() || seen[other])
                continue;
            seen[other] = 1;
            stack.push_back(other);
        }
    }
    return component;
}

int Assembly::brokenLinkCount() const {
    int total = 0;
    for (const auto& joint : _joints)
        total += joint.brokenCount();
    return total;
}

}  // namespace rcs
