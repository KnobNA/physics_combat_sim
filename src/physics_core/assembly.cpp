#include "physics_core/assembly.hpp"

#include "chrono/core/ChMatrix33.h"
#include "chrono/core/ChRotation.h"

#include <algorithm>
#include <cmath>
#include <limits>

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

ContactGains stableGains(double stiffness, double damping, double nodeMass, double dt) {
    if (nodeMass <= 0.0 || dt <= 0.0)
        return {0.0, 0.0};
    const double k = std::min(stiffness, nodeMass / (dt * dt));
    const double criticalDamping = 2.0 * std::sqrt(k * nodeMass);
    const double c = std::min({damping, criticalDamping, 0.5 * nodeMass / dt});
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
    using namespace chrono;

    const int ga = static_cast<int>(globalNode(partA, nodeA));
    const int gb = static_cast<int>(globalNode(partB, nodeB));

    if (!unite(translationSet, ga, gb))
        return;  // already tied together through another link

    const bool lockRotation = wantRotationLock && unite(rotationSet, ga, gb);

    const auto& a = _parts[partA].nodes()[nodeA];
    const auto& b = _parts[partB].nodes()[nodeB];

    auto link = chrono_types::make_shared<ChLinkMateGeneric>(
        true, true, true, lockRotation, lockRotation, lockRotation);
    link->Initialize(a, b, ChFrame<>(0.5 * (a->GetPos() + b->GetPos()), a->GetRot()));
    system.Add(link);

    JointLink entry;
    entry.link = link;
    entry.partA = partA;
    entry.nodeA = nodeA;
    entry.partB = partB;
    entry.nodeB = nodeB;
    entry.breakForce = breakForce;
    entry.breakTorque = breakTorque;
    entry.locksRotation = lockRotation;
    joint.links.push_back(entry);
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
        const ChVector3d pa = a.nodes()[ia]->GetPos();
        if (!insideRegion(pa))
            continue;
        std::size_t best = b.nodes().size();
        double bestDistance = definition.tolerance;
        for (int ib : b.boundaryNodes()) {
            const ChVector3d pb = b.nodes()[ib]->GetPos();
            if (!insideRegion(pb))
                continue;
            const double distance = (pa - pb).Length();
            if (distance <= bestDistance) {
                bestDistance = distance;
                best = static_cast<std::size_t>(ib);
            }
        }
        if (best < b.nodes().size())
            pairs.push_back({static_cast<std::size_t>(ia), best,
                             0.5 * (pa + b.nodes()[best]->GetPos())});
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
        std::size_t best = part.nodes().size();
        double bestDistance = definition.clampRadius;
        for (std::size_t i = 0; i < part.nodes().size(); ++i) {
            const double distance = (part.nodes()[i]->GetPos() - definition.position).Length();
            if (distance <= bestDistance) {
                bestDistance = distance;
                best = i;
            }
        }
        return best;
    };

    const std::size_t nodeA = nearestNode(_parts[definition.partA]);
    const std::size_t nodeB = nearestNode(_parts[definition.partB]);
    if (nodeA >= _parts[definition.partA].nodes().size() ||
        nodeB >= _parts[definition.partB].nodes().size())
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
            if (!entry.link)
                continue;
            const auto& a = _parts[entry.partA].nodes()[entry.nodeA];
            const auto& b = _parts[entry.partB].nodes()[entry.nodeB];
            entry.link->Initialize(
                a, b, ChFrame<>(0.5 * (a->GetPos() + b->GetPos()), a->GetRot()));
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
        for (std::size_t n = 0; n < part.nodes().size(); ++n) {
            if ((part.nodes()[n]->GetPos() - load.point).Length() > load.radius)
                continue;
            targets.push_back({load.partIndex, n, part.nodeMass()[n]});
            totalMass += part.nodeMass()[n];
        }

        // Nothing inside the patch: fall back to the single closest node so a
        // mis-placed load still does something visible.
        if (targets.empty()) {
            std::size_t best = 0;
            double bestDistance = std::numeric_limits<double>::max();
            for (std::size_t n = 0; n < part.nodes().size(); ++n) {
                const double distance = (part.nodes()[n]->GetPos() - load.point).Length();
                if (distance < bestDistance) {
                    bestDistance = distance;
                    best = n;
                }
            }
            if (!part.nodes().empty())
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
    applyUserLoads(loads, time);

    for (auto& part : _parts)
        part.commitLoads();
}

void Assembly::applyGravity(const chrono::ChVector3d& gravity) {
    for (auto& part : _parts) {
        for (std::size_t i = 0; i < part.nodes().size(); ++i)
            part.addLoad(i, part.nodeMass()[i] * gravity);
    }
}

void Assembly::applyAerodynamics(const Environment& environment) {
    const double factor = 0.5 * environment.airDensity * environment.dragCoefficient;
    if (factor <= 0.0)
        return;

    for (auto& part : _parts) {
        const auto& nodes = part.nodes();
        for (const auto& quad : part.quads()) {
            const ChVector3d p0 = nodes[quad[0]]->GetPos();
            const ChVector3d p1 = nodes[quad[1]]->GetPos();
            const ChVector3d p2 = nodes[quad[2]]->GetPos();
            const ChVector3d p3 = nodes[quad[3]]->GetPos();

            // Area vector of the (possibly warped) quad from its diagonals.
            ChVector3d areaVector = 0.5 * (p2 - p0).Cross(p3 - p1);
            const double area = areaVector.Length();
            if (area < 1.0e-12)
                continue;
            const ChVector3d normal = areaVector / area;

            const ChVector3d velocity =
                0.25 * (nodes[quad[0]]->GetPosDt() + nodes[quad[1]]->GetPosDt() +
                        nodes[quad[2]]->GetPosDt() + nodes[quad[3]]->GetPosDt());

            // Pressure drag only responds to the through-panel component; a
            // panel sliding edge-on sees essentially nothing.
            const double vn = velocity.Dot(normal);
            const ChVector3d force = normal * (-factor * area * std::abs(vn) * vn);
            for (int i = 0; i < 4; ++i)
                part.addLoad(quad[i], force * 0.25);
        }
    }
}

void Assembly::applyFloorContact(const Environment& environment, double dt) {
    const double mu = environment.contactFriction;
    const double vStick = std::max(1.0e-4, environment.frictionRegularization);
    const double restitution = std::clamp(environment.floorRestitution, 0.0, 1.0);
    // Energy handed back scales with e². Concrete (e ≈ 0.15) returns about 2%
    // of the impact energy; the panels can still flex, but the floor does not
    // throw the box back up.
    const double reboundScale = restitution * restitution;

    if (_floorPeak.size() != _parts.size())
        _floorPeak.assign(_parts.size(), {});

    for (std::size_t p = 0; p < _parts.size(); ++p) {
        auto& part = _parts[p];
        auto& peak = _floorPeak[p];
        if (peak.size() != part.nodes().size())
            peak.assign(part.nodes().size(), 0.0);

        const double radius = 0.5 * part.thickness();
        for (std::size_t i = 0; i < part.nodes().size(); ++i) {
            const auto& node = part.nodes()[i];
            const double penetration = (environment.floorHeight + radius) - node->GetPos().y();
            if (penetration <= 0.0) {
                peak[i] = 0.0;
                continue;
            }

            const ContactGains gains = stableGains(environment.contactStiffness,
                                                   environment.contactDamping,
                                                   part.nodeMass()[i], dt);
            const ChVector3d v = node->GetPosDt();
            // Full stiffness while the node is digging in or at rest, so the
            // floor holds the box up like concrete. Once the node is actually
            // rising out of a deeper contact, only the restitution fraction of
            // that spring is left, and the stored depth is kept until the node
            // starts descending again. Forgetting it on the way down is what
            // lets a resting box sit on a small elastic sag instead of a
            // leftover impact dent.
            if (penetration > peak[i])
                peak[i] = penetration;
            const bool rebounding = v.y() > 0.0 && penetration < peak[i];
            const double scale = rebounding ? reboundScale : 1.0;
            if (v.y() <= 0.0)
                peak[i] = penetration;
            double fn = gains.stiffness * scale * penetration - gains.damping * v.y();
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
    if (count < 2)
        return;
    if (_contactArmed.size() != count * count)
        _contactArmed.assign(count * count, 0);

    const double mu = environment.contactFriction;
    const double vStick = std::max(1.0e-4, environment.frictionRegularization);

    std::vector<ChVector3d> lo(count), hi(count);
    for (std::size_t i = 0; i < count; ++i)
        _parts[i].bounds(lo[i], hi[i]);

    // Node-against-current-triangle proximity: the contact surface is whatever
    // shape the mesh has right now, deformed or not. Returns true when the two
    // parts overlap at all, which is what arming uses.
    auto resolve = [&](std::size_t nodeSide, std::size_t faceSide, bool apply) {
        bool touching = false;
        Part& nodePart = _parts[nodeSide];
        Part& facePart = _parts[faceSide];
        const double reach = 0.5 * (nodePart.thickness() + facePart.thickness());
        const auto& faceNodes = facePart.nodes();

        for (std::size_t i = 0; i < nodePart.nodes().size(); ++i) {
            const ChVector3d p = nodePart.nodes()[i]->GetPos();
            for (const auto& tri : facePart.triangles()) {
                const ChVector3d a = faceNodes[tri[0]]->GetPos();
                const ChVector3d b = faceNodes[tri[1]]->GetPos();
                const ChVector3d c = faceNodes[tri[2]]->GetPos();
                const ChVector3d closest = closestPointOnTriangle(p, a, b, c);
                ChVector3d delta = p - closest;
                const double distance = delta.Length();
                if (distance >= reach)
                    continue;
                touching = true;
                if (!apply)
                    continue;

                ChVector3d normal;
                if (distance > 1.0e-9) {
                    normal = delta / distance;
                } else {
                    normal = (b - a).Cross(c - a);
                    if (normal.Length() < 1.0e-12)
                        continue;
                    normal.Normalize();
                }

                double wa, wb, wc;
                barycentric(closest, a, b, c, wa, wb, wc);
                const ChVector3d faceVelocity = faceNodes[tri[0]]->GetPosDt() * wa +
                                                faceNodes[tri[1]]->GetPosDt() * wb +
                                                faceNodes[tri[2]]->GetPosDt() * wc;
                const ChVector3d relative = nodePart.nodes()[i]->GetPosDt() - faceVelocity;
                const double vn = relative.Dot(normal);

                // The lighter of the two sides sets the stability limit.
                const double effectiveMass =
                    std::min({nodePart.nodeMass()[i], facePart.nodeMass()[tri[0]],
                              facePart.nodeMass()[tri[1]], facePart.nodeMass()[tri[2]]});
                const ContactGains gains = stableGains(
                    environment.contactStiffness, environment.contactDamping, effectiveMass, dt);

                double fn = gains.stiffness * (reach - distance) - gains.damping * vn;
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
                facePart.addLoad(tri[0], force * -wa);
                facePart.addLoad(tri[1], force * -wb);
                facePart.addLoad(tri[2], force * -wc);
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
            const bool apply = armed != 0;
            const bool touchingA = resolve(i, j, apply);
            const bool touchingB = resolve(j, i, apply);
            if (!apply && !touchingA && !touchingB)
                armed = 1;
        }
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
            if (entry.broken || !entry.link)
                continue;

            const auto reaction = entry.link->GetReaction1();
            double force = reaction.force.Length();
            double torque = reaction.torque.Length();
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
                entry.link->SetBroken(true);
                entry.broken = true;
                changed = true;
            }
        }
    }

    if (changed)
        refreshAdjacency();
    return changed;
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
        entry.nodeCount = static_cast<int>(part.nodes().size());
        entry.quads = part.quads();
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
        out.positions.resize(part.nodes().size());
        for (std::size_t n = 0; n < part.nodes().size(); ++n) {
            const auto& p = part.nodes()[n]->GetPos();
            out.positions[n] = {static_cast<float>(p.x()), static_cast<float>(p.y()),
                                static_cast<float>(p.z())};
        }
        part.evaluateStress(out.stress);
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
        for (const auto& node : _parts[index].nodes()) {
            sum += node->GetPos();
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
