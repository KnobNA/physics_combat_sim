#include "physics_core/assembly.hpp"
#include "physics_core/physics_units.hpp"

#include "chrono/core/ChMatrix33.h"
#include "chrono/core/ChRotation.h"
#include "chrono/fea/ChLinkNodeFrame.h"

#include <algorithm>
#include <array>
#include <cmath>
#include <limits>
#include <unordered_map>

namespace rcs {
namespace {

using chrono::ChVector3d;

/// A link is failed only after it has been over its allowable for this many
/// consecutive steps, so a single solver spike cannot unzip a weld.
constexpr int kOverloadHoldSteps = 8;
/// Time constant of the reaction low-pass filter, seconds.
constexpr double kReactionFilterTau = 8.0e-3;

ContactGains nodeFaceGains(const Part& nodePart, const Part& facePart, std::size_t node, int a, int b,
                           int c, double patchArea) {
    const auto& nm = nodePart.nodeMass();
    const auto& fm = facePart.nodeMass();
    const double mNode = node < nm.size() ? nm[node] : 0.0;
    const double m0 = static_cast<std::size_t>(a) < fm.size() ? fm[static_cast<std::size_t>(a)] : 0.0;
    const double m1 = static_cast<std::size_t>(b) < fm.size() ? fm[static_cast<std::size_t>(b)] : 0.0;
    const double m2 = static_cast<std::size_t>(c) < fm.size() ? fm[static_cast<std::size_t>(c)] : 0.0;
    const double effectiveMass = std::min({mNode, m0, m1, m2});
    const double E = std::min(std::max(1.0e5, nodePart.material().youngsModulus),
                              std::max(1.0e5, facePart.material().youngsModulus));
    return continuumGains(E, patchArea, effectiveMass);
}

bool nearPoints(const std::vector<ChVector3d>& pts, const ChVector3d& p, double radius) {
    if (pts.empty() || radius <= 0.0)
        return false;
    const double r2 = radius * radius;
    for (const ChVector3d& q : pts) {
        if ((p - q).Length2() <= r2)
            return true;
    }
    return false;
}

int hashCell(int x, int y, int z) {
    return (x * 73856093) ^ (y * 19349663) ^ (z * 83492791);
}

bool tetContainsPoint(const Part& part, const std::array<int, 4>& tet, const ChVector3d& p) {
    return rcs::tetContainsPoint(part.nodePosition(static_cast<std::size_t>(tet[0])),
                                 part.nodePosition(static_cast<std::size_t>(tet[1])),
                                 part.nodePosition(static_cast<std::size_t>(tet[2])),
                                 part.nodePosition(static_cast<std::size_t>(tet[3])), p);
}

struct VolumeGrid {
    double cell = 0.04;
    const Part* part = nullptr;
    std::unordered_map<int, std::vector<int>> buckets;

    void build(const Part& host) {
        part = &host;
        buckets.clear();
        if (!host.hasVolume())
            return;
        ChVector3d lo, hi;
        host.bounds(lo, hi);
        const double span = std::max({hi.x() - lo.x(), hi.y() - lo.y(), hi.z() - lo.z(), 1.0e-3});
        cell = std::max(0.02, span / 8.0);
        const auto& tets = host.tets();
        for (int t = 0; t < static_cast<int>(tets.size()); ++t) {
            if (!host.tetLive(static_cast<std::size_t>(t)))
                continue;
            const auto& tet = tets[static_cast<std::size_t>(t)];
            ChVector3d a = host.nodePosition(static_cast<std::size_t>(tet[0]));
            ChVector3d b = host.nodePosition(static_cast<std::size_t>(tet[1]));
            ChVector3d c = host.nodePosition(static_cast<std::size_t>(tet[2]));
            ChVector3d d = host.nodePosition(static_cast<std::size_t>(tet[3]));
            const double minX = std::min({a.x(), b.x(), c.x(), d.x()});
            const double minY = std::min({a.y(), b.y(), c.y(), d.y()});
            const double minZ = std::min({a.z(), b.z(), c.z(), d.z()});
            const double maxX = std::max({a.x(), b.x(), c.x(), d.x()});
            const double maxY = std::max({a.y(), b.y(), c.y(), d.y()});
            const double maxZ = std::max({a.z(), b.z(), c.z(), d.z()});
            const int x0 = static_cast<int>(std::floor(minX / cell));
            const int y0 = static_cast<int>(std::floor(minY / cell));
            const int z0 = static_cast<int>(std::floor(minZ / cell));
            const int x1 = static_cast<int>(std::floor(maxX / cell));
            const int y1 = static_cast<int>(std::floor(maxY / cell));
            const int z1 = static_cast<int>(std::floor(maxZ / cell));
            for (int x = x0; x <= x1 && x < x0 + 8; ++x)
                for (int y = y0; y <= y1 && y < y0 + 8; ++y)
                    for (int z = z0; z <= z1 && z < z0 + 8; ++z)
                        buckets[hashCell(x, y, z)].push_back(t);
        }
    }

    bool contains(const ChVector3d& p) const {
        if (!part || buckets.empty())
            return false;
        const int cx = static_cast<int>(std::floor(p.x() / cell));
        const int cy = static_cast<int>(std::floor(p.y() / cell));
        const int cz = static_cast<int>(std::floor(p.z() / cell));
        const auto bucket = buckets.find(hashCell(cx, cy, cz));
        if (bucket == buckets.end())
            return false;
        const auto& tets = part->tets();
        for (int t : bucket->second) {
            if (t < 0 || static_cast<std::size_t>(t) >= tets.size())
                continue;
            if (tetContainsPoint(*part, tets[static_cast<std::size_t>(t)], p))
                return true;
        }
        return false;
    }
};

bool pointInLiveTets(const Part& part, const ChVector3d& p) {
    if (!part.hasVolume())
        return false;
    const auto& tets = part.tets();
    for (std::size_t t = 0; t < tets.size(); ++t) {
        if (!part.tetLive(t))
            continue;
        if (tetContainsPoint(part, tets[t], p))
            return true;
    }
    return false;
}

bool nodeIsMolten(const Part& part, std::size_t index) {
    return part.nodeTemperature(index) >= part.material().meltTemperature;
}

bool nodeCorrectable(const Part& part, std::size_t index) {
    if (part.nodeRigid(index))
        return false;
    if (part.hasVolume() && !part.nodeOnLiveTet(index))
        return false;
    if (nodeIsMolten(part, index))
        return false;
    return true;
}

bool segmentHitsTriangle(const ChVector3d& p0,
                         const ChVector3d& p1,
                         const ChVector3d& a,
                         const ChVector3d& b,
                         const ChVector3d& c,
                         double& tHit,
                         ChVector3d& hit,
                         ChVector3d& normal) {
    const ChVector3d dir = p1 - p0;
    const ChVector3d ab = b - a;
    const ChVector3d ac = c - a;
    const ChVector3d nrm = ab.Cross(ac);
    const double denom = nrm.Dot(dir);
    if (std::abs(denom) < 1.0e-16)
        return false;
    tHit = nrm.Dot(a - p0) / denom;
    if (tHit < -1.0e-6 || tHit > 1.0 + 1.0e-6)
        return false;
    hit = p0 + dir * tHit;
    double wa, wb, wc;
    barycentric(hit, a, b, c, wa, wb, wc);
    if (wa < -1.0e-3 || wb < -1.0e-3 || wc < -1.0e-3)
        return false;
    normal = nrm;
    const double len = normal.Length();
    if (len < 1.0e-12)
        return false;
    normal /= len;
    return true;
}

constexpr char kSeamNone = 0;
constexpr char kSeamIntact = 1;
constexpr char kSeamDamaged = 2;

ChVector3d aabbSeparate(const ChVector3d& glo,
                        const ChVector3d& ghi,
                        const ChVector3d& hlo,
                        const ChVector3d& hhi,
                        double extra) {
    const double ox = std::min(ghi.x(), hhi.x()) - std::max(glo.x(), hlo.x());
    const double oy = std::min(ghi.y(), hhi.y()) - std::max(glo.y(), hlo.y());
    const double oz = std::min(ghi.z(), hhi.z()) - std::max(glo.z(), hlo.z());
    if (ox <= 0.0 || oy <= 0.0 || oz <= 0.0)
        return {0, 0, 0};
    const ChVector3d gc = 0.5 * (glo + ghi);
    const ChVector3d hc = 0.5 * (hlo + hhi);
    extra = std::max(extra, 1.0e-4);
    if (ox <= oy && ox <= oz) {
        const double s = gc.x() >= hc.x() ? 1.0 : -1.0;
        return {s * (ox + extra), 0, 0};
    }
    if (oy <= oz) {
        const double s = gc.y() >= hc.y() ? 1.0 : -1.0;
        return {0, s * (oy + extra), 0};
    }
    const double s = gc.z() >= hc.z() ? 1.0 : -1.0;
    return {0, 0, s * (oz + extra)};
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
    _continuity.clear();
    _deferred.clear();
    _system = &system;
    _nodeOffset.clear();
    _loadTargets.clear();
    _contactArmed.clear();
    _contacts.clear();
    _contactShape.clear();
    _contactContainer.reset();
    _contactWorld.reset();

    _parts.resize(definition.parts.size());
    std::size_t offset = 0;
    for (std::size_t i = 0; i < definition.parts.size(); ++i) {
        _nodeOffset.push_back(offset);
        if (!definition.parts[i].simulate)
            continue;
        _parts[i].build(definition.parts[i]);
        offset += _parts[i].nodeCount();
        if (_parts[i].mesh())
            system.Add(_parts[i].mesh());
        for (const auto& tie : _parts[i].shellSolidTies()) {
            const int shell = tie.first;
            const int solid = tie.second;
            if (shell < 0 || solid < 0)
                continue;
            if (static_cast<std::size_t>(shell) >= _parts[i].xyzNodes().size() ||
                static_cast<std::size_t>(solid) >= _parts[i].xyzNodes().size())
                continue;
            if (shell == solid)
                continue;
            auto link = chrono_types::make_shared<chrono::fea::ChLinkNodeNode>();
            link->Initialize(_parts[i].xyzNodes()[static_cast<std::size_t>(shell)],
                             _parts[i].xyzNodes()[static_cast<std::size_t>(solid)]);
            system.Add(link);
            _continuity.push_back(link);
        }
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

    _weldDefs = definition.welds;
    std::vector<char> snapped(offset, 0);
    for (const auto& weld : definition.welds)
        buildWeld(system, weld, translationSet, rotationSet, snapped);
    for (const auto& screw : definition.screws)
        buildScrew(system, screw, translationSet, rotationSet);

    // Editor transforms are authored in world metres. A live COM pivot would
    // slide every other part whenever one cube moves, so the drop pose is an
    // offset of the origin, not of the current centre of mass.
    _authoredCenter = ChVector3d(0, 0, 0);
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

    const bool rotA = _parts[partA].rotational();
    const bool rotB = _parts[partB].rotational();
    if (!rotA && !rotB) {
        auto link = chrono_types::make_shared<fea::ChLinkNodeNode>();
        link->Initialize(_parts[partA].xyzNodes()[nodeA], _parts[partB].xyzNodes()[nodeB]);
        const ChVector3d sep =
            _parts[partA].nodePosition(nodeA) - _parts[partB].nodePosition(nodeB);
        // A shell midsurface does not sit on the mating face. Collapsing that
        // gap, or rewriting the rest pose to hide it, strains the mesh.
        const bool shellSeam = !_parts[partA].hasVolume() || !_parts[partB].hasVolume();
        if (shellSeam || sep.Length() > 1.0e-6) {
            link->SetRestOffset(sep);
            entry.offsetLink = link;
        } else {
            entry.pointLink = link;
        }
        system.Add(link);
    } else if (rotA != rotB) {
        const std::size_t xyzPart = rotA ? partB : partA;
        const std::size_t xyzNode = rotA ? nodeB : nodeA;
        const std::size_t framePart = rotA ? partA : partB;
        const std::size_t frameNode = rotA ? nodeA : nodeB;
        auto link = chrono_types::make_shared<fea::ChLinkNodeFrame>();
        link->Initialize(_parts[xyzPart].xyzNodes()[xyzNode],
                         std::static_pointer_cast<ChBodyFrame>(_parts[framePart].nodes()[frameNode]));
        system.Add(link);
        entry.frameLink = link;
    } else {
        const auto& a = _parts[partA].nodes()[nodeA];
        const auto& b = _parts[partB].nodes()[nodeB];
        // A weld bead is rigid in rotation too. One row of position-only
        // ties along an edge leaves a hinge about that edge.
        auto link = chrono_types::make_shared<ChLinkMateGeneric>(true, true, true, lockRotation,
                                                                 lockRotation, lockRotation);
        link->Initialize(a, b, ChFrame<>(0.5 * (a->GetPos() + b->GetPos()), a->GetRot()));
        system.Add(link);
        entry.link = link;
    }
    (void)wantRotationLock;
    joint.links.push_back(entry);
}

bool Assembly::addLink(chrono::ChSystem& system,
                       Joint& joint,
                       std::size_t partA,
                       std::size_t nodeA,
                       std::size_t partB,
                       std::size_t nodeB,
                       double breakForce,
                       double breakTorque,
                       bool wantRotationLock,
                       std::vector<int>& translationSet,
                       std::vector<int>& rotationSet,
                       double breakPeel,
                       double criticalOpening,
                       const chrono::ChVector3d& peelAxis)
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
        skipped.breakPeel = breakPeel;
        skipped.criticalOpening = criticalOpening;
        skipped.peelAxis = peelAxis;
        skipped.wantRotation = wantRotationLock;
        skipped.label = joint.label;
        skipped.kind = joint.kind;
        _deferred.push_back(skipped);
        return false;
    }

    // A rotation lock between nodes already rotationally joined through
    // another path is redundant and leaves the constraint block singular.
    const bool lockRotation = wantRotationLock && unite(rotationSet, ga, gb);
    createLink(system, joint, partA, nodeA, partB, nodeB, breakForce, breakTorque, wantRotationLock,
               lockRotation);
    if (!joint.links.empty()) {
        JointLink& entry = joint.links.back();
        entry.breakPeel = breakPeel;
        entry.criticalOpening = criticalOpening;
        entry.peelAxis = peelAxis;
        const chrono::ChVector3d sep = _parts[partB].nodePosition(nodeB) - _parts[partA].nodePosition(nodeA);
        entry.restGap = sep.Length();
    }
    return true;
}

void Assembly::buildWeld(chrono::ChSystem& system,
                         const WeldDefinition& definition,
                         std::vector<int>& translationSet,
                         std::vector<int>& rotationSet,
                         std::vector<char>& snapped)
{
    if (definition.partA >= _parts.size() || definition.partB >= _parts.size() ||
        definition.partA == definition.partB)
        return;

    Part& a = _parts[definition.partA];
    Part& b = _parts[definition.partB];

    auto insideRegion = [&](const ChVector3d& p) {
        if (!definition.useRegion)
            return true;
        return p.x() >= definition.regionMin.x() && p.x() <= definition.regionMax.x() &&
               p.y() >= definition.regionMin.y() && p.y() <= definition.regionMax.y() &&
               p.z() >= definition.regionMin.z() && p.z() <= definition.regionMax.z();
    };

    const bool haveBead = definition.bead.size() >= 2;
    // Independent tet fills do not share Steiner points, so a bead weld has to
    // pair nodes that sit on the same seam. Capture stays under one wall so
    // the far skin of a thin panel is not pulled onto the bead. A midplane
    // shell (or a patch extruded ±t/2) lives half a thickness from its face.
    auto wallOf = [](const Part& part) {
        const double authored = std::max(part.thickness(), 1.0e-6);
        if (part.nodeCount() == 0)
            return authored;
        ChVector3d lo = part.nodePosition(0);
        ChVector3d hi = lo;
        for (std::size_t i = 1; i < part.nodeCount(); ++i) {
            const ChVector3d p = part.nodePosition(i);
            lo = {std::min(lo.x(), p.x()), std::min(lo.y(), p.y()), std::min(lo.z(), p.z())};
            hi = {std::max(hi.x(), p.x()), std::max(hi.y(), p.y()), std::max(hi.z(), p.z())};
        }
        const ChVector3d span = hi - lo;
        const double shortest =
            std::min({std::abs(span.x()), std::abs(span.y()), std::abs(span.z())});
        // A mid-surface shell is geometrically flat, so the AABB is paper-thin.
        // Using that as the wall made weld snap a micron and the coincident
        // link then yanked the corner 5 mm on the first HHT attempt.
        if (shortest < 0.25 * authored)
            return authored;
        return std::min(authored, std::max(shortest, 1.0e-6));
    };
    const double thick = std::min(wallOf(a), wallOf(b));
    const double midplane = std::max(a.rotational() ? a.thickness() : 0.0,
                                     b.rotational() ? b.thickness() : 0.0);
    const double meshH = std::max({wallOf(a), wallOf(b), thick});
    const double pairCap = thick + midplane;
    // Tet panels already kiss. Capturing a full wall (beadTol = t) includes
    // the far skin and inboard rim Steiners; pairing those with the seam and
    // snapping 0.6 t pulls triangular bites out of the outer silhouette.
    const bool volumeWeld = a.hasVolume() && b.hasVolume();
    const double shellGap = volumeWeld ? 0.0 : std::max(a.thickness(), b.thickness());
    // Only pair nodes that already sit on the kiss. Snap can then close the
    // whole gap, so the coincident link does not yank on the first HHT step.
    const double kiss = std::max(0.05 * thick, 5.0e-4);
    // Midsurface edges of a corner joint sit about half a thickness off the
    // bead, and the two grids rarely share a node. Reach half a boundary
    // pitch past that standoff so every seam is tied; an untied seam is then
    // treated as a collision and the shells never settle.
    auto boundaryPitch = [](const Part& part) {
        double best = 1.0e9;
        const auto& ids = part.boundaryNodes();
        for (std::size_t i = 0; i < ids.size(); ++i) {
            if (ids[i] < 0)
                continue;
            const ChVector3d pi = part.nodePosition(static_cast<std::size_t>(ids[i]));
            for (std::size_t j = i + 1; j < ids.size(); ++j) {
                if (ids[j] < 0)
                    continue;
                const double d = (pi - part.nodePosition(static_cast<std::size_t>(ids[j]))).Length();
                if (d > 1.0e-4 && d < best)
                    best = d;
            }
        }
        return best < 1.0e8 ? best : std::max(part.thickness(), 1.0e-3);
    };
    const double shellReach =
        volumeWeld ? 0.0
                   : std::max({shellGap, 0.55 * std::max(boundaryPitch(a), boundaryPitch(b)), definition.tolerance});
    const double pairTol = volumeWeld && haveBead
                               ? kiss
                               : (haveBead ? std::max({definition.tolerance, 0.4 * meshH, shellGap, shellReach})
                                           : std::max(definition.tolerance, std::max(midplane, std::max(shellGap, shellReach))));
    const double beadTol =
        volumeWeld && haveBead ? kiss : (haveBead ? std::max(pairCap, shellReach) : definition.tolerance);

    struct Candidate {
        std::size_t nodeA = 0;
        std::size_t nodeB = 0;
        double distance = 0.0;
        ChVector3d mid{0, 0, 0};
    };
    std::vector<Candidate> candidates;
    for (int ia : a.boundaryNodes()) {
        if (ia < 0)
            continue;
        const ChVector3d pa = a.nodePosition(static_cast<std::size_t>(ia));
        if (!insideRegion(pa))
            continue;
        if (haveBead && distanceToBead(definition.bead, pa) > beadTol)
            continue;
        for (int ib : b.boundaryNodes()) {
            if (ib < 0)
                continue;
            const ChVector3d pb = b.nodePosition(static_cast<std::size_t>(ib));
            if (!insideRegion(pb))
                continue;
            if (haveBead && distanceToBead(definition.bead, pb) > beadTol)
                continue;
            const double distance = (pa - pb).Length();
            if (distance > pairTol)
                continue;
            candidates.push_back({static_cast<std::size_t>(ia), static_cast<std::size_t>(ib), distance,
                                  0.5 * (pa + pb)});
        }
    }
    std::sort(candidates.begin(), candidates.end(),
              [](const Candidate& l, const Candidate& r) { return l.distance < r.distance; });

    struct Pair {
        std::size_t nodeA;
        std::size_t nodeB;
        ChVector3d position;
    };
    std::vector<Pair> pairs;
    std::vector<char> usedA(a.nodeCount(), 0);
    std::vector<char> usedB(b.nodeCount(), 0);
    for (const Candidate& cand : candidates) {
        if (cand.nodeA >= usedA.size() || cand.nodeB >= usedB.size())
            continue;
        if (usedA[cand.nodeA] || usedB[cand.nodeB])
            continue;
        usedA[cand.nodeA] = 1;
        usedB[cand.nodeB] = 1;
        pairs.push_back({cand.nodeA, cand.nodeB, cand.mid});
    }

    auto tryAddPair = [&](std::size_t nodeA, std::size_t nodeB) {
        if (nodeA >= usedA.size() || nodeB >= usedB.size())
            return;
        if (usedA[nodeA] || usedB[nodeB])
            return;
        usedA[nodeA] = 1;
        usedB[nodeB] = 1;
        const ChVector3d pa = a.nodePosition(nodeA);
        const ChVector3d pb = b.nodePosition(nodeB);
        pairs.push_back({nodeA, nodeB, 0.5 * (pa + pb)});
    };
    auto nearestNode = [&](const Part& part, const std::vector<char>& used, const ChVector3d& target,
                           double maxDist, bool boundaryOnly) -> std::size_t {
        auto consider = [&](std::size_t n, double& bestD, std::size_t& best) {
            if (n >= used.size() || used[n] || part.nodeRigid(n))
                return;
            const ChVector3d p = part.nodePosition(n);
            if (!insideRegion(p))
                return;
            const double d = (p - target).Length();
            if (d <= maxDist && d < bestD) {
                bestD = d;
                best = n;
            }
        };
        std::size_t best = static_cast<std::size_t>(-1);
        double bestD = maxDist + 1.0;
        if (boundaryOnly) {
            for (int i : part.boundaryNodes()) {
                if (i >= 0)
                    consider(static_cast<std::size_t>(i), bestD, best);
            }
        } else {
            for (std::size_t n = 0; n < part.nodeCount(); ++n)
                consider(n, bestD, best);
        }
        return best;
    };

    if (pairs.empty() || !volumeWeld) {
        const double capture = volumeWeld ? kiss : std::max(shellReach, midplane);
        if (haveBead) {
            std::vector<ChVector3d> samples;
            const double sampleStep = std::max(0.5 * thick, 1.0e-3);
            const std::size_t n = definition.bead.size();
            const std::size_t last = n >= 3 ? n : n - 1;
            for (std::size_t i = 0; i < last; ++i) {
                const ChVector3d pa = definition.bead[i];
                const ChVector3d pb = definition.bead[(i + 1) % n];
                const ChVector3d ab = pb - pa;
                const double len = ab.Length();
                if (len < 1.0e-9)
                    continue;
                const int steps = std::max(1, static_cast<int>(std::round(len / sampleStep)));
                for (int s = 0; s < steps; ++s) {
                    const double t = static_cast<double>(s) / static_cast<double>(steps);
                    samples.push_back(pa + ab * t);
                }
            }
            auto pairFromSamples = [&](bool boundaryOnly) {
                for (const ChVector3d& sample : samples) {
                    const std::size_t na = nearestNode(a, usedA, sample, capture, boundaryOnly);
                    const std::size_t nb = nearestNode(b, usedB, sample, capture, boundaryOnly);
                    if (na == static_cast<std::size_t>(-1) || nb == static_cast<std::size_t>(-1))
                        continue;
                    tryAddPair(na, nb);
                }
            };
            pairFromSamples(true);
            if (pairs.empty())
                pairFromSamples(false);
        } else {
            for (int ia : a.boundaryNodes()) {
                if (ia < 0)
                    continue;
                const ChVector3d pa = a.nodePosition(static_cast<std::size_t>(ia));
                if (!insideRegion(pa))
                    continue;
                const std::size_t nb = nearestNode(b, usedB, pa, capture, true);
                if (nb != static_cast<std::size_t>(-1))
                    tryAddPair(static_cast<std::size_t>(ia), nb);
            }
            if (pairs.empty()) {
                for (std::size_t na = 0; na < a.nodeCount(); ++na) {
                    const ChVector3d pa = a.nodePosition(na);
                    if (!insideRegion(pa))
                        continue;
                    const std::size_t nb = nearestNode(b, usedB, pa, capture, false);
                    if (nb != static_cast<std::size_t>(-1))
                        tryAddPair(na, nb);
                }
            }
        }
    }
    if (pairs.empty())
        return;

    (void)snapped;
    (void)kiss;

    const double throat =
        definition.throat > 0.0 ? definition.throat : std::min(a.thickness(), b.thickness());
    const double parentAllowable =
        std::min(a.material().weldShearStrength, b.material().weldShearStrength) *
        definition.strengthScale;
    double allowable = definition.interface.shearStrength > 0.0
                           ? definition.interface.shearStrength * definition.strengthScale
                           : parentAllowable;
    if (definition.interface.boundedByParents)
        allowable = std::min(allowable, parentAllowable);
    const bool closedLoop = definition.bead.size() >= 3;
    const double faceArea = closedLoop ? closedBeadArea(definition.bead) : 0.0;
    const double peri = closedLoop ? beadPerimeter(definition.bead) : 0.0;

    Joint joint;
    joint.kind = JointKind::Weld;
    joint.label = definition.label;
    joint.partA = definition.partA;
    joint.partB = definition.partB;

    for (std::size_t i = 0; i < pairs.size(); ++i) {
        // Tributary bead length: distance to the nearest neighbouring link.
        // Along an edge weld this is the nodal spacing, so N links of spacing
        // L/N recover the fillet capacity throat * L * weldShearStrength.
        double spacing = 0.0;
        if (pairs.size() > 1) {
            double nearest = std::numeric_limits<double>::max();
            for (std::size_t j = 0; j < pairs.size(); ++j) {
                if (j == i)
                    continue;
                nearest = std::min(nearest, (pairs[i].position - pairs[j].position).Length());
            }
            spacing = std::max(nearest, 0.5 * throat);
        } else {
            spacing = throat;
        }

        double area = spacing * throat;
        if (faceArea > 1.0e-12) {
            if (pairs.size() == 1 || peri < 1.0e-9)
                area = faceArea;
            else
                area = faceArea * spacing / peri;
        }
        const double breakForce = allowable * area;
        const double breakTorque = breakForce * throat;
        const double peelPa = definition.interface.peelStrength > 0.0
                                  ? definition.interface.peelStrength * definition.strengthScale
                                  : (definition.peelStrength > 0.0
                                         ? definition.peelStrength * definition.strengthScale
                                         : allowable);
        const double derivedGc =
            std::min(modeIReleaseRate(a.material()), modeIReleaseRate(b.material()));
        const double release = definition.interface.releaseRate > 0.0 ? definition.interface.releaseRate
                                                                      : derivedGc;
        const double opening = cohesiveCriticalOpening(std::max(allowable, peelPa), release, throat);
        chrono::ChVector3d tangent = definition.bead.size() >= 2 ? definition.bead[1] - definition.bead[0]
                                                                 : chrono::ChVector3d{1, 0, 0};
        chrono::ChVector3d axis = tangent.Cross(b.center() - a.center());
        if (axis.Length2() < 1.0e-12)
            axis = {0, 1, 0};
        addLink(system, joint, definition.partA, pairs[i].nodeA, definition.partB, pairs[i].nodeB,
                breakForce, breakTorque, true, translationSet, rotationSet, peelPa * area, opening,
                axis.GetNormalized());
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

    struct Pair {
        std::size_t a = 0;
        std::size_t b = 0;
        double dist = 0.0;
    };
    std::vector<Pair> flange;
    flange.push_back({nodeA, nodeB, 0.0});
    if (definition.locksRotation) {
        const Part& pa = _parts[definition.partA];
        const Part& pb = _parts[definition.partB];
        for (std::size_t ia = 0; ia < pa.nodeCount(); ++ia) {
            if (ia == nodeA)
                continue;
            const double da = (pa.nodePosition(ia) - definition.position).Length();
            if (da > definition.clampRadius)
                continue;
            std::size_t ib = pb.nodeCount();
            double best = definition.clampRadius;
            for (std::size_t j = 0; j < pb.nodeCount(); ++j) {
                if (j == nodeB)
                    continue;
                const double db = (pb.nodePosition(j) - pa.nodePosition(ia)).Length();
                if (db < best) {
                    best = db;
                    ib = j;
                }
            }
            if (ib < pb.nodeCount())
                flange.push_back({ia, ib, da});
        }
        std::sort(flange.begin() + 1, flange.end(),
                  [](const Pair& x, const Pair& y) { return x.dist < y.dist; });
        if (flange.size() > 3)
            flange.resize(3);
    }
    // Extra flange pairs lock rotation; they are not extra bolts. Splitting
    // the proof load made a 0.3× tug pop the first pair.
    addLink(system, joint, definition.partA, flange[0].a, definition.partB, flange[0].b, breakForce,
            breakForce * d * 0.5, definition.locksRotation, translationSet, rotationSet, 0.0, 0.0,
            {0.0, 1.0, 0.0});
    for (std::size_t i = 1; i < flange.size(); ++i)
        createLink(system, joint, definition.partA, flange[i].a, definition.partB, flange[i].b, breakForce,
                   breakForce * d * 0.5, definition.locksRotation, definition.locksRotation);

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
            } else if (entry.offsetLink) {
                entry.offsetLink->Initialize(_parts[entry.partA].xyzNodes()[entry.nodeA],
                                             _parts[entry.partB].xyzNodes()[entry.nodeB]);
                const ChVector3d sep = _parts[entry.partA].nodePosition(entry.nodeA) -
                                       _parts[entry.partB].nodePosition(entry.nodeB);
                entry.offsetLink->SetRestOffset(sep);
                entry.restGap = sep.Length();
            } else if (entry.frameLink) {
                const bool rotA = _parts[entry.partA].rotational();
                const std::size_t xyzPart = rotA ? entry.partB : entry.partA;
                const std::size_t xyzNode = rotA ? entry.nodeB : entry.nodeA;
                const std::size_t framePart = rotA ? entry.partA : entry.partB;
                const std::size_t frameNode = rotA ? entry.nodeA : entry.nodeB;
                entry.frameLink->Initialize(
                    _parts[xyzPart].xyzNodes()[xyzNode],
                    std::static_pointer_cast<ChBodyFrame>(_parts[framePart].nodes()[frameNode]));
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

        auto nodeFixed = [&](std::size_t n) {
            if (n < part.xyzNodes().size() && part.xyzNodes()[n])
                return part.xyzNodes()[n]->IsFixed();
            if (n < part.nodes().size() && part.nodes()[n])
                return part.nodes()[n]->IsFixed();
            return true;
        };

        double totalMass = 0.0;
        std::vector<LoadTarget> targets;
        const bool volume = part.hasVolume();
        const auto& skin = part.boundaryNodes();
        auto onSkin = [&](std::size_t n) {
            if (!volume)
                return true;
            for (int index : skin) {
                if (index >= 0 && static_cast<std::size_t>(index) == n)
                    return true;
            }
            return false;
        };
        ChVector3d loadDir = load.direction;
        if (loadDir.Length() > 1.0e-12)
            loadDir.Normalize();
        else
            loadDir = ChVector3d(0, -1, 0);
        std::vector<LoadTarget> allSkin;
        double allMass = 0.0;
        for (std::size_t n = 0; n < part.nodeCount(); ++n) {
            if (nodeFixed(n) || part.nodeRigid(n))
                continue;
            if (!onSkin(n))
                continue;
            if ((part.nodePosition(n) - load.point).Length() > load.radius)
                continue;
            const LoadTarget entry{load.partIndex, n, part.nodeMass()[n]};
            allSkin.push_back(entry);
            allMass += entry.share;
            // A punch is a face pressure. Nodes on the far face of a thin
            // solid would turn the same force into a body load and cancel.
            if ((part.nodePosition(n) - load.point).Dot(loadDir) <= 0.0) {
                targets.push_back(entry);
                totalMass += entry.share;
            }
        }
        if (targets.empty()) {
            targets = std::move(allSkin);
            totalMass = allMass;
        }

        // Nothing inside the patch: fall back to the single closest free node
        // so a mis-placed load still does something visible. A fixed node
        // would swallow the force without stressing the mesh. Prefer a skin
        // node on a tet body so the 8 MN dent is a face pressure, not a
        // body force through the core.
        if (targets.empty()) {
            std::size_t best = part.nodeCount();
            double bestDistance = std::numeric_limits<double>::max();
            for (int pass = 0; pass < 2 && best >= part.nodeCount(); ++pass) {
                for (std::size_t n = 0; n < part.nodeCount(); ++n) {
                    if (nodeFixed(n) || part.nodeRigid(n))
                        continue;
                    if (pass == 0 && !onSkin(n))
                        continue;
                    const double distance = (part.nodePosition(n) - load.point).Length();
                    if (distance < bestDistance) {
                        bestDistance = distance;
                        best = n;
                    }
                }
            }
            if (best < part.nodeCount()) {
                targets.push_back({load.partIndex, best, 1.0});
                totalMass = 1.0;
            }
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

    const double total = loadMaxStepChange(loads[loadIndex], dt);
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

void Assembly::advanceHeat(double dt, double time, const Environment& environment,
                           const std::vector<HeatSource>& sources)
{
    std::vector<HeatBridge> bridges;
    for (const auto& joint : _joints) {
        for (const auto& link : joint.links) {
            if (link.broken)
                continue;
            bridges.push_back({link.partA, link.nodeA, link.partB, link.nodeB});
        }
    }
    for (const auto& link : _deferred)
        bridges.push_back({link.partA, link.nodeA, link.partB, link.nodeB});
    _heat.advance(dt, time, _parts, environment, sources, bridges);
}

bool Assembly::hasThermalPrestress() const {
    for (const auto& part : _parts) {
        if (part.thermalPrestress())
            return true;
    }
    return false;
}

bool Assembly::thermalLoadPending() const {
    if (_heat.strainMoved())
        return true;
    for (const auto& part : _parts) {
        if (part.thermalLoadChanged())
            return true;
    }
    return false;
}

void Assembly::clearThermalLoadChanged() {
    for (auto& part : _parts)
        part.clearThermalLoadChanged();
}

bool Assembly::applyLoads(const Environment& environment,
                          const std::vector<ImpactLoad>& loads,
                          double time,
                          double dt)
{
    for (auto& part : _parts)
        part.clearLoads();

    applyGravity(environment.gravity);
    if (environment.aerodynamicsEnabled)
        applyAerodynamics(environment);
    if (_contactWorld && dt > 0.0)
        _contactWorld->stepSize = dt;
    // The floor, the walls and the ceiling are implicit loads. They are only
    // refreshed here; the integrator evaluates them inside the step.
    bool contactsChanged = syncContacts(environment);
    ++_penaltyStamp;
    ++_pairStamp;
    _pairsChanged = false;
    // A rebuilt skin or a frozen node can leave a pair pointing at a node
    // the solver no longer integrates. Every pair is found again from scratch.
    if (contactsChanged)
        sweepPairContacts(true);
    if (_contactContainer)
        applyPartContact(environment, dt);
    sweepPairContacts();
    contactsChanged = contactsChanged || _pairsChanged;
    applyShardContact(environment, dt);
    applyUserLoads(loads, time);

    for (auto& part : _parts)
        part.commitLoads();
    return contactsChanged;
}

void Assembly::requestPairContact(std::size_t nodeSide,
                                  std::size_t node,
                                  std::size_t faceSide,
                                  const std::array<int, 3>& tri,
                                  const ChVector3d& normal,
                                  double wa,
                                  double wb,
                                  double wc,
                                  double reach,
                                  double stiffness,
                                  double friction,
                                  double stickSpeed) {
    const std::size_t count = _parts.size();
    const std::uint64_t key = (static_cast<std::uint64_t>(nodeSide * count + faceSide) << 32) |
                              static_cast<std::uint64_t>(node);
    PairContact& slot = _pairContacts[key];
    if (!slot.load || slot.tri != tri) {
        std::vector<std::shared_ptr<chrono::ChLoadable>> loadables;
        std::vector<NodeFaceContact::Corner> corners;
        auto add = [&](Part& part, std::size_t index) {
            NodeFaceContact::Corner corner;
            if (index < part.xyzNodes().size() && part.xyzNodes()[index])
                corner.xyz = part.xyzNodes()[index];
            else if (index < part.nodes().size() && part.nodes()[index])
                corner.rot = part.nodes()[index];
            else
                return false;
            const bool fixed = corner.xyz ? corner.xyz->IsFixed() : corner.rot->IsFixed();
            if (!fixed) {
                if (corner.xyz)
                    loadables.push_back(corner.xyz);
                else
                    loadables.push_back(corner.rot);
            }
            corners.push_back(corner);
            return true;
        };
        bool ok = add(_parts[nodeSide], node);
        for (int k = 0; k < 3 && ok; ++k)
            ok = add(_parts[faceSide], static_cast<std::size_t>(tri[k]));
        if (slot.load) {
            auto& list = _contactContainer->GetLoadList();
            list.erase(std::remove(list.begin(), list.end(), std::static_pointer_cast<chrono::ChLoadBase>(slot.load)),
                       list.end());
            slot.load.reset();
            _pairsChanged = true;
        }
        if (!ok || loadables.empty()) {
            _pairContacts.erase(key);
            return;
        }
        slot.load = chrono_types::make_shared<NodeFaceContact>(loadables, std::move(corners));
        slot.tri = tri;
        _contactContainer->Add(slot.load);
        _pairsChanged = true;
    }
    slot.stamp = _pairStamp;
    slot.load->configure(normal, {1.0, -wa, -wb, -wc}, reach, stiffness, friction, stickSpeed,
                         _contactWorld ? _contactWorld->restitution : 0.15);
}

void Assembly::sweepPairContacts(bool all) {
    for (auto it = _pairContacts.begin(); it != _pairContacts.end();) {
        if (!all && it->second.stamp == _pairStamp) {
            ++it;
            continue;
        }
        if (it->second.load && _contactContainer) {
            auto& list = _contactContainer->GetLoadList();
            list.erase(std::remove(list.begin(), list.end(),
                                   std::static_pointer_cast<chrono::ChLoadBase>(it->second.load)),
                       list.end());
        }
        it = _pairContacts.erase(it);
        _pairsChanged = true;
    }
}

double Assembly::penaltyDamping(std::uint64_t key, double vn) {
    ApproachLatch& latch = _penaltyApproach[key];
    if (latch.stamp + 1 != _penaltyStamp && latch.stamp != _penaltyStamp)
        latch.speed = std::max(0.0, -vn);
    latch.stamp = _penaltyStamp;
    return huntCrossleyFactor(_contactWorld ? _contactWorld->restitution : 0.15, latch.speed);
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
    (void)dt;
    for (auto& owner : _parts) {
        if (owner.shards().empty())
            continue;
        const double ownerReach = owner.contactOffset();
        for (auto& shard : owner.shards()) {
            // A chip that has landed is scenery. It keeps its pose, and it
            // does not push the body it fell off.
            if (shard.asleep || shard.mass <= 0.0)
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
                // A chip that just left this mesh must fall, not restick to
                // the hole it came from. Parent contact was pinning molten
                // drops at the cube COM while their speed still read as g*t.
                if (&other == &owner)
                    continue;
                const double reach = ownerReach + other.contactOffset();
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
                    const double patch = 0.5 * face.Length();
                    const ContactGains gains = continuumGains(owner.material().youngsModulus, patch, mass);
                    const std::uint64_t key =
                        (std::uint64_t{1} << 63) |
                        ((static_cast<std::uint64_t>(&owner - _parts.data()) * 1000003u +
                          static_cast<std::uint64_t>(&shard - owner.shards().data())) * 131u +
                         static_cast<std::uint64_t>(&other - _parts.data())) << 24 ^
                        static_cast<std::uint64_t>(i);
                    const double spring = gains.stiffness * overlap;
                    const double D = penaltyDamping(key, vn);
                    const double fn =
                        spring * huntCrossleyScale(D, vn, _penaltyApproach[key].speed).scale;
                    if (fn <= 0.0)
                        continue;
                    other.addLoad(i, outward * fn);
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
            if (a.asleep || b.asleep)
                continue;
            const double reach = all[i].part->contactOffset() + all[j].part->contactOffset();
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
        part.finishShards(environment.floorEnabled, environment.floorHeight, environment.boundsSize,
                          environment.gravity, environment.floorRestitution,
                          environment.contactFriction);
    sanitizeShards();
}

bool Assembly::syncContacts(const Environment& environment) {
    if (!_system)
        return false;
    if (!_contactWorld)
        _contactWorld = std::make_shared<ContactWorld>();
    bool structureChanged = false;
    if (!_contactContainer) {
        _contactContainer = chrono_types::make_shared<chrono::ChLoadContainer>();
        _system->Add(_contactContainer);
        structureChanged = true;
    }

    _contactWorld->floorEnabled = environment.floorEnabled;
    _contactWorld->floorHeight = environment.floorHeight;
    _contactWorld->friction = environment.contactFriction;
    _contactWorld->stickSpeed = std::max(1.0e-5, environment.frictionRegularization);
    _contactWorld->halfX = std::max(0.5, environment.boundsSize.x() * 0.5);
    _contactWorld->halfZ = std::max(0.5, environment.boundsSize.z() * 0.5);
    _contactWorld->ceiling = environment.floorHeight + std::max(1.0, environment.boundsSize.y());
    _contactWorld->restitution = std::clamp(environment.floorRestitution, 0.0, 1.0);

    _contacts.resize(_parts.size());
    _contactShape.resize(_parts.size(), {static_cast<std::size_t>(-1), 0, 0});

    for (std::size_t p = 0; p < _parts.size(); ++p) {
        Part& part = _parts[p];
        const auto& area = part.contactAreas();
        const std::size_t skin = static_cast<std::size_t>(std::max(0, part.skinTriangleCount()));
        // A fixed node has no solver variables, so a load may not reference
        // it: its Jacobian block would point at a disabled column.
        auto eligible = [&](std::size_t i) {
            if (i >= area.size() || area[i] <= 0.0)
                return false;
            if (i < part.xyzNodes().size() && part.xyzNodes()[i] && !part.xyzNodes()[i]->IsFixed())
                return part.hasVolume() ? part.nodeOnLiveTet(i) : part.nodeOnLiveShell(i);
            if (i < part.nodes().size() && part.nodes()[i] && !part.nodes()[i]->IsFixed())
                return true;
            return false;
        };
        // A digest of which nodes are eligible, not how many. Counting them is
        // not enough: one node freezing in the same step that another thaws
        // leaves the count alone, and the contact on the frozen node survives
        // pointing at a column the solver has switched off. That is an invalid
        // matrix, and Pardiso either rejects it or reads off the end of it.
        std::size_t pattern = 1469598103934665603ull;
        for (std::size_t i = 0; i < part.nodeCount(); ++i) {
            if (!eligible(i))
                continue;
            pattern = (pattern ^ i) * 1099511628211ull;
        }
        // Every free surface node owns an element for as long as the mesh
        // holds still, so the sparsity pattern only moves when the mesh does.
        // A node that is not touching contributes a zero block rather than
        // appearing and disappearing, which would refactorise every step.
        const std::array<std::size_t, 3> shape{part.nodeCount(), skin, pattern};
        if (_contactShape[p] != shape) {
            auto& list = _contactContainer->GetLoadList();
            for (auto& contact : _contacts[p]) {
                if (!contact)
                    continue;
                list.erase(std::remove(list.begin(), list.end(),
                                       std::static_pointer_cast<chrono::ChLoadBase>(contact)),
                           list.end());
            }
            _contacts[p].assign(part.nodeCount(), nullptr);
            _contactShape[p] = shape;
            for (std::size_t i = 0; i < part.nodeCount(); ++i) {
                if (!eligible(i))
                    continue;
                if (i < part.xyzNodes().size() && part.xyzNodes()[i])
                    _contacts[p][i] = chrono_types::make_shared<NodeContact>(part.xyzNodes()[i], _contactWorld);
                else if (i < part.nodes().size() && part.nodes()[i])
                    _contacts[p][i] = chrono_types::make_shared<NodeContact>(part.nodes()[i], _contactWorld);
                else
                    continue;
                _contactContainer->Add(_contacts[p][i]);
            }
            structureChanged = true;
        }

        // A contact patch behaves as an elastic foundation: the material under
        // it compresses over roughly its own width. That gives E*A/L with
        // L = sqrt(A), the real stiffness of metal on concrete, and it is only
        // usable because the integrator differentiates it.
        const double modulus = std::max(1.0e5, part.material().youngsModulus);
        const double offset = part.contactOffset();
        for (std::size_t i = 0; i < part.nodeCount() && i < _contacts[p].size(); ++i) {
            if (!_contacts[p][i])
                continue;
            const double patch = i < area.size() ? area[i] : 0.0;
            // A seam node of a welded volume can sit inside a neighbour's
            // tets. That is not burial: it is still a free surface vs the
            // floor. World contact stays on; part-part is applyPartContact.
            const bool live = eligible(i) && !part.nodeRigid(i);
            _contacts[p][i]->setActive(live);
            if (!live)
                continue;
            const double mass = i < part.nodeMass().size() ? part.nodeMass()[i] : 0.0;
            if (mass <= 1.0e-12 || patch <= 0.0) {
                _contacts[p][i]->setActive(false);
                continue;
            }
            const ContactGains gains = continuumGains(modulus, patch, mass);
            _contacts[p][i]->configure(offset, gains.stiffness, 0.0, mass);
        }
    }
    return structureChanged;
}

void Assembly::applyPartContact(const Environment& environment, double dt) {
    (void)dt;
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

    std::vector<ChVector3d> lo(count), hi(count);
    for (std::size_t i = 0; i < count; ++i)
        _parts[i].bounds(lo[i], hi[i]);

    // Node-against-current-triangle proximity: the contact surface is whatever
    // shape the mesh has right now, deformed or not. Returns true when the two
    // parts overlap at all, which is what arming uses.
    auto resolve = [&](std::size_t nodeSide, std::size_t faceSide, bool apply, bool self, bool fresh,
                       double forceScale, const std::vector<ChVector3d>* weldPts, double weldR,
                       bool skipProximity, const VolumeGrid* hostGrid) {
        bool touching = false;
        Part& nodePart = _parts[nodeSide];
        Part& facePart = _parts[faceSide];
        const double reach = nodePart.contactSkin() + facePart.contactSkin();
        // A pair is built a skin depth before it touches, so the step that
        // closes the gap already carries the contact stiffness.
        const double search = 2.0 * reach;
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
            if (nodePart.hasVolume() && !nodePart.nodeOnLiveTet(i))
                continue;
            const ChVector3d p = nodePart.nodePosition(i);
            const bool insideVolume =
                !self && facePart.hasVolume() &&
                p.x() >= lo[faceSide].x() && p.x() <= hi[faceSide].x() &&
                p.y() >= lo[faceSide].y() && p.y() <= hi[faceSide].y() &&
                p.z() >= lo[faceSide].z() && p.z() <= hi[faceSide].z() &&
                (hostGrid ? hostGrid->contains(p) : pointInLiveTets(facePart, p));
            if (insideVolume) {
                touching = true;
                if (!apply)
                    continue;
                double best = 1.0e30;
                int bestTri = -1;
                ChVector3d bestA, bestB, bestC, bestHit, bestN;
                for (int t = 0; t < static_cast<int>(faces.size()); ++t) {
                    const auto& tri = faces[static_cast<std::size_t>(t)];
                    const ChVector3d a = facePart.nodePosition(tri[0]);
                    const ChVector3d b = facePart.nodePosition(tri[1]);
                    const ChVector3d c = facePart.nodePosition(tri[2]);
                    const ChVector3d hit = closestPointOnTriangle(p, a, b, c);
                    const double d = (p - hit).Length();
                    if (d >= best)
                        continue;
                    ChVector3d face = (b - a).Cross(c - a);
                    if (face.Length() < 1.0e-12)
                        continue;
                    best = d;
                    bestTri = t;
                    bestA = a;
                    bestB = b;
                    bestC = c;
                    bestHit = hit;
                    bestN = face.GetNormalized();
                }
                if (bestTri < 0)
                    continue;
                // An intact weld is a seam, not a burial. Penalty-pushing
                // kissing panels here stole ~0.7 g from a free-falling box.
                if (skipProximity)
                    continue;
                if (weldPts && nearPoints(*weldPts, p, weldR))
                    continue;
                ChVector3d normal = best > 1.0e-9 ? (bestHit - p) / best : bestN;
                if (normal.Dot(bestN) < 0.0)
                    normal = -normal;
                double wa, wb, wc;
                barycentric(bestHit, bestA, bestB, bestC, wa, wb, wc);
                wa = std::max(0.0, wa);
                wb = std::max(0.0, wb);
                wc = std::max(0.0, wc);
                const double wsum = wa + wb + wc;
                if (wsum > 1.0e-12) {
                    wa /= wsum;
                    wb /= wsum;
                    wc /= wsum;
                } else {
                    wa = wb = wc = 1.0 / 3.0;
                }
                const auto& tri = faces[static_cast<std::size_t>(bestTri)];
                const double patch =
                    0.5 * (bestB - bestA).Cross(bestC - bestA).Length();
                const ContactGains gains =
                    nodeFaceGains(nodePart, facePart, i, tri[0], tri[1], tri[2], patch);
                requestPairContact(nodeSide, i, faceSide, {tri[0], tri[1], tri[2]}, normal, wa, wb, wc, reach,
                                   forceScale * gains.stiffness, mu, vStick);
                continue;
            }
            if (skipProximity)
                continue;
            if (!self && weldPts && nearPoints(*weldPts, p, weldR))
                continue;
            const int cx = static_cast<int>(std::floor(p.x() / cell));
            const int cy = static_cast<int>(std::floor(p.y() / cell));
            const int cz = static_cast<int>(std::floor(p.z() / cell));
            ++stamp;
            if (stamp == 0) {
                std::fill(seen.begin(), seen.end(), 0);
                stamp = 1;
            }
            double nearest = search;
            int pick = -1;
            ChVector3d pickNormal;
            std::array<double, 3> pickW{};
            double pickPatch = 0.0;
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
                if (std::abs(planeGap) >= search)
                    continue;
                const ChVector3d projected = p - nHat * planeGap;
                double wa, wb, wc;
                barycentric(projected, a, b, c, wa, wb, wc);
                // On the edge means beside the shell, not through it.
                if (wa <= 1.0e-3 || wb <= 1.0e-3 || wc <= 1.0e-3)
                    continue;
                const double distance = std::abs(planeGap);
                touching = touching || distance < reach;
                if (!apply || distance >= nearest)
                    continue;
                nearest = distance;
                pick = t;
                pickNormal = planeGap >= 0.0 ? nHat : -nHat;
                pickW = {wa, wb, wc};
                pickPatch = 0.5 * faceLen;
                        }
                    }
                }
            }
            if (pick < 0)
                continue;
            ChVector3d normal = pickNormal;
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
            (void)fresh;
            const auto& tri = faces[static_cast<std::size_t>(pick)];
            const ContactGains gains = nodeFaceGains(nodePart, facePart, i, tri[0], tri[1], tri[2], pickPatch);
            requestPairContact(nodeSide, i, faceSide, {tri[0], tri[1], tri[2]}, normal, pickW[0], pickW[1],
                               pickW[2], reach, forceScale * gains.stiffness, mu, vStick);
        }
        return touching;
    };

    std::vector<char> seam;
    classifySeams(seam);
    std::vector<VolumeGrid> grids(count);
    for (std::size_t i = 0; i < count; ++i)
        grids[i].build(_parts[i]);

    std::vector<std::vector<ChVector3d>> liveWeld(count * count);
    for (const auto& joint : _joints) {
        if (joint.partA >= count || joint.partB >= count || joint.partA == joint.partB)
            continue;
        const std::size_t a = std::min(joint.partA, joint.partB);
        const std::size_t b = std::max(joint.partA, joint.partB);
        const std::size_t key = a * count + b;
        if (joint.severed() || seam[key] != kSeamDamaged)
            continue;
        for (const auto& link : joint.links) {
            if (link.broken)
                continue;
            if (link.partA < count && link.nodeA < _parts[link.partA].nodeCount())
                liveWeld[key].push_back(_parts[link.partA].nodePosition(link.nodeA));
            if (link.partB < count && link.nodeB < _parts[link.partB].nodeCount())
                liveWeld[key].push_back(_parts[link.partB].nodePosition(link.nodeB));
        }
    }

    for (std::size_t i = 0; i < count; ++i) {
        for (std::size_t j = i + 1; j < count; ++j) {
            const std::size_t key = i * count + j;
            char& armed = _contactArmed[key];
            if (!aabbOverlap(lo[i], hi[i], lo[j], hi[j])) {
                armed = 1;
                continue;
            }

            if (armed == 0)
                armed = 1;

            const double weldR = std::max({0.008, 3.0 * _parts[i].contactSkin(),
                                           3.0 * _parts[j].contactSkin()});
            const bool skipProximity = seam[key] == kSeamIntact;
            const std::vector<ChVector3d>* weldPts =
                liveWeld[key].empty() ? nullptr : &liveWeld[key];

            const bool apply = armed != 0;
            const bool fresh = armed == 1;
            const bool touchingA =
                resolve(i, j, apply, false, fresh, 0.5, weldPts, weldR, skipProximity, &grids[j]);
            const bool touchingB =
                resolve(j, i, apply, false, fresh, 0.5, weldPts, weldR, skipProximity, &grids[i]);
            if (apply && (touchingA || touchingB))
                armed = 2;
            if (!apply && !touchingA && !touchingB)
                armed = 1;
        }
    }

    for (std::size_t i = 0; i < count; ++i) {
        if (_parts[i].torn())
            resolve(i, i, true, true, false, 1.0, nullptr, 0.0, false, nullptr);
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
            _parts[target.part].addStructuralLoad(target.node, direction * (loadMagnitude(load, time) * target.share));
    }
}

bool Assembly::updateJoints(double dt) {
    const double alpha = dt / (kReactionFilterTau + dt);
    bool changed = false;

    auto release = [&](JointLink& entry) {
        auto drop = [](const auto& link) {
            if (!link)
                return;
            link->SetBroken(true);
            link->SetDisabled(true);
        };
        drop(entry.link);
        drop(entry.pointLink);
        drop(entry.offsetLink);
        drop(entry.frameLink);
        entry.broken = true;
        entry.filteredForce = 0.0;
        entry.filteredForceVec = chrono::ChVector3d(0, 0, 0);
        entry.filteredTorque = 0.0;
        entry.utilisation = 0.0;
        entry.overloadSteps = 0;
        changed = true;
    };

    for (auto& joint : _joints) {
        for (auto& entry : joint.links) {
            if (entry.broken || (!entry.link && !entry.pointLink && !entry.offsetLink && !entry.frameLink))
                continue;

            if (entry.partA < _parts.size() && entry.partB < _parts.size()) {
                Part& sideA = _parts[entry.partA];
                Part& sideB = _parts[entry.partB];
                const bool hot = sideA.heatActive(entry.nodeA) || sideB.heatActive(entry.nodeB);
                auto liveFea = [](const Part& part, std::size_t index) {
                    if (index >= part.nodeCount() || part.nodeRigid(index))
                        return false;
                    if (part.hasVolume())
                        return part.nodeOnLiveTet(index);
                    return part.nodeOnLiveShell(index);
                };
                // A chip on one side of a bead must not pin the partner's live
                // mesh. Break as soon as either end has left the deformable body.
                if (!liveFea(sideA, entry.nodeA) || !liveFea(sideB, entry.nodeB)) {
                    release(entry);
                    continue;
                }
                if (hot) {
                    const double solidus =
                        std::min(sideA.material().meltTemperature, sideB.material().meltTemperature);
                    const double cooler = std::min(sideA.nodeTemperature(entry.nodeA),
                                                   sideB.nodeTemperature(entry.nodeB));
                    if (cooler >= solidus) {
                        release(entry);
                        continue;
                    }
                }
            }

            double force = 0.0;
            double torque = 0.0;
            chrono::ChVector3d forceVec{0, 0, 0};
            if (entry.link) {
                const auto reaction = entry.link->GetReaction1();
                forceVec = reaction.force;
                force = forceVec.Length();
                torque = reaction.torque.Length();
            } else if (entry.pointLink) {
                forceVec = entry.pointLink->GetReactionOnNode1();
                force = forceVec.Length();
            } else if (entry.offsetLink) {
                forceVec = entry.offsetLink->GetReactionOnNode1();
                force = forceVec.Length();
            } else if (entry.frameLink) {
                forceVec = entry.frameLink->GetReactionOnNode();
                force = forceVec.Length();
            }
            if (!std::isfinite(force) || !std::isfinite(torque))
                continue;

            entry.filteredForce += alpha * (force - entry.filteredForce);
            entry.filteredTorque += alpha * (torque - entry.filteredTorque);
            entry.filteredForceVec += (forceVec - entry.filteredForceVec) * alpha;
        }

        double weldCap = 0.0;
        double peelCap = 0.0;
        chrono::ChVector3d weldForce{0, 0, 0};
        chrono::ChVector3d peelAxis{0, 0, 0};
        if (joint.kind == JointKind::Weld) {
            for (const auto& entry : joint.links) {
                if (entry.broken)
                    continue;
                weldCap += std::max(entry.breakForce, 0.0);
                peelCap += std::max(entry.breakPeel, 0.0);
                weldForce += entry.filteredForceVec;
                if (entry.breakPeel > 0.0)
                    peelAxis += entry.peelAxis;
            }
        }
        double weldUtil = 0.0;
        if (peelCap > 1.0e-12 && peelAxis.Length2() > 0.25) {
            peelAxis.Normalize();
            const double peel = std::abs(weldForce.Dot(peelAxis));
            const double shear = (weldForce - peelAxis * weldForce.Dot(peelAxis)).Length();
            const double shearUtil = weldCap > 1.0e-12 ? shear / weldCap : 0.0;
            weldUtil = std::max(shearUtil, peel / peelCap);
        } else {
            const double weldLoad = weldForce.Length();
            weldUtil = weldCap > 1.0e-12 ? weldLoad / weldCap : 0.0;
        }
        if (joint.kind == JointKind::Weld) {
            for (const auto& entry : joint.links) {
                if (entry.broken || entry.criticalOpening <= 0.0)
                    continue;
                if (entry.partA >= _parts.size() || entry.partB >= _parts.size())
                    continue;
                const chrono::ChVector3d sep = _parts[entry.partB].nodePosition(entry.nodeB) -
                                               _parts[entry.partA].nodePosition(entry.nodeA);
                chrono::ChVector3d axis = entry.peelAxis;
                if (axis.Length2() < 0.25)
                    continue;
                axis.Normalize();
                const double extra = std::max(0.0, sep.Length() - entry.restGap);
                const CohesiveDamage damage =
                    cohesiveDamage(0.0, extra, 0.0, 0.0, entry.breakPeel, entry.breakForce, entry.criticalOpening);
                weldUtil = std::max(weldUtil, damage.damage);
            }
        }

        for (auto& entry : joint.links) {
            if (entry.broken || (!entry.link && !entry.pointLink && !entry.offsetLink && !entry.frameLink))
                continue;
            double utilisation =
                entry.breakForce > 0.0 ? entry.filteredForce / entry.breakForce : 0.0;
            if (entry.locksRotation && entry.breakTorque > 0.0) {
                utilisation =
                    std::max(utilisation, entry.filteredTorque / entry.breakTorque);
            }
            // A tension+bend couple overloads the edge of a translation-only
            // bead. The fillet is rated on the whole throat, not one node.
            if (joint.kind == JointKind::Weld)
                utilisation = weldUtil;
            entry.utilisation = utilisation;

            if (utilisation > 1.0)
                ++entry.overloadSteps;
            else
                entry.overloadSteps = 0;

            if (entry.overloadSteps * dt >= 0.04 && entry.overloadSteps >= kOverloadHoldSteps)
                release(entry);
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

bool Assembly::discardAwakeShards() {
    bool any = false;
    for (auto& part : _parts)
        any = part.discardAwakeShards() || any;
    return any;
}

void Assembly::sanitizeShards() {
    for (auto& part : _parts)
        part.sanitizeShards();
}

void Assembly::discardTrial() {
    for (auto& part : _parts)
        part.discardTrial();
}

void Assembly::commitPlastic() {
    for (auto& part : _parts)
        part.commitPlastic();
}

void Assembly::enforceRigidNodes() {
    for (auto& part : _parts)
        part.enforceRigidNodes();
}

void Assembly::applyFloorBounce(const Environment&, bool) {
    // Contact and friction live in NodeContact::ComputeQ. Post-step velocity
    // writes are not a dynamics path.
}

bool Assembly::touchingFloor(const Environment& environment, double pad) const {
    if (!environment.floorEnabled)
        return false;
    const double plane = environment.floorHeight + pad;
    for (const auto& part : _parts) {
        const double offset = part.contactOffset();
        for (std::size_t i = 0; i < part.nodeCount(); ++i) {
            if (part.nodeRigid(i))
                continue;
            if (part.nodePosition(i).y() < plane + offset)
                return true;
        }
    }
    return false;
}

int Assembly::activeContactCount() const {
    if (!_contactWorld)
        return 0;
    int count = 0;
    auto penetrating = [&](const chrono::ChVector3d& p, double offset) {
        const auto gap = [&](const chrono::ChVector3d& n, double plane) {
            return plane + offset - p.Dot(n);
        };
        if (_contactWorld->floorEnabled && gap(chrono::ChVector3d(0, 1, 0), _contactWorld->floorHeight) > 0.0)
            return true;
        if (gap(chrono::ChVector3d(1, 0, 0), -_contactWorld->halfX) > 0.0)
            return true;
        if (gap(chrono::ChVector3d(-1, 0, 0), -_contactWorld->halfX) > 0.0)
            return true;
        if (gap(chrono::ChVector3d(0, 0, 1), -_contactWorld->halfZ) > 0.0)
            return true;
        if (gap(chrono::ChVector3d(0, 0, -1), -_contactWorld->halfZ) > 0.0)
            return true;
        if (gap(chrono::ChVector3d(0, -1, 0), -_contactWorld->ceiling) > 0.0)
            return true;
        return false;
    };
    for (std::size_t p = 0; p < _parts.size(); ++p) {
        const Part& part = _parts[p];
        const double offset = part.contactOffset();
        const std::size_t nMax = p < _contacts.size() ? _contacts[p].size() : 0;
        for (std::size_t i = 0; i < nMax; ++i) {
            if (!_contacts[p][i] || !_contacts[p][i]->active())
                continue;
            if (penetrating(part.nodePosition(i), offset))
                ++count;
        }
    }
    for (const auto& [key, pair] : _pairContacts)
        if (pair.load && pair.load->gap() > 0.0)
            ++count;
    return count;
}

double Assembly::maxApproachSpeed() const {
    if (!_contactWorld)
        return 0.0;
    double fastest = 0.0;
    auto closing = [&](const chrono::ChVector3d& p, const chrono::ChVector3d& v, double offset) {
        const auto into = [&](const chrono::ChVector3d& n, double plane) {
            const double gap = plane + offset - p.Dot(n);
            if (gap <= 0.0)
                return;
            fastest = std::max(fastest, -v.Dot(n));
        };
        if (_contactWorld->floorEnabled)
            into(chrono::ChVector3d(0, 1, 0), _contactWorld->floorHeight);
        into(chrono::ChVector3d(1, 0, 0), -_contactWorld->halfX);
        into(chrono::ChVector3d(-1, 0, 0), -_contactWorld->halfX);
        into(chrono::ChVector3d(0, 0, 1), -_contactWorld->halfZ);
        into(chrono::ChVector3d(0, 0, -1), -_contactWorld->halfZ);
        into(chrono::ChVector3d(0, -1, 0), -_contactWorld->ceiling);
    };
    for (std::size_t p = 0; p < _parts.size(); ++p) {
        const Part& part = _parts[p];
        const double offset = part.contactOffset();
        const std::size_t nMax = p < _contacts.size() ? _contacts[p].size() : 0;
        for (std::size_t i = 0; i < nMax; ++i) {
            if (!_contacts[p][i] || !_contacts[p][i]->active())
                continue;
            closing(part.nodePosition(i), part.nodeVelocity(i), offset);
        }
    }
    return fastest;
}

double Assembly::minTimeToContact(const chrono::ChVector3d& gravity) const {
    if (!_contactWorld)
        return kNoContactTime;
    double soonest = kNoContactTime;
    auto consider = [&](const chrono::ChVector3d& p, const chrono::ChVector3d& v, double offset,
                        const chrono::ChVector3d& n, double plane) {
        const double open = p.Dot(n) - (plane + offset);
        const double closing = -v.Dot(n);
        const double acc = -gravity.Dot(n);
        const double t = timeToClosingPlane(open, closing, acc);
        if (t < soonest)
            soonest = t;
    };
    for (std::size_t p = 0; p < _parts.size(); ++p) {
        const Part& part = _parts[p];
        const double offset = part.contactOffset();
        const std::size_t nMax = p < _contacts.size() ? _contacts[p].size() : 0;
        for (std::size_t i = 0; i < nMax; ++i) {
            if (!_contacts[p][i] || !_contacts[p][i]->active())
                continue;
            const chrono::ChVector3d pos = part.nodePosition(i);
            const chrono::ChVector3d vel = part.nodeVelocity(i);
            if (_contactWorld->floorEnabled)
                consider(pos, vel, offset, chrono::ChVector3d(0, 1, 0), _contactWorld->floorHeight);
            consider(pos, vel, offset, chrono::ChVector3d(1, 0, 0), -_contactWorld->halfX);
            consider(pos, vel, offset, chrono::ChVector3d(-1, 0, 0), -_contactWorld->halfX);
            consider(pos, vel, offset, chrono::ChVector3d(0, 0, 1), -_contactWorld->halfZ);
            consider(pos, vel, offset, chrono::ChVector3d(0, 0, -1), -_contactWorld->halfZ);
            consider(pos, vel, offset, chrono::ChVector3d(0, -1, 0), -_contactWorld->ceiling);
        }
    }
    return soonest;
}

double Assembly::contactHalfPeriod() const {
    if (!_contactWorld)
        return kNoContactTime;
    double shortest = kNoContactTime;
    for (std::size_t p = 0; p < _parts.size(); ++p) {
        const Part& part = _parts[p];
        const double offset = part.contactOffset();
        const std::size_t nMax = p < _contacts.size() ? _contacts[p].size() : 0;
        for (std::size_t i = 0; i < nMax; ++i) {
            const auto& contact = _contacts[p][i];
            if (!contact || !contact->active() || contact->stiffness() <= 0.0)
                continue;
            if (_contactWorld->floorEnabled &&
                part.nodePosition(i).y() >= _contactWorld->floorHeight + offset)
                continue;
            const double m = i < part.nodeMass().size() ? part.nodeMass()[i] : 0.0;
            if (m <= 1.0e-12)
                continue;
            shortest = std::min(shortest, 3.14159265358979323846 * std::sqrt(m / contact->stiffness()));
        }
    }
    return shortest;
}

void Assembly::commitContactApproach() {
    for (auto& contacts : _contacts)
        for (auto& contact : contacts)
            if (contact)
                contact->commitApproach();
    for (auto& [key, pair] : _pairContacts)
        if (pair.load)
            pair.load->commitApproach();
}

double Assembly::maxContactPenetration() const {
    if (!_contactWorld)
        return 0.0;
    double worst = 0.0;
    auto consider = [&](const chrono::ChVector3d& p, double offset) {
        const auto into = [&](const chrono::ChVector3d& n, double plane) {
            worst = std::max(worst, plane + offset - p.Dot(n));
        };
        if (_contactWorld->floorEnabled)
            into(chrono::ChVector3d(0, 1, 0), _contactWorld->floorHeight);
        into(chrono::ChVector3d(1, 0, 0), -_contactWorld->halfX);
        into(chrono::ChVector3d(-1, 0, 0), -_contactWorld->halfX);
        into(chrono::ChVector3d(0, 0, 1), -_contactWorld->halfZ);
        into(chrono::ChVector3d(0, 0, -1), -_contactWorld->halfZ);
        into(chrono::ChVector3d(0, -1, 0), -_contactWorld->ceiling);
    };
    for (std::size_t p = 0; p < _parts.size(); ++p) {
        const Part& part = _parts[p];
        const double offset = part.contactOffset();
        const std::size_t nMax = p < _contacts.size() ? _contacts[p].size() : 0;
        for (std::size_t i = 0; i < nMax; ++i) {
            if (!_contacts[p][i] || !_contacts[p][i]->active())
                continue;
            consider(part.nodePosition(i), offset);
        }
    }
    return std::max(0.0, worst);
}

double Assembly::contactLengthScale() const {
    double scale = 0.0;
    for (const auto& part : _parts) {
        scale = std::max(scale, part.contactOffset());
        scale = std::max(scale, part.contactSkin());
    }
    return std::max(scale, 1.0e-4);
}

bool Assembly::pastBrittleRankine() const {
    for (const auto& part : _parts) {
        if (part.pastBrittleRankine())
            return true;
    }
    return false;
}

double Assembly::maxFloorOvershoot(const Environment& environment) const {
    if (!environment.floorEnabled)
        return 0.0;
    double worst = 0.0;
    for (const auto& part : _parts) {
        const double plane = environment.floorHeight + part.contactOffset();
        for (std::size_t i = 0; i < part.nodeCount(); ++i) {
            if (part.nodeRigid(i))
                continue;
            if (i < part.xyzNodes().size() && part.xyzNodes()[i] && part.xyzNodes()[i]->IsFixed())
                continue;
            if (i < part.nodes().size() && part.nodes()[i] && part.nodes()[i]->IsFixed())
                continue;
            worst = std::max(worst, plane - part.nodePosition(i).y());
        }
    }
    return worst;
}

void Assembly::rigidAdvance(double, const Environment&) {
    // Live path is Chrono HHT. Do not author a rigid pose here.
}

void Assembly::ballisticAdvance(double, const chrono::ChVector3d&) {
    // Live path is Chrono HHT. Free-fall is gravity loads with no contact force.
}

void Assembly::refreshTetTrials() {
    for (auto& part : _parts)
        part.refreshTetTrials();
}

void Assembly::explicitAdvance(double dt) {
    if (dt <= 0.0)
        return;
    auto free = [](const Part& part, std::size_t i) {
        if (part.nodeRigid(i))
            return false;
        if (i < part.xyzNodes().size() && part.xyzNodes()[i] && part.xyzNodes()[i]->IsFixed())
            return false;
        if (i < part.nodes().size() && part.nodes()[i] && part.nodes()[i]->IsFixed())
            return false;
        return true;
    };
    for (std::size_t p = 0; p < _parts.size(); ++p) {
        auto& part = _parts[p];
        for (std::size_t i = 0; i < part.nodeCount(); ++i) {
            if (!free(part, i))
                continue;
            const double mass = i < part.nodeMass().size() ? part.nodeMass()[i] : 0.0;
            if (mass <= 1.0e-12)
                continue;
            chrono::ChVector3d force = part.nodeForce(i);
            if (p < _contacts.size() && i < _contacts[p].size() && _contacts[p][i] &&
                _contacts[p][i]->active()) {
                _contacts[p][i]->ComputeQ(nullptr, nullptr);
                const auto& q = _contacts[p][i]->load_Q;
                if (q.size() >= 3)
                    force += chrono::ChVector3d(q(0), q(1), q(2));
            }
            const chrono::ChVector3d a = force / mass;
            const chrono::ChVector3d v0 = part.nodeVelocity(i);
            const chrono::ChVector3d v1 = v0 + a * dt;
            part.setNodeKinematics(i, part.nodePosition(i) + v0 * dt + a * (0.5 * dt * dt), v1);
        }
    }
}

double Assembly::peakPlastic() const {
    double peak = 0.0;
    for (const auto& part : _parts) {
        std::vector<float> field;
        part.evaluatePlasticStrain(field);
        for (float value : field) {
            if (std::isfinite(value))
                peak = std::max(peak, static_cast<double>(value));
        }
    }
    return peak;
}

void Assembly::captureKinematics() {
    for (auto& part : _parts)
        part.captureKinematics();
}

void Assembly::restoreKinematics() {
    for (auto& part : _parts)
        part.restoreKinematics();
}

double Assembly::kineticEnergy() const {
    double ke = 0.0;
    for (const auto& part : _parts) {
        for (std::size_t n = 0; n < part.nodeCount(); ++n) {
            const double m = n < part.nodeMass().size() ? part.nodeMass()[n] : 0.0;
            ke += 0.5 * m * part.nodeVelocity(n).Length2();
        }
        for (const auto& shard : part.shards())
            ke += 0.5 * shard.mass * shard.velocity.Length2();
    }
    return ke;
}

double Assembly::contactElasticEnergy() const {
    if (!_contactWorld)
        return 0.0;
    double energy = 0.0;
    auto add = [&](const chrono::ChVector3d& p, double offset, double k) {
        if (k <= 0.0)
            return;
        const auto spring = [&](const chrono::ChVector3d& n, double plane) {
            const double gap = plane + offset - p.Dot(n);
            if (gap <= 0.0)
                return;
            const double soft = kContactSoftening;
            energy += gap < soft ? k * gap * gap * gap / (6.0 * soft)
                                 : k * (soft * soft / 6.0 + 0.5 * (gap - soft) * (gap - soft) +
                                        0.5 * soft * (gap - soft));
        };
        if (_contactWorld->floorEnabled)
            spring(chrono::ChVector3d(0, 1, 0), _contactWorld->floorHeight);
        spring(chrono::ChVector3d(1, 0, 0), -_contactWorld->halfX);
        spring(chrono::ChVector3d(-1, 0, 0), -_contactWorld->halfX);
        spring(chrono::ChVector3d(0, 0, 1), -_contactWorld->halfZ);
        spring(chrono::ChVector3d(0, 0, -1), -_contactWorld->halfZ);
        spring(chrono::ChVector3d(0, -1, 0), -_contactWorld->ceiling);
    };
    for (std::size_t p = 0; p < _parts.size(); ++p) {
        const Part& part = _parts[p];
        const double offset = part.contactOffset();
        const std::size_t nMax = p < _contacts.size() ? _contacts[p].size() : 0;
        for (std::size_t i = 0; i < nMax; ++i) {
            if (!_contacts[p][i] || !_contacts[p][i]->active())
                continue;
            add(part.nodePosition(i), offset, _contacts[p][i]->stiffness());
        }
    }
    for (const auto& [key, pair] : _pairContacts)
        if (pair.load)
            energy += pair.load->elasticEnergy();
    return energy;
}

double Assembly::massHeight() const {
    double mass = 0.0;
    double my = 0.0;
    for (const auto& part : _parts) {
        for (std::size_t n = 0; n < part.nodeCount(); ++n) {
            const double m = n < part.nodeMass().size() ? part.nodeMass()[n] : 0.0;
            mass += m;
            my += m * part.nodePosition(n).y();
        }
        for (const auto& shard : part.shards()) {
            mass += shard.mass;
            my += shard.mass * shard.com.y();
        }
    }
    return mass > 0.0 ? my / mass : 0.0;
}

void Assembly::classifySeams(std::vector<char>& state) const {
    const std::size_t count = _parts.size();
    state.assign(count * count, kSeamNone);
    for (const auto& joint : _joints) {
        if (joint.partA >= count || joint.partB >= count || joint.partA == joint.partB)
            continue;
        const std::size_t a = std::min(joint.partA, joint.partB);
        const std::size_t b = std::max(joint.partA, joint.partB);
        const std::size_t key = a * count + b;
        if (joint.severed())
            continue;
        if (joint.brokenCount() == 0) {
            if (state[key] != kSeamDamaged)
                state[key] = kSeamIntact;
        } else {
            state[key] = kSeamDamaged;
        }
    }
}

void Assembly::correctInterpenetration() {
    const std::size_t count = _parts.size();
    if (count < 2)
        return;

    std::vector<char> seam;
    classifySeams(seam);
    std::vector<VolumeGrid> grids(count);
    std::vector<ChVector3d> lo(count), hi(count);
    for (std::size_t i = 0; i < count; ++i) {
        grids[i].build(_parts[i]);
        _parts[i].bounds(lo[i], hi[i]);
    }

    std::vector<char> closingPairs(count * count, 0);
    for (std::size_t i = 0; i < count; ++i) {
        for (std::size_t j = i + 1; j < count; ++j) {
            if (aabbOverlap(lo[i], hi[i], lo[j], hi[j]))
                closingPairs[i * count + j] = 1;
        }
    }

    std::vector<std::vector<ChVector3d>> liveWeld(count * count);
    for (const auto& joint : _joints) {
        if (joint.partA >= count || joint.partB >= count || joint.partA == joint.partB)
            continue;
        const std::size_t key = std::min(joint.partA, joint.partB) * count +
                                std::max(joint.partA, joint.partB);
        if (joint.severed() || seam[key] != kSeamDamaged)
            continue;
        for (const auto& link : joint.links) {
            if (link.broken)
                continue;
            if (link.partA < count && link.nodeA < _parts[link.partA].nodeCount())
                liveWeld[key].push_back(_parts[link.partA].nodePosition(link.nodeA));
            if (link.partB < count && link.nodeB < _parts[link.partB].nodeCount())
                liveWeld[key].push_back(_parts[link.partB].nodePosition(link.nodeB));
        }
    }

    auto closestSkin = [&](const Part& host, const ChVector3d& p, ChVector3d& hit, ChVector3d& normal,
                           int& i0, int& i1, int& i2, double& wa, double& wb, double& wc) {
        const auto& faces = host.triangles();
        double best = 1.0e30;
        i0 = i1 = i2 = -1;
        for (const auto& tri : faces) {
            const ChVector3d a = host.nodePosition(tri[0]);
            const ChVector3d b = host.nodePosition(tri[1]);
            const ChVector3d c = host.nodePosition(tri[2]);
            const ChVector3d q = closestPointOnTriangle(p, a, b, c);
            const double d = (p - q).Length();
            if (d >= best)
                continue;
            ChVector3d n = (b - a).Cross(c - a);
            if (n.Length() < 1.0e-12)
                continue;
            n.Normalize();
            best = d;
            hit = q;
            normal = n;
            i0 = tri[0];
            i1 = tri[1];
            i2 = tri[2];
            barycentric(q, a, b, c, wa, wb, wc);
        }
        if (i0 < 0)
            return false;
        wa = std::max(0.0, wa);
        wb = std::max(0.0, wb);
        wc = std::max(0.0, wc);
        const double sum = wa + wb + wc;
        if (sum > 1.0e-12) {
            wa /= sum;
            wb /= sum;
            wc /= sum;
        } else {
            wa = wb = wc = 1.0 / 3.0;
        }
        if (normal.Dot(hit - p) < 0.0)
            normal = -normal;
        return true;
    };

    auto killClosing = [&](Part& guest, Part& host, const ChVector3d& axis) {
        ChVector3d nHat = axis;
        if (nHat.Length2() < 1.0e-16)
            return;
        nHat.Normalize();
        for (std::size_t n = 0; n < guest.nodeCount(); ++n) {
            ChVector3d vel = guest.nodeVelocity(n);
            const double vn = vel.Dot(nHat);
            if (vn < 0.0)
                guest.setNodeKinematics(n, guest.nodePosition(n), vel - nHat * vn);
        }
        for (std::size_t n = 0; n < host.nodeCount(); ++n) {
            ChVector3d vel = host.nodeVelocity(n);
            const double vn = vel.Dot(nHat);
            if (vn > 0.0)
                host.setNodeKinematics(n, host.nodePosition(n), vel - nHat * vn);
        }
    };

    auto eject = [&](std::size_t guestId, std::size_t hostId) {
        Part& guest = _parts[guestId];
        Part& host = _parts[hostId];
        if (!host.hasVolume())
            return;
        const std::size_t key = std::min(guestId, hostId) * count + std::max(guestId, hostId);
        const double skin = guest.contactSkin() + host.contactSkin();
        const double weldR = std::max({0.008, 3.0 * guest.contactSkin(), 3.0 * host.contactSkin()});
        const double seamBand = std::max(weldR, 2.0 * skin);
        const double maxMove = std::max(0.005, 2.0 * skin);
        const bool welded = seam[key] != kSeamNone;
        if (seam[key] == kSeamIntact)
            return;
        const std::vector<ChVector3d>* weldPts =
            liveWeld[key].empty() ? nullptr : &liveWeld[key];
        const auto& faces = host.triangles();
        const auto& area = guest.contactAreas();

        int buried = 0;
        for (std::size_t n = 0; n < guest.nodeCount(); ++n) {
            if (!nodeCorrectable(guest, n))
                continue;
            const ChVector3d now = guest.nodePosition(n);
            if (now.x() < lo[hostId].x() || now.x() > hi[hostId].x() || now.y() < lo[hostId].y() ||
                now.y() > hi[hostId].y() || now.z() < lo[hostId].z() || now.z() > hi[hostId].z())
                continue;
            if (grids[hostId].contains(now))
                ++buried;
        }
        if (buried == 0)
            return;

        if (seam[key] == kSeamNone && buried * 6 >= static_cast<int>(guest.nodeCount())) {
            const ChVector3d push = aabbSeparate(lo[guestId], hi[guestId], lo[hostId], hi[hostId], skin);
            if (push.Length2() > 1.0e-16) {
                ChVector3d nHat = push;
                nHat.Normalize();
                for (std::size_t n = 0; n < guest.nodeCount(); ++n) {
                    ChVector3d vel = guest.nodeVelocity(n);
                    const double vn = vel.Dot(nHat);
                    if (vn < 0.0)
                        vel -= nHat * vn;
                    const ChVector3d pos = nodeCorrectable(guest, n) ? guest.nodePosition(n) + push
                                                                     : guest.nodePosition(n);
                    guest.setNodeKinematics(n, pos, vel);
                }
                for (std::size_t n = 0; n < host.nodeCount(); ++n) {
                    ChVector3d vel = host.nodeVelocity(n);
                    const double vn = vel.Dot(nHat);
                    if (vn > 0.0)
                        vel -= nHat * vn;
                    host.setNodeKinematics(n, host.nodePosition(n), vel);
                }
                killClosing(guest, host, nHat);
            }
            return;
        }

        auto clampTo = [&](const ChVector3d& from, const ChVector3d& to) {
            const ChVector3d d = to - from;
            const double len = d.Length();
            if (len > maxMove && len > 1.0e-12)
                return from + d * (maxMove / len);
            return to;
        };

        ChVector3d meanDelta(0, 0, 0);
        int meanCount = 0;
        std::vector<char> interior(guest.nodeCount(), 0);

        for (std::size_t n = 0; n < guest.nodeCount(); ++n) {
            if (!nodeCorrectable(guest, n))
                continue;
            const ChVector3d now = guest.nodePosition(n);
            if (now.x() < lo[hostId].x() || now.x() > hi[hostId].x() || now.y() < lo[hostId].y() ||
                now.y() > hi[hostId].y() || now.z() < lo[hostId].z() || now.z() > hi[hostId].z())
                continue;
            if (!grids[hostId].contains(now))
                continue;

            const bool onSkin = !guest.hasVolume() || (n < area.size() && area[n] > 0.0);
            if (!onSkin) {
                interior[n] = 1;
                continue;
            }

            const ChVector3d prev = guest.savedPosition(n);
            const bool prevInside = grids[hostId].contains(prev);
            ChVector3d hit, normal;
            int i0 = -1, i1 = -1, i2 = -1;
            double wa = 1.0 / 3.0, wb = 1.0 / 3.0, wc = 1.0 / 3.0;
            bool found = false;
            if (!prevInside && (prev - now).Length2() > 1.0e-16) {
                double bestT = 2.0;
                for (const auto& tri : faces) {
                    const ChVector3d a = host.nodePosition(tri[0]);
                    const ChVector3d b = host.nodePosition(tri[1]);
                    const ChVector3d c = host.nodePosition(tri[2]);
                    double t = 0.0;
                    ChVector3d q, nrm;
                    if (!segmentHitsTriangle(prev, now, a, b, c, t, q, nrm))
                        continue;
                    if (t >= bestT)
                        continue;
                    bestT = t;
                    hit = q;
                    normal = nrm;
                    i0 = tri[0];
                    i1 = tri[1];
                    i2 = tri[2];
                    barycentric(q, a, b, c, wa, wb, wc);
                    found = true;
                }
                if (found) {
                    if (normal.Dot(prev - hit) < 0.0)
                        normal = -normal;
                    wa = std::max(0.0, wa);
                    wb = std::max(0.0, wb);
                    wc = std::max(0.0, wc);
                    const double sum = wa + wb + wc;
                    if (sum > 1.0e-12) {
                        wa /= sum;
                        wb /= sum;
                        wc /= sum;
                    } else {
                        wa = wb = wc = 1.0 / 3.0;
                    }
                }
            }
            if (!found && !closestSkin(host, now, hit, normal, i0, i1, i2, wa, wb, wc))
                continue;

            const double depth = (now - hit).Length();
            if (welded && depth <= seamBand)
                continue;
            if (weldPts && nearPoints(*weldPts, now, weldR))
                continue;

            const ChVector3d target = hit + normal * std::max(skin, 1.0e-4);
            const ChVector3d pos = clampTo(now, target);
            ChVector3d faceVel = host.nodeVelocity(static_cast<std::size_t>(i0)) * wa +
                                 host.nodeVelocity(static_cast<std::size_t>(i1)) * wb +
                                 host.nodeVelocity(static_cast<std::size_t>(i2)) * wc;
            ChVector3d vel = guest.nodeVelocity(n);
            const ChVector3d rel = vel - faceVel;
            const double vn = rel.Dot(normal);
            if (vn < 0.0 && !welded)
                vel -= normal * vn;

            const double mg = n < guest.nodeMass().size() ? guest.nodeMass()[n] : 0.0;
            if (!welded && vn < 0.0 && mg > 0.0 && i0 >= 0) {
                const ChVector3d impulse = normal * (mg * (-vn));
                auto kick = [&](int idx, double w) {
                    if (idx < 0 || w <= 0.0)
                        return;
                    const std::size_t k = static_cast<std::size_t>(idx);
                    if (k >= host.nodeCount() || !nodeCorrectable(host, k))
                        return;
                    const double mh = k < host.nodeMass().size() ? host.nodeMass()[k] : 0.0;
                    if (mh <= 1.0e-12)
                        return;
                    ChVector3d dv = impulse * (-w / mh);
                    const double cap = std::abs(vn);
                    if (dv.Length() > cap && cap > 0.0)
                        dv *= cap / dv.Length();
                    host.setNodeKinematics(k, host.nodePosition(k), host.nodeVelocity(k) + dv);
                };
                kick(i0, wa);
                kick(i1, wb);
                kick(i2, wc);
            }
            meanDelta += pos - now;
            ++meanCount;
            guest.setNodeKinematics(n, pos, vel);
        }

        ChVector3d shift = meanCount > 0 ? meanDelta / static_cast<double>(meanCount)
                                         : (guest.center() - host.center());
        if (meanCount == 0) {
            if (welded)
                return;
            if (shift.Length() > 1.0e-12)
                shift = shift.GetNormalized() * maxMove;
        } else if (shift.Length() > maxMove && shift.Length() > 1.0e-12)
            shift *= maxMove / shift.Length();

        if (shift.Length2() >= 1.0e-16) {
            for (std::size_t n = 0; n < guest.nodeCount(); ++n) {
                if (!interior[n])
                    continue;
                const ChVector3d now = guest.nodePosition(n);
                guest.setNodeKinematics(n, now + shift, guest.nodeVelocity(n));
            }
        }
        if (!welded)
            killClosing(guest, host, guest.center() - host.center());
    };

    for (int pass = 0; pass < 2; ++pass) {
        for (std::size_t i = 0; i < count; ++i) {
            grids[i].build(_parts[i]);
            _parts[i].bounds(lo[i], hi[i]);
        }
        for (std::size_t i = 0; i < count; ++i) {
            for (std::size_t j = i + 1; j < count; ++j) {
                if (!aabbOverlap(lo[i], hi[i], lo[j], hi[j]))
                    continue;
                eject(i, j);
                eject(j, i);
            }
        }
    }

    auto massCenter = [](const Part& part) {
        ChVector3d c(0, 0, 0);
        double m = 0.0;
        for (std::size_t n = 0; n < part.nodeCount(); ++n) {
            const double mn = n < part.nodeMass().size() ? part.nodeMass()[n] : 0.0;
            c += part.nodePosition(n) * mn;
            m += mn;
        }
        return m > 0.0 ? c / m : part.center();
    };
    for (std::size_t i = 0; i < count; ++i) {
        for (std::size_t j = i + 1; j < count; ++j) {
            if (seam[i * count + j] != kSeamNone || !closingPairs[i * count + j])
                continue;
            killClosing(_parts[i], _parts[j], massCenter(_parts[i]) - massCenter(_parts[j]));
        }
    }
}

bool Assembly::maintainVolume(double dt) {
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
    _volumeWait += dt;
    bool changed = false;
    for (auto& part : _parts)
        changed = part.maintainVolume() || changed;
    if (changed)
        inheritWeldNodes();
    return changed;
}

int Assembly::impactSubsteps(const Environment& environment, double dt) const {
    if (!environment.floorEnabled)
        return 1;
    int steps = 1;
    const double h = std::max(dt, 1.0e-12);
    for (const auto& part : _parts) {
        const double offset = part.contactOffset();
        const auto& area = part.contactAreas();
        for (int index : part.boundaryNodes()) {
            if (index < 0 || static_cast<std::size_t>(index) >= part.nodeCount())
                continue;
            if (static_cast<std::size_t>(index) < part.xyzNodes().size() &&
                part.xyzNodes()[static_cast<std::size_t>(index)] &&
                part.xyzNodes()[static_cast<std::size_t>(index)]->IsFixed())
                continue;
            const double gap = part.nodePosition(static_cast<std::size_t>(index)).y() -
                               (environment.floorHeight + offset);
            const double closing = -part.nodeVelocity(static_cast<std::size_t>(index)).y();
            if (closing <= 0.0)
                continue;
            const double patchArea =
                static_cast<std::size_t>(index) < area.size() ? std::max(0.0, area[static_cast<std::size_t>(index)])
                                                              : 0.0;
            const double patch = std::sqrt(patchArea);
            const double reach = std::max({offset, patch, part.contactSkin(), 1.0e-6});
            if (gap > closing * h + reach)
                continue;
            int need = static_cast<int>(std::ceil(closing * h / reach));
            if (gap < 0.0)
                need = std::max(need, static_cast<int>(std::ceil((-gap) / reach)));
            // Period samples on the outer steps that will reach the floor
            // (this step or the next). A millimetre rest hover with gravity
            // closing does not qualify; a 50 m/s centimetre-per-step hit does.
            if (gap > 0.0 && gap <= 2.0 * closing * h) {
                const double mass = static_cast<std::size_t>(index) < part.nodeMass().size()
                                        ? part.nodeMass()[static_cast<std::size_t>(index)]
                                        : 0.0;
                const double k = continuumGains(part.material().youngsModulus, patchArea, mass).stiffness;
                if (mass > 1.0e-12 && k > 0.0) {
                    const double omega = std::sqrt(k / mass);
                    need = std::max(need, static_cast<int>(std::ceil(omega * h / 3.14159265358979323846)));
                    need = std::min(need, std::max(1, static_cast<int>(std::ceil(h / kHhtMinStep))));
                }
            }
            steps = std::max(steps, std::max(1, need));
        }
    }
    const int cflCap = std::max(1, static_cast<int>(std::ceil(h / kHhtMinStep)));
    return std::min(steps, cflCap);
}

void Assembly::setContactStep(double h) {
    const double step = std::max(h, 1.0e-6);
    if (_contactWorld)
        _contactWorld->stepSize = step;
    for (auto& part : _parts)
        part.setImpactStep(step);
}

void Assembly::inheritWeldNodes() {
    if (!_system)
        return;
    for (const auto& weld : _weldDefs) {
        if (weld.partA >= _parts.size() || weld.partB >= _parts.size())
            continue;
        Joint* joint = nullptr;
        for (auto& existing : _joints) {
            if (existing.kind == JointKind::Weld && existing.partA == weld.partA && existing.partB == weld.partB) {
                joint = &existing;
                break;
            }
        }
        if (!joint || joint->links.empty())
            continue;
        const double breakForce = joint->links.front().breakForce;
        auto linked = [&](std::size_t nodeA, std::size_t nodeB) {
            for (const auto& link : joint->links) {
                if (link.broken)
                    continue;
                if (link.nodeA == nodeA && link.nodeB == nodeB)
                    return true;
            }
            return false;
        };
        const Part& a = _parts[weld.partA];
        const Part& b = _parts[weld.partB];
        std::vector<std::pair<std::size_t, std::size_t>> pairs;
        for (const auto& link : joint->links) {
            if (!link.broken)
                pairs.emplace_back(link.nodeA, link.nodeB);
        }
        auto nearBoundary = [](const Part& part, const ChVector3d& point) {
            std::size_t best = part.nodeCount();
            double bestDistance = 1.0e-3;
            for (int index : part.boundaryNodes()) {
                if (index < 0)
                    continue;
                const double distance = (part.nodePosition(static_cast<std::size_t>(index)) - point).Length();
                if (distance < bestDistance) {
                    bestDistance = distance;
                    best = static_cast<std::size_t>(index);
                }
            }
            return best;
        };
        const std::size_t count = pairs.size();
        for (std::size_t i = 0; i < count; ++i) {
            for (std::size_t j = i + 1; j < count; ++j) {
                const ChVector3d midA =
                    0.5 * (a.nodePosition(pairs[i].first) + a.nodePosition(pairs[j].first));
                const std::size_t midNodeA = nearBoundary(a, midA);
                if (midNodeA >= a.nodeCount() || midNodeA == pairs[i].first || midNodeA == pairs[j].first)
                    continue;
                bool has = false;
                for (const auto& link : joint->links)
                    has = has || (!link.broken && link.nodeA == midNodeA);
                if (has)
                    continue;
                const ChVector3d midB =
                    0.5 * (b.nodePosition(pairs[i].second) + b.nodePosition(pairs[j].second));
                std::size_t midNodeB = nearBoundary(b, midB);
                if (midNodeB >= b.nodeCount() || linked(midNodeA, midNodeB))
                    midNodeB = (midA - b.nodePosition(pairs[i].second)).Length() <=
                                       (midA - b.nodePosition(pairs[j].second)).Length()
                                   ? pairs[i].second
                                   : pairs[j].second;
                if (linked(midNodeA, midNodeB))
                    continue;
                createLink(*_system, *joint, weld.partA, midNodeA, weld.partB, midNodeB, breakForce, 0.0, false,
                           false);
            }
        }
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
    for (auto& joint : _joints) {
        for (auto& entry : joint.links) {
            if (entry.broken)
                continue;
            auto liveFea = [](const Part& part, std::size_t index) {
                if (index >= part.nodeCount() || part.nodeRigid(index))
                    return false;
                if (part.hasVolume())
                    return part.nodeOnLiveTet(index);
                return part.nodeOnLiveShell(index);
            };
            if (entry.partA >= _parts.size() || entry.partB >= _parts.size())
                continue;
            if (!liveFea(_parts[entry.partA], entry.nodeA) || !liveFea(_parts[entry.partB], entry.nodeB)) {
                auto drop = [](const auto& link) {
                    if (!link)
                        return;
                    link->SetBroken(true);
                    link->SetDisabled(true);
                };
                drop(entry.link);
                drop(entry.pointLink);
                drop(entry.offsetLink);
                drop(entry.frameLink);
                entry.broken = true;
                changed = true;
            }
        }
    }
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
        if (!joint->links.empty()) {
            JointLink& created = joint->links.back();
            created.breakPeel = skipped.breakPeel;
            created.criticalOpening = skipped.criticalOpening;
            created.peelAxis = skipped.peelAxis;
            const chrono::ChVector3d sep = _parts[skipped.partB].nodePosition(skipped.nodeB) -
                                           _parts[skipped.partA].nodePosition(skipped.nodeA);
            created.restGap = sep.Length();
        }
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
        entry.volume = part.hasVolume();
        entry.flatTriangles = part.skinTriangleCount();
        entry.tetCount = part.tetCount();
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
            auto p = part.nodePosition(n);
            if (!std::isfinite(p.x()) || !std::isfinite(p.y()) || !std::isfinite(p.z()))
                p = part.restPosition(n);
            if (!std::isfinite(p.x()) || !std::isfinite(p.y()) || !std::isfinite(p.z()))
                p = chrono::ChVector3d(0, 0, 0);
            out.setPosition(n, p.x(), p.y(), p.z());
        }
        if (part.torn()) {
            out.quads = part.quads();
            out.triangles = part.triangles();
        }
        part.evaluateStress(out.stress);
        part.evaluatePlasticStrain(out.plasticStrain);
        part.evaluateTemperature(out.temperature);
        out.nodeTemperature.resize(part.nodeCount());
        for (std::size_t n = 0; n < part.nodeCount(); ++n)
            out.nodeTemperature[n] = static_cast<float>(part.nodeTemperature(n));
        for (float s : out.stress)
            frame.maxStress = std::max(frame.maxStress, s);
    }
    _heat.writeAir(frame.air);

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

int Assembly::invertedTetCount() const {
    int inverted = 0;
    for (const auto& part : _parts)
        inverted += part.invertedTetCount();
    return inverted;
}

bool Assembly::hasFreeSolve() const {
    for (const auto& part : _parts) {
        if (part.hasFreeNodes())
            return true;
    }
    return false;
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
