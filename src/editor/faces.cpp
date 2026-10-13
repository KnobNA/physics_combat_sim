#include "editor/faces.hpp"

#include "chrono/core/ChQuaternion.h"
#include "chrono/core/ChRotation.h"

#include <algorithm>
#include <cmath>
#include <cstdio>
#include <limits>
#include <map>
#include <optional>
#include <set>
#include <utility>

namespace rcs {
namespace {

using chrono::ChVector3d;
using chrono::ChQuaterniond;

constexpr double kPi = 3.14159265358979323846;
constexpr double kCoplanar = 1.0e-6;
constexpr double kNormalDot = 0.999;
constexpr double kHitEps = 1.0e-6;

ChVector3d eulerFromAxes(const ChVector3d& u, const ChVector3d& v);

ChVector3d patchCorner(const PlatePatch& patch, double su, double sv) {
    ChVector3d u = patch.uAxis;
    ChVector3d v = patch.vAxis;
    if (u.Length() > 1.0e-12)
        u.Normalize();
    if (v.Length() > 1.0e-12)
        v.Normalize();
    return patch.origin + u * (su * patch.uSize) + v * (sv * patch.vSize);
}

void sourceFromPatch(PartDefinition& part) {
    if (part.patch.uSize <= 1.0e-9 || part.patch.vSize <= 1.0e-9)
        return;
    ChVector3d u = part.patch.uAxis;
    ChVector3d v = part.patch.vAxis;
    if (u.Length() > 1.0e-12)
        u.Normalize();
    if (v.Length() > 1.0e-12)
        v.Normalize();
    ChVector3d n = u.Cross(v);
    if (n.Length() < 1.0e-12)
        n = {0, 1, 0};
    else
        n.Normalize();
    const double h = 0.5 * std::max(part.thickness, 1.0e-6);
    auto corner = [&](double su, double sv, double side) {
        return patchCorner(part.patch, su, sv) + n * (side * h);
    };
    part.source.vertices = {corner(0.0, 0.0, -1.0), corner(1.0, 0.0, -1.0), corner(1.0, 1.0, -1.0),
                            corner(0.0, 1.0, -1.0), corner(0.0, 0.0, 1.0),  corner(1.0, 0.0, 1.0),
                            corner(1.0, 1.0, 1.0),  corner(0.0, 1.0, 1.0)};
    part.source.faces = {AuthoredFace{{0, 3, 2, 1}}, AuthoredFace{{4, 5, 6, 7}},
                         AuthoredFace{{0, 1, 5, 4}}, AuthoredFace{{1, 2, 6, 5}},
                         AuthoredFace{{2, 3, 7, 6}}, AuthoredFace{{3, 0, 4, 7}}};
    // Leave triangles empty so a plate stays a Reissner sheet until the first splice.
    part.source.triangles.clear();
    part.euler = eulerFromAxes(part.patch.uAxis, part.patch.vAxis);
}

bool plateNeedsShellFaces(const PartDefinition& part) {
    if (part.patch.uSize <= 1.0e-9 || part.patch.vSize <= 1.0e-9)
        return false;
    if (part.thickness <= 1.0e-9)
        return false;
    if (part.source.faces.empty())
        return part.source.triangles.empty();
    return part.source.faces.size() == 1 && part.source.vertices.size() <= 4;
}

ChVector3d faceNormal(const SurfaceMesh& mesh, const AuthoredFace& face) {
    if (face.loop.size() < 3)
        return {0, 1, 0};
    const ChVector3d& a = mesh.vertices[static_cast<std::size_t>(face.loop[0])];
    ChVector3d n{0, 0, 0};
    for (std::size_t i = 1; i + 1 < face.loop.size(); ++i) {
        const ChVector3d& b = mesh.vertices[static_cast<std::size_t>(face.loop[i])];
        const ChVector3d& c = mesh.vertices[static_cast<std::size_t>(face.loop[i + 1])];
        n += (b - a).Cross(c - a);
    }
    if (n.Length() < 1.0e-12)
        return {0, 1, 0};
    n.Normalize();
    return n;
}

ChVector3d triangleNormal(const SurfaceMesh& mesh, const std::array<int, 3>& tri) {
    const ChVector3d& a = mesh.vertices[static_cast<std::size_t>(tri[0])];
    const ChVector3d& b = mesh.vertices[static_cast<std::size_t>(tri[1])];
    const ChVector3d& c = mesh.vertices[static_cast<std::size_t>(tri[2])];
    ChVector3d n = (b - a).Cross(c - a);
    if (n.Length() < 1.0e-12)
        return {0, 0, 0};
    n.Normalize();
    return n;
}

bool extractBoundaryLoop(const SurfaceMesh& mesh, const std::vector<int>& tris, AuthoredFace& out) {
    std::map<std::pair<int, int>, int> directed;
    for (int t : tris) {
        const auto& tri = mesh.triangles[static_cast<std::size_t>(t)];
        for (int k = 0; k < 3; ++k) {
            const int a = tri[static_cast<std::size_t>(k)];
            const int b = tri[static_cast<std::size_t>((k + 1) % 3)];
            directed[{a, b}]++;
        }
    }
    std::map<int, int> next;
    for (const auto& edge : directed) {
        if (directed[{edge.first.second, edge.first.first}] > 0)
            continue;
        if (next.count(edge.first.first))
            return false;
        next[edge.first.first] = edge.first.second;
    }
    if (next.empty())
        return false;
    out.loop.clear();
    int start = next.begin()->first;
    int cursor = start;
    for (int guard = 0; guard < 4096; ++guard) {
        out.loop.push_back(cursor);
        const auto found = next.find(cursor);
        if (found == next.end())
            return false;
        cursor = found->second;
        if (cursor == start)
            break;
    }
    return out.loop.size() >= 3 && cursor == start;
}

void seedFacesFromTriangles(SurfaceMesh& mesh) {
    if (!mesh.faces.empty() || mesh.triangles.empty())
        return;
    const std::size_t count = mesh.triangles.size();
    std::vector<std::vector<int>> adj(count);
    for (std::size_t i = 0; i < count; ++i) {
        for (std::size_t j = i + 1; j < count; ++j) {
            int shared = 0;
            for (int a : mesh.triangles[i])
                for (int b : mesh.triangles[j])
                    shared += a == b ? 1 : 0;
            if (shared == 2) {
                adj[i].push_back(static_cast<int>(j));
                adj[j].push_back(static_cast<int>(i));
            }
        }
    }
    std::vector<char> used(count, 0);
    for (std::size_t seed = 0; seed < count; ++seed) {
        if (used[seed])
            continue;
        const ChVector3d n0 = triangleNormal(mesh, mesh.triangles[seed]);
        if (n0.Length() < 0.5) {
            used[seed] = 1;
            mesh.faces.push_back(AuthoredFace{{mesh.triangles[seed][0], mesh.triangles[seed][1],
                                               mesh.triangles[seed][2]}});
            continue;
        }
        const ChVector3d p0 = mesh.vertices[static_cast<std::size_t>(mesh.triangles[seed][0])];
        std::vector<int> group;
        std::vector<int> stack{static_cast<int>(seed)};
        used[seed] = 1;
        while (!stack.empty()) {
            const int t = stack.back();
            stack.pop_back();
            group.push_back(t);
            const ChVector3d nt = triangleNormal(mesh, mesh.triangles[static_cast<std::size_t>(t)]);
            for (int n : adj[static_cast<std::size_t>(t)]) {
                if (used[static_cast<std::size_t>(n)])
                    continue;
                const ChVector3d nn = triangleNormal(mesh, mesh.triangles[static_cast<std::size_t>(n)]);
                if (nt.Dot(nn) < kNormalDot)
                    continue;
                const ChVector3d p = mesh.vertices[static_cast<std::size_t>(mesh.triangles[static_cast<std::size_t>(n)][0])];
                if (std::abs((p - p0).Dot(n0)) > kCoplanar)
                    continue;
                used[static_cast<std::size_t>(n)] = 1;
                stack.push_back(n);
            }
        }
        AuthoredFace face;
        if (!extractBoundaryLoop(mesh, group, face)) {
            for (int t : group)
                mesh.faces.push_back(AuthoredFace{{mesh.triangles[static_cast<std::size_t>(t)][0],
                                                   mesh.triangles[static_cast<std::size_t>(t)][1],
                                                   mesh.triangles[static_cast<std::size_t>(t)][2]}});
        } else {
            mesh.faces.push_back(std::move(face));
        }
    }
}

int addVertex(SurfaceMesh& mesh, const ChVector3d& point) {
    for (std::size_t i = 0; i < mesh.vertices.size(); ++i) {
        if ((mesh.vertices[i] - point).Length() < 1.0e-9)
            return static_cast<int>(i);
    }
    mesh.vertices.push_back(point);
    return static_cast<int>(mesh.vertices.size() - 1);
}

void shiftMesh(SurfaceMesh& mesh, const ChVector3d& delta) {
    for (auto& v : mesh.vertices)
        v += delta;
}

void rotateMesh(SurfaceMesh& mesh, const ChQuaterniond& rotation, const ChVector3d& pivot) {
    for (auto& v : mesh.vertices)
        v = pivot + rotation.Rotate(v - pivot);
}

ChVector3d eulerFromAxes(const ChVector3d& u, const ChVector3d& v) {
    ChVector3d e1 = u;
    ChVector3d e2 = v;
    if (e1.Length() < 1.0e-12)
        e1 = {1, 0, 0};
    e1.Normalize();
    ChVector3d n = e1.Cross(e2);
    if (n.Length() < 1.0e-12)
        n = {0, 1, 0};
    n.Normalize();
    e2 = n.Cross(e1);
    const double sy = -e1.z();
    const double cy = std::sqrt(std::max(0.0, 1.0 - sy * sy));
    ChVector3d euler;
    if (cy > 1.0e-8) {
        euler.x() = std::atan2(e2.z(), n.z());
        euler.y() = std::atan2(sy, cy);
        euler.z() = std::atan2(e1.y(), e1.x());
    } else {
        euler.x() = std::atan2(-n.y(), e2.y());
        euler.y() = std::atan2(sy, cy);
        euler.z() = 0.0;
    }
    return euler;
}

bool rayHitTriangle(const ChVector3d& origin, const ChVector3d& dir, const ChVector3d& a,
                    const ChVector3d& b, const ChVector3d& c, double& t) {
    const ChVector3d edge1 = b - a;
    const ChVector3d edge2 = c - a;
    const ChVector3d p = dir.Cross(edge2);
    const double det = edge1.Dot(p);
    if (std::abs(det) < 1.0e-12)
        return false;
    const double inv = 1.0 / det;
    const ChVector3d s = origin - a;
    const double u = s.Dot(p) * inv;
    if (u < -1.0e-8 || u > 1.0 + 1.0e-8)
        return false;
    const ChVector3d q = s.Cross(edge1);
    const double v = dir.Dot(q) * inv;
    if (v < -1.0e-8 || u + v > 1.0 + 1.0e-8)
        return false;
    t = edge2.Dot(q) * inv;
    return t > 1.0e-8;
}

bool pointInLoop(const SurfaceMesh& mesh, const AuthoredFace& face, const ChVector3d& point) {
    if (face.loop.size() < 3)
        return false;
    const ChVector3d n = faceNormal(mesh, face);
    ChVector3d e1 = mesh.vertices[static_cast<std::size_t>(face.loop[1])] -
                    mesh.vertices[static_cast<std::size_t>(face.loop[0])];
    if (e1.Length() < 1.0e-12)
        e1 = n.Cross(ChVector3d(0, 1, 0));
    if (e1.Length() < 1.0e-12)
        e1 = {1, 0, 0};
    e1.Normalize();
    const ChVector3d e2 = n.Cross(e1);
    const ChVector3d origin = mesh.vertices[static_cast<std::size_t>(face.loop[0])];
    auto proj = [&](const ChVector3d& p) {
        const ChVector3d d = p - origin;
        return std::pair<double, double>{d.Dot(e1), d.Dot(e2)};
    };
    const auto q = proj(point);
    int cuts = 0;
    for (std::size_t i = 0; i < face.loop.size(); ++i) {
        const auto a = proj(mesh.vertices[static_cast<std::size_t>(face.loop[i])]);
        const auto b = proj(mesh.vertices[static_cast<std::size_t>(face.loop[(i + 1) % face.loop.size()])]);
        if ((a.second > q.second) == (b.second > q.second))
            continue;
        const double t = (q.second - a.second) / (b.second - a.second);
        const double x = a.first + t * (b.first - a.first);
        if (x > q.first)
            ++cuts;
    }
    return (cuts % 2) == 1;
}

int remapPartIndex(int index, const std::vector<char>& gone) {
    if (index < 0)
        return index;
    if (static_cast<std::size_t>(index) >= gone.size() || gone[static_cast<std::size_t>(index)])
        return -1;
    int shift = 0;
    for (int i = 0; i < index; ++i)
        shift += gone[static_cast<std::size_t>(i)] ? 1 : 0;
    return index - shift;
}

int remapGroupIndex(int index, const std::vector<char>& gone) {
    return remapPartIndex(index, gone);
}

}  // namespace

void tessellateAuthoredFaces(SurfaceMesh& mesh) {
    mesh.triangles.clear();
    for (const auto& face : mesh.faces) {
        if (face.loop.size() < 3)
            continue;
        for (std::size_t i = 1; i + 1 < face.loop.size(); ++i)
            mesh.triangles.push_back({face.loop[0], face.loop[i], face.loop[i + 1]});
    }
}

void ensureAuthoredFaces(PartDefinition& part) {
    if (plateNeedsShellFaces(part)) {
        sourceFromPatch(part);
        return;
    }
    if (!part.source.faces.empty())
        return;
    if (part.source.empty())
        sourceFromPatch(part);
    if (part.source.faces.empty() && !part.source.triangles.empty())
        seedFacesFromTriangles(part.source);
}

ChVector3d partCentroid(const PartDefinition& part) {
    ChVector3d sum{0, 0, 0};
    int count = 0;
    auto add = [&](const ChVector3d& p) {
        sum += p;
        ++count;
    };
    if (!part.source.vertices.empty()) {
        for (const auto& v : part.source.vertices)
            add(v);
    } else if (!part.volume.vertices.empty()) {
        for (const auto& v : part.volume.vertices)
            add(v);
    } else if (part.patch.uSize > 1.0e-9 && part.patch.vSize > 1.0e-9) {
        add(patchCorner(part.patch, 0.0, 0.0));
        add(patchCorner(part.patch, 1.0, 0.0));
        add(patchCorner(part.patch, 1.0, 1.0));
        add(patchCorner(part.patch, 0.0, 1.0));
    } else {
        add(part.patch.origin);
    }
    return count ? sum / double(count) : ChVector3d{0, 0, 0};
}

ChVector3d partEuler(const PartDefinition& part) {
    if (!part.source.empty())
        return part.euler;
    return eulerFromAxes(part.patch.uAxis, part.patch.vAxis);
}

void applyPartTranslation(PartDefinition& part, const ChVector3d& delta) {
    if (delta.Length() < 1.0e-12)
        return;
    part.patch.origin += delta;
    shiftMesh(part.source, delta);
    shiftMesh(part.surface, delta);
    for (auto& v : part.volume.vertices)
        v += delta;
}

void applyPartRotation(PartDefinition& part, const ChQuaterniond& rotation, const ChVector3d& pivot) {
    auto spin = [&](ChVector3d& p) { p = pivot + rotation.Rotate(p - pivot); };
    spin(part.patch.origin);
    part.patch.uAxis = rotation.Rotate(part.patch.uAxis);
    part.patch.vAxis = rotation.Rotate(part.patch.vAxis);
    rotateMesh(part.source, rotation, pivot);
    rotateMesh(part.surface, rotation, pivot);
    for (auto& v : part.volume.vertices)
        spin(v);
    const ChQuaterniond old = chrono::QuatFromAngleZ(part.euler.z()) * chrono::QuatFromAngleY(part.euler.y()) *
                              chrono::QuatFromAngleX(part.euler.x());
    const ChQuaterniond next = rotation * old;
    part.euler = eulerFromAxes(next.Rotate(ChVector3d(1, 0, 0)), next.Rotate(ChVector3d(0, 0, 1)));
}

ChVector3d toWorld(const ChVector3d& authored, const ChQuaterniond& rotation, const ChVector3d& origin) {
    return origin + rotation.Rotate(authored);
}

void faceLoopWorld(const PartDefinition& part, std::size_t faceIndex, const ChQuaterniond& rotation,
                   const ChVector3d& origin, std::vector<ChVector3d>& out) {
    out.clear();
    if (faceIndex >= part.source.faces.size())
        return;
    for (int index : part.source.faces[faceIndex].loop) {
        if (index < 0 || static_cast<std::size_t>(index) >= part.source.vertices.size())
            continue;
        out.push_back(toWorld(part.source.vertices[static_cast<std::size_t>(index)], rotation, origin));
    }
}

std::optional<FaceId> pickAuthoredFace(const AssemblyDefinition& definition,
                                       const std::vector<std::size_t>& parts,
                                       const ChVector3d& origin, const ChVector3d& direction,
                                       const ChQuaterniond& worldRotation, const ChVector3d& worldOrigin,
                                       ChVector3d* hitPoint) {
    double best = 1.0e9;
    std::optional<FaceId> hit;
    ChVector3d point{0, 0, 0};
    for (std::size_t partIndex : parts) {
        if (partIndex >= definition.parts.size())
            continue;
        const PartDefinition& part = definition.parts[partIndex];
        for (std::size_t f = 0; f < part.source.faces.size(); ++f) {
            const auto& face = part.source.faces[f];
            if (face.loop.size() < 3)
                continue;
            auto vert = [&](int index, ChVector3d& out) {
                if (index < 0 || static_cast<std::size_t>(index) >= part.source.vertices.size())
                    return false;
                out = toWorld(part.source.vertices[static_cast<std::size_t>(index)], worldRotation,
                              worldOrigin);
                return true;
            };
            ChVector3d a;
            if (!vert(face.loop[0], a))
                continue;
            for (std::size_t i = 1; i + 1 < face.loop.size(); ++i) {
                ChVector3d b, c;
                if (!vert(face.loop[i], b) || !vert(face.loop[i + 1], c))
                    continue;
                double t = 0.0;
                if (rayHitTriangle(origin, direction, a, b, c, t) && t < best) {
                    best = t;
                    hit = FaceId{partIndex, f};
                    point = origin + direction * t;
                }
            }
        }
    }
    if (hit && hitPoint)
        *hitPoint = point;
    return hit;
}

SplicePreview makeSplicePreview(const PartDefinition& part, std::size_t faceIndex, double angle,
                                std::optional<ChVector3d> through) {
    SplicePreview preview;
    preview.face = {0, faceIndex};
    if (faceIndex >= part.source.faces.size())
        return preview;
    const AuthoredFace& face = part.source.faces[faceIndex];
    if (face.loop.size() < 3)
        return preview;
    preview.loop.reserve(face.loop.size());
    ChVector3d centroid{0, 0, 0};
    for (int index : face.loop) {
        preview.loop.push_back(part.source.vertices[static_cast<std::size_t>(index)]);
        centroid += part.source.vertices[static_cast<std::size_t>(index)];
    }
    centroid /= double(face.loop.size());
    const ChVector3d n = faceNormal(part.source, face);
    ChVector3d origin = centroid;
    if (through) {
        origin = *through - n * ((*through - centroid).Dot(n));
    }
    ChVector3d e1 = preview.loop[1] - preview.loop[0];
    e1 = e1 - n * e1.Dot(n);
    if (e1.Length() < 1.0e-12)
        e1 = n.Cross(ChVector3d(0, 1, 0));
    if (e1.Length() < 1.0e-12)
        e1 = {1, 0, 0};
    e1.Normalize();
    const ChVector3d e2 = n.Cross(e1);
    const ChVector3d dir = e1 * std::cos(angle) + e2 * std::sin(angle);
    std::vector<std::pair<double, SpliceHit>> ranked;
    preview.edgeLength.resize(face.loop.size(), 0.0);
    preview.edgeSplitA.resize(face.loop.size(), -1.0);
    preview.edgeSplitB.resize(face.loop.size(), -1.0);
    for (std::size_t i = 0; i < face.loop.size(); ++i) {
        const ChVector3d a = preview.loop[i];
        const ChVector3d b = preview.loop[(i + 1) % face.loop.size()];
        const ChVector3d edge = b - a;
        preview.edgeLength[i] = edge.Length();
        const double denom = n.Dot(dir.Cross(edge));
        if (std::abs(denom) < 1.0e-12)
            continue;
        const ChVector3d r = a - origin;
        const double u = n.Dot(r.Cross(dir)) / denom;
        if (u < -kHitEps || u > 1.0 + kHitEps)
            continue;
        SpliceHit hit;
        hit.edge = static_cast<int>(i);
        hit.along = std::clamp(u, 0.0, 1.0);
        hit.point = a + edge * hit.along;
        const double s = (hit.point - origin).Dot(dir);
        bool dup = false;
        for (const auto& other : ranked) {
            if ((other.second.point - hit.point).Length() < 1.0e-8)
                dup = true;
        }
        if (!dup)
            ranked.push_back({s, hit});
    }
    SpliceHit neg;
    SpliceHit pos;
    double bestNeg = -1.0e30;
    double bestPos = 1.0e30;
    bool haveNeg = false;
    bool havePos = false;
    for (const auto& entry : ranked) {
        if (entry.first < 0.0 && entry.first > bestNeg) {
            bestNeg = entry.first;
            neg = entry.second;
            haveNeg = true;
        }
        if (entry.first >= 0.0 && entry.first < bestPos) {
            bestPos = entry.first;
            pos = entry.second;
            havePos = true;
        }
    }
    if (!haveNeg || !havePos)
        return preview;
    preview.hits[0] = neg;
    preview.hits[1] = pos;
    preview.chordA = neg.point;
    preview.chordB = pos.point;
    preview.chordLength = (preview.chordB - preview.chordA).Length();
    ChVector3d chordDir = preview.chordB - preview.chordA;
    if (chordDir.Length() > 1.0e-12)
        chordDir.Normalize();
    for (int h = 0; h < 2; ++h) {
        const int e = preview.hits[h].edge;
        const ChVector3d a = preview.loop[static_cast<std::size_t>(e)];
        const ChVector3d b = preview.loop[static_cast<std::size_t>((e + 1) % preview.loop.size())];
        ChVector3d edge = b - a;
        const double len = edge.Length();
        preview.edgeSplitA[static_cast<std::size_t>(e)] = len * preview.hits[h].along;
        preview.edgeSplitB[static_cast<std::size_t>(e)] = len * (1.0 - preview.hits[h].along);
        if (edge.Length() > 1.0e-12)
            edge.Normalize();
        const double cosA = std::clamp(std::abs(chordDir.Dot(edge)), 0.0, 1.0);
        const double deg = std::acos(cosA) * 180.0 / kPi;
        preview.interiorDeg[h] = deg;
        preview.square[h] = std::abs(deg - 90.0) < 0.01;
    }
    preview.valid = preview.chordLength > 1.0e-9;
    return preview;
}

bool spliceFace(PartDefinition& part, const SplicePreview& preview, std::size_t& faceA, std::size_t& faceB) {
    ensureAuthoredFaces(part);
    if (!preview.valid || preview.face.face >= part.source.faces.size())
        return false;
    if (part.source.faces[preview.face.face].loop.size() < 3)
        return false;
    const int v0 = addVertex(part.source, preview.hits[0].point);
    const int v1 = addVertex(part.source, preview.hits[1].point);
    const std::vector<int> originalLoop = part.source.faces[preview.face.face].loop;
    const int edgeCount = static_cast<int>(originalLoop.size());
    auto insertOnSharedEdge = [&](int a, int b, int mid) {
        if (a < 0 || b < 0 || mid < 0 || a == b)
            return;
        for (std::size_t f = 0; f < part.source.faces.size(); ++f) {
            if (f == preview.face.face)
                continue;
            auto& other = part.source.faces[f];
            bool already = false;
            for (int index : other.loop)
                already = already || index == mid;
            if (already)
                continue;
            for (std::size_t i = 0; i < other.loop.size(); ++i) {
                const int x = other.loop[i];
                const int y = other.loop[(i + 1) % other.loop.size()];
                if ((x == a && y == b) || (x == b && y == a)) {
                    if (x != mid && y != mid)
                        other.loop.insert(other.loop.begin() + static_cast<int>(i) + 1, mid);
                    break;
                }
            }
        }
    };
    for (int h = 0; h < 2; ++h) {
        const int e = preview.hits[h].edge;
        if (e < 0 || e >= edgeCount)
            continue;
        if (preview.hits[h].along <= 1.0e-6 || preview.hits[h].along >= 1.0 - 1.0e-6)
            continue;
        insertOnSharedEdge(originalLoop[static_cast<std::size_t>(e)],
                           originalLoop[static_cast<std::size_t>((e + 1) % edgeCount)], h == 0 ? v0 : v1);
    }
    std::vector<int> loop;
    for (int i = 0; i < edgeCount; ++i) {
        loop.push_back(originalLoop[static_cast<std::size_t>(i)]);
        for (int h = 0; h < 2; ++h) {
            if (preview.hits[h].edge != i)
                continue;
            const int vertex = h == 0 ? v0 : v1;
            if (preview.hits[h].along > 1.0e-6 && preview.hits[h].along < 1.0 - 1.0e-6)
                loop.push_back(vertex);
            else if (preview.hits[h].along <= 1.0e-6)
                loop.back() = vertex;
            else
                loop.push_back(vertex);
        }
    }
    auto locate = [&](int vertex) {
        for (std::size_t i = 0; i < loop.size(); ++i)
            if (loop[i] == vertex)
                return static_cast<int>(i);
        return -1;
    };
    const int i0 = locate(v0);
    const int i1 = locate(v1);
    if (i0 < 0 || i1 < 0 || i0 == i1)
        return false;
    auto walk = [&](int from, int to) {
        AuthoredFace out;
        int i = from;
        const int count = static_cast<int>(loop.size());
        for (int guard = 0; guard < count + 2; ++guard) {
            out.loop.push_back(loop[static_cast<std::size_t>(i)]);
            if (i == to)
                break;
            i = (i + 1) % count;
        }
        return out;
    };
    AuthoredFace a = walk(i0, i1);
    AuthoredFace b = walk(i1, i0);
    auto clean = [](AuthoredFace face) {
        std::vector<int> loop;
        for (int index : face.loop) {
            if (loop.empty() || loop.back() != index)
                loop.push_back(index);
        }
        if (loop.size() > 1 && loop.front() == loop.back())
            loop.pop_back();
        face.loop = std::move(loop);
        return face;
    };
    a = clean(std::move(a));
    b = clean(std::move(b));
    if (a.loop.size() < 3 || b.loop.size() < 3)
        return false;
    part.source.faces[preview.face.face] = std::move(a);
    part.source.faces.push_back(std::move(b));
    faceA = preview.face.face;
    faceB = part.source.faces.size() - 1;
    return true;
}

namespace {

struct PlaneBasis {
    ChVector3d origin{0, 0, 0};
    ChVector3d e1{1, 0, 0};
    ChVector3d e2{0, 1, 0};
    ChVector3d n{0, 0, 1};
};

bool authoredLoop(const PartDefinition& part, std::size_t faceIndex, std::vector<ChVector3d>& out) {
    out.clear();
    if (faceIndex >= part.source.faces.size())
        return false;
    for (int index : part.source.faces[faceIndex].loop) {
        if (index < 0 || static_cast<std::size_t>(index) >= part.source.vertices.size())
            continue;
        out.push_back(part.source.vertices[static_cast<std::size_t>(index)]);
    }
    return out.size() >= 3;
}

ChVector3d loopCentroid(const std::vector<ChVector3d>& loop) {
    ChVector3d sum{0, 0, 0};
    if (loop.empty())
        return sum;
    for (const auto& p : loop)
        sum += p;
    return sum / double(loop.size());
}

ChVector3d loopNormal(const std::vector<ChVector3d>& loop) {
    if (loop.size() < 3)
        return {0, 1, 0};
    ChVector3d n{0, 0, 0};
    for (std::size_t i = 0; i < loop.size(); ++i) {
        const ChVector3d& a = loop[i];
        const ChVector3d& b = loop[(i + 1) % loop.size()];
        n += a.Cross(b);
    }
    if (n.Length() < 1.0e-12)
        return {0, 1, 0};
    n.Normalize();
    return n;
}

PlaneBasis makeBasis(const std::vector<ChVector3d>& loop) {
    PlaneBasis b;
    b.origin = loop.empty() ? ChVector3d{0, 0, 0} : loop.front();
    b.n = loopNormal(loop);
    ChVector3d e1 = loop.size() > 1 ? loop[1] - loop[0] : ChVector3d{1, 0, 0};
    e1 = e1 - b.n * e1.Dot(b.n);
    if (e1.Length() < 1.0e-12)
        e1 = b.n.Cross(ChVector3d(0, 1, 0));
    if (e1.Length() < 1.0e-12)
        e1 = {1, 0, 0};
    e1.Normalize();
    b.e1 = e1;
    b.e2 = b.n.Cross(e1);
    return b;
}

ChVector3d toPlane(const PlaneBasis& b, const ChVector3d& p) {
    const ChVector3d d = p - b.origin;
    return {d.Dot(b.e1), d.Dot(b.e2), 0.0};
}

ChVector3d fromPlane(const PlaneBasis& b, const ChVector3d& q) {
    return b.origin + b.e1 * q.x() + b.e2 * q.y();
}

double cross2(const ChVector3d& a, const ChVector3d& b) {
    return a.x() * b.y() - a.y() * b.x();
}

bool leftOf(const ChVector3d& p, const ChVector3d& a, const ChVector3d& b) {
    return cross2(b - a, p - a) >= -1.0e-12;
}

ChVector3d intersect2(const ChVector3d& s, const ChVector3d& e, const ChVector3d& a, const ChVector3d& b) {
    const ChVector3d d = e - s;
    const ChVector3d c = b - a;
    const double den = cross2(d, c);
    if (std::abs(den) < 1.0e-18)
        return s;
    const double t = cross2(a - s, c) / den;
    return s + d * t;
}

void ensureCcw(std::vector<ChVector3d>& poly) {
    if (poly.size() < 3)
        return;
    double area = 0.0;
    for (std::size_t i = 0; i < poly.size(); ++i)
        area += cross2(poly[i], poly[(i + 1) % poly.size()]);
    if (area < 0.0)
        std::reverse(poly.begin(), poly.end());
}

std::vector<ChVector3d> clipConvex(std::vector<ChVector3d> subject, std::vector<ChVector3d> clip) {
    if (subject.size() < 3 || clip.size() < 3)
        return {};
    ensureCcw(clip);
    for (std::size_t i = 0; i < clip.size(); ++i) {
        const ChVector3d a = clip[i];
        const ChVector3d b = clip[(i + 1) % clip.size()];
        std::vector<ChVector3d> input = std::move(subject);
        subject.clear();
        if (input.empty())
            break;
        for (std::size_t j = 0; j < input.size(); ++j) {
            const ChVector3d s = input[j];
            const ChVector3d e = input[(j + 1) % input.size()];
            const bool sIn = leftOf(s, a, b);
            const bool eIn = leftOf(e, a, b);
            if (eIn) {
                if (!sIn)
                    subject.push_back(intersect2(s, e, a, b));
                subject.push_back(e);
            } else if (sIn) {
                subject.push_back(intersect2(s, e, a, b));
            }
        }
    }
    return subject;
}

bool segmentHit2(const ChVector3d& a, const ChVector3d& b, const ChVector3d& c, const ChVector3d& d,
                 ChVector3d& at) {
    const ChVector3d r = b - a;
    const ChVector3d s = d - c;
    const double den = cross2(r, s);
    if (std::abs(den) < 1.0e-18)
        return false;
    const double t = cross2(c - a, s) / den;
    const double u = cross2(c - a, r) / den;
    if (t < -kHitEps || t > 1.0 + kHitEps || u < -kHitEps || u > 1.0 + kHitEps)
        return false;
    at = a + r * std::clamp(t, 0.0, 1.0);
    return true;
}

double interiorDeg(const ChVector3d& dirA, const ChVector3d& dirB) {
    ChVector3d a = dirA;
    ChVector3d b = dirB;
    if (a.Length() < 1.0e-12 || b.Length() < 1.0e-12)
        return 90.0;
    a.Normalize();
    b.Normalize();
    return std::acos(std::clamp(std::abs(a.Dot(b)), 0.0, 1.0)) * 180.0 / kPi;
}

ChQuaterniond kissRotation(const ChVector3d& sourceNormal, const ChVector3d& destNormal, double spin) {
    ChVector3d ns = sourceNormal;
    ChVector3d nd = destNormal;
    if (ns.Length() < 1.0e-12)
        ns = {0, 1, 0};
    else
        ns.Normalize();
    if (nd.Length() < 1.0e-12)
        nd = {0, 1, 0};
    else
        nd.Normalize();
    const ChVector3d target = -nd;
    ChQuaterniond align{1, 0, 0, 0};
    const double d = std::clamp(ns.Dot(target), -1.0, 1.0);
    if (d > 0.999999)
        align = ChQuaterniond(1, 0, 0, 0);
    else if (d < -0.999999) {
        ChVector3d axis = ns.Cross(ChVector3d(1, 0, 0));
        if (axis.Length() < 1.0e-8)
            axis = ns.Cross(ChVector3d(0, 1, 0));
        axis.Normalize();
        align = chrono::QuatFromAngleAxis(kPi, axis);
    } else {
        ChVector3d axis = ns.Cross(target);
        axis.Normalize();
        align = chrono::QuatFromAngleAxis(std::acos(d), axis);
    }
    return chrono::QuatFromAngleAxis(spin, nd) * align;
}

struct Key3 {
    long long x = 0;
    long long y = 0;
    long long z = 0;
    bool operator==(const Key3& o) const { return x == o.x && y == o.y && z == o.z; }
    bool operator!=(const Key3& o) const { return !(*this == o); }
    bool operator<(const Key3& o) const {
        if (x != o.x)
            return x < o.x;
        if (y != o.y)
            return y < o.y;
        return z < o.z;
    }
};

Key3 quantize(const ChVector3d& p) {
    return {std::llround(p.x() * 1.0e7), std::llround(p.y() * 1.0e7), std::llround(p.z() * 1.0e7)};
}

std::vector<ChVector3d> perimeterFromPolys(const std::vector<std::vector<ChVector3d>>& polys) {
    std::map<std::pair<Key3, Key3>, std::pair<ChVector3d, ChVector3d>> unique;
    std::map<std::pair<Key3, Key3>, int> count;
    auto undirected = [](Key3 a, Key3 b) {
        return a < b ? std::make_pair(a, b) : std::make_pair(b, a);
    };
    for (const auto& poly : polys) {
        if (poly.size() < 3)
            continue;
        for (std::size_t i = 0; i < poly.size(); ++i) {
            const ChVector3d a = poly[i];
            const ChVector3d b = poly[(i + 1) % poly.size()];
            if ((a - b).Length() < 1.0e-9)
                continue;
            const auto key = undirected(quantize(a), quantize(b));
            ++count[key];
            unique[key] = {a, b};
        }
    }
    std::map<Key3, std::vector<std::pair<Key3, ChVector3d>>> adj;
    std::map<Key3, ChVector3d> point;
    for (const auto& edge : count) {
        if (edge.second != 1)
            continue;
        const auto found = unique.find(edge.first);
        if (found == unique.end())
            continue;
        const Key3 ka = quantize(found->second.first);
        const Key3 kb = quantize(found->second.second);
        point[ka] = found->second.first;
        point[kb] = found->second.second;
        adj[ka].push_back({kb, found->second.second});
        adj[kb].push_back({ka, found->second.first});
    }
    std::vector<ChVector3d> loop;
    if (adj.empty())
        return loop;
    Key3 start = adj.begin()->first;
    Key3 cursor = start;
    Key3 prev{std::numeric_limits<long long>::min(), 0, 0};
    for (int guard = 0; guard < 4096; ++guard) {
        loop.push_back(point[cursor]);
        const auto& nbrs = adj[cursor];
        Key3 next = cursor;
        bool found = false;
        for (const auto& n : nbrs) {
            if (n.first == prev && nbrs.size() > 1)
                continue;
            next = n.first;
            found = true;
            if (n.first != prev)
                break;
        }
        if (!found)
            break;
        prev = cursor;
        cursor = next;
        if (cursor == start)
            break;
    }
    return loop.size() >= 3 ? loop : std::vector<ChVector3d>{};
}

bool raySegment(const ChVector3d& origin, const ChVector3d& dir, const ChVector3d& a,
                const ChVector3d& b, double& t, double& dist) {
    const ChVector3d ab = b - a;
    const ChVector3d n = dir.Cross(ab);
    const double denom = n.Length2();
    if (denom < 1.0e-18)
        return false;
    const ChVector3d ao = origin - a;
    t = ao.Cross(ab).Dot(n) / denom;
    if (t < 1.0e-6)
        return false;
    const double u = ao.Cross(dir).Dot(n) / denom;
    if (u < -1.0e-6 || u > 1.0 + 1.0e-6)
        return false;
    const ChVector3d hit = origin + dir * t;
    const ChVector3d closest = a + ab * std::clamp(u, 0.0, 1.0);
    dist = (hit - closest).Length();
    return true;
}

void fillChord(WeldChord& chord, const std::vector<ChVector3d>& dest, const ChVector3d& mouse,
               const ChVector3d& n, const ChVector3d& dir) {
    chord = {};
    std::vector<std::pair<double, SpliceHit>> ranked;
    for (std::size_t i = 0; i < dest.size(); ++i) {
        const ChVector3d a = dest[i];
        const ChVector3d b = dest[(i + 1) % dest.size()];
        const ChVector3d edge = b - a;
        const double denom = n.Dot(dir.Cross(edge));
        if (std::abs(denom) < 1.0e-12)
            continue;
        const ChVector3d r = a - mouse;
        const double u = n.Dot(r.Cross(dir)) / denom;
        if (u < -kHitEps || u > 1.0 + kHitEps)
            continue;
        SpliceHit hit;
        hit.edge = static_cast<int>(i);
        hit.along = std::clamp(u, 0.0, 1.0);
        hit.point = a + edge * hit.along;
        ranked.push_back({(hit.point - mouse).Dot(dir), hit});
    }
    SpliceHit neg;
    SpliceHit pos;
    double bestNeg = -1.0e30;
    double bestPos = 1.0e30;
    bool haveNeg = false;
    bool havePos = false;
    for (const auto& entry : ranked) {
        if (entry.first < 0.0 && entry.first > bestNeg) {
            bestNeg = entry.first;
            neg = entry.second;
            haveNeg = true;
        }
        if (entry.first >= 0.0 && entry.first < bestPos) {
            bestPos = entry.first;
            pos = entry.second;
            havePos = true;
        }
    }
    if (!haveNeg || !havePos)
        return;
    chord.hits[0] = neg;
    chord.hits[1] = pos;
    chord.a = neg.point;
    chord.b = pos.point;
    chord.length = (chord.b - chord.a).Length();
    ChVector3d chordDir = chord.b - chord.a;
    if (chordDir.Length() > 1.0e-12)
        chordDir.Normalize();
    for (int h = 0; h < 2; ++h) {
        const int e = chord.hits[h].edge;
        if (e < 0)
            continue;
        ChVector3d edge = dest[static_cast<std::size_t>((e + 1) % dest.size())] -
                          dest[static_cast<std::size_t>(e)];
        const double deg = interiorDeg(chordDir, edge);
        chord.interiorDeg[h] = deg;
        chord.square[h] = std::abs(deg - 90.0) < 0.01;
    }
}

void coincidentEdges(const PartDefinition& a, const PartDefinition& b, double tolerance,
                     std::vector<std::vector<ChVector3d>>& out) {
    auto nearOther = [&](const ChVector3d& p) {
        for (const auto& q : b.source.vertices) {
            if ((p - q).Length() <= tolerance)
                return true;
        }
        return false;
    };
    for (const auto& face : a.source.faces) {
        if (face.loop.size() < 2)
            continue;
        for (std::size_t i = 0; i < face.loop.size(); ++i) {
            const int ia = face.loop[i];
            const int ib = face.loop[(i + 1) % face.loop.size()];
            if (ia < 0 || ib < 0)
                continue;
            if (static_cast<std::size_t>(ia) >= a.source.vertices.size() ||
                static_cast<std::size_t>(ib) >= a.source.vertices.size())
                continue;
            const ChVector3d pa = a.source.vertices[static_cast<std::size_t>(ia)];
            const ChVector3d pb = a.source.vertices[static_cast<std::size_t>(ib)];
            if (nearOther(pa) && nearOther(pb))
                out.push_back({pa, pb});
        }
    }
}

}  // namespace

WeldPreview makeWeldPreview(const AssemblyDefinition& definition,
                            const std::vector<FaceId>& destFaces, FaceId source, FaceId hovered,
                            const ChVector3d& mouse, double outlineAngle, double chordAngle) {
    WeldPreview preview;
    preview.source = source;
    preview.dest = hovered;
    if (hovered.part >= definition.parts.size() || source.part >= definition.parts.size())
        return preview;
    if (hovered.part == source.part)
        return preview;
    const PartDefinition& destPart = definition.parts[hovered.part];
    const PartDefinition& sourcePart = definition.parts[source.part];
    std::vector<ChVector3d> destLoop;
    std::vector<ChVector3d> sourceLoop;
    if (!authoredLoop(destPart, hovered.face, destLoop) ||
        !authoredLoop(sourcePart, source.face, sourceLoop))
        return preview;

    const ChVector3d destN = loopNormal(destLoop);
    const ChVector3d sourceN = loopNormal(sourceLoop);
    preview.sourceCentroid = loopCentroid(sourceLoop);
    preview.mouse = mouse - destN * (mouse - loopCentroid(destLoop)).Dot(destN);
    preview.rotation = kissRotation(sourceN, destN, outlineAngle);

    preview.outline.reserve(sourceLoop.size());
    for (const auto& p : sourceLoop)
        preview.outline.push_back(preview.mouse + preview.rotation.Rotate(p - preview.sourceCentroid));

    const PlaneBasis basis = makeBasis(destLoop);
    std::vector<FaceId> canvas;
    for (const FaceId& id : destFaces) {
        if (id.part >= definition.parts.size())
            continue;
        std::vector<ChVector3d> loop;
        if (!authoredLoop(definition.parts[id.part], id.face, loop))
            continue;
        const ChVector3d n = loopNormal(loop);
        if (std::abs(n.Dot(destN)) < kNormalDot)
            continue;
        if (std::abs((loop.front() - destLoop.front()).Dot(destN)) > kCoplanar)
            continue;
        canvas.push_back(id);
        WeldDestLoop drawn;
        drawn.face = id;
        drawn.loop = loop;
        drawn.edgeLength.resize(loop.size(), 0.0);
        for (std::size_t i = 0; i < loop.size(); ++i)
            drawn.edgeLength[i] = (loop[(i + 1) % loop.size()] - loop[i]).Length();
        preview.destLoops.push_back(std::move(drawn));
    }
    if (preview.destLoops.empty()) {
        WeldDestLoop drawn;
        drawn.face = hovered;
        drawn.loop = destLoop;
        drawn.edgeLength.resize(destLoop.size(), 0.0);
        for (std::size_t i = 0; i < destLoop.size(); ++i)
            drawn.edgeLength[i] = (destLoop[(i + 1) % destLoop.size()] - destLoop[i]).Length();
        preview.destLoops.push_back(std::move(drawn));
        canvas.push_back(hovered);
    }

    const bool span = canvas.size() > 1;
    std::vector<ChVector3d> drawnOutline = preview.outline;
    if (!span) {
        std::vector<ChVector3d> dest2;
        std::vector<ChVector3d> src2;
        for (const auto& p : destLoop)
            dest2.push_back(toPlane(basis, p));
        for (const auto& p : preview.outline)
            src2.push_back(toPlane(basis, p));
        const auto clipped = clipConvex(src2, dest2);
        if (clipped.size() >= 3) {
            drawnOutline.clear();
            for (const auto& q : clipped)
                drawnOutline.push_back(fromPlane(basis, q));
            preview.outline = drawnOutline;
        }
    }

    std::vector<std::vector<ChVector3d>> intersections;
    for (const FaceId& id : canvas) {
        std::vector<ChVector3d> loop;
        if (!authoredLoop(definition.parts[id.part], id.face, loop))
            continue;
        std::vector<ChVector3d> dest2;
        std::vector<ChVector3d> src2;
        const PlaneBasis local = makeBasis(loop);
        for (const auto& p : loop)
            dest2.push_back(toPlane(local, p));
        for (const auto& p : preview.outline)
            src2.push_back(toPlane(local, p));
        auto clipped = clipConvex(src2, dest2);
        if (clipped.size() < 3)
            continue;
        std::vector<ChVector3d> world;
        for (const auto& q : clipped)
            world.push_back(fromPlane(local, q));
        intersections.push_back(std::move(world));
    }
    preview.bead = perimeterFromPolys(intersections);
    if (preview.bead.size() < 3 && !intersections.empty())
        preview.bead = intersections.front();

    std::vector<ChVector3d> dest2;
    std::vector<ChVector3d> src2;
    for (const auto& p : destLoop)
        dest2.push_back(toPlane(basis, p));
    for (const auto& p : preview.outline)
        src2.push_back(toPlane(basis, p));
    for (std::size_t i = 0; i < src2.size(); ++i) {
        const ChVector3d a = src2[i];
        const ChVector3d b = src2[(i + 1) % src2.size()];
        for (const FaceId& id : canvas) {
            std::vector<ChVector3d> loop;
            if (!authoredLoop(definition.parts[id.part], id.face, loop))
                continue;
            std::vector<ChVector3d> other2;
            for (const auto& p : loop)
                other2.push_back(toPlane(basis, p));
            for (std::size_t j = 0; j < other2.size(); ++j) {
                ChVector3d at;
                if (!segmentHit2(a, b, other2[j], other2[(j + 1) % other2.size()], at))
                    continue;
                WeldEdgeHit hit;
                hit.point = fromPlane(basis, at);
                const ChVector3d outlineDir = preview.outline[(i + 1) % preview.outline.size()] - preview.outline[i];
                const ChVector3d destDir = loop[(j + 1) % loop.size()] - loop[j];
                hit.interiorDeg = interiorDeg(outlineDir, destDir);
                hit.square = std::abs(hit.interiorDeg - 90.0) < 0.01;
                bool dup = false;
                for (const auto& other : preview.outlineHits) {
                    if ((other.point - hit.point).Length() < 1.0e-6)
                        dup = true;
                }
                if (!dup)
                    preview.outlineHits.push_back(hit);
            }
        }
    }

    ChVector3d e1 = destLoop.size() > 1 ? destLoop[1] - destLoop[0] : ChVector3d{1, 0, 0};
    e1 = e1 - destN * e1.Dot(destN);
    if (e1.Length() < 1.0e-12)
        e1 = destN.Cross(ChVector3d(0, 1, 0));
    if (e1.Length() < 1.0e-12)
        e1 = {1, 0, 0};
    e1.Normalize();
    const ChVector3d e2 = destN.Cross(e1);
    const ChVector3d dirA = e1 * std::cos(chordAngle) + e2 * std::sin(chordAngle);
    const ChVector3d dirB = destN.Cross(dirA);
    fillChord(preview.chords[0], destLoop, preview.mouse, destN, dirA);
    fillChord(preview.chords[1], destLoop, preview.mouse, destN, dirB);

    preview.valid = preview.outline.size() >= 3;
    return preview;
}

bool commitWeldPlacement(AssemblyDefinition& definition, const WeldPreview& preview) {
    if (!preview.valid || preview.source.part >= definition.parts.size() ||
        preview.dest.part >= definition.parts.size() || preview.source.part == preview.dest.part)
        return false;
    PartDefinition& source = definition.parts[preview.source.part];
    PartDefinition& dest = definition.parts[preview.dest.part];
    applyPartRotation(source, preview.rotation, preview.sourceCentroid);
    applyPartTranslation(source, preview.mouse - preview.sourceCentroid);

    WeldDefinition weld;
    weld.partA = preview.dest.part;
    weld.partB = preview.source.part;
    weld.label = dest.name + " - " + source.name;
    weld.throat = std::min(std::max(dest.thickness, 1.0e-4), std::max(source.thickness, 1.0e-4));
    weld.tolerance = std::max(2.0e-3, 0.25 * weld.throat);
    weld.bead = preview.bead;
    if (weld.bead.size() < 3)
        weld.bead = preview.outline;
    definition.welds.push_back(std::move(weld));
    return true;
}

void authoredWeldPolylines(const AssemblyDefinition& definition, std::size_t weldIndex,
                           std::vector<std::vector<ChVector3d>>& out) {
    out.clear();
    if (weldIndex >= definition.welds.size())
        return;
    const WeldDefinition& weld = definition.welds[weldIndex];
    if (weld.bead.size() >= 2) {
        out.push_back(weld.bead);
        return;
    }
    if (weld.partA >= definition.parts.size() || weld.partB >= definition.parts.size())
        return;
    coincidentEdges(definition.parts[weld.partA], definition.parts[weld.partB], weld.tolerance, out);
    coincidentEdges(definition.parts[weld.partB], definition.parts[weld.partA], weld.tolerance, out);
}

std::optional<std::size_t> pickWeld(const AssemblyDefinition& definition,
                                    const std::vector<std::size_t>& parts,
                                    const ChVector3d& origin, const ChVector3d& direction,
                                    const ChQuaterniond& worldRotation, const ChVector3d& worldOrigin,
                                    double maxDistance) {
    std::optional<std::size_t> hit;
    double bestT = 1.0e9;
    auto involves = [&](const WeldDefinition& weld) {
        if (parts.empty())
            return true;
        return std::find(parts.begin(), parts.end(), weld.partA) != parts.end() ||
               std::find(parts.begin(), parts.end(), weld.partB) != parts.end();
    };
    for (std::size_t i = 0; i < definition.welds.size(); ++i) {
        if (!involves(definition.welds[i]))
            continue;
        std::vector<std::vector<ChVector3d>> loops;
        authoredWeldPolylines(definition, i, loops);
        for (const auto& loop : loops) {
            if (loop.size() < 2)
                continue;
            const std::size_t segs = loop.size() >= 3 ? loop.size() : loop.size() - 1;
            for (std::size_t s = 0; s < segs; ++s) {
                const ChVector3d a = toWorld(loop[s], worldRotation, worldOrigin);
                const ChVector3d b = toWorld(loop[(s + 1) % loop.size()], worldRotation, worldOrigin);
                double t = 0.0;
                double dist = 0.0;
                if (!raySegment(origin, direction, a, b, t, dist))
                    continue;
                if (dist > maxDistance || t >= bestT)
                    continue;
                bestT = t;
                hit = i;
            }
        }
    }
    return hit;
}

std::string formatMillimetres(double metres) {
    char buffer[64];
    std::snprintf(buffer, sizeof(buffer), "%.2f mm", metres * 1000.0);
    return buffer;
}

std::string formatDegrees(double degrees) {
    char buffer[64];
    std::snprintf(buffer, sizeof(buffer), "%.2f deg", degrees);
    return buffer;
}

bool partIsDescendant(const AssemblyDefinition& definition, std::size_t candidate, std::size_t ancestor) {
    int cursor = definition.parts[candidate].parentPart;
    std::set<int> seen;
    while (cursor >= 0 && static_cast<std::size_t>(cursor) < definition.parts.size()) {
        if (static_cast<std::size_t>(cursor) == ancestor)
            return true;
        if (!seen.insert(cursor).second)
            return true;
        cursor = definition.parts[static_cast<std::size_t>(cursor)].parentPart;
    }
    return false;
}

bool groupIsDescendant(const AssemblyDefinition& definition, int candidate, int ancestor) {
    int cursor = candidate;
    std::set<int> seen;
    while (cursor >= 0 && static_cast<std::size_t>(cursor) < definition.groups.size()) {
        if (cursor == ancestor)
            return true;
        if (!seen.insert(cursor).second)
            return true;
        cursor = definition.groups[static_cast<std::size_t>(cursor)].parent;
    }
    return false;
}

void renormalizeSiblingOrder(AssemblyDefinition& definition) {
    auto sortParts = [&](int parent, int group) {
        std::vector<std::size_t> kids;
        for (std::size_t i = 0; i < definition.parts.size(); ++i) {
            const auto& p = definition.parts[i];
            if (p.parentPart == parent && p.group == group)
                kids.push_back(i);
        }
        std::sort(kids.begin(), kids.end(), [&](std::size_t a, std::size_t b) {
            if (definition.parts[a].order != definition.parts[b].order)
                return definition.parts[a].order < definition.parts[b].order;
            return a < b;
        });
        for (std::size_t i = 0; i < kids.size(); ++i)
            definition.parts[kids[i]].order = static_cast<int>(i);
    };
    for (std::size_t i = 0; i < definition.parts.size(); ++i)
        sortParts(definition.parts[i].parentPart, definition.parts[i].group);
    sortParts(-1, -1);
    for (int g = 0; g < static_cast<int>(definition.groups.size()); ++g)
        sortParts(-1, g);

    std::vector<int> folders;
    for (int g = 0; g < static_cast<int>(definition.groups.size()); ++g)
        folders.push_back(g);
    std::sort(folders.begin(), folders.end(), [&](int a, int b) {
        if (definition.groups[static_cast<std::size_t>(a)].parent !=
            definition.groups[static_cast<std::size_t>(b)].parent)
            return definition.groups[static_cast<std::size_t>(a)].parent <
                   definition.groups[static_cast<std::size_t>(b)].parent;
        if (definition.groups[static_cast<std::size_t>(a)].order !=
            definition.groups[static_cast<std::size_t>(b)].order)
            return definition.groups[static_cast<std::size_t>(a)].order <
                   definition.groups[static_cast<std::size_t>(b)].order;
        return a < b;
    });
    int lastParent = std::numeric_limits<int>::min();
    int rank = 0;
    for (int g : folders) {
        const int parent = definition.groups[static_cast<std::size_t>(g)].parent;
        if (parent != lastParent) {
            lastParent = parent;
            rank = 0;
        }
        definition.groups[static_cast<std::size_t>(g)].order = rank++;
    }
}

std::vector<std::size_t> collectTransformTargets(const AssemblyDefinition& definition,
                                                 const std::vector<std::size_t>& selectedParts,
                                                 const std::optional<std::size_t>& selectedGroup) {
    std::vector<std::size_t> group;
    const auto& parts = definition.parts;
    const auto& folders = definition.groups;
    std::vector<char> folder(folders.size(), 0);
    if (selectedGroup && *selectedGroup < folders.size()) {
        std::vector<int> stack{static_cast<int>(*selectedGroup)};
        while (!stack.empty()) {
            const int current = stack.back();
            stack.pop_back();
            if (current < 0 || static_cast<std::size_t>(current) >= folder.size() ||
                folder[static_cast<std::size_t>(current)])
                continue;
            folder[static_cast<std::size_t>(current)] = 1;
            for (int child = 0; child < static_cast<int>(folders.size()); ++child)
                if (folders[static_cast<std::size_t>(child)].parent == current)
                    stack.push_back(child);
        }
        for (std::size_t i = 0; i < parts.size(); ++i)
            if (parts[i].group >= 0 && static_cast<std::size_t>(parts[i].group) < folder.size() &&
                folder[static_cast<std::size_t>(parts[i].group)] && parts[i].parentPart < 0)
                group.push_back(i);
    } else {
        group = selectedParts;
    }
    for (std::size_t n = 0; n < group.size(); ++n) {
        for (std::size_t i = 0; i < parts.size(); ++i) {
            if (parts[i].parentPart == static_cast<int>(group[n]) &&
                std::find(group.begin(), group.end(), i) == group.end())
                group.push_back(i);
        }
    }

    std::vector<char> in(parts.size(), 0);
    for (std::size_t index : group)
        if (index < in.size())
            in[index] = 1;
    auto take = [&](std::size_t index) {
        if (index >= in.size() || in[index])
            return false;
        in[index] = 1;
        group.push_back(index);
        return true;
    };
    bool grew = true;
    while (grew) {
        grew = false;
        for (const auto& weld : definition.welds) {
            if (weld.partA >= in.size() || weld.partB >= in.size())
                continue;
            if (in[weld.partA] == in[weld.partB])
                continue;
            grew = take(weld.partA) || grew;
            grew = take(weld.partB) || grew;
        }
        for (const auto& screw : definition.screws) {
            if (screw.partA >= in.size() || screw.partB >= in.size())
                continue;
            if (in[screw.partA] == in[screw.partB])
                continue;
            grew = take(screw.partA) || grew;
            grew = take(screw.partB) || grew;
        }
        for (std::size_t n = 0; n < group.size(); ++n) {
            for (std::size_t i = 0; i < parts.size(); ++i) {
                if (parts[i].parentPart == static_cast<int>(group[n]))
                    grew = take(i) || grew;
            }
        }
    }
    return group;
}

void transformAuthoredJoints(AssemblyDefinition& definition,
                             const std::vector<char>& inGroup,
                             const chrono::ChVector3d& translation) {
    for (auto& screw : definition.screws) {
        if (screw.partA < inGroup.size() && inGroup[screw.partA])
            screw.position += translation;
        else if (screw.partB < inGroup.size() && inGroup[screw.partB])
            screw.position += translation;
    }
    for (auto& weld : definition.welds) {
        const bool incident = (weld.partA < inGroup.size() && inGroup[weld.partA]) ||
                              (weld.partB < inGroup.size() && inGroup[weld.partB]);
        if (!incident)
            continue;
        if (weld.useRegion) {
            weld.regionMin += translation;
            weld.regionMax += translation;
        }
        for (auto& point : weld.bead)
            point += translation;
    }
}

void rotateAuthoredJoints(AssemblyDefinition& definition,
                          const std::vector<char>& inGroup,
                          const chrono::ChQuaterniond& rotation,
                          const chrono::ChVector3d& pivot) {
    auto spin = [&](chrono::ChVector3d& p) { p = pivot + rotation.Rotate(p - pivot); };
    for (auto& screw : definition.screws) {
        if (screw.partA < inGroup.size() && inGroup[screw.partA])
            spin(screw.position);
        else if (screw.partB < inGroup.size() && inGroup[screw.partB])
            spin(screw.position);
    }
    for (auto& weld : definition.welds) {
        const bool incident = (weld.partA < inGroup.size() && inGroup[weld.partA]) ||
                              (weld.partB < inGroup.size() && inGroup[weld.partB]);
        if (!incident)
            continue;
        if (weld.useRegion) {
            spin(weld.regionMin);
            spin(weld.regionMax);
        }
        for (auto& point : weld.bead)
            spin(point);
    }
}

void deleteSelectedParts(AssemblyDefinition& definition, std::vector<std::size_t>& selected,
                         std::vector<ImpactLoad>& loads) {
    if (selected.empty())
        return;
    std::vector<char> gone(definition.parts.size(), 0);
    for (std::size_t i : selected) {
        if (i < gone.size())
            gone[i] = 1;
    }
    for (auto& part : definition.parts) {
        if (part.parentPart >= 0 && static_cast<std::size_t>(part.parentPart) < gone.size() &&
            gone[static_cast<std::size_t>(part.parentPart)])
            part.parentPart = -1;
    }
    std::vector<WeldDefinition> welds;
    for (auto& weld : definition.welds) {
        if (weld.partA < gone.size() && weld.partB < gone.size() && !gone[weld.partA] && !gone[weld.partB]) {
            weld.partA = static_cast<std::size_t>(remapPartIndex(static_cast<int>(weld.partA), gone));
            weld.partB = static_cast<std::size_t>(remapPartIndex(static_cast<int>(weld.partB), gone));
            welds.push_back(weld);
        }
    }
    definition.welds.swap(welds);
    std::vector<ScrewDefinition> screws;
    for (auto& screw : definition.screws) {
        if (screw.partA < gone.size() && screw.partB < gone.size() && !gone[screw.partA] && !gone[screw.partB]) {
            screw.partA = static_cast<std::size_t>(remapPartIndex(static_cast<int>(screw.partA), gone));
            screw.partB = static_cast<std::size_t>(remapPartIndex(static_cast<int>(screw.partB), gone));
            screws.push_back(screw);
        }
    }
    definition.screws.swap(screws);
    std::vector<ImpactLoad> keptLoads;
    for (auto& load : loads) {
        if (load.partIndex < gone.size() && gone[load.partIndex])
            continue;
        const int mapped = remapPartIndex(static_cast<int>(load.partIndex), gone);
        if (mapped < 0)
            continue;
        load.partIndex = static_cast<std::size_t>(mapped);
        keptLoads.push_back(load);
    }
    loads.swap(keptLoads);

    std::vector<PartDefinition> parts;
    parts.reserve(definition.parts.size());
    for (std::size_t i = 0; i < definition.parts.size(); ++i) {
        if (gone[i])
            continue;
        auto part = definition.parts[i];
        part.parentPart = remapPartIndex(part.parentPart, gone);
        parts.push_back(std::move(part));
    }
    definition.parts.swap(parts);
    selected.clear();
}

void deleteSelectedGroup(AssemblyDefinition& definition, std::size_t groupIndex,
                         std::vector<std::size_t>& selectedParts) {
    if (groupIndex >= definition.groups.size())
        return;
    std::vector<char> gone(definition.groups.size(), 0);
    gone[groupIndex] = 1;
    const int parent = definition.groups[groupIndex].parent;
    for (auto& group : definition.groups) {
        if (group.parent == static_cast<int>(groupIndex))
            group.parent = parent;
    }
    for (auto& part : definition.parts) {
        if (part.group == static_cast<int>(groupIndex))
            part.group = parent;
    }
    std::vector<GroupDefinition> groups;
    for (std::size_t i = 0; i < definition.groups.size(); ++i) {
        if (gone[i])
            continue;
        auto group = definition.groups[i];
        group.parent = remapGroupIndex(group.parent, gone);
        groups.push_back(std::move(group));
    }
    definition.groups.swap(groups);
    for (auto& part : definition.parts)
        part.group = remapGroupIndex(part.group, gone);
    selectedParts.clear();
}

}  // namespace rcs
