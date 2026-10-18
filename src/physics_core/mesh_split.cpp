#include "physics_core/mesh_split.hpp"
#include "physics_core/assembly.hpp"
#include "physics_core/element_policy.hpp"
#include "physics_core/physics_units.hpp"

#define TETLIBRARY
#include "tetgen.h"

#include <algorithm>
#include <cmath>
#include <cstdio>
#include <map>
#include <set>
#include <unordered_map>
#include <utility>
#include <vector>

namespace rcs {
namespace {

using chrono::ChVector3d;

struct RayHit {
    double distance = 1.0e30;
    int face = -1;
};

bool rayTriangle(const ChVector3d& origin, const ChVector3d& direction,
                 const ChVector3d& v0, const ChVector3d& v1, const ChVector3d& v2,
                 double& t) {
    const ChVector3d e1 = v1 - v0;
    const ChVector3d e2 = v2 - v0;
    const ChVector3d p = direction.Cross(e2);
    const double det = e1.Dot(p);
    if (std::abs(det) < 1.0e-12)
        return false;
    const double inv = 1.0 / det;
    const ChVector3d s = origin - v0;
    const double u = s.Dot(p) * inv;
    if (u < -1.0e-6 || u > 1.0 + 1.0e-6)
        return false;
    const ChVector3d q = s.Cross(e1);
    const double v = direction.Dot(q) * inv;
    if (v < -1.0e-6 || u + v > 1.0 + 1.0e-6)
        return false;
    t = e2.Dot(q) * inv;
    return t > 1.0e-5;
}

ChVector3d faceNormal(const SurfaceMesh& mesh, const std::array<int, 3>& face) {
    const ChVector3d& a = mesh.vertices[static_cast<std::size_t>(face[0])];
    const ChVector3d& b = mesh.vertices[static_cast<std::size_t>(face[1])];
    const ChVector3d& c = mesh.vertices[static_cast<std::size_t>(face[2])];
    ChVector3d n = (b - a).Cross(c - a);
    const double len = n.Length();
    if (len < 1.0e-16)
        return {0, 0, 0};
    return n / len;
}

bool sharesVertex(const std::array<int, 3>& a, const std::array<int, 3>& b) {
    for (int i : a)
        for (int j : b)
            if (i == j)
                return true;
    return false;
}

struct Key {
    long long x = 0;
    long long y = 0;
    long long z = 0;
    bool operator==(const Key& other) const { return x == other.x && y == other.y && z == other.z; }
};

struct KeyHash {
    std::size_t operator()(const Key& key) const {
        return static_cast<std::size_t>(key.x * 73856093 ^ key.y * 19349663 ^ key.z * 83492791);
    }
};

Key quantize(const ChVector3d& p, double tol) {
    return {std::llround(p.x() / tol), std::llround(p.y() / tol), std::llround(p.z() / tol)};
}

int weldVertex(std::vector<ChVector3d>& vertices, std::unordered_map<Key, int, KeyHash>& index,
               const ChVector3d& p, double tol) {
    const Key key = quantize(p, tol);
    const auto found = index.find(key);
    if (found != index.end())
        return found->second;
    const int id = static_cast<int>(vertices.size());
    vertices.push_back(p);
    index.emplace(key, id);
    return id;
}

double clampSize(double value, const MeshSizing& sizing) {
    const double lo = std::max(sizing.minElement, 1.0e-6);
    const double hi = std::max(lo, sizing.maxElement > lo ? sizing.maxElement : lo);
    return std::clamp(value, lo, hi);
}

bool fillSolid(const SurfaceMesh& skin, const std::vector<double>& sizeAtVertex, VolumeMesh& volume,
               MeshSizing& sizing) {
    if (skin.triangles.size() < 4 || skin.vertices.size() < 4)
        return false;

    SurfaceMesh plc;
    std::vector<double> size;
    std::vector<int> marker;
    plc.vertices.reserve(skin.vertices.size());
    size.reserve(skin.vertices.size());
    plc.triangles.reserve(skin.triangles.size());
    marker.reserve(skin.triangles.size());
    std::vector<int> remap(skin.vertices.size(), -1);
    auto take = [&](int index) {
        if (index < 0 || static_cast<std::size_t>(index) >= remap.size())
            return -1;
        if (remap[static_cast<std::size_t>(index)] < 0) {
            remap[static_cast<std::size_t>(index)] = static_cast<int>(plc.vertices.size());
            plc.vertices.push_back(skin.vertices[static_cast<std::size_t>(index)]);
            const double local = static_cast<std::size_t>(index) < sizeAtVertex.size()
                                     ? sizeAtVertex[static_cast<std::size_t>(index)]
                                     : 0.0;
            size.push_back(local);
        }
        return remap[static_cast<std::size_t>(index)];
    };
    for (const auto& tri : skin.triangles) {
        const int a = take(tri[0]);
        const int b = take(tri[1]);
        const int c = take(tri[2]);
        if (a < 0 || b < 0 || c < 0 || a == b || b == c || c == a)
            continue;
        plc.triangles.push_back({a, b, c});
        const std::size_t slot = static_cast<std::size_t>(&tri - skin.triangles.data());
        marker.push_back(slot < skin.markers.size() ? skin.markers[slot] : 1);
    }
    if (plc.triangles.size() < 4 || plc.vertices.size() < 4)
        return false;

    const bool flip = meshVolume(plc) < 0.0;
    double requested = 0.0;
    int sized = 0;
    for (double h : size) {
        if (h > 1.0e-12) {
            requested += h;
            ++sized;
        }
    }
    if (sized > 0)
        requested /= static_cast<double>(sized);
    else
        requested = budgetEdge(std::abs(meshVolume(plc)), sizing);
    // The caller's sizing field is a request, not a floor. Raising only the
    // ceiling cannot bring a mesh back under budget when the field itself is
    // what is refining it, which is the case for anything with a thin wall:
    // the wall asks for elements a third of its thickness no matter how coarse
    // the ceiling gets. Both have to be relaxed together.
    double relax = 1.0;
    const int kBudget = std::max(400, sizing.tetBudget);
    {
        double finest = requested;
        for (double h : size)
            if (h > 1.0e-12)
                finest = std::min(finest, h);
        const double guess = std::abs(meshVolume(plc)) / std::max(finest * finest * finest / 6.0, 1.0e-18);
        if (guess > static_cast<double>(kBudget))
            relax = std::clamp(std::cbrt(guess / static_cast<double>(kBudget)), 1.0, 8.0);
    }
    const bool allowProtect = relax <= 1.0;
    for (int attempt = 0; attempt < 4; ++attempt) {
        tetgenio in;
        tetgenio out;
        const int n = static_cast<int>(plc.vertices.size());
        in.firstnumber = 0;
        in.numberofpoints = n;
        in.pointlist = new REAL[static_cast<std::size_t>(n) * 3];
        in.numberofpointmtrs = 1;
        in.pointmtrlist = new REAL[static_cast<std::size_t>(n)];
        for (int i = 0; i < n; ++i) {
            const ChVector3d& p = plc.vertices[static_cast<std::size_t>(i)];
            in.pointlist[3 * i] = p.x();
            in.pointlist[3 * i + 1] = p.y();
            in.pointlist[3 * i + 2] = p.z();
            const double local = static_cast<std::size_t>(i) < size.size() && size[static_cast<std::size_t>(i)] > 1.0e-12
                                     ? size[static_cast<std::size_t>(i)]
                                     : requested;
            // A fine request is a thin feature. Coarsen the stocky field on the
            // early attempts. The last attempt coarsens everything so a slender
            // authored solid still comes back under budget instead of exploding.
            const bool protect = allowProtect && attempt + 1 < 4 && local < 0.5 * requested;
            in.pointmtrlist[i] = protect ? local : local * relax;
        }
        const int faceCount = static_cast<int>(plc.triangles.size());
        in.numberoffacets = faceCount;
        in.facetlist = new tetgenio::facet[static_cast<std::size_t>(faceCount)];
        in.facetmarkerlist = new int[static_cast<std::size_t>(faceCount)];
        for (int f = 0; f < faceCount; ++f) {
            tetgenio::facet& facet = in.facetlist[f];
            facet.numberofpolygons = 1;
            facet.polygonlist = new tetgenio::polygon[1];
            facet.numberofholes = 0;
            facet.holelist = nullptr;
            facet.polygonlist[0].numberofvertices = 3;
            facet.polygonlist[0].vertexlist = new int[3];
            const auto& tri = plc.triangles[static_cast<std::size_t>(f)];
            facet.polygonlist[0].vertexlist[0] = tri[0];
            facet.polygonlist[0].vertexlist[1] = tri[flip ? 2 : 1];
            facet.polygonlist[0].vertexlist[2] = tri[flip ? 1 : 2];
            in.facetmarkerlist[f] = f < static_cast<int>(marker.size()) ? marker[static_cast<std::size_t>(f)] : 1;
        }

        // Volume cap from the stocky part of the field. The -m metric still
        // refines a fine vertex below this cap, so a web is not coarsened
        // into one flat tet just to pay for the boss.
        double coarse = 0.0;
        int coarseCount = 0;
        for (double h : size) {
            if (h >= 0.5 * requested) {
                coarse += h;
                ++coarseCount;
            }
        }
        const double capBase = coarseCount > 0 ? coarse / static_cast<double>(coarseCount) : requested;
        const double capH = std::max(capBase * relax, 1.0e-6);
        const double maxVolume = std::max(1.0e-12, capH * capH * capH / 6.0);
        char switches[80];
        // The radius-edge bound alone still permits slivers: four nodes can be
        // almost coplanar and pass it. A sliver has a near-singular shape
        // function, which drives the element stiffness up without adding any
        // real stiffness to the part, and the integrator then pays for it with
        // a collapsed step size. Bounding the minimum dihedral angle as well
        // makes tetgen split them out. Sharp features such as a cone apex are
        // where this matters, because that is where slivers are unavoidable
        // without the extra constraint.
        std::snprintf(switches, sizeof(switches), "pq1.4/12mQO9/7a%.8g", maxVolume);
        bool ok = true;
        try {
            tetrahedralize(switches, &in, &out);
        } catch (int code) {
            std::fprintf(stderr, "tetgen failed (%d) on %d faces\n", code, faceCount);
            ok = false;
        } catch (...) {
            std::fprintf(stderr, "tetgen failed on %d faces\n", faceCount);
            ok = false;
        }
        if (!ok || out.numberoftetrahedra <= 0 || out.numberofpoints <= 0 || out.numberofcorners < 4) {
            requested *= 1.6;
            relax *= 1.6;
            continue;
        }
        if (out.numberoftetrahedra > kBudget && attempt + 1 < 4) {
            // Element count goes as the cube of the size, so scale by the
            // cube root of the overshoot and land near the budget in one go
            // instead of creeping up on it.
            const double over = static_cast<double>(out.numberoftetrahedra) / double(kBudget);
            const double scale = std::clamp(std::cbrt(over), 1.2, 3.0);
            requested *= scale;
            relax *= scale;
            continue;
        }

        volume.vertices.clear();
        volume.tets.clear();
        volume.boundary.clear();
        volume.vertices.reserve(static_cast<std::size_t>(out.numberofpoints));
        // TetGen keeps PLC vertices as the first n output points. Write those
        // back from the input so a rounding nudge cannot dent the authored
        // skin. Steiner points (n..) lie on PLC facets or in the interior.
        for (int i = 0; i < out.numberofpoints; ++i) {
            if (i < n)
                volume.vertices.push_back(plc.vertices[static_cast<std::size_t>(i)]);
            else
                volume.vertices.push_back(
                    {out.pointlist[3 * i], out.pointlist[3 * i + 1], out.pointlist[3 * i + 2]});
        }
        volume.tets.reserve(static_cast<std::size_t>(out.numberoftetrahedra));
        const int corners = out.numberofcorners;
        for (int t = 0; t < out.numberoftetrahedra; ++t) {
            std::array<int, 4> tet = {
                out.tetrahedronlist[t * corners + 0],
                out.tetrahedronlist[t * corners + 1],
                out.tetrahedronlist[t * corners + 2],
                out.tetrahedronlist[t * corners + 3]};
            if (tet[0] == tet[1] || tet[0] == tet[2] || tet[0] == tet[3] || tet[1] == tet[2] ||
                tet[1] == tet[3] || tet[2] == tet[3])
                continue;
            volume.tets.push_back(tet);
        }

        struct Oriented {
            int a = 0;
            int b = 0;
            int c = 0;
            int count = 0;
            ChVector3d inside;
        };
        std::map<std::array<int, 3>, Oriented> boundaryFaces;
        auto addFace = [&](int a, int b, int c, const ChVector3d& inside) {
            std::array<int, 3> key = {a, b, c};
            std::sort(key.begin(), key.end());
            Oriented& slot = boundaryFaces[key];
            if (slot.count == 0) {
                slot.a = a;
                slot.b = b;
                slot.c = c;
                slot.inside = inside;
            }
            ++slot.count;
        };
        for (const auto& tet : volume.tets) {
            const ChVector3d inside = 0.25 * (volume.vertices[static_cast<std::size_t>(tet[0])] +
                                              volume.vertices[static_cast<std::size_t>(tet[1])] +
                                              volume.vertices[static_cast<std::size_t>(tet[2])] +
                                              volume.vertices[static_cast<std::size_t>(tet[3])]);
            addFace(tet[0], tet[2], tet[1], inside);
            addFace(tet[0], tet[1], tet[3], inside);
            addFace(tet[0], tet[3], tet[2], inside);
            addFace(tet[1], tet[2], tet[3], inside);
        }
        for (const auto& face : boundaryFaces) {
            if (face.second.count != 1)
                continue;
            int a = face.second.a;
            int b = face.second.b;
            int c = face.second.c;
            const ChVector3d& pa = volume.vertices[static_cast<std::size_t>(a)];
            const ChVector3d& pb = volume.vertices[static_cast<std::size_t>(b)];
            const ChVector3d& pc = volume.vertices[static_cast<std::size_t>(c)];
            const ChVector3d mid = (pa + pb + pc) / 3.0;
            const ChVector3d outward = (pb - pa).Cross(pc - pa);
            if (outward.Dot(mid - face.second.inside) < 0.0)
                std::swap(b, c);
            volume.boundary.push_back({a, b, c});
            int tag = 1;
            if (a < n && b < n && c < n) {
                std::array<int, 3> key = {a, b, c};
                std::sort(key.begin(), key.end());
                for (std::size_t f = 0; f < plc.triangles.size() && f < marker.size(); ++f) {
                    std::array<int, 3> plcKey = plc.triangles[f];
                    std::sort(plcKey.begin(), plcKey.end());
                    if (plcKey == key && marker[f] != 1) {
                        tag = marker[f];
                        break;
                    }
                }
            }
            volume.boundaryMarker.push_back(tag);
        }
        volume.elementSize = requested;
        if (worstTetQuality(volume) < 1.0e-5 && attempt + 1 < 4) {
            requested *= 1.4;
            relax *= 1.4;
            continue;
        }
        return !volume.tets.empty();
    }
    return false;
}

void snapToSolid(SurfaceMesh& shell, const VolumeMesh& volume, double tolerance) {
    if (volume.vertices.empty() || tolerance <= 0.0)
        return;
    const double tol2 = tolerance * tolerance;
    for (ChVector3d& vertex : shell.vertices) {
        double best = tol2;
        const ChVector3d* match = nullptr;
        for (const ChVector3d& solid : volume.vertices) {
            const double d = (vertex - solid).Length2();
            if (d < best) {
                best = d;
                match = &solid;
            }
        }
        if (match)
            vertex = *match;
    }
}

void refineShellEdges(SurfaceMesh& shell, const MeshSizing& sizing) {
    if (shell.triangles.empty())
        return;
    for (int pass = 0; pass < 4; ++pass) {
        if (shell.triangles.size() > 800)
            return;
        std::map<std::pair<int, int>, int> midpoints;
        auto midpoint = [&](int a, int b) {
            if (a > b)
                std::swap(a, b);
            const auto key = std::make_pair(a, b);
            const auto found = midpoints.find(key);
            if (found != midpoints.end())
                return found->second;
            const ChVector3d mid = 0.5 * (shell.vertices[static_cast<std::size_t>(a)] +
                                          shell.vertices[static_cast<std::size_t>(b)]);
            const int id = static_cast<int>(shell.vertices.size());
            shell.vertices.push_back(mid);
            midpoints.emplace(key, id);
            return id;
        };
        std::vector<std::array<int, 3>> next;
        std::vector<double> nextThickness;
        bool split = false;
        for (std::size_t t = 0; t < shell.triangles.size(); ++t) {
            const auto face = shell.triangles[t];
            const double thick = t < shell.thickness.size() ? shell.thickness[t] : 0.0;
            const double target = std::max(sizing.refineLimit > 1.0e-12 ? sizing.refineLimit : sizing.maxElement,
                                           sizing.minElement);
            double longest = 0.0;
            int edge = 0;
            for (int k = 0; k < 3; ++k) {
                const double len = (shell.vertices[static_cast<std::size_t>(face[k])] -
                                    shell.vertices[static_cast<std::size_t>(face[(k + 1) % 3])])
                                       .Length();
                if (len > longest) {
                    longest = len;
                    edge = k;
                }
            }
            if (longest <= target || longest < sizing.minElement * 2.0) {
                next.push_back(face);
                nextThickness.push_back(thick);
                continue;
            }
            const int i0 = face[edge];
            const int i1 = face[(edge + 1) % 3];
            const int i2 = face[(edge + 2) % 3];
            const int mid = midpoint(i0, i1);
            next.push_back({i0, mid, i2});
            next.push_back({mid, i1, i2});
            nextThickness.push_back(thick);
            nextThickness.push_back(thick);
            split = true;
        }
        shell.triangles.swap(next);
        shell.thickness.swap(nextThickness);
        if (!split)
            return;
    }
}

void sealBoundary(VolumeMesh& volume) {
    struct Oriented {
        int a = 0;
        int b = 0;
        int c = 0;
        int count = 0;
        ChVector3d inside;
    };
    std::map<std::array<int, 3>, Oriented> faces;
    auto addFace = [&](int a, int b, int c, const ChVector3d& inside) {
        std::array<int, 3> key = {a, b, c};
        std::sort(key.begin(), key.end());
        Oriented& slot = faces[key];
        if (slot.count == 0) {
            slot.a = a;
            slot.b = b;
            slot.c = c;
            slot.inside = inside;
        }
        ++slot.count;
    };
    for (const auto& tet : volume.tets) {
        const ChVector3d inside = 0.25 * (volume.vertices[static_cast<std::size_t>(tet[0])] +
                                          volume.vertices[static_cast<std::size_t>(tet[1])] +
                                          volume.vertices[static_cast<std::size_t>(tet[2])] +
                                          volume.vertices[static_cast<std::size_t>(tet[3])]);
        addFace(tet[0], tet[2], tet[1], inside);
        addFace(tet[0], tet[1], tet[3], inside);
        addFace(tet[0], tet[3], tet[2], inside);
        addFace(tet[1], tet[2], tet[3], inside);
    }
    volume.boundary.clear();
    for (const auto& face : faces) {
        if (face.second.count != 1)
            continue;
        int a = face.second.a;
        int b = face.second.b;
        int c = face.second.c;
        const ChVector3d& pa = volume.vertices[static_cast<std::size_t>(a)];
        const ChVector3d& pb = volume.vertices[static_cast<std::size_t>(b)];
        const ChVector3d& pc = volume.vertices[static_cast<std::size_t>(c)];
        const ChVector3d mid = (pa + pb + pc) / 3.0;
        const ChVector3d outward = (pb - pa).Cross(pc - pa);
        if (outward.Dot(mid - face.second.inside) < 0.0)
            std::swap(b, c);
        volume.boundary.push_back({a, b, c});
    }
}

void addOrientedTet(VolumeMesh& volume, int a, int b, int c, int d) {
    const ChVector3d& pa = volume.vertices[static_cast<std::size_t>(a)];
    const ChVector3d& pb = volume.vertices[static_cast<std::size_t>(b)];
    const ChVector3d& pc = volume.vertices[static_cast<std::size_t>(c)];
    const ChVector3d& pd = volume.vertices[static_cast<std::size_t>(d)];
    const double signedVolume = (pb - pa).Dot((pc - pa).Cross(pd - pa));
    if (std::abs(signedVolume) < 1.0e-18)
        return;
    if (signedVolume < 0.0)
        std::swap(c, d);
    volume.tets.push_back({a, b, c, d});
}

/// Tetrahedra through `thickness`. The base triangles stay put so a welded
/// seam that already coincides keeps its nodes. `layers` is the number of
/// hex-equivalent slices (each becomes 3 tets per triangle).
bool extrudeSheet(const SurfaceMesh& mid, double thickness, VolumeMesh& volume, int layers) {
    if (mid.triangles.empty() || thickness <= 1.0e-6)
        return false;
    layers = std::clamp(layers, 1, 16);
    volume = {};
    std::unordered_map<Key, int, KeyHash> weld;
    std::vector<ChVector3d> normals(mid.vertices.size(), ChVector3d(0, 0, 0));
    for (const auto& face : mid.triangles) {
        const ChVector3d n = faceNormal(mid, face) * faceArea(mid, face);
        for (int index : face)
            if (index >= 0)
                normals[static_cast<std::size_t>(index)] += n;
    }
    std::vector<std::vector<int>> sheet(static_cast<std::size_t>(layers + 1), std::vector<int>(mid.vertices.size(), -1));
    for (std::size_t i = 0; i < mid.vertices.size(); ++i) {
        ChVector3d n = normals[i];
        if (n.Length() < 1.0e-12)
            n = {0, 1, 0};
        n.Normalize();
        for (int layer = 0; layer <= layers; ++layer) {
            const double offset = thickness * static_cast<double>(layer) / layers;
            sheet[static_cast<std::size_t>(layer)][i] =
                weldVertex(volume.vertices, weld, mid.vertices[i] + n * offset, 1.0e-6);
        }
    }
    for (const auto& face : mid.triangles) {
        for (int layer = 0; layer < layers; ++layer) {
            const auto& bot = sheet[static_cast<std::size_t>(layer)];
            const auto& top = sheet[static_cast<std::size_t>(layer + 1)];
            const int a = bot[static_cast<std::size_t>(face[0])];
            const int b = bot[static_cast<std::size_t>(face[1])];
            const int c = bot[static_cast<std::size_t>(face[2])];
            const int d = top[static_cast<std::size_t>(face[0])];
            const int e = top[static_cast<std::size_t>(face[1])];
            const int f = top[static_cast<std::size_t>(face[2])];
            if (a == b || b == c || c == a)
                continue;
            addOrientedTet(volume, a, b, c, f);
            addOrientedTet(volume, a, b, f, e);
            addOrientedTet(volume, a, e, f, d);
        }
    }
    sealBoundary(volume);
    volume.elementSize = thickness;
    return !volume.tets.empty();
}

SurfaceMesh sheetFromPatch(const PlatePatch& patch) {
    SurfaceMesh mid;
    const int nu = std::max(1, patch.divU);
    const int nv = std::max(1, patch.divV);
    ChVector3d u = patch.uAxis;
    ChVector3d v = patch.vAxis;
    if (u.Length() < 1.0e-12 || v.Length() < 1.0e-12)
        return mid;
    u.Normalize();
    v.Normalize();
    const int stride = nv + 1;
    mid.vertices.reserve(static_cast<std::size_t>(nu + 1) * stride);
    for (int iu = 0; iu <= nu; ++iu) {
        for (int iv = 0; iv <= nv; ++iv) {
            const double su = static_cast<double>(iu) / nu;
            const double sv = static_cast<double>(iv) / nv;
            mid.vertices.push_back(patch.origin + u * (su * patch.uSize) + v * (sv * patch.vSize));
        }
    }
    for (int iu = 0; iu < nu; ++iu) {
        for (int iv = 0; iv < nv; ++iv) {
            const int i00 = iu * stride + iv;
            const int i10 = (iu + 1) * stride + iv;
            const int i11 = (iu + 1) * stride + iv + 1;
            const int i01 = iu * stride + iv + 1;
            mid.triangles.push_back({i00, i10, i11});
            mid.triangles.push_back({i00, i11, i01});
        }
    }
    return mid;
}

int extrusionLayers(double thickness, double h) {
    const double edge = std::max(h, 1.0e-9);
    const int layers = std::max(1, static_cast<int>(std::ceil(thickness / edge)));
    return std::min(layers, 16);
}

int sharedCount(const std::array<int, 3>& a, const std::array<int, 3>& b) {
    int hits = 0;
    for (int i : a)
        for (int j : b)
            hits += i == j ? 1 : 0;
    return hits;
}

double faceLongestEdge(const SurfaceMesh& mesh, const std::array<int, 3>& face) {
    double longest = 0.0;
    for (int k = 0; k < 3; ++k) {
        const ChVector3d& a = mesh.vertices[static_cast<std::size_t>(face[k])];
        const ChVector3d& b = mesh.vertices[static_cast<std::size_t>(face[(k + 1) % 3])];
        longest = std::max(longest, (b - a).Length());
    }
    return longest;
}

/// Characteristic in-plane width: triangle altitude, not the long span of a
/// bar-side face. A 0.8 m square rod is not a 12 cm plate.
double faceInplaneWidth(const SurfaceMesh& mesh, const std::array<int, 3>& face) {
    const double longest = std::max(faceLongestEdge(mesh, face), 1.0e-12);
    return 2.0 * std::max(faceArea(mesh, face), 0.0) / longest;
}

void measureFaceFields(const SurfaceMesh& source, std::vector<double>& thickness,
                       std::vector<double>& curvature, std::vector<ChVector3d>& normals,
                       std::vector<int>& opposite) {
    const double inward = meshVolume(source) >= 0.0 ? -1.0 : 1.0;
    thickness.assign(source.triangles.size(), 1.0e30);
    curvature.assign(source.triangles.size(), 0.0);
    opposite.assign(source.triangles.size(), -1);
    normals.resize(source.triangles.size());
    for (std::size_t i = 0; i < source.triangles.size(); ++i)
        normals[i] = faceNormal(source, source.triangles[i]);
    for (std::size_t i = 0; i < source.triangles.size(); ++i) {
        const auto& face = source.triangles[i];
        if (normals[i].Length() < 0.5)
            continue;
        const ChVector3d origin = (source.vertices[static_cast<std::size_t>(face[0])] +
                                   source.vertices[static_cast<std::size_t>(face[1])] +
                                   source.vertices[static_cast<std::size_t>(face[2])]) /
                                  3.0;
        const ChVector3d direction = normals[i] * inward;
        RayHit best;
        for (std::size_t j = 0; j < source.triangles.size(); ++j) {
            if (i == j)
                continue;
            if (normals[i].Dot(normals[j]) > -0.2)
                continue;
            const auto& other = source.triangles[j];
            double t = 0.0;
            if (!rayTriangle(origin, direction, source.vertices[static_cast<std::size_t>(other[0])],
                             source.vertices[static_cast<std::size_t>(other[1])],
                             source.vertices[static_cast<std::size_t>(other[2])], t))
                continue;
            if (t < best.distance) {
                best.distance = t;
                best.face = static_cast<int>(j);
            }
        }
        thickness[i] = best.distance;
        opposite[i] = best.face;
    }
    for (std::size_t i = 0; i < source.triangles.size(); ++i) {
        for (std::size_t j = i + 1; j < source.triangles.size(); ++j) {
            if (sharedCount(source.triangles[i], source.triangles[j]) != 2)
                continue;
            double edge = 1.0e-9;
            const auto& fa = source.triangles[i];
            const auto& fb = source.triangles[j];
            for (int k = 0; k < 3; ++k) {
                const int u = fa[k];
                const int v = fa[(k + 1) % 3];
                for (int m = 0; m < 3; ++m) {
                    const int p = fb[m];
                    const int q = fb[(m + 1) % 3];
                    if ((u == p && v == q) || (u == q && v == p)) {
                        edge = std::max(edge, (source.vertices[static_cast<std::size_t>(u)] -
                                               source.vertices[static_cast<std::size_t>(v)])
                                                  .Length());
                    }
                }
            }
            const double kappa = (normals[i] - normals[j]).Length() / edge;
            curvature[i] = std::max(curvature[i], kappa);
            curvature[j] = std::max(curvature[j], kappa);
        }
    }
}

void gradateSizeField(const SurfaceMesh& mesh, std::vector<double>& size);

std::vector<double> vertexSizeField(const SurfaceMesh& source, const std::vector<double>& thickness,
                                    const std::vector<double>& curvature, double hTarget,
                                    const MeshSizing& sizing) {
    std::vector<double> size(source.vertices.size(), 0.0);
    std::vector<double> weight(source.vertices.size(), 0.0);
    for (std::size_t i = 0; i < source.triangles.size(); ++i) {
        const double t = i < thickness.size() ? thickness[i] : 1.0e30;
        const double kappa = i < curvature.size() ? curvature[i] : 0.0;
        double h = geometricEdgeRequest(t < 1.0e20 ? t : 0.0, kappa);
        if (h > 1.0e29)
            h = hTarget;
        const double area = std::max(faceArea(source, source.triangles[i]), 1.0e-18);
        for (int index : source.triangles[i]) {
            if (index < 0)
                continue;
            size[static_cast<std::size_t>(index)] += h * area;
            weight[static_cast<std::size_t>(index)] += area;
        }
    }
    for (std::size_t i = 0; i < size.size(); ++i)
        size[i] = weight[i] > 1.0e-18 ? size[i] / weight[i] : hTarget;
    clampSizeField(size, hTarget, sizing.minElement,
                   sizing.maxElement > sizing.minElement ? sizing.maxElement : hTarget);
    gradateSizeField(source, size);
    return size;
}

/// One thin-wall test. Span is the larger in-plane size, not the diagonal, and
/// not a mesh-size knob. The same cutoff is used for patches, closed bodies,
/// open sheets, and weld seeding.
constexpr double kThinWallSlenderness = 0.15;

bool isSlenderWall(double thickness, double span) {
    if (thickness >= 1.0e20 || thickness <= 1.0e-12)
        return false;
    return thickness / std::max(span, 1.0e-9) < kThinWallSlenderness;
}

bool faceIsThinWall(double thickness, double inplane) {
    return isSlenderWall(thickness, inplane);
}

/// Small-strain shells are legal for this material. Hyperelastic laws are not:
/// a patch is refused and a closed body is tet-filled instead.
bool linearShellAllowed(const Material& material) {
    return compatibility(material, ElementFamily::ReissnerShell).allowed &&
           compatibility(material, ElementFamily::KirchhoffTriShell).allowed;
}

/// In-plane extent of each face's connected region of similar ray thickness.
/// A curved CAD skin is classified by that patch, not by one triangle's altitude.
std::vector<double> regionSpans(const SurfaceMesh& mesh, const std::vector<double>& thickness) {
    const std::size_t count = mesh.triangles.size();
    std::vector<double> span(count, 0.0);
    if (count == 0)
        return span;
    using Edge = std::pair<int, int>;
    auto canon = [](int a, int b) { return a < b ? Edge{a, b} : Edge{b, a}; };
    std::map<Edge, std::vector<int>> edgeFaces;
    for (std::size_t i = 0; i < count; ++i) {
        const auto& face = mesh.triangles[i];
        for (int k = 0; k < 3; ++k)
            edgeFaces[canon(face[k], face[(k + 1) % 3])].push_back(static_cast<int>(i));
    }
    std::vector<std::vector<int>> adj(count);
    for (const auto& entry : edgeFaces) {
        const auto& faces = entry.second;
        for (std::size_t a = 0; a < faces.size(); ++a) {
            for (std::size_t b = a + 1; b < faces.size(); ++b) {
                adj[static_cast<std::size_t>(faces[a])].push_back(faces[b]);
                adj[static_cast<std::size_t>(faces[b])].push_back(faces[a]);
            }
        }
    }
    auto similar = [&](std::size_t a, std::size_t b) {
        const double ta = a < thickness.size() ? thickness[a] : 1.0e30;
        const double tb = b < thickness.size() ? thickness[b] : 1.0e30;
        if (ta >= 1.0e20 || tb >= 1.0e20 || ta <= 1.0e-12 || tb <= 1.0e-12)
            return false;
        const double ratio = ta > tb ? ta / tb : tb / ta;
        return ratio <= 1.5;
    };
    std::vector<int> region(count, -1);
    for (std::size_t seed = 0; seed < count; ++seed) {
        if (region[seed] >= 0)
            continue;
        const int id = static_cast<int>(seed);
        std::vector<int> stack{static_cast<int>(seed)};
        std::vector<int> faces;
        region[seed] = id;
        while (!stack.empty()) {
            const int face = stack.back();
            stack.pop_back();
            faces.push_back(face);
            for (int next : adj[static_cast<std::size_t>(face)]) {
                if (next < 0 || region[static_cast<std::size_t>(next)] >= 0)
                    continue;
                if (!similar(static_cast<std::size_t>(face), static_cast<std::size_t>(next)))
                    continue;
                region[static_cast<std::size_t>(next)] = id;
                stack.push_back(next);
            }
        }
        ChVector3d normal{0, 0, 0};
        std::vector<int> used;
        for (int face : faces) {
            const auto& tri = mesh.triangles[static_cast<std::size_t>(face)];
            normal += faceNormal(mesh, tri) * std::max(faceArea(mesh, tri), 0.0);
            for (int index : tri)
                used.push_back(index);
        }
        std::sort(used.begin(), used.end());
        used.erase(std::unique(used.begin(), used.end()), used.end());
        if (normal.Length2() < 1.0e-16)
            normal = {0, 1, 0};
        else
            normal.Normalize();
        ChVector3d e1 = std::abs(normal.x()) < 0.9 ? normal.Cross(ChVector3d{1, 0, 0}) : normal.Cross(ChVector3d{0, 1, 0});
        if (e1.Length2() < 1.0e-16)
            e1 = {1, 0, 0};
        else
            e1.Normalize();
        const ChVector3d e2 = normal.Cross(e1);
        double uMin = 1.0e30;
        double uMax = -1.0e30;
        double vMin = 1.0e30;
        double vMax = -1.0e30;
        for (int index : used) {
            if (index < 0 || static_cast<std::size_t>(index) >= mesh.vertices.size())
                continue;
            const ChVector3d& p = mesh.vertices[static_cast<std::size_t>(index)];
            const double u = p.Dot(e1);
            const double v = p.Dot(e2);
            uMin = std::min(uMin, u);
            uMax = std::max(uMax, u);
            vMin = std::min(vMin, v);
            vMax = std::max(vMax, v);
        }
        const double du = uMax >= uMin ? uMax - uMin : 0.0;
        const double dv = vMax >= vMin ? vMax - vMin : 0.0;
        const double extent = std::max(du, dv);
        for (int face : faces)
            span[static_cast<std::size_t>(face)] = std::max(extent, 1.0e-9);
    }
    return span;
}

/// Slender faces become shells. Faces within a few thicknesses of a stocky
/// neighbour stay tets so the shell ties into a transition band.
void classifyWalls(const SurfaceMesh& mesh, const std::vector<double>& thickness, double minElement,
                   std::vector<char>& shell, std::vector<char>& band) {
    const std::size_t count = mesh.triangles.size();
    const std::vector<double> span = regionSpans(mesh, thickness);
    std::vector<char> slender(count, 0);
    std::vector<char> stocky(count, 0);
    shell.assign(count, 0);
    band.assign(count, 0);
    for (std::size_t i = 0; i < count; ++i) {
        const double t = i < thickness.size() ? thickness[i] : 1.0e30;
        const double width = i < span.size() ? span[i] : 0.0;
        if (faceIsThinWall(t, width))
            slender[i] = 1;
        else if (t < 1.0e20)
            stocky[i] = 1;
    }
    using Edge = std::pair<int, int>;
    auto canon = [](int a, int b) { return a < b ? Edge{a, b} : Edge{b, a}; };
    std::map<Edge, std::vector<int>> edgeFaces;
    for (std::size_t i = 0; i < count; ++i) {
        const auto& face = mesh.triangles[i];
        for (int k = 0; k < 3; ++k)
            edgeFaces[canon(face[k], face[(k + 1) % 3])].push_back(static_cast<int>(i));
    }
    std::vector<std::vector<int>> adj(count);
    for (const auto& entry : edgeFaces) {
        for (std::size_t a = 0; a < entry.second.size(); ++a) {
            for (std::size_t b = a + 1; b < entry.second.size(); ++b) {
                adj[static_cast<std::size_t>(entry.second[a])].push_back(entry.second[b]);
                adj[static_cast<std::size_t>(entry.second[b])].push_back(entry.second[a]);
            }
        }
    }
    std::vector<int> stack;
    std::vector<double> reach(count, 1.0e30);
    for (std::size_t i = 0; i < count; ++i) {
        if (!slender[i])
            continue;
        for (int next : adj[i]) {
            if (next >= 0 && stocky[static_cast<std::size_t>(next)]) {
                band[i] = 1;
                reach[i] = 0.0;
                stack.push_back(static_cast<int>(i));
                break;
            }
        }
    }
    for (std::size_t cursor = 0; cursor < stack.size(); ++cursor) {
        const int face = stack[cursor];
        const double here = reach[static_cast<std::size_t>(face)];
        const double limit = 3.0 * std::max(face < static_cast<int>(thickness.size()) ? thickness[static_cast<std::size_t>(face)] : 0.0, minElement);
        for (int next : adj[static_cast<std::size_t>(face)]) {
            if (next < 0 || !slender[static_cast<std::size_t>(next)])
                continue;
            const auto& tri = mesh.triangles[static_cast<std::size_t>(face)];
            const double step = faceLongestEdge(mesh, tri);
            if (here + step > limit || here + step >= reach[static_cast<std::size_t>(next)])
                continue;
            reach[static_cast<std::size_t>(next)] = here + step;
            band[static_cast<std::size_t>(next)] = 1;
            stack.push_back(next);
        }
    }
    for (std::size_t i = 0; i < count; ++i)
        shell[i] = slender[i] && !band[i] ? 1 : 0;
}

/// Neighboring vertices may differ by at most this factor, so a fine web
/// cannot sit against a tet orders of magnitude larger.
void gradateSizeField(const SurfaceMesh& mesh, std::vector<double>& size) {
    if (size.size() < 2)
        return;
    // About three per edge. Tighter than this pulls a boss down to the web and
    // the resolution ceiling can no longer shrink the boss.
    constexpr double kRatio = 3.0;
    std::vector<std::vector<int>> adj(size.size());
    auto link = [&](int a, int b) {
        if (a == b || a < 0 || b < 0)
            return;
        if (static_cast<std::size_t>(a) >= size.size() || static_cast<std::size_t>(b) >= size.size())
            return;
        adj[static_cast<std::size_t>(a)].push_back(b);
        adj[static_cast<std::size_t>(b)].push_back(a);
    };
    for (const auto& tri : mesh.triangles) {
        link(tri[0], tri[1]);
        link(tri[1], tri[2]);
        link(tri[2], tri[0]);
    }
    for (int pass = 0; pass < 48; ++pass) {
        bool changed = false;
        for (std::size_t i = 0; i < size.size(); ++i) {
            for (int next : adj[i]) {
                const double cap = size[i] * kRatio;
                if (size[static_cast<std::size_t>(next)] > cap) {
                    size[static_cast<std::size_t>(next)] = cap;
                    changed = true;
                }
            }
        }
        if (!changed)
            break;
    }
}

SurfaceMesh midSurfaceOf(const SurfaceMesh& source, const std::vector<double>& thickness,
                         const std::vector<char>& keep, const std::vector<ChVector3d>& normals) {
    SurfaceMesh mid;
    const double inward = meshVolume(source) >= 0.0 ? -1.0 : 1.0;
    std::vector<char> interface(source.vertices.size(), 0);
    std::vector<char> used(source.vertices.size(), 0);
    for (std::size_t i = 0; i < source.triangles.size(); ++i) {
        const bool thin = i < keep.size() && keep[i];
        for (int index : source.triangles[i]) {
            if (index < 0)
                continue;
            if (thin)
                used[static_cast<std::size_t>(index)] = 1;
        }
    }
    for (std::size_t i = 0; i < source.triangles.size(); ++i) {
        const bool thin = i < keep.size() && keep[i];
        if (thin)
            continue;
        for (int index : source.triangles[i]) {
            if (index < 0)
                continue;
            if (used[static_cast<std::size_t>(index)])
                interface[static_cast<std::size_t>(index)] = 1;
        }
    }
    std::vector<ChVector3d> offset(source.vertices.size(), ChVector3d(0, 0, 0));
    std::vector<double> weight(source.vertices.size(), 0.0);
    for (std::size_t i = 0; i < source.triangles.size(); ++i) {
        if (i >= keep.size() || !keep[i])
            continue;
        if (normals[i].Length() < 0.5)
            continue;
        const double t = i < thickness.size() && thickness[i] < 1.0e20 ? thickness[i] : 0.0;
        const ChVector3d n = normals[i] * inward;
        for (int index : source.triangles[i]) {
            if (index < 0)
                continue;
            offset[static_cast<std::size_t>(index)] += n * (0.5 * t);
            weight[static_cast<std::size_t>(index)] += 1.0;
        }
    }
    std::vector<int> remap(source.vertices.size(), -1);
    auto take = [&](int index) {
        if (index < 0)
            return -1;
        if (remap[static_cast<std::size_t>(index)] >= 0)
            return remap[static_cast<std::size_t>(index)];
        ChVector3d p = source.vertices[static_cast<std::size_t>(index)];
        if (!interface[static_cast<std::size_t>(index)] &&
            weight[static_cast<std::size_t>(index)] > 0.5)
            p += offset[static_cast<std::size_t>(index)] / weight[static_cast<std::size_t>(index)];
        remap[static_cast<std::size_t>(index)] = static_cast<int>(mid.vertices.size());
        mid.vertices.push_back(p);
        return remap[static_cast<std::size_t>(index)];
    };
    for (std::size_t i = 0; i < source.triangles.size(); ++i) {
        if (i >= keep.size() || !keep[i])
            continue;
        const auto& face = source.triangles[i];
        const int a = take(face[0]);
        const int b = take(face[1]);
        const int c = take(face[2]);
        if (a < 0 || b < 0 || c < 0 || a == b || b == c || c == a)
            continue;
        mid.triangles.push_back({a, b, c});
        mid.thickness.push_back(i < thickness.size() && thickness[i] < 1.0e20 ? thickness[i] : 0.0);
    }
    return mid;
}

bool capInterface(const SurfaceMesh& source, const std::vector<char>& thin, SurfaceMesh& bulk) {
    bulk.vertices = source.vertices;
    bulk.triangles.clear();
    using Edge = std::pair<int, int>;
    auto canon = [](int a, int b) { return a < b ? Edge{a, b} : Edge{b, a}; };
    std::map<Edge, std::vector<int>> edgeFaces;
    for (std::size_t i = 0; i < source.triangles.size(); ++i) {
        const auto& face = source.triangles[i];
        if (!(i < thin.size() && thin[i]))
            bulk.triangles.push_back(face);
        for (int k = 0; k < 3; ++k)
            edgeFaces[canon(face[k], face[(k + 1) % 3])].push_back(static_cast<int>(i));
    }
    std::map<int, std::vector<int>> adj;
    int interfaceEdges = 0;
    for (const auto& entry : edgeFaces) {
        bool hasThin = false;
        bool hasBulk = false;
        for (int f : entry.second) {
            if (static_cast<std::size_t>(f) < thin.size() && thin[static_cast<std::size_t>(f)])
                hasThin = true;
            else
                hasBulk = true;
        }
        if (!hasThin || !hasBulk)
            continue;
        ++interfaceEdges;
        adj[entry.first.first].push_back(entry.first.second);
        adj[entry.first.second].push_back(entry.first.first);
    }
    if (interfaceEdges == 0)
        return !bulk.triangles.empty();
    std::set<Edge> used;
    int loops = 0;
    for (const auto& start : adj) {
        for (int neighbor : start.second) {
            const Edge first = canon(start.first, neighbor);
            if (used.count(first))
                continue;
            std::vector<int> loop = {start.first};
            int prev = start.first;
            int cur = neighbor;
            used.insert(first);
            bool closed = false;
            for (int step = 0; step < static_cast<int>(source.vertices.size()) + 2; ++step) {
                if (cur == start.first) {
                    closed = true;
                    break;
                }
                loop.push_back(cur);
                int next = -1;
                for (int n : adj[cur]) {
                    if (n == prev)
                        continue;
                    if (!used.count(canon(cur, n))) {
                        next = n;
                        break;
                    }
                }
                if (next < 0)
                    break;
                used.insert(canon(cur, next));
                prev = cur;
                cur = next;
            }
            if (!closed || loop.size() < 3)
                continue;
            ++loops;
            for (std::size_t i = 1; i + 1 < loop.size(); ++i) {
                const int a = loop[0];
                const int b = loop[i];
                const int c = loop[i + 1];
                if (a == b || b == c || c == a)
                    continue;
                bulk.triangles.push_back({a, b, c});
            }
        }
    }
    return loops > 0 && bulk.triangles.size() >= 4;
}

bool meshClosedBody(PartDefinition& part, const SurfaceMesh& source, double h) {
    std::vector<double> thickness;
    std::vector<double> curvature;
    std::vector<ChVector3d> normals;
    std::vector<int> opposite;
    measureFaceFields(source, thickness, curvature, normals, opposite);
    const std::vector<double> size = vertexSizeField(source, thickness, curvature, h, part.sizing);
    std::vector<char> thin(source.triangles.size(), 0);
    std::vector<char> band;
    std::vector<char> dropped(source.triangles.size(), 0);
    int thinCount = 0;
    int bulkCount = 0;
    double thinArea = 0.0;
    double bulkArea = 0.0;
    // Shell versus solid is geometry only. throughThickness is how many layers
    // a solid extrusion asked for; it does not skip this test, and the
    // material does not either. A finite-strain material is stripped back to
    // a tet fill after the geometric decision, below.
    {
        std::vector<char> shellFace;
        classifyWalls(source, thickness, part.sizing.minElement, shellFace, band);
        (void)band;
        for (std::size_t i = 0; i < source.triangles.size(); ++i) {
            const double area = std::max(faceArea(source, source.triangles[i]), 0.0);
            const bool wall = i < shellFace.size() && shellFace[i];
            if (wall) {
                thin[i] = 1;
                ++thinCount;
                thinArea += area;
            } else {
                ++bulkCount;
                bulkArea += area;
            }
        }
        // Opposite skins of the same wall both offset to the mid-plane and
        // become two coincident BST layers. Keep one skin per pair; the
        // discarded skin is not a solid to tet-fill either.
        for (std::size_t i = 0; i < source.triangles.size(); ++i) {
            if (!thin[i])
                continue;
            const int hit = i < opposite.size() ? opposite[i] : -1;
            if (hit < 0 || static_cast<std::size_t>(hit) >= thin.size() || !thin[static_cast<std::size_t>(hit)])
                continue;
            if (static_cast<int>(i) > hit) {
                thin[i] = 0;
                dropped[i] = 1;
                --thinCount;
                thinArea -= std::max(faceArea(source, source.triangles[i]), 0.0);
            }
        }
        bulkCount = 0;
        bulkArea = 0.0;
        for (std::size_t i = 0; i < source.triangles.size(); ++i) {
            if (thin[i] || dropped[i])
                continue;
            ++bulkCount;
            bulkArea += std::max(faceArea(source, source.triangles[i]), 0.0);
        }
        // Rim faces of a plate are bulk by the local test but they are not a
        // solid to tet. If almost all area is thin-wall, the part is a shell.
        if (thinCount > 0 && thinArea > 3.0 * std::max(bulkArea, 1.0e-18))
            bulkCount = 0;
    }

    auto finishVolume = [&]() {
        part.thickness = part.volume.elementSize;
        part.sizing.cached = true;
        ++part.sizing.generation;
        return true;
    };

    if (!source.vertices.empty()) {
        ChVector3d lo = source.vertices.front();
        ChVector3d hi = lo;
        for (const ChVector3d& p : source.vertices) {
            lo = {std::min(lo.x(), p.x()), std::min(lo.y(), p.y()), std::min(lo.z(), p.z())};
            hi = {std::max(hi.x(), p.x()), std::max(hi.y(), p.y()), std::max(hi.z(), p.z())};
        }
        double extent[3] = {hi.x() - lo.x(), hi.y() - lo.y(), hi.z() - lo.z()};
        std::sort(extent, extent + 3);
        // A plate has one small dimension and two comparable large ones. A bar
        // has two small dimensions, so it stays a tet body. Span is the middle
        // extent, the same larger in-plane size a square patch uses.
        if (isSlenderWall(extent[0], extent[1]) && extent[1] / std::max(extent[2], 1.0e-9) > 0.5) {
            thin.assign(source.triangles.size(), 1);
            dropped.assign(source.triangles.size(), 0);
            thinCount = static_cast<int>(source.triangles.size());
            bulkCount = 0;
            for (double& t : thickness) {
                if (t > 2.0 * extent[0])
                    t = extent[0];
            }
            part.thickness = extent[0];
        }
    }

    // Geometry may call faces shells. A finite-strain material has no linear
    // shell, so a closed body is tet-filled instead of mid-surfaced.
    if (!linearShellAllowed(part.material)) {
        thin.assign(source.triangles.size(), 0);
        dropped.assign(source.triangles.size(), 0);
        thinCount = 0;
        bulkCount = static_cast<int>(source.triangles.size());
    }

    if (thinCount > 0 && bulkCount == 0) {
        part.volume = {};
        double thick = part.thickness;
        int counted = 0;
        double sum = 0.0;
        for (std::size_t i = 0; i < thickness.size(); ++i) {
            if (i < thin.size() && !thin[i])
                continue;
            if (thickness[i] > 1.0e-12 && thickness[i] < 1.0e20) {
                sum += thickness[i];
                ++counted;
            }
        }
        if (counted > 0)
            thick = sum / static_cast<double>(counted);
        thick = std::max(thick, part.sizing.minElement);
        // A rectangular plate is one Reissner midsurface. Triangle BST shells
        // tied by point constraints do not converge under HHT: Newton exits
        // and the panel never meets the floor. The default 1 m patch is not
        // an authored plate; only a patch that covers this solid is.
        bool usePatch = part.patch.uSize > 1.0e-6 && part.patch.vSize > 1.0e-6 && !source.vertices.empty();
        if (usePatch) {
            ChVector3d boxLo = source.vertices.front();
            ChVector3d boxHi = boxLo;
            for (const ChVector3d& p : source.vertices) {
                boxLo = {std::min(boxLo.x(), p.x()), std::min(boxLo.y(), p.y()), std::min(boxLo.z(), p.z())};
                boxHi = {std::max(boxHi.x(), p.x()), std::max(boxHi.y(), p.y()), std::max(boxHi.z(), p.z())};
            }
            double span[3] = {boxHi.x() - boxLo.x(), boxHi.y() - boxLo.y(), boxHi.z() - boxLo.z()};
            std::sort(span, span + 3);
            const double plateArea = std::max(span[1] * span[2], 1.0e-12);
            const double patchArea = part.patch.uSize * part.patch.vSize;
            const double pad = std::max({part.patch.uSize, part.patch.vSize, 1.0e-3});
            const ChVector3d& origin = part.patch.origin;
            const bool areaOk = patchArea > 0.25 * plateArea && patchArea < 4.0 * plateArea;
            const bool originOk = origin.x() >= boxLo.x() - pad && origin.x() <= boxHi.x() + pad &&
                                  origin.y() >= boxLo.y() - pad && origin.y() <= boxHi.y() + pad &&
                                  origin.z() >= boxLo.z() - pad && origin.z() <= boxHi.z() + pad;
            usePatch = areaOk && originOk;
        }
        if (usePatch) {
            double localH = std::max(h, part.sizing.minElement);
            int divU = std::max(1, static_cast<int>(std::ceil(part.patch.uSize / localH)));
            int divV = std::max(1, static_cast<int>(std::ceil(part.patch.vSize / localH)));
            const int cap = tetBudgetCap(part.sizing, part.material);
            while (divU * divV > cap && (divU > 1 || divV > 1)) {
                localH *= 1.25;
                divU = std::max(1, static_cast<int>(std::ceil(part.patch.uSize / localH)));
                divV = std::max(1, static_cast<int>(std::ceil(part.patch.vSize / localH)));
            }
            part.patch.divU = divU;
            part.patch.divV = divV;
            part.surface = {};
            part.volume = {};
            part.thickness = thick;
            part.sizing.throughThickness = 0;
            part.sizing.cached = true;
            ++part.sizing.generation;
            return true;
        } else {
            part.surface = midSurfaceOf(source, thickness, thin, normals);
            if (part.surface.triangles.empty())
                part.surface = source;
            refineShellEdges(part.surface, part.sizing);
            counted = 0;
            sum = 0.0;
            for (double t : part.surface.thickness) {
                if (t > 1.0e-12 && t < 1.0e20) {
                    sum += t;
                    ++counted;
                }
            }
            if (counted > 0)
                thick = sum / static_cast<double>(counted);
            thick = std::max(thick, part.sizing.minElement);
        }
        part.thickness = thick;
        part.sizing.throughThickness = 0;
        part.sizing.cached = true;
        ++part.sizing.generation;
        return !part.surface.empty();
    }

    if (thinCount > 0 && bulkCount > 0 && thinArea <= bulkArea) {
        SurfaceMesh bulk;
        SurfaceMesh shell = midSurfaceOf(source, thickness, thin, normals);
        std::vector<char> capKeep = thin;
        for (std::size_t i = 0; i < dropped.size(); ++i) {
            if (dropped[i])
                capKeep[i] = 1;
        }
        if (capInterface(source, capKeep, bulk) && fillSolid(bulk, size, part.volume, part.sizing)) {
            refineShellEdges(shell, part.sizing);
            snapToSolid(shell, part.volume, 1.0e-8);
            part.surface = std::move(shell);
            return finishVolume();
        }
        part.surface = {};
        part.volume = {};
    } else if (thinCount > 0 && thinArea > bulkArea) {
        bulkCount = 0;
    }

    if (thinCount > 0 && bulkCount == 0) {
        part.volume = {};
        double thick = part.thickness;
        int counted = 0;
        double sum = 0.0;
        for (std::size_t i = 0; i < thickness.size(); ++i) {
            if (i < thin.size() && !thin[i] && !(i < dropped.size() && dropped[i]))
                continue;
            if (thickness[i] > 1.0e-12 && thickness[i] < 1.0e20) {
                sum += thickness[i];
                ++counted;
            }
        }
        if (counted > 0)
            thick = sum / static_cast<double>(counted);
        thick = std::max(thick, part.sizing.minElement);
        part.surface = midSurfaceOf(source, thickness, thin, normals);
        if (part.surface.triangles.empty())
            part.surface = source;
        if (part.surface.thickness.size() != part.surface.triangles.size()) {
            part.surface.thickness.assign(part.surface.triangles.size(), thick);
        }
        part.thickness = thick;
        part.sizing.throughThickness = 0;
        part.sizing.cached = true;
        ++part.sizing.generation;
        return !part.surface.empty();
    }

    if (!fillSolid(source, size, part.volume, part.sizing))
        return false;
    part.surface = {};
    return finishVolume();
}

}  // namespace

void applyResolutionSizing(PartDefinition& part, double resolution) {
    const double scale = std::max(resolution, 1.0e-6);
    const int base = std::max(200, part.sizing.tetBudget);
    part.sizing.tetBudget = std::max(200, static_cast<int>(std::lround(static_cast<double>(base) * scale)));
    const double h = targetEdge(part);
    part.sizing.maxElement = std::max(h, part.sizing.minElement);
    part.sizing.refineLimit = h;
    part.sizing.cached = false;
    part.volume = {};
    part.surface = {};
}

bool ensurePartMesh(PartDefinition& part) {
    if (part.sizing.cached)
        return true;
    const double h = std::max(targetEdge(part), 1.0e-6);
    auto keepAsShell = [&](double thickness, double localH, int divU, int divV) {
        part.patch.divU = divU;
        part.patch.divV = divV;
        part.sizing.maxElement = std::max(localH, part.sizing.minElement);
        part.volume = {};
        part.surface = {};
        part.thickness = thickness;
        part.sizing.throughThickness = 0;
        part.sizing.cached = true;
        ++part.sizing.generation;
        return true;
    };
    auto fillFromPatch = [&](bool centerOnPatch) {
        if (part.patch.uSize <= 1.0e-6 || part.patch.vSize <= 1.0e-6)
            return false;
        const double thickness = std::max(part.thickness, part.sizing.minElement);
        double localH = h;
        int layers = extrusionLayers(thickness, h);
        int divU = std::max(1, static_cast<int>(std::ceil(part.patch.uSize / localH)));
        int divV = std::max(1, static_cast<int>(std::ceil(part.patch.vSize / localH)));
        const int cap = tetBudgetCap(part.sizing, part.material);
        const double span = std::max(part.patch.uSize, part.patch.vSize);
        const bool slender = isSlenderWall(thickness, span);
        // A slender patch is a Reissner shell when that element is legal.
        // Finite-strain materials have no such shell: refuse the patch rather
        // than extrude locking tets. A bulky patch stays a solid even when the
        // tet budget is tight; coarsen it, do not flip it back to a shell.
        if (slender) {
            if (!linearShellAllowed(part.material))
                return false;
            return keepAsShell(thickness, localH, divU, divV);
        }
        layers = std::max(layers, 3);
        while (6 * divU * divV * layers > cap && (divU > 1 || divV > 1)) {
            localH *= 1.25;
            divU = std::max(1, static_cast<int>(std::ceil(part.patch.uSize / localH)));
            divV = std::max(1, static_cast<int>(std::ceil(part.patch.vSize / localH)));
        }
        while (6 * divU * divV * layers > cap && layers > 1)
            --layers;
        part.patch.divU = divU;
        part.patch.divV = divV;
        part.sizing.maxElement = std::max(localH, part.sizing.minElement);
        SurfaceMesh mid = sheetFromPatch(part.patch);
        if (centerOnPatch && !mid.vertices.empty()) {
            chrono::ChVector3d u = part.patch.uAxis;
            chrono::ChVector3d v = part.patch.vAxis;
            if (u.Length() > 1.0e-12)
                u.Normalize();
            if (v.Length() > 1.0e-12)
                v.Normalize();
            chrono::ChVector3d n = u.Cross(v);
            if (n.Length() < 1.0e-12)
                n = {0, 1, 0};
            else
                n.Normalize();
            const chrono::ChVector3d shift = n * (-0.5 * thickness);
            for (auto& p : mid.vertices)
                p += shift;
        }
        if (!extrudeSheet(mid, thickness, part.volume, layers))
            return false;
        part.surface = {};
        part.thickness = thickness;
        part.sizing.throughThickness = layers;
        part.sizing.cached = true;
        ++part.sizing.generation;
        return true;
    };

    if (part.source.empty()) {
        if (part.patch.uSize > 1.0e-6 && part.patch.vSize > 1.0e-6 && !fillFromPatch(true))
            return false;
        if (!part.sizing.cached) {
            part.sizing.cached = true;
            ++part.sizing.generation;
        }
        return true;
    }

    const SurfaceMesh& source = part.source;
    const bool closed = closedEnough(source) && std::abs(meshVolume(source)) > 1.0e-12;
    part.surface = {};
    part.volume = {};
    if (closed)
        return meshClosedBody(part, source, h);

    std::vector<double> thickness;
    std::vector<double> curvature;
    std::vector<ChVector3d> normals;
    std::vector<int> opposite;
    measureFaceFields(source, thickness, curvature, normals, opposite);
    (void)opposite;
    double thick = part.thickness;
    double thickSum = 0.0;
    int thickCount = 0;
    for (double t : thickness) {
        if (t < 1.0e20 && t > 1.0e-6) {
            thickSum += t;
            ++thickCount;
        }
    }
    if (thickCount > 0)
        thick = thickSum / static_cast<double>(thickCount);
    thick = std::max(thick, part.sizing.minElement);
    double span = std::max(part.patch.uSize, part.patch.vSize);
    if (!source.vertices.empty()) {
        ChVector3d lo = source.vertices.front();
        ChVector3d hi = lo;
        for (const ChVector3d& p : source.vertices) {
            lo = {std::min(lo.x(), p.x()), std::min(lo.y(), p.y()), std::min(lo.z(), p.z())};
            hi = {std::max(hi.x(), p.x()), std::max(hi.y(), p.y()), std::max(hi.z(), p.z())};
        }
        span = std::max({span, hi.x() - lo.x(), hi.y() - lo.y(), hi.z() - lo.z()});
    }
    int layers = extrusionLayers(thick, h);
    const bool slender = isSlenderWall(thick, span);
    // One linear tet through a thin open sheet has no plate bending. A welded
    // glass pane then rode into the floor as a rigid lid. A finite-strain
    // sheet cannot use that shell, so it falls through to the extrusion.
    if (slender && linearShellAllowed(part.material)) {
        part.surface = source;
        part.surface.thickness.resize(part.surface.triangles.size());
        for (std::size_t i = 0; i < part.surface.triangles.size(); ++i) {
            const double t = i < thickness.size() && thickness[i] < 1.0e20 ? thickness[i] : thick;
            part.surface.thickness[i] = std::max(t, part.sizing.minElement);
        }
        part.volume = {};
        part.thickness = thick;
        part.sizing.throughThickness = 0;
        part.sizing.cached = true;
        ++part.sizing.generation;
        return true;
    }
    const int cap = tetBudgetCap(part.sizing, part.material);
    const int faces = std::max(1, static_cast<int>(source.triangles.size()));
    while (3 * faces * layers > cap && layers > 1)
        --layers;
    if (!extrudeSheet(source, thick, part.volume, layers))
        return false;
    part.thickness = thick;
    part.sizing.throughThickness = layers;
    part.sizing.cached = true;
    ++part.sizing.generation;
    return true;
}

ChVector3d closestOnTriangle(const ChVector3d& p, const ChVector3d& a, const ChVector3d& b,
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
    if (vc <= 0.0 && d1 >= 0.0 && d3 <= 0.0) {
        const double v = d1 / (d1 - d3);
        return a + ab * v;
    }
    const ChVector3d cp = p - c;
    const double d5 = ab.Dot(cp);
    const double d6 = ac.Dot(cp);
    if (d6 >= 0.0 && d5 <= d6)
        return c;
    const double vb = d5 * d2 - d1 * d6;
    if (vb <= 0.0 && d2 >= 0.0 && d6 <= 0.0) {
        const double w = d2 / (d2 - d6);
        return a + ac * w;
    }
    const double va = d3 * d6 - d5 * d4;
    if (va <= 0.0 && (d4 - d3) >= 0.0 && (d5 - d6) >= 0.0) {
        const double w = (d4 - d3) / ((d4 - d3) + (d5 - d6));
        return b + (c - b) * w;
    }
    const double denom = 1.0 / (va + vb + vc);
    return a + ab * (vb * denom) + ac * (vc * denom);
}

double pointSegmentDistance(const ChVector3d& p, const ChVector3d& a, const ChVector3d& b) {
    const ChVector3d ab = b - a;
    const double len2 = ab.Length2();
    if (len2 < 1.0e-24)
        return (p - a).Length();
    const double t = std::clamp((p - a).Dot(ab) / len2, 0.0, 1.0);
    return (p - (a + ab * t)).Length();
}

void splitSkinEdge(SurfaceMesh& mesh, int ia, int ib, int ip) {
    std::vector<std::array<int, 3>> next;
    std::vector<double> nextT;
    const bool thick = !mesh.thickness.empty();
    next.reserve(mesh.triangles.size() + 2);
    for (std::size_t t = 0; t < mesh.triangles.size(); ++t) {
        const auto tri = mesh.triangles[t];
        const int hit = (tri[0] == ia || tri[0] == ib ? 1 : 0) + (tri[1] == ia || tri[1] == ib ? 1 : 0) +
                        (tri[2] == ia || tri[2] == ib ? 1 : 0);
        const double th = thick && t < mesh.thickness.size() ? mesh.thickness[t] : 0.0;
        if (hit < 2) {
            next.push_back(tri);
            if (thick)
                nextT.push_back(th);
            continue;
        }
        bool split = false;
        for (int k = 0; k < 3; ++k) {
            const int u = tri[k];
            const int v = tri[(k + 1) % 3];
            const int w = tri[(k + 2) % 3];
            if ((u != ia || v != ib) && (u != ib || v != ia))
                continue;
            next.push_back({u, ip, w});
            next.push_back({ip, v, w});
            if (thick) {
                nextT.push_back(th);
                nextT.push_back(th);
            }
            split = true;
            break;
        }
        if (!split) {
            next.push_back(tri);
            if (thick)
                nextT.push_back(th);
        }
    }
    mesh.triangles = std::move(next);
    if (thick)
        mesh.thickness = std::move(nextT);
}

void splitSkinFace(SurfaceMesh& mesh, std::size_t t, int ip) {
    if (t >= mesh.triangles.size())
        return;
    const auto tri = mesh.triangles[t];
    mesh.triangles[t] = {tri[0], tri[1], ip};
    mesh.triangles.push_back({tri[1], tri[2], ip});
    mesh.triangles.push_back({tri[2], tri[0], ip});
    if (!mesh.thickness.empty() && t < mesh.thickness.size()) {
        const double th = mesh.thickness[t];
        mesh.thickness.push_back(th);
        mesh.thickness.push_back(th);
    }
}

bool insertSkinPoint(SurfaceMesh& mesh, const ChVector3d& p, double maxDist, bool allowInterior) {
    if (mesh.triangles.empty() || mesh.vertices.empty())
        return false;

    constexpr double kSnap = 1.0e-4;
    int existing = -1;
    double bestVert = 1.0e30;
    for (std::size_t i = 0; i < mesh.vertices.size(); ++i) {
        const double d = (mesh.vertices[i] - p).Length();
        if (d < bestVert) {
            bestVert = d;
            existing = static_cast<int>(i);
        }
    }
    if (existing >= 0 && bestVert <= kSnap)
        return true;

    int bestTri = -1;
    double bestD = maxDist;
    for (std::size_t t = 0; t < mesh.triangles.size(); ++t) {
        const auto& tri = mesh.triangles[t];
        if (tri[0] < 0 || tri[1] < 0 || tri[2] < 0)
            continue;
        if (static_cast<std::size_t>(tri[0]) >= mesh.vertices.size() ||
            static_cast<std::size_t>(tri[1]) >= mesh.vertices.size() ||
            static_cast<std::size_t>(tri[2]) >= mesh.vertices.size())
            continue;
        const ChVector3d q = closestOnTriangle(p, mesh.vertices[static_cast<std::size_t>(tri[0])],
                                               mesh.vertices[static_cast<std::size_t>(tri[1])],
                                               mesh.vertices[static_cast<std::size_t>(tri[2])]);
        const double d = (q - p).Length();
        if (d <= bestD) {
            bestD = d;
            bestTri = static_cast<int>(t);
        }
    }
    if (bestTri < 0)
        return false;

    const auto tri = mesh.triangles[static_cast<std::size_t>(bestTri)];
    const ChVector3d va = mesh.vertices[static_cast<std::size_t>(tri[0])];
    const ChVector3d vb = mesh.vertices[static_cast<std::size_t>(tri[1])];
    const ChVector3d vc = mesh.vertices[static_cast<std::size_t>(tri[2])];
    const ChVector3d onSkin = closestOnTriangle(p, va, vb, vc);
    if ((onSkin - va).Length() <= kSnap || (onSkin - vb).Length() <= kSnap ||
        (onSkin - vc).Length() <= kSnap)
        return true;

    const double da = pointSegmentDistance(onSkin, vb, vc);
    const double db = pointSegmentDistance(onSkin, vc, va);
    const double dc = pointSegmentDistance(onSkin, va, vb);
    // Interior Steiner points on a large CAD face become sliver triangles and
    // tetgen then emits inverted tets. Only stamp the bead onto existing edges.
    constexpr double kEdge = 1.0e-3;
    int ea = -1;
    int eb = -1;
    ChVector3d onEdge = onSkin;
    if (dc <= kEdge) {
        ea = tri[0];
        eb = tri[1];
        onEdge = va + (vb - va) * (std::clamp((onSkin - va).Dot(vb - va) / std::max((vb - va).Length2(), 1.0e-24), 0.0, 1.0));
    } else if (da <= kEdge) {
        ea = tri[1];
        eb = tri[2];
        onEdge = vb + (vc - vb) * (std::clamp((onSkin - vb).Dot(vc - vb) / std::max((vc - vb).Length2(), 1.0e-24), 0.0, 1.0));
    } else if (db <= kEdge) {
        ea = tri[2];
        eb = tri[0];
        onEdge = vc + (va - vc) * (std::clamp((onSkin - vc).Dot(va - vc) / std::max((va - vc).Length2(), 1.0e-24), 0.0, 1.0));
    } else if (!allowInterior) {
        return false;
    } else {
        const int ip = static_cast<int>(mesh.vertices.size());
        mesh.vertices.push_back(onSkin);
        splitSkinFace(mesh, static_cast<std::size_t>(bestTri), ip);
        return true;
    }

    const int ip = static_cast<int>(mesh.vertices.size());
    mesh.vertices.push_back(onEdge);
    splitSkinEdge(mesh, ea, eb, ip);
    return true;
}

void seedWeldBeads(AssemblyDefinition& assembly) {
    auto hydratePatch = [](PartDefinition& part) {
        if (!part.source.triangles.empty())
            return;
        if (part.patch.uSize <= 1.0e-6 || part.patch.vSize <= 1.0e-6)
            return;
        const double thickness = std::max(part.thickness, part.sizing.minElement);
        const double span = std::max(part.patch.uSize, part.patch.vSize);
        // A slender patch stays a Reissner grid when the shell is legal, and
        // stays a patch (so meshing can refuse it) when the shell is not.
        // Baking either one into `source` would pick a different family.
        if (isSlenderWall(thickness, span))
            return;
        SurfaceMesh mid = sheetFromPatch(part.patch);
        chrono::ChVector3d u = part.patch.uAxis;
        chrono::ChVector3d v = part.patch.vAxis;
        if (u.Length() > 1.0e-12)
            u.Normalize();
        if (v.Length() > 1.0e-12)
            v.Normalize();
        chrono::ChVector3d n = u.Cross(v);
        if (n.Length() < 1.0e-12)
            n = {0, 1, 0};
        else
            n.Normalize();
        const chrono::ChVector3d shift = n * (-0.5 * thickness);
        for (auto& p : mid.vertices)
            p += shift;
        part.source = std::move(mid);
        part.sizing.cached = false;
        part.volume = {};
        part.surface = {};
    };

    for (const auto& weld : assembly.welds) {
        if (weld.bead.size() < 2 || weld.partA >= assembly.parts.size() ||
            weld.partB >= assembly.parts.size() || weld.partA == weld.partB)
            continue;
        PartDefinition& a = assembly.parts[weld.partA];
        PartDefinition& b = assembly.parts[weld.partB];
        hydratePatch(a);
        hydratePatch(b);
        if (a.source.triangles.empty() && b.source.triangles.empty())
            continue;

        double spacing = targetEdge(a);
        spacing = std::max(spacing, targetEdge(b));

        const double capture =
            std::max({weld.tolerance, 1.5 * std::min(std::max(a.thickness, 1.0e-4), std::max(b.thickness, 1.0e-4)),
                      2.0e-3});

        std::vector<ChVector3d> samples;
        const std::size_t n = weld.bead.size();
        for (std::size_t i = 0; i < n; ++i) {
            const ChVector3d pa = weld.bead[i];
            const ChVector3d pb = weld.bead[(i + 1) % n];
            const ChVector3d ab = pb - pa;
            const double len = ab.Length();
            if (len < 1.0e-9)
                continue;
            samples.push_back(pa);
            // Short thickness-direction bead edges are ~8 mm. Interior stamps
            // on those edges force a sliver row through an 8 mm plate. The
            // T-joint lives on the long CAD edges.
            if (len < 2.5 * std::max({a.thickness, b.thickness, 1.0e-3}))
                continue;
            const int steps = std::max(1, static_cast<int>(std::round(len / spacing)));
            for (int s = 1; s < steps; ++s) {
                const double t = static_cast<double>(s) / static_cast<double>(steps);
                samples.push_back(pa + ab * t);
            }
        }

        // Interior Steiner points belong on parts that will be solids. A
        // slender wall is a shell (or a refused patch); stamping the skin of
        // a 400 mm panel made tetgen blow the budget. throughThickness is not
        // the family switch.
        auto willBeSolid = [](const PartDefinition& part) {
            const double span = std::max(part.patch.uSize, part.patch.vSize);
            const double thickness = std::max(part.thickness, part.sizing.minElement);
            return part.patch.uSize > 1.0e-6 && part.patch.vSize > 1.0e-6 && !isSlenderWall(thickness, span);
        };
        const bool interiorA = willBeSolid(a);
        const bool interiorB = willBeSolid(b);
        for (const ChVector3d& sample : samples) {
            if (!a.source.triangles.empty())
                insertSkinPoint(a.source, sample, capture, interiorA);
            if (!b.source.triangles.empty())
                insertSkinPoint(b.source, sample, capture, interiorB);
        }
        auto markBead = [&](SurfaceMesh& mesh) {
            if (mesh.triangles.empty())
                return;
            if (mesh.markers.size() != mesh.triangles.size())
                mesh.markers.assign(mesh.triangles.size(), 1);
            for (std::size_t t = 0; t < mesh.triangles.size(); ++t) {
                const auto& tri = mesh.triangles[t];
                if (tri[0] < 0 || tri[1] < 0 || tri[2] < 0)
                    continue;
                const ChVector3d mid = (mesh.vertices[static_cast<std::size_t>(tri[0])] +
                                        mesh.vertices[static_cast<std::size_t>(tri[1])] +
                                        mesh.vertices[static_cast<std::size_t>(tri[2])]) /
                                       3.0;
                for (const ChVector3d& sample : samples) {
                    double nearest = (mid - sample).Length();
                    for (int index : tri) {
                        if (index < 0 || static_cast<std::size_t>(index) >= mesh.vertices.size())
                            continue;
                        nearest = std::min(nearest, (mesh.vertices[static_cast<std::size_t>(index)] - sample).Length());
                    }
                    if (nearest <= capture) {
                        mesh.markers[t] = 2;
                        break;
                    }
                }
            }
        };
        markBead(a.source);
        markBead(b.source);
    }
}

PartDefinition makeSteppedSolid(const ChVector3d& origin) {
    PartDefinition part;
    part.name = "Stepped solid";
    part.material = aluminum6061();
    part.initialTemperature = 293.15;
    const double boss = 0.12;
    const double web = 0.08;
    const double tall = 0.08;
    const double thin = 0.008;
    const double depth = 0.08;
    const auto at = [&](double x, double y, double z) { return origin + ChVector3d{x, y, z}; };
    part.source.vertices = {
        at(0.0, 0.0, 0.0),            at(boss + web, 0.0, 0.0), at(boss + web, 0.0, depth),
        at(0.0, 0.0, depth),          at(0.0, tall, 0.0),       at(boss, tall, 0.0),
        at(boss, tall, depth),        at(0.0, tall, depth),     at(boss, thin, 0.0),
        at(boss + web, thin, 0.0),    at(boss + web, thin, depth), at(boss, thin, depth)};
    const int tris[][3] = {
        {0, 3, 2}, {0, 2, 1}, {4, 5, 6}, {4, 6, 7}, {8, 9, 10}, {8, 10, 11},
        {0, 1, 9}, {0, 9, 8}, {0, 8, 5}, {0, 5, 4}, {3, 7, 6}, {3, 6, 11},
        {3, 11, 10}, {3, 10, 2}, {0, 4, 7}, {0, 7, 3}, {1, 2, 10}, {1, 10, 9},
        {5, 8, 11}, {5, 11, 6}};
    for (const auto& tri : tris)
        part.source.triangles.push_back({tri[0], tri[1], tri[2]});
    part.sizing.shellCutoff = 0.008;
    part.sizing.minElement = 0.002;
    part.sizing.maxElement = 0.025;
    part.sizing.tetBudget = 2500;
    return part;
}

PartDefinition makeSolidCube(double outerSize, const ChVector3d& center) {
    PartDefinition part;
    part.name = "Solid cube";
    part.material = aluminum6061();
    part.initialTemperature = 293.15;
    const double h = 0.5 * std::max(outerSize, 1.0e-4);
    const ChVector3d corner[8] = {
        center + ChVector3d{-h, -h, -h}, center + ChVector3d{h, -h, -h},
        center + ChVector3d{h, h, -h},   center + ChVector3d{-h, h, -h},
        center + ChVector3d{-h, -h, h},  center + ChVector3d{h, -h, h},
        center + ChVector3d{h, h, h},    center + ChVector3d{-h, h, h}};
    part.source.vertices.assign(corner, corner + 8);
    const int quads[6][4] = {
        {0, 3, 2, 1}, {4, 5, 6, 7}, {0, 1, 5, 4}, {3, 7, 6, 2}, {0, 4, 7, 3}, {1, 2, 6, 5}};
    for (const auto& face : quads) {
        part.source.faces.push_back(AuthoredFace{{face[0], face[1], face[2], face[3]}});
        part.source.triangles.push_back({face[0], face[1], face[2]});
        part.source.triangles.push_back({face[0], face[2], face[3]});
    }
    part.sizing.shellCutoff = 0.008;
    part.sizing.minElement = 0.002;
    part.sizing.maxElement = 0.025;
    part.sizing.throughThickness = 3;
    return part;
}

PartDefinition makeSolidSphere(double radius, const ChVector3d& center) {
    PartDefinition part;
    part.name = "Solid sphere";
    part.material = aluminum6061();
    part.initialTemperature = 293.15;
    const double r = std::max(radius, 1.0e-4);
    constexpr int kSegments = 16;
    constexpr int kRings = 8;
    constexpr double kPi = 3.14159265358979323846;
    for (int ring = 1; ring < kRings; ++ring) {
        const double phi = kPi * static_cast<double>(ring) / static_cast<double>(kRings);
        const double y = r * std::cos(phi);
        const double ringR = r * std::sin(phi);
        for (int i = 0; i < kSegments; ++i) {
            const double a = 2.0 * kPi * static_cast<double>(i) / static_cast<double>(kSegments);
            part.source.vertices.push_back(center + ChVector3d{ringR * std::cos(a), y, ringR * std::sin(a)});
        }
    }
    auto ringStart = [&](int ring) { return ring * kSegments; };
    for (int ring = 0; ring + 2 < kRings; ++ring) {
        const int a = ringStart(ring);
        const int b = ringStart(ring + 1);
        for (int i = 0; i < kSegments; ++i) {
            const int j = (i + 1) % kSegments;
            part.source.triangles.push_back({a + i, a + j, b + j});
            part.source.triangles.push_back({a + i, b + j, b + i});
        }
    }
    const int north = static_cast<int>(part.source.vertices.size());
    part.source.vertices.push_back(center + ChVector3d{0.0, r, 0.0});
    const int south = static_cast<int>(part.source.vertices.size());
    part.source.vertices.push_back(center + ChVector3d{0.0, -r, 0.0});
    for (int i = 0; i < kSegments; ++i) {
        const int j = (i + 1) % kSegments;
        part.source.triangles.push_back({ringStart(0) + i, north, ringStart(0) + j});
        const int last = ringStart(kRings - 2);
        part.source.triangles.push_back({last + i, last + j, south});
    }
    part.thickness = r;
    part.sizing.shellCutoff = 0.008;
    part.sizing.minElement = 0.002;
    part.sizing.maxElement = 0.025;
    part.sizing.throughThickness = 3;
    part.sizing.tetBudget = 800;
    return part;
}

}  // namespace rcs
