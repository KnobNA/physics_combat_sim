#include "physics_core/mesh_split.hpp"

#define TETLIBRARY
#include "tetgen.h"

#include <algorithm>
#include <cmath>
#include <cstdio>
#include <map>
#include <unordered_map>
#include <utility>

namespace rcs {
namespace {

using chrono::ChVector3d;

double budgetEdge(double volume, const MeshSizing& sizing);

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

double faceArea(const SurfaceMesh& mesh, const std::array<int, 3>& face) {
    const ChVector3d& a = mesh.vertices[static_cast<std::size_t>(face[0])];
    const ChVector3d& b = mesh.vertices[static_cast<std::size_t>(face[1])];
    const ChVector3d& c = mesh.vertices[static_cast<std::size_t>(face[2])];
    return 0.5 * (b - a).Cross(c - a).Length();
}

bool sharesVertex(const std::array<int, 3>& a, const std::array<int, 3>& b) {
    for (int i : a)
        for (int j : b)
            if (i == j)
                return true;
    return false;
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
    return std::clamp(value, sizing.minElement, sizing.maxElement);
}

bool fillSolid(const SurfaceMesh& skin, const std::vector<double>& sizeAtVertex, VolumeMesh& volume,
               MeshSizing& sizing) {
    if (skin.triangles.size() < 4 || skin.vertices.size() < 4)
        return false;

    SurfaceMesh plc;
    std::vector<double> size;
    plc.vertices.reserve(skin.vertices.size());
    size.reserve(skin.vertices.size());
    plc.triangles.reserve(skin.triangles.size());
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
    }
    if (plc.triangles.size() < 4 || plc.vertices.size() < 4)
        return false;

    const bool flip = meshVolume(plc) < 0.0;
    double requested = budgetEdge(std::abs(meshVolume(plc)), sizing);
    // The caller's sizing field is a request, not a floor. Raising only the
    // ceiling cannot bring a mesh back under budget when the field itself is
    // what is refining it, which is the case for anything with a thin wall:
    // the wall asks for elements a third of its thickness no matter how coarse
    // the ceiling gets. Both have to be relaxed together.
    double relax = 1.0;
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
            const double local = i < static_cast<int>(size.size()) && size[static_cast<std::size_t>(i)] > 1.0e-12
                                     ? size[static_cast<std::size_t>(i)]
                                     : requested;
            // maxElement is the finest the initial mesh may use. Without this
            // floor a thin wall's thickness/N field (a few millimetres) wins
            // over the budget and tetgen emits tens of thousands of tets.
            const double minH = std::min(requested, std::max(sizing.minElement, sizing.maxElement));
            const double h = std::clamp(std::min(local * relax, requested), minH, requested);
            in.pointmtrlist[i] = h;
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
            in.facetmarkerlist[f] = 1;
        }

        const double maxVolume = std::max(1.0e-12, requested * requested * requested / 6.0);
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
            continue;
        }
        const int kBudget = std::max(4000, sizing.tetBudget * 4);
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
        for (int i = 0; i < out.numberofpoints; ++i) {
            volume.vertices.push_back({out.pointlist[3 * i], out.pointlist[3 * i + 1], out.pointlist[3 * i + 2]});
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
        }
        volume.elementSize = requested;
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
        if (shell.triangles.size() > 20000)
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
            const double target = clampSize(std::max(4.0 * thick, sizing.minElement), sizing);
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
            if (longest <= target * 1.4 || longest < sizing.minElement * 2.0) {
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

double budgetEdge(double volume, const MeshSizing& sizing) {
    const int budget = std::max(200, sizing.tetBudget);
    const double h = std::cbrt(std::max(volume, 1.0e-12) * 6.0 / static_cast<double>(budget));
    return std::max(h, sizing.minElement);
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
    layers = std::clamp(layers, 2, 16);
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

void appendVolume(VolumeMesh& dst, const VolumeMesh& src) {
    if (src.tets.empty())
        return;
    std::unordered_map<Key, int, KeyHash> weld;
    for (std::size_t i = 0; i < dst.vertices.size(); ++i)
        weld.emplace(quantize(dst.vertices[i], 1.0e-5), static_cast<int>(i));
    std::vector<int> map(src.vertices.size(), -1);
    for (std::size_t i = 0; i < src.vertices.size(); ++i) {
        int nearest = -1;
        double best = 0.012 * 0.012;
        for (std::size_t v = 0; v < dst.vertices.size(); ++v) {
            const double d = (dst.vertices[v] - src.vertices[i]).Length2();
            if (d < best) {
                best = d;
                nearest = static_cast<int>(v);
            }
        }
        map[i] = nearest >= 0 ? nearest : weldVertex(dst.vertices, weld, src.vertices[i], 1.0e-5);
    }
    for (const auto& tet : src.tets) {
        std::array<int, 4> mapped = {map[static_cast<std::size_t>(tet[0])], map[static_cast<std::size_t>(tet[1])],
                                     map[static_cast<std::size_t>(tet[2])], map[static_cast<std::size_t>(tet[3])]};
        addOrientedTet(dst, mapped[0], mapped[1], mapped[2], mapped[3]);
    }
    sealBoundary(dst);
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

}  // namespace

bool ensurePartMesh(PartDefinition& part) {
    if (part.sizing.cached)
        return true;
    if (part.source.empty()) {
        if (part.patch.uSize > 1.0e-6 && part.patch.vSize > 1.0e-6 &&
            part.sizing.throughThickness >= 2) {
            const double thickness = std::max(part.thickness, part.sizing.minElement);
            const int layers = std::clamp(part.sizing.throughThickness, 2, 16);
            if (!extrudeSheet(sheetFromPatch(part.patch), thickness, part.volume, layers))
                return false;
            part.surface = {};
            part.thickness = thickness;
        }
        part.sizing.cached = true;
        ++part.sizing.generation;
        return true;
    }

    const SurfaceMesh& source = part.source;
    const double inward = meshVolume(source) >= 0.0 ? -1.0 : 1.0;
    std::vector<double> thickness(source.triangles.size(), 1.0e30);
    std::vector<int> partner(source.triangles.size(), -1);
    std::vector<ChVector3d> normals(source.triangles.size());
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
        partner[i] = best.face;
    }

    const double cutoff = std::max(part.sizing.shellCutoff, part.sizing.minElement);
    std::vector<char> shell(source.triangles.size(), 0);
    int shellCount = 0;
    int solidCount = 0;
    for (std::size_t i = 0; i < source.triangles.size(); ++i) {
        if (thickness[i] <= cutoff) {
            shell[i] = 1;
            ++shellCount;
        } else if (thickness[i] < 1.0e20) {
            ++solidCount;
        }
    }
    // A watertight surface already says exactly what the solid is, so it gets
    // filled with tetrahedra whatever its proportions: the sizing field below
    // puts enough elements across a thin wall. Splitting a closed part into
    // "shell" and "solid" regions by wall thickness is what used to produce
    // two kinds of object out of one piece of CAD, and it fails outright on
    // shapes like a tube, whose wall reads as thin while its end faces read as
    // thick. Only a genuinely open surface, which encloses nothing, still has
    // to be given a thickness and extruded.
    const bool chunky = closedEnough(source) && std::abs(meshVolume(source)) > 1.0e-12;

    SurfaceMesh shells;
    part.surface = {};
    part.volume = {};

    auto emitShells = [&](bool onlyPaired) {
        std::vector<char> used(source.triangles.size(), 0);
        std::vector<char> onSolid(source.vertices.size(), 0);
        for (std::size_t i = 0; i < source.triangles.size(); ++i) {
            if (shell[i])
                continue;
            for (int index : source.triangles[i])
                if (index >= 0)
                    onSolid[static_cast<std::size_t>(index)] = 1;
        }
        std::unordered_map<Key, int, KeyHash> weld;
        for (std::size_t i = 0; i < source.triangles.size(); ++i) {
            if (used[i] || !shell[i])
                continue;
            if (onlyPaired && partner[i] < 0)
                continue;
            used[i] = 1;
            if (partner[i] >= 0 && static_cast<std::size_t>(partner[i]) < used.size())
                used[static_cast<std::size_t>(partner[i])] = 1;
            const double gap = std::min(thickness[i], cutoff);
            const ChVector3d shift = normals[i] * inward * (0.5 * gap);
            std::array<int, 3> tri{};
            for (int k = 0; k < 3; ++k) {
                const int index = source.triangles[i][static_cast<std::size_t>(k)];
                ChVector3d p = source.vertices[static_cast<std::size_t>(index)];
                if (!onSolid[static_cast<std::size_t>(index)])
                    p += shift;
                tri[k] = weldVertex(shells.vertices, weld, p, 1.0e-5);
            }
            if (tri[0] == tri[1] || tri[1] == tri[2] || tri[2] == tri[0])
                continue;
            shells.triangles.push_back(tri);
            shells.thickness.push_back(std::max(gap, part.sizing.minElement * 0.5));
        }
    };

    if (chunky) {
        const double coarse = budgetEdge(std::abs(meshVolume(source)), part.sizing);
        std::vector<double> size(source.vertices.size(), coarse);
        for (std::size_t v = 0; v < source.vertices.size(); ++v) {
            double h = coarse;
            double angle = 0.0;
            double edge = coarse;
            for (std::size_t i = 0; i < source.triangles.size(); ++i) {
                bool incident = false;
                for (int index : source.triangles[i])
                    incident = incident || index == static_cast<int>(v);
                if (!incident)
                    continue;
                const double finest = std::max(part.sizing.minElement, part.sizing.maxElement);
                const double local = thickness[i] < 1.0e20
                                         ? thickness[i] / std::max(1, part.sizing.throughThickness)
                                         : coarse;
                const double lo = std::min(finest, coarse);
                const double hi = std::max(finest, coarse);
                h = std::min(h, std::clamp(local, lo, hi));
                for (std::size_t j = i + 1; j < source.triangles.size(); ++j) {
                    if (!sharesVertex(source.triangles[i], source.triangles[j]))
                        continue;
                    const double turn = std::acos(std::clamp(normals[i].Dot(normals[j]), -1.0, 1.0));
                    if (turn > angle)
                        angle = turn;
                }
                for (int k = 0; k < 3; ++k) {
                    if (source.triangles[i][k] != static_cast<int>(v))
                        continue;
                    const int next = source.triangles[i][(k + 1) % 3];
                    edge = std::min(edge, (source.vertices[v] - source.vertices[static_cast<std::size_t>(next)]).Length());
                }
            }
            if (angle > 1.0e-3 && edge > 1.0e-6) {
                // Below a right angle this is curvature: the local radius is
                // edge/turn, and half of it resolves the curve. A rim or a box
                // corner is left alone, because refining every hard edge costs
                // elements everywhere and buys nothing. Past 115 degrees the
                // surface folds back on itself, as it does at a cone apex, and
                // there tetgen has no room to build a well-shaped element
                // unless the edges meeting at the point are shortened.
                // Curvature only: the local radius is edge/turn, and half of
                // it resolves the curve. Hard features such as a rim or a cone
                // apex are left coarse on purpose. Refining them multiplies
                // the element count without improving element shape, because
                // the corner stays just as sharp at any size.
                if (angle < 1.0)
                    h = std::min(h, std::clamp(edge / angle * 0.5, part.sizing.minElement, coarse));
            }
            size[v] = h;
        }
        if (!fillSolid(source, size, part.volume, part.sizing))
            return false;
        part.thickness = part.volume.elementSize;
    } else if (solidCount == 0) {
        if (shellCount > 0)
            emitShells(true);
        if (shells.empty())
            shells = source;
        double thick = part.thickness;
        if (!shells.thickness.empty()) {
            double sum = 0.0;
            for (double t : shells.thickness)
                sum += t;
            thick = sum / static_cast<double>(shells.thickness.size());
        }
        thick = std::max(thick, part.sizing.minElement);
        if (!extrudeSheet(shells, thick, part.volume, std::max(2, part.sizing.throughThickness)))
            return false;
        part.thickness = thick;
    } else {
        SurfaceMesh solidSkin;
        std::unordered_map<Key, int, KeyHash> weld;
        std::vector<int> remap(source.vertices.size(), -1);
        for (std::size_t i = 0; i < source.triangles.size(); ++i) {
            if (shell[i])
                continue;
            std::array<int, 3> tri{};
            for (int k = 0; k < 3; ++k) {
                const int index = source.triangles[i][static_cast<std::size_t>(k)];
                if (remap[static_cast<std::size_t>(index)] < 0)
                    remap[static_cast<std::size_t>(index)] =
                        weldVertex(solidSkin.vertices, weld, source.vertices[static_cast<std::size_t>(index)], 1.0e-6);
                tri[k] = remap[static_cast<std::size_t>(index)];
            }
            solidSkin.triangles.push_back(tri);
        }
        std::map<std::pair<int, int>, int> boundary;
        if (!closedEnough(solidSkin)) {
        for (const auto& face : solidSkin.triangles) {
            for (int k = 0; k < 3; ++k) {
                const int a = face[k];
                const int b = face[(k + 1) % 3];
                if (boundary[{b, a}])
                    boundary[{b, a}] = 0;
                else
                    boundary[{a, b}] = 1;
            }
        }
        }
        std::vector<std::pair<int, int>> open;
        for (const auto& edge : boundary)
            if (edge.second)
                open.push_back(edge.first);
        std::map<int, std::vector<int>> next;
        for (const auto& edge : open)
            next[edge.first].push_back(edge.second);
        std::map<std::pair<int, int>, char> consumed;
        for (const auto& edge : open) {
            if (consumed[edge])
                continue;
            std::vector<int> loop;
            int cursor = edge.first;
            int guard = 0;
            while (guard++ < 10000) {
                loop.push_back(cursor);
                const auto found = next.find(cursor);
                if (found == next.end())
                    break;
                int step = -1;
                for (int candidate : found->second) {
                    if (!consumed[{cursor, candidate}]) {
                        step = candidate;
                        break;
                    }
                }
                if (step < 0)
                    break;
                consumed[{cursor, step}] = 1;
                cursor = step;
                if (cursor == edge.first)
                    break;
            }
            if (loop.size() < 3)
                continue;
            for (std::size_t k = 1; k + 1 < loop.size(); ++k)
                solidSkin.triangles.push_back({loop[0], loop[k], loop[k + 1]});
        }

        std::vector<double> size(solidSkin.vertices.size(), part.sizing.maxElement);
        for (double& h : size)
            h = clampSize(part.sizing.maxElement, part.sizing);
        if (!fillSolid(solidSkin, size, part.volume, part.sizing)) {
            std::fprintf(stderr, "solid region %zu faces, volume %.6g, shells %d\n",
                         solidSkin.triangles.size(), meshVolume(solidSkin), shellCount);
            return false;
        }
        emitShells(true);
        if (!shells.empty()) {
            double thick = part.thickness;
            if (!shells.thickness.empty()) {
                double sum = 0.0;
                for (double t : shells.thickness)
                    sum += t;
                thick = sum / static_cast<double>(shells.thickness.size());
            }
            VolumeMesh fin;
            if (extrudeSheet(shells, std::max(thick, part.sizing.minElement), fin,
                             std::max(2, part.sizing.throughThickness)))
                appendVolume(part.volume, fin);
        }
    }

    part.sizing.cached = true;
    ++part.sizing.generation;
    return true;
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

}  // namespace rcs
