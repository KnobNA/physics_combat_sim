#include "physics_core/cad.hpp"
#include "physics_core/mesh_split.hpp"

#include <algorithm>
#include <cctype>
#include <cmath>
#include <cstdint>
#include <cstdlib>
#include <fstream>
#include <map>
#include <sstream>

namespace rcs {
namespace {

using chrono::ChVector3d;

PartDefinition panel(const std::string& name,
                     const Material& material,
                     double thickness,
                     int divisions,
                     const ChVector3d& origin,
                     const ChVector3d& uAxis,
                     const ChVector3d& vAxis,
                     double uSize,
                     double vSize)
{
    PartDefinition part;
    part.name = name;
    part.material = material;
    part.thickness = thickness;
    part.patch.origin = origin;
    part.patch.uAxis = uAxis;
    part.patch.vAxis = vAxis;
    part.patch.uSize = uSize;
    part.patch.vSize = vSize;
    part.patch.divU = std::max(1, divisions);
    part.patch.divV = std::max(1, divisions);

    ChVector3d u = uAxis;
    ChVector3d v = vAxis;
    if (u.Length() > 1.0e-12)
        u.Normalize();
    if (v.Length() > 1.0e-12)
        v.Normalize();
    ChVector3d n = u.Cross(v);
    if (n.Length() < 1.0e-12)
        n = {0, 1, 0};
    else
        n.Normalize();
    const double h = 0.5 * std::max(thickness, 1.0e-6);
    auto corner = [&](double su, double sv, double side) {
        return origin + u * (su * uSize) + v * (sv * vSize) + n * (side * h);
    };
    // Keep the 8-vertex solid. An inset 8 mm frame on the large skins is a
    // sliver PLC: tetgen fills it with needle tets, weld pairing then puts
    // two coincident rows 8 mm apart, and HHT Newton stalls.
    part.source.vertices = {corner(0.0, 0.0, -1.0), corner(1.0, 0.0, -1.0), corner(1.0, 1.0, -1.0),
                            corner(0.0, 1.0, -1.0), corner(0.0, 0.0, 1.0),  corner(1.0, 0.0, 1.0),
                            corner(1.0, 1.0, 1.0),  corner(0.0, 1.0, 1.0)};
    const int quads[6][4] = {
        {0, 3, 2, 1}, {4, 5, 6, 7}, {0, 1, 5, 4}, {1, 2, 6, 5}, {2, 3, 7, 6}, {3, 0, 4, 7}};
    for (const auto& face : quads) {
        part.source.faces.push_back(AuthoredFace{{face[0], face[1], face[2], face[3]}});
        part.source.triangles.push_back({face[0], face[1], face[2]});
        part.source.triangles.push_back({face[0], face[2], face[3]});
    }
    part.sizing.shellCutoff = 0.008;
    part.sizing.minElement = 0.008;
    part.sizing.maxElement = 0.025;
    part.sizing.throughThickness = 2;
    part.sizing.tetBudget = 800;
    return part;
}

struct PanelBox {
    ChVector3d lo;
    ChVector3d hi;
};

PanelBox panelBox(const PartDefinition& part) {
    PanelBox box{part.source.vertices.front(), part.source.vertices.front()};
    for (const auto& p : part.source.vertices) {
        box.lo = {std::min(box.lo.x(), p.x()), std::min(box.lo.y(), p.y()), std::min(box.lo.z(), p.z())};
        box.hi = {std::max(box.hi.x(), p.x()), std::max(box.hi.y(), p.y()), std::max(box.hi.z(), p.z())};
    }
    return box;
}

/// Contact patch between two kissing panels: the rectangle where their AABBs
/// meet. Hand-authored Back-Right/Left beads used to sit on the outer back
/// face, 8 mm from the real seam, which stretched weld links through air.
std::vector<ChVector3d> kissBead(const PanelBox& a, const PanelBox& b) {
    const ChVector3d lo{std::max(a.lo.x(), b.lo.x()), std::max(a.lo.y(), b.lo.y()),
                        std::max(a.lo.z(), b.lo.z())};
    const ChVector3d hi{std::min(a.hi.x(), b.hi.x()), std::min(a.hi.y(), b.hi.y()),
                        std::min(a.hi.z(), b.hi.z())};
    const double span[3] = {hi.x() - lo.x(), hi.y() - lo.y(), hi.z() - lo.z()};
    int axis = 0;
    for (int i = 1; i < 3; ++i) {
        if (span[i] < span[axis])
            axis = i;
    }
    if (span[axis] < -1.0e-6)
        return {};
    const int u = (axis + 1) % 3;
    const int v = (axis + 2) % 3;
    const double mid = 0.5 * ((axis == 0 ? lo.x() : axis == 1 ? lo.y() : lo.z()) +
                              (axis == 0 ? hi.x() : axis == 1 ? hi.y() : hi.z()));
    const double ulo = u == 0 ? lo.x() : u == 1 ? lo.y() : lo.z();
    const double uhi = u == 0 ? hi.x() : u == 1 ? hi.y() : hi.z();
    const double vlo = v == 0 ? lo.x() : v == 1 ? lo.y() : lo.z();
    const double vhi = v == 0 ? hi.x() : v == 1 ? hi.y() : hi.z();
    if (uhi - ulo < 1.0e-9 || vhi - vlo < 1.0e-9)
        return {};
    auto at = [&](double uu, double vv) {
        ChVector3d p(0, 0, 0);
        p[axis] = mid;
        p[u] = uu;
        p[v] = vv;
        return p;
    };
    return {at(ulo, vlo), at(uhi, vlo), at(uhi, vhi), at(ulo, vhi)};
}

std::string extensionOf(const std::string& path) {
    const auto dot = path.find_last_of('.');
    if (dot == std::string::npos)
        return {};
    std::string ext = path.substr(dot + 1);
    for (char& c : ext)
        c = static_cast<char>(std::tolower(static_cast<unsigned char>(c)));
    return ext;
}

std::string fileStem(const std::string& path) {
    auto slash = path.find_last_of("/\\");
    std::string name = slash == std::string::npos ? path : path.substr(slash + 1);
    auto dot = name.find_last_of('.');
    if (dot != std::string::npos)
        name.resize(dot);
    return name.empty() ? "Imported" : name;
}

int parseIndex(const std::string& token, int count) {
    const auto slash = token.find('/');
    const int raw = std::atoi(token.substr(0, slash).c_str());
    if (raw > 0)
        return raw - 1;
    if (raw < 0)
        return count + raw;
    return -1;
}

void addFace(SurfaceMesh& mesh, const std::vector<int>& corners) {
    if (corners.size() < 3)
        return;
    mesh.faces.push_back(AuthoredFace{corners});
    for (std::size_t i = 1; i + 1 < corners.size(); ++i)
        mesh.triangles.push_back({corners[0], corners[i], corners[i + 1]});
}

bool loadObj(const std::string& path, std::vector<SurfaceMesh>& meshes, std::vector<std::string>& names,
             std::string& error)
{
    std::ifstream in(path);
    if (!in) {
        error = "Could not open " + path;
        return false;
    }
    std::vector<ChVector3d> vertices;
    SurfaceMesh* current = nullptr;
    std::string line;
    while (std::getline(in, line)) {
        if (line.size() >= 2 && line[0] == 'v' && line[1] == ' ') {
            std::istringstream stream(line.substr(2));
            double x = 0, y = 0, z = 0;
            stream >> x >> y >> z;
            vertices.push_back({x, y, z});
            continue;
        }
        if (line.size() >= 2 && (line[0] == 'o' || line[0] == 'g') && line[1] == ' ') {
            names.push_back(line.substr(2));
            meshes.push_back({});
            current = &meshes.back();
            continue;
        }
        if (line.size() >= 2 && line[0] == 'f' && line[1] == ' ') {
            if (!current) {
                names.push_back(fileStem(path));
                meshes.push_back({});
                current = &meshes.back();
            }
            std::istringstream stream(line.substr(2));
            std::string token;
            std::vector<int> corners;
            while (stream >> token) {
                const int index = parseIndex(token, static_cast<int>(vertices.size()));
                if (index >= 0)
                    corners.push_back(index);
            }
            addFace(*current, corners);
        }
    }
    for (auto& mesh : meshes) {
        mesh.vertices = vertices;
    }
    if (meshes.empty()) {
        error = "No faces in " + path;
        return false;
    }
    return true;
}

bool loadStl(const std::string& path, SurfaceMesh& mesh, std::string& error) {
    std::ifstream in(path, std::ios::binary);
    if (!in) {
        error = "Could not open " + path;
        return false;
    }
    in.seekg(0, std::ios::end);
    const auto size = in.tellg();
    in.seekg(0);
    if (size < 84) {
        error = "STL file is too small.";
        return false;
    }
    char header[80] = {};
    in.read(header, 80);
    const bool ascii = std::string(header, header + 5) == "solid" &&
                       std::string(header, header + 80).find('\0') == std::string::npos;
    in.seekg(0);
    if (ascii) {
        std::string line;
        ChVector3d corners[3];
        int corner = 0;
        bool inFacet = false;
        while (std::getline(in, line)) {
            std::istringstream stream(line);
            std::string word;
            stream >> word;
            if (word == "facet")
                inFacet = true;
            else if (word == "vertex" && inFacet && corner < 3) {
                double x = 0, y = 0, z = 0;
                stream >> x >> y >> z;
                corners[corner++] = {x, y, z};
            } else if (word == "endfacet" && corner == 3) {
                const int base = static_cast<int>(mesh.vertices.size());
                mesh.vertices.insert(mesh.vertices.end(), corners, corners + 3);
                mesh.triangles.push_back({base, base + 1, base + 2});
                corner = 0;
                inFacet = false;
            }
        }
    } else {
        in.seekg(80);
        std::uint32_t count = 0;
        in.read(reinterpret_cast<char*>(&count), 4);
        mesh.vertices.reserve(static_cast<std::size_t>(count) * 3);
        for (std::uint32_t i = 0; i < count; ++i) {
            float values[12] = {};
            in.read(reinterpret_cast<char*>(values), 48);
            char pad[2];
            in.read(pad, 2);
            if (!in)
                break;
            const int base = static_cast<int>(mesh.vertices.size());
            for (int v = 0; v < 3; ++v)
                mesh.vertices.push_back({values[3 + v * 3], values[4 + v * 3], values[5 + v * 3]});
            mesh.triangles.push_back({base, base + 1, base + 2});
        }
    }
    if (mesh.triangles.empty()) {
        error = "No triangles in " + path;
        return false;
    }
    return true;
}

struct Hole {
    std::size_t part = 0;
    ChVector3d center{0, 0, 0};
    double radius = 0.0;
};

void findHoles(const SurfaceMesh& mesh, std::size_t part, std::vector<Hole>& holes) {
    std::map<std::pair<int, int>, int> edges;
    for (const auto& face : mesh.triangles) {
        const int id[3] = {face[0], face[1], face[2]};
        for (int i = 0; i < 3; ++i) {
            const int a = id[i];
            const int b = id[(i + 1) % 3];
            edges[{a, b}]++;
        }
    }
    std::map<int, std::vector<int>> boundary;
    for (const auto& edge : edges) {
        const auto reverse = edges.find({edge.first.second, edge.first.first});
        if (reverse != edges.end())
            continue;
        boundary[edge.first.first].push_back(edge.first.second);
    }
    std::map<int, bool> used;
    for (const auto& start : boundary) {
        if (used[start.first])
            continue;
        std::vector<int> loop;
        int cursor = start.first;
        for (int guard = 0; guard < 64 && boundary.count(cursor) && !used[cursor]; ++guard) {
            used[cursor] = true;
            loop.push_back(cursor);
            cursor = boundary[cursor].front();
            if (cursor == start.first)
                break;
        }
        if (loop.size() < 6 || loop.size() > 32)
            continue;
        ChVector3d center{0, 0, 0};
        for (int index : loop)
            center += mesh.vertices[static_cast<std::size_t>(index)];
        center /= static_cast<double>(loop.size());
        double radius = 0.0;
        for (int index : loop)
            radius += (mesh.vertices[static_cast<std::size_t>(index)] - center).Length();
        radius /= static_cast<double>(loop.size());
        if (radius < 0.001 || radius > 0.02)
            continue;
        double variance = 0.0;
        for (int index : loop) {
            const double delta = (mesh.vertices[static_cast<std::size_t>(index)] - center).Length() - radius;
            variance += delta * delta;
        }
        if (std::sqrt(variance / static_cast<double>(loop.size())) > 0.2 * radius)
            continue;
        holes.push_back({part, center, radius});
    }
}

bool isChunkySolid(const SurfaceMesh& mesh) {
    std::map<std::pair<int, int>, int> edges;
    double area = 0.0;
    double volume = 0.0;
    for (const auto& face : mesh.triangles) {
        const ChVector3d& a = mesh.vertices[static_cast<std::size_t>(face[0])];
        const ChVector3d& b = mesh.vertices[static_cast<std::size_t>(face[1])];
        const ChVector3d& c = mesh.vertices[static_cast<std::size_t>(face[2])];
        area += 0.5 * (b - a).Cross(c - a).Length();
        volume += a.Dot(b.Cross(c)) / 6.0;
        for (int i = 0; i < 3; ++i)
            edges[{face[i], face[(i + 1) % 3]}]++;
    }
    bool closed = true;
    for (const auto& edge : edges) {
        if (edges.find({edge.first.second, edge.first.first}) == edges.end()) {
            closed = false;
            break;
        }
    }
    if (!closed || area < 1.0e-8)
        return false;
    const double thickness = std::abs(volume) / (0.5 * area);
    return thickness > 0.02;
}

}  // namespace

AssemblyDefinition makeTestCube(const TestCubeOptions& options) {
    const double t = options.thickness;
    const int n = std::max(1, options.divisions);
    const double outer = options.outerSize;
    const double ho = 0.5 * outer;
    const double hi = ho - t;
    const double inner = outer - 2.0 * t;

    const ChVector3d x{1, 0, 0};
    const ChVector3d y{0, 1, 0};
    const ChVector3d z{0, 0, 1};

    AssemblyDefinition assembly;
    assembly.name = "Welded test box";
    assembly.groups.push_back(GroupDefinition{"Hollow cube", -1});

    // Thickness-edge faces kiss at the inner envelope so tet skins share a
    // contact strip instead of overlapping through each other. Back uses the
    // same (x, y) UV as Front so the two walls are mirrors, not a swapped
    // parametrisation that hid a weld-bead copy-paste on the outer back face.
    assembly.parts.push_back(panel("Bottom", options.material, t, n, {-ho, -ho + 0.5 * t, -ho}, x, z, outer, outer));
    assembly.parts.push_back(panel("Top", options.material, t, n, {-ho, ho - 0.5 * t, -ho}, z, x, outer, outer));
    assembly.parts.push_back(panel("Front", options.material, t, n, {-ho, -hi, ho - 0.5 * t}, x, y, outer, inner));
    assembly.parts.push_back(panel("Back", options.material, t, n, {-ho, -hi, -ho + 0.5 * t}, x, y, outer, inner));
    assembly.parts.push_back(panel("Right", options.material, t, n, {ho - 0.5 * t, -hi, -hi}, y, z, inner, inner));
    assembly.parts.push_back(panel("Left", options.material, t, n, {-ho + 0.5 * t, -hi, -hi}, z, y, inner, inner));
    for (auto& part : assembly.parts)
        part.group = 0;

    enum { kBottom = 0, kTop, kFront, kBack, kRight, kLeft };
    const std::pair<int, int> seams[12] = {
        {kBottom, kFront}, {kBottom, kBack}, {kBottom, kRight}, {kBottom, kLeft},
        {kTop, kFront},    {kTop, kBack},    {kTop, kRight},    {kTop, kLeft},
        {kFront, kRight},  {kFront, kLeft},  {kBack, kRight},   {kBack, kLeft},
    };

    std::vector<PanelBox> boxes;
    boxes.reserve(assembly.parts.size());
    for (const auto& part : assembly.parts)
        boxes.push_back(panelBox(part));

    for (int i = 0; i < 12; ++i) {
        WeldDefinition weld;
        weld.partA = static_cast<std::size_t>(seams[i].first);
        weld.partB = static_cast<std::size_t>(seams[i].second);
        weld.label = assembly.parts[weld.partA].name + " - " + assembly.parts[weld.partB].name;
        weld.tolerance = 0.5 * t;
        weld.throat = t;
        weld.bead = kissBead(boxes[weld.partA], boxes[weld.partB]);
        assembly.welds.push_back(weld);
    }

    return assembly;
}

bool importAssembly(const std::string& path, AssemblyDefinition& out, std::string& error) {
    error.clear();
    const std::string ext = extensionOf(path);
    std::vector<SurfaceMesh> meshes;
    std::vector<std::string> names;
    if (ext == "obj") {
        if (!loadObj(path, meshes, names, error))
            return false;
    } else if (ext == "stl") {
        SurfaceMesh mesh;
        if (!loadStl(path, mesh, error))
            return false;
        meshes.push_back(std::move(mesh));
        names.push_back(fileStem(path));
    } else {
        error = "Import reads OBJ and STL surfaces. STEP solids need a mesher this build does not include.";
        return false;
    }

    AssemblyDefinition assembly;
    assembly.name = fileStem(path);
    std::string skipped;
    std::vector<SurfaceMesh> kept;
    std::vector<std::size_t> keptIndex;
    for (std::size_t i = 0; i < meshes.size(); ++i) {
        if (meshes[i].triangles.empty())
            continue;
        PartDefinition part;
        part.name = names[i].empty() ? ("Part " + std::to_string(assembly.parts.size() + 1)) : names[i];
        part.material = aluminum6061();
        part.thickness = 0.008;
        part.source = std::move(meshes[i]);
        keptIndex.push_back(assembly.parts.size());
        kept.push_back(part.source);
        assembly.parts.push_back(std::move(part));
    }
    if (assembly.parts.empty()) {
        error = skipped.empty() ? "The file had no surfaces."
                                : "Could not mesh: " + skipped;
        return false;
    }

    std::vector<Hole> holes;
    for (std::size_t i = 0; i < kept.size(); ++i)
        findHoles(kept[i], keptIndex[i], holes);
    for (std::size_t i = 0; i < holes.size(); ++i) {
        for (std::size_t j = i + 1; j < holes.size(); ++j) {
            if (holes[i].part == holes[j].part)
                continue;
            if ((holes[i].center - holes[j].center).Length() > 0.008)
                continue;
            ScrewDefinition screw;
            screw.partA = holes[i].part;
            screw.partB = holes[j].part;
            screw.position = 0.5 * (holes[i].center + holes[j].center);
            screw.shankDiameter = std::min(holes[i].radius, holes[j].radius);
            screw.clampRadius = std::max(0.012, screw.shankDiameter * 2.0);
            screw.label = "Screw " + std::to_string(assembly.screws.size() + 1);
            assembly.screws.push_back(screw);
        }
    }

    if (!skipped.empty())
        error = "Imported what could be meshed. Failed: " + skipped;
    out = std::move(assembly);
    return true;
}

}  // namespace rcs
