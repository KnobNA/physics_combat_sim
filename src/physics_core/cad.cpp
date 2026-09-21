#include "physics_core/cad.hpp"

#include <algorithm>

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
                     double size)
{
    PartDefinition part;
    part.name = name;
    part.material = material;
    part.thickness = thickness;
    part.patch.origin = origin;
    part.patch.uAxis = uAxis;
    part.patch.vAxis = vAxis;
    part.patch.uSize = size;
    part.patch.vSize = size;
    part.patch.divU = divisions;
    part.patch.divV = divisions;
    return part;
}

}  // namespace

AssemblyDefinition makeTestCube(const TestCubeOptions& options) {
    const double t = options.thickness;
    const int n = std::max(1, options.divisions);
    // Shell elements live on the panel mid-surface, so the six mid-surfaces
    // form a cube one thickness smaller than the outside dimension. Laying
    // them out this way makes the twelve seams geometrically coincident, which
    // is what the weld node-pairing looks for.
    const double size = options.outerSize - t;
    const double h = 0.5 * size;

    const ChVector3d x{1, 0, 0};
    const ChVector3d y{0, 1, 0};
    const ChVector3d z{0, 0, 1};

    AssemblyDefinition assembly;
    assembly.name = "Welded test box";

    // Axes are ordered so that u x v points out of the box.
    assembly.parts.push_back(panel("Bottom", options.material, t, n, {-h, -h, -h}, x, z, size));
    assembly.parts.push_back(panel("Top", options.material, t, n, {-h, h, -h}, z, x, size));
    assembly.parts.push_back(panel("Front", options.material, t, n, {-h, -h, h}, x, y, size));
    assembly.parts.push_back(panel("Back", options.material, t, n, {-h, -h, -h}, y, x, size));
    assembly.parts.push_back(panel("Right", options.material, t, n, {h, -h, -h}, y, z, size));
    assembly.parts.push_back(panel("Left", options.material, t, n, {-h, -h, -h}, z, y, size));

    enum { kBottom = 0, kTop, kFront, kBack, kRight, kLeft };
    const std::pair<int, int> seams[12] = {
        {kBottom, kFront}, {kBottom, kBack}, {kBottom, kRight}, {kBottom, kLeft},
        {kTop, kFront},    {kTop, kBack},    {kTop, kRight},    {kTop, kLeft},
        {kFront, kRight},  {kFront, kLeft},  {kBack, kRight},   {kBack, kLeft},
    };

    for (const auto& seam : seams) {
        WeldDefinition weld;
        weld.partA = static_cast<std::size_t>(seam.first);
        weld.partB = static_cast<std::size_t>(seam.second);
        weld.label = assembly.parts[weld.partA].name + " - " + assembly.parts[weld.partB].name;
        weld.tolerance = 0.25 * t;
        weld.throat = t;
        assembly.welds.push_back(weld);
    }

    return assembly;
}

bool importAssembly(const std::string& path, AssemblyDefinition& out, std::string& error) {
    (void)path;
    (void)out;
    error =
        "CAD import is not wired up yet. Each body in the file will become one Part "
        "and each fastener feature one ScrewDefinition; see the deformation plan.";
    return false;
}

}  // namespace rcs
