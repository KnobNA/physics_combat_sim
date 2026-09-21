#pragma once

#include "physics_core/assembly.hpp"
#include "physics_core/material.hpp"

#include <string>

namespace rcs {

/// Parameters for the stand-in assembly used until real CAD import lands:
/// six finite-thickness panels welded along the twelve shared edges.
struct TestCubeOptions {
    double outerSize = 0.40;   ///< m, outside dimension of the box
    double thickness = 0.008;  ///< m, panel thickness
    int divisions = 4;         ///< shell elements per panel edge
    Material material = aluminum6061();
};

/// Build the six-panel welded test box. Every panel is an independent part and
/// every edge is an explicit manual weld, exactly as a CAD import would produce.
AssemblyDefinition makeTestCube(const TestCubeOptions& options = {});

/// Load an assembly from a CAD file, one part per solid body in the file.
/// Not implemented yet; returns false and fills `error`.
bool importAssembly(const std::string& path, AssemblyDefinition& out, std::string& error);

}  // namespace rcs
