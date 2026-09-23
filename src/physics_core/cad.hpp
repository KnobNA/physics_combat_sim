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

/// Load an OBJ or STL surface. One part per object. Coincident circular holes
/// become screws. Watertight chunky solids are skipped and named in `error`.
/// Welds are not inferred.
bool importAssembly(const std::string& path, AssemblyDefinition& out, std::string& error);

}  // namespace rcs
