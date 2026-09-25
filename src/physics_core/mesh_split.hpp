#pragma once

#include "physics_core/part.hpp"

#include "chrono/core/ChVector3.h"

namespace rcs {

/// Measure local thickness, keep thin walls as shells, and fill thick regions
/// once. A later call does nothing until `sizing.cached` is cleared, so
/// simulate does not remesh.
bool ensurePartMesh(PartDefinition& part);

/// Solid aluminum cube of the same outer size as the shell box. The mesh is
/// cached on the part. `center` is the cube centre in assembly coordinates.
PartDefinition makeSolidCube(double outerSize = 0.40,
                             const chrono::ChVector3d& center = chrono::ChVector3d(0, 0, 0));

}  // namespace rcs
