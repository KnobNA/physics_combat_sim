#pragma once

#include "physics_core/part.hpp"

#include "chrono/core/ChVector3.h"

namespace rcs {

struct AssemblyDefinition;

/// Measure local thickness, keep thin walls as shells, and fill thick regions
/// once. A later call does nothing until `sizing.cached` is cleared, so
/// simulate does not remesh.
bool ensurePartMesh(PartDefinition& part);

/// Scale tetBudget by scenario resolution and set refineLimit/maxElement to
/// the uniform target edge for this part.
void applyResolutionSizing(PartDefinition& part, double resolution);

/// Stamp each weld bead onto both parts' authored skins so tetgen keeps the
/// same nodes along the seam. Independent fills otherwise never share vertices,
/// and a T-joint then only pins at the CAD corners.
void seedWeldBeads(AssemblyDefinition& assembly);

/// Solid aluminum cube of the same outer size as the shell box. The mesh is
/// cached on the part. `center` is the cube centre in assembly coordinates.
PartDefinition makeSolidCube(double outerSize = 0.40,
                             const chrono::ChVector3d& center = chrono::ChVector3d(0, 0, 0));

/// Solid aluminum sphere. The mesh is cached on the part. `center` is the
/// sphere centre in assembly coordinates.
PartDefinition makeSolidSphere(double radius = 0.10,
                               const chrono::ChVector3d& center = chrono::ChVector3d(0, 0, 0));

/// Closed L-step: 12 cm × 8 cm × 8 cm boss plus an 8 cm × 8 mm × 8 cm web.
PartDefinition makeSteppedSolid(const chrono::ChVector3d& origin = chrono::ChVector3d(0, 0, 0));

}  // namespace rcs
