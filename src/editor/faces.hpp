#pragma once

#include "physics_core/assembly.hpp"
#include "physics_core/part.hpp"
#include "physics_core/scenario.hpp"

#include "chrono/core/ChQuaternion.h"
#include "chrono/core/ChVector3.h"

#include <cstddef>
#include <optional>
#include <string>
#include <vector>

namespace rcs {

struct FaceId {
    std::size_t part = 0;
    std::size_t face = 0;

    bool operator==(const FaceId& other) const {
        return part == other.part && face == other.face;
    }
};

struct SpliceHit {
    int edge = -1;
    double along = 0.0;
    chrono::ChVector3d point{0, 0, 0};
};

struct SplicePreview {
    bool valid = false;
    FaceId face;
    chrono::ChVector3d chordA{0, 0, 0};
    chrono::ChVector3d chordB{0, 0, 0};
    double chordLength = 0.0;
    SpliceHit hits[2];
    double interiorDeg[2] = {90.0, 90.0};
    bool square[2] = {false, false};
    std::vector<chrono::ChVector3d> loop;
    std::vector<double> edgeLength;
    std::vector<double> edgeSplitA;
    std::vector<double> edgeSplitB;
};

/// Fill `source.faces` from a plate patch, OBJ loops, or coplanar triangle
/// groups. Safe to call often; does nothing when faces already exist.
void ensureAuthoredFaces(PartDefinition& part);

/// Rebuild `source.triangles` as a fan of each authored face.
void tessellateAuthoredFaces(SurfaceMesh& mesh);

chrono::ChVector3d partCentroid(const PartDefinition& part);
chrono::ChVector3d partEuler(const PartDefinition& part);
void applyPartTranslation(PartDefinition& part, const chrono::ChVector3d& delta);
void applyPartRotation(PartDefinition& part, const chrono::ChQuaterniond& rotation,
                       const chrono::ChVector3d& pivot);

std::optional<FaceId> pickAuthoredFace(const AssemblyDefinition& definition,
                                       const std::vector<std::size_t>& parts,
                                       const chrono::ChVector3d& origin,
                                       const chrono::ChVector3d& direction,
                                       const chrono::ChQuaterniond& worldRotation,
                                       const chrono::ChVector3d& worldOrigin,
                                       chrono::ChVector3d* hitPoint = nullptr);

SplicePreview makeSplicePreview(const PartDefinition& part, std::size_t faceIndex, double angle,
                                std::optional<chrono::ChVector3d> through = {});
/// Split the face along the preview chord. Returns the two new face indices.
bool spliceFace(PartDefinition& part, const SplicePreview& preview, std::size_t& faceA,
                std::size_t& faceB);

chrono::ChVector3d toWorld(const chrono::ChVector3d& authored, const chrono::ChQuaterniond& rotation,
                           const chrono::ChVector3d& origin);
void faceLoopWorld(const PartDefinition& part, std::size_t faceIndex,
                   const chrono::ChQuaterniond& rotation, const chrono::ChVector3d& origin,
                   std::vector<chrono::ChVector3d>& out);

std::string formatMillimetres(double metres);
std::string formatDegrees(double degrees);

bool partIsDescendant(const AssemblyDefinition& definition, std::size_t candidate, std::size_t ancestor);
bool groupIsDescendant(const AssemblyDefinition& definition, int candidate, int ancestor);
void renormalizeSiblingOrder(AssemblyDefinition& definition);

std::vector<std::size_t> collectTransformTargets(const AssemblyDefinition& definition,
                                                 const std::vector<std::size_t>& selectedParts,
                                                 const std::optional<std::size_t>& selectedGroup);

void deleteSelectedParts(AssemblyDefinition& definition, std::vector<std::size_t>& selected,
                         std::vector<ImpactLoad>& loads);
void deleteSelectedGroup(AssemblyDefinition& definition, std::size_t groupIndex,
                         std::vector<std::size_t>& selectedParts);

}  // namespace rcs
