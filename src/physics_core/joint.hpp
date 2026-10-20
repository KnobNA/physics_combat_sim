#pragma once

#include "physics_core/interface.hpp"

#include "chrono/core/ChVector3.h"
#include "chrono/fea/ChLinkNodeFrame.h"
#include "chrono/fea/ChLinkNodeNode.h"
#include "chrono/physics/ChLinkMate.h"

#include <cstddef>
#include <memory>
#include <string>
#include <vector>

namespace rcs {

/// A manual weld. The user states which two parts are joined; nodes along
/// the (optional) bead become one weld link each, so a seam can fail
/// progressively instead of all at once. Break force is
/// weldShearStrength times tributary area: spacing*throat on an edge, or
/// a share of the enclosed face on a closed bead. Uses the weaker parent.
struct WeldDefinition {
    std::string label = "Weld";
    std::size_t partA = 0;
    std::size_t partB = 1;
    double tolerance = 2.0e-3;   ///< node pairing distance, m
    double throat = 0.0;         ///< bead throat, m; 0 means min part thickness
    double strengthScale = 1.0;  ///< derating factor for a poor bead
    /// Peel strength, Pa. Zero uses the shear allowable as the peel strength.
    double peelStrength = 0.0;
    /// Explicit seam rating. Zero fields are derived from the two parents.
    InterfaceDefinition interface;
    bool useRegion = false;
    chrono::ChVector3d regionMin{0, 0, 0};
    chrono::ChVector3d regionMax{0, 0, 0};
    /// Closed perimeter in authored metres. Empty keeps the old behaviour:
    /// every coincident boundary pair inside the optional AABB. A nonempty
    /// bead keeps only pairs whose midpoint sits on this loop, so a flush
    /// face weld does not glue the whole skin.
    std::vector<chrono::ChVector3d> bead;
};

/// A manual fastener at a point. A bolt pattern is several of these; CAD
/// import can emit one per hole it finds.
struct ScrewDefinition {
    std::string label = "Screw";
    std::size_t partA = 0;
    std::size_t partB = 1;
    chrono::ChVector3d position{0, 0, 0};
    double shankDiameter = 0.005;  ///< m
    double clampRadius = 0.012;    ///< how far from `position` nodes may be grabbed
    double proofStress = 6.0e8;    ///< Pa, class 8.8 fastener
    /// A single screw through a clearance hole is a pin; a torqued flange bolt
    /// with a large washer resists rotation too.
    bool locksRotation = false;
};

enum class JointKind { Weld, Screw };

/// One constraint inside a joint, with its own failure budget.
struct JointLink {
    std::shared_ptr<chrono::ChLinkMateGeneric> link;
    /// Translation-only joint used by triangle shells, which have no rotation.
    std::shared_ptr<chrono::fea::ChLinkNodeNode> pointLink;
    /// Position tie between a shell node (xyzrot) and a solid node (xyz).
    std::shared_ptr<chrono::fea::ChLinkNodeFrame> frameLink;
    /// xyz-xyz tie that keeps the built separation instead of collapsing it.
    /// Not a coincident point link: the nodes stay where the mesh put them.
    std::shared_ptr<chrono::fea::ChLinkNodeNode> offsetLink;
    std::size_t partA = 0;
    std::size_t nodeA = 0;
    std::size_t partB = 0;
    std::size_t nodeB = 0;
    double breakForce = 0.0;
    double breakTorque = 0.0;
    /// Peel capacity, N. Zero means the link is rated on |force| alone.
    double breakPeel = 0.0;
    /// Opening already present when the link was tied. Later damage is only
    /// growth past this, so a midsurface seam is not born broken.
    double restGap = 0.0;
    /// Opening at which the cohesive law is exhausted, metres.
    double criticalOpening = 0.0;
    chrono::ChVector3d peelAxis{0.0, 1.0, 0.0};
    double filteredForce = 0.0;
    chrono::ChVector3d filteredForceVec{0, 0, 0};
    double filteredTorque = 0.0;
    double utilisation = 0.0;
    int overloadSteps = 0;
    bool locksRotation = false;
    bool broken = false;
};

struct Joint {
    JointKind kind = JointKind::Weld;
    std::string label;
    std::size_t partA = 0;
    std::size_t partB = 0;
    std::vector<JointLink> links;

    int brokenCount() const {
        int n = 0;
        for (const auto& l : links)
            n += l.broken ? 1 : 0;
        return n;
    }

    bool severed() const {
        return links.empty() || brokenCount() == static_cast<int>(links.size());
    }

    double peakUtilisation() const {
        double u = 0.0;
        for (const auto& l : links)
            if (!l.broken && l.utilisation > u)
                u = l.utilisation;
        return u;
    }
};

}  // namespace rcs
