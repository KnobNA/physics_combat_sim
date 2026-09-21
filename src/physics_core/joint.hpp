#pragma once

#include "chrono/core/ChVector3.h"
#include "chrono/physics/ChLinkMate.h"

#include <cstddef>
#include <memory>
#include <string>
#include <vector>

namespace rcs {

/// A manual weld. The user states which two parts are joined; every pair of
/// coincident boundary nodes inside the (optional) region becomes one weld
/// link, so a bead can fail progressively instead of all at once.
struct WeldDefinition {
    std::string label = "Weld";
    std::size_t partA = 0;
    std::size_t partB = 1;
    double tolerance = 2.0e-3;   ///< node pairing distance, m
    double throat = 0.0;         ///< bead throat, m; 0 means min part thickness
    double strengthScale = 1.0;  ///< derating factor for a poor bead
    bool useRegion = false;
    chrono::ChVector3d regionMin{0, 0, 0};
    chrono::ChVector3d regionMax{0, 0, 0};
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
    std::size_t partA = 0;
    std::size_t nodeA = 0;
    std::size_t partB = 0;
    std::size_t nodeB = 0;
    double breakForce = 0.0;
    double breakTorque = 0.0;
    double filteredForce = 0.0;
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
