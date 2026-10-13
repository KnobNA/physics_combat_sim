#pragma once

#include "physics_core/part.hpp"
#include "physics_core/recording.hpp"
#include "physics_core/scenario.hpp"

#include <vector>

namespace rcs {

/// One thermally tied pair: a live weld or screw, or two nodes in contact.
struct HeatBridge {
    std::size_t partA = 0;
    std::size_t nodeA = 0;
    std::size_t partB = 0;
    std::size_t nodeB = 0;
};

/// Still-air grid plus the slow clock that advances node and cell temperatures.
/// An idle scene, with no source on and every node at its initial temperature,
/// returns without touching the mesh.
class HeatField {
public:
    /// Advance conduction, convection, and in-sphere exchange. `dt` is the
    /// mechanics step; heat itself moves on a coarser clock. Returns false
    /// when nothing changed.
    bool advance(double dt,
                 double time,
                 std::vector<Part>& parts,
                 const Environment& environment,
                 const std::vector<HeatSource>& sources,
                 const std::vector<HeatBridge>& bridges);

    void writeAir(AirFrame& out) const;

    /// True when the last advance moved lagged strain temperature, so the
    /// mechanical step still has an unbalanced thermal residual.
    bool strainMoved() const { return _strainMoved; }
    /// True while strain temperature is still catching node temperature.
    bool easing() const { return _easing; }

private:
    struct Cache {
        bool on = false;
        int nodeStamp = -1;
        std::vector<int> part;
        std::vector<int> node;
        std::vector<int> cell;
    };

    void configure(const Environment& environment);
    int cellAt(const chrono::ChVector3d& point) const;
    chrono::ChVector3d cellCenter(int index) const;
    void rebuildCache(Cache& cache, const HeatSource& source, const std::vector<Part>& parts, int stamp);
    void exchangeEdge(Part& part, int a, int b, double hdt, bool volume = false);
    void exchangePair(Part& a, std::size_t ia, Part& b, std::size_t ib, double conductance, double hdt);
    void exchangeContact(std::vector<Part>& parts, double hdt);
    void exchangeRadiation(std::vector<Part>& parts, const HeatSource& source, double hdt);
    void wake(std::vector<Part>& parts);
    bool easeStrain(std::vector<Part>& parts, double dt);

    int _nx = 0;
    int _ny = 0;
    int _nz = 0;
    double _originX = 0.0;
    double _originY = 0.0;
    double _originZ = 0.0;
    double _dx = 0.25;
    double _dy = 0.25;
    double _dz = 0.25;
    double _ambient = 293.15;
    double _spanX = 0.0;
    double _spanY = 0.0;
    double _spanZ = 0.0;
    double _floor = 0.0;
    std::vector<float> _air;
    std::vector<char> _airActive;
    std::vector<int> _airList;
    std::vector<Cache> _caches;
    double _wait = 0.0;
    double _airWait = 0.0;
    bool _sourceOn = false;
    bool _dirty = false;
    bool _easing = false;
    bool _strainMoved = false;
};

}  // namespace rcs
