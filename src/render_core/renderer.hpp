#pragma once

#include "physics_core/recording.hpp"
#include "physics_core/scenario.hpp"

#include <Magnum/GL/Buffer.h>
#include <Magnum/GL/Mesh.h>
#include <Magnum/Math/Matrix4.h>
#include <Magnum/Math/Range.h>
#include <Magnum/Math/Vector3.h>
#include <Magnum/Shaders/FlatGL.h>
#include <Magnum/Shaders/PhongGL.h>
#include <Magnum/Tags.h>

#include <cstddef>
#include <memory>
#include <optional>
#include <vector>

namespace rcs {

/// How part surfaces are tinted.
enum class ShadeMode {
    Material,     ///< flat metal colour
    Stress,       ///< von Mises colour ramp
    Temperature   ///< absolute kelvin ramp, 0 K is always dark blue
};

class SceneRenderer {
public:
    SceneRenderer();

    /// Rebuild GPU buffers for a new topology. Call whenever the part list or
    /// mesh resolution changes.
    void setLayout(const RecordingLayout& layout);

    /// Push one frame of node positions and stresses to the GPU.
    void update(const RecordingLayout& layout,
                const Frame& frame,
                ShadeMode mode,
                float stressScale,
                const std::vector<std::size_t>& highlighted,
                bool showAir);

    void draw(const Magnum::Matrix4& view,
              const Magnum::Matrix4& projection,
              const Environment& environment,
              const std::vector<HeatSource>& sources,
              double time);

    bool partBounds(std::size_t partIndex, Magnum::Range3D& out) const;
    std::size_t partCount() const { return _parts.size(); }

    /// Colour used for a normalised stress value in [0, 1].
    static Magnum::Color3 stressColor(float normalised);
    /// Absolute temperature. 0 K is dark blue; the stops do not move with the scene.
    static Magnum::Color3 temperatureColor(float kelvin);

private:
    struct PartMesh {
        Magnum::GL::Buffer vertices{Magnum::NoCreate};
        Magnum::GL::Mesh mesh{Magnum::NoCreate};
        Magnum::Range3D bounds;
    };

    void rebuildArena(const Environment& environment);

    Magnum::Shaders::PhongGL _flat;
    Magnum::Shaders::PhongGL _shaded;
    Magnum::Shaders::FlatGL3D _points;
    Magnum::GL::Mesh _floorMesh;
    Magnum::GL::Buffer _gridVertices{Magnum::NoCreate};
    Magnum::GL::Mesh _gridMesh{Magnum::NoCreate};
    Magnum::GL::Buffer _edgeVertices{Magnum::NoCreate};
    Magnum::GL::Mesh _edgeMesh{Magnum::NoCreate};
    Magnum::GL::Buffer _airVertices{Magnum::NoCreate};
    Magnum::GL::Mesh _airMesh{Magnum::NoCreate};
    Magnum::GL::Mesh _heatSphere{Magnum::NoCreate};
    double _arenaWidth = -1.0;
    double _arenaHeight = -1.0;
    double _arenaDepth = -1.0;
    double _arenaFloor = -1.0e9;
    std::vector<std::unique_ptr<PartMesh>> _parts;
};

}  // namespace rcs
