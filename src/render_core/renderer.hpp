#pragma once

#include "physics_core/recording.hpp"
#include "physics_core/scenario.hpp"

#include <Magnum/GL/Buffer.h>
#include <Magnum/GL/Mesh.h>
#include <Magnum/Math/Matrix4.h>
#include <Magnum/Math/Range.h>
#include <Magnum/Math/Vector3.h>
#include <Magnum/Shaders/PhongGL.h>
#include <Magnum/Tags.h>

#include <cstddef>
#include <memory>
#include <optional>
#include <vector>

namespace rcs {

/// How part surfaces are tinted.
enum class ShadeMode {
    Material,  ///< flat metal colour
    Stress     ///< von Mises colour ramp
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
                std::optional<std::size_t> selectedPart);

    void draw(const Magnum::Matrix4& view,
              const Magnum::Matrix4& projection,
              const Environment& environment);

    bool partBounds(std::size_t partIndex, Magnum::Range3D& out) const;
    std::size_t partCount() const { return _parts.size(); }

    /// Colour used for a normalised stress value in [0, 1].
    static Magnum::Color3 stressColor(float normalised);

private:
    struct PartMesh {
        Magnum::GL::Buffer vertices{Magnum::NoCreate};
        Magnum::GL::Mesh mesh{Magnum::NoCreate};
        Magnum::Range3D bounds;
    };

    Magnum::Shaders::PhongGL _flat;
    Magnum::Shaders::PhongGL _shaded;
    Magnum::GL::Mesh _floorMesh;
    std::vector<std::unique_ptr<PartMesh>> _parts;
};

}  // namespace rcs
