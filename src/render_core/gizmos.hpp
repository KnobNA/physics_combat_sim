#pragma once

#include "editor/tools.hpp"

#include <Magnum/GL/Mesh.h>
#include <Magnum/Math/Matrix4.h>
#include <Magnum/Math/Vector3.h>
#include <Magnum/Shaders/PhongGL.h>

namespace rcs {

enum class GizmoHandle {
    None,
    AxisX,
    AxisY,
    AxisZ
};

class TransformGizmos {
public:
    TransformGizmos();

    void draw(const Magnum::Matrix4& view,
              const Magnum::Matrix4& projection,
              const Magnum::Vector3& origin,
              float scale,
              Tool tool,
              GizmoHandle active);

    GizmoHandle pickMove(const Magnum::Vector3& rayOrigin,
                         const Magnum::Vector3& rayDir,
                         const Magnum::Vector3& origin,
                         float scale) const;

    GizmoHandle pickRotate(const Magnum::Vector3& rayOrigin,
                           const Magnum::Vector3& rayDir,
                           const Magnum::Vector3& origin,
                           float scale) const;

    static Magnum::Vector3 axisVector(GizmoHandle handle);

private:
    Magnum::Shaders::PhongGL _phong;
    Magnum::GL::Mesh _shaft;
    Magnum::GL::Mesh _head;
    Magnum::GL::Mesh _ring;
};

}  // namespace rcs
