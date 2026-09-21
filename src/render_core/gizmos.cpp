#include "render_core/gizmos.hpp"

#include <Corrade/Containers/ArrayViewStl.h>
#include <Magnum/GL/Buffer.h>
#include <Magnum/Math/Angle.h>
#include <Magnum/Math/Color.h>
#include <Magnum/Math/Constants.h>
#include <Magnum/Math/Functions.h>
#include <Magnum/Math/Matrix4.h>
#include <Magnum/MeshTools/Compile.h>
#include <Magnum/Primitives/Cone.h>
#include <Magnum/Primitives/Cylinder.h>
#include <Magnum/Shaders/GenericGL.h>
#include <Magnum/Trade/MeshData.h>

#include <cmath>
#include <vector>

namespace rcs {

using namespace Magnum;
using namespace Magnum::Math::Literals;

namespace {

GL::Mesh makeRingMesh() {
    constexpr UnsignedInt radial = 48;
    constexpr UnsignedInt tube = 8;
    constexpr float R = 1.0f;
    constexpr float r = 0.045f;

    struct Vertex {
        Vector3 position;
        Vector3 normal;
    };
    std::vector<Vertex> vertices;
    vertices.reserve(radial * tube);
    std::vector<UnsignedInt> indices;
    indices.reserve(radial * tube * 6);

    for (UnsignedInt i = 0; i < radial; ++i) {
        const float u = Float(i) / Float(radial) * Math::Constants<Float>::pi() * 2.0f;
        const Vector3 center{Math::cos(Rad{u}) * R, Math::sin(Rad{u}) * R, 0.0f};
        const Vector3 tangent{-Math::sin(Rad{u}), Math::cos(Rad{u}), 0.0f};
        const Vector3 normal = center.normalized();
        const Vector3 binormal = Math::cross(tangent, normal);
        for (UnsignedInt j = 0; j < tube; ++j) {
            const float v = Float(j) / Float(tube) * Math::Constants<Float>::pi() * 2.0f;
            const Vector3 n = normal * Math::cos(Rad{v}) + binormal * Math::sin(Rad{v});
            vertices.push_back({center + n * r, n});
        }
    }

    for (UnsignedInt i = 0; i < radial; ++i) {
        const UnsignedInt i0 = i * tube;
        const UnsignedInt i1 = ((i + 1) % radial) * tube;
        for (UnsignedInt j = 0; j < tube; ++j) {
            const UnsignedInt j1 = (j + 1) % tube;
            indices.push_back(i0 + j);
            indices.push_back(i1 + j);
            indices.push_back(i1 + j1);
            indices.push_back(i0 + j);
            indices.push_back(i1 + j1);
            indices.push_back(i0 + j1);
        }
    }

    GL::Buffer vertexBuffer;
    vertexBuffer.setData(vertices);
    GL::Buffer indexBuffer;
    indexBuffer.setData(indices);
    GL::Mesh mesh;
    mesh.setPrimitive(GL::MeshPrimitive::Triangles)
        .setCount(Int(indices.size()))
        .addVertexBuffer(std::move(vertexBuffer), 0,
                         Shaders::PhongGL::Position{},
                         Shaders::PhongGL::Normal{})
        .setIndexBuffer(std::move(indexBuffer), 0, GL::MeshIndexType::UnsignedInt);
    return mesh;
}

bool closestOnRayAndSegment(const Vector3& rayOrigin,
                            const Vector3& rayDir,
                            const Vector3& segOrigin,
                            const Vector3& segDir,
                            float segLen,
                            float& outRayT,
                            float& outSegT,
                            float& outDist)
{
    const Vector3 w0 = rayOrigin - segOrigin;
    const float a = Math::dot(rayDir, rayDir);
    const float b = Math::dot(rayDir, segDir);
    const float c = Math::dot(segDir, segDir);
    const float d = Math::dot(rayDir, w0);
    const float e = Math::dot(segDir, w0);
    const float denom = a * c - b * b;
    if (std::abs(denom) < 1.0e-8f)
        return false;
    outRayT = (b * e - c * d) / denom;
    outSegT = (a * e - b * d) / denom;
    outSegT = Math::clamp(outSegT, 0.0f, segLen);
    if (outRayT < 0.0f)
        return false;
    const Vector3 pRay = rayOrigin + rayDir * outRayT;
    const Vector3 pSeg = segOrigin + segDir * outSegT;
    outDist = (pRay - pSeg).length();
    return true;
}

}  // namespace

TransformGizmos::TransformGizmos():
    _phong{Shaders::PhongGL::Configuration{}.setLightCount(1)},
    _shaft{MeshTools::compile(Primitives::cylinderSolid(1, 16, 0.5f, Primitives::CylinderFlag::CapEnds))},
    _head{MeshTools::compile(Primitives::coneSolid(1, 16, 0.5f, Primitives::ConeFlag::CapEnd))},
    _ring{makeRingMesh()}
{
    _phong.setLightPositions({Vector4{4.0f, 8.0f, 4.0f, 0.0f}})
        .setSpecularColor(0x111111_rgbf)
        .setShininess(20.0f);
}

Vector3 TransformGizmos::axisVector(GizmoHandle handle) {
    switch (handle) {
        case GizmoHandle::AxisX: return Vector3::xAxis();
        case GizmoHandle::AxisY: return Vector3::yAxis();
        case GizmoHandle::AxisZ: return Vector3::zAxis();
        default: return {};
    }
}

void TransformGizmos::draw(const Matrix4& view,
                           const Matrix4& projection,
                           const Vector3& origin,
                           float scale,
                           Tool tool,
                           GizmoHandle active)
{
    const Color3 colors[] = {0xe74c3c_rgbf, 0x27ae60_rgbf, 0x3498db_rgbf};
    const Vector3 axes[] = {Vector3::xAxis(), Vector3::yAxis(), Vector3::zAxis()};
    _phong.setProjectionMatrix(projection);

    auto alignYTo = [](const Vector3& axis) {
        const Vector3 y = Vector3::yAxis();
        const Vector3 a = axis.normalized();
        if (Math::dot(y, a) > 0.999f)
            return Matrix4{};
        if (Math::dot(y, a) < -0.999f)
            return Matrix4::rotationX(180.0_degf);
        const Vector3 rotAxis = Math::cross(y, a).normalized();
        const Rad angle{std::acos(Math::clamp(Math::dot(y, a), -1.0f, 1.0f))};
        return Matrix4::rotation(angle, rotAxis);
    };

    if (tool == Tool::Move) {
        for (int i = 0; i < 3; ++i) {
            const bool hot = active == GizmoHandle(int(GizmoHandle::AxisX) + i);
            const Matrix4 align = alignYTo(axes[i]);
            const Matrix4 shaftTf = Matrix4::translation(origin) * align *
                Matrix4::translation({0.0f, scale * 0.28f, 0.0f}) *
                Matrix4::scaling({scale * 0.025f, scale * 0.28f, scale * 0.025f});
            const Matrix4 headTf = Matrix4::translation(origin) * align *
                Matrix4::translation({0.0f, scale * 0.62f, 0.0f}) *
                Matrix4::scaling({scale * 0.055f, scale * 0.12f, scale * 0.055f});
            const Color3 color = hot ? Color3{1.0f} : colors[i];
            _phong.setDiffuseColor(color)
                .setAmbientColor(color * 0.3f)
                .setTransformationMatrix(view * shaftTf)
                .setNormalMatrix((view * shaftTf).normalMatrix())
                .draw(_shaft);
            _phong.setTransformationMatrix(view * headTf)
                .setNormalMatrix((view * headTf).normalMatrix())
                .draw(_head);
        }
    } else if (tool == Tool::Rotate) {
        const Matrix4 ringAlign[] = {
            Matrix4::rotationY(90.0_degf),
            Matrix4::rotationX(90.0_degf),
            Matrix4{}
        };
        for (int i = 0; i < 3; ++i) {
            const bool hot = active == GizmoHandle(int(GizmoHandle::AxisX) + i);
            const Matrix4 tf = Matrix4::translation(origin) * ringAlign[i] * Matrix4::scaling(Vector3{scale});
            const Color3 color = hot ? Color3{1.0f} : colors[i];
            _phong.setDiffuseColor(color)
                .setAmbientColor(color * 0.25f)
                .setTransformationMatrix(view * tf)
                .setNormalMatrix((view * tf).normalMatrix())
                .draw(_ring);
        }
    }
}

GizmoHandle TransformGizmos::pickMove(const Vector3& rayOrigin,
                                      const Vector3& rayDir,
                                      const Vector3& origin,
                                      float scale) const
{
    GizmoHandle best = GizmoHandle::None;
    float bestDist = scale * 0.08f;
    const Vector3 axes[] = {Vector3::xAxis(), Vector3::yAxis(), Vector3::zAxis()};
    for (int i = 0; i < 3; ++i) {
        float rayT = 0, segT = 0, dist = 0;
        if (!closestOnRayAndSegment(rayOrigin, rayDir, origin, axes[i], scale, rayT, segT, dist))
            continue;
        if (dist < bestDist) {
            bestDist = dist;
            best = GizmoHandle(int(GizmoHandle::AxisX) + i);
        }
    }
    return best;
}

GizmoHandle TransformGizmos::pickRotate(const Vector3& rayOrigin,
                                        const Vector3& rayDir,
                                        const Vector3& origin,
                                        float scale) const
{
    GizmoHandle best = GizmoHandle::None;
    float bestErr = scale * 0.12f;
    const Vector3 axes[] = {Vector3::xAxis(), Vector3::yAxis(), Vector3::zAxis()};
    for (int i = 0; i < 3; ++i) {
        const float denom = Math::dot(rayDir, axes[i]);
        if (std::abs(denom) < 1.0e-4f)
            continue;
        const float t = Math::dot(origin - rayOrigin, axes[i]) / denom;
        if (t < 0.0f)
            continue;
        const Vector3 hit = rayOrigin + rayDir * t;
        const float err = std::abs((hit - origin).length() - scale);
        if (err < bestErr) {
            bestErr = err;
            best = GizmoHandle(int(GizmoHandle::AxisX) + i);
        }
    }
    return best;
}

}  // namespace rcs
