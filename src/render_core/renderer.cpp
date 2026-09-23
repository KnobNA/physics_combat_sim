#include "render_core/renderer.hpp"

#include <Corrade/Containers/ArrayViewStl.h>
#include <Magnum/GL/Renderer.h>
#include <Magnum/Math/Color.h>
#include <Magnum/Math/Functions.h>
#include <Magnum/MeshTools/Compile.h>
#include <Magnum/Primitives/Cube.h>
#include <Magnum/Trade/MeshData.h>

#include <algorithm>
#include <limits>

namespace rcs {

using namespace Magnum;
using namespace Magnum::Math::Literals;

namespace {

struct Vertex {
    Vector3 position;
    Vector3 normal;
    Color3 color;
};

Vector3 toMagnum(const Vec3f& v) {
    return {v.x, v.y, v.z};
}

/// Expand a mid-surface quad into a closed box of the part's real thickness,
/// so panels look and pick like solid plate rather than zero-width sheets.
void appendPlateBox(std::vector<Vertex>& out,
                    const Vector3& a,
                    const Vector3& b,
                    const Vector3& c,
                    const Vector3& d,
                    float halfThickness,
                    const Color3& color)
{
    Vector3 n = Math::cross(b - a, d - a);
    if (n.length() < 1.0e-12f)
        return;
    n = n.normalized();
    const Vector3 offset = n * halfThickness;
    const Vector3 p[8] = {
        a - offset, b - offset, c - offset, d - offset,
        a + offset, b + offset, c + offset, d + offset
    };

    const int faces[6][4] = {
        {0, 3, 2, 1}, {4, 5, 6, 7}, {0, 1, 5, 4},
        {1, 2, 6, 5}, {2, 3, 7, 6}, {3, 0, 4, 7}
    };

    for (const auto& face : faces) {
        const Vector3& q0 = p[face[0]];
        const Vector3& q1 = p[face[1]];
        const Vector3& q2 = p[face[2]];
        const Vector3& q3 = p[face[3]];
        Vector3 fn = Math::cross(q1 - q0, q3 - q0);
        if (fn.length() < 1.0e-12f)
            continue;
        fn = fn.normalized();
        out.push_back({q0, fn, color});
        out.push_back({q1, fn, color});
        out.push_back({q2, fn, color});
        out.push_back({q0, fn, color});
        out.push_back({q2, fn, color});
        out.push_back({q3, fn, color});
    }
}

}  // namespace

SceneRenderer::SceneRenderer():
    _flat{Shaders::PhongGL::Configuration{}.setLightCount(1)},
    _shaded{Shaders::PhongGL::Configuration{}
                .setLightCount(1)
                .setFlags(Shaders::PhongGL::Flag::VertexColor)},
    _floorMesh{MeshTools::compile(Primitives::cubeSolid())}
{
    for (Shaders::PhongGL* shader : {&_flat, &_shaded}) {
        shader->setLightPositions({Vector4{8.0f, 14.0f, 6.0f, 0.0f}})
            .setLightColors({0xffffff_rgbf})
            .setSpecularColor(0x222222_rgbf)
            .setShininess(40.0f);
    }
}

Color3 SceneRenderer::stressColor(float normalised) {
    const float t = Math::clamp(normalised, 0.0f, 1.0f);
    // Cool-to-hot ramp with distinguishable bands: navy, cyan, green, amber, red.
    const Color3 stops[5] = {
        0x203a8f_rgbf, 0x00b7c2_rgbf, 0x3ec46d_rgbf, 0xf2b134_rgbf, 0xd64545_rgbf
    };
    const float scaled = t * 4.0f;
    const int index = Math::min(int(scaled), 3);
    return Math::lerp(stops[index], stops[index + 1], scaled - float(index));
}

void SceneRenderer::setLayout(const RecordingLayout& layout) {
    _parts.clear();
    _parts.reserve(layout.parts.size());
    for (std::size_t i = 0; i < layout.parts.size(); ++i) {
        auto part = std::make_unique<PartMesh>();
        part->vertices = GL::Buffer{};
        part->mesh = GL::Mesh{};
        part->mesh.setPrimitive(GL::MeshPrimitive::Triangles)
            .addVertexBuffer(part->vertices, 0,
                             Shaders::PhongGL::Position{},
                             Shaders::PhongGL::Normal{},
                             Shaders::PhongGL::Color3{});
        _parts.push_back(std::move(part));
    }
}

void SceneRenderer::update(const RecordingLayout& layout,
                           const Frame& frame,
                           ShadeMode mode,
                           float stressScale,
                           std::optional<std::size_t> selectedPart)
{
    if (_parts.size() != layout.parts.size())
        setLayout(layout);

    const float scale = stressScale > 1.0f ? stressScale : 1.0f;

    for (std::size_t i = 0; i < layout.parts.size() && i < frame.parts.size(); ++i) {
        const PartLayout& partLayout = layout.parts[i];
        const PartFrame& partFrame = frame.parts[i];
        const bool selected = selectedPart && *selectedPart == i;

        std::vector<Vertex> vertices;
        vertices.reserve(partLayout.quads.size() * 36);

        Vector3 lo{std::numeric_limits<float>::max()};
        Vector3 hi{std::numeric_limits<float>::lowest()};

        const auto& quads = partFrame.quads.empty() ? partLayout.quads : partFrame.quads;
        const auto& triangles = partFrame.triangles.empty() ? partLayout.triangles : partFrame.triangles;
        const std::size_t nodes = partFrame.nodeCount();

        auto drawQuad = [&](int i0, int i1, int i2, int i3, std::size_t colorIndex) {
            if (i0 < 0 || i1 < 0 || i2 < 0 || i3 < 0)
                return;
            if (static_cast<std::size_t>(i0) >= nodes || static_cast<std::size_t>(i1) >= nodes ||
                static_cast<std::size_t>(i2) >= nodes || static_cast<std::size_t>(i3) >= nodes)
                return;
            const Vector3 a = toMagnum(partFrame.position(static_cast<std::size_t>(i0)));
            const Vector3 b = toMagnum(partFrame.position(static_cast<std::size_t>(i1)));
            const Vector3 c = toMagnum(partFrame.position(static_cast<std::size_t>(i2)));
            const Vector3 d = toMagnum(partFrame.position(static_cast<std::size_t>(i3)));

            Color3 color = selected ? 0xd4af37_rgbf : 0xb0b7be_rgbf;
            if (mode == ShadeMode::Stress && colorIndex < partFrame.stress.size()) {
                color = stressColor(partFrame.stress[colorIndex] / scale);
                if (selected)
                    color = Math::lerp(color, 0xffffff_rgbf, 0.25f);
            } else if (mode == ShadeMode::Stress && colorIndex < partFrame.plasticStrain.size() &&
                       partFrame.plasticStrain[colorIndex] > 0.0f) {
                color = Math::lerp(color, 0x8d6e63_rgbf, 0.35f);
            }

            appendPlateBox(vertices, a, b, c, d, partLayout.thickness * 0.5f, color);
            for (const Vector3& p : {a, b, c, d}) {
                lo = Math::min(lo, p);
                hi = Math::max(hi, p);
            }
        };

        if (!quads.empty()) {
            for (std::size_t q = 0; q < quads.size(); ++q)
                drawQuad(quads[q][0], quads[q][1], quads[q][2], quads[q][3], q);
        } else {
            for (std::size_t t = 0; t < triangles.size(); ++t)
                drawQuad(triangles[t][0], triangles[t][1], triangles[t][2], triangles[t][2], t);
        }

        PartMesh& gpu = *_parts[i];
        gpu.vertices.setData(vertices, GL::BufferUsage::DynamicDraw);
        gpu.mesh.setCount(Int(vertices.size()));
        const Vector3 pad{partLayout.thickness};
        gpu.bounds = vertices.empty() ? Range3D{} : Range3D{lo - pad, hi + pad};
    }
}

void SceneRenderer::draw(const Matrix4& view,
                         const Matrix4& projection,
                         const Environment& environment)
{
    if (environment.floorEnabled) {
        const float half = Float(environment.floorHalfSize);
        const float thickness = 0.05f;
        const Matrix4 floorTransform =
            Matrix4::translation({0.0f, Float(environment.floorHeight) - thickness, 0.0f}) *
            Matrix4::scaling({half, thickness, half});
        _flat.setProjectionMatrix(projection)
            .setTransformationMatrix(view * floorTransform)
            .setNormalMatrix((view * floorTransform).normalMatrix())
            .setDiffuseColor(0x6d6d6d_rgbf)
            .setAmbientColor(0x2a2a2a_rgbf)
            .draw(_floorMesh);
    }

    _shaded.setProjectionMatrix(projection)
        .setTransformationMatrix(view)
        .setNormalMatrix(view.normalMatrix())
        .setDiffuseColor(0xffffff_rgbf)
        .setAmbientColor(0x3a3d42_rgbf);

    for (const auto& part : _parts) {
        if (part->mesh.count() == 0)
            continue;
        _shaded.draw(part->mesh);
    }
}

bool SceneRenderer::partBounds(std::size_t partIndex, Range3D& out) const {
    if (partIndex >= _parts.size())
        return false;
    out = _parts[partIndex]->bounds;
    return true;
}

}  // namespace rcs
