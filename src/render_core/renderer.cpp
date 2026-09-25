#include "render_core/renderer.hpp"

#include <Corrade/Containers/ArrayViewStl.h>
#include <Magnum/GL/Renderer.h>
#include <Magnum/Math/Color.h>
#include <Magnum/Math/Functions.h>
#include <Magnum/MeshTools/Compile.h>
#include <Magnum/Primitives/Cube.h>
#include <Magnum/Primitives/UVSphere.h>
#include <Magnum/Trade/MeshData.h>

#include <algorithm>
#include <cmath>
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

struct AirPoint {
    Vector3 position;
    Color4 color;
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

void addDoubleSidedQuad(std::vector<Vertex>& out,
                        const Vector3& q0,
                        const Vector3& q1,
                        const Vector3& q2,
                        const Vector3& q3,
                        const Color3& color)
{
    Vector3 n = Math::cross(q1 - q0, q3 - q0);
    if (n.length() < 1.0e-12f)
        return;
    n = n.normalized();
    const Vector3 tri[6] = {q0, q1, q2, q0, q2, q3};
    for (int i = 0; i < 6; ++i)
        out.push_back({tri[i], n, color});
    for (int i = 5; i >= 0; --i)
        out.push_back({tri[i], -n, color});
}

void addRod(std::vector<Vertex>& out, const Vector3& a, const Vector3& b, float half, const Color3& color)
{
    const Vector3 dir = b - a;
    if (dir.length() < 1.0e-6f)
        return;
    const Vector3 axis = dir.normalized();
    const Vector3 helper = Math::abs(axis.y()) < 0.9f ? Vector3::yAxis() : Vector3::xAxis();
    const Vector3 side = Math::cross(axis, helper).normalized() * half;
    const Vector3 up = Math::cross(axis, side).normalized() * half;
    const Vector3 p[8] = {
        a - side - up, a + side - up, a + side + up, a - side + up,
        b - side - up, b + side - up, b + side + up, b - side + up
    };
    const int faces[6][4] = {
        {0, 1, 2, 3}, {4, 7, 6, 5}, {0, 4, 5, 1},
        {1, 5, 6, 2}, {2, 6, 7, 3}, {3, 7, 4, 0}
    };
    for (const auto& face : faces)
        addDoubleSidedQuad(out, p[face[0]], p[face[1]], p[face[2]], p[face[3]], color);
}

}  // namespace

SceneRenderer::SceneRenderer():
    _flat{Shaders::PhongGL::Configuration{}.setLightCount(1)},
    _shaded{Shaders::PhongGL::Configuration{}
                .setLightCount(1)
                .setFlags(Shaders::PhongGL::Flag::VertexColor)},
    _points{Shaders::FlatGL3D::Configuration{}.setFlags(Shaders::FlatGL3D::Flag::VertexColor)},
    _floorMesh{MeshTools::compile(Primitives::cubeSolid())},
    _gridVertices{},
    _gridMesh{},
    _edgeVertices{},
    _edgeMesh{},
    _airVertices{},
    _airMesh{},
    _heatSphere{MeshTools::compile(Primitives::uvSphereSolid(10, 16))}
{
    _gridMesh.setPrimitive(GL::MeshPrimitive::Triangles)
        .addVertexBuffer(_gridVertices, 0, Shaders::PhongGL::Position{},
                         Shaders::PhongGL::Normal{}, Shaders::PhongGL::Color3{});
    _edgeMesh.setPrimitive(GL::MeshPrimitive::Triangles)
        .addVertexBuffer(_edgeVertices, 0, Shaders::PhongGL::Position{},
                         Shaders::PhongGL::Normal{}, Shaders::PhongGL::Color3{});
    _airMesh.setPrimitive(GL::MeshPrimitive::Points)
        .addVertexBuffer(_airVertices, 0, Shaders::FlatGL3D::Position{}, Shaders::FlatGL3D::Color4{});
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

Color3 SceneRenderer::temperatureColor(float kelvin) {
    const float stopsK[7] = {0.0f, 200.0f, 293.0f, 450.0f, 700.0f, 1100.0f, 2000.0f};
    const Color3 stopsC[7] = {
        Color3{0.02f, 0.02f, 0.16f},
        Color3{0.05f, 0.18f, 0.78f},
        Color3{0.20f, 0.62f, 0.88f},
        Color3{0.25f, 0.78f, 0.38f},
        Color3{0.96f, 0.84f, 0.16f},
        Color3{0.95f, 0.28f, 0.06f},
        Color3{1.0f, 0.96f, 0.92f}
    };
    if (kelvin <= stopsK[0])
        return stopsC[0];
    if (kelvin >= stopsK[6])
        return stopsC[6];
    int index = 0;
    while (index < 5 && kelvin > stopsK[index + 1])
        ++index;
    const float span = stopsK[index + 1] - stopsK[index];
    const float t = span > 0.0f ? (kelvin - stopsK[index]) / span : 0.0f;
    return Math::lerp(stopsC[index], stopsC[index + 1], t);
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
                           const std::vector<std::size_t>& highlighted,
                           bool showAir)
{
    if (_parts.size() != layout.parts.size())
        setLayout(layout);

    const float scale = stressScale > 1.0f ? stressScale : 1.0f;

    for (std::size_t i = 0; i < layout.parts.size() && i < frame.parts.size(); ++i) {
        const PartLayout& partLayout = layout.parts[i];
        const PartFrame& partFrame = frame.parts[i];
        const bool selected = std::find(highlighted.begin(), highlighted.end(), i) != highlighted.end();

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
            if (mode == ShadeMode::Temperature && colorIndex < partFrame.temperature.size()) {
                color = temperatureColor(partFrame.temperature[colorIndex]);
                if (selected)
                    color = Math::lerp(color, 0xffffff_rgbf, 0.25f);
            } else if (mode == ShadeMode::Stress && colorIndex < partFrame.stress.size()) {
                color = stressColor(partFrame.stress[colorIndex] / scale);
                if (selected)
                    color = Math::lerp(color, 0xffffff_rgbf, 0.25f);
            } else if (mode == ShadeMode::Stress && colorIndex < partFrame.plasticStrain.size() &&
                       partFrame.plasticStrain[colorIndex] > 0.0f) {
                color = Math::lerp(color, 0x8d6e63_rgbf, 0.35f);
            }

            const bool flat = i2 == i3 && (partLayout.volume || static_cast<int>(colorIndex) < partLayout.flatTriangles);
            if (flat) {
                Vector3 n = Math::cross(b - a, c - a);
                if (n.length() > 1.0e-12f) {
                    n = n.normalized();
                    vertices.push_back({a, n, color});
                    vertices.push_back({b, n, color});
                    vertices.push_back({c, n, color});
                    vertices.push_back({a, -n, color});
                    vertices.push_back({c, -n, color});
                    vertices.push_back({b, -n, color});
                }
            } else {
                appendPlateBox(vertices, a, b, c, d, partLayout.thickness * 0.5f, color);
            }
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

    std::vector<AirPoint> air;
    const AirFrame& sample = frame.air;
    if (showAir && !sample.temperature.empty() && sample.nx > 0 && sample.ny > 0 && sample.nz > 0) {
        const int count = sample.nx * sample.ny * sample.nz;
        // A few samples across each cell, so the field reads as a grid of
        // points rather than one box per cell.
        constexpr int kSub = 3;
        for (int index = 0; index < count && index < int(sample.temperature.size()); ++index) {
            const float kelvin = sample.temperature[static_cast<std::size_t>(index)];
            const float warmth = std::abs(kelvin - sample.ambient);
            if (warmth < 0.5f)
                continue;
            const float alpha = Math::clamp(warmth / 280.0f, 0.08f, 0.92f);
            const Color3 rgb = temperatureColor(kelvin);
            const Color4 color{rgb.r(), rgb.g(), rgb.b(), alpha};
            const int ix = index % sample.nx;
            const int iy = (index / sample.nx) % sample.ny;
            const int iz = index / (sample.nx * sample.ny);
            const Vector3 origin{sample.originX + ix * sample.dx,
                                 sample.originY + iy * sample.dy,
                                 sample.originZ + iz * sample.dz};
            const Vector3 step{sample.dx / kSub, sample.dy / kSub, sample.dz / kSub};
            for (int sx = 0; sx < kSub; ++sx) {
                for (int sy = 0; sy < kSub; ++sy) {
                    for (int sz = 0; sz < kSub; ++sz) {
                        const Vector3 position = origin + Vector3{(sx + 0.5f) * step.x(),
                                                                  (sy + 0.5f) * step.y(),
                                                                  (sz + 0.5f) * step.z()};
                        air.push_back({position, color});
                    }
                }
            }
        }
    }
    _airVertices.setData(air, GL::BufferUsage::DynamicDraw);
    _airMesh.setCount(Int(air.size()));
}

void SceneRenderer::rebuildArena(const Environment& environment)
{
    const double width = std::max(1.0, environment.boundsSize.x());
    const double height = std::max(1.0, environment.boundsSize.y());
    const double depth = std::max(1.0, environment.boundsSize.z());
    if (width == _arenaWidth && height == _arenaHeight && depth == _arenaDepth &&
        environment.floorHeight == _arenaFloor && _gridMesh.count() > 0)
        return;
    _arenaWidth = width;
    _arenaHeight = height;
    _arenaDepth = depth;
    _arenaFloor = environment.floorHeight;

    const float hx = Float(width * 0.5);
    const float hz = Float(depth * 0.5);
    const float y0 = Float(environment.floorHeight);
    const float y1 = y0 + Float(height);
    const float yGrid = y0 + 0.02f;
    const Color3 gridColor = 0xc5c8cc_rgbf;
    const Color3 edgeColor = 0xffffff_rgbf;

    std::vector<Vertex> grid;
    const auto span = [&](float limit, float step, auto emit) {
        const float start = std::ceil((-limit) / step - 1.0e-4f) * step;
        for (float t = start; t <= limit + 1.0e-3f; t += step)
            emit(t);
    };
    const float stepX = Float(Environment::gridStep(width));
    const float stepZ = Float(Environment::gridStep(depth));
    span(hx, stepX, [&](float x) {
        addRod(grid, {x, yGrid, -hz}, {x, yGrid, hz}, 0.004f, gridColor);
    });
    span(hz, stepZ, [&](float z) {
        addRod(grid, {-hx, yGrid, z}, {hx, yGrid, z}, 0.004f, gridColor);
    });
    _gridVertices.setData(grid, GL::BufferUsage::StaticDraw);
    _gridMesh.setCount(Int(grid.size()));

    std::vector<Vertex> edges;
    const Vector3 corner[8] = {
        {-hx, y0, -hz}, {hx, y0, -hz}, {hx, y0, hz}, {-hx, y0, hz},
        {-hx, y1, -hz}, {hx, y1, -hz}, {hx, y1, hz}, {-hx, y1, hz}
    };
    const int index[12][2] = {
        {0, 1}, {1, 2}, {2, 3}, {3, 0},
        {4, 5}, {5, 6}, {6, 7}, {7, 4},
        {0, 4}, {1, 5}, {2, 6}, {3, 7}
    };
    for (const auto& edge : index)
        addRod(edges, corner[edge[0]], corner[edge[1]], 0.012f, edgeColor);
    _edgeVertices.setData(edges, GL::BufferUsage::StaticDraw);
    _edgeMesh.setCount(Int(edges.size()));
}

void SceneRenderer::draw(const Matrix4& view,
                         const Matrix4& projection,
                         const Environment& environment,
                         const std::vector<HeatSource>& sources,
                         double time)
{
    if (environment.floorEnabled) {
        const float hx = Float(std::max(1.0, environment.boundsSize.x()) * 0.5);
        const float hz = Float(std::max(1.0, environment.boundsSize.z()) * 0.5);
        const float thickness = 0.05f;
        const Matrix4 floorTransform =
            Matrix4::translation({0.0f, Float(environment.floorHeight) - thickness, 0.0f}) *
            Matrix4::scaling({hx, thickness, hz});
        _flat.setProjectionMatrix(projection)
            .setTransformationMatrix(view * floorTransform)
            .setNormalMatrix((view * floorTransform).normalMatrix())
            .setDiffuseColor(0x6d6d6d_rgbf)
            .setAmbientColor(0x2a2a2a_rgbf)
            .draw(_floorMesh);

        rebuildArena(environment);
        _shaded.setProjectionMatrix(projection)
            .setTransformationMatrix(view)
            .setNormalMatrix(view.normalMatrix())
            .setDiffuseColor(0xffffff_rgbf)
            .setAmbientColor(0xffffff_rgbf);
        if (_gridMesh.count() > 0)
            _shaded.draw(_gridMesh);
        if (_edgeMesh.count() > 0)
            _shaded.draw(_edgeMesh);
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

    const bool air = _airMesh.count() > 0;
    bool sphere = false;
    for (const HeatSource& source : sources) {
        if (source.enabled && time >= source.startTime && time < source.startTime + source.duration) {
            sphere = true;
            break;
        }
    }
    if (!air && !sphere)
        return;

    GL::Renderer::enable(GL::Renderer::Feature::Blending);
    GL::Renderer::setBlendFunction(GL::Renderer::BlendFunction::SourceAlpha,
                                   GL::Renderer::BlendFunction::OneMinusSourceAlpha);
    GL::Renderer::setDepthMask(false);
    if (air) {
        GL::Renderer::setPointSize(7.0f);
        _points.setTransformationProjectionMatrix(projection * view)
            .setColor(Color4{1.0f})
            .draw(_airMesh);
    }
    if (sphere) {
        _flat.setProjectionMatrix(projection)
            .setDiffuseColor(Color4{1.0f, 0.55f, 0.2f, 0.28f})
            .setAmbientColor(Color4{1.0f, 0.45f, 0.15f, 0.28f});
        for (const HeatSource& source : sources) {
            if (!source.enabled || time < source.startTime || time >= source.startTime + source.duration)
                continue;
            const float radius = Float(std::max(source.radius, 1.0e-3));
            const Matrix4 transform =
                Matrix4::translation({Float(source.point.x()), Float(source.point.y()), Float(source.point.z())}) *
                Matrix4::scaling(Vector3{radius});
            _flat.setTransformationMatrix(view * transform)
                .setNormalMatrix((view * transform).normalMatrix())
                .draw(_heatSphere);
        }
    }
    GL::Renderer::setDepthMask(true);
    GL::Renderer::disable(GL::Renderer::Feature::Blending);
}

bool SceneRenderer::partBounds(std::size_t partIndex, Range3D& out) const {
    if (partIndex >= _parts.size())
        return false;
    out = _parts[partIndex]->bounds;
    return true;
}

}  // namespace rcs
