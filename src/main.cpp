#include "editor/ui.hpp"
#include "physics_core/cad.hpp"
#include "physics_core/mesh_split.hpp"
#include "physics_core/simulator.hpp"
#include "render_core/gizmos.hpp"
#include "render_core/orbit_camera.hpp"
#include "render_core/renderer.hpp"

#include <Magnum/GL/DefaultFramebuffer.h>
#include <Magnum/GL/Renderer.h>
#include <Magnum/ImGuiIntegration/Context.hpp>
#include <Magnum/Math/Color.h>
#include <Magnum/Math/Functions.h>
#include <Magnum/Math/Range.h>
#include <Magnum/Platform/Sdl2Application.h>
#include <Magnum/Timeline.h>

#include "chrono/core/ChQuaternion.h"
#include "chrono/core/ChRotation.h"

#include <algorithm>
#include <cmath>
#include <optional>
#include <vector>

namespace rcs {

using namespace Magnum;
using namespace Magnum::Math::Literals;

namespace {

bool rayTriangle(const Vector3& origin, const Vector3& dir,
                const Vector3& v0, const Vector3& v1, const Vector3& v2,
                float& t, float& u, float& v)
{
    const Vector3 edge1 = v1 - v0;
    const Vector3 edge2 = v2 - v0;
    const Vector3 p = Math::cross(dir, edge2);
    const float det = Math::dot(edge1, p);
    if (std::abs(det) < 1.0e-8f)
        return false;
    const float inv = 1.0f / det;
    const Vector3 s = origin - v0;
    u = Math::dot(s, p) * inv;
    if (u < 0.0f || u > 1.0f)
        return false;
    const Vector3 q = Math::cross(s, edge1);
    v = Math::dot(dir, q) * inv;
    if (v < 0.0f || u + v > 1.0f)
        return false;
    t = Math::dot(edge2, q) * inv;
    return t > 1.0e-4f;
}

Vector3 frameNode(const PartFrame& part, int index) {
    if (index < 0)
        return {};
    const Vec3f p = part.position(static_cast<std::size_t>(index));
    return {p.x, p.y, p.z};
}

float nodeKelvin(const PartFrame& part, int index, float fallback) {
    if (index < 0 || static_cast<std::size_t>(index) >= part.nodeTemperature.size())
        return fallback;
    return part.nodeTemperature[static_cast<std::size_t>(index)];
}

/// Closest hit on the visible shells. Temperature is interpolated from the
/// triangle corners. Stress is the quad the point lands on.
bool sampleSurface(const Frame& frame, const RecordingLayout& layout,
                   const Vector3& origin, const Vector3& dir, ShadeMode shade, float& value)
{
    float bestT = 1.0e8f;
    bool hit = false;
    const std::size_t count = std::min(frame.parts.size(), layout.parts.size());
    for (std::size_t i = 0; i < count; ++i) {
        const PartFrame& part = frame.parts[i];
        const PartLayout& partLayout = layout.parts[i];
        const auto& quads = part.quads.empty() ? partLayout.quads : part.quads;
        const auto& triangles = part.triangles.empty() ? partLayout.triangles : part.triangles;
        auto consider = [&](int i0, int i1, int i2, std::size_t fieldIndex) {
            float t = 0.0f, u = 0.0f, v = 0.0f;
            const Vector3 a = frameNode(part, i0);
            const Vector3 b = frameNode(part, i1);
            const Vector3 c = frameNode(part, i2);
            if (!rayTriangle(origin, dir, a, b, c, t, u, v) || t >= bestT)
                return;
            bestT = t;
            hit = true;
            if (shade == ShadeMode::Stress) {
                value = fieldIndex < part.stress.size() ? part.stress[fieldIndex] : 0.0f;
                return;
            }
            const float fallback = fieldIndex < part.temperature.size() ? part.temperature[fieldIndex] : 293.15f;
            const float w = 1.0f - u - v;
            value = w * nodeKelvin(part, i0, fallback) + u * nodeKelvin(part, i1, fallback) +
                    v * nodeKelvin(part, i2, fallback);
        };
        if (!quads.empty()) {
            for (std::size_t q = 0; q < quads.size(); ++q) {
                consider(quads[q][0], quads[q][1], quads[q][2], q);
                consider(quads[q][0], quads[q][2], quads[q][3], q);
            }
        } else {
            for (std::size_t t = 0; t < triangles.size(); ++t)
                consider(triangles[t][0], triangles[t][1], triangles[t][2], t);
        }
    }
    return hit;
}

bool rayAabb(const Vector3& origin, const Vector3& dir, const Range3D& box, float& tHit) {
    float tmin = 0.0f;
    float tmax = 1.0e8f;
    for (int i = 0; i < 3; ++i) {
        if (std::abs(dir[i]) < 1.0e-8f) {
            if (origin[i] < box.min()[i] || origin[i] > box.max()[i])
                return false;
            continue;
        }
        float t1 = (box.min()[i] - origin[i]) / dir[i];
        float t2 = (box.max()[i] - origin[i]) / dir[i];
        if (t1 > t2)
            std::swap(t1, t2);
        tmin = std::max(tmin, t1);
        tmax = std::min(tmax, t2);
        if (tmax < tmin)
            return false;
    }
    tHit = tmin;
    return tHit <= tmax;
}

chrono::ChVector3d toChrono(const Vector3& v) {
    return {double(v.x()), double(v.y()), double(v.z())};
}

Vector3 toMagnum(const chrono::ChVector3d& v) {
    return {Float(v.x()), Float(v.y()), Float(v.z())};
}

float axisClosestParam(const Vector3& rayOrigin,
                       const Vector3& rayDir,
                       const Vector3& axisOrigin,
                       const Vector3& axis)
{
    const Vector3 w0 = rayOrigin - axisOrigin;
    const float a = Math::dot(rayDir, rayDir);
    const float b = Math::dot(rayDir, axis);
    const float c = Math::dot(axis, axis);
    const float d = Math::dot(rayDir, w0);
    const float e = Math::dot(axis, w0);
    const float denom = a * c - b * b;
    if (std::abs(denom) < 1.0e-8f)
        return 0.0f;
    return (a * e - b * d) / denom;
}

}  // namespace

class CombatApp: public Platform::Sdl2Application {
public:
    explicit CombatApp(const Arguments& arguments);

private:
    void drawEvent() override;
    void viewportEvent(ViewportEvent& event) override;
    void pointerPressEvent(PointerEvent& event) override;
    void pointerReleaseEvent(PointerEvent& event) override;
    void pointerMoveEvent(PointerMoveEvent& event) override;
    void scrollEvent(ScrollEvent& event) override;
    void keyPressEvent(KeyEvent& event) override;
    void keyReleaseEvent(KeyEvent& event) override;
    void textInputEvent(TextInputEvent& event) override;

    void rebuildPreview();
    void handleEditorRequests();
    void advancePlayback(float realSeconds);
    void startSimulation();
    const Recording& activeRecording() const;
    const Frame* activeFrame() const;

    Vector2 framebufferPosition(const Vector2& windowPosition) const;
    void rayAt(const Vector2& windowPosition, Vector3& origin, Vector3& direction) const;
    void selectAt(const Vector2& windowPosition);
    void applyTranslation(const Vector3& delta);
    void applyRotation(const Vector3& axis, Rad angle);
    float gizmoScale() const;
    std::vector<std::size_t> selectionGroup() const;
    Vector3 selectionCenter() const;
    bool editingAllowed() const;

    ImGuiIntegration::Context _imgui{NoCreate};
    Timeline _timeline;

    AssemblyDefinition _definition;
    Simulator _preview;
    Recording _pose;
    Recording _run;
    SimulationJob _job;

    OrbitCamera _camera;
    SceneRenderer _renderer;
    TransformGizmos _gizmos;
    EditorState _editor;

    GizmoHandle _hoverHandle = GizmoHandle::None;
    GizmoHandle _activeHandle = GizmoHandle::None;
    Vector3 _dragOrigin;
    Vector3 _dragAxis;
    Vector3 _dragLastPoint;
    float _dragGrabU = 0.0f;
    bool _orbiting = false;
    bool _panning = false;
};

CombatApp::CombatApp(const Arguments& arguments):
    Sdl2Application{arguments, Configuration{}
        .setTitle("Robot Combat Sim")
        .setWindowFlags(Configuration::WindowFlag::Resizable)
        .setSize({1440, 810}),
        GLConfiguration{}.setSampleCount(4)}
{
    GL::Renderer::enable(GL::Renderer::Feature::DepthTest);
    GL::Renderer::enable(GL::Renderer::Feature::FaceCulling);

    _imgui = ImGuiIntegration::Context(Vector2{windowSize()}/dpiScaling(),
                                       windowSize(), framebufferSize());
    _camera.setViewport(framebufferSize());

    _definition = makeTestCube();

    // A sensible default experiment: drop the box on a corner and punch the
    // front panel once it has landed.
    _editor.scenario.initial.position = {0.0, 0.55, 0.0};
    _editor.scenario.initial.orientation = {0.35, 0.0, 0.45};
    _editor.stressScale = float(_definition.parts.front().material.yieldStress);

    rebuildPreview();
    setSwapInterval(1);
    _timeline.start();
}

void CombatApp::rebuildPreview() {
    _preview.reset(_definition, _editor.scenario);
    _pose = _preview.capture();
    // The old run no longer describes this assembly or this starting state.
    _run = Recording{};
    _editor.playback.time = 0.0;
    _editor.playback.playing = false;
    if (_editor.selectedPart && *_editor.selectedPart >= _pose.layout.parts.size())
        _editor.selectedPart.reset();
    if (_editor.selectedGroup && *_editor.selectedGroup >= _definition.groups.size())
        _editor.selectedGroup.reset();
    _renderer.setLayout(_pose.layout);
}

const Recording& CombatApp::activeRecording() const {
    return _run.frames.empty() ? _pose : _run;
}

const Frame* CombatApp::activeFrame() const {
    return activeRecording().sample(_editor.playback.time);
}

bool CombatApp::editingAllowed() const {
    return !_job.running();
}

void CombatApp::startSimulation() {
    _editor.status = "Meshing...";
    AssemblyDefinition meshed = _definition;
    const double resolution = std::clamp(_editor.scenario.settings.meshResolution, 0.25, 8.0);
    for (auto& part : meshed.parts) {
        if (!part.simulate)
            continue;
        if (part.source.empty() && (part.patch.uSize <= 1.0e-6 || part.patch.vSize <= 1.0e-6))
            continue;
        part.sizing.tetBudget = std::clamp(static_cast<int>(2500.0 * resolution), 400, 8000);
        part.sizing.refineLimit = part.sizing.minElement;
        part.sizing.cached = false;
        part.volume = {};
        part.surface = {};
        if (!ensurePartMesh(part))
            _editor.status = "Could not mesh " + part.name + ".";
    }
    _editor.status = "Simulating...";
    _editor.playback.playing = false;
    _editor.playback.time = 0.0;
    _job.start(meshed, _editor.scenario);
}

void CombatApp::handleEditorRequests() {
    if (_editor.addLoadRequested) {
        _editor.addLoadRequested = false;
        ImpactLoad load;
        load.name = "Load " + std::to_string(_editor.scenario.loads.size() + 1);
        load.partIndex = _editor.selectedPart.value_or(0);
        load.initialForce = 0.0;
        load.finalForce = 1000.0;
        load.rampTime = 0.1;
        load.duration = 1.0;
        load.ramp = LoadRamp::Linear;
        load.polynomialDegree = 2;
        load.point = toChrono(selectionGroup().empty()
                                  ? toMagnum(_preview.assembly().center())
                                  : selectionCenter());
        _editor.scenario.loads.push_back(load);
        _editor.selectedLoad = int(_editor.scenario.loads.size()) - 1;
        _editor.rebuildRequested = true;
    }

    if (_editor.addHeatRequested) {
        _editor.addHeatRequested = false;
        HeatSource source;
        source.name = "Heat " + std::to_string(_editor.scenario.heatSources.size() + 1);
        source.point = toChrono(selectionGroup().empty()
                                    ? toMagnum(_preview.assembly().center())
                                    : selectionCenter());
        _editor.scenario.heatSources.push_back(source);
        _editor.selectedHeat = int(_editor.scenario.heatSources.size()) - 1;
    }

    if (_editor.removeHeatIndex >= 0) {
        const int index = _editor.removeHeatIndex;
        _editor.removeHeatIndex = -1;
        if (index < int(_editor.scenario.heatSources.size())) {
            _editor.scenario.heatSources.erase(_editor.scenario.heatSources.begin() + index);
            _editor.selectedHeat = -1;
        }
    }

    if (_editor.removeLoadIndex >= 0) {
        const int index = _editor.removeLoadIndex;
        _editor.removeLoadIndex = -1;
        if (index < int(_editor.scenario.loads.size())) {
            _editor.scenario.loads.erase(_editor.scenario.loads.begin() + index);
            _editor.selectedLoad = -1;
            _editor.rebuildRequested = true;
        }
    }

    if (_editor.rebuildRequested) {
        _editor.rebuildRequested = false;
        if (editingAllowed())
            rebuildPreview();
    }

    if (_editor.cancelRequested) {
        _editor.cancelRequested = false;
        _job.cancel();
    }

    if (_editor.addGroupRequested) {
        _editor.addGroupRequested = false;
        if (editingAllowed()) {
            GroupDefinition folder;
            folder.name = "Folder " + std::to_string(_definition.groups.size() + 1);
            if (_editor.selectedGroup)
                folder.parent = static_cast<int>(*_editor.selectedGroup);
            _definition.groups.push_back(folder);
            _editor.selectedGroup = _definition.groups.size() - 1;
            _editor.selectedPart.reset();
            _editor.status = "Added " + folder.name + ".";
        }
    }

    if (_editor.addSolidCubeRequested) {
        _editor.addSolidCubeRequested = false;
        if (editingAllowed()) {
            _definition.parts.push_back(makeSolidCube(0.40, {0.50, 0.0, 0.0}));
            _editor.selectedPart = _definition.parts.size() - 1;
            rebuildPreview();
            _editor.status = "Added a 400 mm solid aluminum cube.";
        }
    }

    if (_editor.importRequested) {
        _editor.importRequested = false;
        if (editingAllowed()) {
            AssemblyDefinition loaded;
            std::string error;
            if (importAssembly(_editor.importPath, loaded, error)) {
                _definition = std::move(loaded);
                _editor.selectedPart.reset();
                _editor.selectedGroup.reset();
                rebuildPreview();
                _editor.status = error.empty() ? "Imported " + _definition.name : error;
            } else {
                _editor.status = error.empty() ? "Import failed." : error;
            }
        }
    }

    if (_editor.simulateRequested) {
        _editor.simulateRequested = false;
        if (editingAllowed())
            startSimulation();
    }

    Recording finished;
    if (_job.collect(finished)) {
        _run = std::move(finished);
        _editor.status = _run.message;
        _editor.playback.time = 0.0;
        _editor.playback.playing = !_run.frames.empty();
        _renderer.setLayout(_run.layout);
    }
}

void CombatApp::advancePlayback(float realSeconds) {
    if (!_editor.playback.playing || _run.frames.size() < 2)
        return;
    _editor.playback.time += double(realSeconds) * double(_editor.playback.speed);
    if (_editor.playback.time >= _run.duration) {
        if (_editor.playback.loop)
            _editor.playback.time = 0.0;
        else {
            _editor.playback.time = _run.duration;
            _editor.playback.playing = false;
        }
    }
}

Vector2 CombatApp::framebufferPosition(const Vector2& windowPosition) const {
    return windowPosition * Vector2{framebufferSize()} / Vector2{windowSize()};
}

void CombatApp::rayAt(const Vector2& windowPosition, Vector3& origin, Vector3& direction) const {
    _camera.rayFromFramebuffer(framebufferPosition(windowPosition), origin, direction);
}

float CombatApp::gizmoScale() const {
    return std::max(0.12f, _camera.distance() * 0.12f);
}

std::vector<std::size_t> CombatApp::selectionGroup() const {
    std::vector<std::size_t> group;
    const auto& parts = _definition.parts;
    const auto& folders = _definition.groups;
    std::vector<char> folder(_definition.groups.size(), 0);
    if (_editor.selectedGroup && *_editor.selectedGroup < folders.size()) {
        std::vector<int> stack{static_cast<int>(*_editor.selectedGroup)};
        while (!stack.empty()) {
            const int current = stack.back();
            stack.pop_back();
            if (current < 0 || static_cast<std::size_t>(current) >= folder.size() || folder[static_cast<std::size_t>(current)])
                continue;
            folder[static_cast<std::size_t>(current)] = 1;
            for (int child = 0; child < static_cast<int>(folders.size()); ++child)
                if (folders[static_cast<std::size_t>(child)].parent == current)
                    stack.push_back(child);
        }
        for (std::size_t i = 0; i < parts.size(); ++i)
            if (parts[i].group >= 0 && static_cast<std::size_t>(parts[i].group) < folder.size() &&
                folder[static_cast<std::size_t>(parts[i].group)] && parts[i].parentPart < 0)
                group.push_back(i);
    } else if (_editor.selectedPart && *_editor.selectedPart < parts.size()) {
        group.push_back(*_editor.selectedPart);
    }
    for (std::size_t n = 0; n < group.size(); ++n) {
        for (std::size_t i = 0; i < parts.size(); ++i) {
            if (parts[i].parentPart == static_cast<int>(group[n]) &&
                std::find(group.begin(), group.end(), i) == group.end())
                group.push_back(i);
        }
    }
    return group;
}

Vector3 CombatApp::selectionCenter() const {
    return toMagnum(_preview.assembly().centerOf(selectionGroup()));
}

void CombatApp::selectAt(const Vector2& windowPosition) {
    Vector3 origin, dir;
    rayAt(windowPosition, origin, dir);

    std::optional<std::size_t> hit;
    float bestT = 1.0e8f;
    for (std::size_t i = 0; i < _renderer.partCount(); ++i) {
        Range3D bounds;
        if (!_renderer.partBounds(i, bounds))
            continue;
        float t = 0.0f;
        if (rayAabb(origin, dir, bounds, t) && t < bestT) {
            bestT = t;
            hit = i;
        }
    }
    _editor.selectedPart = hit;
    if (hit)
        _editor.selectedGroup.reset();
}

void CombatApp::applyTranslation(const Vector3& delta) {
    const auto group = selectionGroup();
    if (group.empty() || delta.length() < 1.0e-8f)
        return;

    const chrono::ChVector3d world = toChrono(delta);

    // place() drops the mean of the authored nodes onto InitialState::position.
    // Shifting every patch by the same vector moves that mean by the same
    // vector, and the drop cancels it, so dragging the whole welded assembly
    // has to edit the drop position itself.
    if (group.size() == _definition.parts.size()) {
        _editor.scenario.initial.position += world;
        rebuildPreview();
        return;
    }

    const Assembly& assembly = _preview.assembly();
    const chrono::ChVector3d local = assembly.placementRotation().RotateBack(world);

    std::vector<char> inGroup(_definition.parts.size(), 0);
    std::size_t movedNodes = 0;
    std::size_t totalNodes = 0;
    for (const auto& part : assembly.parts())
        totalNodes += part.nodeCount();

    auto shiftMesh = [&](SurfaceMesh& mesh) {
        for (auto& vertex : mesh.vertices)
            vertex += local;
    };
    auto shiftVolume = [&](VolumeMesh& mesh) {
        for (auto& vertex : mesh.vertices)
            vertex += local;
    };
    for (std::size_t index : group) {
        if (index >= _definition.parts.size())
            continue;
        inGroup[index] = 1;
        PartDefinition& part = _definition.parts[index];
        part.patch.origin += local;
        shiftMesh(part.source);
        shiftMesh(part.surface);
        shiftVolume(part.volume);
        if (index < assembly.parts().size())
            movedNodes += assembly.parts()[index].nodeCount();
    }

    // Screws and optional weld regions live in authoring space. A joint only
    // ties parts in the same connected group, so both ends move together.
    for (auto& screw : _definition.screws) {
        if (screw.partA < inGroup.size() && screw.partB < inGroup.size() &&
            inGroup[screw.partA] && inGroup[screw.partB])
            screw.position += local;
    }
    for (auto& weld : _definition.welds) {
        if (!weld.useRegion)
            continue;
        if (weld.partA < inGroup.size() && weld.partB < inGroup.size() &&
            inGroup[weld.partA] && inGroup[weld.partB]) {
            weld.regionMin += local;
            weld.regionMax += local;
        }
    }

    // Recentring onto the old drop position would slide the parts that were
    // not dragged the other way. Shift the drop by the fraction the mean moved.
    if (totalNodes > 0) {
        const double fraction = double(movedNodes) / double(totalNodes);
        _editor.scenario.initial.position += world * fraction;
    }
    rebuildPreview();
}

void CombatApp::applyRotation(const Vector3& axis, Rad angle) {
    if (axis.length() < 1.0e-6f)
        return;
    const auto group = selectionGroup();
    if (group.empty())
        return;

    const Assembly& assembly = _preview.assembly();
    const chrono::ChVector3d localAxis =
        assembly.placementRotation().RotateBack(toChrono(axis.normalized()));
    const chrono::ChQuaterniond rotation =
        chrono::QuatFromAngleAxis(double(float(angle)), localAxis);
    const chrono::ChVector3d worldPivot = toChrono(selectionCenter());
    const chrono::ChVector3d pivot =
        assembly.placementRotation().RotateBack(worldPivot - assembly.placementPosition()) +
        assembly.authoredCenter();

    for (std::size_t index : group) {
        if (index >= _definition.parts.size())
            continue;
        PartDefinition& part = _definition.parts[index];
        auto spin = [&](chrono::ChVector3d& point) {
            point = pivot + rotation.Rotate(point - pivot);
        };
        spin(part.patch.origin);
        part.patch.uAxis = rotation.Rotate(part.patch.uAxis);
        part.patch.vAxis = rotation.Rotate(part.patch.vAxis);
        for (auto* mesh : {&part.source, &part.surface})
            for (auto& vertex : mesh->vertices)
                spin(vertex);
        for (auto& vertex : part.volume.vertices)
            spin(vertex);
    }
    rebuildPreview();
}

void CombatApp::drawEvent() {
    _timeline.nextFrame();

    handleEditorRequests();
    advancePlayback(_timeline.previousFrameDuration());

    const Recording& recording = activeRecording();
    const Frame* frame = activeFrame();

    if (_editor.autoStressScale)
        _editor.stressScale = std::max(recording.maxStress, 1.0e6f);

    if (frame) {
        _renderer.update(recording.layout, *frame, _editor.shade, _editor.stressScale,
                         selectionGroup(), _editor.showAir);
    }

    GL::Renderer::setClearColor(0x1b1e23_rgbf);
    GL::defaultFramebuffer.clear(GL::FramebufferClear::Color|GL::FramebufferClear::Depth);

    const Matrix4 view = _camera.viewMatrix();
    const Matrix4 projection = _camera.projectionMatrix();
    _renderer.draw(view, projection, _editor.scenario.environment, _editor.scenario.heatSources,
                   _editor.playback.time);

    const bool showGizmo = editingAllowed() && _run.frames.empty() &&
                           !selectionGroup().empty() && _editor.tool != Tool::Select;
    if (showGizmo) {
        _gizmos.draw(view, projection, selectionCenter(), gizmoScale(),
                     _editor.tool, _activeHandle != GizmoHandle::None ? _activeHandle : _hoverHandle);
    }

    _imgui.newFrame();
    if (ImGui::GetIO().WantTextInput && !isTextInputActive())
        startTextInput();
    else if (!ImGui::GetIO().WantTextInput && isTextInputActive())
        stopTextInput();

    UiContext context;
    context.definition = &_definition;
    context.layout = &recording.layout;
    context.recording = &recording;
    context.frame = frame;
    context.simulating = _job.running();
    context.progress = _job.progress();
    context.assemblyMass = _preview.assembly().mass();
    context.loadAdvice.reserve(_editor.scenario.loads.size());
    for (std::size_t i = 0; i < _editor.scenario.loads.size(); ++i) {
        context.loadAdvice.push_back(_preview.assembly().adviseLoad(
            i, _editor.scenario.loads, _editor.scenario.settings.timestep));
    }
    drawEditorUi(_editor, context);
    drawArenaOverlay(_editor.scenario.environment, view, projection);
    _imgui.updateApplicationCursor(*this);

    GL::Renderer::enable(GL::Renderer::Feature::Blending);
    GL::Renderer::enable(GL::Renderer::Feature::ScissorTest);
    GL::Renderer::disable(GL::Renderer::Feature::FaceCulling);
    GL::Renderer::disable(GL::Renderer::Feature::DepthTest);
    GL::Renderer::setBlendEquation(GL::Renderer::BlendEquation::Add,
                                   GL::Renderer::BlendEquation::Add);
    GL::Renderer::setBlendFunction(GL::Renderer::BlendFunction::SourceAlpha,
                                   GL::Renderer::BlendFunction::OneMinusSourceAlpha);
    _imgui.drawFrame();
    GL::Renderer::enable(GL::Renderer::Feature::DepthTest);
    GL::Renderer::enable(GL::Renderer::Feature::FaceCulling);
    GL::Renderer::disable(GL::Renderer::Feature::ScissorTest);
    GL::Renderer::disable(GL::Renderer::Feature::Blending);

    swapBuffers();
    redraw();
}

void CombatApp::viewportEvent(ViewportEvent& event) {
    GL::defaultFramebuffer.setViewport({{}, event.framebufferSize()});
    _imgui.relayout(Vector2{event.windowSize()}/event.dpiScaling(),
                    event.windowSize(), event.framebufferSize());
    _camera.setViewport(event.framebufferSize());
}

void CombatApp::pointerPressEvent(PointerEvent& event) {
    if (_imgui.handlePointerPressEvent(event)) {
        event.setAccepted();
        return;
    }
    if (!event.isPrimary())
        return;

    if (event.pointer() == Pointer::MouseRight) {
        _orbiting = true;
        event.setAccepted();
        return;
    }
    if (event.pointer() == Pointer::MouseMiddle) {
        _panning = true;
        event.setAccepted();
        return;
    }
    if (event.pointer() != Pointer::MouseLeft)
        return;

    const bool canEdit = editingAllowed() && _run.frames.empty();
    if (canEdit && !selectionGroup().empty() && _editor.tool != Tool::Select) {
        Vector3 origin, dir;
        rayAt(event.position(), origin, dir);
        const Vector3 center = selectionCenter();
        const float scale = gizmoScale();
        _activeHandle = (_editor.tool == Tool::Move)
            ? _gizmos.pickMove(origin, dir, center, scale)
            : _gizmos.pickRotate(origin, dir, center, scale);
        if (_activeHandle != GizmoHandle::None) {
            _dragOrigin = center;
            _dragAxis = TransformGizmos::axisVector(_activeHandle);
            if (_editor.tool == Tool::Move) {
                _dragGrabU = axisClosestParam(origin, dir, center, _dragAxis);
            } else {
                const float denom = Math::dot(dir, _dragAxis);
                float t = 0.0f;
                if (std::abs(denom) > 1.0e-5f)
                    t = Math::dot(_dragOrigin - origin, _dragAxis) / denom;
                _dragLastPoint = origin + dir * t;
            }
            event.setAccepted();
            return;
        }
    }

    selectAt(event.position());
    event.setAccepted();
}

void CombatApp::pointerReleaseEvent(PointerEvent& event) {
    if (_imgui.handlePointerReleaseEvent(event)) {
        event.setAccepted();
        return;
    }
    if (event.pointer() == Pointer::MouseRight)
        _orbiting = false;
    if (event.pointer() == Pointer::MouseMiddle)
        _panning = false;
    if (event.pointer() == Pointer::MouseLeft)
        _activeHandle = GizmoHandle::None;
    event.setAccepted();
}

void CombatApp::pointerMoveEvent(PointerMoveEvent& event) {
    if (_imgui.handlePointerMoveEvent(event)) {
        _editor.probeVisible = false;
        event.setAccepted();
        return;
    }

    Vector3 origin, dir;
    rayAt(event.position(), origin, dir);
    _editor.probeVisible = false;
    if (_editor.shade == ShadeMode::Temperature || _editor.shade == ShadeMode::Stress) {
        if (const Frame* frame = activeFrame()) {
            _editor.probeVisible = sampleSurface(*frame, activeRecording().layout, origin, dir, _editor.shade,
                                                 _editor.probeValue);
        }
    }

    if (_orbiting && (event.pointers() & Pointer::MouseRight)) {
        _camera.orbit(event.relativePosition());
        event.setAccepted();
        return;
    }
    if (_panning && (event.pointers() & Pointer::MouseMiddle)) {
        _camera.pan(event.relativePosition());
        event.setAccepted();
        return;
    }

    const bool canEdit = editingAllowed() && _run.frames.empty();
    if (canEdit && !selectionGroup().empty() && _editor.tool != Tool::Select &&
        _activeHandle == GizmoHandle::None) {
        const Vector3 center = selectionCenter();
        _hoverHandle = (_editor.tool == Tool::Move)
            ? _gizmos.pickMove(origin, dir, center, gizmoScale())
            : _gizmos.pickRotate(origin, dir, center, gizmoScale());
    }

    if (canEdit && _activeHandle != GizmoHandle::None && (event.pointers() & Pointer::MouseLeft)) {
        if (_editor.tool == Tool::Move) {
            const Vector3 center = selectionCenter();
            const float u = axisClosestParam(origin, dir, center, _dragAxis);
            applyTranslation(_dragAxis * (u - _dragGrabU));
        } else if (_editor.tool == Tool::Rotate) {
            const float denom = Math::dot(dir, _dragAxis);
            if (std::abs(denom) > 1.0e-5f) {
                const float t = Math::dot(_dragOrigin - origin, _dragAxis) / denom;
                const Vector3 hit = origin + dir * t;
                Vector3 from = _dragLastPoint - _dragOrigin;
                Vector3 to = hit - _dragOrigin;
                from -= _dragAxis * Math::dot(from, _dragAxis);
                to -= _dragAxis * Math::dot(to, _dragAxis);
                if (from.length() > 1.0e-5f && to.length() > 1.0e-5f) {
                    from = from.normalized();
                    to = to.normalized();
                    const float sinA = Math::dot(_dragAxis, Math::cross(from, to));
                    const float cosA = Math::dot(from, to);
                    applyRotation(_dragAxis, Rad{std::atan2(sinA, cosA)});
                    _dragLastPoint = hit;
                    _dragOrigin = selectionCenter();
                }
            }
        }
        event.setAccepted();
        return;
    }

    event.setAccepted();
}

void CombatApp::scrollEvent(ScrollEvent& event) {
    if (_imgui.handleScrollEvent(event)) {
        event.setAccepted();
        return;
    }
    _camera.zoom(event.offset().y());
    event.setAccepted();
}

void CombatApp::keyPressEvent(KeyEvent& event) {
    if (_imgui.handleKeyPressEvent(event)) {
        event.setAccepted();
        return;
    }
    if (event.key() == Key::Space && !_run.frames.empty()) {
        _editor.playback.playing = !_editor.playback.playing;
        event.setAccepted();
    }
}

void CombatApp::keyReleaseEvent(KeyEvent& event) {
    if (_imgui.handleKeyReleaseEvent(event))
        event.setAccepted();
}

void CombatApp::textInputEvent(TextInputEvent& event) {
    if (_imgui.handleTextInputEvent(event))
        event.setAccepted();
}

}  // namespace rcs

MAGNUM_SDL2APPLICATION_MAIN(rcs::CombatApp)
