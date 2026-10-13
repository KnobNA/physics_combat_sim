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
#include <cstddef>
#include <optional>
#include <vector>

#include <imgui.h>

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
    void selectAt(const Vector2& windowPosition, bool additive);
    void applyTranslation(const Vector3& delta);
    void applyRotation(const Vector3& axis, Rad angle);
    float gizmoScale() const;
    std::vector<std::size_t> selectionGroup() const;
    Vector3 selectionCenter() const;
    bool editingAllowed() const;
    void pruneEditorSelection();
    void seedSelectedPartFaces();
    void refreshFaceVisuals();
    void updateSplicePreview(const Vector3& origin, const Vector3& dir);
    bool commitSplice();
    void updateWeldPreview(const Vector3& origin, const Vector3& dir);
    bool commitWeld();
    std::vector<FaceHighlight> gpuFaceHighlights() const;
    chrono::ChQuaterniond worldRotation() const;
    chrono::ChVector3d worldOrigin() const;

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
    Vector2 _lastPointer{};
    bool _havePointer = false;
    bool _spliceHovering = false;
    SplicePreview _authoredSplice;
    bool _weldHovering = false;
    WeldPreview _authoredWeld;
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
    _preview.reset(_definition, _editor.scenario, false);
    _pose = _preview.capture();
    // The old run no longer describes this assembly or this starting state.
    _run = Recording{};
    _editor.playback.time = 0.0;
    _editor.playback.playing = false;
    pruneEditorSelection();
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
        applyResolutionSizing(part, resolution);
    }
    seedWeldBeads(meshed);
    for (auto& part : meshed.parts) {
        if (!part.simulate)
            continue;
        if (part.source.empty() && (part.patch.uSize <= 1.0e-6 || part.patch.vSize <= 1.0e-6))
            continue;
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
        load.partIndex = _editor.selectedParts.empty() ? 0 : _editor.selectedParts.front();
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
        _editor.status = "Stopping...";
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
            _editor.selectedParts.clear();
            _editor.selectedFaces.clear();
            _editor.selectFaces = false;
            _editor.status = "Added " + folder.name + ".";
        }
    }

    if (_editor.addSolidCubeRequested) {
        _editor.addSolidCubeRequested = false;
        if (editingAllowed()) {
            _definition.parts.push_back(makeSolidCube(0.40, {0.50, 0.0, 0.0}));
            _editor.selectedParts = {_definition.parts.size() - 1};
            _editor.selectedGroup.reset();
            _editor.selectedFaces.clear();
            rebuildPreview();
            _editor.status = "Added a 400 mm solid aluminum cube.";
        }
    }

    if (_editor.addSolidSphereRequested) {
        _editor.addSolidSphereRequested = false;
        if (editingAllowed()) {
            _definition.parts.push_back(makeSolidSphere(0.10, {-0.50, 0.0, 0.0}));
            _editor.selectedParts = {_definition.parts.size() - 1};
            _editor.selectedGroup.reset();
            _editor.selectedFaces.clear();
            rebuildPreview();
            _editor.status = "Added a 200 mm solid aluminum sphere.";
        }
    }

    if (_editor.importRequested) {
        _editor.importRequested = false;
        if (editingAllowed()) {
            AssemblyDefinition loaded;
            std::string error;
            if (importAssembly(_editor.importPath, loaded, error)) {
                _definition = std::move(loaded);
                _editor.selectedParts.clear();
                _editor.selectedGroup.reset();
                _editor.selectedFaces.clear();
                _editor.selectFaces = false;
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

    if (_editor.deleteSelectionRequested) {
        _editor.deleteSelectionRequested = false;
        if (editingAllowed() && _run.frames.empty()) {
            if (_editor.selectedWeld && *_editor.selectedWeld < _definition.welds.size()) {
                _definition.welds.erase(_definition.welds.begin() +
                                        static_cast<std::ptrdiff_t>(*_editor.selectedWeld));
                _editor.selectedWeld.reset();
                _editor.status = "Removed weld.";
                rebuildPreview();
            } else if (_editor.selectedGroup && *_editor.selectedGroup < _definition.groups.size()) {
                deleteSelectedGroup(_definition, *_editor.selectedGroup, _editor.selectedParts);
                _editor.selectedGroup.reset();
                _editor.selectedFaces.clear();
                _editor.selectFaces = false;
                _editor.splicePreview = {};
                _editor.weldSource.reset();
                _editor.weldPreview = {};
                rebuildPreview();
            } else {
                deleteSelectedParts(_definition, _editor.selectedParts, _editor.scenario.loads);
                _editor.selectedFaces.clear();
                _editor.selectFaces = false;
                _editor.splicePreview = {};
                _editor.weldSource.reset();
                _editor.weldPreview = {};
                rebuildPreview();
            }
        }
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

void CombatApp::pruneEditorSelection() {
    _editor.selectedParts.erase(
        std::remove_if(_editor.selectedParts.begin(), _editor.selectedParts.end(),
                       [&](std::size_t i) { return i >= _definition.parts.size(); }),
        _editor.selectedParts.end());
    if (_editor.selectedGroup && *_editor.selectedGroup >= _definition.groups.size())
        _editor.selectedGroup.reset();
    _editor.selectedFaces.erase(
        std::remove_if(_editor.selectedFaces.begin(), _editor.selectedFaces.end(),
                       [&](const FaceId& id) {
                           if (id.part >= _definition.parts.size() || !partIsSelected(_editor, id.part))
                               return true;
                           return id.face >= _definition.parts[id.part].source.faces.size();
                       }),
        _editor.selectedFaces.end());
    if (_editor.selectedWeld && *_editor.selectedWeld >= _definition.welds.size())
        _editor.selectedWeld.reset();
    if (_editor.weldSource) {
        const FaceId id = *_editor.weldSource;
        if (id.part >= _definition.parts.size() ||
            id.face >= _definition.parts[id.part].source.faces.size())
            _editor.weldSource.reset();
    }
    if (_editor.selectedParts.empty())
        _editor.selectFaces = false;
    if (_editor.selectedFaces.empty() &&
        (_editor.tool == Tool::Splice || _editor.tool == Tool::Weld))
        _editor.tool = Tool::Select;
    if (_editor.tool != Tool::Weld) {
        _editor.weldPreview = {};
        _editor.weldSource.reset();
    }
}

void CombatApp::seedSelectedPartFaces() {
    for (std::size_t index : _editor.selectedParts) {
        if (index < _definition.parts.size())
            ensureAuthoredFaces(_definition.parts[index]);
    }
}

chrono::ChQuaterniond CombatApp::worldRotation() const {
    return _preview.assembly().placementRotation();
}

chrono::ChVector3d CombatApp::worldOrigin() const {
    return _preview.assembly().placementPosition();
}

std::vector<std::size_t> CombatApp::selectionGroup() const {
    return collectTransformTargets(_definition, _editor.selectedParts, _editor.selectedGroup);
}

Vector3 CombatApp::selectionCenter() const {
    const auto group = selectionGroup();
    if (group.empty())
        return {};
    chrono::ChVector3d sum{0, 0, 0};
    std::size_t count = 0;
    for (std::size_t index : group) {
        if (index >= _definition.parts.size())
            continue;
        sum += partCentroid(_definition.parts[index]);
        ++count;
    }
    if (!count)
        return {};
    return toMagnum(toWorld(sum / double(count), worldRotation(), worldOrigin()));
}

void CombatApp::refreshFaceVisuals() {
    _editor.faceOutlines.clear();
    _editor.weldOutlines.clear();
    const chrono::ChQuaterniond rotation = worldRotation();
    const chrono::ChVector3d origin = worldOrigin();
    for (const FaceId& id : _editor.selectedFaces) {
        if (id.part >= _definition.parts.size())
            continue;
        std::vector<chrono::ChVector3d> loop;
        faceLoopWorld(_definition.parts[id.part], id.face, rotation, origin, loop);
        if (!loop.empty())
            _editor.faceOutlines.push_back(std::move(loop));
    }
    auto pushWeld = [&](std::size_t index) {
        std::vector<std::vector<chrono::ChVector3d>> loops;
        authoredWeldPolylines(_definition, index, loops);
        for (auto& loop : loops) {
            for (auto& point : loop)
                point = toWorld(point, rotation, origin);
            if (loop.size() >= 2)
                _editor.weldOutlines.push_back(std::move(loop));
        }
    };
    if (_editor.selectedWeld && *_editor.selectedWeld < _definition.welds.size()) {
        pushWeld(*_editor.selectedWeld);
    } else if (_editor.selectedFaces.empty()) {
        for (std::size_t i = 0; i < _definition.welds.size(); ++i) {
            const WeldDefinition& weld = _definition.welds[i];
            if (!partIsSelected(_editor, weld.partA) && !partIsSelected(_editor, weld.partB))
                continue;
            pushWeld(i);
        }
    }
}

std::vector<FaceHighlight> CombatApp::gpuFaceHighlights() const {
    std::vector<FaceHighlight> faces;
    const chrono::ChQuaterniond rotation = worldRotation();
    const chrono::ChVector3d origin = worldOrigin();
    for (const FaceId& id : _editor.selectedFaces) {
        if (id.part >= _definition.parts.size())
            continue;
        std::vector<chrono::ChVector3d> loop;
        faceLoopWorld(_definition.parts[id.part], id.face, rotation, origin, loop);
        FaceHighlight highlight;
        highlight.part = id.part;
        highlight.loop.reserve(loop.size());
        for (const auto& point : loop)
            highlight.loop.push_back(toMagnum(point));
        faces.push_back(std::move(highlight));
    }
    return faces;
}

void CombatApp::updateSplicePreview(const Vector3& origin, const Vector3& dir) {
    _editor.splicePreview = {};
    _authoredSplice = {};
    _spliceHovering = false;
    if (_editor.tool != Tool::Splice || _editor.selectedFaces.empty())
        return;
    seedSelectedPartFaces();
    std::vector<std::size_t> parts;
    for (const FaceId& id : _editor.selectedFaces) {
        if (std::find(parts.begin(), parts.end(), id.part) == parts.end())
            parts.push_back(id.part);
    }
    chrono::ChVector3d worldHit{0, 0, 0};
    const auto hit = pickAuthoredFace(_definition, parts, toChrono(origin), toChrono(dir),
                                      worldRotation(), worldOrigin(), &worldHit);
    FaceId target;
    std::optional<chrono::ChVector3d> through;
    if (hit && faceIsSelected(_editor, *hit)) {
        target = *hit;
        _spliceHovering = true;
        through = worldRotation().RotateBack(worldHit - worldOrigin());
    } else if (_editor.selectedFaces.size() == 1) {
        target = _editor.selectedFaces.front();
    } else {
        return;
    }
    if (target.part >= _definition.parts.size())
        return;
    SplicePreview preview =
        makeSplicePreview(_definition.parts[target.part], target.face, _editor.spliceAngle, through);
    if (!preview.valid && through)
        preview = makeSplicePreview(_definition.parts[target.part], target.face, _editor.spliceAngle, {});
    preview.face = target;
    _authoredSplice = preview;
    auto toPlaced = [&](chrono::ChVector3d point) {
        return toWorld(point, worldRotation(), worldOrigin());
    };
    preview.chordA = toPlaced(preview.chordA);
    preview.chordB = toPlaced(preview.chordB);
    for (auto& point : preview.loop)
        point = toPlaced(point);
    preview.hits[0].point = toPlaced(preview.hits[0].point);
    preview.hits[1].point = toPlaced(preview.hits[1].point);
    _editor.splicePreview = std::move(preview);
}

bool CombatApp::commitSplice() {
    if (!_spliceHovering || !_authoredSplice.valid ||
        _authoredSplice.face.part >= _definition.parts.size())
        return false;
    PartDefinition& part = _definition.parts[_authoredSplice.face.part];
    std::size_t faceA = 0;
    std::size_t faceB = 0;
    if (!spliceFace(part, _authoredSplice, faceA, faceB))
        return false;
    const FaceId old = _authoredSplice.face;
    _editor.selectedFaces.erase(
        std::remove(_editor.selectedFaces.begin(), _editor.selectedFaces.end(), old),
        _editor.selectedFaces.end());
    _editor.selectedFaces.push_back(FaceId{old.part, faceA});
    _editor.selectedFaces.push_back(FaceId{old.part, faceB});
    _editor.splicePreview = {};
    _authoredSplice = {};
    _spliceHovering = false;
    rebuildPreview();
    return true;
}

void CombatApp::updateWeldPreview(const Vector3& origin, const Vector3& dir) {
    _editor.weldPreview = {};
    _authoredWeld = {};
    _weldHovering = false;
    if (_editor.tool != Tool::Weld || !_editor.weldSource || _editor.selectedFaces.empty())
        return;
    seedSelectedPartFaces();
    if (_editor.weldSource->part < _definition.parts.size())
        ensureAuthoredFaces(_definition.parts[_editor.weldSource->part]);
    std::vector<std::size_t> destParts;
    for (const FaceId& id : _editor.selectedFaces) {
        if (std::find(destParts.begin(), destParts.end(), id.part) == destParts.end())
            destParts.push_back(id.part);
    }
    chrono::ChVector3d worldHit{0, 0, 0};
    const auto hit = pickAuthoredFace(_definition, destParts, toChrono(origin), toChrono(dir),
                                      worldRotation(), worldOrigin(), &worldHit);
    if (!hit || !faceIsSelected(_editor, *hit))
        return;
    _weldHovering = true;
    const chrono::ChVector3d mouse = worldRotation().RotateBack(worldHit - worldOrigin());
    WeldPreview preview = makeWeldPreview(_definition, _editor.selectedFaces, *_editor.weldSource,
                                          *hit, mouse, _editor.weldAngle, _editor.weldChordAngle);
    _authoredWeld = preview;
    auto toPlaced = [&](chrono::ChVector3d point) {
        return toWorld(point, worldRotation(), worldOrigin());
    };
    preview.mouse = toPlaced(preview.mouse);
    for (auto& point : preview.outline)
        point = toPlaced(point);
    for (auto& dest : preview.destLoops)
        for (auto& point : dest.loop)
            point = toPlaced(point);
    for (auto& hitPt : preview.outlineHits)
        hitPt.point = toPlaced(hitPt.point);
    for (int c = 0; c < 2; ++c) {
        preview.chords[c].a = toPlaced(preview.chords[c].a);
        preview.chords[c].b = toPlaced(preview.chords[c].b);
        preview.chords[c].hits[0].point = toPlaced(preview.chords[c].hits[0].point);
        preview.chords[c].hits[1].point = toPlaced(preview.chords[c].hits[1].point);
    }
    _editor.weldPreview = std::move(preview);
}

bool CombatApp::commitWeld() {
    if (!_weldHovering || !_authoredWeld.valid)
        return false;
    if (!commitWeldPlacement(_definition, _authoredWeld))
        return false;
    _editor.weldSource.reset();
    _editor.weldPreview = {};
    _authoredWeld = {};
    _weldHovering = false;
    _editor.tool = Tool::Select;
    _editor.status = "Welded.";
    rebuildPreview();
    return true;
}

void CombatApp::selectAt(const Vector2& windowPosition, bool additive) {
    Vector3 origin, dir;
    rayAt(windowPosition, origin, dir);

    if (!_editor.selectFaces && _editor.selectedFaces.empty()) {
        const double maxDist = std::max(0.008, double(gizmoScale()) * 0.15);
        std::vector<std::size_t> parts = _editor.selectedParts;
        if (_editor.selectedWeld && *_editor.selectedWeld < _definition.welds.size()) {
            parts.push_back(_definition.welds[*_editor.selectedWeld].partA);
            parts.push_back(_definition.welds[*_editor.selectedWeld].partB);
        }
        const auto weld = pickWeld(_definition, parts, toChrono(origin), toChrono(dir),
                                   worldRotation(), worldOrigin(), maxDist);
        if (weld) {
            _editor.selectedWeld = weld;
            _editor.selectedParts.clear();
            _editor.selectedGroup.reset();
            _editor.selectedFaces.clear();
            _editor.selectFaces = false;
            return;
        }
    }

    if (_editor.selectFaces && !_editor.selectedParts.empty()) {
        seedSelectedPartFaces();
        const auto face = pickAuthoredFace(_definition, _editor.selectedParts, toChrono(origin),
                                           toChrono(dir), worldRotation(), worldOrigin());
        if (face) {
            if (additive) {
                auto found = std::find(_editor.selectedFaces.begin(), _editor.selectedFaces.end(), *face);
                if (found == _editor.selectedFaces.end())
                    _editor.selectedFaces.push_back(*face);
                else
                    _editor.selectedFaces.erase(found);
            } else {
                _editor.selectedFaces = {*face};
            }
            return;
        }
    }

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
    if (!hit) {
        if (!additive) {
            _editor.selectedParts.clear();
            _editor.selectedGroup.reset();
            _editor.selectedFaces.clear();
            _editor.selectFaces = false;
            _editor.selectedWeld.reset();
        }
        return;
    }
    _editor.selectedGroup.reset();
    _editor.selectedWeld.reset();
    if (additive) {
        auto found = std::find(_editor.selectedParts.begin(), _editor.selectedParts.end(), *hit);
        if (found == _editor.selectedParts.end())
            _editor.selectedParts.push_back(*hit);
        else {
            _editor.selectedParts.erase(found);
            _editor.selectedFaces.erase(
                std::remove_if(_editor.selectedFaces.begin(), _editor.selectedFaces.end(),
                               [&](const FaceId& id) { return id.part == *hit; }),
                _editor.selectedFaces.end());
        }
    } else {
        _editor.selectedParts = {*hit};
        _editor.selectedFaces.clear();
    }
    if (_editor.selectedParts.empty())
        _editor.selectFaces = false;
}

void CombatApp::applyTranslation(const Vector3& delta) {
    if (!_editor.selectedFaces.empty())
        return;
    const auto group = selectionGroup();
    if (group.empty() || delta.length() < 1.0e-8f)
        return;

    const chrono::ChVector3d world = toChrono(delta);
    const chrono::ChVector3d local = worldRotation().RotateBack(world);

    std::vector<char> inGroup(_definition.parts.size(), 0);
    for (std::size_t index : group) {
        if (index >= _definition.parts.size())
            continue;
        inGroup[index] = 1;
        applyPartTranslation(_definition.parts[index], local);
    }

    transformAuthoredJoints(_definition, inGroup, local);
    rebuildPreview();
}

void CombatApp::applyRotation(const Vector3& axis, Rad angle) {
    if (!_editor.selectedFaces.empty())
        return;
    if (axis.length() < 1.0e-6f)
        return;
    const auto group = selectionGroup();
    if (group.empty())
        return;

    const chrono::ChVector3d localAxis = worldRotation().RotateBack(toChrono(axis.normalized()));
    const chrono::ChQuaterniond rotation = chrono::QuatFromAngleAxis(double(float(angle)), localAxis);
    const chrono::ChVector3d worldPivot = toChrono(selectionCenter());
    const chrono::ChVector3d pivot = worldRotation().RotateBack(worldPivot - worldOrigin()) +
                                     _preview.assembly().authoredCenter();

    for (std::size_t index : group) {
        if (index >= _definition.parts.size())
            continue;
        applyPartRotation(_definition.parts[index], rotation, pivot);
    }
    std::vector<char> inGroup(_definition.parts.size(), 0);
    for (std::size_t index : group) {
        if (index < inGroup.size())
            inGroup[index] = 1;
    }
    rotateAuthoredJoints(_definition, inGroup, rotation, pivot);
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

    pruneEditorSelection();
    if (editingAllowed() && _run.frames.empty()) {
        Vector3 origin, dir;
        rayAt(_havePointer ? _lastPointer : Vector2{windowSize()} * 0.5f, origin, dir);
        if (_editor.tool == Tool::Splice)
            updateSplicePreview(origin, dir);
        if (_editor.tool == Tool::Weld)
            updateWeldPreview(origin, dir);
    }
    refreshFaceVisuals();

    if (frame) {
        _renderer.update(recording.layout, *frame, _editor.shade, _editor.stressScale,
                         selectionGroup(), gpuFaceHighlights(), _editor.showAir);
    }

    GL::Renderer::setClearColor(0x1b1e23_rgbf);
    GL::defaultFramebuffer.clear(GL::FramebufferClear::Color|GL::FramebufferClear::Depth);

    const Matrix4 view = _camera.viewMatrix();
    const Matrix4 projection = _camera.projectionMatrix();
    _renderer.draw(view, projection, _editor.scenario.environment, _editor.scenario.heatSources,
                   _editor.playback.time);

    const bool showGizmo = editingAllowed() && _run.frames.empty() &&
                           _editor.selectedFaces.empty() &&
                           (_editor.tool == Tool::Move || _editor.tool == Tool::Rotate) &&
                           !selectionGroup().empty();
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
    context.cancelling = _job.cancelling();
    context.progress = _job.progress();
    context.simulatedTime = _job.simulatedTime();
    context.runDuration = _job.runDuration();
    context.simStep = _job.stepIndex();
    context.simStepCount = _job.stepCount();
    context.assemblyMass = _preview.assembly().mass();
    context.loadAdvice.reserve(_editor.scenario.loads.size());
    for (std::size_t i = 0; i < _editor.scenario.loads.size(); ++i) {
        context.loadAdvice.push_back(_preview.assembly().adviseLoad(
            i, _editor.scenario.loads, _editor.scenario.settings.timestep));
    }
    drawEditorUi(_editor, context);
    drawArenaOverlay(_editor.scenario.environment, view, projection);
    drawFaceOverlays(_editor, view, projection);
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

    _lastPointer = event.position();
    _havePointer = true;

    const bool canEdit = editingAllowed() && _run.frames.empty();
    const bool gizmoTool = _editor.selectedFaces.empty() &&
                           (_editor.tool == Tool::Move || _editor.tool == Tool::Rotate);
    if (canEdit && !selectionGroup().empty() && gizmoTool) {
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

    if (canEdit && _editor.tool == Tool::Splice) {
        Vector3 origin, dir;
        rayAt(event.position(), origin, dir);
        updateSplicePreview(origin, dir);
        if (commitSplice()) {
            event.setAccepted();
            return;
        }
    }

    if (canEdit && _editor.tool == Tool::Weld) {
        Vector3 origin, dir;
        rayAt(event.position(), origin, dir);
        if (_editor.weldSource) {
            updateWeldPreview(origin, dir);
            if (commitWeld()) {
                event.setAccepted();
                return;
            }
            event.setAccepted();
            return;
        }
        for (std::size_t i = 0; i < _definition.parts.size(); ++i)
            ensureAuthoredFaces(_definition.parts[i]);
        std::vector<std::size_t> others;
        for (std::size_t i = 0; i < _definition.parts.size(); ++i) {
            bool dest = false;
            for (const FaceId& id : _editor.selectedFaces)
                dest = dest || id.part == i;
            if (!dest)
                others.push_back(i);
        }
        const auto face = pickAuthoredFace(_definition, others, toChrono(origin), toChrono(dir),
                                           worldRotation(), worldOrigin());
        if (face) {
            _editor.weldSource = *face;
            _editor.weldAngle = 0.0;
            _editor.weldChordAngle = 0.0;
            _editor.status = "Place the weld on the selected faces.";
        }
        event.setAccepted();
        return;
    }

    const bool additive = bool(event.modifiers() & Modifier::Ctrl);
    selectAt(event.position(), additive);
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
    _lastPointer = event.position();
    _havePointer = true;
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
    const bool gizmoTool = _editor.selectedFaces.empty() &&
                           (_editor.tool == Tool::Move || _editor.tool == Tool::Rotate);
    if (canEdit && !selectionGroup().empty() && gizmoTool &&
        _activeHandle == GizmoHandle::None) {
        const Vector3 center = selectionCenter();
        _hoverHandle = (_editor.tool == Tool::Move)
            ? _gizmos.pickMove(origin, dir, center, gizmoScale())
            : _gizmos.pickRotate(origin, dir, center, gizmoScale());
    } else if (!gizmoTool) {
        _hoverHandle = GizmoHandle::None;
    }

    if (canEdit && _editor.tool == Tool::Splice)
        updateSplicePreview(origin, dir);
    if (canEdit && _editor.tool == Tool::Weld)
        updateWeldPreview(origin, dir);

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
    if (editingAllowed() && _run.frames.empty() && _editor.tool == Tool::Splice &&
        !_editor.selectedFaces.empty()) {
        Vector3 origin, dir;
        rayAt(event.position(), origin, dir);
        updateSplicePreview(origin, dir);
        if (_editor.splicePreview.valid) {
            _editor.spliceAngle += double(event.offset().y()) * (5.0 * 3.14159265358979323846 / 180.0);
            updateSplicePreview(origin, dir);
            event.setAccepted();
            return;
        }
    }
    if (editingAllowed() && _run.frames.empty() && _editor.tool == Tool::Weld &&
        _editor.weldSource && !_editor.selectedFaces.empty()) {
        Vector3 origin, dir;
        rayAt(event.position(), origin, dir);
        updateWeldPreview(origin, dir);
        if (_editor.weldPreview.valid) {
            const double step = double(event.offset().y()) * (5.0 * 3.14159265358979323846 / 180.0);
            if (event.modifiers() & Modifier::Ctrl)
                _editor.weldChordAngle += step;
            else
                _editor.weldAngle += step;
            updateWeldPreview(origin, dir);
            event.setAccepted();
            return;
        }
    }
    _camera.zoom(event.offset().y());
    event.setAccepted();
}

void CombatApp::keyPressEvent(KeyEvent& event) {
    const bool typing = ImGui::GetIO().WantTextInput;
    if (!typing && event.key() == Key::Delete && editingAllowed() && _run.frames.empty()) {
        if (_editor.selectedWeld || !_editor.selectedParts.empty() || _editor.selectedGroup) {
            _editor.deleteSelectionRequested = true;
            event.setAccepted();
            return;
        }
    }
    if (!typing && editingAllowed() && _run.frames.empty() &&
        event.key() == Key::W && (event.modifiers() & Modifier::Ctrl) &&
        !_editor.selectedFaces.empty()) {
        _editor.tool = Tool::Weld;
        _editor.weldSource.reset();
        _editor.weldPreview = {};
        _editor.status = "Pick a face on another part to weld.";
        event.setAccepted();
        return;
    }
    if (!typing && editingAllowed() && _run.frames.empty() && event.key() == Key::Esc &&
        _editor.tool == Tool::Weld) {
        if (_editor.weldSource) {
            _editor.weldSource.reset();
            _editor.weldPreview = {};
            _editor.status = "Pick a face on another part to weld.";
        } else {
            _editor.tool = Tool::Select;
        }
        event.setAccepted();
        return;
    }
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
