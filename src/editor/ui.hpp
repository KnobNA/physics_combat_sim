#pragma once

#include "editor/faces.hpp"
#include "editor/tools.hpp"
#include "physics_core/assembly.hpp"
#include "physics_core/recording.hpp"
#include "physics_core/scenario.hpp"
#include "render_core/renderer.hpp"

#include <algorithm>
#include <cstddef>
#include <optional>
#include <string>
#include <vector>

namespace rcs {

struct PlaybackState {
    bool playing = false;
    double time = 0.0;
    float speed = 1.0f;
    bool loop = true;
};

/// All editor-owned state. The app treats the `*Requested` members as
/// one-shot commands and clears them once handled.
struct EditorState {
    Tool tool = Tool::Select;
    Scenario scenario;
    PlaybackState playback;

    ShadeMode shade = ShadeMode::Stress;
    bool autoStressScale = true;
    bool showAir = false;
    /// Hover readout for the temperature or stress view. `probeValue` is
    /// kelvin or pascals, matching the active shade mode.
    bool probeVisible = false;
    float probeValue = 0.0f;
    float stressScale = 2.76e8f;  ///< Pa mapped to the top of the colour ramp

    std::vector<std::size_t> selectedParts;
    /// Empty folder selection. Mutually exclusive with a part selection.
    std::optional<std::size_t> selectedGroup;
    bool selectFaces = false;
    std::vector<FaceId> selectedFaces;
    double spliceAngle = 0.0;
    SplicePreview splicePreview;
    std::optional<FaceId> weldSource;
    double weldAngle = 0.0;
    double weldChordAngle = 0.0;
    WeldPreview weldPreview;
    std::optional<std::size_t> selectedWeld;
    std::vector<std::vector<chrono::ChVector3d>> faceOutlines;
    std::vector<std::vector<chrono::ChVector3d>> weldOutlines;

    bool addGroupRequested = false;
    bool deleteSelectionRequested = false;
    int selectedLoad = -1;
    int selectedHeat = -1;

    bool simulateRequested = false;
    bool cancelRequested = false;
    bool rebuildRequested = false;
    bool addLoadRequested = false;
    int removeLoadIndex = -1;
    bool addHeatRequested = false;
    int removeHeatIndex = -1;
    bool importRequested = false;
    bool addSolidCubeRequested = false;
    bool addSolidSphereRequested = false;
    std::string importPath;

    std::string status;
};

inline bool partIsSelected(const EditorState& state, std::size_t index) {
    return std::find(state.selectedParts.begin(), state.selectedParts.end(), index) !=
           state.selectedParts.end();
}

inline bool faceIsSelected(const EditorState& state, FaceId id) {
    return std::find(state.selectedFaces.begin(), state.selectedFaces.end(), id) != state.selectedFaces.end();
}

/// View of everything the panels display. The assembly definition is mutable
/// so a selected part's material can be edited in place.
struct UiContext {
    AssemblyDefinition* definition = nullptr;
    const RecordingLayout* layout = nullptr;
    const Recording* recording = nullptr;
    const Frame* frame = nullptr;
    bool simulating = false;
    bool cancelling = false;
    float progress = 0.0f;
    double simulatedTime = 0.0;
    double runDuration = 0.0;
    long long simStep = 0;
    long long simStepCount = 0;
    double assemblyMass = 0.0;
    /// One entry per scenario load, in the same order.
    std::vector<Assembly::LoadAdvisory> loadAdvice;
};

void drawEditorUi(EditorState& state, const UiContext& context);

/// Metre ticks on the floor grid, and the length of the three axis edges.
void drawArenaOverlay(const Environment& environment,
                      const Magnum::Matrix4& view,
                      const Magnum::Matrix4& projection);

void drawFaceOverlays(const EditorState& state,
                      const Magnum::Matrix4& view,
                      const Magnum::Matrix4& projection);

}  // namespace rcs
