#pragma once

#include "editor/tools.hpp"
#include "physics_core/assembly.hpp"
#include "physics_core/recording.hpp"
#include "physics_core/scenario.hpp"
#include "render_core/renderer.hpp"

#include <cstddef>
#include <optional>
#include <string>

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
    float stressScale = 2.76e8f;  ///< Pa mapped to the top of the colour ramp

    std::optional<std::size_t> selectedPart;
    int selectedLoad = -1;

    bool simulateRequested = false;
    bool cancelRequested = false;
    bool rebuildRequested = false;
    bool addLoadRequested = false;
    int removeLoadIndex = -1;
    bool importRequested = false;
    std::string importPath;

    std::string status;
};

/// View of everything the panels display. The assembly definition is mutable
/// so a selected part's material can be edited in place.
struct UiContext {
    AssemblyDefinition* definition = nullptr;
    const RecordingLayout* layout = nullptr;
    const Recording* recording = nullptr;
    const Frame* frame = nullptr;
    bool simulating = false;
    float progress = 0.0f;
    double assemblyMass = 0.0;
    /// One entry per scenario load, in the same order.
    std::vector<Assembly::LoadAdvisory> loadAdvice;
};

void drawEditorUi(EditorState& state, const UiContext& context);

}  // namespace rcs
