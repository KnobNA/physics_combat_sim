#include "editor/ui.hpp"
#include "physics_core/material.hpp"

#include <Magnum/Math/Color.h>
#include <imgui.h>

#include <algorithm>
#include <cmath>
#include <cstdio>

namespace rcs {
namespace {

bool vectorControl(const char* label, chrono::ChVector3d& value, float speed,
                   const char* format = "%.3f")
{
    float temp[3] = {float(value.x()), float(value.y()), float(value.z())};
    if (!ImGui::DragFloat3(label, temp, speed, 0.0f, 0.0f, format))
        return false;
    value = {temp[0], temp[1], temp[2]};
    return true;
}

bool scalarControl(const char* label, double& value, float speed, float low, float high,
                   const char* format = "%.4f")
{
    float temp = float(value);
    if (!ImGui::DragFloat(label, &temp, speed, low, high, format))
        return false;
    value = temp;
    return true;
}

const char* presetName(void* data, int index) {
    const auto* list = static_cast<const std::vector<Material>*>(data);
    if (!list || index < 0 || index >= static_cast<int>(list->size()))
        return "";
    return (*list)[static_cast<std::size_t>(index)].name.c_str();
}

std::string formatStress(float pascals) {
    char buffer[64];
    if (pascals >= 1.0e9f)
        std::snprintf(buffer, sizeof(buffer), "%.2f GPa", pascals * 1.0e-9f);
    else if (pascals >= 1.0e6f)
        std::snprintf(buffer, sizeof(buffer), "%.1f MPa", pascals * 1.0e-6f);
    else if (pascals >= 1.0e3f)
        std::snprintf(buffer, sizeof(buffer), "%.1f kPa", pascals * 1.0e-3f);
    else
        std::snprintf(buffer, sizeof(buffer), "%.0f Pa", pascals);
    return buffer;
}

void drawStressLegend(float scale) {
    ImDrawList* draw = ImGui::GetWindowDrawList();
    const ImVec2 origin = ImGui::GetCursorScreenPos();
    const float width = ImGui::GetContentRegionAvail().x;
    const float height = 14.0f;
    constexpr int kSegments = 48;

    for (int i = 0; i < kSegments; ++i) {
        const float t0 = float(i) / kSegments;
        const float t1 = float(i + 1) / kSegments;
        const Magnum::Color3 c = SceneRenderer::stressColor(t0);
        const ImVec4 rgba{c.r(), c.g(), c.b(), 1.0f};
        draw->AddRectFilled({origin.x + width * t0, origin.y},
                            {origin.x + width * t1, origin.y + height},
                            ImGui::GetColorU32(rgba));
    }
    ImGui::Dummy({width, height + 2.0f});
    ImGui::TextDisabled("0");
    ImGui::SameLine(width * 0.5f - 20.0f);
    ImGui::TextDisabled("%s", formatStress(scale * 0.5f).c_str());
    ImGui::SameLine(width - 60.0f);
    ImGui::TextDisabled("%s", formatStress(scale).c_str());
}

void drawAssemblyPanel(EditorState& state, const UiContext& context) {
    if (!ImGui::CollapsingHeader("Assembly", ImGuiTreeNodeFlags_DefaultOpen))
        return;
    if (!context.definition || !context.layout)
        return;

    ImGui::Text("%s", context.definition->name.c_str());
    ImGui::TextDisabled("%zu parts, %zu joints, %.2f kg",
                        context.layout->parts.size(), context.layout->joints.size(),
                        context.assemblyMass);

    if (ImGui::BeginListBox("##parts", {-1.0f, 110.0f})) {
        for (std::size_t i = 0; i < context.layout->parts.size(); ++i) {
            const PartLayout& part = context.layout->parts[i];
            const bool selected = state.selectedPart && *state.selectedPart == i;
            char label[128];
            std::snprintf(label, sizeof(label), "%s  (%s, %.1f mm)", part.name.c_str(),
                          part.materialName.c_str(), part.thickness * 1000.0f);
            if (ImGui::Selectable(label, selected))
                state.selectedPart = i;
        }
        ImGui::EndListBox();
    }

    if (state.selectedPart && context.definition &&
        *state.selectedPart < context.definition->parts.size()) {
        PartDefinition& part = context.definition->parts[*state.selectedPart];
        const auto presets = materialPresets();
        int preset = presetIndex(part.material);
        if (ImGui::Combo("Material##preset", &preset, presetName,
                         const_cast<std::vector<Material>*>(&presets),
                         static_cast<int>(presets.size()))) {
            part.material = presets[static_cast<std::size_t>(preset)];
            state.rebuildRequested = true;
        }
        Material& material = part.material;
        bool edited = false;
        int mode = static_cast<int>(material.fracture);
        const char* modes[] = {"Ductile", "Brittle", "Compliant"};
        if (ImGui::Combo("Fracture", &mode, modes, 3)) {
            material.fracture = static_cast<FractureMode>(mode);
            edited = true;
        }
        edited |= ImGui::Checkbox("Dice on first crack", &material.shatterToDust);

        if (ImGui::TreeNode("Elastic")) {
            edited |= scalarControl("Density (kg/m3)", material.density, 10.0f, 1.0f, 2.0e4f, "%.0f");
            edited |= scalarControl("Young's modulus", material.youngsModulus, 1.0e8f, 1.0e5f, 1.0e12f, "%.4g");
            edited |= scalarControl("Poisson", material.poissonRatio, 0.005f, 0.0f, 0.49f, "%.3f");
            edited |= scalarControl("Rayleigh beta", material.dampingBeta, 1.0e-5f, 0.0f, 0.05f, "%.5f");
            ImGui::TreePop();
        }
        if (ImGui::TreeNode("Failure")) {
            edited |= scalarControl("Yield (Pa)", material.yieldStress, 1.0e6f, 0.0f, 5.0e9f, "%.4g");
            edited |= scalarControl("Ultimate (Pa)", material.ultimateStress, 1.0e6f, 0.0f, 5.0e9f, "%.4g");
            edited |= scalarControl("Weld shear (Pa)", material.weldShearStrength, 1.0e6f, 0.0f, 5.0e9f, "%.4g");
            edited |= scalarControl("Tear stretch", material.tearStretch, 0.01f, 1.0f, 8.0f, "%.3f");
            edited |= scalarControl("Plastic fail", material.failPlasticStrain, 0.005f, 0.0f, 2.0f, "%.3f");
            edited |= scalarControl("Fracture strength", material.fractureStrength, 1.0e6f, 0.0f, 5.0e9f, "%.4g");
            edited |= scalarControl("Toughness Kic", material.fractureToughness, 1.0e4f, 0.0f, 1.0e8f, "%.4g");
            edited |= scalarControl("Min shard (m)", material.minimumShardSize, 0.001f, 0.005f, 0.25f, "%.3f");
            ImGui::TextDisabled("Pieces smaller than this become rigid bodies.");
            ImGui::TreePop();
        }
        if (ImGui::TreeNode("Thermal")) {
            edited |= scalarControl("Expansion (1/K)", material.thermalExpansion, 1.0e-7f, 0.0f, 5.0e-4f, "%.3g");
            edited |= scalarControl("Specific heat", material.specificHeat, 5.0f, 0.0f, 5000.0f, "%.0f");
            edited |= scalarControl("Conductivity", material.conductivity, 0.5f, 0.0f, 500.0f, "%.2f");
            ImGui::TextDisabled("Temperature stays at the reference until heat is simulated.");
            ImGui::TreePop();
        }
        if (edited) {
            material.allowPlasticity =
                material.fracture == FractureMode::Ductile && material.yieldStress > 0.0;
            state.rebuildRequested = true;
        }
        ImGui::TextDisabled("Contact stiffness and floor restitution are under Scenario.");
        if (context.frame && *state.selectedPart < context.frame->parts.size()) {
            const auto& strains = context.frame->parts[*state.selectedPart].plasticStrain;
            float peak = 0.0f;
            for (float value : strains)
                peak = std::max(peak, value);
            if (peak > 0.0f)
                ImGui::Text("Plastic strain %.4f", peak);
        }
    }

    if (ImGui::TreeNode("Joints")) {
        for (std::size_t i = 0; i < context.layout->joints.size(); ++i) {
            const JointLayout& joint = context.layout->joints[i];
            int broken = 0;
            float peak = 0.0f;
            if (context.frame && i < context.frame->joints.size()) {
                const JointFrame& jf = context.frame->joints[i];
                for (std::size_t l = 0; l < jf.broken.size(); ++l) {
                    broken += jf.broken[l] ? 1 : 0;
                    peak = std::max(peak, jf.utilisation[l]);
                }
            }
            const ImVec4 colour = broken > 0 ? ImVec4{0.9f, 0.35f, 0.35f, 1.0f}
                                             : ImVec4{0.75f, 0.78f, 0.8f, 1.0f};
            ImGui::TextColored(colour, "%s  %d/%d links, peak %.0f%%", joint.label.c_str(),
                               joint.linkCount - broken, joint.linkCount, peak * 100.0f);
        }
        ImGui::TreePop();
    }

    char path[512] = {};
    std::snprintf(path, sizeof(path), "%s", state.importPath.c_str());
    if (ImGui::InputText("OBJ / STL", path, sizeof(path)))
        state.importPath = path;
    ImGui::BeginDisabled(state.importPath.empty() || context.simulating);
    if (ImGui::Button("Import surface", {-1.0f, 0.0f}))
        state.importRequested = true;
    ImGui::EndDisabled();
    ImGui::TextDisabled("One part per object. Welds stay manual. Solids are skipped.");
}

void drawScenarioPanel(EditorState& state) {
    if (!ImGui::CollapsingHeader("Scenario", ImGuiTreeNodeFlags_DefaultOpen))
        return;

    SimulationSettings& settings = state.scenario.settings;
    bool dirty = false;

    scalarControl("Duration (s)##run", settings.duration, 0.05f, 0.05f, 60.0f, "%.2f");
    scalarControl("Timestep (s)", settings.timestep, 1.0e-5f, 1.0e-5f, 1.0e-2f, "%.5f");
    scalarControl("Capture (fps)", settings.captureRate, 1.0f, 10.0f, 1000.0f, "%.0f");
    ImGui::TextDisabled("%.0f steps, %.0f frames",
                        settings.duration / std::max(1.0e-5, settings.timestep),
                        settings.duration * settings.captureRate);

    ImGui::SeparatorText("Initial state");
    InitialState& initial = state.scenario.initial;
    dirty |= vectorControl("Position (m)", initial.position, 0.01f);
    dirty |= vectorControl("Rotation (rad)", initial.orientation, 0.01f);
    dirty |= vectorControl("Velocity (m/s)", initial.linearVelocity, 0.05f);
    dirty |= vectorControl("Spin (rad/s)", initial.angularVelocity, 0.05f);

    ImGui::SeparatorText("Environment");
    Environment& environment = state.scenario.environment;
    dirty |= vectorControl("Gravity (m/s2)", environment.gravity, 0.05f, "%.2f");
    ImGui::Checkbox("Aerodynamic drag", &environment.aerodynamicsEnabled);
    if (environment.aerodynamicsEnabled) {
        scalarControl("Air density", environment.airDensity, 0.005f, 0.0f, 10.0f, "%.3f");
        scalarControl("Drag coefficient", environment.dragCoefficient, 0.01f, 0.0f, 3.0f, "%.2f");
        ImGui::TextDisabled("1.225 kg/m3 is sea level standard.");
    }
    dirty |= ImGui::Checkbox("Floor", &environment.floorEnabled);
    if (environment.floorEnabled)
        dirty |= scalarControl("Floor height (m)", environment.floorHeight, 0.01f, -10.0f, 10.0f, "%.3f");

    if (ImGui::TreeNode("Contact")) {
        scalarControl("Stiffness (N/m)", environment.contactStiffness, 1.0e4f, 1.0e3f, 1.0e9f, "%.3g");
        scalarControl("Damping (Ns/m)", environment.contactDamping, 10.0f, 0.0f, 1.0e6f, "%.3g");
        scalarControl("Floor restitution", environment.floorRestitution, 0.01f, 0.0f, 1.0f, "%.2f");
        scalarControl("Friction", environment.contactFriction, 0.01f, 0.0f, 2.0f, "%.2f");
        ImGui::TextDisabled("Concrete is about 0.15. 1.00 gives the impact energy back.");
        ImGui::TreePop();
    }

    if (dirty)
        state.rebuildRequested = true;
}

void drawLoadsPanel(EditorState& state, const UiContext& context) {
    if (!ImGui::CollapsingHeader("Applied loads", ImGuiTreeNodeFlags_DefaultOpen))
        return;

    auto& loads = state.scenario.loads;

    if (ImGui::Button("Add load"))
        state.addLoadRequested = true;
    ImGui::SameLine();
    ImGui::BeginDisabled(state.selectedLoad < 0);
    if (ImGui::Button("Remove"))
        state.removeLoadIndex = state.selectedLoad;
    ImGui::EndDisabled();

    if (ImGui::BeginListBox("##loads", {-1.0f, 80.0f})) {
        for (int i = 0; i < int(loads.size()); ++i) {
            char label[160];
            std::snprintf(label, sizeof(label), "%s  %.0f N @ %.3f s%s", loads[i].name.c_str(),
                          loads[i].force, loads[i].startTime,
                          loads[i].enabled ? "" : "  (off)");
            if (ImGui::Selectable(label, state.selectedLoad == i))
                state.selectedLoad = i;
        }
        ImGui::EndListBox();
    }

    if (state.selectedLoad < 0 || state.selectedLoad >= int(loads.size())) {
        ImGui::TextDisabled("Select a load to edit it.");
        return;
    }

    ImpactLoad& load = loads[std::size_t(state.selectedLoad)];
    char nameBuffer[64];
    std::snprintf(nameBuffer, sizeof(nameBuffer), "%s", load.name.c_str());
    if (ImGui::InputText("Name", nameBuffer, sizeof(nameBuffer)))
        load.name = nameBuffer;

    ImGui::Checkbox("Enabled", &load.enabled);

    int partIndex = int(load.partIndex);
    const std::size_t partCount = context.layout ? context.layout->parts.size() : 0;
    if (ImGui::BeginCombo("Target part",
                          partCount && partIndex < int(partCount)
                              ? context.layout->parts[std::size_t(partIndex)].name.c_str()
                              : "-")) {
        for (std::size_t i = 0; i < partCount; ++i) {
            if (ImGui::Selectable(context.layout->parts[i].name.c_str(), partIndex == int(i)))
                load.partIndex = i;
        }
        ImGui::EndCombo();
    }

    vectorControl("Hit point (m)", load.point, 0.005f);
    vectorControl("Direction", load.direction, 0.01f, "%.2f");
    scalarControl("Force (N)", load.force, 50.0f, 0.0f, 1.0e7f, "%.0f");
    scalarControl("Patch radius (m)", load.radius, 0.002f, 0.001f, 1.0f, "%.3f");
    scalarControl("Start (s)", load.startTime, 0.005f, 0.0f, 60.0f, "%.3f");
    scalarControl("Duration (s)##load", load.duration, 0.001f, 0.0001f, 10.0f, "%.4f");
    ImGui::TextDisabled("Impulse %.1f N*s", load.force * load.duration);

    if (std::size_t(state.selectedLoad) < context.loadAdvice.size()) {
        const auto& advice = context.loadAdvice[std::size_t(state.selectedLoad)];
        ImGui::TextDisabled("Spread over %zu nodes, %.0f N on the worst one",
                            advice.nodeCount, advice.peakNodalForce);
        if (advice.deltaVPerStep > 1.0) {
            ImGui::TextColored({0.95f, 0.6f, 0.25f, 1.0f},
                               "Hit is sharper than the timestep can resolve: %.1f m/s per "
                               "step. Shorten the timestep or widen the patch.",
                               advice.deltaVPerStep);
        } else {
            ImGui::TextDisabled("%.2f m/s per step, well inside the timestep.",
                                advice.deltaVPerStep);
        }
    }
}

void drawRunPanel(EditorState& state, const UiContext& context) {
    ImGui::SeparatorText("Run");

    if (context.simulating) {
        ImGui::ProgressBar(context.progress, {-1.0f, 0.0f});
        if (ImGui::Button("Cancel", {-1.0f, 32.0f}))
            state.cancelRequested = true;
    } else {
        if (ImGui::Button("Simulate", {-1.0f, 32.0f}))
            state.simulateRequested = true;
    }

    if (!state.status.empty())
        ImGui::TextWrapped("%s", state.status.c_str());

    scalarControl("Tolerance", state.scenario.settings.solverTolerance, 1.0e-5f, 1.0e-12f, 1.0e6f,
                  "%.4g");
    ImGui::TextDisabled("Looser lets a very large load finish instead of aborting the step.");
}

void drawPlaybackWindow(EditorState& state, const UiContext& context) {
    const ImGuiViewport* viewport = ImGui::GetMainViewport();
    ImGui::SetNextWindowPos({viewport->WorkPos.x + 12.0f,
                             viewport->WorkPos.y + viewport->WorkSize.y - 12.0f},
                            ImGuiCond_Always, {0.0f, 1.0f});
    ImGui::SetNextWindowSize({viewport->WorkSize.x - 24.0f, 0.0f}, ImGuiCond_Always);
    ImGui::Begin("Playback", nullptr,
                 ImGuiWindowFlags_NoResize | ImGuiWindowFlags_NoMove |
                     ImGuiWindowFlags_NoCollapse | ImGuiWindowFlags_AlwaysAutoResize);

    const bool hasRun = context.recording && context.recording->frames.size() > 1;
    ImGui::BeginDisabled(!hasRun);

    if (ImGui::Button(state.playback.playing ? "Pause" : "Play", {90.0f, 0.0f}))
        state.playback.playing = !state.playback.playing;
    ImGui::SameLine();
    if (ImGui::Button("Restart", {90.0f, 0.0f})) {
        state.playback.time = 0.0;
        state.playback.playing = true;
    }
    ImGui::SameLine();
    ImGui::SetNextItemWidth(120.0f);
    ImGui::DragFloat("Speed", &state.playback.speed, 0.05f, 0.05f, 4.0f, "%.2fx");
    ImGui::SameLine();
    ImGui::Checkbox("Loop", &state.playback.loop);

    const double duration = context.recording ? context.recording->duration : 0.0;
    float time = float(state.playback.time);
    ImGui::SetNextItemWidth(-260.0f);
    if (ImGui::SliderFloat("##time", &time, 0.0f, float(std::max(duration, 1.0e-3)), "t = %.4f s")) {
        state.playback.time = time;
        state.playback.playing = false;
    }
    ImGui::SameLine();
    if (context.frame) {
        ImGui::Text("peak %s  |  broken links %d", formatStress(context.frame->maxStress).c_str(),
                    context.frame->brokenLinks);
    }
    ImGui::EndDisabled();

    if (!hasRun)
        ImGui::TextDisabled("Press Simulate to generate a run, then scrub it here.");

    ImGui::End();
}

}  // namespace

void drawEditorUi(EditorState& state, const UiContext& context) {
    ImGui::SetNextWindowPos({12.0f, 12.0f}, ImGuiCond_FirstUseEver);
    ImGui::SetNextWindowSize({330.0f, 720.0f}, ImGuiCond_FirstUseEver);
    ImGui::Begin("Robot Combat Sim");

    drawAssemblyPanel(state, context);
    drawScenarioPanel(state);
    drawLoadsPanel(state, context);
    drawRunPanel(state, context);

    ImGui::SeparatorText("View");
    int shade = state.shade == ShadeMode::Stress ? 1 : 0;
    if (ImGui::RadioButton("Material##shade", shade == 0))
        state.shade = ShadeMode::Material;
    ImGui::SameLine();
    if (ImGui::RadioButton("Stress", shade == 1))
        state.shade = ShadeMode::Stress;

    if (state.shade == ShadeMode::Stress) {
        ImGui::Checkbox("Auto scale", &state.autoStressScale);
        ImGui::BeginDisabled(state.autoStressScale);
        float megapascals = state.stressScale * 1.0e-6f;
        if (ImGui::DragFloat("Range", &megapascals, 1.0f, 0.1f, 5.0e3f, "%.0f MPa"))
            state.stressScale = megapascals * 1.0e6f;
        ImGui::EndDisabled();
        drawStressLegend(state.stressScale);
        if (context.recording)
            ImGui::TextDisabled("Run peak: %s", formatStress(context.recording->maxStress).c_str());
    }

    ImGui::SeparatorText("Tools");
    ImGui::BeginDisabled(context.simulating);
    auto radio = [&](const char* label, Tool tool) {
        if (ImGui::RadioButton(label, state.tool == tool))
            state.tool = tool;
    };
    radio("Select", Tool::Select);
    ImGui::SameLine();
    radio("Move", Tool::Move);
    ImGui::SameLine();
    radio("Rotate", Tool::Rotate);
    ImGui::EndDisabled();
    ImGui::TextDisabled("Right drag orbits, middle drag pans, wheel zooms.");

    ImGui::End();

    drawPlaybackWindow(state, context);
}

}  // namespace rcs
