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

void drawTemperatureLegend() {
    ImDrawList* draw = ImGui::GetWindowDrawList();
    const ImVec2 origin = ImGui::GetCursorScreenPos();
    const float width = 16.0f;
    const float height = 120.0f;
    constexpr int kSegments = 32;
    for (int i = 0; i < kSegments; ++i) {
        const float t0 = float(i) / kSegments;
        const float t1 = float(i + 1) / kSegments;
        const Magnum::Color3 c = SceneRenderer::temperatureColor(2000.0f * (1.0f - t0));
        const ImVec4 rgba{c.r(), c.g(), c.b(), 1.0f};
        draw->AddRectFilled({origin.x, origin.y + height * t0},
                            {origin.x + width, origin.y + height * t1},
                            ImGui::GetColorU32(rgba));
    }
    ImGui::Dummy({width, height});
    ImGui::SameLine();
    ImGui::BeginGroup();
    ImGui::TextDisabled("2000 K");
    ImGui::Dummy({1.0f, height * 0.72f});
    ImGui::TextDisabled("293 K");
    ImGui::Dummy({1.0f, height * 0.08f});
    ImGui::TextDisabled("0 K");
    ImGui::EndGroup();
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

    const std::size_t partCount = std::min(context.definition->parts.size(), context.layout->parts.size());
    auto setPartBranch = [&](std::size_t index, bool on, auto& self) -> void {
        if (index >= context.definition->parts.size())
            return;
        context.definition->parts[index].simulate = on;
        for (std::size_t i = 0; i < context.definition->parts.size(); ++i)
            if (context.definition->parts[i].parentPart == static_cast<int>(index))
                self(i, on, self);
    };
    auto setFolder = [&](int groupIndex, bool on, auto& self) -> void {
        if (groupIndex < 0 || static_cast<std::size_t>(groupIndex) >= context.definition->groups.size())
            return;
        context.definition->groups[static_cast<std::size_t>(groupIndex)].simulate = on;
        for (int child = 0; child < static_cast<int>(context.definition->groups.size()); ++child)
            if (context.definition->groups[static_cast<std::size_t>(child)].parent == groupIndex)
                self(child, on, self);
        for (std::size_t i = 0; i < context.definition->parts.size(); ++i)
            if (context.definition->parts[i].group == groupIndex && context.definition->parts[i].parentPart < 0)
                setPartBranch(i, on, setPartBranch);
        state.rebuildRequested = true;
    };
    auto showPart = [&](std::size_t i, auto& drawChildren) -> void {
        const PartDefinition& defined = context.definition->parts[i];
        const PartLayout& layout = context.layout->parts[i];
        char label[160];
        if (defined.simulate && layout.volume)
            std::snprintf(label, sizeof(label), "%s  (%s, solid, %d tets)", defined.name.c_str(),
                          defined.material.name.c_str(), layout.tetCount);
        else if (defined.simulate)
            std::snprintf(label, sizeof(label), "%s  (%s, %.1f mm)", defined.name.c_str(),
                          defined.material.name.c_str(), defined.thickness * 1000.0);
        else
            std::snprintf(label, sizeof(label), "%s  (off)", defined.name.c_str());
        const bool partSelected = state.selectedPart && *state.selectedPart == i;
        ImGui::PushID(static_cast<int>(i) + 1000);
        ImGui::BeginDisabled(context.simulating);
        bool on = context.definition->parts[i].simulate;
        if (ImGui::Checkbox("##sim", &on)) {
            setPartBranch(i, on, setPartBranch);
            state.rebuildRequested = true;
        }
        ImGui::EndDisabled();
        if (ImGui::IsItemHovered())
            ImGui::SetTooltip("Include this part in the view and the simulation");
        ImGui::SameLine();
        if (ImGui::Selectable(label, partSelected)) {
            state.selectedPart = i;
            state.selectedGroup.reset();
        }
        drawChildren(i);
        ImGui::PopID();
    };
    auto drawChildren = [&](std::size_t parent, auto& self) -> void {
        for (std::size_t i = 0; i < partCount; ++i) {
            if (context.definition->parts[i].parentPart != static_cast<int>(parent))
                continue;
            ImGui::Indent();
            showPart(i, [&](std::size_t child) { self(child, self); });
            ImGui::Unindent();
        }
    };
    auto drawGroup = [&](int groupIndex, auto& self) -> void {
        const GroupDefinition& group = context.definition->groups[static_cast<std::size_t>(groupIndex)];
        const bool selected = state.selectedGroup && *state.selectedGroup == static_cast<std::size_t>(groupIndex);
        ImGui::PushID(groupIndex);
        ImGui::BeginDisabled(context.simulating);
        bool on = group.simulate;
        if (ImGui::Checkbox("##sim", &on))
            setFolder(groupIndex, on, setFolder);
        ImGui::EndDisabled();
        if (ImGui::IsItemHovered())
            ImGui::SetTooltip("Turn every child of this folder on or off");
        ImGui::SameLine();
        const ImGuiTreeNodeFlags flags = ImGuiTreeNodeFlags_DefaultOpen | ImGuiTreeNodeFlags_OpenOnArrow |
                                         (selected ? ImGuiTreeNodeFlags_Selected : 0);
        const bool open = ImGui::TreeNodeEx(group.name.c_str(), flags);
        if (ImGui::IsItemClicked() && !ImGui::IsItemToggledOpen()) {
            state.selectedGroup = static_cast<std::size_t>(groupIndex);
            state.selectedPart.reset();
        }
        if (open) {
            for (int child = 0; child < static_cast<int>(context.definition->groups.size()); ++child) {
                if (context.definition->groups[static_cast<std::size_t>(child)].parent == groupIndex)
                    self(child, self);
            }
            for (std::size_t i = 0; i < partCount; ++i) {
                const PartDefinition& part = context.definition->parts[i];
                if (part.group != groupIndex || part.parentPart >= 0)
                    continue;
                showPart(i, [&](std::size_t child) { drawChildren(child, drawChildren); });
            }
            ImGui::TreePop();
        }
        ImGui::PopID();
    };
    for (int g = 0; g < static_cast<int>(context.definition->groups.size()); ++g) {
        if (context.definition->groups[static_cast<std::size_t>(g)].parent < 0)
            drawGroup(g, drawGroup);
    }
    for (std::size_t i = 0; i < partCount; ++i) {
        const PartDefinition& part = context.definition->parts[i];
        if (part.group >= 0 || part.parentPart >= 0)
            continue;
        showPart(i, [&](std::size_t child) { drawChildren(child, drawChildren); });
    }
    ImGui::BeginDisabled(context.simulating);
    if (ImGui::Button("New folder", {-1.0f, 0.0f}))
        state.addGroupRequested = true;
    ImGui::EndDisabled();

    if (state.selectedPart && context.definition &&
        *state.selectedPart < context.definition->parts.size()) {
        PartDefinition& part = context.definition->parts[*state.selectedPart];
        const char* folderLabel = "(none)";
        if (part.group >= 0 && static_cast<std::size_t>(part.group) < context.definition->groups.size())
            folderLabel = context.definition->groups[static_cast<std::size_t>(part.group)].name.c_str();
        if (ImGui::BeginCombo("Folder", folderLabel)) {
            if (ImGui::Selectable("(none)", part.group < 0)) {
                part.group = -1;
                state.rebuildRequested = true;
            }
            for (int g = 0; g < static_cast<int>(context.definition->groups.size()); ++g) {
                const bool chosen = part.group == g;
                if (ImGui::Selectable(context.definition->groups[static_cast<std::size_t>(g)].name.c_str(), chosen)) {
                    part.group = g;
                    part.parentPart = -1;
                    state.rebuildRequested = true;
                }
            }
            ImGui::EndCombo();
        }
        ImGui::TextDisabled("Moving a folder moves its children and leaves everything else where it is.");
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
            edited |= scalarControl("Emissivity", material.emissivity, 0.01f, 0.0f, 1.0f, "%.2f");
            ImGui::TextDisabled("Stiffness and yield follow the handbook curve and keep falling past the solidus.");
            edited |= scalarControl("Initial temperature (K)", part.initialTemperature, 1.0f, 0.0f, 4000.0f, "%.2f");
            ImGui::TextDisabled("The mesh is unstressed at this temperature. Room is 293.15 K.");
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
    ImGui::BeginDisabled(context.simulating);
    if (ImGui::Button("Add solid aluminum cube", {-1.0f, 0.0f}))
        state.addSolidCubeRequested = true;
    ImGui::EndDisabled();
    ImGui::TextDisabled("Thin walls stay shells. Thick regions are tetrahedra. Welds stay manual.");

    if (state.selectedPart && context.definition &&
        *state.selectedPart < context.definition->parts.size()) {
        PartDefinition& sized = context.definition->parts[*state.selectedPart];
        if (!sized.source.empty() || !sized.volume.empty()) {
            if (ImGui::TreeNode("Mesh size")) {
                bool edited = false;
                edited |= ImGui::InputDouble("Shell cutoff (m)", &sized.sizing.shellCutoff, 0.0, 0.0, "%.4f");
                edited |= ImGui::InputDouble("Min element (m)", &sized.sizing.minElement, 0.0, 0.0, "%.4f");
                edited |= ImGui::InputDouble("Max element (m)", &sized.sizing.maxElement, 0.0, 0.0, "%.4f");
                int through = sized.sizing.throughThickness;
                if (ImGui::InputInt("Elements through thickness", &through)) {
                    sized.sizing.throughThickness = std::clamp(through, 1, 8);
                    edited = true;
                }
                if (!sized.volume.empty())
                    ImGui::TextDisabled("%zu tets, element %.1f mm", sized.volume.tets.size(),
                                        sized.volume.elementSize * 1000.0);
                if (!sized.surface.empty())
                    ImGui::TextDisabled("%zu shell triangles", sized.surface.triangles.size());
                if (edited) {
                    sized.sizing.shellCutoff = std::clamp(sized.sizing.shellCutoff, 0.0005, 0.05);
                    sized.sizing.minElement = std::clamp(sized.sizing.minElement, 0.0005, sized.sizing.maxElement);
                    sized.sizing.maxElement = std::clamp(sized.sizing.maxElement, sized.sizing.minElement, 0.1);
                    sized.sizing.cached = false;
                    sized.volume = {};
                    sized.surface = {};
                    state.rebuildRequested = true;
                }
                ImGui::TreePop();
            }
        }
    }
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
    scalarControl("Arena width (m)", environment.boundsSize.x(), 0.1f, 1.0f, 200.0f, "%.2f");
    scalarControl("Arena height (m)", environment.boundsSize.y(), 0.1f, 1.0f, 200.0f, "%.2f");
    scalarControl("Arena depth (m)", environment.boundsSize.z(), 0.1f, 1.0f, 200.0f, "%.2f");
    scalarControl("Ambient (K)", environment.ambientTemperature, 0.5f, 0.0f, 4000.0f, "%.2f");
    ImGui::TextDisabled("White edges mark the arena. Grid lines are %g m. Parts stay inside.",
                        Environment::gridStep(std::max(environment.boundsSize.x(),
                                                       environment.boundsSize.z())));

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
            if (loads[i].rampTime > 1.0e-9) {
                std::snprintf(label, sizeof(label), "%s  %.0f to %.0f N @ %.3f s%s", loads[i].name.c_str(),
                              loads[i].initialForce, loads[i].finalForce, loads[i].startTime,
                              loads[i].enabled ? "" : "  (off)");
            } else {
                std::snprintf(label, sizeof(label), "%s  %.0f N @ %.3f s%s", loads[i].name.c_str(),
                              loads[i].finalForce, loads[i].startTime,
                              loads[i].enabled ? "" : "  (off)");
            }
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
    scalarControl("Initial force (N)", load.initialForce, 50.0f, -1.0e7f, 1.0e7f, "%.0f");
    scalarControl("Final force (N)", load.finalForce, 50.0f, -1.0e7f, 1.0e7f, "%.0f");
    scalarControl("Ramp time (s)", load.rampTime, 0.005f, 0.0f, 30.0f, "%.3f");
    const char* rampNames[] = {"Linear", "Logarithmic", "Exponential", "Polynomial"};
    int ramp = static_cast<int>(load.ramp);
    if (ImGui::Combo("Ramp", &ramp, rampNames, 4))
        load.ramp = static_cast<LoadRamp>(ramp);
    if (load.ramp == LoadRamp::Polynomial)
        ImGui::DragInt("Degree", &load.polynomialDegree, 0.1f, 1, 8);
    scalarControl("Patch radius (m)", load.radius, 0.002f, 0.001f, 1.0f, "%.3f");
    scalarControl("Start (s)", load.startTime, 0.005f, 0.0f, 60.0f, "%.3f");
    scalarControl("Duration (s)##load", load.duration, 0.001f, 0.0001f, 30.0f, "%.4f");
    double impulse = 0.0;
    const int samples = 48;
    for (int i = 0; i < samples; ++i) {
        const double t = load.startTime + (static_cast<double>(i) + 0.5) / samples * load.duration;
        impulse += loadMagnitude(load, t);
    }
    impulse *= load.duration / samples;
    ImGui::TextDisabled("Impulse %.1f N*s", impulse);
    ImGui::TextDisabled(
        "Force moves from initial to final over the ramp, then holds until the duration ends. A ramp time of zero applies the final force immediately.");
    if (load.rampTime > load.duration + 1.0e-6) {
        const double reached =
            loadMagnitude(load, load.startTime + std::max(0.0, load.duration - 1.0e-6));
        ImGui::TextColored({0.95f, 0.6f, 0.25f, 1.0f},
                           "Duration is shorter than the ramp, so this hit only reaches %.0f N "
                           "(impulse %.1f N*s), not the final force.",
                           reached, impulse);
    }

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

void drawHeatPanel(EditorState& state) {
    if (!ImGui::CollapsingHeader("Heat sources"))
        return;

    auto& sources = state.scenario.heatSources;
    if (ImGui::Button("Add heat"))
        state.addHeatRequested = true;
    ImGui::SameLine();
    ImGui::BeginDisabled(state.selectedHeat < 0);
    if (ImGui::Button("Remove##heat"))
        state.removeHeatIndex = state.selectedHeat;
    ImGui::EndDisabled();

    if (ImGui::BeginListBox("##heat", {-1.0f, 72.0f})) {
        for (int i = 0; i < int(sources.size()); ++i) {
            char label[160];
            std::snprintf(label, sizeof(label), "%s  %.0f K @ %.2f s%s", sources[i].name.c_str(),
                          sources[i].temperature, sources[i].startTime,
                          sources[i].enabled ? "" : "  (off)");
            if (ImGui::Selectable(label, state.selectedHeat == i))
                state.selectedHeat = i;
        }
        ImGui::EndListBox();
    }

    if (state.selectedHeat < 0 || state.selectedHeat >= int(sources.size())) {
        ImGui::TextDisabled("A sphere. No collision. Hotter than the air adds energy.");
        return;
    }

    HeatSource& source = sources[std::size_t(state.selectedHeat)];
    char nameBuffer[64];
    std::snprintf(nameBuffer, sizeof(nameBuffer), "%s", source.name.c_str());
    if (ImGui::InputText("Name##heat", nameBuffer, sizeof(nameBuffer)))
        source.name = nameBuffer;
    ImGui::Checkbox("Enabled##heat", &source.enabled);
    vectorControl("Centre (m)", source.point, 0.01f);
    scalarControl("Radius (m)##heat", source.radius, 0.005f, 0.01f, 5.0f, "%.3f");
    scalarControl("Temperature (K)", source.temperature, 5.0f, 0.0f, 4000.0f, "%.1f");
    scalarControl("Start (s)##heat", source.startTime, 0.01f, 0.0f, 60.0f, "%.3f");
    scalarControl("Duration (s)##heat", source.duration, 0.01f, 0.0f, 60.0f, "%.3f");
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
        scalarControl("Resolution", state.scenario.settings.meshResolution, 0.05f, 0.25f, 8.0f, "%.2f");
        ImGui::TextDisabled("1 is the part's own element size. Higher adds tetrahedra everywhere.");
        ImGui::Checkbox("Dynamic Resolution", &state.scenario.settings.dynamicResolution);
        ImGui::TextDisabled("When on, a hot or loaded spot is split finer than this mesh.");
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
    drawHeatPanel(state);
    drawRunPanel(state, context);

    ImGui::SeparatorText("View");
    int shade = state.shade == ShadeMode::Temperature ? 2 : state.shade == ShadeMode::Stress ? 1 : 0;
    if (ImGui::RadioButton("Material##shade", shade == 0))
        state.shade = ShadeMode::Material;
    ImGui::SameLine();
    if (ImGui::RadioButton("Stress", shade == 1))
        state.shade = ShadeMode::Stress;
    ImGui::SameLine();
    if (ImGui::RadioButton("Temperature", shade == 2))
        state.shade = ShadeMode::Temperature;
    ImGui::Checkbox("Show air temperature", &state.showAir);

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
    if (state.shade == ShadeMode::Temperature)
        drawTemperatureLegend();

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

    if (state.probeVisible &&
        (state.shade == ShadeMode::Temperature || state.shade == ShadeMode::Stress)) {
        const ImVec2 mouse = ImGui::GetMousePos();
        ImGui::SetNextWindowPos({mouse.x + 16.0f, mouse.y + 16.0f}, ImGuiCond_Always);
        ImGui::SetNextWindowBgAlpha(0.92f);
        ImGui::Begin("##probe", nullptr,
                     ImGuiWindowFlags_NoDecoration | ImGuiWindowFlags_AlwaysAutoResize |
                         ImGuiWindowFlags_NoInputs | ImGuiWindowFlags_NoNav | ImGuiWindowFlags_NoMove);
        if (state.shade == ShadeMode::Temperature)
            ImGui::Text("%.1f K", state.probeValue);
        else
            ImGui::Text("%s", formatStress(state.probeValue).c_str());
        ImGui::End();
    }

    drawPlaybackWindow(state, context);
}

void drawArenaOverlay(const Environment& environment,
                      const Magnum::Matrix4& view,
                      const Magnum::Matrix4& projection)
{
    if (!environment.floorEnabled)
        return;

    const Magnum::Matrix4 clipFromWorld = projection * view;
    const ImVec2 display = ImGui::GetIO().DisplaySize;
    auto project = [&](const Magnum::Vector3& world, ImVec2& screen) {
        const Magnum::Vector4 clip = clipFromWorld * Magnum::Vector4{world, 1.0f};
        if (clip.w() <= 0.05f)
            return false;
        const float x = clip.x() / clip.w();
        const float y = clip.y() / clip.w();
        if (x < -1.05f || x > 1.05f || y < -1.05f || y > 1.05f)
            return false;
        screen.x = (x * 0.5f + 0.5f) * display.x;
        screen.y = (1.0f - (y * 0.5f + 0.5f)) * display.y;
        return true;
    };
    auto label = [](ImDrawList* draw, const ImVec2& at, const char* text, bool strong) {
        const ImU32 fill = strong ? IM_COL32(255, 255, 255, 255) : IM_COL32(220, 224, 230, 230);
        const ImVec2 size = ImGui::CalcTextSize(text);
        const ImVec2 pos{at.x - size.x * 0.5f, at.y - size.y * 0.5f};
        draw->AddText({pos.x + 1.0f, pos.y + 1.0f}, IM_COL32(0, 0, 0, 200), text);
        draw->AddText(pos, fill, text);
    };
    auto formatLength = [](double metres, char* buffer, std::size_t count) {
        if (std::abs(metres - std::round(metres)) < 0.05)
            std::snprintf(buffer, count, "%.0f m", metres);
        else
            std::snprintf(buffer, count, "%.2f m", metres);
    };

    const double hx = std::max(0.5, environment.boundsSize.x() * 0.5);
    const double hz = std::max(0.5, environment.boundsSize.z() * 0.5);
    const double height = std::max(1.0, environment.boundsSize.y());
    const float y0 = float(environment.floorHeight);
    ImDrawList* draw = ImGui::GetBackgroundDrawList();

    const auto ticks = [&](double limit, double step, auto place, auto textFor) {
        const double start = std::ceil((-limit) / step - 1.0e-4) * step;
        for (double t = start; t <= limit + 1.0e-3; t += step) {
            ImVec2 screen;
            if (!project(place(float(t)), screen))
                continue;
            char buffer[32];
            textFor(t, buffer, sizeof(buffer));
            label(draw, screen, buffer, false);
        }
    };
    const double stepX = Environment::gridStep(environment.boundsSize.x());
    const double stepZ = Environment::gridStep(environment.boundsSize.z());
    ticks(hx, stepX,
          [&](float x) { return Magnum::Vector3{x, y0, float(-hz) - 0.08f}; },
          [](double t, char* buffer, std::size_t count) {
              if (std::abs(t) < 1.0e-3)
                  std::snprintf(buffer, count, "0 m");
              else
                  std::snprintf(buffer, count, "%.0f", t);
          });
    ticks(hz, stepZ,
          [&](float z) { return Magnum::Vector3{float(-hx) - 0.08f, y0, z}; },
          [](double t, char* buffer, std::size_t count) {
              if (std::abs(t) < 1.0e-3)
                  std::snprintf(buffer, count, "0 m");
              else
                  std::snprintf(buffer, count, "%.0f", t);
          });

    char length[32];
    const Magnum::Vector3 axisEdge[3] = {
        {0.0f, y0 + 0.35f, float(-hz) - 0.45f},
        {float(-hx) - 0.45f, y0 + float(height) * 0.5f, float(-hz) - 0.25f},
        {float(-hx) - 0.45f, y0 + 0.35f, 0.0f}
    };
    const double axisLength[3] = {environment.boundsSize.x(), height, environment.boundsSize.z()};
    for (int i = 0; i < 3; ++i) {
        ImVec2 screen;
        if (!project(axisEdge[i], screen))
            continue;
        formatLength(axisLength[i], length, sizeof(length));
        label(draw, {screen.x, screen.y - 14.0f}, length, true);
    }
}

}  // namespace rcs
