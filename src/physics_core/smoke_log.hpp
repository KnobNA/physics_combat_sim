#pragma once

#include "physics_core/mesh_split.hpp"
#include "physics_core/simulator.hpp"

#include <algorithm>
#include <chrono>
#include <cmath>
#include <cstdio>
#include <cstring>
#include <string>
#include <vector>

namespace {

const char* currentGroup = "smoke";
int failures = 0;
int passes = 0;
bool keepEditorDefaultsFlag = false;
bool groupKeepsEditorDefaults = false;
std::chrono::steady_clock::time_point groupStart = std::chrono::steady_clock::now();

void keepDynamicResolution() {}
void keepEditorDefaults() { keepEditorDefaultsFlag = true; }

void beginGroup(const char* name) {
    currentGroup = name ? name : "smoke";
    groupStart = std::chrono::steady_clock::now();
    groupKeepsEditorDefaults = std::strcmp(currentGroup, "editor") == 0;
    keepEditorDefaultsFlag = false;
    std::printf("\n== group %s ==\n", currentGroup);
    if (groupKeepsEditorDefaults)
        std::printf("SMOKE group=%s meshResolution=1.0\n", currentGroup);
    else
        std::printf("SMOKE group=%s meshResolution=%s\n", currentGroup, "0.25");
    std::fflush(stdout);
}

void applySmokeDefaults(rcs::Scenario& scenario) {
    if (keepEditorDefaultsFlag || groupKeepsEditorDefaults) {
        scenario.settings.meshResolution = 1.0;
        std::printf("SMOKE group=%s meshResolution=1.0 (editor)\n", currentGroup);
        return;
    }
    scenario.settings.meshResolution = 0.25;
}

void applyCoarseMesh(rcs::PartDefinition& part) {
    if (part.sizing.cached && !part.volume.empty())
        return;
    const int base = std::max(200, part.sizing.tetBudget);
    part.sizing.tetBudget = std::clamp(static_cast<int>(base * 0.25), 400, 8000);
    part.sizing.cached = false;
    part.volume = {};
    part.surface = {};
}

void applyCoarseAssembly(rcs::AssemblyDefinition& definition) {
    if (keepEditorDefaultsFlag || groupKeepsEditorDefaults) {
        for (auto& part : definition.parts)
            rcs::applyResolutionSizing(part, 1.0);
        keepEditorDefaultsFlag = false;
        return;
    }
    for (auto& part : definition.parts)
        applyCoarseMesh(part);
}

bool meshPart(rcs::PartDefinition& part) {
    applyCoarseMesh(part);
    return rcs::ensurePartMesh(part);
}

double dumpElapsed() {
    const auto now = std::chrono::steady_clock::now();
    return std::chrono::duration<double>(now - groupStart).count();
}

void dumpSim(const char* tag, const rcs::Simulator& simulator, double t = 0.0, bool finite = true) {
    int liveTets = 0;
    int totalTets = 0;
    std::size_t shards = 0;
    int torn = 0;
    double mass = simulator.assembly().mass();
    chrono::ChVector3d com{0, 0, 0};
    double comMass = 0.0;
    double lowest = 1.0e9;
    double ke = 0.0;
    double peakStress = 0.0;
    double peakPlastic = 0.0;
    double peakTension = 0.0;
    double worstAspect = 1.0;
    for (const auto& part : simulator.assembly().parts()) {
        liveTets += part.liveElementCount();
        totalTets += part.totalElementCount();
        shards += part.shards().size();
        torn += part.torn() ? 1 : 0;
        peakTension = std::max(peakTension, part.peakTension());
        std::vector<float> stress;
        part.evaluateStress(stress);
        for (float value : stress)
            peakStress = std::max(peakStress, static_cast<double>(value));
        std::vector<float> plastic;
        part.evaluatePlasticStrain(plastic);
        for (float value : plastic) {
            if (std::isfinite(value))
                peakPlastic = std::max(peakPlastic, static_cast<double>(value));
        }
        for (std::size_t n = 0; n < part.nodeCount(); ++n) {
            const double m = n < part.nodeMass().size() ? part.nodeMass()[n] : 0.0;
            const auto p = part.nodePosition(n);
            const auto v = part.nodeVelocity(n);
            com += m * p;
            comMass += m;
            lowest = std::min(lowest, p.y());
            ke += 0.5 * m * v.Length2();
        }
        for (const auto& shard : part.shards())
            ke += 0.5 * shard.mass * shard.velocity.Length2();
        for (int ti = 0; ti < part.tetCount(); ++ti) {
            if (!part.tetLive(static_cast<std::size_t>(ti)))
                continue;
            const auto& tet = part.tets()[static_cast<std::size_t>(ti)];
            double longest = 0.0;
            double shortest = 1.0e9;
            for (int a = 0; a < 4; ++a) {
                const int b = (a + 1) % 4;
                if (tet[a] < 0 || tet[b] < 0)
                    continue;
                const double edge = (part.nodePosition(static_cast<std::size_t>(tet[a])) -
                                     part.nodePosition(static_cast<std::size_t>(tet[b])))
                                        .Length();
                longest = std::max(longest, edge);
                shortest = std::min(shortest, edge);
            }
            if (shortest > 1.0e-9)
                worstAspect = std::max(worstAspect, longest / shortest);
        }
    }
    if (comMass > 1.0e-12)
        com /= comMass;
    std::printf("DUMP tag=%s group=%s t=%.4f wall=%.3fs finite=%d liveTets=%d totalTets=%d shards=%zu "
                "torn=%d mass=%.4f com=(%.4f,%.4f,%.4f) lowest=%.4f ke=%.4f stress=%.3e plastic=%.5f "
                "tension=%.3e broken=%d aspect=%.2f newtonRejects=%d inverted=%d\n",
                tag ? tag : "", currentGroup, t, dumpElapsed(), finite ? 1 : 0, liveTets, totalTets, shards,
                torn, mass, com.x(), com.y(), com.z(), lowest, ke, peakStress, peakPlastic, peakTension,
                simulator.assembly().brokenLinkCount(), worstAspect, simulator.newtonRejectCount(),
                simulator.assembly().invertedTetCount());
    std::fflush(stdout);
}

void emitCheck(const char* what, bool ok, double got = 0.0, double lo = 0.0, double hi = 0.0,
               bool ranged = false) {
    if (ranged) {
        std::printf("CHECK group=%s name=\"%s\" result=%s got=%.6g want=[%.6g,%.6g]\n", currentGroup, what,
                    ok ? "ok" : "FAIL", got, lo, hi);
    } else {
        std::printf("CHECK group=%s name=\"%s\" result=%s\n", currentGroup, what, ok ? "ok" : "FAIL");
    }
    std::fflush(stdout);
}

void check(bool condition, const char* what) {
    std::printf("%-46s %s\n", what, condition ? "ok" : "FAILED");
    emitCheck(what, condition);
    std::fflush(stdout);
    if (condition)
        ++passes;
    else
        ++failures;
}

void expectFinite(bool finite, const char* what) { check(finite, what); }

void expectRange(double value, double lo, double hi, const char* what) {
    const bool ok = std::isfinite(value) && std::isfinite(lo) && std::isfinite(hi) && value >= lo &&
                    value <= hi;
    if (ok) {
        std::printf("%-46s ok\n", what);
        ++passes;
    } else {
        std::printf("%-46s FAILED  (got %.6g, want [%.6g, %.6g])\n", what, value, lo, hi);
        ++failures;
    }
    emitCheck(what, ok, value, lo, hi, true);
    std::fflush(stdout);
}

void expectRange(double value, double lo, double hi, const std::string& what) {
    expectRange(value, lo, hi, what.c_str());
}

class SmokeSim : public rcs::Simulator {
public:
    void reset(const rcs::AssemblyDefinition& definition, rcs::Scenario scenario, bool mesh = true) {
        applySmokeDefaults(scenario);
        rcs::AssemblyDefinition meshed = definition;
        applyCoarseAssembly(meshed);
        rcs::Simulator::reset(meshed, scenario, mesh);
        dumpSim("reset", *this, 0.0, true);
    }

    rcs::Recording run(const rcs::AssemblyDefinition& definition, rcs::Scenario scenario,
                       rcs::RunControl* control) {
        applySmokeDefaults(scenario);
        rcs::AssemblyDefinition meshed = definition;
        applyCoarseAssembly(meshed);
        rcs::Recording recording = rcs::Simulator::run(meshed, scenario, control);
        dumpSim("run", *this, recording.duration, !recording.diverged);
        return recording;
    }
};

}  // namespace
