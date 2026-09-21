#pragma once

#include "physics_core/assembly.hpp"
#include "physics_core/recording.hpp"
#include "physics_core/scenario.hpp"

#include "chrono/physics/ChSystemSMC.h"

#include <atomic>
#include <mutex>
#include <thread>

namespace rcs {

/// Shared state a running simulation exposes to the UI thread.
struct RunControl {
    std::atomic<float> progress{0.0f};
    std::atomic<bool> cancel{false};
};

/// Owns one Chrono system and the assembly inside it. A Simulator is always
/// rebuilt from the authoring definition, so a run always starts from t = 0.
class Simulator {
public:
    Simulator();

    /// Tear down and rebuild the system from scratch at the scenario's
    /// initial state.
    void reset(const AssemblyDefinition& definition, const Scenario& scenario);

    /// One integrator step at simulated time `time`. False if the state blew up.
    bool step(double dt, double time, const Environment& environment,
              const std::vector<ImpactLoad>& loads);

    /// Snapshot the current state as a single-frame recording, for the
    /// authoring viewport.
    Recording capture() const;

    /// Run the whole scenario and return the playback recording.
    Recording run(const AssemblyDefinition& definition,
                  const Scenario& scenario,
                  RunControl* control);

    Assembly& assembly() { return _assembly; }
    const Assembly& assembly() const { return _assembly; }
    chrono::ChSystemSMC& system() { return _system; }

private:
    void refreshLinearSolver();

    chrono::ChSystemSMC _system;
    Assembly _assembly;
};

/// Runs a scenario on a worker thread so a multi-second batch never stalls
/// the viewport.
class SimulationJob {
public:
    ~SimulationJob();

    void start(const AssemblyDefinition& definition, const Scenario& scenario);
    void cancel();
    /// Moves the finished recording into `out`. True exactly once per run.
    bool collect(Recording& out);

    bool running() const { return _running; }
    float progress() const { return _control.progress.load(); }

private:
    void join();

    std::thread _thread;
    RunControl _control;
    std::mutex _mutex;
    Recording _result;
    bool _hasResult = false;
    std::atomic<bool> _done{false};
    bool _running = false;
};

}  // namespace rcs
