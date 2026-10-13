#pragma once

#include "physics_core/assembly.hpp"
#include "physics_core/recording.hpp"
#include "physics_core/scenario.hpp"

#include "chrono/physics/ChSystemSMC.h"
#include "chrono/timestepper/ChTimestepperExplicit.h"

#include <atomic>
#include <memory>
#include <mutex>
#include <thread>

namespace rcs {

/// Shared state a running simulation exposes to the UI thread.
struct RunControl {
    std::atomic<float> progress{0.0f};
    std::atomic<double> simulatedTime{0.0};
    std::atomic<double> duration{0.0};
    std::atomic<long long> step{0};
    std::atomic<long long> stepCount{0};
    std::atomic<bool> cancel{false};
};

/// Owns one Chrono system and the assembly inside it. A Simulator is always
/// rebuilt from the authoring definition, so a run always starts from t = 0.
class Simulator {
public:
    Simulator();

    /// Tear down and rebuild the system from scratch at the scenario's
    /// initial state. `mesh` runs tetgen when a part still needs a volume;
    /// the editor preview passes false so opening a file does not stall.
    void reset(const AssemblyDefinition& definition, const Scenario& scenario,
               bool mesh = true);

    /// One integrator step at simulated time `time`. False if the state blew up.
    /// When `cancel` is set, HHT stops starting new increments after the
    /// current Newton attempt returns, instead of draining the rest of `dt`.
    bool step(double dt, double time, const Environment& environment,
              const std::vector<ImpactLoad>& loads,
              const std::vector<HeatSource>& heat = {},
              std::atomic<bool>* cancel = nullptr);

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
    int newtonRejectCount() const { return _newtonRejects; }
    /// HHT substeps that kept a finite, non-inverted iterate.
    int hhtCommitCount() const { return _hhtCommits; }

private:
    void refreshLinearSolver();
    void ensureExplicit(bool lumped);
    bool tryExplicitStep(double h, bool lumped);

    chrono::ChSystemSMC _system;
    Assembly _assembly;
    std::unique_ptr<chrono::ChTimestepperEulerExplicitIIorder> _explicit;
    bool _explicitLumped = false;
    bool _contactStepFailed = false;
    int _newtonRejects = 0;
    int _hhtCommits = 0;
    int _prevContactCount = -1;
    int _newtonSuccessStreak = 0;
    double _localH = 0.0;
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
    bool cancelling() const { return _running && _control.cancel.load(); }
    float progress() const { return _control.progress.load(); }
    double simulatedTime() const { return _control.simulatedTime.load(); }
    double runDuration() const { return _control.duration.load(); }
    long long stepIndex() const { return _control.step.load(); }
    long long stepCount() const { return _control.stepCount.load(); }

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
