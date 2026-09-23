#include "physics_core/simulator.hpp"

#include "chrono/ChConfig.h"
#include "chrono/solver/ChDirectSolverLS.h"
#include "chrono/timestepper/ChTimestepperHHT.h"

#include <exception>
#include "chrono_pardisomkl/ChSolverPardisoMKL.h"

#ifndef CHRONO_PARDISO_MKL
#error "combat_sim requires Chrono built with CH_ENABLE_MODULE_PARDISO_MKL (Intel oneMKL Pardiso)."
#endif

#include <algorithm>
#include <cmath>
#include <thread>

namespace rcs {
namespace {

int solverThreadCount() {
    const int hardware = static_cast<int>(std::thread::hardware_concurrency());
    return std::clamp(hardware, 1, 8);
}

}  // namespace

Simulator::Simulator() {
    using namespace chrono;

    const int threads = solverThreadCount();
    _system.SetNumThreads(threads, threads, threads);

    auto solver = chrono_types::make_shared<ChSolverPardisoMKL>(static_cast<unsigned int>(threads));
    solver->UseSparsityPatternLearner(true);
    solver->LockSparsityPattern(true);
    _system.SetSolver(solver);

    // HHT rather than the single-iteration linearized Euler. That stepper folds
    // the entire weld-constraint drift into one velocity correction (C/dt) and
    // accepts it. On a tumbling welded shell the drift grows quietly and then
    // arrives as a shove while nothing is touching the box. HHT iterates the
    // correction out, and if a step still will not converge it cuts the step
    // instead of keeping the bad one.
    _system.SetTimestepperType(ChTimestepper::Type::HHT);
    if (auto hht = std::dynamic_pointer_cast<ChTimestepperHHT>(_system.GetTimestepper())) {
        hht->SetAlpha(-0.05);
        hht->SetMaxIters(10);
        hht->SetRelTolerance(1e-4);
        hht->SetAbsTolerances(1e-6, 1e-4);
        hht->SetStepControl(true);
        hht->SetMinStepSize(1e-6);
        hht->AcceptTerminatedStep(false);
    }

    // Every external force (gravity, drag, contact, user loads) is assembled
    // per node in Assembly::applyLoads, so Chrono must not add its own.
    _system.SetGravitationalAcceleration(ChVector3d(0, 0, 0));
    _system.SetCollisionSystemType(ChCollisionSystem::Type::BULLET);
}

void Simulator::reset(const AssemblyDefinition& definition, const Scenario& scenario) {
    _system.Clear();
    _assembly = Assembly{};
    _assembly.build(_system, definition);
    _assembly.place(scenario.initial);
    _assembly.bindLoads(scenario.loads);

    // Relative tolerance, plus the two absolute ones kept in the same ratio as
    // the constructor defaults (states 100x tighter, multipliers equal).
    if (auto hht = std::dynamic_pointer_cast<chrono::ChTimestepperHHT>(_system.GetTimestepper())) {
        const double tol = std::max(scenario.settings.solverTolerance, 1.0e-12);
        hht->SetRelTolerance(tol);
        hht->SetAbsTolerances(tol * 1.0e-2, tol);
    }

    _system.Setup();
    _system.Update(chrono::UpdateFlags::UPDATE_ALL);
    refreshLinearSolver();
}

bool Simulator::step(double dt, double time, const Environment& environment,
                     const std::vector<ImpactLoad>& loads)
{
    _assembly.applyLoads(environment, loads, time, dt);
    try {
        _system.DoStepDynamics(dt);
    } catch (const std::exception&) {
        // HHT throws when it has already cut the step to the minimum and the
        // Newton solve still will not converge. That is a diverged run.
        return false;
    }
    _assembly.integrateShards(dt, environment);
    const bool joints = _assembly.updateJoints(dt);
    const bool torn = _assembly.tearOverstretched();
    const bool spanned = (joints || torn) && _assembly.respan();
    if (joints || torn || spanned) {
        // New nodes change the coordinate count. Constraint-only edits (a weld
        // breaking) already reanalyse from ForceUpdate, which is what the
        // weak-weld run depends on.
        if (torn)
            _system.Setup();
        refreshLinearSolver();
    }
    return _assembly.finite();
}

void Simulator::refreshLinearSolver() {
    // Breaking a link changes the constraint block, so the cached sparsity
    // pattern has to be rebuilt before the next factorisation. HHT only
    // re-runs that analysis when the system reports itself modified.
    _system.ForceUpdate();
    if (auto solver = std::dynamic_pointer_cast<chrono::ChDirectSolverLS>(_system.GetSolver()))
        solver->ForceSparsityPatternUpdate();
}

Recording Simulator::capture() const {
    Recording recording;
    _assembly.writeLayout(recording.layout);
    Frame frame;
    frame.time = 0.0;
    _assembly.writeFrame(frame);
    recording.maxStress = frame.maxStress;
    recording.brokenLinks = frame.brokenLinks;
    recording.frames.push_back(std::move(frame));
    return recording;
}

Recording Simulator::run(const AssemblyDefinition& definition,
                         const Scenario& scenario,
                         RunControl* control)
{
    reset(definition, scenario);

    Recording recording;
    _assembly.writeLayout(recording.layout);

    const double dt = std::max(1.0e-5, scenario.settings.timestep);
    const double duration = std::max(dt, scenario.settings.duration);
    const double captureInterval =
        1.0 / std::max(1.0, scenario.settings.captureRate);
    const auto totalSteps = static_cast<long long>(std::ceil(duration / dt));

    auto captureFrame = [&](double time) {
        Frame frame;
        frame.time = time;
        _assembly.writeFrame(frame);
        recording.maxStress = std::max(recording.maxStress, frame.maxStress);
        recording.brokenLinks = std::max(recording.brokenLinks, frame.brokenLinks);
        recording.frames.push_back(std::move(frame));
    };

    double time = 0.0;
    double nextCapture = 0.0;
    captureFrame(time);
    nextCapture += captureInterval;

    for (long long i = 0; i < totalSteps; ++i) {
        if (control && control->cancel.load()) {
            recording.message = "Cancelled";
            break;
        }

        if (!step(dt, time, scenario.environment, scenario.loads)) {
            recording.diverged = true;
            recording.message = "Solver diverged at t = " + std::to_string(time) + " s";
            break;
        }
        time += dt;

        if (time + 1.0e-12 >= nextCapture) {
            captureFrame(time);
            nextCapture += captureInterval;
        }
        if (control)
            control->progress.store(static_cast<float>(time / duration));
    }

    if (recording.frames.empty() || recording.frames.back().time < time)
        captureFrame(time);

    recording.duration = time;
    if (recording.message.empty()) {
        recording.message = "Simulated " + std::to_string(recording.frames.size()) +
                            " frames over " + std::to_string(time) + " s";
    }
    if (control)
        control->progress.store(1.0f);
    return recording;
}

SimulationJob::~SimulationJob() {
    cancel();
    join();
}

void SimulationJob::join() {
    if (_thread.joinable())
        _thread.join();
    _running = false;
}

void SimulationJob::start(const AssemblyDefinition& definition, const Scenario& scenario) {
    cancel();
    join();

    _control.cancel.store(false);
    _control.progress.store(0.0f);
    _done.store(false);
    _hasResult = false;
    _running = true;

    _thread = std::thread([this, definition, scenario]() {
        // The worker owns its own Chrono system; nothing is shared with the
        // viewport except the control block and the finished recording.
        Simulator simulator;
        Recording recording = simulator.run(definition, scenario, &_control);
        {
            std::lock_guard<std::mutex> lock(_mutex);
            _result = std::move(recording);
            _hasResult = true;
        }
        _done.store(true);
    });
}

void SimulationJob::cancel() {
    _control.cancel.store(true);
}

bool SimulationJob::collect(Recording& out) {
    if (!_done.load())
        return false;
    join();

    std::lock_guard<std::mutex> lock(_mutex);
    if (!_hasResult)
        return false;
    out = std::move(_result);
    _hasResult = false;
    _done.store(false);
    return true;
}

}  // namespace rcs
