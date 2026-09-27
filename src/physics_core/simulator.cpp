#include "physics_core/simulator.hpp"
#include "physics_core/mesh_split.hpp"

#include "chrono/ChConfig.h"
#include "chrono/solver/ChDirectSolverLS.h"
#include "chrono/timestepper/ChTimestepperHHT.h"

#include <exception>
#include <stdexcept>
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

void configureHht(chrono::ChTimestepperHHT& hht) {
    hht.SetAlpha(-0.05);
    hht.SetMaxIters(30);
    hht.SetRelTolerance(1e-4);
    hht.SetAbsTolerances(1e-6, 1e-4);
    hht.SetStepControl(true);
    hht.SetMinStepSize(1e-6);
    hht.SetRequiredSuccessfulSteps(5);
    hht.SetStepIncreaseFactor(2.0);
    hht.AcceptTerminatedStep(false);
    hht.SetJacobianUpdateMethod(
        chrono::ChTimestepperImplicit::JacobianUpdate::EVERY_STEP);
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
    if (auto hht = std::dynamic_pointer_cast<ChTimestepperHHT>(_system.GetTimestepper()))
        configureHht(*hht);

    // Every external force (gravity, drag, contact, user loads) is assembled
    // per node in Assembly::applyLoads, so Chrono must not add its own.
    _system.SetGravitationalAcceleration(ChVector3d(0, 0, 0));
    _system.SetCollisionSystemType(ChCollisionSystem::Type::BULLET);
}

void Simulator::reset(const AssemblyDefinition& definition, const Scenario& scenario, bool mesh) {
    _system.Clear();
    _assembly = Assembly{};
    AssemblyDefinition meshed = definition;
    if (mesh) {
        for (auto& part : meshed.parts) {
            if (part.simulate)
                ensurePartMesh(part);
        }
    }
    _assembly.build(_system, meshed);
    _assembly.setDynamicResolution(scenario.settings.dynamicResolution);
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
    _contactStepFailed = false;
}

bool Simulator::step(double dt, double time, const Environment& environment,
                     const std::vector<ImpactLoad>& loads,
                     const std::vector<HeatSource>& heat)
{
    _assembly.advanceHeat(dt, time, environment, heat);
    const int impact = _assembly.impactSubsteps(environment, dt);
    bool alreadyTorn = false;
    for (const auto& part : _assembly.parts())
        alreadyTorn = alreadyTorn || part.torn();
    const bool allowRefine = impact < 8 && !alreadyTorn;
    const bool volume = _assembly.maintainVolume(dt, allowRefine);
    const bool joints = _assembly.updateJoints(dt);
    if (volume || joints) {
        if (volume)
            _system.Setup();
        refreshLinearSolver();
    }
    if (_assembly.applyLoads(environment, loads, time, dt)) {
        // Contact elements carry a Jacobian, so adding one changes the matrix.
        _system.Setup();
        refreshLinearSolver();
    }
    bool rejected = false;
    try {
        if (_assembly.hasFreeSolve()) {
            const int subs = std::max(1, _assembly.impactSubsteps(environment, dt));
            auto hht = std::dynamic_pointer_cast<chrono::ChTimestepperHHT>(_system.GetTimestepper());
            auto restoreHht = [&]() {
                if (hht)
                    configureHht(*hht);
            };
            // High-speed floor hits are already split to ~2 mm. Cutting
            // further on thousands of tets is the glass-test stall only when
            // the Advance dt is still 0.5 ms. Inside a 2 mm substep, a cut
            // to 1e-5 is a handful of Newton loops and is how J2/Rankine
            // actually get to commit. A miss keeps Chrono's last accepted
            // state; rewinding retries the same hit forever.
            const bool bulkyContact = environment.floorEnabled && subs >= 8;
            if (hht && environment.floorEnabled) {
                hht->SetJacobianUpdateMethod(
                    chrono::ChTimestepperImplicit::JacobianUpdate::EVERY_ITERATION);
                if (bulkyContact)
                    hht->SetMinStepSize(1.0e-5);
            }
            if (hht && _contactStepFailed) {
                // After a min-stepsize throw Chrono leaves h below h_min.
                // One Advance with StepControl off sets h to this
                // substep so the next cut starts from the user dt, not 1e-5.
                hht->SetStepControl(false);
            }
            auto runSubs = [&](int n) {
                const double sub = dt / static_cast<double>(n);
                _assembly.setContactStep(sub);
                for (int i = 0; i < n; ++i)
                    _system.DoStepDynamics(sub);
            };
            _assembly.captureKinematics();
            try {
                runSubs(subs);
                _contactStepFailed = false;
            } catch (...) {
                _contactStepFailed = true;
                restoreHht();
                _assembly.discardTrial();
                rejected = true;
                // A coupon Newton miss at the user dt must not spin at
                // h_min (tempered/rubber hang). Crawl a few fixed slices
                // with cutting off so the outer clock can unload.
                if (hht && !environment.floorEnabled) {
                    hht->SetStepControl(false);
                    const int slices = 8;
                    const double crawl = dt / static_cast<double>(slices);
                    _assembly.setContactStep(crawl);
                    int got = 0;
                    for (int i = 0; i < slices; ++i) {
                        try {
                            _system.DoStepDynamics(crawl);
                            ++got;
                        } catch (...) {
                            break;
                        }
                    }
                    if (got > 0)
                        rejected = false;
                }
            }
            if (!_assembly.finite()) {
                _assembly.restoreKinematics();
                _assembly.discardTrial();
                _assembly.discardAwakeShards();
                rejected = true;
            }
            restoreHht();
        }
    } catch (const std::exception&) {
        // A rejected HHT step keeps the last accepted positions. Discard the
        // tet trial state so a failed 50 m/s hit cannot store a 1e8 plastic
        // strain and then explode the mesh. Do not peel or tear from that
        // iterate; the next user step tries again on the intact core.
        _assembly.restoreKinematics();
        _assembly.discardTrial();
        _assembly.discardAwakeShards();
        rejected = true;
        if (!_assembly.finite())
            return false;
    }
    _assembly.integrateShards(dt, environment);
    bool torn = false;
    if (!rejected)
        torn = _assembly.tearOverstretched();
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
    if (control) {
        control->duration.store(duration);
        control->stepCount.store(totalSteps);
        control->simulatedTime.store(0.0);
        control->step.store(0);
    }

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

        if (!step(dt, time, scenario.environment, scenario.loads, scenario.heatSources)) {
            recording.diverged = true;
            recording.message = "Solver diverged at t = " + std::to_string(time) + " s";
            break;
        }
        time += dt;

        if (time + 1.0e-12 >= nextCapture) {
            captureFrame(time);
            nextCapture += captureInterval;
        }
        if (control) {
            control->simulatedTime.store(time);
            control->step.store(i + 1);
            control->progress.store(static_cast<float>(time / duration));
        }
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
    _control.simulatedTime.store(0.0);
    _control.duration.store(std::max(scenario.settings.timestep, scenario.settings.duration));
    _control.step.store(0);
    _control.stepCount.store(0);
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
