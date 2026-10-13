#include "physics_core/simulator.hpp"
#include "physics_core/mesh_split.hpp"
#include "physics_core/physics_units.hpp"

#include "chrono/ChConfig.h"
#include "chrono/solver/ChDirectSolverLS.h"
#include "chrono/timestepper/ChTimestepperExplicit.h"
#include "chrono/timestepper/ChTimestepperHHT.h"

#include <exception>
#include <stdexcept>
#include "chrono_pardisomkl/ChSolverPardisoMKL.h"

#ifndef CHRONO_PARDISO_MKL
#error "combat_sim requires Chrono built with CH_ENABLE_MODULE_PARDISO_MKL (Intel oneMKL Pardiso)."
#endif

#include <algorithm>
#include <cmath>
#include <iostream>
#include <sstream>
#include <thread>

namespace rcs {
namespace {

int solverThreadCount() {
    const int hardware = static_cast<int>(std::thread::hardware_concurrency());
    return std::clamp(hardware, 1, 8);
}

void configureHht(chrono::ChTimestepperHHT& hht) {
    hht.SetAlpha(-0.05);
    hht.SetMaxIters(12);
    hht.SetRelTolerance(1e-4);
    hht.SetAbsTolerances(1e-6, 1e-4);
    // We own incrementation: reject a terminated Newton, then cut h.
    // Extra iters only at h_min or a seated increment. Explicit is the
    // h_min fallback in contact.
    hht.SetStepControl(false);
    hht.SetMinStepSize(kHhtMinStep);
    hht.SetRequiredSuccessfulSteps(2);
    hht.SetStepIncreaseFactor(2.0);
    hht.AcceptTerminatedStep(false);
    hht.SetJacobianUpdateMethod(
        chrono::ChTimestepperImplicit::JacobianUpdate::AUTOMATIC);
}

bool tryHhtStep(chrono::ChSystemSMC& system, double h) {
    // Chrono writes "[ERROR] Newton did not converge" on the throw we
    // catch and cut. That is incrementation, not a failed run.
    std::ostringstream sink;
    auto* prev = std::cerr.rdbuf(sink.rdbuf());
    bool threw = false;
    try {
        system.DoStepDynamics(h);
    } catch (...) {
        threw = true;
    }
    std::cerr.rdbuf(prev);
    return !threw;
}

bool assemblyHasRotNodes(const Assembly& assembly) {
    for (const auto& part : assembly.parts()) {
        if (!part.nodes().empty())
            return true;
    }
    return false;
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
        seedWeldBeads(meshed);
        for (auto& part : meshed.parts) {
            if (part.simulate)
                ensurePartMesh(part);
        }
    }
    _assembly.build(_system, meshed);
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
    _newtonRejects = 0;
    _hhtCommits = 0;
    _prevContactCount = -1;
    _newtonSuccessStreak = 0;
    _localH = 0.0;
    _explicit.reset();
    _explicitLumped = false;
}

bool Simulator::step(double dt, double time, const Environment& environment,
                     const std::vector<ImpactLoad>& loads,
                     const std::vector<HeatSource>& heat,
                     std::atomic<bool>* cancel)
{
    _assembly.advanceHeat(dt, time, environment, heat);
    const int impact = _assembly.impactSubsteps(environment, dt);
    const bool joints = _assembly.updateJoints(dt);
    if (joints) {
        _system.Setup();
        _assembly.enforceRigidNodes();
        refreshLinearSolver();
    }
    if (_assembly.applyLoads(environment, loads, time, dt)) {
        // Contact elements carry a Jacobian, so adding one changes the matrix.
        _system.Setup();
        _assembly.enforceRigidNodes();
        refreshLinearSolver();
    }
    bool rejected = false;
    try {
        if (_assembly.hasFreeSolve()) {
            auto hht = std::dynamic_pointer_cast<chrono::ChTimestepperHHT>(_system.GetTimestepper());
            int contacts = _assembly.activeContactCount();
            const bool contactOnset = _prevContactCount <= 0 && contacts > 0;
            _prevContactCount = contacts;
            const int slices = std::max(1, impact);
            // Predictor-limited first contact increment: CFL / period from
            // impactSubsteps, not 8 h_min discovered by failing Newton.
            const double contactH = std::min(
                dt, std::max(kHhtMinStep, dt / static_cast<double>(slices)));
            const double g = environment.gravity.Length();
            const double mass = std::max(_assembly.mass(), 1.0e-12);
            const double mgd =
                mass * g * std::max(_assembly.contactLengthScale(), 1.0e-4);
            const double bounceHeight =
                g > 1.0e-12 ? _assembly.kineticEnergy() / (mass * g) : 0.0;
            // Quasi-static seating: KE is tiny vs mgδ and nothing is still
            // slamming the plane. Same residual, full dt, no extra cutbacks.
            const bool seated = contacts > 0 &&
                                _assembly.kineticEnergy() <= 0.05 * std::max(mgd, 1.0e-12) &&
                                _assembly.maxApproachSpeed() <= 2.0e-3;
            // Leftover hops under 5 cm: same residual as the last hit, try
            // the user dt and reuse K.
            const bool ringDown = contacts > 0 && bounceHeight <= 0.05 && !seated;
            const bool hardOnset = contactOnset && !ringDown && !seated;
            auto capFreeFlight = [&](double proposed) {
                const double tau = _assembly.minTimeToContact(environment.gravity);
                if (tau >= kNoContactTime)
                    return proposed;
                return std::min(proposed, std::max(kHhtMinStep, tau));
            };
            // Contact loads already occupy Jacobian slots (zero when open).
            // Do not ForceUpdate: sparsity did not change. Reform K numerically
            // on the increment that reaches the plane.
            if (contacts == 0 || seated || (contactOnset && ringDown))
                _localH = dt;
            else if (hardOnset)
                _localH = contactH;
            if (hht) {
                hht->SetMinStepSize(kHhtMinStep);
                hht->SetStepControl(false);
                hht->AcceptTerminatedStep(false);
                hht->SetMaxIters(12);
                hht->SetJacobianUpdateMethod(
                    chrono::ChTimestepperImplicit::JacobianUpdate::AUTOMATIC);
            }
            double remaining = dt;
            double h = _localH > 0.0 ? std::min(_localH, dt) : (contacts > 0 ? contactH : dt);
            h = std::min(std::max(h, kHhtMinStep), dt);
            if (contacts == 0)
                h = capFreeFlight(h);
            int maxIters = seated ? 40 : 12;
            int sameHFails = 0;
            bool reformK = hardOnset || _newtonSuccessStreak < 1;
            auto commitIncrement = [&](double taken, bool keepPlastic) {
                if (keepPlastic)
                    _assembly.commitPlastic();
                else
                    _assembly.discardTrial();
                ++_newtonSuccessStreak;
                ++_hhtCommits;
                remaining -= taken;
                sameHFails = 0;
                maxIters = seated ? 40 : 12;
                reformK = false;
                contacts = _assembly.activeContactCount();
                if (contacts == 0) {
                    _localH = dt;
                    if (remaining > 1.0e-15) {
                        const double tau = _assembly.minTimeToContact(environment.gravity);
                        h = remaining;
                        if (tau < kNoContactTime) {
                            h = std::min(remaining, std::max(kHhtMinStep, tau));
                            if (tau <= h)
                                reformK = true;
                        }
                    }
                } else if (seated) {
                    _localH = dt;
                    if (remaining > 1.0e-15)
                        h = remaining;
                } else {
                    const double grown = std::min(std::max(taken, kHhtMinStep) * 2.0, dt);
                    _localH = grown;
                    if (remaining > 1.0e-15)
                        h = std::min(grown, remaining);
                }
                _prevContactCount = contacts;
                if (hht)
                    hht->SetJacobianUpdateMethod(
                        chrono::ChTimestepperImplicit::JacobianUpdate::AUTOMATIC);
            };
            while (remaining > 1.0e-15) {
                if (cancel && cancel->load()) {
                    rejected = true;
                    break;
                }
                h = std::min(h, remaining);
                if (contacts == 0) {
                    const double tau = _assembly.minTimeToContact(environment.gravity);
                    if (tau < kNoContactTime) {
                        if (tau <= h)
                            reformK = true;
                        h = std::min(h, std::max(kHhtMinStep, tau));
                    }
                }
                if (remaining >= kHhtMinStep)
                    h = std::max(h, kHhtMinStep);
                const bool atMin = h <= kHhtMinStep * 1.01;
                _assembly.setContactStep(h);
                _assembly.captureKinematics();
                const unsigned terminated0 = hht ? hht->GetNumTerminated() : 0;
                const double time0 = _system.GetChTime();
                const double ke0 = _assembly.kineticEnergy();
                const double height0 = _assembly.massHeight();
                const double contact0 = _assembly.contactElasticEnergy();
                const double pen0 = _assembly.maxContactPenetration();
                const double close0 = std::max(0.0, _assembly.maxApproachSpeed());
                if (hht) {
                    hht->SetMaxIters(maxIters);
                    // Last cutback at h_min, or a seated static increment:
                    // keep a finite iterate so the mesh advances, then drop
                    // it if mechanical energy or penetration jumped.
                    hht->AcceptTerminatedStep(atMin || seated);
                    hht->SetJacobianUpdateMethod(
                        (reformK || sameHFails > 0)
                            ? chrono::ChTimestepperImplicit::JacobianUpdate::EVERY_STEP
                            : chrono::ChTimestepperImplicit::JacobianUpdate::AUTOMATIC);
                }
                const bool threw = !tryHhtStep(_system, h);
                const bool finite = !threw && _assembly.finite();
                const bool inverted = finite && _assembly.invertedTetCount() > 0;
                const bool terminated = hht && !threw && hht->GetNumTerminated() > terminated0;
                // Slack is 2% of mgδ, not of current KE. A 2% of KE allowance
                // compounds across last-cutbacks and was launching later hops.
                // Contact ½kδ² is included so a real rebound (spring → KE) is
                // not mistaken for injection.
                const double e0 = ke0 + mass * g * height0 + contact0;
                const double e1 = _assembly.kineticEnergy() + mass * g * _assembly.massHeight() +
                                  _assembly.contactElasticEnergy();
                const bool energyOk = e1 <= e0 + 0.02 * std::max(mgd, 1.0e-12);
                const bool penetrationOk =
                    _assembly.maxContactPenetration() <=
                    pen0 + close0 * h + _assembly.contactLengthScale();
                const bool lastCutback = finite && !inverted && energyOk && penetrationOk &&
                                         (atMin || (seated && terminated));
                const bool converged = finite && !inverted && !terminated;
                if (converged) {
                    commitIncrement(h, true);
                    continue;
                }
                if (lastCutback) {
                    commitIncrement(h, false);
                    continue;
                }
                _assembly.restoreKinematics();
                _assembly.discardTrial();
                _system.SetChTime(time0);
                if (cancel && cancel->load()) {
                    rejected = true;
                    break;
                }
                _newtonSuccessStreak = 0;
                contacts = _assembly.activeContactCount();
                ++sameHFails;
                // Seated: extra Newton on the same h. A slam whose predictor
                // already missed the plane cuts immediately; one more attempt
                // only at h_min before the explicit fallback.
                const int extra = seated ? 3 : (atMin ? 2 : 1);
                if (sameHFails < extra) {
                    maxIters = std::min(40, (seated ? 24 : 12) * (1 << sameHFails));
                    reformK = true;
                    continue;
                }
                sameHFails = 0;
                maxIters = seated ? 40 : 12;
                if (seated) {
                    rejected = true;
                    remaining -= h;
                    _localH = dt;
                    h = remaining > 1.0e-15 ? remaining : dt;
                    reformK = true;
                    continue;
                }
                if (!atMin) {
                    h *= 0.5;
                    if (h < kHhtMinStep)
                        h = kHhtMinStep;
                    _localH = h;
                    reformK = true;
                    continue;
                }
                // Solver threw or the state is non-finite at h_min. One
                // explicit slice, then HHT again — never drain remaining
                // on restored poses (that was the 0.2 s freeze).
                if (contacts > 0) {
                    const bool lumped = _system.GetNumConstraints() == 0 &&
                                        !assemblyHasRotNodes(_assembly);
                    const double hExp = std::min(kHhtMinStep, remaining);
                    _assembly.setContactStep(hExp);
                    _assembly.captureKinematics();
                    const double expT0 = _system.GetChTime();
                    const double expKe0 = _assembly.kineticEnergy();
                    const double expH0 = _assembly.massHeight();
                    const double expC0 = _assembly.contactElasticEnergy();
                    const double expPen0 = _assembly.maxContactPenetration();
                    const double expClose0 = std::max(0.0, _assembly.maxApproachSpeed());
                    const bool expThrew = !tryExplicitStep(hExp, lumped);
                    const bool expOk = !expThrew && _assembly.finite() &&
                                       _assembly.invertedTetCount() == 0;
                    const double expE0 = expKe0 + mass * g * expH0 + expC0;
                    const double expE1 = _assembly.kineticEnergy() +
                                         mass * g * _assembly.massHeight() +
                                         _assembly.contactElasticEnergy();
                    const bool expEnergy = expOk &&
                        expE1 <= expE0 + 0.02 * std::max(mgd, 1.0e-12);
                    const bool expPen = expOk &&
                        _assembly.maxContactPenetration() <=
                            expPen0 + expClose0 * hExp + _assembly.contactLengthScale();
                    if (expEnergy && expPen) {
                        if (hht)
                            hht->SetTime(_system.GetChTime());
                        refreshLinearSolver();
                        commitIncrement(hExp, false);
                        continue;
                    }
                    _assembly.restoreKinematics();
                    _assembly.discardTrial();
                    _system.SetChTime(expT0);
                }
                rejected = true;
                remaining -= std::min(remaining, kHhtMinStep);
                _localH = kHhtMinStep;
                h = kHhtMinStep;
                reformK = true;
            }
            if (!_assembly.finite()) {
                _assembly.restoreKinematics();
                _assembly.discardTrial();
                _assembly.sanitizeShards();
                rejected = true;
            } else {
                _assembly.refreshTetTrials();
            }
            _contactStepFailed = rejected;
            if (hht)
                configureHht(*hht);
        }
    } catch (const std::exception&) {
        _assembly.restoreKinematics();
        _assembly.discardTrial();
        _assembly.sanitizeShards();
        rejected = true;
        ++_newtonRejects;
        if (!_assembly.finite())
            return false;
    }
    if (_assembly.finite())
        _assembly.captureKinematics();
    _assembly.integrateShards(dt, environment);
    const bool volume = _assembly.finite() && _assembly.maintainVolume(dt);
    bool torn = false;
    if (_assembly.finite())
        torn = _assembly.tearOverstretched();
    const bool spanned = (joints || torn) && _assembly.respan();
    if (joints || torn || spanned || volume) {
        // Breaking or promoting a constraint changes the coordinate count and
        // the sparsity pattern. ForceUpdate alone left a stale HHT Jacobian
        // on the free bar after a screw or tree weld let go.
        _system.Setup();
        _assembly.enforceRigidNodes();
        refreshLinearSolver();
    }
    if (rejected)
        ++_newtonRejects;
    if (!rejected)
        _assembly.clearThermalLoadChanged();
    return _assembly.finite();
}

void Simulator::ensureExplicit(bool lumped) {
    if (!_explicit) {
        _explicit = std::make_unique<chrono::ChTimestepperEulerExplicitIIorder>(&_system);
        _explicitLumped = false;
    }
    // Chrono's SetDiagonalLumpingON leaks if called twice; OFF leaves a
    // dangling pointer. Only arm lumping once.
    if (lumped && !_explicitLumped) {
        _explicit->SetDiagonalLumpingON();
        _explicitLumped = true;
    }
}

bool Simulator::tryExplicitStep(double h, bool lumped) {
    ensureExplicit(lumped);
    std::ostringstream sink;
    auto* prev = std::cerr.rdbuf(sink.rdbuf());
    bool threw = false;
    try {
        _explicit->Advance(h);
    } catch (...) {
        threw = true;
    }
    std::cerr.rdbuf(prev);
    return !threw;
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

    bool sharpLoad = false;
    for (std::size_t i = 0; i < scenario.loads.size(); ++i) {
        if (_assembly.adviseLoad(i, scenario.loads, dt).deltaVPerStep > 80.0)
            sharpLoad = true;
    }
    if (sharpLoad) {
        recording.diverged = true;
        recording.duration = 0.0;
        recording.message = "Load step is too sharp for the timestep";
        if (control)
            control->progress.store(1.0f);
        return recording;
    }

    for (long long i = 0; i < totalSteps; ++i) {
        if (control && control->cancel.load()) {
            recording.message = "Cancelled";
            break;
        }

        if (!step(dt, time, scenario.environment, scenario.loads, scenario.heatSources,
                  control ? &control->cancel : nullptr)) {
            recording.diverged = true;
            recording.message = "Solver diverged at t = " + std::to_string(time) + " s";
            break;
        }
        if (control && control->cancel.load()) {
            const double now = std::max(time, _system.GetChTime());
            captureFrame(now);
            time = now;
            recording.message = "Cancelled";
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
    if (sharpLoad)
        recording.diverged = true;
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
