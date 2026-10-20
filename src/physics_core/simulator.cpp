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
#include <cstdlib>
#include <iostream>
#include <sstream>
#include <thread>

namespace rcs {
namespace {

/// One MKL thread. Pardiso's parallel reductions do not repeat bit for bit,
/// and two identical drops then land a fraction of a millimetre apart.
int solverThreadCount() {
    _putenv("MKL_NUM_THREADS=1");
    _putenv("MKL_DYNAMIC=FALSE");
    return 1;
}

/// Chrono's HHT with a backtracking line search on the Newton update, as
/// Abaqus/Standard and Adams do. Element tangents here are approximate (the
/// Reissner shell's assumed-strain shear, the corotated tet), and without a
/// line search an inexact tangent can lock Newton into a two-iterate cycle
/// that no step cut breaks. The equations, predictor and convergence test
/// are Chrono's; only the length of each update changes, and only when the
/// full update raised the residual.
class LineSearchHht : public chrono::ChTimestepperHHT {
public:
    explicit LineSearchHht(chrono::ChIntegrableIIorder* integrable) : chrono::ChTimestepperHHT(integrable) {}

private:
    void coefficients() {
        _alpha = GetAlpha();
        _gamma = (1.0 - 2.0 * _alpha) / 2.0;
        _beta = std::pow(1.0 - _alpha, 2) / 4.0;
    }

    void InitializeStep() override {
        coefficients();
        integrable->StateSetup(X, V, A);
        Ds.setZero(integrable->GetNumCoordsVelLevel(), integrable);
        Dl.setZero(integrable->GetNumConstraints());
        _xNew.setZero(integrable->GetNumCoordsPosLevel(), integrable);
        _vNew.setZero(integrable->GetNumCoordsVelLevel(), integrable);
        _aNew.setZero(integrable->GetNumCoordsVelLevel(), integrable);
        R.setZero(integrable->GetNumCoordsVelLevel());
        _rOld.setZero(integrable->GetNumCoordsVelLevel());
        Qc.setZero(integrable->GetNumConstraints());
        L.setZero(integrable->GetNumConstraints());
        _lNew.setZero(integrable->GetNumConstraints());
        integrable->StateGather(X, V, T);
        integrable->StateGatherAcceleration(A);
        integrable->StateGatherReactions(L);
    }

    void PrepareStep() override {
        coefficients();
        if (step_control)
            _aNew = A;
        else
            _aNew.setZero(integrable->GetNumCoordsVelLevel(), integrable);
        _vNew = V + _aNew * h;
        _xNew = X + _vNew * h + _aNew * (h * h);
        integrable->LoadResidual_F(_rOld, -_alpha / (1.0 + _alpha));
        integrable->LoadResidual_CqL(_rOld, L, -_alpha / (1.0 + _alpha));
        CalcErrorWeights(A, reltol, abstolS, ewtS);
        _lNew = L;
        CalcErrorWeights(L, reltol, abstolL, ewtL);
        _havePrevious = false;
    }

    void place(const chrono::ChStateDelta& a) {
        _xNew = X + V * h + A * (h * h * (0.5 - _beta)) + a * (h * h * _beta);
        _vNew = V + A * (h * (1.0 - _gamma)) + a * (h * _gamma);
    }

    /// HHT residual at the current iterate, into R and Qc. Returns |R|² + |Qc|².
    double residual() {
        integrable->StateScatter(_xNew, _vNew, T + h, chrono::UpdateFlags::UPDATE_ALL_NO_VISUAL);
        R = _rOld;
        Qc.setZero();
        integrable->LoadResidual_F(R, 1.0);
        integrable->LoadResidual_CqL(R, _lNew, 1.0);
        integrable->LoadResidual_Mv(R, _aNew, -1.0 / (1.0 + _alpha));
        integrable->LoadConstraint_C(Qc, 1.0 / (_beta * h * h), 1.0 / (_gamma * h), Qc_do_clamp, Qc_clamping);
        return R.squaredNorm() + Qc.squaredNorm();
    }

    void Increment() override {
        double merit = residual();
        // The last full update made things worse: walk back along it.
        if (_havePrevious && merit > _previousMerit) {
            double step = 1.0;
            double bestStep = 1.0;
            double bestMerit = merit;
            for (int cut = 0; cut < kCuts; ++cut) {
                step *= 0.5;
                _aNew = _previousA + _previousDs * step;
                _lNew = _previousL + _previousDl * step;
                place(_aNew);
                const double trial = residual();
                if (trial < bestMerit) {
                    bestMerit = trial;
                    bestStep = step;
                }
                if (trial < _previousMerit)
                    break;
            }
            if (bestStep != step) {
                _aNew = _previousA + _previousDs * bestStep;
                _lNew = _previousL + _previousDl * bestStep;
                place(_aNew);
                bestMerit = residual();
            }
            merit = bestMerit;
        }

        integrable->StateSolveCorrection(Ds, Dl, R, Qc, 1.0 / (1.0 + _alpha), -h * _gamma, -h * h * _beta, _xNew,
                                         _vNew, T + h, false, chrono::UpdateFlags::UPDATE_ALL_NO_VISUAL, call_setup,
                                         call_analyze);
        _previousA = _aNew;
        _previousL = _lNew;
        _previousDs = Ds;
        _previousDl = Dl;
        _previousMerit = merit;
        _havePrevious = true;
        if (verbose) {
            Eigen::Index imax = 0;
            Ds.cwiseAbs().maxCoeff(&imax);
            Eigen::Index rmax = 0;
            R.cwiseAbs().maxCoeff(&rmax);
            std::cout << "    dsMax idx " << imax << " ds " << Ds(imax) << " R " << R(imax) << " | rMax idx " << rmax
                      << " R " << R(rmax) << " ds " << Ds(rmax) << " n " << Ds.size() << "\n";
        }
        _lNew += Dl;
        _aNew += Ds;
        place(_aNew);
    }

    void ResetStep() override {
        integrable->StateScatter(X, V, T, chrono::UpdateFlags::UPDATE_ALL_NO_VISUAL);
        _rOld.setZero();
        _aNew.setZero(integrable->GetNumCoordsVelLevel(), integrable);
        _havePrevious = false;
    }

    void AcceptStep() override {
        X = _xNew;
        V = _vNew;
        A = _aNew;
        L = _lNew;
    }

    void FinalizeStep() override {
        integrable->StateScatter(X, V, T, chrono::UpdateFlags::UPDATE_ALL);
        integrable->StateScatterAcceleration(A);
        integrable->StateScatterReactions(L);
        integrable->StateOnEndStep(T);
    }

    static constexpr int kCuts = 5;
    double _alpha = -0.05;
    double _gamma = 0.55;
    double _beta = 0.275625;
    chrono::ChState _xNew;
    chrono::ChStateDelta _vNew;
    chrono::ChStateDelta _aNew;
    chrono::ChVectorDynamic<> _lNew;
    chrono::ChVectorDynamic<> _rOld;
    chrono::ChStateDelta _previousA;
    chrono::ChStateDelta _previousDs;
    chrono::ChVectorDynamic<> _previousL;
    chrono::ChVectorDynamic<> _previousDl;
    double _previousMerit = 0.0;
    bool _havePrevious = false;
};

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

/// Velocities and element internal forces (elastic, plastic, damping) at the
/// current state, in system velocity coordinates. Node loads are left out:
/// gravity and contact enter the balance through their exact potentials.
struct InternalState {
    chrono::ChStateDelta v;
    chrono::ChVectorDynamic<> f;
    double kinetic = 0.0;  ///< ½ vᵀ M v with the integrator's own mass matrix
};

InternalState internalState(chrono::ChSystem& system, const Assembly& assembly) {
    InternalState s;
    const auto nx = system.GetNumCoordsPosLevel();
    const auto nv = system.GetNumCoordsVelLevel();
    chrono::ChState x(nx, &system);
    s.v.setZero(nv, &system);
    double t = 0.0;
    system.StateGather(x, s.v, t);
    s.f.setZero(nv);
    for (const auto& part : assembly.parts()) {
        if (!part.mesh())
            continue;
        for (const auto& element : part.mesh()->GetElements())
            element->EleIntLoadResidual_F(s.f, 1.0);
    }
    // Joint reactions: a weld holding back thermal growth carries a large
    // force that only cancels the element force when both are counted.
    if (system.GetNumConstraints() > 0) {
        chrono::ChVectorDynamic<> l;
        l.setZero(system.GetNumConstraints());
        system.StateGatherReactions(l);
        system.LoadResidual_CqL(s.f, l, 1.0);
    }
    chrono::ChVectorDynamic<> mv;
    mv.setZero(nv);
    system.LoadResidual_Mv(mv, s.v, 1.0);
    s.kinetic = 0.5 * s.v.dot(mv);
    return s;
}

/// Work done on the nodes by internal forces over one increment, trapezoid
/// rule, the way an implicit code accumulates internal energy. Positive when
/// stored strain energy is released into motion.
double internalWork(const InternalState& a, const InternalState& b, double h) {
    if (a.f.size() != b.f.size() || a.v.size() != b.v.size() || a.f.size() != a.v.size())
        return 0.0;
    return 0.25 * h * (a.f + b.f).dot(a.v + b.v);
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
    auto hht = chrono_types::make_shared<LineSearchHht>(&_system);
    configureHht(*hht);
    _system.SetTimestepper(hht);

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
    bool advanced = false;
    bool incomplete = false;
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
            // A millimetre hover after a bounce is not free flight. MixedJoints
            // seating left nC=0 then took h=dt, Newton ran gap ~ 1e23, and
            // HHT threw. Keep contact-sized h while any node is within the
            // penalty gap cap of the floor.
            const bool nearFloor =
                _assembly.touchingFloor(environment, kContactGapCap);
            // Quasi-static seating: KE is tiny vs mgδ and nothing is still
            // slamming the plane. Same residual, full dt, no extra cutbacks.
            // After a 1 mm bounce nC drops to 0 while still on the floor;
            // treating that as seated avoids EVERY_ITERATION on 121 nodes.
            const bool seated = (contacts > 0 || nearFloor) &&
                                _assembly.kineticEnergy() <= 0.05 * std::max(mgd, 1.0e-12) &&
                                _assembly.maxApproachSpeed() <= 2.0e-3;
            const bool hardOnset = contactOnset && !seated;
            // Stop an increment at first touch only when the approach would
            // overshoot the contact skin within that increment. A body
            // settling a fraction of a millimetre is an ordinary implicit
            // contact increment, as in an implicit FE code; capping it to
            // the time-to-contact of a hovering skin drives h to h_min.
            auto capFreeFlight = [&](double proposed) {
                const double tau = _assembly.minTimeToContact(environment.gravity);
                if (tau >= proposed)
                    return proposed;
                // Quasi-static approach only (the seated test): an impact
                // integrated inside one long increment lands at an arbitrary
                // point of the contact period, and the rebound depends on it.
                const double speed =
                    std::sqrt(2.0 * _assembly.kineticEnergy() / mass) + g * proposed;
                const double ke = 0.5 * mass * speed * speed;
                if (ke <= 0.05 * std::max(mgd, 1.0e-12) &&
                    speed * proposed <= _assembly.contactLengthScale())
                    return proposed;
                return std::min(proposed, std::max(kHhtMinStep, tau));
            };
            // Contact loads already occupy Jacobian slots (zero when open).
            // Do not ForceUpdate: sparsity did not change. Reform K numerically
            // on the increment that reaches the plane.
            if (seated)
                _localH = dt;
            else if (contacts == 0 && !nearFloor)
                _localH = dt;
            else if (hardOnset) {
                const double halfPeriod = _assembly.contactHalfPeriod();
                _localH = halfPeriod < kNoContactTime ? std::min(dt, std::max(kHhtMinStep, halfPeriod)) : contactH;
            }
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
            auto commitIncrement = [&](double taken, bool keepPlastic, int kind, double energyBefore) {
                if (keepPlastic)
                    _assembly.commitPlastic();
                else
                    _assembly.discardTrial();
                _assembly.commitContactApproach();
                ++_newtonSuccessStreak;
                ++_hhtCommits;
                advanced = true;
                if (kind == 0)
                    ++_stats.converged;
                else if (kind == 1)
                    ++_stats.terminatedAccepts;
                else
                    ++_stats.explicitSlices;
                if (_stats.minStep <= 0.0 || taken < _stats.minStep)
                    _stats.minStep = taken;
                const double gained =
                    (_assembly.kineticEnergy() + mass * g * _assembly.massHeight() +
                     _assembly.contactElasticEnergy()) -
                    energyBefore;
                _stats.maxEnergyGain = std::max(_stats.maxEnergyGain, gained);
                _stats.maxPenetration =
                    std::max(_stats.maxPenetration, _assembly.maxContactPenetration());
                remaining -= taken;
                sameHFails = 0;
                maxIters = seated ? 40 : 12;
                reformK = false;
                const int contactsBefore = contacts;
                contacts = _assembly.activeContactCount();
                if (contactsBefore == 0 && contacts > 0 && !seated) {
                    // First touch: size the next increment from the contact
                    // period, not from how much of the step happened to remain
                    // before touch, so nearby drops integrate the same impact.
                    // Skip when KE is already quasi-static: MixedJoints seating
                    // chattered off the floor, then this reset drained at 25 µs.
                    const double keNow = _assembly.kineticEnergy();
                    const bool impact = keNow > 0.05 * std::max(mgd, 1.0e-12);
                    const double halfPeriod = _assembly.contactHalfPeriod();
                    const double onset = halfPeriod < kNoContactTime ? halfPeriod : taken;
                    if (impact)
                        _localH = std::min(dt, std::max(kHhtMinStep, onset));
                    if (remaining > 1.0e-15)
                        h = std::min(_localH, remaining);
                } else if (contacts == 0) {
                    if (_assembly.touchingFloor(environment, kContactGapCap)) {
                        const double grown = std::min(std::max(taken, kHhtMinStep) * 2.0, dt);
                        _localH = grown;
                        if (remaining > 1.0e-15)
                            h = std::min(grown, remaining);
                    } else {
                        _localH = dt;
                        if (remaining > 1.0e-15)
                            h = capFreeFlight(remaining);
                    }
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
                // A sliver smaller than h_min is leftover rounding, not a
                // physical increment. Factoring it threw on the hollow box
                // (h = 7.84e-7 at t = 0.668 s).
                if (remaining < kHhtMinStep)
                    break;
                h = std::min(h, remaining);
                if (contacts == 0)
                    h = capFreeFlight(h);
                const double tau = _assembly.minTimeToContact(environment.gravity);
                if (remaining >= kHhtMinStep)
                    h = std::max(h, kHhtMinStep);
                const bool atMin = h <= kHhtMinStep * 1.01;
                _assembly.setContactStep(h);
                _assembly.captureKinematics();
                const unsigned terminated0 = hht ? hht->GetNumTerminated() : 0;
                const double time0 = _system.GetChTime();
                const double height0 = _assembly.massHeight();
                const double contact0 = _assembly.contactElasticEnergy();
                const double pen0 = _assembly.maxContactPenetration();
                const double close0 = std::max(0.0, _assembly.maxApproachSpeed());
                if (hht) {
                    // Contact status changes inside this increment. Like an
                    // implicit FE contact search, re-form the tangent at every
                    // iterate so the active set and its stiffness agree, and
                    // allow the extra iterations an active-set change needs.
                    const bool crossing = !seated && (contacts > 0 || tau <= h * 1.01);
                    hht->SetMaxIters(crossing ? std::max(maxIters, 30) : maxIters);
                    hht->AcceptTerminatedStep(false);
                    hht->SetJacobianUpdateMethod(
                        (crossing || reformK || sameHFails > 0)
                            ? chrono::ChTimestepperImplicit::JacobianUpdate::EVERY_ITERATION
                            : chrono::ChTimestepperImplicit::JacobianUpdate::AUTOMATIC);
                }
                const InternalState internal0 = internalState(_system, _assembly);
                const bool threw = !tryHhtStep(_system, h);
                const bool finite = !threw && _assembly.finite();
                const bool inverted = finite && _assembly.invertedTetCount() > 0;
                const bool terminated = hht && !threw && hht->GetNumTerminated() > terminated0;
                // Energy balance: KE + gravity + contact potential may change
                // only by the work internal forces did (strain released or
                // stored, damping and plasticity removed). Slack is 2% of mgδ
                // per increment. Without the internal term a panel springing
                // back from its sag reads as injected energy.
                const InternalState internal1 = finite ? internalState(_system, _assembly) : internal0;
                const double wInternal = finite ? internalWork(internal0, internal1, h) : 0.0;
                const double e0 = internal0.kinetic + mass * g * height0 + contact0 + wInternal;
                const double e1 = internal1.kinetic + mass * g * _assembly.massHeight() +
                                  _assembly.contactElasticEnergy();
                // Slack is 2% of the larger of the gravity scale and the
                // energy actually moving this increment, so a run without
                // gravity still has a meaningful tolerance.
                const double energyScale = std::max(
                    {mgd, internal0.kinetic, internal1.kinetic, std::abs(wInternal), 1.0e-9});
                // HHT with alpha < 0 is dissipative on its own; penalty contact
                // is what can inject energy. Without contact the f.v estimate
                // also misreads a thermal-strain jump that the integrator
                // damps out within one increment.
                const bool contactIncrement =
                    contacts > 0 || contact0 > 0.0 || _assembly.activeContactCount() > 0;
                const bool energyOk = !contactIncrement || e1 <= e0 + 0.02 * energyScale;
                const bool penetrationOk =
                    _assembly.maxContactPenetration() <=
                    pen0 + close0 * h + _assembly.contactLengthScale();
                const bool lastCutback = finite && !inverted && energyOk && penetrationOk &&
                                         (atMin || (seated && terminated));
                const bool converged = finite && !inverted && !terminated;
                if (converged && energyOk && penetrationOk) {
                    commitIncrement(h, true, 0, e0);
                    continue;
                }
                if (lastCutback) {
                    commitIncrement(h, false, 1, e0);
                    continue;
                }
                _assembly.restoreKinematics();
                _assembly.discardTrial();
                _system.SetChTime(time0);
                // Loads and constraint Jacobians are cached from the last
                // Newton iterate; HHT reads them before its first scatter.
                _system.Update(chrono::UpdateFlags::UPDATE_ALL_NO_VISUAL);
                if (cancel && cancel->load()) {
                    rejected = true;
                    break;
                }
                _newtonSuccessStreak = 0;
                contacts = _assembly.activeContactCount();
                ++sameHFails;
                ++_stats.cutbacks;
                ++_stats.rejects;
                // A Newton failure may converge with more iterations on the
                // same h (seated: a few, otherwise one more). A converged step
                // that failed the energy or penetration check gives the same
                // answer again, so it cuts h at once. Rejected time is never
                // skipped: every second of a run is a committed increment.
                const bool newtonFailed = threw || terminated || !finite;
                const int extra = newtonFailed ? (seated ? 3 : (atMin ? 2 : 1)) : 0;
                if (sameHFails < extra) {
                    maxIters = std::min(40, (seated ? 24 : 12) * (1 << sameHFails));
                    reformK = true;
                    continue;
                }
                sameHFails = 0;
                maxIters = seated ? 40 : 12;
                if (!atMin) {
                    h *= 0.5;
                    if (h < kHhtMinStep)
                        h = kHhtMinStep;
                    _localH = h;
                    reformK = true;
                    continue;
                }
                // One explicit slice, then stop. Walking the rest of the user
                // step at h_min on a restored pose repeats the same rejected
                // factorization thousands of times.
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
                        commitIncrement(hExp, false, 2, expE0);
                        continue;
                    }
                    _assembly.restoreKinematics();
                    _assembly.discardTrial();
                    _system.SetChTime(expT0);
                    _system.Update(chrono::UpdateFlags::UPDATE_ALL_NO_VISUAL);
                }
                rejected = true;
                incomplete = true;
                break;
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
    // A finite mesh that could not accept all of this user step has not
    // advanced to t + dt. Reporting success would silently drop the rest.
    if (_assembly.hasFreeSolve() && (!advanced || incomplete))
        return false;
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
