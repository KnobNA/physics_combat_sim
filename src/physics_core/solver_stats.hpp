#pragma once

namespace rcs {

/// One run's incrementation record. Every kept HHT commit is exactly one of
/// `converged`, `terminatedAccepts`, or `explicitSlices`.
struct SolverStats {
    int converged = 0;
    int terminatedAccepts = 0;
    int explicitSlices = 0;
    int cutbacks = 0;
    int rejects = 0;
    double minStep = 0.0;
    double maxEnergyGain = 0.0;  ///< J, max of E(t) - E(t-1) over commits
    double maxPenetration = 0.0;  ///< m
};

}  // namespace rcs
