# Physics unit inventory

Each exported calculation has a GoogleTest whose name contains the function name.
Chrono orchestrators stay on `combat_sim_physics_smoke` and are not unit-tested.

## Units (`src/physics_core/physics_units.hpp`)

| Function | TEST |
|---|---|
| `vonMisesStress` | `vonMisesStress.UniaxialTensionEqualsAbsSigma`, `vonMisesStress.ChVectorDynamicUniaxial` |
| `vonMisesPlaneStress` | `vonMisesPlaneStress.UniaxialEqualsAbsSigma` |
| `equivalentPlastic` | `equivalentPlastic.UniaxialDeviators` |
| `projectRankine` | `projectRankine.ClipsOnlyOverLimitPrincipal`, `projectRankine.UntouchedInsideLimits`, `projectRankine.ChVectorDynamicWrapper` |
| `brittleShouldOpen` | `brittleShouldOpen.TensionBelowRankineIsFalse`, `brittleShouldOpen.TensionAboveButEnergyBelowGriffithIsFalse`, `brittleShouldOpen.BothGatesOpen`, `brittleShouldOpen.TemperedUsesTenthGc`, `brittleShouldOpen.DuctileNeverOpens`, `brittleShouldOpen.NoToughnessFallsBackToStrain` |
| `findRoot` | `findRoot.PathCompressedToRepresentative` |
| `unite` | `unite.AlreadyUnitedIsFalse` |
| `distanceToBead` | `distanceToBead.PointOnEdgeIsZero`, `distanceToBead.CentreOfUnitSquare` |
| `beadPerimeter` | `beadPerimeter.UnitSquare` |
| `closedBeadArea` | `closedBeadArea.UnitSquare` |
| `closestOnBead` | `closestOnBead.ProjectsOntoEdge` |
| `closestPointOnTriangle` | `closestPointOnTriangle.ProjectsToInterior` |
| `barycentric` | `barycentric.CentroidIsThirds` |
| `tetContainsPoint` | `tetContainsPoint.CentreInsideFarOutside` |
| `stableGains` | `stableGains.ClampsStiffnessAndDamping`, `stableGains.ZeroMassOrDt` |
| `restGains` | `restGains.OverdampedCap`, `restGains.ZeroMassOrDt` |
| `faceArea` | `faceArea.RightTriangleHalf` |
| `meshVolume` | `meshVolume.UnitCube` |
| `surfaceVolume` | `surfaceVolume.UnitCubeAbsolute` |
| `closedEnough` | `closedEnough.ClosedCubeAndOpenTriangle` |
| `budgetEdge` | `budgetEdge.CubeWithDefaultBudget` |
| `worstTetQuality` | `worstTetQuality.RegularTetNearPointOneOneSeven`, `worstTetQuality.SliverNearZero` |

## Already-public inlines

| Function | TEST |
|---|---|
| `mapThermalAxis` | `mapThermalAxis.RoomTemperatureIdentity`, `mapThermalAxis.PastMeltScales` |
| `sampleThermalKnots` | `sampleThermalKnots.InterpolatesMidSpan`, `sampleThermalKnots.ClampsBelowFirstKnot` |
| `thermalStiffnessRetention` | `thermalStiffnessRetention.AluminumRoomNearUnity`, `thermalStiffnessRetention.AluminumAboveSolidusNearFloor` |
| `thermalStrengthRetention` | `thermalStrengthRetention.AluminumRoomNearUnity` |
| `thermalDuctilityRetention` | `thermalDuctilityRetention.AluminumRoomIsOne` |
| `kThermalFloor` | `kThermalFloor.Residual` |
| `rankineTensionLimit` | `rankineTensionLimit.DuctileIsZero`, `rankineTensionLimit.GlassUsesFractureStrength` |
| `rankineCompressionLimit` | `rankineCompressionLimit.GlassUsesExplicitCrush`, `rankineCompressionLimit.DerivedFifteenTimesTension` |
| `aluminum6061` | `aluminum6061.DefaultNameAndDensity` |
| `mildSteel` | `mildSteel.MeltTemperature` |
| `ar500` | `ar500.YieldStress` |
| `glass` | `glass.BrittleFracture` |
| `temperedGlass` | `temperedGlass.ShatterToDust` |
| `rubber` | `rubber.CompliantFracture` |
| `materialPresets` | `materialPresets.LibraryOrder` |
| `presetIndex` | `presetIndex.MatchesLibraryOrder` |
| `loadRampBlend` | `loadRampBlend.LinearMidpoint`, `loadRampBlend.ExponentialEndpoints`, `loadRampBlend.LogarithmicEndpoints`, `loadRampBlend.MidRampLogGreaterThanLinearGreaterThanExp`, `loadRampBlend.PolynomialDegreeTwo` |
| `loadRampMaxSlope` | `loadRampMaxSlope.KnownShapes` |
| `loadMagnitude` | `loadMagnitude.ZeroOutsideWindow`, `loadMagnitude.LinearRampThenHold`, `loadMagnitude.DisabledOrZeroDuration` |
| `loadMaxStepChange` | `loadMaxStepChange.StepAndRamp` |
| `Environment::gridStep` | `gridStep.ArenaSpans` |
| `packHalf` | `packHalf.MillimetreRoundTrip` |
| `unpackHalf` | `unpackHalf.MillimetreScalePositions` |
| `Joint::brokenCount` | `brokenCount.CountsFailedLinks` |
| `Joint::severed` | `severed.EmptyOrAllBroken` |
| `Joint::peakUtilisation` | `peakUtilisation.IgnoresBrokenLinks` |

Smoke-group fidelity is in `tests/unit/smoke_inventory.md`.

## Smoke-only (not unit-tested)

These are Chrono/HHT orchestration loops. Known I/O lives in `combat_sim_physics_smoke`.

- `Simulator::step`
- `Simulator::run`
- `Assembly::build`
- `Assembly::applyContact` / nodal contact force assembly
- `Part::build`
- `Part::maintainVolume`
- `Part::tear`
- `Part::releaseMoltenTets`
- `HeatField::advance`
- `NodeContact::ComputeQ`
- `ThermalTet4::ComputeInternalForces` / `ComputeKRMmatricesGlobal`
- `ensurePartMesh` / tetgen fill
- ImGui, OpenGL, SDL
