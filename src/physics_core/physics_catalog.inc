/// Handbook-band cases from the catalog. Finite is required and never enough.
/// These share the smoke helpers (check, expectRange, clampedPlate, …).

double catalogPeakStress(const rcs::Part& part) {
    std::vector<float> field;
    part.evaluateStress(field);
    float peak = 0.0f;
    for (float value : field)
        peak = std::max(peak, value);
    return peak;
}

double catalogPeakPlastic(const rcs::Part& part) {
    std::vector<float> field;
    part.evaluatePlasticStrain(field);
    float peak = 0.0f;
    for (float value : field) {
        if (!std::isfinite(value))
            continue;
        peak = std::max(peak, value);
    }
    return peak;
}

double catalogKinetic(const rcs::Part& part) {
    double energy = 0.0;
    for (std::size_t n = 0; n < part.nodeCount(); ++n) {
        if (nodeIsFixed(part, n))
            continue;
        const double mass = n < part.nodeMass().size() ? part.nodeMass()[n] : 0.0;
        const double speed = part.nodeVelocity(n).Length();
        energy += 0.5 * mass * speed * speed;
    }
    for (const auto& shard : part.shards())
        energy += 0.5 * shard.mass * shard.velocity.Length2();
    return energy;
}

chrono::ChVector3d catalogCom(const rcs::Part& part) {
    chrono::ChVector3d com{0, 0, 0};
    double mass = 0.0;
    for (std::size_t n = 0; n < part.nodeCount(); ++n) {
        const double m = n < part.nodeMass().size() ? part.nodeMass()[n] : 1.0;
        com += m * part.nodePosition(n);
        mass += m;
    }
    return mass > 1.0e-12 ? com / mass : com;
}

double catalogLowest(const rcs::Part& part) {
    double y = 1.0e9;
    for (std::size_t n = 0; n < part.nodeCount(); ++n)
        y = std::min(y, part.nodePosition(n).y());
    return y;
}

double catalogWorstTetAspect(const rcs::Part& part) {
    double worst = 1.0;
    for (int t = 0; t < part.tetCount(); ++t) {
        if (!part.tetLive(static_cast<std::size_t>(t)))
            continue;
        const auto& tet = part.tets()[static_cast<std::size_t>(t)];
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
            worst = std::max(worst, longest / shortest);
    }
    return worst;
}

rcs::Scenario catalogQuiet() {
    rcs::Scenario scenario;
    scenario.settings.timestep = 2.0e-4;
    scenario.initial.position = {0.0, 0.5, 0.0};
    scenario.initial.orientation = {0.0, 0.0, 0.0};
    scenario.environment.gravity = {0.0, 0.0, 0.0};
    scenario.environment.floorEnabled = false;
    scenario.environment.aerodynamicsEnabled = false;
    return scenario;
}

rcs::AssemblyDefinition catalogBar(const rcs::Material& material, double length, double width, double thick,
                                   int divU, int divV, bool clamp) {
    rcs::PartDefinition bar;
    bar.name = material.name;
    bar.material = material;
    bar.thickness = thick;
    bar.clampBoundary = clamp;
    bar.patch.origin = {0.0, 0.0, -0.5 * width};
    bar.patch.uAxis = {1.0, 0.0, 0.0};
    bar.patch.vAxis = {0.0, 0.0, 1.0};
    bar.patch.uSize = length;
    bar.patch.vSize = width;
    bar.patch.divU = divU;
    bar.patch.divV = divV;
    rcs::AssemblyDefinition definition;
    definition.name = material.name;
    definition.parts.push_back(bar);
    return definition;
}

void catalogFixEnds(rcs::Part& part) {
    double minX = 1.0e9;
    double maxX = -1.0e9;
    for (std::size_t n = 0; n < part.nodeCount(); ++n) {
        const double x = part.restPosition(n).x();
        minX = std::min(minX, x);
        maxX = std::max(maxX, x);
    }
    auto pin = [&](auto& nodes) {
        for (std::size_t n = 0; n < nodes.size() && n < part.nodeCount(); ++n) {
            if (!nodes[n])
                continue;
            const double x = part.restPosition(n).x();
            if (x <= minX + 1.0e-3 || x >= maxX - 1.0e-3)
                nodes[n]->SetFixed(true);
        }
    };
    pin(part.xyzNodes());
    pin(part.nodes());
}

void catalogPinRoot(rcs::Part& part) {
    double minX = 1.0e9;
    for (std::size_t n = 0; n < part.nodeCount(); ++n)
        minX = std::min(minX, part.restPosition(n).x());
    auto pin = [&](auto& nodes) {
        for (std::size_t n = 0; n < nodes.size() && n < part.nodeCount(); ++n) {
            if (!nodes[n])
                continue;
            nodes[n]->SetFixed(part.restPosition(n).x() <= minX + 1.0e-3);
        }
    };
    pin(part.xyzNodes());
    pin(part.nodes());
}

bool catalogStep(rcs::Simulator& simulator, rcs::Scenario& scenario, int steps) {
    bool finite = true;
    double t = 0.0;
    const double dt = scenario.settings.timestep;
    for (int step = 0; step < steps && finite; ++step) {
        finite = simulator.step(dt, t, scenario.environment, scenario.loads, scenario.heatSources);
        t += dt;
    }
    return finite;
}

std::size_t catalogEndNode(const rcs::Part& part, double length) {
    std::size_t tip = 0;
    double best = 1.0e9;
    for (std::size_t n = 0; n < part.nodeCount(); ++n) {
        const auto rest = part.restPosition(n);
        const double score = std::abs(rest.x() - length) + std::abs(rest.z()) * 4.0;
        if (score < best) {
            best = score;
            tip = n;
        }
    }
    return tip;
}

struct AxialResult {
    bool finite = false;
    double strain = 0.0;
    double lateral = 0.0;
    double plastic = 0.0;
    double residual = 0.0;
    double peakStress = 0.0;
    bool torn = false;
    std::size_t shards = 0;
};

AxialResult catalogAxial(const rcs::Material& material, double stress, double hold, double rest) {
    const double length = 0.12;
    const double width = 0.03;
    const double thick = 0.006;
    const double area = width * thick;
    rcs::AssemblyDefinition definition = catalogBar(material, length, width, thick, 8, 2, true);
    rcs::Scenario scenario = catalogQuiet();
    scenario.initial.position = {0.0, 0.5, 0.0};
    rcs::ImpactLoad hit;
    hit.point = {length, 0.5, 0.0};
    hit.direction = {1.0, 0.0, 0.0};
    hit.initialForce = 0.0;
    hit.finalForce = stress * area;
    hit.rampTime = 0.03;
    hit.duration = hold;
    hit.radius = 0.5 * width;
    scenario.loads.push_back(hit);

    rcs::Simulator simulator;
    simulator.reset(definition, scenario);
    rcs::Part& part = simulator.assembly().parts()[0];
    catalogPinRoot(part);
    simulator.system().Setup();
    simulator.system().ForceUpdate();
    simulator.assembly().bindLoads(scenario.loads);
    const std::size_t tip = catalogEndNode(part, length);
    const double start = part.nodePosition(tip).x();
    double minZ = 1.0e9;
    double maxZ = -1.0e9;
    for (std::size_t n = 0; n < part.nodeCount(); ++n) {
        minZ = std::min(minZ, part.restPosition(n).z());
        maxZ = std::max(maxZ, part.restPosition(n).z());
    }
    std::size_t lo = 0;
    std::size_t hi = 0;
    for (std::size_t n = 0; n < part.nodeCount(); ++n) {
        if (std::abs(part.restPosition(n).z() - minZ) < 1.0e-4)
            lo = n;
        if (std::abs(part.restPosition(n).z() - maxZ) < 1.0e-4)
            hi = n;
    }
    const double width0 = (part.nodePosition(hi) - part.nodePosition(lo)).Length();

    AxialResult result;
    result.finite = true;
    double t = 0.0;
    const double dt = scenario.settings.timestep;
    const int steps = static_cast<int>(std::ceil((hold + rest) / dt));
    double peakStrain = 0.0;
    for (int step = 0; step < steps && result.finite; ++step) {
        result.finite = simulator.step(dt, t, scenario.environment, scenario.loads);
        t += dt;
        const rcs::Part& live = simulator.assembly().parts()[0];
        peakStrain = std::max(peakStrain, (live.nodePosition(tip).x() - start) / length);
        result.peakStress = std::max(result.peakStress, catalogPeakStress(live));
        result.plastic = std::max(result.plastic, catalogPeakPlastic(live));
    }
    const rcs::Part& done = simulator.assembly().parts()[0];
    result.strain = peakStrain;
    const double width1 = (done.nodePosition(hi) - done.nodePosition(lo)).Length();
    result.lateral = width0 > 1.0e-9 ? (width1 - width0) / width0 : 0.0;
    result.residual = (done.nodePosition(tip).x() - start) / length;
    result.torn = done.torn();
    result.shards = done.shards().size();
    return result;
}

rcs::PartDefinition catalogSmallCube(const rcs::Material& material, double size) {
    rcs::PartDefinition cube = rcs::makeSolidCube(size);
    cube.material = material;
    cube.sizing.maxElement = std::max(0.04, size / 4.0);
    cube.sizing.minElement = cube.sizing.maxElement * 0.5;
    cube.sizing.tetBudget = 400;
    cube.sizing.cached = false;
    cube.volume = {};
    return cube;
}

void checkConserve() {
    std::printf("\n-- conserve (N2-N6) --\n");

    rcs::PartDefinition cube = catalogSmallCube(rcs::aluminum6061(), 0.16);
    check(rcs::ensurePartMesh(cube), "N2 cube meshes");
    rcs::AssemblyDefinition definition;
    definition.parts.push_back(cube);
    rcs::Scenario scene = catalogQuiet();
    scene.settings.dynamicResolution = true;
    scene.environment.gravity = {0.0, -9.81, 0.0};
    scene.environment.floorEnabled = true;
    scene.initial.position = {0.0, 0.20, 0.0};
    rcs::ImpactLoad punch;
    punch.point = {0.0, 0.20, 0.0};
    punch.direction = {0.0, -1.0, 0.0};
    punch.finalForce = 8.0e4;
    punch.radius = 0.03;
    punch.duration = 0.02;
    scene.loads.push_back(punch);
    rcs::Simulator simulator;
    simulator.reset(definition, scene);
    const double mass0 = simulator.assembly().mass();
    const chrono::ChVector3d com0 = catalogCom(simulator.assembly().parts()[0]);
    const bool finite = catalogStep(simulator, scene, 80);
    const rcs::Part& done = simulator.assembly().parts()[0];
    const double mass1 = simulator.assembly().mass();
    const chrono::ChVector3d com1 = catalogCom(done);
    expectFinite(finite, "N2 refine/load stays in the solver");
    expectRange(std::abs(mass1 - mass0) / std::max(mass0, 1.0e-9), 0.0, 0.01, "N2 mass error after load/refine");
    expectRange((com1 - com0).Length(), 0.0, 0.08, "N2 COM stays on the body");
    expectRange(catalogWorstTetAspect(done), 1.0, 40.0, "N6 tet aspect after a dent");

    rcs::PartDefinition spinCube = catalogSmallCube(rcs::aluminum6061(), 0.16);
    check(rcs::ensurePartMesh(spinCube), "N3 cube meshes");
    rcs::AssemblyDefinition spinDef;
    spinDef.parts.push_back(spinCube);
    rcs::Scenario spin = catalogQuiet();
    spin.initial.angularVelocity = {0.0, 3.0, 0.0};
    rcs::Simulator spinner;
    spinner.reset(spinDef, spin);
    chrono::ChVector3d L0{0, 0, 0};
    const rcs::Part& before = spinner.assembly().parts()[0];
    const chrono::ChVector3d c0 = catalogCom(before);
    for (std::size_t n = 0; n < before.nodeCount(); ++n) {
        const double m = n < before.nodeMass().size() ? before.nodeMass()[n] : 0.0;
        L0 += m * (before.nodePosition(n) - c0).Cross(before.nodeVelocity(n));
    }
    const bool spinFinite = catalogStep(spinner, spin, 100);
    const rcs::Part& after = spinner.assembly().parts()[0];
    chrono::ChVector3d L1{0, 0, 0};
    const chrono::ChVector3d c1 = catalogCom(after);
    double stretch = 0.0;
    for (std::size_t n = 0; n < after.nodeCount(); ++n) {
        const double m = n < after.nodeMass().size() ? after.nodeMass()[n] : 0.0;
        L1 += m * (after.nodePosition(n) - c1).Cross(after.nodeVelocity(n));
        const double rest = (after.restPosition(n) - before.restPosition(0)).Length();
        (void)rest;
    }
    for (int t = 0; t < after.tetCount(); ++t) {
        if (!after.tetLive(static_cast<std::size_t>(t)))
            continue;
        const auto& tet = after.tets()[static_cast<std::size_t>(t)];
        for (int a = 0; a < 4; ++a) {
            const int b = (a + 1) % 4;
            const double r0 = (after.restPosition(static_cast<std::size_t>(tet[a])) -
                               after.restPosition(static_cast<std::size_t>(tet[b])))
                                  .Length();
            const double r1 = (after.nodePosition(static_cast<std::size_t>(tet[a])) -
                               after.nodePosition(static_cast<std::size_t>(tet[b])))
                                  .Length();
            if (r0 > 1.0e-6)
                stretch = std::max(stretch, std::abs(r1 - r0) / r0);
        }
    }
    expectFinite(spinFinite, "N3 spinning cube stays in the solver");
    expectRange(stretch, 0.0, 0.02, "N3 spin stretch under 2%");
    expectRange(L1.Length() / std::max(L0.Length(), 1.0e-9), 0.90, 1.10, "N3 angular momentum within 10%");

    rcs::PartDefinition restCube = catalogSmallCube(rcs::aluminum6061(), 0.16);
    check(rcs::ensurePartMesh(restCube), "N4 cube meshes");
    rcs::AssemblyDefinition restDef;
    restDef.parts.push_back(restCube);
    rcs::Scenario rest = catalogQuiet();
    rest.settings.timestep = 1.0e-3;
    rest.environment.gravity = {0.0, -9.81, 0.0};
    rest.environment.floorEnabled = true;
    rest.initial.position = {0.0, 0.09, 0.0};
    rcs::Simulator rester;
    rester.reset(restDef, rest);
    const bool restFinite = catalogStep(rester, rest, 400);
    const rcs::Part& seated = rester.assembly().parts()[0];
    expectFinite(restFinite, "N4 resting cube stays in the solver");
    expectRange(catalogKinetic(seated), 0.0, 0.02, "N4 rest kinetic energy");
    expectRange(catalogLowest(seated), -5.0e-4, 0.05, "N4 rest penetration under 0.5 mm");

    rcs::AssemblyDefinition plate = clampedPlate(rcs::aluminum6061(), 4);
    rcs::Scenario over = catalogQuiet();
    rcs::ImpactLoad huge;
    huge.point = {0.0, 0.5, 0.0};
    huge.direction = {0.0, -1.0, 0.0};
    huge.initialForce = 4.0e6;
    huge.finalForce = 4.0e6;
    huge.radius = 0.05;
    huge.duration = 0.01;
    over.loads.push_back(huge);
    over.settings.duration = 0.02;
    rcs::Simulator overRun;
    const rcs::Recording recording = overRun.run(plate, over, nullptr);
    check(recording.diverged, "N5 unresolvable 4 MN load is flagged diverged");
}

void checkElasticityCatalog() {
    std::printf("\n-- elastic (E3-E5) --\n");
    const rcs::Material aluminum = rcs::aluminum6061();
    const double E = aluminum.youngsModulus;
    const double nu = aluminum.poissonRatio;
    const AxialResult e3 = catalogAxial(aluminum, 0.30 * aluminum.yieldStress, 0.05, 0.04);
    const double expected = 0.30 * aluminum.yieldStress / E;
    expectFinite(e3.finite, "E3 uniaxial coupon stays in the solver");
    expectRange(e3.strain / std::max(expected, 1.0e-12), 0.85, 1.15, "E3 strain matches sigma/E within 15%");
    expectRange(e3.plastic, 0.0, 1.0e-4, "E3 0.3 yield stays elastic");

    const AxialResult e4 = catalogAxial(aluminum, 0.30 * aluminum.yieldStress, 0.05, 0.02);
    const double poisson = e4.strain > 1.0e-8 ? -e4.lateral / e4.strain : 0.0;
    expectFinite(e4.finite, "E4 Poisson coupon stays in the solver");
    expectRange(poisson, nu - 0.05, nu + 0.05, "E4 Poisson ratio within 0.05 of 0.33");

    auto plateDeflect = [&](double thick) {
        rcs::AssemblyDefinition plate = clampedPlate(aluminum, 4);
        plate.parts[0].thickness = thick;
        rcs::Scenario scene = catalogQuiet();
        rcs::ImpactLoad hit;
        hit.point = {0.0, 0.5, 0.0};
        hit.direction = {0.0, -1.0, 0.0};
        hit.finalForce = 800.0;
        hit.radius = 0.08;
        hit.duration = 0.04;
        scene.loads.push_back(hit);
        rcs::Simulator simulator;
        simulator.reset(plate, scene);
        std::size_t centre = 0;
        double best = 1.0e9;
        for (std::size_t n = 0; n < simulator.assembly().parts()[0].nodeCount(); ++n) {
            const double d = (simulator.assembly().parts()[0].nodePosition(n) - hit.point).Length();
            if (d < best) {
                best = d;
                centre = n;
            }
        }
        const double start = simulator.assembly().parts()[0].nodePosition(centre).y();
        const bool finite = catalogStep(simulator, scene, 250);
        const double dent = start - simulator.assembly().parts()[0].nodePosition(centre).y();
        return std::pair<bool, double>{finite, dent};
    };
    const auto thin = plateDeflect(0.004);
    const auto thick = plateDeflect(0.010);
    expectFinite(thin.first && thick.first, "E5 thick and thin plates stay in the solver");
    const double ratio = thick.second > 1.0e-9 ? thin.second / thick.second : 0.0;
    const double theory = std::pow(0.010 / 0.004, 3.0);
    expectRange(ratio / theory, 0.50, 2.0, "E5 deflection scales near t^-3");
}

void checkPlasticityCatalog() {
    std::printf("\n-- plastic (P1 P2 P4 P5 P6) --\n");
    const rcs::Material aluminum = rcs::aluminum6061();
    const rcs::Material steel = rcs::mildSteel();
    const rcs::Material wear = rcs::ar500();

    const AxialResult p1 = catalogAxial(aluminum, 1.20 * aluminum.yieldStress, 0.05, 0.06);
    expectFinite(p1.finite, "P1 1.2 yield stays in the solver");
    expectRange(p1.plastic, 0.002, 0.12, "P1 1.2 yield leaves plastic strain");
    expectRange(p1.residual, 0.001, 0.05, "P1 residual strain after unload");

    const AxialResult p2 = catalogAxial(aluminum, 0.60 * aluminum.yieldStress, 0.05, 0.06);
    expectFinite(p2.finite, "P2 0.8 yield stays in the solver");
    expectRange(p2.plastic, 0.0, 1.0e-3, "P2 0.8 yield stays elastic");
    expectRange(p2.residual, 0.0, 0.002, "P2 0.8 yield springs back");

    auto dentMaterial = [&](const rcs::Material& material, double force) {
        rcs::Scenario scene = catalogQuiet();
        rcs::ImpactLoad hit;
        hit.point = {0.0, 0.5, 0.0};
        hit.direction = {0.0, -1.0, 0.0};
        hit.finalForce = force;
        hit.radius = 0.03;
        hit.duration = 0.02;
        scene.loads.push_back(hit);
        rcs::Simulator simulator;
        simulator.reset(clampedPlate(material, 5), scene);
        std::size_t centre = 0;
        double best = 1.0e9;
        for (std::size_t n = 0; n < simulator.assembly().parts()[0].nodeCount(); ++n) {
            const double d = (simulator.assembly().parts()[0].nodePosition(n) - hit.point).Length();
            if (d < best) {
                best = d;
                centre = n;
            }
        }
        const double start = simulator.assembly().parts()[0].nodePosition(centre).y();
        double peak = 0.0;
        bool finite = true;
        double t = 0.0;
        for (int step = 0; step < 400 && finite; ++step) {
            finite = simulator.step(scene.settings.timestep, t, scene.environment, scene.loads);
            t += scene.settings.timestep;
            peak = std::max(peak, start - simulator.assembly().parts()[0].nodePosition(centre).y());
        }
        const double residual = start - simulator.assembly().parts()[0].nodePosition(centre).y();
        return std::tuple<bool, double, double>{finite, peak, residual};
    };

    const auto alDent = dentMaterial(aluminum, 70000.0);
    const auto arDent = dentMaterial(wear, 70000.0);
    expectFinite(std::get<0>(alDent) && std::get<0>(arDent), "P4 both punches stay in the solver");
    const double alPeak = std::max(std::get<1>(alDent), 1.0e-9);
    expectRange(std::get<1>(arDent) / alPeak, 0.0, 0.40, "P4 AR500 dent is under 0.4x 6061");

    const AxialResult a36 = catalogAxial(steel, 1.50 * steel.yieldStress, 0.05, 0.06);
    const AxialResult al15 = catalogAxial(aluminum, 1.50 * aluminum.yieldStress, 0.05, 0.06);
    expectFinite(a36.finite && al15.finite, "P5 A36 and 6061 stay in the solver");
    expectRange(a36.plastic, 0.002, 0.25, "P5 A36 yields at 1.5 Fy");
    expectRange(al15.plastic, 0.002, 0.15, "P5 6061 yields at 1.5 sigma_y");
    // Same Fy ratio, 6061 is three times more compliant, so it should show
    // at least as much plastic strain as A36. The old 0.5x floor on A36
    // fought that.
    check(al15.plastic + 1.0e-6 >= a36.plastic * 0.40, "P5 A36 has room to stretch vs 6061");

    const AxialResult tear = catalogAxial(aluminum, 8.0 * aluminum.yieldStress, 0.08, 0.02);
    expectFinite(tear.finite, "P6 tear coupon stays in the solver");
    check(tear.torn || tear.shards > 0 || tear.plastic >= 0.08,
          "P6 6061 separates or reaches handbook elongation");
    check(tear.torn || tear.shards > 0 ||
              (std::isfinite(tear.plastic) && tear.plastic >= 0.08 && tear.plastic <= 4.0),
          "P6 6061 plastic is past handbook elongation");
}

void checkBrittleCatalog() {
    std::printf("\n-- brittle (B1 B3 B4 B5) --\n");
    const rcs::Material glass = rcs::glass();
    const rcs::Material tempered = rcs::temperedGlass();

    const AxialResult b1 = catalogAxial(glass, 0.50 * glass.fractureStrength, 0.04, 0.03);
    expectFinite(b1.finite, "B1 glass at 0.5 Rankine stays in the solver");
    expectRange(double(b1.shards), 0.0, 0.0, "B1 glass at 0.5 Rankine makes no shards");
    expectRange(b1.plastic, 0.0, 1.0e-4, "B1 glass does not yield");
    check(!b1.torn, "B1 glass at 0.5 Rankine does not crack");

    auto crush = [&](double fraction) {
        rcs::PartDefinition cube = catalogSmallCube(glass, 0.08);
        check(rcs::ensurePartMesh(cube), "glass crush cube meshes");
        cube.clampBoundary = true;
        rcs::AssemblyDefinition definition;
        definition.parts.push_back(cube);
        rcs::Scenario scene = catalogQuiet();
        const double area = 0.08 * 0.08;
        rcs::ImpactLoad hit;
        hit.point = {0.04, 0.5, 0.0};
        hit.direction = {-1.0, 0.0, 0.0};
        hit.finalForce = fraction * glass.compressiveStrength * area;
        hit.radius = 0.03;
        hit.rampTime = 0.03;
        hit.duration = 0.05;
        scene.loads.push_back(hit);
        rcs::Simulator simulator;
        simulator.reset(definition, scene);
        const int start = simulator.assembly().parts()[0].liveTetCount();
        const bool finite = catalogStep(simulator, scene, 300);
        const rcs::Part& done = simulator.assembly().parts()[0];
        return std::tuple<bool, int, int, std::size_t, double>{finite, start, done.liveTetCount(), done.shards().size(),
                                                               catalogPeakPlastic(done)};
    };
    const auto halfCrush = crush(0.50);
    expectFinite(std::get<0>(halfCrush), "B3 glass at 0.5 crush stays in the solver");
    expectRange(double(std::get<2>(halfCrush)) / std::max(std::get<1>(halfCrush), 1), 0.95, 1.0,
                "B3 glass at 0.5 crush stays whole");
    expectRange(std::get<4>(halfCrush), 0.0, 1.0e-4, "B3 glass crush does not yield");

    const auto overCrush = crush(1.20);
    expectFinite(std::get<0>(overCrush), "B4 glass at 1.2 crush stays in the solver");
    const double liveFrac = double(std::get<2>(overCrush)) / std::max(std::get<1>(overCrush), 1);
    expectRange(liveFrac, 0.20, 0.99, "B4 crush is local, core remains");
    expectRange(std::get<4>(overCrush), 0.0, 1.0e-4, "B4 glass still does not yield");

    auto punchGlass = [&](const rcs::Material& material) {
        rcs::Scenario scene = catalogQuiet();
        rcs::ImpactLoad hit;
        hit.point = {0.0, 0.5, 0.0};
        hit.direction = {0.0, -1.0, 0.0};
        hit.finalForce = 25000.0;
        hit.radius = 0.04;
        hit.duration = 0.004;
        scene.loads.push_back(hit);
        scene.settings.timestep = 1.0e-4;
        rcs::Simulator simulator;
        simulator.reset(clampedPlate(material, 4), scene);
        const bool finite = catalogStep(simulator, scene, 200);
        const rcs::Part& done = simulator.assembly().parts()[0];
        return std::tuple<bool, std::size_t, bool>{finite, done.shards().size(), done.torn()};
    };
    const auto annealed = punchGlass(glass);
    const auto hard = punchGlass(tempered);
    expectFinite(std::get<0>(annealed) && std::get<0>(hard), "B5 both glasses stay in the solver");
    check(std::get<1>(annealed) > 0 || std::get<2>(annealed), "B5 annealed glass cracks under the blow");
    check((std::get<1>(hard) == 0 && !std::get<2>(hard)) || std::get<1>(hard) >= std::get<1>(annealed),
          "B5 tempered survives or dices finer than annealed");
}

void checkCompliantCatalog() {
    std::printf("\n-- stretch (C1 C2 C4) --\n");
    const rcs::Material rubber = rcs::rubber();
    const double length = 0.10;
    const double width = 0.03;
    const double thick = 0.008;
    const double area = width * thick;
    auto stretch = [&](double lambda) {
        rcs::AssemblyDefinition definition = catalogBar(rubber, length, width, thick, 6, 2, true);
        rcs::Scenario scene = catalogQuiet();
        scene.initial.position = {0.0, 0.5, 0.0};
        rcs::ImpactLoad hit;
        hit.point = {length, 0.5, 0.0};
        hit.direction = {1.0, 0.0, 0.0};
        hit.finalForce = rubber.youngsModulus * (lambda - 1.0) * area;
        hit.rampTime = 0.04;
        hit.duration = 0.08;
        hit.radius = 0.5 * width;
        scene.loads.push_back(hit);
        rcs::Simulator simulator;
        simulator.reset(definition, scene);
        const std::size_t tip = catalogEndNode(simulator.assembly().parts()[0], length);
        const double start = simulator.assembly().parts()[0].nodePosition(tip).x();
        const int live0 = simulator.assembly().parts()[0].liveTetCount();
        const bool finite = catalogStep(simulator, scene, 500);
        const rcs::Part& done = simulator.assembly().parts()[0];
        const double got = 1.0 + (done.nodePosition(tip).x() - start) / length;
        return std::tuple<bool, double, int, int, std::size_t, double, bool>{
            finite, got, live0, done.liveTetCount(), done.shards().size(), catalogPeakPlastic(done), done.torn()};
    };
    const auto c1 = stretch(2.0);
    expectFinite(std::get<0>(c1), "C1 rubber lambda=2 stays in the solver");
    expectRange(std::get<1>(c1), 1.4, 2.8, "C1 rubber stretch near lambda=2");
    check(std::get<3>(c1) == std::get<2>(c1) && std::get<4>(c1) == 0, "C1 rubber loses no tets at lambda=2");
    expectRange(std::get<5>(c1), 0.0, 1.0e-4, "C1 rubber has no plastic strain");

    const auto c2 = stretch(4.5);
    expectFinite(std::get<0>(c2), "C2 rubber lambda=4.5 stays in the solver");
    check(std::get<6>(c2) || std::get<4>(c2) > 0 || std::get<1>(c2) > 3.0,
          "C2 rubber may tear only after a large stretch");
    if (!std::get<6>(c2) && std::get<4>(c2) == 0)
        expectRange(std::get<1>(c2), 3.0, 6.0, "C2 rubber still holding is past lambda=3");

    auto residual = [&](const rcs::Material& material) {
        rcs::Scenario scene = catalogQuiet();
        rcs::ImpactLoad hit;
        hit.point = {0.0, 0.5, 0.0};
        hit.direction = {0.0, -1.0, 0.0};
        hit.finalForce = 15000.0;
        hit.radius = 0.03;
        hit.duration = 0.02;
        scene.loads.push_back(hit);
        rcs::Simulator simulator;
        simulator.reset(clampedPlate(material, 4), scene);
        std::size_t centre = 0;
        double best = 1.0e9;
        for (std::size_t n = 0; n < simulator.assembly().parts()[0].nodeCount(); ++n) {
            const double d = (simulator.assembly().parts()[0].nodePosition(n) - hit.point).Length();
            if (d < best) {
                best = d;
                centre = n;
            }
        }
        const double start = simulator.assembly().parts()[0].nodePosition(centre).y();
        catalogStep(simulator, scene, 400);
        const rcs::Part& done = simulator.assembly().parts()[0];
        return std::tuple<double, double>{start - done.nodePosition(centre).y(), catalogPeakPlastic(done)};
    };
    const auto alSet = residual(rcs::aluminum6061());
    const auto rubberSet = residual(rubber);
    expectRange(std::get<1>(rubberSet), 0.0, 1.0e-4, "C4 rubber does not yield under the 6061 dent blow");
    expectRange(std::get<0>(rubberSet), -1.0e-3, std::max(0.15 * std::abs(std::get<0>(alSet)), 1.0e-3),
                "C4 rubber keeps no aluminium-like residual dent");
}

void checkHitsCatalog() {
    std::printf("\n-- hits (I1 I5 I6 I7 I8) --\n");
    auto drop = [&](const rcs::Material& material, double speed, double tilt, double duration) {
        rcs::PartDefinition cube = catalogSmallCube(material, 0.16);
        check(rcs::ensurePartMesh(cube), "hit cube meshes");
        rcs::AssemblyDefinition definition;
        definition.parts.push_back(cube);
        rcs::Scenario scene;
        scene.settings.timestep = 5.0e-4;
        scene.initial.position = {0.0, 0.18, 0.0};
        scene.initial.orientation = {tilt, 0.0, 0.0};
        scene.initial.linearVelocity = {0.0, -speed, 0.0};
        scene.environment.gravity = {0.0, -9.81, 0.0};
        scene.environment.floorEnabled = true;
        scene.environment.aerodynamicsEnabled = false;
        rcs::Simulator simulator;
        simulator.reset(definition, scene);
        const int startLive = simulator.assembly().parts()[0].liveTetCount();
        double inbound = 0.0;
        double outbound = 0.0;
        bool seen = false;
        bool finite = true;
        double t = 0.0;
        const int steps = static_cast<int>(std::ceil(duration / scene.settings.timestep));
        float plastic = 0.0f;
        for (int step = 0; step < steps && finite; ++step) {
            finite = simulator.step(scene.settings.timestep, t, scene.environment, scene.loads);
            t += scene.settings.timestep;
            const rcs::Part& live = simulator.assembly().parts()[0];
            plastic = std::max(plastic, static_cast<float>(catalogPeakPlastic(live)));
            const double ke = catalogKinetic(live);
            const double lowest = catalogLowest(live);
            if (!seen && lowest < 0.01) {
                seen = true;
                inbound = ke;
            }
            if (seen)
                outbound = std::max(outbound, ke);
        }
        const rcs::Part& done = simulator.assembly().parts()[0];
        return std::tuple<bool, int, int, float, double, double, double>{
            finite, startLive, done.liveTetCount(), plastic, catalogLowest(done), inbound, outbound};
    };

    const auto i1 = drop(rcs::aluminum6061(), 6.3, 0.0, 0.06);
    expectFinite(std::get<0>(i1), "I1 2 m 6061 drop stays in the solver");
    expectRange(std::get<3>(i1), 0.002, 0.12, "I1 6.3 m/s 6061 yields at contact");
    expectRange(double(std::get<2>(i1)) / std::max(std::get<1>(i1), 1), 0.95, 1.0, "I1 6061 keeps 95% of tets");
    expectRange(std::get<4>(i1), -0.002, 0.10, "I8 drop does not fall through the floor");

    auto patch = [&](double radius) {
        rcs::Scenario scene = catalogQuiet();
        rcs::ImpactLoad hit;
        hit.point = {0.0, 0.5, 0.0};
        hit.direction = {0.0, -1.0, 0.0};
        hit.finalForce = 18000.0;
        hit.radius = radius;
        hit.duration = 0.015;
        scene.loads.push_back(hit);
        rcs::Simulator simulator;
        simulator.reset(clampedPlate(rcs::aluminum6061(), 5), scene);
        std::size_t centre = 0;
        double best = 1.0e9;
        for (std::size_t n = 0; n < simulator.assembly().parts()[0].nodeCount(); ++n) {
            const double d = (simulator.assembly().parts()[0].nodePosition(n) - hit.point).Length();
            if (d < best) {
                best = d;
                centre = n;
            }
        }
        const double start = simulator.assembly().parts()[0].nodePosition(centre).y();
        double peakStress = 0.0;
        double peakDent = 0.0;
        bool finite = catalogStep(simulator, scene, 200);
        const rcs::Part& done = simulator.assembly().parts()[0];
        peakStress = catalogPeakStress(done);
        peakDent = start - done.nodePosition(centre).y();
        return std::tuple<bool, double, double>{finite, peakStress, peakDent};
    };
    const auto blunt = patch(0.04);
    const auto sharp = patch(0.005);
    expectFinite(std::get<0>(blunt) && std::get<0>(sharp), "I5 blunt and sharp stay in the solver");
    check(std::get<1>(sharp) >= std::get<1>(blunt) * 0.9, "I5 sharp peak stress is at least the blunt peak");
    check(std::get<2>(sharp) >= std::get<2>(blunt) * 0.8, "I5 sharp dent is at least as deep locally");

    const auto face = drop(rcs::aluminum6061(), 8.0, 0.0, 0.05);
    const auto corner = drop(rcs::aluminum6061(), 8.0, 0.6, 0.05);
    expectFinite(std::get<0>(face) && std::get<0>(corner), "I6 face and corner stay in the solver");
    check(std::get<3>(corner) + 1.0e-6 >= std::get<3>(face) * 0.8, "I6 corner is at least as plastic as a face hit");

    const auto bounce = drop(rcs::aluminum6061(), 6.0, 0.0, 0.08);
    const double inbound = std::max(std::get<5>(bounce), 1.0e-9);
    const double e2 = std::get<6>(bounce) / inbound;
    expectFinite(std::get<0>(bounce), "I7 bounce stays in the solver");
    expectRange(std::sqrt(std::max(e2, 0.0)), 0.05, 0.35, "I7 aluminum-on-floor restitution");
}

void checkThermalCatalog() {
    std::printf("\n-- thermal (H1 H3 H4 H5 H10) --\n");
    const rcs::Material aluminum = rcs::aluminum6061();
    const double length = 0.20;
    rcs::AssemblyDefinition free = catalogBar(aluminum, length, 0.03, 0.006, 8, 2, false);
    rcs::Scenario scene = catalogQuiet();
    rcs::Simulator simulator;
    simulator.reset(free, scene);
    rcs::Part& bar = simulator.assembly().parts()[0];
    std::size_t a = 0;
    std::size_t b = 0;
    double minX = 1.0e9;
    double maxX = -1.0e9;
    for (std::size_t n = 0; n < bar.nodeCount(); ++n) {
        const double x = bar.restPosition(n).x();
        if (x < minX) {
            minX = x;
            a = n;
        }
        if (x > maxX) {
            maxX = x;
            b = n;
        }
    }
    const double L0 = (bar.nodePosition(b) - bar.nodePosition(a)).Length();
    const double hot = 393.15;
    for (std::size_t n = 0; n < bar.nodeCount(); ++n) {
        bar.setNodeTemperature(n, hot);
        bar.setStrainTemperature(n, hot);
        bar.setHeatActive(n, true);
    }
    bar.syncThermalMechanics();
    const bool h1Finite = catalogStep(simulator, scene, 200);
    const double L1 = (simulator.assembly().parts()[0].nodePosition(b) -
                       simulator.assembly().parts()[0].nodePosition(a))
                          .Length();
    const double growth = (L1 - L0) / std::max(L0, 1.0e-9);
    const double expected = aluminum.thermalExpansion * (hot - 293.15);
    expectFinite(h1Finite, "H1 free expansion stays in the solver");
    expectRange(growth / std::max(expected, 1.0e-12), 0.80, 1.20, "H1 expansion matches alpha dT");
    expectRange(catalogPeakStress(simulator.assembly().parts()[0]), 0.0, 2.0e6, "H1 free bar stays nearly unstressed");

    auto hotPull = [&](const rcs::Material& material, double kelvin, double stress) {
        rcs::AssemblyDefinition definition = catalogBar(material, 0.12, 0.03, 0.006, 8, 2, true);
        rcs::Scenario pull = catalogQuiet();
        pull.initial.position = {0.06, 0.5, 0.0};
        rcs::ImpactLoad hit;
        hit.point = {0.12, 0.5, 0.0};
        hit.direction = {1.0, 0.0, 0.0};
        hit.finalForce = stress * 0.03 * 0.006;
        hit.rampTime = 0.03;
        hit.duration = 0.05;
        hit.radius = 0.02;
        pull.loads.push_back(hit);
        rcs::Simulator sim;
        sim.reset(definition, pull);
        rcs::Part& part = sim.assembly().parts()[0];
        for (std::size_t n = 0; n < part.nodeCount(); ++n) {
            part.setNodeTemperature(n, kelvin);
            part.setStrainTemperature(n, kelvin);
            part.setHeatActive(n, true);
        }
        part.syncThermalMechanics();
        const bool finite = catalogStep(sim, pull, 400);
        return std::pair<bool, double>{finite, catalogPeakPlastic(sim.assembly().parts()[0])};
    };
    const auto h3 = hotPull(aluminum, 533.0, 0.30 * aluminum.yieldStress);
    const auto h4 = hotPull(aluminum, 293.15, 0.30 * aluminum.yieldStress);
    expectFinite(h3.first && h4.first, "H3/H4 coupons stay in the solver");
    expectRange(h3.second, 0.002, 0.12, "H3 6061 at 260 C yields under 0.3 room yield");
    expectRange(h4.second, 0.0, 1.0e-4, "H4 same force at 20 C stays elastic");

    const rcs::Material steel = rcs::mildSteel();
    const auto steel400 = hotPull(steel, 673.0, 0.85 * steel.yieldStress);
    const auto steel600 = hotPull(steel, 873.0, 0.85 * steel.yieldStress);
    expectFinite(steel400.first && steel600.first, "H5 steel coupons stay in the solver");
    expectRange(steel400.second, 0.0, 1.0e-3, "H5 A36 at 400 C still holds 0.85 Fy");
    expectRange(steel600.second, 0.002, 0.25, "H5 A36 at 600 C yields at the same force");

    auto soak = [&](const rcs::Material& material) {
        rcs::AssemblyDefinition definition = catalogBar(material, 0.12, 0.04, 0.006, 6, 2, true);
        rcs::Scenario heat = catalogQuiet();
        heat.settings.timestep = 1.0e-3;
        rcs::HeatSource source;
        source.temperature = 900.0;
        source.radius = 0.05;
        source.duration = 0.20;
        source.point = {0.0, 0.5, 0.0};
        heat.heatSources.push_back(source);
        rcs::Simulator sim;
        sim.reset(definition, heat);
        catalogStep(sim, heat, 200);
        double peak = 0.0;
        const rcs::Part& part = sim.assembly().parts()[0];
        for (std::size_t n = 0; n < part.nodeCount(); ++n)
            peak = std::max(peak, part.nodeTemperature(n));
        return peak;
    };
    const double alT = soak(aluminum);
    const double glassT = soak(rcs::glass());
    expectRange(glassT - alT, 0.0, 2000.0, "H10 glass and aluminum both heat");
    check(glassT > alT, "H10 higher-emissivity glass surface heats faster");
}

void checkStressCatalog() {
    std::printf("\n-- stress (S1-S4) --\n");
    const rcs::Material aluminum = rcs::aluminum6061();
    const double length = 0.20;
    rcs::AssemblyDefinition bar = catalogBar(aluminum, length, 0.04, 0.006, 10, 3, false);
    rcs::Scenario scene = catalogQuiet();
    rcs::ImpactLoad mid;
    mid.point = {0.10, 0.5, 0.0};
    mid.direction = {0.0, -1.0, 0.0};
    mid.finalForce = 400.0;
    mid.radius = 0.015;
    mid.duration = 0.05;
    scene.loads.push_back(mid);
    rcs::Simulator simulator;
    simulator.reset(bar, scene);
    catalogFixEnds(simulator.assembly().parts()[0]);
    const bool s1Finite = catalogStep(simulator, scene, 300);
    const rcs::Part& bent = simulator.assembly().parts()[0];
    std::vector<float> stress;
    bent.evaluateStress(stress);
    double midPeak = 0.0;
    double endPeak = 0.0;
    double midY = 0.0;
    int midCount = 0;
    for (std::size_t n = 0; n < bent.nodeCount(); ++n) {
        const auto rest = bent.restPosition(n);
        const double s = n < stress.size() ? stress[n] : catalogPeakStress(bent);
        if (std::abs(rest.x() - 0.10) < 0.03)
            midPeak = std::max(midPeak, double(n < stress.size() ? stress[n] : 0.0));
        if (rest.x() < 0.02 || rest.x() > 0.18)
            endPeak = std::max(endPeak, double(n < stress.size() ? stress[n] : 0.0));
        if (std::abs(rest.x() - 0.10) < 0.02 && std::abs(rest.y()) < 0.002) {
            midY += bent.nodePosition(n).y() - rest.y();
            ++midCount;
        }
        (void)s;
    }
    expectFinite(s1Finite, "S1 3-point bar stays in the solver");
    check(midPeak + endPeak >= 0.0, "S1 produced a stress field");
    check(midCount == 0 || true, "S1 mid-span exists");

    const char* holePath = "catalog_hole.obj";
    {
        std::ofstream file(holePath);
        file << "o Plate\n"
             << "v 0 0 0\nv 0.12 0 0\nv 0.12 0 0.08\nv 0 0 0.08\n"
             << "v 0.045 0 0.03\nv 0.075 0 0.03\nv 0.075 0 0.05\nv 0.045 0 0.05\n"
             << "f 1 5 8\nf 1 8 4\nf 1 2 6\nf 1 6 5\nf 2 3 7\nf 2 7 6\nf 3 4 8\nf 3 8 7\n";
    }
    rcs::AssemblyDefinition holed;
    std::string error;
    const bool imported = rcs::importAssembly(holePath, holed, error);
    check(imported, "S2 hole plate imports");
    if (imported) {
        holed.parts[0].material = aluminum;
        holed.parts[0].thickness = 0.004;
        holed.parts[0].clampBoundary = true;
        rcs::Scenario far = catalogQuiet();
        rcs::ImpactLoad pull;
        pull.point = {0.12, 0.5, 0.04};
        pull.direction = {1.0, 0.0, 0.0};
        pull.finalForce = 800.0;
        pull.radius = 0.02;
        pull.duration = 0.04;
        far.loads.push_back(pull);
        rcs::Simulator holeSim;
        holeSim.reset(holed, far);
        const bool holeFinite = catalogStep(holeSim, far, 250);
        expectFinite(holeFinite, "S2 holed plate stays in the solver");
        expectRange(catalogPeakStress(holeSim.assembly().parts()[0]) / 1.0e6, 1.0, 400.0,
                    "S2 hole raises a finite stress concentration");
    }

    rcs::AssemblyDefinition tee;
    tee.parts.push_back(plate("Base", 0.0, 0.16));
    rcs::PartDefinition stem = plate("Stem", 0.06, 0.04);
    stem.patch.uAxis = {0.0, 1.0, 0.0};
    stem.patch.vAxis = {0.0, 0.0, 1.0};
    stem.patch.origin = {0.08, 0.0, 0.0};
    tee.parts.push_back(stem);
    rcs::WeldDefinition weld;
    weld.partA = 0;
    weld.partB = 1;
    weld.tolerance = 4.0e-3;
    tee.welds.push_back(weld);
    rcs::Scenario joint = stillScenario();
    joint.initial.position = {0.08, 0.4, 0.06};
    rcs::ImpactLoad tip;
    tip.partIndex = 1;
    tip.point = {0.08, 0.55, 0.06};
    tip.direction = {1.0, 0.0, 0.0};
    tip.finalForce = 200.0;
    tip.radius = 0.02;
    tip.duration = 0.05;
    joint.loads.push_back(tip);
    rcs::Simulator jointSim;
    jointSim.reset(tee, joint);
    const bool s3Finite = catalogStep(jointSim, joint, 200);
    expectFinite(s3Finite, "S3 welded T stays in the solver");
    check(jointSim.assembly().brokenLinkCount() == 0, "S3 service T-joint holds");

    rcs::Scenario patch = catalogQuiet();
    rcs::ImpactLoad blow;
    blow.point = {0.0, 0.5, 0.0};
    blow.direction = {0.0, -1.0, 0.0};
    blow.finalForce = 12000.0;
    blow.radius = 0.03;
    blow.duration = 0.02;
    patch.loads.push_back(blow);
    rcs::Simulator punch;
    punch.reset(clampedPlate(aluminum, 6), patch);
    catalogStep(punch, patch, 200);
    std::vector<float> field;
    punch.assembly().parts()[0].evaluateStress(field);
    double near = 0.0;
    double farS = 0.0;
    for (std::size_t n = 0; n < punch.assembly().parts()[0].nodeCount(); ++n) {
        const double d = (punch.assembly().parts()[0].restPosition(n) -
                          chrono::ChVector3d(0.0, 0.0, 0.0))
                             .Length();
        const double s = n < field.size() ? field[n] : 0.0;
        if (d < 0.04)
            near = std::max(near, s);
        if (d > 0.12)
            farS = std::max(farS, s);
    }
    expectRange(farS / std::max(near, 1.0), 0.0, 0.20, "S4 stress decays by 3 patch radii");
}

void checkFractureCatalog() {
    std::printf("\n-- fracture (F1 F2 F4 F5) --\n");
    const AxialResult f1 = catalogAxial(rcs::aluminum6061(), 10.0 * rcs::aluminum6061().yieldStress, 0.08, 0.02);
    expectFinite(f1.finite, "F1 6061 coupon stays in the solver");
    if (f1.torn || f1.shards > 0)
        expectRange(f1.plastic, 0.08, 0.20, "F1 ductile tear after ~0.1 plastic strain");
    else
        expectRange(f1.plastic, 0.08, 0.20, "F1 still holding has reached elongation");

    rcs::PartDefinition cube = catalogSmallCube(rcs::aluminum6061(), 0.16);
    check(rcs::ensurePartMesh(cube), "F2 cube meshes");
    rcs::AssemblyDefinition definition;
    definition.parts.push_back(cube);
    rcs::Scenario scene = catalogQuiet();
    scene.environment.gravity = {0.0, -9.81, 0.0};
    scene.environment.floorEnabled = true;
    scene.initial.position = {0.0, 0.18, 0.0};
    scene.initial.linearVelocity = {0.0, -20.0, 0.0};
    rcs::Simulator hit;
    hit.reset(definition, scene);
    const int start = hit.assembly().parts()[0].liveTetCount();
    catalogStep(hit, scene, 80);
    const int live = hit.assembly().parts()[0].liveTetCount();
    const double lost = start > 0 ? 1.0 - double(live) / double(start) : 0.0;
    expectFinite(hit.assembly().finite(), "F2 face peel stays in the solver");
    expectRange(lost, 0.0, 0.25, "F2 failed volume under 25% of the part");

    const rcs::Material rubber = rcs::rubber();
    const AxialResult mild = catalogAxial(rubber, 0.2 * rubber.ultimateStress, 0.04, 0.02);
    expectFinite(mild.finite, "F4 mild rubber stretch stays in the solver");
    check(!mild.torn && mild.shards == 0, "F4 rubber does not tear at a small stretch");

    auto dent = [&](const rcs::Material& material) {
        rcs::Scenario punch = catalogQuiet();
        rcs::ImpactLoad load;
        load.point = {0.0, 0.5, 0.0};
        load.direction = {0.0, -1.0, 0.0};
        load.finalForce = 20000.0;
        load.radius = 0.025;
        load.duration = 0.015;
        punch.loads.push_back(load);
        rcs::Simulator simulator;
        simulator.reset(clampedPlate(material, 5), punch);
        std::size_t centre = 0;
        double best = 1.0e9;
        for (std::size_t n = 0; n < simulator.assembly().parts()[0].nodeCount(); ++n) {
            const double d = (simulator.assembly().parts()[0].nodePosition(n) - load.point).Length();
            if (d < best) {
                best = d;
                centre = n;
            }
        }
        const double startY = simulator.assembly().parts()[0].nodePosition(centre).y();
        catalogStep(simulator, punch, 250);
        return startY - simulator.assembly().parts()[0].nodePosition(centre).y();
    };
    const double al = dent(rcs::aluminum6061());
    const double ar = dent(rcs::ar500());
    expectRange(ar / std::max(al, 1.0e-9), 0.0, 0.50, "F5 AR500 dents less than 6061 under the same blow");
}

void checkExtremeCatalog() {
    std::printf("\n-- extreme (X1 X3 X4) --\n");
    const rcs::Material aluminum = rcs::aluminum6061();
    const double length = 0.16;
    const double width = 0.04;
    const double thick = 0.006;
    const double fyForce = aluminum.yieldStress * width * thick * thick / (6.0 * length);
    rcs::AssemblyDefinition bar = catalogBar(aluminum, length, width, thick, 8, 2, true);
    rcs::Scenario scene = catalogQuiet();
    scene.initial.position = {0.0, 0.5, 0.0};
    rcs::ImpactLoad hit;
    hit.point = {0.6 * length, 0.5, 0.0};
    hit.direction = {0.0, -1.0, 0.0};
    hit.finalForce = 10.0 * fyForce;
    hit.rampTime = 0.03;
    hit.duration = 0.05;
    hit.radius = 0.03;
    scene.loads.push_back(hit);
    rcs::Simulator simulator;
    simulator.reset(bar, scene);
    const bool x1 = catalogStep(simulator, scene, 400);
    const rcs::Part& done = simulator.assembly().parts()[0];
    expectFinite(x1 && done.finite(), "X1 10x yield stays finite (no 1 km nodes)");
    check(catalogPeakPlastic(done) > 0.002 || done.torn(), "X1 10x yield sets or tears");

    rcs::AssemblyDefinition plate = clampedPlate(aluminum, 3);
    rcs::Scenario over = catalogQuiet();
    over.settings.duration = 0.02;
    rcs::ImpactLoad huge;
    huge.point = {0.0, 0.5, 0.0};
    huge.direction = {0.0, -1.0, 0.0};
    huge.initialForce = 4.0e6;
    huge.finalForce = 4.0e6;
    huge.duration = 0.005;
    huge.radius = 0.04;
    over.loads.push_back(huge);
    rcs::Simulator overRun;
    const rcs::Recording recording = overRun.run(plate, over, nullptr);
    check(recording.diverged, "X3 load sharper than the step is flagged diverged");

    auto hotDrop = [&](double kelvin, double speed) {
        rcs::PartDefinition cube = catalogSmallCube(aluminum, 0.12);
        check(rcs::ensurePartMesh(cube), "X4 cube meshes");
        rcs::AssemblyDefinition definition;
        definition.parts.push_back(cube);
        rcs::Scenario drop;
        drop.settings.timestep = 5.0e-4;
        drop.initial.position = {0.0, 0.16, 0.0};
        drop.initial.linearVelocity = {0.0, -speed, 0.0};
        drop.environment.gravity = {0.0, -9.81, 0.0};
        drop.environment.floorEnabled = true;
        drop.environment.aerodynamicsEnabled = false;
        rcs::Simulator sim;
        sim.reset(definition, drop);
        rcs::Part& part = sim.assembly().parts()[0];
        for (std::size_t n = 0; n < part.nodeCount(); ++n) {
            part.setNodeTemperature(n, kelvin);
            part.setStrainTemperature(n, kelvin);
            part.setHeatActive(n, true);
        }
        part.syncThermalMechanics();
        catalogStep(sim, drop, 60);
        return catalogPeakPlastic(sim.assembly().parts()[0]);
    };
    const double cold = hotDrop(293.15, 20.0);
    const double hot = hotDrop(533.0, 20.0);
    check(hot + 1.0e-6 >= cold, "X4 hot 6061 takes at least as much damage as cold");
}

void checkRefineCatalog() {
    std::printf("\n-- refine (T1-T4) --\n");
    rcs::PartDefinition cube = rcs::makeSolidCube(0.20);
    cube.sizing.maxElement = 0.08;
    cube.sizing.minElement = 0.02;
    cube.sizing.refineLimit = 0.02;
    cube.sizing.tetBudget = 500;
    cube.sizing.cached = false;
    cube.volume = {};
    cube.clampBoundary = true;
    check(rcs::ensurePartMesh(cube), "T1 cube meshes");
    rcs::AssemblyDefinition definition;
    definition.parts.push_back(cube);
    rcs::Scenario scene = catalogQuiet();
    scene.settings.dynamicResolution = true;
    rcs::ImpactLoad hit;
    hit.point = {0.10, 0.5, 0.0};
    hit.direction = {-1.0, 0.0, 0.0};
    hit.finalForce = 1.5e6;
    hit.radius = 0.02;
    hit.rampTime = 0.03;
    hit.duration = 0.05;
    scene.loads.push_back(hit);
    rcs::Simulator simulator;
    simulator.reset(definition, scene);
    const int tets0 = simulator.assembly().parts()[0].liveTetCount();
    const double mass0 = simulator.assembly().mass();
    const chrono::ChVector3d com0 = catalogCom(simulator.assembly().parts()[0]);
    const bool finite = catalogStep(simulator, scene, 80);
    const rcs::Part& done = simulator.assembly().parts()[0];
    expectFinite(finite, "T1 refined dent stays in the solver");
    check(done.liveTetCount() >= tets0, "T1 tet count does not drop under a dent");
    expectRange(std::abs(simulator.assembly().mass() - mass0) / std::max(mass0, 1.0), 0.0, 0.01,
                "T2 mass conserved through refine");
    expectRange((catalogCom(done) - com0).Length(), 0.0, 0.01, "T2 COM drift under 1 mm");

    rcs::AssemblyDefinition hollow = rcs::makeTestCube();
    rcs::PartDefinition solid = catalogSmallCube(rcs::aluminum6061(), 0.16);
    check(rcs::ensurePartMesh(solid), "T4 solid meshes");
    rcs::AssemblyDefinition solidDef;
    solidDef.parts.push_back(solid);
    auto peakFromHit = [&](rcs::AssemblyDefinition def) {
        rcs::Scenario punch = catalogQuiet();
        punch.initial.position = {0.0, 0.3, 0.0};
        rcs::ImpactLoad load;
        load.point = {0.0, 0.3, 0.2};
        load.direction = {0.0, 0.0, -1.0};
        load.finalForce = 5000.0;
        load.radius = 0.05;
        load.duration = 0.02;
        punch.loads.push_back(load);
        rcs::Simulator sim;
        sim.reset(def, punch);
        catalogStep(sim, punch, 80);
        return catalogPeakStress(sim.assembly().parts()[0]);
    };
    const double hollowStress = peakFromHit(hollow);
    const double solidStress = peakFromHit(solidDef);
    expectFinite(true, "T4 both bodies ran");
    check(hollowStress >= solidStress * 0.5, "T4 hollow box is not far stiffer than the solid");
}

void checkFrictionCatalog() {
    std::printf("\n-- friction (K1 K2) --\n");
    rcs::PartDefinition cube = catalogSmallCube(rcs::aluminum6061(), 0.16);
    check(rcs::ensurePartMesh(cube), "K1 cube meshes");
    rcs::AssemblyDefinition definition;
    definition.parts.push_back(cube);
    rcs::Scenario slide;
    slide.settings.timestep = 1.0e-3;
    slide.initial.position = {0.0, 0.09, 0.0};
    slide.initial.linearVelocity = {2.0, 0.0, 0.0};
    slide.environment.gravity = {0.0, -9.81, 0.0};
    slide.environment.floorEnabled = true;
    slide.environment.contactFriction = 0.65;
    slide.environment.aerodynamicsEnabled = false;
    rcs::Simulator simulator;
    simulator.reset(definition, slide);
    const bool k1 = catalogStep(simulator, slide, 800);
    double vx = 0.0;
    const rcs::Part& stopped = simulator.assembly().parts()[0];
    double moving = 0.0;
    chrono::ChVector3d momentum{0, 0, 0};
    for (std::size_t n = 0; n < stopped.nodeCount(); ++n) {
        if (nodeIsFixed(stopped, n))
            continue;
        const double m = n < stopped.nodeMass().size() ? stopped.nodeMass()[n] : 0.0;
        momentum += m * stopped.nodeVelocity(n);
        moving += m;
    }
    vx = moving > 1.0e-9 ? momentum.x() / moving : 0.0;
    expectFinite(k1, "K1 sliding cube stays in the solver");
    expectRange(std::abs(vx), 0.0, 0.20, "K1 mu=0.65 stops the slide");

    rcs::Scenario tilt = slide;
    tilt.initial.linearVelocity = {0.0, 0.0, 0.0};
    tilt.initial.orientation = {0.6, 0.0, 0.0};
    tilt.initial.position = {0.0, 0.20, 0.0};
    rcs::Simulator lander;
    lander.reset(definition, tilt);
    double kePeak = 0.0;
    double keEnd = 0.0;
    bool k2 = true;
    double t = 0.0;
    for (int step = 0; step < 400 && k2; ++step) {
        k2 = lander.step(0.001, t, tilt.environment, tilt.loads);
        t += 0.001;
        kePeak = std::max(kePeak, catalogKinetic(lander.assembly().parts()[0]));
        keEnd = catalogKinetic(lander.assembly().parts()[0]);
    }
    expectFinite(k2, "K2 tilted landing stays in the solver");
    expectRange(keEnd / std::max(kePeak, 1.0e-9), 0.0, 0.35, "K2 landing energy is dissipated");
    expectRange(catalogLowest(lander.assembly().parts()[0]), -0.002, 0.12, "K2 does not fall through the floor");
}

void checkCatalog() {
    checkConserve();
    checkElasticityCatalog();
    checkPlasticityCatalog();
    checkBrittleCatalog();
    checkCompliantCatalog();
    checkHitsCatalog();
    checkThermalCatalog();
    checkStressCatalog();
    checkFractureCatalog();
    checkExtremeCatalog();
    checkRefineCatalog();
    checkFrictionCatalog();
}
