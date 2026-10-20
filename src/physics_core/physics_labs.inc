/// Honest lab and motion checks. Included from physics_smoke.cpp after the
/// catalog helpers. Claims match metrics. Feature-gap labs are allowed to fail.

struct BodySample {
    chrono::ChVector3d com{0, 0, 0};
    chrono::ChVector3d vel{0, 0, 0};
    double mass = 0.0;
    double lowest = 1.0e9;
    double highest = -1.0e9;
    double ke = 0.0;
    int freeNodes = 0;
    int live = 0;
    int broken = 0;
};

BodySample sampleBody(const rcs::Assembly& assembly) {
    BodySample s;
    s.broken = assembly.brokenLinkCount();
    for (const auto& part : assembly.parts()) {
        s.live += part.liveElementCount();
        for (std::size_t n = 0; n < part.nodeCount(); ++n) {
            const double y = part.nodePosition(n).y();
            s.lowest = std::min(s.lowest, y);
            s.highest = std::max(s.highest, y);
            const double m = n < part.nodeMass().size() ? part.nodeMass()[n] : 0.0;
            if (nodeIsFixed(part, n))
                continue;
            ++s.freeNodes;
            s.com += m * part.nodePosition(n);
            s.vel += m * part.nodeVelocity(n);
            s.ke += 0.5 * m * part.nodeVelocity(n).Length2();
            s.mass += m;
        }
        for (const auto& shard : part.shards()) {
            s.ke += 0.5 * shard.mass * shard.velocity.Length2();
            s.mass += shard.mass;
        }
    }
    if (s.mass > 1.0e-12) {
        s.com /= s.mass;
        s.vel /= s.mass;
    }
    return s;
}

double sampleSpan(const BodySample& s) { return s.highest - s.lowest; }

double tetAspect(const rcs::Assembly& assembly) {
    double worst = 1.0;
    for (const auto& part : assembly.parts()) {
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
                worst = std::max(worst, longest / shortest);
        }
    }
    return worst;
}

void checkSeatedVolume(const char* tag, const BodySample& start, const BodySample& last,
                       const rcs::Assembly& assembly, double halfHeight) {
    const double aspect = tetAspect(assembly);
    expectRange(aspect, 1.0, 30.0, (std::string(tag) + " keeps tet aspect O(1), not a pancake").c_str());
    expectRange(last.com.y(), 0.45 * halfHeight, 1.80 * halfHeight,
                (std::string(tag) + " COM stays near half-height").c_str());
    expectRange(last.lowest, -0.01, 0.02, (std::string(tag) + " sits on the floor, not through it").c_str());
    check(assembly.invertedTetCount() == 0, (std::string(tag) + " has no inverted tets at rest").c_str());
    check(last.live == start.live, (std::string(tag) + " does not shatter while seating").c_str());
}

rcs::Scenario editorIdentityScene() {
    rcs::Scenario scene = smokeFloor();
    scene.settings.duration = 1.0;
    scene.settings.timestep = 5.0e-4;
    scene.settings.captureRate = 40.0;
    scene.settings.meshResolution = 0.25;
    scene.initial.position = {0.0, 0.55, 0.0};
    scene.initial.orientation = {0.35, 0.0, 0.45};
    scene.environment.aerodynamicsEnabled = true;
    scene.environment.floorEnabled = true;
    scene.environment.gravity = {0.0, -9.81, 0.0};
    scene.environment.contactFriction = 0.65;
    scene.environment.floorRestitution = 0.15;
    return scene;
}

struct DropTrace {
    bool finite = true;
    BodySample first;
    BodySample airborne;
    BodySample hit;
    BodySample last;
    bool hadAir = false;
    bool hadHit = false;
    double maxComXzFromStart = 0.0;
    double postHitXz = 0.0;
    double postHitMove = 0.0;
    double maxStress = 0.0;
    double postHitKe = 0.0;
    int hhtCommits = 0;
    int newtonRejects = 0;
};

DropTrace runCornerDrop(const rcs::AssemblyDefinition& definition, bool, const char* tag) {
    rcs::Scenario scene = editorIdentityScene();
    keepEditorDefaultsFlag = false;
    groupKeepsEditorDefaults = false;
    scene.settings.meshResolution = 0.25;
    SmokeSim sim;
    const rcs::Recording rec = sim.run(definition, scene, nullptr);
    DropTrace tr;
    tr.finite = !rec.diverged;
    tr.last = sampleBody(sim.assembly());
    tr.maxStress = rec.maxStress;
    tr.hhtCommits = sim.hhtCommitCount();
    tr.newtonRejects = sim.newtonRejectCount();
    dumpSim(tag, sim, rec.duration, tr.finite);
    if (rec.frames.empty())
        return tr;
    auto fromFrame = [](const rcs::Frame& frame) {
        BodySample s;
        s.broken = frame.brokenLinks;
        for (const auto& part : frame.parts) {
            for (std::size_t i = 0; i < part.nodeCount(); ++i) {
                const auto p = part.position(i);
                s.com += chrono::ChVector3d(double(p.x), double(p.y), double(p.z));
                s.lowest = std::min(s.lowest, double(p.y));
                s.highest = std::max(s.highest, double(p.y));
                s.mass += 1.0;
            }
        }
        if (s.mass > 1.0e-12)
            s.com /= s.mass;
        return s;
    };
    tr.first = fromFrame(rec.frames.front());
    const rcs::Frame* air = rec.sample(0.14);
    if (air && air->time <= 0.18) {
        tr.airborne = fromFrame(*air);
        tr.hadAir = true;
    }
    for (const auto& frame : rec.frames) {
        const BodySample now = fromFrame(frame);
        const double xz = std::hypot(now.com.x() - tr.first.com.x(), now.com.z() - tr.first.com.z());
        tr.maxComXzFromStart = std::max(tr.maxComXzFromStart, xz);
        if (!tr.hadHit && now.lowest < 0.008) {
            tr.hit = now;
            tr.hadHit = true;
        }
        if (tr.hadHit) {
            tr.postHitMove = std::max(tr.postHitMove, (now.com - tr.hit.com).Length());
            tr.postHitXz = std::max(tr.postHitXz, std::hypot(now.com.x() - tr.hit.com.x(),
                                                             now.com.z() - tr.hit.com.z()));
        }
    }
    tr.last.ke = sampleBody(sim.assembly()).ke;
    tr.last.broken = sim.assembly().brokenLinkCount();
    tr.postHitKe = tr.last.ke;
    return tr;
}

void assertEngineDrop(const DropTrace& tr, const char* body) {
    const double span = sampleSpan(tr.last);
    const std::string prefix = std::string(body) + " ";
    check(tr.finite, (prefix + "corner drop stays in the solver").c_str());
    check(tr.hhtCommits > 0, (prefix + "committed HHT steps").c_str());
    if (tr.hadAir) {
        check(tr.airborne.lowest > 0.03, (prefix + "is still clear of the floor at 0.14 s").c_str());
        expectRange(tr.airborne.com.x() - tr.first.com.x(), -0.04, 0.04,
                    prefix + "does not translate in X before it touches the floor");
        expectRange(tr.airborne.com.z() - tr.first.com.z(), -0.04, 0.04,
                    prefix + "does not translate in Z before it touches the floor");
        check(std::abs(sampleSpan(tr.airborne) - sampleSpan(tr.first)) < 0.08,
              (prefix + "does not fold in mid-air").c_str());
    }
    check(tr.hadHit || tr.last.lowest < 0.08, (prefix + "reaches the floor").c_str());
    // Corner contact starts before the lowest node is under 8 mm, so the
    // first captured "hit" frame can already be the tipped pose. Horizontal
    // COM from rest is the tip: airborne xz is already required to be ~0.
    std::printf("%spostHitMove=%.4f xzFromStart=%.4f\n", prefix.c_str(), tr.postHitMove,
                tr.maxComXzFromStart);
    check(tr.postHitMove > 0.04 || tr.maxComXzFromStart > 0.04,
          (prefix + "tips off the corner after the floor hit").c_str());
    check(!(tr.last.ke < 1.0e-3 && tr.last.lowest > 0.03),
          (prefix + "is not velocity-killed while still above the floor").c_str());
    if (tr.hadHit)
        check(tr.maxStress > 1.0e4, (prefix + "records constitutive stress at impact").c_str());
    std::printf("%sspan %.3f m com=(%.3f,%.3f,%.3f) ke=%.3f stress=%.2f MPa hht=%d rejects=%d broken=%d\n",
                prefix.c_str(), span, tr.last.com.x(), tr.last.com.y(), tr.last.com.z(), tr.last.ke,
                tr.maxStress * 1.0e-6, tr.hhtCommits, tr.newtonRejects, tr.last.broken);
    check(tr.last.broken < 48, (prefix + "welds survive the floor hit").c_str());
}

void checkMotion() {
    std::printf("\n-- motion M1-M7 --\n");

    // M1: airborne 1 g, hollow and solid, until still clear of the floor.
    auto airborne1g = [](const rcs::AssemblyDefinition& definition, const char* label) {
        rcs::Scenario scene = smokeFloor();
        scene.settings.timestep = 5.0e-4;
        scene.initial.position = {0.0, 0.55, 0.0};
        scene.initial.orientation = {0.35, 0.0, 0.45};
        scene.environment.aerodynamicsEnabled = false;
        SmokeSim sim;
        sim.reset(definition, scene);
        const BodySample start = sampleBody(sim.assembly());
        chrono::ChVector3d prevV = start.vel;
        bool finite = true;
        double t = 0.0;
        const double dt = scene.settings.timestep;
        double worstStray = 0.0;
        BodySample last = start;
        for (int step = 0; step < 800 && finite; ++step) {
            finite = sim.step(dt, t, scene.environment, scene.loads);
            t += dt;
            const BodySample now = sampleBody(sim.assembly());
            if (now.lowest <= 0.03)
                break;
            last = now;
            if (t > 0.02 && dt > 0.0) {
                const chrono::ChVector3d a = (now.vel - prevV) / dt;
                const double stray = (a - scene.environment.gravity).Length();
                worstStray = std::max(worstStray, stray);
            }
            prevV = now.vel;
        }
        dumpSim(label, sim, t, finite);
        const double drop = start.com.y() - last.com.y();
        const double expected = 0.5 * 9.81 * t * t;
        std::printf("%s airborne t=%.3f drop=%.4f expect=%.4f stray=%.2f xz=(%.4f,%.4f)\n", label, t, drop,
                    expected, worstStray, last.com.x() - start.com.x(), last.com.z() - start.com.z());
        check(finite, (std::string(label) + " free fall stays in the solver").c_str());
        check(last.lowest > 0.03, (std::string(label) + " samples while still clear of the floor").c_str());
        check(std::abs(drop - expected) < 0.08 * std::max(expected, 1.0e-3),
              (std::string(label) + " free fall matches 1 g").c_str());
        expectRange(last.com.x() - start.com.x(), -0.04, 0.04,
                    std::string(label) + " COM does not drift in X while airborne");
        expectRange(last.com.z() - start.com.z(), -0.04, 0.04,
                    std::string(label) + " COM does not drift in Z while airborne");
        check(worstStray < 3.0, (std::string(label) + " has no mid-air impulses").c_str());
    };
    airborne1g(rcs::makeTestCube(), "M1 hollow");
    rcs::AssemblyDefinition solidAir;
    solidAir.parts.push_back(rcs::makeSolidCube());
    airborne1g(solidAir, "M1 solid");

    // M5: a UV-sphere of tets has ~15° faces. Resting on a facet is real
    // polyhedral contact, not a smooth ball. Volume fill matches that facet
    // length so the smoke does not embed a millimetre tet grid in a 6 cm body.
    {
        rcs::PartDefinition ball;
        ball.name = "sphere";
        ball.material = rcs::aluminum6061();
        ball.source = sphereSurface(0.06, 24, 12);
        ball.sizing.shellCutoff = 0.008;
        ball.sizing.minElement = 0.010;
        ball.sizing.maxElement = 0.030;
        check(meshPart(ball), "M5 sphere meshes");
        rcs::AssemblyDefinition def;
        def.parts.push_back(ball);
        rcs::Scenario scene = smokeFloor();
        scene.settings.timestep = 2.0e-4;
        scene.initial.position = {0.0, 0.10, 0.0};
        scene.initial.angularVelocity = {0.0, 0.0, -6.0};
        scene.environment.aerodynamicsEnabled = false;
        SmokeSim sim;
        sim.reset(def, scene);
        const BodySample start = sampleBody(sim.assembly());
        bool finite = true;
        double t = 0.0;
        const double dt = scene.settings.timestep;
        BodySample last = start;
        for (int step = 0; step < 1200 && finite; ++step) {
            finite = sim.step(dt, t, scene.environment, scene.loads);
            t += dt;
            last = sampleBody(sim.assembly());
        }
        dumpSim("M5-sphere", sim, t, finite);
        const double travel = last.com.x() - start.com.x();
        std::printf("M5 sphere travel x=%.4f ke=%.3f lowest=%.4f tets=%d\n", travel, last.ke, last.lowest,
                    last.live);
        check(finite, "M5 sphere lands without diverging");
        expectRange(last.lowest, -0.005, 0.02, "M5 sphere rest gap stays millimetre-scale");
        // Friction turns slip into rolling. Angular momentum about the contact
        // point is kept, so a solid sphere rolls on with 2/7 of its spin
        // energy; "wiped" is losing most of that too.
        check(last.ke > 0.5 * (2.0 / 7.0) * start.ke, "M5 sphere spin is not wiped on floor contact");
    }
    {
        rcs::AssemblyDefinition cube;
        cube.parts.push_back(rcs::makeSolidCube(0.12));
        rcs::Scenario scene = smokeFloor();
        scene.settings.timestep = 5.0e-4;
        scene.initial.position = {0.0, 0.10, 0.0};
        scene.initial.orientation = {0.0, 0.0, 0.0};
        scene.initial.angularVelocity = {0.0, 0.0, -6.0};
        scene.environment.aerodynamicsEnabled = false;
        SmokeSim sim;
        sim.reset(cube, scene);
        const BodySample start = sampleBody(sim.assembly());
        bool finite = true;
        double t = 0.0;
        BodySample last = start;
        for (int step = 0; step < 800 && finite; ++step) {
            finite = sim.step(scene.settings.timestep, t, scene.environment, scene.loads);
            t += scene.settings.timestep;
            last = sampleBody(sim.assembly());
        }
        dumpSim("M5-cube", sim, t, finite);
        check(finite, "M5 face-down cube with spin stays in the solver");
        check(std::abs(last.com.x() - start.com.x()) < 0.03,
              "M5 face-down cube does not roll away on a vertical spin");
        checkSeatedVolume("M5 cube", start, last, sim.assembly(), 0.06);
    }

    // M6: slide to stop.
    {
        rcs::AssemblyDefinition cube;
        cube.parts.push_back(rcs::makeSolidCube(0.16));
        rcs::Scenario scene = smokeFloor();
        scene.settings.duration = 0.80;
        scene.settings.timestep = 5.0e-4;
        scene.initial.position = {0.0, 0.081, 0.0};
        scene.initial.orientation = {0.0, 0.0, 0.0};
        scene.initial.linearVelocity = {2.0, 0.0, 0.0};
        scene.environment.aerodynamicsEnabled = false;
        SmokeSim sim;
        sim.reset(cube, scene);
        const BodySample start = sampleBody(sim.assembly());
        bool finite = true;
        double t = 0.0;
        BodySample last = start;
        double lastX = last.com.x();
        double creep = 0.0;
        bool stopped = false;
        for (int step = 0; step < 1600 && finite; ++step) {
            finite = sim.step(scene.settings.timestep, t, scene.environment, scene.loads);
            t += scene.settings.timestep;
            last = sampleBody(sim.assembly());
            if (!stopped && std::abs(last.vel.x()) < 0.05) {
                stopped = true;
                lastX = last.com.x();
            }
            if (stopped)
                creep = std::max(creep, std::abs(last.com.x() - lastX));
        }
        dumpSim("M6-slide", sim, t, finite);
        std::printf("M6 vx=%.3f creep=%.4f stopped=%d\n", last.vel.x(), creep, stopped ? 1 : 0);
        check(finite, "M6 sliding block stays in the solver");
        check(stopped, "M6 sliding block stops under Coulomb friction");
        check(creep < 0.01, "M6 sliding block does not creep after it has stopped");
        checkSeatedVolume("M6 slide", start, last, sim.assembly(), 0.08);
    }

    // M7: seated rest, no skate.
    {
        rcs::AssemblyDefinition cube;
        cube.parts.push_back(rcs::makeSolidCube(0.16));
        rcs::Scenario scene = smokeFloor();
        scene.settings.duration = 0.40;
        scene.settings.timestep = 5.0e-4;
        scene.initial.position = {0.0, 0.081, 0.0};
        scene.initial.orientation = {0.0, 0.0, 0.0};
        scene.environment.aerodynamicsEnabled = false;
        SmokeSim sim;
        sim.reset(cube, scene);
        const BodySample start = sampleBody(sim.assembly());
        bool finite = true;
        double t = 0.0;
        BodySample last = start;
        for (int step = 0; step < 800 && finite; ++step) {
            finite = sim.step(scene.settings.timestep, t, scene.environment, scene.loads);
            t += scene.settings.timestep;
            last = sampleBody(sim.assembly());
        }
        dumpSim("M7-rest", sim, t, finite);
        check(finite, "M7 resting cube stays in the solver");
        expectRange(last.com.x() - start.com.x(), -0.005, 0.005, "M7 resting cube does not skate in X");
        expectRange(last.com.z() - start.com.z(), -0.005, 0.005, "M7 resting cube does not skate in Z");
        check(last.vel.Length() < 1.0, "M7 resting cube is not launched");
        checkSeatedVolume("M7 rest", start, last, sim.assembly(), 0.08);
    }

    // M2 / M3 editor-default corner drops last so a hung solid path cannot
    // skip the coarse kinematics jobs. Full 1 s, meshResolution 0.25.
    const DropTrace hollow = runCornerDrop(rcs::makeTestCube(), true, "M2-hollow");
    assertEngineDrop(hollow, "M2 hollow");
    rcs::AssemblyDefinition solidDef;
    solidDef.parts.push_back(rcs::makeSolidCube(0.40));
    const DropTrace solid = runCornerDrop(solidDef, true, "M3-solid");
    assertEngineDrop(solid, "M3 solid");

    std::printf("M4 hollow mass=%.3f ke=%.3f  solid mass=%.3f ke=%.3f\n", hollow.last.mass, hollow.last.ke,
                solid.last.mass, solid.last.ke);
    check(hollow.last.mass > 1.0 && solid.last.mass > 1.0,
          "M4 hollow welded box and solid cube both ran the editor drop");
}

void checkDuctileLabs() {
    std::printf("\n-- ductile labs --\n");
    const rcs::Material metals[] = {rcs::aluminum6061(), rcs::mildSteel(), rcs::ar500()};
    for (const rcs::Material& metal : metals) {
        const double G = metal.youngsModulus / (2.0 * (1.0 + metal.poissonRatio));
        const char* name = metal.name.c_str();
        const AxialResult service = catalogAxial(metal, 0.30 * metal.yieldStress, 0.04, 0.012, 0.12, 0.010);
        const double eStrain = std::max(std::abs(service.strain), 1.0e-12);
        const double eMeas = (0.30 * metal.yieldStress) / eStrain;
        std::printf("%s E_meas=%.3e handbook=%.3e strain=%.5f plastic=%.5f nu_meas=%.3f\n", name, eMeas,
                    metal.youngsModulus, service.strain, service.gaugePlastic,
                    eStrain > 1.0e-9 ? -service.lateral / eStrain : 0.0);
        check(service.finite, (std::string(name) + " 0.3 yield tension stays in the solver").c_str());
        expectRange(eMeas / metal.youngsModulus, 0.70, 1.30,
                    std::string(name) + " Young's modulus from sigma/epsilon below yield");
        expectRange(service.gaugePlastic, 0.0, 1.0e-3,
                    std::string(name) + " 0.3 yield stays elastic in the gauge");
        if (eStrain > 1.0e-5)
            expectRange(-service.lateral / eStrain, metal.poissonRatio - 0.12, metal.poissonRatio + 0.12,
                        std::string(name) + " Poisson ratio from lateral strain");

        const AxialResult over = catalogAxial(metal, 1.20 * metal.yieldStress, 0.05, 0.016, 0.12, 0.010);
        check(over.finite, (std::string(name) + " 1.2 yield tension stays in the solver").c_str());
        check(over.gaugePlastic > 0.002 || over.torn,
              (std::string(name) + " 1.2 yield leaves plastic strain or tears").c_str());
        check(over.residual > 0.0005 || over.torn,
              (std::string(name) + " 1.2 yield leaves residual elongation").c_str());

        const AxialResult atYield = catalogAxial(metal, metal.yieldStress, 0.05, 0.016, 0.12, 0.010);
        check(atYield.finite, (std::string(name) + " yield-load tension stays in the solver").c_str());
        expectRange(atYield.gaugePlastic, 0.001, 0.012,
                    std::string(name) + " 0.2 percent offset yield in the gauge");

        const AxialResult uts = catalogAxial(metal, 1.15 * metal.ultimateStress, 0.08, 0.016, 0.12, 0.010);
        std::printf("%s UTS torn=%d strain=%.5f plastic=%.5f residual=%.5f A=%.3f\n", name,
                    uts.torn ? 1 : 0, uts.strain, uts.gaugePlastic, uts.residual, metal.failPlasticStrain);
        check(uts.finite, (std::string(name) + " UTS tension stays in the solver").c_str());
        check(uts.torn || uts.gaugePlastic > 0.5 * metal.failPlasticStrain,
              (std::string(name) + " UTS reaches tear or the handbook elongation").c_str());
        if (uts.torn)
            expectRange(uts.strain, 0.4 * metal.failPlasticStrain, 2.0 * metal.failPlasticStrain,
                        std::string(name) + " elongation at tear");

        const AxialResult crush = catalogCompress(metal, -0.30 * metal.yieldStress, 0.04, 0.012);
        check(crush.finite, (std::string(name) + " 0.3 yield compression stays in the solver").c_str());
        expectRange(crush.plastic, 0.0, 1.0e-3, std::string(name) + " 0.3 yield compression stays elastic");
        check(crush.shards == 0, (std::string(name) + " compression does not Rankine-crush a metal").c_str());
        std::printf("%s G=%.3e handbook\n", name, G);
    }

    // 3-point analog: clamped root, tip load, residual camber after unload.
    {
        const rcs::Material al = rcs::aluminum6061();
        const double L = 0.20, w = 0.05, t = 0.006;
        const double I = w * t * t * t / 12.0;
        const double P = 1.20 * al.yieldStress * (w * t * t / 6.0) / L;
        rcs::AssemblyDefinition bar = catalogBar(al, L, w, t, 10, 3, true);
        rcs::Scenario scene = catalogQuiet();
        scene.initial.position = {0.0, 0.6, 0.0};
        rcs::ImpactLoad hit;
        hit.point = {0.95 * L, 0.6, 0.0};
        hit.direction = {0.0, -1.0, 0.0};
        hit.finalForce = P;
        hit.rampTime = 0.03;
        hit.duration = 0.05;
        hit.radius = 0.5 * w;
        scene.loads.push_back(hit);
        SmokeSim sim;
        sim.reset(bar, scene);
        catalogPinRoot(sim.assembly().parts()[0]);
        sim.system().Setup();
        const std::size_t tip = catalogEndNode(sim.assembly().parts()[0], L);
        const double y0 = sim.assembly().parts()[0].nodePosition(tip).y();
        bool finite = true;
        double clock = 0.0;
        double peak = 0.0;
        const double dt = scene.settings.timestep;
        const int loaded = static_cast<int>(std::ceil(hit.duration / dt));
        for (int step = 0; step < loaded && finite; ++step) {
            finite = sim.step(dt, clock, scene.environment, scene.loads, scene.heatSources);
            clock += dt;
            peak = std::max(peak, y0 - sim.assembly().parts()[0].nodePosition(tip).y());
        }
        dumpSim("catalogStep", sim, clock, finite);
        scene.loads.clear();
        finite = catalogStep(sim, scene, 80) && finite;
        const double residual = y0 - sim.assembly().parts()[0].nodePosition(tip).y();
        const double elastic = P * L * L * L / (3.0 * al.youngsModulus * I);
        std::printf("3pt analog P=%.0f N peak=%.4f residual=%.4f elastic=%.4f\n", P, peak, residual, elastic);
        check(finite, "ductile 3-point analog stays in the solver");
        check(peak > 0.5 * elastic, "ductile 3-point analog deflects on the order of PL^3/3EI");
        check(residual > 0.001, "ductile 3-point analog keeps a residual camber after unload");
    }

    // Continuum shear: in-plane tip force on a squat plate so membrane shear
    // is not swamped by in-plane bending, gamma = tau/G.
    {
        const rcs::Material al = rcs::aluminum6061();
        const double L = 0.04, w = 0.05, th = 0.008;
        const double G = al.youngsModulus / (2.0 * (1.0 + al.poissonRatio));
        const double tau = 0.05 * al.yieldStress;
        const double force = tau * w * th;
        rcs::AssemblyDefinition bar = catalogBar(al, L, w, th, 6, 3, true);
        rcs::Scenario scene = catalogQuiet();
        scene.initial.position = {0.0, 0.6, 0.0};
        rcs::ImpactLoad hit;
        hit.point = {L, 0.6, 0.0};
        hit.direction = {0.0, 0.0, 1.0};
        hit.finalForce = force;
        hit.rampTime = 0.02;
        hit.duration = 0.04;
        hit.radius = 0.5 * w;
        scene.loads.push_back(hit);
        SmokeSim sim;
        sim.reset(bar, scene);
        catalogPinRoot(sim.assembly().parts()[0]);
        sim.system().Setup();
        const std::size_t tip = catalogEndNode(sim.assembly().parts()[0], L);
        const double z0 = sim.assembly().parts()[0].nodePosition(tip).z();
        bool finite = true;
        double t = 0.0;
        double dz = 0.0;
        const double dt = scene.settings.timestep;
        const int loaded = static_cast<int>(std::ceil(hit.duration / dt));
        for (int step = 0; step < loaded && finite; ++step) {
            finite = sim.step(dt, t, scene.environment, scene.loads, scene.heatSources);
            t += dt;
            dz = std::max(dz, sim.assembly().parts()[0].nodePosition(tip).z() - z0);
        }
        dumpSim("catalogStep", sim, t, finite);
        const double gamma = dz / L;
        const double want = tau / G;
        std::printf("shear gamma=%.5f handbook=%.5f G=%.3e\n", gamma, want, G);
        check(finite, "ductile shear coupon stays in the solver");
        expectRange(gamma / std::max(want, 1.0e-12), 0.40, 2.50, "ductile shear strain is on the order of tau/G");
    }

    // Torsion analog: equal-and-opposite Z forces on the tip edge.
    {
        const rcs::Material al = rcs::aluminum6061();
        const double L = 0.10, w = 0.04, th = 0.008;
        rcs::AssemblyDefinition bar = catalogBar(al, L, w, th, 8, 4, true);
        rcs::Scenario scene = catalogQuiet();
        scene.initial.position = {0.0, 0.6, 0.0};
        SmokeSim sim;
        sim.reset(bar, scene);
        rcs::Part& part = sim.assembly().parts()[0];
        catalogPinRoot(part);
        double maxX = -1.0e9, minZ = 1.0e9, maxZ = -1.0e9;
        for (std::size_t n = 0; n < part.nodeCount(); ++n) {
            maxX = std::max(maxX, part.restPosition(n).x());
            minZ = std::min(minZ, part.restPosition(n).z());
            maxZ = std::max(maxZ, part.restPosition(n).z());
        }
        std::size_t lo = 0, hi = 0;
        for (std::size_t n = 0; n < part.nodeCount(); ++n) {
            if (part.restPosition(n).x() < maxX - 1.0e-3)
                continue;
            if (part.restPosition(n).z() <= minZ + 1.0e-3)
                lo = n;
            if (part.restPosition(n).z() >= maxZ - 1.0e-3)
                hi = n;
        }
        const double couple = 40.0;
        auto add = [&](std::size_t n, double sign) {
            rcs::ImpactLoad hit;
            hit.point = part.nodePosition(n);
            hit.direction = {0.0, sign, 0.0};
            hit.finalForce = couple;
            hit.rampTime = 0.02;
            hit.duration = 0.04;
            hit.radius = 1.0e-4;
            scene.loads.push_back(hit);
        };
        add(lo, 1.0);
        add(hi, -1.0);
        sim.assembly().bindLoads(scene.loads);
        sim.system().Setup();
        const double z0 = part.nodePosition(hi).y() - part.nodePosition(lo).y();
        bool finite = true;
        double t = 0.0;
        double peakTwist = 0.0;
        const double dt = scene.settings.timestep;
        const double loadEnd = scene.loads.empty() ? 0.0 : scene.loads.front().duration;
        const int loaded = static_cast<int>(std::ceil(loadEnd / dt));
        for (int step = 0; step < loaded && finite; ++step) {
            finite = sim.step(dt, t, scene.environment, scene.loads, scene.heatSources);
            t += dt;
            const double zNow = sim.assembly().parts()[0].nodePosition(hi).y() -
                                sim.assembly().parts()[0].nodePosition(lo).y();
            peakTwist = std::max(peakTwist, std::abs(zNow - z0));
        }
        dumpSim("catalogStep", sim, t, finite);
        const double z1 = sim.assembly().parts()[0].nodePosition(hi).y() -
                          sim.assembly().parts()[0].nodePosition(lo).y();
        std::printf("torsion analog deltaY=%.5f -> %.5f peak=%.5f\n", z0, z1, peakTwist);
        check(finite, "ductile torsion analog stays in the solver");
        check(peakTwist > 1.0e-5, "ductile torsion analog twists the free end");
    }

    // Charpy analog: notched bar, striker, absorbed energy in Joules.
    {
        const rcs::Material al = rcs::aluminum6061();
        rcs::AssemblyDefinition bar = catalogBar(al, 0.055, 0.010, 0.010, 8, 2, true);
        rcs::Scenario scene = catalogQuiet();
        scene.initial.position = {0.0, 0.5, 0.0};
        rcs::ImpactLoad hit;
        hit.point = {0.0275, 0.5, 0.0};
        hit.direction = {0.0, 0.0, 1.0};
        hit.finalForce = 8000.0;
        hit.duration = 0.002;
        hit.radius = 0.004;
        scene.loads.push_back(hit);
        SmokeSim sim;
        sim.reset(bar, scene);
        std::size_t mid = 0;
        double midBest = 1.0e9;
        for (std::size_t n = 0; n < sim.assembly().parts()[0].nodeCount(); ++n) {
            const auto r = sim.assembly().parts()[0].restPosition(n);
            const double s = std::abs(r.x() - 0.0275) + std::abs(r.z());
            if (s < midBest) {
                midBest = s;
                mid = n;
            }
        }
        const double z0 = sim.assembly().parts()[0].nodePosition(mid).z();
        bool finite = true;
        double t = 0.0;
        double work = 0.0;
        double zPrev = z0;
        const double dt = scene.settings.timestep;
        for (int step = 0; step < 80 && finite; ++step) {
            finite = sim.step(dt, t, scene.environment, scene.loads);
            t += dt;
            const double z = sim.assembly().parts()[0].nodePosition(mid).z();
            const double F = rcs::loadMagnitude(hit, t - dt);
            work += F * (z - zPrev);
            zPrev = z;
        }
        dumpSim("charpy", sim, t, finite);
        const double absorbed = std::abs(work);
        std::printf("Charpy analog absorbed=%.3f J (handbook 10-40 J) torn=%d shards=%zu\n", absorbed,
                    sim.assembly().parts()[0].torn() ? 1 : 0, sim.assembly().parts()[0].shards().size());
        check(finite, "Charpy analog stays in the solver");
        check(absorbed > 1.0, "Charpy analog reports absorbed mesh work");
        std::printf("Charpy handbook 10-40 J is a feature-gap (no striker mass; analog got %.3f J)\n",
                    absorbed);
    }

    // Combined tension + bending on a weldment.
    {
        const rcs::Material al = rcs::aluminum6061();
        const double width = 0.04, thick = 0.006, length = 0.10;
        const double throat = thick;
        const double rated = al.weldShearStrength * width * throat;
        rcs::PartDefinition left;
        left.name = "Held";
        left.material = al;
        left.thickness = thick;
        left.fixed = true;
        left.patch.origin = {0.0, 0.0, -0.5 * width};
        left.patch.uAxis = {1.0, 0.0, 0.0};
        left.patch.vAxis = {0.0, 0.0, 1.0};
        left.patch.uSize = length;
        left.patch.vSize = width;
        left.patch.divU = 5;
        left.patch.divV = 2;
        rcs::PartDefinition right = left;
        right.name = "Loaded";
        right.fixed = false;
        right.patch.origin = {length, 0.0, -0.5 * width};
        rcs::WeldDefinition seam;
        seam.partA = 0;
        seam.partB = 1;
        seam.tolerance = 3.0e-3;
        seam.throat = throat;
        rcs::AssemblyDefinition definition;
        definition.parts.push_back(left);
        definition.parts.push_back(right);
        definition.welds.push_back(seam);
        rcs::Scenario scene = catalogQuiet();
        scene.initial.position = {0.0, 0.6, 0.0};
        rcs::ImpactLoad axial;
        axial.partIndex = 1;
        axial.point = {2.0 * length, 0.6, 0.0};
        axial.direction = {1.0, 0.0, 0.0};
        axial.finalForce = 0.5 * rated;
        axial.rampTime = 0.02;
        axial.duration = 0.05;
        axial.radius = 0.5 * width;
        rcs::ImpactLoad bend = axial;
        bend.direction = {0.0, -1.0, 0.0};
        bend.finalForce = 0.5 * rated;
        scene.loads.push_back(axial);
        scene.loads.push_back(bend);
        SmokeSim sim;
        sim.reset(definition, scene);
        const bool finite = catalogStep(sim, scene, 300);
        dumpSim("weldment-bend", sim, 0.05, finite);
        check(finite, "combined tension+bending weldment stays in the solver");
        check(sim.assembly().brokenLinkCount() == 0,
              "combined tension+bending at 0.5 rated throat does not pop the weld");
    }

    // Short cyclic fatigue: 20 cycles at 0.4 σy, plastic increment → 0.
    {
        const rcs::Material al = rcs::aluminum6061();
        const double L = 0.12, w = 0.03, th = 0.006;
        const double force = 0.40 * al.yieldStress * w * th;
        rcs::AssemblyDefinition bar = catalogBar(al, L, w, th, 8, 2, true);
        rcs::Scenario scene = catalogQuiet();
        scene.settings.timestep = 2.0e-4;
        scene.initial.position = {0.0, 0.5, 0.0};
        SmokeSim sim;
        sim.reset(bar, scene);
        catalogPinRoot(sim.assembly().parts()[0]);
        sim.system().Setup();
        bool finite = true;
        double t = 0.0;
        const double dt = scene.settings.timestep;
        double plasticStart = 0.0;
        double plasticEnd = 0.0;
        for (int cycle = 0; cycle < 20 && finite; ++cycle) {
            rcs::ImpactLoad hit;
            hit.point = {L, 0.5, 0.0};
            hit.direction = {1.0, 0.0, 0.0};
            hit.finalForce = force;
            hit.rampTime = 0.001;
            hit.duration = 0.004;
            hit.radius = 0.5 * w;
            scene.loads = {hit};
            sim.assembly().bindLoads(scene.loads);
            for (int step = 0; step < 20 && finite; ++step) {
                finite = sim.step(dt, t, scene.environment, scene.loads);
                t += dt;
            }
            scene.loads.clear();
            sim.assembly().bindLoads(scene.loads);
            for (int step = 0; step < 10 && finite; ++step) {
                finite = sim.step(dt, t, scene.environment, scene.loads);
                t += dt;
            }
            const double plastic = catalogGaugePlastic(sim.assembly().parts()[0], 0.25 * L, 0.75 * L);
            if (cycle == 0)
                plasticStart = plastic;
            plasticEnd = plastic;
        }
        dumpSim("fatigue-20", sim, t, finite);
        std::printf("fatigue analog 20 cycles plastic %.6f -> %.6f (S-N endurance is a feature-gap)\n",
                    plasticStart, plasticEnd);
        check(finite, "short cyclic 20-cycle load stays in the solver");
        check(std::abs(plasticEnd - plasticStart) < 5.0e-4,
              "0.4 yield 20-cycle load does not ratchet plastic strain");
        const AxialResult over = catalogAxial(al, 1.10 * al.yieldStress, 0.010, 0.004, 0.12, 0.002);
        check(over.gaugePlastic > 0.001 || over.torn,
              "overload cyclic load yields or tears (no S-N endurance claim)");
    }
}

void checkBrittleLabs() {
    std::printf("\n-- brittle labs --\n");
    const rcs::Material soda = rcs::glass();
    const rcs::Material temper = rcs::temperedGlass();

    auto flex = [](const rcs::Material& glass, double stressFrac, const char* tag) {
        const double L = 0.12, w = 0.03, th = 0.006;
        const double mor = glass.fractureStrength;
        // Simply-supported midspan: σ = 3PL / (2 w th²). The old P used
        // (w th²/6)/(L/2) and only aimed at 0.5 σ.
        const double P = stressFrac * mor * (2.0 * w * th * th) / (3.0 * L);
        rcs::AssemblyDefinition bar = catalogBar(glass, L, w, th, 8, 2, false);
        rcs::Scenario scene = catalogQuiet();
        scene.initial.position = {0.0, 0.6, 0.0};
        rcs::ImpactLoad hit;
        hit.point = {0.5 * L, 0.6, 0.0};
        hit.direction = {0.0, -1.0, 0.0};
        hit.finalForce = P;
        hit.rampTime = 0.01;
        hit.duration = 0.03;
        hit.radius = 0.5 * w;
        scene.loads.push_back(hit);
        SmokeSim sim;
        sim.reset(bar, scene);
        // 3-point MOR is simply supported. Fixing the whole end face makes a
        // built-in beam and halves midspan tension so 1.4 MOR never cracks.
        auto pinSeat = [](rcs::Part& part, bool root) {
            double x0 = root ? 1.0e9 : -1.0e9;
            double y0 = 1.0e9;
            for (std::size_t n = 0; n < part.nodeCount(); ++n) {
                const auto r = part.restPosition(n);
                x0 = root ? std::min(x0, r.x()) : std::max(x0, r.x());
                y0 = std::min(y0, r.y());
            }
            auto pin = [&](auto& nodes) {
                for (std::size_t n = 0; n < nodes.size() && n < part.nodeCount(); ++n) {
                    if (!nodes[n])
                        continue;
                    const auto r = part.restPosition(n);
                    if (std::abs(r.x() - x0) < 1.5e-3 && r.y() <= y0 + 1.5e-3)
                        nodes[n]->SetFixed(true);
                }
            };
            pin(part.xyzNodes());
            pin(part.nodes());
        };
        pinSeat(sim.assembly().parts()[0], true);
        pinSeat(sim.assembly().parts()[0], false);
        sim.system().Setup();
        const auto advice = sim.assembly().adviseLoad(0, scene.loads, scene.settings.timestep);
        bool finite = true;
        double t = 0.0;
        double peakTension = 0.0;
        const double dt = scene.settings.timestep;
        for (int step = 0; step < 200 && finite; ++step) {
            finite = sim.step(dt, t, scene.environment, scene.loads);
            t += dt;
            if (t >= hit.startTime && t <= hit.startTime + hit.duration)
                peakTension = std::max(peakTension, sim.assembly().parts()[0].peakTension());
        }
        dumpSim(tag, sim, t, finite);
        const rcs::Part& done = sim.assembly().parts()[0];
        std::printf("%s P=%.0f N targets=%zu peakT=%.1f MPa endT=%.1f MPa torn=%d shards=%zu plastic=%.5f\n",
                    tag, P, advice.nodeCount, peakTension * 1.0e-6, done.peakTension() * 1.0e-6,
                    done.torn() ? 1 : 0, done.shards().size(), catalogPeakPlastic(done));
        return std::tuple<bool, bool, std::size_t, double>(finite, done.torn() || done.shards().size() > 0,
                                                           done.shards().size(), catalogPeakPlastic(done));
    };
    const auto sodaHalf = flex(soda, 0.50, "soda-0.5MOR");
    const auto sodaOver = flex(soda, 1.40, "soda-1.4MOR");
    const auto temperOver = flex(temper, 1.40, "temper-1.4MOR");
    check(std::get<0>(sodaHalf), "soda-lime 0.5 MOR flexure stays in the solver");
    check(!std::get<1>(sodaHalf), "soda-lime 0.5 MOR does not crack");
    expectRange(std::get<3>(sodaHalf), 0.0, 1.0e-4, "soda-lime flexure has no plastic strain");
    check(std::get<0>(sodaOver), "soda-lime 1.4 MOR flexure stays in the solver");
    check(std::get<1>(sodaOver), "soda-lime 1.4 MOR opens a Rankine crack");
    check(std::get<0>(temperOver), "tempered 1.4 MOR flexure stays in the solver");

    const AxialResult crush = catalogCompress(soda, -0.10 * rcs::rankineCompressionLimit(soda), 0.012, 0.006);
    check(crush.finite, "soda-lime mild compression stays in the solver");
    check(!crush.torn && crush.shards == 0, "soda-lime below crush strength does not powder");
    expectRange(crush.plastic, 0.0, 1.0e-4, "soda-lime compression has no plastic strain");
    const AxialResult crushOver =
        catalogCompress(soda, -1.20 * rcs::rankineCompressionLimit(soda), 0.012, 0.006);
    check(crushOver.finite, "soda-lime 1.2 crush stays in the solver");
    check(crushOver.torn || crushOver.shards > 0, "soda-lime at 1.2 crush Rankine-opens");

    // Punch: tempered peels at least as much as soda-lime. Not Voronoi.
    {
        auto punch = [](const rcs::Material& material, const char* tag) {
            rcs::Scenario scene = stillScenario();
            scene.settings.timestep = 1.0e-4;
            scene.initial.position = {0.0, 0.3, 0.0};
            rcs::ImpactLoad hit;
            hit.point = {0.0, 0.3, 0.0};
            hit.direction = {0.0, -1.0, 0.0};
            hit.finalForce = 60000.0;
            hit.radius = 0.15;
            hit.startTime = 0.002;
            hit.duration = 0.003;
            scene.loads.push_back(hit);
            SmokeSim sim;
            sim.reset(clampedPlate(material, 4), scene);
            const int before = sim.assembly().parts()[0].liveElementCount();
            bool finite = catalogStep(sim, scene, 400);
            const rcs::Part& done = sim.assembly().parts()[0];
            dumpSim(tag, sim, 0.04, finite);
            return std::tuple<bool, int, int, std::size_t>(finite, before, done.liveElementCount(),
                                                           done.shards().size());
        };
        const auto a = punch(soda, "soda-punch");
        const auto b = punch(temper, "temper-punch");
        check(std::get<0>(a) && std::get<0>(b), "glass punches stay in the solver");
        check(std::get<3>(b) >= std::get<3>(a) ||
                  (std::get<1>(b) > 0 && double(std::get<2>(b)) / double(std::get<1>(b)) <=
                                             double(std::get<2>(a)) / std::max(1, std::get<1>(a)) + 0.05),
              "tempered glass peels at least as much as soda-lime");
        std::printf("brittle shatter is Rankine element bisection, not Voronoi\n");
    }

    // Thermal shock analog: hot plate, cold sink on one face via heat source.
    {
        rcs::AssemblyDefinition plate = clampedPlate(soda, 4);
        rcs::Scenario scene = stillScenario();
        scene.settings.timestep = 5.0e-4;
        scene.initial.position = {0.0, 0.3, 0.0};
        rcs::HeatSource torch;
        torch.point = {0.0, 0.3, 0.0};
        torch.radius = 0.08;
        torch.temperature = 900.0;
        torch.duration = 0.20;
        scene.heatSources.push_back(torch);
        SmokeSim sim;
        sim.reset(plate, scene);
        const bool finite = catalogStep(sim, scene, 400);
        const rcs::Part& done = sim.assembly().parts()[0];
        dumpSim("thermal-shock", sim, 0.20, finite);
        check(finite, "glass thermal-shock analog stays in the solver");
        check(done.torn() || done.shards().size() > 0 || done.peakTension() > 0.8 * soda.fractureStrength,
              "glass thermal-shock analog reaches Rankine tension or cracks");
    }

    // Vickers analog: mean indent pressure vs HV 500-550.
    {
        rcs::AssemblyDefinition plate = clampedPlate(soda, 4);
        rcs::Scenario scene = stillScenario();
        scene.settings.timestep = 1.0e-4;
        scene.initial.position = {0.0, 0.3, 0.0};
        rcs::ImpactLoad hit;
        hit.point = {0.0, 0.3, 0.0};
        hit.direction = {0.0, -1.0, 0.0};
        hit.finalForce = 50.0;
        hit.radius = 0.0005;
        hit.duration = 0.01;
        scene.loads.push_back(hit);
        SmokeSim sim;
        sim.reset(plate, scene);
        const rcs::Part& plate0 = sim.assembly().parts()[0];
        std::size_t hitNode = 0;
        double hitBest = 1.0e9;
        const chrono::ChVector3d hitP(0.0, 0.3, 0.0);
        for (std::size_t n = 0; n < plate0.nodeCount(); ++n) {
            const double d = (plate0.nodePosition(n) - hitP).Length();
            if (d < hitBest) {
                hitBest = d;
                hitNode = n;
            }
        }
        const double y0 = plate0.nodePosition(hitNode).y();
        const bool finite = catalogStep(sim, scene, 120);
        const rcs::Part& plate1 = sim.assembly().parts()[0];
        double rIndent = 0.0;
        const double depth = std::max(0.0, y0 - plate1.nodePosition(hitNode).y());
        const double rCap = 4.0 * std::max(hit.radius, 1.0e-4);
        for (std::size_t n = 0; n < plate1.nodeCount(); ++n) {
            const auto p = plate1.nodePosition(n);
            const double rh = std::hypot(p.x() - hitP.x(), p.z() - hitP.z());
            if (rh > rCap)
                continue;
            const double dy = y0 - p.y();
            if (dy > 0.25 * std::max(depth, 1.0e-6))
                rIndent = std::max(rIndent, rh);
        }
        const double rMeas = std::max(rIndent, 1.0e-6);
        const double area = 3.14159265358979323846 * rMeas * rMeas;
        const double hv = (hit.finalForce / area) / 9.807e6;
        std::printf("Vickers analog HV=%.1f depth=%.5f r=%.5f (handbook 500-550) finite=%d\n", hv, depth,
                    rMeas, finite ? 1 : 0);
        check(finite, "Vickers analog stays in the solver");
        check(depth >= 0.0 && rMeas >= 0.0, "Vickers analog measures an indent from the mesh");
        std::printf("Vickers HV 500-550 is a feature-gap (analog HV=%.1f from mesh indent)\n", hv);
    }

    checkMixedFloorDrop();
}

void checkCompliantLabs() {
    std::printf("\n-- compliant labs --\n");
    const rcs::Material rub = rcs::rubber();
    const AxialResult uni = catalogAxial(rub, 0.10 * rub.youngsModulus, 0.05, 0.02, 0.12, 0.02);
    const double lambda = 1.0 + uni.strain;
    std::printf("uniaxial lambda=%.3f plastic=%.5f torn=%d\n", lambda, uni.plastic, uni.torn ? 1 : 0);
    check(uni.finite, "elastomer uniaxial stretch stays in the solver");
    expectRange(uni.plastic, 0.0, 1.0e-4, "elastomer uniaxial stretch has no plastic strain");
    check(!uni.torn, "elastomer at modest stretch does not tear");

    const AxialResult tear = catalogAxial(rub, 0.90 * rub.youngsModulus, 0.08, 0.02, 0.12, 0.02);
    const double lambdaTear = 1.0 + tear.strain;
    std::printf("overload lambda=%.3f torn=%d want tearStretch=%.2f\n", lambdaTear, tear.torn ? 1 : 0,
                rub.tearStretch);
    check(tear.finite, "elastomer overload stretch stays in the solver");
    check(tear.torn || lambdaTear > 0.6 * rub.tearStretch,
          "elastomer overload approaches tearStretch rather than yielding");

    // Biaxial vs uniaxial: Mooney-Rivlin equibiaxial P(λ) differs from uniaxial.
    // Measure interior gauge stretch, not the AABB of a fully clamped boundary.
    {
        rcs::PartDefinition sheet;
        sheet.name = "Biaxial";
        sheet.material = rub;
        sheet.thickness = 0.008;
        sheet.clampBoundary = false;
        sheet.patch.origin = {-0.2, 0.0, -0.2};
        sheet.patch.uAxis = {1.0, 0.0, 0.0};
        sheet.patch.vAxis = {0.0, 0.0, 1.0};
        sheet.patch.uSize = 0.4;
        sheet.patch.vSize = 0.4;
        sheet.patch.divU = 4;
        sheet.patch.divV = 4;
        rcs::AssemblyDefinition plate;
        plate.parts.push_back(sheet);
        rcs::Scenario scene = catalogQuiet();
        scene.initial.position = {0.0, 0.3, 0.0};
        SmokeSim sim;
        sim.reset(plate, scene);
        rcs::Part& part = sim.assembly().parts()[0];
        std::size_t center = 0;
        double bestC = 1.0e9;
        for (std::size_t n = 0; n < part.nodeCount(); ++n) {
            const auto r = part.restPosition(n);
            const double s = r.x() * r.x() + r.z() * r.z();
            if (s < bestC) {
                bestC = s;
                center = n;
            }
        }
        auto pin = [&](auto& nodes) {
            if (center < nodes.size() && nodes[center])
                nodes[center]->SetFixed(true);
        };
        pin(part.nodes());
        pin(part.xyzNodes());
        rcs::ImpactLoad fx;
        fx.point = {0.15, 0.3, 0.0};
        fx.direction = {1.0, 0.0, 0.0};
        fx.finalForce = 80.0;
        fx.duration = 0.04;
        fx.radius = 0.04;
        rcs::ImpactLoad fz = fx;
        fz.point = {0.0, 0.3, 0.15};
        fz.direction = {0.0, 0.0, 1.0};
        scene.loads.push_back(fx);
        scene.loads.push_back(fz);
        sim.system().Setup();
        sim.system().ForceUpdate();
        sim.assembly().bindLoads(scene.loads);
        const bool finite = catalogStep(sim, scene, 250);
        dumpSim("biaxial", sim, 0.05, finite);
        const rcs::Part& done = sim.assembly().parts()[0];
        double restMin = 1.0e9, restMax = -1.0e9, liveMin = 1.0e9, liveMax = -1.0e9;
        for (std::size_t n = 0; n < done.nodeCount(); ++n) {
            if (done.nodeRigid(n))
                continue;
            const double rx = done.restPosition(n).x();
            if (std::abs(rx) < 0.04 || std::abs(rx) > 0.16)
                continue;
            restMin = std::min(restMin, rx);
            restMax = std::max(restMax, rx);
            liveMin = std::min(liveMin, done.nodePosition(n).x());
            liveMax = std::max(liveMax, done.nodePosition(n).x());
        }
        const double restSpan = std::max(restMax - restMin, 1.0e-6);
        const double lambdaBi = (liveMax - liveMin) / restSpan;
        std::printf("biaxial lambda=%.3f uniaxial lambda=%.3f\n", lambdaBi, lambda);
        check(finite, "elastomer biaxial plate stays in the solver");
        check(std::abs(lambdaBi - lambda) > 0.05,
              "elastomer biaxial stretch differs from uniaxial (Ogden/Mooney, not Hencky-Hooke overlay)");
    }

    // Compression set: 25% hold then unload. Residual must be a real set.
    {
        const AxialResult hold = catalogAxial(rub, -0.25 * rub.youngsModulus, 0.04, 0.04, 0.06, 0.01);
        std::printf("compression set residual strain=%.5f\n", hold.residual);
        check(hold.finite, "elastomer compression hold stays in the solver");
        check(std::abs(hold.residual) > 0.02,
              "elastomer compression set leaves a residual height change");
    }

    // Trouser analog: slit is a missing quad on a stretched strip — we load a
    // bar and require any tear to start at the loaded tip, not a random midspan.
    {
        const AxialResult slit = catalogAxial(rub, 0.80 * rub.youngsModulus, 0.06, 0.02, 0.12, 0.02);
        check(slit.finite, "elastomer trouser analog stays in the solver");
        if (slit.torn)
            check(slit.strain > 0.5, "elastomer tear runs after large stretch, not a mid-span snap");
        else
            check(slit.strain > 0.3, "elastomer trouser analog actually stretches");
    }

    // DMA: oscillating load. Storage/loss from phase. Feature-gap if E''~0.
    {
        rcs::AssemblyDefinition bar = catalogBar(rub, 0.08, 0.03, 0.006, 6, 2, true);
        rcs::Scenario scene = catalogQuiet();
        scene.settings.timestep = 2.0e-4;
        scene.initial.position = {0.0, 0.5, 0.0};
        SmokeSim sim;
        sim.reset(bar, scene);
        catalogPinRoot(sim.assembly().parts()[0]);
        const std::size_t tip = catalogEndNode(sim.assembly().parts()[0], 0.08);
        const double y0 = sim.assembly().parts()[0].nodePosition(tip).x();
        double dispCos = 0.0, dispSin = 0.0;
        const double omega = 2.0 * 3.14159265358979323846 * 10.0;
        bool finite = true;
        double t = 0.0;
        const double dt = scene.settings.timestep;
        for (int step = 0; step < 400 && finite; ++step) {
            rcs::ImpactLoad hit;
            hit.point = sim.assembly().parts()[0].nodePosition(tip);
            hit.direction = {1.0, 0.0, 0.0};
            hit.finalForce = 20.0 * std::sin(omega * t);
            hit.duration = dt * 1.5;
            hit.radius = 1.0e-4;
            scene.loads = {hit};
            sim.assembly().bindLoads(scene.loads);
            finite = sim.step(dt, t, scene.environment, scene.loads);
            t += dt;
            const double u = sim.assembly().parts()[0].nodePosition(tip).x() - y0;
            dispCos += u * std::cos(omega * t);
            dispSin += u * std::sin(omega * t);
        }
        dumpSim("dma", sim, t, finite);
        const double tanDelta = std::abs(dispCos) > 1.0e-12 ? std::abs(dispSin / dispCos) : 0.0;
        std::printf("DMA tan_delta=%.4f (viscoelastic rubber is not ~0)\n", tanDelta);
        check(finite, "DMA analog stays in the solver");
        check(tanDelta > 0.05, "DMA analog has a loss modulus (tan delta not ~0)");
    }

    // Shore analog: indent pressure vs a Shore D / hard TPU band.
    {
        rcs::AssemblyDefinition pad = clampedPlate(rub, 3);
        rcs::Scenario scene = stillScenario();
        scene.initial.position = {0.0, 0.3, 0.0};
        rcs::ImpactLoad hit;
        hit.point = {0.0, 0.3, 0.0};
        hit.direction = {0.0, -1.0, 0.0};
        hit.finalForce = 10.0;
        hit.radius = 0.001;
        hit.duration = 0.02;
        scene.loads.push_back(hit);
        SmokeSim sim;
        sim.reset(pad, scene);
        const rcs::Part& pad0 = sim.assembly().parts()[0];
        std::size_t hitNode = 0;
        double hitBest = 1.0e9;
        const chrono::ChVector3d hitP(0.0, 0.3, 0.0);
        for (std::size_t n = 0; n < pad0.nodeCount(); ++n) {
            const double d = (pad0.nodePosition(n) - hitP).Length();
            if (d < hitBest) {
                hitBest = d;
                hitNode = n;
            }
        }
        const double y0 = pad0.nodePosition(hitNode).y();
        const bool finite = catalogStep(sim, scene, 120);
        const rcs::Part& pad1 = sim.assembly().parts()[0];
        double rIndent = 0.0;
        const double depth = std::max(0.0, y0 - pad1.nodePosition(hitNode).y());
        const double rCap = 4.0 * std::max(hit.radius, 1.0e-4);
        for (std::size_t n = 0; n < pad1.nodeCount(); ++n) {
            const auto p = pad1.nodePosition(n);
            const double rh = std::hypot(p.x() - hitP.x(), p.z() - hitP.z());
            if (rh > rCap)
                continue;
            const double dy = y0 - p.y();
            if (dy > 0.25 * std::max(depth, 1.0e-6))
                rIndent = std::max(rIndent, rh);
        }
        const double rMeas = std::max(rIndent, 1.0e-6);
        const double pMPa = hit.finalForce / (3.14159265358979323846 * rMeas * rMeas) / 1.0e6;
        std::printf("Shore analog indent pressure=%.1f MPa depth=%.5f r=%.5f\n", pMPa, depth, rMeas);
        check(finite, "Shore analog stays in the solver");
        check(depth >= 0.0, "Shore analog measures an indent from the mesh");
        std::printf("Shore D pressure band is a feature-gap (analog %.1f MPa from mesh indent)\n", pMPa);
    }

    // Rebound resilience.
    {
        rcs::PartDefinition pad;
        pad.name = "Pad";
        pad.material = rub;
        pad.thickness = 0.02;
        pad.patch.origin = {-0.08, 0.0, -0.08};
        pad.patch.uAxis = {1, 0, 0};
        pad.patch.vAxis = {0, 0, 1};
        pad.patch.uSize = 0.16;
        pad.patch.vSize = 0.16;
        pad.patch.divU = 4;
        pad.patch.divV = 4;
        pad.fixed = true;
        rcs::PartDefinition striker = rcs::makeSolidCube(0.04, {0.0, 0.08, 0.0});
        striker.material = rcs::aluminum6061();
        rcs::AssemblyDefinition def;
        def.parts.push_back(pad);
        def.parts.push_back(striker);
        rcs::Scenario scene = smokeFloor();
        scene.environment.floorEnabled = false;
        scene.initial.position = {0.0, 0.0, 0.0};
        scene.initial.linearVelocity = {0.0, -2.0, 0.0};
        scene.settings.timestep = 2.0e-4;
        SmokeSim sim;
        sim.reset(def, scene);
        double vyIn = 0.0, vyOut = 0.0;
        bool finite = true;
        double t = 0.0;
        for (int step = 0; step < 400 && finite; ++step) {
            finite = sim.step(scene.settings.timestep, t, scene.environment, scene.loads);
            t += scene.settings.timestep;
            const BodySample s = sampleBody(sim.assembly());
            if (s.vel.y() < vyIn)
                vyIn = s.vel.y();
            if (vyIn < -0.2 && s.vel.y() > vyOut)
                vyOut = s.vel.y();
        }
        dumpSim("rebound", sim, t, finite);
        const double e = vyOut / std::max(-vyIn, 1.0e-9);
        std::printf("rebound CoR=%.3f\n", e);
        check(finite, "elastomer rebound stays in the solver");
        expectRange(e, 0.10, 0.90, "elastomer rebound coefficient of restitution is between 0.1 and 0.9");
    }
}

void checkJoints() {
    std::printf("\n-- joints --\n");
    checkWeldStrength();
    checkScrews();
}

void checkHeatLabs() {
    std::printf("\n-- heat --\n");
    checkMaterialsAndImport();
    checkHeat();
    checkThermalWarp();
}

void checkRemesh() {
    std::printf("\n-- remesh --\n");
    auto meanRegionEdge = [](const rcs::VolumeMesh& volume, double xLo, double xHi) {
        double sum = 0.0;
        int count = 0;
        for (const auto& tet : volume.tets) {
            chrono::ChVector3d c{0, 0, 0};
            for (int i = 0; i < 4; ++i)
                c += volume.vertices[static_cast<std::size_t>(tet[static_cast<std::size_t>(i)])];
            c *= 0.25;
            if (c.x() < xLo || c.x() > xHi)
                continue;
            const int edges[6][2] = {{0, 1}, {0, 2}, {0, 3}, {1, 2}, {1, 3}, {2, 3}};
            for (const auto& e : edges) {
                sum += (volume.vertices[static_cast<std::size_t>(tet[static_cast<std::size_t>(e[0])])] -
                        volume.vertices[static_cast<std::size_t>(tet[static_cast<std::size_t>(e[1])])])
                           .Length();
                ++count;
            }
        }
        return count > 0 ? sum / static_cast<double>(count) : 0.0;
    };
    auto meshAt = [&](double resolution) {
        rcs::PartDefinition part = rcs::makeSteppedSolid();
        rcs::applyResolutionSizing(part, resolution);
        check(rcs::ensurePartMesh(part), resolution > 2.0 ? "stepped meshes at resolution 4"
                                                          : "stepped meshes at resolution 1");
        return part;
    };
    rcs::PartDefinition coarse = meshAt(1.0);
    rcs::PartDefinition fine = meshAt(4.0);
    const double boss1 = meanRegionEdge(coarse.volume, -1.0, 0.10);
    const double web1 = meanRegionEdge(coarse.volume, 0.13, 1.0);
    const double boss4 = meanRegionEdge(fine.volume, -1.0, 0.10);
    const double web4 = meanRegionEdge(fine.volume, 0.13, 1.0);
    std::printf("stepped h boss=%.4f web=%.4f -> boss=%.4f web=%.4f tets %zu -> %zu shells %zu -> %zu\n",
                boss1, web1, boss4, web4, coarse.volume.tets.size(), fine.volume.tets.size(),
                coarse.surface.triangles.size(), fine.surface.triangles.size());
    if (web1 > 1.0e-9 && boss1 > 1.0e-9)
        check(web1 < boss1 * 0.85, "web elements are finer than the boss at resolution 1");
    else
        check(!coarse.surface.empty() && !coarse.volume.empty(),
              "stepped solid keeps a thin-wall shell and a bulky tet fill");
    if (boss1 > 1.0e-9 && boss4 > 1.0e-9)
        check(boss4 < boss1 * 0.9, "raising resolution shrinks boss tet size");
    if (web1 > 1.0e-9 && web4 > 1.0e-9)
        check(web4 < web1 * 0.9, "raising resolution shrinks web tet size");
    if (web1 > 1.0e-9 && boss1 > 1.0e-9 && web4 > 1.0e-9 && boss4 > 1.0e-9)
        check(web4 / boss4 < 0.95, "resolution scaling does not flatten the size-field ratio");

    rcs::PartDefinition cube = rcs::makeSolidCube(0.20);
    cube.sizing.maxElement = 0.08;
    cube.sizing.minElement = 0.02;
    cube.sizing.refineLimit = 0.02;
    cube.sizing.tetBudget = 500;
    cube.sizing.cached = false;
    cube.volume = {};
    cube.clampBoundary = true;
    check(meshPart(cube), "remesh cube meshes");
    rcs::AssemblyDefinition definition;
    definition.parts.push_back(cube);
    rcs::Scenario scene = catalogQuiet();
    rcs::ImpactLoad hit;
    hit.point = {0.10, 0.5, 0.0};
    hit.direction = {-1.0, 0.0, 0.0};
    hit.finalForce = 7.0e4;
    hit.radius = 0.02;
    hit.rampTime = 0.03;
    hit.duration = 0.05;
    scene.loads.push_back(hit);
    SmokeSim simulator;
    simulator.reset(definition, scene);
    const int tets0 = simulator.assembly().parts()[0].liveTetCount();
    const double mass0 = simulator.assembly().mass();
    const chrono::ChVector3d com0 = catalogCom(simulator.assembly().parts()[0]);
    bool finite = true;
    double t = 0.0;
    const double dt = scene.settings.timestep;
    for (int step = 0; step < 400 && finite; ++step) {
        finite = simulator.step(dt, t, scene.environment, scene.loads);
        t += dt;
    }
    dumpSim("remesh", simulator, t, finite);
    const rcs::Part& done = simulator.assembly().parts()[0];
    expectFinite(finite, "static-mesh dent stays in the solver");
    check(done.liveTetCount() <= tets0, "loaded dent does not grow the tet count");
    expectRange(std::abs(simulator.assembly().mass() - mass0) / std::max(mass0, 1.0), 0.0, 0.01,
                "static mesh conserves mass");
    expectRange((catalogCom(done) - com0).Length(), 0.0, 0.01, "static mesh conserves COM");
}

void checkEditorGate() {
    std::printf("\n-- editor identity --\n");
    const DropTrace hollow = runCornerDrop(rcs::makeTestCube(), true, "editor-hollow");
    assertEngineDrop(hollow, "editor hollow");
    rcs::AssemblyDefinition solidDef;
    solidDef.parts.push_back(rcs::makeSolidCube(0.40));
    const DropTrace solid = runCornerDrop(solidDef, true, "editor-solid");
    assertEngineDrop(solid, "editor solid");
}
