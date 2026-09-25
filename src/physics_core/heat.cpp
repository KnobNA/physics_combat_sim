#include "physics_core/heat.hpp"

#include <algorithm>
#include <cmath>
#include <unordered_map>

namespace rcs {
namespace {

constexpr double kHeatInterval = 0.02;
constexpr double kAirInterval = 0.05;
constexpr double kWake = 0.05;
constexpr double kAirRho = 1.20;
constexpr double kAirCp = 1006.0;
constexpr double kAirK = 0.026;
constexpr double kConvection = 8.0;
constexpr double kSourceTau = 0.08;
constexpr double kMaxSourceStep = 80.0;

bool sourceActive(const HeatSource& source, double time) {
    return source.enabled && source.duration > 0.0 && time + 1.0e-12 >= source.startTime &&
           time < source.startTime + source.duration;
}

double relaxToward(double current, double target, double hdt) {
    const double gain = 1.0 - std::exp(-hdt / kSourceTau);
    double delta = (target - current) * gain;
    if (delta > kMaxSourceStep)
        delta = kMaxSourceStep;
    else if (delta < -kMaxSourceStep)
        delta = -kMaxSourceStep;
    return current + delta;
}

/// Closest point on a triangle. A sphere resting on a face often misses every
/// corner node while still covering the shell between them.
chrono::ChVector3d closestOnTriangle(const chrono::ChVector3d& p,
                                    const chrono::ChVector3d& a,
                                    const chrono::ChVector3d& b,
                                    const chrono::ChVector3d& c)
{
    const chrono::ChVector3d ab = b - a;
    const chrono::ChVector3d ac = c - a;
    const chrono::ChVector3d ap = p - a;
    const double d1 = ab.Dot(ap);
    const double d2 = ac.Dot(ap);
    if (d1 <= 0.0 && d2 <= 0.0)
        return a;
    const chrono::ChVector3d bp = p - b;
    const double d3 = ab.Dot(bp);
    const double d4 = ac.Dot(bp);
    if (d3 >= 0.0 && d4 <= d3)
        return b;
    const double vc = d1 * d4 - d3 * d2;
    if (vc <= 0.0 && d1 >= 0.0 && d3 <= 0.0) {
        const double v = d1 / (d1 - d3);
        return a + v * ab;
    }
    const chrono::ChVector3d cp = p - c;
    const double d5 = ab.Dot(cp);
    const double d6 = ac.Dot(cp);
    if (d6 >= 0.0 && d5 <= d6)
        return c;
    const double vb = d5 * d2 - d1 * d6;
    if (vb <= 0.0 && d2 >= 0.0 && d6 <= 0.0) {
        const double w = d2 / (d2 - d6);
        return a + w * ac;
    }
    const double va = d3 * d6 - d5 * d4;
    if (va <= 0.0 && (d4 - d3) >= 0.0 && (d5 - d6) >= 0.0) {
        const double w = (d4 - d3) / ((d4 - d3) + (d5 - d6));
        return b + w * (c - b);
    }
    const double denom = 1.0 / (va + vb + vc);
    return a + ab * (vb * denom) + ac * (vc * denom);
}

}  // namespace

void HeatField::configure(const Environment& environment) {
    const double spanX = std::max(1.0, environment.boundsSize.x());
    const double spanY = std::max(1.0, environment.boundsSize.y());
    const double spanZ = std::max(1.0, environment.boundsSize.z());
    const double floor = environment.floorHeight;
    const double ambient = environment.ambientTemperature;
    if (_nx > 0 && spanX == _spanX && spanY == _spanY && spanZ == _spanZ && floor == _floor &&
        std::abs(ambient - _ambient) < 1.0e-6)
        return;

    // A torch is tens of millimetres across. At quarter-metre cells its heat
    // lands in one cell that also contains most of the part, and the whole
    // part soaks instead of developing a gradient. Only active cells are
    // stepped, so a finer grid costs memory rather than time.
    auto count = [](double span) {
        const int n = static_cast<int>(std::ceil(span / 0.05));
        return std::clamp(n, 8, 48);
    };
    _spanX = spanX;
    _spanY = spanY;
    _spanZ = spanZ;
    _floor = floor;
    _ambient = ambient;
    _nx = count(spanX);
    _ny = count(spanY);
    _nz = count(spanZ);
    _dx = spanX / _nx;
    _dy = spanY / _ny;
    _dz = spanZ / _nz;
    _originX = -0.5 * spanX;
    _originY = floor;
    _originZ = -0.5 * spanZ;
    _air.assign(static_cast<std::size_t>(_nx * _ny * _nz), static_cast<float>(ambient));
    _airActive.assign(_air.size(), 0);
    _airList.clear();
    _dirty = false;
}

int HeatField::cellAt(const chrono::ChVector3d& point) const {
    if (_nx <= 0)
        return -1;
    const int ix = static_cast<int>(std::floor((point.x() - _originX) / _dx));
    const int iy = static_cast<int>(std::floor((point.y() - _originY) / _dy));
    const int iz = static_cast<int>(std::floor((point.z() - _originZ) / _dz));
    if (ix < 0 || iy < 0 || iz < 0 || ix >= _nx || iy >= _ny || iz >= _nz)
        return -1;
    return ix + _nx * (iy + _ny * iz);
}

chrono::ChVector3d HeatField::cellCenter(int index) const {
    const int ix = index % _nx;
    const int iy = (index / _nx) % _ny;
    const int iz = index / (_nx * _ny);
    return {_originX + (ix + 0.5) * _dx, _originY + (iy + 0.5) * _dy, _originZ + (iz + 0.5) * _dz};
}

void HeatField::rebuildCache(Cache& cache, const HeatSource& source, const std::vector<Part>& parts, int stamp) {
    cache.part.clear();
    cache.node.clear();
    cache.cell.clear();
    cache.nodeStamp = stamp;
    const double radius2 = source.radius * source.radius;
    std::vector<std::vector<char>> seen(parts.size());
    auto addNode = [&](int partIndex, int nodeIndex) {
        if (partIndex < 0 || nodeIndex < 0)
            return;
        auto& marks = seen[static_cast<std::size_t>(partIndex)];
        const Part& part = parts[static_cast<std::size_t>(partIndex)];
        if (static_cast<std::size_t>(nodeIndex) >= part.nodeCount())
            return;
        if (marks.size() < part.nodeCount())
            marks.assign(part.nodeCount(), 0);
        if (marks[static_cast<std::size_t>(nodeIndex)])
            return;
        marks[static_cast<std::size_t>(nodeIndex)] = 1;
        cache.part.push_back(partIndex);
        cache.node.push_back(nodeIndex);
    };
    for (int p = 0; p < static_cast<int>(parts.size()); ++p) {
        const Part& part = parts[static_cast<std::size_t>(p)];
        // The drawn plate stands off the mid-surface by half the thickness, so
        // a sphere that is touching that surface still covers the shell.
        const double reach = source.radius + 0.5 * part.thickness() + 0.002;
        const double reach2 = reach * reach;
        for (int n = 0; n < static_cast<int>(part.nodeCount()); ++n) {
            if ((part.nodePosition(static_cast<std::size_t>(n)) - source.point).Length2() <= reach2)
                addNode(p, n);
        }
        auto faceHit = [&](int i0, int i1, int i2) {
            if (i0 < 0 || i1 < 0 || i2 < 0)
                return false;
            const chrono::ChVector3d a = part.nodePosition(static_cast<std::size_t>(i0));
            const chrono::ChVector3d b = part.nodePosition(static_cast<std::size_t>(i1));
            const chrono::ChVector3d c = part.nodePosition(static_cast<std::size_t>(i2));
            return (closestOnTriangle(source.point, a, b, c) - source.point).Length2() <= reach2;
        };
        for (const auto& quad : part.quads()) {
            if (!faceHit(quad[0], quad[1], quad[2]) && !faceHit(quad[0], quad[2], quad[3]))
                continue;
            for (int node : quad)
                addNode(p, node);
        }
        if (!part.quads().empty())
            continue;
        for (const auto& tri : part.triangles()) {
            if (!faceHit(tri[0], tri[1], tri[2]))
                continue;
            for (int node : tri)
                addNode(p, node);
        }
    }
    const double hx = 0.5 * _dx;
    const double hy = 0.5 * _dy;
    const double hz = 0.5 * _dz;
    for (int c = 0; c < static_cast<int>(_air.size()); ++c) {
        const chrono::ChVector3d center = cellCenter(c);
        const double cx = std::clamp(source.point.x(), center.x() - hx, center.x() + hx);
        const double cy = std::clamp(source.point.y(), center.y() - hy, center.y() + hy);
        const double cz = std::clamp(source.point.z(), center.z() - hz, center.z() + hz);
        const double dx = cx - source.point.x();
        const double dy = cy - source.point.y();
        const double dz = cz - source.point.z();
        if (dx * dx + dy * dy + dz * dz <= radius2)
            cache.cell.push_back(c);
    }
}

void HeatField::exchangeEdge(Part& part, int a, int b, double hdt, bool volume) {
    if (a == b || a < 0 || b < 0)
        return;
    const std::size_t ia = static_cast<std::size_t>(a);
    const std::size_t ib = static_cast<std::size_t>(b);
    if (ia >= part.nodeCount() || ib >= part.nodeCount())
        return;
    if (!part.heatActive(ia) && !part.heatActive(ib))
        return;
    const double length = (part.nodePosition(ia) - part.nodePosition(ib)).Length();
    if (length < 1.0e-6)
        return;
    const double massA = ia < part.nodeMass().size() ? part.nodeMass()[ia] : 0.0;
    const double massB = ib < part.nodeMass().size() ? part.nodeMass()[ib] : 0.0;
    const double cp = std::max(part.material().specificHeat, 1.0);
    const double cA = std::max(massA, 1.0e-8) * cp;
    const double cB = std::max(massB, 1.0e-8) * cp;
    const double density = std::max(part.material().density, 1.0);
    const double section = volume ? 0.5 * (massA + massB) / density / std::max(length, 1.0e-6)
                                  : std::max(part.thickness(), 1.0e-5);
    const double conductance = part.material().conductivity * section / std::max(length, 1.0e-6);
    const double ti = part.nodeTemperature(ia);
    const double tj = part.nodeTemperature(ib);
    double flow = conductance * hdt * (tj - ti);
    const double cap = (volume ? 0.02 : 0.25) * std::min(cA, cB) * std::abs(tj - ti);
    if (std::abs(flow) > cap)
        flow = std::copysign(cap, flow);
    part.setNodeTemperature(ia, ti + flow / cA);
    part.setNodeTemperature(ib, tj - flow / cB);
}

void HeatField::exchangePair(Part& a, std::size_t ia, Part& b, std::size_t ib, double conductance, double hdt) {
    if (conductance <= 0.0 || ia >= a.nodeCount() || ib >= b.nodeCount())
        return;
    const double cpA = std::max(a.material().specificHeat, 1.0);
    const double cpB = std::max(b.material().specificHeat, 1.0);
    const double massA = ia < a.nodeMass().size() ? a.nodeMass()[ia] : 0.0;
    const double massB = ib < b.nodeMass().size() ? b.nodeMass()[ib] : 0.0;
    const double cA = std::max(massA, 1.0e-8) * cpA;
    const double cB = std::max(massB, 1.0e-8) * cpB;
    const double ti = a.nodeTemperature(ia);
    const double tj = b.nodeTemperature(ib);
    if (!a.heatActive(ia) && !b.heatActive(ib) && std::abs(ti - tj) <= kWake)
        return;
    double flow = conductance * hdt * (tj - ti);
    const double cap = 0.25 * std::min(cA, cB) * std::abs(tj - ti);
    if (std::abs(flow) > cap)
        flow = std::copysign(cap, flow);
    if (std::abs(flow) < 1.0e-12)
        return;
    a.setNodeTemperature(ia, ti + flow / cA);
    b.setNodeTemperature(ib, tj - flow / cB);
    if (std::abs(a.nodeTemperature(ia) - a.initialTemperature()) > kWake)
        a.setHeatActive(ia, true);
    if (std::abs(b.nodeTemperature(ib) - b.initialTemperature()) > kWake)
        b.setHeatActive(ib, true);
    _dirty = true;
}

void HeatField::exchangeRadiation(std::vector<Part>& parts, const HeatSource& source, double hdt) {
    // A few-thousand-kelvin body dumps far more heat by radiation than still
    // air can carry by natural convection. Facing shell within a few radii
    // warms; the back of a panel does not.
    constexpr double kSigma = 5.670374419e-8;
    const double radius = std::max(source.radius, 1.0e-4);
    const double inside2 = radius * radius;
    const double far = radius * 8.0;
    const double far2 = far * far;
    const double sourceT4 = std::pow(std::max(source.temperature, 0.0), 4.0);

    for (Part& part : parts) {
        const double emissivity = std::clamp(part.material().emissivity, 0.0, 1.0);
        if (emissivity <= 0.0 || part.nodeCount() == 0)
            continue;
        const std::size_t count = part.nodeCount();
        std::vector<chrono::ChVector3d> normal(count, chrono::ChVector3d(0.0, 0.0, 0.0));
        auto addFace = [&](int i0, int i1, int i2) {
            if (i0 < 0 || i1 < 0 || i2 < 0)
                return;
            const chrono::ChVector3d a = part.nodePosition(static_cast<std::size_t>(i0));
            const chrono::ChVector3d b = part.nodePosition(static_cast<std::size_t>(i1));
            const chrono::ChVector3d c = part.nodePosition(static_cast<std::size_t>(i2));
            chrono::ChVector3d n = (b - a).Cross(c - a);
            if (n.Length2() < 1.0e-16)
                return;
            n.Normalize();
            normal[static_cast<std::size_t>(i0)] += n;
            normal[static_cast<std::size_t>(i1)] += n;
            normal[static_cast<std::size_t>(i2)] += n;
        };
        if (!part.quads().empty()) {
            for (const auto& quad : part.quads()) {
                addFace(quad[0], quad[1], quad[2]);
                addFace(quad[0], quad[2], quad[3]);
            }
        } else {
            for (const auto& tri : part.triangles())
                addFace(tri[0], tri[1], tri[2]);
        }

        const double density = std::max(part.material().density, 1.0);
        const double thickness = std::max(part.thickness(), 1.0e-5);
        const double cp = std::max(part.material().specificHeat, 1.0);
        for (std::size_t n = 0; n < count; ++n) {
            if (normal[n].Length2() < 1.0e-12)
                continue;
            normal[n].Normalize();
            const chrono::ChVector3d offset = source.point - part.nodePosition(n);
            const double dist2 = offset.Length2();
            if (dist2 <= inside2 || dist2 > far2)
                continue;
            const double dist = std::sqrt(dist2);
            const double facing = normal[n].Dot(offset) / dist;
            if (facing <= 0.0)
                continue;
            const double view = (radius / dist) * (radius / dist) * facing;
            const double mass = n < part.nodeMass().size() ? part.nodeMass()[n] : 0.0;
            if (mass <= 0.0)
                continue;
            const double area = mass / (density * thickness);
            const double nodeT = part.nodeTemperature(n);
            const double flux = emissivity * kSigma * view * (sourceT4 - std::pow(std::max(nodeT, 0.0), 4.0));
            double delta = flux * area * hdt / (mass * cp);
            if (delta > kMaxSourceStep)
                delta = kMaxSourceStep;
            else if (delta < -kMaxSourceStep)
                delta = -kMaxSourceStep;
            if (std::abs(delta) < 1.0e-4)
                continue;
            part.setNodeTemperature(n, nodeT + delta);
            if (std::abs(part.nodeTemperature(n) - part.initialTemperature()) > kWake)
                part.setHeatActive(n, true);
            _dirty = true;
        }
    }
}

void HeatField::exchangeContact(std::vector<Part>& parts, double hdt) {
    struct Ref {
        int part = 0;
        int node = 0;
    };
    std::unordered_map<int, std::vector<Ref>> buckets;
    for (int p = 0; p < static_cast<int>(parts.size()); ++p) {
        Part& part = parts[static_cast<std::size_t>(p)];
        for (int n = 0; n < static_cast<int>(part.nodeCount()); ++n) {
            const std::size_t index = static_cast<std::size_t>(n);
            const int cell = cellAt(part.nodePosition(index));
            if (cell < 0)
                continue;
            const bool cellHot =
                std::abs(_air[static_cast<std::size_t>(cell)] - static_cast<float>(_ambient)) > kWake;
            if (!part.heatActive(index) && !cellHot)
                continue;
            buckets[cell].push_back({p, n});
        }
    }

    for (const auto& bucket : buckets) {
        const auto& refs = bucket.second;
        for (std::size_t i = 0; i < refs.size(); ++i) {
            for (std::size_t j = i + 1; j < refs.size(); ++j) {
                if (refs[i].part == refs[j].part)
                    continue;
                Part& a = parts[static_cast<std::size_t>(refs[i].part)];
                Part& b = parts[static_cast<std::size_t>(refs[j].part)];
                const std::size_t ia = static_cast<std::size_t>(refs[i].node);
                const std::size_t ib = static_cast<std::size_t>(refs[j].node);
                const double distance = (a.nodePosition(ia) - b.nodePosition(ib)).Length();
                const double gap = 0.5 * (a.thickness() + b.thickness()) + 0.004;
                // Coincident weld nodes are tied separately. A real touch is
                // the two mid-surfaces a plate-thickness apart.
                if (distance < 1.0e-3 || distance > gap)
                    continue;
                const double densityA = std::max(a.material().density, 1.0);
                const double densityB = std::max(b.material().density, 1.0);
                const double massA = ia < a.nodeMass().size() ? a.nodeMass()[ia] : 0.0;
                const double massB = ib < b.nodeMass().size() ? b.nodeMass()[ib] : 0.0;
                const double areaA = massA / (densityA * std::max(a.thickness(), 1.0e-5));
                const double areaB = massB / (densityB * std::max(b.thickness(), 1.0e-5));
                const double area = std::min(areaA, areaB);
                // Pressed metal-to-metal contact, W/(m^2*K), well below a weld.
                exchangePair(a, ia, b, ib, 4000.0 * area, hdt);
            }
        }
    }
}

void HeatField::wake(std::vector<Part>& parts) {
    _dirty = false;
    for (Part& part : parts) {
        const std::size_t count = part.nodeCount();
        std::vector<char> next(count, 0);
        auto mark = [&](int index) {
            if (index >= 0 && static_cast<std::size_t>(index) < count)
                next[static_cast<std::size_t>(index)] = 1;
        };
        auto consider = [&](int a, int b) {
            if (a < 0 || b < 0)
                return;
            const std::size_t ia = static_cast<std::size_t>(a);
            const std::size_t ib = static_cast<std::size_t>(b);
            if (ia >= count || ib >= count)
                return;
            const double ta = part.nodeTemperature(ia);
            const double tb = part.nodeTemperature(ib);
            const bool drifted = std::abs(ta - part.initialTemperature()) > kWake ||
                                 std::abs(tb - part.initialTemperature()) > kWake ||
                                 std::abs(ta - tb) > kWake;
            if (!drifted)
                return;
            if (part.heatActive(ia) || part.heatActive(ib)) {
                mark(a);
                mark(b);
            }
        };
        for (const auto& quad : part.quads()) {
            consider(quad[0], quad[1]);
            consider(quad[1], quad[2]);
            consider(quad[2], quad[3]);
            consider(quad[3], quad[0]);
        }
        for (const auto& tri : part.triangles()) {
            consider(tri[0], tri[1]);
            consider(tri[1], tri[2]);
            consider(tri[2], tri[0]);
        }
        for (std::size_t t = 0; t < part.tets().size(); ++t) {
            if (!part.tetLive(t))
                continue;
            const auto& tet = part.tets()[t];
            for (int i = 0; i < 4; ++i)
                for (int j = i + 1; j < 4; ++j)
                    consider(tet[i], tet[j]);
        }
        for (std::size_t n = 0; n < count; ++n) {
            part.setHeatActive(n, next[n] != 0);
            if (next[n])
                _dirty = true;
        }
    }

    if (_air.empty())
        return;
    std::vector<char> next = _airActive;
    auto wakeCell = [&](int index) {
        if (index < 0 || static_cast<std::size_t>(index) >= next.size())
            return;
        if (std::abs(_air[static_cast<std::size_t>(index)] - _ambient) <= kWake && !_airActive[static_cast<std::size_t>(index)])
            return;
        next[static_cast<std::size_t>(index)] = 1;
    };
    for (int index : _airList) {
        if (index < 0 || static_cast<std::size_t>(index) >= _air.size())
            continue;
        const int ix = index % _nx;
        const int iy = (index / _nx) % _ny;
        const int iz = index / (_nx * _ny);
        const float here = _air[static_cast<std::size_t>(index)];
        auto neighbour = [&](int x, int y, int z) {
            if (x < 0 || y < 0 || z < 0 || x >= _nx || y >= _ny || z >= _nz)
                return;
            const int other = x + _nx * (y + _ny * z);
            if (std::abs(here - _air[static_cast<std::size_t>(other)]) > kWake ||
                std::abs(here - static_cast<float>(_ambient)) > kWake)
                wakeCell(other);
        };
        neighbour(ix - 1, iy, iz);
        neighbour(ix + 1, iy, iz);
        neighbour(ix, iy - 1, iz);
        neighbour(ix, iy + 1, iz);
        neighbour(ix, iy, iz - 1);
        neighbour(ix, iy, iz + 1);
        if (std::abs(here - static_cast<float>(_ambient)) <= kWake)
            next[static_cast<std::size_t>(index)] = 0;
        else
            _dirty = true;
    }
    _airActive.swap(next);
    _airList.clear();
    for (int i = 0; i < static_cast<int>(_airActive.size()); ++i) {
        if (_airActive[static_cast<std::size_t>(i)])
            _airList.push_back(i);
    }
}

bool HeatField::easeStrain(std::vector<Part>& parts, double dt) {
    const double cap = 2000.0 * std::max(dt, 0.0);
    bool moved = false;
    for (Part& part : parts) {
        bool partMoved = false;
        for (std::size_t n = 0; n < part.nodeCount(); ++n) {
            const double target = part.nodeTemperature(n);
            const double applied = part.strainTemperature(n);
            const double delta = target - applied;
            if (std::abs(delta) <= 0.02)
                continue;
            const double step = std::clamp(delta, -cap, cap);
            part.setStrainTemperature(n, applied + step);
            partMoved = true;
        }
        // Node temperature can already be at the colour the user sees while
        // strain is still catching up. Scale has to follow the metal now.
        part.syncThermalMechanics();
        if (partMoved)
            moved = true;
    }
    return moved;
}

bool HeatField::advance(double dt,
                        double time,
                        std::vector<Part>& parts,
                        const Environment& environment,
                        const std::vector<HeatSource>& sources,
                        const std::vector<HeatBridge>& bridges)
{
    if (dt <= 0.0)
        return false;
    configure(environment);

    bool sourceOn = false;
    for (const HeatSource& source : sources) {
        if (sourceActive(source, time)) {
            sourceOn = true;
            break;
        }
    }
    if (!sourceOn && !_dirty && !_easing) {
        _wait = 0.0;
        _sourceOn = false;
        return false;
    }

    _wait += dt;
    const bool edge = sourceOn != _sourceOn;
    if (_wait < kHeatInterval && !edge) {
        _easing = easeStrain(parts, dt);
        return _easing;
    }
    const double hdt = std::min(std::max(_wait, dt), 0.25);
    _wait = 0.0;
    _sourceOn = sourceOn;

    int stamp = 0;
    for (const Part& part : parts)
        stamp += static_cast<int>(part.nodeCount());
    if (_caches.size() != sources.size())
        _caches.assign(sources.size(), {});

    for (std::size_t s = 0; s < sources.size(); ++s) {
        Cache& cache = _caches[s];
        const bool on = sourceActive(sources[s], time);
        if (!on) {
            cache.on = false;
            continue;
        }
        // Positions move, and a sphere can cover a face without covering a
        // corner, so the overlap list is rebuilt on each heat step.
        rebuildCache(cache, sources[s], parts, stamp);
        cache.on = true;
        for (std::size_t i = 0; i < cache.node.size(); ++i) {
            Part& part = parts[static_cast<std::size_t>(cache.part[i])];
            const std::size_t node = static_cast<std::size_t>(cache.node[i]);
            part.setNodeTemperature(node, relaxToward(part.nodeTemperature(node), sources[s].temperature, hdt));
            part.setHeatActive(node, true);
            _dirty = true;
        }
        for (int cell : cache.cell) {
            if (cell < 0 || static_cast<std::size_t>(cell) >= _air.size())
                continue;
            const double updated = relaxToward(_air[static_cast<std::size_t>(cell)], sources[s].temperature, hdt);
            _air[static_cast<std::size_t>(cell)] = static_cast<float>(updated);
            if (!_airActive[static_cast<std::size_t>(cell)]) {
                _airActive[static_cast<std::size_t>(cell)] = 1;
                _airList.push_back(cell);
            }
            _dirty = true;
        }
        exchangeRadiation(parts, sources[s], hdt);
    }

    for (const HeatBridge& bridge : bridges) {
        if (bridge.partA >= parts.size() || bridge.partB >= parts.size())
            continue;
        Part& a = parts[bridge.partA];
        Part& b = parts[bridge.partB];
        if (bridge.nodeA >= a.nodeCount() || bridge.nodeB >= b.nodeCount())
            continue;
        // A weld or screw is the same metal at one point, so the two nodes
        // share heat the way a continuous plate would. The cap keeps the
        // exchange stable; the conductance itself is not the limit.
        exchangePair(a, bridge.nodeA, b, bridge.nodeB, 1.0e9, hdt);
    }
    exchangeContact(parts, hdt);

    for (Part& part : parts) {
        for (const auto& quad : part.quads()) {
            const bool live = quad[0] >= 0 && part.heatActive(static_cast<std::size_t>(quad[0])) ||
                              quad[1] >= 0 && part.heatActive(static_cast<std::size_t>(quad[1])) ||
                              quad[2] >= 0 && part.heatActive(static_cast<std::size_t>(quad[2])) ||
                              quad[3] >= 0 && part.heatActive(static_cast<std::size_t>(quad[3]));
            if (!live)
                continue;
            exchangeEdge(part, quad[0], quad[1], hdt, false);
            exchangeEdge(part, quad[1], quad[2], hdt, false);
            exchangeEdge(part, quad[2], quad[3], hdt, false);
            exchangeEdge(part, quad[3], quad[0], hdt, false);
        }
        if (part.hasVolume()) {
            for (std::size_t t = 0; t < part.tets().size(); ++t) {
                if (!part.tetLive(t))
                    continue;
                const auto& tet = part.tets()[t];
                for (int i = 0; i < 4; ++i)
                    for (int j = i + 1; j < 4; ++j)
                        exchangeEdge(part, tet[i], tet[j], hdt, true);
            }
            const auto& tris = part.triangles();
            for (std::size_t t = static_cast<std::size_t>(part.skinTriangleCount()); t < tris.size(); ++t) {
                exchangeEdge(part, tris[t][0], tris[t][1], hdt, false);
                exchangeEdge(part, tris[t][1], tris[t][2], hdt, false);
                exchangeEdge(part, tris[t][2], tris[t][0], hdt, false);
            }
            continue;
        }
        if (!part.quads().empty())
            continue;
        for (const auto& tri : part.triangles()) {
            const bool live = tri[0] >= 0 && part.heatActive(static_cast<std::size_t>(tri[0])) ||
                              tri[1] >= 0 && part.heatActive(static_cast<std::size_t>(tri[1])) ||
                              tri[2] >= 0 && part.heatActive(static_cast<std::size_t>(tri[2]));
            if (!live)
                continue;
            exchangeEdge(part, tri[0], tri[1], hdt);
            exchangeEdge(part, tri[1], tri[2], hdt);
            exchangeEdge(part, tri[2], tri[0], hdt);
        }
    }

    const double cellVolume = _dx * _dy * _dz;
    const double airCapacity = kAirRho * kAirCp * cellVolume;
    for (Part& part : parts) {
        const double cp = std::max(part.material().specificHeat, 1.0);
        // Air touches the outside of a part only. A node buried inside a
        // solid has no wetted area and must not exchange with the room.
        // A shell is a midsurface, so both of its faces are wetted.
        const auto& wetted = part.contactAreas();
        const double faces = part.hasVolume() ? 1.0 : 2.0;
        for (std::size_t n = 0; n < part.nodeCount(); ++n) {
            const int cell = cellAt(part.nodePosition(n));
            if (cell < 0)
                continue;
            const double area = faces * (n < wetted.size() ? wetted[n] : 0.0);
            if (area <= 0.0)
                continue;
            const double airT = _air[static_cast<std::size_t>(cell)];
            const double nodeT = part.nodeTemperature(n);
            const bool cellHot = std::abs(airT - _ambient) > kWake;
            const bool nodeHot = part.heatActive(n) || std::abs(nodeT - part.initialTemperature()) > kWake;
            if (!cellHot && !nodeHot)
                continue;
            const double mass = n < part.nodeMass().size() ? part.nodeMass()[n] : 0.0;
            if (mass <= 0.0)
                continue;
            const double nodeCapacity = mass * cp;
            const double conductance = kConvection * area;
            double flow = conductance * hdt * (airT - nodeT);
            const double cap = 0.25 * std::min(nodeCapacity, airCapacity) * std::abs(airT - nodeT);
            if (std::abs(flow) > cap)
                flow = std::copysign(cap, flow);
            part.setNodeTemperature(n, nodeT + flow / nodeCapacity);
            if (std::abs(part.nodeTemperature(n) - part.initialTemperature()) > kWake)
                part.setHeatActive(n, true);
            _air[static_cast<std::size_t>(cell)] = static_cast<float>(airT - flow / airCapacity);
            if (!_airActive[static_cast<std::size_t>(cell)]) {
                _airActive[static_cast<std::size_t>(cell)] = 1;
                _airList.push_back(cell);
            }
        }
    }

    _airWait += hdt;
    if (_airWait >= kAirInterval || edge) {
        const double airDt = std::min(std::max(_airWait, hdt), 0.25);
        _airWait = 0.0;
        const double alpha = kAirK / (kAirRho * kAirCp);
        auto coeff = [&](double spacing) {
            const double value = alpha * airDt / (spacing * spacing);
            return std::min(value, 0.08);
        };
        const double cx = coeff(_dx);
        const double cy = coeff(_dy);
        const double cz = coeff(_dz);
        std::vector<float> next = _air;
        for (int index : _airList) {
            if (index < 0 || static_cast<std::size_t>(index) >= _air.size())
                continue;
            const int ix = index % _nx;
            const int iy = (index / _nx) % _ny;
            const int iz = index / (_nx * _ny);
            const float here = _air[static_cast<std::size_t>(index)];
            double delta = 0.0;
            auto pull = [&](int x, int y, int z, double c) {
                if (x < 0 || y < 0 || z < 0 || x >= _nx || y >= _ny || z >= _nz)
                    return;
                const int other = x + _nx * (y + _ny * z);
                delta += c * (_air[static_cast<std::size_t>(other)] - here);
            };
            pull(ix - 1, iy, iz, cx);
            pull(ix + 1, iy, iz, cx);
            pull(ix, iy - 1, iz, cy);
            pull(ix, iy + 1, iz, cy);
            pull(ix, iy, iz - 1, cz);
            pull(ix, iy, iz + 1, cz);
            next[static_cast<std::size_t>(index)] = static_cast<float>(here + delta);
        }
        _air.swap(next);
    }

    wake(parts);
    _easing = true;
    easeStrain(parts, dt);
    return true;
}

void HeatField::writeAir(AirFrame& out) const {
    out = {};
    out.ambient = static_cast<float>(_ambient);
    bool off = false;
    for (float value : _air) {
        if (std::abs(value - static_cast<float>(_ambient)) > 0.5f) {
            off = true;
            break;
        }
    }
    if (!off || _air.empty())
        return;
    out.nx = _nx;
    out.ny = _ny;
    out.nz = _nz;
    out.originX = static_cast<float>(_originX);
    out.originY = static_cast<float>(_originY);
    out.originZ = static_cast<float>(_originZ);
    out.dx = static_cast<float>(_dx);
    out.dy = static_cast<float>(_dy);
    out.dz = static_cast<float>(_dz);
    out.temperature = _air;
}

}  // namespace rcs
