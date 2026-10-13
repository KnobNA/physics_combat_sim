#pragma once

#include <array>
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <cstring>
#include <string>
#include <vector>

namespace rcs {

struct Vec3f {
    float x = 0.0f;
    float y = 0.0f;
    float z = 0.0f;
};

/// IEEE binary16. Playback only needs millimetre-level positions, and a
/// multi-second recording is otherwise dominated by three floats per node.
inline std::uint16_t packHalf(float value) {
    if (!std::isfinite(value))
        return 0;
    std::uint32_t bits = 0;
    static_assert(sizeof(float) == 4, "float is not 32-bit");
    std::memcpy(&bits, &value, sizeof(bits));
    const std::uint32_t sign = (bits >> 16) & 0x8000u;
    const std::int32_t exponent = static_cast<std::int32_t>((bits >> 23) & 0xffu) - 127 + 15;
    const std::uint32_t mantissa = bits & 0x7fffffu;
    if (exponent <= 0)
        return static_cast<std::uint16_t>(sign);
    if (exponent >= 31)
        return static_cast<std::uint16_t>(sign | (30u << 10) | 0x3ffu);
    return static_cast<std::uint16_t>(sign | (static_cast<std::uint32_t>(exponent) << 10) | (mantissa >> 13));
}

inline float unpackHalf(std::uint16_t value) {
    const std::uint32_t sign = static_cast<std::uint32_t>(value & 0x8000u) << 16;
    const std::uint32_t exponent = (value >> 10) & 0x1fu;
    const std::uint32_t mantissa = value & 0x3ffu;
    std::uint32_t bits = 0;
    if (exponent == 0) {
        bits = sign;
    } else if (exponent == 31) {
        bits = sign | 0x7f800000u | (mantissa << 13);
    } else {
        bits = sign | ((exponent + 127 - 15) << 23) | (mantissa << 13);
    }
    float out = 0.0f;
    std::memcpy(&out, &bits, sizeof(out));
    return out;
}

/// Static topology shared by every frame, so frames only carry state.
struct PartLayout {
    std::string name;
    std::string materialName;
    float thickness = 0.0f;
    int nodeCount = 0;
    bool volume = false;
    /// Leading triangles are the solid boundary and are drawn flat.
    int flatTriangles = 0;
    int tetCount = 0;
    std::vector<std::array<int, 4>> quads;
    std::vector<std::array<int, 3>> triangles;
};

struct JointLayout {
    std::string label;
    int linkCount = 0;
};

struct RecordingLayout {
    std::vector<PartLayout> parts;
    std::vector<JointLayout> joints;
};

struct PartFrame {
    /// Packed node positions, three half-floats each.
    std::vector<std::uint16_t> xyz;
    std::vector<float> stress;          ///< von Mises per quad (or triangle), Pa
    std::vector<float> plasticStrain;   ///< accumulated equivalent plastic strain
    std::vector<float> temperature;     ///< mean corner temperature, K
    std::vector<float> nodeTemperature; ///< kelvin at each node, for a point readout
    /// Filled once tearing changes connectivity. Empty means "use the layout".
    std::vector<std::array<int, 4>> quads;
    std::vector<std::array<int, 3>> triangles;

    std::size_t nodeCount() const { return xyz.size() / 3; }

    Vec3f position(std::size_t index) const {
        Vec3f out;
        const std::size_t base = index * 3;
        if (base + 2 >= xyz.size())
            return out;
        out.x = unpackHalf(xyz[base]);
        out.y = unpackHalf(xyz[base + 1]);
        out.z = unpackHalf(xyz[base + 2]);
        return out;
    }

    void setPosition(std::size_t index, double x, double y, double z) {
        const std::size_t base = index * 3;
        if (xyz.size() < base + 3)
            xyz.resize(base + 3);
        xyz[base] = packHalf(static_cast<float>(x));
        xyz[base + 1] = packHalf(static_cast<float>(y));
        xyz[base + 2] = packHalf(static_cast<float>(z));
    }
};

struct JointFrame {
    /// Per-link utilisation: reaction magnitude / allowable, 1.0 == at limit.
    std::vector<float> utilisation;
    std::vector<std::uint8_t> broken;
};

/// Coarse still-air sample. Empty `temperature` means the whole arena is at
/// ambient, so a cold run stores nothing extra.
struct AirFrame {
    int nx = 0;
    int ny = 0;
    int nz = 0;
    float originX = 0.0f;
    float originY = 0.0f;
    float originZ = 0.0f;
    float dx = 0.0f;
    float dy = 0.0f;
    float dz = 0.0f;
    float ambient = 293.15f;
    std::vector<float> temperature;
};

struct Frame {
    double time = 0.0;
    std::vector<PartFrame> parts;
    std::vector<JointFrame> joints;
    AirFrame air;
    float maxStress = 0.0f;
    int brokenLinks = 0;
};

struct Recording {
    RecordingLayout layout;
    std::vector<Frame> frames;
    double duration = 0.0;      ///< simulated time actually reached
    float maxStress = 0.0f;     ///< peak over the whole run
    int brokenLinks = 0;
    bool diverged = false;
    std::string message;

    bool empty() const { return frames.empty(); }

    /// Nearest recorded frame at or before `t`, clamped to the range.
    const Frame* sample(double t) const {
        if (frames.empty())
            return nullptr;
        if (t <= frames.front().time)
            return &frames.front();
        if (t >= frames.back().time)
            return &frames.back();
        std::size_t lo = 0;
        std::size_t hi = frames.size() - 1;
        while (lo + 1 < hi) {
            const std::size_t mid = (lo + hi) / 2;
            if (frames[mid].time <= t)
                lo = mid;
            else
                hi = mid;
        }
        return &frames[lo];
    }
};

}  // namespace rcs
