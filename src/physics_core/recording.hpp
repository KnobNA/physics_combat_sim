#pragma once

#include <array>
#include <cstddef>
#include <cstdint>
#include <string>
#include <vector>

namespace rcs {

struct Vec3f {
    float x = 0.0f;
    float y = 0.0f;
    float z = 0.0f;
};

/// Static topology shared by every frame, so frames only carry state.
struct PartLayout {
    std::string name;
    std::string materialName;
    float thickness = 0.0f;
    int nodeCount = 0;
    std::vector<std::array<int, 4>> quads;
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
    std::vector<Vec3f> positions;     ///< one per node
    std::vector<float> stress;        ///< von Mises per quad, Pa
};

struct JointFrame {
    /// Per-link utilisation: reaction magnitude / allowable, 1.0 == at limit.
    std::vector<float> utilisation;
    std::vector<std::uint8_t> broken;
};

struct Frame {
    double time = 0.0;
    std::vector<PartFrame> parts;
    std::vector<JointFrame> joints;
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
