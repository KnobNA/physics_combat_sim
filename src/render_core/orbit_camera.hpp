#pragma once

#include <Magnum/Magnum.h>
#include <Magnum/Math/Matrix4.h>
#include <Magnum/Math/Range.h>
#include <Magnum/Math/Vector2.h>
#include <Magnum/Math/Vector3.h>

namespace rcs {

class OrbitCamera {
public:
    void setViewport(const Magnum::Vector2i& size);
    void orbit(const Magnum::Vector2& deltaPixels);
    void pan(const Magnum::Vector2& deltaPixels);
    void zoom(float scroll);

    Magnum::Matrix4 viewMatrix() const;
    Magnum::Matrix4 projectionMatrix() const;
    Magnum::Vector3 position() const;
    Magnum::Vector3 target() const { return _target; }
    float distance() const { return _distance; }

    void rayFromFramebuffer(const Magnum::Vector2& framebufferPos,
                            Magnum::Vector3& origin,
                            Magnum::Vector3& direction) const;

private:
    Magnum::Vector3 _target{0.0f, 0.7f, 0.0f};
    float _yaw = -0.7f;
    float _pitch = 0.45f;
    float _distance = 3.6f;
    float _aspect = 16.0f / 9.0f;
    Magnum::Vector2i _viewport{1280, 720};
};

}  // namespace rcs
