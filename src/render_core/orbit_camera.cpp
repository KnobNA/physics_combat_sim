#include "render_core/orbit_camera.hpp"

#include <Magnum/Math/Constants.h>
#include <Magnum/Math/Functions.h>
#include <Magnum/Math/Angle.h>

namespace rcs {

using namespace Magnum;
using namespace Magnum::Math::Literals;

void OrbitCamera::setViewport(const Vector2i& size) {
    _viewport = size;
    _aspect = size.y() != 0 ? Float(size.x()) / Float(size.y()) : 1.0f;
}

void OrbitCamera::orbit(const Vector2& deltaPixels) {
    _yaw -= deltaPixels.x() * 0.005f;
    _pitch += deltaPixels.y() * 0.005f;
    _pitch = Math::clamp(_pitch, -1.4f, 1.4f);
}

void OrbitCamera::pan(const Vector2& deltaPixels) {
    const Vector3 eye = position();
    const Vector3 forward = (_target - eye).normalized();
    const Vector3 right = Math::cross(forward, Vector3::yAxis()).normalized();
    const Vector3 up = Math::cross(right, forward).normalized();
    const float scale = _distance * 0.0025f;
    _target -= right * deltaPixels.x() * scale;
    _target += up * deltaPixels.y() * scale;
}

void OrbitCamera::zoom(float scroll) {
    _distance *= Math::pow(0.9f, scroll);
    _distance = Math::clamp(_distance, 0.4f, 20.0f);
}

Vector3 OrbitCamera::position() const {
    const float cy = Math::cos(Rad{_yaw});
    const float sy = Math::sin(Rad{_yaw});
    const float cp = Math::cos(Rad{_pitch});
    const float sp = Math::sin(Rad{_pitch});
    return _target + Vector3{sy * cp, sp, cy * cp} * _distance;
}

Matrix4 OrbitCamera::viewMatrix() const {
    return Matrix4::lookAt(position(), _target, Vector3::yAxis()).invertedRigid();
}

Matrix4 OrbitCamera::projectionMatrix() const {
    return Matrix4::perspectiveProjection(35.0_degf, _aspect, 0.05f, 80.0f);
}

void OrbitCamera::rayFromFramebuffer(const Vector2& framebufferPos,
                                     Vector3& origin,
                                     Vector3& direction) const
{
    const Vector2 ndc{
        (2.0f * framebufferPos.x()) / Float(_viewport.x()) - 1.0f,
        1.0f - (2.0f * framebufferPos.y()) / Float(_viewport.y())
    };
    const Matrix4 inv = (projectionMatrix() * viewMatrix()).inverted();
    const Vector4 near = inv * Vector4{ndc.x(), ndc.y(), -1.0f, 1.0f};
    const Vector4 far = inv * Vector4{ndc.x(), ndc.y(), 1.0f, 1.0f};
    origin = near.xyz() / near.w();
    const Vector3 farPoint = far.xyz() / far.w();
    direction = (farPoint - origin).normalized();
}

}  // namespace rcs
