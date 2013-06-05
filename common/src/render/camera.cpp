// +-------------------------------------------------------------------------
// | camera.cpp
// | 
// | Author: Gilbert Bernstein
// +-------------------------------------------------------------------------
// | COPYRIGHT:
// |    Copyright Gilbert Bernstein 2012
// |    See the included COPYRIGHT file for further details.
// |    
// |    This file is part of the TopTop library.
// |
// |    TopTop is free software: you can redistribute it and/or modify
// |    it under the terms of the GNU Lesser General Public License as
// |    published by the Free Software Foundation, either version 3 of
// |    the License, or (at your option) any later version.
// |
// |    TopTop is distributed in the hope that it will be useful,
// |    but WITHOUT ANY WARRANTY; without even the implied warranty of
// |    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// |    GNU Lesser General Public License for more details.
// |
// |    You should have received a copy 
// |    of the GNU Lesser General Public License
// |    along with TopTop.  If not, see <http://www.gnu.org/licenses/>.
// +-------------------------------------------------------------------------
#include "camera.h"

// minimum separation distance between pos and lookAt
float Camera::zoomLimit = 1.0e-3;

Camera::Camera() :
    yaw(0), pitch(0), focus(0,0,0), position(0,0,1)
{
    updateTransform();
}

CameraPtr Camera::create()
{
    CameraPtr camera(new Camera);
    return camera;
}
CameraPtr Camera::create(const Vec3f &pos)
{
    CameraPtr camera(new Camera);
    camera->trackTo(pos);
    return camera;
}
CameraPtr Camera::create(const Vec3f &pos, float yaw, float pitch)
{
    CameraPtr camera(new Camera);
    camera->trackTo(pos);
    camera->panLeft(yaw);
    camera->panUp(pitch);
    return camera;
}
CameraPtr Camera::create(const Vec3f &pos, const Vec3f &lookAt)
{
    CameraPtr camera(new Camera);
    camera->trackTo(pos);
    camera->lookAt(lookAt);
    return camera;
}

CameraPtr Camera::clone() const
{
    CameraPtr camera = create(position, focus);
    return camera;
}



float Camera::zoomDist() const
{
    return len(focus-position);
}




void Camera::updateTransform()
{
    Mat3f rot = rotByYawPitch(yaw, pitch);
    
    transform = Transform::translation(position) *
                Transform::rotation(rot, transpose(rot));
}

void Camera::recomputeYawPitch()
{
    Vec3f looking = normalized(focus - position);
    Vec3f xzproj = Vec3f(looking.x, 0, looking.z);
    float horizontal = len(xzproj);
    if(horizontal < 1.0e-6) { // near gimble lock singularity...
        if(looking.y > 0)
            pitch = 90.0;
        else
            pitch = -90.0;
    } else {
        normalize(xzproj);
        // no value wrap necessary b/c atan2 returns val in range (-pi,pi)
        yaw = rad2deg(std::atan2(-xzproj.x, -xzproj.z));
        pitch = rad2deg(std::atan2(looking.y, horizontal));
        // note that 0 yaw and 0 pitch corresponds to direction (0,0,-1)
    }
}
Vec3f lookingFromYawPitch(float yaw, float pitch) {
    yaw = deg2rad(yaw);
    pitch = deg2rad(pitch);
    Vec3f horizontal(-std::sin(yaw), 0, -std::cos(yaw));
    Vec3f vertical(0,1,0);
    return std::cos(pitch)*horizontal + std::sin(pitch)*vertical;
}
void Camera::recomputeFocus()
{
    float dist = len(focus - position);
    Vec3f looking = lookingFromYawPitch(yaw,pitch);
    focus = position + dist * looking;
}
void Camera::recomputePosition()
{
    float dist = len(focus - position);
    Vec3f looking = lookingFromYawPitch(yaw,pitch);
    position = focus - dist * looking;
}

void Camera::lookAt(const Vec3f &point)
{
    if(len(point - position) < zoomLimit)   return;
    focus = point;
    recomputeYawPitch(); // TODO
    updateTransform();
}

void Camera::panLeft(float degrees)
{
    yaw = wrap(yaw + degrees, -180.0f, 180.0f);
    recomputeFocus(); // TODO
    updateTransform();
}
void Camera::panUp(float degrees)
{
    pitch = clamp(pitch + degrees, -90.0f, 90.0f);
    recomputeFocus(); // TODO
    updateTransform();
}

void Camera::orbitLeft(float degrees)
{
    yaw = wrap(yaw - degrees, -180.0f, 180.0f);
    recomputePosition(); // TODO
    updateTransform();
}
void Camera::orbitUp(float degrees)
{
    pitch = clamp(pitch - degrees, -90.0f, 90.0f);
    recomputePosition(); // TODO
    updateTransform();
}

void Camera::dollyLeft(float dist)
{
    Vec3f left(-1,0,0);
    Vec3f offset = dist * transform.rot() * left;
    
    position += offset;
    focus += offset;
    updateTransform();
}
void Camera::dollyUp(float dist)
{
    Vec3f up(0,1,0);
    Vec3f offset = dist * transform.rot() * up;
    
    position += offset;
    focus += offset;
    updateTransform();
}
void Camera::dollyForward(float dist)
{
    Vec3f forward(0,0,-1);
    Vec3f offset = dist * transform.rot() * forward;
    
    position += offset;
    focus += offset;
    updateTransform();
}

void Camera::trackTo(const Vec3f &pos)
{
    Vec3f offset = pos - position;
    
    position = pos;
    focus += offset;
    updateTransform();
}

void Camera::zoomIn(float dist)
{
    Vec3f looking = focus - position;
    if(len(looking) - dist < zoomLimit) return; // do not zoom in too much!
    
    position += dist * normalized(looking);
    updateTransform();
}


