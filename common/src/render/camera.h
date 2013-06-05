// +-------------------------------------------------------------------------
// | camera.h
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
#pragma once

#include "scene.h"

class Camera;
typedef std::shared_ptr<Camera> CameraPtr;

class Camera : public SceneNode
{
public:
    Camera();
    virtual ~Camera() {}
    static CameraPtr create();
    static CameraPtr create(const Vec3f &pos);
    static CameraPtr create(const Vec3f &pos, float yaw, float pitch);
    static CameraPtr create(const Vec3f &pos, const Vec3f &lookAt);
    
    CameraPtr clone() const;
    
    // minimum separation distance between pos and lookAt
    static float zoomLimit;
    
public:
    float zoomDist() const; // get distance from camera to focus point
        
public:
    void lookAt(const Vec3f &point);
    // pan -- rotate camera in an ego-centric manner
    void panLeft(float degrees);
    void panUp(float degrees);
    // orbit -- rotate camera around the look-at point
    void orbitLeft(float degrees);
    void orbitUp(float degrees);
    // dolly -- translate the camera ego-centrically
    void dollyLeft(float dist);
    void dollyUp(float dist);
    void dollyForward(float dist);
    // track -- move the camera to a specified position
    void trackTo(const Vec3f &pos);
    // zoom -- move the camera closer to or farther away from the lookAt point
    //          will have no effect if zoomLimit constraint is violated
    void zoomIn(float dist);

private:
    void recomputeYawPitch(); // from position and focus data
    void recomputeFocus(); // from other data
    void recomputePosition();
    
    void updateTransform();

private: // maintain alternate form of transformation...
    float yaw, pitch;
    Vec3f focus;
    Vec3f position;
};



